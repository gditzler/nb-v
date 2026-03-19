// pipeline provides the train and classify orchestrators that drive the full
// NBV workflow. Both pipelines support single-threaded and multithreaded execution
// paths selected by config.Config.threads. This file covers training; classify.v
// covers classification. The two files share the pipeline module.
module pipeline

import os
import sync
import src.config
import src.model
import src.io as nbio
import src.kmer as kmod

// TrainJob carries one (class, file) pair queued for k-mer counting.
struct TrainJob {
	class_id string
	path     string
}

// TrainResult carries the k-mer counts produced by one training file.
struct TrainResult {
	class_id    string
	kmer_counts map[int]int
}

// scan_training_dir enumerates source_dir for per-class subdirectories and
// returns one TrainJob per file matching extension found inside each subdirectory.
fn scan_training_dir(source_dir string, extension string) ![]TrainJob {
	mut jobs := []TrainJob{}
	entries := os.ls(source_dir)!
	for entry in entries {
		subdir := '${source_dir}/${entry}'
		if !os.is_dir(subdir) {
			continue
		}
		class_id := entry
		files := os.ls(subdir)!
		for file in files {
			if file.ends_with(extension) {
				jobs << TrainJob{
					class_id: class_id
					path:     '${subdir}/${file}'
				}
			}
		}
	}
	return jobs
}

// load_kmer_counts reads k-mer counts from a single file, dispatching to the
// FASTA reader or the .kmr reader depending on input_type.
fn load_kmer_counts(path string, input_type config.InputType, k int) !map[int]int {
	if input_type == .fasta {
		records := nbio.read_fasta(path)!
		mut merged := map[int]int{}
		for record in records {
			counts := kmod.count_from_buffer(record.sequence, k)
			for km, count in counts {
				merged[km] = merged[km] + count
			}
		}
		return merged
	} else {
		return nbio.read_kmer_file(path, k)!
	}
}

// train_worker drains job_ch, computes k-mer counts for each file, and sends
// TrainResult values to result_ch. Signals completion to wg when job_ch is closed.
fn train_worker(job_ch chan TrainJob, result_ch chan TrainResult, input_type config.InputType, k int, mut wg sync.WaitGroup) {
	defer {
		wg.done()
	}
	for {
		job := <-job_ch or { break }
		kmer_counts := load_kmer_counts(job.path, input_type, k) or {
			eprintln('Warning: skipping ${job.path}: ${err}')
			continue
		}
		result_ch <- TrainResult{
			class_id:    job.class_id
			kmer_counts: kmer_counts
		}
	}
}

// accumulate_results folds a batch of TrainResult values into the classes map,
// creating a new NbClass for any class_id seen for the first time.
fn accumulate_results(mut classes map[string]model.NbClass, results []TrainResult, c config.Config) {
	for result in results {
		if result.class_id !in classes {
			savefile := '${c.save_dir}/${result.class_id}.nbv'
			classes[result.class_id] = model.NbClass.new(result.class_id, c.kmer_size, savefile)
		}
		mut cls := classes[result.class_id]
		cls.add_genome(result.kmer_counts)
		classes[result.class_id] = cls
	}
}

// train_single processes jobs sequentially on the calling goroutine.
// When batch_size > 0 it checkpoints class models to disk after every batch_size
// files, resuming from the saved state on the next run.
fn train_single(c config.Config, jobs []TrainJob) !map[string]model.NbClass {
	mut classes := map[string]model.NbClass{}
	mut processed := 0

	for job in jobs {
		if job.class_id !in classes {
			savefile := '${c.save_dir}/${job.class_id}.nbv'
			if c.batch_size > 0 && os.exists(savefile) {
				classes[job.class_id] = nbio.load_class(savefile)!
			} else {
				classes[job.class_id] = model.NbClass.new(job.class_id, c.kmer_size, savefile)
			}
		}

		kmer_counts := load_kmer_counts(job.path, c.input_type, c.kmer_size) or {
			eprintln('Warning: skipping ${job.path}: ${err}')
			continue
		}

		mut cls := classes[job.class_id]
		cls.add_genome(kmer_counts)
		classes[job.class_id] = cls
		processed++

		if c.batch_size > 0 && processed % c.batch_size == 0 {
			os.mkdir_all(c.save_dir)!
			for _, cl in classes {
				nbio.save_class(cl, cl.savefile)!
			}
		}
	}

	return classes
}

// train_multi distributes jobs across c.threads worker goroutines and accumulates
// results on the main thread. Delegates to train_multi_batched when batch_size > 0.
fn train_multi(c config.Config, jobs []TrainJob) !map[string]model.NbClass {
	if c.batch_size > 0 {
		return train_multi_batched(c, jobs)
	}

	n_workers := if c.threads > jobs.len { jobs.len } else { c.threads }

	job_ch := chan TrainJob{cap: jobs.len}
	result_ch := chan TrainResult{cap: jobs.len}

	mut wg := sync.new_waitgroup()
	wg.add(n_workers)

	// Spawn worker threads
	for _ in 0 .. n_workers {
		spawn train_worker(job_ch, result_ch, c.input_type, c.kmer_size, mut wg)
	}

	// Feed all jobs into the channel, then close it
	for job in jobs {
		job_ch <- job
	}
	job_ch.close()

	// Wait for all workers to finish, then close result channel
	wg.wait()
	result_ch.close()

	// Drain results from the channel
	mut results := []TrainResult{}
	for {
		result := <-result_ch or { break }
		results << result
	}

	// Accumulate into classes in the main thread
	mut classes := map[string]model.NbClass{}
	accumulate_results(mut classes, results, c)

	return classes
}

// train_multi_batched processes jobs in fixed-size chunks using c.threads workers,
// flushing all class models to disk after each chunk so training can be resumed if
// interrupted. Existing savefiles are loaded at the start of each chunk.
fn train_multi_batched(c config.Config, jobs []TrainJob) !map[string]model.NbClass {
	mut classes := map[string]model.NbClass{}
	mut batch_start := 0

	for batch_start < jobs.len {
		batch_end := if batch_start + c.batch_size > jobs.len {
			jobs.len
		} else {
			batch_start + c.batch_size
		}
		chunk := jobs[batch_start..batch_end]

		n_workers := if c.threads > chunk.len { chunk.len } else { c.threads }

		job_ch := chan TrainJob{cap: chunk.len}
		result_ch := chan TrainResult{cap: chunk.len}

		mut wg := sync.new_waitgroup()
		wg.add(n_workers)

		for _ in 0 .. n_workers {
			spawn train_worker(job_ch, result_ch, c.input_type, c.kmer_size, mut wg)
		}

		for job in chunk {
			job_ch <- job
		}
		job_ch.close()

		wg.wait()
		result_ch.close()

		mut results := []TrainResult{}
		for {
			result := <-result_ch or { break }
			results << result
		}

		// Load existing classes from disk if not yet in memory
		for result in results {
			if result.class_id !in classes {
				savefile := '${c.save_dir}/${result.class_id}.nbv'
				if os.exists(savefile) {
					classes[result.class_id] = nbio.load_class(savefile)!
				} else {
					classes[result.class_id] = model.NbClass.new(result.class_id, c.kmer_size, savefile)
				}
			}
		}
		accumulate_results(mut classes, results, c)

		// Flush batch to disk
		os.mkdir_all(c.save_dir)!
		for _, cl in classes {
			nbio.save_class(cl, cl.savefile)!
		}

		batch_start = batch_end
	}

	return classes
}

// train runs the full training pipeline described by c: scans the source directory
// for per-class training files, builds one NbClass model per class, and writes
// all class savefiles plus a meta.nbv index to c.save_dir.
pub fn train(c config.Config) ! {
	jobs := scan_training_dir(c.source_dir, c.extension)!
	if jobs.len == 0 {
		return error('no training files found in ${c.source_dir}')
	}

	classes := if c.threads > 1 {
		train_multi(c, jobs)!
	} else {
		train_single(c, jobs)!
	}

	os.mkdir_all(c.save_dir)!
	for _, cls in classes {
		nbio.save_class(cls, cls.savefile)!
	}
	nbio.save_meta(c.save_dir, c.kmer_size)!
}
