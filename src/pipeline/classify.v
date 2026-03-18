module pipeline

import os
import math
import sync
import src.config
import src.model
import src.io as nbio
import src.kmer as kmod

struct SeqJob {
	seq_id      string
	kmer_counts map[int]int
}

struct ClassifyResult {
	seq_id     string
	best_class string
	best_score f64
	all_scores map[string]f64
}

pub fn classify(c config.Config) ! {
	// Verify kmer size matches training (only when meta.nbv exists, i.e. native savefiles)
	meta_path := '${c.save_dir}/meta.nbv'
	if os.exists(meta_path) {
		trained_k := nbio.load_meta(c.save_dir)!
		if trained_k != c.kmer_size {
			return error('kmer_size mismatch: config has ${c.kmer_size}, training used ${trained_k}')
		}
	}

	classes := load_classes(c.save_dir, c.kmer_size)!
	if classes.len == 0 {
		return error('no class savefiles found in ${c.save_dir}')
	}

	class_ids := classes.map(it.id)

	input_files := find_input_files(c.source_dir, c.extension)!

	if c.threads > 1 {
		classify_multi(c, classes, class_ids, input_files)!
	} else {
		classify_single(c, classes, class_ids, input_files)!
	}
}

fn classify_single(c config.Config, classes []model.NbClass, class_ids []string, input_files []string) ! {
	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!

	if c.full_result {
		writer.write_header(class_ids)!
	}

	for input_file in input_files {
		if c.input_type == .fasta {
			records := nbio.read_fasta(input_file)!
			for record in records {
				seq_id := record.header.split(' ')[0]
				kmer_counts := kmod.count_from_buffer(record.sequence, c.kmer_size)

				if kmer_counts.len == 0 {
					writer.write_no_valid_kmers(seq_id)!
					continue
				}

				result := classify_read(seq_id, kmer_counts, classes, c.full_result)
				write_classify_result(mut writer, result, class_ids, c.full_result)!
			}
		} else {
			kmer_counts := nbio.read_kmer_file(input_file, c.kmer_size)!
			seq_id := os.file_name(input_file).replace(c.extension, '')
			if kmer_counts.len == 0 {
				writer.write_no_valid_kmers(seq_id)!
				continue
			}
			result := classify_read(seq_id, kmer_counts, classes, c.full_result)
			write_classify_result(mut writer, result, class_ids, c.full_result)!
		}
	}

	writer.close()!
}

fn classify_worker(seq_ch chan SeqJob, result_ch chan ClassifyResult, classes []model.NbClass, full_result bool, mut wg sync.WaitGroup) {
	defer {
		wg.done()
	}
	for {
		job := <-seq_ch or { break }
		result := classify_read(job.seq_id, job.kmer_counts, classes, full_result)
		result_ch <- result
	}
}

fn classify_multi(c config.Config, classes []model.NbClass, class_ids []string, input_files []string) ! {
	// First, build all sequence jobs so we know total count for channel capacity
	mut jobs := []SeqJob{}
	mut no_kmer_ids := []string{}

	for input_file in input_files {
		if c.input_type == .fasta {
			records := nbio.read_fasta(input_file)!
			for record in records {
				seq_id := record.header.split(' ')[0]
				kmer_counts := kmod.count_from_buffer(record.sequence, c.kmer_size)

				if kmer_counts.len == 0 {
					no_kmer_ids << seq_id
					continue
				}

				jobs << SeqJob{
					seq_id:      seq_id
					kmer_counts: kmer_counts
				}
			}
		} else {
			kmer_counts := nbio.read_kmer_file(input_file, c.kmer_size)!
			seq_id := os.file_name(input_file).replace(c.extension, '')
			if kmer_counts.len == 0 {
				no_kmer_ids << seq_id
				continue
			}
			jobs << SeqJob{
				seq_id:      seq_id
				kmer_counts: kmer_counts
			}
		}
	}

	n_workers := if c.threads > jobs.len { jobs.len } else { c.threads }

	// Handle edge case: no valid jobs
	if n_workers == 0 {
		output_path := nbio.output_filename(c.prefix, c.format)
		mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!
		if c.full_result {
			writer.write_header(class_ids)!
		}
		for sid in no_kmer_ids {
			writer.write_no_valid_kmers(sid)!
		}
		writer.close()!
		return
	}

	seq_ch := chan SeqJob{cap: jobs.len}
	result_ch := chan ClassifyResult{cap: jobs.len}

	mut wg := sync.new_waitgroup()
	wg.add(n_workers)

	// Spawn worker threads
	for _ in 0 .. n_workers {
		spawn classify_worker(seq_ch, result_ch, classes, c.full_result, mut wg)
	}

	// Feed all jobs into the channel, then close it
	for job in jobs {
		seq_ch <- job
	}
	seq_ch.close()

	// Wait for all workers to finish, then close result channel
	wg.wait()
	result_ch.close()

	// Drain results
	mut results := []ClassifyResult{}
	for {
		result := <-result_ch or { break }
		results << result
	}

	// Write output
	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!

	if c.full_result {
		writer.write_header(class_ids)!
	}

	// Write sequences with no valid kmers first
	for sid in no_kmer_ids {
		writer.write_no_valid_kmers(sid)!
	}

	// Write classification results
	for result in results {
		write_classify_result(mut writer, result, class_ids, c.full_result)!
	}

	writer.close()!
}

fn classify_read(seq_id string, kmer_counts map[int]int, classes []model.NbClass, full_result bool) ClassifyResult {
	mut best_class := ''
	mut best_score := -math.max_f64
	mut all_scores := map[string]f64{}

	for cls in classes {
		score := cls.compute_log_likelihood(kmer_counts)
		if full_result {
			all_scores[cls.id] = score
		}
		if score > best_score {
			best_score = score
			best_class = cls.id
		}
	}

	return ClassifyResult{
		seq_id:     seq_id
		best_class: best_class
		best_score: best_score
		all_scores: all_scores
	}
}

fn write_classify_result(mut writer nbio.Writer, result ClassifyResult, class_ids []string, full_result bool) ! {
	if full_result {
		writer.write_full_result(result.seq_id, result.all_scores, class_ids)!
	} else {
		writer.write_result(result.seq_id, result.best_class, result.best_score)!
	}
}

fn load_classes(save_dir string, kmer_size int) ![]model.NbClass {
	mut classes := []model.NbClass{}
	entries := os.ls(save_dir)!
	for entry in entries {
		path := '${save_dir}/${entry}'
		if entry.ends_with('.nbv') && entry != 'meta.nbv' {
			cls := nbio.load_class(path)!
			classes << cls
		} else if entry.ends_with('-save.dat') {
			cls := nbio.load_legacy_class(path, kmer_size)!
			classes << cls
		}
	}
	return classes
}

fn find_input_files(source_dir string, extension string) ![]string {
	mut files := []string{}
	entries := os.ls(source_dir)!
	for entry in entries {
		if entry.ends_with(extension) {
			files << '${source_dir}/${entry}'
		}
	}
	if files.len == 0 {
		return error('no input files with extension ${extension} found in ${source_dir}')
	}
	return files
}
