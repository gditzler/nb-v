module pipeline

import os
import math
import sync
import src.config
import src.model
import src.io as nbio
import src.kmer as kmod

// SeqJob carries one read (sequence ID, k-mer counts, and original input order index)
// queued for classification.
struct SeqJob {
	seq_id      string
	kmer_counts map[int]int
	index       int
}

// ClassifyResult holds the classification outcome for one read, including its
// original input index so results can be sorted back into input order.
struct ClassifyResult {
	seq_id     string
	best_class string
	best_score f64
	all_scores map[string]f64
	index      int
}

// BestResult records the highest-scoring class and its log-likelihood for one read.
// Used by classify_multiround to track the global best across multiple class chunks.
struct BestResult {
	best_class string
	best_score f64
}

// ClassPath holds the path and type of a class savefile for deferred loading.
struct ClassPath {
	path    string
	legacy  bool
}

// classify runs the full classification pipeline described by c: discovers class
// savefiles in c.save_dir, reads input sequences from c.source_dir, scores each
// read against all classes, and writes results to a file named by c.prefix and c.format.
// When c.limit_mb > 0 it uses multi-round classification to stay within the memory
// budget; otherwise all classes are loaded once and all reads are scored in a single pass.
pub fn classify(c config.Config) ! {
	// Verify kmer size matches training (only when meta.nbv exists, i.e. native savefiles)
	meta_path := '${c.save_dir}/meta.nbv'
	if os.exists(meta_path) {
		trained_k := nbio.load_meta(c.save_dir)!
		if trained_k != c.kmer_size {
			return error('kmer_size mismatch: config has ${c.kmer_size}, training used ${trained_k}')
		}
	}

	all_class_paths := discover_class_paths(c.save_dir, c.max_cols)!
	if all_class_paths.len == 0 {
		return error('no class savefiles found in ${c.save_dir}')
	}

	input_files := find_input_files(c.source_dir, c.extension)!

	limit_bytes := i64(c.limit_mb) * 1024 * 1024

	if limit_bytes <= 0 {
		// No memory limit -- load all classes and classify normally
		classes := load_all_classes(all_class_paths, c.kmer_size)!
		class_ids := classes.map(it.id)

		if c.threads > 1 {
			classify_multi(c, classes, class_ids, input_files)!
		} else {
			classify_single(c, classes, class_ids, input_files)!
		}
	} else {
		// Multi-round classification with memory limit
		classify_multiround(c, all_class_paths, input_files, limit_bytes)!
	}
}

// classify_single scores reads sequentially on the calling goroutine, writing
// results immediately as each read is processed.
fn classify_single(c config.Config, classes []model.NbClass, class_ids []string, input_files []string) ! {
	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!

	if c.full_result {
		writer.write_header(class_ids)!
	}

	mut rows_written := 0

	for input_file in input_files {
		if c.input_type == .fasta {
			records := nbio.read_fasta(input_file)!
			for record in records {
				if c.max_rows > 0 && rows_written >= c.max_rows {
					break
				}
				seq_id := record.header.split(' ')[0]
				kmer_counts := kmod.count_from_buffer(record.sequence, c.kmer_size)

				if kmer_counts.len == 0 {
					writer.write_no_valid_kmers(seq_id)!
					rows_written++
					continue
				}

				result := classify_read(seq_id, kmer_counts, classes, c.full_result, rows_written)
				write_classify_result(mut writer, result, class_ids, c.full_result)!
				rows_written++
			}
		} else {
			if c.max_rows > 0 && rows_written >= c.max_rows {
				break
			}
			kmer_counts := nbio.read_kmer_file(input_file, c.kmer_size)!
			seq_id := os.file_name(input_file).replace(c.extension, '')
			if kmer_counts.len == 0 {
				writer.write_no_valid_kmers(seq_id)!
				rows_written++
				continue
			}
			result := classify_read(seq_id, kmer_counts, classes, c.full_result, rows_written)
			write_classify_result(mut writer, result, class_ids, c.full_result)!
			rows_written++
		}
		if c.max_rows > 0 && rows_written >= c.max_rows {
			break
		}
	}

	writer.close()!
}

// classify_worker drains seq_ch, scores each read against all classes, and sends
// ClassifyResult values to result_ch. Signals completion to wg when seq_ch is closed.
fn classify_worker(seq_ch chan SeqJob, result_ch chan ClassifyResult, classes []model.NbClass, full_result bool, mut wg sync.WaitGroup) {
	defer {
		wg.done()
	}
	for {
		job := <-seq_ch or { break }
		result := classify_read(job.seq_id, job.kmer_counts, classes, full_result, job.index)
		result_ch <- result
	}
}

// classify_multi distributes reads across c.threads worker goroutines, collects
// results into a map keyed by input index, then writes output in original input order.
fn classify_multi(c config.Config, classes []model.NbClass, class_ids []string, input_files []string) ! {
	// Build all sequence jobs, tracking input order via index
	mut jobs := []SeqJob{}
	mut no_kmer_entries := map[int]string{}
	mut seq_index := 0

	for input_file in input_files {
		if c.input_type == .fasta {
			records := nbio.read_fasta(input_file)!
			for record in records {
				if c.max_rows > 0 && seq_index >= c.max_rows {
					break
				}
				seq_id := record.header.split(' ')[0]
				kmer_counts := kmod.count_from_buffer(record.sequence, c.kmer_size)

				if kmer_counts.len == 0 {
					no_kmer_entries[seq_index] = seq_id
					seq_index++
					continue
				}

				jobs << SeqJob{
					seq_id:      seq_id
					kmer_counts: kmer_counts
					index:       seq_index
				}
				seq_index++
			}
		} else {
			if c.max_rows > 0 && seq_index >= c.max_rows {
				break
			}
			kmer_counts := nbio.read_kmer_file(input_file, c.kmer_size)!
			seq_id := os.file_name(input_file).replace(c.extension, '')
			if kmer_counts.len == 0 {
				no_kmer_entries[seq_index] = seq_id
				seq_index++
				continue
			}
			jobs << SeqJob{
				seq_id:      seq_id
				kmer_counts: kmer_counts
				index:       seq_index
			}
			seq_index++
		}
		if c.max_rows > 0 && seq_index >= c.max_rows {
			break
		}
	}

	total_seqs := seq_index
	n_workers := if c.threads > jobs.len { jobs.len } else { c.threads }

	// Handle edge case: no valid jobs
	if n_workers == 0 {
		output_path := nbio.output_filename(c.prefix, c.format)
		mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!
		if c.full_result {
			writer.write_header(class_ids)!
		}
		for i in 0 .. total_seqs {
			if i in no_kmer_entries {
				writer.write_no_valid_kmers(no_kmer_entries[i])!
			}
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

	// Drain results into a map keyed by input index
	mut result_map := map[int]ClassifyResult{}
	for {
		result := <-result_ch or { break }
		result_map[result.index] = result
	}

	// Write output in original input order
	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!

	if c.full_result {
		writer.write_header(class_ids)!
	}

	for i in 0 .. total_seqs {
		if i in no_kmer_entries {
			writer.write_no_valid_kmers(no_kmer_entries[i])!
		} else if i in result_map {
			write_classify_result(mut writer, result_map[i], class_ids, c.full_result)!
		}
	}

	writer.close()!
}

// classify_multiround implements a memory-bounded classification strategy. All reads
// are loaded into memory up front (reads are small relative to class models), then the
// full set of class savefiles is split into chunks that each fit within limit_bytes.
// The reads are scored against each chunk in turn, keeping a running best (class, score)
// per read. This allows classifying against arbitrarily many classes with a fixed
// memory footprint determined by limit_mb in the config.
fn classify_multiround(c config.Config, all_class_paths []ClassPath, input_files []string, limit_bytes i64) ! {
	// Build all read jobs up front (reads are small relative to class models)
	mut jobs := []SeqJob{}
	mut no_kmer_ids := []string{}

	for input_file in input_files {
		if c.input_type == .fasta {
			records := nbio.read_fasta(input_file)!
			for record in records {
				if c.max_rows > 0 && (jobs.len + no_kmer_ids.len) >= c.max_rows {
					break
				}
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
			if c.max_rows > 0 && (jobs.len + no_kmer_ids.len) >= c.max_rows {
				break
			}
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
		if c.max_rows > 0 && (jobs.len + no_kmer_ids.len) >= c.max_rows {
			break
		}
	}

	// Split class paths into chunks that fit within the memory limit
	chunks := split_classes_by_memory(all_class_paths, c.kmer_size, limit_bytes)!

	// Track the best result for each read across all rounds
	mut best_results := map[string]BestResult{}

	for chunk in chunks {
		classes := load_all_classes(chunk, c.kmer_size)!

		// Score every read against this chunk of classes
		for job in jobs {
			for cls in classes {
				score := cls.compute_log_likelihood(job.kmer_counts)
				if job.seq_id in best_results {
					existing := best_results[job.seq_id]
					if score > existing.best_score {
						best_results[job.seq_id] = BestResult{
							best_class: cls.id
							best_score: score
						}
					}
				} else {
					best_results[job.seq_id] = BestResult{
						best_class: cls.id
						best_score: score
					}
				}
			}
		}
		// classes go out of scope here, freeing memory
	}

	// Write final output
	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!

	for sid in no_kmer_ids {
		writer.write_no_valid_kmers(sid)!
	}

	// Write results in the same order as jobs
	for job in jobs {
		if job.seq_id in best_results {
			br := best_results[job.seq_id]
			writer.write_result(job.seq_id, br.best_class, br.best_score)!
		}
	}

	writer.close()!
}

// classify_read scores one read against all classes and returns the best class and
// score. When full_result is true all per-class scores are also recorded.
fn classify_read(seq_id string, kmer_counts map[int]int, classes []model.NbClass, full_result bool, index int) ClassifyResult {
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
		index:      index
	}
}

// write_classify_result dispatches to write_full_result or write_result on the writer
// depending on whether full per-class scores were requested.
fn write_classify_result(mut writer nbio.Writer, result ClassifyResult, class_ids []string, full_result bool) ! {
	if full_result {
		writer.write_full_result(result.seq_id, result.all_scores, class_ids)!
	} else {
		writer.write_result(result.seq_id, result.best_class, result.best_score)!
	}
}

// discover_class_paths returns all class savefile paths found in save_dir, up to
// max_cols entries (0 means unlimited). Recognises both native .nbv files and
// legacy NBC++ -save.dat files; meta.nbv is excluded.
fn discover_class_paths(save_dir string, max_cols int) ![]ClassPath {
	mut paths := []ClassPath{}
	entries := os.ls(save_dir)!
	for entry in entries {
		if max_cols > 0 && paths.len >= max_cols {
			break
		}
		path := '${save_dir}/${entry}'
		if entry.ends_with('.nbv') && entry != 'meta.nbv' {
			paths << ClassPath{
				path:   path
				legacy: false
			}
		} else if entry.ends_with('-save.dat') {
			paths << ClassPath{
				path:   path
				legacy: true
			}
		}
	}
	return paths
}

// load_all_classes loads every ClassPath entry into an NbClass slice, using the
// native loader for .nbv files and the legacy NBC++ loader for -save.dat files.
fn load_all_classes(paths []ClassPath, kmer_size int) ![]model.NbClass {
	mut classes := []model.NbClass{}
	for cp in paths {
		if cp.legacy {
			cls := nbio.load_legacy_class(cp.path, kmer_size)!
			classes << cls
		} else {
			cls := nbio.load_class(cp.path)!
			classes << cls
		}
	}
	return classes
}

// split_classes_by_memory partitions paths into chunks whose estimated memory
// footprint stays within limit_bytes. Each class file's on-disk size is multiplied
// by 3 as a conservative estimate of its in-memory size (map overhead). A class that
// alone would exceed the limit still forms its own single-element chunk.
fn split_classes_by_memory(paths []ClassPath, kmer_size int, limit_bytes i64) ![][]ClassPath {
	mut chunks := [][]ClassPath{}
	mut current_chunk := []ClassPath{}
	mut current_bytes := i64(0)

	for cp in paths {
		// Estimate in-memory size from file size (loaded class is typically larger
		// than on-disk due to map overhead, so multiply by 3 as a conservative estimate)
		file_size := os.file_size(cp.path)
		estimated_mem := i64(file_size) * 3

		// If adding this class would exceed the limit and we already have classes
		// in the current chunk, start a new chunk
		if current_chunk.len > 0 && current_bytes + estimated_mem > limit_bytes {
			chunks << current_chunk
			current_chunk = []ClassPath{}
			current_bytes = 0
		}

		current_chunk << cp
		current_bytes += estimated_mem
	}

	if current_chunk.len > 0 {
		chunks << current_chunk
	}

	return chunks
}

// find_input_files returns all files in source_dir whose name ends with extension.
// Returns an error if no matching files are found.
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
