module pipeline

import os
import math
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
	// Verify kmer size matches training
	trained_k := nbio.load_meta(c.save_dir)!
	if trained_k != c.kmer_size {
		return error('kmer_size mismatch: config has ${c.kmer_size}, training used ${trained_k}')
	}

	classes := load_classes(c.save_dir, c.kmer_size)!
	if classes.len == 0 {
		return error('no class savefiles found in ${c.save_dir}')
	}

	class_ids := classes.map(it.id)

	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!

	if c.full_result {
		writer.write_header(class_ids)!
	}

	input_files := find_input_files(c.source_dir, c.extension)!

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
