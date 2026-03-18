module pipeline

import os
import src.config
import src.model
import src.io as nbio
import src.kmer as kmod

struct TrainJob {
	class_id string
	path     string
}

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

pub fn train(c config.Config) ! {
	jobs := scan_training_dir(c.source_dir, c.extension)!
	if jobs.len == 0 {
		return error('no training files found in ${c.source_dir}')
	}

	mut classes := map[string]model.NbClass{}

	for job in jobs {
		if job.class_id !in classes {
			savefile := '${c.save_dir}/${job.class_id}.nbv'
			classes[job.class_id] = model.NbClass.new(job.class_id, c.kmer_size, savefile)
		}

		kmer_counts := load_kmer_counts(job.path, c.input_type, c.kmer_size) or {
			eprintln('Warning: skipping ${job.path}: ${err}')
			continue
		}

		mut cls := classes[job.class_id]
		cls.add_genome(kmer_counts)
		classes[job.class_id] = cls
	}

	os.mkdir_all(c.save_dir)!
	for _, cls in classes {
		nbio.save_class(cls, cls.savefile)!
	}
	nbio.save_meta(c.save_dir, c.kmer_size)!
}
