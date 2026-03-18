module pipeline

import os
import src.config

fn test_train_creates_savefiles() {
	out_dir := '/tmp/nbv_test_train'
	os.rmdir_all(out_dir) or {}

	c := config.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   out_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    1
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 0
	}

	train(c)!

	assert os.exists('${out_dir}/class_a.nbv')
	assert os.exists('${out_dir}/class_b.nbv')
	assert os.exists('${out_dir}/meta.nbv')

	meta_content := os.read_file('${out_dir}/meta.nbv')!
	assert meta_content.trim_space() == '4'

	os.rmdir_all(out_dir) or {}
}
