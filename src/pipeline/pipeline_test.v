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

fn test_train_then_classify() {
	train_dir := '/tmp/nbv_test_e2e_train'
	os.rmdir_all(train_dir) or {}

	train_cfg := config.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   train_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    1
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 0
	}
	train(train_cfg)!

	// Create test read from class_a's genome
	os.mkdir_all('/tmp/nbv_test_e2e_reads') or {}
	os.write_file('/tmp/nbv_test_e2e_reads/test_read.fasta', '>read1\nACGTACGTACGTACGT\n')!

	classify_cfg := config.Config{
		mode:        .classify
		kmer_size:   4
		save_dir:    train_dir
		source_dir:  '/tmp/nbv_test_e2e_reads'
		threads:     1
		input_type:  .fasta
		extension:   '.fasta'
		format:      .csv
		prefix:      '/tmp/nbv_test_e2e_output'
		full_result: false
		temp_dir:    '/tmp'
		limit_mb:    0
		max_rows:    1000
		max_cols:    100
	}
	classify(classify_cfg)!

	output := os.read_file('/tmp/nbv_test_e2e_output.csv')!
	assert output.contains('read1')

	os.rmdir_all(train_dir) or {}
	os.rmdir_all('/tmp/nbv_test_e2e_reads') or {}
	os.rm('/tmp/nbv_test_e2e_output.csv') or {}
}

fn test_train_multithreaded() {
	out_dir := '/tmp/nbv_test_train_mt'
	os.rmdir_all(out_dir) or {}

	c := config.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   out_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    2
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 0
	}

	train(c)!

	assert os.exists('${out_dir}/class_a.nbv')
	assert os.exists('${out_dir}/class_b.nbv')

	os.rmdir_all(out_dir) or {}
}
