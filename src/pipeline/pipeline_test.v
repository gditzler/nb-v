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

fn test_classify_multithreaded() {
	train_dir := '/tmp/nbv_test_e2e_mt_train'
	os.rmdir_all(train_dir) or {}

	train(config.Config{
		mode: .train, kmer_size: 4, save_dir: train_dir,
		source_dir: 'src/pipeline/testdata/training',
		threads: 1, input_type: .fasta, extension: '.fasta',
	})!

	os.mkdir_all('/tmp/nbv_test_e2e_mt_reads') or {}
	os.write_file('/tmp/nbv_test_e2e_mt_reads/test.fasta', '>r1\nACGTACGTACGTACGT\n>r2\nGGGGCCCCAAAATTTT\n')!

	classify(config.Config{
		mode: .classify, kmer_size: 4, save_dir: train_dir,
		source_dir: '/tmp/nbv_test_e2e_mt_reads',
		threads: 2, input_type: .fasta, extension: '.fasta',
		format: .csv, prefix: '/tmp/nbv_test_e2e_mt_out',
		full_result: false, temp_dir: '/tmp',
	})!

	output := os.read_file('/tmp/nbv_test_e2e_mt_out.csv')!
	assert output.contains('r1')
	assert output.contains('r2')

	os.rmdir_all(train_dir) or {}
	os.rmdir_all('/tmp/nbv_test_e2e_mt_reads') or {}
	os.rm('/tmp/nbv_test_e2e_mt_out.csv') or {}
}

fn test_classify_with_legacy_savefiles() {
	c := config.Config{
		mode:        .classify
		kmer_size:   9
		save_dir:    'example/training_classes'
		source_dir:  'example/reads'
		threads:     1
		input_type:  .fasta
		extension:   '.fna'
		format:      .csv
		prefix:      '/tmp/nbv_legacy_test_output'
		full_result: false
		temp_dir:    '/tmp'
		limit_mb:    0
		max_rows:    450000
		max_cols:    20000
	}

	classify(c)!

	output := os.read_file('/tmp/nbv_legacy_test_output.csv')!
	expected := os.read_file('example/results_max_1.csv')!

	output_lines := output.trim_space().split('\n')
	expected_lines := expected.trim_space().split('\n')

	// Check that class assignments match for as many reads as possible
	mut matches := 0
	mut total := 0
	for i, exp_line in expected_lines {
		if i >= output_lines.len {
			break
		}
		exp_parts := exp_line.split(',')
		out_parts := output_lines[i].split(',')
		if exp_parts.len >= 2 && out_parts.len >= 2 {
			total++
			if exp_parts[0] == out_parts[0] && exp_parts[1] == out_parts[1] {
				matches++
			}
		}
	}

	// At least 95% of class assignments should match
	assert total > 0
	match_pct := f64(matches) / f64(total)
	assert match_pct > 0.95

	os.rm('/tmp/nbv_legacy_test_output.csv') or {}
}
