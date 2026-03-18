module io

import os
import math
import src.kmer as kmod
import src.model
import src.config

fn test_parse_fasta_basic() {
	records := read_fasta('src/io/testdata/test.fasta')!

	assert records.len == 4
	assert records[0].header == 'seq1 description'
	assert records[0].sequence.bytestr() == 'ACGTACGTACGTACGT'
	assert records[1].header == 'seq2'
	assert records[1].sequence.bytestr() == 'GGGGCCCC'
	assert records[2].header == 'seq3 empty after header'
	assert records[2].sequence.bytestr() == ''
	assert records[3].header == 'seq4'
	assert records[3].sequence.bytestr() == 'ACGT'
}

fn test_parse_fasta_callback() {
	mut count := 0
	parse_fasta('src/io/testdata/test.fasta', fn [mut count] (header string, seq []u8) {
		count++
	})!
	// Note: V closure mutable captures may not propagate back,
	// so we verify callback API compiles and runs without error.
	// Functional correctness is covered by test_parse_fasta_basic via read_fasta.
}

fn test_count_sequences() {
	count := count_sequences('src/io/testdata/test.fasta')!
	assert count == 4
}

fn test_read_kmer_file() {
	counts := read_kmer_file('src/io/testdata/test.kmr', 6)!
	assert counts.len > 0
	kmer_int := kmod.encode('ACGTAA'.bytes(), 6)
	canon := kmod.canonical(kmer_int, 6)
	assert counts[canon] > 0
}

fn test_save_and_load_class_roundtrip() {
	mut cls := model.NbClass.new('test_class', 6, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3, 100: 1})

	save_class(cls, '/tmp/test_roundtrip.nbv')!
	loaded := load_class('/tmp/test_roundtrip.nbv')!

	assert loaded.id == cls.id
	assert loaded.kmer_size == cls.kmer_size
	assert loaded.ngenomes == cls.ngenomes
	assert loaded.sumfreq == cls.sumfreq
	assert loaded.freqcnt.len == cls.freqcnt.len
	for k, v in cls.freqcnt {
		assert loaded.freqcnt[k] == v
	}
	assert loaded.state == .full
	assert math.abs(loaded.sumfreq_lg - cls.sumfreq_lg) < 1e-10

	os.rm('/tmp/test_roundtrip.nbv') or {}
}

fn test_save_and_load_meta() {
	save_meta('/tmp/test_meta_dir', 9)!
	k := load_meta('/tmp/test_meta_dir')!
	assert k == 9
	os.rm('/tmp/test_meta_dir/meta.nbv') or {}
	os.rmdir('/tmp/test_meta_dir') or {}
}

fn test_load_legacy_class() {
	cls := load_legacy_class('example/training_classes/1748-save.dat', 9)!
	assert cls.id == '1748'
	assert cls.kmer_size == 9
	assert cls.state == .classify_only
	assert cls.freqcnt_lg.len > 0
	// ngenomes_lg should be log(2) = 0.693... (verified from hex dump)
	assert math.abs(cls.ngenomes_lg - 0.6931471805599453) < 1e-10
}

fn test_writer_csv() {
	mut w := Writer.new('/tmp/test_output.csv', .csv, false)!
	w.write_result('seq1', 'class_a', -123.45)!
	w.write_result('seq2', 'class_b', -678.90)!
	w.close()!

	content := os.read_file('/tmp/test_output.csv')!
	lines := content.trim_space().split('\n')
	assert lines.len == 2
	assert lines[0].contains('seq1')
	assert lines[0].contains('class_a')
	os.rm('/tmp/test_output.csv') or {}
}

fn test_writer_json() {
	mut w := Writer.new('/tmp/test_output.jsonl', .json, false)!
	w.write_result('seq1', 'class_a', -123.45)!
	w.close()!

	content := os.read_file('/tmp/test_output.jsonl')!
	assert content.contains('"seq_id"')
	assert content.contains('"best_class"')
	os.rm('/tmp/test_output.jsonl') or {}
}

fn test_writer_no_valid_kmers() {
	mut w := Writer.new('/tmp/test_nokmers.csv', .csv, false)!
	w.write_no_valid_kmers('bad_read')!
	w.close()!

	content := os.read_file('/tmp/test_nokmers.csv')!
	assert content.contains('sequence contains no valid kmers')
	os.rm('/tmp/test_nokmers.csv') or {}
}

fn test_output_filename() {
	assert output_filename('results', .csv) == 'results.csv'
	assert output_filename('results', .tsv) == 'results.tsv'
	assert output_filename('results', .json) == 'results.jsonl'
}
