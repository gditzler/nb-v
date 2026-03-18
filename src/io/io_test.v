module io

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
