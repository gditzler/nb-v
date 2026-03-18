module kmer

fn test_encode_single_bases() {
	assert encode('A'.bytes(), 1) == 0
	assert encode('C'.bytes(), 1) == 1
	assert encode('G'.bytes(), 1) == 2
	assert encode('T'.bytes(), 1) == 3
}

fn test_encode_kmer() {
	assert encode('AC'.bytes(), 2) == 1
	assert encode('GT'.bytes(), 2) == 11
	assert encode('ACG'.bytes(), 3) == 6
}

fn test_reverse_complement() {
	assert reverse_complement(0, 1) == 3
	assert reverse_complement(1, 1) == 2
	assert reverse_complement(1, 2) == 11
	assert reverse_complement(11, 2) == 1
}

fn test_canonical() {
	assert canonical(0, 1) == 0
	assert canonical(3, 1) == 0
	assert canonical(1, 2) == 1
	assert canonical(11, 2) == 1
}

fn test_count_from_buffer() {
	counts := count_from_buffer('ACGT'.bytes(), 2)
	assert counts[canonical(encode('AC'.bytes(), 2), 2)] == 2
	assert counts[canonical(encode('CG'.bytes(), 2), 2)] == 1
}

fn test_count_from_buffer_skips_invalid() {
	counts := count_from_buffer('ACNGT'.bytes(), 2)
	assert counts[1] == 2
}

fn test_count_from_buffer_skips_newlines() {
	counts := count_from_buffer('AC\nGT'.bytes(), 2)
	assert counts[1] == 2
}

fn test_num_canonical_kmers() {
	assert num_canonical_kmers(1) == 2
	assert num_canonical_kmers(2) == 10
	assert num_canonical_kmers(3) == 32
	assert num_canonical_kmers(6) == 2080
}
