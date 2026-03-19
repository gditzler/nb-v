module io

import os
import src.kmer as kmod

// read_kmer_file reads a tab-separated k-mer count file at path and returns a
// map from canonical k-mer integer encoding to total count. Each line must
// contain a k-mer string and an integer count separated by a tab. K-mers are
// encoded and collapsed to their canonical (lexicographically smaller) form
// before accumulation. Returns an error if the file cannot be read, a line is
// malformed, or a k-mer's length does not match k.
pub fn read_kmer_file(path string, k int) !map[int]int {
	lines := os.read_lines(path)!
	mut counts := map[int]int{}

	for i, line in lines {
		trimmed := line.trim_space()
		if trimmed.len == 0 {
			continue
		}
		parts := trimmed.split('\t')
		if parts.len != 2 {
			return error('malformed line ${i + 1} in ${path}: expected <kmer>\\t<count>')
		}
		kmer_str := parts[0]
		count := parts[1].int()
		if kmer_str.len != k {
			return error('kmer length mismatch on line ${i + 1}: expected ${k}, got ${kmer_str.len}')
		}
		kmer_int := kmod.encode(kmer_str.bytes(), k)
		canon := kmod.canonical(kmer_int, k)
		counts[canon] = counts[canon] + count
	}

	return counts
}
