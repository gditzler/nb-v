module io

import os
import src.kmer as kmod

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
