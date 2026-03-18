module kmer

fn base_to_int(b u8) int {
	return match b {
		`A`, `a` { 0 }
		`C`, `c` { 1 }
		`G`, `g` { 2 }
		`T`, `t` { 3 }
		else { -1 }
	}
}

fn complement(val int) int {
	return 3 - val
}

pub fn encode(kmer []u8, k int) int {
	mut result := 0
	for i in 0 .. k {
		result = result * 4 + base_to_int(kmer[i])
	}
	return result
}

pub fn reverse_complement(kmer_int int, k int) int {
	mut result := 0
	mut val := kmer_int
	for _ in 0 .. k {
		result = result * 4 + complement(val & 3)
		val >>= 2
	}
	return result
}

pub fn canonical(kmer_int int, k int) int {
	rc := reverse_complement(kmer_int, k)
	if kmer_int <= rc {
		return kmer_int
	}
	return rc
}

pub fn count_from_buffer(buf []u8, k int) map[int]int {
	mut counts := map[int]int{}
	mut window := 0
	mut valid_len := 0

	for i in 0 .. buf.len {
		val := base_to_int(buf[i])
		if val < 0 {
			valid_len = 0
			window = 0
			continue
		}
		window = (window * 4 + val) & ((1 << (2 * k)) - 1)
		valid_len++
		if valid_len >= k {
			canon := canonical(window, k)
			counts[canon] = counts[canon] + 1
		}
	}
	return counts
}

pub fn num_canonical_kmers(k int) i64 {
	total := i64(u64(1) << (2 * k))
	if k % 2 == 1 {
		return total / 2
	}
	palindromes := i64(u64(1) << k)
	return (total + palindromes) / 2
}
