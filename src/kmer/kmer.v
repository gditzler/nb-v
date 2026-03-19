module kmer

// kmer provides pure functions for DNA k-mer encoding, reverse complementation,
// canonicalization, and frequency counting from byte buffers. All functions are
// thread-safe and operate without shared mutable state. Supports k-mer sizes up
// to k=15 using base-4 integer encoding in a 32-bit int.

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

// encode converts a DNA k-mer byte slice into its base-4 integer representation.
// Each base is mapped as A=0, C=1, G=2, T=3. The slice must contain at least k
// bytes and each byte must be a valid ACGT character (upper or lower case).
pub fn encode(kmer []u8, k int) int {
	mut result := 0
	for i in 0 .. k {
		result = result * 4 + base_to_int(kmer[i])
	}
	return result
}

// reverse_complement returns the base-4 integer encoding of the reverse complement
// of a k-mer given its integer encoding and length k. Complement values are computed
// as 3-val (A<->T, C<->G) and the bases are reversed by iterating from LSB to MSB.
pub fn reverse_complement(kmer_int int, k int) int {
	mut result := 0
	mut val := kmer_int
	for _ in 0 .. k {
		result = result * 4 + complement(val & 3)
		val >>= 2
	}
	return result
}

// canonical returns the canonical form of a k-mer, defined as the lexicographically
// smaller of the k-mer and its reverse complement. Using canonical k-mers reduces
// the feature space by treating a strand and its complement as the same feature.
pub fn canonical(kmer_int int, k int) int {
	rc := reverse_complement(kmer_int, k)
	if kmer_int <= rc {
		return kmer_int
	}
	return rc
}

// count_from_buffer counts canonical k-mer occurrences in a raw byte buffer using
// a sliding window. Non-ACGT characters reset the window, so k-mers spanning
// ambiguous bases or sequence boundaries are excluded. Returns a map from canonical
// k-mer integer encoding to observation count.
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

// num_canonical_kmers returns the number of distinct canonical k-mers for a given k.
// For odd k, no palindromic k-mers exist, so the count is exactly 4^k / 2. For even
// k, palindromic k-mers (which are their own reverse complement) are counted once,
// yielding (4^k + 2^k) / 2. This value is used to initialize Laplace smoothing.
pub fn num_canonical_kmers(k int) i64 {
	total := i64(u64(1) << (2 * k))
	if k % 2 == 1 {
		return total / 2
	}
	palindromes := i64(u64(1) << k)
	return (total + palindromes) / 2
}
