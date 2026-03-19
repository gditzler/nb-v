// model provides the per-class Naive Bayes model used for metagenomic sequence
// classification. It implements k-mer frequency counting, Laplace-smoothed
// log-likelihood scoring, and numerically stable accumulation via Kahan summation.
module model

import math
import src.kmer as kmod

// KahanAccumulator holds the running sum and compensation term for Kahan
// compensated summation, reducing floating-point error when accumulating many
// log-probability values.
pub struct KahanAccumulator {
pub mut:
	sum  f64
	comp f64
}

// kahan_add returns a new KahanAccumulator that incorporates val into acc using
// Kahan compensated summation to reduce floating-point rounding error.
pub fn kahan_add(acc KahanAccumulator, val f64) KahanAccumulator {
	y := val - acc.comp
	t := acc.sum + y
	return KahanAccumulator{
		sum:  t
		comp: (t - acc.sum) - y
	}
}

// LoadState represents the loading status of an NbClass model, indicating
// whether parameter data has been read from disk and at what level of detail.
pub enum LoadState {
	// unloaded indicates no model data has been loaded from disk.
	unloaded
	// full indicates the model was trained or loaded with all parameters,
	// including raw counts and log-space values.
	full
	// classify_only indicates the model was loaded with only the log-space
	// parameters needed for classification, without raw counts.
	classify_only
}

// NbClass represents a single Naive Bayes class model for one taxonomic group.
// It accumulates k-mer counts across training genomes and maintains both raw
// counts and precomputed log-space values for efficient classification scoring.
// Laplace smoothing is applied via sumfreq, which is initialized to the number
// of canonical k-mers for the given k.
pub struct NbClass {
pub mut:
	id          string
	kmer_size   int
	savefile    string
	ngenomes_lg f64
	sumfreq_lg  f64
	freqcnt_lg  map[int]f64
	ngenomes    int
	sumfreq     i64
	freqcnt     map[int]int
	state       LoadState
}

// NbClass.new returns a new, unloaded NbClass for the given taxonomic id,
// k-mer size, and save file path. sumfreq is initialized to the number of
// canonical k-mers for kmer_size to provide Laplace smoothing from the start.
pub fn NbClass.new(id string, kmer_size int, savefile string) NbClass {
	v := kmod.num_canonical_kmers(kmer_size)
	return NbClass{
		id:        id
		kmer_size: kmer_size
		savefile:  savefile
		sumfreq:   v
		state:     .unloaded
	}
}

// add_genome incorporates the k-mer counts from one genome into the model,
// updating raw counts, the genome counter, and all precomputed log-space
// parameters. Log-space values are recomputed eagerly after each call.
// Sets state to .full.
pub fn (mut self NbClass) add_genome(kmer_counts map[int]int) {
	self.ngenomes += 1
	mut total := i64(0)
	for km, count in kmer_counts {
		self.freqcnt[km] = self.freqcnt[km] + count
		total += count
	}
	self.sumfreq += total

	self.ngenomes_lg = math.log(f64(self.ngenomes))
	self.sumfreq_lg = math.log(f64(self.sumfreq))
	for km, count in self.freqcnt {
		self.freqcnt_lg[km] = math.log(f64(count + 1))
	}

	self.state = .full
}

// get_freq_count_lg returns the precomputed log-space frequency count for the
// given canonical k-mer index. Returns 0.0 for k-mers not present in the model.
pub fn (self &NbClass) get_freq_count_lg(km int) f64 {
	return self.freqcnt_lg[km]
}

// compute_log_likelihood returns the Naive Bayes log-likelihood score for a
// read represented by kmer_counts against this class model. Uses Kahan
// compensated summation over log-space frequency counts, normalized by
// sumfreq_lg. The class prior is not included, matching NBC++ behavior.
pub fn (self &NbClass) compute_log_likelihood(kmer_counts map[int]int) f64 {
	mut acc := KahanAccumulator{}
	mut total_count := i64(0)

	for km, freq in kmer_counts {
		fcl := self.get_freq_count_lg(km)
		acc = kahan_add(acc, f64(freq) * fcl)
		total_count += freq
	}

	return acc.sum - f64(total_count) * self.sumfreq_lg
}

// size_bytes returns an estimate of the memory used by this NbClass instance
// in bytes, including the base struct and both the raw and log-space k-mer maps.
pub fn (self &NbClass) size_bytes() u64 {
	base := u64(sizeof(NbClass))
	freq_plain := u64(self.freqcnt.len) * u64(sizeof(int) + sizeof(int))
	freq_log := u64(self.freqcnt_lg.len) * u64(sizeof(int) + sizeof(f64))
	return base + freq_plain + freq_log
}
