module model

import math
import src.kmer as kmod

pub struct KahanAccumulator {
pub mut:
	sum  f64
	comp f64
}

pub fn kahan_add(acc KahanAccumulator, val f64) KahanAccumulator {
	y := val - acc.comp
	t := acc.sum + y
	return KahanAccumulator{
		sum:  t
		comp: (t - acc.sum) - y
	}
}

pub enum LoadState {
	unloaded
	full
	classify_only
}

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

pub fn (self &NbClass) get_freq_count_lg(km int) f64 {
	return self.freqcnt_lg[km]
}

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

pub fn (self &NbClass) size_bytes() u64 {
	base := u64(sizeof(NbClass))
	freq_plain := u64(self.freqcnt.len) * u64(sizeof(int) + sizeof(int))
	freq_log := u64(self.freqcnt_lg.len) * u64(sizeof(int) + sizeof(f64))
	return base + freq_plain + freq_log
}
