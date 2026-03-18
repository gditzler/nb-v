module model

import math

fn test_kahan_add_basic() {
	mut acc := KahanAccumulator{}
	acc = kahan_add(acc, 1.0)
	acc = kahan_add(acc, 2.0)
	acc = kahan_add(acc, 3.0)
	assert acc.sum == 6.0
}

fn test_kahan_add_precision() {
	mut acc := KahanAccumulator{}
	mut naive := f64(0)
	acc = kahan_add(acc, 1.0)
	naive += 1.0
	for _ in 0 .. 10000 {
		acc = kahan_add(acc, 1e-16)
		naive += 1e-16
	}
	expected := 1.0 + 10000.0 * 1e-16
	kahan_err := math.abs(acc.sum - expected)
	naive_err := math.abs(naive - expected)
	assert kahan_err <= naive_err
}

fn test_nbclass_new() {
	cls := NbClass.new('test_class', 6, '/tmp/test.nbv')
	assert cls.id == 'test_class'
	assert cls.kmer_size == 6
	assert cls.ngenomes == 0
	assert cls.sumfreq == 2080 // num_canonical_kmers(6)
	assert cls.state == .unloaded
}

fn test_nbclass_add_genome() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	kmer_counts := {1: 5, 6: 3}
	cls.add_genome(kmer_counts)

	assert cls.ngenomes == 1
	assert cls.sumfreq == 10 + 8
	assert cls.freqcnt[1] == 5
	assert cls.freqcnt[6] == 3
	assert cls.state == .full

	assert math.abs(cls.ngenomes_lg - math.log(1.0)) < 1e-10
	assert math.abs(cls.sumfreq_lg - math.log(18.0)) < 1e-10
	assert math.abs(cls.freqcnt_lg[1] - math.log(6.0)) < 1e-10
	assert math.abs(cls.freqcnt_lg[6] - math.log(4.0)) < 1e-10
}

fn test_nbclass_add_genome_twice() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3})
	cls.add_genome({1: 2, 9: 1})

	assert cls.ngenomes == 2
	assert cls.sumfreq == 10 + 8 + 3
	assert cls.freqcnt[1] == 7
	assert cls.freqcnt[6] == 3
	assert cls.freqcnt[9] == 1
}

fn test_get_freq_count_lg_seen() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5})
	assert math.abs(cls.get_freq_count_lg(1) - math.log(6.0)) < 1e-10
}

fn test_get_freq_count_lg_unseen() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5})
	assert cls.get_freq_count_lg(999) == 0.0
}

fn test_compute_log_likelihood() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3})
	read_counts := {1: 2, 6: 1}
	ll := cls.compute_log_likelihood(read_counts)
	expected := 2.0 * math.log(6.0) + 1.0 * math.log(4.0) - 3.0 * math.log(18.0)
	assert math.abs(ll - expected) < 1e-10
}

fn test_size_bytes() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3})
	bytes := cls.size_bytes()
	assert bytes > 0
}
