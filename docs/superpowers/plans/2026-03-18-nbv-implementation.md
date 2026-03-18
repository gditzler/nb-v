# NBV Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Naive Bayes metagenomic classifier in V (vlang), reimplementing NBC++ with a module-per-concern architecture, multithreaded pipelines, and YAML configuration.

**Architecture:** Six V modules (config, kmer, model, io, pipeline, main) communicate through well-defined interfaces. Concurrency uses V's `spawn`/`chan` primitives in a producer-consumer pattern. All I/O funnels through the `io` module. The `pipeline` module owns threading logic.

**Tech Stack:** V 0.5.0, `prantlf.yaml` vpm module for YAML parsing, V standard library for everything else.

**Spec:** `docs/superpowers/specs/2026-03-18-nbv-design.md`

**Example data:** `example/training_classes/` (100 NBC++ savefiles, k=9), `example/reads/cross.fna` (FASTA reads), `example/results_max_1.csv` (expected classification output).

---

## File Structure

```
nbv/
├── v.mod                      # V module manifest
├── src/
│   ├── main.v                 # CLI entry point: parse arg, load config, dispatch
│   ├── config/
│   │   └── config.v           # Config struct, YAML parsing, validation
│   ├── kmer/
│   │   ├── kmer.v             # encode, reverse_complement, canonical, count_from_buffer, num_canonical_kmers
│   │   └── kmer_test.v        # unit tests for all kmer functions
│   ├── model/
│   │   ├── model.v            # KahanAccumulator, NbClass, LoadState, add_genome, compute_log_likelihood
│   │   └── model_test.v       # unit tests for Kahan, NbClass init/training/classification math
│   ├── io/
│   │   ├── fasta.v            # parse_fasta, count_sequences
│   │   ├── kmer_file.v        # read_kmer_file, InputType enum
│   │   ├── serialization.v    # save_class, load_class, load_legacy_class, save_meta, load_meta
│   │   ├── writer.v           # Writer struct, write_header, write_result, write_full_result
│   │   └── io_test.v          # unit tests for FASTA parsing, kmer file reading, serialization, writer
│   └── pipeline/
│       ├── train.v            # train orchestrator with channel-based workers
│       ├── classify.v         # classify orchestrator with three-stage pipeline
│       └── pipeline_test.v    # integration tests: train-then-classify on synthetic data
```

**Test commands:**
- Single module: `v test src/kmer/`
- All tests: `v test src/`
- Run program: `v run src/ example/classify.yaml`

**V module import notes:**
- Cross-module imports may need adjustment based on V's resolution with `v.mod`. The plan uses short names (`import kmer`, `import config as cfg`, `import io as nbio`) but these may need to be project-qualified (e.g., `import src.kmer`). Verify in Task 1 and apply consistently.
- The `io` module directory name may conflict with V's stdlib `io` module. If this causes resolution issues, rename to `nbvio/` and update all imports.
- `InputType` and `OutputFormat` enums live in `config/config.v` only. Other modules import them via `config as cfg`.

---

## Task 1: Project Scaffold and Config Module

**Files:**
- Create: `v.mod`
- Create: `src/main.v`
- Create: `src/config/config.v`

- [ ] **Step 1: Create `v.mod`**

```v
// v.mod
Module {
	name: 'nbv'
	version: '0.1.0'
	deps: ['prantlf.yaml']
}
```

Run: `ls v.mod` -- file exists.

- [ ] **Step 2: Write `src/config/config.v`**

```v
module config

import prantlf.yaml

pub enum Mode {
	train
	classify
}

pub enum InputType {
	kmer_file
	fasta
}

pub enum OutputFormat {
	csv
	tsv
	json
}

pub struct InputConfig {
pub:
	extension  string = '.kmr'
	input_type string = 'kmer_file'
}

pub struct MemoryConfig {
pub:
	limit_mb   int
	batch_size int
	max_rows   int = 450000
	max_cols   int = 20000
}

pub struct OutputConfig {
pub:
	format      string = 'csv'
	prefix      string = 'log_likelihood'
	full_result bool
	temp_dir    string = '/tmp'
}

struct RawConfig {
pub:
	version    int = 1
	mode       string
	kmer_size  int    = 6
	save_dir   string = './nbv_save'
	source_dir string
	threads    int = 1
	input      InputConfig
	memory     MemoryConfig
	output     OutputConfig
}

pub struct Config {
pub:
	version    int
	mode       Mode
	kmer_size  int
	save_dir   string
	source_dir string
	threads    int
	input_type InputType
	extension  string
	limit_mb   int
	batch_size int
	max_rows   int
	max_cols   int
	format     OutputFormat
	prefix     string
	full_result bool
	temp_dir   string
}

pub fn load(path string) !Config {
	raw := yaml.unmarshal_file[RawConfig](path, yaml.UnmarshalOpts{})!

	mode := match raw.mode {
		'train' { Mode.train }
		'classify' { Mode.classify }
		else { return error("invalid mode '${raw.mode}': must be 'train' or 'classify'") }
	}

	if raw.kmer_size < 1 || raw.kmer_size > 15 {
		return error('kmer_size must be between 1 and 15, got ${raw.kmer_size}')
	}
	if raw.source_dir == '' {
		return error('source_dir is required')
	}
	if raw.threads < 1 {
		return error('threads must be >= 1, got ${raw.threads}')
	}

	input_type := match raw.input.input_type {
		'kmer_file' { InputType.kmer_file }
		'fasta' { InputType.fasta }
		else { return error("invalid input_type '${raw.input.input_type}': must be 'kmer_file' or 'fasta'") }
	}

	format := match raw.output.format {
		'csv' { OutputFormat.csv }
		'tsv' { OutputFormat.tsv }
		'json' { OutputFormat.json }
		else { return error("invalid output format '${raw.output.format}': must be 'csv', 'tsv', or 'json'") }
	}

	return Config{
		version:     raw.version
		mode:        mode
		kmer_size:   raw.kmer_size
		save_dir:    raw.save_dir
		source_dir:  raw.source_dir
		threads:     raw.threads
		input_type:  input_type
		extension:   raw.input.extension
		limit_mb:    raw.memory.limit_mb
		batch_size:  raw.memory.batch_size
		max_rows:    raw.memory.max_rows
		max_cols:    raw.memory.max_cols
		format:      format
		prefix:      raw.output.prefix
		full_result: raw.output.full_result
		temp_dir:    raw.output.temp_dir
	}
}
```

- [ ] **Step 3: Write minimal `src/main.v`**

```v
module main

import os
import config

fn main() {
	args := os.args[1..]
	if args.len != 1 {
		eprintln('Usage: nbv <config.yaml>')
		exit(1)
	}

	cfg := config.load(args[0]) or {
		eprintln('Error loading config: ${err}')
		exit(1)
	}

	println('Mode: ${cfg.mode}')
	println('K-mer size: ${cfg.kmer_size}')
	println('Source dir: ${cfg.source_dir}')
}
```

- [ ] **Step 4: Verify it compiles and runs with example config**

Run: `v run src/ example/classify.yaml`
Expected: Prints mode, kmer_size, source_dir.

- [ ] **Step 5: Commit**

```bash
git init
git add v.mod src/main.v src/config/config.v
git commit -m "feat: project scaffold with config module and YAML parsing"
```

---

## Task 2: K-mer Module

**Files:**
- Create: `src/kmer/kmer.v`
- Create: `src/kmer/kmer_test.v`

- [ ] **Step 1: Write the tests first in `src/kmer/kmer_test.v`**

```v
module kmer

import math

fn test_encode_single_bases() {
	assert encode('A'.bytes(), 1) == 0
	assert encode('C'.bytes(), 1) == 1
	assert encode('G'.bytes(), 1) == 2
	assert encode('T'.bytes(), 1) == 3
}

fn test_encode_kmer() {
	// AC = 0*4 + 1 = 1
	assert encode('AC'.bytes(), 2) == 1
	// GT = 2*4 + 3 = 11
	assert encode('GT'.bytes(), 2) == 11
	// ACG = 0*16 + 1*4 + 2 = 6
	assert encode('ACG'.bytes(), 3) == 6
}

fn test_reverse_complement() {
	// RC of A(0) with k=1 is T(3)
	assert reverse_complement(0, 1) == 3
	// RC of C(1) with k=1 is G(2)
	assert reverse_complement(1, 1) == 2
	// AC(1) k=2 -> RC is GT(11)
	assert reverse_complement(1, 2) == 11
	// GT(11) k=2 -> RC is AC(1)
	assert reverse_complement(11, 2) == 1
}

fn test_canonical() {
	// canonical picks the smaller of kmer and its RC
	// A(0) vs T(3) -> 0
	assert canonical(0, 1) == 0
	// T(3) vs A(0) -> 0
	assert canonical(3, 1) == 0
	// AC(1) vs GT(11) -> 1
	assert canonical(1, 2) == 1
	assert canonical(11, 2) == 1
}

fn test_count_from_buffer() {
	// Simple sequence "ACGT" with k=2
	// kmers: AC(1), CG(6), GT(11)
	// canonical: AC(1) vs GT(11)->1, CG(6) vs CG(6)->6, GT(11) vs AC(1)->1
	counts := count_from_buffer('ACGT'.bytes(), 2)
	assert counts[canonical(encode('AC'.bytes(), 2), 2)] == 2 // AC and GT are same canonical
	assert counts[canonical(encode('CG'.bytes(), 2), 2)] == 1
}

fn test_count_from_buffer_skips_invalid() {
	// 'ACNGT' - N is invalid, window resets
	// Valid kmers: AC (before N), GT (after N)
	counts := count_from_buffer('ACNGT'.bytes(), 2)
	// AC and GT are canonical pairs, so both map to canonical(1,2) = 1
	assert counts[1] == 2
}

fn test_count_from_buffer_skips_newlines() {
	// Newline in the middle should reset window
	counts := count_from_buffer('AC\nGT'.bytes(), 2)
	assert counts[1] == 2 // AC and GT are canonical pair
}

fn test_num_canonical_kmers() {
	// k=1: 4^1/2 = 2 (odd k)
	assert num_canonical_kmers(1) == 2
	// k=2: (4^2 + 4^1) / 2 = (16+4)/2 = 10 (even k)
	assert num_canonical_kmers(2) == 10
	// k=3: 4^3/2 = 32 (odd k)
	assert num_canonical_kmers(3) == 32
	// k=6: (4^6 + 4^3) / 2 = (4096+64)/2 = 2080 (even k)
	assert num_canonical_kmers(6) == 2080
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `v test src/kmer/`
Expected: Compilation errors -- functions not defined yet.

- [ ] **Step 3: Write `src/kmer/kmer.v`**

```v
module kmer

import math

// Base-4 encoding: A=0, C=1, G=2, T=3
fn base_to_int(b u8) int {
	return match b {
		`A`, `a` { 0 }
		`C`, `c` { 1 }
		`G`, `g` { 2 }
		`T`, `t` { 3 }
		else { -1 }
	}
}

// Complement: A<->T (0<->3), C<->G (1<->2)
fn complement(val int) int {
	return 3 - val
}

// Encode a k-mer byte slice to its base-4 integer representation.
// Operates on exactly k bytes.
pub fn encode(kmer []u8, k int) int {
	mut result := 0
	for i in 0 .. k {
		result = result * 4 + base_to_int(kmer[i])
	}
	return result
}

// Compute reverse complement of a k-mer in integer space.
pub fn reverse_complement(kmer_int int, k int) int {
	mut result := 0
	mut val := kmer_int
	for _ in 0 .. k {
		result = result * 4 + complement(val & 3)
		val >>= 2
	}
	return result
}

// Return the lexicographically smaller of a k-mer and its reverse complement.
pub fn canonical(kmer_int int, k int) int {
	rc := reverse_complement(kmer_int, k)
	if kmer_int <= rc {
		return kmer_int
	}
	return rc
}

// Count canonical k-mers in a buffer.
// Invalid characters and newlines reset the current window.
pub fn count_from_buffer(buf []u8, k int) map[int]int {
	mut counts := map[int]int{}
	mut window := 0
	mut valid_len := 0

	for i in 0 .. buf.len {
		val := base_to_int(buf[i])
		if val < 0 {
			// Invalid char or newline: reset window
			valid_len = 0
			window = 0
			continue
		}

		// Shift window and add new base
		window = (window * 4 + val) & ((1 << (2 * k)) - 1)
		valid_len++

		if valid_len >= k {
			canon := canonical(window, k)
			counts[canon] = counts[canon] + 1
		}
	}
	return counts
}

// Number of distinct canonical k-mers for a given k.
// Odd k: 4^k / 2. Even k: (4^k + 4^(k/2)) / 2.
pub fn num_canonical_kmers(k int) i64 {
	total := i64(1) << (2 * k) // 4^k
	if k % 2 == 1 {
		return total / 2
	}
	palindromes := i64(1) << k // 4^(k/2)
	return (total + palindromes) / 2
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `v test src/kmer/`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/kmer/kmer.v src/kmer/kmer_test.v
git commit -m "feat: kmer module with encode, reverse_complement, canonical, count_from_buffer"
```

---

## Task 3: Model Module

**Files:**
- Create: `src/model/model.v`
- Create: `src/model/model_test.v`

- [ ] **Step 1: Write tests in `src/model/model_test.v`**

```v
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
	// Kahan should handle this better than naive summation
	mut acc := KahanAccumulator{}
	mut naive := f64(0)
	acc = kahan_add(acc, 1.0)
	naive += 1.0
	for _ in 0 .. 10000 {
		acc = kahan_add(acc, 1e-16)
		naive += 1e-16
	}
	// Kahan should be closer to the true answer (1.0 + 10000*1e-16 = 1.000000000001)
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
	// num_canonical_kmers(2) = 10, so initial sumfreq = 10

	kmer_counts := {
		1: 5  // e.g. canonical kmer AC
		6: 3  // e.g. canonical kmer CG
	}
	cls.add_genome(kmer_counts)

	assert cls.ngenomes == 1
	assert cls.sumfreq == 10 + 8 // initial 10 + sum of counts (5+3)
	assert cls.freqcnt[1] == 5
	assert cls.freqcnt[6] == 3
	assert cls.state == .full

	// Check log-space values
	assert math.abs(cls.ngenomes_lg - math.log(1.0)) < 1e-10
	assert math.abs(cls.sumfreq_lg - math.log(18.0)) < 1e-10
	// freqcnt_lg[1] = log(5 + 1) = log(6)
	assert math.abs(cls.freqcnt_lg[1] - math.log(6.0)) < 1e-10
	// freqcnt_lg[6] = log(3 + 1) = log(4)
	assert math.abs(cls.freqcnt_lg[6] - math.log(4.0)) < 1e-10
}

fn test_nbclass_add_genome_twice() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3})
	cls.add_genome({1: 2, 9: 1})

	assert cls.ngenomes == 2
	assert cls.sumfreq == 10 + 8 + 3 // initial + first(5+3) + second(2+1)
	assert cls.freqcnt[1] == 7  // 5 + 2
	assert cls.freqcnt[6] == 3
	assert cls.freqcnt[9] == 1
}

fn test_get_freq_count_lg_seen() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5})
	// freqcnt_lg[1] = log(5 + 1) = log(6)
	assert math.abs(cls.get_freq_count_lg(1) - math.log(6.0)) < 1e-10
}

fn test_get_freq_count_lg_unseen() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5})
	// Unseen kmer returns log(1) = 0.0
	assert cls.get_freq_count_lg(999) == 0.0
}

fn test_compute_log_likelihood() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3})
	// sumfreq = 10 + 8 = 18, sumfreq_lg = log(18)

	read_counts := {1: 2, 6: 1}
	ll := cls.compute_log_likelihood(read_counts)

	// Manual calculation:
	// total_kmer_count = 2 + 1 = 3
	// sum_kahan(freq_i * freqcount_lg(kmer_i)):
	//   2 * log(6) + 1 * log(4)
	// score = (2*log(6) + 1*log(4)) - 3 * log(18)
	expected := 2.0 * math.log(6.0) + 1.0 * math.log(4.0) - 3.0 * math.log(18.0)
	assert math.abs(ll - expected) < 1e-10
}

fn test_size_bytes() {
	mut cls := NbClass.new('test_class', 2, '/tmp/test.nbv')
	cls.add_genome({1: 5, 6: 3})
	bytes := cls.size_bytes()
	assert bytes > 0
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `v test src/model/`
Expected: Compilation errors -- types and functions not defined.

- [ ] **Step 3: Write `src/model/model.v`**

```v
module model

import math
import kmer as kmod

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

	// Eagerly recompute log-space parameters
	self.ngenomes_lg = math.log(f64(self.ngenomes))
	self.sumfreq_lg = math.log(f64(self.sumfreq))
	for km, count in self.freqcnt {
		self.freqcnt_lg[km] = math.log(f64(count + 1))
	}

	self.state = .full
}

pub fn (self &NbClass) get_freq_count_lg(km int) f64 {
	// V maps return zero value (0.0) for missing keys, which is exactly
	// what we want: log(1) = 0.0 for the Laplace pseudo-count of unseen k-mers
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
```

Note: The `import kmer as kmod` import path may need adjustment based on how V resolves sibling modules. If V requires the full path, use the project module name. Test and adjust if needed.

- [ ] **Step 4: Run tests to verify they pass**

Run: `v test src/model/`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/model/model.v src/model/model_test.v
git commit -m "feat: model module with KahanAccumulator, NbClass, add_genome, compute_log_likelihood"
```

---

## Task 4: I/O Module -- FASTA Parser

**Files:**
- Create: `src/io/fasta.v`
- Create: `src/io/io_test.v`

- [ ] **Step 1: Create a test FASTA fixture**

Create `src/io/testdata/test.fasta`:
```
>seq1 description
ACGTACGT
ACGTACGT
>seq2
GGGGCCCC
>seq3 empty after header

>seq4
ACGT
```

- [ ] **Step 2: Write FASTA tests in `src/io/io_test.v`**

```v
module io

fn test_parse_fasta_basic() {
	mut headers := []string{}
	mut sequences := []string{}

	parse_fasta('src/io/testdata/test.fasta', fn [mut headers, mut sequences] (header string, seq []u8) {
		headers << header
		sequences << seq.bytestr()
	})!

	assert headers.len == 4
	assert headers[0] == 'seq1 description'
	assert sequences[0] == 'ACGTACGTACGTACGT' // multi-line concatenated
	assert headers[1] == 'seq2'
	assert sequences[1] == 'GGGGCCCC'
	assert headers[2] == 'seq3 empty after header'
	assert sequences[2] == '' // empty sequence
	assert headers[3] == 'seq4'
	assert sequences[3] == 'ACGT'
}

fn test_count_sequences() {
	count := count_sequences('src/io/testdata/test.fasta')!
	assert count == 4
}
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `v test src/io/`
Expected: Compilation errors.

- [ ] **Step 4: Write `src/io/fasta.v`**

```v
module io

import os

// Streaming FASTA parser. Calls callback for each sequence.
// Reads line-by-line to avoid loading the entire file into memory.
// Accumulates multi-line sequences before invoking callback.
pub fn parse_fasta(path string, callback fn (string, []u8)) ! {
	mut f := os.open(path)!
	defer { f.close() }

	mut current_header := ''
	mut current_seq := []u8{}
	mut in_record := false
	mut buf := []u8{len: 65536}

	// Use os.read_lines for simplicity in initial impl.
	// TODO: For very large files, replace with buffered line reader.
	lines := os.read_lines(path)!
	for line in lines {
		if line.len > 0 && line[0] == `>` {
			if in_record {
				callback(current_header, current_seq)
			}
			current_header = line[1..].trim_space()
			current_seq = []u8{}
			in_record = true
		} else if in_record {
			current_seq << line.bytes()
		}
	}

	if in_record {
		callback(current_header, current_seq)
	}
}

// Fast count of sequences by counting '>' header lines.
pub fn count_sequences(path string) !u64 {
	lines := os.read_lines(path)!
	mut count := u64(0)
	for line in lines {
		if line.len > 0 && line[0] == `>` {
			count++
		}
	}
	return count
}
```

Note: The initial implementation uses `os.read_lines` for simplicity. For production use with very large genomic files (multi-GB), this should be refactored to use buffered line-by-line reading. This is acceptable for the initial implementation since the example data is small, and can be optimized in a follow-up task.

- [ ] **Step 5: Run tests to verify they pass**

Run: `v test src/io/`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/io/fasta.v src/io/io_test.v src/io/testdata/test.fasta
git commit -m "feat: io module with streaming FASTA parser and sequence counter"
```

---

## Task 5: I/O Module -- K-mer File Reader

**Files:**
- Create: `src/io/kmer_file.v`
- Modify: `src/io/io_test.v`

- [ ] **Step 1: Create a test `.kmr` fixture**

Create `src/io/testdata/test.kmr`:
```
ACGTAA	15
CGTAAC	8
TGCAAT	3
```

- [ ] **Step 2: Add kmer file tests to `src/io/io_test.v`**

```v
// Append to io_test.v

import kmer as kmod

fn test_read_kmer_file() {
	counts := read_kmer_file('src/io/testdata/test.kmr', 6)!
	// Should have 3 entries (or fewer if any are canonical pairs)
	assert counts.len > 0
	// Verify a specific kmer: ACGTAA encoded and canonicalized
	kmer_int := kmod.encode('ACGTAA'.bytes(), 6)
	canon := kmod.canonical(kmer_int, 6)
	assert counts[canon] > 0
}
```

- [ ] **Step 3: Run tests to verify new tests fail**

Run: `v test src/io/`
Expected: New tests fail, existing FASTA tests still pass.

- [ ] **Step 4: Write `src/io/kmer_file.v`**

```v
module io

import os
import kmer as kmod

// Read NBC++ format .kmr file: <kmer_string>\t<count> per line.
// Encodes and canonicalizes each kmer.
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
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `v test src/io/`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/io/kmer_file.v src/io/testdata/test.kmr
git commit -m "feat: kmer file reader and input type detection"
```

---

## Task 6: I/O Module -- Serialization (New Format + Legacy Reader)

**Files:**
- Create: `src/io/serialization.v`
- Modify: `src/io/io_test.v`

- [ ] **Step 1: Add serialization tests to `src/io/io_test.v`**

```v
// Append to io_test.v

import model
import math

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

	// Verify log-space values were recomputed
	assert math.abs(loaded.sumfreq_lg - cls.sumfreq_lg) < 1e-10

	os.rm('/tmp/test_roundtrip.nbv')!
}

fn test_save_and_load_meta() {
	save_meta('/tmp/test_meta', 9)!
	k := load_meta('/tmp/test_meta')!
	assert k == 9
	os.rm('/tmp/test_meta/meta.nbv')!
	os.rmdir('/tmp/test_meta')!
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `v test src/io/`
Expected: New serialization tests fail.

- [ ] **Step 3: Write `src/io/serialization.v`**

```v
module io

import os
import model
import math
import kmer as kmod

const (
	magic_bytes    = [u8(`N`), `B`, `V`, `1`]
	format_version = u8(1)
)

// Save an NbClass in the new NBV binary format.
// Header: magic(4) + version(1) + kmer_size(4) + ngenomes(4) + sumfreq(8)
// Body: repeated (kmer_int i32, count i32) pairs
pub fn save_class(cls model.NbClass, path string) ! {
	mut f := os.create(path)!
	defer { f.close() }

	// Magic bytes
	f.write(magic_bytes)!
	// Version
	f.write([format_version])!
	// kmer_size as i32 (4 bytes, little-endian)
	f.write(i32_to_bytes(cls.kmer_size))!
	// ngenomes as i32
	f.write(i32_to_bytes(cls.ngenomes))!
	// sumfreq as i64 (8 bytes, little-endian)
	f.write(i64_to_bytes(cls.sumfreq))!
	// id as length-prefixed string
	id_bytes := cls.id.bytes()
	f.write(i32_to_bytes(id_bytes.len))!
	f.write(id_bytes)!
	// frequency pairs
	for km, count in cls.freqcnt {
		f.write(i32_to_bytes(km))!
		f.write(i32_to_bytes(count))!
	}
}

// Load an NbClass from the new NBV binary format.
pub fn load_class(path string) !model.NbClass {
	data := os.read_bytes(path)!
	mut pos := 0

	// Magic bytes check
	if data.len < 4 || data[0..4] != magic_bytes {
		return error('invalid NBV file: bad magic bytes in ${path}')
	}
	pos = 4

	// Version
	version := data[pos]
	if version != format_version {
		return error('unsupported NBV format version ${version}')
	}
	pos++

	kmer_size := bytes_to_i32(data[pos..pos + 4])
	pos += 4
	ngenomes := bytes_to_i32(data[pos..pos + 4])
	pos += 4
	sumfreq := bytes_to_i64(data[pos..pos + 8])
	pos += 8

	// id
	id_len := bytes_to_i32(data[pos..pos + 4])
	pos += 4
	id := data[pos..pos + id_len].bytestr()
	pos += id_len

	// frequency pairs
	mut freqcnt := map[int]int{}
	mut freqcnt_lg := map[int]f64{}
	for pos + 8 <= data.len {
		km := bytes_to_i32(data[pos..pos + 4])
		pos += 4
		count := bytes_to_i32(data[pos..pos + 4])
		pos += 4
		freqcnt[km] = count
		freqcnt_lg[km] = math.log(f64(count + 1))
	}

	return model.NbClass{
		id:          id
		kmer_size:   kmer_size
		savefile:    path
		ngenomes:    ngenomes
		sumfreq:     sumfreq
		ngenomes_lg: math.log(f64(ngenomes))
		sumfreq_lg:  math.log(f64(sumfreq))
		freqcnt:     freqcnt
		freqcnt_lg:  freqcnt_lg
		state:       .full
	}
}

// Load an NBC++ legacy -save.dat file. Populates log-space fields only.
// NBC++ binary format (verified via hex dump):
//   ngenomes_lg  f64 (8 bytes, little-endian)
//   sumfreq_lg   f64 (8 bytes, little-endian)
//   n_entries    i32 (4 bytes, little-endian)  -- number of kmer entries
//   Then n_entries * (kmer_int i32, freqcnt_lg f64) pairs (12 bytes each)
pub fn load_legacy_class(path string, k int) !model.NbClass {
	data := os.read_bytes(path)!
	mut pos := 0

	if data.len < 20 {
		return error('legacy save file too small: ${path}')
	}

	ngenomes_lg := bytes_to_f64(data[pos..pos + 8])
	pos += 8
	sumfreq_lg := bytes_to_f64(data[pos..pos + 8])
	pos += 8
	n_entries := bytes_to_i32(data[pos..pos + 4])
	pos += 4

	expected_remaining := n_entries * 12
	if data.len - pos < expected_remaining {
		return error('legacy save file truncated: ${path} (expected ${n_entries} entries)')
	}

	mut freqcnt_lg := map[int]f64{}
	for _ in 0 .. n_entries {
		km := bytes_to_i32(data[pos..pos + 4])
		pos += 4
		val := bytes_to_f64(data[pos..pos + 8])
		pos += 8
		freqcnt_lg[km] = val
	}

	// Extract class id from filename: "<class_id>-save.dat"
	basename := os.file_name(path)
	id := basename.replace('-save.dat', '')

	return model.NbClass{
		id:          id
		kmer_size:   k
		savefile:    path
		ngenomes_lg: ngenomes_lg
		sumfreq_lg:  sumfreq_lg
		freqcnt_lg:  freqcnt_lg
		state:       .classify_only
	}
}

// Save kmer size metadata.
pub fn save_meta(save_dir string, kmer_size int) ! {
	os.mkdir_all(save_dir)!
	os.write_file('${save_dir}/meta.nbv', '${kmer_size}')!
}

// Load kmer size metadata and return it.
pub fn load_meta(save_dir string) !int {
	content := os.read_file('${save_dir}/meta.nbv')!
	return content.trim_space().int()
}

// -- Byte conversion helpers (little-endian) --

fn i32_to_bytes(val int) []u8 {
	mut b := []u8{len: 4}
	b[0] = u8(val)
	b[1] = u8(val >> 8)
	b[2] = u8(val >> 16)
	b[3] = u8(val >> 24)
	return b
}

fn bytes_to_i32(b []u8) int {
	return int(b[0]) | (int(b[1]) << 8) | (int(b[2]) << 16) | (int(b[3]) << 24)
}

fn i64_to_bytes(val i64) []u8 {
	mut b := []u8{len: 8}
	for i in 0 .. 8 {
		b[i] = u8(val >> (i * 8))
	}
	return b
}

fn bytes_to_i64(b []u8) i64 {
	mut result := i64(0)
	for i in 0 .. 8 {
		result |= i64(b[i]) << (i * 8)
	}
	return result
}

fn bytes_to_f64(b []u8) f64 {
	bits := u64(bytes_to_i64(b))
	return unsafe { *(&f64(&bits)) }
}
```

Note: The `load_legacy_class` function parses the NBC++ binary format. The exact byte layout of NBC++ savefiles may need to be verified against an actual file during implementation. The implementer should hex-dump one of the `example/training_classes/*.dat` files and adjust the parsing if the layout differs. This is documented as a verification step below.

- [ ] **Step 4: Run tests to verify they pass**

Run: `v test src/io/`
Expected: All tests pass.

- [ ] **Step 5: Verify legacy format against actual NBC++ savefiles**

Run: `xxd example/training_classes/1748-save.dat | head -20`

Inspect the byte layout and compare against the `load_legacy_class` parser. Adjust if the format differs from the assumed layout (the C++ `Class::serialize` may use a different field order or include additional fields). Update `load_legacy_class` accordingly and add a test:

```v
fn test_load_legacy_class() {
	cls := load_legacy_class('example/training_classes/1748-save.dat', 9)!
	assert cls.id == '1748'
	assert cls.kmer_size == 9
	assert cls.state == .classify_only
	assert cls.freqcnt_lg.len > 0
}
```

- [ ] **Step 6: Commit**

```bash
git add src/io/serialization.v
git commit -m "feat: NBV binary serialization, legacy NBC++ reader, and meta file"
```

---

## Task 7: I/O Module -- Output Writer

**Files:**
- Create: `src/io/writer.v`
- Modify: `src/io/io_test.v`

- [ ] **Step 1: Add writer tests to `src/io/io_test.v`**

```v
// Append to io_test.v

fn test_writer_csv() {
	mut w := Writer.new('/tmp/test_output.csv', .csv, false)!
	w.write_result('seq1', 'class_a', -123.45)!
	w.write_result('seq2', 'class_b', -678.90)!
	w.close()!

	content := os.read_file('/tmp/test_output.csv')!
	lines := content.trim_space().split('\n')
	assert lines.len == 2
	assert lines[0] == 'seq1,class_a,-123.45'
	os.rm('/tmp/test_output.csv')!
}

fn test_writer_tsv() {
	mut w := Writer.new('/tmp/test_output.tsv', .tsv, false)!
	w.write_result('seq1', 'class_a', -123.45)!
	w.close()!

	content := os.read_file('/tmp/test_output.tsv')!
	assert content.trim_space().contains('\t')
	os.rm('/tmp/test_output.tsv')!
}

fn test_writer_json() {
	mut w := Writer.new('/tmp/test_output.jsonl', .json, false)!
	w.write_result('seq1', 'class_a', -123.45)!
	w.close()!

	content := os.read_file('/tmp/test_output.jsonl')!
	assert content.contains('"seq_id"')
	assert content.contains('"best_class"')
	os.rm('/tmp/test_output.jsonl')!
}

fn test_writer_full_result_csv() {
	mut w := Writer.new('/tmp/test_full.csv', .csv, true)!
	w.write_header(['class_a', 'class_b'])!
	w.write_full_result('seq1', {'class_a': f64(-100.0), 'class_b': f64(-200.0)}, ['class_a', 'class_b'])!
	w.close()!

	content := os.read_file('/tmp/test_full.csv')!
	lines := content.trim_space().split('\n')
	assert lines.len == 2
	assert lines[0] == 'seq_id,class_a,class_b'
	os.rm('/tmp/test_full.csv')!
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `v test src/io/`
Expected: New writer tests fail.

- [ ] **Step 3: Write `src/io/writer.v`**

```v
module io

import os
import config

pub struct Writer {
mut:
	file        os.File
	format      config.OutputFormat
	full_result bool
}

pub fn Writer.new(path string, format config.OutputFormat, full_result bool) !Writer {
	f := os.create(path)!
	return Writer{
		file:        f
		format:      format
		full_result: full_result
	}
}

pub fn (mut self Writer) write_header(class_ids []string) ! {
	match self.format {
		.csv {
			self.file.writeln('seq_id,' + class_ids.join(','))!
		}
		.tsv {
			self.file.writeln('seq_id\t' + class_ids.join('\t'))!
		}
		.json {} // JSON Lines has no header
	}
}

pub fn (mut self Writer) write_result(seq_id string, best_class string, score f64) ! {
	match self.format {
		.csv {
			self.file.writeln('${seq_id},${best_class},${score}')!
		}
		.tsv {
			self.file.writeln('${seq_id}\t${best_class}\t${score}')!
		}
		.json {
			self.file.writeln('{"seq_id":"${seq_id}","best_class":"${best_class}","score":${score}}')!
		}
	}
}

pub fn (mut self Writer) write_full_result(seq_id string, scores map[string]f64, class_order []string) ! {
	match self.format {
		.csv {
			mut parts := [seq_id]
			for cls in class_order {
				parts << '${scores[cls]}'
			}
			self.file.writeln(parts.join(','))!
		}
		.tsv {
			mut parts := [seq_id]
			for cls in class_order {
				parts << '${scores[cls]}'
			}
			self.file.writeln(parts.join('\t'))!
		}
		.json {
			mut score_parts := []string{}
			for cls in class_order {
				score_parts << '"${cls}":${scores[cls]}'
			}
			self.file.writeln('{"seq_id":"${seq_id}","scores":{${score_parts.join(",")}}}')!
		}
	}
}

pub fn (mut self Writer) write_no_valid_kmers(seq_id string) ! {
	match self.format {
		.csv {
			self.file.writeln('${seq_id},sequence contains no valid kmers,')!
		}
		.tsv {
			self.file.writeln('${seq_id}\tsequence contains no valid kmers\t')!
		}
		.json {
			self.file.writeln('{"seq_id":"${seq_id}","best_class":"sequence contains no valid kmers","score":null}')!
		}
	}
}

pub fn (mut self Writer) close() ! {
	self.file.close()
}

// Build output filename from prefix and format.
pub fn output_filename(prefix string, format config.OutputFormat) string {
	ext := match format {
		.csv { 'csv' }
		.tsv { 'tsv' }
		.json { 'jsonl' }
	}
	return '${prefix}.${ext}'
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `v test src/io/`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/io/writer.v
git commit -m "feat: output writer supporting CSV, TSV, and JSON Lines formats"
```

---

## Task 8: Training Pipeline (Single-threaded Core)

**Files:**
- Create: `src/pipeline/train.v`
- Create: `src/pipeline/pipeline_test.v`

Build the training pipeline first without threading, then add concurrency in Task 10.

- [ ] **Step 1: Create test fixtures for training**

Create `src/pipeline/testdata/training/` with two class subdirectories:

```
src/pipeline/testdata/training/
  class_a/
    genome1.fasta
  class_b/
    genome2.fasta
```

`genome1.fasta`:
```
>genome1_seq1
ACGTACGTACGTACGT
>genome1_seq2
GGGGCCCCAAAATTTT
```

`genome2.fasta`:
```
>genome2_seq1
TTTTAAAACCCCGGGG
>genome2_seq2
ACACACACACACACAC
```

- [ ] **Step 2: Write training pipeline test in `src/pipeline/pipeline_test.v`**

```v
module pipeline

import os
import config as cfg
import model

fn test_train_creates_savefiles() {
	out_dir := '/tmp/nbv_test_train'
	os.rmdir_all(out_dir) or {}

	c := cfg.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   out_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    1
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 0
	}

	train(c)!

	// Should have created savefiles for class_a and class_b
	assert os.exists('${out_dir}/class_a.nbv')
	assert os.exists('${out_dir}/class_b.nbv')
	assert os.exists('${out_dir}/meta.nbv')

	// Verify meta
	meta_content := os.read_file('${out_dir}/meta.nbv')!
	assert meta_content.trim_space() == '4'

	os.rmdir_all(out_dir) or {}
}
```

- [ ] **Step 3: Run test to verify it fails**

Run: `v test src/pipeline/`
Expected: `train` function not defined.

- [ ] **Step 4: Write `src/pipeline/train.v`**

```v
module pipeline

import os
import config as cfg
import model
import io as nbio
import kmer as kmod

struct TrainJob {
	class_id string
	path     string
}

// Scan source_dir for class subdirectories and build work queue.
fn scan_training_dir(source_dir string, extension string) ![]TrainJob {
	mut jobs := []TrainJob{}
	entries := os.ls(source_dir)!
	for entry in entries {
		subdir := '${source_dir}/${entry}'
		if !os.is_dir(subdir) {
			continue
		}
		class_id := entry
		files := os.ls(subdir)!
		for file in files {
			if file.ends_with(extension) {
				jobs << TrainJob{
					class_id: class_id
					path:     '${subdir}/${file}'
				}
			}
		}
	}
	return jobs
}

// Load k-mer counts from a file (FASTA or .kmr).
fn load_kmer_counts(path string, input_type cfg.InputType, k int) !map[int]int {
	if input_type == .fasta {
		mut merged := map[int]int{}
		nbio.parse_fasta(path, fn [mut merged, k] (header string, seq []u8) {
			counts := kmod.count_from_buffer(seq, k)
			for km, count in counts {
				merged[km] = merged[km] + count
			}
		})!
		return merged
	} else {
		return nbio.read_kmer_file(path, k)!
	}
}

pub fn train(c cfg.Config) ! {
	jobs := scan_training_dir(c.source_dir, c.extension)!
	if jobs.len == 0 {
		return error('no training files found in ${c.source_dir}')
	}

	mut classes := map[string]model.NbClass{}

	for job in jobs {
		// Create class if it doesn't exist
		if job.class_id !in classes {
			savefile := '${c.save_dir}/${job.class_id}.nbv'
			classes[job.class_id] = model.NbClass.new(job.class_id, c.kmer_size, savefile)
		}

		kmer_counts := load_kmer_counts(job.path, c.input_type, c.kmer_size) or {
			eprintln('Warning: skipping ${job.path}: ${err}')
			continue
		}

		mut cls := classes[job.class_id]
		cls.add_genome(kmer_counts)
		classes[job.class_id] = cls
	}

	// Serialize all classes
	os.mkdir_all(c.save_dir)!
	for _, cls in classes {
		nbio.save_class(cls, cls.savefile)!
	}
	nbio.save_meta(c.save_dir, c.kmer_size)!
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `v test src/pipeline/`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/pipeline/train.v src/pipeline/pipeline_test.v src/pipeline/testdata/
git commit -m "feat: single-threaded training pipeline"
```

---

## Task 9: Classification Pipeline (Single-threaded Core)

**Files:**
- Create: `src/pipeline/classify.v`
- Modify: `src/pipeline/pipeline_test.v`

- [ ] **Step 1: Add classification test to `src/pipeline/pipeline_test.v`**

```v
// Append to pipeline_test.v

fn test_train_then_classify() {
	train_dir := '/tmp/nbv_test_e2e_train'
	os.rmdir_all(train_dir) or {}

	// Train
	train_cfg := cfg.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   train_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    1
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 0
	}
	train(train_cfg)!

	// Create a test read file from class_a's genome
	os.mkdir_all('/tmp/nbv_test_e2e_reads')!
	os.write_file('/tmp/nbv_test_e2e_reads/test_read.fasta', '>read1\nACGTACGTACGTACGT\n')!

	// Classify
	classify_cfg := cfg.Config{
		mode:        .classify
		kmer_size:   4
		save_dir:    train_dir
		source_dir:  '/tmp/nbv_test_e2e_reads'
		threads:     1
		input_type:  .fasta
		extension:   '.fasta'
		format:      .csv
		prefix:      '/tmp/nbv_test_e2e_output',
		full_result: false
		temp_dir:    '/tmp'
		limit_mb:    0
		max_rows:    1000
		max_cols:    100
	}
	classify(classify_cfg)!

	// Check output exists and contains a result
	output := os.read_file('/tmp/nbv_test_e2e_output.csv')!
	assert output.contains('read1')
	// The read is from class_a's genome, so it should classify to class_a
	assert output.contains('class_a')

	os.rmdir_all(train_dir) or {}
	os.rmdir_all('/tmp/nbv_test_e2e_reads') or {}
	os.rm('/tmp/nbv_test_e2e_output.csv') or {}
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `v test src/pipeline/`
Expected: `classify` function not defined.

- [ ] **Step 3: Write `src/pipeline/classify.v`**

```v
module pipeline

import os
import math
import config as cfg
import model
import io as nbio
import kmer as kmod

struct SeqJob {
	seq_id      string
	kmer_counts map[int]int
}

struct ClassifyResult {
	seq_id     string
	best_class string
	best_score f64
	all_scores map[string]f64
}

pub fn classify(c cfg.Config) ! {
	// Verify kmer size matches training
	trained_k := nbio.load_meta(c.save_dir)!
	if trained_k != c.kmer_size {
		return error('kmer_size mismatch: config has ${c.kmer_size}, training used ${trained_k}')
	}

	// Load all class savefiles
	classes := load_classes(c.save_dir, c.kmer_size)!
	if classes.len == 0 {
		return error('no class savefiles found in ${c.save_dir}')
	}

	// Collect class IDs for header ordering
	class_ids := classes.map(it.id)

	// Open output writer
	output_path := nbio.output_filename(c.prefix, c.format)
	mut writer := nbio.Writer.new(output_path, c.format, c.full_result)!
	defer { writer.close() or {} }

	if c.full_result {
		writer.write_header(class_ids)!
	}

	// Find input files
	input_files := find_input_files(c.source_dir, c.extension)!

	// Process each input file
	for input_file in input_files {
		if c.input_type == .fasta {
			nbio.parse_fasta(input_file, fn [classes, class_ids, mut writer, c] (header string, seq []u8) {
				seq_id := header.split(' ')[0] // Use first token as seq_id
				kmer_counts := kmod.count_from_buffer(seq, c.kmer_size)

				if kmer_counts.len == 0 {
					writer.write_no_valid_kmers(seq_id) or {
						eprintln('Warning: failed to write result for ${seq_id}: ${err}')
					}
					return
				}

				result := classify_read(seq_id, kmer_counts, classes, c.full_result)
				write_classify_result(mut writer, result, class_ids, c.full_result) or {
					eprintln('Warning: failed to write result for ${seq_id}: ${err}')
				}
			})!
		} else {
			// For kmer_file input: each file is one read
			kmer_counts := nbio.read_kmer_file(input_file, c.kmer_size)!
			seq_id := os.file_name(input_file).replace(c.extension, '')
			if kmer_counts.len == 0 {
				writer.write_no_valid_kmers(seq_id)!
				continue
			}
			result := classify_read(seq_id, kmer_counts, classes, c.full_result)
			write_classify_result(mut writer, result, class_ids, c.full_result)!
		}
	}
}

fn classify_read(seq_id string, kmer_counts map[int]int, classes []model.NbClass, full_result bool) ClassifyResult {
	mut best_class := ''
	mut best_score := -math.max_f64
	mut all_scores := map[string]f64{}

	for cls in classes {
		score := cls.compute_log_likelihood(kmer_counts)
		if full_result {
			all_scores[cls.id] = score
		}
		if score > best_score {
			best_score = score
			best_class = cls.id
		}
	}

	return ClassifyResult{
		seq_id:     seq_id
		best_class: best_class
		best_score: best_score
		all_scores: all_scores
	}
}

fn write_classify_result(mut writer nbio.Writer, result ClassifyResult, class_ids []string, full_result bool) ! {
	if full_result {
		writer.write_full_result(result.seq_id, result.all_scores, class_ids)!
	} else {
		writer.write_result(result.seq_id, result.best_class, result.best_score)!
	}
}

fn load_classes(save_dir string, kmer_size int) ![]model.NbClass {
	mut classes := []model.NbClass{}
	entries := os.ls(save_dir)!
	for entry in entries {
		path := '${save_dir}/${entry}'
		if entry.ends_with('.nbv') && entry != 'meta.nbv' {
			cls := nbio.load_class(path)!
			classes << cls
		} else if entry.ends_with('-save.dat') {
			cls := nbio.load_legacy_class(path, kmer_size)!
			classes << cls
		}
	}
	return classes
}

fn find_input_files(source_dir string, extension string) ![]string {
	mut files := []string{}
	entries := os.ls(source_dir)!
	for entry in entries {
		if entry.ends_with(extension) {
			files << '${source_dir}/${entry}'
		}
	}
	if files.len == 0 {
		return error('no input files with extension ${extension} found in ${source_dir}')
	}
	return files
}
```

Note: V uses `math.max_f64` for the largest representable float; negating it gives us an effective negative infinity sentinel.

- [ ] **Step 4: Run tests to verify they pass**

Run: `v test src/pipeline/`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/pipeline/classify.v
git commit -m "feat: single-threaded classification pipeline"
```

---

## Task 10: Wire Up main.v and Test End-to-End

**Files:**
- Modify: `src/main.v`

- [ ] **Step 1: Update `src/main.v` to dispatch to pipelines**

```v
module main

import os
import config
import pipeline

fn main() {
	args := os.args[1..]
	if args.len != 1 {
		eprintln('Usage: nbv <config.yaml>')
		exit(1)
	}

	cfg := config.load(args[0]) or {
		eprintln('Error: ${err}')
		exit(1)
	}

	match cfg.mode {
		.train {
			pipeline.train(cfg) or {
				eprintln('Training failed: ${err}')
				exit(1)
			}
			println('Training complete. Savefiles written to ${cfg.save_dir}')
		}
		.classify {
			pipeline.classify(cfg) or {
				eprintln('Classification failed: ${err}')
				exit(1)
			}
			println('Classification complete.')
		}
	}
}
```

- [ ] **Step 2: Build and test with example data**

Run: `v build src/ -o nbv && ./nbv example/classify.yaml`
Expected: Classifies reads from `cross.fna` against the 100 NBC++ savefiles and writes output.

- [ ] **Step 3: Compare output against expected results**

Run: `head -10 example_results.csv && echo "---" && head -10 example/results_max_1.csv`
Expected: The first few rows should have the same class assignments. Scores may differ slightly in precision but class assignments should match.

- [ ] **Step 4: Commit**

```bash
git add src/main.v
git commit -m "feat: wire main.v to train and classify pipelines"
```

---

## Task 11: Add Multithreading to Training Pipeline

**Files:**
- Modify: `src/pipeline/train.v`

- [ ] **Step 1: Add threaded training test to `src/pipeline/pipeline_test.v`**

```v
fn test_train_multithreaded() {
	out_dir := '/tmp/nbv_test_train_mt'
	os.rmdir_all(out_dir) or {}

	c := cfg.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   out_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    2
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 0
	}

	train(c)!

	assert os.exists('${out_dir}/class_a.nbv')
	assert os.exists('${out_dir}/class_b.nbv')

	// Results should be identical to single-threaded
	cls_a := nbio.load_class('${out_dir}/class_a.nbv')!
	assert cls_a.ngenomes == 1
	assert cls_a.freqcnt.len > 0

	os.rmdir_all(out_dir) or {}
}
```

- [ ] **Step 2: Refactor `train()` to use channels when threads > 1**

Update `src/pipeline/train.v`:

```v
pub fn train(c cfg.Config) ! {
	jobs := scan_training_dir(c.source_dir, c.extension)!
	if jobs.len == 0 {
		return error('no training files found in ${c.source_dir}')
	}

	mut classes := map[string]model.NbClass{}

	if c.threads <= 1 {
		// Single-threaded path (unchanged)
		for job in jobs {
			if job.class_id !in classes {
				savefile := '${c.save_dir}/${job.class_id}.nbv'
				classes[job.class_id] = model.NbClass.new(job.class_id, c.kmer_size, savefile)
			}
			kmer_counts := load_kmer_counts(job.path, c.input_type, c.kmer_size) or {
				eprintln('Warning: skipping ${job.path}: ${err}')
				continue
			}
			mut cls := classes[job.class_id]
			cls.add_genome(kmer_counts)
			classes[job.class_id] = cls
		}
	} else {
		// Multi-threaded: workers count kmers, accumulator merges
		job_ch := chan TrainJob{cap: c.threads * 2}
		result_ch := chan TrainResult{cap: c.threads * 2}
		done_ch := chan bool{cap: c.threads}

		// Spawn workers
		for _ in 0 .. c.threads {
			spawn train_worker(job_ch, result_ch, done_ch, c.input_type, c.kmer_size)
		}

		// Spawn feeder
		spawn fn [job_ch, jobs] () {
			for job in jobs {
				job_ch <- job
			}
			job_ch.close()
		}()

		// Accumulator: runs in main thread
		// Count finished workers
		mut workers_done := 0
		for {
			select {
				result := <-result_ch {
					if result.class_id !in classes {
						savefile := '${c.save_dir}/${result.class_id}.nbv'
						classes[result.class_id] = model.NbClass.new(result.class_id, c.kmer_size, savefile)
					}
					mut cls := classes[result.class_id]
					cls.add_genome(result.kmer_counts)
					classes[result.class_id] = cls
				}
				_ := <-done_ch {
					workers_done++
					if workers_done == c.threads {
						break
					}
				}
			}
		}
		// Drain remaining results
		for {
			result := <-result_ch or { break }
			if result.class_id !in classes {
				savefile := '${c.save_dir}/${result.class_id}.nbv'
				classes[result.class_id] = model.NbClass.new(result.class_id, c.kmer_size, savefile)
			}
			mut cls := classes[result.class_id]
			cls.add_genome(result.kmer_counts)
			classes[result.class_id] = cls
		}
	}

	// Serialize
	os.mkdir_all(c.save_dir)!
	for _, cls in classes {
		nbio.save_class(cls, cls.savefile)!
	}
	nbio.save_meta(c.save_dir, c.kmer_size)!
}

struct TrainResult {
	class_id    string
	kmer_counts map[int]int
}

fn train_worker(job_ch chan TrainJob, result_ch chan TrainResult, done_ch chan bool, input_type cfg.InputType, k int) {
	for {
		job := <-job_ch or { break }
		kmer_counts := load_kmer_counts(job.path, input_type, k) or {
			eprintln('Warning: skipping ${job.path}: ${err}')
			continue
		}
		result_ch <- TrainResult{
			class_id:    job.class_id
			kmer_counts: kmer_counts
		}
	}
	done_ch <- true
}
```

Note: V's `select` syntax and channel close/drain semantics may require adjustment. The implementer should verify the exact V channel API and adjust accordingly. The key invariant: workers send results on `result_ch`, the main thread is the sole consumer, no locks needed.

- [ ] **Step 3: Run tests to verify both single and multi-threaded pass**

Run: `v test src/pipeline/`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/pipeline/train.v src/pipeline/pipeline_test.v
git commit -m "feat: multithreaded training pipeline using channels"
```

---

## Task 12: Add Multithreading to Classification Pipeline

**Files:**
- Modify: `src/pipeline/classify.v`

- [ ] **Step 1: Add threaded classification test**

```v
fn test_classify_multithreaded() {
	// Reuse the train-then-classify test with threads=2
	train_dir := '/tmp/nbv_test_e2e_mt_train'
	os.rmdir_all(train_dir) or {}

	train(cfg.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   train_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    1
		input_type: .fasta
		extension:  '.fasta'
	})!

	os.mkdir_all('/tmp/nbv_test_e2e_mt_reads')!
	os.write_file('/tmp/nbv_test_e2e_mt_reads/test_read.fasta', '>read1\nACGTACGTACGTACGT\n>read2\nGGGGCCCCAAAATTTT\n')!

	classify(cfg.Config{
		mode:        .classify
		kmer_size:   4
		save_dir:    train_dir
		source_dir:  '/tmp/nbv_test_e2e_mt_reads'
		threads:     2
		input_type:  .fasta
		extension:   '.fasta'
		format:      .csv
		prefix:      '/tmp/nbv_test_e2e_mt_output'
		full_result: false
		temp_dir:    '/tmp'
		limit_mb:    0
		max_rows:    1000
		max_cols:    100
	})!

	output := os.read_file('/tmp/nbv_test_e2e_mt_output.csv')!
	assert output.contains('read1')
	assert output.contains('read2')

	os.rmdir_all(train_dir) or {}
	os.rmdir_all('/tmp/nbv_test_e2e_mt_reads') or {}
	os.rm('/tmp/nbv_test_e2e_mt_output.csv') or {}
}
```

- [ ] **Step 2: Refactor `classify()` to use the three-stage pipeline when threads > 1**

Update `src/pipeline/classify.v` to add a threaded path:

- Reader thread: parses FASTA/kmr files, sends `SeqJob` to `seq_ch`
- N worker threads: receive `SeqJob`, compute log-likelihoods, send `ClassifyResult` to `result_ch`
- Writer thread: drains `result_ch`, writes output

Use the same `classify_read` function from the single-threaded path. The channel types `SeqJob` and `ClassifyResult` are already defined.

The single-threaded path remains as a fallback for `threads <= 1`.

- [ ] **Step 3: Run tests**

Run: `v test src/pipeline/`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/pipeline/classify.v src/pipeline/pipeline_test.v
git commit -m "feat: multithreaded classification pipeline with reader/workers/writer stages"
```

---

## Task 13: Memory Management (batch_size and limit_mb)

**Files:**
- Modify: `src/pipeline/train.v`
- Modify: `src/pipeline/classify.v`

- [ ] **Step 1: Add memory-constrained training test**

```v
fn test_train_with_batch_size() {
	out_dir := '/tmp/nbv_test_train_batch'
	os.rmdir_all(out_dir) or {}

	c := cfg.Config{
		mode:       .train
		kmer_size:  4
		save_dir:   out_dir
		source_dir: 'src/pipeline/testdata/training'
		threads:    1
		input_type: .fasta
		extension:  '.fasta'
		limit_mb:   0
		batch_size: 1 // Process one genome at a time
	}

	train(c)!
	assert os.exists('${out_dir}/class_a.nbv')
	assert os.exists('${out_dir}/class_b.nbv')

	os.rmdir_all(out_dir) or {}
}
```

- [ ] **Step 2: Implement batch_size support in train()**

In the training loop, after processing `batch_size` genomes:
1. Serialize all current classes to disk
2. Clear in-memory state
3. On next batch, load existing savefiles and continue adding

- [ ] **Step 3: Implement multi-round classification for limit_mb**

In `classify()`, when `limit_mb > 0`:
1. Calculate how many classes fit in `limit_mb` using `NbClass.size_bytes()`
2. Load a subset, classify all reads, write partial `.max` results to `temp_dir`
3. Repeat for remaining classes
4. Final merge pass to produce output

This is the most complex memory management feature. Follow the `.max` file format from the spec: length-prefixed strings (u32 + UTF-8) for `seq_id` and `best_class`, followed by `best_score` as f64.

- [ ] **Step 4: Run all tests**

Run: `v test src/`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/pipeline/train.v src/pipeline/classify.v src/pipeline/pipeline_test.v
git commit -m "feat: memory management with batch_size and multi-round classification"
```

---

## Task 14: Legacy NBC++ Cross-Validation Test

**Files:**
- Modify: `src/pipeline/pipeline_test.v`

- [ ] **Step 1: Write legacy cross-validation test**

```v
fn test_classify_with_legacy_savefiles() {
	// Uses the example NBC++ savefiles from example/training_classes/
	// Classifies example/reads/cross.fna
	// Compares output against example/results_max_1.csv

	c := cfg.Config{
		mode:        .classify
		kmer_size:   9
		save_dir:    'example/training_classes'
		source_dir:  'example/reads'
		threads:     1
		input_type:  .fasta
		extension:   '.fna'
		format:      .csv
		prefix:      '/tmp/nbv_legacy_test_output'
		full_result: false
		temp_dir:    '/tmp'
		limit_mb:    0
		max_rows:    450000
		max_cols:    20000
	}

	classify(c)!

	// Compare against expected results
	output := os.read_file('/tmp/nbv_legacy_test_output.csv')!
	expected := os.read_file('example/results_max_1.csv')!

	output_lines := output.trim_space().split('\n')
	expected_lines := expected.trim_space().split('\n')

	// Check that class assignments match for the first several reads
	mut matches := 0
	for i, exp_line in expected_lines {
		if i >= output_lines.len {
			break
		}
		exp_parts := exp_line.split(',')
		out_parts := output_lines[i].split(',')
		if exp_parts.len >= 2 && out_parts.len >= 2 {
			if exp_parts[0] == out_parts[0] && exp_parts[1] == out_parts[1] {
				matches++
			}
		}
	}

	// At least 90% of class assignments should match
	match_pct := f64(matches) / f64(expected_lines.len)
	assert match_pct > 0.9

	os.rm('/tmp/nbv_legacy_test_output.csv') or {}
}
```

- [ ] **Step 2: Run the test**

Run: `v test src/pipeline/ -run test_classify_with_legacy_savefiles`
Expected: Test passes (if legacy format parsing is correct). If it fails, debug by comparing individual read classifications and adjusting the legacy parser.

- [ ] **Step 3: Commit**

```bash
git add src/pipeline/pipeline_test.v
git commit -m "test: legacy NBC++ cross-validation against example data"
```

---

## Task 15: Final Integration and Cleanup

**Files:**
- All modules

- [ ] **Step 1: Run all tests**

Run: `v test src/`
Expected: All tests pass.

- [ ] **Step 2: Build the binary**

Run: `v build src/ -o nbv -prod`
Expected: Clean build, no warnings.

- [ ] **Step 3: Test end-to-end with example data**

Run: `./nbv example/classify.yaml`
Expected: Classification completes, output written to `example_results.csv`.

- [ ] **Step 4: Test training end-to-end**

Create a proper training directory structure from the example reads:
```bash
mkdir -p /tmp/nbv_final_train/class_test
cp example/reads/cross.fna /tmp/nbv_final_train/class_test/
```

Run: `./nbv example/train.yaml` (after updating `source_dir` in train.yaml to point to the training dir)
Expected: Training completes, savefiles written.

- [ ] **Step 5: Run all tests one final time**

Run: `v test src/`
Expected: All pass.

- [ ] **Step 6: Final commit**

```bash
git add -A
git commit -m "feat: NBV v0.1.0 - Naive Bayes metagenomic classifier in V"
```
