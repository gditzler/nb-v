# NBV: Naive Bayes Classifier for Metagenomic Data in V

**Date:** 2026-03-18
**Status:** Approved
**Source reference:** [EESI/Naive_Bayes (NBC++)](https://github.com/EESI/Naive_Bayes)

## Overview

NBV is a reimplementation of NBC++ in V (vlang). It is a Naive Bayes classifier for metagenomic sequence classification using k-mer frequencies. It supports training models from genomic data and classifying unknown reads against trained models.

This is not a line-by-line port. The goal is to leverage V's strengths (channels, modules, simplicity) and fix architectural weaknesses in the original C++ codebase, while preserving algorithmic correctness.

## Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Modes | Train + Classify (full pipeline) | User needs complete workflow |
| Concurrency | Multithreaded from day one via V's `spawn`/channels | Metagenomic datasets are large; single-threaded is too slow |
| External deps | Pure V + `yaml` vpm module | Grinder is test-only; Jellyfish replaced by internal k-mer counting; Boost replaced by V stdlib + small Kahan impl. The only external dep is the `yaml` vpm module for config parsing. |
| CLI interface | Single YAML config file argument | User preference; replaces 12+ CLI flags |
| Input formats | Raw FASTA or pre-computed k-mer files | Flexible; config-specified via `input_type` |
| Serialization | New binary format (primary) + legacy NBC++ reader (verification) | Clean design, but cross-validation with C++ version is valuable |
| Memory management | Fine-grained controls (limit_mb, batch_size, max_rows, max_cols) | Necessary for large metagenomic databases on constrained machines |
| Output formats | CSV, TSV, JSON Lines (configurable) | Flexibility for downstream tools |
| Numerical precision | Kahan summation for log-likelihood accumulation | Prevents floating-point drift on reads with thousands of k-mers |
| Max k-mer size | k <= 15 (signed 32-bit int encoding) | Matches NBC++ limitation. Base-4 encoding of a k-mer requires 2k bits; `int` (i32) supports up to k=15. Use `i64` in the future if larger k needed. |

## YAML Configuration

The program accepts a single CLI argument: `./nbv config.yaml`

```yaml
version: 1             # config format version

mode: train            # or "classify"
kmer_size: 6
save_dir: ./nbv_save
source_dir: ./data/training_genomes  # for train: dir of class subdirs; for classify: dir or file of reads
threads: 4

input:
  extension: .kmr        # or .fasta for raw sequences
  input_type: kmer_file  # or "fasta" (triggers internal counting)

memory:
  limit_mb: 16000        # 0 = unlimited
  batch_size: 0          # genomes per batch, 0 = all
  max_rows: 450000
  max_cols: 20000

output:
  format: csv            # csv, tsv, or json
  prefix: log_likelihood
  full_result: false     # true = all class scores per read
  temp_dir: /tmp
```

**Field semantics:**

- `source_dir` is dual-purpose: in train mode, it points to a directory of class subdirectories (see Training Directory Structure below). In classify mode, it points to the directory containing the input file(s) to classify.
- `version` allows forward-compatible config parsing.

Validation happens at startup with clear error messages for missing or invalid fields. Functions that can fail return V's `!T` (result type) for structured error handling.

## Training Directory Structure

The expected layout for `source_dir` in train mode:

```
source_dir/
  class_name_1/
    genome1.kmr       # or genome1.fasta if input_type is "fasta"
    genome2.kmr
  class_name_2/
    genome3.kmr
    genome4.kmr
```

The subdirectory name becomes the class ID. Each file within a subdirectory is one genome belonging to that class.

## K-mer File Format (.kmr)

Pre-computed `.kmr` files use the NBC++ text format:

```
ACGTAA	15
CGTAAC	8
TGCAAT	3
```

Each line contains a k-mer DNA string, a tab separator, and an integer count. The reader encodes the k-mer string to its integer representation using base-4 encoding (A=0, C=1, G=2, T=3) and canonicalizes it (takes the lexicographically smaller of the k-mer and its reverse complement).

## Architecture

Module-per-concern design with six modules:

```
nbv/
├── main.v
├── config/
│   └── config.v
├── kmer/
│   └── kmer.v
├── model/
│   └── model.v
├── io/
│   ├── fasta.v
│   ├── kmer_file.v
│   ├── serialization.v
│   └── writer.v
├── pipeline/
│   ├── train.v
│   └── classify.v
└── tests/
    ├── kmer_test.v
    ├── model_test.v
    ├── io_test.v
    ├── pipeline_test.v
    └── legacy_test.v
```

### config/

Parses the YAML configuration file into a `Config` struct. Validates all required fields at startup and errors early with clear messages. Replaces Boost `program_options`.

### kmer/

All k-mer logic, extracted from C++ `Diskutil::countKmer()`:

- `encode(kmer []u8, k int) int` -- encodes a single k-mer-length subsequence to its base-4 integer representation (A=0, C=1, G=2, T=3). Operates on exactly `k` bytes, not an entire buffer.
- `reverse_complement(kmer_int int, k int) int` -- computes reverse complement in integer space (no string allocation)
- `canonical(kmer_int int, k int) int` -- returns the lexicographically smaller of a k-mer and its reverse complement. This works because the base-4 encoding preserves lexicographic order (A(0) < C(1) < G(2) < T(3) matches alphabetical order).
- `count_from_buffer(buf []u8, k int) map[int]int` -- slides a window of size `k` across the buffer, calls `encode` + `canonical` on each window, counts occurrences. When a newline or invalid character (any character not in {A, C, G, T}) is encountered, the current k-mer window is invalidated and a new window starts from the next valid character.
- `num_canonical_kmers(k int) i64` -- computes the number of distinct canonical k-mers for a given k. Formula: for odd k, `4^k / 2`; for even k, `(4^k + 4^(k/2)) / 2` (accounts for palindromic reverse complements).

All functions are pure (no shared state), so workers can call them concurrently without synchronization.

**Improvement over C++:** K-mer logic was buried inside `Diskutil` alongside unrelated filesystem code. Isolation enables independent unit testing and potential future optimization.

### model/

The Naive Bayes math and per-class state:

```
KahanAccumulator {
    sum   f64
    comp  f64    // running compensation
}

LoadState = .unloaded | .full | .classify_only

NbClass {
    id          string
    kmer_size   int
    savefile    string    // path used by io/ for save/load; set at construction or load time
    // Log-space parameters (classification)
    ngenomes_lg  f64
    sumfreq_lg   f64
    freqcnt_lg   map[int]f64
    // Plain parameters (training updates)
    ngenomes     int
    sumfreq      i64
    freqcnt      map[int]int
    state        LoadState  // .unloaded, .full (both plain+log), .classify_only (log-space only)
}
```

#### Initialization

When creating a new `NbClass` with kmer_size `k`:

- `ngenomes = 0`
- `sumfreq` is initialized to the number of canonical k-mers for size `k` (via `kmer.num_canonical_kmers(k)`, call this value `V`). This serves as the Laplace pseudo-count base. After training, `sumfreq = V + total_training_counts`, which gives the standard Laplace smoothing denominator: `P(kmer | class) = (count + 1) / (V + total_training_counts)`.
- `freqcnt` starts empty

#### Training: `add_genome(kmer_counts map[int]int)`

1. Increment `ngenomes` by 1
2. For each `(kmer, count)` in `kmer_counts`: add `count` to `freqcnt[kmer]` (insert with value `count` if absent)
3. Add the sum of all counts in `kmer_counts` to `sumfreq`
4. Eagerly recompute log-space parameters:
   - `ngenomes_lg = log(ngenomes)`
   - `sumfreq_lg = log(sumfreq)`
   - For each kmer in `freqcnt`: `freqcnt_lg[kmer] = log(freqcnt[kmer] + 1)` (the +1 is Laplace smoothing)

Note: The C++ version uses lazy recomputation via dirty flags (`double_wflag`). We use eager recomputation for simplicity. This is equivalent in correctness -- the log values are always consistent after `add_genome` returns.

#### Laplace Smoothing

The smoothed log-probability of a k-mer given a class is:

```
log P(kmer | class) = log(count(kmer) + 1) - log(sumfreq)
```

Where:
- `count(kmer)` is the raw count (0 if unseen)
- `sumfreq` includes the canonical k-mer count as its initial value, so the denominator naturally accounts for the pseudo-counts

For unseen k-mers: `get_freq_count_lg(kmer)` returns `log(1) = 0` (just the Laplace pseudo-count in the numerator), and the normalization comes from `sumfreq_lg`.

#### Classification Log-Likelihood Formula

`compute_log_likelihood(kmer_counts map[int]int) f64` computes:

```
score = sum_kahan( freq_i * freqcount_lg(kmer_i) ) - total_kmer_count * sumfreq_lg
```

Where:
- `freq_i` is the count of k-mer `i` in the read
- `freqcount_lg(kmer_i)` is `log(count(kmer_i) + 1)` (from `freqcnt_lg`, or `log(1)` if unseen)
- `total_kmer_count` is the sum of all `freq_i` values
- `sumfreq_lg` is `log(sumfreq)` for the class
- The summation uses Kahan accumulation for numerical stability

The prior (`ngenomes_lg`) is **not** included in the score, matching NBC++ behavior (the C++ code comments out the prior term).

#### Other functions

- `NbClass.get_freq_count_lg(kmer int) f64` -- returns `freqcnt_lg[kmer]` if present, else `0.0` (log of Laplace pseudo-count of 1)
- `NbClass.size_bytes() u64` -- reports memory footprint: `sizeof(NbClass) + len(freqcnt) * (sizeof(int) + sizeof(int)) + len(freqcnt_lg) * (sizeof(int) + sizeof(f64))`
- `kahan_add(acc KahanAccumulator, val f64) KahanAccumulator` -- numerically stable addition

**Improvement over C++:** The C++ `Class<T>` mixes model logic with genome queuing, load/unload lifecycle, and disk I/O. Here, `NbClass` is purely the mathematical model. Serialization lives in `io/`, lifecycle management lives in `pipeline/`.

### io/

All file reading, writing, and serialization:

**FASTA Parser:**
- `parse_fasta(path string, callback fn(header string, sequence []u8))` -- streaming parser, callback per sequence, never loads the entire file. Accumulates multi-line sequences before invoking the callback. For classification (short reads), this is memory-efficient. For training on large genomes, the full sequence is held in memory temporarily during counting; this is acceptable because k-mer counting is fast and the sequence is released after the callback returns.
- `count_sequences(path string) u64` -- fast header scan counting `>` lines for memory allocation planning

**K-mer File Reader:**
- `read_kmer_file(path string, k int) !map[int]int` -- reads `.kmr` text files (see K-mer File Format above). Parses each line, encodes the k-mer string via `kmer.encode()`, canonicalizes via `kmer.canonical()`, and builds the count map. Returns error on malformed lines.
- `detect_input_type(extension string) InputType` -- determines FASTA vs k-mer file based on the configured extension. `.fasta`, `.fa`, `.fna` map to FASTA; `.kmr` maps to k-mer file. No magic-byte sniffing; relies on the config-specified extension.

**Serialization (new format):**
- `save_class(cls NbClass, path string) !` -- binary format: header (magic bytes `NBV1`, version u8, kmer_size i32, ngenomes i32, sumfreq i64) + packed `(i32, i32)` frequency pairs (kmer_int, count). Stores plain-space data for continued training.
- `load_class(path string) !NbClass` -- reads new format, populates both plain and log-space fields
- `load_legacy_class(path string, k int) !NbClass` -- reads NBC++ `-save.dat` files. These store log-space values. Populates `freqcnt_lg`, `ngenomes_lg`, `sumfreq_lg`. Plain-space fields (`freqcnt`, `ngenomes`, `sumfreq`) are **not** populated (set to 0/empty) and `state` is set to `.classify_only`. This is sufficient for cross-validation of classification results.

**Save file naming convention:**
- New format: `<save_dir>/<class_id>.nbv` (one file per class)
- Legacy format: `<save_dir>/<class_id>-save.dat` (NBC++ convention)
- `meta.nbv`: single-line text file containing the integer k-mer size (e.g., `6`). Written during training. On classify, parsed and checked against the config's `kmer_size`; an error is raised on mismatch.

**Output Formatter:**
- `Writer` struct with buffered file handle and `format` field (csv/tsv/json)
- `Writer.write_header(class_ids []string)` -- column headers for full_result mode
- `Writer.write_result(seq_id string, best_class string, score f64)` -- top-hit result
- `Writer.write_full_result(seq_id string, scores map[string]f64)` -- all class scores
- JSON mode: one JSON object per line (JSON Lines)
- Output filename: `<prefix>.<ext>` where `ext` is `csv`, `tsv`, or `jsonl` based on the configured format (e.g., `log_likelihood.csv`)

**Improvement over C++:** The C++ version scatters I/O across `Diskutil`, `Class::serialize/deserialize`, `NB::writeToCSV`, `NB::concatenateCSVByColumns`, and memory-mapped reads. Here, all I/O goes through one module with a consistent interface.

### pipeline/

Concurrency orchestration using V's `spawn` and `chan` primitives:

**Training Pipeline (`train`):**

1. Scan `source_dir` for class subdirectories containing k-mer or FASTA files
2. Build work queue of `(class_id, file_path)` tuples
3. Spawn `n` worker threads pulling from an input channel (buffered, capacity = `2 * threads`):
   - Read/count k-mers per genome
   - Send `(class_id, kmer_counts)` to accumulator channel
4. **One** accumulator thread receives all `(class_id, kmer_counts)` messages and calls `NbClass.add_genome()` on the appropriate class. No locks needed because only this single thread mutates any `NbClass`.
5. Respect `batch_size`: flush and serialize when batch limit reached
6. Respect `limit_mb`: track memory, serialize and unload when approaching cap
7. After all genomes processed, serialize all classes to `save_dir` and write `meta.nbv`

**FASTA-mode training:** When `input_type` is `fasta`, each FASTA file is treated as one genome. The worker thread calls `io.parse_fasta()` on the file, and for each sequence in the file, counts k-mers via `kmer.count_from_buffer()`. If the file contains multiple sequences, their k-mer counts are merged (summed) into a single `map[int]int` for the genome. This merged map is then sent to the accumulator as one genome.

**Classification Pipeline (`classify`):**

1. Load class savefiles from `save_dir`, respecting `limit_mb` (chunk if needed)
2. Open output `Writer`
3. Three-stage pipeline with typed channels:
   - **Sequence channel** `chan SeqJob` where `SeqJob { seq_id string; kmer_counts map[int]int }` -- buffered, capacity = `2 * threads`
   - **Result channel** `chan ClassifyResult` where `ClassifyResult { seq_id string; best_class string; best_score f64; all_scores map[string]f64 }` -- buffered, capacity = `2 * threads`
   - **Reader thread**: streams FASTA/k-mer input, counts k-mers if needed, feeds sequence channel
   - **N worker threads**: receive `SeqJob`, compute log-likelihood against all loaded classes, send `ClassifyResult` to result channel
   - **Writer thread**: drains result channel, writes results via `io.Writer`
4. **Multi-round classification** (when not all classes fit in memory):
   - Load a subset of classes that fits within `limit_mb`
   - Score all reads against the loaded subset
   - For each read, write `(seq_id, best_class, best_score)` to a temporary `.max` file on disk (in `temp_dir`). Format: length-prefixed strings (u32 length + UTF-8 bytes) for `seq_id` and `best_class`, followed by `best_score` as f64. One record per read, sequential.
   - Unload current classes, load next subset
   - Re-read the input and score again; for each read, compare the new best score against the stored best in the `.max` file, keeping whichever is higher
   - Repeat until all classes have been scored
   - Final pass: read the `.max` file and write final output via `io.Writer`
   - For `full_result` mode: each round appends per-class scores to temporary column files; a final merge pass concatenates them
5. Respect `max_rows`/`max_cols` for output dimensioning

**Improvement over C++:** The C++ version uses raw pthreads, manual mutexes, and condition variables scattered across `NB.cpp` and `Genome.hpp`. V's channels provide the same producer-consumer pattern with less boilerplate and no deadlock risk from lock ordering. The reader/workers/writer form a clean three-stage pipeline rather than interleaving buffer management with thread coordination.

### main.v

Entry point. Parses CLI argument, reads and validates YAML config, dispatches to `pipeline.train()` or `pipeline.classify()`. Minimal -- all logic lives in modules.

## Error Handling

V's `!T` result type is used for all operations that can fail. The strategy by layer:

- **config/**: Validation errors at startup cause an immediate exit with a descriptive message (file not found, missing required field, invalid value).
- **io/**: File I/O errors, malformed FASTA, corrupt save files return `!T`. Callers in `pipeline/` handle these -- log the error and skip the file for training (with a warning), or abort for classification (partial results are unreliable).
- **model/**: Pure math, does not fail. Only `size_bytes` and accessors.
- **pipeline/**: Catches `io` errors. For training, a single bad genome file logs a warning and continues. For classification, a corrupt save file or unreadable input file is fatal. Individual malformed reads mid-stream during classification are skipped with a warning (logged to stderr) and do not abort the pipeline.

## Testing Strategy

- **Unit tests (`kmer_test.v`):** encoding with known sequences, canonical k-mer verification, count correctness against hand-counted examples, `num_canonical_kmers` against known values
- **Unit tests (`model_test.v`):** log-likelihood against hand-computed values, Kahan accumulator accuracy vs naive summation on adversarial inputs, `add_genome` verifying ngenomes/sumfreq/freqcnt increments
- **Unit tests (`io_test.v`):** FASTA parsing edge cases (multi-line sequences, empty headers, invalid characters), serialization round-trip, output format correctness, `.kmr` file parsing
- **Integration test (`pipeline_test.v`):** train on small synthetic dataset (few short sequences across 2-3 classes), classify held-out reads, verify top prediction matches expected class
- **Legacy test (`legacy_test.v`):** load NBC++ savefile, load same data through new format, assert log-likelihoods match within floating-point tolerance

## Dependencies

Pure V with one external module: `yaml` from vpm for config parsing. All other functionality (k-mer counting, Naive Bayes math, Kahan summation, serialization, concurrency) is implemented directly.
