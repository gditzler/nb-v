# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

NBV is a Naive Bayes classifier for metagenomic sequence classification using k-mer frequencies, written in V (vlang). It is a reimplementation of [NBC++](https://github.com/EESI/Naive_Bayes) with improved architecture.

## Build & Test Commands

```bash
v build src/ -o nbv          # Build binary
v test src/                   # Run all tests
v test src/kmer/              # Test single module
./nbv example/classify.yaml  # Run classification with example data
./nbv example/train.yaml     # Run training with example data
```

## Architecture

Module-per-concern design under `src/`. No circular dependencies.

```
main.v → config/ → prantlf.yaml (external)
              ↑
         pipeline/ → model/ → kmer/
              ↓
             io/ → kmer/
```

- **config/** — YAML config parsing into flat `Config` struct. All validation at load time.
- **kmer/** — Pure functions for base-4 encoding, reverse complement, canonical k-mers, buffer counting. Thread-safe (no shared state).
- **model/** — `NbClass` (Naive Bayes per-class model), `KahanAccumulator` for numerically stable log-likelihood. Laplace smoothing with `sumfreq` initialized to `num_canonical_kmers(k)`.
- **io/** — FASTA parsing (`read_fasta` returns `[]FastaRecord`), `.kmr` file reading, NBV binary serialization, legacy NBC++ savefile reader, CSV/TSV/JSON output writer.
- **pipeline/** — Train and classify orchestrators with channel-based concurrency (`spawn`/`chan`).

## V-Specific Patterns

- **Imports:** Sibling modules require `import src.module_name` (e.g., `import src.kmer as kmod`), not `import module_name`.
- **YAML:** `prantlf.yaml` uses `yaml.unmarshal_file[T](path)` (one arg). The two-arg variant is `unmarshal_file_opt`.
- **Mutable closures:** V 0.5.0 closure captures with `fn [mut var] (...)` don't propagate mutations back. Use return values or `read_fasta()` instead of callback+capture pattern.
- **Error handling:** Result types `!T` with `or { }` blocks. No `try/catch`.
- **Map access:** Missing keys return zero value (no option type), which is relied on for `get_freq_count_lg` returning 0.0 for unseen k-mers.

## Key Design Decisions

- Max k-mer size: k <= 15 (32-bit int encoding)
- Prior (`ngenomes_lg`) is NOT included in classification score (matches NBC++)
- Log-space parameters recomputed eagerly after each `add_genome` (simpler than C++ lazy/dirty-flag approach)
- Legacy NBC++ save format: `f64(ngenomes_lg) + f64(sumfreq_lg) + i32(n_entries) + n_entries*(i32 kmer, f64 freqcnt_lg)`

## Documentation

- **Design spec:** `docs/superpowers/specs/2026-03-18-nbv-design.md`
- **Implementation plan:** `docs/superpowers/plans/2026-03-18-nbv-implementation.md`
- **Example data:** `example/training_classes/` (100 NBC++ savefiles, k=9), `example/reads/cross.fna`, `example/results_max_1.csv`
