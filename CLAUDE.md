# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

NBV is a Naive Bayes classifier for metagenomic sequence classification using k-mer frequencies, written in V (vlang). It reimplements [NBC++](https://github.com/EESI/Naive_Bayes) with a module-per-concern architecture, multithreaded pipelines, and YAML-based configuration.

## Build & Test Commands

```bash
v build src/ -o nbv           # Build binary
v -prod -o nbv src/            # Production build (optimized)
v test src/                    # Run all tests (~3 min, dominated by legacy cross-validation)
v test src/kmer/               # Test single module
v test src/model/              # Test model module
v test src/io/                 # Test I/O module
v test src/pipeline/           # Test pipelines (~3 min, includes legacy cross-validation)
./nbv example/classify.yaml   # Classify reads against example data
./nbv example/train.yaml      # Train on example data
```

Note: Pipeline tests take ~3 minutes because `test_classify_with_legacy_savefiles` loads 100 NBC++ class models and classifies 94K reads.

## Architecture

Module-per-concern design under `src/`. No circular dependencies.

```
main.v → config/ → prantlf.yaml (external)
              ↑
         pipeline/ → model/ → kmer/
              ↓
             io/ → kmer/
```

- **config/** — YAML config parsing into flat `Config` struct. Enums: `Mode`, `InputType`, `OutputFormat`. All validation at load time via `config.load(path)`.
- **kmer/** — Pure functions for base-4 encoding, reverse complement, canonical k-mers, buffer counting. Thread-safe (no shared state). Max k-mer size: k <= 15 (32-bit int encoding).
- **model/** — `NbClass` (per-class Naive Bayes model), `KahanAccumulator` for numerically stable log-likelihood. Laplace smoothing with `sumfreq` initialized to `num_canonical_kmers(k)`.
- **io/** — `FastaRecord` + `read_fasta()`, `.kmr` file reading, NBV binary serialization (`save_class`/`load_class`), legacy NBC++ reader (`load_legacy_class`), output `Writer` (CSV/TSV/JSON Lines).
- **pipeline/** — `train()` and `classify()` orchestrators. Single-threaded and multithreaded paths (channels + `sync.WaitGroup`). Memory management: `batch_size` for training, `limit_mb` for multi-round classification, `max_rows`/`max_cols` limits.

## V-Specific Patterns

- **Imports:** Sibling modules require `import src.module_name` (e.g., `import src.kmer as kmod`), not `import module_name`.
- **YAML:** `prantlf.yaml` uses `yaml.unmarshal_file[T](path)` (one arg). The two-arg variant is `unmarshal_file_opt`.
- **Mutable closures:** V 0.5.0 closure captures with `fn [mut var] (...)` don't propagate mutations back. Use return values or `read_fasta()` instead of callback+capture pattern.
- **Error handling:** Result types `!T` with `or { }` blocks. No `try/catch`.
- **Map access:** Missing keys return zero value (no option type), which is relied on for `get_freq_count_lg` returning 0.0 for unseen k-mers.
- **Concurrency:** Use `sync.WaitGroup` for worker completion tracking (V 0.5.0 thread `.wait()` panics). Buffered channels with `cap:` sized to expected items.
- **Unsafe:** `bytes_to_f64` in `io/serialization.v` uses `unsafe { *(&f64(&bits)) }` for bit-level reinterpretation.

## Key Design Decisions

- Prior (`ngenomes_lg`) is NOT included in classification score (matches NBC++)
- Log-space parameters recomputed eagerly after each `add_genome` (simpler than C++ lazy/dirty-flag approach)
- Legacy NBC++ save format: `f64(ngenomes_lg) + f64(sumfreq_lg) + i32(n_entries) + n_entries*(i32 kmer, f64 freqcnt_lg)` — verified via hex dump
- Multithreaded classification output may differ in order from single-threaded (content is identical)

## Documentation

- **Design spec:** `docs/superpowers/specs/2026-03-18-nbv-design.md`
- **Implementation plan:** `docs/superpowers/plans/2026-03-18-nbv-implementation.md`
- **Example data:** `example/training_classes/` (100 NBC++ savefiles, k=9), `example/reads/cross.fna`, `example/results_max_1.csv`


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:1105d646 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/core-concepts/sync-concepts.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->
