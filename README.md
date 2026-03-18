# NBV

A Naive Bayes classifier for metagenomic sequence classification using k-mer frequencies, written in [V](https://vlang.io/).

NBV reimplements [NBC++](https://github.com/EESI/Naive_Bayes) with a cleaner module-per-concern architecture, multithreaded train/classify pipelines, and a single YAML configuration file instead of command-line flags.

## Features

- **Train** models from FASTA or pre-computed k-mer files
- **Classify** metagenomic reads against trained models
- **Multithreaded** training and classification using V's channels
- **Memory management** with configurable limits (`limit_mb`, `batch_size`, `max_rows`, `max_cols`)
- **Multiple output formats:** CSV, TSV, JSON Lines
- **Legacy compatibility:** reads NBC++ `-save.dat` savefiles for cross-validation
- **Numerically stable:** Kahan summation for log-likelihood accumulation
- **Pure V** with no external dependencies beyond YAML parsing

## Requirements

- [V](https://vlang.io/) 0.5.0+
- `prantlf.yaml` module: `v install prantlf.yaml`

## Building

```bash
v build src/ -o nbv
```

For an optimized build:

```bash
v -prod -o nbv src/
```

## Usage

NBV takes a single argument: a YAML configuration file.

```bash
./nbv config.yaml
```

### Configuration

```yaml
version: 1

mode: classify           # "train" or "classify"
kmer_size: 9
save_dir: ./model_dir    # where trained models are saved/loaded
source_dir: ./reads      # training: dir of class subdirs; classify: dir of input files
threads: 4

input:
  extension: .fna        # file extension to look for
  input_type: fasta      # "fasta" or "kmer_file"

memory:
  limit_mb: 0            # 0 = unlimited; >0 = multi-round classification
  batch_size: 0          # 0 = all at once; >0 = flush to disk every N genomes
  max_rows: 0            # 0 = unlimited; >0 = stop after N reads
  max_cols: 0            # 0 = unlimited; >0 = load at most N classes

output:
  format: csv            # "csv", "tsv", or "json"
  prefix: results        # output filename prefix (e.g., results.csv)
  full_result: false     # true = log-likelihoods for all classes per read
  temp_dir: /tmp
```

### Training

Organize training data as class subdirectories, each containing genome files:

```
training_data/
  species_a/
    genome1.fasta
    genome2.fasta
  species_b/
    genome3.fasta
```

```yaml
mode: train
source_dir: ./training_data
save_dir: ./trained_model
kmer_size: 9
input:
  extension: .fasta
  input_type: fasta
```

```bash
./nbv train_config.yaml
# Output: trained_model/species_a.nbv, trained_model/species_b.nbv, trained_model/meta.nbv
```

### Classification

```yaml
mode: classify
source_dir: ./reads_to_classify
save_dir: ./trained_model
kmer_size: 9
input:
  extension: .fasta
  input_type: fasta
output:
  format: csv
  prefix: classification_results
```

```bash
./nbv classify_config.yaml
# Output: classification_results.csv
```

Output format (CSV):
```
read_id,best_class,log_likelihood
NZ_CP031447.1_0_0/1,370777,-1281.47
NZ_CP031447.1_1_0/1,1748,-1274.63
```

### Legacy NBC++ Savefiles

NBV can classify against existing NBC++ trained models (`.dat` savefiles) without retraining:

```yaml
mode: classify
save_dir: ./nbc_savefiles   # directory containing *-save.dat files
```

## Example

The `example/` directory contains NBC++ training data and reads for testing:

```bash
# Classify example reads against 100 pre-trained NBC++ models (k=9)
./nbv example/classify.yaml

# Compare output against expected results
diff <(cut -d, -f1,2 example_results.csv | sort) \
     <(cut -d, -f1,2 example/results_max_1.csv | sort)
```

## Testing

```bash
v test src/           # all tests
v test src/kmer/      # single module
```

## Architecture

```
src/
  main.v              # CLI entry point
  config/config.v     # YAML config parsing and validation
  kmer/kmer.v         # k-mer encoding, reverse complement, canonical form, counting
  model/model.v       # Naive Bayes model, Kahan summation, Laplace smoothing
  io/
    fasta.v           # FASTA parser
    kmer_file.v       # NBC++ .kmr file reader
    serialization.v   # NBV binary format + legacy NBC++ reader
    writer.v          # CSV/TSV/JSON output
  pipeline/
    train.v           # Training orchestrator (single/multithreaded)
    classify.v        # Classification orchestrator (single/multithreaded/multi-round)
```

## References

- [NBC++ (EESI/Naive_Bayes)](https://github.com/EESI/Naive_Bayes) -- the original C++ implementation
- Rosen, G., Garbarine, E., Caseiro, D., Polikar, R., & Sokhansanj, B. (2008). Metagenome fragment classification using N-mer frequency profiles. *Advances in Bioinformatics*.

## License

MIT License. See [LICENSE](LICENSE).
