// config parses and validates YAML configuration files into a flat Config struct.
// All validation is performed at load time; callers receive either a ready-to-use
// Config or an error describing the first constraint violation found.
module config

import prantlf.yaml

// Mode specifies whether the run will train new class models or classify reads.
pub enum Mode {
	train
	classify
}

// InputType specifies the format of the input sequence data.
// kmer_file expects pre-computed k-mer count files (.kmr); fasta expects raw FASTA files.
pub enum InputType {
	kmer_file
	fasta
}

// OutputFormat specifies the format of the classification result file.
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

// Config holds the fully parsed and validated configuration for a single NBV run.
// All string-typed enum fields from the YAML source are resolved to their typed enum
// equivalents; nested input/memory/output sections are flattened into a single struct.
pub struct Config {
pub:
	version     int
	mode        Mode
	kmer_size   int
	save_dir    string
	source_dir  string
	threads     int
	input_type  InputType
	extension   string
	limit_mb    int
	batch_size  int
	max_rows    int
	max_cols    int
	format      OutputFormat
	prefix      string
	full_result bool
	temp_dir    string
}

// load parses the YAML file at path, validates all fields, and returns a Config.
// Returns an error if the file cannot be read, any enum value is unrecognised,
// kmer_size is outside [1, 15], source_dir is empty, or threads is less than 1.
pub fn load(path string) !Config {
	raw := yaml.unmarshal_file[RawConfig](path)!

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
