module io

import os
import src.config

// Writer writes classification results to a file in CSV, TSV, or JSON Lines
// format. Use Writer.new to construct, then write_header (for CSV/TSV),
// write_result or write_full_result per sequence, and close when done.
pub struct Writer {
mut:
	file        os.File
	format      config.OutputFormat
	full_result bool
}

// Writer.new creates a new Writer that writes to path in the given format.
// When full_result is true, callers should use write_full_result to emit all
// per-class scores; otherwise use write_result to emit only the best class.
// Returns an error if the output file cannot be created.
pub fn Writer.new(path string, format config.OutputFormat, full_result bool) !Writer {
	f := os.create(path)!
	return Writer{
		file:        f
		format:      format
		full_result: full_result
	}
}

// write_header writes the column header row using class_ids as the class
// column names. Has no effect for JSON Lines output.
pub fn (mut self Writer) write_header(class_ids []string) ! {
	match self.format {
		.csv { self.file.writeln('seq_id,' + class_ids.join(','))! }
		.tsv { self.file.writeln('seq_id\t' + class_ids.join('\t'))! }
		.json {}
	}
}

// write_result writes a single classification result containing only the
// winning class and its log-likelihood score. Use this when full_result is false.
pub fn (mut self Writer) write_result(seq_id string, best_class string, score f64) ! {
	match self.format {
		.csv { self.file.writeln('${seq_id},${best_class},${score}')! }
		.tsv { self.file.writeln('${seq_id}\t${best_class}\t${score}')! }
		.json { self.file.writeln('{"seq_id":"${seq_id}","best_class":"${best_class}","score":${score}}')! }
	}
}

// write_full_result writes one row containing the log-likelihood score for
// every class. class_order determines column order and must be consistent with
// the header written by write_header. Use this when full_result is true.
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

// write_no_valid_kmers writes a sentinel row for sequences that contain no
// valid k-mers and therefore cannot be classified.
pub fn (mut self Writer) write_no_valid_kmers(seq_id string) ! {
	match self.format {
		.csv { self.file.writeln('${seq_id},sequence contains no valid kmers,')! }
		.tsv { self.file.writeln('${seq_id}\tsequence contains no valid kmers\t')! }
		.json { self.file.writeln('{"seq_id":"${seq_id}","best_class":"sequence contains no valid kmers","score":null}')! }
	}
}

// close flushes and closes the underlying output file.
pub fn (mut self Writer) close() ! {
	self.file.close()
}

// output_filename returns the full output file path for a given prefix and
// format, appending the appropriate extension (.csv, .tsv, or .jsonl).
pub fn output_filename(prefix string, format config.OutputFormat) string {
	ext := match format {
		.csv { 'csv' }
		.tsv { 'tsv' }
		.json { 'jsonl' }
	}
	return '${prefix}.${ext}'
}
