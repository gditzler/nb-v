module io

import os
import src.config

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
		.csv { self.file.writeln('seq_id,' + class_ids.join(','))! }
		.tsv { self.file.writeln('seq_id\t' + class_ids.join('\t'))! }
		.json {}
	}
}

pub fn (mut self Writer) write_result(seq_id string, best_class string, score f64) ! {
	match self.format {
		.csv { self.file.writeln('${seq_id},${best_class},${score}')! }
		.tsv { self.file.writeln('${seq_id}\t${best_class}\t${score}')! }
		.json { self.file.writeln('{"seq_id":"${seq_id}","best_class":"${best_class}","score":${score}}')! }
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
		.csv { self.file.writeln('${seq_id},sequence contains no valid kmers,')! }
		.tsv { self.file.writeln('${seq_id}\tsequence contains no valid kmers\t')! }
		.json { self.file.writeln('{"seq_id":"${seq_id}","best_class":"sequence contains no valid kmers","score":null}')! }
	}
}

pub fn (mut self Writer) close() ! {
	self.file.close()
}

pub fn output_filename(prefix string, format config.OutputFormat) string {
	ext := match format {
		.csv { 'csv' }
		.tsv { 'tsv' }
		.json { 'jsonl' }
	}
	return '${prefix}.${ext}'
}
