module io

// io is the I/O layer for NBV. It handles FASTA parsing, k-mer file reading,
// model serialization (NBV binary format and legacy NBC++ format), and
// formatted output writing (CSV, TSV, JSON Lines).

import os

// FastaRecord holds a single parsed FASTA entry, consisting of the header
// line (without the leading '>') and the concatenated sequence bytes.
pub struct FastaRecord {
pub:
	header   string
	sequence []u8
}

// parse_fasta reads the FASTA file at path and invokes callback once per
// record, passing the header string and sequence bytes. Records are delivered
// in file order. Returns an error if the file cannot be read.
pub fn parse_fasta(path string, callback fn (string, []u8)) ! {
	lines := os.read_lines(path)!
	mut current_header := ''
	mut current_seq := []u8{}
	mut in_record := false

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

// read_fasta reads the FASTA file at path and returns all records as a slice
// of FastaRecord. Prefer this over parse_fasta when the full record list is
// needed, since V 0.5.0 mutable closure captures do not propagate mutations.
// Returns an error if the file cannot be read.
pub fn read_fasta(path string) ![]FastaRecord {
	lines := os.read_lines(path)!
	mut records := []FastaRecord{}
	mut current_header := ''
	mut current_seq := []u8{}
	mut in_record := false

	for line in lines {
		if line.len > 0 && line[0] == `>` {
			if in_record {
				records << FastaRecord{
					header:   current_header
					sequence: current_seq
				}
			}
			current_header = line[1..].trim_space()
			current_seq = []u8{}
			in_record = true
		} else if in_record {
			current_seq << line.bytes()
		}
	}

	if in_record {
		records << FastaRecord{
			header:   current_header
			sequence: current_seq
		}
	}

	return records
}

// count_sequences counts the number of sequence records in the FASTA file at
// path by counting '>' header lines. Returns an error if the file cannot be read.
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
