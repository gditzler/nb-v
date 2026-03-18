module io

import os

pub struct FastaRecord {
pub:
	header   string
	sequence []u8
}

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
