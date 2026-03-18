module main

import os
import src.config

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
