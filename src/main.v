module main

import os
import src.config
import src.pipeline

fn main() {
	args := os.args[1..]
	if args.len != 1 {
		eprintln('Usage: nbv <config.yaml>')
		exit(1)
	}

	cfg := config.load(args[0]) or {
		eprintln('Error: ${err}')
		exit(1)
	}

	match cfg.mode {
		.train {
			pipeline.train(cfg) or {
				eprintln('Training failed: ${err}')
				exit(1)
			}
			println('Training complete. Savefiles written to ${cfg.save_dir}')
		}
		.classify {
			pipeline.classify(cfg) or {
				eprintln('Classification failed: ${err}')
				exit(1)
			}
			println('Classification complete.')
		}
	}
}
