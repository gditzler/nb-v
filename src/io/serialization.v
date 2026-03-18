module io

import os
import src.model
import math

const magic_bytes = [u8(`N`), `B`, `V`, `1`]
const format_version = u8(1)

// Save an NbClass in NBV binary format.
// Header: magic(4) + version(1) + kmer_size(4) + ngenomes(4) + sumfreq(8) + id_len(4) + id(id_len)
// Body: repeated (kmer_int i32, count i32) pairs
pub fn save_class(cls model.NbClass, path string) ! {
	mut f := os.create(path)!
	defer {
		f.close()
	}

	f.write(magic_bytes)!
	f.write([format_version])!
	f.write(i32_to_bytes(cls.kmer_size))!
	f.write(i32_to_bytes(cls.ngenomes))!
	f.write(i64_to_bytes(cls.sumfreq))!
	id_bytes := cls.id.bytes()
	f.write(i32_to_bytes(id_bytes.len))!
	f.write(id_bytes)!
	for km, count in cls.freqcnt {
		f.write(i32_to_bytes(km))!
		f.write(i32_to_bytes(count))!
	}
}

pub fn load_class(path string) !model.NbClass {
	data := os.read_bytes(path)!
	mut pos := 0

	if data.len < 4 || data[0..4] != magic_bytes {
		return error('invalid NBV file: bad magic bytes in ${path}')
	}
	pos = 4

	version := data[pos]
	if version != format_version {
		return error('unsupported NBV format version ${version}')
	}
	pos++

	kmer_size := bytes_to_i32(data[pos..pos + 4])
	pos += 4
	ngenomes := bytes_to_i32(data[pos..pos + 4])
	pos += 4
	sumfreq := bytes_to_i64(data[pos..pos + 8])
	pos += 8

	id_len := bytes_to_i32(data[pos..pos + 4])
	pos += 4
	id := data[pos..pos + id_len].bytestr()
	pos += id_len

	mut freqcnt := map[int]int{}
	mut freqcnt_lg := map[int]f64{}
	for pos + 8 <= data.len {
		km := bytes_to_i32(data[pos..pos + 4])
		pos += 4
		count := bytes_to_i32(data[pos..pos + 4])
		pos += 4
		freqcnt[km] = count
		freqcnt_lg[km] = math.log(f64(count + 1))
	}

	return model.NbClass{
		id:          id
		kmer_size:   kmer_size
		savefile:    path
		ngenomes:    ngenomes
		sumfreq:     sumfreq
		ngenomes_lg: math.log(f64(ngenomes))
		sumfreq_lg:  math.log(f64(sumfreq))
		freqcnt:     freqcnt
		freqcnt_lg:  freqcnt_lg
		state:       .full
	}
}

// Load NBC++ legacy -save.dat file. Log-space fields only.
// Format: f64(ngenomes_lg) + f64(sumfreq_lg) + i32(n_entries) + n_entries*(i32 kmer, f64 freqcnt_lg)
pub fn load_legacy_class(path string, k int) !model.NbClass {
	data := os.read_bytes(path)!
	mut pos := 0

	if data.len < 20 {
		return error('legacy save file too small: ${path}')
	}

	ngenomes_lg := bytes_to_f64(data[pos..pos + 8])
	pos += 8
	sumfreq_lg := bytes_to_f64(data[pos..pos + 8])
	pos += 8
	n_entries := bytes_to_i32(data[pos..pos + 4])
	pos += 4

	expected_remaining := n_entries * 12
	if data.len - pos < expected_remaining {
		return error('legacy save file truncated: ${path} (expected ${n_entries} entries)')
	}

	mut freqcnt_lg := map[int]f64{}
	for _ in 0 .. n_entries {
		km := bytes_to_i32(data[pos..pos + 4])
		pos += 4
		val := bytes_to_f64(data[pos..pos + 8])
		pos += 8
		freqcnt_lg[km] = val
	}

	basename := os.file_name(path)
	id := basename.replace('-save.dat', '')

	return model.NbClass{
		id:          id
		kmer_size:   k
		savefile:    path
		ngenomes_lg: ngenomes_lg
		sumfreq_lg:  sumfreq_lg
		freqcnt_lg:  freqcnt_lg
		state:       .classify_only
	}
}

pub fn save_meta(save_dir string, kmer_size int) ! {
	os.mkdir_all(save_dir)!
	os.write_file('${save_dir}/meta.nbv', '${kmer_size}')!
}

pub fn load_meta(save_dir string) !int {
	content := os.read_file('${save_dir}/meta.nbv')!
	return content.trim_space().int()
}

// -- Byte conversion helpers (little-endian) --

fn i32_to_bytes(val int) []u8 {
	mut b := []u8{len: 4}
	b[0] = u8(val)
	b[1] = u8(val >> 8)
	b[2] = u8(val >> 16)
	b[3] = u8(val >> 24)
	return b
}

fn bytes_to_i32(b []u8) int {
	return int(b[0]) | (int(b[1]) << 8) | (int(b[2]) << 16) | (int(b[3]) << 24)
}

fn i64_to_bytes(val i64) []u8 {
	mut b := []u8{len: 8}
	for i in 0 .. 8 {
		b[i] = u8(val >> (i * 8))
	}
	return b
}

fn bytes_to_i64(b []u8) i64 {
	mut result := i64(0)
	for i in 0 .. 8 {
		result |= i64(b[i]) << (i * 8)
	}
	return result
}

fn bytes_to_f64(b []u8) f64 {
	bits := u64(bytes_to_i64(b))
	return unsafe { *(&f64(&bits)) }
}
