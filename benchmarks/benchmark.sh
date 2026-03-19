#!/usr/bin/env bash
set -euo pipefail

# Benchmark: NBV (V) vs NBC++ (C++) classification
#
# Prerequisites:
#   - V compiler (https://vlang.io)
#   - Docker (NBC++ is built and run in a container since it requires Linux/GCC)
#
# Usage:
#   ./benchmarks/benchmark.sh [threads] [runs]

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
NBCPP_DIR="$SCRIPT_DIR/Naive_Bayes"
THREADS="${1:-4}"
RUNS="${2:-3}"

EXAMPLE_DIR="$PROJECT_DIR/example"
READS_DIR="$EXAMPLE_DIR/reads"
SAVE_DIR="$EXAMPLE_DIR/training_classes"
TMP_DIR="/tmp/nbv_benchmark"
DOCKER_IMAGE="nbcpp-bench"

mkdir -p "$TMP_DIR"

echo "============================================"
echo "  NBV vs NBC++ Classification Benchmark"
echo "============================================"
echo "Threads:    $THREADS"
echo "Runs:       $RUNS"
echo "Reads:      $READS_DIR"
echo "Classes:    $SAVE_DIR"
echo ""

# --- Build NBV ---
echo "--- Building NBV (V) ---"
cd "$PROJECT_DIR"
v -prod -o "$SCRIPT_DIR/nbv_bench" src/
echo "Built: $SCRIPT_DIR/nbv_bench"

# --- Build NBC++ via Docker ---
echo ""
echo "--- Building NBC++ (C++) via Docker ---"
NBCPP_AVAILABLE=false

if ! command -v docker &>/dev/null; then
    echo "Docker not found. NBC++ benchmark will be skipped."
elif ! docker info &>/dev/null; then
    echo "Docker daemon not running. NBC++ benchmark will be skipped."
else
    if [ ! -d "$NBCPP_DIR" ]; then
        echo "Cloning NBC++ from https://github.com/EESI/Naive_Bayes ..."
        git clone https://github.com/EESI/Naive_Bayes.git "$NBCPP_DIR"
    fi

    # Build Docker image if needed
    if ! docker image inspect "$DOCKER_IMAGE" &>/dev/null; then
        echo "Building Docker image ..."
        cat > "$NBCPP_DIR/Dockerfile.bench" <<'DOCKERFILE'
FROM ubuntu:22.04
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ make libboost-all-dev && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /nbc
COPY *.cpp *.hpp Makefile ./
RUN make clean 2>/dev/null; make
DOCKERFILE
        docker build -t "$DOCKER_IMAGE" -f "$NBCPP_DIR/Dockerfile.bench" "$NBCPP_DIR"
    fi

    if docker image inspect "$DOCKER_IMAGE" &>/dev/null; then
        echo "Built: Docker image $DOCKER_IMAGE"
        NBCPP_AVAILABLE=true
    else
        echo "ERROR: Docker build failed."
    fi
fi

# --- Write NBV benchmark config ---
NBV_CONFIG="$TMP_DIR/bench_classify.yaml"
cat > "$NBV_CONFIG" <<YAML
version: 1
mode: classify
kmer_size: 9
save_dir: $SAVE_DIR
source_dir: $READS_DIR
threads: $THREADS

input:
  extension: .fna
  input_type: fasta

memory:
  limit_mb: 0
  batch_size: 0
  max_rows: 0
  max_cols: 0

output:
  format: csv
  prefix: $TMP_DIR/nbv_result
  full_result: false
  temp_dir: $TMP_DIR
YAML

# --- Benchmark function ---
run_benchmark() {
    local label="$1"
    shift
    local cmd=("$@")
    local times=()

    echo ""
    echo "--- $label ---"

    for i in $(seq 1 "$RUNS"); do
        start=$(python3 -c 'import time; print(time.time())')
        "${cmd[@]}" > /dev/null 2>&1
        end=$(python3 -c 'import time; print(time.time())')
        elapsed=$(python3 -c "print(f'{${end} - ${start}:.3f}')")
        times+=("$elapsed")
        echo "  Run $i: ${elapsed}s"
    done

    local IFS=,
    mean=$(python3 -c "
t = [${times[*]}]
print(f'{sum(t)/len(t):.3f}')
")
    echo "  Mean:  ${mean}s"
}

# --- Run benchmarks ---
echo ""
echo "============================================"
echo "  Running benchmarks ($RUNS runs each)"
echo "============================================"

run_benchmark "NBV (V, $THREADS threads)" \
    "$SCRIPT_DIR/nbv_bench" "$NBV_CONFIG"

if [ "$NBCPP_AVAILABLE" = true ]; then
    run_benchmark "NBC++ (C++, $THREADS threads, Docker)" \
        docker run --rm \
        -v "$READS_DIR:/data/reads:ro" \
        -v "$SAVE_DIR:/data/classes:ro" \
        -v "$TMP_DIR:/data/tmp" \
        "$DOCKER_IMAGE" \
        ./NB.run classify /data/reads \
        -s /data/classes -k 9 -m 30000 -t "$THREADS" \
        -d /data/tmp -o /data/tmp/nbcpp_result
else
    echo ""
    echo "--- NBC++ (C++) ---"
    echo "  SKIPPED (Docker required; install from https://docs.docker.com/get-docker/)"
fi

# --- Cleanup ---
rm -f "$TMP_DIR"/nbv_result* "$TMP_DIR"/nbcpp_result*

echo ""
echo "============================================"
echo "  Benchmark complete"
echo "============================================"
echo ""
echo "Note: NBC++ runs inside Docker (Linux container)."
echo "Container overhead may add a small constant to each run."
