# Expected Output

```
% bash benchmarks/benchmark.sh

============================================
  NBV vs NBC++ Classification Benchmark
============================================
Threads:    4
Runs:       3
Reads:      /Users/gditzler/git/nbv/example/reads
Classes:    /Users/gditzler/git/nbv/example/training_classes

--- Building NBV (V) ---
Built: /Users/gditzler/git/nbv/benchmarks/nbv_bench

--- Building NBC++ (C++) via Docker ---
Building Docker image ...
[+] Building 48.5s (11/11) FINISHED                                docker:desktop-linux
 => [internal] load build definition from Dockerfile.bench                         0.0s
 => => transferring dockerfile: 282B                                               0.0s
 => [internal] load metadata for docker.io/library/ubuntu:22.04                    0.9s
 => [auth] library/ubuntu:pull token for registry-1.docker.io                      0.0s
 => [internal] load .dockerignore                                                  0.0s
 => => transferring context: 2B                                                    0.0s
 => [1/5] FROM docker.io/library/ubuntu:22.04@sha256:445586e41c1de7dfda82d2637f5f  2.6s
 => => resolve docker.io/library/ubuntu:22.04@sha256:445586e41c1de7dfda82d2637f5f  0.0s
 => => sha256:cf67f3f0b7b3a837aac5c0be2974a3574a6b600345d9528de 27.39MB / 27.39MB  2.3s
 => => extracting sha256:cf67f3f0b7b3a837aac5c0be2974a3574a6b600345d9528def747c7e  0.3s
 => [internal] load build context                                                  0.0s
 => => transferring context: 126.99kB                                              0.0s
 => [2/5] RUN apt-get update &&     apt-get install -y --no-install-recommends    22.9s
 => [3/5] WORKDIR /nbc                                                             0.1s
 => [4/5] COPY *.cpp *.hpp Makefile ./                                             0.0s
 => [5/5] RUN make clean 2>/dev/null; make                                         9.6s
 => exporting to image                                                            12.4s
 => => exporting layers                                                           10.1s
 => => exporting manifest sha256:5fdb715c719e7b0d5b71ddc6da1bb30c9f8ce0fb07176327  0.0s
 => => exporting config sha256:50adc7ee00c1fb61f1ba66b93e280dd21ad0b9c41d2a0f8420  0.0s
 => => exporting attestation manifest sha256:66c10f34c90386806bbf43c146fcd013a61c  0.0s
 => => exporting manifest list sha256:c75f5b0c5d9ec4830fa0c975fae23e4532413e9b1dd  0.0s
 => => naming to docker.io/library/nbcpp-bench:latest                              0.0s
 => => unpacking to docker.io/library/nbcpp-bench:latest                           2.3s
Built: Docker image nbcpp-bench

============================================
  Running benchmarks (3 runs each)
============================================

--- NBV (V, 4 threads) ---
  Run 1: 6.428s
  Run 2: 6.031s
  Run 3: 6.061s
  Mean:  6.173s

--- NBC++ (C++, 4 threads, Docker) ---
(base) gditzler@Gregs-MacBook-Pro nbv % bash benchmarks/benchmark.sh
============================================
  NBV vs NBC++ Classification Benchmark
============================================
Threads:    4
Runs:       3
Reads:      /Users/gditzler/git/nbv/example/reads
Classes:    /Users/gditzler/git/nbv/example/training_classes

--- Building NBV (V) ---
Built: /Users/gditzler/git/nbv/benchmarks/nbv_bench

--- Building NBC++ (C++) via Docker ---
Built: Docker image nbcpp-bench

============================================
  Running benchmarks (3 runs each)
============================================

--- NBV (V, 4 threads) ---
  Run 1: 6.306s
  Run 2: 6.012s
  Run 3: 6.010s
  Mean:  6.109s

--- NBC++ (C++, 4 threads, Docker) ---
  Run 1: 24.883s
  Run 2: 24.932s
  Run 3: 24.831s
  Mean:  24.882s

============================================
  Benchmark complete
============================================

Note: NBC++ runs inside Docker (Linux container).
Container overhead may add a small constant to each run.
```

