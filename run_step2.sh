#!/usr/bin/env bash
# Runs the k-mer counting benchmark (exact dict vs Count-Min Sketch)
# against whatever reads are currently in data/reads/.
#
# Run ./run_step1.sh first if data/reads/ is empty.

set -euo pipefail

cd "$(dirname "$0")/src"
python3 benchmark_kmer_structures.py \
  --fastq ../data/reads/normal_.fq ../data/reads/tumor_.fq \
  --k 21 \
  --cms-width 2000 \
  --cms-depth 4
