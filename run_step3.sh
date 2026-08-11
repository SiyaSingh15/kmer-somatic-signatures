#!/usr/bin/env bash
# Builds the colored de Bruijn graph from data/reads/ and detects bubbles,
# validating against the known mutations in data/mutations/.
#
# Run ./run_step1.sh first if data/reads/ is empty.

set -euo pipefail

GENE="${1:-PLACEHOLDER_synthetic_test_gene}"

cd "$(dirname "$0")/src"
python3 call_variants_from_bubbles.py \
  --normal-fastq ../data/reads/normal_.fq \
  --tumor-fastq ../data/reads/tumor_.fq \
  --reference "../data/reference/${GENE}.fasta" \
  --mutations ../data/mutations/example_mutations.csv \
  --gene "$GENE" \
  --k 21 --min-support 3 --min-color-support 2
