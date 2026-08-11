#!/usr/bin/env bash
# Runs the head-to-head: BWA+bcftools alignment-based calling vs
# reference-free bubble calling, both checked against known mutations.
#
# Run ./run_step1.sh and ./run_step3.sh first, and pass the bubble
# positions that step 3 printed as the argument here.

set -euo pipefail

GENE="${1:-PLACEHOLDER_synthetic_test_gene}"
BUBBLE_POSITIONS="${2:?Usage: ./run_step4_compare.sh <gene> <bubble_positions_comma_separated>}"

python3 src/compare_to_alignment_based_calling.py \
  --reference "data/reference/${GENE}.fasta" \
  --tumor-fastq data/reads/tumor_.fq \
  --normal-fastq data/reads/normal_.fq \
  --mutations data/mutations/example_mutations.csv \
  --gene "$GENE" \
  --bubble-positions "$BUBBLE_POSITIONS"
