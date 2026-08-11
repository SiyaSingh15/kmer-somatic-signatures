#!/usr/bin/env bash
# Runs the full "grounded read simulation" pipeline end to end.
#
# By default this uses the synthetic placeholder gene so the pipeline is
# runnable immediately. Once you've run src/download_reference.py locally
# and updated data/mutations/ with your real IDR mutation calls, change
# GENE / REFERENCE / TUMOR_REFERENCE below (or pass them as args).

set -euo pipefail

GENE="${1:-PLACEHOLDER_synthetic_test_gene}"
REFERENCE="data/reference/${GENE}.fasta"
MUTATIONS="data/mutations/example_mutations.csv"
TUMOR_REFERENCE="data/reference/${GENE}_tumor.fasta"

echo "== Step 1: inject mutations =="
python3 src/mutate_reference.py \
  --reference "$REFERENCE" \
  --mutations "$MUTATIONS" \
  --gene "$GENE" \
  --out "$TUMOR_REFERENCE"

echo
echo "== Step 2: simulate normal + tumor reads =="
python3 src/simulate_reads.py \
  --normal-ref "$REFERENCE" \
  --tumor-ref "$TUMOR_REFERENCE" \
  --out-dir data/reads \
  --read-length 100 \
  --coverage 30

echo
echo "Done. Reads in data/reads/normal_.fq and data/reads/tumor_.fq"
