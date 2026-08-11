"""
compare_to_alignment_based_calling.py

The honest head-to-head this project has been building toward: run the
same tumor reads through a standard alignment-based pipeline and
compare its recall/precision against reference-free bubble calling,
both checked against the known ground-truth mutations.

Two backends:

- BWA-MEM + samtools + bcftools, if all three are on PATH. Real,
  established tools (same as KRAS-targeted-NGS-pipeline). No native
  Windows build exists for any of them, and bioconda doesn't ship
  win-64 builds either — WSL is the real fix if you want this path.
- A pure-Python fallback, used automatically otherwise: brute-force
  best-offset alignment (try every reference position, both strand
  orientations, pick the lowest Hamming distance) plus a naive
  frequency-threshold pileup caller. This works because our simulated
  reads only carry substitution errors, never indels — a brute-force
  ungapped aligner is a legitimate, honest approach for that case, not
  a shortcut pretending to be BWA. It would NOT be adequate for real
  sequencing data with indels or a larger, repetitive reference; that
  limitation is real and worth stating rather than glossing over.

This does NOT assume the reference-free approach wins — the honest
expectation, going in, is that alignment-based calling should do at
least as well on a single small, non-repetitive reference like this.
The interesting result is *how close* bubble calling gets without ever
touching a reference coordinate, and where exactly it falls short.

Usage:
    python src/compare_to_alignment_based_calling.py \
        --reference data/reference/PLACEHOLDER_synthetic_test_gene.fasta \
        --tumor-fastq data/reads/tumor_.fq \
        --normal-fastq data/reads/normal_.fq \
        --mutations data/mutations/example_mutations.csv \
        --gene PLACEHOLDER_synthetic_test_gene \
        --bubble-positions 53,158,335,408
"""

import argparse
import csv
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from kmer_utils import read_fastq, reverse_complement


def run(cmd, **kwargs):
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{result.stderr}")
    return result


def require_tools():
    return all(shutil.which(t) is not None for t in ("bwa", "samtools", "bcftools"))


def read_fasta(path: Path) -> str:
    seq_chunks = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def best_offset_alignment(read: str, reference: str, max_mismatch_fraction: float = 0.2):
    """
    Brute-force ungapped alignment: try every valid start offset in the
    reference, in both orientations, and return whichever gives the
    fewest mismatches. O(len(reference) * len(read)) per read — fine at
    gene-panel scale, would not scale to a real genome.

    Returns (offset, oriented_read, n_mismatches) or None if nothing
    beats max_mismatch_fraction of the read length.
    """
    best = None  # (mismatches, offset, oriented_read)
    for orientation in (read, reverse_complement(read)):
        for offset in range(len(reference) - len(orientation) + 1):
            window = reference[offset:offset + len(orientation)]
            mismatches = sum(1 for a, b in zip(orientation, window) if a != b)
            if best is None or mismatches < best[0]:
                best = (mismatches, offset, orientation)

    if best is None:
        return None
    mismatches, offset, oriented_read = best
    if mismatches > max_mismatch_fraction * len(oriented_read):
        return None
    return offset, oriented_read, mismatches


def build_pileup_pure_python(fastq_path: Path, reference: str):
    """Align every read to `reference` and accumulate per-position base
    counts: pileup[position_0based][base] = count."""
    pileup = defaultdict(lambda: defaultdict(int))
    n_aligned, n_unaligned = 0, 0
    for _header, seq, _qual in read_fastq(fastq_path):
        result = best_offset_alignment(seq, reference)
        if result is None:
            n_unaligned += 1
            continue
        n_aligned += 1
        offset, oriented_read, _mismatches = result
        for i, base in enumerate(oriented_read):
            pileup[offset + i][base] += 1
    return pileup, n_aligned, n_unaligned


def call_variants_naive(pileup, reference: str, min_alt_reads: int = 3,
                         min_alt_fraction: float = 0.2):
    """
    Simple frequency-threshold caller: at each reference position, call
    a variant if a non-reference base has at least min_alt_reads
    supporting reads AND makes up at least min_alt_fraction of total
    coverage there. This is a real simplification versus bcftools'
    genotype-likelihood model (no base quality weighting, no proper
    statistical model) — stated here rather than dressed up as
    equivalent.
    """
    called = set()
    for pos0, base_counts in pileup.items():
        if pos0 >= len(reference):
            continue
        ref_base = reference[pos0]
        total = sum(base_counts.values())
        for base, count in base_counts.items():
            if base == ref_base or total == 0:
                continue
            if count >= min_alt_reads and (count / total) >= min_alt_fraction:
                called.add(pos0 + 1)  # convert to 1-based
                break
    return called


def align_and_call_pure_python(reference_path: Path, tumor_fastq: Path):
    reference = read_fasta(reference_path)
    pileup, n_aligned, n_unaligned = build_pileup_pure_python(tumor_fastq, reference)
    print(f"[pure-python aligner] {n_aligned} reads aligned, {n_unaligned} discarded "
          f"(exceeded mismatch threshold)")
    return call_variants_naive(pileup, reference)


def align_and_call(reference: Path, fastq: Path, work_dir: Path, label: str):
    """BWA-MEM -> sort -> index -> bcftools mpileup+call, same tools and
    ordering as the KRAS-targeted-NGS-pipeline project."""
    work_dir.mkdir(parents=True, exist_ok=True)
    sam = work_dir / f"{label}.sam"
    bam = work_dir / f"{label}.sorted.bam"
    vcf = work_dir / f"{label}.vcf"

    if not (reference.with_suffix(reference.suffix + ".bwt")).exists():
        run(["bwa", "index", str(reference)])

    with open(sam, "w") as f:
        result = subprocess.run(["bwa", "mem", str(reference), str(fastq)],
                                 stdout=f, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bwa mem failed: {result.stderr}")

    run(["samtools", "sort", "-o", str(bam), str(sam)])
    run(["samtools", "index", str(bam)])

    mpileup = subprocess.run(["bcftools", "mpileup", "-f", str(reference), str(bam)],
                              capture_output=True, text=True)
    if mpileup.returncode != 0:
        raise RuntimeError(f"bcftools mpileup failed: {mpileup.stderr}")
    call = subprocess.run(["bcftools", "call", "-mv", "-Ov"],
                           input=mpileup.stdout, capture_output=True, text=True)
    vcf.write_text(call.stdout)
    return vcf


def parse_vcf_positions(vcf_path: Path):
    """Return the set of 1-based positions bcftools called a variant at."""
    positions = set()
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            positions.add(int(fields[1]))
    return positions


def load_ground_truth(mutations_csv: Path, gene: str):
    with open(mutations_csv) as f:
        rows = [r for r in csv.DictReader(f) if r["gene"] == gene]
    return set(int(r["cds_position"]) for r in rows)


def precision_recall(called: set, truth: set):
    tp = called & truth
    fp = called - truth
    fn = truth - called
    precision = len(tp) / len(called) if called else float("nan")
    recall = len(tp) / len(truth) if truth else float("nan")
    return precision, recall, tp, fp, fn


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--tumor-fastq", required=True, type=Path)
    parser.add_argument("--normal-fastq", required=True, type=Path)
    parser.add_argument("--mutations", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--bubble-positions", required=True,
                         help="Comma-separated positions recovered by run_step3.sh, "
                              "e.g. 53,158,335,408")
    parser.add_argument("--work-dir", type=Path, default=Path("data/alignment_calling"))
    args = parser.parse_args()

    truth = load_ground_truth(args.mutations, args.gene)
    bubble_called = set(int(p) for p in args.bubble_positions.split(","))

    use_bwa = require_tools()
    backend = "BWA-MEM + bcftools" if use_bwa else "pure-Python fallback"
    print(f"[backend] {backend}\n")

    if use_bwa:
        print("Aligning tumor reads to the NORMAL reference (standard somatic calling setup) ...")
        tumor_vcf = align_and_call(args.reference, args.tumor_fastq, args.work_dir, "tumor_vs_normal_ref")
        alignment_called = parse_vcf_positions(tumor_vcf)
    else:
        print("Aligning tumor reads to the NORMAL reference (standard somatic calling setup) ...")
        alignment_called = align_and_call_pure_python(args.reference, args.tumor_fastq)

    print(f"[caller] called {len(alignment_called)} variant position(s): {sorted(alignment_called)}")

    print(f"\n=== Ground truth ===")
    print(f"Injected mutations: {sorted(truth)}")

    print(f"\n=== Alignment-based calling ===")
    p, r, tp, fp, fn = precision_recall(alignment_called, truth)
    print(f"precision={p:.2f}  recall={r:.2f}  "
          f"true_positive={sorted(tp)}  false_positive={sorted(fp)}  false_negative={sorted(fn)}")

    print(f"\n=== Reference-free bubble calling ===")
    p2, r2, tp2, fp2, fn2 = precision_recall(bubble_called, truth)
    print(f"precision={p2:.2f}  recall={r2:.2f}  "
          f"true_positive={sorted(tp2)}  false_positive={sorted(fp2)}  false_negative={sorted(fn2)}")

    print(f"\n=== Head-to-head ===")
    print(f"Alignment-based: {len(tp)}/{len(truth)} recovered, {len(fp)} false positive(s)")
    print(f"Bubble calling:  {len(tp2)}/{len(truth)} recovered, {len(fp2)} false positive(s)")


if __name__ == "__main__":
    main()