"""
simulate_reads.py

Simulates Illumina short reads from a "normal" reference and its mutated
"tumor" counterpart using ART (art_illumina), so the rest of the pipeline
can work with realistic FASTQ data (sequencing errors and all) instead of
clean reference sequence.

We reuse ART here deliberately rather than writing a simulator from scratch —
it's an established, well-validated tool (you've already used it in your
KRAS-targeted-NGS-pipeline project), and reimplementing Illumina error
profiles from scratch wouldn't demonstrate anything useful. The parts of
this project worth hand-implementing are the k-mer/de-Bruijn-graph steps
downstream, not read simulation.

Usage:
    python src/simulate_reads.py \
        --normal-ref data/reference/PLACEHOLDER_synthetic_test_gene.fasta \
        --tumor-ref data/reference/PLACEHOLDER_synthetic_test_gene_tumor.fasta \
        --out-dir data/reads \
        --read-length 100 \
        --coverage 30
"""

import argparse
import shutil
import subprocess
from pathlib import Path


def run_art(reference: Path, out_prefix: Path, read_length: int,
            coverage: float, seq_system: str, seed: int):
    if shutil.which("art_illumina") is None:
        raise RuntimeError(
            "art_illumina not found on PATH. Install it, e.g.:\n"
            "  sudo apt-get install art-nextgen-simulation-tools"
        )

    cmd = [
        "art_illumina",
        "-ss", seq_system,       # sequencing system error/quality profile
        "-i", str(reference),
        "-l", str(read_length),
        "-f", str(coverage),
        "-rs", str(seed),
        "-o", str(out_prefix),
        "-na",                   # skip the .aln alignment file, we don't need it
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"art_illumina failed:\n{result.stderr}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normal-ref", required=True, type=Path)
    parser.add_argument("--tumor-ref", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--coverage", type=float, default=30.0)
    parser.add_argument("--seq-system", default="HS25",
                         help="ART sequencing system profile, e.g. HS25 = HiSeq 2500")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for label, ref in [("normal", args.normal_ref), ("tumor", args.tumor_ref)]:
        out_prefix = args.out_dir / f"{label}_"
        print(f"[simulate] {label}: {ref} -> {out_prefix}*.fq "
              f"(len={args.read_length}, cov={args.coverage}x)")
        run_art(ref, out_prefix, args.read_length, args.coverage,
                args.seq_system, args.seed)

    for fq in sorted(args.out_dir.glob("*.fq")):
        n_reads = sum(1 for _ in open(fq)) // 4
        print(f"[ok] {fq.name}: {n_reads} reads")


if __name__ == "__main__":
    main()
