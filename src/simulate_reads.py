"""
simulate_reads.py

Simulates Illumina-style short reads from a "normal" reference and its
mutated "tumor" counterpart, so the rest of the pipeline can work with
realistic FASTQ data (sequencing errors and all) instead of clean
reference sequence.

Two backends:

- ART (art_illumina), if it's on PATH. This is an established,
  well-validated tool (you've already used it in your
  KRAS-targeted-NGS-pipeline project) with real empirically-derived
  Illumina error/quality profiles — reimplementing that from scratch
  wouldn't demonstrate anything useful.
- A built-in pure-Python fallback, used automatically when ART isn't
  available (e.g. native Windows, where ART has no straightforward
  binary). It's a simpler uniform-error model, not ART's empirical
  profile — good enough for what this pipeline actually needs
  (plausible substitution errors, both strands represented), but
  worth being explicit that it's the simpler of the two, not a silent
  equivalent.
"""

import argparse
import random
import shutil
import subprocess
from pathlib import Path

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq))


def read_reference(path: Path) -> str:
    seq_chunks = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def simulate_reads_pure_python(reference: Path, out_fastq: Path, read_length: int,
                                coverage: float, error_rate: float, seed: int,
                                label: str):
    """
    Uniform-error fallback simulator: samples random start positions
    across the reference (both strands, 50/50) at the requested
    coverage, then injects independent per-base substitution errors at
    a fixed rate. Quality scores are a constant Q30 ('?' in Phred+33)
    since nothing downstream in this pipeline reads quality values —
    only the sequence matters for k-mer counting and the de Bruijn
    graph, so a realistic quality *curve* isn't worth the complexity
    here. That's a deliberate scope decision, not an oversight.
    """
    rng = random.Random(seed)
    ref_seq = read_reference(reference)
    ref_len = len(ref_seq)
    if ref_len < read_length:
        raise ValueError(f"Reference ({ref_len} bp) shorter than read length ({read_length} bp)")

    n_reads = max(1, round(coverage * ref_len / read_length))
    bases = "ACGT"
    qual_char = chr(30 + 33)  # constant Q30

    with open(out_fastq, "w") as f:
        for i in range(n_reads):
            start = rng.randint(0, ref_len - read_length)
            fragment = ref_seq[start:start + read_length]
            if rng.random() < 0.5:
                fragment = reverse_complement(fragment)

            read_chars = list(fragment)
            for pos in range(len(read_chars)):
                if rng.random() < error_rate:
                    true_base = read_chars[pos]
                    read_chars[pos] = rng.choice([b for b in bases if b != true_base])
            read_seq = "".join(read_chars)

            f.write(f"@{label}_read{i}_pos{start}\n")
            f.write(f"{read_seq}\n")
            f.write("+\n")
            f.write(f"{qual_char * len(read_seq)}\n")

    return n_reads


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
    parser.add_argument("--error-rate", type=float, default=0.001,
                         help="Per-base substitution error rate for the pure-Python fallback")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force-pure-python", action="store_true",
                         help="Skip ART even if it's on PATH, and use the fallback simulator")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    use_art = shutil.which("art_illumina") is not None and not args.force_pure_python
    backend = "ART (art_illumina)" if use_art else "pure-Python fallback"
    print(f"[backend] {backend}\n")

    for label, ref in [("normal", args.normal_ref), ("tumor", args.tumor_ref)]:
        if use_art:
            out_prefix = args.out_dir / f"{label}_"
            print(f"[simulate] {label}: {ref} -> {out_prefix}*.fq "
                  f"(len={args.read_length}, cov={args.coverage}x)")
            run_art(ref, out_prefix, args.read_length, args.coverage,
                    args.seq_system, args.seed)
        else:
            out_fastq = args.out_dir / f"{label}_.fq"
            print(f"[simulate] {label}: {ref} -> {out_fastq} "
                  f"(len={args.read_length}, cov={args.coverage}x, "
                  f"error_rate={args.error_rate})")
            n_reads = simulate_reads_pure_python(
                ref, out_fastq, args.read_length, args.coverage,
                args.error_rate, args.seed, label)
            print(f"[ok] {out_fastq.name}: {n_reads} reads")

    if use_art:
        for fq in sorted(args.out_dir.glob("*.fq")):
            n_reads = sum(1 for _ in open(fq)) // 4
            print(f"[ok] {fq.name}: {n_reads} reads")


if __name__ == "__main__":
    main()
