"""
benchmark_kmer_structures.py

Runs exact vs. Count-Min Sketch k-mer counting on the same FASTQ file(s)
and reports memory, speed, and accuracy trade-offs.

Usage:
    python src/benchmark_kmer_structures.py \
        --fastq data/reads/normal_.fq data/reads/tumor_.fq \
        --k 21 --cms-width 2000 --cms-depth 4
"""

import argparse
import time
import tracemalloc
from pathlib import Path

from kmer_utils import kmers_from_fastq
from kmer_structures import ExactCounter, CountMinSketch


def build_exact(kmer_iter):
    tracemalloc.start()
    t0 = time.perf_counter()
    counter = ExactCounter()
    for kmer in kmer_iter:
        counter.add(kmer)
    elapsed = time.perf_counter() - t0
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return counter, elapsed, peak


def build_cms(kmer_iter, width, depth):
    tracemalloc.start()
    t0 = time.perf_counter()
    cms = CountMinSketch(width=width, depth=depth)
    for kmer in kmer_iter:
        cms.add(kmer)
    elapsed = time.perf_counter() - t0
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return cms, elapsed, peak


def evaluate_accuracy(exact: ExactCounter, cms: CountMinSketch, sample_size: int = 50):
    """
    Compare CMS estimates against ground truth for a sample of k-mers:
    the most frequent ones (where collisions matter most) plus a random
    sample of the rest. CMS overestimates only — errors should never
    be negative.
    """
    import random
    all_kmers = list(exact._counts.keys())
    if not all_kmers:
        return []

    top_kmers = [k for k, _ in exact.top(min(sample_size // 2, len(all_kmers)))]
    remaining = [k for k in all_kmers if k not in set(top_kmers)]
    random.seed(0)
    sample_rest = random.sample(remaining, min(sample_size - len(top_kmers), len(remaining)))
    sample = top_kmers + sample_rest

    rows = []
    for kmer in sample:
        true_count = exact.query(kmer)
        est_count = cms.query(kmer)
        rows.append((kmer, true_count, est_count, est_count - true_count))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", nargs="+", required=True, type=Path)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--cms-width", type=int, default=2000)
    parser.add_argument("--cms-depth", type=int, default=4)
    args = parser.parse_args()

    def all_kmers():
        for fq in args.fastq:
            yield from kmers_from_fastq(fq, args.k)

    print(f"Reading k-mers (k={args.k}) from: {', '.join(str(f) for f in args.fastq)}\n")

    exact, exact_time, exact_mem = build_exact(all_kmers())
    print(f"[exact dict]        distinct k-mers = {exact.n_distinct():>8}   "
          f"time = {exact_time*1000:6.2f} ms   peak memory = {exact_mem:,} bytes "
          f"({exact.approx_memory_bytes():,} bytes by object accounting)")

    cms, cms_time, cms_mem = build_cms(all_kmers(), args.cms_width, args.cms_depth)
    print(f"[count-min sketch]  fixed table     = {args.cms_width}x{args.cms_depth}  "
          f"time = {cms_time*1000:6.2f} ms   peak memory = {cms_mem:,} bytes "
          f"({cms.approx_memory_bytes():,} bytes, table itself)")

    print(f"\nMemory ratio (exact / CMS, object-accounting basis): "
          f"{exact.approx_memory_bytes() / cms.approx_memory_bytes():.2f}x")

    print(f"\nAccuracy check on {min(50, exact.n_distinct())} sampled k-mers "
          f"(top-count k-mers + a random sample):")
    print(f"{'k-mer':<25}{'true':>6}{'estimated':>11}{'error':>8}")
    rows = evaluate_accuracy(exact, cms)
    n_exact_matches = 0
    max_error = 0
    for kmer, true_c, est_c, err in rows:
        print(f"{kmer:<25}{true_c:>6}{est_c:>11}{err:>8}")
        if err == 0:
            n_exact_matches += 1
        max_error = max(max_error, err)
        assert err >= 0, "CMS estimate came in BELOW true count — that should never happen"

    print(f"\n{n_exact_matches}/{len(rows)} sampled k-mers matched exactly, "
          f"max overestimation = {max_error}")
    print("\nNote: at this toy gene-panel scale, the exact dict is already small, "
          "so the memory gap is modest — the point is the methodology and error "
          "bounds, which are what actually matter once this runs against a real "
          "sequencing-depth dataset.")
    if cms_time > exact_time:
        print(f"\nNote: CMS was slower here ({cms_time*1000:.0f}ms vs {exact_time*1000:.0f}ms) "
              f"despite using less memory — that's the real tradeoff, not a bug. "
              f"blake2b (cryptographic, used here for determinism with zero extra "
              f"dependencies) is heavier per call than the non-cryptographic hashes "
              f"(xxHash, MurmurHash) that production k-mer tools use for this exact "
              f"reason. Swapping the hash function is the natural next optimization.")


if __name__ == "__main__":
    main()
