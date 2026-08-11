"""
call_variants_from_bubbles.py

Builds a colored de Bruijn graph from normal + tumor reads, detects
bubbles, and — since this project works from simulated data with known
injected mutations — validates the detected bubbles against ground
truth: how many injected mutations produced a bubble, and how many
bubbles don't correspond to any injected mutation (false positives,
most likely from sequencing-error noise that survived pruning).

This validation step only works because the ground truth is available
here. It's the honest way to test whether reference-free bubble calling
actually works before trusting it on data where you don't already know
the answer.

Usage:
    python src/call_variants_from_bubbles.py \
        --normal-fastq data/reads/normal_.fq \
        --tumor-fastq data/reads/tumor_.fq \
        --reference data/reference/PLACEHOLDER_synthetic_test_gene.fasta \
        --mutations data/mutations/example_mutations.csv \
        --gene PLACEHOLDER_synthetic_test_gene \
        --k 21 --min-support 3
"""

import argparse
import csv
from pathlib import Path

from debruijn_graph import DeBruijnGraph
from kmer_utils import canonical_kmers


def read_fasta(path):
    seq_chunks = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def load_ground_truth(mutations_csv, gene):
    with open(mutations_csv) as f:
        rows = [r for r in csv.DictReader(f) if r["gene"] == gene]
    return [int(r["cds_position"]) for r in rows]


def anchors_flank_position(graph: DeBruijnGraph, anchor_nodes, reference_seq, position_1based, window=30):
    """
    Cross-check: does either anchor's (k-1)-mer sequence appear in the
    reference within `window` bases of the known mutation position?
    This is only possible because we have the reference here for
    validation — the bubble-calling itself never uses it.
    """
    idx0 = position_1based - 1
    region_start = max(0, idx0 - window)
    region_end = min(len(reference_seq), idx0 + window)
    region = reference_seq[region_start:region_end]

    for node in anchor_nodes:
        # a canonical (k-1)-mer might be the reverse complement of what's
        # in the reference — check both orientations
        rc = "".join({"A": "T", "T": "A", "C": "G", "G": "C"}.get(b, b) for b in reversed(node))
        if node in region or rc in region:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normal-fastq", required=True, type=Path)
    parser.add_argument("--tumor-fastq", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--mutations", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--min-support", type=int, default=3)
    parser.add_argument("--min-color-support", type=int, default=2,
                         help="Minimum reads from a sample before an edge counts as "
                              "'supported by' that sample — guards against a single "
                              "stray error-read fragmenting a clean bubble")
    args = parser.parse_args()

    print(f"Building de Bruijn graph (k={args.k}) from {args.normal_fastq} + {args.tumor_fastq} ...")
    graph = DeBruijnGraph(k=args.k)
    graph.add_reads_from_fastq(args.normal_fastq, "normal")
    graph.add_reads_from_fastq(args.tumor_fastq, "tumor")
    n_edges_before = len(graph.edge_colors)
    print(f"[graph] {n_edges_before} distinct edges, {len(graph.neighbors)} nodes before pruning")

    n_pruned = graph.prune_low_support(args.min_support)
    print(f"[prune] removed {n_pruned} edges with support < {args.min_support} "
          f"({len(graph.edge_colors)} edges remain)")

    normal_only, tumor_only, shared = graph.classify_edges(args.min_color_support)
    print(f"[classify] normal-only={len(normal_only)}  tumor-only={len(tumor_only)}  shared={len(shared)}")

    bubbles = graph.find_bubbles(args.min_color_support)
    print(f"\n[bubbles] found {len(bubbles)} candidate variant site(s)\n")

    reference_seq = read_fasta(args.reference)
    ground_truth_positions = load_ground_truth(args.mutations, args.gene)

    matched_positions = set()
    print(f"{'#':<4}{'tumor nodes':>12}{'normal nodes':>14}   anchors")
    for i, bubble in enumerate(bubbles):
        print(f"{i:<4}{len(bubble['tumor_nodes']):>12}{len(bubble['normal_nodes']):>14}   "
              f"{bubble['anchors'][0][:12]}.. / {bubble['anchors'][1][:12]}..")
        for pos in ground_truth_positions:
            if anchors_flank_position(graph, bubble["anchors"], reference_seq, pos):
                matched_positions.add(pos)

    print(f"\n=== Validation against known ground truth ===")
    print(f"Injected mutations: {ground_truth_positions}")
    print(f"Recovered via bubble detection: {sorted(matched_positions)}")
    n_missed = len(set(ground_truth_positions) - matched_positions)
    n_unexplained_bubbles = len(bubbles) - len(matched_positions)
    print(f"{len(matched_positions)}/{len(ground_truth_positions)} injected mutations recovered, "
          f"{n_missed} missed, {n_unexplained_bubbles} bubble(s) not matched to any known mutation "
          f"(candidate false positives, or sequencing-error noise that survived pruning)")


if __name__ == "__main__":
    main()
