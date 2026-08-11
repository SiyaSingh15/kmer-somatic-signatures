"""
mutate_reference.py

Applies point mutations from a CSV table onto a reference CDS FASTA to
produce a "tumor" reference, and reports the resulting amino-acid changes.

Deliberately does not use Biopython's Seq.translate() — the codon table
and translation logic are written out explicitly here, since the whole
point of this project is to demonstrate hand-rolled sequence logic rather
than library calls.

Mutation table format (CSV):
    gene,cds_position,ref_base,alt_base,source
    SNAI1,53,G,A,illustrative_example

cds_position is 1-based, indexed into the CDS FASTA (start codon = position 1).

Usage:
    python src/mutate_reference.py \
        --reference data/reference/PLACEHOLDER_synthetic_test_gene.fasta \
        --mutations data/mutations/example_mutations.csv \
        --gene PLACEHOLDER_synthetic_test_gene \
        --out data/reference/PLACEHOLDER_synthetic_test_gene_tumor.fasta
"""

import argparse
import csv
from pathlib import Path

# Standard genetic code, written out explicitly (not imported).
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def read_fasta(path: Path):
    """Return (header, sequence) for a single-record FASTA file."""
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
            elif line:
                seq_chunks.append(line)
    return header, "".join(seq_chunks)


def write_fasta(path: Path, header: str, seq: str, width: int = 70):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")  # X = unknown/ambiguous


def codon_at_position(seq: str, cds_position_1based: int):
    """
    Given a 1-based nucleotide position in the CDS, return
    (codon_number_1based, codon_seq, position_within_codon_0to2).
    """
    idx0 = cds_position_1based - 1
    codon_number = idx0 // 3 + 1
    pos_in_codon = idx0 % 3
    codon_start = idx0 - pos_in_codon
    codon_seq = seq[codon_start:codon_start + 3]
    return codon_number, codon_seq, pos_in_codon


def apply_mutations(seq: str, mutations: list, gene: str):
    """
    Apply a list of mutation dicts (from the CSV) to seq.
    Returns (mutated_seq, report_rows). Raises ValueError on any
    ref-base mismatch, since a silent mismatch would mean the mutation
    table doesn't actually correspond to this reference.
    """
    seq_chars = list(seq)
    report_rows = []

    for m in mutations:
        if m["gene"] != gene:
            continue

        pos = int(m["cds_position"])
        ref_base = m["ref_base"].upper()
        alt_base = m["alt_base"].upper()

        actual_base = seq_chars[pos - 1]
        if actual_base != ref_base:
            raise ValueError(
                f"Ref-base mismatch at {gene}:{pos} — table says {ref_base}, "
                f"reference has {actual_base}. Check your coordinates "
                f"(are you sure this mutation table matches this FASTA?)."
            )

        codon_number, codon_before, pos_in_codon = codon_at_position(seq, pos)
        aa_before = translate_codon(codon_before)

        seq_chars[pos - 1] = alt_base
        mutated_seq_so_far = "".join(seq_chars)
        _, codon_after, _ = codon_at_position(mutated_seq_so_far, pos)
        aa_after = translate_codon(codon_after)

        if aa_before == aa_after:
            consequence = "synonymous"
        elif aa_after == "*":
            consequence = "nonsense"
        else:
            consequence = "missense"

        report_rows.append({
            "gene": gene,
            "cds_position": pos,
            "ref_base": ref_base,
            "alt_base": alt_base,
            "codon_number": codon_number,
            "codon_before": codon_before,
            "codon_after": codon_after,
            "aa_before": aa_before,
            "aa_after": aa_after,
            "consequence": consequence,
        })

    return "".join(seq_chars), report_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--mutations", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    header, seq = read_fasta(args.reference)
    with open(args.mutations) as f:
        mutations = list(csv.DictReader(f))

    mutated_seq, report_rows = apply_mutations(seq, mutations, args.gene)

    if not report_rows:
        print(f"[warn] no mutations found for gene={args.gene} in {args.mutations}")

    write_fasta(args.out, f"{header} | TUMOR (mutated)", mutated_seq)
    print(f"Wrote mutated reference: {args.out}\n")

    print(f"{'pos':>5} {'ref>alt':>8} {'codon#':>7} {'codon':>12} {'aa':>6} {'consequence'}")
    for r in report_rows:
        print(f"{r['cds_position']:>5} "
              f"{r['ref_base']}>{r['alt_base']:>6} "
              f"{r['codon_number']:>7} "
              f"{r['codon_before']}>{r['codon_after']:>9} "
              f"{r['aa_before']}>{r['aa_after']:>4} "
              f"{r['consequence']}")


if __name__ == "__main__":
    main()
