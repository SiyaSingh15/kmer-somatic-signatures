"""
kmer_utils.py

Shared, hand-written utilities for k-mer extraction. Used by both the
counting engine (step 2) and the de Bruijn graph builder (step 3), so
the canonicalization logic only lives in one place.

No Biopython, no external k-mer libraries — this is the part of the
project meant to demonstrate implementation, not library usage.
"""

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA sequence. Non-ACGT bases pass through
    unchanged (they'll get filtered out upstream by canonical_kmers)."""
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq))


def canonical(kmer: str) -> str:
    """
    A k-mer and its reverse complement represent the same underlying
    genomic locus when you don't know which strand a read came from.
    The canonical form is whichever of the two sorts first
    lexicographically — this is the standard convention used by
    k-mer counters like KMC and Jellyfish, so that a k-mer and its
    reverse complement always collapse to the same count bucket.
    """
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def canonical_kmers(seq: str, k: int):
    """
    Yield canonical k-mers from a sequence via a sliding window.
    Windows containing any non-ACGT character (N, or sequencing
    adapter artifacts) are skipped entirely, rather than guessed at —
    an ambiguous base shouldn't silently become a wrong one.
    """
    seq = seq.upper()
    valid = set("ACGT")
    n = len(seq)
    for i in range(n - k + 1):
        window = seq[i:i + k]
        if all(b in valid for b in window):
            yield canonical(window)


def read_fastq(path):
    """
    Minimal FASTQ parser: yields (header, sequence, quality) tuples.
    FASTQ is a strict 4-line-per-record format (@header / seq / +
    / quality) so this doesn't need to handle anything fancier than
    reading four lines at a time.
    """
    with open(path) as f:
        while True:
            header = f.readline().rstrip("\n")
            if not header:
                break
            seq = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")  # noqa: F841 (unused, format spacer)
            qual = f.readline().rstrip("\n")
            if not qual:
                break
            yield header, seq, qual


def kmers_from_fastq(path, k: int):
    """Yield every canonical k-mer across all reads in a FASTQ file."""
    for _header, seq, _qual in read_fastq(path):
        yield from canonical_kmers(seq, k)
