"""
kmer_structures.py

Two ways of counting k-mer frequencies:

1. ExactCounter — a plain dict. Exact, but memory grows linearly with
   the number of *distinct* k-mers, and each Python dict entry carries
   real object overhead (~100 bytes, not just the bytes of the string
   and an int).

2. CountMinSketch — a fixed-size 2D array of counters. Memory is
   constant regardless of how many distinct k-mers you see, at the
   cost of counts that can be *overestimated* (never underestimated)
   due to hash collisions. This is the standard sketch-based approach
   to approximate frequency counting, and it's the data structure that
   makes "count k-mers across an entire sequencing run without loading
   every distinct k-mer into memory" possible at real scale.

Hashing uses hashlib (blake2b) with a per-row salt rather than Python's
built-in hash(), because str hashing in CPython is randomized per
process (PYTHONHASHSEED) — a Count-Min Sketch needs the same k-mer to
hash to the same slot every run, or querying it later would be
meaningless.
"""

import hashlib
import sys
from array import array


class ExactCounter:
    """A thin wrapper around dict, just so it shares an interface with
    CountMinSketch for the benchmark script."""

    def __init__(self):
        self._counts = {}

    def add(self, kmer: str):
        self._counts[kmer] = self._counts.get(kmer, 0) + 1

    def query(self, kmer: str) -> int:
        return self._counts.get(kmer, 0)

    def n_distinct(self) -> int:
        return len(self._counts)

    def top(self, n: int):
        return sorted(self._counts.items(), key=lambda kv: -kv[1])[:n]

    def approx_memory_bytes(self) -> int:
        """
        Rough accounting of what the dict actually costs: the dict's
        own bucket-table overhead, plus each key string and each int
        value as separate Python objects (dicts store references, not
        inline values). This undercounts slightly (doesn't walk every
        object with full recursion) but is far more honest than
        sys.getsizeof(self._counts) alone, which only reports the
        bucket table, not the objects it points to.
        """
        total = sys.getsizeof(self._counts)
        for k, v in self._counts.items():
            total += sys.getsizeof(k) + sys.getsizeof(v)
        return total


class CountMinSketch:
    def __init__(self, width: int, depth: int):
        """
        width  — counters per row. Larger width -> fewer collisions ->
                 lower overestimation error. Error bound: expected
                 overestimate <= (total items added) / width.
        depth  — number of independent rows. Larger depth -> lower
                 probability that *all* rows collide for a given item.
                 Failure probability roughly e^(-depth) for the error
                 bound above to be exceeded.
        """
        self.width = width
        self.depth = depth
        # unsigned 32-bit counters, one flat array of width*depth
        self.table = array("I", [0] * (width * depth))
        self.n_added = 0

    def _row_indices(self, item: str):
        """
        Derive `depth` independent hash values from a single blake2b
        digest by salting with the row index, rather than needing
        `depth` genuinely different hash functions. Each row gets its
        own slice of hash output, then reduced mod width.
        """
        for row in range(self.depth):
            h = hashlib.blake2b(item.encode(), person=str(row).encode(), digest_size=8)
            idx = int.from_bytes(h.digest(), "big") % self.width
            yield row, idx

    def add(self, kmer: str):
        for row, idx in self._row_indices(kmer):
            self.table[row * self.width + idx] += 1
        self.n_added += 1

    def query(self, kmer: str) -> int:
        return min(self.table[row * self.width + idx]
                    for row, idx in self._row_indices(kmer))

    def approx_memory_bytes(self) -> int:
        """Exact, not approximate: a CMS's whole point is that its
        memory footprint is known and fixed ahead of time."""
        return self.table.itemsize * len(self.table) + sys.getsizeof(self)
