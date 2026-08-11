"""
debruijn_graph.py

Builds a colored de Bruijn graph from normal + tumor reads, and detects
"bubbles" — the graph signature of a point mutation — without ever
aligning a read to a reference.

GRAPH MODEL (and its honest limitation):
Nodes are canonical (k-1)-mers. Each k-mer observed in a read contributes
an edge between canonical(kmer[:-1]) and canonical(kmer[1:]), tagged with
which sample(s) support it (colors: "normal", "tumor", or both — this is
a *colored* de Bruijn graph, the same concept behind tools like Cortex
and the compact colored-graph representations Rayan Chikhi's group has
published on).

This is built as an UNDIRECTED graph. Canonicalizing each k-mer's two
(k-1)-mer endpoints independently correctly merges the same genomic
locus regardless of which strand a read came from (a reverse-complement
read produces the mirror-image edge, which canonicalizes to the same
undirected edge) — but it means the graph has no notion of "assembly
direction" or explicit orientation at each node. Production assemblers
track this via bidirected graphs so they can walk a path and read off
an actual sequence. This implementation deliberately does NOT do that:
it detects bubble topology (which is sufficient to flag a candidate
variant site) but does not reconstruct the mutant allele's sequence
from the graph. That's a real, substantial next step — named here
rather than glossed over.

WHY A BUBBLE MEANS A MUTATION:
A single-base substitution changes every k-mer that overlaps it (k of
them). Upstream and downstream of the mutation, tumor and normal reads
produce identical k-mers (shared, "both"-colored edges). Across the
mutation itself, they diverge — producing a short run of tumor-only
edges and a parallel run of normal-only edges, both anchored at the
same two shared nodes (the unmutated sequence flanking the mutation on
either side). That divergent-then-reconvergent shape is a "bubble".
"""

from collections import defaultdict, deque

from kmer_utils import canonical, read_fastq


class DeBruijnGraph:
    def __init__(self, k: int):
        self.k = k
        # edge key: sorted tuple of the two canonical (k-1)-mer endpoints
        # edge value: {"normal": count, "tumor": count}
        self.edge_colors = defaultdict(lambda: {"normal": 0, "tumor": 0})
        self.neighbors = defaultdict(set)  # node -> set of neighbor nodes

    def add_kmer(self, kmer: str, label: str):
        prefix, suffix = kmer[:-1], kmer[1:]
        node_a, node_b = canonical(prefix), canonical(suffix)
        if node_a == node_b:
            return  # degenerate self-loop (low-complexity repeat) — skip
        edge_key = (node_a, node_b) if node_a < node_b else (node_b, node_a)
        self.edge_colors[edge_key][label] += 1
        self.neighbors[node_a].add(node_b)
        self.neighbors[node_b].add(node_a)

    def add_reads_from_fastq(self, path, label: str):
        for _header, seq, _qual in read_fastq(path):
            seq = seq.upper()
            valid = set("ACGT")
            for i in range(len(seq) - self.k + 1):
                window = seq[i:i + self.k]
                if all(b in valid for b in window):
                    self.add_kmer(window, label)

    def prune_low_support(self, min_support: int):
        """
        Drop edges with total support below min_support. Sequencing
        errors create spurious low-count k-mers that would otherwise
        show up as fake tiny "bubbles" — this is the standard
        error-correction step real assemblers apply before bubble
        calling (usually called "tip clipping" for dead-end paths;
        this is the analogous idea applied directly to edge support).
        """
        to_remove = [
            edge for edge, colors in self.edge_colors.items()
            if (colors["normal"] + colors["tumor"]) < min_support
        ]
        for edge in to_remove:
            a, b = edge
            del self.edge_colors[edge]
            self.neighbors[a].discard(b)
            self.neighbors[b].discard(a)
        return len(to_remove)

    def classify_edges(self, min_color_support: int = 1):
        """
        Split edges into normal-only, tumor-only, and shared, based on
        which sample(s) support each (post-pruning).

        min_color_support matters: a single stray read (most likely a
        sequencing error) can produce one coincidental count for the
        "wrong" sample on an otherwise clean private edge. Treating any
        nonzero count as "supported" lets that one noisy read
        reclassify a real private edge as shared, fragmenting an
        otherwise clean bubble into a messier shape with the wrong
        number of anchors. Requiring at least min_color_support reads
        per color is a small but real correction for that.
        """
        normal_only, tumor_only, shared = set(), set(), set()
        for edge, colors in self.edge_colors.items():
            has_normal = colors["normal"] >= min_color_support
            has_tumor = colors["tumor"] >= min_color_support
            if has_normal and has_tumor:
                shared.add(edge)
            elif has_tumor:
                tumor_only.add(edge)
            elif has_normal:
                normal_only.add(edge)
        return normal_only, tumor_only, shared

    def connected_components(self, edge_subset):
        """Connected components of the subgraph induced by edge_subset,
        as node sets. Plain BFS — no need for anything fancier at this
        graph size."""
        adj = defaultdict(set)
        for a, b in edge_subset:
            adj[a].add(b)
            adj[b].add(a)

        seen = set()
        components = []
        for start in adj:
            if start in seen:
                continue
            component = set()
            queue = deque([start])
            seen.add(start)
            while queue:
                node = queue.popleft()
                component.add(node)
                for nbr in adj[node]:
                    if nbr not in seen:
                        seen.add(nbr)
                        queue.append(nbr)
            components.append(component)
        return components

    def boundary_nodes(self, component: set, min_color_support: int = 1):
        """
        Nodes in `component` that also have a SHARED edge (present in
        both samples) connecting them to a node outside the component
        — i.e. the anchor points where this private region reconnects
        to the common backbone graph.
        """
        _normal_only, _tumor_only, shared = self.classify_edges(min_color_support)
        boundary = set()
        for edge in shared:
            a, b = edge
            if a in component and b not in component:
                boundary.add(a)
            elif b in component and a not in component:
                boundary.add(b)
        return boundary

    def find_bubbles(self, min_color_support: int = 1):
        """
        Returns a list of candidate variant sites: dicts with the
        shared anchor nodes and the tumor-only / normal-only node sets
        that diverge between them. A bubble is a tumor-only component
        and a normal-only component that share the exact same pair of
        boundary (anchor) nodes.
        """
        normal_only, tumor_only, _shared = self.classify_edges(min_color_support)

        tumor_components = self.connected_components(tumor_only)
        normal_components = self.connected_components(normal_only)

        # index normal components by their frozenset of boundary anchors
        normal_by_boundary = {}
        for comp in normal_components:
            boundary = frozenset(self.boundary_nodes(comp, min_color_support))
            if len(boundary) == 2:  # a clean bubble has exactly 2 anchors
                normal_by_boundary[boundary] = comp

        bubbles = []
        for t_comp in tumor_components:
            t_boundary = frozenset(self.boundary_nodes(t_comp, min_color_support))
            if len(t_boundary) != 2:
                continue  # not a clean two-anchor divergence — skip
            if t_boundary in normal_by_boundary:
                bubbles.append({
                    "anchors": tuple(t_boundary),
                    "tumor_nodes": t_comp,
                    "normal_nodes": normal_by_boundary[t_boundary],
                })
        return bubbles
