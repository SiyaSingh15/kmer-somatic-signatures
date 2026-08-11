# Reference-Free Somatic Mutation Detection from k-mer Spectra

A prototype pipeline that detects somatic point mutations directly from
raw sequencing reads via k-mer/de Bruijn graph structure — no reference
alignment step — then layers a machine-learning classifier on top.

## Why this project

Standard somatic variant calling aligns reads to a reference genome first
(BWA, samtools, bcftools — see my [KRAS-targeted-NGS-pipeline](https://github.com/SiyaSingh15/KRAS-targeted-NGS-pipeline)
for that approach). This project instead asks: how much can be recovered
about a mutation from k-mer counts and de Bruijn graph topology alone,
without ever mapping a read to a reference coordinate?

This is not a novel idea — it's a simplified, from-scratch implementation
of the general approach behind tools like **DiscoSnp++** (bubble-calling
in de Bruijn graphs) and k-mer-based association methods like **HAWK**.
The contribution here is implementing the core mechanics myself and
applying them to a gene panel and mutation set grounded in my own
pan-cancer IDR mutation work ([IDR-EMT-PanCancer](https://github.com/SiyaSingh15/IDR-EMT-PanCancer)),
not inventing a new method.

## What's implemented from scratch vs. what uses existing tools

| Component | From scratch | Existing tool |
|---|---|---|
| Reference fetching | | NCBI E-utilities |
| Mutation injection + codon translation | ✅ hand-written codon table & logic | |
| Read simulation | | ART, with a pure-Python fallback for platforms without it |
| Canonical k-mer extraction | ✅ | |
| Exact k-mer counter | ✅ | |
| Count-Min Sketch (approximate frequency counting) | ✅ | |
| Colored de Bruijn graph + bubble calling | ✅ | |
| Alignment-based comparison baseline | | BWA-MEM, samtools, bcftools |
| ML classification layer | | XGBoost, SHAP (Step 5) |

## Project status

**Step 4: alignment-based comparison — harness done, real-data validation pending.**

- `src/compare_to_alignment_based_calling.py` — runs the same tumor
  reads through a standard BWA-MEM + bcftools pipeline (same tools as
  [KRAS-targeted-NGS-pipeline](https://github.com/SiyaSingh15/KRAS-targeted-NGS-pipeline))
  and reports precision/recall against ground truth, alongside the
  same metrics for bubble calling — the direct head-to-head.
- `run_step4_compare.sh` — runs it end to end.

**Result on placeholder data:** both approaches tie, 4/4 recovered,
0 false positives each. That's the expected, honest outcome on a
single small non-repetitive reference — it doesn't show whether
either approach has a real edge. The interesting comparison needs
real gene-scale data (repeats, more complex local sequence context),
which is the other half of this step, still pending:

**Still needed — real data swap (requires your local machine + your
own mutation catalog, not something I can do from here):**
1. Run `python src/download_reference.py` locally to fetch real
   SNAI1/ZEB1 CDS sequences (my build sandbox can't reach NCBI).
2. Replace `example_mutations.csv` with real calls from your
   IDR-EMT-PanCancer results, in the same
   `gene,cds_position,ref_base,alt_base,source` format.
3. Rerun `run_step1.sh` → `run_step3.sh` → `run_step4_compare.sh`
   against the real gene(s) instead of the placeholder.

**Step 3: de Bruijn graph + bubble calling — done, validated against ground truth.**

- `src/debruijn_graph.py` — builds an undirected, colored de Bruijn
  graph (nodes = canonical (k-1)-mers, edges = k-mers tagged with which
  sample — normal, tumor, or both — supports them), prunes low-support
  edges (sequencing-error noise), and detects bubbles: a tumor-only
  component and a normal-only component that share the same two
  "anchor" nodes, which is the graph signature a single point mutation
  leaves behind.
- `src/call_variants_from_bubbles.py` — orchestrates graph construction
  and bubble detection, then validates the result against the known
  injected mutations (available here because this is simulated data —
  a real reference-free caller wouldn't have this luxury, which is
  exactly why it's worth checking hard against ground truth now).
- `run_step3.sh` — runs it end to end.

**Result:** 4/4 injected mutations recovered as bubbles, 0 false
positives, on both the ART and pure-Python read-simulation backends.

**A real bug turned up during validation, worth keeping in the
writeup:** the first version only recovered 3/4 mutations. A single
stray read — almost certainly a sequencing error — coincidentally
produced a "wrong-sample" k-mer right at one mutation site, which was
enough to make `classify_edges` treat an otherwise-clean tumor-only
region as "shared," fragmenting the bubble. Fixed by requiring a
minimum read count per sample (not just nonzero presence) before
calling an edge sample-supported (`--min-color-support`, default 2).
This is the kind of failure that only shows up when you actually check
against ground truth rather than trust a plausible-looking result —
worth being upfront about rather than only showing the clean 4/4.

**Known limitation, stated rather than hidden:** this builds an
*undirected* graph. Canonicalizing each k-mer's two (k-1)-mer endpoints
correctly merges the same locus regardless of which strand a read came
from, but the graph has no explicit orientation at each node — so it
detects bubble *topology* (sufficient to flag a candidate variant site)
but doesn't reconstruct the mutant allele's actual sequence from a
graph walk. Production assemblers handle this via bidirected graphs,
which is real complexity (and, not coincidentally, close to what
Rayan Chikhi's group has published on for compact colored de Bruijn
graph representations) — named here as a real next step, not glossed
over.

**Step 2: k-mer counting engine — done.**

- `src/kmer_utils.py` — canonical k-mer extraction (sliding window +
  reverse-complement canonicalization, so a k-mer and its reverse
  complement always collapse to the same bucket) and a minimal FASTQ
  parser. Shared by the counter and, in step 3, the de Bruijn graph.
- `src/kmer_structures.py` — `ExactCounter` (a dict wrapper, for
  comparison) and `CountMinSketch`, a fixed-memory approximate
  frequency counter built from a flat array of counters and
  `hashlib.blake2b`-derived hash rows (chosen for determinism across
  runs, since Python's built-in `hash()` is randomized per process).
- `src/benchmark_kmer_structures.py` — runs both structures on the same
  reads and reports memory, speed, and per-k-mer accuracy.
- `run_step2.sh` — runs the benchmark against whatever's in `data/reads/`.

**Result on the current placeholder data** (k=21, CMS width=2000,
depth=4, 1,513 distinct k-mers across both FASTQs): CMS matched the
exact count on 47/50 sampled k-mers, with a max overestimation of 2 —
and by construction, it never *under*estimates. Memory footprint was
~5.9x smaller by object accounting, though at this toy gene-panel
scale that gap is modest; it becomes the whole point at real
sequencing-depth scale, where a dict of billions of distinct k-mers
just isn't an option. One honest tradeoff surfaced by the benchmark:
CMS was actually *slower* here than the dict, because `blake2b` (a
cryptographic hash, used for determinism) costs more per call than the
non-cryptographic hashes (xxHash, MurmurHash) production k-mer tools
use for this reason — a natural next optimization, not a flaw to
paper over.

**Step 1: grounded read simulation — done.**

- `src/download_reference.py` — fetches real RefSeq CDS FASTA for a panel
  of EMT/IDR genes (SNAI1, ZEB1 so far) via NCBI E-utilities. **Run this on
  your own machine** — it needs open internet access that a locked-down
  build sandbox won't have.
- `src/mutate_reference.py` — injects point mutations from a CSV table
  into a reference CDS and reports codon-level amino acid consequences
  (missense/nonsense/synonymous), using a manually written codon table.
- `src/simulate_reads.py` — wraps ART to simulate normal and tumor
  short-read FASTQ sets from the unmutated and mutated references.
  **Falls back automatically to a built-in pure-Python simulator if
  `art_illumina` isn't on PATH** (e.g. native Windows, where ART has no
  straightforward binary — WSL/Linux/Mac all have it via `apt`/`conda`).
  The fallback uses a simpler uniform substitution-error model rather
  than ART's empirical Illumina profiles; documented as the simpler
  option, not treated as equivalent. Force it explicitly with
  `--force-pure-python` if you have ART installed but want to compare.
- `run_step1.sh` — runs the above two scripts end to end.

**Currently using placeholder data.** `data/reference/PLACEHOLDER_synthetic_test_gene.fasta`
is a synthetic 603 bp ORF (clearly labeled, not a real gene) used to
validate the pipeline runs correctly. `data/mutations/example_mutations.csv`
has 4 illustrative mutations positioned against that placeholder.

To switch to real data:
1. Run `python src/download_reference.py` locally to fetch real SNAI1/ZEB1
   CDS sequences into `data/reference/`.
2. Replace `example_mutations.csv` with mutations from your own
   IDR-EMT-PanCancer results, in the same `gene,cds_position,ref_base,alt_base,source`
   format (CDS-nucleotide coordinates, 1-based from the ATG).

**Next: finish Step 4 with real data** (see above), then **Step 5** —
layer XGBoost + SHAP on top of k-mer/bubble features for tumor vs
normal classification.

## Setup

```bash
pip install -r requirements.txt   # currently empty — steps 1-3 use only the standard library
sudo apt-get install art-nextgen-simulation-tools   # or your platform's equivalent
```

## Repo layout

```
kmer-somatic-signatures/
├── data/
│   ├── reference/     # reference + mutated ("tumor") CDS FASTA files
│   ├── mutations/     # mutation tables (CSV)
│   └── reads/         # simulated FASTQ (git-ignored, regenerate via run_step1.sh)
├── src/
│   ├── download_reference.py
│   ├── mutate_reference.py
│   ├── simulate_reads.py
│   ├── kmer_utils.py
│   ├── kmer_structures.py
│   ├── benchmark_kmer_structures.py
│   ├── debruijn_graph.py
│   ├── call_variants_from_bubbles.py
│   └── compare_to_alignment_based_calling.py
├── run_step1.sh
├── run_step2.sh
├── run_step3.sh
├── run_step4_compare.sh
└── README.md
```

## Author

Siya Singh | IISER Tirupati
