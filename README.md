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
| Read simulation | | ART (`art_illumina`) |
| k-mer counting engine | ✅ (Step 2) | |
| Compact k-mer data structure (Bloom filter / Count-Min Sketch) | ✅ (Step 2) | |
| de Bruijn graph + bubble calling | ✅ (Step 3) | |
| ML classification layer | | XGBoost, SHAP (Step 5) |

## Project status

**Step 1 (this commit): grounded read simulation — done.**

- `src/download_reference.py` — fetches real RefSeq CDS FASTA for a panel
  of EMT/IDR genes (SNAI1, ZEB1 so far) via NCBI E-utilities. **Run this on
  your own machine** — it needs open internet access that a locked-down
  build sandbox won't have.
- `src/mutate_reference.py` — injects point mutations from a CSV table
  into a reference CDS and reports codon-level amino acid consequences
  (missense/nonsense/synonymous), using a manually written codon table.
- `src/simulate_reads.py` — wraps ART to simulate normal and tumor
  short-read FASTQ sets from the unmutated and mutated references.
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

**Next: Step 2** — a from-scratch k-mer counting engine (canonical k-mers,
hash-based exact counting) plus a compact approximate structure (Bloom
filter or Count-Min Sketch), benchmarked against each other for memory
and speed.

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
│   └── simulate_reads.py
├── run_step1.sh
└── README.md
```

## Author

Siya Singh | IISER Tirupati
