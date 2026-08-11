"""
download_reference.py

Fetches real RefSeq CDS (coding sequence) FASTA files from NCBI for a small
panel of EMT/IDR genes, using NCBI E-utilities (efetch).

NOTE ON WHERE TO RUN THIS:
This script needs open internet access to eutils.ncbi.nlm.nih.gov.
Run it on your own laptop, not inside a locked-down sandbox/CI environment.

Usage:
    python src/download_reference.py
    python src/download_reference.py --genes SNAI1 ZEB1 TWIST1

To add a gene, look up its RefSeq mRNA accession (format NM_XXXXXX) on
https://www.ncbi.nlm.nih.gov/gene/ and add it to GENE_PANEL below.
"""

import argparse
import json
import time
import urllib.request
import urllib.parse
from pathlib import Path

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Verified RefSeq mRNA accessions for a starter panel of EMT/IDR genes.
# Extend this as you pull more genes from your pan-cancer IDR catalog.
GENE_PANEL = {
    "SNAI1": "NM_005985",   # Snail — canonical EMT-TF, well-characterized IDR
    "ZEB1": "NM_030751",    # Zinc finger E-box binding homeobox 1
    # "TWIST1": "NM_000474",  # add once you've verified the accession
}

REFERENCE_DIR = Path(__file__).resolve().parent.parent / "data" / "reference"


def fetch_cds_fasta(accession: str, tool: str = "kmer-somatic-signatures",
                     email: str = "your_email@example.com") -> str:
    """
    Fetch the CDS-only FASTA for a RefSeq mRNA accession via NCBI efetch.
    rettype=fasta_cds_na returns just the coding sequence (no 5'/3' UTR),
    which is what we want for codon-level mutation injection.
    """
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta_cds_na",
        "retmode": "text",
        "tool": tool,
        "email": email,
    }
    url = f"{EUTILS_BASE}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genes", nargs="+", default=list(GENE_PANEL.keys()),
        help="Gene symbols to fetch (must be keys in GENE_PANEL)",
    )
    args = parser.parse_args()

    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    manifest = {}

    for gene in args.genes:
        if gene not in GENE_PANEL:
            print(f"[skip] {gene}: no accession in GENE_PANEL, add it first")
            continue

        accession = GENE_PANEL[gene]
        print(f"[fetch] {gene} ({accession}) ...")
        try:
            fasta_text = fetch_cds_fasta(accession)
        except Exception as e:
            print(f"[error] {gene}: {e}")
            continue

        if not fasta_text.strip().startswith(">"):
            print(f"[error] {gene}: unexpected response, skipping. "
                  f"First 200 chars: {fasta_text[:200]!r}")
            continue

        out_path = REFERENCE_DIR / f"{gene}_cds.fasta"
        out_path.write_text(fasta_text)
        seq_len = len("".join(
            line.strip() for line in fasta_text.splitlines()
            if not line.startswith(">")
        ))
        print(f"[ok] {gene}: wrote {out_path} ({seq_len} bp)")
        manifest[gene] = {"accession": accession, "length_bp": seq_len,
                           "path": str(out_path.relative_to(REFERENCE_DIR.parent.parent))}

        time.sleep(0.4)  # be polite to NCBI's rate limits (max ~3 req/sec)

    manifest_path = REFERENCE_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nWrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
