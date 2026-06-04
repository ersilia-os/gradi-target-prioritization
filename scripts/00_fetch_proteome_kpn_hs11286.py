"""Fetch the K. pneumoniae HS11286 reference proteome from UniProt (FASTA + TSV).

HS11286 (UniProt proteome UP000007841) is the only K. pneumoniae proteome flagged
by UniProt as a "Reference and representative proteome" (5,728 proteins). It is the
project anchor strain: all orthology and identifier mappings are resolved onto it
(see CLAUDE.md). UniProt accessions are the canonical identifier throughout this repo.

Two files are written, both keyed by UniProt accession:
  - <PROTEOME_ID>_HS11286.fasta : sequences (accession in each header)
  - <PROTEOME_ID>_HS11286.tsv   : accession, gene_names, sequence
The TSV schema matches the pre-staged data/raw/klebsiella_pneumoniae_proteome.tsv
so the two are interchangeable.

Both downloads use the UniProt REST stream endpoint:
  https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000007841&format=<fmt>
Output lands under data/raw/proteome/, which is tracked by eosvc (DVC + S3), not Git.
"""

from __future__ import annotations

from pathlib import Path

import requests
from tenacity import retry, stop_after_attempt, wait_exponential

PROTEOME_ID = "UP000007841"  # K. pneumoniae subsp. pneumoniae HS11286
STREAM = "https://rest.uniprot.org/uniprotkb/stream"
TSV_FIELDS = "accession,gene_names,sequence"

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "data" / "raw" / "proteome"
FASTA_PATH = OUT_DIR / f"{PROTEOME_ID}_HS11286.fasta"
TSV_PATH = OUT_DIR / f"{PROTEOME_ID}_HS11286.tsv"


def stream_url(fmt: str, fields: str | None = None) -> str:
    url = f"{STREAM}?compressed=false&format={fmt}&query=proteome:{PROTEOME_ID}"
    if fields:
        url += f"&fields={fields}"
    return url


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def download(url: str) -> str:
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    return resp.text


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {PROTEOME_ID} (HS11286) FASTA from UniProt ...")
    fasta = download(stream_url("fasta"))
    FASTA_PATH.write_text(fasta)
    n_seqs = fasta.count("\n>") + (1 if fasta.startswith(">") else 0)
    print(f"  wrote {n_seqs} sequences to {FASTA_PATH.relative_to(REPO_ROOT)}")

    print(f"Downloading {PROTEOME_ID} (HS11286) TSV ({TSV_FIELDS}) from UniProt ...")
    tsv = download(stream_url("tsv", TSV_FIELDS))
    TSV_PATH.write_text(tsv)
    n_rows = tsv.count("\n") - 1  # minus header line
    print(f"  wrote {n_rows} rows to {TSV_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
