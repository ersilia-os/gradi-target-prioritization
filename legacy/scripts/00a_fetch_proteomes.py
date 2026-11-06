"""Fetch the project's reference proteomes from UniProt (FASTA + TSV), keyed by UniProt accession.

A registry-driven generalization of the original Kp-only fetcher: the only things that differ
between proteomes are the UniProt proteome ID and an output label, so they live in one REGISTRY
rather than three near-identical scripts.

  - HS11286  (UP000007841) — K. pneumoniae anchor strain, 5,728 proteins (~all unreviewed, so
               NOT filtered to reviewed; see CLAUDE.md). Every downstream track derives from it.
  - EcoliK12 (UP000000625) — E. coli K-12 MG1655, 4,403 proteins (100% Swiss-Prot reviewed); the
               curation hub used for ortholog / literature transfer.
  - Human    (UP000005640) — fetched as reviewed canonical (~20,416, one sequence per gene), NOT
               the ~83k set that is bloated with unreviewed TrEMBL fragments; used for the
               selectivity-vs-human axis and ortholog detection.

Each proteome yields two files under data/raw/<organism>/proteome/ (eosvc-tracked, not Git):
  <PROTEOME_ID>_<LABEL>.fasta : sequences (accession in each header)
  <PROTEOME_ID>_<LABEL>.tsv   : accession, gene_names, sequence

Downloads use the UniProt REST stream endpoint:
  https://rest.uniprot.org/uniprotkb/stream?query=proteome:<ID>[ AND reviewed:true]&format=<fmt>

Run with the `gradi` conda env interpreter, e.g.:
  python scripts/00_fetch_proteomes.py                 # all three
  python scripts/00_fetch_proteomes.py --only EcoliK12 Human
"""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.parse import quote

import requests
from tenacity import retry, stop_after_attempt, wait_exponential

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
TSV_FIELDS = "accession,gene_names,sequence"

# label -> (UniProt proteome ID, restrict to reviewed entries?, organism folder)
REGISTRY: dict[str, dict] = {
    "HS11286": {
        "id": "UP000007841",
        "reviewed_only": False,
        "organism": "kpneumoniae",
    },  # Kp ~all unreviewed — must NOT filter
    "EcoliK12": {
        "id": "UP000000625",
        "reviewed_only": False,
        "organism": "ecoli",
    },  # 100% reviewed anyway
    "Human": {
        "id": "UP000005640",
        "reviewed_only": True,
        "organism": "human",
    },  # 20,416 canonical; avoid 83k bloat
}

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = (
    REPO_ROOT / "data" / "raw"
)  # each proteome lands in data/raw/<organism>/proteome/


def stream_url(
    proteome_id: str, fmt: str, reviewed_only: bool, fields: str | None = None
) -> str:
    query = f"proteome:{proteome_id}"
    if reviewed_only:
        query += " AND reviewed:true"
    url = f"{STREAM}?compressed=false&format={fmt}&query={quote(query, safe=':')}"
    if fields:
        url += f"&fields={fields}"
    return url


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def download(url: str) -> str:
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    return resp.text


def fetch_one(label: str) -> None:
    spec = REGISTRY[label]
    pid, reviewed = spec["id"], spec["reviewed_only"]
    out_dir = RAW_DIR / spec["organism"] / "proteome"
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{label} ({pid}{', reviewed' if reviewed else ''})"
    fasta_path = out_dir / f"{pid}_{label}.fasta"
    tsv_path = out_dir / f"{pid}_{label}.tsv"

    print(f"Downloading {tag} FASTA from UniProt ...")
    fasta = download(stream_url(pid, "fasta", reviewed))
    fasta_path.write_text(fasta)
    n_seqs = fasta.count("\n>") + (1 if fasta.startswith(">") else 0)
    print(f"  wrote {n_seqs} sequences to {fasta_path.relative_to(REPO_ROOT)}")

    print(f"Downloading {tag} TSV ({TSV_FIELDS}) from UniProt ...")
    tsv = download(stream_url(pid, "tsv", reviewed, TSV_FIELDS))
    tsv_path.write_text(tsv)
    n_rows = tsv.count("\n") - 1  # minus header line
    print(f"  wrote {n_rows} rows to {tsv_path.relative_to(REPO_ROOT)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--only",
        nargs="+",
        choices=list(REGISTRY),
        metavar="LABEL",
        help="fetch only these proteomes (default: all)",
    )
    args = ap.parse_args()

    for label in args.only or list(REGISTRY):
        fetch_one(label)


if __name__ == "__main__":
    main()
