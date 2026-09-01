"""Stage 09f — experimental localization from STEPdb, and its ortholog transfer onto K. pneumoniae.

The gap this closes is stark. 09a shows that of HS11286's 5,728 proteins exactly **one** carries an
experimentally-evidenced UniProt localization, against **803** for E. coli K-12. Every other Kp call
in the axis is therefore curated-by-similarity or predicted. E. coli, by contrast, has been mapped
exhaustively: STEPdb 2.0 (Loos et al.) assigns all ~4,300 K-12 proteins to subcellular classes,
built from proteomics, biochemistry and topology data re-examined by hand across 426 papers.

So this track does two things:

  * `--organism ecoli` — ingest STEPdb and normalise its classes onto our vocabulary. STEPdb is the
    only source here that resolves *which side* a peripheral membrane protein faces, which matters
    directly: a peripheral inner-membrane protein facing the cytoplasm is as reachable by ClpXP as a
    soluble cytoplasmic one, and is scored as such.
  * `--organism kpneumoniae` — carry E. coli's experimental calls (UniProt-experimental ∪ STEPdb)
    onto HS11286 through the 03a ortholog table. That table is OrthoFinder orthogroup membership
    with no identity scores, so multi-ortholog anchors are resolved by donor consensus (majority
    label, ties abstain) rather than by identity; donors and agreement are kept as provenance.

That is the difference between Kp localization being entirely predicted and a large slice of it
resting on E. coli bench work — the same evidence-transfer logic the essentiality axis already uses
(07c), applied to compartments, which are better conserved across Enterobacteriaceae than most
properties.

Output: output/results/<organism>/<prefix>_loc_transfer.csv
Run with the `gradi` conda env. Network fetch for E. coli (stepdb.eu), cached on disk.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

STEPDB_URL = ("http://stepdb.eu/static/downs/"
              "E%20coli%20K-12%20strain%20MG1655%20basic%20proteome.csv")
STEPDB_ACC_COL = "Accession (UniProt)"
STEPDB_LOC_COL = "STEPdb Sub-cellular Location (Full Name)"
STEPDB_REF_COL = "Annotation References"


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=3, max=60))
def fetch_stepdb(dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  using cached {dest.name}", flush=True)
        return dest
    r = requests.get(STEPDB_URL, timeout=300)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


def build_ecoli() -> pd.DataFrame:
    raw = LOC.localization_raw_dir("ecoli", "stepdb")
    src = fetch_stepdb(raw / "stepdb_k12_basic_proteome.csv")

    # Semicolon-delimited, and some free-text cells contain semicolons of their own — rows that
    # shift are handled by classify_stepdb() returning `unknown` rather than a wrong compartment.
    with open(src, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter=";"))

    canon = set(LOC.load_accessions("ecoli"))
    recs = []
    for r in rows:
        acc = (r.get(STEPDB_ACC_COL) or "").strip()
        if acc not in canon:
            continue
        label = LOC.classify_stepdb(r.get(STEPDB_LOC_COL) or "")
        if label == "unknown":
            continue
        # STEPdb mixes hand-curated and inferred assignments; a cited reference is the best
        # available proxy for the assignment resting on published bench work.
        has_ref = bool((r.get(STEPDB_REF_COL) or "").strip())
        recs.append({
            "uniprot_accession": acc,
            "transfer_localization": label,
            "transfer_evidence": "experimental" if has_ref else "curated",
            "transfer_source": "stepdb",
            "transfer_donor": "",
            "transfer_n_donors": None,
            "transfer_agreement": None,
            "stepdb_raw": (r.get(STEPDB_LOC_COL) or "").strip(),
        })
    return pd.DataFrame(recs).drop_duplicates(subset="uniprot_accession")


def build_kpneumoniae() -> pd.DataFrame:
    rdir = LOC.results_dir("ecoli")

    ec_uni = pd.read_csv(rdir / "ec_loc_uniprot.csv")
    ec_exp = ec_uni[(ec_uni["uniprot_evidence"] == "experimental")
                    & (ec_uni["uniprot_localization"] != "unknown")]
    donors: dict[str, str] = dict(zip(ec_exp["uniprot_accession"], ec_exp["uniprot_localization"]))
    print(f"  {len(donors)} E. coli UniProt-experimental donors", flush=True)

    step_path = rdir / "ec_loc_transfer.csv"
    if step_path.exists():
        step = pd.read_csv(step_path)
        step = step[step["transfer_evidence"] == "experimental"]
        # UniProt-experimental wins where both speak; STEPdb fills the rest.
        added = 0
        for acc, label in zip(step["uniprot_accession"], step["transfer_localization"]):
            if acc not in donors:
                donors[acc] = label
                added += 1
        print(f"  +{added} STEPdb-experimental donors -> {len(donors)} total", flush=True)
    else:
        print("  [warn] ec_loc_transfer.csv missing — run `--organism ecoli` first for STEPdb "
              "donors; transferring from UniProt-experimental only", flush=True)

    transferred = LOC.transfer_categorical_ecoli_to_kp(donors)
    recs = [{
        "uniprot_accession": acc,
        "transfer_localization": rec["label"],
        "transfer_evidence": "experimental",
        "transfer_source": "ortholog_transfer",
        "transfer_donor": rec["donors"],
        "transfer_n_donors": rec["n_donors"],
        "transfer_agreement": rec["agreement"],
        "stepdb_raw": "",
    } for acc, rec in transferred.items()]
    return pd.DataFrame(recs)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    print(f"[{org}] building experimental-localization transfer table", flush=True)

    df = build_ecoli() if org == "ecoli" else build_kpneumoniae()

    out = LOC.results_dir(org) / f"{prefix}_loc_transfer.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(df)} proteins)", flush=True)
    if len(df):
        print(f"[{org}] classes: {df['transfer_localization'].value_counts().to_dict()}", flush=True)
        print(f"[{org}] evidence: {df['transfer_evidence'].value_counts().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
