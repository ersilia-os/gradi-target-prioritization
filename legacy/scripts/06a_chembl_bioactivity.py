"""ChEMBL bioactivity by SEQUENCE mapping (docs §2.1b).

For each protein: how many molecules were tested and how many are potent (pChEMBL >= 6, i.e.
<= 1 uM) — both directly and across homologs. Crucially we map to ChEMBL targets by SEQUENCE
(DIAMOND), NOT by exact UniProt accession: HS11286 is a dark TrEMBL proteome whose accessions
rarely equal ChEMBL's (reviewed / other-strain) accessions, so accession matching silently
misses direct K. pneumoniae targets (e.g. the SHV / OXA-48 / CTX-M / NDM beta-lactamases).

Two phases:
  1. extract  — one pass over the local ChEMBL SQLite dump aggregating per target accession
     (n_compounds, n_potent, best_pchembl) + organism/tax_id + sequence. Writes a summary CSV
     and a target FASTA under data/processed/other/chembl/ (cached, organism-agnostic).
  2. map+aggregate — DIAMOND blastp the proteome vs the ChEMBL target FASTA, then per protein
     take the best hit in each bucket: direct (pident >= 95, ~same protein), bacterial/non-human,
     human (tax 9606, kept separate — family druggability, not selective). -> <prefix>_chembl.csv

ChEMBL path auto-detected under data/raw/other/chembl/. DIAMOND runs in `gradi-ortho`; this
script runs in `gradi`.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

CHEMBL_DIR = L.REPO_ROOT / "data" / "processed" / "other" / "chembl"
SUMMARY = CHEMBL_DIR / "chembl_target_summary.csv"
TARGET_FAA = CHEMBL_DIR / "chembl_targets.faa"
HUMAN_TAX = 9606
DIRECT_PIDENT = 95.0
TRANSFER_FLOOR = 40.0  # min %identity to transfer bioactivity from a homolog (below = too distant)

EXTRACT_SQL = """
SELECT cs.accession                                                       AS accession,
       MAX(td.organism)                                                   AS organism,
       MAX(td.tax_id)                                                     AS tax_id,
       MAX(oc.l1)                                                         AS superkingdom,
       cs.sequence                                                        AS sequence,
       COUNT(DISTINCT act.molregno)                                       AS n_compounds,
       COUNT(DISTINCT CASE WHEN act.pchembl_value >= 6 THEN act.molregno END) AS n_potent,
       MAX(act.pchembl_value)                                             AS best_pchembl
FROM component_sequences cs
JOIN target_components   tc ON tc.component_id = cs.component_id
JOIN target_dictionary   td ON td.tid = tc.tid
LEFT JOIN organism_class oc ON oc.tax_id = td.tax_id
JOIN assays              a  ON a.tid = td.tid
JOIN activities          act ON act.assay_id = a.assay_id
WHERE a.assay_type IN ('B', 'F') AND cs.accession IS NOT NULL AND cs.sequence IS NOT NULL
GROUP BY cs.accession;
"""


def find_db(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    hits = sorted((L.REPO_ROOT / "data" / "raw" / "other" / "chembl").rglob("chembl_*.db"))
    if not hits:
        raise FileNotFoundError("No chembl_*.db under data/raw/other/chembl/; extract the dump.")
    return hits[0]


def extract_targets(db_path: Path) -> pd.DataFrame:
    print(f"  extracting ChEMBL targets from {db_path.name} ...", flush=True)
    con = sqlite3.connect(str(db_path))
    try:
        df = pd.read_sql_query(EXTRACT_SQL, con)
    finally:
        con.close()
    CHEMBL_DIR.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["sequence"]).to_csv(SUMMARY, index=False)
    with open(TARGET_FAA, "w") as fh:
        for r in df.itertuples():
            seq = (r.sequence or "").strip().replace("*", "")
            if seq:
                fh.write(f">{r.accession}\n{seq}\n")
    print(f"  wrote {SUMMARY.name} + {TARGET_FAA.name} ({len(df)} targets)", flush=True)
    return df


def best_hit(rows: list[dict]):
    """Pick the most-informative hit: max (n_potent, best_pchembl)."""
    best = None
    for h in rows:
        key = (h["n_potent"], h["best_pchembl"] if h["best_pchembl"] is not None else 0.0)
        if best is None or key > best[0]:
            best = (key, h)
    return best[1] if best else None


def aggregate(org: str, smap: dict) -> pd.DataFrame:
    _, prefix = L.ORGANISMS[org]
    hits_tsv = L.processed_dir(org, "chembl") / "diamond_hits.tsv"
    L.run_diamond_blastp(L.proteome_fasta(org), TARGET_FAA, hits_tsv)
    hits = L.load_diamond_hits(hits_tsv)
    hits["q"] = hits["qseqid"].map(L.acc_from_header)

    by_q: dict[str, list[dict]] = {}
    for r in hits.itertuples():
        acc = r.sseqid
        s = smap.get(acc)
        if not s:
            continue
        org_name, tax_id, superkingdom, nc, npot, bp = s
        if nc <= 0:
            continue
        by_q.setdefault(r.q, []).append(
            {
                "acc": acc, "organism": org_name, "tax_id": tax_id, "superkingdom": superkingdom,
                "n_compounds": nc, "n_potent": npot, "best_pchembl": bp,
                "pident": float(r.pident),
            }
        )

    rows = []
    for acc in L.load_accessions(org):
        row = {"uniprot_accession": acc}
        hs = by_q.get(acc, [])
        # human bucket: family druggability (kept separate, not selective)
        human = [h for h in hs if h["tax_id"] == HUMAN_TAX and h["pident"] >= TRANSFER_FLOOR]
        # bacterial bucket: TRUE bacteria only (organism_class), above the transfer-identity floor
        # — excludes coincidental eukaryote / model-organism matches at low identity
        bact = [h for h in hs if h["superkingdom"] == "Bacteria" and h["pident"] >= TRANSFER_FLOOR]
        b = best_hit(bact)
        h = best_hit(human)
        direct = max((x for x in bact if x["pident"] >= DIRECT_PIDENT),
                     key=lambda x: (x["n_potent"], x["best_pchembl"] or 0), default=None)
        # direct
        row["chembl_direct_has"] = direct is not None
        row["chembl_direct_pident"] = round(direct["pident"], 1) if direct else None
        row["chembl_direct_acc"] = direct["acc"] if direct else ""
        row["chembl_direct_organism"] = direct["organism"] if direct else ""
        row["chembl_direct_n_compounds"] = direct["n_compounds"] if direct else 0
        row["chembl_direct_n_potent"] = direct["n_potent"] if direct else 0
        row["chembl_direct_best_pchembl"] = direct["best_pchembl"] if direct else None
        # bacterial / non-human (antibacterial-relevant)
        row["chembl_bact_n_compounds"] = b["n_compounds"] if b else 0
        row["chembl_bact_n_potent"] = b["n_potent"] if b else 0
        row["chembl_bact_best_pchembl"] = b["best_pchembl"] if b else None
        row["chembl_bact_best_acc"] = b["acc"] if b else ""
        row["chembl_bact_best_organism"] = b["organism"] if b else ""
        row["chembl_bact_best_pident"] = round(b["pident"], 1) if b else None
        # human (separate; family druggability, not selective)
        row["chembl_human_n_compounds"] = h["n_compounds"] if h else 0
        row["chembl_human_n_potent"] = h["n_potent"] if h else 0
        row["chembl_human_best_pchembl"] = h["best_pchembl"] if h else None
        row["chembl_human_best_acc"] = h["acc"] if h else ""
        row["chembl_human_best_pident"] = round(h["pident"], 1) if h else None
        # any = bacterial bucket
        row["chembl_any_n_compounds"] = row["chembl_bact_n_compounds"]
        row["chembl_any_n_potent"] = row["chembl_bact_n_potent"]
        row["chembl_any_best_pchembl"] = row["chembl_bact_best_pchembl"]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--db", default=None)
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    if args.refresh or not SUMMARY.exists() or not TARGET_FAA.exists():
        extract_targets(find_db(args.db))
    s = pd.read_csv(SUMMARY)
    smap = {
        r.accession: (r.organism, int(r.tax_id) if pd.notna(r.tax_id) else -1,
                      r.superkingdom if pd.notna(r.superkingdom) else "",
                      int(r.n_compounds), int(r.n_potent),
                      float(r.best_pchembl) if pd.notna(r.best_pchembl) else None)
        for r in s.itertuples()
    }

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    df = aggregate(org, smap)
    out = L.results_dir(org) / f"{prefix}_chembl.csv"
    df.to_csv(out, index=False)
    print(
        f"[{org}] wrote {out} ({len(df)} proteins; "
        f"direct(≥95%id)={int(df['chembl_direct_has'].sum())}, "
        f"bacterial potent={int((df['chembl_bact_n_potent'] > 0).sum())}, "
        f"any potent={int((df['chembl_any_n_potent'] > 0).sum())}, "
        f"human potent={int((df['chembl_human_n_potent'] > 0).sum())})",
        flush=True,
    )


if __name__ == "__main__":
    main()
