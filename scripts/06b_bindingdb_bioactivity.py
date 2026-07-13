"""BindingDB bioactivity by SEQUENCE mapping (docs §2.1b).

Like 06a, but for BindingDB and mapped by target SEQUENCE rather than accession, so direct
K. pneumoniae entries deposited under other-strain / reviewed accessions are still found
(BindingDB has thousands of rows for Kp beta-lactamases under accessions absent from HS11286).

Phases (cached):
  1. build  — stream BindingDB_All.tsv once; collect each unique target-chain sequence -> a FASTA
              + meta (organism, representative accession) under data/processed/other/bindingdb/.
  2. map    — DIAMOND blastp each proteome vs that FASTA (gradi-ortho).
  3. tally  — stream the TSV again, aggregating distinct ligand InChIKeys tested / potent (<=1 uM,
              pAff = -log10(M) >= 6) / best pAff, ONLY for target sequences that matched a proteome.
  4. aggregate per protein -> output/results/<org>/<prefix>_bindingdb.csv (direct / bacterial /
     human buckets, %identity reported). Run in `gradi` (DIAMOND in `gradi-ortho`).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

csv.field_size_limit(50_000_000)

BDB_DIR = L.REPO_ROOT / "data" / "raw" / "other" / "bindingdb"
OUT_DIR = L.REPO_ROOT / "data" / "processed" / "other" / "bindingdb"
TARGET_FAA = OUT_DIR / "bindingdb_targets.faa"
META = OUT_DIR / "bindingdb_target_meta.csv"
SUMMARY = OUT_DIR / "bindingdb_target_summary.csv"

SEQ_COL = "BindingDB Target Chain Sequence"
ORG_COL = "Target Source Organism According to Curator or DataSource"
IK_COL = "Ligand InChI Key"
ACC_COL = "UniProt (SwissProt) Primary ID of Target Chain"
AFF_COLS = ["Ki (nM)", "IC50 (nM)", "Kd (nM)", "EC50 (nM)"]
DIRECT_PIDENT = 95.0
TRANSFER_FLOOR = 40.0
CHEMBL_SUMMARY = L.REPO_ROOT / "data" / "processed" / "other" / "chembl" / "chembl_target_summary.csv"


def bacterial_genera() -> set[str]:
    """Genus names classified as Bacteria by ChEMBL's organism_class (built by 06a). Used to
    classify BindingDB targets (BindingDB has no taxonomy id, only an organism string)."""
    if not CHEMBL_SUMMARY.exists():
        print("  [warn] chembl_target_summary.csv missing; run 06a first for bacterial classification", flush=True)
        return set()
    s = pd.read_csv(CHEMBL_SUMMARY)
    if "superkingdom" not in s.columns:
        return set()
    orgs = s.loc[s["superkingdom"] == "Bacteria", "organism"].dropna().astype(str)
    return {o.split()[0].lower() for o in orgs if o.split()}


def find_tsv(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    hits = sorted(BDB_DIR.glob("BindingDB_All*.tsv"))
    if not hits:
        raise FileNotFoundError(f"No BindingDB_All*.tsv under {BDB_DIR}")
    return hits[0]


def seq_id(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()[:16]


def parse_nm(val: str):
    if not val:
        return None, False
    v = val.strip()
    qual = ""
    while v and v[0] in "<>=~ ":
        qual += v[0]
        v = v[1:]
    try:
        num = float(v)
    except ValueError:
        return None, False
    if num <= 0:
        return None, False
    return num, (">" not in qual)


def col_index(header: list[str], name: str) -> int | None:
    return header.index(name) if name in header else None


def build_targets(tsv: Path) -> None:
    """Pass 1: one FASTA + meta row per unique target-chain sequence."""
    print(f"  [pass 1] scanning {tsv.name} for unique target sequences ...", flush=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    seen: dict[str, tuple[str, str]] = {}
    with open(tsv, newline="", encoding="utf-8", errors="replace") as fh, open(TARGET_FAA, "w") as faa:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        si, oi, ai = col_index(header, SEQ_COL), col_index(header, ORG_COL), col_index(header, ACC_COL)
        n = 0
        for row in reader:
            n += 1
            if n % 1_000_000 == 0:
                print(f"    {n:,} rows, {len(seen)} targets", flush=True)
            if si is None or si >= len(row):
                continue
            seq = row[si].strip().replace("*", "").upper()
            if len(seq) < 30 or not seq.isalpha():
                continue
            sid = seq_id(seq)
            if sid not in seen:
                org = row[oi].strip() if (oi is not None and oi < len(row)) else ""
                acc = row[ai].strip() if (ai is not None and ai < len(row)) else ""
                seen[sid] = (org, acc)
                faa.write(f">{sid}\n{seq}\n")
    pd.DataFrame(
        [{"sid": k, "organism": v[0], "accession": v[1]} for k, v in seen.items()]
    ).to_csv(META, index=False)
    print(f"  [pass 1] {len(seen)} unique targets -> {TARGET_FAA.name}", flush=True)


def tally(tsv: Path, matched: set[str]) -> None:
    """Pass 2: distinct tested / potent InChIKeys + best pAff per matched target sequence."""
    print(f"  [pass 2] tallying bioactivity for {len(matched)} matched targets ...", flush=True)
    tested: dict[str, set[str]] = {}
    potent: dict[str, set[str]] = {}
    best: dict[str, float] = {}
    with open(tsv, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = next(reader)
        si = col_index(header, SEQ_COL)
        iki = col_index(header, IK_COL)
        affi = [col_index(header, c) for c in AFF_COLS]
        n = 0
        for row in reader:
            n += 1
            if n % 1_000_000 == 0:
                print(f"    {n:,} rows", flush=True)
            if si is None or si >= len(row):
                continue
            seq = row[si].strip().replace("*", "").upper()
            if len(seq) < 30:
                continue
            sid = seq_id(seq)
            if sid not in matched:
                continue
            ik = row[iki].strip() if (iki is not None and iki < len(row)) else ""
            if not ik:
                continue
            best_nm, can_potent = None, False
            for ci in affi:
                if ci is not None and ci < len(row):
                    nm, ok = parse_nm(row[ci])
                    if nm is not None and (best_nm is None or nm < best_nm):
                        best_nm, can_potent = nm, ok
            tested.setdefault(sid, set()).add(ik)
            if best_nm is not None:
                paff = 9.0 - math.log10(best_nm)
                if paff > best.get(sid, -99):
                    best[sid] = paff
                if can_potent and best_nm <= 1000.0:
                    potent.setdefault(sid, set()).add(ik)
    rows = [
        {"sid": sid, "n_compounds": len(tested[sid]),
         "n_potent": len(potent.get(sid, set())), "best_paff": round(best.get(sid, float("nan")), 3)}
        for sid in tested
    ]
    pd.DataFrame(rows).to_csv(SUMMARY, index=False)
    print(f"  [pass 2] wrote {SUMMARY.name} ({len(rows)} targets)", flush=True)


def best_hit(rows: list[dict]):
    best = None
    for h in rows:
        key = (h["n_potent"], h["best_paff"] if h["best_paff"] is not None else 0.0)
        if best is None or key > best[0]:
            best = (key, h)
    return best[1] if best else None


def aggregate(org: str, smap: dict, meta: dict, bact_genera: set[str]) -> pd.DataFrame:
    _, prefix = L.ORGANISMS[org]
    hits = L.load_diamond_hits(L.processed_dir(org, "bindingdb") / "diamond_hits.tsv")
    hits["q"] = hits["qseqid"].map(L.acc_from_header)
    by_q: dict[str, list[dict]] = {}
    for r in hits.itertuples():
        sid = r.sseqid
        s = smap.get(sid)
        if not s or s[0] <= 0:
            continue
        nc, npot, bp = s
        organism, acc = meta.get(sid, ("", ""))
        genus = organism.split()[0].lower() if organism.split() else ""
        by_q.setdefault(r.q, []).append(
            {"sid": sid, "organism": organism, "acc": acc,
             "is_human": "homo sapiens" in organism.lower(),
             "is_bact": genus in bact_genera,
             "n_compounds": nc, "n_potent": npot, "best_paff": bp, "pident": float(r.pident)}
        )

    rows = []
    for acc in L.load_accessions(org):
        row = {"uniprot_accession": acc}
        hs = by_q.get(acc, [])
        # bacterial = true bacteria (genus per ChEMBL organism_class) above the transfer floor
        bact = [h for h in hs if h["is_bact"] and h["pident"] >= TRANSFER_FLOOR]
        human = [h for h in hs if h["is_human"] and h["pident"] >= TRANSFER_FLOOR]
        b, h = best_hit(bact), best_hit(human)
        direct = max((x for x in bact if x["pident"] >= DIRECT_PIDENT),
                     key=lambda x: (x["n_potent"], x["best_paff"] or 0), default=None)
        row["bindingdb_direct_has"] = direct is not None
        row["bindingdb_direct_pident"] = round(direct["pident"], 1) if direct else None
        row["bindingdb_direct_acc"] = direct["acc"] if direct else ""
        row["bindingdb_direct_organism"] = direct["organism"] if direct else ""
        row["bindingdb_direct_n_compounds"] = direct["n_compounds"] if direct else 0
        row["bindingdb_direct_n_potent"] = direct["n_potent"] if direct else 0
        row["bindingdb_direct_best_paff"] = direct["best_paff"] if direct else None
        row["bindingdb_bact_n_compounds"] = b["n_compounds"] if b else 0
        row["bindingdb_bact_n_potent"] = b["n_potent"] if b else 0
        row["bindingdb_bact_best_paff"] = b["best_paff"] if b else None
        row["bindingdb_bact_best_acc"] = b["acc"] if b else ""
        row["bindingdb_bact_best_organism"] = b["organism"] if b else ""
        row["bindingdb_bact_best_pident"] = round(b["pident"], 1) if b else None
        row["bindingdb_human_n_compounds"] = h["n_compounds"] if h else 0
        row["bindingdb_human_n_potent"] = h["n_potent"] if h else 0
        row["bindingdb_human_best_paff"] = h["best_paff"] if h else None
        row["bindingdb_human_best_acc"] = h["acc"] if h else ""
        row["bindingdb_human_best_pident"] = round(h["pident"], 1) if h else None
        row["bindingdb_any_n_compounds"] = row["bindingdb_bact_n_compounds"]
        row["bindingdb_any_n_potent"] = row["bindingdb_bact_n_potent"]
        row["bindingdb_any_best_paff"] = row["bindingdb_bact_best_paff"]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--tsv", default=None)
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    tsv = find_tsv(args.tsv)

    if args.refresh or not TARGET_FAA.exists() or not META.exists():
        build_targets(tsv)
    # DIAMOND for BOTH organisms (so the shared summary covers either), then tally once
    if args.refresh or not SUMMARY.exists():
        matched: set[str] = set()
        for o in L.ORGANISMS:
            ht = L.processed_dir(o, "bindingdb") / "diamond_hits.tsv"
            L.run_diamond_blastp(L.proteome_fasta(o), TARGET_FAA, ht)
            matched.update(L.load_diamond_hits(ht)["sseqid"].astype(str))
        tally(tsv, matched)
    else:
        L.run_diamond_blastp(L.proteome_fasta(args.organism), TARGET_FAA,
                             L.processed_dir(args.organism, "bindingdb") / "diamond_hits.tsv")

    s = pd.read_csv(SUMMARY)
    smap = {r.sid: (int(r.n_compounds), int(r.n_potent),
                    float(r.best_paff) if pd.notna(r.best_paff) else None) for r in s.itertuples()}
    m = pd.read_csv(META).fillna("")
    meta = {r.sid: (str(r.organism), str(r.accession)) for r in m.itertuples()}

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    df = aggregate(org, smap, meta, bacterial_genera())
    out = L.results_dir(org) / f"{prefix}_bindingdb.csv"
    df.to_csv(out, index=False)
    print(
        f"[{org}] wrote {out} ({len(df)} proteins; "
        f"direct(≥95%id)={int(df['bindingdb_direct_has'].sum())}, "
        f"bacterial potent={int((df['bindingdb_bact_n_potent'] > 0).sum())}, "
        f"any potent={int((df['bindingdb_any_n_potent'] > 0).sum())}, "
        f"human potent={int((df['bindingdb_human_n_potent'] > 0).sum())})",
        flush=True,
    )


if __name__ == "__main__":
    main()
