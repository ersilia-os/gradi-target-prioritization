"""E. coli essentiality transfer + broad-spectrum conservation (docs §4.2a + §4.2c).

Lifts the well-characterised *E. coli* essential-gene set onto each anchor protein via ortholog,
and attaches a graded cross-Enterobacteriaceae "how broadly essential is this gene" score.

Data source: the Gardner-BinfLab **Enterobacteriaceae-TraDIS** compendium (Goodall/Gardner,
PMID 39207104), staged by 07a at data/raw/other/essentiality/enterobacteriaceae_tradis/. Its
`giant-tab_final.tsv` is a fully ortholog-cluster-aligned essentiality matrix keyed by the E. coli
Keio b-number (`Locus: Escherichia coli BW25113 (Keio)`), carrying:
  * `EcoGene Essentiality: Escherichia coli BW25113` — the curated binary essential call (299 genes;
    the Keio/EcoGene reference set) -> our strict §4.2a `ec_transfer_essential`;
  * `Enterobacteriaceae %essential` / `Bacteria %essential` — the fraction of genomes in which the
    gene is essential -> a graded broad-spectrum vulnerability score (§4.2c).

Mapping onto the anchor:
  * `--organism ecoli`  — the anchor IS E. coli K-12; map each accession to its b-number directly
    (proteome TSV Gene Names) and read the call. This is a *direct* essentiality call, not a transfer.
  * `--organism kpneumoniae` — map each Kp protein to its *E. coli K-12 ortholog* (the 03a orthology
    table, species Ecoli_K12_MG1655), then that ortholog's b-number -> the compendium row. Where a Kp
    protein has several E. coli orthologs, keep the most-essential one (essential > higher %essential).

Output: output/results/<org>/<prefix>_ess_ecoli.csv  (keyed by uniprot_accession).
Run with the `gradi` conda env interpreter. No network.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

TRADIS_DIR = E.REPO_ROOT / "data" / "raw" / "other" / "essentiality" / "enterobacteriaceae_tradis"
GIANT_TAB = TRADIS_DIR / "giant-tab_final.tsv"

COL_BNUM = "Locus: Escherichia coli BW25113 (Keio)"
COL_ECOGENE = "EcoGene Essentiality: Escherichia coli BW25113"
COL_ENTERO = "Enterobacteriaceae %essential"
COL_BACTERIA = "Bacteria %essential"
BNUM_RE = re.compile(r"^b\d{4}$")


def load_bnumber_essentiality() -> pd.DataFrame:
    """Per E. coli b-number: ecogene_essential (bool), entero_pct, bacteria_pct (both [0,1])."""
    if not GIANT_TAB.exists():
        raise SystemExit(
            f"Missing {GIANT_TAB}. Run scripts/07a_fetch_essentiality.py first "
            "(it downloads the Enterobacteriaceae-TraDIS compendium)."
        )
    df = pd.read_csv(GIANT_TAB, sep="\t", low_memory=False,
                     usecols=[COL_BNUM, COL_ECOGENE, COL_ENTERO, COL_BACTERIA])
    df = df.rename(columns={COL_BNUM: "bnumber"})
    df["bnumber"] = df["bnumber"].astype(str).str.strip()
    df = df[df["bnumber"].str.match(BNUM_RE)]
    df["ecogene_essential"] = pd.to_numeric(df[COL_ECOGENE], errors="coerce").fillna(0).astype(int).eq(1)
    df["entero_pct"] = (pd.to_numeric(df[COL_ENTERO], errors="coerce") / 100.0).clip(0, 1)
    df["bacteria_pct"] = (pd.to_numeric(df[COL_BACTERIA], errors="coerce") / 100.0).clip(0, 1)
    df = df.groupby("bnumber", as_index=False).agg(
        ecogene_essential=("ecogene_essential", "max"),
        entero_pct=("entero_pct", "max"),
        bacteria_pct=("bacteria_pct", "max"),
    )
    return df.set_index("bnumber")


def ecoli_acc_to_bnumber() -> dict[str, str]:
    """E. coli K-12 UniProt accession -> its Blattner b-number (from the proteome TSV Gene Names)."""
    tsv = E.proteome_tsv("ecoli")
    d = pd.read_csv(tsv, sep="\t", usecols=["Entry", "Gene Names"])
    out: dict[str, str] = {}
    for acc, gn in zip(d["Entry"], d["Gene Names"].fillna("")):
        for tok in str(gn).split():
            if BNUM_RE.match(tok):
                out[str(acc)] = tok
                break
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]

    bn_ess = load_bnumber_essentiality()
    print(f"compendium: {len(bn_ess)} E. coli b-numbers; "
          f"{int(bn_ess['ecogene_essential'].sum())} EcoGene-essential", flush=True)

    accs = E.load_accessions(org)
    rows = []

    if org == "ecoli":
        # direct: anchor accession -> its own b-number -> call
        acc2b = ecoli_acc_to_bnumber()
        for acc in accs:
            b = acc2b.get(acc)
            rec = bn_ess.loc[b] if (b is not None and b in bn_ess.index) else None
            rows.append({
                "uniprot_accession": acc,
                "ec_transfer_essential": bool(rec["ecogene_essential"]) if rec is not None else False,
                "ec_transfer_via": b or "",
                "ec_transfer_pident": 100.0 if rec is not None else np.nan,
                "n_ecoli_orthologs": 1 if b else 0,
                "entero_pct_essential": float(rec["entero_pct"]) if rec is not None else np.nan,
                "bacteria_pct_essential": float(rec["bacteria_pct"]) if rec is not None else np.nan,
            })
    else:
        # transfer: Kp anchor -> E. coli K-12 ortholog(s) -> b-number -> call; keep the best ortholog
        acc2b = ecoli_acc_to_bnumber()
        orth = E.load_orthologs(org)
        ec = orth[orth["species"] == "Ecoli_K12_MG1655"].copy()
        ec["bnumber"] = ec["target_uniprot"].map(acc2b)
        ec = ec[ec["bnumber"].notna()]
        ec = ec.join(bn_ess, on="bnumber")
        ec = ec.dropna(subset=["entero_pct"])
        # rank orthologs per anchor: essential first, then highest broad-spectrum %essential
        ec["_rank"] = ec["ecogene_essential"].astype(int) * 2 + ec["entero_pct"]
        ec = ec.sort_values("_rank", ascending=False)
        best = ec.drop_duplicates("anchor_uniprot").set_index("anchor_uniprot")
        n_orth = ec.groupby("anchor_uniprot").size()
        for acc in accs:
            b = best.loc[acc] if acc in best.index else None
            rows.append({
                "uniprot_accession": acc,
                "ec_transfer_essential": bool(b["ecogene_essential"]) if b is not None else False,
                "ec_transfer_via": (b["bnumber"] if b is not None else ""),
                "ec_transfer_pident": float(b["pident"]) if (b is not None and not pd.isna(b["pident"])) else np.nan,
                "n_ecoli_orthologs": int(n_orth.get(acc, 0)),
                "entero_pct_essential": float(b["entero_pct"]) if b is not None else np.nan,
                "bacteria_pct_essential": float(b["bacteria_pct"]) if b is not None else np.nan,
            })

    df = pd.DataFrame(rows)

    # For Kp: also transfer the E. coli experimental SCREEN evidence (07n: Keio/Goodall/CRISPRi/
    # vulnerability) onto each Kp anchor via the E. coli K-12 orthology — richer than the EcoGene call.
    if org != "ecoli":
        ec_exp = E.results_dir("ecoli") / "ec_ess_experimental.csv"
        if ec_exp.exists():
            e = pd.read_csv(ec_exp).set_index("uniprot_accession")
            ess_map = e["ecoli_experimental_essential"].map({True: 1.0, False: 0.0}).dropna().to_dict()
            vuln_map = pd.to_numeric(e["ecoli_vulnerability_score"], errors="coerce").dropna().to_dict()
            ess_t = E.transfer_ecoli_to_kp(ess_map, reduce="max")
            vuln_t = E.transfer_ecoli_to_kp(vuln_map, reduce="max")
            df["ec_screens_essential_transfer"] = df["uniprot_accession"].map(ess_t).fillna(0.0) >= 1.0
            df["ec_screens_vulnerability_transfer"] = df["uniprot_accession"].map(vuln_t)
            print(f"  E. coli-screen transfer onto Kp: essential {int(df.ec_screens_essential_transfer.sum())} | "
                  f"with a vulnerability value {int(df.ec_screens_vulnerability_transfer.notna().sum())}", flush=True)
        else:
            df["ec_screens_essential_transfer"] = False
            df["ec_screens_vulnerability_transfer"] = np.nan

    out = E.results_dir(org) / f"{prefix}_ess_ecoli.csv"
    df.to_csv(out, index=False)

    n = len(df)
    print(f"[{org}] {n} proteins | E. coli-essential (transfer): {int(df.ec_transfer_essential.sum())} "
          f"({100*df.ec_transfer_essential.mean():.1f}%) | with an E. coli ortholog: "
          f"{int((df.n_ecoli_orthologs>0).sum())}", flush=True)
    print(f"  broad-spectrum (Enterobacteriaceae %essential >= 0.8): "
          f"{int((df.entero_pct_essential>=0.8).sum())}", flush=True)
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
