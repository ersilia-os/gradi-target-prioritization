"""Direct E. coli K-12 experimental essentiality (docs §4.1/§4.2) — the E. coli analog of 07b.

E. coli is the most functionally-screened bacterium; this makes it a first-class experimental target
(not just a transfer reference). Ingests the major published screens and lands each per-gene call on
the E. coli proteome (b-number / Keio-JW / gene-symbol keyed — near-complete mapping):

  * **Keio / PEC** (Baba 2006) — single-gene KO essential set (`Class==1` in PECData.dat).
  * **Goodall 2018** (mBio) — BW25113 TraDIS essential call (Table S4 `Essential`).
  * **Rousset 2018** (PLoS Genet) — genome-wide CRISPRi gene-level median log2FC depletion (S2 table).
  * **Wang 2018** (Nat Commun) — pooled CRISPRi gene fitness (essential-genes sheet; essential < -6).
  * **Rousset 2021** (Nat Microbiol) — CRISPRi across 18 E. coli strains → fraction-essential.
  * **Hawkins 2020** (Cell Sys) — mismatch-CRISPRi relative fitness → quantitative vulnerability.
  * **RB-TnSeq / Fitness Browser** (Price 2018) — per-gene fitness across ~3,500 conditions
    (min-fitness vulnerability + the 280-condition stress/antibiotic matrix for the 07o condition view).

Output: output/results/ecoli/ec_ess_experimental.csv (keyed by uniprot_accession) + a condition matrix
cached under data/processed/ecoli/essentiality/conditions/. Run with the `gradi` env. No network.
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

B = E.essentiality_raw_dir("ecoli")
ROUSSET18_CUT = -2.0     # CRISPRi gene-level median log2FC <= this -> essential/depleted
WANG18_CUT = -6.0        # Wang gene fitness < this -> essential (their definition)
MULTISTRAIN_CUT = -3.0   # Rousset 2021 gene score <= this -> essential in that strain/condition
BNUM_RE = re.compile(r"\bb\d{4}\b")


def _norm(g):
    return str(g).strip().lower()


def keio_pec(sym, bnum):
    """Keio/PEC single-KO essential set (Class==1)."""
    f = B / "pec" / "PECData.dat"
    if not f.exists():
        return {}
    d = pd.read_csv(f, sep="\t", dtype=str)
    cls_col = [c for c in d.columns if c.startswith("Class")][0]
    orf_col = "Orf"
    alt_col = "Alternative name"
    out = {}
    for _, r in d.iterrows():
        if str(r[cls_col]).strip() != "1":
            continue
        acc = sym.get(_norm(r[orf_col]))
        if not acc:
            m = BNUM_RE.search(str(r.get(alt_col, "")))
            if m:
                acc = bnum.get(m.group(0).upper())
        if acc:
            out[acc] = True
    return out


def goodall(sym):
    f = B / "goodall2018_Ec_BW25113" / "mbo001183726st4.xlsx"
    if not f.exists():
        return {}
    xl = pd.ExcelFile(f)
    d = xl.parse(xl.sheet_names[0], header=1)
    out = {}
    for _, r in d.iterrows():
        acc = sym.get(_norm(r.get("Gene")))
        if acc and str(r.get("Essential")).strip().lower() in ("true", "1", "yes"):
            out[acc] = True
    return out


def rousset18(sym):
    f = B / "rousset2018_crispri" / "pgen.1007749.s012.csv"
    if not f.exists():
        return {}, {}
    d = pd.read_csv(f)
    lfc, ess = {}, {}
    for _, r in d.iterrows():
        acc = sym.get(_norm(r.get("gene")))
        if not acc:
            continue
        v = pd.to_numeric(r.get("median_coding"), errors="coerce")
        if pd.notna(v):
            if acc not in lfc or v < lfc[acc]:
                lfc[acc] = float(v)
            ess[acc] = ess.get(acc, False) or (v <= ROUSSET18_CUT)
    return lfc, ess


def wang18(sym):
    f = B / "wang2018_crispri" / "41467_2018_4899_MOESM8_ESM.xlsx"
    if not f.exists():
        return {}, {}
    d = pd.read_excel(f, sheet_name="essential genes")
    fit, ess = {}, {}
    fcol = [c for c in d.columns if "fitness" in c.lower()][0]
    for _, r in d.iterrows():
        acc = sym.get(_norm(r.get("gene")))
        v = pd.to_numeric(r.get(fcol), errors="coerce")
        if acc and pd.notna(v):
            if acc not in fit or v < fit[acc]:
                fit[acc] = float(v)
            ess[acc] = ess.get(acc, False) or (v < WANG18_CUT)
    return fit, ess


def rousset21(sym):
    """Fraction of the 18-strain × condition panel where the gene scores as essential."""
    f = B / "rousset2021_crispri" / "MOESM4.xlsx"
    if not f.exists():
        return {}
    d = pd.read_excel(f, sheet_name="Table S4", header=1)
    gcol = d.columns[0]
    score_cols = [c for c in d.columns[1:] if pd.api.types.is_numeric_dtype(d[c])]
    if not score_cols:
        # coerce
        for c in d.columns[1:]:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        score_cols = list(d.columns[1:])
    frac = {}
    for _, r in d.iterrows():
        acc = sym.get(_norm(r.get(gcol)))
        if not acc:
            continue
        vals = pd.to_numeric(r[score_cols], errors="coerce")
        n = vals.notna().sum()
        if n:
            frac[acc] = float((vals <= MULTISTRAIN_CUT).sum() / n)
    return frac


def hawkins(sym):
    """Per-gene vulnerability from mismatch-CRISPRi relative fitness (min = most severe knockdown)."""
    f = B / "hawkins2020_crispri" / "mmc4.xlsx"
    if not f.exists():
        return {}
    xl = pd.ExcelFile(f)
    sh = [s for s in xl.sheet_names if "eco" in s.lower() and "fitness" in s.lower()]
    if not sh:
        return {}
    d = xl.parse(sh[0])
    fcol = [c for c in d.columns if "relative fitness (mean)" in c.lower()]
    if not fcol or "gene" not in d.columns:
        return {}
    fcol = fcol[0]
    vmin = {}
    for _, r in d.iterrows():
        acc = sym.get(_norm(r.get("gene")))
        v = pd.to_numeric(r.get(fcol), errors="coerce")
        if acc and pd.notna(v):
            vmin[acc] = min(vmin.get(acc, np.inf), float(v))
    return vmin


def rbtnseq(bnum):
    """RB-TnSeq per-gene min-fitness (strongest defect across ~3,500 conditions)."""
    f = B / "fitness_browser" / "ec_rbtnseq_gene_summary.csv"
    if not f.exists():
        return {}
    d = pd.read_csv(f)
    out = {}
    for _, r in d.iterrows():
        acc = bnum.get(str(r.get("sysName")).strip().upper())
        v = pd.to_numeric(r.get("rbtnseq_min_fitness"), errors="coerce")
        if acc and pd.notna(v):
            out[acc] = float(v)
    return out


def build_condition_matrix(bnum):
    """Map the RB-TnSeq stress/antibiotic condition matrix onto E. coli accessions (for 07o)."""
    f = B / "fitness_browser" / "ec_rbtnseq_stress_conditions.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    d["uniprot_accession"] = d["sysName"].astype(str).str.upper().map(bnum)
    d = d.dropna(subset=["uniprot_accession"])
    out = E.essentiality_processed_dir("ecoli", "conditions") / "ec_rbtnseq_conditions.csv"
    cond_cols = [c for c in d.columns if c not in ("locusId", "sysName", "desc", "uniprot_accession")]
    d[["uniprot_accession", "desc"] + cond_cols].to_csv(out, index=False)
    print(f"  condition matrix: {len(d)} genes x {len(cond_cols)} conditions -> {out.relative_to(E.REPO_ROOT)}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=["ecoli"], default="ecoli")
    args = ap.parse_args()
    org = "ecoli"
    accs = E.load_accessions(org)
    sym = E.gene_aliases_to_uniprot(org)
    bnum = E.locus_to_uniprot(org)  # {B####: acc}

    keio = keio_pec(sym, bnum)
    good = goodall(sym)
    r18_lfc, r18_ess = rousset18(sym)
    w18_fit, w18_ess = wang18(sym)
    r21 = rousset21(sym)
    haw = hawkins(sym)
    rbt = rbtnseq(bnum)
    build_condition_matrix(bnum)

    df = pd.DataFrame({"uniprot_accession": accs})
    df["ecoli_keio_essential"] = df["uniprot_accession"].isin(keio)
    df["ecoli_goodall_essential"] = df["uniprot_accession"].isin(good)
    df["ecoli_crispri_rousset18_log2fc"] = df["uniprot_accession"].map(r18_lfc)
    df["ecoli_crispri_rousset18_essential"] = df["uniprot_accession"].map(lambda a: r18_ess.get(a, False))
    df["ecoli_crispri_wang18_fitness"] = df["uniprot_accession"].map(w18_fit)
    df["ecoli_crispri_wang18_essential"] = df["uniprot_accession"].map(lambda a: w18_ess.get(a, False))
    df["ecoli_crispri_multistrain_frac"] = df["uniprot_accession"].map(r21)
    df["ecoli_vulnerability_hawkins"] = df["uniprot_accession"].map(haw)
    df["ecoli_rbtnseq_min_fitness"] = df["uniprot_accession"].map(rbt)

    # consensus essential across the 4 binary/thresholded screens
    ess_cols = ["ecoli_keio_essential", "ecoli_goodall_essential",
                "ecoli_crispri_rousset18_essential", "ecoli_crispri_wang18_essential"]
    df["n_ecoli_screens_essential"] = df[ess_cols].sum(axis=1)
    # first-class experimental essential: the gold KO/TraDIS binary, or >=2 screens agree
    df["ecoli_experimental_essential"] = (df["ecoli_keio_essential"] | df["ecoli_goodall_essential"]
                                          | (df["n_ecoli_screens_essential"] >= 2))

    # graded GENERAL (standard-growth) vulnerability score [0,1]: strongest depletion across the
    # standard-condition CRISPRi screens (Rousset18 LB, Wang18) + Hawkins mismatch-CRISPRi. RB-TnSeq
    # min-fitness is deliberately EXCLUDED here — it spans ~3,500 conditions and captures
    # condition-specific essentiality (e.g. lacZ on lactose), which belongs in the 07o condition view,
    # not the general vulnerability. `ecoli_rbtnseq_min_fitness` is kept as its own column for 07o.
    def scale(v, lo, hi):  # lo (most vulnerable) -> 1, hi (neutral) -> 0
        return np.clip((hi - v) / (hi - lo), 0, 1)
    vuln = pd.concat([
        scale(pd.to_numeric(df["ecoli_crispri_rousset18_log2fc"], errors="coerce"), -6, 0),
        scale(pd.to_numeric(df["ecoli_crispri_wang18_fitness"], errors="coerce"), -10, 0),
        scale(pd.to_numeric(df["ecoli_vulnerability_hawkins"], errors="coerce"), 0, 1),
    ], axis=1).max(axis=1)
    # essentials without a graded score still count as maximally vulnerable
    df["ecoli_vulnerability_score"] = np.where(df["ecoli_experimental_essential"] & vuln.isna(), 1.0, vuln).round(4)

    out = E.results_dir(org) / "ec_ess_experimental.csv"
    df.to_csv(out, index=False)

    print(f"[ecoli] screens mapped — Keio/PEC:{len(keio)} Goodall:{len(good)} "
          f"Rousset18:{len(r18_lfc)} Wang18:{len(w18_fit)} Rousset21:{len(r21)} "
          f"Hawkins:{len(haw)} RB-TnSeq:{len(rbt)}", flush=True)
    print(f"[ecoli] experimentally essential: {int(df.ecoli_experimental_essential.sum())} "
          f"(Keio {int(df.ecoli_keio_essential.sum())} / Goodall {int(df.ecoli_goodall_essential.sum())} / "
          f">=2 screens {int((df.n_ecoli_screens_essential>=2).sum())})", flush=True)
    print(f"[ecoli] wrote {out.relative_to(E.REPO_ROOT)} ({len(df)} proteins)", flush=True)


if __name__ == "__main__":
    main()
