"""Merge essentiality tracks + graded composite score/tier (docs §4, result).

Joins every per-track table on uniprot_accession and computes a transparent, graded [0,1]
`essentiality_score` (the spec asks for "a graded 0–1 vulnerability, not a binary call") as a
weighted ensemble of three independently-kept sub-scores:

  * evidence_experimental — direct Kp Tn-seq/CRISPRi (07b): ECL8 in-vitro essential is strongest;
    urine/serum niche fitness, in-vivo defect and CRISPRi vulnerability are supporting.
  * evidence_transfer     — E. coli essential ortholog + graded broad-spectrum conservation (07c:
    EcoGene essential -> 1.0, else Enterobacteriaceae %essential).
  * evidence_predictor    — consensus of the computational predictors (ProteomeLM 07d, Geptop 07e,
    FBA 07f), a weighted mean over whichever predictors are available for the protein.

Crucially, a sub-score that is entirely ABSENT for a protein (no Kp screen coverage, no E. coli
ortholog) is DROPPED and the remaining weights are renormalised — a missing track lowers confidence,
it does not push the score toward zero. All sub-scores/weights are columns/constants (auditable).

Tier (evidence-driven, not a pure threshold): `essential` if there is a direct Kp essential call, or
a strong predictor+transfer consensus, or score >= 0.60; `likely_essential` if score >= 0.35 or any
partial signal; else `non_essential`.

Outputs (output/results/<org>/):
  <prefix>_essentiality.csv           — full per-protein table
  <prefix>_essentiality_shortlist.csv — broad-spectrum-selective (03c) AND essential, ranked
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

W_EXPERIMENTAL = E.W_EXPERIMENTAL   # 0.40
W_TRANSFER = E.W_ECOLI_TRANSFER     # 0.20
W_PREDICTOR = E.W_PREDICTOR         # 0.40
PRED_W = E.PRED_WEIGHTS             # {proteomelm:.5, geptop:.3, fba:.2}
TIER_ESSENTIAL = E.TIER_ESSENTIAL   # 0.60
TIER_LIKELY = E.TIER_LIKELY         # 0.35
CATEGORIES = E.REPO_ROOT / "data" / "processed" / "other" / "orthology" / "three_way_protein_categories.tsv"


def _read(org: str, suffix: str) -> pd.DataFrame | None:
    _, prefix = E.ORGANISMS[org]
    p = E.results_dir(org) / f"{prefix}_ess_{suffix}.csv"
    if not p.exists():
        print(f"  [warn] missing track table {p.name}; its columns will default", flush=True)
        return None
    return pd.read_csv(p).drop_duplicates("uniprot_accession")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]

    base = pd.DataFrame({"uniprot_accession": E.load_accessions(org)})
    base["gene"] = base["uniprot_accession"].map(E.load_genes(org))
    for suffix in ("kp", "ecoli", "proteomelm", "geptop", "fba", "deeplyessential",
                   "publications", "experimental"):
        t = _read(org, suffix)
        if t is not None:
            # keep only the publication columns the composite consumes (avoid clobbering shared names)
            if suffix == "publications":
                keep = [c for c in ("uniprot_accession", "crispri_ce_library", "kpnih1_essential",
                                    "crispri_invivo_hit") if c in t.columns]
                t = t[keep]
            base = base.merge(t, on="uniprot_accession", how="left")
    df = base

    # ---- evidence_experimental (direct Kp screens; NaN where the anchor has no Kp measurement)
    # Combines the genome-wide ECL8 Tn-seq call (essential/unclear/non), the Kp Mobile-CRISPRi-seq
    # conditionally-essential library + KPNIH1 Tn-seq essential set (both positive-only essential
    # lists, 07l), and the niche / in-vivo / vulnerability supporting signals. `has_any` is true when
    # the protein carries ANY of these measurements; the CRISPRi/KPNIH1 lists only ADD essential calls
    # (absence is not evidence of non-essentiality). Conditionally-essential-in-CRISPRi beats
    # non-essential-in-LB via max().
    def _true(v):
        return v is True or str(v) == "True"

    def experimental_row(r):
        iv = str(r.get("kp_ess_in_vitro_call"))
        crispri_lib = _true(r.get("crispri_ce_library"))
        kpnih1 = _true(r.get("kpnih1_essential"))
        crispri_vivo = _true(r.get("crispri_invivo_hit"))
        ecl8_measured = pd.notna(r.get("kp_ess_in_vitro_call")) and iv not in ("nan", "")
        niche = (str(r.get("kp_ess_urine_call")) == "conditional"
                 or str(r.get("kp_ess_serum_call")) == "conditional"
                 or str(r.get("kp_ess_in_vivo_call")) == "in_vivo_defect")
        vuln = pd.notna(r.get("kp_ess_vulnerability_call")) and str(r.get("kp_ess_vulnerability_call")) not in ("nan", "")
        # no direct Kp measurement at all -> drop (renormalised in the composite)
        if not (ecl8_measured or crispri_lib or kpnih1 or crispri_vivo or niche or vuln):
            return np.nan
        # rigorous genome-wide essential call (bimodal Tn-seq): ECL8-essential or KPNIH1-essential -> 1.0
        rigorous = 1.0 if (E.essentiality_ladder(iv if iv != "nan" else None) >= 1.0 or kpnih1) else 0.0
        # CRISPRi conditionally-essential library membership: strong but softer than a Tn-seq call
        crispri = 0.85 if crispri_lib else 0.0
        support = 0.0
        if niche:
            support = max(support, 0.6)
        if crispri_vivo or vuln:
            support = max(support, 0.7)
        unclear = E.essentiality_ladder(iv if iv != "nan" else None)  # 0.4 for ECL8 'unclear', else 0
        return float(max(rigorous, crispri, support, unclear))

    def experimental_row_ecoli(r):
        """E. coli's own screens (07n) — the direct 0.40 experimental axis, symmetric with Kp."""
        rigorous = _true(r.get("ecoli_keio_essential")) or _true(r.get("ecoli_goodall_essential"))
        n_scr = sum(_true(r.get(c)) for c in ("ecoli_keio_essential", "ecoli_goodall_essential",
                                              "ecoli_crispri_rousset18_essential", "ecoli_crispri_wang18_essential"))
        v = pd.to_numeric(r.get("ecoli_vulnerability_score"), errors="coerce")
        measured = rigorous or n_scr > 0 or pd.notna(v) or pd.notna(pd.to_numeric(r.get("ecoli_crispri_rousset18_log2fc"), errors="coerce"))
        if not measured:
            return np.nan
        if rigorous:
            return 1.0
        # consensus CRISPRi essential (>=2 screens) is strong; else the graded vulnerability
        cons = 0.9 if n_scr >= 2 else (0.75 if n_scr == 1 else 0.0)
        return float(max(cons, v if pd.notna(v) else 0.0))

    df["evidence_experimental"] = df.apply(
        experimental_row_ecoli if org == "ecoli" else experimental_row, axis=1)

    # ---- evidence_transfer
    # For Kp: best of the EcoGene ortholog call, the graded broad-spectrum %essential, and the richer
    # E. coli-SCREEN transfer (07c: Keio/Goodall/CRISPRi/vulnerability lifted via the E. coli ortholog).
    # For E. coli: broad-spectrum conservation only (its own essentiality is the experimental axis).
    ec_ess = df.get("ec_transfer_essential", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    entero = pd.to_numeric(df.get("entero_pct_essential"), errors="coerce")
    transfer = np.where(ec_ess, 1.0, entero).astype(float)
    if org != "ecoli":
        scr_ess = df.get("ec_screens_essential_transfer", pd.Series(False, index=df.index)).fillna(False).astype(bool)
        scr_vuln = pd.to_numeric(df.get("ec_screens_vulnerability_transfer"), errors="coerce")
        transfer = np.fmax(transfer, np.where(scr_ess, 1.0, scr_vuln))
    df["evidence_transfer"] = transfer  # NaN stays NaN (no ortholog / no data)

    # ---- evidence_predictor (weighted mean over available predictors)
    plm = pd.to_numeric(df.get("proteomelm_ess_score"), errors="coerce")
    gep = pd.to_numeric(df.get("geptop_score"), errors="coerce")
    fba = df.get("fba_essential")
    fba = fba.map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0}) if fba is not None else pd.Series(np.nan, index=df.index)
    fba = pd.to_numeric(fba, errors="coerce")
    preds = {"proteomelm": plm, "geptop": gep, "fba": fba}
    num = pd.Series(0.0, index=df.index)
    den = pd.Series(0.0, index=df.index)
    for name, s in preds.items():
        w = PRED_W[name]
        avail = s.notna()
        num = num.add((s.fillna(0.0) * w).where(avail, 0.0), fill_value=0.0)
        den = den.add(pd.Series(np.where(avail, w, 0.0), index=df.index), fill_value=0.0)
    df["evidence_predictor"] = (num / den).where(den > 0, np.nan)

    # ---- composite: renormalised weighted sum over AVAILABLE sub-scores
    subs = {"evidence_experimental": W_EXPERIMENTAL,
            "evidence_transfer": W_TRANSFER,
            "evidence_predictor": W_PREDICTOR}
    cnum = pd.Series(0.0, index=df.index)
    cden = pd.Series(0.0, index=df.index)
    for col, w in subs.items():
        s = pd.to_numeric(df[col], errors="coerce")
        avail = s.notna()
        cnum = cnum.add((s.fillna(0.0) * w).where(avail, 0.0), fill_value=0.0)
        cden = cden.add(pd.Series(np.where(avail, w, 0.0), index=df.index), fill_value=0.0)
    df["essentiality_score"] = (cnum / cden).where(cden > 0, 0.0).round(4)

    # ---- hard flags + tier
    kp_essential_call = df.get("kp_ess_in_vitro_call", pd.Series(index=df.index)).astype(str).eq("essential")
    kpnih1_ess = df.get("kpnih1_essential", pd.Series(False, index=df.index)).map(_true)
    # E. coli's own rigorous KO/TraDIS essential call (07n), for the E. coli organism run
    ec_own_ess = (df.get("ecoli_keio_essential", pd.Series(False, index=df.index)).map(_true)
                  | df.get("ecoli_goodall_essential", pd.Series(False, index=df.index)).map(_true))
    # experimentally essential = a rigorous genome-wide essential call: Kp ECL8/KPNIH1 Tn-seq, E. coli
    # Keio/Goodall, or the E. coli EcoGene ortholog. (CRISPRi CE-library membership strongly boosts the
    # SCORE via evidence_experimental, but is not itself a hard essential-tier override.)
    df["experimentally_essential"] = kp_essential_call | kpnih1_ess | ec_ess | ec_own_ess
    strong_consensus = (df["evidence_predictor"].fillna(0) >= 0.7) & (df["evidence_transfer"].fillna(0) >= 0.7)

    def tier_row(r):
        if r["experimentally_essential"] or strong_consensus.loc[r.name] or r["essentiality_score"] >= TIER_ESSENTIAL:
            return "essential"
        partial = (r["essentiality_score"] >= TIER_LIKELY
                   or (pd.notna(r["evidence_experimental"]) and r["evidence_experimental"] >= 0.5)
                   or (pd.notna(r["evidence_predictor"]) and r["evidence_predictor"] >= 0.5))
        return "likely_essential" if partial else "non_essential"

    df["essentiality_tier"] = df.apply(tier_row, axis=1)

    # ---- provenance string
    def evidence_sources(r):
        s = []
        if kp_essential_call.loc[r.name]: s.append("Kp_ECL8")
        if _true(r.get("crispri_ce_library")): s.append("CRISPRi")
        if _true(r.get("kpnih1_essential")): s.append("KPNIH1")
        if _true(r.get("crispri_invivo_hit")): s.append("CRISPRi_invivo")
        if str(r.get("kp_ess_sources")) not in ("nan", ""): s.append(str(r.get("kp_ess_sources")))
        # E. coli's own screens (07n) — provenance for the E. coli organism run
        if _true(r.get("ecoli_keio_essential")): s.append("Keio")
        if _true(r.get("ecoli_goodall_essential")): s.append("Goodall")
        if _true(r.get("ecoli_crispri_rousset18_essential")): s.append("Rousset18")
        if _true(r.get("ecoli_crispri_wang18_essential")): s.append("Wang18")
        if ec_ess.loc[r.name]: s.append("Ecoli")
        if _true(r.get("ec_screens_essential_transfer")): s.append("Ecoli_screens")
        if pd.notna(r.get("proteomelm_ess_score")) and r.get("proteomelm_ess_score", 0) >= 0.5: s.append("ProteomeLM")
        if pd.notna(r.get("geptop_score")) and r.get("geptop_score", 0) >= E.GEPTOP_CUTOFF: s.append("Geptop")
        if r.get("fba_essential") in (True, "True"): s.append("FBA")
        return ";".join([x for x in s if x])

    df["evidence_sources"] = df.apply(evidence_sources, axis=1)

    # ---- selectivity (03c) for the shortlist
    if CATEGORIES.exists():
        cat = pd.read_csv(CATEGORIES, sep="\t")
        cat = cat[cat["organism"] == org][["uniprot_accession", "selectivity"]]
        df = df.merge(cat, on="uniprot_accession", how="left")
    else:
        df["selectivity"] = pd.NA

    out = E.results_dir(org) / f"{prefix}_essentiality.csv"
    df.to_csv(out, index=False)

    # ---- shortlist: broad-spectrum-selective AND essential, ranked by score then predictor
    short = df[(df["selectivity"] == "broad_selective") & (df["essentiality_tier"] == "essential")].copy()
    short = short.sort_values(["experimentally_essential", "essentiality_score",
                               "evidence_predictor"], ascending=False)
    keep = ["uniprot_accession", "gene", "essentiality_score", "essentiality_tier",
            "experimentally_essential", "evidence_experimental", "evidence_transfer",
            "evidence_predictor", "entero_pct_essential", "selectivity", "evidence_sources"]
    keep = [c for c in keep if c in short.columns]
    short_out = E.results_dir(org) / f"{prefix}_essentiality_shortlist.csv"
    short[keep].to_csv(short_out, index=False)

    # ---- summary
    tiers = df["essentiality_tier"].value_counts().to_dict()
    print(f"[{org}] {len(df)} proteins | tiers: {tiers}", flush=True)
    print(f"  experimentally essential: {int(df.experimentally_essential.sum())} | "
          f"score>=0.6: {int((df.essentiality_score>=0.6).sum())}", flush=True)
    print(f"  shortlist (broad-selective & essential): {len(short)} -> {short_out.name}", flush=True)
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
