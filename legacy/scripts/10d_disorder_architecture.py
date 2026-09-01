"""Disorder architecture: terminal, internal and global (docs §3.1; the Regime B feature layer).

Why this table exists. The 2026-08-06 spec revision splits the composite into two bars, and the
**validation bar** (`partnerless_clpP` — the one matching the assay that actually gates progression)
is built from global disorder, `two_domain_architecture`, stability and assembly state. Activated
ClpP has no sequence degron at all, so once the unfoldase is removed, disorder *is* the recognition
story. This script is the authoritative disorder table for both organisms.

It supersedes `10b`'s structural columns, which used 10/20-residue windows; the spec asks for 15/30/50
and for the measured ClpXP thresholds 5/20/37. `10b` keeps the sequence motifs. Downstream, `10i`
should read structure from here.

What it adds beyond 10b:

  * **The spec's windows** — mean pLDDT over the terminal 15/30/50 residues, per terminus. 30 matters
    specifically for parity with Won 2024, whose top feature is N-terminal-30 disorder.
  * **Initiation regions at two strictness levels** (pLDDT < 70 and < 50), bucketed against the
    measured thresholds: ~5 aa closed-channel recognition, ~20 aa open-channel engagement, ~37 aa
    pore-to-active-site reach.
  * **Internal disordered runs** — degrons need not be terminal. Internally placed ssrA tags are
    degraded by ClpAP/ClpXP (Hoskins 2002), and E. coli FtsZ carries an internal ClpXP site in its
    disordered linker at residues 349-358 (Camberg 2014).
  * **MobiDB** — 8 independent sequence-side disorder predictors, a consensus, LIP regions and
    Pfam/Gene3D domain boundaries, at 100% of both proteomes. This is the guard against the known
    pLDDT failure mode: an obligate complex subunit is ordered *in situ* but predicted disordered in
    isolation, so a structure-only disorder call reads assembly context as intrinsic disorder.

**And it scores every one of those features against measured turnover.** `10b`'s validation table
covers only the sequence motifs, so the axis's central claim — that disorder predicts degradability —
has never actually been tested here. The `<prefix>_disorder_validation.csv` written alongside is that
test (spec §8 check #5), using Nagar 2021's classes: direct for E. coli, ortholog-transferred for Kp.

Outputs:
  output/results/<org>/<prefix>_disorder.csv
  data/processed/<org>/degradability/<prefix>_disorder_validation.csv
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

# An "engageable" terminus is one clearing the open-channel threshold (~20 aa): enough unstructured
# sequence for five pore-1 loops to grip. Below that, recognition may still occur but processive
# unfolding is not licensed by the measurements.
ENGAGE_AA = D.INIT_THRESHOLDS[1]


def _architecture(row: pd.Series) -> str:
    """Coarse disorder architecture — the shape of the disorder, not just its amount.

    Ordering matters: `mostly_disordered` wins outright (there is no folded core to speak of), then a
    terminal IDR (the canonical engagement route), then an internal linker (the FtsZ route).
    """
    g = row.get("global_disorder_frac_plddt")
    if g is None or pd.isna(g):
        return "unknown"
    if g >= 0.50:
        return "mostly_disordered"
    term = max(row.get("nterm_init_len_70") or 0, row.get("cterm_init_len_70") or 0)
    internal = row.get("internal_disorder_max_run") or 0
    if term >= ENGAGE_AA and internal >= D.MIN_INTERNAL_RUN:
        return "multi_idr"
    if term >= ENGAGE_AA:
        return "terminal_idr"
    if internal >= D.MIN_INTERNAL_RUN:
        return "internal_linker"
    return "ordered"


def _row(acc: str, organism: str) -> dict:
    plddt = D.plddt_series(organism, acc)
    r: dict = {"uniprot_accession": acc, "af_available": bool(plddt),
               "modeled_len": len(plddt) or None}
    if not plddt:
        return r

    for w in D.TERM_WINDOWS_SPEC:
        r[f"nterm_plddt_{w}"] = D.terminal_plddt(plddt, w, "n")
        r[f"cterm_plddt_{w}"] = D.terminal_plddt(plddt, w, "c")
    for cut, tag in ((D.PLDDT_DISORDER_CUT, 70), (D.PLDDT_DISORDER_CUT_STRICT, 50)):
        for t in ("n", "c"):
            r[f"{t}term_init_len_{tag}"] = D.init_region_len(plddt, t, cut)
    for t in ("n", "c"):
        r[f"{t}term_init_bucket"] = D.init_bucket(r[f"{t}term_init_len_70"])
        r[f"{t}term_exposure"] = D.exposure(r[f"{t}term_plddt_15"])

    r.update(D.internal_disorder_runs(plddt))
    r["internal_degron_candidate"] = bool(
        (r["internal_disorder_max_run"] or 0) >= D.MIN_INTERNAL_RUN)
    r["global_disorder_frac_plddt"] = D.global_disorder_frac(plddt)
    r["global_disorder_frac_plddt50"] = D.global_disorder_frac(plddt, D.PLDDT_DISORDER_CUT_STRICT)
    r["plddt_mean"] = round(sum(plddt) / len(plddt), 2)

    n_eng = (r["nterm_init_len_70"] or 0) >= ENGAGE_AA
    c_eng = (r["cterm_init_len_70"] or 0) >= ENGAGE_AA
    r["engageable_terminus"] = ("both" if n_eng and c_eng else
                                "n" if n_eng else "c" if c_eng else "none")
    r["max_init_len_70"] = max(r["nterm_init_len_70"] or 0, r["cterm_init_len_70"] or 0)
    return r


# --------------------------------------------------------------------------- validation
def _turnover_labels(organism: str) -> pd.DataFrame:
    """Measured turnover class per accession — direct for E. coli, orthology-transferred for Kp."""
    nagar = D.load_nagar()
    if nagar.empty:
        return pd.DataFrame(columns=["uniprot_accession", "halflife_class", "halflife_min"])
    keep = ["uniprot_accession", "halflife_class", "halflife_min"]
    if organism == "ecoli":
        return nagar[keep]
    orth = D.load_orthologs(organism)
    sub = orth[orth["species"] == D.SPECIES_ECOLI][["anchor_uniprot", "target_uniprot"]]
    ec = nagar[keep].rename(columns={"uniprot_accession": "ec_uniprot"})
    m = sub.merge(ec, left_on="target_uniprot", right_on="ec_uniprot", how="inner")
    rank = {"fast": 0, "intermediate": 1, "stable": 2}
    m["r"] = m["halflife_class"].map(rank)
    m = m.sort_values("r").drop_duplicates("anchor_uniprot")
    return m.rename(columns={"anchor_uniprot": "uniprot_accession"})[keep]


# Every feature that claims to measure disorder, so the ranking is a fair fight between the
# structure-derived, the sequence-derived (MobiDB) and Nagar's own descriptors.
def _feature_list(df: pd.DataFrame) -> list[tuple[str, str, int]]:
    """(column, human label, sign) — sign +1 when HIGHER should mean faster turnover."""
    feats: list[tuple[str, str, int]] = [
        ("global_disorder_frac_plddt", "global disorder (pLDDT<70)", +1),
        ("global_disorder_frac_plddt50", "global disorder (pLDDT<50)", +1),
        ("plddt_mean", "mean pLDDT", -1),
        ("internal_disorder_max_run", "internal max run", +1),
        ("internal_disorder_frac", "internal disorder frac", +1),
        ("max_init_len_70", "best terminus init-region", +1),
    ]
    for w in D.TERM_WINDOWS_SPEC:
        feats.append((f"nterm_plddt_{w}", f"N-term pLDDT ({w} aa)", -1))
        feats.append((f"cterm_plddt_{w}", f"C-term pLDDT ({w} aa)", -1))
    for t, lab in (("n", "N"), ("c", "C")):
        feats.append((f"{t}term_init_len_70", f"{lab}-term init-region (<70)", +1))
        feats.append((f"{t}term_init_len_50", f"{lab}-term init-region (<50)", +1))
    for c, lab in (("mobidb_disorder_priority", "MobiDB consensus"),
                   ("mobidb_dis_dis465", "MobiDB dis465"),
                   ("mobidb_dis_iupl", "MobiDB IUPred-long"),
                   ("mobidb_disorder_alphafold", "MobiDB AF-disorder"),
                   ("mobidb_lip_frac", "MobiDB LIP fraction"),
                   ("mobidb_n_domains_merge", "n domains")):
        if c in df.columns:
            feats.append((c, lab, +1))
    return [(c, lab, s) for c, lab, s in feats if c in df.columns]


def _validate(df: pd.DataFrame, organism: str) -> pd.DataFrame:
    lab = _turnover_labels(organism)
    if lab.empty:
        print("  [validate] no Nagar turnover data — run scripts/10a_fetch_degradability.py", flush=True)
        return pd.DataFrame()
    j = df.merge(lab, on="uniprot_accession", how="inner")
    j = j[j["halflife_class"].isin(["fast", "intermediate", "stable"])]
    # fast vs stable — the cleanest contrast; `intermediate` is the noisy middle and is excluded from
    # the AUROC so the null stays interpretable.
    fs = j[j["halflife_class"].isin(["fast", "stable"])].copy()
    fs["is_fast"] = fs["halflife_class"] == "fast"
    print(f"\n  [validate] {len(j)} proteins with a measured class "
          f"({int((j.halflife_class == 'fast').sum())} fast); AUROC computed on "
          f"{len(fs)} fast-vs-stable", flush=True)

    rows = []
    for col, label, sign in _feature_list(df):
        v = fs[col].astype(float) * sign     # orient so higher = expected faster
        res = D.disorder_auroc(v, fs["is_fast"])
        rows.append({"feature": col, "label": label, "oriented_sign": sign, **res})
    tab = pd.DataFrame(rows).sort_values("auroc", ascending=False, na_position="last")
    show = tab.copy()
    show["95% CI"] = [f"[{lo:.2f}, {hi:.2f}]" if pd.notna(lo) else ""
                      for lo, hi in zip(show["auroc_lo"], show["auroc_hi"])]
    print(show[["label", "auroc", "95% CI", "spearman", "n", "n_pos"]].to_string(index=False),
          flush=True)

    best = tab.iloc[0]
    # With ~60-75 positives the CI is roughly +/-0.07 wide, so the ranking order is not meaningful on
    # its own. What IS meaningful is whether the interval excludes the 0.5 null.
    sig = tab[tab["auroc_lo"] > 0.5]
    print(f"\n  [validate] {len(sig)}/{len(tab)} features have a 95% CI excluding the 0.5 null",
          flush=True)
    if pd.isna(best["auroc"]) or best["auroc"] < 0.55 or len(sig) == 0:
        print(f"  [WARN] best disorder feature reaches AUROC {best['auroc']} "
              f"[{best['auroc_lo']}, {best['auroc_hi']}] — at or near the 0.5 null. Disorder does "
              f"NOT usefully predict measured turnover here; say so in the log and the figure "
              f"rather than weighting it.", flush=True)
    else:
        overlap = tab[(tab["auroc_lo"] <= best["auroc"]) & (tab["auroc"] >= 0.5)]
        print(f"  [validate] best = {best['label']} at AUROC {best['auroc']} "
              f"[{best['auroc_lo']}, {best['auroc_hi']}]; {len(overlap)} features are statistically "
              f"indistinguishable from it — do not read the ranking order as a result.", flush=True)

    # The Won 2024 asymmetry, stated explicitly: N-terminal-30 was his top feature and C-terminal-30
    # was null — but that is ClpC1 in an Actinobacterium, and ClpX reads both termini.
    got = {r["feature"]: r["auroc"] for _, r in tab.iterrows()}
    n30, c30 = got.get("nterm_plddt_30"), got.get("cterm_plddt_30")
    if n30 and c30:
        print(f"  [validate] Won asymmetry check — N-30 AUROC {n30:.3f} vs C-30 AUROC {c30:.3f}: "
              f"{'N-dominant (reproduces Won)' if n30 - c30 > 0.03 else 'C-dominant (opposite to Won)' if c30 - n30 > 0.03 else 'comparable — asymmetry does NOT transfer'}",
              flush=True)
    return tab


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", choices=list(D.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = D.ORGANISMS[org]

    accs = D.load_accessions(org)
    print(f"[{org}] {len(accs)} accessions; reading cached AlphaFold pLDDT profiles…", flush=True)
    df = pd.DataFrame([_row(a, org) for a in accs])

    mob = D.load_mobidb(org)
    if mob.empty:
        print(f"[{org}] [warn] no MobiDB bulk — sequence-side disorder and domain architecture "
              f"will be absent (run: python scripts/10a_fetch_degradability.py --only mobidb)",
              flush=True)
    else:
        df = df.merge(mob, on="uniprot_accession", how="left")
        # exactly two structural domains = the Regime B `two_domain_architecture` feature
        df["two_domain_architecture"] = df["mobidb_n_domains_merge"].fillna(-1).eq(2)
        print(f"[{org}] MobiDB joined: {df['mobidb_dis_dis465'].notna().sum()}/{len(df)} accessions",
              flush=True)

    df["disorder_architecture"] = df.apply(_architecture, axis=1)

    n_af = int(df["af_available"].sum())
    print(f"[{org}] AlphaFold models {n_af}/{len(df)} ({n_af/len(df):.1%}); "
          f"median global disorder {df['global_disorder_frac_plddt'].median():.3f}", flush=True)
    for t, lab in (("n", "N"), ("c", "C")):
        cnt = df[f"{t}term_init_bucket"].value_counts()
        print(f"[{org}]   {lab}-term buckets: "
              + "  ".join(f"{b}={int(cnt.get(b, 0))}" for b in D.INIT_BUCKETS), flush=True)
    eng = df["engageable_terminus"].value_counts()
    print(f"[{org}]   engageable (>={ENGAGE_AA} aa): "
          + "  ".join(f"{k}={int(v)}" for k, v in eng.items()), flush=True)
    print(f"[{org}]   architecture: "
          + "  ".join(f"{k}={int(v)}" for k, v in df['disorder_architecture'].value_counts().items()),
          flush=True)
    print(f"[{org}]   internal degron candidates (>= {D.MIN_INTERNAL_RUN} aa): "
          f"{int(df['internal_degron_candidate'].sum())}", flush=True)

    # Cross-check against the ligandability disorder_frac — same definition, so a mismatch is a bug.
    lig = D.results_dir(org) / f"{prefix}_ligandability.csv"
    if lig.exists():
        L = pd.read_csv(lig, usecols=["uniprot_accession", "disorder_frac"])
        chk = df.merge(L, on="uniprot_accession", how="inner").dropna(
            subset=["global_disorder_frac_plddt", "disorder_frac"])
        delta = (chk["global_disorder_frac_plddt"] - chk["disorder_frac"]).abs()
        print(f"[{org}] cross-check vs ligandability disorder_frac: n={len(chk)} "
              f"max|Δ|={delta.max():.4f} mean|Δ|={delta.mean():.5f}"
              + ("  OK" if delta.max() < 0.02 else "  <-- MISMATCH, investigate"), flush=True)

    tab = _validate(df, org)

    out = D.results_dir(org) / f"{prefix}_disorder.csv"
    df.to_csv(out, index=False)
    print(f"\n[{org}] wrote {out.relative_to(D.REPO_ROOT)}  ({len(df)} rows, {df.shape[1]} cols)",
          flush=True)
    if not tab.empty:
        vout = D.degradability_processed_dir(org) / f"{prefix}_disorder_validation.csv"
        tab.to_csv(vout, index=False)
        print(f"[{org}] wrote {vout.relative_to(D.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
