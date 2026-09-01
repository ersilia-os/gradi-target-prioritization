"""Stage 09g — merge the localization tracks into the axis table (docs §5.1).

Combines 09a (UniProt), 09b (PSORTb 3.0 via PSORTdb), 09c (DeepLocPro), 09d (TMbed topology),
09e (signal-peptide type / lipoprotein sorting) and 09f (STEPdb + ortholog transfer) into one call
per protein, with an explicit account of where that call came from.

**Precedence** — evidence beats curation, curation beats prediction:

    1. uniprot_experimental   UniProt CC with an experimental ECO code
    2. ortholog_transfer      E. coli experimental call carried over an OrthoFinder orthogroup (Kp)
       stepdb                 STEPdb 2.0 assignment with a cited reference (Ec)
    3. stepdb_curated         STEPdb assignment without a reference
    4. uniprot_curated        UniProt CC without experimental evidence
    5. deeplocpro             the ML predictor
    6. psortb                 the rule-based predictor (tiebreak / corroboration)
    7. unknown

**Topology overrides.** TMbed and SignalP fix two specific, well-documented failure modes of
compartment predictors, so they are allowed to overrule a *predicted* call (never an experimental
or curated one):

  * a transmembrane β-barrel is an outer-membrane protein — mis-calling one as cytoplasmic is the
    worst-case error for this axis, because it would hand a surface protein a Clp-accessibility
    of 1.0;
  * a Sec/SPII lipoprotein is sorted by the Lol "+2 rule", which is what PSORTb characteristically
    gets wrong (it bins outer-membrane lipoproteins as periplasmic).

The provisional `membrane` label (UniProt "Membrane" with no side stated) is resolved from topology
regardless of source, since it is not a final answer in the first place.

**`clp_accessibility`** is then the graded ladder from `src/localization.py` — the key refinement
over the previous flat 1.0/0.5/0.0 stub being that an inner-membrane protein is scored on how much
of it actually faces the cytoplasm, since that is what ClpXP can engage.

Outputs: output/results/<organism>/<prefix>_localization.csv
         output/results/<organism>/<prefix>_localization_shortlist.csv
Run with the `gradi` conda env.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

# Confidence assigned to a call by the tier of source that produced it. Predictor calls carry
# their own probability instead (DeepLocPro's softmax, PSORTb's final score / 10).
TIER_CONFIDENCE = {"experimental": 1.0, "curated": 0.8}
PSORTB_SCORE_MAX = 10.0

# The independent opinions compared to compute `predictor_agreement`.
AGREEMENT_TRACKS = ("uniprot_localization", "dlp_localization", "psortb_localization")


def _read(path: Path) -> pd.DataFrame:
    """Defensive read — a missing track degrades the merge, it does not break it."""
    if not path.exists():
        print(f"  [warn] missing {path.name} — that track will not contribute", flush=True)
        return pd.DataFrame(columns=["uniprot_accession"])
    return pd.read_csv(path)


def build(org: str) -> pd.DataFrame:
    _, prefix = LOC.ORGANISMS[org]
    rdir = LOC.results_dir(org)

    df = pd.DataFrame({"uniprot_accession": LOC.load_accessions(org)})
    for name in ("uniprot", "psortdb", "deeplocpro", "topology", "signalp", "transfer"):
        part = _read(rdir / f"{prefix}_loc_{name}.csv")
        if len(part.columns) > 1:
            df = df.merge(part, on="uniprot_accession", how="left")

    for col, default in (
        ("uniprot_localization", "unknown"), ("uniprot_evidence", "none"),
        ("dlp_localization", np.nan), ("psortb_localization", np.nan),
        ("transfer_localization", np.nan), ("transfer_evidence", np.nan),
        ("is_beta_barrel", False), ("is_lipoprotein", False),
        ("cyto_residue_fraction", np.nan), ("lipoprotein_sorting", np.nan),
    ):
        if col not in df.columns:
            df[col] = default

    df["is_beta_barrel"] = df["is_beta_barrel"].fillna(False).astype(bool)
    df["is_lipoprotein"] = df["is_lipoprotein"].fillna(False).astype(bool)

    records = []
    for r in df.itertuples():
        candidates: list[tuple[str, str]] = []          # (source, label)

        uni = getattr(r, "uniprot_localization", "unknown")
        uni_ev = getattr(r, "uniprot_evidence", "none")
        if isinstance(uni, str) and uni != "unknown":
            candidates.append(
                ("uniprot_experimental" if uni_ev == "experimental" else "uniprot_curated", uni))

        tr = getattr(r, "transfer_localization", None)
        tr_ev = getattr(r, "transfer_evidence", None)
        if isinstance(tr, str) and tr not in ("unknown", ""):
            src = getattr(r, "transfer_source", "stepdb")
            candidates.append((src if tr_ev == "experimental" else "stepdb_curated", tr))

        dlp = getattr(r, "dlp_localization", None)
        if isinstance(dlp, str) and dlp != "unknown":
            candidates.append(("deeplocpro", dlp))

        ps = getattr(r, "psortb_localization", None)
        if isinstance(ps, str) and ps != "unknown":
            candidates.append(("psortb", ps))

        ranked = {src: lbl for src, lbl in candidates}
        source, label = "none", "unknown"
        for src in LOC.SOURCE_PRECEDENCE:
            if src in ranked:
                source, label = src, ranked[src]
                break

        evidence = LOC.SOURCE_EVIDENCE.get(source, "predicted")
        beta = bool(r.is_beta_barrel)
        cyto = getattr(r, "cyto_residue_fraction", np.nan)
        cyto = None if (cyto is None or (isinstance(cyto, float) and np.isnan(cyto))) else float(cyto)

        # `membrane` states no side; resolve it from topology whatever the source.
        override = ""
        if label == "membrane":
            label = "outer_membrane" if beta else "inner_membrane"
            override = "membrane_resolved"

        # Topology may overrule a *predicted* call, never an experimental or curated one.
        if evidence == "predicted":
            if beta and label != "outer_membrane":
                label, override = "outer_membrane", "beta_barrel"
            elif r.is_lipoprotein and isinstance(getattr(r, "lipoprotein_sorting", None), str):
                sorted_to = r.lipoprotein_sorting
                if sorted_to != label:
                    label, override = sorted_to, "lipoprotein_lol_rule"

        if source == "deeplocpro":
            conf = getattr(r, "dlp_confidence", np.nan)
        elif source == "psortb":
            score = getattr(r, "psortb_score", np.nan)
            conf = (float(score) / PSORTB_SCORE_MAX
                    if score is not None and not (isinstance(score, float) and np.isnan(score))
                    else np.nan)
        else:
            conf = TIER_CONFIDENCE.get(evidence, np.nan)

        opinions = [getattr(r, t, None) for t in AGREEMENT_TRACKS]
        opinions = [o for o in opinions if isinstance(o, str) and o != "unknown"]
        agreement = (round(sum(o == label for o in opinions) / len(opinions), 3)
                     if opinions else None)

        records.append({
            "uniprot_accession": r.uniprot_accession,
            "localization": label,
            "localization_confidence": None if conf is None or (isinstance(conf, float) and np.isnan(conf)) else round(float(conf), 3),
            "localization_source": source,
            "localization_evidence": evidence,
            "localization_override": override,
            "predictor_agreement": agreement,
            "n_localization_sources": len(candidates),
            "clp_accessibility": LOC.clp_accessibility(label, cyto, beta),
        })

    merged = pd.DataFrame(records).merge(df, on="uniprot_accession", how="left")
    return merged


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    print(f"[{org}] merging localization tracks", flush=True)

    df = build(org)
    out = LOC.results_dir(org) / f"{prefix}_localization.csv"
    df.to_csv(out, index=False)

    short = df[(df["clp_accessibility"] >= LOC.MIN_CLP_ACCESSIBILITY)
               & (df["localization"] != "unknown")]
    short_out = LOC.results_dir(org) / f"{prefix}_localization_shortlist.csv"
    short.to_csv(short_out, index=False)

    n = len(df)
    known = int((df["localization"] != "unknown").sum())
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({n} proteins)", flush=True)
    print(f"[{org}] localized: {known}/{n} ({known / n:.1%})", flush=True)
    print(f"[{org}] classes: {df['localization'].value_counts().to_dict()}", flush=True)
    print(f"[{org}] evidence: {df['localization_evidence'].value_counts().to_dict()}", flush=True)
    print(f"[{org}] source: {df['localization_source'].value_counts().to_dict()}", flush=True)
    ov = df[df["localization_override"] != ""]["localization_override"].value_counts().to_dict()
    print(f"[{org}] overrides: {ov}", flush=True)
    print(f"[{org}] clp-accessible shortlist (>= {LOC.MIN_CLP_ACCESSIBILITY}): "
          f"{len(short)} -> {short_out.relative_to(LOC.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
