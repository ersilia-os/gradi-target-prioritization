"""Stage 10e — measured protein turnover and protease attribution (docs §3.3).

Ingests the three E. coli datasets that actually *measure* degradation, as opposed to predicting it
from sequence or structure. Everything else in the degradability axis so far (`10b` degrons, `10d`
disorder) is computed; this is the layer that can falsify it.

  * **Nagar 2021** (mSystems, pulsed-SILAC) — 1,149 measured half-lives, 3 classes. `src/degradability.py`
    already parses this via `load_nagar()`; it has only ever been used in-memory as a validation
    label, never materialised, so this script writes it out for the first time.
  * **Gupta 2024** (Nat Commun) — half-lives across 13 growth conditions **plus a protease-knockout
    panel** (ΔclpP / Δlon / ΔhslV / triple / ΔsmpB at N-lim6). The only dataset in the repo that can
    say a protein is degraded *by ClpP* rather than merely degraded.
  * **Niwa 2022** (Molecules) — GroE-depletion proteomics run in WT, Δlon, ΔclpPX and ΔhslVU
    backgrounds: the only head-to-head of the three cytoplasmic protease systems.

**Two derivations here are ours, not the papers'**, and are named as such in the output columns:

  * `gupta_log2_<protease>` = log2(KO half-life / WT half-life) at N-lim6. A protein stabilised in a
    knockout was being degraded by that protease. Gupta's published substrate categories are not a
    column in Table S1, so the attribution is recomputed from the panel with an explicit threshold.
  * `niwa_rescue_<protease>` = Foldchange(KO) − Foldchange(WT). Niwa's fold-changes are GroE-on vs
    GroE-off *within* each background, so a protein whose depletion is blunted in a protease KO was
    being cleared by that protease once GroE stopped protecting it.

K. pneumoniae has **no measured turnover data of its own** — not one protein. Its values are
transferred from E. coli across the 03a orthogroups, which caps coverage at 3,179/5,728 (55.5%), and
every transferred row carries `measured_source` / `measured_donor` so the distinction survives into
the plots.

Output: output/results/<organism>/<prefix>_deg_measured.csv
Run with the `gradi` conda env.
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402
from src import localization as LOC  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

GUPTA = "gupta2024_turnover/gupta2024_TableS1_turnover.xlsx"
NIWA = "niwa2022_lon_clpxp/Supplementary_Dataset_S1.xlsx"

# Gupta Table S1 carries two header rows: long descriptions then short names. The short names are
# the ones worth keying on.
GUPTA_HEADER_ROW = 4
GUPTA_WT_COL = "N-lim6_mean"          # wild-type reference the KO panel was run against
GUPTA_KO = {                          # output suffix -> Table S1 column
    "clpP": "clpPN-lim6_mean",
    "lon": "lonN-lim6_mean",
    "hslV": "hslVN-lim6_mean",
    "triple": "TripleN-lim6_mean",
    "smpB": "smpBN-lim6_mean",
}
# A knockout must lengthen the half-life by this much (log2) before we call the protein its
# substrate. 0.5 = 1.41x stabilised.
#
# NOTE this is more permissive than Gupta's own published categories (ClpP 64 / Lon 14 / HslV 1 /
# additive 82 / redundant 41); at 0.5 we get ClpP 189 / Lon 130 / HslV 29 / additive 74 /
# redundant 224. The paper applied significance testing across replicates that Table S1 does not
# expose, so the sets are not directly comparable — treat these as a recomputed, wider net, and use
# the continuous `gupta_log2_*` columns rather than the category when the distinction matters.
STABILISATION_LOG2 = 0.5

NIWA_WT = "Foldchange_MGM100"
NIWA_KO = {"lon": "Foldchange_MGM100dlon",
           "clpxp": "Foldchange_MGM100dPX",
           "hslvu": "Foldchange_MGM100dVU"}
NIWA_RESCUE_LOG2 = 0.5

_SP = re.compile(r"^[a-z]{2}\|([A-Z0-9]+)\|")


def _acc(value: str) -> str | None:
    """`sp|A5A614|YCIZ_ECOLI` -> `A5A614`; pass through a bare accession."""
    s = str(value).strip()
    m = _SP.match(s)
    if m:
        return m.group(1)
    return s.split(";")[0].strip() or None


def load_gupta(canon: set[str]) -> pd.DataFrame:
    path = D.degradability_raw_dir("ecoli") / GUPTA
    if not path.exists():
        print(f"  [warn] missing {path.name} — Gupta track skipped", flush=True)
        return pd.DataFrame(columns=["uniprot_accession"])
    g = pd.read_excel(path, sheet_name="TableS1", header=GUPTA_HEADER_ROW)
    g["uniprot_accession"] = g["Protein ID"].map(_acc)
    g = g[g["uniprot_accession"].isin(canon)].copy()

    out = pd.DataFrame({"uniprot_accession": g["uniprot_accession"]})
    out["gupta_halflife_hrs"] = pd.to_numeric(g.get("MinimalMedia_mean"), errors="coerce")
    wt = pd.to_numeric(g.get(GUPTA_WT_COL), errors="coerce")
    out["gupta_halflife_hrs_nlim6"] = wt

    # How many of the 13 growth conditions this protein was measured in.
    mean_cols = [c for c in g.columns if str(c).endswith("_mean") and not any(
        str(c).startswith(k) for k in ("clpP", "lon", "hslV", "Triple", "smpB"))]
    out["gupta_n_conditions"] = g[mean_cols].notna().sum(axis=1).to_numpy()

    for name, col in GUPTA_KO.items():
        ko = pd.to_numeric(g.get(col), errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"gupta_log2_{name}"] = np.log2(ko / wt)

    # Attribution: which single-protease knockout stabilised it. `redundant` = the triple knockout
    # stabilises it while no single one does, which is exactly what redundancy looks like.
    singles = ["clpP", "lon", "hslV"]
    mat = out[[f"gupta_log2_{s}" for s in singles]]
    measured = mat.notna().any(axis=1)
    hit = mat >= STABILISATION_LOG2
    n_hit = hit.sum(axis=1)
    best = pd.Series(index=out.index, dtype=object)
    best[measured] = mat[measured].idxmax(axis=1).str.replace("gupta_log2_", "", regex=False)
    triple_hit = out["gupta_log2_triple"] >= STABILISATION_LOG2

    attribution = pd.Series(
        np.where(n_hit == 1, best,
                 np.where(n_hit > 1, "additive",
                          np.where(triple_hit, "redundant", "unattributed"))),
        index=out.index)
    attribution[~measured] = None
    out["gupta_protease_attribution"] = attribution
    return out.dropna(subset=["uniprot_accession"]).drop_duplicates("uniprot_accession")


def load_niwa(canon: set[str]) -> pd.DataFrame:
    path = D.degradability_raw_dir("ecoli") / NIWA
    if not path.exists():
        print(f"  [warn] missing {path.name} — Niwa track skipped", flush=True)
        return pd.DataFrame(columns=["uniprot_accession"])
    n = pd.read_excel(path, sheet_name=0)
    n["uniprot_accession"] = n["UniProt_ID"].map(_acc)
    n = n[n["uniprot_accession"].isin(canon)].copy()

    out = pd.DataFrame({"uniprot_accession": n["uniprot_accession"]})
    wt = pd.to_numeric(n.get(NIWA_WT), errors="coerce")
    out["niwa_log2fc_wt"] = wt
    for name, col in NIWA_KO.items():
        ko = pd.to_numeric(n.get(col), errors="coerce")
        out[f"niwa_rescue_{name}"] = ko - wt

    rescue = out[[f"niwa_rescue_{k}" for k in NIWA_KO]]
    measured = rescue.notna().any(axis=1)
    hit = rescue >= NIWA_RESCUE_LOG2
    n_hit = hit.sum(axis=1)
    best = pd.Series(index=out.index, dtype=object)
    best[measured] = rescue[measured].idxmax(axis=1).str.replace("niwa_rescue_", "", regex=False)
    attribution = pd.Series(np.where(n_hit == 1, best,
                                     np.where(n_hit > 1, "multiple", "none")), index=out.index)
    attribution[~measured] = None
    out["niwa_attribution"] = attribution
    return out.dropna(subset=["uniprot_accession"]).drop_duplicates("uniprot_accession")


def load_nagar_frame(canon: set[str]) -> pd.DataFrame:
    n = D.load_nagar().copy()
    n = n[n["uniprot_accession"].isin(canon)]
    keep = {"halflife_min": "nagar_halflife_min",
            "halflife_class": "nagar_halflife_class",
            "kdeg_growth_corrected": "nagar_kdeg_growth_corrected",
            "halflife_proteolytic_min": "nagar_halflife_proteolytic_min"}
    out = n[["uniprot_accession"] + [c for c in keep if c in n.columns]].rename(columns=keep)
    return out.drop_duplicates("uniprot_accession")


NUMERIC_PREFIXES = ("nagar_halflife", "nagar_kdeg", "gupta_halflife", "gupta_log2",
                    "gupta_n_conditions", "niwa_log2fc", "niwa_rescue")
CATEGORICAL = ("nagar_halflife_class", "gupta_protease_attribution", "niwa_attribution")


def build_ecoli() -> pd.DataFrame:
    canon = set(LOC.load_accessions("ecoli"))
    frames = [load_nagar_frame(canon), load_gupta(canon), load_niwa(canon)]
    for name, f in zip(("Nagar", "Gupta", "Niwa"), frames):
        print(f"  {name}: {len(f)} proteins", flush=True)
    d = pd.DataFrame({"uniprot_accession": LOC.load_accessions("ecoli")})
    for f in frames:
        if len(f.columns) > 1:
            d = d.merge(f, on="uniprot_accession", how="left")
    d["measured_source"] = np.where(
        d[[c for c in d.columns if c != "uniprot_accession"]].notna().any(axis=1), "native", None)
    d["measured_donor"] = None
    return d


def build_kpneumoniae() -> pd.DataFrame:
    ec = pd.read_csv(LOC.results_dir("ecoli") / "ec_deg_measured.csv")
    accs = LOC.load_accessions("kpneumoniae")
    d = pd.DataFrame({"uniprot_accession": accs})

    value_cols = [c for c in ec.columns
                  if c not in ("uniprot_accession", "measured_source", "measured_donor")]
    categorical = [c for c in value_cols if c in CATEGORICAL]
    # CATEGORICAL must be subtracted explicitly: `nagar_halflife_class` also matches the
    # `nagar_halflife` numeric prefix, and averaging a string column fails at groupby time.
    numeric = [c for c in value_cols if c.startswith(NUMERIC_PREFIXES) and c not in CATEGORICAL]

    for c in numeric:
        series = pd.to_numeric(ec.set_index("uniprot_accession")[c], errors="coerce").dropna()
        series = series[np.isfinite(series)]
        d[c] = d["uniprot_accession"].map(D.transfer_ecoli_to_kp(series.to_dict(), reduce="mean"))

    donors: dict[str, str] = {}
    for c in categorical:
        labels = ec.dropna(subset=[c]).set_index("uniprot_accession")[c].astype(str).to_dict()
        transferred = LOC.transfer_categorical_ecoli_to_kp(labels)
        d[c] = d["uniprot_accession"].map(lambda a: transferred.get(a, {}).get("label"))
        for a, rec in transferred.items():
            donors.setdefault(a, rec["donors"])

    has_any = d[value_cols].notna().any(axis=1)
    d["measured_source"] = np.where(has_any, "ortholog_transfer", None)
    d["measured_donor"] = d["uniprot_accession"].map(donors)
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    print(f"[{org}] ingesting measured turnover", flush=True)

    d = build_ecoli() if org == "ecoli" else build_kpneumoniae()

    out = LOC.results_dir(org) / f"{prefix}_deg_measured.csv"
    d.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(d)} rows)", flush=True)

    n_any = int(d["measured_source"].notna().sum())
    print(f"[{org}] any measured evidence: {n_any}/{len(d)} ({n_any / len(d):.1%})", flush=True)
    for c in CATEGORICAL:
        if c in d.columns and d[c].notna().any():
            print(f"[{org}] {c}: {d[c].value_counts().to_dict()}", flush=True)
    for c in ("nagar_halflife_min", "nagar_halflife_proteolytic_min", "gupta_halflife_hrs"):
        if c in d.columns:
            v = d[c].replace([np.inf, -np.inf], np.nan)
            print(f"[{org}] {c}: {int(v.notna().sum())} finite", flush=True)


if __name__ == "__main__":
    main()
