"""2x3 overview of the essentiality PREDICTORS (07d ProteomeLM, 07e Geptop, 07f FBA).

Same house style as the ligandability slides (06i etc.): stylia "slide" format, NPG palette. Every
panel is specific to the single `--organism` of the slide.

  1  predictor score distributions   ProteomeLM vs Geptop score histograms
  2  ProteomeLM vs Geptop            per-protein agreement scatter, coloured by E. coli-essential
  3  FBA metabolic knockouts         in-model essential / non-essential / not-in-model
  4  predictor agreement             how many of the 3 predictors call each protein essential
  5  ProteomeLM calibration          predicted score vs observed E. coli-essential fraction (+CV AUROC)
  6  predictors vs experiment        concordance with the direct call (Kp ECL8 / E. coli EcoGene)

Reads output/results/<org>/<prefix>_ess_{proteomelm,geptop,fba,ecoli,kp}.csv.
Output: output/plots/07i_predictors_<prefix>.png (one slide per --organism). Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

REPO_ROOT = E.REPO_ROOT
NPG = stylia.CategoricalPalette("npg").colors
ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}
PLM_C, GEP_C, FBA_C = "#00A087", "#3C5488", "#F39B7F"  # NPG green / navy / salmon
SS = stylia.SLIDE_FONTSIZE_SMALL


def load(org: str) -> pd.DataFrame:
    _, prefix = E.ORGANISMS[org]
    d = pd.DataFrame({"uniprot_accession": E.load_accessions(org)})
    for suffix in ("proteomelm", "geptop", "fba", "ecoli", "kp"):
        p = E.results_dir(org) / f"{prefix}_ess_{suffix}.csv"
        if p.exists():
            d = d.merge(pd.read_csv(p).drop_duplicates("uniprot_accession"), on="uniprot_accession", how="left")
    d["plm"] = pd.to_numeric(d.get("proteomelm_ess_score"), errors="coerce")
    d["gep"] = pd.to_numeric(d.get("geptop_score"), errors="coerce")
    d["fba_ess"] = d.get("fba_essential").map({True: 1, False: 0, "True": 1, "False": 0})
    d["ec_ess"] = d.get("ec_transfer_essential", pd.Series(False, index=d.index)).fillna(False).astype(bool)
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    orgname = E.ORG_DISPLAY[org]
    d = load(org)

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: predictor score distributions ----
    ax = axs.next()
    bins = np.linspace(0, 1, 26)
    ax.hist(d["plm"].dropna(), bins=bins, color=PLM_C, alpha=0.6, label="ProteomeLM")
    ax.hist(d["gep"].dropna(), bins=bins, color=GEP_C, alpha=0.5, label="Geptop 2.0")
    ax.axvline(E.GEPTOP_CUTOFF, color=GEP_C, ls="--", lw=1)
    ax.axvline(0.5, color=PLM_C, ls="--", lw=1)
    stylia.label(ax, xlabel="essentiality score", ylabel="proteins", title=f"Predictor scores — {orgname}")
    ax.set_yscale("log")
    ax.legend(fontsize=SS, frameon=False)

    # ---- panel 2: ProteomeLM vs Geptop agreement ----
    ax = axs.next()
    for lab, mask, c in [("E. coli-essential", d["ec_ess"], ORG_COLOR[prefix]),
                         ("other", ~d["ec_ess"], "#C9C9C7")]:
        s = d[mask]
        ax.scatter(s["gep"], s["plm"], s=6, alpha=0.5 if lab == "other" else 0.8,
                   color=c, linewidths=0, rasterized=True, label=lab)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="Geptop score", ylabel="ProteomeLM score", title=f"Predictor agreement — {orgname}")
    r = d[["gep", "plm"]].corr().iloc[0, 1]
    ax.legend(fontsize=SS, frameon=False, loc="lower right")

    # ---- panel 3: FBA metabolic knockouts ----
    ax = axs.next()
    fs = d.get("fba_status", pd.Series("not_in_model", index=d.index)).fillna("not_in_model")
    cats = ["FBA-essential", "in-model,\nnon-essential", "not in\nmodel"]
    counts = [int((d["fba_ess"] == 1).sum()),
              int(((fs == "in_model") & (d["fba_ess"] == 0)).sum()),
              int((fs == "not_in_model").sum())]
    bars = ax.bar(range(3), counts, color=[FBA_C, "#F7C6B4", "#D8D8D6"])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels(cats, fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"FBA single-gene KO — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 4: predictor agreement (# predictors calling essential) ----
    ax = axs.next()
    plm_call = (d["plm"] >= 0.5).astype(int)
    gep_call = (d["gep"] >= E.GEPTOP_CUTOFF).astype(int)
    fba_call = (d["fba_ess"] == 1).astype(int)
    nvote = plm_call.fillna(0) + gep_call.fillna(0) + fba_call.fillna(0)
    vc = nvote.value_counts().reindex([0, 1, 2, 3]).fillna(0).astype(int)
    bars = ax.bar([0, 1, 2, 3], vc.to_numpy(), color=["#D8D8D6", NPG[4], NPG[1], NPG[2]])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["0", "1", "2", "3"])
    stylia.label(ax, xlabel="# predictors calling essential", ylabel="proteins",
                 title=f"Predictor consensus — {orgname}")
    ax.set_yscale("log"); ax.margins(y=0.2)

    # ---- panel 5: ProteomeLM calibration vs E. coli labels ----
    ax = axs.next()
    bins = np.linspace(0, 1, 11)
    mid = (bins[:-1] + bins[1:]) / 2
    frac = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (d["plm"] >= lo) & (d["plm"] < hi)
        frac.append(d.loc[m, "ec_ess"].mean() if m.sum() else np.nan)
    ax.plot([0, 1], [0, 1], ls=":", color="#999")
    ax.plot(mid, frac, "-o", color=PLM_C, ms=4)
    stylia.label(ax, xlabel="ProteomeLM score", ylabel="observed E. coli-essential frac",
                 title=f"ProteomeLM calibration — {orgname}")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # ---- panel 6: predictors vs the direct experimental call ----
    ax = axs.next()
    if org == "kpneumoniae":
        direct = d.get("kp_ess_in_vitro_call", pd.Series(index=d.index)).astype(str).eq("essential")
        dlabel = "Kp ECL8-essential"
    else:
        direct = d["ec_ess"]
        dlabel = "E. coli EcoGene-essential"
    grp = d[direct]
    labels = ["ProteomeLM", "Geptop", "FBA"]
    recovered = [
        float((grp["plm"] >= 0.5).mean()),
        float((grp["gep"] >= E.GEPTOP_CUTOFF).mean()),
        float((grp["fba_ess"] == 1).mean()) if grp["fba_ess"].notna().any() else 0.0,
    ]
    bars = ax.bar(range(3), recovered, color=[PLM_C, GEP_C, FBA_C])
    ax.bar_label(bars, labels=[f"{v:.0%}" for v in recovered], padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=SS)
    ax.set_ylim(0, 1.15)  # headroom so a 100% bar + its label clear the top frame
    stylia.label(ax, xlabel="", ylabel=f"frac of {dlabel} recovered",
                 title=f"Predictor sensitivity — {orgname}")

    out = REPO_ROOT / "output" / "plots" / f"07i_predictors_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
