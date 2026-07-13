"""2x3 summary of the composite essentiality annotation (07h).

Same house style as the ligandability summary (06h): stylia "slide" format, NPG palette. Every panel
is specific to the single `--organism` of the slide.

  1  evidence sources          # proteins supported by each track (ECL8 / E. coli / ProteomeLM / Geptop / FBA)
  2  essentiality tiers        essential / likely_essential / non_essential
  3  sub-score coverage        proteins with an experimental / transfer / predictor sub-score
  4  composite score           essentiality_score histogram, split by tier
  5  evidence agreement        # independent evidence sources per essential-tier protein
  6  top essential targets     highest-scoring proteins, shaded by evidence breadth

Reads output/results/<org>/<prefix>_essentiality.csv (07h).
Output: output/plots/07j_essentiality_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

REPO_ROOT = E.REPO_ROOT
NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
TIER_COLOR = {"essential": NPG[2], "likely_essential": NPG[4], "non_essential": "#C9C9C7"}
TOP_N = 15


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    orgname = E.ORG_DISPLAY[org]
    d = pd.read_csv(E.results_dir(org) / f"{prefix}_essentiality.csv")

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: evidence sources ----
    ax = axs.next()
    src = d["evidence_sources"].fillna("")
    tracks = [("Kp ECL8", src.str.contains("Kp_ECL8")),
              ("E. coli", src.str.contains("Ecoli")),
              ("ProteomeLM", src.str.contains("ProteomeLM")),
              ("Geptop", src.str.contains("Geptop")),
              ("FBA", src.str.contains(r"\bFBA\b"))]
    counts = [int(m.sum()) for _, m in tracks]
    colors = [NPG[0], NPG[1], "#00A087", "#3C5488", "#F39B7F"]
    bars = ax.bar(range(len(tracks)), counts, color=colors)
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(len(tracks))); ax.set_xticklabels([t for t, _ in tracks], fontsize=SS, rotation=20)
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Essentiality evidence — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 2: tiers ----
    ax = axs.next()
    order = ["essential", "likely_essential", "non_essential"]
    vc = d["essentiality_tier"].value_counts().reindex(order).fillna(0).astype(int)
    bars = ax.bar(range(3), vc.to_numpy(), color=[TIER_COLOR[t] for t in order])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels(["essential", "likely", "non-essential"], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Essentiality tiers — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 3: sub-score coverage ----
    ax = axs.next()
    subs = [("experimental", "evidence_experimental"), ("E. coli\ntransfer", "evidence_transfer"),
            ("predictor", "evidence_predictor")]
    cov = [int(pd.to_numeric(d[c], errors="coerce").notna().sum()) for _, c in subs]
    bars = ax.bar(range(3), cov, color=[NPG[0], NPG[1], "#00A087"])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels([t for t, _ in subs], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="proteins with sub-score", title=f"Sub-score coverage — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 4: composite score histogram by tier ----
    ax = axs.next()
    bins = np.linspace(0, 1, 26)
    bottom = np.zeros(len(bins) - 1)
    for t in order:
        h, _ = np.histogram(d.loc[d["essentiality_tier"] == t, "essentiality_score"], bins=bins)
        ax.bar(bins[:-1], h, width=np.diff(bins), bottom=bottom, align="edge",
               color=TIER_COLOR[t], label=t.replace("_", " "))
        bottom += h
    ax.set_yscale("log")
    stylia.label(ax, xlabel="essentiality score", ylabel="proteins", title=f"Composite score — {orgname}")
    ax.legend(fontsize=SS, frameon=False)

    # ---- panel 5: evidence breadth for essential-tier proteins ----
    ax = axs.next()
    ess = d[d["essentiality_tier"] == "essential"].copy()
    nsrc = ess["evidence_sources"].fillna("").apply(
        lambda s: len({x for x in s.replace("ECL8", "Kp_ECL8").split(";")
                       if x in {"Kp_ECL8", "Ecoli", "ProteomeLM", "Geptop", "FBA"}}))
    vc = nsrc.value_counts().reindex([1, 2, 3, 4, 5]).fillna(0).astype(int)
    bars = ax.bar([1, 2, 3, 4, 5], vc.to_numpy(), color=NPG[2])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks([1, 2, 3, 4, 5])
    stylia.label(ax, xlabel="# independent evidence sources", ylabel="essential-tier proteins",
                 title=f"Evidence convergence — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 6: top essential targets ----
    ax = axs.next()
    top = d.sort_values("essentiality_score", ascending=False).head(TOP_N).iloc[::-1]
    labels = [g if isinstance(g, str) and g else a for g, a in zip(top["gene"], top["uniprot_accession"])]
    nsrc_top = top["evidence_sources"].fillna("").apply(
        lambda s: len({x for x in s.replace("ECL8", "Kp_ECL8").split(";")
                       if x in {"Kp_ECL8", "Ecoli", "ProteomeLM", "Geptop", "FBA"}}))
    import matplotlib.cm as cm
    colors = [cm.get_cmap("YlGn")(0.35 + 0.13 * n) for n in nsrc_top]
    ax.barh(range(len(top)), top["essentiality_score"], color=colors)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(labels, fontsize=SS)
    ax.set_xlim(0.9, 1.001)
    stylia.label(ax, xlabel="essentiality score", ylabel="", title=f"Top essential targets — {orgname}")

    out = REPO_ROOT / "output" / "plots" / f"07j_essentiality_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
