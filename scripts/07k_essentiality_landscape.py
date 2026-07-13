"""Capstone: the essentiality landscape & cross-axis target prioritization (docs §4, result).

Ties the 07h composite to the other axes: the ESM-C protein-universe map, the 06g ligandability
score, the 05a bibliometric studiedness, and the 03c broad-spectrum-selectivity call. This is where
essentiality becomes a target-prioritization statement: broad-spectrum, human-selective, essential AND
ligandable = a prime antibacterial target; essential AND understudied = a neglected opportunity.

  1  essentiality landscape     ESM-C map coloured by essentiality tier
  2  essential & ligandable     essentiality vs ligandability; prime (broad-selective) highlighted
  3  broad-spectrum essentials  Enterobacteriaceae %essential for essential-tier proteins
  4  prime target set           essential ∧ ligandable ∧ broad-selective, on the ESM-C map
  5  neglected & essential      studiedness (05a) vs essentiality; understudied + essential = opportunity
  6  top prime targets          highest essential∧ligandable prime targets, shaded by studiedness

Reads output/results/<org>/<prefix>_{essentiality,ligandability,esmc600m_projection}.csv and
data/processed/<org>/bibliometric/<prefix>_popularity.csv.
Output: output/plots/07k_landscape_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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
TIER_COLOR = {"essential": NPG[2], "likely_essential": NPG[4], "non_essential": "#D8D8D6"}
PRIME_C = "#E64B35"
TOP_N = 12


def load(org: str) -> pd.DataFrame:
    _, prefix = E.ORGANISMS[org]
    d = pd.read_csv(E.results_dir(org) / f"{prefix}_essentiality.csv")
    lg = E.results_dir(org) / f"{prefix}_ligandability.csv"
    if lg.exists():
        d = d.merge(pd.read_csv(lg)[["uniprot_accession", "ligandability_score", "ligandability_tier"]],
                    on="uniprot_accession", how="left")
    em = E.results_dir(org) / f"{prefix}_esmc600m_projection.csv"
    if em.exists():
        d = d.merge(pd.read_csv(em)[["uniprot_accession", "tsne_x", "tsne_y"]],
                    on="uniprot_accession", how="left")
    pop = REPO_ROOT / "data" / "processed" / org / "bibliometric" / f"{prefix}_popularity.csv"
    if pop.exists():
        d = d.merge(pd.read_csv(pop)[["uniprot_accession", "popularity_score"]],
                    on="uniprot_accession", how="left")
    for c in ("ligandability_score", "popularity_score", "tsne_x", "tsne_y"):
        if c not in d:
            d[c] = np.nan
    d["is_essential"] = d["essentiality_tier"] == "essential"
    d["is_broad"] = d.get("selectivity") == "broad_selective"
    d["is_ligandable"] = d.get("ligandability_tier").isin(["tractable"]) if "ligandability_tier" in d else False
    d["is_prime"] = d["is_essential"] & d["is_broad"] & d["is_ligandable"]
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    orgname = E.ORG_DISPLAY[org]
    d = load(org)
    has_xy = d["tsne_x"].notna() & d["tsne_y"].notna()

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: essentiality on the ESM-C map ----
    ax = axs.next()
    for t in ["non_essential", "likely_essential", "essential"]:
        s = d[has_xy & (d["essentiality_tier"] == t)]
        ax.scatter(s["tsne_x"], s["tsne_y"], s=4, alpha=0.35 if t == "non_essential" else 0.8,
                   color=TIER_COLOR[t], linewidths=0, rasterized=True, label=t.replace("_", " "))
    stylia.label(ax, xlabel="ESM-C dim 1", ylabel="ESM-C dim 2", title=f"Essentiality landscape — {orgname}")
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(fontsize=SS, frameon=False, markerscale=2, loc="best")

    # ---- panel 2: essential vs ligandable ----
    ax = axs.next()
    bg = d[~d["is_prime"]]
    pr = d[d["is_prime"]]
    ax.scatter(bg["essentiality_score"], bg["ligandability_score"], s=4, alpha=0.25,
               color="#C9C9C7", linewidths=0, rasterized=True)
    ax.scatter(pr["essentiality_score"], pr["ligandability_score"], s=10, alpha=0.8,
               color=PRIME_C, linewidths=0, rasterized=True, label=f"prime (n={len(pr)})")
    ax.axvline(0.6, color="#999", ls=":", lw=1); ax.axhline(0.6, color="#999", ls=":", lw=1)
    stylia.label(ax, xlabel="essentiality score", ylabel="ligandability score",
                 title=f"Essential & ligandable — {orgname}")
    ax.legend(fontsize=SS, frameon=False, loc="lower left")

    # ---- panel 3: broad-spectrum conservation of essentials ----
    ax = axs.next()
    ep = pd.to_numeric(d.loc[d["is_essential"], "entero_pct_essential"], errors="coerce").dropna()
    if len(ep):
        ax.hist(ep, bins=np.linspace(0, 1, 21), color=NPG[1])
        ax.axvline(0.8, color="#555", ls="--", lw=1)
    stylia.label(ax, xlabel="Enterobacteriaceae %essential", ylabel="essential-tier proteins",
                 title=f"Broad-spectrum essentials — {orgname}")

    # ---- panel 4: prime target set on the map ----
    ax = axs.next()
    ax.scatter(d.loc[has_xy, "tsne_x"], d.loc[has_xy, "tsne_y"], s=4, alpha=0.3,
               color="#D8D8D6", linewidths=0, rasterized=True)
    p = d[has_xy & d["is_prime"]]
    ax.scatter(p["tsne_x"], p["tsne_y"], s=12, alpha=0.9, color=PRIME_C, linewidths=0, rasterized=True)
    stylia.label(ax, xlabel="ESM-C dim 1", ylabel="ESM-C dim 2",
                 title=f"Prime targets (essential & ligandable & broad) — {orgname}")
    ax.set_xticks([]); ax.set_yticks([])

    # ---- panel 5: neglected & essential ----
    ax = axs.next()
    sub = d[d["popularity_score"].notna()]
    ax.scatter(sub["popularity_score"], sub["essentiality_score"], s=4, alpha=0.25,
               color="#C9C9C7", linewidths=0, rasterized=True)
    prs = sub[sub["is_prime"]]
    ax.scatter(prs["popularity_score"], prs["essentiality_score"], s=8, alpha=0.75,
               color=PRIME_C, linewidths=0, rasterized=True)
    px = round(float(sub["popularity_score"].quantile(0.33)), 2)
    ax.axvline(px, color="#999", ls=":", lw=1); ax.axhline(0.6, color="#999", ls=":", lw=1)
    opp = prs[(prs["popularity_score"] <= px) & (prs["essentiality_score"] >= 0.6)] \
        .sort_values("essentiality_score", ascending=False)
    for _, r in opp.head(6).iterrows():
        lab = r["gene"] if isinstance(r["gene"], str) and r["gene"] else r["uniprot_accession"]
        ax.annotate(lab, (r["popularity_score"], r["essentiality_score"]),
                    fontsize=7, color="#7a2417", xytext=(2, 2), textcoords="offset points")
    stylia.label(ax, xlabel="studiedness (05a)", ylabel="essentiality score",
                 title=f"Neglected & essential — {orgname}")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)

    # ---- panel 6: top prime targets ----
    ax = axs.next()
    pr2 = d[d["is_prime"]].copy()
    pr2["combined"] = pr2["essentiality_score"].fillna(0) * pr2["ligandability_score"].fillna(0)
    top = pr2.sort_values("combined", ascending=False).head(TOP_N).iloc[::-1]
    if len(top):
        labels = [g if isinstance(g, str) and g else a for g, a in zip(top["gene"], top["uniprot_accession"])]
        import matplotlib.cm as cm
        pops = top["popularity_score"].fillna(0.0)
        colors = [cm.get_cmap("YlOrRd_r")(0.2 + 0.6 * p) for p in pops]
        ax.barh(range(len(top)), top["combined"], color=colors)
        ax.set_yticks(range(len(top))); ax.set_yticklabels(labels, fontsize=SS)
        stylia.label(ax, xlabel="essentiality × ligandability", ylabel="",
                     title=f"Top prime targets — {orgname}")

    out = REPO_ROOT / "output" / "plots" / f"07k_landscape_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)} (prime targets: {int(d['is_prime'].sum())})")


if __name__ == "__main__":
    main()
