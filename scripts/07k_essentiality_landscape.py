"""Capstone: essentiality target-prioritization slide (docs §4, result).

A single 6-panel synthesis that turns the graded essentiality score into a prioritization statement,
crossed with the ESM-C protein-universe map, the 06g ligandability score, the 03c broad-spectrum
selectivity call, and the 05a bibliometric studiedness. Replaces the former composite-summary +
landscape slides.

  1  global essentiality score   reverse-cumulative curve; top targets labelled in the sparse tail
  2  prioritization map          ESM-C 2-D map: all / essential / prime (essential∧ligandable∧broad)
  3  essential vs ligandable     essentiality vs ligandability; prime highlighted
  4  essential vs studiedness    essentiality vs studiedness; neglected-yet-essential = opportunity
  5  target-profile scorecard    top prime targets × 6 evidence axes (heatmap)
  6  prioritization funnel       all → essential → broad-selective → ligandable → prime

Reads output/results/<org>/<prefix>_{essentiality,ligandability,esmc600m_projection}.csv and
data/processed/<org>/bibliometric/<prefix>_popularity.csv.
Output: output/plots/07k_prioritization_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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
ESS_C = NPG[2]        # essential tier (green)
PRIME_C = "#E64B35"   # prime targets (NPG red)


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


def _label(row):
    g = row.get("gene")
    return g if isinstance(g, str) and g and str(g) != "nan" else row["uniprot_accession"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    orgname = E.ORG_DISPLAY[org]
    d = load(org)
    d["gene_lab"] = d.apply(_label, axis=1)
    has_xy = d["tsne_x"].notna() & d["tsne_y"].notna()

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: global essentiality score, reverse-cumulative + top targets labelled ----
    ax = axs.next()
    sc = d["essentiality_score"].to_numpy()
    xs = np.linspace(0, 1, 200)
    ys = np.array([(sc >= x).sum() for x in xs])
    ax.plot(xs, ys, color=ESS_C, lw=2)
    ax.fill_between(xs, ys, where=(xs >= 0.6), color=ESS_C, alpha=0.18)
    ax.set_yscale("log")
    ax.set_ylim(1, len(d))
    ax.axvline(0.6, color="#999", ls=":", lw=1)
    stylia.label(ax, xlabel="essentiality score", ylabel="proteins with score ≥ x",
                 title=f"Global essentiality score — {orgname}")
    # list the top targets as a clean ordered column in the shaded tail (no leader lines -> no clutter);
    # prime targets in red, others dark. They all sit at score≈1, so a vertical list reads cleanly.
    top = d.sort_values("essentiality_score", ascending=False).head(10).reset_index(drop=True)
    ys_lab = np.geomspace(len(d) * 0.5, 3, len(top))
    for i, r in top.iterrows():
        ax.text(0.635, ys_lab[i], r["gene_lab"], fontsize=SS, va="center", ha="left",
                color=PRIME_C if r["is_prime"] else "#333")

    # ---- panel 2: prioritization on the ESM-C map ----
    ax = axs.next()
    ax.scatter(d.loc[has_xy, "tsne_x"], d.loc[has_xy, "tsne_y"], s=4, alpha=0.28,
               color="#D8D8D6", linewidths=0, rasterized=True, label="all proteins")
    es = d[has_xy & d["is_essential"] & ~d["is_prime"]]
    ax.scatter(es["tsne_x"], es["tsne_y"], s=7, alpha=0.7, color=ESS_C, linewidths=0,
               rasterized=True, label="essential")
    pr = d[has_xy & d["is_prime"]]
    ax.scatter(pr["tsne_x"], pr["tsne_y"], s=26, color=PRIME_C, marker="*", linewidths=0,
               rasterized=True, label=f"prime ({int(d['is_prime'].sum())})")
    ax.set_xticks([]); ax.set_yticks([])
    stylia.label(ax, xlabel="ESM-C dim 1", ylabel="ESM-C dim 2",
                 title=f"Prioritization map — {orgname}")
    ax.legend(fontsize=SS, frameon=False, markerscale=1.4, loc="best")

    # ---- panel 3: essential vs ligandable ----
    ax = axs.next()
    bg = d[~d["is_prime"]]
    ax.scatter(bg["essentiality_score"], bg["ligandability_score"], s=4, alpha=0.22,
               color="#C9C9C7", linewidths=0, rasterized=True)
    ax.scatter(pr["essentiality_score"], pr["ligandability_score"], s=12, alpha=0.85,
               color=PRIME_C, linewidths=0, rasterized=True, label=f"prime ({len(pr)})")
    ax.axvline(0.6, color="#999", ls=":", lw=1); ax.axhline(0.6, color="#999", ls=":", lw=1)
    stylia.label(ax, xlabel="essentiality score", ylabel="ligandability score",
                 title=f"Essential & ligandable — {orgname}")
    ax.legend(fontsize=SS, frameon=False, loc="lower left")

    # ---- panel 4: essential vs studiedness (neglected-yet-essential = opportunity) ----
    ax = axs.next()
    sub = d[d["popularity_score"].notna()]
    ax.scatter(sub["essentiality_score"], sub["popularity_score"], s=4, alpha=0.22,
               color="#C9C9C7", linewidths=0, rasterized=True)
    prs = sub[sub["is_prime"]]
    ax.scatter(prs["essentiality_score"], prs["popularity_score"], s=12, alpha=0.85,
               color=PRIME_C, linewidths=0, rasterized=True)
    px = round(float(sub["popularity_score"].quantile(0.33)), 2)
    ax.axhline(px, color="#999", ls=":", lw=1); ax.axvline(0.6, color="#999", ls=":", lw=1)
    opp = prs[(prs["popularity_score"] <= px) & (prs["essentiality_score"] >= 0.6)] \
        .sort_values("essentiality_score", ascending=False)
    for _, r in opp.head(6).iterrows():
        ax.annotate(r["gene_lab"], (r["essentiality_score"], r["popularity_score"]),
                    fontsize=SS, color="#7a2417", xytext=(2, 2), textcoords="offset points")
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.0)
    stylia.label(ax, xlabel="essentiality score", ylabel="studiedness",
                 title=f"Neglected & essential — {orgname}")

    # ---- panel 5: target-profile scorecard (top prime targets × evidence axes) ----
    ax = axs.next()
    axes_cols = [("experimental", "evidence_experimental"), ("conservation", "entero_pct_essential"),
                 ("predictor", "evidence_predictor"), ("essentiality", "essentiality_score"),
                 ("ligandability", "ligandability_score"), ("studiedness", "popularity_score")]
    pr2 = d[d["is_prime"]].copy()
    pr2["rank"] = pr2["essentiality_score"].fillna(0) * pr2["ligandability_score"].fillna(0)
    pr2 = pr2.sort_values("rank", ascending=False).head(15)
    M = np.vstack([pd.to_numeric(pr2[c], errors="coerce").fillna(0).clip(0, 1).to_numpy()
                   for _, c in axes_cols]).T
    ax.imshow(M, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(axes_cols))); ax.set_xticklabels([a for a, _ in axes_cols], rotation=90, fontsize=SS)
    ax.set_yticks(range(len(pr2))); ax.set_yticklabels(pr2["gene_lab"], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="", title=f"Top prime-target scorecard — {orgname}")

    # ---- panel 6: prioritization funnel ----
    ax = axs.next()
    n_all = len(d)
    n_ess = int(d["is_essential"].sum())
    n_broad = int((d["is_essential"] & d["is_broad"]).sum())
    n_prime = int(d["is_prime"].sum())
    stages = [("all proteins", n_all, "#C9C9C7"),
              ("essential", n_ess, ESS_C),
              ("+ broad-spectrum selective", n_broad, NPG[1]),
              ("+ ligandable = prime", n_prime, PRIME_C)]
    counts = np.array([n for _, n, _ in stages], dtype=float)
    widths = np.sqrt(counts / counts[0])  # sqrt scaling keeps the small stages visible
    for i, (lab, n, c) in enumerate(stages):
        w = widths[i]
        ax.barh(i, w, height=0.74, left=(1 - w) / 2, color=c)
        # dark label centred; it overflows narrow bars onto the white ground and stays legible
        ax.text(0.5, i, f"{lab}: {n:,}", ha="center", va="center", fontsize=SS, color="#2B2333")
    ax.set_ylim(-0.6, len(stages) - 0.4); ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    stylia.label(ax, xlabel="", ylabel="", title=f"Prioritization funnel — {orgname}")

    out = REPO_ROOT / "output" / "plots" / f"07k_prioritization_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)} "
          f"(essential {n_ess}, broad-selective {n_broad}, prime {n_prime})")


if __name__ == "__main__":
    main()
