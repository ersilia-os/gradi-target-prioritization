"""Capstone synthesis: the ligandability landscape & target prioritization (docs §2, result).

Ties the whole assessment together on the ESM-C protein-universe map (01b): where in protein space
do druggable, human-selective targets live, and which are the prime opportunities? Combines the
06g composite (ligandability_score / tier / selectivity / per-track evidence), the ESM-C 2D
embedding, and the bibliometric "studiedness" score (05a) to surface neglected-but-druggable prime
targets — the actionable payoff of the pipeline.

Prime target := broad-spectrum (conserved across bacteria) + human-selective (no human ortholog) +
tractable tier. "Human-selective" = selectivity in {broad_selective, narrow_selective} (no human
homolog); "broad-spectrum" = broad_selective.

  1  druggability landscape       ESM-C map coloured by ligandability tier
  2  prime targets in space       ESM-C map, prime (broad-selective + tractable) highlighted
  3  prioritization funnel        proteome → tractable → human-selective → broad-spectrum (prime)
  4  evidence basis of prime      which evidence supports the prime set (binding / structure / pocket)
  5  neglected & druggable        studiedness vs ligandability; understudied + druggable = opportunity
  6  top prime targets            highest-ligandability prime targets, shaded by studiedness

Reads output/results/<org>/<prefix>_{ligandability,esmc600m_projection}.csv and
data/processed/<org>/bibliometric/<prefix>_popularity.csv. Output: output/plots/06o_landscape_<prefix>.png
(one slide per --organism). Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)  # stylia rmtree's this on import; ensure it exists
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

REPO_ROOT = L.REPO_ROOT
TOP_N = 10
NPG = stylia.CategoricalPalette("npg").colors
ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}
TIER_COLOR = {"intractable": "#D2D2D0", "partial": NPG[4], "tractable": NPG[2]}
DARK = "#E64B35"   # understudied highlight (opportunity)
ORGNAME = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli K-12"}
HUMAN_SELECTIVE = {"broad_selective", "narrow_selective"}  # no human ortholog


def load(org: str, prefix: str) -> pd.DataFrame:
    lg = pd.read_csv(L.results_dir(org) / f"{prefix}_ligandability.csv")
    em = pd.read_csv(L.results_dir(org) / f"{prefix}_esmc600m_projection.csv")
    pop = pd.read_csv(REPO_ROOT / "data" / "processed" / org / "bibliometric" / f"{prefix}_popularity.csv",
                      usecols=["uniprot_accession", "popularity_score", "popularity_tier"])
    d = lg.merge(em[["uniprot_accession", "tsne_x", "tsne_y"]], on="uniprot_accession", how="left")
    d = d.merge(pop, on="uniprot_accession", how="left")
    d["is_tractable"] = d["ligandability_tier"] == "tractable"
    d["is_hsel"] = d["selectivity"].isin(HUMAN_SELECTIVE)
    d["is_prime"] = d["is_tractable"] & (d["selectivity"] == "broad_selective")
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = L.ORGANISMS[org]
    orgname = ORGNAME[org]
    d = load(org, prefix)
    col = ORG_COLOR.get(prefix, NPG[0])
    has_xy = d["tsne_x"].notna() & d["tsne_y"].notna()

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)  # 16:9 slide, 2x3 panels

    # ---- panel 1: druggability landscape — ESM-C map coloured by tier ----
    ax = axs.next()
    for tier in ["intractable", "partial", "tractable"]:  # draw tractable on top
        s = d[has_xy & (d["ligandability_tier"] == tier)]
        ax.scatter(s["tsne_x"], s["tsne_y"], s=4, alpha=0.55, linewidths=0,
                   color=TIER_COLOR[tier], rasterized=True, label=tier)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("ESM-C tSNE-1"); ax.set_ylabel("ESM-C tSNE-2")
    ax.set_title(f"Druggability landscape — {orgname}")
    ax.legend(loc="upper right", markerscale=2.2, fontsize="small", title="tier")

    # ---- panel 2: prime targets in protein space ----
    ax = axs.next()
    bg = d[has_xy & ~d["is_prime"]]
    pr = d[has_xy & d["is_prime"]]
    ax.scatter(bg["tsne_x"], bg["tsne_y"], s=4, alpha=0.35, linewidths=0, color="#D8D8D6", rasterized=True)
    ax.scatter(pr["tsne_x"], pr["tsne_y"], s=7, alpha=0.8, linewidths=0, color=col, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("ESM-C tSNE-1"); ax.set_ylabel("ESM-C tSNE-2")
    ax.set_title(f"Prime targets in protein space — {orgname}")
    ax.legend(handles=[Patch(color=col, label=f"prime (n={len(pr)})"),
                       Patch(color="#D8D8D6", label="rest of proteome")],
              loc="upper right", fontsize="small")

    # ---- panel 3: prioritization funnel ----
    ax = axs.next()
    n_all = len(d)
    n_tract = int(d["is_tractable"].sum())
    n_tract_hsel = int((d["is_tractable"] & d["is_hsel"]).sum())
    n_prime = int(d["is_prime"].sum())
    labels = ["Proteome", "Tractable", "+ human-\nselective", "+ broad-\nspectrum"]
    vals = [n_all, n_tract, n_tract_hsel, n_prime]
    colors = ["#D2D2D0", TIER_COLOR["tractable"], NPG[0], col]
    bars = ax.bar(range(4), vals, color=colors)
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(4)); ax.set_xticklabels(labels); ax.set_xlabel("")
    ax.set_ylabel("Number of proteins"); ax.margins(y=0.16)
    ax.set_title(f"Prioritization funnel — {orgname}")
    ax.text(0.97, 0.95, f"prime shortlist: {n_prime}", transform=ax.transAxes,
            ha="right", va="top", color=col, fontsize="small", fontweight="bold")

    # ---- panel 4: evidence basis of the prime set ----
    ax = axs.next()
    pr = d[d["is_prime"]]
    exp_bind = (pr["has_hard_evidence"].fillna(False).astype(bool)) | (pr["evidence_binding"].fillna(0) > 0)
    struct = pr["evidence_structural"].fillna(0) > 0
    pocket = pr["evidence_pocket"].fillna(0) >= 0.5
    struct_or_pocket = struct | pocket
    struct_only = struct_or_pocket & ~exp_bind  # only findable via the structure-first route
    cats = ["Experimental\nbinding", "Structural\n(PDB/AlphaFill)", "Druggable\npocket", "Structure-only\n(no exp. binding)"]
    vals = [int(exp_bind.sum()), int(struct.sum()), int(pocket.sum()), int(struct_only.sum())]
    bars = ax.bar(range(4), vals, color=[NPG[0], NPG[3], NPG[2], col])
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(4)); ax.set_xticklabels(cats); ax.set_xlabel("")
    ax.set_ylabel(f"Prime targets (of {len(pr)})"); ax.margins(y=0.16)
    ax.set_title(f"What makes prime targets ligandable — {orgname}")

    # ---- panel 5: neglected & druggable — studiedness vs ligandability ----
    ax = axs.next()
    sub = d[d["popularity_score"].notna()]
    ax.scatter(sub["popularity_score"], sub["ligandability_score"], s=4, alpha=0.3,
               linewidths=0, color="#D8D8D6", rasterized=True)
    pr = sub[sub["is_prime"]]
    ax.scatter(pr["popularity_score"], pr["ligandability_score"], s=8, alpha=0.7,
               linewidths=0, color=col, rasterized=True)
    # "neglected & druggable" opportunity box: understudied (low studiedness) + above the prime
    # median ligandability — i.e. druggable, selective, yet uncharacterised.
    px = 0.33
    py = round(float(pr["ligandability_score"].median()), 2) if len(pr) else 0.4
    ax.add_patch(Rectangle((0, py), px, 1.05 - py, facecolor=DARK, alpha=0.08, zorder=0))
    ax.axvline(px, color="#999", ls="--", lw=0.7); ax.axhline(py, color="#999", ls="--", lw=0.7)
    opp = pr[(pr["popularity_score"] <= px) & (pr["ligandability_score"] >= py)] \
        .sort_values("ligandability_score", ascending=False)
    genes = L.load_genes(org)
    for _, r in opp.head(6).iterrows():
        lab = (r.get("gene") if isinstance(r.get("gene"), str) and r.get("gene") else None) \
            or genes.get(r["uniprot_accession"]) or r["uniprot_accession"]
        ax.annotate(lab, (r["popularity_score"], r["ligandability_score"]),
                    fontsize=7, color="#7a2417", xytext=(2, 2), textcoords="offset points")
    ax.set_xlim(0, 1); ax.set_ylim(0, max(0.95, float(sub["ligandability_score"].max()) * 1.02))
    ax.set_xlabel("Studiedness (bibliometric)"); ax.set_ylabel("Ligandability score")
    ax.set_title(f"Neglected & druggable — {orgname}")
    ax.text(px / 2, py + (ax.get_ylim()[1] - py) * 0.93, f"opportunity: {len(opp)} prime",
            ha="center", va="top", color="#7a2417", fontsize="small", fontweight="bold")

    # ---- panel 6: top prime targets, shaded by studiedness ----
    ax = axs.next()
    top = d[d["is_prime"]].sort_values("ligandability_score", ascending=False).head(TOP_N).iloc[::-1]
    ys = range(len(top))
    bar_cols = [DARK if t == "dark" else (NPG[4] if t == "studied" else "#B9B9B7")
                for t in top["popularity_tier"].fillna("dark")]
    bars = ax.barh(list(ys), top["ligandability_score"].to_numpy(), color=bar_cols)
    labels = [(g if isinstance(g, str) and g else a) for g, a in zip(top["gene"], top["uniprot_accession"])]
    ax.set_yticks(list(ys)); ax.set_yticklabels(labels)
    ax.bar_label(bars, labels=[f"{v:.2f}" for v in top["ligandability_score"]], label_type="edge", padding=2)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Ligandability score"); ax.set_ylabel("")
    ax.set_title(f"Top prime targets — {orgname}")
    ax.legend(handles=[Patch(color=DARK, label="dark"), Patch(color=NPG[4], label="studied"),
                       Patch(color="#B9B9B7", label="well-studied")],
              loc="lower right", fontsize="small", title="studiedness")

    out = REPO_ROOT / "output" / "plots" / f"06o_landscape_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
