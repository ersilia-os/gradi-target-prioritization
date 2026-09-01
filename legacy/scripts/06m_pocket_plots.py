"""2x3 overview of structure-based druggability / binding-site detection on AlphaFold models (06e).

The user's favourite ligandability route: instead of needing experimental data, detect candidate
binding pockets directly on the AlphaFold model. Two orthogonal detectors are run (06e_pockets.py):
fpocket (geometric / alpha-sphere druggability) and P2Rank (ML ligand-site predictor, AlphaFold-
tuned). Each pocket's score is pLDDT-weighted — a high score in a low-confidence region is not
trusted — and combined into `pocket_consensus_score`. A consensus >= 0.50 is a confident druggable
site (matches POCKET_STRONG in 06g_ligandability_merge.py).

Same house style as 06i/06j/06k/06l: stylia "slide" format, NPG palette, white/dark in-bar labels.
Every panel is specific to the single `--organism` of the slide.

  1  AlphaFold model confidence    af_mean_plddt distribution (are the models good enough?)
  2  fpocket × P2Rank agreement    geometric vs ML detector, coloured by pLDDT
  3  why pLDDT-weight pockets       druggability vs best-pocket pLDDT (spurious low-pLDDT pockets)
  4  druggability funnel            AF model → ≥1 pocket → confident druggable site (consensus≥0.5)
  5  top druggable targets          best proteins by pocket_consensus_score
  6  consensus spread               pocket_consensus_score distribution across the proteome

Reads output/results/<org>/<prefix>_pockets.csv (06e) and _alphafold_structure.csv (04a).
Output: output/plots/06m_pocket_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

REPO_ROOT = L.REPO_ROOT
TOP_N = 8
NPG = stylia.CategoricalPalette("npg").colors
ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}                 # NPG per organism (matches 06i/06j/06k/06l)
DRUGGABLE = 0.50  # pocket_consensus_score >= this == confident druggable site (POCKET_STRONG in 06g)
PLDDT_CONF = 70.0  # AlphaFold "confident" threshold
PLDDT_HIGH = 90.0  # AlphaFold "very high" threshold
HIT = "#00A087"    # NPG green — passes consensus
MISS = "#C9C9C7"   # grey — fails consensus
ORGANISMS = [("kpneumoniae", "kp", "K. pneumoniae"), ("ecoli", "ec", "E. coli K-12")]


def load(prefix: str, org: str) -> pd.DataFrame:
    d = pd.read_csv(L.results_dir(org) / f"{prefix}_pockets.csv")
    af = pd.read_csv(L.results_dir(org) / f"{prefix}_alphafold_structure.csv",
                     usecols=["uniprot_accession", "af_available"])
    d = d.merge(af, on="uniprot_accession", how="left")
    d["af_available"] = d["af_available"].fillna(False).astype(bool)
    d["consensus"] = d["pocket_consensus_score"].fillna(0.0)
    d["is_druggable"] = d["consensus"] >= DRUGGABLE
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    orgname = {p: n for _, p, n in ORGANISMS}[prefix]
    d = load(prefix, org)
    genes = L.load_genes(org)
    col = ORG_COLOR.get(prefix, NPG[0])
    n_total = len(L.load_accessions(org))

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)  # 16:9 slide, 2x3 panels

    # ---- panel 1: AlphaFold model confidence — af_mean_plddt distribution ----
    ax = axs.next()
    plddt = d["af_mean_plddt"].dropna()
    ax.hist(plddt, bins=range(40, 101, 5), histtype="stepfilled", alpha=0.85, color=col)
    for thr in (PLDDT_CONF, PLDDT_HIGH):
        ax.axvline(thr, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Mean pLDDT"); ax.set_ylabel("Number of proteins")
    ax.set_title(f"AlphaFold model confidence — {orgname}")
    # left side of the axis is empty (mass piles at high pLDDT) — put the read-out there
    pct70 = 100.0 * (plddt >= PLDDT_CONF).mean() if len(plddt) else 0.0
    pct90 = 100.0 * (plddt >= PLDDT_HIGH).mean() if len(plddt) else 0.0
    ax.text(0.03, 0.96,
            f"median {plddt.median():.0f}\n{pct70:.0f}% ≥70 (confident)\n{pct90:.0f}% ≥90 (very high)",
            transform=ax.transAxes, ha="left", va="top", color="#555555", fontsize="small")

    # ---- panel 2: two detectors → consensus (per-method druggable-site counts) ----
    # fpocket (geometric) and P2Rank (ML) are only weakly correlated, so a raw scatter is an
    # uninformative blob; the useful read is how many sites each flags and how the pLDDT-weighted
    # consensus integrates them.
    ax = axs.next()
    fp = d["fpocket_max_drug_score"].fillna(0) >= DRUGGABLE
    p2 = d["p2rank_top_prob"].fillna(0) >= DRUGGABLE
    cats = ["fpocket\n(geometric)", "P2Rank\n(ML)", "Both\nagree"]
    vals = [int(fp.sum()), int(p2.sum()), int((fp & p2).sum())]
    bars = ax.bar(range(len(cats)), vals, color=["#3C5488", "#8491B4", HIT])
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats)
    ax.set_ylabel("Proteins with a druggable site (≥0.5)"); ax.set_xlabel("")
    ax.margins(y=0.16)
    ax.set_title(f"Two pocket detectors — {orgname}")

    # ---- panel 3: why we pLDDT-weight — raw druggability vs weighted consensus, by pocket pLDDT ----
    # binned so the effect is legible: a confident-looking pocket sitting in a low-pLDDT (unreliable)
    # region has its score pulled down by the weighting; the gap closes as pocket pLDDT rises.
    ax = axs.next()
    sub = d[d["fpocket_best_pocket_plddt"].notna()].copy()
    edges = [0, 50, 70, 90, 100]
    labs = ["<50", "50–70", "70–90", "90+"]
    sub["pb"] = pd.cut(sub["fpocket_best_pocket_plddt"], bins=edges, labels=labs, include_lowest=True)
    g = (sub.groupby("pb", observed=False)
            .agg(raw=("fpocket_max_drug_score", "mean"), cons=("consensus", "mean"),
                 n=("consensus", "size"))
            .reindex(labs))
    xi = np.arange(len(labs)); w = 0.38
    b1 = ax.bar(xi - w / 2, g["raw"].to_numpy(), w, color=MISS, label="raw fpocket druggability")
    b2 = ax.bar(xi + w / 2, g["cons"].to_numpy(), w, color=HIT, label="pLDDT-weighted consensus")
    ax.set_xticks(xi)
    ax.set_xticklabels([f"{lab}\n(n={int(v)})" for lab, v in zip(labs, g["n"].fillna(0))])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Best-pocket mean pLDDT"); ax.set_ylabel("Mean score (0–1)")
    ax.set_title(f"Why we pLDDT-weight — {orgname}")
    ax.legend(loc="upper left", fontsize="small")

    # ---- panel 4: druggability funnel ----
    ax = axs.next()
    labels = ["AlphaFold\nmodel", "≥1 detected\npocket", f"Druggable site\n(consensus≥{DRUGGABLE:g})"]
    has_pocket = (d["fpocket_n_pockets"].fillna(0) > 0) | (d["p2rank_n_pockets"].fillna(0) > 0)
    counts = [int(d["af_available"].sum()), int(has_pocket.sum()), int(d["is_druggable"].sum())]
    bars = ax.bar(range(len(labels)), counts, color=[NPG[5], MISS, HIT])
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
    ax.set_ylabel("Number of proteins"); ax.set_xlabel("")
    ax.margins(y=0.18)
    ax.set_title(f"Druggability funnel — {orgname}")
    ax.text(0.5, -0.22, f"of {n_total:,} proteins in the proteome",
            transform=ax.transAxes, ha="center", va="top", color="#777777", fontsize="small")

    # ---- panel 5: top druggable targets by consensus score ----
    ax = axs.next()
    top = d.sort_values("consensus", ascending=False).head(TOP_N).iloc[::-1]
    ys = range(len(top))
    bars = ax.barh(list(ys), top["consensus"].to_numpy(), color=col)
    ax.set_yticks(list(ys)); ax.set_yticklabels([genes.get(a) or a for a in top["uniprot_accession"]])
    ax.bar_label(bars, labels=[f"{v:.2f}" for v in top["consensus"]], label_type="edge", padding=2)
    ax.set_xlim(0, 1.08)
    ax.set_ylabel(""); ax.set_xlabel("pocket consensus score")
    ax.set_title(f"Top druggable targets — {orgname}")

    # ---- panel 6: consensus druggability across the proteome ----
    ax = axs.next()
    cons = d["consensus"]
    bins = np.linspace(0, 1, 21)
    n_arr, edges, patches = ax.hist(cons, bins=bins, color=MISS)
    for patch, left in zip(patches, edges[:-1]):  # shade the druggable region
        if left >= DRUGGABLE:
            patch.set_facecolor(HIT)
    ax.axvline(DRUGGABLE, color="#555555", linestyle="--", linewidth=0.8)
    n_drug = int(d["is_druggable"].sum())
    ax.set_xlabel("pocket consensus score"); ax.set_ylabel("Number of proteins")
    ax.set_title(f"Druggability across proteome — {orgname}")
    ax.legend(handles=[Patch(color=HIT, label=f"druggable (≥{DRUGGABLE:g}): {n_drug:,}")],
              loc="upper right", fontsize="small")

    out = REPO_ROOT / "output" / "plots" / f"06m_pocket_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
