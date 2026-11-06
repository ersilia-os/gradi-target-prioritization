"""Slide figure for the focused 3-way orthology analysis (03c outputs).

One combined stylia slide (npg palette) with four panels, for the target-prioritization story
(broad-spectrum + human-selective = prime targets):
  (a) UpSet of orthogroup membership across kp / ecoli / human
  (b) per-protein selectivity categories (kp vs ecoli)
  (c) RBH %identity distributions per pair (orthology depth)
  (d) proteome composition by selectivity class (fractions)

Reads the 03c tables in data/processed/other/orthology/. Run with the `gradi` env.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import stylia

REPO_ROOT = Path(__file__).resolve().parents[1]
ORTHO_DIR = REPO_ROOT / "data" / "processed" / "other" / "orthology"
OUT_PATH = REPO_ROOT / "output" / "plots" / "03d_orthology_focused.png"

NPG = stylia.CategoricalPalette("npg").colors
KP, EC, HS = NPG[0], NPG[4], NPG[-1]  # red / teal / grey(anti-target)
ORG_COLOR = {"kpneumoniae": KP, "ecoli": EC}
SETS = [("kp", "K. pneumoniae", KP), ("ec", "E. coli", EC), ("hs", "human", HS)]
SEL_ORDER = [
    "broad_selective",
    "narrow_selective",
    "broad_human_homolog",
    "narrow_human_homolog",
]
SEL_COLOR = {
    "broad_selective": NPG[3],
    "narrow_selective": NPG[2],
    "broad_human_homolog": NPG[5],
    "narrow_human_homolog": NPG[7],
}


def panel_upset(ax, og: pd.DataFrame) -> None:
    """Hand-rolled UpSet: intersection-size bars (top) + kp/ec/human dot-matrix (bottom)."""
    counts = og["region"].value_counts().to_dict()
    inter = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)  # (region, size)
    sizes = [s for _, s in inter]
    members = [set(region.split("+")) for region, _ in inter]
    x = np.arange(len(inter))

    ax.axis("off")
    bars_ax = ax.inset_axes([0.18, 0.34, 0.82, 0.66])
    mat_ax = ax.inset_axes([0.18, 0.0, 0.82, 0.30])

    bars_ax.bar(x, sizes, color=NPG[5], width=0.7)
    for xi, s in zip(x, sizes):
        bars_ax.text(
            xi,
            s,
            f"{s}",
            ha="center",
            va="bottom",
            fontsize=stylia.SLIDE_FONTSIZE_SMALL,
        )
    bars_ax.set_xlim(-0.6, len(inter) - 0.4)
    bars_ax.set_ylim(0, max(sizes) * 1.18)
    bars_ax.set_xticks([])
    bars_ax.set_ylabel("orthogroups")
    for sp in ("top", "right"):
        bars_ax.spines[sp].set_visible(False)

    for row, (key, label, color) in enumerate(SETS):
        y = len(SETS) - 1 - row
        for xi, mem in zip(x, members):
            on = key in mem
            mat_ax.scatter(xi, y, s=70, color=color if on else "#DDDDDD", zorder=3)
    # connect filled dots within each column
    for xi, mem in zip(x, members):
        ys = [len(SETS) - 1 - r for r, (k, _, _) in enumerate(SETS) if k in mem]
        if len(ys) > 1:
            mat_ax.plot([xi, xi], [min(ys), max(ys)], color="#888888", lw=1.2, zorder=2)
    mat_ax.set_xlim(-0.6, len(inter) - 0.4)
    mat_ax.set_ylim(-0.5, len(SETS) - 0.5)
    mat_ax.set_xticks([])
    mat_ax.set_yticks(range(len(SETS)))
    mat_ax.set_yticklabels(
        [lab for _, lab, _ in SETS][::-1], fontsize=stylia.SLIDE_FONTSIZE_SMALL
    )
    for sp in mat_ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Orthogroup membership (UpSet)")


def panel_selectivity(ax, cats: pd.DataFrame) -> None:
    ct = (
        cats.groupby(["organism", "selectivity"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=SEL_ORDER, fill_value=0)
    )
    x = np.arange(len(SEL_ORDER))
    w = 0.38
    for i, org in enumerate(["kpneumoniae", "ecoli"]):
        if org in ct.index:
            ax.bar(
                x + (i - 0.5) * w,
                ct.loc[org].values,
                width=w,
                color=ORG_COLOR[org],
                label="K. pneumoniae" if org == "kpneumoniae" else "E. coli",
            )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [
            "broad\nselective",
            "narrow\nselective",
            "broad\nhuman-hom.",
            "narrow\nhuman-hom.",
        ],
        fontsize=stylia.SLIDE_FONTSIZE_SMALL,
    )
    stylia.label(ax, xlabel="", ylabel="proteins", title="Selectivity categories")
    ax.legend(fontsize=stylia.SLIDE_FONTSIZE_SMALL, frameon=False)


def panel_identity(ax, rbh: pd.DataFrame) -> None:
    pair_color = {"kp-ec": NPG[3], "kp-human": KP, "ec-human": EC}
    for pair, color in pair_color.items():
        v = rbh.loc[rbh.pair == pair, "pident"].to_numpy()
        if len(v):
            ax.hist(
                v,
                bins=40,
                range=(0, 100),
                histtype="step",
                linewidth=1.8,
                color=color,
                label=f"{pair} (med {np.median(v):.0f}%)",
            )
    stylia.label(
        ax, xlabel="RBH % identity", ylabel="ortholog pairs", title="Orthology depth"
    )
    ax.legend(fontsize=stylia.SLIDE_FONTSIZE_SMALL, frameon=False, loc="upper left")


def panel_composition(ax, cats: pd.DataFrame) -> None:
    frac = (
        cats.groupby(["organism", "selectivity"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=SEL_ORDER, fill_value=0)
    )
    frac = frac.div(frac.sum(axis=1), axis=0)
    orgs = [o for o in ["kpneumoniae", "ecoli"] if o in frac.index]
    x = np.arange(len(orgs))
    bottom = np.zeros(len(orgs))
    for sel in SEL_ORDER:
        vals = frac.loc[orgs, sel].values
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=SEL_COLOR[sel],
            width=0.6,
            label=sel.replace("_", " "),
        )
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["K. pneumoniae", "E. coli"], fontsize=stylia.SLIDE_FONTSIZE_SMALL
    )
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel="fraction of proteome", title="Proteome composition")
    # park the legend in the empty gap between the two proteome bars
    ax.legend(
        fontsize=stylia.SLIDE_FONTSIZE_SMALL,
        frameon=False,
        loc="center",
        handlelength=1.2,
    )


def main() -> None:
    og = pd.read_csv(ORTHO_DIR / "three_way_orthogroups.tsv", sep="\t")
    cats = pd.read_csv(ORTHO_DIR / "three_way_protein_categories.tsv", sep="\t")
    rbh = pd.read_csv(ORTHO_DIR / "three_way_rbh_identity.tsv", sep="\t")

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 2, width=1.0, height=0.66)
    panel_upset(axs.next(), og)
    panel_selectivity(axs.next(), cats)
    panel_identity(axs.next(), rbh)
    panel_composition(axs.next(), cats)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(OUT_PATH))
    print(f"Wrote {OUT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
