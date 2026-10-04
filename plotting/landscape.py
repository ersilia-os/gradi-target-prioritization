"""The K. pneumoniae proteome as one picture: 5,728 proteins, embedded and annotated.

    landscape.png   A  the t-SNE of the proteome, coloured by predicted localization
                    B  the same map, coloured by COG functional group
                    C  how much of the proteome each COG group accounts for

The opening slide. It argues one thing: the unit of analysis is a COMPLETE reference proteome, not
a gene list -- every protein has a vector, and the axes that follow are columns on these same rows.

Coordinates are openTSNE over 1,152-d ESM-C, inherited from v1 and NOT re-derived. They are
PER-SPECIES with no shared frame: never measure a distance across species on these columns.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/landscape.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import stylia  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import function as FN  # noqa: E402
from src import localization as LOC  # noqa: E402
from src import projections as PROJ  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: ersilia. The STYLE sets typography and grid; the COLOURS come from
# `plotting/palette.py`, which is stylia's **npg** palette and not stylia's ersilia palette
# (owner's instruction). Do not reintroduce `stylia.NamedColors()` here.
stylia.set_format("slide")
stylia.set_style("ersilia")

SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

# COG letters -> the three top-level groups, collapsed for a legible legend. 'R'/'S' are the
# uninformative pair and are kept visible rather than dropped: 'poorly characterised' is 20-26% of
# a bacterial proteome and hiding it would overstate how much is understood.
COG_GROUP_SHORT = {
    "INFORMATION STORAGE AND PROCESSING": "information",
    "CELLULAR PROCESSES AND SIGNALING": "cellular processes",
    "METABOLISM": "metabolism",
    "POORLY CHARACTERIZED": "poorly characterised",
}
COG_GROUP_COLOR = PAL.COG_GROUP_COLOR


def _scatter(ax, xy, colors, s=1.6):
    ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=s, linewidths=0, rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])


def plot_localization(ax, df, abc: str) -> None:
    """The map by compartment. Colours are pinned in `src/localization.py` so figures cannot drift."""
    colors = [LOC.LOC_CLASS_COLOR.get(v, "#DDDDDD") for v in df["localization"]]
    _scatter(ax, df[["tsne_x", "tsne_y"]].to_numpy(), colors)
    order = [c for c in LOC.LOC_CLASS_COLOR if (df["localization"] == c).any()]
    handles = [
        matplotlib.lines.Line2D(
            [], [], marker="o", ls="", markersize=4, color=LOC.LOC_CLASS_COLOR[c],
            label=f"{LOC.LOC_CLASS_ABBREV[c]}  {(df['localization'] == c).sum():,}",
        )
        for c in order
    ]
    ax.legend(handles=handles, fontsize=SS * 0.85, frameon=False, loc="upper left",
              handletextpad=0.2, labelspacing=0.25)
    stylia.label(ax, xlabel="", ylabel="", title="Where the protein goes", abc=abc)


def plot_cog(ax, df, abc: str) -> None:
    """The same map by COG group. Deliberately the SAME coordinates -- the point is that function
    and compartment partition the one embedding differently."""
    colors = [COG_GROUP_COLOR[g] for g in df["cog_group_short"]]
    _scatter(ax, df[["tsne_x", "tsne_y"]].to_numpy(), colors)
    handles = [
        matplotlib.lines.Line2D([], [], marker="o", ls="", markersize=4, color=v, label=k)
        for k, v in COG_GROUP_COLOR.items()
        if (df["cog_group_short"] == k).any()
    ]
    ax.legend(handles=handles, fontsize=SS * 0.85, frameon=False, loc="upper left",
              handletextpad=0.2, labelspacing=0.25)
    stylia.label(ax, xlabel="", ylabel="", title="What the protein does", abc=abc)


def plot_cog_bars(ax, df, abc: str) -> None:
    """How much of the proteome each group accounts for, so the map's colours have a scale."""
    counts = df["cog_group_short"].value_counts()
    order = [k for k in COG_GROUP_COLOR if k in counts.index]
    vals = [counts[k] / len(df) * 100 for k in order]
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, color=[COG_GROUP_COLOR[k] for k in order], height=0.6)
    for yi, v in zip(y, vals):
        ax.text(v + 1, yi, f"{v:.0f}%", va="center", fontsize=SS, color=PAL.INK)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=SS)
    ax.set_xlim(0, max(vals) * 1.22)
    stylia.label(ax, xlabel="% of proteome", ylabel="", title="How much is understood", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  landscape -- {LABELS[sp]}")

    df = PROJ.load(sp).merge(LOC.load(sp)[["uniprot_ac", "localization"]], on="uniprot_ac")
    cog = FN.load_cog(sp)[["uniprot_ac", "cog_category", "cog_group"]]
    df = df.merge(cog, on="uniprot_ac", how="left")
    df["cog_group_short"] = (
        df["cog_group"].map(COG_GROUP_SHORT).fillna("not classified")
    )
    say(f"  {len(df):,} proteins with coordinates, localization and COG")

    fig, axs = stylia.create_figure(1, 3, width_ratios=[3, 3, 2.2], width=1.0, height=0.38)
    plot_localization(axs.next(), df, abc="A")
    plot_cog(axs.next(), df, abc="B")
    plot_cog_bars(axs.next(), df, abc="C")
    out = OUT_DIR / "landscape.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  LOCALIZATION")
    for c, n in df["localization"].value_counts().items():
        say(f"    {c:<22} {n:>6,}  {n / len(df) * 100:>5.1f}%")
    say("\n  COG GROUP")
    for c, n in df["cog_group_short"].value_counts().items():
        say(f"    {c:<22} {n:>6,}  {n / len(df) * 100:>5.1f}%")
    say("\n  CAVEATS")
    say("    - Coordinates are per-species with NO shared frame: never compute a cross-species")
    say("      distance on them. Use the 1,152-d space.")
    say("    - DeepLocPro always returns a call, so 100% localization coverage is a property of")
    say("      the method, not evidence. 12-15% of calls sit below 0.7 confidence.")


if __name__ == "__main__":
    main()
