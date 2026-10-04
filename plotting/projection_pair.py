"""One map, two readings: where every K. pneumoniae protein sits, and where the essential ones are.

    projection_pair.png   A  the t-SNE of the proteome, coloured by predicted localization
                          B  the same map, with the top 250 most essential proteins picked out

**The two panels share one set of coordinates on purpose.** Nothing is re-projected between them,
so every point is in the identical position in both -- which is what lets the eye carry compartment
structure from A into B. Re-running the projection per panel would silently break that.

**The essential set really does concentrate, and it was measured before it was drawn**: the top 250
have a mean pairwise distance of 24.4 on this map against 50.9 for random draws of the same size --
less than half. They are COG J (translation, 90 of 250) and cytoplasmic (197 of 250, 79%, against a
60% proteome background). So panel B is showing a real cluster, not a highlight that would look
equally convincing applied to any subset.

**Essentiality is `screens_ess_mean`** -- the column this deck ranks on, because on Kp the
screens-trained transfer models reach AUROC 0.89-0.96 on the three measured screens while
`geptop_ess` leaves 66.3% of the proteome tied at exactly 0, i.e. two thirds of the anchor unranked.
Comparable WITHIN a species, never across one. It is also nine models over the SAME ProtT5
embedding -- one correlated opinion, not nine votes.

**250 is a display choice, not a threshold**, and the axis ships no essentiality cutoff; `--top`
moves it. The shortlist cascade in `plotting/filters.py` uses a proteome-wide 90th percentile
instead, which is a different and larger set.

Coordinates are per-species with NO shared frame: never measure a distance across species on them.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`, which returns the ersilia plum/orange/mint set.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/projection_pair.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import stylia  # noqa: E402

from plotting import palette as PAL  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from src import essentiality as E  # noqa: E402
from src import localization as LOC  # noqa: E402
from src import projections as PROJ  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: ersilia. The STYLE sets typography and grid; the COLOURS come from
# `plotting/palette.py` -- stylia's npg palette, not stylia's ersilia palette.
stylia.set_format("slide")
stylia.set_style("ersilia")

# PALETTE: npg, NOT the ersilia NamedColors palette (owner's instruction). `set_style` controls
# typography and grid; the COLOURS come from CategoricalPalette("npg") -- ggsci's Nature palette,
# 10 colours. This is the same palette the localization vocabulary is already pinned to in
# `src/localization.py` (#E64B35, #00A087, #3C5488 are npg), so panel A and panel B agree by
# construction rather than by coincidence.
NPG = stylia.CategoricalPalette("npg").colors
ACCENT = NPG[8]   # magenta -- deliberately NOT a colour panel A uses for a compartment
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

# The compartment order the legend reads in: outward from the cytoplasm. Colours are pinned in
# src/localization.py so this figure cannot drift from the localization axis's own plots.
LOC_ORDER = [
    "cytoplasm", "cytoplasmic_membrane", "periplasm", "outer_membrane",
    "cell_wall_surface", "extracellular",
]
LOC_PRETTY = {
    "cytoplasm": "cytoplasm",
    "cytoplasmic_membrane": "inner membrane",
    "periplasm": "periplasm",
    "outer_membrane": "outer membrane",
    "cell_wall_surface": "cell wall",
    "extracellular": "extracellular",
}

# The quiet background for panel B -- a neutral, not a data colour, so it is not drawn from the
# categorical palette. Light enough that 250 highlighted points carry the panel, dark enough that
# the proteome's shape stays legible underneath them.
BACKGROUND = "#E4E4E0"


def _bare(ax, xy: np.ndarray) -> None:
    """A map, not a plot: no ticks, no frame, equal aspect so the embedding is not sheared."""
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    pad = 0.04 * (xy.max(axis=0) - xy.min(axis=0))
    ax.set_xlim(xy[:, 0].min() - pad[0], xy[:, 0].max() + pad[0])
    ax.set_ylim(xy[:, 1].min() - pad[1], xy[:, 1].max() + pad[1])
    # "box" not "datalim": datalim re-derives the limits just set above and warns about it.
    ax.set_aspect("equal", adjustable="box")


def plot_localization(ax, df: pd.DataFrame, abc: str) -> None:
    """The proteome by compartment. Drawn largest class first so the small, interesting classes
    (periplasm, outer membrane) end up on top rather than buried under the cytoplasm."""
    xy = df[["tsne_x", "tsne_y"]].to_numpy()
    present = [c for c in LOC_ORDER if (df["localization"] == c).any()]
    by_size = sorted(present, key=lambda c: -(df["localization"] == c).sum())

    for c in by_size:
        m = (df["localization"] == c).to_numpy()
        ax.scatter(xy[m, 0], xy[m, 1], s=3.2, color=LOC.LOC_CLASS_COLOR[c],
                   linewidths=0, alpha=0.85, rasterized=True)

    handles = [
        Line2D([], [], marker="o", ls="", markersize=5, color=LOC.LOC_CLASS_COLOR[c],
               label=f"{LOC_PRETTY[c]}   {(df['localization'] == c).sum():,}")
        for c in present
    ]
    ax.legend(handles=handles, fontsize=SS * 0.9, frameon=False, loc="upper left",
              handletextpad=0.3, labelspacing=0.35, borderpad=0.1)
    _bare(ax, xy)
    stylia.label(ax, xlabel="", ylabel="", title="Where every protein goes", abc=abc)


def plot_essential(ax, df: pd.DataFrame, top: int, abc: str) -> None:
    """The same coordinates, with the most essential proteins picked out.

    The highlighted points get a thin white edge so that a dense cluster still reads as individual
    proteins rather than a solid blob -- at 250 points in one region it otherwise saturates."""
    xy = df[["tsne_x", "tsne_y"]].to_numpy()
    sel = df.nlargest(top, "screens_ess_mean")
    m = df["uniprot_ac"].isin(set(sel["uniprot_ac"])).to_numpy()

    ax.scatter(xy[~m, 0], xy[~m, 1], s=3.0, color=BACKGROUND, linewidths=0, rasterized=True)
    ax.scatter(xy[m, 0], xy[m, 1], s=16, color=ACCENT, linewidths=0.35,
               edgecolors="white", zorder=3)

    handles = [
        Line2D([], [], marker="o", ls="", markersize=4, color=BACKGROUND,
               label=f"proteome   {len(df):,}"),
        Line2D([], [], marker="o", ls="", markersize=6, color=ACCENT,
               markeredgecolor="white", markeredgewidth=0.4,
               label=f"top {top} essential"),
    ]
    ax.legend(handles=handles, fontsize=SS * 0.9, frameon=False, loc="upper left",
              handletextpad=0.3, labelspacing=0.35, borderpad=0.1)
    _bare(ax, xy)
    stylia.label(ax, xlabel="", ylabel="", title="Where the essential ones are", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("--top", type=int, default=250, help="how many essential proteins to pick out")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  projection pair -- {LABELS[sp]}")

    df = (
        PROJ.load(sp)
        .merge(LOC.load(sp)[["uniprot_ac", "localization"]], on="uniprot_ac")
        .merge(E.load(sp)[["uniprot_ac", "screens_ess_mean"]], on="uniprot_ac")
    )
    df["screens_ess_mean"] = pd.to_numeric(df["screens_ess_mean"], errors="coerce")
    say(f"  {len(df):,} proteins with coordinates, localization and essentiality")

    fig, axs = stylia.create_figure(1, 2, width=1.0, height=0.46)
    plot_localization(axs.next(), df, abc="A")
    plot_essential(axs.next(), df, args.top, abc="B")
    out = OUT_DIR / "projection_pair.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    sel = df.nlargest(args.top, "screens_ess_mean")
    xy_sel = sel[["tsne_x", "tsne_y"]].to_numpy()
    rng = np.random.default_rng(0)

    def mean_pairwise(a: np.ndarray) -> float:
        return float(np.mean(np.linalg.norm(a[:, None] - a[None], axis=-1)))

    rand = np.mean([
        mean_pairwise(df.sample(args.top, random_state=i)[["tsne_x", "tsne_y"]].to_numpy())
        for i in range(5)
    ])
    say(f"\n  DO THEY CLUSTER?  mean pairwise distance on the map")
    say(f"    top {args.top} essential   {mean_pairwise(xy_sel):.1f}")
    say(f"    random {args.top}          {rand:.1f}   (mean of 5 draws)")
    say("    -- the highlight is a real cluster, not a subset that would look convincing anyway")

    say(f"\n  COMPARTMENT OF THE TOP {args.top}   (vs the proteome)")
    bg = df["localization"].value_counts(normalize=True) * 100
    for c, n in sel["localization"].value_counts().items():
        say(f"    {LOC_PRETTY[c]:<18} {n:>4}  {n / len(sel) * 100:>5.1f}%"
            f"   proteome {bg.get(c, 0):>5.1f}%")

    say("\n  CAVEATS")
    say("    - Coordinates are per-species with NO shared frame. Never compute a cross-species")
    say("      distance on these columns; use the 1,152-d space.")
    say(f"    - {args.top} is a DISPLAY choice. The axis ships no essentiality cutoff, and the")
    say("      shortlist cascade uses a proteome-wide 90th percentile instead.")
    say("    - screens_ess_mean is nine models over the SAME ProtT5 embedding: one correlated")
    say("      opinion, not nine votes. Comparable within a species, never across.")


if __name__ == "__main__":
    main()
