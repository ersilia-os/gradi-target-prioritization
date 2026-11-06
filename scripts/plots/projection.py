"""Figures for the stage-01 2D projections: one map per species, plus a 3-panel comparison.

Plotting is kept out of `embeddings/projection.py` -- matplotlib and stylia are a presentation
concern, not a dimensionality-reduction one, and the maps are expensive enough that re-styling a
figure must not mean re-projecting a proteome.

Uncoloured by design. Each map is a single stylia hue faded by 2D point density, so dense family
cores read saturated and sparse outliers read pale. That is an aesthetic for an overplotted scatter,
not an encoding of any variable -- there is no legend because there is nothing to look up.

To colour by something real, join it on `uniprot_ac` and pass `c=`; the projection table carries
coordinates only:

    proj = projections.load("kpneumoniae")
    cog  = function.load_cog("kpneumoniae")
    df   = proj.merge(cog[["uniprot_ac", "cog_category"]], on="uniprot_ac", how="left")

Each panel keeps its own axis limits on purpose. The three species are independent embeddings with
no shared frame, so a common `xlim`/`ylim` would imply a comparability that does not exist.

Outputs:
    output/plots/embeddings/projection_<species>.png
    output/plots/embeddings/projection_all.png

Run with the `gradi` env:

    ~/miniconda3/envs/gradi/bin/python scripts/plots/projection.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# stylia rmtree's the matplotlib cache dir on import; recreate it first or the run can crash.
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)

import numpy as np  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import stylia  # noqa: E402

from src import projections as PROJ  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "embeddings"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
# Kp cobalt and Ec crimson are v1's hues, carried over so figures do not drift between versions.
HUES = {"kpneumoniae": "cobalt", "ecoli": "crimson", "saureus": "turquoise"}


def fade_colors(xy: np.ndarray, hue: str) -> np.ndarray:
    """A single stylia hue, faded by 2D point density: sparse = light, dense = saturated.

    The floor of 0.2 keeps isolated points visible rather than fading them into the page.
    """
    d = gaussian_kde(xy.T)(xy.T)
    d = (d - d.min()) / (d.max() - d.min() + 1e-12)
    return stylia.FadingColormap(hue).cmap(0.2 + 0.8 * d)


def draw(ax, species: str, abc: str | None = None) -> int:
    """Scatter one species' map onto `ax`. Returns the number of points drawn."""
    df = PROJ.load(species)
    xy = df[list(PROJ.COORD_COLS)].to_numpy(dtype=float)
    ax.scatter(
        xy[:, 0], xy[:, 1],
        c=fade_colors(xy, HUES[species]),
        linewidths=0, alpha=0.85, rasterized=True,
    )
    # `adjustable="datalim"` keeps every axes box the same size and pads the limits instead of
    # shrinking the box to the data. With the default ("box") each panel's box shrinks to its own
    # data aspect, which leaves the three titles sitting at three different heights.
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    # The axes carry no units worth labelling -- t-SNE coordinates are not a measurement.
    stylia.label(ax, xlabel="", ylabel="",
                 title=f"{LABELS[species]}  ({len(df):,})", abc=abc)
    return len(df)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(PROJ.SPECIES), choices=list(PROJ.SPECIES))
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # One map per species: the figure you actually look at to find families. Square data space,
    # hence width == height.
    for sp in args.species:
        fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
        n = draw(axs.next(), sp)
        out = OUT_DIR / f"projection_{sp}.png"
        stylia.save_figure(str(out))
        print(f"wrote {out.relative_to(REPO_ROOT)}  ({n:,} proteins)")

    # And a 3-panel sheet, for the deck. Independent axes -- see the module docstring.
    fig, axs = stylia.create_figure(1, 3)
    for abc, sp in zip("ABC", PROJ.SPECIES):
        draw(axs.next(), sp, abc=abc)
    out = OUT_DIR / "projection_all.png"
    stylia.save_figure(str(out))
    print(f"wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
