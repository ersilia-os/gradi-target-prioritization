"""Lightweight 2x2 plot of relative vs absolute ESM projections, for both organisms.

Reads only precomputed coordinates (no heavy compute) and draws a 2x2 map:

              relative (openTSNE)        absolute (ESM Atlas)
  K. pneumoniae   01b CSV                    01c CSV
  E. coli K-12    01b CSV                    01c CSV

Style: one NPG (Nature Publishing Group) colour per organism, alpha-blended so density
reads as colour intensity. No frame, no chart junk; square panels.

Inputs (per organism <org>/<prefix>):
  output/results/<org>/<prefix>_esmc600m_projection.csv   (uniprot_accession, tsne_x, tsne_y, family)
  output/results/<org>/<prefix>_esmatlas_coords.csv       (uniprot_accession, atlas_x, atlas_y, matched)

Output: output/plots/01d_esm_projections.png
Styling via stylia (ersilia-os/stylia). Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_PATH = REPO_ROOT / "output" / "plots" / "01d_esm_projections.png"

# Nature Publishing Group (npg) palette — one colour per organism
NPG = [
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148",
    "#B09C85",
]

# (organism key, file prefix, display name, npg colour)
ORGANISMS = [
    ("kpneumoniae", "kp", "K. pneumoniae HS11286", NPG[0]),  # red
    ("ecoli", "ec", "E. coli K-12 MG1655", NPG[1]),  # blue
]


def draw(ax, xv, yv, color, size, alpha, invert=False):
    """Single-colour, alpha-blended density scatter into a clean square panel."""
    ax.scatter(xv, yv, s=size, c=color, alpha=alpha, edgecolors="none", rasterized=True)
    cx, cy = (xv.min() + xv.max()) / 2, (yv.min() + yv.max()) / 2
    r = max(xv.max() - xv.min(), yv.max() - yv.min()) / 2 * 1.05
    ax.set_xlim(cx - r, cx + r)
    ax.set_ylim(cy + r, cy - r) if invert else ax.set_ylim(cy - r, cy + r)
    ax.set_box_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> None:
    fig, axs = stylia.create_figure(2, 2, width=1.35, height=1.3)

    for i, (org, prefix, name, color) in enumerate(ORGANISMS):
        rel = pd.read_csv(
            REPO_ROOT / "output" / "results" / org / f"{prefix}_esmc600m_projection.csv"
        )
        atlas = pd.read_csv(
            REPO_ROOT / "output" / "results" / org / f"{prefix}_esmatlas_coords.csv"
        )
        atlas = atlas[atlas["matched"]]

        # relative map is dense (one tight cloud) -> low alpha; atlas is sparse/spread -> higher alpha
        ax = axs.next()
        draw(ax, rel["tsne_x"], rel["tsne_y"], color, size=6, alpha=0.22)
        ax.set_ylabel(name)
        if i == 0:
            ax.set_title("Relative (openTSNE)")

        ax = axs.next()
        draw(
            ax,
            atlas["atlas_x"],
            atlas["atlas_y"],
            color,
            size=6,
            alpha=0.45,
            invert=True,
        )
        if i == 0:
            ax.set_title("Absolute (ESM Atlas)")

        print(f"{name}: relative {len(rel)} pts, absolute {len(atlas)} matched pts")

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(PLOT_PATH))
    print(f"\nWrote {PLOT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
