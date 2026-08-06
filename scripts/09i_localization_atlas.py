"""Compartment atlas — where each subcellular compartment sits in protein space (docs §5.1).

Small multiples over the ESM-C protein-universe map (01b): one panel per compartment, that
compartment highlighted against the faded proteome. The question it answers is whether subcellular
localization is *legible in sequence space at all* — if compartments occupy distinct territory, the
axis is describing real structure rather than a label we bolted on; if they smear evenly across the
map, that is worth knowing and saying.

Panels are ordered inward → outward through the cell envelope, which is also the order Clp
accessibility falls in, so the slide reads as a journey from the reachable interior to the
unreachable surface.

  1  cytoplasm            2  inner membrane      3  periplasm
  4  outer membrane       5  extracellular       6  cell wall & surface

Each panel carries its protein count, its share of the proteome, and the mean Clp accessibility of
that compartment. A convex-hull-free density hint (the highlighted points themselves) is deliberate:
no smoothing, no interpolation, just the proteins.

Reads output/results/<org>/<prefix>_localization.csv + <prefix>_esmc600m_projection.csv.
Output: output/plots/09i_atlas_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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
from src import localization as LOC  # noqa: E402

SS = stylia.SLIDE_FONTSIZE_SMALL
BG = "#D8D8D6"

# The six real compartments, inward -> outward. `membrane` is provisional (09g resolves it) and
# `unknown` is empty at 100% coverage, so neither earns a panel.
PANEL_ORDER = ["cytoplasm", "inner_membrane", "periplasm",
               "outer_membrane", "extracellular", "cell_wall_surface"]

TITLE = {
    "cytoplasm": "Cytoplasm",
    "inner_membrane": "Inner membrane",
    "periplasm": "Periplasm",
    "outer_membrane": "Outer membrane",
    "extracellular": "Extracellular",
    "cell_wall_surface": "Cell wall & surface",
}


def load(org: str) -> pd.DataFrame:
    _, prefix = LOC.ORGANISMS[org]
    r = LOC.results_dir(org)
    d = pd.read_csv(r / f"{prefix}_localization.csv",
                    usecols=["uniprot_accession", "localization", "clp_accessibility",
                             "localization_evidence"])
    proj = r / f"{prefix}_esmc600m_projection.csv"
    if proj.exists():
        d = d.merge(pd.read_csv(proj)[["uniprot_accession", "tsne_x", "tsne_y"]],
                    on="uniprot_accession", how="left")
    for c in ("tsne_x", "tsne_y"):
        if c not in d:
            d[c] = np.nan
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]

    d = load(org)
    has_xy = d["tsne_x"].notna() & d["tsne_y"].notna()
    if not has_xy.any():
        raise SystemExit(f"[{org}] no ESM-C map coordinates — run 01b_esmc_projections.py first")
    xy = d[has_xy]
    n_total = len(xy)
    print(f"[{org}] {n_total}/{len(d)} proteins have map coordinates", flush=True)

    # Shared limits so all six panels are directly comparable.
    xlim = (xy["tsne_x"].min(), xy["tsne_x"].max())
    ylim = (xy["tsne_y"].min(), xy["tsne_y"].max())
    pad_x = 0.03 * (xlim[1] - xlim[0])
    pad_y = 0.03 * (ylim[1] - ylim[0])

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    for cls in PANEL_ORDER:
        ax = axs.next()
        sub = xy[xy["localization"] == cls]
        colour = LOC.LOC_CLASS_COLOR[cls]

        ax.scatter(xy["tsne_x"], xy["tsne_y"], s=3.5, alpha=0.30, color=BG,
                   linewidths=0, rasterized=True)
        if len(sub):
            ax.scatter(sub["tsne_x"], sub["tsne_y"], s=7, alpha=0.85, color=colour,
                       linewidths=0, rasterized=True)

        n = len(sub)
        pct = 100 * n / n_total if n_total else 0.0
        mean_clp = sub["clp_accessibility"].mean() if len(sub) else float("nan")
        n_exp = int((sub["localization_evidence"] == "experimental").sum())
        ax.text(0.03, 0.97,
                f"{n:,} proteins  ({pct:.1f}%)\n"
                f"mean Clp access. {mean_clp:.2f}\n"
                f"{n_exp:,} experimental",
                transform=ax.transAxes, fontsize=SS, va="top", ha="left", color="#2B2333",
                linespacing=1.45)

        ax.set_xlim(xlim[0] - pad_x, xlim[1] + pad_x)
        ax.set_ylim(ylim[0] - pad_y, ylim[1] + pad_y)
        ax.set_xticks([]); ax.set_yticks([])
        stylia.label(ax, xlabel="ESM-C tSNE-1", ylabel="ESM-C tSNE-2",
                     title=f"{TITLE[cls]} — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"09i_atlas_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")


if __name__ == "__main__":
    main()
