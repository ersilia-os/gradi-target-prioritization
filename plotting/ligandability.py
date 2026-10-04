"""Has anyone ever made a molecule against this protein -- and did it work? Two questions, not one.

    ligandability.png   -  effort against potency: the quadrant that separates the two questions
                        -  how little of the proteome has any potent ligand at all
                        -  where the evidence comes from -- this protein, or a bacterial homolog

**A single ligand count conflates two different states, and the distinction is the whole axis.**
`n_ligands_*` is POTENT (pChEMBL >= 6); `n_assayed_*` is "has anyone looked". A 0 against 158
assayed compounds is a MEASURED DISCOURAGEMENT; a 0 against 0 assayed is an OPEN QUESTION. Kp
`pyrH` is the case to remember: 158 compounds assayed against a 98.3%-identical target, not one
potent -- previously indistinguishable from "nobody looked".

`n_assayed_*` is NA, never 0, where no effort extract exists, because a 0 there would claim nobody
ever assayed the protein -- the opposite evidence from "we do not know".

**This axis does NOT say "has an antibiotic".** `pchembl_value` exists only for `=` relations on
IC50/EC50/Ki/Kd/Potency in nM, so MIC and %-inhibition are absent by construction. That is why the
ribosome is largely missing here while leading the degradability ranking -- the two are not in
conflict, they are measuring different things.

~2% of each proteome has a potent ligand, and that must not be forced upward: it is a fact about
how little of the bacterial proteome anyone has screened, not about druggability. The transfer
bands (95/60/40%) control COVERAGE, not reliability -- they could not be calibrated, because
P(potent | neighbour potent) is FLAT from 25% to 100% identity.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`, which returns the ersilia plum/orange/mint set.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/ligandability.py
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

from src import ligandability as LG  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
SPECIES = ["kpneumoniae", "ecoli", "saureus"]


def plot_quadrant(ax, lg: pd.DataFrame) -> None:
    """Effort (x) against potency (y), both log1p. The bottom-right corner is the one the project
    built this axis to make visible: screened hard, nothing potent."""
    x = np.log1p(pd.to_numeric(lg["n_assayed_bacterial"], errors="coerce").fillna(0))
    y = np.log1p(pd.to_numeric(lg["n_ligands_bacterial"], errors="coerce").fillna(0))
    jitter = np.random.default_rng(0).normal(0, 0.045, len(x))
    ax.scatter(x + jitter, y + jitter, s=5, color=PAL.PRIMARY, alpha=0.45, linewidths=0,
               rasterized=True)

    discouraged = int(((pd.to_numeric(lg["n_assayed_bacterial"], errors="coerce").fillna(0) >= 10)
                       & (pd.to_numeric(lg["n_ligands_bacterial"], errors="coerce").fillna(0) == 0)).sum())
    ax.text(0.97, 0.13, f"screened, nothing potent: {discouraged:,}", transform=ax.transAxes,
            fontsize=SS, color=PAL.INK, ha="right")
    ticks = [0, 1, 10, 100, 1000]
    ax.set_xticks(np.log1p(ticks))
    ax.set_xticklabels([str(t) for t in ticks], fontsize=SS)
    ax.set_yticks(np.log1p(ticks))
    ax.set_yticklabels([str(t) for t in ticks], fontsize=SS)
    stylia.label(ax, xlabel="Compounds assayed (bacterial)", ylabel="Potent compounds",
                 title="A 0 is not always a 0")


def plot_fraction(ax, data: dict) -> None:
    """How little of each proteome carries any potent ligand. ~2%, and that is the real number."""
    x = np.arange(len(SPECIES))
    vals, ns = [], []
    for sp in SPECIES:
        lg = data[sp]
        n = int((pd.to_numeric(lg["n_ligands_bacterial"], errors="coerce").fillna(0) > 0).sum())
        ns.append(n)
        vals.append(n / len(lg) * 100)
    ax.bar(x, vals, color=PAL.PRIMARY, width=0.55)
    for xi, v, n in zip(x, vals, ns):
        ax.text(xi, v + 0.05, f"{n}\n{v:.1f}%", ha="center", fontsize=SS, color=PAL.INK)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS, style="italic")
    ax.set_ylim(0, max(vals) * 1.35)
    stylia.label(ax, xlabel="", ylabel="% of proteome", title="Chemical precedent is rare")


def plot_scope(ax, data: dict) -> None:
    """Own vs bacterial-homolog vs human. Human is a LIABILITY column and is never summed into the
    bacterial count -- clpP carries 106 human compounds against 61 bacterial."""
    x = np.arange(len(SPECIES))
    width = 0.26
    series = [("this protein", "n_ligands_own", PAL.PRIMARY),
              ("bacterial homolog", "n_ligands_bacterial", PAL.SECONDARY),
              ("human (liability)", "n_ligands_human", PAL.TERTIARY)]
    for i, (name, col, color) in enumerate(series):
        vals = [int((pd.to_numeric(data[sp][col], errors="coerce").fillna(0) > 0).sum())
                for sp in SPECIES]
        ax.bar(x + (i - 1) * width, vals, width=width, color=color, label=name)
        for xi, v in zip(x + (i - 1) * width, vals):
            ax.text(xi, v + 2, f"{v}", ha="center", fontsize=SS * 0.8, color=PAL.INK)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS, style="italic")
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="", ylabel="Proteins with a potent compound",
                 title="Homology transfer is the whole game")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  ligandability -- chemical precedent")
    data = {sp: LG.load(sp) for sp in SPECIES}

    fig, axs = stylia.create_figure(1, 3, width=1.0, height=0.38)
    plot_quadrant(axs.next(), data[args.species])
    plot_fraction(axs.next(), data)
    plot_scope(axs.next(), data)
    out = OUT_DIR / "ligandability.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  THE TWO QUESTIONS")
    for sp in SPECIES:
        lg = data[sp]
        assayed = pd.to_numeric(lg["n_assayed_bacterial"], errors="coerce").fillna(0)
        potent = pd.to_numeric(lg["n_ligands_bacterial"], errors="coerce").fillna(0)
        say(f"    {LABELS[sp]:<16} potent {int((potent > 0).sum()):>4,}"
            f"   assayed {int((assayed > 0).sum()):>5,}"
            f"   screened-but-nothing {int(((assayed >= 10) & (potent == 0)).sum()):>4,}"
            f"   never looked {int((assayed == 0).sum()):>5,}")
    say("\n  CAVEATS")
    say("    - This axis does NOT say 'has an antibiotic': pchembl_value excludes MIC and")
    say("      %-inhibition, which is why the ribosome looks empty here.")
    say("    - The identity bands control COVERAGE, not reliability. P(potent | neighbour")
    say("      potent) is FLAT from 25% to 100% identity -- they could not be calibrated.")
    say("    - ~2% is a fact about screening effort, not about druggability. Do not force it up.")


if __name__ == "__main__":
    main()
