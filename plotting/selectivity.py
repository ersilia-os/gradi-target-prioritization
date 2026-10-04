"""Would a degrader hit the patient too? Human orthology as the selectivity filter.

    selectivity.png   A  how much of each proteome has a human ortholog, by method
                      B  the identity distribution of the human calls
                      C  what widening the comparator panel did -- in BOTH directions

**Two methods are kept and never merged, and this figure is the argument for that.** Kp `clpP` is a
documented ortholog of human CLPP at 56.3% identity; after the panel expansion OrthoFinder calls it
`of=0` and only RBH still finds it. At >=50% identity RBH rescues 9 Kp proteins OrthoFinder misses.
The union holds where either single method would not.

**The error direction is why this is conservative.** Under-detecting human homology makes a target
look MORE selective than it is -- the expensive direction -- so the union, plus `--very-sensitive`,
plus identity and coverage as COLUMNS rather than filters.

**Widening the panel moved both directions at once**, which an earlier plan got wrong: bacterial
orthogroup membership went UP and human orthologs went DOWN (Kp 1,184 -> 951). That is refinement,
not breakage -- the panel adds 26 bacteria and zero eukaryotes, so more bacterial resolution lets
OrthoFinder separate orthologs from out-paralogs and withdraw marginal human calls. At >=60%
identity 13/13 Kp human relationships survive; the losses sit at 40-50%, where ortholog-vs-paralog
is genuinely ambiguous.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
`stylia.set_style("ersilia")` MUST precede `NamedColors()`, or NC.plum raises AttributeError.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/selectivity.py
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

from src import orthology as O  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
SPECIES = ["kpneumoniae", "ecoli", "saureus"]


def plot_counts(ax, dense: dict, abc: str) -> None:
    """Per species: OrthoFinder alone, RBH alone, and the union that the deliverable ships."""
    x = np.arange(len(SPECIES))
    width = 0.26
    series = [
        ("OrthoFinder", "n_orthologs_of_human", NC.blue),
        ("RBH", "n_orthologs_rbh_human", NC.orange),
        ("union (shipped)", "has_human_ortholog", NC.plum),
    ]
    for i, (name, col, color) in enumerate(series):
        vals = []
        for sp in SPECIES:
            d = dense[sp]
            v = (d[col].astype(bool) if col == "has_human_ortholog"
                 else pd.to_numeric(d[col], errors="coerce").fillna(0) > 0)
            vals.append(v.mean() * 100)
        ax.bar(x + (i - 1) * width, vals, width=width, color=color, label=name)
        for xi, v in zip(x + (i - 1) * width, vals):
            ax.text(xi, v + 0.3, f"{v:.0f}", ha="center", fontsize=SS * 0.8, color=NC.black)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS, style="italic")
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="", ylabel="% with a human ortholog",
                 title="Neither method alone", abc=abc)


def plot_identity(ax, dense: dict, abc: str) -> None:
    """Identity of the human calls. The 40-50% band is where the two methods disagree and where
    ortholog-vs-paralog is genuinely ambiguous -- so identity ships as a column, never a filter."""
    for sp, color in zip(SPECIES, [NC.plum, NC.orange, NC.blue]):
        d = dense[sp]
        v = pd.to_numeric(d.loc[d["has_human_ortholog"].astype(bool), "best_identity_human"],
                          errors="coerce").dropna()
        if v.empty:
            continue
        ax.hist(v, bins=40, histtype="step", lw=2, color=color,
                label=f"{LABELS[sp]}  median {v.median():.1f}%")
    ax.axvspan(40, 50, color=NC.gray, alpha=0.35)
    ax.text(45, ax.get_ylim()[1] * 0.93, "ambiguous", fontsize=SS, color=NC.black, ha="center")
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="% identity to the human ortholog", ylabel="proteins",
                 title="Identity is a column, not a filter", abc=abc)


def plot_selectivity_gain(ax, dense: dict, abc: str) -> None:
    """What the filter buys: the proteome before and after removing human-orthologous proteins."""
    x = np.arange(len(SPECIES))
    width = 0.38
    before = [len(dense[sp]) for sp in SPECIES]
    after = [int((~dense[sp]["has_human_ortholog"].astype(bool)).sum()) for sp in SPECIES]
    ax.bar(x - width / 2, before, width=width, color=NC.gray, label="proteome")
    ax.bar(x + width / 2, after, width=width, color=NC.plum, label="no human ortholog")
    for xi, b, a in zip(x, before, after):
        ax.text(xi + width / 2, a + 60, f"{a / b * 100:.0f}%", ha="center", fontsize=SS,
                color=NC.black)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS, style="italic")
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="", ylabel="proteins", title="What the filter keeps", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  selectivity -- human orthology")
    dense = {sp: O.load_dense(sp) for sp in SPECIES}

    fig, axs = stylia.create_figure(1, 3, width=1.0, height=0.38)
    plot_counts(axs.next(), dense, abc="A")
    plot_identity(axs.next(), dense, abc="B")
    plot_selectivity_gain(axs.next(), dense, abc="C")
    out = OUT_DIR / "selectivity.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  HUMAN ORTHOLOGS")
    for sp in SPECIES:
        d = dense[sp]
        of = (pd.to_numeric(d["n_orthologs_of_human"], errors="coerce").fillna(0) > 0).sum()
        rbh = (pd.to_numeric(d["n_orthologs_rbh_human"], errors="coerce").fillna(0) > 0).sum()
        uni = d["has_human_ortholog"].astype(bool).sum()
        say(f"    {LABELS[sp]:<16} OrthoFinder {of:>5,}   RBH {rbh:>5,}   union {uni:>5,}"
            f"   ({uni / len(d) * 100:.1f}%)   RBH-only {uni - of:>4,}")
    say("\n  CAVEATS")
    say("    - Under-detecting human homology makes a target look MORE selective than it is, so")
    say("      the UNION of both methods is the conservative choice. Never merge them.")
    say("    - A human ortholog is not automatically disqualifying -- ClpP itself is conserved.")
    say("      It is a liability to price, which is why identity ships as a column.")


if __name__ == "__main__":
    main()
