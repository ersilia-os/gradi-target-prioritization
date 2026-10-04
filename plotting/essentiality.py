"""Three independent predictors of essentiality, side by side, and deliberately NO verdict.

    essentiality.png   A  the three predictors' distributions on K. pneumoniae
                       B  how much they agree -- pairwise Spearman
                       C  how little their shortlists overlap -- shared members of each top 500

**The missing verdict column is a design decision, not an omission.** It was removed because it
MIXED UNITS -- a measured 1.0 pinned against a predicted 1.0, so every measured essential outranked
every prediction by construction -- and because on Kp it was a verbatim copy of `screens_ess_mean`
for all 5,728 rows. Which column to rank on is now an explicit choice.

**This deck ranks on `screens_ess_mean`**, and the reason is measured on Kp, the only anchor where
every predictor is honest: the screens-trained transfer models reach AUROC 0.89-0.96 on the three
measured Kp screens against Geptop's validated 0.59-0.81, and `geptop_ess` leaves 66.3% of the Kp
proteome tied at exactly 0 -- two-thirds of the anchor unranked, on an axis consumed by ranking.
Do NOT re-derive this from an E. coli comparison: `geptop_ess` is 53.2% self-derived there.

Panel C is the panel that justifies keeping all three: they are genuinely different opinions, not a
consensus. Contrast degradability's two activator columns at rho 0.89.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
`stylia.set_style("ersilia")` MUST precede `NamedColors()`, or NC.plum raises AttributeError.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/essentiality.py
"""

from __future__ import annotations

import argparse
import itertools
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

from src import essentiality as E  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

PREDICTORS = ["screens_ess_mean", "proteomelm_ess", "geptop_ess"]
PRED_SHORT = {"screens_ess_mean": "screens", "proteomelm_ess": "ProteomeLM", "geptop_ess": "Geptop"}
PRED_COLOR = {"screens_ess_mean": NC.plum, "proteomelm_ess": NC.blue, "geptop_ess": NC.orange}


def plot_distributions(ax, ess: pd.DataFrame, abc: str) -> None:
    """Three distributions. Geptop's spike at 0 is the point: 66.3% of Kp ties there, which is why
    it is not the column this deck ranks on."""
    for col in PREDICTORS:
        v = pd.to_numeric(ess[col], errors="coerce").dropna()
        ax.hist(v, bins=60, histtype="step", lw=2, color=PRED_COLOR[col],
                label=f"{PRED_SHORT[col]}  {(v == 0).mean() * 100:.0f}% tied at 0")
    ax.set_yscale("log")
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="predicted essentiality", ylabel="proteins (log)",
                 title="Three predictors, three shapes", abc=abc)


def plot_agreement(ax, ess: pd.DataFrame, abc: str) -> None:
    """Pairwise Spearman. Rank correlation, not Pearson, because the axis is consumed by ranking."""
    pairs = list(itertools.combinations(PREDICTORS, 2))
    vals = [ess[a].corr(ess[b], method="spearman") for a, b in pairs]
    y = np.arange(len(pairs))[::-1]
    ax.barh(y, vals, color=NC.plum, height=0.55)
    for yi, v in zip(y, vals):
        ax.text(v + 0.012, yi, f"{v:.2f}", va="center", fontsize=SS, color=NC.black)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{PRED_SHORT[a]} / {PRED_SHORT[b]}" for a, b in pairs], fontsize=SS)
    ax.set_xlim(0, 1)
    stylia.label(ax, xlabel="Spearman rho", ylabel="", title="They do not agree", abc=abc)


def plot_overlap(ax, ess: pd.DataFrame, top: int, abc: str) -> None:
    """Shared members of each pair's top-N. The operational question: would swapping the column
    change the shortlist? It would."""
    pairs = list(itertools.combinations(PREDICTORS, 2))
    tops = {c: set(ess.nlargest(top, c)["uniprot_ac"]) for c in PREDICTORS}
    vals = [len(tops[a] & tops[b]) for a, b in pairs]
    y = np.arange(len(pairs))[::-1]
    ax.barh(y, vals, color=NC.orange, height=0.55)
    for yi, v in zip(y, vals):
        ax.text(v + top * 0.015, yi, f"{v}/{top}", va="center", fontsize=SS, color=NC.black)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{PRED_SHORT[a]} / {PRED_SHORT[b]}" for a, b in pairs], fontsize=SS)
    ax.set_xlim(0, top * 1.15)
    stylia.label(ax, xlabel=f"shared members of each top {top}", ylabel="",
                 title="Different shortlists", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("--top", type=int, default=500)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  essentiality -- {LABELS[sp]}")
    ess = E.load(sp)
    for c in PREDICTORS:
        ess[c] = pd.to_numeric(ess[c], errors="coerce")

    fig, axs = stylia.create_figure(1, 3, width=1.0, height=0.38)
    plot_distributions(axs.next(), ess, abc="A")
    plot_agreement(axs.next(), ess, abc="B")
    plot_overlap(axs.next(), ess, args.top, abc="C")
    out = OUT_DIR / "essentiality.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  TIED AT ZERO (why the ranking column matters)")
    for c in PREDICTORS:
        say(f"    {PRED_SHORT[c]:<12} {(ess[c] == 0).sum():>6,}  ({(ess[c] == 0).mean() * 100:.1f}%)")
    say("\n  PAIRWISE AGREEMENT")
    for a, b in itertools.combinations(PREDICTORS, 2):
        rho = ess[a].corr(ess[b], method="spearman")
        shared = len(set(ess.nlargest(args.top, a)["uniprot_ac"])
                     & set(ess.nlargest(args.top, b)["uniprot_ac"]))
        say(f"    {PRED_SHORT[a]:<12} / {PRED_SHORT[b]:<12} rho {rho:>5.2f}"
            f"   top-{args.top} shared {shared}")
    say("\n  CAVEATS")
    say("    - No verdict column, on purpose: it mixed measured with predicted units.")
    say("    - This deck ranks on screens_ess_mean (AUROC 0.89-0.96 on measured Kp screens).")
    say("    - screens_ess_mean is the mean of NINE models that all read the SAME ProtT5")
    say("      embedding -- one correlated opinion, not nine votes.")
    say("    - Comparable WITHIN a species, never across one.")


if __name__ == "__main__":
    main()
