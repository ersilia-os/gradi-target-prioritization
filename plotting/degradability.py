"""Is this protein a substrate of activated, partnerless ClpP? The axis the whole project turns on.

    degradability.png   -  the Kp probability distribution, with the cut that actually applies
                        -  measured hit rate by compartment -- why "cytoplasmic only" is wrong
                        -  how far the S. aureus -> K. pneumoniae extrapolation is being pushed

**The compartment panel is the mechanistic finding worth the slide.** On the measured
S. aureus labels the hit rate is cytoplasm 0.180/0.275, membrane 0.031/0.105 and
EXTRACELLULAR 0.049/0.346. Membrane
proteins are protected -- inserted co-translationally, never a soluble cytoplasmic chain. Secreted
proteins are not: they transit the cytoplasm unfolded and ARE reachable. So the right filter is
"not membrane", not "cytoplasm only", and the obvious filter would discard the compartment with the
highest measured ONC212 rate.

**0.5 is the wrong threshold.** Probabilities top out at 0.868 on Kp, and the right cut MOVED when
the estimator became TabPFN-3.5: `BASE_RATE_THRESHOLD` is 0.328/0.313, re-derived as empirical
quantiles of the OOF distribution. Re-derive again on any estimator change.

**The two activator columns are not independent evidence**: they correlate at rho 0.89 while the
labels agree at only rho 0.52, so "both activators agree" is one opinion wearing two hats.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`, which returns the ersilia plum/orange/mint set.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/degradability.py
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

from src import degradability as D  # noqa: E402
from src import localization as LOC  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
ACT_COLOR = {"adep4": PAL.PRIMARY, "onc212": PAL.SECONDARY}


def plot_distribution(ax, deg: pd.DataFrame) -> None:
    """Both activators on Kp, with 0.5 drawn so the room can see why it is the wrong cut."""
    for act in D.ACTIVATORS:
        v = pd.to_numeric(deg[f"{act}_prob"], errors="coerce").dropna()
        ax.hist(v, bins=60, histtype="step", lw=2, color=ACT_COLOR[act],
                label=f"{act.upper()}  max {v.max():.3f}")
        ax.axvline(D.BASE_RATE_THRESHOLD[act], color=ACT_COLOR[act], lw=1.2, ls="--")
    ax.axvline(0.5, color=PAL.INK, lw=1.4, ls=":")
    ax.text(0.5, ax.get_ylim()[1] * 0.55, " 0.5 (wrong)", fontsize=SS, color=PAL.INK)
    ax.legend(fontsize=SS, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.98))
    stylia.label(ax, xlabel="Predicted probability", ylabel="Proteins",
                 title="Dashed = the cut that applies")


def plot_hit_rate(ax) -> None:
    """Measured hit rate by compartment, from the S. aureus labels -- the only measurements here.

    Counts are printed on the bars because two of these compartments are small, and a rate without
    its denominator is how an 8-member category becomes a headline."""
    lab = D.load_labels()
    loc = LOC.load("saureus")[["uniprot_ac", "localization"]]
    df = lab.merge(loc, on="uniprot_ac", how="inner")

    order = [c for c in LOC.LOC_CLASS_COLOR if (df["localization"] == c).sum() >= 20]
    x = np.arange(len(order))
    width = 0.38
    for i, act in enumerate(D.ACTIVATORS):
        rates, ns = [], []
        for c in order:
            sub = df[df["localization"] == c][act].dropna()
            rates.append(sub.mean() if len(sub) else np.nan)
            ns.append(len(sub))
        ax.bar(x + (i - 0.5) * width, rates, width=width, color=ACT_COLOR[act],
               label=act.upper())
        for xi, r, n in zip(x + (i - 0.5) * width, rates, ns):
            if not np.isnan(r):
                ax.text(xi, r + 0.008, f"{n}", ha="center", fontsize=SS * 0.8, color=PAL.INK)

    ax.set_xticks(x)
    ax.set_xticklabels([LOC.LOC_CLASS_ABBREV[c] for c in order], fontsize=SS)
    ax.legend(fontsize=SS, frameon=False, loc="upper left")
    stylia.label(ax, xlabel="Compartment (measured, S. aureus)", ylabel="Hit rate",
                 title="Secreted proteins ARE reachable")


def plot_extrapolation(ax, species: list[str]) -> None:
    """How similar each proteome is to the nearest labelled S. aureus protein.

    This is the price of the extrapolation, and it is the panel to read before quoting any Kp
    number: S. aureus rows are near their own labels by construction; Kp rows are not."""
    for sp in species:
        v = pd.to_numeric(D.load(sp)["nn_similarity"], errors="coerce").dropna()
        ax.hist(v, bins=60, histtype="step", lw=2,
                label=f"{LABELS[sp]}  median {v.median():.3f}")
    ax.legend(fontsize=SS, frameon=False, loc="upper left")
    stylia.label(ax, xlabel="Cosine to nearest labelled S. aureus protein", ylabel="Proteins",
                 title="What the model is extrapolating across")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  degradability -- {LABELS[sp]}")
    deg = D.load(sp)

    fig, axs = stylia.create_figure(1, 3, width=1.0, height=0.38)
    plot_distribution(axs.next(), deg)
    plot_hit_rate(axs.next())
    plot_extrapolation(axs.next(), ["kpneumoniae", "ecoli", "saureus"])
    out = OUT_DIR / "degradability.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  WHAT THE CUT SELECTS")
    for act in D.ACTIVATORS:
        v = pd.to_numeric(deg[f"{act}_prob"], errors="coerce")
        cut = D.BASE_RATE_THRESHOLD[act]
        say(f"    {act.upper():<8} max {v.max():.3f}   cut {cut}"
            f"   selects {(v >= cut).sum():,} ({(v >= cut).mean() * 100:.1f}%)"
            f"   at 0.5: {(v >= 0.5).sum():,}")
    rho = deg["adep4_prob"].corr(deg["onc212_prob"], method="spearman")
    say(f"\n    the two columns correlate at rho {rho:.3f} -- ONE opinion, not two")
    say(f"    cross-assay yardstick (NOT a ceiling): {D.CROSS_ASSAY_AUROC}")
    say("\n  CAVEATS")
    say("    - Kp and Ec rows are RANKING HYPOTHESES, not measurements. The premise that ESM-C")
    say("      cosine measures transferability is unvalidated without Gram-negative labels.")
    say("    - Rank on the probability; the cut reproduces the base rate on the LABELED set,")
    say("      not on a proteome.")


if __name__ == "__main__":
    main()
