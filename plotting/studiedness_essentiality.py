"""How dark K. pneumoniae is, and why "unexplored" and "essential" are not independent.

    studiedness_essentiality.png   1  how much literature exists per protein, Kp against Ec
                                   2  why a zero means different things in the two organisms
                                   3  studiedness against essentiality, within each species

**THE FINDING: the two are correlated, and dropping the ambiguous zeros makes it STRONGER.**
Spearman between `studiedness_consensus` and `essentiality_consensus`:

    all proteins              Kp +0.332 (n=5,728)    Ec +0.440 (n=4,403)
    studiedness_evidence >=2  Kp +0.433 (n=3,767)    Ec +0.448 (n=4,374)
    studiedness_evidence ==3  Kp +0.639 (n=364)      Ec +0.456 (n=3,855)

On Kp it nearly doubles as the literature evidence improves, so this is not an artifact of the
34% of the proteome that scores zero. Controlling for protein length changes nothing either
(+0.332 -> +0.364 on Kp), so it is not a size effect. **This figure is the honest companion to
`novelty_vs_essentiality.png`**: the "unexplored AND essential" quadrant that the deck sells is
genuinely thinner than it looks, because people have already studied the essential genes.

**It is an association, and it is not causal in either direction.** People study essential genes
because they matter; and genes that are easy to study get called essential more often, because an
essentiality screen needs a gene you can disrupt and detect. Nothing here separates those.

**WHY THE CROSS-SPECIES PANELS DO NOT USE THE CONSENSUS SCORE.** `studiedness_consensus` is a mean
of WITHIN-SPECIES percentile ranks, so its distribution is uniform by construction: measured mean
**0.500** and median 0.506 in BOTH Kp and Ec. Two overlaid consensus histograms would be identical
and would say nothing. `src/consensus.py` states the rule -- *rank within a species, never across* --
and this is the first figure in the deck where obeying it visibly changes the design. So panels 1
and 2 use the RAW paper counts and the evidence tiers; only panel 3, which stays inside a species,
uses the consensus.

**A Kp zero and an Ec zero are not the same claim**, which is what panel 2 is for. Kp sits at
1,961 / 3,403 / 364 across evidence tiers 1/2/3 against Ec's 29 / 519 / 3,855 -- Ec is 87.6% at
tier 3, its own curated literature, where Kp is 6.4%. A Kp zero is usually `no_hit` or
`below_floor`, and the deliverable cannot tell those apart; only `load_transfer()` can.

Everything is read from the stages' own tables through `src/` loaders.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/studiedness_essentiality.py
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
from matplotlib.lines import Line2D  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import essentiality as E  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import studiedness as ST  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
SPECIES = ["kpneumoniae", "ecoli"]
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli"}

#: The column the project ranks studiedness on. NOT the consensus -- see the module docstring for
#: why the cross-species panels cannot use that.
PAPERS = "n_papers_uniprot_prokaryotic"

#: Evidence tiers, weakest first, so the stacked bar reads bottom-up as "how sure are we".
TIER_LABEL = {1: "1  no usable donor", 2: "2  one donor", 3: "3  its own literature"}
TIER_COLOR = {1: PAL.NPG[0], 2: PAL.NPG[1], 3: PAL.NPG[4]}


def plot_papers(ax, data: dict) -> None:
    """How much literature each protein has, as a violin per species.

    **A VIOLIN CANNOT REPRESENT A POINT MASS, and 34% of K. pneumoniae sits at exactly zero.** The
    kernel smears that spike into a smooth hump and, unclipped, puts density below zero where no
    protein can be. So the body is clipped at zero and the zero fraction is written on the panel as
    a number -- the one part of this distribution the shape genuinely cannot carry.

    Counts are log1p because the range is 0-58 with a median of 4: on a linear axis both proteomes
    collapse against the bottom."""
    for i, sp in enumerate(SPECIES):
        v = data[sp][PAPERS].to_numpy(float)
        parts = ax.violinplot([np.log1p(v)], positions=[i], widths=0.8, showextrema=False)
        for b in parts["bodies"]:
            # clip the kernel's tail at zero: it is an artifact of the smoothing, not data
            pv = b.get_paths()[0].vertices
            pv[:, 1] = np.clip(pv[:, 1], 0.0, None)
            b.set_alpha(1.0)
            b.set_facecolor(PAL.SPECIES_COLOR[sp])
            b.set_edgecolor("white")
            b.set_linewidth(1.0)

        med = float(np.median(v))
        ax.scatter([i], [np.log1p(med)], s=26, color=PAL.INK, zorder=5)
        ax.text(i, np.log1p(med) + 0.17, f"median {med:.0f}", ha="center", fontsize=SS * 0.8,
                color=PAL.INK)

    ticks = [0, 1, 2, 5, 10, 20, 50]
    ax.set_yticks(np.log1p(ticks))
    ax.set_yticklabels([str(t) for t in ticks], fontsize=SS)
    ax.set_ylim(-0.1, np.log1p(70))
    ax.set_xticks(range(len(SPECIES)))
    # The zero fraction rides in the tick label: as free-floating text under the violin it
    # collided with the species name, and it is the number the violin shape cannot carry.
    ax.set_xticklabels(
        [f"{LABELS[sp]}\n{(data[sp][PAPERS] == 0).mean() * 100:.1f}% at zero" for sp in SPECIES],
        fontsize=SS * 0.9,
    )
    ax.set_xlim(-0.65, len(SPECIES) - 0.35)
    stylia.label(ax, xlabel="", ylabel="Curated papers (prokaryotic donor)",
                 title="K. pneumoniae is a dark proteome")


def plot_tiers(ax, data: dict) -> None:
    """The evidence tiers, as a proportion bar per species.

    This panel licenses the one beside it: it says a zero in Kp and a zero in Ec are not the same
    claim. Proportions, with counts annotated, because the proteomes differ in size and the point
    is compositional."""
    x = np.arange(len(SPECIES))
    bottom = np.zeros(len(SPECIES))
    for tier in (1, 2, 3):
        vals, counts = [], []
        for sp in SPECIES:
            n = int((data[sp]["studiedness_evidence"] == tier).sum())
            counts.append(n)
            vals.append(n / len(data[sp]) * 100)
        ax.bar(x, vals, bottom=bottom, width=0.62, color=TIER_COLOR[tier],
               label=TIER_LABEL[tier], edgecolor="white", linewidth=0.8)
        for xi, v, b, n in zip(x, vals, bottom, counts):
            if v >= 6:
                ax.text(xi, b + v / 2, f"{n:,}", ha="center", va="center", fontsize=SS * 0.85,
                        color="white")
        bottom += np.array(vals)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[sp] for sp in SPECIES], fontsize=SS, style="italic")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=SS * 0.78, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.155),
              ncol=1, handletextpad=0.4, labelspacing=0.25)
    stylia.label(ax, xlabel="", ylabel="% of proteome", title="A zero is not the same claim")


def plot_relationship(ax, data: dict, n_bins: int) -> tuple[dict, dict]:
    """Studiedness against essentiality, within each species.

    Both axes are `<axis>_consensus`, which is a within-species percentile rank -- comparing the
    NUMBER across the two lines would be a mistake, but comparing the SHAPE is exactly what this
    panel is for. The faint scatter is the proteins; the line is the median per studiedness bin,
    which is what makes a monotone trend visible through 10,000 points."""
    curves, rhos = {}, {}
    for sp in SPECIES:
        d = data[sp]
        ax.scatter(d["studiedness_consensus"], d["essentiality_consensus"], s=2.0,
                   color=PAL.SPECIES_COLOR[sp], alpha=0.10, linewidths=0, rasterized=True)

        # duplicates="drop": Kp's 34% zero block is ONE tie, so 10 requested bins yield 8.
        d = d.assign(_bin=pd.qcut(d["studiedness_consensus"], n_bins, labels=False,
                                  duplicates="drop"))
        g = d.groupby("_bin").agg(x=("studiedness_consensus", "median"),
                                  y=("essentiality_consensus", "median"))
        curves[sp] = g
        ev = d[d["studiedness_evidence"] >= 2]
        rhos[sp] = (
            d["studiedness_consensus"].corr(d["essentiality_consensus"], method="spearman"),
            ev["studiedness_consensus"].corr(ev["essentiality_consensus"], method="spearman"),
            len(g),
        )
        ax.plot(g["x"], g["y"], color=PAL.SPECIES_COLOR[sp], lw=2.4, marker="o", markersize=5,
                markeredgecolor="white", markeredgewidth=0.6, zorder=4)

    handles = [
        Line2D([], [], color=PAL.SPECIES_COLOR[sp], lw=2.4, marker="o", markersize=5,
               label=f"{LABELS[sp]}   ρ {rhos[sp][0]:+.2f}   (ev≥2: {rhos[sp][1]:+.2f})")
        for sp in SPECIES
    ]
    ax.legend(handles=handles, fontsize=SS * 0.82, frameon=False, loc="upper left",
              handletextpad=0.5, labelspacing=0.35)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="Studiedness (within-species percentile)",
                 ylabel="Essentiality (within-species percentile)",
                 title="Better studied, more essential")
    return curves, rhos


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bins", type=int, default=10,
                    help="studiedness bins for the median curve (Kp yields fewer; see --help)")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  studiedness vs essentiality -- K. pneumoniae against E. coli")
    data = {}
    for sp in SPECIES:
        d = ST.load(sp).merge(E.load(sp), on="uniprot_ac")
        for c in ("studiedness_consensus", "essentiality_consensus", PAPERS):
            d[c] = pd.to_numeric(d[c], errors="coerce")
        pr = P.load(sp)[["uniprot_ac", "sequence"]].copy()
        pr["length"] = pr["sequence"].str.len()
        data[sp] = d.merge(pr[["uniprot_ac", "length"]], on="uniprot_ac")

    fig, axs = stylia.create_figure(1, 3, width_ratios=[1.25, 0.8, 1.25], width=1.0, height=0.42)
    plot_papers(axs.next(), data)
    plot_tiers(axs.next(), data)
    curves, rhos = plot_relationship(axs.next(), data, args.bins)
    out = OUT_DIR / "studiedness_essentiality.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  HOW MUCH LITERATURE")
    for sp in SPECIES:
        d = data[sp]
        say(f"    {LABELS[sp]:<16} n={len(d):>5,}  median {d[PAPERS].median():>3.0f}  "
            f"zero {100 * (d[PAPERS] == 0).mean():>5.1f}%  "
            f"own-median {pd.to_numeric(d['n_papers_uniprot_own'], errors='coerce').median():>3.0f}")

    say("\n  EVIDENCE TIERS")
    for sp in SPECIES:
        c = data[sp]["studiedness_evidence"].value_counts().sort_index()
        tot = len(data[sp])
        say(f"    {LABELS[sp]:<16} " + "  ".join(
            f"{t}:{int(c.get(t, 0)):>5,} ({100 * c.get(t, 0) / tot:>4.1f}%)" for t in (1, 2, 3)))

    say("\n  THE CORRELATION  rho(studiedness_consensus, essentiality_consensus)")
    say(f"    {'':<16} {'all':>18} {'evidence>=2':>18} {'evidence==3':>18}")
    for sp in SPECIES:
        d = data[sp]
        out_cells = []
        for mask, name in ((d.index == d.index, "all"), (d["studiedness_evidence"] >= 2, "ev2"),
                           (d["studiedness_evidence"] == 3, "ev3")):
            s = d[mask]
            r = s["studiedness_consensus"].corr(s["essentiality_consensus"], method="spearman")
            out_cells.append(f"{r:+.3f} (n={len(s):,})")
        say(f"    {LABELS[sp]:<16} " + " ".join(f"{c:>18}" for c in out_cells))
    say("    It STRENGTHENS as the evidence improves -- so it is not an artifact of the zero block.")

    say("\n  LENGTH CONTROL  (partial rho given protein length)")
    for sp in SPECIES:
        d = data[sp]
        r = d[["studiedness_consensus", "essentiality_consensus", "length"]].rank()
        X = np.c_[np.ones(len(r)), r["length"]]
        res = {}
        for c in ("studiedness_consensus", "essentiality_consensus"):
            beta = np.linalg.lstsq(X, r[c].to_numpy(), rcond=None)[0]
            res[c] = r[c].to_numpy() - X @ beta
        par = float(np.corrcoef(res["studiedness_consensus"], res["essentiality_consensus"])[0, 1])
        raw = d["studiedness_consensus"].corr(d["essentiality_consensus"], method="spearman")
        say(f"    {LABELS[sp]:<16} raw {raw:+.3f}   controlling for length {par:+.3f}")
    say("    Unchanged -- not a protein-size effect either.")

    say("\n  MEDIAN ESSENTIALITY BY STUDIEDNESS BIN")
    for sp in SPECIES:
        g = curves[sp]
        say(f"    {LABELS[sp]:<16} {len(g)} bins   " + " ".join(f"{v:.2f}" for v in g["y"]))
    say("    Kp yields fewer bins than requested: its 34% zero block is a single tie.")

    say("\n  CAVEATS")
    say("    - ASSOCIATION, NOT CAUSATION, in either direction. People study essential genes; and")
    say("      a gene that is easy to study is easier to call essential. This separates neither.")
    say("    - `studiedness_consensus` is a WITHIN-SPECIES percentile rank (mean 0.500 in both")
    say("      species by construction). Compare the shapes of the two curves, never the numbers.")
    say("    - A Kp zero is `no_hit` OR `below_floor` and the deliverable cannot tell them apart;")
    say("      only `load_transfer()` can. That is 1,961 Kp proteins.")
    say("    - FOR THE DECK: the 'unexplored and essential' quadrant in")
    say("      novelty_vs_essentiality.png is genuinely thinner than it looks. This is its")
    say("      honest companion, not a contradiction of it.")


if __name__ == "__main__":
    main()
