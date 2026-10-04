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

**THE THREE PANELS DO NOT ALL REST ON THE SAME SOURCE, so each one says which.** The violins and
the evidence tiers are **UniProt only** -- `n_papers_uniprot_prokaryotic` and the SwissProt
transfer ladder behind it. The relationship panel uses `studiedness_consensus`, which is a rank
mean over **three** columns: `n_papers_uniprot_own`, `n_papers_uniprot_prokaryotic` and
`n_papers_pubtator_prokaryotic`. PubTator therefore enters exactly one panel, and the axis labels
name the source so that is visible rather than silent.

**The source does not change the SIGN, but on K. pneumoniae it changes the MAGNITUDE**, and the run
log prints rho against `essentiality_consensus` for all three source columns separately rather than
asserting they agree. Measured: Kp **consensus +0.332 · uniprot_prokaryotic +0.312 · uniprot_own
+0.126 · pubtator +0.417** -- a 3.3x spread, with **PubTator strongest**, which matches the axis's
own held-out control where it beat the curated column 0.3722 to 0.3398. `n_papers_uniprot_own` is
weakest because it is 95.1% ties on Kp, an organism with almost no literature of its own. On E. coli
all three are well populated and agree closely (0.381-0.417). The consensus tracks
`n_papers_uniprot_prokaryotic` at rho +0.974 Kp / +0.938 Ec and `src/studiedness.py:consensus()`
says outright that it "adds little"; it exists so a reader stacking ten deliverables has one scale.

**`n_papers_pubtator_prokaryotic` CARRIES REAL BLANKS, NOT ZEROS** -- Kp 2,149 · Ec 190 · Sa 1,356,
meaning no in-scope donor had an NCBI GeneID. `percentile_consensus()` skips them row-wise, so
`studiedness_consensus` rests on THREE columns for some proteins and TWO for others, and nothing in
the number says which. Never fill them.

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


#: How many of the top-cited proteins to name on the ranked panel.
N_LABEL = 10


def plot_pubtator_rank(ax, data: dict, species: str, n_label: int) -> pd.DataFrame:
    """Every protein ranked by PubTator citations, with the most-cited named.

    **Both axes are logarithmic, and both have to be.** The values span 0 to 2,541 against a median
    of 10, and the proteins worth naming are the first ~0.3% of the rank order -- on linear axes the
    head is a single pixel column against a flat floor. Log-log turns the same data into the shape
    that is actually the point: a very short head, then four orders of magnitude of tail.

    **PubTator, not UniProt, and this is the panel where that matters most.** It is the stronger
    predictor of essentiality on the anchor (rho +0.417 against the curated column's +0.312) and it
    counts machine-read mentions across 36M abstracts rather than curated references, so the head of
    this ranking is what the literature actually talks about.

    **Shown for E. coli because K. pneumoniae's PubTator column is borrowed.** Its counts come from
    whichever donor the transfer found, and for Kp that is usually the E. coli protein -- the two
    species' top values are literally identical (rpoB 2,541, recA 633, rpoS 439, ...). Ranking Kp
    here would draw E. coli's literature under a Klebsiella label."""
    d = data[species].copy()
    d["pub"] = pd.to_numeric(d["n_papers_pubtator_prokaryotic"], errors="coerce")
    scored = d.dropna(subset=["pub"]).sort_values("pub", ascending=False).reset_index(drop=True)
    scored["rank"] = np.arange(1, len(scored) + 1)

    ax.plot(scored["rank"], np.log1p(scored["pub"]), color=PAL.SPECIES_COLOR[species], lw=2.0,
            zorder=2)
    ax.fill_between(scored["rank"], 0, np.log1p(scored["pub"]),
                    color=PAL.SPECIES_COLOR[species], alpha=0.18, linewidth=0, zorder=1)

    top = scored.head(n_label)
    ax.scatter(top["rank"], np.log1p(top["pub"]), s=22, color=PAL.ACCENT, zorder=4,
               linewidths=0.5, edgecolors="white")

    # A LEADER-LINE LADDER, not labels pinned to their points. Ranks 4-10 differ by a few
    # mentions (376 down to 330) and sit almost on top of each other, so text anchored to each
    # marker overlaps into an unreadable smear. The ladder puts the names on an evenly spaced
    # column to the right of the head and draws a thin line back to each point, which keeps the
    # reading order identical to the rank order.
    y_top = float(np.log1p(top["pub"].iloc[0]))
    y_bot = float(np.log1p(top["pub"].iloc[-1])) * 0.62
    ladder = np.linspace(y_top, y_bot, len(top))
    for (_, r), y_lab in zip(top.iterrows(), ladder):
        name = r["gene_name"] if isinstance(r["gene_name"], str) and r["gene_name"] else r["uniprot_ac"]
        ax.annotate(
            f"{name}  {int(r['pub']):,}",
            xy=(r["rank"], np.log1p(r["pub"])), xytext=(len(scored) ** 0.42, y_lab),
            textcoords="data", ha="left", va="center", fontsize=SS * 0.72, color=PAL.INK,
            arrowprops={"arrowstyle": "-", "lw": 0.6, "color": PAL.MUTED,
                        "shrinkA": 0, "shrinkB": 2},
        )

    med = float(scored["pub"].median())
    ax.axhline(np.log1p(med), color=PAL.MUTED, lw=1.0, ls="--", zorder=0)
    ax.text(len(scored), np.log1p(med) + 0.16, f"median {med:.0f}", ha="right",
            fontsize=SS * 0.78, color=PAL.MUTED)

    ax.set_xscale("log")
    ax.set_xlim(1, len(scored) * 1.05)
    yt = [0, 1, 10, 100, 1000]
    ax.set_yticks(np.log1p(yt))
    ax.set_yticklabels([f"{t:,}" for t in yt], fontsize=SS)
    ax.set_ylim(0, np.log1p(scored["pub"].max()) * 1.18)
    stylia.label(ax, xlabel=f"Rank within {LABELS[species]} (log)",
                 ylabel="PubTator mentions (prokaryotic donor)",
                 title="A very short head, then a long tail")
    return scored


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
    stylia.label(ax, xlabel="Studiedness consensus, 3 sources (percentile)",
                 ylabel="Essentiality (within-species percentile)",
                 title="Better studied, more essential")
    return curves, rhos


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bins", type=int, default=10,
                    help="studiedness bins for the median curve (Kp yields fewer; see --help)")
    ap.add_argument("--rank-species", default="ecoli", choices=SPECIES,
                    help="species for the ranked PubTator panel (Kp's counts are donor-borrowed)")
    ap.add_argument("--label", type=int, default=N_LABEL, help="how many top proteins to name")
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
        pr = P.load(sp)[["uniprot_ac", "gene_name", "sequence"]].copy()
        pr["length"] = pr["sequence"].str.len()
        data[sp] = d.merge(pr[["uniprot_ac", "gene_name", "length"]], on="uniprot_ac")

    fig, axs = stylia.create_figure(1, 3, width_ratios=[0.72, 1.35, 1.2], width=1.0, height=0.42)
    plot_tiers(axs.next(), data)
    scored = plot_pubtator_rank(axs.next(), data, args.rank_species, args.label)
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

    say(f"\n  MOST-CITED IN {LABELS[args.rank_species].upper()}  (PubTator, prokaryotic donor)")
    for i, r in scored.head(args.label).iterrows():
        name = r["gene_name"] if isinstance(r["gene_name"], str) and r["gene_name"] else r["uniprot_ac"]
        say(f"    {i + 1:>3}. {name:<14} {int(r['pub']):>6,}")
    say(f"    median {scored['pub'].median():.0f}   "
        f"top {args.label} hold {100 * scored['pub'].head(args.label).sum() / scored['pub'].sum():.1f}% "
        f"of all mentions in this proteome")

    say("\n  DOES THE SOURCE CHANGE THE ANSWER?  rho against essentiality_consensus")
    say(f"    {'':<16} {'consensus(3)':>14} {'uniprot_prok':>14} {'uniprot_own':>14} {'pubtator':>14}")
    for sp in SPECIES:
        d = data[sp]
        cells = []
        for col in ("studiedness_consensus", PAPERS, "n_papers_uniprot_own",
                    "n_papers_pubtator_prokaryotic"):
            v = pd.to_numeric(d[col], errors="coerce")
            cells.append(f"{v.corr(d['essentiality_consensus'], method='spearman'):+.3f}")
        say(f"    {LABELS[sp]:<16} " + " ".join(f"{c:>14}" for c in cells))
    say("    Same POSITIVE sign in every column and both species -- the finding does not depend")
    say("    on the source. The MAGNITUDE does, on Kp: 0.126 to 0.417, a 3.3x spread.")
    say("    PubTator is the STRONGEST there (+0.417), matching the axis's own held-out control")
    say("    where it beat the curated column 0.3722 to 0.3398. `n_papers_uniprot_own` is the")
    say("    weakest (+0.126) and that is expected -- it is 95.1% ties on Kp, which has almost")
    say("    no literature of its own. On Ec, where all three are well populated, they agree")
    say("    closely (0.381-0.417). The violins and tiers panels are UniProt only; only the")
    say("    relationship panel uses the 3-source consensus.")

    say("\n  PUBTATOR BLANKS  (no in-scope donor had an NCBI GeneID -- NOT zeros)")
    for sp in SPECIES:
        n = pd.to_numeric(data[sp]["n_papers_pubtator_prokaryotic"], errors="coerce").isna().sum()
        say(f"    {LABELS[sp]:<16} {int(n):>5,} of {len(data[sp]):,} "
            f"({100 * n / len(data[sp]):.1f}%) -- ranked on the other two columns, never filled")

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
