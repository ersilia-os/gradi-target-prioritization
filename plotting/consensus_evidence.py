"""What `<axis>_consensus` and `<axis>_evidence` are, and why an axis ships both.

    consensus_evidence.png   1  the consensus is a blend, not a copy of its best input
                             2  where the evidence comes from -- same ladder, opposite provenance
                             3  the two columns are near-orthogonal, and mildly ANTI-correlated

Every axis in this project ships the same pair: `<axis>_consensus`, a mean of WITHIN-SPECIES
percentile ranks, and `<axis>_evidence`, an integer 1-3 grading how well corroborated that value is.
**Essentiality is the worked example because it is the only axis with several independent predictors
AND several independent experiments.** The convention is defined in `src/consensus.py`; these
columns are built in `scripts/essentiality/merge.py:consensus_and_evidence()`. This figure only
describes what ships -- it reimplements none of that logic.

**THE CONSENSUS IS PREDICTORS ONLY. No measurement enters it**, so the evidence stays an independent
statement rather than a restatement of the score, and the label-transfer this axis rejected stays
rejected.

**The ladder is count + concordance**: 3 = two or more independent sources, unanimous, AND agreeing
with the consensus · 2 = one source, or several that conflict · 1 = no measurement. Level 0 cannot
occur, because every protein has a prediction. **"One source" means one DATASET OR STRAIN, not one
family.** Concordance is judged against a **base-rate cut, never 0.5** -- essentials are 11-17% of a
proteome -- recorded per species in `evidence/consensus_audit.tsv` (Kp 0.824, Ec 0.738).

**THE EVIDENCE IS NOT PURELY EXPERIMENTAL, and that is the subtlety the third panel exists for.**
A protein measured twice, unanimously, whose consensus contradicts it lands at **2**, not 3 -- never
read a 2 as "the experiment was weak". And concordance is tested with equality, so a unanimous
NEGATIVE sitting below the cut is concordant too. Tier 3 therefore fills with corroborated
NON-essentials, which are the overwhelming majority of any proteome, while a protein the models rank
highly but the screens called dispensable is demoted. The measured consequence is that consensus and
evidence run mildly NEGATIVE together: rho -0.216 on Kp, -0.349 on E. coli.

**Essentiality does NOT deviate from the generic ladder.** The documented deviations in
`src/consensus.py` belong to `degradability` (reads 3 as ONE measurement) and `ligands` (pins its
zero block to 0 instead of the average rank). Do not attribute either to this axis.

Everything is read from the stage's own tables; the two audit tables have no loader and are read by
path, which is the documented exception rather than the rule.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/consensus_evidence.py
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
from matplotlib.lines import Line2D  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import essentiality as E  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# The two audit tables have NO loader in `src/`. Reading them by path is the documented exception
# -- `scripts/plots/proteomes.py` does the same for `name_audit.tsv` -- not a licence to bypass the
# loaders for anything that has one.
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "essentiality" / "evidence"
CONSENSUS_AUDIT = EVIDENCE_DIR / "consensus_audit.tsv"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
SPECIES = ["kpneumoniae", "ecoli"]
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli"}

#: The three columns `percentile_consensus()` averages. Named here only to label the panel; the
#: list that actually builds the column is CONSENSUS_INPUTS in scripts/essentiality/merge.py.
PREDICTORS = ["geptop_ess", "proteomelm_ess", "screens_ess_mean"]
PRED_SHORT = {"geptop_ess": "Geptop", "proteomelm_ess": "ProteomeLM",
              "screens_ess_mean": "Screens"}

#: Evidence tiers. Same ramp as `studiedness_essentiality.png`, so a reader who has seen one
#: evidence bar in this deck reads this one the same way.
TIERS = {
    1: ("1  no measurement", PAL.NPG[0]),
    2: ("2  one source, or they conflict", PAL.NPG[1]),
    3: ("3  >=2 sources, unanimous, concordant", PAL.NPG[4]),
}


def plot_blend(ax, data: dict) -> dict:
    """Each predictor against the consensus, and against the other two predictors.

    The gap is the panel. Inputs that agree with each other at rho 0.32-0.44 each land at
    0.69-0.78 against the blend -- which is what rank-averaging buys, and why the consensus is not
    simply a copy of whichever input is best. Drawn for the anchor only; E. coli behaves the same
    and the second species would double the rows without adding an argument."""
    sp = "kpneumoniae"
    d = data[sp]
    cons = pd.to_numeric(d["essentiality_consensus"], errors="coerce")
    stats = {}
    for col in PREDICTORS:
        v = pd.to_numeric(d[col], errors="coerce")
        others = [c for c in PREDICTORS if c != col]
        stats[col] = {
            "to_consensus": v.corr(cons, method="spearman"),
            "to_others": float(np.mean([
                v.corr(pd.to_numeric(d[o], errors="coerce"), method="spearman") for o in others])),
            "ties": (v == 0).mean() * 100,
        }

    y = np.arange(len(PREDICTORS))[::-1]
    height = 0.36
    ax.barh(y + height / 2, [stats[c]["to_consensus"] for c in PREDICTORS], height=height,
            color=PAL.PRIMARY, label="to the consensus")
    ax.barh(y - height / 2, [stats[c]["to_others"] for c in PREDICTORS], height=height,
            color=PAL.MUTED, label="to the other two predictors")
    for yi, col in zip(y, PREDICTORS):
        ax.text(stats[col]["to_consensus"] + 0.015, yi + height / 2,
                f"{stats[col]['to_consensus']:.2f}", va="center", fontsize=SS * 0.8, color=PAL.INK)
        ax.text(stats[col]["to_others"] + 0.015, yi - height / 2,
                f"{stats[col]['to_others']:.2f}", va="center", fontsize=SS * 0.8, color=PAL.INK)

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{PRED_SHORT[c]}" + (f"\n{stats[c]['ties']:.0f}% tied at 0" if stats[c]["ties"] > 1 else "")
         for c in PREDICTORS], fontsize=SS * 0.85)
    ax.set_xlim(0, 1.0)
    ax.legend(fontsize=SS * 0.8, frameon=False, loc="lower right", handletextpad=0.5)
    stylia.label(ax, xlabel="Spearman rho", ylabel="",
                 title="Three predictors, one score")
    return stats


def plot_provenance(ax, data: dict, audit: pd.DataFrame) -> None:
    """The ladder per species, with the provenance of the experimental tier written on it.

    The bars are nearly the same height in the two species and mean opposite things, which is the
    only reason this panel exists: K. pneumoniae's entire experimental tier is a >=95%-identity
    counterpart in ANOTHER Klebsiella strain, because DEG contains no Klebsiella at all, while
    E. coli's is measured on the anchor strain itself."""
    x = np.arange(len(SPECIES))
    bottom = np.zeros(len(SPECIES))
    for tier, (label, color) in TIERS.items():
        pct, counts = [], []
        for sp in SPECIES:
            n = int((data[sp]["essentiality_evidence"] == tier).sum())
            counts.append(n)
            pct.append(n / len(data[sp]) * 100)
        ax.bar(x, pct, bottom=bottom, width=0.58, color=color, label=label,
               edgecolor="white", linewidth=0.9)
        for xi, v, b, n in zip(x, pct, bottom, counts):
            if v >= 5:
                ax.text(xi, b + v / 2, f"{n:,}", ha="center", va="center",
                        fontsize=SS * 0.88, color="white")
        bottom += np.array(pct)

    rows = {r["species"]: r for _, r in audit.iterrows()}
    ticks = []
    for sp in SPECIES:
        r = rows[sp]
        prox = int(r["n_proxy_only"]) / max(int(r["n_covered"]), 1) * 100
        ticks.append(f"{LABELS[sp]}\n{int(r['n_covered']):,} covered · "
                     f"{prox:.0f}% proxy")
    ax.set_xticks(x)
    ax.set_xticklabels(ticks, fontsize=SS * 0.82)
    ax.set_xlim(-0.6, len(SPECIES) - 0.4)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=SS * 0.74, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.235), handletextpad=0.4, labelspacing=0.25)
    stylia.label(ax, xlabel="", ylabel="% of proteome",
                 title="Same ladder, opposite provenance")


def plot_orthogonality(ax, data: dict) -> dict:
    """Consensus distribution by evidence tier -- the panel that justifies shipping both columns.

    Box plots rather than violins: the comparison is between six medians and their spread, and six
    kernels on one axis is more ink than information. Whiskers at 5-95% so the tails do not drag
    the boxes flat."""
    rhos = {}
    positions, box_data, colors = [], [], []
    for i, sp in enumerate(SPECIES):
        d = data[sp]
        cons = pd.to_numeric(d["essentiality_consensus"], errors="coerce")
        ev = pd.to_numeric(d["essentiality_evidence"], errors="coerce")
        rhos[sp] = cons.corr(ev, method="spearman")
        for j, tier in enumerate((1, 2, 3)):
            positions.append(i * 3.6 + j)
            box_data.append(cons[ev == tier].dropna().to_numpy())
            colors.append(TIERS[tier][1])

    bp = ax.boxplot(box_data, positions=positions, widths=0.72, patch_artist=True,
                    showfliers=False, whis=(5, 95))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_edgecolor(PAL.INK)
        patch.set_linewidth(0.8)
    for el in ("medians", "whiskers", "caps"):
        for ln in bp[el]:
            ln.set_color(PAL.INK)
            ln.set_linewidth(0.9)

    ax.set_xticks(positions)
    ax.set_xticklabels([str(t) for _sp in SPECIES for t in (1, 2, 3)], fontsize=SS)
    for i, sp in enumerate(SPECIES):
        ax.text(i * 3.6 + 1, -0.10, f"{LABELS[sp]}\nρ {rhos[sp]:+.2f}", ha="center", va="top",
                fontsize=SS * 0.85, color=PAL.SPECIES_COLOR[sp],
                transform=matplotlib.transforms.blended_transform_factory(ax.transData,
                                                                         ax.transAxes))
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.8, 3.6 + 2 + 0.8)
    stylia.label(ax, xlabel="", ylabel="Essentiality consensus",
                 title="A high score is not a well-evidenced score")
    return rhos


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  consensus + evidence, worked on the essentiality axis")
    data = {sp: E.load(sp) for sp in SPECIES}
    audit = pd.read_csv(CONSENSUS_AUDIT, sep="\t")
    audit = audit[audit["species"] != "species"]          # the file repeats its header once

    fig, axs = stylia.create_figure(1, 3, width_ratios=[1.15, 0.85, 1.25], width=1.0, height=0.40)
    stats = plot_blend(axs.next(), data)
    plot_provenance(axs.next(), data, audit)
    rhos = plot_orthogonality(axs.next(), data)
    out = OUT_DIR / "consensus_evidence.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  THE CONSENSUS IS A BLEND  (K. pneumoniae)")
    say(f"    {'predictor':<14} {'rho to consensus':>18} {'rho to other two':>18} {'tied at 0':>11}")
    for col in PREDICTORS:
        s = stats[col]
        say(f"    {PRED_SHORT[col]:<14} {s['to_consensus']:>18.3f} {s['to_others']:>18.3f}"
            f" {s['ties']:>10.1f}%")
    say("    Pairwise between predictors:")
    d = data["kpneumoniae"]
    for a, b in itertools.combinations(PREDICTORS, 2):
        r = pd.to_numeric(d[a], errors="coerce").corr(pd.to_numeric(d[b], errors="coerce"),
                                                      method="spearman")
        say(f"      {PRED_SHORT[a]:<12} ~ {PRED_SHORT[b]:<12} {r:+.3f}")

    say("\n  THE EVIDENCE LADDER")
    for _, r in audit.iterrows():
        if r["species"] not in LABELS:
            continue
        say(f"    {LABELS[r['species']]:<16} 1:{int(r['level_1']):>5,}  2:{int(r['level_2']):>5,}  "
            f"3:{int(r['level_3']):>5,}   covered {int(r['n_covered']):>5,}   "
            f"proxy-only {int(r['n_proxy_only']):>5,}   "
            f"base rate {float(r['base_rate']):.4f}   concordance cut {float(r['consensus_cut']):.3f}")
    say("    Kp's experimental tier is 100% PROXY: DEG has no Klebsiella, so every call is a")
    say("    >=95%-identity counterpart in ATCC 43816, ECL8 or RH201207 -- and NOT ONE of the")
    say("    3,973 level-3 proteins is measured on HS11286 itself. Ec's is 0% proxy: two DEG")
    say("    datasets on the anchor strain plus the OGEE label.")

    say("\n  THE TWO COLUMNS ARE NOT THE SAME THING")
    for sp in SPECIES:
        d = data[sp]
        cons = pd.to_numeric(d["essentiality_consensus"], errors="coerce")
        ev = pd.to_numeric(d["essentiality_evidence"], errors="coerce")
        say(f"    {LABELS[sp]:<16} rho(consensus, evidence) = {rhos[sp]:+.3f}")
        for tier in (1, 2, 3):
            s = cons[ev == tier]
            say(f"      tier {tier}: n={len(s):>5,}  median {s.median():.3f}  "
                f"below 0.5 {100 * (s < 0.5).mean():>5.1f}%  above 0.9 {100 * (s >= 0.9).mean():>5.1f}%")
        top = d[cons >= cons.quantile(0.9)]
        say(f"      top consensus decile by tier: "
            f"{top['essentiality_evidence'].value_counts().sort_index().to_dict()}")
    say("    NEGATIVE, because concordance is tested with equality: a unanimous NON-essential")
    say("    below the cut is concordant too, so tier 3 fills with corroborated non-essentials.")

    say("\n  CAVEATS")
    say("    - `essentiality_evidence` is NOT purely experimental. A protein measured twice,")
    say("      unanimously, whose consensus contradicts it lands at 2. Never read 2 as 'the")
    say("      experiment was weak'.")
    say("    - Concordance uses a BASE-RATE cut, never 0.5, because essentials are 11-17% of a")
    say("      proteome. The cut is per species and is in evidence/consensus_audit.tsv.")
    say("    - The measured call NEVER enters the consensus -- the evidence is a flag, so the")
    say("      label transfer this axis rejected stays rejected.")
    say("    - `essentiality_consensus` is WITHIN-SPECIES. Never compare the number across one.")
    say("    - A `geptop_ess` of 0 is `orthologs_none_essential` (a confident non-essential call)")
    say("      or `no_orthologs`; read geptop_<sp>.tsv before reading that zero.")


if __name__ == "__main__":
    main()
