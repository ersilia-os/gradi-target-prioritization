"""Several ways of asking whether a protein is druggable, and the two consensus columns.

    druggability_consensus.png   1  ligands: three views, one column
                                 2  pockets: four views, one column
                                 3  what the pocket consensus is NOT -- it is largely length

The two axes each ship the project's standard pair, `<axis>_consensus` (a mean of WITHIN-SPECIES
percentile ranks) and `<axis>_evidence` (1-3). They are built in `src/consensus.py` from

    src/ligandability.py  CONSENSUS_COLUMNS = n_ligands_own, n_ligands_bacterial,
                                              best_pactivity_bacterial
    src/pockets.py        CONSENSUS_COLUMNS = p2rank_score, fpocket_score,
                                              n_ligands_pdb, n_ligands_alphafill

**This figure describes what already ships and combines nothing.** `CLAUDE.md` forbids a composite
score across axes; the two consensus columns here are each WITHIN one axis, and the project keeps
them separate deliberately -- `src/pockets.py`: *"Combine at prioritisation time, across axes, where
the weighting is an explicit choice someone owns."*

**PANEL 3 IS NOT OPTIONAL.** `pockets_consensus` is substantially a ranking by protein length
(rho ~0.66), and within length deciles P2Rank sits at chance on the anchor. A figure that showed
four measures converging on a consensus, without that, would be selling a length proxy as
structural insight.

**NEVER score a consensus against one of its own inputs.** Checking `pockets_consensus` against
`n_ligands_pdb` returns AUROC 0.958 and is circular -- that column is one of the four averaged.
The honest external comparison is protein length, which is what panel 3 draws.

Every number this figure asserts is computed here, from the stage deliverables, and printed to
stdout. The values in `docs/pockets.md` are the expected result, used to check the figure, not
copied into it.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/`. Copying
`parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError names `src`, not the path,
so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/druggability_consensus.py
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

from src import ligandability as L  # noqa: E402
from src import pockets as PK  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
SPECIES = ["kpneumoniae", "ecoli", "saureus"]
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
ANCHOR = "kpneumoniae"

#: Short labels for the panel. The lists that BUILD the columns are CONSENSUS_COLUMNS in
#: src/ligandability.py and src/pockets.py; these are for the y axis only.
LIGAND_INPUTS = {
    "n_ligands_own": "potent, this protein",
    "n_ligands_bacterial": "potent, bacterial pool",
    "best_pactivity_bacterial": "best potency",
}
POCKET_INPUTS = {
    "p2rank_score": "P2Rank (predicted)",
    "fpocket_score": "fpocket (predicted)",
    "n_ligands_pdb": "ligands in own PDB",
    "n_ligands_alphafill": "AlphaFill transplants",
}


def auroc(y_true: np.ndarray, score: np.ndarray) -> float:
    """Rank AUROC, ties averaged. Local because the deck must not import a stage's internals."""
    ok = ~pd.isna(score)
    y, s = np.asarray(y_true)[ok].astype(bool), np.asarray(score, dtype=float)[ok]
    n_pos, n_neg = int(y.sum()), int((~y).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = pd.Series(s).rank(method="average").to_numpy()
    return float((r[y].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auroc_within_deciles(d: pd.DataFrame, score: str, label: str, length: str) -> float:
    """AUROC pooled over protein-length deciles -- the only honest read of a pocket score.

    Big proteins are both crystallised more often and offer more surface, so a raw AUROC rewards
    size. Scoring inside each decile removes that and is where P2Rank drops to chance on Kp.
    """
    d = d.dropna(subset=[score, label, length]).copy()
    d["_dec"] = pd.qcut(d[length].rank(method="first"), 10, labels=False)
    num = den = 0.0
    for _, g in d.groupby("_dec"):
        y = g[label].to_numpy().astype(bool)
        if y.sum() == 0 or (~y).sum() == 0:
            continue
        w = y.sum() * (~y).sum()
        num += auroc(y, g[score].to_numpy()) * w
        den += w
    return num / den if den else float("nan")


def load() -> dict:
    """Both deliverables plus sequence length, per species."""
    out = {}
    for sp in SPECIES:
        lg = L.load(sp)
        pk = PK.load(sp)
        pr = P.load(sp)[["uniprot_ac", "sequence"]].copy()
        pr["length"] = pr["sequence"].str.len()
        d = lg.merge(pk, on="uniprot_ac", how="inner").merge(
            pr[["uniprot_ac", "length"]], on="uniprot_ac", how="left")
        out[sp] = d
    return out


def _input_panel(ax, d: pd.DataFrame, inputs: dict, cons_col: str, title: str) -> dict:
    """Each input's rank correlation to the consensus, and to the other inputs.

    The GAP between the two bars is the panel's argument: inputs that agree only loosely with each
    other all land high against the blend, which is what rank-averaging buys and why the consensus
    is not a copy of whichever input happens to be best.
    """
    cons = pd.to_numeric(d[cons_col], errors="coerce")
    stats = {}
    for col in inputs:
        v = pd.to_numeric(d[col], errors="coerce")
        others = [c for c in inputs if c != col]
        stats[col] = {
            "to_consensus": float(v.corr(cons, method="spearman")),
            "to_others": float(np.mean([v.corr(pd.to_numeric(d[o], errors="coerce"),
                                               method="spearman") for o in others])),
        }

    cols = list(inputs)
    y = np.arange(len(cols))[::-1]
    h = 0.36
    ax.barh(y + h / 2, [stats[c]["to_consensus"] for c in cols], height=h,
            color=PAL.PRIMARY, label="to the consensus")
    ax.barh(y - h / 2, [stats[c]["to_others"] for c in cols], height=h,
            color=PAL.MUTED, label="to the other inputs")
    ax.set_yticks(y)
    ax.set_yticklabels([inputs[c] for c in cols], fontsize=SS * 0.85)
    ax.set_xlim(0, 1)
    ax.axvline(0, color=PAL.INK, lw=0.6)
    ax.legend(fontsize=SS * 0.8, loc="lower right", frameon=False)
    stylia.label(ax, xlabel="Spearman rho", ylabel="", title=title)
    return stats


def plot_ligands(ax, data: dict) -> dict:
    """Ligands: three views of 'has somebody found a compound', and their blend."""
    return _input_panel(ax, data[ANCHOR], LIGAND_INPUTS, "ligands_consensus",
                        f"Ligand evidence, {LABELS[ANCHOR]}")


def plot_pockets(ax, data: dict) -> dict:
    """Pockets: two predictions and two measurements, and their blend."""
    return _input_panel(ax, data[ANCHOR], POCKET_INPUTS, "pockets_consensus",
                        f"Pocket evidence, {LABELS[ANCHOR]}")


def plot_length(ax, data: dict) -> dict:
    """The caveat panel: the pocket consensus is largely a ranking by protein length.

    Two things per species. The DOT is rho(pockets_consensus, length) -- how much of the column is
    size. The BARS are the external check: AUROC for separating proteins that have a drug-like
    ligand in their OWN PDB structures, for protein length alone against P2Rank scored WITHIN
    length deciles. Where the bar sits at 0.5 the pocket score has added nothing over size.

    `n_ligands_pdb` is used as the LABEL here and is also one of the consensus inputs, so the
    consensus itself is deliberately not scored against it -- that comparison is circular and
    returns 0.958. Only `p2rank_score`, which is a different input, is scored.
    """
    stats = {}
    for sp in SPECIES:
        d = data[sp]
        lab = (pd.to_numeric(d["n_ligands_pdb"], errors="coerce").fillna(0) > 0).to_numpy()
        stats[sp] = {
            "rho_cons_length": float(pd.to_numeric(d["pockets_consensus"], errors="coerce").corr(
                d["length"], method="spearman")),
            "auroc_length": auroc(lab, d["length"].to_numpy()),
            "auroc_p2rank_within": auroc_within_deciles(
                d.assign(_lab=lab), "p2rank_score", "_lab", "length"),
            "base": float(lab.mean()),
        }

    x = np.arange(len(SPECIES))
    w = 0.36
    ax.bar(x - w / 2, [stats[s]["auroc_length"] for s in SPECIES], width=w,
           color=PAL.MUTED, label="protein length alone")
    ax.bar(x + w / 2, [stats[s]["auroc_p2rank_within"] for s in SPECIES], width=w,
           color=PAL.PRIMARY, label="P2Rank, within length deciles")
    ax.axhline(0.5, color=PAL.INK, lw=0.8, ls="--")
    ax.text(len(SPECIES) - 0.5, 0.505, "chance", fontsize=SS * 0.8, ha="right", va="bottom",
            color=PAL.INK)

    # rho is NOT drawn on this axis. It is a correlation, not an AUROC, and a dot at 0.66 sitting
    # beside a 0.67 bar reads as if the two were the same quantity. It goes under the species name
    # as text instead.
    for i, sp in enumerate(SPECIES):
        top = max(stats[sp]["auroc_length"], stats[sp]["auroc_p2rank_within"])
        ax.text(i, top + 0.02, rf"$\rho_{{len}}$ {stats[sp]['rho_cons_length']:.2f}",
                fontsize=SS * 0.8, ha="center", va="bottom", color=PAL.ACCENT)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS * 0.85)
    ax.set_ylim(0.3, 1.0)
    ax.legend(fontsize=SS * 0.8, loc="upper left", frameon=False)
    stylia.label(ax, xlabel="", ylabel="AUROC vs a measured PDB ligand",
                 title="Pocket scores add little over size")
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    data = load()
    fig, axs = stylia.create_figure(1, 3, width_ratios=[3, 3, 2.6], width=1.0, height=0.38)
    lig = plot_ligands(axs.next(), data)
    pkt = plot_pockets(axs.next(), data)
    length = plot_length(axs.next(), data)
    out = OUT_DIR / "druggability_consensus.png"
    stylia.save_figure(str(out))

    say(f"  -> {out.relative_to(REPO_ROOT)}")
    say(f"\n  LIGANDS, {LABELS[ANCHOR]} -- rho to consensus / to the other inputs")
    for c, s in lig.items():
        say(f"    {LIGAND_INPUTS[c]:<26} {s['to_consensus']:>6.3f}  {s['to_others']:>6.3f}")
    say(f"\n  POCKETS, {LABELS[ANCHOR]} -- rho to consensus / to the other inputs")
    for c, s in pkt.items():
        say(f"    {POCKET_INPUTS[c]:<26} {s['to_consensus']:>6.3f}  {s['to_others']:>6.3f}")
    say("\n  THE LENGTH CONFOUND -- expected from docs/pockets.md:")
    say("    rho(consensus, length) ~0.663 Kp / 0.638 Ec / 0.654 Sa;"
        " P2Rank within deciles ~0.494 / 0.561 / 0.619")
    for sp in SPECIES:
        s = length[sp]
        say(f"    {LABELS[sp]:<16} rho {s['rho_cons_length']:>6.3f}   "
            f"length alone {s['auroc_length']:>6.3f}   "
            f"P2Rank within length {s['auroc_p2rank_within']:>6.3f}   "
            f"base {100 * s['base']:.1f}%")
    say("\n  The consensus is NOT scored against n_ligands_pdb: that column is one of its four"
        " inputs, so the comparison is circular (it returns 0.958).")


if __name__ == "__main__":
    main()
