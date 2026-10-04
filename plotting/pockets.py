"""Could a small molecule bind this fold? Structure for nearly everything, measurement for very few.

    pockets.png   A  structural coverage: AlphaFold models vs experimental PDB structures
                  B  measured vs modelled ligands -- two different kinds of evidence
                  C  the honest control: pocket scores against protein LENGTH

**Panel C is the caveat that has to be on the slide.** Length alone predicts a measured ligand at
AUROC 0.65-0.67 -- big proteins are crystallised more and have more surface. Within length deciles,
against `n_ligands_pdb > 0`, P2Rank scores 0.494 on K. pneumoniae (0.561 Ec, 0.619 Sa) and fpocket
0.435. So **on the anchor the pocket scores add nothing over protein size**, and they are a soft
prior at best. Always quote the length-controlled number, never the raw AUROC.

**NEVER SUM THE TWO LIGAND COUNTS.** A co-crystal of this protein (88 Kp proteins) and a transplant
from a ~30%-identity homolog (1,533 Kp) are different evidence. Both are non-redundant distinct
Bemis-Murcko scaffolds.

**An NA in the pocket columns means "could not look", not "looked and found nothing"** -- a protein
WITH a model and no admitted pocket gets a 0. Never `fillna(0)`; that was the v1 mistake.

A pocket counts only if its residues average pLDDT >= 70, and confidence enters exactly ONCE (v1
applied it twice). "Drug-like" is built from published sources -- BioLiP, PLINDER, PDBe cofactor
classes, a nucleotide SMARTS and ECMDB -- because QED >= 0.2 and Ro3 were measured and REJECTED:
both delete antibiotics (novobiocin 0.184, rifampicin 0.109).

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
`stylia.set_style("ersilia")` MUST precede `NamedColors()`, or NC.plum raises AttributeError.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/pockets.py
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

from src import pockets as PK  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
SHORT = {"kpneumoniae": "Kp", "ecoli": "Ec", "saureus": "Sa"}
SPECIES = ["kpneumoniae", "ecoli", "saureus"]


def _auroc(score: np.ndarray, label: np.ndarray) -> float:
    """Rank-based AUROC. No sklearn import for three lines of arithmetic."""
    ok = ~np.isnan(score)
    s, y = score[ok], label[ok]
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan")
    r = pd.Series(s).rank().to_numpy()
    n1, n0 = y.sum(), (~y.astype(bool)).sum()
    return float((r[y.astype(bool)].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def plot_coverage(ax, data: dict, abc: str) -> None:
    """Modelled vs experimentally determined. The gap is the point: structure is nearly universal,
    experimental structure is not."""
    x = np.arange(len(SPECIES))
    width = 0.38
    modelled, solved = [], []
    for sp in SPECIES:
        pk = data[sp]
        modelled.append(pd.to_numeric(pk["af_plddt"], errors="coerce").notna().mean() * 100)
        solved.append((pd.to_numeric(pk["n_pdb_structures"], errors="coerce").fillna(0) > 0).mean() * 100)
    ax.bar(x - width / 2, modelled, width=width, color=NC.plum, label="AlphaFold model")
    ax.bar(x + width / 2, solved, width=width, color=NC.orange, label="experimental PDB")
    for xi, m, s in zip(x, modelled, solved):
        ax.text(xi - width / 2, m + 1, f"{m:.0f}%", ha="center", fontsize=SS, color=NC.black)
        ax.text(xi + width / 2, s + 1, f"{s:.0f}%", ha="center", fontsize=SS, color=NC.black)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS, style="italic")
    ax.set_ylim(0, 112)
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="", ylabel="% of proteome", title="Structure is not the bottleneck",
                 abc=abc)


def plot_ligands(ax, data: dict, abc: str) -> None:
    """Measured co-crystal ligands vs AlphaFill transplants. Side by side and NEVER summed."""
    x = np.arange(len(SPECIES))
    width = 0.38
    pdb = [int((pd.to_numeric(data[sp]["n_ligands_pdb"], errors="coerce").fillna(0) > 0).sum())
           for sp in SPECIES]
    af = [int((pd.to_numeric(data[sp]["n_ligands_alphafill"], errors="coerce").fillna(0) > 0).sum())
          for sp in SPECIES]
    ax.bar(x - width / 2, pdb, width=width, color=NC.orange, label="measured (own PDB)")
    ax.bar(x + width / 2, af, width=width, color=NC.blue, label="modelled (AlphaFill)")
    for xi, a, b in zip(x, pdb, af):
        ax.text(xi - width / 2, a + 18, f"{a}", ha="center", fontsize=SS, color=NC.black)
        ax.text(xi + width / 2, b + 18, f"{b:,}", ha="center", fontsize=SS, color=NC.black)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SPECIES], fontsize=SS, style="italic")
    ax.set_ylim(0, max(af) * 1.28)
    ax.legend(fontsize=SS, frameon=False, loc="upper center", ncol=2)
    stylia.label(ax, xlabel="", ylabel="proteins with a drug-like ligand",
                 title="Never sum these two", abc=abc)


def plot_length_control(ax, data: dict, lengths: dict, abc: str) -> None:
    """Raw AUROC against length-controlled AUROC, for both pocket scores.

    The honest comparison. Raw says the scores work; within length deciles they do not, on the
    anchor. Drawn side by side so the gap cannot be quoted away."""
    rows = []
    for sp in SPECIES:
        pk = data[sp].merge(lengths[sp], on="uniprot_ac", how="left")
        y = (pd.to_numeric(pk["n_ligands_pdb"], errors="coerce").fillna(0) > 0).to_numpy()
        pk["_decile"] = pd.qcut(pk["length"], 10, labels=False, duplicates="drop")
        for col, name in (("p2rank_score", "P2Rank"), ("fpocket_score", "fpocket")):
            s = pd.to_numeric(pk[col], errors="coerce").to_numpy()
            raw = _auroc(s, y)
            within = [
                _auroc(s[m], y[m]) for d in pk["_decile"].dropna().unique()
                for m in [(pk["_decile"] == d).to_numpy()]
            ]
            within = np.nanmean([w for w in within if not np.isnan(w)])
            rows.append((sp, name, raw, within))

    x = np.arange(len(rows))
    ax.bar(x - 0.19, [r[2] for r in rows], width=0.38, color=NC.gray, label="raw AUROC")
    ax.bar(x + 0.19, [r[3] for r in rows], width=0.38, color=NC.plum, label="within length deciles")
    ax.axhline(0.5, color=NC.black, lw=1.1, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r[1]}\n{SHORT[r[0]]}" for r in rows], fontsize=SS * 0.78)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="", ylabel="AUROC vs a measured ligand",
                 title="Control for length before believing this", abc=abc)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  pockets -- structural ligandability")
    data = {sp: PK.load(sp) for sp in SPECIES}
    lengths = {}
    for sp in SPECIES:
        pr = P.load(sp)[["uniprot_ac", "sequence"]].copy()
        pr["length"] = pr["sequence"].str.len()
        lengths[sp] = pr[["uniprot_ac", "length"]]

    fig, axs = stylia.create_figure(1, 3, width_ratios=[2.4, 2.4, 3.2], width=1.0, height=0.38)
    plot_coverage(axs.next(), data, abc="A")
    plot_ligands(axs.next(), data, abc="B")
    rows = plot_length_control(axs.next(), data, lengths, abc="C")
    out = OUT_DIR / "pockets.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  COVERAGE")
    for sp in SPECIES:
        pk = data[sp]
        say(f"    {LABELS[sp]:<16} model {pd.to_numeric(pk['af_plddt'], errors='coerce').notna().mean() * 100:>5.1f}%"
            f"   PDB {int((pd.to_numeric(pk['n_pdb_structures'], errors='coerce').fillna(0) > 0).sum()):>5,}"
            f"   own ligands {int((pd.to_numeric(pk['n_ligands_pdb'], errors='coerce').fillna(0) > 0).sum()):>4,}"
            f"   AlphaFill {int((pd.to_numeric(pk['n_ligands_alphafill'], errors='coerce').fillna(0) > 0).sum()):>5,}")

    say("\n  THE LENGTH CONTROL (quote the second number, never the first)")
    for sp, name, raw, within in rows:
        say(f"    {LABELS[sp]:<16} {name:<8} raw {raw:.3f}   within length deciles {within:.3f}")
    say("\n  CAVEATS")
    say("    - On the anchor the pocket scores add nothing over protein size. A soft prior.")
    say("    - NA means 'could not look'. A modelled protein with no admitted pocket gets 0.")
    say("      Never fillna(0) -- that was the v1 mistake.")
    say("    - No pocket COUNT is shipped: without a probability cutoff it is just size.")


if __name__ == "__main__":
    main()
