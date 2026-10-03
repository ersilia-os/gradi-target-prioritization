"""Studiedness figures: what the axis measured, and whether the transfer can be believed.

    python scripts/plots/studiedness.py

    output/plots/studiedness/studiedness.png     the distributions and the own-vs-family gap
    output/plots/studiedness/control.png         the E. coli held-out control + donor organisms
    output/plots/studiedness/interest_panel.png  where the consortium's own targets sit

**Run plot scripts ONE AT A TIME.** `stylia/__init__.py` calls `shutil.rmtree` on matplotlib's
cache directory at import time, so two plot scripts running concurrently delete each other's font
cache and die inside `shutil` naming `~/.matplotlib`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import stylia  # noqa: E402
from src import interest as I  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import studiedness as S  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "studiedness"
NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL

LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
TIER_COLOUR = {"swissprot_direct": NPG[0], "swissprot_close": NPG[1], "swissprot_homolog": NPG[2],
               "below_floor": "#bdbdbd", "no_hit": "#eeeeee"}
TIERS = ("swissprot_direct", "swissprot_close", "swissprot_homolog", "below_floor", "no_hit")


def fig_distributions(species: list[str]) -> Path:
    """Own vs family per species, the gap between them, and the evidence tiers."""
    stylia.set_format("slide")
    fig, axs = stylia.create_figure(1, 3, width=1.0, height=0.55)

    # Proteins with no donor score exactly 0 and would otherwise be a spike that flattens the
    # rest of the distribution. They are counted in the legend and shown as their own tiers in
    # panel 3, so nothing is hidden by plotting the scored proteins here.
    # A SURVIVAL CURVE, not a histogram. The counts are small integers, so log-spaced bins
    # alias -- some bins catch two integers and some one, producing alternating heights that
    # read as structure and are pure binning artifact. "What fraction of the proteome has at
    # least N papers" is monotone, bin-free and directly readable.
    ax = axs[0]
    for i, sp in enumerate(species):
        v = S.load(sp)["n_papers_uniprot_prokaryotic"].astype(float).to_numpy()
        grid = np.arange(0, 61)
        frac = [(v >= n).mean() * 100 for n in grid]
        ax.step(grid, frac, where="post", lw=1.8, color=NPG[i],
                label=f"{LABELS[sp]}  ({(v == 0).mean() * 100:.0f}% zero)")
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=SS, frameon=False)
    stylia.label(ax, xlabel="papers on the best-studied homolog (N)",
                 ylabel="% of proteome with at least N",
                 title="What is known about the family")

    # +1 so the many zero-paper accessions stay visible on a log axis.
    ax = axs[1]
    for i, sp in enumerate(species):
        d = S.load(sp)
        ax.scatter(d["n_papers_uniprot_own"].astype(float) + 1, d["n_papers_uniprot_prokaryotic"].astype(float) + 1,
                   s=3, alpha=0.25, color=NPG[i], edgecolors="none", label=LABELS[sp])
    lim = [1, 400]
    ax.plot(lim, lim, ls="--", lw=0.8, color="#888")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    leg = ax.legend(fontsize=SS, frameon=False, markerscale=4)
    for h in leg.legend_handles:
        h.set_alpha(1.0)
    stylia.label(ax, xlabel="papers on this accession (+1)",
                 ylabel="papers on the family (+1)",
                 title="Dark accession, known family")

    ax = axs[2]
    tiers = list(TIERS)
    bottoms = [0.0] * len(species)
    xs = range(len(species))
    for tier in tiers:
        vals = []
        for sp in species:
            d = S.load(sp)
            vals.append(100 * float((d["evidence"] == tier).mean()))
        ax.bar(xs, vals, bottom=bottoms, color=TIER_COLOUR[tier], label=tier, width=0.6)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax.set_xticks(list(xs))
    ax.set_xticklabels([LABELS[s] for s in species], fontsize=SS, rotation=20, ha="right")
    ax.legend(fontsize=SS, frameon=False, loc="lower right")
    stylia.label(ax, xlabel="", ylabel="% of proteome", title="Where the number came from")

    path = OUT_DIR / "studiedness.png"
    stylia.save_figure(str(path))
    return path


def fig_control(species: list[str]) -> Path:
    """The held-out control, and who the donors actually are."""
    stylia.set_format("slide")
    fig, axs = stylia.create_figure(1, 2, width=1.0, height=0.58)

    ax = axs[0]
    try:
        ctrl = S.control()
    except FileNotFoundError:
        ctrl = None
    if ctrl is not None and len(ctrl):
        for c in ("n_papers_uniprot_own", "n_papers_uniprot_prokaryotic"):
            ctrl[c] = ctrl[c].astype(float)
        scored = ctrl[ctrl["donor_ac"].astype(str) != ""]
        ax.scatter(scored["n_papers_uniprot_own"] + 1, scored["n_papers_uniprot_prokaryotic"] + 1, s=4, alpha=0.3,
                   color=NPG[0], edgecolors="none")
        rho = scored["n_papers_uniprot_prokaryotic"].corr(scored["n_papers_uniprot_own"], method="spearman")
        ax.plot([1, 400], [1, 400], ls="--", lw=0.8, color="#888")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(1, 400); ax.set_ylim(1, 400)
        ax.text(0.04, 0.94, f"rho {rho:.3f}   n {len(scored):,}", transform=ax.transAxes,
                fontsize=SS, va="top")
    stylia.label(ax, xlabel="E. coli measured papers (+1)",
                 ylabel="papers from non-Escherichia donors (+1)",
                 title="Held-out transfer control")

    ax = axs[1]
    counts: dict[str, int] = {}
    for sp in species:
        tr = S.load_transfer(sp)
        for org in tr.loc[tr["donor_organism"].astype(str) != "", "donor_organism"]:
            short = str(org).split(" (")[0]
            counts[short] = counts.get(short, 0) + 1
    top = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:12][::-1]
    ax.barh(range(len(top)), [v for _, v in top], color=NPG[2], height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([k for k, _ in top], fontsize=SS)
    stylia.label(ax, xlabel="proteins served", ylabel="", title="Who the donors are")

    path = OUT_DIR / "control.png"
    stylia.save_figure(str(path))
    return path


def fig_interest(species: list[str]) -> Path:
    """Are the consortium's own targets novel, or already well studied?"""
    stylia.set_format("slide")
    fig, axs = stylia.create_figure(1, len(species), width=1.0, height=0.55)
    axs = [axs] if len(species) == 1 else list(axs)

    for i, sp in enumerate(species):
        ax = axs[i]
        d = S.load(sp).merge(P.load(sp)[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left")
        d = I.annotate(d)
        panel = d[d["is_interest"]]
        rest = d[~d["is_interest"]]
        grid = np.arange(0, 61)
        rest_v = rest["n_papers_uniprot_prokaryotic"].astype(float).to_numpy()
        ax.step(grid, [(rest_v >= n).mean() * 100 for n in grid], where="post",
                lw=1.6, color="#999",
                label=f"proteome ({(rest_v == 0).mean() * 100:.0f}% zero)")
        ax.set_xlim(0, 60)
        ax.set_ylim(0, 100)
        if len(panel):
            # Percentile against the WHOLE proteome, zeros included -- the panel's standing is a
            # claim about the proteome, not about the scored subset.
            pct = [100 * float((d["n_papers_uniprot_prokaryotic"] < v).mean())
                   for v in panel["n_papers_uniprot_prokaryotic"]]
            pv = panel["n_papers_uniprot_prokaryotic"].astype(float).to_numpy()
            ax.scatter(pv, [(rest_v >= n).mean() * 100 for n in pv],
                       s=20, color=NPG[3], zorder=5,
                       label=f"panel ({len(panel)}), median {np.median(pct):.0f}th pct")
        ax.legend(fontsize=SS, frameon=False, loc="upper left")
        stylia.label(ax, xlabel="papers on the family (N)",
                     ylabel="% with at least N" if i == 0 else "", title=LABELS[sp])

    path = OUT_DIR / "interest_panel.png"
    stylia.save_figure(str(path))
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(S.SPECIES), choices=list(S.SPECIES))
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for fn in (fig_distributions, fig_control, fig_interest):
        path = fn(args.species)
        print(f"wrote {path.relative_to(REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
