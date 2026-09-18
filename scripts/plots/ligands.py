"""Figures for stage 06 part 1 -- what the ChEMBL evidence looks like.

Three questions, three files, all into output/plots/ligands/:

    chembl_coverage.png     how much of each proteome has ligand evidence, and at what distance
    chembl_redundancy.png   does collapsing analog series to scaffolds change the picture
    chembl_selectivity.png  bacterial evidence against human evidence -- the liability plane

Run with the `gradi` env:
    python scripts/plots/ligands.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import interest as I  # noqa: E402
from src import ligandability as L  # noqa: E402
from src import proteomes as P  # noqa: E402

PLOT_DIR = REPO_ROOT / "output" / "plots" / "ligands"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

LABEL = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}


def plot_buckets(ax, frames: dict[str, pd.DataFrame]) -> None:
    """Proteins with a potent ligand, per bucket, as a share of the proteome."""
    nc = stylia.NamedColors()
    buckets = ["species", "close", "remote"]
    colors = [nc.plum, nc.purple, nc.blue]
    width = 0.26
    xs = np.arange(len(frames))
    for i, (b, color) in enumerate(zip(buckets, colors)):
        vals = [100 * (f[f"{b}_n_compounds"] > 0).sum() / len(f) for f in frames.values()]
        ax.bar(xs + (i - 1) * width, vals, width, color=color, label=b)
    ax.set_xticks(xs)
    ax.set_xticklabels([LABEL[s] for s in frames], rotation=20, ha="right")
    ax.legend()
    stylia.label(ax, xlabel="", ylabel="% of proteome", title="Evidence by distance", abc="A")


def plot_cutoff(ax) -> None:
    """The pChEMBL sweep. Reported, never optimised."""
    nc = stylia.NamedColors()
    sens = L.load_cutoff_sensitivity()
    for sp, color in zip(L.SPECIES, [nc.plum, nc.mint, nc.orange]):
        s = sens[sens.species == sp].sort_values("pchembl")
        ax.plot(s.pchembl, s.n_remote, marker="o", color=color, label=LABEL[sp])
    ax.legend()
    stylia.label(ax, xlabel="pChEMBL cutoff", ylabel="proteins with evidence",
                 title="Cutoff sensitivity", abc="B")


def plot_identity(ax) -> None:
    """Where in identity space the evidence actually comes from."""
    nc = stylia.NamedColors()
    hits = L.load_hits()
    tgt = L.load_targets()
    lig = L.load_ligands()
    potent = set(lig.loc[lig.pchembl >= L.PCHEMBL_HEADLINE, "tid"])
    bact = tgt[(tgt.superkingdom == "Bacteria") & tgt.target_type.isin(L.SINGLE_TYPES)]
    j = hits.merge(bact[["component_id", "tid"]], on="component_id")
    j = j[(j.qcov >= L.MIN_QCOV) & (j.scov >= L.MIN_SCOV) & j.tid.isin(potent)]
    bins = np.arange(40, 101, 5)
    ax.hist(j.pident.clip(40, 100), bins=bins, color=nc.plum)
    for cut, color in ((L.CLOSE_PIDENT, nc.orange), (L.DIRECT_PIDENT, nc.mint)):
        ax.axvline(cut, color=color, linestyle="--")
    stylia.label(ax, xlabel="% identity to liganded target", ylabel="hits",
                 title="Where evidence comes from", abc="C")


def plot_compounds_vs_scaffolds(ax, all_df: pd.DataFrame) -> None:
    nc = stylia.NamedColors()
    d = all_df[all_df.remote_n_compounds > 0]
    ax.scatter(d.remote_n_compounds, d.remote_n_scaffolds, color=nc.plum, alpha=0.6)
    lim = [1, max(d.remote_n_compounds.max(), 1) * 1.2]
    ax.plot(lim, lim, color=nc.gray, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    stylia.label(ax, xlabel="compounds", ylabel="scaffolds",
                 title="Analog series collapse", abc="A")


def plot_ratio(ax, all_df: pd.DataFrame) -> None:
    nc = stylia.NamedColors()
    d = all_df[all_df.remote_n_compounds >= 5].copy()
    d["ratio"] = d.remote_n_compounds / d.remote_n_scaffolds.clip(lower=1)
    ax.hist(d.ratio, bins=30, color=nc.purple)
    stylia.label(ax, xlabel="compounds per scaffold", ylabel="proteins",
                 title="Redundancy per protein", abc="B")


def plot_selectivity(ax, all_df: pd.DataFrame) -> None:
    """Bacterial evidence against human evidence -- higher and to the right is worse."""
    nc = stylia.NamedColors()
    d = all_df[(all_df.remote_n_scaffolds > 0) | (all_df.human_n_scaffolds > 0)].copy()
    panel = d.gene_name.isin(I.genes()) if "gene_name" in d else pd.Series(False, index=d.index)
    # one generator for both axes: re-seeding per axis gives x and y identical offsets and draws
    # spurious diagonal streaks through the cloud
    rng = np.random.default_rng(0)
    bx = np.log10(d.remote_n_scaffolds + 1)
    by = np.log10(d.human_n_scaffolds + 1)
    x = bx + rng.uniform(-0.04, 0.04, len(d))
    y = by + rng.uniform(-0.04, 0.04, len(d))
    ax.scatter(x[~panel], y[~panel], color=nc.gray, alpha=0.4)
    ax.scatter(x[panel], y[panel], color=nc.plum)
    lim = [0, max(bx.max(), by.max()) * 1.05]
    ax.plot(lim, lim, color=nc.gray, linestyle="--")
    seen: list[tuple[float, float]] = []
    for _, r in d[d.gene_name.isin(["clpP", "gyrB", "lpxC", "folA", "ampC"])].iterrows():
        px, py = np.log10(r.remote_n_scaffolds + 1), np.log10(r.human_n_scaffolds + 1)
        if any(abs(px - sx) < 0.25 and abs(py - sy) < 0.12 for sx, sy in seen):
            continue
        seen.append((px, py))
        ax.annotate(r.gene_name, (px, py), xytext=(4, 5), textcoords="offset points",
                    fontsize=stylia.FONTSIZE_SMALL)
    stylia.label(ax, xlabel="log10 bacterial scaffolds + 1", ylabel="log10 human scaffolds + 1",
                 title="Selectivity liability")


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    frames = {sp: L.load(sp) for sp in L.SPECIES}
    named = []
    for sp, f in frames.items():
        m = f.merge(P.load(sp)[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left")
        m.insert(0, "species", sp)
        named.append(m)
    all_df = pd.concat(named, ignore_index=True)

    fig, axs = stylia.create_figure(1, 3)
    plot_buckets(axs.next(), frames)
    plot_cutoff(axs.next())
    plot_identity(axs.next())
    stylia.save_figure(str(PLOT_DIR / "chembl_coverage.png"))
    print(f"wrote {PLOT_DIR / 'chembl_coverage.png'}")

    fig, axs = stylia.create_figure(1, 2)
    plot_compounds_vs_scaffolds(axs.next(), all_df)
    plot_ratio(axs.next(), all_df)
    stylia.save_figure(str(PLOT_DIR / "chembl_redundancy.png"))
    print(f"wrote {PLOT_DIR / 'chembl_redundancy.png'}")

    fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
    plot_selectivity(axs.next(), all_df)
    stylia.save_figure(str(PLOT_DIR / "chembl_selectivity.png"))
    print(f"wrote {PLOT_DIR / 'chembl_selectivity.png'}")


if __name__ == "__main__":
    main()
