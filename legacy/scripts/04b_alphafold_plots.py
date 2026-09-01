"""Comparative AlphaFold analysis plots: K. pneumoniae vs E. coli (docs §1.2b).

Lightweight reader over the per-protein metrics produced by 04a -- no network, no structure
parsing. Draws a 2x3 slide figure comparing the two proteomes across structural-confidence and
domain-organisation axes:

  1. AFDB coverage (% modeled)        2. mean pLDDT distribution     3. pLDDT band composition
  4. domain-count distribution        5. inter-domain PAE (multidom) 6. disorder ECDF

K. pneumoniae = NPG red, E. coli = NPG blue (matching 01d/02c). Output:
  output/plots/04b_alphafold_plots.png
Styling via stylia (ersilia-os/stylia), default slide parameters. Run with the `gradi` env.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_PATH = REPO_ROOT / "output" / "plots" / "04b_alphafold_plots.png"

NPG = {
    "kp": "#E64B35",
    "ec": "#4DBBD5",
}  # one NPG colour per organism (matches 01d/02c)
ORGANISMS = [
    ("kpneumoniae", "kp", "K. pneumoniae HS11286"),
    ("ecoli", "ec", "E. coli K-12 MG1655"),
]
# pLDDT confidence bands: green (good) -> red (bad)
BANDS = [
    ("af_frac_very_high_plddt", "very high (≥90)", "#006d2c"),
    ("af_frac_confident_plddt", "confident (70–90)", "#74c476"),
    ("af_frac_low_plddt", "low (50–70)", "#fdae6b"),
    ("af_frac_very_low_plddt", "very low (<50)", "#d7301f"),
]


def load() -> dict[str, pd.DataFrame]:
    data = {}
    for org, prefix, _ in ORGANISMS:
        path = (
            REPO_ROOT / "output" / "results" / org / f"{prefix}_alphafold_structure.csv"
        )
        data[prefix] = pd.read_csv(path)
    return data


def main() -> None:
    data = load()
    avail = {p: df[df["af_available"]] for p, df in data.items()}

    stylia.set_format("slide")
    # slide width is full SIZE (default); slide's default height ratio (0.3) is per single row,
    # so double it for the two-row grid to keep each row at slide proportions.
    fig, axs = stylia.create_figure(2, 3, height=0.6)

    # --- 1. AFDB coverage (% of proteome modeled) ---
    ax = axs.next()
    xs = list(range(len(ORGANISMS)))
    pcts = [100 * data[p]["af_available"].mean() for _, p, _ in ORGANISMS]
    ax.bar(xs, pcts, color=[NPG[p] for _, p, _ in ORGANISMS])
    ax.set_xticks(xs)
    ax.set_xticklabels([p for _, p, _ in ORGANISMS])
    ax.set_xlabel("")
    ax.set_ylim(90, 100.5)
    ax.set_ylabel("% proteome modeled")
    ax.set_title("AFDB coverage")
    for x, (_, p, _) in zip(xs, ORGANISMS):
        n_miss = int((~data[p]["af_available"]).sum())
        ax.text(x, pcts[x], f"−{n_miss}", ha="center", va="bottom")

    # --- 2. mean pLDDT distribution ---
    ax = axs.next()
    bins = np.linspace(30, 100, 40)
    for _, p, name in ORGANISMS:
        ax.hist(
            avail[p]["af_mean_plddt"],
            bins=bins,
            density=True,
            alpha=0.55,
            color=NPG[p],
            label=name,
        )
    ax.set_xlabel("mean pLDDT")
    ax.set_ylabel("density")
    ax.set_title("Per-protein confidence")
    ax.legend()

    # --- 3. pLDDT band composition (proteome-mean residue fractions) ---
    ax = axs.next()
    xs = np.arange(len(ORGANISMS))
    bottom = np.zeros(len(ORGANISMS))
    for col, label, color in BANDS:
        vals = np.array([avail[p][col].mean() for _, p, _ in ORGANISMS])
        ax.bar(xs, vals, bottom=bottom, color=color, label=label, width=0.6)
        bottom += vals
    ax.set_xticks(xs)
    ax.set_xticklabels([p for _, p, _ in ORGANISMS])
    ax.set_xlabel("")
    ax.set_ylabel("mean residue fraction")
    ax.set_title("pLDDT band composition")
    ax.legend()

    # --- 4. disorder ECDF (fraction of pLDDT<50 residues) -- bottom-left ---
    ax = axs.next()
    for _, p, name in ORGANISMS:
        v = np.sort(avail[p]["af_frac_very_low_plddt"].values)
        y = np.arange(1, len(v) + 1) / len(v)
        ax.step(v, y, where="post", color=NPG[p], label=name)
    ax.set_xlim(0, 0.5)
    ax.set_xlabel("fraction pLDDT<50 (disorder)")
    ax.set_ylabel("cumulative fraction")
    ax.set_title("Disorder content")

    # --- 5. domain-count distribution (% of modeled) ---
    ax = axs.next()
    cats = ["1", "2", "3", "4", "5+"]
    xs = np.arange(len(cats))
    w = 0.38
    for k, (_, p, name) in enumerate(ORGANISMS):
        nd = avail[p]["af_n_domains"].clip(upper=5)
        frac = [
            100 * (nd == 1).mean(),
            100 * (nd == 2).mean(),
            100 * (nd == 3).mean(),
            100 * (nd == 4).mean(),
            100 * (nd == 5).mean(),
        ]
        ax.bar(xs + (k - 0.5) * w, frac, width=w, color=NPG[p], label=name)
    ax.set_xticks(xs)
    ax.set_xticklabels(cats)
    ax.set_xlabel("PAE domains")
    ax.set_ylabel("% of modeled")
    ax.set_title("Domain count")

    # --- 6. inter-domain PAE (multidomain subset) -- bottom-right ---
    ax = axs.next()
    bins = np.linspace(0, 32, 33)
    for _, p, name in ORGANISMS:
        vals = (
            avail[p]
            .loc[avail[p]["af_is_multidomain"], "af_interdomain_pae_mean"]
            .dropna()
        )
        ax.hist(vals, bins=bins, density=True, alpha=0.55, color=NPG[p], label=name)
    ax.set_xlabel("mean inter-domain PAE (Å)")
    ax.set_ylabel("density")
    ax.set_title("Domain-arrangement uncertainty")

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(PLOT_PATH))

    # --- compact per-organism summary ---
    for _, p, name in ORGANISMS:
        df, a = data[p], avail[p]
        print(
            f"{name}: modeled {int(df.af_available.sum())}/{len(df)} "
            f"({100 * df.af_available.mean():.1f}%), median pLDDT {a.af_mean_plddt.median():.1f}, "
            f"high-conf(≥90) {100 * (a.af_mean_plddt >= 90).mean():.0f}%, "
            f"multidomain {100 * a.af_is_multidomain.mean():.0f}%, "
            f"median disorder {a.af_frac_very_low_plddt.median():.3f}"
        )
    print(f"\nWrote {PLOT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
