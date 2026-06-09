"""Comparative experimental-PDB analysis plots: K. pneumoniae vs E. coli (docs §1.2a).

Lightweight reader over the per-protein PDB-coverage tables produced by 04c -- no network.
Draws a 2x3 slide figure comparing how much experimental structural data each proteome has:

  1. PDB availability (% with a structure)   2. sequence coverage ECDF   3. apo vs holo
  4. best resolution distribution             5. experimental method     6. structures per protein

K. pneumoniae = NPG red, E. coli = NPG blue (matching 04b/01d/02c). Panels 2-6 use the subset
of proteins that have at least one PDB structure. Output:
  output/plots/04d_pdb_plots.png
Styling via stylia (ersilia-os/stylia), default slide format. Run with the `gradi` env.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_PATH = REPO_ROOT / "output" / "plots" / "04d_pdb_plots.png"

NPG = {
    "kp": "#E64B35",
    "ec": "#4DBBD5",
}  # one NPG colour per organism (matches 04b/01d/02c)
ORGANISMS = [
    ("kpneumoniae", "kp", "K. pneumoniae HS11286"),
    ("ecoli", "ec", "E. coli K-12 MG1655"),
]


def short_method(m: str) -> str:
    m = (m or "").lower()
    if "x-ray" in m or "xray" in m:
        return "X-ray"
    if "microscop" in m or "em" == m or "cryo" in m:
        return "EM"
    if "nmr" in m:
        return "NMR"
    return "other"


def load() -> dict[str, pd.DataFrame]:
    out = {}
    for org, prefix, _ in ORGANISMS:
        out[prefix] = pd.read_csv(
            REPO_ROOT / "output" / "results" / org / f"{prefix}_pdb_coverage.csv"
        )
    return out


def main() -> None:
    data = load()
    pdb = {
        p: df[df["pdb_has_structure"]] for p, df in data.items()
    }  # PDB-having subset

    stylia.set_format("slide")
    # full slide width (default); double the per-row slide height ratio for the 2-row grid
    fig, axs = stylia.create_figure(2, 3, height=0.6)

    # --- 1. PDB availability (% of proteome with any structure) ---
    ax = axs.next()
    xs = list(range(len(ORGANISMS)))
    pct = [100 * data[p]["pdb_has_structure"].mean() for _, p, _ in ORGANISMS]
    ax.bar(xs, pct, color=[NPG[p] for _, p, _ in ORGANISMS])
    ax.set_xticks(xs)
    ax.set_xticklabels([p for _, p, _ in ORGANISMS])
    ax.set_xlabel("")
    ax.set_ylabel("% proteome with PDB")
    ax.set_title("PDB availability")
    for x, (_, p, _) in zip(xs, ORGANISMS):
        n = int(data[p]["pdb_has_structure"].sum())
        tot = len(data[p])
        ax.text(x, pct[x], f"{n}/{tot}", ha="center", va="bottom")

    # --- 2. sequence coverage ECDF (union, PDB-having subset) ---
    ax = axs.next()
    for _, p, name in ORGANISMS:
        v = np.sort(pdb[p]["pdb_coverage_union"].values)
        y = np.arange(1, len(v) + 1) / len(v)
        ax.step(v, y, where="post", color=NPG[p], label=name)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("union sequence coverage")
    ax.set_ylabel("cumulative fraction")
    ax.set_title("Coverage (PDB subset)")
    ax.legend()

    # --- 3. apo vs holo (% of PDB-having proteins) ---
    ax = axs.next()
    xs = np.arange(len(ORGANISMS))
    holo = np.array([100 * pdb[p]["pdb_has_holo"].mean() for _, p, _ in ORGANISMS])
    ax.bar(xs, holo, color="#3C5488", label="holo (ligand-bound)", width=0.6)
    ax.bar(xs, 100 - holo, bottom=holo, color="#B0BEC5", label="apo", width=0.6)
    ax.set_xticks(xs)
    ax.set_xticklabels([p for _, p, _ in ORGANISMS])
    ax.set_xlabel("")
    ax.set_ylabel("% of PDB-having")
    ax.set_title("Apo vs holo")
    ax.legend()

    # --- 4. best resolution distribution (PDB-having subset) ---
    ax = axs.next()
    bins = np.linspace(1, 5, 33)
    for _, p, name in ORGANISMS:
        v = pdb[p]["pdb_best_resolution_A"].dropna()
        ax.hist(v, bins=bins, density=True, alpha=0.55, color=NPG[p], label=name)
    ax.set_xlabel("best resolution (Å)")
    ax.set_ylabel("density")
    ax.set_title("Resolution")

    # --- 5. experimental method composition (% of PDB-having) ---
    ax = axs.next()
    cats = ["X-ray", "EM", "NMR", "other"]
    xs = np.arange(len(cats))
    w = 0.38
    for k, (_, p, name) in enumerate(ORGANISMS):
        meth = pdb[p]["pdb_best_method"].map(short_method)
        frac = [100 * (meth == c).mean() for c in cats]
        ax.bar(xs + (k - 0.5) * w, frac, width=w, color=NPG[p], label=name)
    ax.set_xticks(xs)
    ax.set_xticklabels(cats)
    ax.set_xlabel("")
    ax.set_ylabel("% of PDB-having")
    ax.set_title("Experimental method")

    # --- 6. structures per protein (% of PDB-having) ---
    ax = axs.next()
    cats = ["1", "2–5", "6–20", "21–100", ">100"]

    def binned(n: pd.Series) -> list[float]:
        b = [(n == 1), n.between(2, 5), n.between(6, 20), n.between(21, 100), (n > 100)]
        return [100 * x.mean() for x in b]

    xs = np.arange(len(cats))
    w = 0.38
    for k, (_, p, name) in enumerate(ORGANISMS):
        ax.bar(
            xs + (k - 0.5) * w,
            binned(pdb[p]["pdb_n_structures"]),
            width=w,
            color=NPG[p],
            label=name,
        )
    ax.set_xticks(xs)
    ax.set_xticklabels(cats)
    ax.set_xlabel("PDB entries per protein")
    ax.set_ylabel("% of PDB-having")
    ax.set_title("Structures per protein")

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(PLOT_PATH))

    for _, p, name in ORGANISMS:
        df, s = data[p], pdb[p]
        print(
            f"{name}: PDB {int(df.pdb_has_structure.sum())}/{len(df)} "
            f"({100 * df.pdb_has_structure.mean():.1f}%), holo {int(s.pdb_has_holo.sum())} "
            f"({100 * s.pdb_has_holo.mean():.0f}% of PDB-having), "
            f"median coverage {s.pdb_coverage_union.median():.2f}, "
            f"median resolution {s.pdb_best_resolution_A.median():.2f} Å"
        )
    print(f"\nWrote {PLOT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
