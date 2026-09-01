"""Slide figure for the ligandability assessment (06g output), docs §2.

One stylia npg slide, four panels, telling the ligandability story for presentations:
  (a) evidence sources  — # proteins with each independent ligandability signal (kp vs ec)
  (b) ligandability tiers — proteome composition by tier (kp vs ec)
  (c) pocket quality     — fpocket druggability vs pocket pLDDT (why we pLDDT-weight pockets)
  (d) prime targets      — tier breakdown of the whole proteome vs the broad-spectrum,
                           human-selective shortlist (the prioritization payoff)

Reads output/results/<org>/<prefix>_ligandability.csv + the 03c selectivity table.
Run with the `gradi` conda env. Output: output/plots/06h_ligandability.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import stylia  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

OUT_PATH = L.REPO_ROOT / "output" / "plots" / "06h_ligandability.png"
CATS = L.REPO_ROOT / "data" / "processed" / "other" / "orthology" / "three_way_protein_categories.tsv"

NPG = stylia.CategoricalPalette("npg").colors
KP, EC = NPG[0], NPG[1]
ORG_COLOR = {"kpneumoniae": KP, "ecoli": EC}
TIER_ORDER = ["intractable", "partial", "tractable"]
TIER_COLOR = {"intractable": "#C9C9C7", "partial": NPG[4], "tractable": NPG[2]}

EVIDENCE = [
    ("ChEMBL\n≤1µM", lambda d: d["chembl_any_n_potent"].fillna(0) > 0),
    ("BindingDB\n≤1µM", lambda d: d["bindingdb_any_n_potent"].fillna(0) > 0),
    ("PDB\nco-crystal", lambda d: d["pdb_lig_any_has_druglike"].fillna(False).astype(bool)),
    ("AlphaFill\ntransplant", lambda d: d["alphafill_n_druglike"].fillna(0) > 0),
    ("Druggable\npocket", lambda d: d["evidence_pocket"].fillna(0) >= 0.5),
]


def load(org: str) -> pd.DataFrame:
    _, prefix = L.ORGANISMS[org]
    return pd.read_csv(L.results_dir(org) / f"{prefix}_ligandability.csv")


def panel_evidence(ax, data: dict[str, pd.DataFrame]) -> None:
    labels = [e[0] for e in EVIDENCE]
    x = np.arange(len(EVIDENCE))
    w = 0.38
    for i, org in enumerate(["kpneumoniae", "ecoli"]):
        d = data.get(org)
        if d is None:
            continue
        counts = [int(fn(d).sum()) for _, fn in EVIDENCE]
        ax.bar(x + (i - 0.5) * w, counts, width=w, color=ORG_COLOR[org],
               label="K. pneumoniae" if org == "kpneumoniae" else "E. coli")
        for xi, c in zip(x + (i - 0.5) * w, counts):
            ax.text(xi, c, str(c), ha="center", va="bottom", fontsize=stylia.SLIDE_FONTSIZE_SMALL)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=stylia.SLIDE_FONTSIZE_SMALL)
    stylia.label(ax, xlabel="", ylabel="proteins", title="Ligandability evidence sources")
    ax.legend(fontsize=stylia.SLIDE_FONTSIZE_SMALL, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def panel_tiers(ax, data: dict[str, pd.DataFrame]) -> None:
    orgs = [o for o in ["kpneumoniae", "ecoli"] if o in data]
    x = np.arange(len(orgs))
    bottom = np.zeros(len(orgs))
    for tier in TIER_ORDER:
        fracs = []
        for org in orgs:
            d = data[org]
            fracs.append((d["ligandability_tier"] == tier).mean())
        ax.bar(x, fracs, bottom=bottom, width=0.6, color=TIER_COLOR[tier], label=tier)
        bottom += np.array(fracs)
    ax.set_xticks(x)
    ax.set_xticklabels(["K. pneumoniae", "E. coli"][: len(orgs)], fontsize=stylia.SLIDE_FONTSIZE_SMALL)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel="fraction of proteome", title="Ligandability tiers")
    ax.legend(fontsize=stylia.SLIDE_FONTSIZE_SMALL, frameon=False, loc="upper center", ncol=3)


def panel_pocket_quality(ax, data: dict[str, pd.DataFrame]) -> None:
    d = data["kpneumoniae"].copy()
    d = d[d["fpocket_max_drug_score"].notna() & d["fpocket_best_pocket_plddt"].notna()]
    if len(d) > 4000:
        d = d.sample(4000, random_state=0)
    colors = d["ligandability_tier"].map(TIER_COLOR).fillna("#C9C9C7")
    ax.scatter(d["fpocket_max_drug_score"], d["fpocket_best_pocket_plddt"],
               s=6, c=colors, alpha=0.45, linewidths=0)
    ax.axhline(70, color="#888888", lw=1.0, ls="--")
    ax.text(0.02, 71, "pLDDT 70", fontsize=stylia.SLIDE_FONTSIZE_SMALL, color="#666666")
    stylia.label(ax, xlabel="fpocket druggability score", ylabel="best-pocket mean pLDDT",
                 title="Pocket quality (K. pneumoniae)")


def panel_prime(ax, data: dict[str, pd.DataFrame]) -> None:
    d = data["kpneumoniae"]
    cats = pd.read_csv(CATS, sep="\t")
    prime = set(cats[(cats.organism == "kpneumoniae") & (cats.selectivity == "broad_selective")]["uniprot_accession"])
    groups = [("All Kp", d), ("Prime\n(broad-selective)", d[d["uniprot_accession"].isin(prime)])]
    x = np.arange(len(groups))
    bottom = np.zeros(len(groups))
    for tier in TIER_ORDER:
        fracs = [(g["ligandability_tier"] == tier).mean() if len(g) else 0 for _, g in groups]
        ax.bar(x, fracs, bottom=bottom, width=0.55, color=TIER_COLOR[tier], label=tier)
        bottom += np.array(fracs)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{n}\n(n={len(g)})" for n, g in groups], fontsize=stylia.SLIDE_FONTSIZE_SMALL)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel="fraction", title="Prime antibacterial targets")
    ax.legend(fontsize=stylia.SLIDE_FONTSIZE_SMALL, frameon=False, loc="upper center", ncol=3)


def main() -> None:
    data = {}
    for org in L.ORGANISMS:
        _, prefix = L.ORGANISMS[org]
        p = L.results_dir(org) / f"{prefix}_ligandability.csv"
        if p.exists():
            data[org] = load(org)
    if "kpneumoniae" not in data:
        raise SystemExit("kp ligandability table missing; run 06g first.")

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 2, width=1.0, height=0.66)
    panel_evidence(axs.next(), data)
    panel_tiers(axs.next(), data)
    panel_pocket_quality(axs.next(), data)
    panel_prime(axs.next(), data)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(OUT_PATH))
    print(f"Wrote {OUT_PATH.relative_to(L.REPO_ROOT)}")


if __name__ == "__main__":
    main()
