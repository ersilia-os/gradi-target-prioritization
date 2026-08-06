"""2x3 overview: everything worth knowing about where this proteome sits (docs §5.1).

Six panels, chosen to be *informative about localization itself* rather than about how well the
axis was built — no predictor concordance, no evidence tiers, no coverage accounting. Those live in
09h/09j. This is the slide for someone who wants to understand the proteome's architecture.

  1  the envelope, inward to outward  compartment sizes with their Clp accessibility
  2  the protein universe             ESM-C map coloured by compartment, all classes at once
  3  Clp accessibility ladder         the graded score and what each level means
  4  what each compartment is made of β-barrel / TM helix / signal peptide / soluble
  5  which biology lives where        functional class x compartment
  6  the membrane proteome            TM helices vs how much of the chain faces the cytoplasm

Reads output/results/<org>/<prefix>_localization.csv, the ESM-C projection, and — for panel 5 only —
`functional_class` from app/data/<prefix>.json, which is where 08a materialises it. Panel 5 degrades
to a message if that file is absent; nothing else depends on it.

Output: output/plots/09l_overview_<prefix>.png. Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
BG = "#D8D8D6"
GREY = "#C9C9C7"

# Inward -> outward through the envelope, which is also the order Clp accessibility falls in.
ORDER = ["cytoplasm", "inner_membrane", "periplasm", "outer_membrane",
         "extracellular", "cell_wall_surface"]
NICE = {"cytoplasm": "Cytoplasm", "inner_membrane": "Inner membrane", "periplasm": "Periplasm",
        "outer_membrane": "Outer membrane", "extracellular": "Extracellular",
        "cell_wall_surface": "Cell wall & surface"}

FC_NICE = {
    "ribosomal_translation": "Ribosome / translation", "dna_replication_repair": "DNA replication",
    "transcription_regulation": "Transcription", "signaling": "Signalling",
    "transport": "Transport", "cell_envelope": "Cell envelope",
    "oxidoreductase": "Oxidoreductase", "transferase": "Transferase",
    "hydrolase_protease": "Hydrolase / protease", "lyase_isomerase_ligase": "Lyase / ligase",
    "uncharacterized": "Uncharacterized", "other": "Other",
}

LADDER = [
    (1.0, "Cytoplasm", "freely reachable"),
    (0.6, "IM, cytoplasm-facing", "≥30% of the chain inside"),
    (0.2, "IM buried / periplasm", "across or inside the bilayer"),
    (0.0, "OM, β-barrel, surface", "out of reach"),
]


def load(org: str) -> pd.DataFrame:
    _, prefix = LOC.ORGANISMS[org]
    r = LOC.results_dir(org)
    d = pd.read_csv(r / f"{prefix}_localization.csv")
    proj = r / f"{prefix}_esmc600m_projection.csv"
    if proj.exists():
        d = d.merge(pd.read_csv(proj)[["uniprot_accession", "tsne_x", "tsne_y"]],
                    on="uniprot_accession", how="left")
    for c in ("tsne_x", "tsne_y"):
        if c not in d:
            d[c] = np.nan
    # functional_class only exists in the webapp payload
    app = LOC.REPO_ROOT / "app" / "data" / f"{prefix}.json"
    if app.exists():
        j = json.load(open(app))
        fc = pd.DataFrame(j["rows"], columns=j["columns"])
        if "functional_class" in fc.columns:
            d = d.merge(fc[["uniprot_accession", "functional_class"]],
                        on="uniprot_accession", how="left")
    if "functional_class" not in d:
        d["functional_class"] = None
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]

    d = load(org)
    n = len(d)
    present = [c for c in ORDER if (d["localization"] == c).sum() > 0]

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: the envelope, inward to outward ----
    ax = axs.next()
    counts = [int((d["localization"] == c).sum()) for c in present]
    clp = [d.loc[d["localization"] == c, "clp_accessibility"].mean() for c in present]
    y = np.arange(len(present))
    bars = ax.barh(y, counts, color=[LOC.LOC_CLASS_COLOR[c] for c in present], height=0.68)
    for i, (cnt, a) in enumerate(zip(counts, clp)):
        ax.text(cnt + max(counts) * 0.015, i,
                f"{cnt:,}  ({cnt / n:.0%})   Clp {a:.2f}",
                va="center", ha="left", fontsize=SS, color="#2B2333")
    ax.set_yticks(y)
    ax.set_yticklabels([NICE[c] for c in present], fontsize=SS)
    ax.invert_yaxis()
    ax.margins(x=0.42)
    stylia.label(ax, xlabel="proteins", ylabel="",
                 title=f"The envelope, inward → outward — {orgname}")

    # ---- panel 2: the protein universe ----
    ax = axs.next()
    has_xy = d["tsne_x"].notna() & d["tsne_y"].notna()
    ax.scatter(d.loc[has_xy, "tsne_x"], d.loc[has_xy, "tsne_y"], s=3.5, alpha=0.25, color=BG,
               linewidths=0, rasterized=True)
    for c in present:
        m = has_xy & (d["localization"] == c)
        ax.scatter(d.loc[m, "tsne_x"], d.loc[m, "tsne_y"], s=6, alpha=0.8,
                   color=LOC.LOC_CLASS_COLOR[c], linewidths=0, rasterized=True,
                   label=f"{LOC.LOC_ABBREV[c]} ({int(m.sum()):,})")
    # The map is full edge to edge, so the legend needs its own ground.
    ax.legend(fontsize=SS - 1, loc="upper left", ncol=2, handletextpad=0.3, columnspacing=0.8,
              borderpad=0.35, frameon=True, facecolor="white", framealpha=0.9, edgecolor="none")
    ax.set_xticks([]); ax.set_yticks([])
    stylia.label(ax, xlabel="ESM-C tSNE-1", ylabel="ESM-C tSNE-2",
                 title=f"Compartments in sequence space — {orgname}")

    # ---- panel 3: the Clp accessibility ladder ----
    ax = axs.next()
    vals = [int((d["clp_accessibility"] == lv).sum()) for lv, _, _ in LADDER]
    y = np.arange(len(LADDER))
    cmap = stylia.FadingColormap("turquoise").cmap
    bars = ax.barh(y, vals, color=[cmap(lv) for lv, _, _ in LADDER], height=0.66)
    for i, ((lv, lab, sub), v) in enumerate(zip(LADDER, vals)):
        ax.text(v + max(vals) * 0.015, i, f"{v:,}", va="center", ha="left", fontsize=SS,
                color="#2B2333")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{lv:g}   {lab}\n        {sub}" for lv, lab, sub in LADDER],
                       fontsize=SS - 1)
    ax.invert_yaxis()
    ax.margins(x=0.22)
    n_reach = int((d["clp_accessibility"] >= LOC.MIN_CLP_ACCESSIBILITY).sum())
    ax.text(0.97, 0.05, f"{n_reach:,} reachable ({n_reach / n:.0%})", transform=ax.transAxes,
            fontsize=SS, va="bottom", ha="right", color="#555555")
    stylia.label(ax, xlabel="proteins", ylabel="",
                 title=f"Can ClpXP reach it? — {orgname}")

    # ---- panel 4: what each compartment is made of ----
    ax = axs.next()
    beta = d["is_beta_barrel"].fillna(False).astype(bool)
    helix = (~beta) & (d["n_tm_helix"].fillna(0) > 0)
    sp = (~beta) & (~helix) & (d["has_signal_peptide_tmbed"].fillna(False).astype(bool))
    sol = (~beta) & (~helix) & (~sp)
    feats = [("β-barrel", "#8491B4", beta), ("TM helix", "#F39B7F", helix),
             ("signal peptide", "#FCBF49", sp), ("soluble", GREY, sol)]
    bottom = np.zeros(len(present))
    for name, colour, mask in feats:
        vals = np.array([100 * (mask & (d["localization"] == c)).sum()
                         / max((d["localization"] == c).sum(), 1) for c in present])
        ax.bar(range(len(present)), vals, bottom=bottom, color=colour, label=name, width=0.7)
        bottom += vals
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in present], fontsize=SS)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=SS - 1, ncol=2, loc="lower center", frameon=True, facecolor="white",
              framealpha=0.93, edgecolor="none")
    stylia.label(ax, xlabel="", ylabel="% of compartment",
                 title=f"What each compartment is made of — {orgname}")

    # ---- panel 5: which biology lives where ----
    ax = axs.next()
    if d["functional_class"].notna().any():
        classes = [c for c in FC_NICE if (d["functional_class"] == c).sum() >= 10]
        mat = np.zeros((len(classes), len(present)))
        for i, fc in enumerate(classes):
            sub = d[d["functional_class"] == fc]
            for j, c in enumerate(present):
                mat[i, j] = 100 * (sub["localization"] == c).sum() / max(len(sub), 1)
        ax.imshow(mat, cmap="Greens", aspect="auto", vmin=0, vmax=100)
        for i in range(len(classes)):
            for j in range(len(present)):
                if mat[i, j] >= 1:
                    ax.text(j, i, f"{mat[i, j]:.0f}", ha="center", va="center", fontsize=SS - 2,
                            color="white" if mat[i, j] > 55 else "#222222")
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in present], fontsize=SS)
        ax.set_yticks(range(len(classes)))
        ax.set_yticklabels([f"{FC_NICE[c]} ({int((d['functional_class'] == c).sum()):,})"
                            for c in classes], fontsize=SS - 2)
    else:
        ax.text(0.5, 0.5, "functional_class unavailable\n(run 08a_webapp_export.py)",
                ha="center", va="center", transform=ax.transAxes, fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="% of each functional class", ylabel="",
                 title=f"Which biology lives where — {orgname}")

    # ---- panel 6: the membrane proteome ----
    ax = axs.next()
    mem = d[(d["n_tm_helix"].fillna(0) > 0) | beta].copy()
    rng = np.random.default_rng(0)
    jitter = rng.uniform(-0.28, 0.28, len(mem))
    for c in present:
        m = (mem["localization"] == c) & (~mem["is_beta_barrel"].fillna(False).astype(bool))
        if m.any():
            ax.scatter(mem.loc[m, "n_tm_helix"] + jitter[m.to_numpy()],
                       mem.loc[m, "cyto_residue_fraction"], s=11, alpha=0.65,
                       color=LOC.LOC_CLASS_COLOR[c], linewidths=0, rasterized=True)
    b = mem["is_beta_barrel"].fillna(False).astype(bool)
    if b.any():
        ax.scatter(mem.loc[b, "n_tm_helix"] + jitter[b.to_numpy()],
                   mem.loc[b, "cyto_residue_fraction"], s=30, marker="*", color="#8491B4",
                   linewidths=0, rasterized=True, label=f"β-barrel ({int(b.sum())})")
    ax.axhline(LOC.CYTO_DOMAIN_FRACTION_MIN, ls="--", color="#555555", lw=1.1)
    ax.text(0.98, LOC.CYTO_DOMAIN_FRACTION_MIN + 0.02, "≥ 0.30 → Clp-reachable (0.6)",
            transform=ax.get_yaxis_transform(), fontsize=SS, va="bottom", ha="right",
            color="#555555")
    if b.any():
        ax.legend(fontsize=SS, frameon=False, loc="upper right")
    ax.set_ylim(-0.03, 1.03)
    # Bottom-left is where the β-barrels pile up (no helices, little cytoplasmic face), so the
    # count goes bottom-right instead.
    ax.text(0.98, 0.03, f"{len(mem):,} membrane proteins", transform=ax.transAxes, fontsize=SS,
            va="bottom", ha="right", color="#555555")
    stylia.label(ax, xlabel="transmembrane α-helices", ylabel="fraction facing the cytoplasm",
                 title=f"The membrane proteome — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"09l_overview_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")


if __name__ == "__main__":
    main()
