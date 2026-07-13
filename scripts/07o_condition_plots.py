"""Condition / stress-specific essentiality (docs §4.1 conditional) — organism-symmetric slide.

Essentiality is not absolute — genes can be required only under a stress. This slide surfaces the
condition-dependent axis, from the richest condition resource each organism has:
  * E. coli — RB-TnSeq / Fitness Browser (Price 2018): gene fitness across ~280 antibiotic/stress
    conditions. Shows antibiotic-conditional essential genes, the drug × gene sensitivity landscape,
    and constitutive- vs condition-specific vulnerability.
  * K. pneumoniae — the ECL8 niche screens (Eichelberger 2024: urine, serum) + KPPR1 in-vivo: which
    genes become required in a host-relevant niche.

Reads (ec) data/processed/ecoli/essentiality/conditions/ec_rbtnseq_conditions.csv and
output/results/ecoli/ec_ess_experimental.csv; (kp) output/results/kpneumoniae/kp_ess_kp.csv.
Output: output/plots/07o_conditions_<prefix>.png. Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

REPO_ROOT = E.REPO_ROOT
NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
DEFECT = -2.0  # RB-TnSeq fitness below this = strong sensitization / conditional requirement

# antibiotic (and key stress) classes -> regex over the condition label
DRUGS = {
    "ciprofloxacin": r"ciprofloxacin", "nalidixic": r"nalidixic", "norfloxacin": r"norfloxacin",
    "trimethoprim": r"trimethoprim", "nitrofurantoin": r"nitrofurantoin",
    "carbenicillin": r"carbenicillin", "ampicillin": r"ampicillin", "mecillinam": r"mecillinam",
    "cephalothin": r"cephalothin", "ceftriaxone": r"ceftriaxone", "cefoxitin": r"cefoxitin",
    "tetracycline": r"tetracycline|doxycycline|oxytetracycline", "chloramphenicol": r"chloramphenicol",
    "gentamicin": r"gentamicin", "spectinomycin": r"spectinomycin", "polymyxin": r"polymyxin",
    "bacitracin": r"bacitracin", "vancomycin": r"vancomycin", "novobiocin": r"novobiocin",
    "A22": r"\bA22\b",
}


def _cond_label(col):  # "set8IT031 LB Aerobic with Ciprofloxacin ... mM" -> the descriptive tail
    return re.sub(r"^set\w+?\s+", "", col)


def ecoli_panels(axs):
    orgname = E.ORG_DISPLAY["ecoli"]
    f = E.essentiality_processed_dir("ecoli", "conditions") / "ec_rbtnseq_conditions.csv"
    d = pd.read_csv(f)
    gene = E.load_genes("ecoli")
    d["g"] = d["uniprot_accession"].map(gene).fillna(d["uniprot_accession"])
    cond_cols = [c for c in d.columns if c not in ("uniprot_accession", "desc", "g")]
    F = d[cond_cols].apply(pd.to_numeric, errors="coerce")

    # drug -> its condition columns
    drug_cols = {drug: [c for c in cond_cols if re.search(pat, _cond_label(c), re.I)]
                 for drug, pat in DRUGS.items()}
    drug_cols = {k: v for k, v in drug_cols.items() if v}
    # per-gene min fitness within each drug (strongest sensitization)
    drug_min = pd.DataFrame({drug: F[cols].min(axis=1) for drug, cols in drug_cols.items()})

    # panel 1: # genes sensitized per antibiotic
    ax = axs.next()
    nsens = (drug_min < DEFECT).sum(axis=0).sort_values()
    bars = ax.barh(range(len(nsens)), nsens.to_numpy(), color=NPG[0])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_yticks(range(len(nsens))); ax.set_yticklabels(nsens.index, fontsize=SS)
    stylia.label(ax, xlabel=f"# genes sensitized (fitness < {DEFECT})", ylabel="",
                 title=f"Antibiotic-conditional essentials — {orgname}")
    ax.margins(x=0.12)

    # panel 2: drug × gene sensitivity heatmap (top sensitized genes)
    ax = axs.next()
    top_genes = drug_min.min(axis=1).nsmallest(28).index
    M = drug_min.loc[top_genes]
    order = nsens.index[::-1]
    M = M[order]
    im = ax.imshow(M.to_numpy(), aspect="auto", cmap="RdBu", vmin=-6, vmax=6, interpolation="nearest")
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=5)
    ax.set_yticks(range(len(top_genes))); ax.set_yticklabels(d.loc[top_genes, "g"], fontsize=5)
    stylia.label(ax, xlabel="", ylabel="", title=f"Drug × gene sensitivity (RB-TnSeq) — {orgname}")

    # panel 3: most condition-variable genes (widest fitness range)
    ax = axs.next()
    rng = (F.max(axis=1) - F.min(axis=1))
    top = rng.nlargest(15).index[::-1]
    ax.barh(range(len(top)), rng.loc[top].to_numpy(), color=NPG[4])
    ax.set_yticks(range(len(top))); ax.set_yticklabels(d.loc[top, "g"], fontsize=SS)
    stylia.label(ax, xlabel="fitness range across conditions", ylabel="",
                 title=f"Most condition-variable genes — {orgname}")

    # panel 4: constitutive vs condition-specific (min fitness vs # conditions with defect)
    ax = axs.next()
    minf = F.min(axis=1); ndef = (F < DEFECT).sum(axis=1)
    ax.scatter(ndef, minf, s=6, alpha=0.4, color="#C9C9C7", linewidths=0, rasterized=True)
    ax.set_xscale("symlog")
    stylia.label(ax, xlabel="# conditions with strong defect", ylabel="strongest fitness (min)",
                 title=f"Constitutive vs conditional — {orgname}")

    # panel 5: fluoroquinolone (ciprofloxacin) example — top sensitized genes
    ax = axs.next()
    cip = [c for c in cond_cols if re.search(r"ciprofloxacin", _cond_label(c), re.I)]
    if cip:
        cipmin = F[cip].min(axis=1)
        top = cipmin.nsmallest(12).index[::-1]
        ax.barh(range(len(top)), cipmin.loc[top].to_numpy(), color=NPG[3])
        ax.set_yticks(range(len(top))); ax.set_yticklabels(d.loc[top, "g"], fontsize=SS)
        stylia.label(ax, xlabel="ciprofloxacin fitness (min)", ylabel="",
                     title=f"Ciprofloxacin-sensitizing genes — {orgname}")
    # panel 6: trimethoprim example
    ax = axs.next()
    tmp = [c for c in cond_cols if re.search(r"trimethoprim", _cond_label(c), re.I)]
    if tmp:
        tmin = F[tmp].min(axis=1)
        top = tmin.nsmallest(12).index[::-1]
        ax.barh(range(len(top)), tmin.loc[top].to_numpy(), color=NPG[1])
        ax.set_yticks(range(len(top))); ax.set_yticklabels(d.loc[top, "g"], fontsize=SS)
        stylia.label(ax, xlabel="trimethoprim fitness (min)", ylabel="",
                     title=f"Trimethoprim-sensitizing genes — {orgname}")


def kp_panels(axs):
    orgname = E.ORG_DISPLAY["kpneumoniae"]
    kp = pd.read_csv(E.results_dir("kpneumoniae") / "kp_ess_kp.csv")
    gene = E.load_genes("kpneumoniae")
    kp["g"] = kp["uniprot_accession"].map(gene).fillna(kp["uniprot_accession"])
    urine = kp["kp_ess_urine_call"] == "conditional"
    serum = kp["kp_ess_serum_call"] == "conditional"
    invivo = kp["kp_ess_in_vivo_call"] == "in_vivo_defect"

    # panel 1: niche-required counts
    ax = axs.next()
    cats = ["urine\n(ECL8)", "serum\n(ECL8)", "in-vivo\n(KPPR1)"]
    vals = [int(urine.sum()), int(serum.sum()), int(invivo.sum())]
    bars = ax.bar(range(3), vals, color=[NPG[0], NPG[1], NPG[3]])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels(cats, fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="genes required in niche", title=f"Niche-conditional genes — {orgname}")
    ax.margins(y=0.18)

    # panel 2: urine vs serum score scatter
    ax = axs.next()
    us = pd.to_numeric(kp["kp_ess_urine_score"], errors="coerce").fillna(0)
    ss = pd.to_numeric(kp["kp_ess_serum_score"], errors="coerce").fillna(0)
    ax.scatter(us, ss, s=8, alpha=0.5, color="#C9C9C7", linewidths=0, rasterized=True)
    both = (us > 0) & (ss > 0)
    ax.scatter(us[both], ss[both], s=12, alpha=0.85, color=NPG[2], linewidths=0, label="required in both")
    stylia.label(ax, xlabel="urine requirement score", ylabel="serum requirement score",
                 title=f"Urine vs serum requirement — {orgname}")
    ax.legend(fontsize=SS, frameon=False, loc="upper right")

    # panels 3/4: top urine / serum required genes
    for call, score, title, color in [("kp_ess_urine_call", "kp_ess_urine_score", "urine", NPG[0]),
                                       ("kp_ess_serum_call", "kp_ess_serum_score", "serum", NPG[1])]:
        ax = axs.next()
        sub = kp[kp[call] == "conditional"].copy()
        sub["s"] = pd.to_numeric(sub[score], errors="coerce")
        top = sub.nlargest(12, "s").iloc[::-1]
        ax.barh(range(len(top)), top["s"], color=color)
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top["g"], fontsize=SS)
        stylia.label(ax, xlabel=f"{title} requirement score", ylabel="",
                     title=f"Top {title}-required genes — {orgname}")

    # panel 5: niche overlap
    ax = axs.next()
    only_u = int((urine & ~serum).sum()); only_s = int((serum & ~urine).sum()); both_n = int((urine & serum).sum())
    bars = ax.bar(range(3), [only_u, both_n, only_s], color=[NPG[0], NPG[2], NPG[1]])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels(["urine\nonly", "both", "serum\nonly"], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="genes", title=f"Niche overlap — {orgname}")
    ax.margins(y=0.18)

    # panel 6: genes required in BOTH niches (shared host-adaptation core)
    ax = axs.next()
    both_df = kp[urine & serum].copy()
    both_df["s"] = pd.to_numeric(both_df["kp_ess_urine_score"], errors="coerce").fillna(0) \
        + pd.to_numeric(both_df["kp_ess_serum_score"], errors="coerce").fillna(0)
    top = both_df.nlargest(12, "s").iloc[::-1]
    if len(top):
        ax.barh(range(len(top)), top["s"], color=NPG[2])
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top["g"], fontsize=SS)
    stylia.label(ax, xlabel="urine + serum requirement", ylabel="",
                 title=f"Required in both niches — {orgname}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)
    if org == "ecoli":
        ecoli_panels(axs)
    else:
        kp_panels(axs)

    out = REPO_ROOT / "output" / "plots" / f"07o_conditions_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
