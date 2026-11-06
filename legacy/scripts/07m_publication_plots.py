"""2x3 slide of PUBLICATION / experimental essentiality — no predictions (docs §4.1/§4.2).

Organism-symmetric: each tab shows THAT organism's real published screens (no placeholders).
  * K. pneumoniae — Mobile-CRISPRi-seq (Jana 2023), KPNIH1 & ECL8 Tn-seq, cross-species compendium.
  * E. coli K-12 — Keio KO, Goodall TraDIS, Rousset 2018 & Wang 2018 CRISPRi, cross-species compendium.

Panels (both organisms):
  1  experimental screen coverage        # proteins essential in each published screen
  2  CRISPRi depletion / in-vivo hits     the headline CRISPRi signal for that organism
  3  screen agreement                     # independent screens calling each gene essential
  4  cross-species conservation           # of 12 Enterobacteriaceae genomes essential in
  5  core essential-genome heatmap        genes x the 12 experimental genomes
  6  method concordance / library         KO vs TraDIS vs CRISPRi overlap (ec) / CRISPRi-by-function (kp)

Reads output/results/<org>/<prefix>_ess_publications.csv (07l), <prefix>_ess_experimental.csv (07n, ec),
and the Jana s0001 sheet (kp). Output: output/plots/07m_publications_<prefix>.png. Run with `gradi`.
"""

from __future__ import annotations

import argparse
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
from src import essentiality as E  # noqa: E402

REPO_ROOT = E.REPO_ROOT
NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
CRISPRI_C = "#E64B35"
JANA = E.essentiality_raw_dir("kpneumoniae", "jana2023_crispri") / "aem.00956-23-s0001.xlsx"
MB = E.essentiality_raw_dir("kpneumoniae", "mikebachman2023_KPPR1") / "ppat.1011233.s011.xlsx"

# compact genome tick labels for the cross-species heatmap (the full names are too long rotated)
GENOME_ABBR = {
    "K. pneumoniae ECL8": "Kp ECL8", "K. pneumoniae RH201207": "Kp RH201207",
    "E. coli BW25113": "Ec BW25113", "E. coli EC958": "Ec EC958",
    "E. coli NCTC13441": "Ec NCTC13441", "C. rodentium ICC168": "Cr ICC168",
    "S. Typhi Ty2": "S.Ty Ty2", "S. Tm A130": "S.Tm A130", "S. Tm D23580": "S.Tm D23580",
    "S. Tm SL3261": "S.Tm SL3261", "S. Tm SL1344": "S.Tm SL1344",
    "S. Enteritidis P125109": "S.Ent P125109",
}


def _note(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", color="#999", fontsize=SS)
    ax.set_xticks([]); ax.set_yticks([])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    orgname = E.ORG_DISPLAY[org]
    d = pd.read_csv(E.results_dir(org) / f"{prefix}_ess_publications.csv")
    genomes = [c.replace("pub_ess__", "") for c in d.columns if c.startswith("pub_ess__")]
    exp = None
    if org == "ecoli":
        p = E.results_dir(org) / "ec_ess_experimental.csv"
        if p.exists():
            exp = pd.read_csv(p)

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    def cnt(col, val=True):
        return int((d[col] == val).sum()) if col in d else 0

    # ---- panel 1: experimental screen coverage ----
    ax = axs.next()
    if org == "kpneumoniae":
        items = [("CRISPRi\nlibrary", cnt("crispri_ce_library"), CRISPRI_C),
                 ("KPNIH1\nTn-seq", cnt("kpnih1_essential"), NPG[0]),
                 ("ECL8\nTn-seq", cnt("ecl8_essential"), NPG[1]),
                 ("cross-sp.\ncore", cnt("pub_core_essential"), NPG[4])]
    else:
        items = [("Keio\nKO", cnt("ecoli_keio_essential"), NPG[3]),
                 ("Goodall\nTraDIS", cnt("ecoli_goodall_essential"), NPG[1]),
                 ("Rousset18\nCRISPRi", cnt("ecoli_crispri_rousset18_essential"), CRISPRI_C),
                 ("Wang18\nCRISPRi", cnt("ecoli_crispri_wang18_essential"), NPG[6] if len(NPG) > 6 else "#8491B4"),
                 ("cross-sp.\ncore", cnt("pub_core_essential"), NPG[4])]
    bars = ax.bar(range(len(items)), [v for _, v, _ in items], color=[c for _, _, c in items])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(len(items))); ax.set_xticklabels([a for a, _, _ in items], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="proteins (experimentally essential)",
                 title=f"Experimental screen coverage — {orgname}")
    ax.margins(y=0.2)

    # ---- panel 2: the headline CRISPRi signal ----
    ax = axs.next()
    if org == "kpneumoniae" and JANA.exists():
        iv = pd.read_excel(JANA, sheet_name="in vivo screening-KPPR1", header=1)
        iv["ratio"] = pd.to_numeric(iv["ratio"], errors="coerce")
        v2s = {}
        if MB.exists():
            mb = pd.read_excel(MB, sheet_name="S1 Table. TnSeq Results")
            v2s = {str(i).strip(): str(g).strip() for i, g in zip(mb["ID"], mb["Gene_name"]) if str(g).lower() != "nan"}
        def lbl(row):
            t = str(row.iloc[0]); t = t if t.startswith("VK055_") else f"VK055_{t.strip()}"
            return v2s.get(t, str(row["Gene"]))
        iv["glabel"] = iv.apply(lbl, axis=1)
        top = iv.sort_values("ratio", ascending=False).head(15).iloc[::-1]
        vals = np.log2(top["ratio"].clip(lower=1e-3)).to_numpy()
        ax.barh(range(len(top)), vals, color=CRISPRI_C)
        # gene/product labels INSIDE the bars (left-aligned, white) — no long y-tick labels
        ax.set_yticks([])
        for i, lab in enumerate(top["glabel"]):
            lab = lab if len(lab) <= 34 else lab[:32] + "…"
            ax.text(vals.max() * 0.015, i, lab, ha="left", va="center", fontsize=SS, color="white")
        ax.set_xlim(0, vals.max() * 1.02)
        stylia.label(ax, xlabel="log2 in-vivo depletion ratio", ylabel="",
                     title=f"Top CRISPRi in-vivo hits — {orgname}")
    elif org == "ecoli" and exp is not None:
        r = pd.to_numeric(exp["ecoli_crispri_rousset18_log2fc"], errors="coerce")
        w = pd.to_numeric(exp["ecoli_crispri_wang18_fitness"], errors="coerce")
        keio = exp["ecoli_keio_essential"] == True  # noqa: E712
        m = r.notna() & w.notna()
        ax.scatter(r[m & ~keio], w[m & ~keio], s=6, alpha=0.4, color="#C9C9C7", linewidths=0, rasterized=True, label="non-essential")
        ax.scatter(r[m & keio], w[m & keio], s=9, alpha=0.8, color=CRISPRI_C, linewidths=0, rasterized=True, label="Keio-essential")
        ax.axvline(-2, color="#999", ls=":", lw=1); ax.axhline(-6, color="#999", ls=":", lw=1)
        stylia.label(ax, xlabel="Rousset 2018 gene log2FC", ylabel="Wang 2018 gene fitness",
                     title=f"Genome-wide CRISPRi depletion — {orgname}")
        ax.legend(fontsize=SS, frameon=False, loc="upper left")
    else:
        _note(ax, "no CRISPRi data")
        stylia.label(ax, xlabel="", ylabel="", title=f"CRISPRi — {orgname}")

    # ---- panel 3: screen agreement ----
    ax = axs.next()
    if "n_experimental_screens" in d and (d["experimental_essential"] == True).any():  # noqa: E712
        ess = d[d["experimental_essential"] == True]  # noqa: E712
        maxn = 4
        vc = ess["n_experimental_screens"].clip(upper=maxn).value_counts().reindex(range(1, maxn + 1)).fillna(0).astype(int)
        bars = ax.bar(range(1, maxn + 1), vc.to_numpy(), color=[NPG[4], NPG[1], NPG[0], NPG[2]])
        ax.bar_label(bars, padding=2, fontsize=SS)
        ax.set_xticks(range(1, maxn + 1)); ax.set_xticklabels(["1", "2", "3", "4"])
        lab = "Kp screens" if org == "kpneumoniae" else "E. coli screens"
        stylia.label(ax, xlabel=f"# independent {lab}", ylabel="experimentally-essential genes",
                     title=f"Screen agreement — {orgname}")
        ax.margins(y=0.18)
    else:
        _note(ax, "screen agreement")
        stylia.label(ax, xlabel="", ylabel="", title=f"Screen agreement — {orgname}")

    # ---- panel 4: cross-species conservation (both) ----
    ax = axs.next()
    cov = d[d["pub_covered"] == True]  # noqa: E712
    ax.hist(cov["pub_n_species_essential"], bins=np.arange(-0.5, len(genomes) + 1.5, 1), color=NPG[3], rwidth=0.9)
    ax.axvline(len(genomes) * 0.8, color="#555", ls="--", lw=1)
    ax.set_yscale("log")
    stylia.label(ax, xlabel=f"# of {len(genomes)} genomes essential in", ylabel="genes",
                 title=f"Conservation of essentiality — {orgname}")

    # ---- panel 5: cross-species essentiality heatmap (both) ----
    # Keep tick labels at the deck-standard font size (SS); to do that we bound the number of rows to
    # what fits legibly (~15) rather than shrinking the font. Rows are sampled evenly across the
    # conservation gradient so the core (all-red) -> accessory (patchy) structure still reads.
    ax = axs.next()
    N_HEATMAP = 15
    ess_cols = [f"pub_ess__{g}" for g in genomes]
    var = cov[cov["pub_n_species_essential"] >= 1].copy()
    gene = E.load_genes(org)
    var["g"] = var["uniprot_accession"].map(gene).fillna(var["uniprot_accession"])
    var = var.sort_values("pub_n_species_essential", ascending=False)
    idx = (np.linspace(0, len(var) - 1, N_HEATMAP).round().astype(int)
           if len(var) > N_HEATMAP else range(len(var)))
    sel = var.iloc[idx]
    ax.imshow(sel[ess_cols].astype(float).to_numpy(), aspect="auto", cmap="Reds", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(genomes)))
    ax.set_xticklabels([GENOME_ABBR.get(g, g) for g in genomes], rotation=90, fontsize=SS)
    ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel["g"], fontsize=SS)
    # title states this is a representative slice of ALL genes essential in >=1 genome
    stylia.label(ax, xlabel="", ylabel="",
                 title=f"Cross-species essentiality — representative sample "
                       f"({len(sel)} of {len(var):,}) — {orgname}")

    # ---- panel 6: method concordance (ec) / CRISPRi-by-function (kp) ----
    ax = axs.next()
    if org == "ecoli" and exp is not None:
        ko = exp["ecoli_keio_essential"] == True  # noqa: E712
        tn = exp["ecoli_goodall_essential"] == True  # noqa: E712
        cr = (exp["ecoli_crispri_rousset18_essential"] == True) | (exp["ecoli_crispri_wang18_essential"] == True)  # noqa: E712
        cats = ["KO\n(Keio)", "TraDIS\n(Goodall)", "CRISPRi\n(R+W)", "all 3\nmethods"]
        vals = [int(ko.sum()), int(tn.sum()), int(cr.sum()), int((ko & tn & cr).sum())]
        bars = ax.bar(range(4), vals, color=[NPG[3], NPG[1], CRISPRI_C, NPG[2]])
        ax.bar_label(bars, padding=2, fontsize=SS)
        ax.set_xticks(range(4)); ax.set_xticklabels(cats, fontsize=SS)
        stylia.label(ax, xlabel="", ylabel="essential genes", title=f"Method concordance — {orgname}")
        ax.margins(y=0.18)
    elif org == "kpneumoniae" and JANA.exists():
        lib = pd.read_excel(JANA, sheet_name="870 selected essential genes", header=1)
        fn = lib["Function"].dropna().astype(str).str.strip()
        fn = fn[~fn.str.lower().isin(["", "nan", "-"])]
        vc = fn.value_counts().head(8).iloc[::-1]
        yy = np.arange(len(vc))
        # thin bars with the full COG category name placed ABOVE each bar (inline) — the categories are
        # too long for y-tick labels; count sits at the bar end.
        bars = ax.barh(yy, vc.to_numpy(), height=0.42, color=CRISPRI_C)
        ax.bar_label(bars, padding=3, fontsize=SS)
        ax.set_yticks([])
        for i, cat in enumerate(vc.index):
            ax.text(0, i + 0.28, cat, ha="left", va="bottom", fontsize=SS, color="#333")
        ax.set_ylim(-0.6, len(vc) - 0.1)
        ax.set_xlim(0, vc.max() * 1.15)
        stylia.label(ax, xlabel="genes in 870-gene CRISPRi library", ylabel="",
                     title=f"CRISPRi library by function — {orgname}")
    else:
        _note(ax, "method concordance")
        stylia.label(ax, xlabel="", ylabel="", title=f"Method concordance — {orgname}")

    out = REPO_ROOT / "output" / "plots" / f"07m_publications_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
