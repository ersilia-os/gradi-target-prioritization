"""2x3 slide: the measured activated-ClpP evidence (docs §3.3c).

This is the only dataset in the repo that measures the Gr-ADI modality itself. ADEP4 (Conlon 2013)
and ONC212 (Jacques 2020) both open ClpP **without an unfoldase partner**, so a protein depleted or
cleaved under them is direct evidence for partnerless ClpP degradation — not a proxy, not a
sequence rule.

  1  do the two activators agree?   ADEP4 vs ONC212 abundance log2FC — independent chemistries
  2  do the two readouts agree?     abundance depletion vs cleavage-peptide evidence (ADEP4)
  3  volcano                        ADEP4 log2FC vs adjusted p, substrates named
  4  evidence composition           depleted / cleaved / both / neither, and the evidence score
  5  transfer quality               DIAMOND RBH identity vs coverage — the honesty panel
  6  named targets                  AcpP, DnaK, GyrA, GyrB against the measured distribution

Everything here was measured in *S. aureus* and mapped onto the anchor organism by reciprocal best
hit, reaching only ~11-14% of the proteome. Panel 5 exists so that caveat is never off-screen.

Reads output/results/<org>/<prefix>_clpp_activator.csv (10c).
Output: output/plots/10l_activator_<prefix>.png. Run with `gradi`.
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
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL

DEPLETED_C = "#E64B35"     # depleted under activator = evidence of degradation
CLEAVED_C = "#00A087"
BOTH_C = "#3C5488"
NONE_C = "#C9C9C7"

# Thresholds 10c scored with, redrawn here so the panels and the table agree.
ABUNDANCE_LOG2 = -1.0
CLEAVAGE_LOG2 = 1.0
PADJ = 0.05

# The proposal's declared targets. The docs pre-register their activator evidence, and note the
# proposal's own quantitative wording does not survive the data — so plot them as measured.
NAMED = ["acpP", "dnaK", "gyrA", "gyrB"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]

    d = pd.read_csv(LOC.results_dir(org) / f"{prefix}_clpp_activator.csv")
    anchor_genes = LOC.load_genes(org)
    d["anchor_gene"] = d["uniprot_accession"].map(anchor_genes).fillna("")
    n_proteome = len(LOC.load_accessions(org))

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: ADEP4 vs ONC212 abundance ----
    ax = axs.next()
    both = d.dropna(subset=["adep4_abundance_log2fc", "onc212_abundance_log2fc"])
    x, y = both["adep4_abundance_log2fc"], both["onc212_abundance_log2fc"]
    ax.axhline(0, color="#DDDDDD", lw=0.8); ax.axvline(0, color="#DDDDDD", lw=0.8)
    ax.axvline(ABUNDANCE_LOG2, ls="--", color="#999999", lw=0.9)
    ax.axhline(ABUNDANCE_LOG2, ls="--", color="#999999", lw=0.9)
    agree = (x <= ABUNDANCE_LOG2) & (y <= ABUNDANCE_LOG2)
    ax.scatter(x[~agree], y[~agree], s=9, alpha=0.5, color=NONE_C, linewidths=0, rasterized=True)
    ax.scatter(x[agree], y[agree], s=14, alpha=0.85, color=DEPLETED_C, linewidths=0, rasterized=True)
    if len(both) > 3:
        rho, p = stats.spearmanr(x, y)
        ax.text(0.03, 0.03, f"n = {len(both)}\nρ = {rho:.2f}\nboth depleted: {int(agree.sum())}",
                transform=ax.transAxes, fontsize=SS, va="bottom", ha="left", color="#2B2333",
                linespacing=1.4)
    stylia.label(ax, xlabel="ADEP4 abundance log2FC", ylabel="ONC212 abundance log2FC",
                 title=f"Two activators, one protease — {orgname}")

    # ---- panel 2: abundance vs cleavage (ADEP4) ----
    ax = axs.next()
    rd = d.dropna(subset=["adep4_abundance_log2fc", "adep4_cleavage_log2fc_max"])
    x, y = rd["adep4_abundance_log2fc"], rd["adep4_cleavage_log2fc_max"]
    ax.axhline(0, color="#DDDDDD", lw=0.8); ax.axvline(0, color="#DDDDDD", lw=0.8)
    ax.axvline(ABUNDANCE_LOG2, ls="--", color="#999999", lw=0.9)
    ax.axhline(CLEAVAGE_LOG2, ls="--", color="#999999", lw=0.9)
    hit = (x <= ABUNDANCE_LOG2) & (y >= CLEAVAGE_LOG2)
    ax.scatter(x[~hit], y[~hit], s=9, alpha=0.5, color=NONE_C, linewidths=0, rasterized=True)
    ax.scatter(x[hit], y[hit], s=16, alpha=0.9, color=BOTH_C, linewidths=0, rasterized=True)
    if len(rd) > 3:
        rho, _ = stats.spearmanr(x, y)
        ax.text(0.03, 0.97, f"n = {len(rd)}\nρ = {rho:.2f}\ndepleted + cleaved: {int(hit.sum())}",
                transform=ax.transAxes, fontsize=SS, va="top", ha="left", color="#2B2333",
                linespacing=1.4)
    stylia.label(ax, xlabel="ADEP4 abundance log2FC", ylabel="ADEP4 cleavage log2FC (max peptide)",
                 title=f"Two readouts, one activator — {orgname}")

    # ---- panel 3: volcano ----
    ax = axs.next()
    v = d.dropna(subset=["adep4_abundance_log2fc", "adep4_abundance_padj"])
    lp = -np.log10(v["adep4_abundance_padj"].clip(lower=1e-12))
    sig = (v["adep4_abundance_log2fc"] <= ABUNDANCE_LOG2) & (v["adep4_abundance_padj"] < PADJ)
    ax.scatter(v.loc[~sig, "adep4_abundance_log2fc"], lp[~sig], s=9, alpha=0.5, color=NONE_C,
               linewidths=0, rasterized=True)
    ax.scatter(v.loc[sig, "adep4_abundance_log2fc"], lp[sig], s=16, alpha=0.9, color=DEPLETED_C,
               linewidths=0, rasterized=True)
    ax.axvline(ABUNDANCE_LOG2, ls="--", color="#999999", lw=0.9)
    ax.axhline(-np.log10(PADJ), ls="--", color="#999999", lw=0.9)
    # Leader lines into a vertical stack, the 07k pattern — the most-depleted proteins cluster
    # tightly at the top-left and direct labels collide.
    top = v[sig].nsmallest(6, "adep4_abundance_log2fc")
    if len(top):
        y0, y1 = ax.get_ylim()
        xs = ax.get_xlim()
        lab_x = xs[0] + 0.34 * (xs[1] - xs[0])
        ys_lab = np.linspace(y0 + 0.42 * (y1 - y0), y0 + 0.94 * (y1 - y0), len(top))[::-1]
        for i, (_, r) in enumerate(top.iterrows()):
            lab = r["anchor_gene"] or r.get("gene") or r["uniprot_accession"]
            ax.annotate(str(lab),
                        xy=(r["adep4_abundance_log2fc"],
                            -np.log10(max(r["adep4_abundance_padj"], 1e-12))),
                        xytext=(lab_x, ys_lab[i]), fontsize=SS - 1, color="#2B2333",
                        va="center", ha="right",
                        arrowprops=dict(arrowstyle="-", color="#CCCCCC", lw=0.5))
    ax.text(0.97, 0.97, f"n = {len(v)}\n{int(sig.sum())} significant", transform=ax.transAxes,
            fontsize=SS, va="top", ha="right", color="#555555", linespacing=1.4)
    stylia.label(ax, xlabel="ADEP4 abundance log2FC", ylabel="−log10 adjusted p",
                 title=f"ADEP4 depletion — {orgname}")

    # ---- panel 4: evidence composition ----
    ax = axs.next()
    dep = d["activator_depleted"].fillna(False).astype(bool)
    cle = d["activator_cleaved"].fillna(False).astype(bool)
    cats = [("both", int((dep & cle).sum()), BOTH_C),
            ("cleaved only", int((~dep & cle).sum()), CLEAVED_C),
            ("depleted only", int((dep & ~cle).sum()), DEPLETED_C),
            ("neither", int((~dep & ~cle).sum()), NONE_C)]
    bars = ax.bar(range(len(cats)), [n for _, n, _ in cats], color=[c for _, _, c in cats], width=0.62)
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([lab for lab, _, _ in cats], fontsize=SS, rotation=12)
    ax.margins(y=0.2)
    ax.text(0.97, 0.97,
            f"{len(d):,} of {n_proteome:,} proteins\nreached ({len(d) / n_proteome:.1%})",
            transform=ax.transAxes, fontsize=SS, va="top", ha="right", color="#555555",
            linespacing=1.4)
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Activator evidence — {orgname}")

    # ---- panel 5: transfer quality ----
    ax = axs.next()
    t = d.dropna(subset=["pident", "coverage"])
    has_ev = t["activator_evidence"].fillna(0) > 0
    ax.scatter(t.loc[~has_ev, "pident"], t.loc[~has_ev, "coverage"], s=9, alpha=0.45, color=NONE_C,
               linewidths=0, rasterized=True, label="no evidence")
    ax.scatter(t.loc[has_ev, "pident"], t.loc[has_ev, "coverage"], s=13, alpha=0.8, color=BOTH_C,
               linewidths=0, rasterized=True, label="has evidence")
    ax.axvline(40, ls="--", color="#999999", lw=0.9)
    ax.legend(fontsize=SS, frameon=False, loc="lower left")
    ax.text(0.97, 0.03,
            f"median identity {t['pident'].median():.0f}%\n"
            f"S. aureus → {orgname.split()[0]}.\nreciprocal best hits",
            transform=ax.transAxes, fontsize=SS, va="bottom", ha="right", color="#555555",
            linespacing=1.4)
    stylia.label(ax, xlabel="DIAMOND identity (%)", ylabel="coverage (%)",
                 title=f"All of this is transferred — {orgname}")

    # ---- panel 6: named targets ----
    ax = axs.next()
    ev = d["activator_evidence"].dropna()
    ax.hist(ev, bins=np.linspace(0, 1, 26), color=NONE_C)
    ymax = ax.get_ylim()[1]
    found = []
    for i, g in enumerate(NAMED):
        r = d[d["anchor_gene"].str.lower() == g.lower()]
        if r.empty:
            continue
        val = float(r.iloc[0]["activator_evidence"])
        found.append((g, val))
        ax.axvline(val, color=NPG[i % len(NPG)], lw=1.6)
        ax.annotate(f"{g}  {val:.2f}", xy=(val, ymax * (0.92 - 0.12 * i)),
                    fontsize=SS, color=NPG[i % len(NPG)], ha="left" if val < 0.6 else "right",
                    xytext=(4 if val < 0.6 else -4, 0), textcoords="offset points")
    if not found:
        ax.text(0.5, 0.5, "named targets not in the transferred set", ha="center", va="center",
                transform=ax.transAxes, fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="activator_evidence", ylabel="proteins",
                 title=f"The proposal's targets, as measured — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"10l_activator_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")
    if found:
        print(f"[{org}] named targets: {dict(found)}")


if __name__ == "__main__":
    main()
