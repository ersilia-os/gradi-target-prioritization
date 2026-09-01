"""2x3 slide: why degrons are the wrong instrument for THIS modality (docs §3.2) — E. coli K-12.

Companion to `10k_degron_plots.py`. 10k asks "do degron motifs work at all?" and answers: barely.
This slide asks the question that decides the axis — "do they work for **partnerless activated
ClpP**, the handle Gr-ADI funds?" — and answers: they cannot, by construction.

  1  the mechanism      a degron is read by ClpX/ClpA; the activator occupies that docking cleft
  2  cleavage rate      by motif status, against baseline
  3  what predicts it   AUROC for activated-ClpP cleavage
  4  the overlap        the degron set and the activated-ClpP set barely intersect
  5  machine inventory  what E. coli has, and what survives partnerless activation
  6  score distribution degron_score of cleaved vs uncleaved proteins

PROVENANCE, which every activator panel states on its own face: **the ADEP4 (Conlon 2013) and
ONC212 (Jacques 2020) proteomics were measured in *S. aureus*.** No equivalent dataset exists for
E. coli or K. pneumoniae — the Gr-ADI SoW notes that gap itself. The E. coli values here are
ortholog transfers by DIAMOND RBH (`10c`), covering 609/4,403 at a median ~42% identity. Panels 2,
3, 4 and 6 are therefore S. aureus measurements attributed to E. coli orthologs, and are labelled
that way rather than presented as native E. coli data.

Reads <prefix>_deg_degrons.csv (10b), <prefix>_clpp_activator.csv (10c), <prefix>_disorder.csv (10d),
and the proteome TSV (00a) for the machine inventory in panel 5.
Output: output/plots/10o_degron_relevance_<prefix>.png. Run with `gradi`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# See the note in 10j/10k: stylia rmtree's the matplotlib cache dir at import time.
_MPL_CACHE = Path(__file__).resolve().parents[1] / "tmp" / "mplcache"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL

NTERM_C, CTERM_C = "#3C5488", "#E64B35"
HIT, MISS, GREY, BG = "#00A087", "#C9C9C7", "#555555", "#D8D8D6"
INK = "#2B2333"

PRETTY = {
    "cterm_cm1_strict": "CM1 strict",
    "cterm_cm1_broad": "CM1 broad",
    "cterm_cm2": "CM2 (MuA)",
    "nterm_nm2": "NM2",
    "nterm_nm1": "NM1",
    "nterm_nm3": "NM3",
    "nend_primary_destabilizing": "N-end destab.",
    "nend_imet_cleaved": "iMet cleaved",
}

# The degradation machinery, by role. Presence is MEASURED against the proteome, not asserted.
# `usable` = does this subunit still act once a partnerless activator occupies the docking cleft.
MACHINES = [
    ("ClpP", "clpP", "protease core", True),
    ("ClpX", "clpX", "unfoldase", False),
    ("ClpA", "clpA", "unfoldase", False),
    ("ClpS", "clpS", "adaptor", False),
    ("SspB", "sspB", "adaptor", False),
    ("Lon", "lon", "unfoldase+protease", False),
    ("HslU", "hslU", "unfoldase", False),
    ("FtsH", "ftsH", "unfoldase+protease", False),
    ("ClpC", "clpC", "unfoldase", False),
    ("McsB", "mcsB", "adaptor", False),
]
# Every degron in 10b is read by one of these, and not one of them is ClpP.
SOURCE = "ADEP4 / ONC212 measured in S. aureus, transferred by RBH"


def _save(out: Path) -> None:
    """stylia.save_figure ignores its `pad` argument; see the note in 10j."""
    import matplotlib.pyplot as plt
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(str(out), dpi=600, transparent=False, bbox_inches="tight")


def _note(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", color="#999",
            fontsize=SS)
    ax.set_xticks([]); ax.set_yticks([])


def _wilson(k: int, n: int, z: float = 1.96):
    """Wilson score interval for a proportion. Behaves at k=0 and tiny n, unlike Wald."""
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def _box(ax, x, y, w, h, text, edge, fontsize, textcolor=INK):
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                facecolor="#FFFFFF", edgecolor=edge, linewidth=1.1, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            color=textcolor, zorder=3)


def _sub(ax, text):
    """One-line provenance note under the axes, in place of a caption paragraph."""
    ax.set_xlabel(ax.get_xlabel() + f"\n{text}")
    ax.xaxis.label.set_fontsize(SS - 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(D.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = D.ORGANISMS[org]
    orgname = D.ORG_DISPLAY[org]

    d = pd.read_csv(D.results_dir(org) / f"{prefix}_deg_degrons.csv")
    n = len(d)
    act_f = D.results_dir(org) / f"{prefix}_clpp_activator.csv"
    act = pd.read_csv(act_f) if act_f.exists() else pd.DataFrame()
    dis_f = D.results_dir(org) / f"{prefix}_disorder.csv"
    dis = pd.read_csv(dis_f) if dis_f.exists() else pd.DataFrame()

    motifs = [m for m in D.MOTIFS if m in d.columns]
    feats = motifs + [f for f in ("nend_primary_destabilizing", "nend_imet_cleaved")
                      if f in d.columns]
    wcols = [f for f in feats if D.MOTIF_WEIGHTS.get(f, 0.0) > 0]
    d["has_weighted_motif"] = d[wcols].any(axis=1) if wcols else False

    j = d.merge(act, on="uniprot_accession", how="inner") if not act.empty else pd.DataFrame()
    pid = float(j["pident"].median()) if "pident" in j.columns and len(j) else float("nan")
    print(f"[{org}] {n:,} proteins; {len(act):,} with S. aureus activated-ClpP evidence "
          f"transferred by RBH (median {pid:.0f}% id); {len(j):,} in both", flush=True)

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: the mechanism ----
    # The only drawn panel on the slide. It is first because nothing measured below is interpretable
    # until the reader knows the activator displaces the subunit that reads degrons.
    ax = axs.next()
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    fs = SS - 1
    ax.text(0.1, 9.4, "NATIVE ClpXP", fontsize=SS, color=INK, fontweight="bold", va="center")
    _box(ax, 0.1, 7.0, 2.5, 1.6, "degron\n-LAA", CTERM_C, fs, CTERM_C)
    _box(ax, 3.3, 7.0, 3.3, 1.6, "ClpX / ClpA\nreads it, PULLS", NTERM_C, fs, NTERM_C)
    _box(ax, 7.3, 7.0, 2.5, 1.6, "ClpP\ncleaves", INK, fs)
    for x0, x1 in ((2.7, 3.2), (6.7, 7.2)):
        ax.annotate("", xy=(x1, 7.8), xytext=(x0, 7.8),
                    arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.3))

    ax.plot([0.1, 9.9], [5.9, 5.9], color=BG, lw=1.2)
    ax.text(0.1, 5.0, "ACTIVATED ClpP  (no partner)", fontsize=SS, color=HIT,
            fontweight="bold", va="center")
    _box(ax, 0.1, 2.6, 2.5, 1.6, "degron\n-LAA", MISS, fs, "#9A9A98")
    _box(ax, 3.3, 2.6, 3.3, 1.6, "activator fills the\nClpX docking cleft", HIT, fs, HIT)
    _box(ax, 7.3, 2.6, 2.5, 1.6, "ClpP\npore open", INK, fs)
    ax.annotate("", xy=(3.2, 3.4), xytext=(2.7, 3.4),
                arrowprops=dict(arrowstyle="-|>", color=MISS, lw=1.3))
    # Drawn, not typed: Arial has no U+2715 and matplotlib drops it silently.
    ax.scatter([2.95], [3.4], marker="x", s=110, color=CTERM_C, linewidths=2.4, zorder=5)
    ax.text(5.0, 1.5, "the degron is a handle for the unfoldase, not for ClpP",
            fontsize=fs, color=GREY, ha="center", va="center")
    ax.set_title(f"Who reads a degron? — {orgname}", fontsize=stylia.SLIDE_FONTSIZE, color=INK)

    # ---- panel 2: cleavage rate by motif status ----
    ax = axs.next()
    if j.empty or "activator_cleaved" not in j.columns:
        _note(ax, "no activated-ClpP evidence")
    else:
        base_k, base_n = int(j["activator_cleaved"].sum()), len(j)
        rows = []
        for f in feats:
            s = j[j[f] == True]  # noqa: E712
            if len(s) < 5:       # below this the interval spans the axis and says nothing
                print(f"[{org}] panel 2 omits {f}: only {len(s)} transferred proteins carry it",
                      flush=True)
                continue
            rows.append((f, int(s["activator_cleaved"].sum()), len(s)))
        if not rows:
            _note(ax, "no motif has >=5 proteins\nwith activator evidence")
        else:
            ys = np.arange(len(rows))
            for y, (f, k, m) in zip(ys, rows):
                p, lo, hi = _wilson(k, m)
                c = HIT if D.MOTIF_WEIGHTS.get(f, 0.0) > 0 else MISS
                ax.plot([100 * lo, 100 * hi], [y, y], color=c, lw=1.4, solid_capstyle="round")
                ax.scatter([100 * p], [y], s=26, color=c, zorder=3, linewidths=0)
                ax.text(101, y, f"{k}/{m}", fontsize=SS - 1, color=GREY, va="center")
            bp, blo, bhi = _wilson(base_k, base_n)
            ax.axvspan(100 * blo, 100 * bhi, color=BG, alpha=0.5, zorder=0)
            ax.axvline(100 * bp, color=INK, lw=1.0, ls="--", zorder=1)
            ax.text(100 * bp + 1.5, len(rows) - 0.35, f"baseline {100 * bp:.0f}%", fontsize=SS - 1,
                    color=INK, va="center")
            ax.set_yticks(ys)
            ax.set_yticklabels([PRETTY.get(f, f) for f, _, _ in rows], fontsize=SS - 1)
            ax.set_xlim(0, 118); ax.set_ylim(-0.6, len(rows) + 0.1)
            ax.set_xticks([0, 25, 50, 75, 100])
            stylia.label(ax, xlabel="% cleaved", ylabel="",
                         title=f"Cleavage rate by motif — {orgname}")
            _sub(ax, SOURCE)

    # ---- panel 3: what predicts it ----
    ax = axs.next()
    jj = j.merge(dis, on="uniprot_accession", how="left", suffixes=("", "_dis")) if not j.empty \
        else pd.DataFrame()
    cand = [("degron_motif_score", "degron motif score"),
            ("degron_score", "degron_score"),
            ("cterm_exposure", "C-terminal exposure"),
            ("nterm_exposure", "N-terminal exposure"),
            ("global_disorder_frac_plddt", "global disorder"),
            ("internal_disorder_max_run", "longest internal IDR"),
            ("max_init_len_70", "terminal initiation region"),
            ("seq_len", "sequence length")]
    if jj.empty or "activator_cleaved" not in jj.columns:
        _note(ax, "no activated-ClpP evidence")
    else:
        res = []
        for c, lab in cand:
            if c not in jj.columns:
                continue
            r = D.disorder_auroc(jj[c], jj["activator_cleaved"], n_boot=1000)
            if r["auroc"] is not None:
                res.append((lab, r, c))
        res.sort(key=lambda t: t[1]["auroc"])
        ys = np.arange(len(res))
        for y, (lab, r, c) in zip(ys, res):
            col = CTERM_C if c.startswith("degron") else NTERM_C
            ax.plot([r["auroc_lo"], r["auroc_hi"]], [y, y], color=col, lw=1.4,
                    solid_capstyle="round")
            ax.scatter([r["auroc"]], [y], s=26, color=col, zorder=3, linewidths=0)
        ax.axvline(0.5, color=INK, lw=1.0, ls="--")
        ax.set_yticks(ys)
        ax.set_yticklabels([lab for lab, _, _ in res], fontsize=SS - 1)
        for y, (lab, r, c) in zip(ys, res):
            if c.startswith("degron"):
                ax.get_yticklabels()[int(y)].set_color(CTERM_C)
        ax.set_xlim(0.30, 0.72); ax.set_ylim(-0.6, len(res) - 0.4)
        ax.text(0.03, 0.96, f"n = {len(jj):,}  ·  {int(jj['activator_cleaved'].sum()):,} cleaved",
                transform=ax.transAxes, ha="left", va="top", fontsize=SS - 1, color=GREY)
        stylia.label(ax, xlabel="AUROC for cleavage (95% CI)", ylabel="",
                     title=f"What predicts cleavage? — {orgname}")
        _sub(ax, SOURCE)

    # ---- panel 4: the overlap ----
    ax = axs.next()
    if j.empty:
        _note(ax, "no activated-ClpP evidence")
    else:
        n_mot = int(d["has_weighted_motif"].sum())
        n_both = int(j["has_weighted_motif"].sum())
        n_cleaved = int(j["activator_cleaved"].sum()) if "activator_cleaved" in j else 0
        n_both_cleaved = int(j.loc[j["has_weighted_motif"], "activator_cleaved"].sum()) \
            if "activator_cleaved" in j else 0
        steps = [("proteome", n, "#BFBFBD"),
                 ("activator evidence", len(act), HIT),
                 ("...cleaved", n_cleaved, HIT),
                 ("weighted degron motif", n_mot, CTERM_C),
                 ("motif AND evidence", n_both, CTERM_C),
                 ("motif AND cleaved", n_both_cleaved, CTERM_C)]
        ys = np.arange(len(steps))[::-1]
        bars = ax.barh(ys, [max(s[1], 0) for s in steps], color=[s[2] for s in steps], height=0.66)
        ax.bar_label(bars, labels=[f"{s[1]:,}" for s in steps], padding=3, fontsize=SS - 1)
        ax.set_yticks(ys)
        ax.set_yticklabels([s[0] for s in steps], fontsize=SS - 1)
        ax.set_xscale("symlog", linthresh=10)
        ax.set_xlim(0, n * 6); ax.set_ylim(-0.7, len(steps) - 0.3)
        stylia.label(ax, xlabel="proteins (log scale)", ylabel="",
                     title=f"Do the two layers meet? — {orgname}")
        _sub(ax, SOURCE)

    # ---- panel 5: machine inventory ----
    # Measured, not asserted: gene symbols matched against the proteome TSV this pipeline uses.
    # Two dot columns replace what used to be two columns of repeated words.
    ax = axs.next()
    prot = pd.read_csv(D.proteome_tsv(org), sep="\t", usecols=["Entry", "Gene Names"])
    syms = set()
    for g in prot["Gene Names"].fillna("").astype(str):
        syms.update(x.lower() for x in g.split())
    present = {name: (gene.lower() in syms) for name, gene, _, _ in MACHINES}
    ys = np.arange(len(MACHINES))[::-1]
    for y, (name, gene, role, usable) in zip(ys, MACHINES):
        ok = present[name]
        ax.scatter([0], [y], s=95, color=HIT if ok else MISS, marker="o" if ok else "X",
                   linewidths=0, zorder=3)
        ax.scatter([1], [y], s=95, color=HIT if (ok and usable) else MISS,
                   marker="o" if (ok and usable) else "X", linewidths=0, zorder=3)
        ax.text(-0.28, y, name, fontsize=SS, color=INK if ok else "#9A9A98", va="center",
                ha="right", fontweight="bold")
        ax.text(1.34, y, role, fontsize=SS - 1, color=GREY, va="center")
    ax.text(0, len(MACHINES) - 0.5, "in the\nproteome", fontsize=SS - 1, color=GREY, ha="center",
            va="bottom", fontweight="bold")
    ax.text(1, len(MACHINES) - 0.5, "acts under\nactivated ClpP", fontsize=SS - 1, color=GREY,
            ha="center", va="bottom", fontweight="bold")
    ax.set_xlim(-1.05, 2.9); ax.set_ylim(-0.7, len(MACHINES) + 1.3)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(""); ax.set_ylabel("")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"The machinery present — {orgname}", fontsize=stylia.SLIDE_FONTSIZE, color=INK)

    # ---- panel 6: degron_score of cleaved vs uncleaved ----
    ax = axs.next()
    if j.empty or "activator_cleaved" not in j.columns:
        _note(ax, "no activated-ClpP evidence")
    else:
        for flag, lab, c in ((True, "cleaved", HIT), (False, "not cleaved", MISS)):
            v = j.loc[j["activator_cleaved"] == flag, "degron_score"].dropna().sort_values()
            if not len(v):
                continue
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=c, lw=1.8,
                    label=f"{lab} ({len(v)})")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        stylia.label(ax, xlabel="degron_score (10b)", ylabel="cumulative fraction",
                     title=f"Same score, either way — {orgname}")
        _sub(ax, SOURCE)
        ax.legend(fontsize=SS, frameon=False, loc="lower right")

    out = D.REPO_ROOT / "output" / "plots" / f"10o_degron_relevance_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    _save(out)
    print(f"[{org}] wrote {out.relative_to(D.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
