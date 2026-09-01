"""2x3 degron/motif slide for the degradability axis (docs §3.2) — E. coli K-12.

House style: stylia "slide" format, NPG palette, single organism per slide.

  1  motif prevalence      how rare a sequence degron actually is, per motif
  2  enrichment            odds ratio for fast turnover, WITH intervals — the load-bearing panel
  3  motif x accessibility what survives once the motif must also be exposed (docs §3, "never alone")
  4  the exposure map      N vs C terminal exposure, motif hits overlaid
  5  what drives the score degron_score is exposure, not motif — the slide's thesis
  6  the ssrA-like set     the CM1-strict proteins, named, with their real C-terminal residues

This slide argues a NEGATIVE result, and argues it honestly. Measured against the Nagar 2021
turnover classes, only the ssrA-like C-terminal motif beats chance — and it does so on a handful of
labelled proteins, which is why panel 2 plots confidence intervals rather than bare point estimates.
`src/degradability.py` states the conclusion in a comment ("sequence degron motifs cannot carry this
axis"); this figure is where it becomes legible, and where the zeros in `MOTIF_WEIGHTS` are justified.

Motif columns are read from whatever the CSV actually carries, not from `D.MOTIFS`: the current
`*_deg_degrons.csv` predates the 2026-08-06 addition of Flynn's real N-M1/N-M3 consensuses, so those
two motifs are absent from the data and panel 1 says so. Re-running 10b picks them up here with no
code change.

Reads <prefix>_deg_degrons.csv (10b) + <prefix>_degron_validation.csv (10b) + Nagar 2021 (10a).
Output: output/plots/10k_degrons_<prefix>.png. Run with `gradi`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# stylia does an unconditional `shutil.rmtree(mpl.get_cachedir())` at import time. The usual
# `makedirs(get_cachedir())` guard is not enough: matplotlib silently falls back to a temporary
# cache dir when the default is unwritable, so the path we create and the path stylia deletes can
# differ, and the import dies with FileNotFoundError. Pin the cache to a stable repo-local path
# before matplotlib is imported, then re-create it after stylia has wiped it.
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

ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}
NTERM_C, CTERM_C = "#3C5488", "#E64B35"
HIT, MISS, GREY, BG = "#00A087", "#C9C9C7", "#555555", "#D8D8D6"

# The two N-end-rule provenance columns are not motifs — they are carried by 10b so the legacy sign
# error stays auditable — but they belong on panels 1 and 2 precisely because their prevalence is
# what disqualifies them.
NEND_FEATURES = ["nend_primary_destabilizing", "nend_imet_cleaved"]

# Short labels; the regexes themselves are in src/degradability.py MOTIFS[*]["desc"].
PRETTY = {
    "cterm_cm1_strict": "CM1 strict   φ-A-A$",
    "cterm_cm1_broad": "CM1 broad   -A-A$",
    "cterm_cm2": "CM2 (MuA)   basic+aliphatic$",
    "nterm_nm2": "NM2   ^M-[KR]-φ-φ",
    "nterm_nm1": "NM1   T-X-K-[ILV]",
    "nterm_nm3": "NM3   φ-X-pol-X-pol-X-bas-pol",
    # No set-membership glyph here: Arial has no U+2208, and matplotlib silently drops it.
    "nend_primary_destabilizing": "N-end destabilizing   pos2 = LFYW",
    "nend_imet_cleaved": "iMet cleaved   pos2 = ACGPSTV",
}
# Exposure at which a terminus is called usable. `exposure()` is a logistic on terminal pLDDT, so the
# distribution is strongly bimodal and 0.5 is the natural cut rather than a tuned one.
EXPOSED = 0.5


def _save(out: Path) -> None:
    """stylia.save_figure with real padding.

    `stylia.save_figure(pad=...)` accepts a pad argument but ignores it — both branches call a bare
    `plt.tight_layout()` — so the 2x3 grid ends up with the matplotlib default 1.08 inter-panel pad,
    which is a lot of white at slide scale. Same dpi/bbox/facecolor as stylia so the output is
    otherwise identical.
    """
    import matplotlib.pyplot as plt
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(str(out), dpi=600, transparent=False, bbox_inches="tight")


def _note(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", color="#999",
            fontsize=SS)
    ax.set_xticks([]); ax.set_yticks([])


def _terminus(feature: str) -> str:
    """Which terminus' exposure modulates this feature. The N-end priors read the N-terminus."""
    return "cterm" if feature.startswith("cterm") else "nterm"


def _nagar_labels_for(organism: str) -> pd.DataFrame:
    """Measured turnover class per accession — direct for E. coli, orthology for Kp.

    Same derivation as scripts/10b_degrons.py:116-134, repeated here (rather than imported, since
    10b is a script) so panel 2 can rebuild the 2x2 tables it needs for intervals.
    """
    nagar = D.load_nagar()
    if nagar.empty:
        return pd.DataFrame(columns=["uniprot_accession", "halflife_class"])
    if organism == "ecoli":
        return nagar[["uniprot_accession", "halflife_class"]]
    orth = D.load_orthologs(organism)
    sub = orth[orth["species"] == D.SPECIES_ECOLI][["anchor_uniprot", "target_uniprot"]]
    ec = nagar[["uniprot_accession", "halflife_class"]].rename(
        columns={"uniprot_accession": "ec_uniprot"})
    m = sub.merge(ec, left_on="target_uniprot", right_on="ec_uniprot", how="inner")
    rank = {"fast": 0, "intermediate": 1, "stable": 2}
    m["r"] = m["halflife_class"].map(rank)
    m = m.sort_values("r").drop_duplicates("anchor_uniprot")
    return m.rename(columns={"anchor_uniprot": "uniprot_accession"})[
        ["uniprot_accession", "halflife_class"]]


def _or_ci(pos_with, n_with, pos_without, n_without):
    """Haldane-Anscombe odds ratio + 95% Wald interval on log-OR.

    10b reports the point estimate only. At n_with = 4 that estimate is nearly uninformative, and
    saying so is the whole point of panel 2, so the interval is computed here rather than stored.
    """
    a, b = pos_with, n_with - pos_with
    c, d = pos_without, n_without - pos_without
    if min(n_with, n_without) == 0:
        return None, None, None
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a / b) / (c / d)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orr, orr * np.exp(-1.96 * se), orr * np.exp(1.96 * se)


def _enrichment(df: pd.DataFrame, feats: list[str], organism: str):
    """Per-feature OR with intervals against the Nagar fast/slow labels. None if unavailable."""
    lab = _nagar_labels_for(organism)
    if lab.empty:
        return None, None, None
    j = df.merge(lab, on="uniprot_accession", how="inner")
    if j.empty:
        return None, None, None
    j["is_fast"] = j["halflife_class"] == "fast"
    base = float(j["is_fast"].mean())
    rows = []
    for f in feats:
        with_f = j[j[f] == True]      # noqa: E712
        without = j[j[f] != True]     # noqa: E712
        orr, lo, hi = _or_ci(int(with_f["is_fast"].sum()), len(with_f),
                             int(without["is_fast"].sum()), len(without))
        if orr is None:
            continue
        rows.append({"feature": f, "odds_ratio": orr, "lo": lo, "hi": hi,
                     "n_labelled_hits": len(with_f), "n_fast_hits": int(with_f["is_fast"].sum()),
                     "weight": D.MOTIF_WEIGHTS.get(f, 0.0)})
    return pd.DataFrame(rows), base, len(j)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(D.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = D.ORGANISMS[org]
    orgname = D.ORG_DISPLAY[org]
    col = ORG_COLOR[prefix]

    src = D.results_dir(org) / f"{prefix}_deg_degrons.csv"
    d = pd.read_csv(src)
    n = len(d)
    print(f"[{org}] {n:,} proteins from {src.relative_to(D.REPO_ROOT)}", flush=True)

    # Motif columns come from the FILE, not from D.MOTIFS — see the module docstring.
    motifs = [m for m in D.MOTIFS if m in d.columns]
    missing = [m for m in D.MOTIFS if m not in d.columns]
    feats = motifs + [f for f in NEND_FEATURES if f in d.columns]
    if missing:
        print(f"[{org}] NOT plotted (absent from this CSV vintage): {', '.join(missing)} "
              f"— re-run scripts/10b_degrons.py to include them", flush=True)
    for f in feats:
        k = int((d[f] == True).sum())  # noqa: E712
        print(f"[{org}]   {f:<28s} {k:5d}  ({k / n:6.2%})  weight "
              f"{D.MOTIF_WEIGHTS.get(f, 0.0):.2f}", flush=True)

    d["gene"] = (d["gene_names"].fillna("").astype(str).str.split().str[0]
                 .replace("", np.nan).fillna(d["uniprot_accession"]))

    enr, base, n_lab = _enrichment(d, feats, org)

    # Cross-check the recomputed point estimates against what 10b stored. A disagreement would mean
    # 10b is not reproducible from its own inputs, which is a finding, not a rounding detail.
    vfile = D.degradability_processed_dir(org) / f"{prefix}_degron_validation.csv"
    if enr is not None and vfile.exists():
        stored = pd.read_csv(vfile).set_index("feature")["odds_ratio"]
        for _, r in enr.iterrows():
            if r["feature"] in stored.index:
                delta = abs(round(r["odds_ratio"], 2) - stored[r["feature"]])
                if delta > 0.01:
                    print(f"[{org}] [WARN] OR mismatch vs stored validation for "
                          f"{r['feature']}: recomputed {r['odds_ratio']:.2f} vs stored "
                          f"{stored[r['feature']]:.2f}", flush=True)

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: how rare is a degron, really ----
    ax = axs.next()
    counts = {f: int((d[f] == True).sum()) for f in feats}  # noqa: E712
    order = sorted(feats, key=lambda f: counts[f])
    ys = np.arange(len(order))
    pcts = [100 * counts[f] / n for f in order]
    cols = [HIT if D.MOTIF_WEIGHTS.get(f, 0.0) > 0 else MISS for f in order]
    bars = ax.barh(ys, pcts, color=cols, height=0.68)
    ax.bar_label(bars, labels=[f"{counts[f]:,}" for f in order], padding=3, fontsize=SS)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{PRETTY.get(f, f)}\nweight {D.MOTIF_WEIGHTS.get(f, 0.0):.2f}"
                        for f in order], fontsize=SS - 1)
    ax.set_xlim(0, max(pcts) * 1.28)
    # Union, not sum: CM1-strict is a strict subset of CM1-broad (every phi-AA ending also ends -AA),
    # so adding the two counts double-counts all 11 strict hits.
    wcols = [f for f in feats if D.MOTIF_WEIGHTS.get(f, 0.0) > 0]
    weighted = int(d[wcols].any(axis=1).sum()) if wcols else 0
    ax.text(0.97, 0.36, f"{weighted:,} of {n:,} ({weighted / n:.1%}) carry a weighted motif",
            transform=ax.transAxes, ha="right", va="top", fontsize=SS, color=HIT)
    if missing:
        ax.text(0.97, 0.22, ", ".join(m.split("_")[-1].upper() for m in missing)
                            + " absent — re-run 10b",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=SS - 1, color=GREY)
    stylia.label(ax, xlabel="% of proteome carrying the motif",
                 ylabel="", title=f"How rare is a degron motif? — {orgname}")

    # ---- panel 2: does the motif actually predict fast turnover ----
    ax = axs.next()
    if enr is None or enr.empty:
        _note(ax, "no Nagar turnover labels\n(run scripts/10a_fetch_degradability.py)")
    else:
        e = enr.sort_values("odds_ratio").reset_index(drop=True)
        ys = np.arange(len(e))
        for y, r in zip(ys, e.itertuples()):
            c = HIT if r.weight > 0 else MISS
            ax.plot([r.lo, r.hi], [y, y], color=c, lw=1.4, solid_capstyle="round", zorder=2)
            ax.scatter([r.odds_ratio], [y], s=18 + 2.2 * r.n_labelled_hits, color=c, zorder=3,
                       linewidths=0)
        ax.axvline(1.0, color="#1A1A1A", lw=0.9, ls="--", zorder=1)
        ax.set_xscale("log")
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{PRETTY.get(r.feature, r.feature).split('  ')[0]}   "
                            f"(n={r.n_labelled_hits})" for r in e.itertuples()], fontsize=SS - 1)
        # Extra headroom below the lowest row so the caveat sits in an empty band rather than on top
        # of CM2's (very long) interval.
        ax.set_ylim(-0.8, len(e) + 0.3)
        ax.text(0.03, 0.97, f"{n_lab:,} labelled  ·  baseline P(fast) = {base:.3f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=SS, color=GREY)
        stylia.label(ax, xlabel="odds ratio for fast turnover (log scale, 95% CI)", ylabel="",
                     title=f"Do the motifs predict turnover? — {orgname}")

    # ---- panel 3: motif x accessibility ----
    # docs/03_degradability.md: "Motif channels must never be scored alone ... always motif x
    # accessibility". Drawn as the accounting it implies: of the already-small motif-positive sets,
    # how many sit on a terminus that is actually exposed.
    ax = axs.next()
    rng = np.random.default_rng(0)
    rows = [(f, d.loc[d[f] == True, f"{_terminus(f)}_exposure"].dropna())  # noqa: E712
            for f in order]
    bg = pd.concat([d["nterm_exposure"], d["cterm_exposure"]]).dropna()
    ys = np.arange(len(rows) + 1)
    sample = bg.sample(min(len(bg), 1500), random_state=0)
    ax.scatter(sample, np.zeros(len(sample)) + rng.normal(0, 0.11, len(sample)),
               s=1.2, color=BG, linewidths=0, zorder=1)
    for i, (f, v) in enumerate(rows, start=1):
        if not len(v):
            continue
        cs = [HIT if x >= EXPOSED else MISS for x in v]
        ax.scatter(v, np.full(len(v), i) + rng.normal(0, 0.13, len(v)), s=7, c=cs,
                   linewidths=0, alpha=0.85, zorder=3)
        k = int((v >= EXPOSED).sum())
        ax.text(1.04, i, f"{k}  {k / len(v):.0%}", fontsize=SS, color=HIT if k else CTERM_C,
                va="center", fontweight="bold")
    ax.axvline(EXPOSED, color=GREY, lw=0.8, ls=":")
    # The background row is per-TERMINUS (both termini of every protein pooled), which is the
    # like-for-like comparator for a motif that sits on one specific terminus.
    ax.text(1.04, 0, f"{int((bg >= EXPOSED).sum()):,}  {(bg >= EXPOSED).mean():.0%}", fontsize=SS,
            color=GREY, va="center")
    ax.set_yticks(ys)
    ax.set_yticklabels([f"all termini  n={len(bg):,}"]
                       + [f"{PRETTY.get(f, f).split('  ')[0]}  n={len(v):,}" for f, v in rows],
                       fontsize=SS - 1)
    ax.set_xlim(-0.04, 1.48); ax.set_ylim(-0.6, len(rows) + 0.6)
    ax.set_xticks([0, 0.5, 1.0])
    stylia.label(ax, xlabel="exposure of the motif's own terminus",
                 ylabel="", title=f"A motif is not enough — {orgname}")
    ax.xaxis.label.set_fontsize(SS)

    # ---- panel 4: the exposure map ----
    ax = axs.next()
    s = d.dropna(subset=["nterm_exposure", "cterm_exposure"])
    from matplotlib.colors import LinearSegmentedColormap
    # Density in grey, not the organism hue: the motif overlays are the content here, and a blue
    # background hides the dark-blue N-terminal points.
    cmap = LinearSegmentedColormap.from_list("dens", ["#FFFFFF", "#9E9E9C"])
    ax.hexbin(s["nterm_exposure"], s["cterm_exposure"], gridsize=26, bins="log", cmap=cmap,
              mincnt=1, extent=(0, 1, 0, 1), zorder=1)
    cpos = s[s[[m for m in motifs if m.startswith("cterm")]].any(axis=1)]
    npos = s[s[[m for m in motifs if m.startswith("nterm")]].any(axis=1)]
    ax.scatter(npos["nterm_exposure"], npos["cterm_exposure"], s=6, color=NTERM_C, linewidths=0,
               alpha=0.75, zorder=3, label=f"N-terminal motif ({len(npos)})")
    ax.scatter(cpos["nterm_exposure"], cpos["cterm_exposure"], s=6, color=CTERM_C, linewidths=0,
               alpha=0.75, zorder=3, label=f"C-terminal motif ({len(cpos)})")
    ax.axvline(EXPOSED, color=GREY, lw=0.8, ls="--"); ax.axhline(EXPOSED, color=GREY, lw=0.8,
                                                                ls="--")
    both = int(((s["nterm_exposure"] >= EXPOSED) & (s["cterm_exposure"] >= EXPOSED)).sum())
    neither = int(((s["nterm_exposure"] < EXPOSED) & (s["cterm_exposure"] < EXPOSED)).sum())
    # Quadrant counts go in the two free corners; the legend takes the third. Putting either count
    # at top-right collides with the legend at slide scale.
    # Every annotation in this panel sits on top of data, so each gets a translucent white plate.
    plate = dict(facecolor="white", alpha=0.78, edgecolor="none", pad=1.6)
    ax.text(0.97, 0.03, f"both {both:,} ({both / len(s):.0%})", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=SS, color=GREY, bbox=plate)
    ax.text(0.03, 0.03, f"neither {neither:,} ({neither / len(s):.0%})", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=SS, color=GREY, bbox=plate)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="N-terminal exposure", ylabel="C-terminal exposure",
                 title=f"Motif hits vs accessible termini — {orgname}")
    ax.legend(fontsize=SS, loc="upper left", markerscale=2.2, handletextpad=0.4,
              borderaxespad=0.3, frameon=True, framealpha=0.78, edgecolor="none",
              facecolor="white")

    # ---- panel 5: what actually drives degron_score ----
    ax = axs.next()
    v = d.dropna(subset=["degron_exposure_score", "degron_score"])
    has_motif = v["degron_motif_score"] > 0
    ax.scatter(v.loc[~has_motif, "degron_exposure_score"], v.loc[~has_motif, "degron_score"],
               s=3, color=MISS, linewidths=0, alpha=0.5, zorder=2,
               label=f"no motif contribution ({int((~has_motif).sum()):,})")
    ax.scatter(v.loc[has_motif, "degron_exposure_score"], v.loc[has_motif, "degron_score"],
               s=11, color=CTERM_C, linewidths=0, alpha=0.9, zorder=3,
               label=f"motif contributes ({int(has_motif.sum()):,})")
    rho = v["degron_exposure_score"].corr(v["degron_score"], method="spearman")
    ax.text(0.03, 0.97, f"ρ = {rho:.3f}", transform=ax.transAxes, ha="left", va="top",
            fontsize=SS, color=GREY)
    stylia.label(ax, xlabel="terminal exposure score", ylabel="degron_score",
                 title=f"What actually drives the score? — {orgname}")
    # Upper-left is the empty half of the plot: the cloud is a near-perfect diagonal.
    ax.legend(fontsize=SS, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.82),
              markerscale=2.2)

    # ---- panel 6: the ssrA-like set ----
    # Not a ranking by degron_score, which would flatter the set. The question a shortlist has to
    # answer is whether the C-terminus carrying the motif is actually reachable, so both termini are
    # plotted per protein against the same bar. In E. coli the answer is that none of them clear it,
    # and the two highest-scoring members score on the OPPOSITE terminus from their motif.
    ax = axs.next()
    key = "cterm_cm1_strict"
    sub = d[d[key] == True].sort_values("cterm_exposure")  # noqa: E712
    if sub.empty:
        _note(ax, f"no {key} hits")
    else:
        MAX_BARS = 20
        dropped = max(0, len(sub) - MAX_BARS)
        if dropped:
            print(f"[{org}] panel 6 shows the {MAX_BARS} most C-terminally exposed of "
                  f"{len(sub)} {key} proteins", flush=True)
            sub = sub.tail(MAX_BARS)
        ys = np.arange(len(sub))
        h = 0.38
        ax.barh(ys + h / 2, sub["cterm_exposure"], height=h, color=CTERM_C)
        ax.barh(ys - h / 2, sub["nterm_exposure"], height=h, color=NTERM_C)
        ax.axvline(EXPOSED, color="#1A1A1A", lw=1.0, ls="--", zorder=4)
        # Direct labels on the top row instead of a legend box: with one long N-terminal bar in the
        # middle of the panel there is no rectangle a legend can occupy without covering data.
        topi = len(sub) - 1
        ax.text(float(sub["cterm_exposure"].iloc[-1]) + 0.015, topi + h / 2, "C-term (has the motif)",
                fontsize=SS - 1, color=CTERM_C, va="center", fontweight="bold")
        ax.text(float(sub["nterm_exposure"].iloc[-1]) + 0.015, topi - h / 2, "N-term",
                fontsize=SS - 1, color=NTERM_C, va="center", fontweight="bold")
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{g:<6s} ···{c}" for g, c in zip(sub["gene"], sub["cterm_last5"])],
                           fontsize=SS - 1, fontfamily="monospace")
        ax.set_xlim(0, 1.0); ax.set_ylim(-0.8, len(sub) + 0.1)
        usable = int((sub["cterm_exposure"] >= EXPOSED).sum())
        best = sub.iloc[-1]
        ax.text(0.55, 0.04,
                f"{usable}/{len(sub)} clear the C-terminal bar"
                + (f"  ·  {dropped} not shown" if dropped else ""),
                transform=ax.transAxes, ha="left", va="bottom", fontsize=SS,
                color=CTERM_C if not usable else HIT)
        stylia.label(ax, xlabel="terminal exposure  ·  dashed line = usable-handle threshold",
                     ylabel="", title=f"The ssrA-like set — {orgname}")
        ax.xaxis.label.set_fontsize(SS)

    out = D.REPO_ROOT / "output" / "plots" / f"10k_degrons_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    _save(out)
    print(f"[{org}] wrote {out.relative_to(D.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
