"""Six-panel slide for the activated-ClpP sequence-feature analysis (docs §3.3c, from 10e).

The story the six panels tell, in order:

  1. What labels exist at all — four label sets, and how few positives some carry.
  2. What correlates with ADEP4 depletion — forest plot with cluster-bootstrap CIs. Disorder looks
     good (0.697)... but length looks better (0.225 = 0.775 inverted).
  3. Why panel 2 is misleading — length and disorder are correlated, so disorder's apparent signal
     may just be size.
  4. THE GATE — nested models. Disorder adds ~0.01 over length. It is not the driver.
  5. SHAP over all features. Disorder collects ~45% of attribution yet adds nothing (panel 4) —
     SHAP splits credit among correlated features, the nested models measure what is ADDED. And a
     gradient-boosted model on everything still fails to beat logistic length-only, so the
     conclusion is not an artifact of model form. Motifs attribute 0%, exactly as a
     recognition-bypassing mechanism predicts.
  6. What does hold — abundance susceptibility generalises across two distinct activators
     (~0.69-0.77), while cleavage sits at chance.

The MS-depth control (is the size effect an artifact of peptide count? no — it survives in all four
quartiles) lives in 10e's report and `activator_peptide_strata.csv`, not here: it is a control, not a
result, and this figure is capped at six panels.

House style matches 06i-06m: stylia "slide", NPG palette, 2x3 panels, 16:9.
Run with the `gradi` env, after `scripts/10e_activator_features.py`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)  # stylia rmtree's this on import
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

RES = D.REPO_ROOT / "output" / "results" / "other"
OUT = D.REPO_ROOT / "output" / "plots" / "10l_activator_features.png"

NPG = stylia.CategoricalPalette("npg").colors
SIGNAL = "#00A087"   # NPG green — CI excludes the 0.5 null
NOSIG = "#C9C9C7"    # grey — indistinguishable from chance
WARN = "#E64B35"     # NPG red — at or below chance
BLUE = "#4DBBD5"     # NPG blue — the "other experiment" / secondary series
DARK = "#3C5488"     # NPG dark blue

TARGET_LABEL = {
    "conlon_abundance": "ADEP4\nabundance", "conlon_cleavage": "ADEP4\ncleavage",
    "jacques_abundance": "ONC212\nabundance", "jacques_cleavage": "ONC212\ncleavage",
    "conlon_abundance_bin": "ADEP4\nabundance", "conlon_cleavage_bin": "ADEP4\ncleavage",
    "jacques_abundance_bin": "ONC212\nabundance", "jacques_cleavage_bin": "ONC212\ncleavage",
}
FEAT_LABEL = {
    "length": "length", "log_length": "log length", "disorder_mean": "disorder (mean)",
    "disorder_frac": "disorder fraction", "disorder_longest_run": "longest disordered run",
    "disorder_nterm30": "disorder, N-term 30", "disorder_cterm30": "disorder, C-term 30",
    "disorder_nterm50": "disorder, N-term 50", "disorder_cterm50": "disorder, C-term 50",
    "frac_charged": "charged fraction", "gravy": "GRAVY (hydrophobicity)",
    "frac_hydrophobic": "hydrophobic fraction", "frac_PEST": "PEST fraction",
    "frac_GS": "Gly+Ser fraction", "nend_imet_cleaved": "iMet cleavable",
    "nend_primary_destabilizing": "N-end destabilizing",
    "motif_cterm_cm1_strict": "motif CM1 strict (ssrA)", "motif_cterm_cm1_broad": "motif CM1 broad",
    "motif_cterm_cm2": "motif CM2 (MuA)", "motif_nterm_nm1": "motif NM1 (Flynn)",
    "motif_nterm_nm2": "motif NM2", "motif_nterm_nm3": "motif NM3 (Flynn)",
}


def _auroc(y: np.ndarray, x: np.ndarray) -> float:
    n1 = int(y.sum()); n0 = int(len(y) - n1)
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(x)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _panel_title(ax, n: str, text: str) -> None:
    ax.set_title(f"{n}  {text}", loc="left", fontsize="medium")


def _clear_labels(ax) -> None:
    """stylia pre-fills every axis with 'X-axis / Units' placeholders; clear them.

    Called first in each panel, so a panel that genuinely wants a label sets it afterwards and a panel
    that does not is left blank rather than shipping the placeholder.
    """
    ax.set_xlabel(""); ax.set_ylabel("")


def _legend_above(ax, handles=None, ncol=1) -> None:
    """Legend inside the axes at loc='best'.

    Two earlier attempts failed: fixed corners collided with bars, and parking it above the axes
    collided with the panel title. `loc="best"` scores candidate positions against the artists and is
    the right tool; a translucent frame keeps it readable if the panel is genuinely full.
    """
    kw = dict(loc="best", fontsize="small", frameon=True, framealpha=0.88, edgecolor="none",
              ncol=ncol, handlelength=1.1, columnspacing=1.0, borderpad=0.4)
    ax.legend(handles=handles, **kw) if handles else ax.legend(**kw)


# --------------------------------------------------------------------------- panels
def p1_labels(ax, feats: pd.DataFrame) -> None:
    """How much labelled data exists, and how imbalanced it is."""
    _clear_labels(ax)
    cols = ["conlon_abundance_bin", "jacques_abundance_bin",
            "conlon_cleavage_bin", "jacques_cleavage_bin"]
    rows = []
    for c in cols:
        v = feats[c].dropna()
        rows.append((TARGET_LABEL[c].replace("\n", " "), len(v), int(v.sum())))
    y = np.arange(len(rows))[::-1]
    tot = [r[1] for r in rows]; pos = [r[2] for r in rows]
    ax.barh(y, tot, height=0.55, color=NOSIG, label="measured", zorder=2)
    ax.barh(y, pos, height=0.55, color=SIGNAL, label="positive (susceptible)", zorder=3)
    for yy, t, p in zip(y, tot, pos):
        ax.text(t + 30, yy, f"{p:,} / {t:,}", va="center", fontsize="small", color="#555555")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize="small")
    ax.set_xlabel("proteins")
    ax.set_xlim(0, max(tot) * 1.34)
    _legend_above(ax)
    _panel_title(ax, "1", "The labels that exist")


def p2_forest(ax, stats: pd.DataFrame) -> None:
    """Single-feature AUROC vs ADEP4 abundance, with cluster-bootstrap CIs."""
    _clear_labels(ax)
    s = stats[stats["target"] == "conlon_abundance"].copy()
    s["dev"] = (s["auroc"] - 0.5).abs()
    s = s.sort_values("dev", ascending=True).tail(11)
    y = np.arange(len(s))
    colors = [SIGNAL if b else NOSIG for b in s["beats_null"]]
    ax.hlines(y, s["ci_lo"], s["ci_hi"], color=colors, linewidth=2.2, zorder=2)
    ax.scatter(s["auroc"], y, s=34, color=colors, zorder=3, edgecolor="white", linewidth=0.6)
    ax.axvline(0.5, color="#777777", linestyle="--", linewidth=1, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([FEAT_LABEL.get(f, f) for f in s["feature"]], fontsize="small")
    ax.set_xlabel("AUROC  (0.5 = chance; < 0.5 = inverted)")
    ax.set_xlim(0.12, 0.88)
    # name the two that matter
    for lbl, col in (("length", WARN), ("disorder_mean", DARK)):
        if lbl in list(s["feature"]):
            i = list(s["feature"]).index(lbl)
            ax.annotate(f"{s['auroc'].iloc[i]:.3f}", (s["auroc"].iloc[i], y[i]),
                        textcoords="offset points", xytext=(0, 8), ha="center",
                        fontsize="small", color=col, fontweight="bold")
    _legend_above(ax, handles=[Patch(facecolor=SIGNAL, label="CI excludes 0.5"),
                               Patch(facecolor=NOSIG, label="indistinguishable")])
    _panel_title(ax, "2", "What correlates with ADEP4 depletion")


def p3_confound(ax, feats: pd.DataFrame) -> None:
    """Length and disorder are correlated — so panel 2 cannot separate them."""
    _clear_labels(ax)
    d = feats.dropna(subset=["conlon_abundance_bin", "log_length", "disorder_mean"])
    dep = d["conlon_abundance_bin"] == 1
    ax.scatter(d.loc[~dep, "log_length"], d.loc[~dep, "disorder_mean"], s=5, alpha=0.3,
               color=NOSIG, linewidth=0, label="not depleted", zorder=2)
    ax.scatter(d.loc[dep, "log_length"], d.loc[dep, "disorder_mean"], s=8, alpha=0.75,
               color=WARN, linewidth=0, label="depleted ≥2×", zorder=3)
    rho = spearmanr(d["log_length"], d["disorder_mean"]).statistic
    ax.text(0.03, 0.04, f"Spearman(length, disorder) = {rho:+.2f}", transform=ax.transAxes,
            fontsize="small", color="#555555")
    ax.set_xlabel("log10 protein length"); ax.set_ylabel("mean predicted disorder")
    _legend_above(ax)
    _panel_title(ax, "3", "…but the two features are entangled")


def p5_gate(ax, inc: pd.DataFrame) -> None:
    """THE GATE — nested models. Does disorder add anything beyond length?"""
    _clear_labels(ax)
    order = ["conlon_abundance_bin", "jacques_abundance_bin",
             "conlon_cleavage_bin", "jacques_cleavage_bin"]
    inc = inc.set_index("target").reindex([t for t in order if t in set(inc["target"])])
    series = [("length", WARN), ("disorder", DARK),
              ("length+disorder", SIGNAL), ("length+disorder+composition", BLUE)]
    x = np.arange(len(inc)); w = 0.2
    for i, (name, col) in enumerate(series):
        ax.bar(x + (i - 1.5) * w, inc[name], w, color=col, label=name, zorder=3)
    ax.axhline(0.5, color="#777777", linestyle="--", linewidth=1, zorder=2)
    for xi, (t, r) in enumerate(inc.iterrows()):
        g = r["gain_disorder_over_length"]
        ax.text(xi, max(r[n] for n, _ in series) + 0.015, f"gain {g:+.3f}",
                ha="center", fontsize="small",
                color=SIGNAL if g >= 0.02 else "#555555")
    ax.set_xticks(x)
    short = {"conlon_abundance_bin": "ADEP4\nabund.", "jacques_abundance_bin": "ONC212\nabund.",
             "conlon_cleavage_bin": "ADEP4\ncleav.", "jacques_cleavage_bin": "ONC212\ncleav."}
    ax.set_xticklabels([short.get(t, t) for t in inc.index], fontsize="small")
    ax.set_xlabel("")
    ax.set_ylabel("AUROC (cluster-grouped CV)")
    ax.set_ylim(0.40, 0.94)
    _legend_above(ax, ncol=2)
    _panel_title(ax, "4", "The gate: disorder adds ~nothing over length")


def p6_cross(ax, cx: pd.DataFrame) -> None:
    """Cross-activator generalisation — the only non-circular test available."""
    _clear_labels(ax)
    if cx.empty:
        ax.set_axis_off(); return
    lbl = [f"{TARGET_LABEL.get(r['train_on'], r['train_on'])}".replace("\n", " ")
           + "\n→ " + f"{TARGET_LABEL.get(r['evaluate_on'], r['evaluate_on'])}".replace("\n", " ")
           for _, r in cx.iterrows()]
    y = np.arange(len(cx))[::-1]; h = 0.34
    ax.barh(y + h / 2, cx["auroc_same_experiment"], h, color=NOSIG,
            label="same experiment (in-sample label)", zorder=3)
    colors = [SIGNAL if v >= 0.60 else WARN for v in cx["auroc_other_experiment"]]
    ax.barh(y - h / 2, cx["auroc_other_experiment"], h, color=colors,
            label="OTHER activator (held-out label)", zorder=3)
    for yy, v in zip(y, cx["auroc_other_experiment"]):
        ax.text(v + 0.008, yy - h / 2, f"{v:.3f}", va="center", fontsize="small", color="#333333")
    ax.axvline(0.5, color="#777777", linestyle="--", linewidth=1, zorder=2)
    ax.set_yticks(y); ax.set_yticklabels(lbl, fontsize="small")
    ax.set_xlabel("AUROC")
    ax.set_xlim(0.40, 0.90)
    _legend_above(ax, ncol=1)
    _panel_title(ax, "6", "Does it generalise to the other activator?")


def p7_shap(ax, summ: pd.DataFrame) -> None:
    """SHAP attribution, coloured by feature family.

    The subtlety worth reading off this panel: disorder collects ~45% of attribution, yet adds ~0.01
    AUROC over length (panel 5). Both are true — SHAP SPLITS credit among correlated features, while
    the nested models measure what is ADDED. Disorder carries information about the label; it does not
    carry information that length does not already carry.
    """
    _clear_labels(ax)
    if summ.empty:
        ax.set_axis_off(); return
    fam = lambda f: ("length" if f.startswith("log_length")
                     else "disorder" if f.startswith("disorder_")
                     else "motif" if f.startswith(("motif_", "nend_"))
                     else "composition")
    col = {"length": WARN, "disorder": DARK, "composition": BLUE, "motif": NOSIG}
    s = summ.head(11).iloc[::-1]
    y = np.arange(len(s))
    ax.barh(y, s["mean_abs_shap"], height=0.62,
            color=[col[fam(f)] for f in s["feature"]], zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([FEAT_LABEL.get(f, f) for f in s["feature"]], fontsize="small")
    ax.set_xlabel("mean |SHAP|")
    ax.set_xlim(0, float(s["mean_abs_shap"].max()) * 1.30)
    # family totals — the number that actually matters
    tot = summ.assign(fam=summ["feature"].map(fam)).groupby("fam")["share_pct"].sum()
    txt = "\n".join(f"{k} {tot.get(k, 0):.0f}%" for k in ("length", "disorder", "composition", "motif"))
    auc = summ.attrs.get("cv_auroc")
    sub = f"family share\n{txt}"
    if auc:
        sub += f"\nGBM on all features: AUROC {auc:.3f}"
    ax.set_xlabel("mean |SHAP|")
    _legend_above(ax, handles=[Patch(facecolor=col[k], label=k) for k in
                               ("length", "disorder", "composition", "motif")], ncol=4)
    _panel_title(ax, "5", "SHAP attribution (all features)")
    ax.text(0.97, 0.05, sub, transform=ax.transAxes, ha="right", va="bottom",
            fontsize="small", color="#555555", linespacing=1.5)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    need = ["activator_features.csv", "activator_feature_stats.csv",
            "activator_incremental.csv", "activator_cross_activator.csv",
            "activator_shap_summary.csv"]
    missing = [n for n in need if not (RES / n).exists()]
    if missing:
        raise SystemExit(f"missing {missing} — run scripts/10e_activator_features.py first")
    feats = pd.read_csv(RES / "activator_features.csv")
    stats = pd.read_csv(RES / "activator_feature_stats.csv")
    inc = pd.read_csv(RES / "activator_incremental.csv")
    cx = pd.read_csv(RES / "activator_cross_activator.csv")
    shsum = pd.read_csv(RES / "activator_shap_summary.csv")
    # 10e writes the run-level numbers as columns because .attrs is lost on a CSV round-trip
    if "cv_auroc" in shsum.columns:
        shsum.attrs["cv_auroc"] = float(shsum["cv_auroc"].iloc[0])
        shsum = shsum.drop(columns=[c for c in ("cv_auroc", "target") if c in shsum.columns])

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)   # 16:9, 2x3
    # stylia hands back an AxisManager, consumed with .next() — same idiom as 06i-06m
    p1_labels(axs.next(), feats)
    p2_forest(axs.next(), stats)
    p3_confound(axs.next(), feats)
    p5_gate(axs.next(), inc)
    p7_shap(axs.next(), shsum)
    p6_cross(axs.next(), cx)
    # ONE short line only. A long suptitle is laid out as a single unwrapped line, and with
    # bbox_inches="tight" that WIDENS the whole figure — which visibly redistributed the panels when
    # the verdict text was crammed in here. The conclusions live in the docstring, 10e's report and
    # docs/degradability_datasets.md instead.
    fig.suptitle("Activated-ClpP susceptibility from sequence — S. aureus (ADEP4 + ONC212), "
                 f"n={len(feats):,} proteins", fontsize="large", x=0.01, ha="left")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(args.out))
    print(f"[out] {args.out.relative_to(D.REPO_ROOT)}")


if __name__ == "__main__":
    main()
