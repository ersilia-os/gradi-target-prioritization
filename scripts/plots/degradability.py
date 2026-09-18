"""ROC curves for the stage-04 degradability classifier.

Reads `evidence/oof_<activator>_saureus.tsv` -- the out-of-fold probabilities stage 04 already
wrote -- and draws one panel per activator: the five cluster-grouped folds of the FIRST CV seed
individually, plus their mean. Nothing is refitted here.

The AUROC in each panel title is the canonical one from `cv_<activator>.tsv`: the mean +/- SD over
`N_CV_SEEDS` repeats of the split. The drawn curves come from one seed, because five overlaid seeds
of five folds each is 25 lines and reads as noise.

Two reference lines are drawn on purpose:

  * the diagonal, i.e. no skill;
  * the **cross-assay reference**, as its own curve: the OTHER activator's measured log2FC used to
    rank this activator's calls, on the proteins both screens measured. That is the best a perfect
    model of one chemistry could do at predicting the other, and it is what a predictor should be
    read against -- not 1.0. Drawn as a curve rather than a horizontal line because an AUROC is an
    area, not a true-positive rate.

Run with the `gradi` env:
    python scripts/plots/degradability.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import stylia
from sklearn.metrics import roc_auc_score, roc_curve

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "degradability"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

# A common grid for averaging: each fold's TPR is interpolated onto it before the mean is taken.
FPR_GRID = np.linspace(0.0, 1.0, 201)


def _estimator_name() -> str:
    """Which estimator produced the shipped predictions, read from the model artifact.

    NOT hard-coded: this legend said "forest" for the whole life of the forest, and would have kept
    saying it after the TabPFN retrofit -- a shipped figure asserting the wrong model.
    """
    import numpy as _np

    for act in D.ACTIVATORS:
        p = D.EVIDENCE_DIR / f"model_{act}.npz"
        if p.exists():
            try:
                return str(_np.load(p, allow_pickle=True)["estimator"])
            except Exception:
                pass
    return "model"


def cross_assay_roc(activator: str):
    """The other activator's measured readout, ranking THIS activator's calls.

    A real predictor with a real ROC curve, on the proteins both screens measured -- the best a
    perfect model of one chemistry could manage on the other. `-log2FC` because depletion is
    negative, so more-negative must rank as more-likely-a-hit.
    """
    other = "onc212" if activator == "adep4" else "adep4"
    lab = D.load_labels()
    m = lab[lab[activator].notna() & lab[f"{other}_log2fc"].notna()]
    y = m[activator].to_numpy(dtype=float).astype(int)
    x = -m[f"{other}_log2fc"].to_numpy(dtype=float)
    fpr, tpr, _ = roc_curve(y, x)
    return fpr, tpr, float(roc_auc_score(y, x)), len(m), other


def plot_roc(ax, oof: dict, activator: str, abc: str) -> float:
    """Per-fold ROC curves, their mean, and the cross-assay reference, for one activator."""
    nc = stylia.NamedColors()
    folds, hit, prob = oof["fold"], oof["hit"], oof["prob"]

    tprs = []
    for k in sorted(set(folds)):
        m = folds == k
        if hit[m].sum() < 5:
            continue
        fpr, tpr, _ = roc_curve(hit[m], prob[m])
        tprs.append(np.interp(FPR_GRID, fpr, tpr))
        ax.plot(fpr, tpr, color=nc.get("plum", lighten=0.6))

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[0], mean_tpr[-1] = 0.0, 1.0
    ax.plot(FPR_GRID, mean_tpr, color=nc.plum, label=f"ESM-C {_estimator_name()}")

    c_fpr, c_tpr, c_auc, c_n, other = cross_assay_roc(activator)
    ax.plot(c_fpr, c_tpr, color=nc.orange, linestyle=":",
            label=f"{other} measured (n={c_n})")
    ax.plot([0, 1], [0, 1], color=nc.gray, linestyle="--")

    # The headline number comes from `cv_<activator>.tsv`, not from this pooled curve. They differ
    # slightly -- pooling the seed-averaged probability scores a little higher than the mean of the
    # per-seed AUROCs (measured 0.867 vs 0.865 on ADEP4 under the forest) -- and having two headline
    # numbers in circulation is worse than picking the canonical one.
    row = D.load_cv(activator)
    row = row[row["block"] == "esmc"].iloc[0]
    auc, sd = float(row["roc_auc_clustered"]), float(row["roc_auc_clustered_sd"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right")
    stylia.label(ax, xlabel="False positive rate", ylabel="True positive rate",
                 title=f"{activator}  AUROC {auc:.3f}±{sd:.3f}  ·  cross-assay {c_auc:.3f}",
                 abc=abc)
    return auc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--activator", nargs="+", default=list(D.ACTIVATORS),
                    choices=list(D.ACTIVATORS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    data = {}
    for a in args.activator:
        t = D.load_oof(a)
        data[a] = {"fold": t["fold"].to_numpy(), "hit": t["hit"].to_numpy(),
                   "prob": t["oof_esmc"].to_numpy()}

    fig, axs = stylia.create_figure(1, len(data))
    for abc, (a, oof) in zip("ABCD", data.items()):
        auc = plot_roc(axs.next(), oof, a, abc)
        say(f"  {a:<8} n={len(oof['hit']):>5}  positives={int(oof['hit'].sum()):>4}  "
            f"pooled out-of-fold AUROC {auc:.3f}  folds={len(set(oof['fold']))}")

    out = OUT_DIR / "roc_curves.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
