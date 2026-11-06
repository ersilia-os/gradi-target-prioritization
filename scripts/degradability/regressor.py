"""Does TabPFN-3.5 help as a REGRESSOR too, or only as a classifier?

CLAUDE.md makes TabPFN-3.5 the project's default for supervised ML, classification AND regression.
The classifier half was measured (6/6 arms better on PR-AUC); **the regressor half was adopted by
instruction, not evidence.** This closes that gap.

The target is the continuous abundance log2FC, and the baseline is the one stage 04 already
measured and superseded -- an **ESM-C ridge** scoring rho 0.515 / AUROC 0.874 on ADEP4 and
rho 0.363 / 0.753 on ONC212 (`docs/degradability.md`, "Superseded: the regression variant").

**Two metrics, because they answer different questions.** Spearman rho against the continuous
log2FC is what a regressor is actually fitting. AUROC against the BINARY hit is what makes the
number comparable to the shipped classifier and to v1 -- and ranking is what a shortlist consumes.
A regressor can win on one and lose on the other.

Sign convention: depletion is NEGATIVE log2FC, so a more-negative prediction means more-likely-a-hit
and the AUROC is computed on **-prediction**. Getting this backwards yields ~1-AUROC, which looks
like catastrophic failure rather than a sign error.

Same folds as everything else: stage 04's cluster-grouped splitters, materialised once and handed
to both estimators, so the comparison is PAIRED.

Run with the `gradi` env (dispatches to gradi-tabpfn). ~50 hosted calls = ~500,000 credits.
  python scripts/degradability/regressor.py
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402

OUT = REPO_ROOT / "output" / "results" / "degradability" / "regressor_benchmark.tsv"


def stage04():
    spec = importlib.util.spec_from_file_location("s4", REPO_ROOT / "scripts" / "degradability" / "predict.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.VERBOSE = False
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--activator", nargs="+", default=list(D.ACTIVATORS), choices=list(D.ACTIVATORS))
    ap.add_argument("--estimator", nargs="+", default=["ridge", "tabpfn"],
                    choices=["ridge", "tabpfn"])
    args = ap.parse_args()

    from scipy.stats import spearmanr
    from sklearn.linear_model import RidgeCV
    from sklearn.metrics import roc_auc_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    m = stage04()
    feats = m.feature_frame("saureus")
    cols = D.block_columns("esmc", feats)
    lab = D.load_labels()

    print("=" * 92, flush=True)
    print("TabPFN as a REGRESSOR on continuous log2FC -- the half adopted by instruction, not evidence")
    print("=" * 92, flush=True)
    print("  baseline  ESM-C ridge, as measured in the superseded regression variant:")
    print("            ADEP4 rho 0.515 / AUROC 0.874  |  ONC212 rho 0.363 / 0.753", flush=True)
    print(f"  cv        cluster-grouped, {D.N_SPLITS} folds, {D.N_CV_SEEDS} seeds, PAIRED")
    print("=" * 92, flush=True)

    rows = []
    for act in args.activator:
        sub = lab[lab[f"{act}_log2fc"].notna() & lab["uniprot_ac"].isin(feats.index)].copy()
        X = feats.loc[sub["uniprot_ac"], cols].to_numpy(dtype=float)
        y = sub[f"{act}_log2fc"].to_numpy(dtype=float)
        hit = sub[act].to_numpy(dtype=float)
        groups = sub["cluster"].to_numpy()
        ok = ~np.isnan(hit)          # AUROC needs the binary call; rho does not
        print(f"\n  {act}: n={len(sub)}  with binary call={int(ok.sum())}  "
              f"log2FC median {np.median(y):+.2f}")

        for est in args.estimator:
            rhos, aucs = [], []
            for si in range(D.N_CV_SEEDS):
                cv = dict(m._splitters(D.N_SPLITS, D.SEED + si))["clustered"]
                splits = list(cv.split(X, (hit == 1).astype(int), groups=groups))
                pred = np.full(len(y), np.nan)
                for k, (tr, te) in enumerate(splits):
                    if est == "ridge":
                        pipe = make_pipeline(StandardScaler(),
                                             RidgeCV(alphas=np.logspace(-2, 4, 25)))
                        pred[te] = pipe.fit(X[tr], y[tr]).predict(X[te])
                    else:
                        pred[te] = m.predict_fold_regression(
                            X[tr], y[tr], X[te], tag=f"reg/{act}/s{D.SEED + si}/f{k}")
                rhos.append(spearmanr(pred, y).statistic)
                # depletion is negative, so rank on -prediction
                aucs.append(roc_auc_score(hit[ok], -pred[ok]))
            rows.append({"activator": act, "estimator": est, "n": len(sub),
                         "n_binary": int(ok.sum()), "n_cv_seeds": D.N_CV_SEEDS,
                         "spearman": round(float(np.mean(rhos)), 4),
                         "spearman_sd": round(float(np.std(rhos)), 4),
                         "roc_auc": round(float(np.mean(aucs)), 4),
                         "roc_auc_sd": round(float(np.std(aucs)), 4),
                         "spearman_seeds": ";".join(f"{v:.4f}" for v in rhos),
                         "roc_auc_seeds": ";".join(f"{v:.4f}" for v in aucs)})
            print(f"    {est:8} rho {rows[-1]['spearman']:.4f} +/- {rows[-1]['spearman_sd']:.4f}   "
                  f"AUROC {rows[-1]['roc_auc']:.4f} +/- {rows[-1]['roc_auc_sd']:.4f}")

    res = pd.DataFrame(rows)
    print("\n" + "-" * 92, flush=True)
    print("VERDICT - paired per-seed difference, TabPFN minus ridge")
    print("-" * 92, flush=True)
    for act in args.activator:
        r = res[res.activator == act].set_index("estimator")
        if not {"ridge", "tabpfn"} <= set(r.index):
            continue
        for col, lab_ in (("spearman", "rho  "), ("roc_auc", "AUROC")):
            a = np.array([float(v) for v in r.loc["tabpfn", f"{col}_seeds"].split(";")])
            b = np.array([float(v) for v in r.loc["ridge", f"{col}_seeds"].split(";")])
            d = a - b
            se = float(np.std(d, ddof=1) / np.sqrt(len(d)))
            call = "BETTER" if d.mean() > 2 * se else ("WORSE" if d.mean() < -2 * se
                                                       else "no difference")
            print(f"  {act:7} {lab_} {d.mean():+.4f} +/- {se:.4f}  "
                  f"[{int((d > 0).sum())}/{len(d)} seeds]  {call}")
    print("\n  A regressor can win on rho and lose on AUROC: rho scores the continuous fit, AUROC", flush=True)
    print("  scores the ranking a shortlist actually consumes. Report both; prefer AUROC for")
    print("  comparability with the shipped classifier and with v1.", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT, sep="\t", index=False)
    print(f"\n  wrote {OUT.relative_to(REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
