"""Score a feature matrix with lazy-qsar's LazyClassifier on FOLDS DECIDED ELSEWHERE.

**A worker, not a stage.** It runs under the `gradi-lazyqsar` interpreter and is never imported --
lazy-qsar pins `numpy==2.1.3` and `scikit-learn==1.6.1` while `gradi` runs numpy 2.4.6 / sklearn
1.9.0 under torch 2.12, so installing it into `gradi` would downgrade both and take ESM-C (stage 01)
and the stage-04 forest down with it. The same trap CLAUDE.md records for `blast`, `eggnog-mapper`
and `rdkit`. Hence the process boundary.

WHY THE FOLDS COME IN FROM OUTSIDE
----------------------------------
`LazyClassifier.fit()` reports an `oof_auc_` of its own. **Do not use it.** It comes from lazy-qsar's
internal random splits, which know nothing about sequence homology -- exactly the leakage stage 04
guards against with `StratifiedGroupKFold` over MMseqs2 clusters, and the reason a plain `cv=5` here
would read high and mean nothing (v1's `07d` made that mistake). So the caller computes the folds
with stage 04's own splitters and passes them in; this worker only fits and predicts. It is reported
in the audit beside the honest number so the gap is visible.

PROTOCOL
--------
For each seed and each fold: fit on that seed's training rows, predict the held-out rows. One model
per fold, so every prediction comes from a model that never saw the protein or its homologs.
lazy-qsar selects its portfolio and calibrates *inside* `fit`, i.e. inside the training fold, which
makes this a correct nested protocol rather than a tuned-on-test one.

IN   a .npz bundle: X (n, d) float32 | y (n,) int | folds (n_seeds, n) int, -1 never allowed
OUT  a .npz: oof (n_seeds, n) float64 -- the held-out probability of class 1
            portfolio (n_seeds, n_folds) str | fit_seconds (n_seeds, n_folds) float
            internal_oof_auc (n_seeds, n_folds) float -- lazy-qsar's own leaky estimate, for contrast
            lazyqsar_version, calibrated

Run with the `gradi-lazyqsar` env:
  ~/miniconda3/envs/gradi-lazyqsar/bin/python scripts/degradability/workers/lazyqsar_cv.py --in b.npz --out o.npz
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True, help="the .npz bundle: X, y, folds")
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--tag", default="", help="printed with every progress line")
    ap.add_argument("--uncalibrated", action="store_true",
                    help="fit with calibrated=False (default is lazy-qsar's own default, True)")
    args = ap.parse_args()

    import importlib.metadata as md

    from lazyqsar.agnostic import LazyClassifier
    from sklearn.metrics import roc_auc_score

    version = md.version("lazyqsar")
    z = np.load(args.inp, allow_pickle=False)
    X, y, folds = z["X"].astype("float32"), z["y"].astype(int), z["folds"].astype(int)
    n_seeds, n = folds.shape
    n_folds = int(folds.max()) + 1

    # A bad fold matrix would silently turn this into a train-on-test run, so check it rather than
    # trust it: every protein assigned in every seed, and no fold empty.
    if X.shape[0] != n or len(y) != n:
        sys.exit(f"FATAL shape mismatch: X {X.shape}, y {len(y)}, folds {folds.shape}")
    if (folds < 0).any():
        sys.exit("FATAL the fold matrix holds negative entries -- some protein was never assigned")
    for si in range(n_seeds):
        counts = np.bincount(folds[si], minlength=n_folds)
        if (counts == 0).any():
            sys.exit(f"FATAL seed {si} has an empty fold: {counts.tolist()}")

    print(f"lazy-qsar {version} | X {X.shape} | {int(y.sum())}/{n} positive "
          f"({100 * y.mean():.1f}%) | {n_seeds} seeds x {n_folds} folds", flush=True)

    oof = np.full((n_seeds, n), np.nan)
    portfolio = np.empty((n_seeds, n_folds), dtype=object)
    secs = np.zeros((n_seeds, n_folds))
    internal = np.full((n_seeds, n_folds), np.nan)

    for si in range(n_seeds):
        for k in range(n_folds):
            te = folds[si] == k
            tr = ~te
            t0 = time.time()
            clf = LazyClassifier(calibrated=not args.uncalibrated)
            clf.fit(X[tr], y[tr])
            oof[si, te] = clf.predict_proba(X[te])[:, 1]
            secs[si, k] = time.time() - t0
            try:
                portfolio[si, k] = ",".join(clf._model.portfolio)
                internal[si, k] = float(clf.oof_auc_)
            except Exception:                      # never let the audit fields fail the run
                portfolio[si, k] = "?"
            print(f"  {args.tag} seed {si} fold {k}: n_train={int(tr.sum())} "
                  f"n_test={int(te.sum())} {secs[si, k]:.0f}s "
                  f"portfolio={portfolio[si, k]} internal_oof_auc={internal[si, k]:.4f}", flush=True)
        # Seed complete: the honest, cluster-grouped, held-out AUROC for this seed alone. Printed as
        # a progress signal -- the reportable number is the mean over all seeds, never one of these.
        done = si + 1
        auc_si = roc_auc_score(y, oof[si])
        eta = secs[secs > 0].mean() * (n_seeds - done) * n_folds / 60
        print(f"  {args.tag} SEED {si} DONE ({done}/{n_seeds}): held-out AUROC {auc_si:.4f}  "
              f"[running mean {np.mean([roc_auc_score(y, oof[j]) for j in range(done)]):.4f}]  "
              f"~{eta:.0f} min left in this arm", flush=True)

    if np.isnan(oof).any():
        sys.exit(f"FATAL {int(np.isnan(oof).sum())} predictions are NaN -- some protein was never "
                 "in a test fold, so the out-of-fold matrix is incomplete")

    per_seed = [roc_auc_score(y, oof[si]) for si in range(n_seeds)]
    print(f"  {args.tag} per-seed held-out AUROC: "
          f"{', '.join(f'{a:.4f}' for a in per_seed)}", flush=True)
    print(f"  {args.tag} lazy-qsar's own internal (LEAKY, ungrouped) oof_auc_ mean: "
          f"{np.nanmean(internal):.4f}", flush=True)

    np.savez_compressed(
        args.out, oof=oof, portfolio=portfolio.astype(str), fit_seconds=secs,
        internal_oof_auc=internal, lazyqsar_version=version,
        calibrated=not args.uncalibrated,
    )
    print(f"  {args.tag} wrote {args.out} ({secs.sum() / 60:.1f} min of fitting)", flush=True)


if __name__ == "__main__":
    main()
