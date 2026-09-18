"""Score a feature matrix with TabPFN on FOLDS DECIDED ELSEWHERE. A worker, not a stage.

TabPFN is a *tabular foundation model*: a transformer pre-trained on millions of synthetic tabular
tasks that performs supervised learning **in context**. `fit` does not train -- it stores the
training rows; `predict_proba` runs one forward pass that attends over them. So there are no
hyperparameters to tune and nothing to overfit in the usual sense, which is exactly why it is worth
testing against a forest that was hand-set and an AutoML portfolio that tunes itself.

Runs under `gradi-tabpfn` and is never imported: it pulls torch 2.14 + mlx, and `gradi` is on torch
2.12 under ESM-C. The same process-boundary rule as the lazy-qsar and ProtT5 workers.

WHY THE FOLDS COME IN FROM OUTSIDE
----------------------------------
Identical reasoning to `lazyqsar_cv.py`: the caller computes cluster-grouped folds with stage 04's
own splitters so no homology cluster spans train and test, and every estimator is scored on
byte-identical partitions. TabPFN is *transductive* -- its prediction for a row attends over the
training rows it was given -- so feeding it a leaky split would be especially flattering.

MODEL, LICENCE AND WHERE IT RUNS
--------------------------------
TabPFN-3.5 accepts up to 1M rows / 20k features, so our 1152-d x 1677-row matrix needs no PCA -- it
is scored on exactly the features the other heads get, which is what makes the comparison clean.

Two routes, both keyed by `TABPFN_TOKEN`:
  local (default)  downloads the weights. Needs the licence ACCEPTED on the Licenses tab at
                   ux.priorlabs.ai -- a valid token alone is not enough (measured: token valid,
                   `accepted: False`, download refused).
  --hosted         runs on Prior Labs' servers via `tabpfn_client`. Measured 9.3 s per fold against
                   ~95 s for lazy-qsar locally, at 10,000 credits per fold. NOTE this UPLOADS the
                   feature matrix and labels to a third party; fine for embeddings of published
                   proteomes, a decision to make consciously for anything else.

The weights are **non-commercial licensed** either way. That is an OPEN question for GraDi's
deliverable, not a settled one -- it must be resolved against how the outputs are used before this
ships externally. It does not block methods work, and stage 04 now runs on it by decision.

The hosted route fails TRANSIENTLY (a GCS 500 killed a 40-minute run at call 64 of ~160), so every
call retries with backoff. A failed attempt is not billed -- the charge lands on a completed predict
-- so retrying costs wall time, not credits.

IN   a .npz bundle: X (n, d) float32 | y (n,) int | folds (n_seeds, n) int
OUT  a .npz: oof (n_seeds, n) float64 | fit_seconds | device | model_version | tabpfn_version

  ~/miniconda3/envs/gradi-tabpfn/bin/python scripts/workers/tabpfn_cv.py --in b.npz --out o.npz
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time
from pathlib import Path

import numpy as np


MAX_ATTEMPTS = 5          # 1 try + 4 retries at 5/10/20/40 s
BACKOFF_BASE = 5

# ---------------------------------------------------------------------------------------------
# CACHE KEY FORMAT -- READ BEFORE TOUCHING `_key`.
#
# **Any change to what `_key` hashes silently invalidates the ENTIRE cache.** Not an error, not a
# warning: every lookup simply misses and every call is paid for again. Measured cost of learning
# this: adding `task` to the hashed string (`"{model}|{seed}"` -> `"{model}|{seed}|{task}"`) threw
# away ~160 entries and cost **1,540,000 credits, 7.7% of a monthly quota**, for numbers that
# reproduced identically.
#
# If you must change the format: bump KEY_VERSION, and expect to re-pay for everything. Prefer
# ADDING a new field with a default that reproduces the old digest for existing tasks.
#
# To check whether a call is cached WITHOUT spending anything, compute the key and test for the
# file. Do not "verify the cache" by making a call with made-up inputs -- a novel input always
# misses, which looks identical to a broken cache and proves nothing (this mistake was also made).
KEY_VERSION = 2           # v1: "{model}|{seed}"   v2: "{model}|{seed}|{task}"
# ---------------------------------------------------------------------------------------------


def _key(Xtr, ytr, Xte, model: str, seed: int, task: str = "classification") -> str:
    """Content address for one predict call. Hashes the DATA, not a label.

    Keying on anything smaller (a fold index, an arm name) would collide the moment the features
    change -- and silently serve a stale probability vector, which is worse than paying again.
    """
    # y's dtype is part of the bytes hashed, so it must stay STABLE per task: switching
    # classification from int64 to float64 would silently invalidate every cached call (measured:
    # 160 entries, 1.6M credits of rework) while changing nothing about the data.
    y_dtype = "float64" if task == "regression" else "int64"
    h = hashlib.sha256()
    for a in (np.ascontiguousarray(Xtr, dtype="float32"),
              np.ascontiguousarray(ytr, dtype=y_dtype),
              np.ascontiguousarray(Xte, dtype="float32")):
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    h.update(f"{model}|{seed}|{task}".encode())
    return h.hexdigest()[:32]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--model", default="auto", help="e.g. v3.5, v3.5-fast; 'auto' = the default")
    ap.add_argument("--n-estimators", default="auto")
    ap.add_argument("--hosted", action="store_true",
                    help="run on Prior Labs' API instead of local weights (uploads the data)")
    ap.add_argument("--mode", default="cv", choices=["cv", "fit_predict"])
    ap.add_argument("--cache-dir", default=None,
                    help="content-addressed predict cache; omit to disable (and pay every time)")
    ap.add_argument("--task", default="classification",
                    choices=["classification", "regression"],
                    help="regression uses TabPFNRegressor and returns point predictions")
    ap.add_argument("--credit-budget", type=int, default=None,
                    help="abort before exceeding this many credits in THIS run")
    args = ap.parse_args()

    if not os.environ.get("TABPFN_TOKEN"):
        sys.exit("FATAL TABPFN_TOKEN is not set. TabPFN requires a one-time licence acceptance "
                 "before it will download weights:\n"
                 "  1. log in at https://ux.priorlabs.ai\n"
                 "  2. accept the licence on the Licenses tab\n"
                 "  3. copy the API key from https://ux.priorlabs.ai/account\n"
                 "  4. export TABPFN_TOKEN=\"<key>\"\n"
                 "This is a licence gate, not a technical one -- do not route around it.")

    import importlib.metadata as md

    from sklearn.metrics import roc_auc_score

    credits_used = None
    if args.hosted:
        import re as _re

        import tabpfn_client
        from tabpfn_client import TabPFNClassifier, TabPFNRegressor
        tabpfn_client.set_access_token(os.environ["TABPFN_TOKEN"])
        version = f"client-{md.version('tabpfn-client')}"

        def credits_used():                                   # noqa: F811
            m = _re.search(r"used ([\d,]+) of", tabpfn_client.get_api_usage())
            return int(m.group(1).replace(",", "")) if m else None
    else:
        from tabpfn import TabPFNClassifier, TabPFNRegressor
        version = md.version("tabpfn")
    z = np.load(args.inp, allow_pickle=False)
    if args.mode == "cv":
        X, y, folds = z["X"].astype("float32"), z["y"].astype(int), z["folds"].astype(int)
        n_seeds, n = folds.shape
        n_folds = int(folds.max()) + 1
    else:
        X = y = folds = None
        n_seeds = n = n_folds = 0

    if args.mode == "cv":
        if X.shape[0] != n or len(y) != n:
            sys.exit(f"FATAL shape mismatch: X {X.shape}, y {len(y)}, folds {folds.shape}")
        if (folds < 0).any():
            sys.exit("FATAL the fold matrix holds negative entries -- some protein was never "
                     "assigned")

    if args.mode == "cv":
        print(f"tabpfn {version} | X {X.shape} | {int(y.sum())}/{n} positive "
              f"({100 * y.mean():.1f}%) | {n_seeds} seeds x {n_folds} folds | "
              f"device={args.device} model={args.model}", flush=True)
    else:
        print(f"tabpfn {version} | mode=fit_predict | device={args.device} model={args.model}",
              flush=True)

    # The hosted client takes neither `device` (it has its own) nor a string n_estimators.
    kw = {} if args.hosted else {"device": args.device, "n_estimators": args.n_estimators}
    if args.model != "auto":
        kw["model_path"] = args.model
    c0 = credits_used() if credits_used else None
    if c0 is not None:
        print(f"  credits before: {c0:,}", flush=True)

    cache_dir = Path(args.cache_dir) if args.cache_dir else None
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
    stats = {"calls": 0, "cache_hits": 0}

    def call(Xtr, ytr, Xte, seed: int) -> np.ndarray:
        """One fit+predict, served from the content-addressed cache when possible.

        The cache is checked BEFORE the budget: a cached answer costs nothing, so a run that is
        fully cached must never be refused for lack of credits.
        """
        if cache_dir:
            f = cache_dir / f"{_key(Xtr, ytr, Xte, args.model, seed, args.task)}.npy"
            if f.exists():
                stats["cache_hits"] += 1
                return np.load(f)
        if args.credit_budget is not None and credits_used:
            spent = credits_used() - (c0 or 0)
            if spent + 10_000 > args.credit_budget:
                sys.exit(f"FATAL credit budget exhausted: {spent:,} spent of "
                         f"{args.credit_budget:,} allowed, and the next call costs ~10,000. "
                         "Raise --credit-budget or accept the partial result; do NOT silently "
                         "continue spending.")
        # RETRY, because the hosted route fails transiently and a bare call makes the whole stage
        # hostage to it: a single GCS 500 ("We encountered an internal error. Please try again.")
        # killed a 40-minute run at call 64 of ~160. Server-side, nothing to do with our data, and
        # the next attempt succeeds. A failed attempt is not billed -- the charge is on a completed
        # predict -- so retrying costs time, not credits.
        last = None
        for attempt in range(MAX_ATTEMPTS):
            try:
                if args.task == "regression":
                    est = TabPFNRegressor(random_state=seed, **kw)
                    est.fit(Xtr, ytr)
                    p = np.asarray(est.predict(Xte), dtype=float)
                else:
                    est = TabPFNClassifier(random_state=seed, **kw)
                    est.fit(Xtr, ytr)
                    p = est.predict_proba(Xte)[:, 1]
                break
            except Exception as exc:                      # noqa: BLE001 -- any transport failure
                last = exc
                if attempt == MAX_ATTEMPTS - 1:
                    raise
                wait = BACKOFF_BASE * (2 ** attempt)
                print(f"    [retry {attempt + 1}/{MAX_ATTEMPTS - 1}] {type(exc).__name__}: "
                      f"{str(exc)[:160]} -- waiting {wait}s", flush=True)
                time.sleep(wait)
        stats["calls"] += 1
        if cache_dir:
            np.save(cache_dir / f"{_key(Xtr, ytr, Xte, args.model, seed, args.task)}.npy", p)
        return p

    # ---- fit_predict: the calls that are not k-fold (cross-activator, full-proteome scoring)
    if args.mode == "fit_predict":
        Xtr = z["X_train"].astype("float32")
        ytr = z["y_train"].astype(float if args.task == "regression" else int)
        Xte = z["X_test"].astype("float32")
        print(f"fit_predict: train {Xtr.shape} ({int(ytr.sum())} positive) -> test {Xte.shape}",
              flush=True)
        prob = call(Xtr, ytr, Xte, seed=0)
        c1 = credits_used() if credits_used else None
        spent = (c1 - c0) if (c0 is not None and c1 is not None) else -1
        print(f"  {args.tag} {stats['calls']} call(s), {stats['cache_hits']} cache hit(s), "
              f"spent {spent:,}", flush=True)
        np.savez_compressed(args.out, prob=prob, tabpfn_version=version,
                            device=("hosted" if args.hosted else args.device),
                            model_version=args.model, credits_spent=spent,
                            n_calls=stats["calls"], n_cache_hits=stats["cache_hits"])
        print(f"  {args.tag} wrote {args.out}", flush=True)
        return

    oof = np.full((n_seeds, n), np.nan)
    secs = np.zeros((n_seeds, n_folds))

    for si in range(n_seeds):
        for k in range(n_folds):
            te = folds[si] == k
            tr = ~te
            t0 = time.time()
            # random_state varies the model's internal ensembling, not the split -- the split is
            # fixed by the caller. Tied to the seed so the arms stay comparable.
            oof[si, te] = call(X[tr], y[tr], X[te], seed=si)
            secs[si, k] = time.time() - t0
            print(f"  {args.tag} seed {si} fold {k}: n_train={int(tr.sum())} "
                  f"n_test={int(te.sum())} {secs[si, k]:.0f}s", flush=True)
        auc_si = roc_auc_score(y, oof[si])
        eta = secs[secs > 0].mean() * (n_seeds - si - 1) * n_folds / 60
        print(f"  {args.tag} SEED {si} DONE ({si + 1}/{n_seeds}): held-out AUROC {auc_si:.4f}  "
              f"[running mean "
              f"{np.mean([roc_auc_score(y, oof[j]) for j in range(si + 1)]):.4f}]  "
              f"~{eta:.0f} min left in this arm", flush=True)

    if np.isnan(oof).any():
        sys.exit(f"FATAL {int(np.isnan(oof).sum())} predictions are NaN")

    per_seed = [roc_auc_score(y, oof[si]) for si in range(n_seeds)]
    print(f"  {args.tag} per-seed held-out AUROC: {', '.join(f'{a:.4f}' for a in per_seed)}",
          flush=True)
    c1 = credits_used() if credits_used else None
    spent = (c1 - c0) if (c0 is not None and c1 is not None) else -1
    if spent >= 0:
        print(f"  {args.tag} credits: {c0:,} -> {c1:,} (spent {spent:,}) | "
              f"{stats['calls']} call(s), {stats['cache_hits']} cache hit(s))", flush=True)
    np.savez_compressed(args.out, oof=oof, fit_seconds=secs, tabpfn_version=version,
                        device=("hosted" if args.hosted else args.device),
                        model_version=args.model, credits_spent=spent,
                        n_calls=stats["calls"], n_cache_hits=stats["cache_hits"])
    print(f"  {args.tag} wrote {args.out} ({secs.sum() / 60:.1f} min of fitting)", flush=True)


if __name__ == "__main__":
    main()
