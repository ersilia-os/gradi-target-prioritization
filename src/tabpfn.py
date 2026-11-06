"""TabPFN-3.5, the project's default learner, callable from any stage.

CLAUDE.md, *Supervised ML: always TabPFN-3.5*: whenever this project needs a classifier or a
regressor, it uses TabPFN-3.5. That makes it transversal, so the dispatch machinery lives here
rather than inside the stage that happened to need it first (degradability). Stages call
`predict_fold` / `predict_fold_regression` / `oof_predict`; nothing else should shell out to the
worker directly.

**It cannot run in `gradi`.** TabPFN pulls torch 2.14 + mlx against gradi's 2.12, which would take
stage 01's ESM-C down with it. So every call is an out-of-process dispatch to
`scripts/workers/tabpfn_cv.py` under the `gradi-tabpfn` interpreter. That is why this module is a
subprocess wrapper and not a scikit-learn estimator.

**Hosted TabPFN bills ~10,000 credits per CALL, flat** -- measured at 350, 1,000 and 2,889 test
rows, all the same price; `fit` is free, because fit only stores the rows and predict is the single
forward pass. The consequence is that cost scales with the NUMBER of calls, not their size: batch
your test rows into one call wherever the question allows it.

Everything is served from a content-addressed cache keyed on the data itself, so a re-run that
recomputes nothing also spends nothing. `cache_key` and `is_cached` let a caller ask "would this
cost anything?" without spending anything -- which is the ONLY correct way to check the cache. The
key function is imported from the worker rather than copied, so the two cannot drift apart; a
drift would not raise, it would just silently miss on every lookup.

Licence: TabPFN's weights are **non-commercial**. That is a question for the collaboration, not a
thing this module can decide.
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKER = REPO_ROOT / "scripts" / "workers" / "tabpfn_cv.py"
DEFAULT_BIN = Path.home() / "miniconda3" / "envs" / "gradi-tabpfn" / "bin" / "python"

# Shared across axes, deliberately. The cache is content-addressed, so two stages asking the same
# question of the same data SHOULD hit the same entry -- that is the saving. Filenames are the keys,
# so this directory can be moved without invalidating anything.
CACHE_DIR = REPO_ROOT / "data" / "processed" / "tabpfn" / "cache"

MODEL = "v3.5"
CREDITS_PER_CALL = 10_000     # measured, flat, independent of test-set size


def _worker_module():
    """The worker, imported for its pure helpers (`_key`). Never for its TabPFN dependencies."""
    spec = importlib.util.spec_from_file_location("tabpfn_worker", WORKER)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def cache_key(X_train, y_train, X_test, *, model: str = "auto", seed: int = 0,
              task: str = "classification") -> str:
    """The content address this call would use. Computing it costs nothing.

    Imported from the worker so the two definitions cannot diverge -- see the CACHE KEY FORMAT
    block in `scripts/workers/tabpfn_cv.py` before changing anything it hashes.
    """
    return _worker_module()._key(X_train, y_train, X_test, model, seed, task)


def is_cached(X_train, y_train, X_test, *, model: str = "auto", seed: int = 0,
              task: str = "classification", cache_dir: Path | None = None) -> bool:
    """Would this call be free? Checks for the key's file; makes no API call."""
    d = Path(cache_dir) if cache_dir else CACHE_DIR
    return (d / f"{cache_key(X_train, y_train, X_test, model=model, seed=seed, task=task)}.npy").exists()


def binary() -> str:
    """The `gradi-tabpfn` interpreter, plus the licence gate, checked before any expensive work."""
    p = Path(os.environ.get("GRADI_TABPFN_BIN", DEFAULT_BIN))
    if not p.exists():
        sys.exit(f"FATAL no tabpfn interpreter at {p}. See install.sh. Do NOT install tabpfn into "
                 "`gradi`: it pulls torch 2.14 and would break stage 01's ESM-C.")
    if not os.environ.get("TABPFN_TOKEN"):
        sys.exit("FATAL TABPFN_TOKEN is not set. TabPFN needs a one-time licence acceptance at "
                 "https://ux.priorlabs.ai before it will serve predictions.")
    return str(p)


def _run(X_train, y_train, X_test, *, task: str, seed: int, tag: str,
         cache_dir: Path | None, credit_budget: int | None) -> np.ndarray:
    d = Path(cache_dir) if cache_dir else CACHE_DIR
    d.mkdir(parents=True, exist_ok=True)

    # `fit_predict` in the worker always keys on seed 0 and model "auto", so this lookup reproduces
    # the worker's own and lets the whole call be skipped when it would be free anyway.
    hit = d / f"{cache_key(X_train, y_train, X_test, task=task)}.npy"
    if hit.exists():
        return np.load(hit)

    # GRADI_TABPFN_CACHE_ONLY=1 -- refuse to spend. The point is to be able to prove a run is fully
    # cached WITHOUT the proof itself costing 10,000 credits a call. Use it after any refactor that
    # could have disturbed what the key hashes (see CACHE KEY FORMAT in the worker): if the stage
    # still reproduces its numbers under this flag, the cache survived.
    if os.environ.get("GRADI_TABPFN_CACHE_ONLY") == "1":
        sys.exit(f"FATAL cache-only mode: {tag} is NOT cached ({hit.name}). This call would cost "
                 f"~{CREDITS_PER_CALL:,} credits. Unset GRADI_TABPFN_CACHE_ONLY to pay for it.")
    y_dtype = float if task == "regression" else int
    with tempfile.TemporaryDirectory() as td:
        bundle, out = Path(td) / "in.npz", Path(td) / "out.npz"
        np.savez_compressed(bundle,
                            X_train=np.asarray(X_train, dtype="float32"),
                            y_train=np.asarray(y_train, dtype=y_dtype),
                            X_test=np.asarray(X_test, dtype="float32"))
        cmd = [binary(), str(WORKER), "--mode", "fit_predict", "--task", task,
               "--in", str(bundle), "--out", str(out), "--tag", tag, "--cache-dir", str(d)]
        if credit_budget is not None:
            cmd += ["--credit-budget", str(credit_budget)]
        if os.environ.get("GRADI_TABPFN_HOSTED", "1") != "0":
            cmd.append("--hosted")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or not out.exists():
            sys.exit(f"FATAL the tabpfn worker failed for {tag} (exit {r.returncode}):\n"
                     f"{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        return np.load(out, allow_pickle=False)["prob"]


def predict_fold(X_train, y_train, X_test, *, seed: int = 0, tag: str = "",
                 cache_dir: Path | None = None, credit_budget: int | None = None) -> np.ndarray:
    """Fit on (X_train, y_train); return P(class 1) for X_test."""
    return _run(X_train, y_train, X_test, task="classification", seed=seed, tag=tag,
                cache_dir=cache_dir, credit_budget=credit_budget)


def predict_fold_regression(X_train, y_train, X_test, *, seed: int = 0, tag: str = "",
                            cache_dir: Path | None = None,
                            credit_budget: int | None = None) -> np.ndarray:
    """The same on a continuous target.

    Separate from `predict_fold` because the cache key includes the task: a regression call and a
    classification call on the same arrays are different questions and must not share an entry.
    """
    return _run(X_train, y_train, X_test, task="regression", seed=seed, tag=tag,
                cache_dir=cache_dir, credit_budget=credit_budget)


def oof_predict(X, y, splits, *, seed: int = 0, tag: str = "", task: str = "classification",
                cache_dir: Path | None = None, credit_budget: int | None = None) -> np.ndarray:
    """Out-of-fold predictions over a MATERIALISED split list.

    Splits are passed in rather than generated here so that several estimators can be compared on
    byte-identical folds -- the difference between two heads is smaller than the difference between
    two partitions, so a comparison on differently-generated folds measures the wrong thing.
    """
    pred = np.full(len(y), np.nan)
    for k, (tr, te) in enumerate(splits):
        pred[te] = _run(X[tr], y[tr], X[te], task=task, seed=seed, tag=f"{tag}/fold{k}",
                        cache_dir=cache_dir, credit_budget=credit_budget)
    if np.isnan(pred).any():
        sys.exit(f"FATAL {tag}: {int(np.isnan(pred).sum())} rows were in no test fold")
    return pred


def would_cost(calls, *, cache_dir: Path | None = None) -> tuple[int, int]:
    """(n_uncached, credits) for an iterable of (X_train, y_train, X_test[, task]) tuples.

    Ask this BEFORE a long run. Costs nothing to answer.
    """
    n = 0
    for c in calls:
        Xtr, ytr, Xte = c[0], c[1], c[2]
        task = c[3] if len(c) > 3 else "classification"
        if not is_cached(Xtr, ytr, Xte, task=task, cache_dir=cache_dir):
            n += 1
    return n, n * CREDITS_PER_CALL
