"""Does fusing ESM-C and ProtT5 beat either alone?  FREE: reads cached out-of-fold TabPFN predictions.

Reads   data/processed/degradability/scratch/lazyqsar/{bundle,oof_tabpfn}_<activator>_<esmc|prott5>.npz
        output/results/degradability/proteomelm_vs_esmc.tsv      (the grid the arms are checked against)
Writes  output/results/degradability/fusion_eval.tsv

NO TabPFN CALL IS MADE and no credit is spent -- `head_comparison.py` already fitted both single
arms on identical cluster-grouped folds and kept the per-seed out-of-fold probabilities.

THE FUSION IS LATE: the mean of the two arms' out-of-fold probabilities, per seed. Concatenating
ESM-C + ProtT5 is 2,176 features, above the installed `tabpfn` limit of 2,000, and would need a
reduction step chosen for the purpose.

THREE GUARDS before any number is trusted:
  1. the two arms' labels `y` and fold matrices are byte-identical, or the pairing is meaningless;
  2. each single arm REPRODUCES its row in `proteomelm_vs_esmc.tsv` (roc_auc and pr_auc to 1e-4);
  3. every out-of-fold probability is finite and in [0, 1].

HOW TO READ IT (house rules): PR-AUC beside ROC-AUC, never ROC alone; a PAIRED comparison -- per-seed
differences on identical folds, mean +/- SE, and `[k/5 seeds]` won. Two independent SDs are badly
under-powered here because partition difficulty is shared.

Run with the `gradi` env:  python scripts/degradability/fusion_eval.py [-q]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from sklearn.metrics import roc_auc_score  # noqa: E402

from src import degradability as D  # noqa: E402

WORK = REPO_ROOT / "data/processed/degradability/scratch/lazyqsar"
GRID = REPO_ROOT / "output/results/degradability/proteomelm_vs_esmc.tsv"
OUT = REPO_ROOT / "output/results/degradability/fusion_eval.tsv"
ACTIVATORS = ("adep4", "onc212")
ARMS = ("esmc", "prott5")
TOL = 1e-4


def say(msg: str = "", quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def load_arm(act: str, fs: str):
    b = np.load(WORK / f"bundle_{act}_{fs}.npz", allow_pickle=False)
    o = np.load(WORK / f"oof_tabpfn_{act}_{fs}.npz", allow_pickle=False)
    return b["y"].astype(int), b["folds"], o["oof"]


def per_seed(y: np.ndarray, oof: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (np.array([roc_auc_score(y, p) for p in oof]),
            np.array([D.pr_auc(y, p) for p in oof]))


def paired(a: np.ndarray, b: np.ndarray) -> str:
    """a - b per seed: mean +/- SE and how many seeds a wins."""
    d = a - b
    se = d.std(ddof=1) / np.sqrt(len(d))
    return f"{d.mean():+.4f} +/- {se:.4f}  [{int((d > 0).sum())}/{len(d)} seeds]"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    q = ap.parse_args().quiet

    say("fusion_eval.py  ESM-C + ProtT5 late fusion, from cached out-of-fold TabPFN predictions", q)
    say(f"  inputs : {WORK.relative_to(REPO_ROOT)}/ (cached, read-only)   spends: 0 credits", q)
    say(f"  output : {OUT.relative_to(REPO_ROOT)}", q)
    grid = pd.read_csv(GRID, sep="\t")
    rows, failures = [], []

    for act in ACTIVATORS:
        data = {fs: load_arm(act, fs) for fs in ARMS}
        y, folds, _ = data["esmc"]
        # guard 1: the pairing
        if not (np.array_equal(y, data["prott5"][0]) and np.array_equal(folds, data["prott5"][1])):
            sys.exit(f"FATAL {act}: the arms do not share labels and folds; a paired test is meaningless")
        oofs = {fs: data[fs][2] for fs in ARMS}
        # guard 3
        for fs, o in oofs.items():
            if not (np.isfinite(o).all() and (o >= 0).all() and (o <= 1).all()):
                sys.exit(f"FATAL {act} {fs}: out-of-fold probabilities are not finite probabilities")
        oofs["fused"] = (oofs["esmc"] + oofs["prott5"]) / 2

        say(f"\n  {act.upper()}: n={len(y):,}  positives={int(y.sum())} ({100 * y.mean():.1f}%)  "
            f"seeds={oofs['esmc'].shape[0]}  folds={int(folds.max()) + 1}", q)
        stats = {k: per_seed(y, o) for k, o in oofs.items()}

        # guard 2: reproduce the grid
        for fs in ARMS:
            g = grid[(grid.activator == act) & (grid.estimator == "tabpfn") & (grid.feature_set == fs)].iloc[0]
            roc, pr = stats[fs][0].mean(), stats[fs][1].mean()
            ok = abs(roc - g.roc_auc) < TOL and abs(pr - g.pr_auc) < TOL
            say(f"    reproduces grid {fs:<7} roc {roc:.4f} vs {g.roc_auc:.4f}   pr {pr:.4f} vs {g.pr_auc:.4f}"
                f"   {'OK' if ok else 'MISMATCH'}", q)
            if not ok:
                failures.append(f"{act}/{fs}")

        say("", q)
        say(f"    {'arm':<8}{'ROC-AUC':>16}{'PR-AUC':>16}", q)
        for k in ("esmc", "prott5", "fused"):
            r, p = stats[k]
            say(f"    {k:<8}{r.mean():>9.4f} +/-{r.std():.4f}{p.mean():>9.4f} +/-{p.std():.4f}", q)
            rows.append({"activator": act, "arm": k, "n": len(y), "n_pos": int(y.sum()),
                         "roc_auc": round(r.mean(), 4), "roc_auc_sd": round(r.std(), 4),
                         "pr_auc": round(p.mean(), 4), "pr_auc_sd": round(p.std(), 4),
                         "roc_auc_seeds": ";".join(f"{x:.4f}" for x in r),
                         "pr_auc_seeds": ";".join(f"{x:.4f}" for x in p)})

        best = max(("esmc", "prott5"), key=lambda k: stats[k][1].mean())  # best single by PR-AUC
        say(f"\n    PAIRED, fused minus ...   (best single arm by PR-AUC: {best})", q)
        for ref in ("esmc", "prott5"):
            say(f"      vs {ref:<7} ROC {paired(stats['fused'][0], stats[ref][0])}", q)
            say(f"      {'':<10} PR  {paired(stats['fused'][1], stats[ref][1])}", q)
        rho = np.corrcoef(oofs["esmc"].mean(0), oofs["prott5"].mean(0))[0, 1]
        say(f"    correlation of the two arms' seed-averaged OOF probabilities: {rho:.3f}", q)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    if failures:
        print(f"\nFAIL: arms did not reproduce the grid: {failures}", file=sys.stderr)
        return 1
    say("\n  guards: pairing OK, grid reproduced, probabilities finite -- the numbers above can be trusted", q)
    return 0


if __name__ == "__main__":
    sys.exit(main())
