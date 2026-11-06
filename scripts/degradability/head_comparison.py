"""ProteomeLM vs ESM-C, and lazy-qsar vs the stage-04 forest, for activated-ClpP degradability.

Stage 04 predicts whether a protein is a substrate of activated partnerless ClpP, from ESM-C 600M
embeddings, trained on the two S. aureus activator screens. Measured there, under cluster-grouped CV
repeated over 5 seeds: **ADEP4 0.865 +/- 0.002** against a length-only baseline of 0.776, and
**ONC212 0.752 +/- 0.007** against 0.687.

This script varies the two things stage 04 held fixed, one at a time and then together:

  FEATURES    three protein language models over the same proteins --
                esmc        ESM-C 600M, 1152-d, per-protein (stage 01, what stage 04 ships)
                proteomelm  ProteomeLM-L layer 8, 1152-d, proteome-CONTEXTUALISED (stage 07)
                prott5      ProtT5-XL-U50, 1024-d, per-protein (embeddings/prott5.py), the one
                            arm validated against an external reference: our E. coli vectors
                            reproduce UniProt's published ProtT5 at median cosine 1.000000
  ESTIMATOR   three heads on identical folds --
                forest      stage 04's hand-set RandomForest(500) -- the FORMER shipped
                            configuration, kept as the baseline this comparison is measured against
                            and regenerable with `degradability/predict.py --estimator forest`
                lazyqsar    lazy-qsar's LazyClassifier: an lr/xgb/rf/svc portfolio with internal
                            calibration and a gating pooler, selected inside each training fold
                tabpfn      TabPFN-3.5, a tabular FOUNDATION MODEL -- pre-trained on synthetic
                            tasks, learns in context in one forward pass, nothing tuned. **This
                            comparison is why it is now stage 04's shipped estimator.** Its weights
                            are non-commercial licensed, which is an open question for the
                            deliverable, not a blocker for methods work.

WHAT IS HELD FIXED
------------------
The labels, the homology clusters, the fold count, the seeds -- and, critically, **the folds
themselves**: the split list is materialised ONCE per activator per seed and handed to both
estimators, so the two arms are scored on byte-identical partitions rather than on two independent
draws that merely used the same seed. The estimator, splitters and fold-disjointness guard are
IMPORTED FROM `scripts/degradability/predict.py`, so the forest arm cannot drift from the shipped stage.

Every comparison is therefore PAIRED, and the seed SD is a fair error bar on the difference.

READ THE SD, NOT THE MEAN ALONE
-------------------------------
A single-seed estimate carries ~+/-0.004 of pure arbitrariness here, larger than several effects
stage 04 tested and rejected. Every number is mean +/- SD over `N_CV_SEEDS` seeds, and a difference
smaller than the SDs is not a result.

TWO NUMBERS THAT LOOK COMPARABLE AND ARE NOT
--------------------------------------------
1. The `length` row is scored with the SAME estimator as everything else. Stage 04 measured that a
   forest is the wrong estimator for a one-feature baseline (500 trees bin one variable into steps):
   RF 0.682 vs logistic 0.776 on ADEP4. So the printed `length` number reads ~0.09 LOW and is not
   the honest length baseline -- stage 04's logistic 0.776 / 0.687 is. It is kept because holding
   the estimator fixed is what makes the real arms comparable.
2. lazy-qsar's own `oof_auc_` is reported in the audit and **must not be quoted as a result**. It
   comes from its internal ungrouped splits, so homologs sit on both sides of them -- the exact
   leakage the cluster-grouped folds exist to prevent. The worker prints both so the gap is visible.

Outputs
  output/results/degradability/proteomelm_vs_esmc.tsv          the result table
  data/processed/degradability/scratch/lazyqsar/             per-arm OOF probabilities + audit

Run with the `gradi` env; the lazy-qsar arm is dispatched to `gradi-lazyqsar` across a process
boundary (`GRADI_LAZYQSAR_BIN` overrides). ~10 min for the forest arm alone; ~2.5 h with lazy-qsar.
  python scripts/degradability/head_comparison.py --estimator forest
  python scripts/degradability/head_comparison.py                      # both
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402
from src import embeddings as E  # noqa: E402
from src import proteomelm as M  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "results" / "degradability"
WORK_DIR = REPO_ROOT / "data" / "processed" / "degradability" / "scratch" / "lazyqsar"
WORKER = REPO_ROOT / "scripts" / "degradability" / "workers" / "lazyqsar_cv.py"
TABPFN_WORKER = REPO_ROOT / "scripts" / "workers" / "tabpfn_cv.py"
DEFAULT_LAZYQSAR_BIN = Path.home() / "miniconda3" / "envs" / "gradi-lazyqsar" / "bin" / "python"
DEFAULT_TABPFN_BIN = Path.home() / "miniconda3" / "envs" / "gradi-tabpfn" / "bin" / "python"

SPECIES = "saureus"
FEATURE_SETS = ("length", "esmc", "proteomelm", "prott5", "both")
DEFAULT_FEATURE_SETS = ("esmc", "proteomelm", "prott5")
ESTIMATORS = ("forest", "lazyqsar", "tabpfn")
DEFAULT_ESTIMATORS = ("forest", "lazyqsar")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 100) -> None:
    say(char * width)


def tabpfn_bin() -> tuple[str, str]:
    """The `gradi-tabpfn` interpreter, and the licence gate checked BEFORE any expensive work."""
    p = Path(os.environ.get("GRADI_TABPFN_BIN", DEFAULT_TABPFN_BIN))
    if not p.exists():
        sys.exit(f"FATAL no tabpfn interpreter at {p}.\n"
                 "  conda create -y -n gradi-tabpfn python=3.11\n"
                 "  ~/miniconda3/envs/gradi-tabpfn/bin/pip install 'tabpfn==9.0.0'")
    if not os.environ.get("TABPFN_TOKEN"):
        sys.exit("FATAL TABPFN_TOKEN is not set -- TabPFN needs a one-time licence acceptance "
                 "before it downloads weights. Accept at https://ux.priorlabs.ai (Licenses tab), "
                 "copy the key from /account, then `export TABPFN_TOKEN=...`.")
    hosted = os.environ.get("GRADI_TABPFN_HOSTED", "1") != "0"
    mod = "tabpfn_client" if hosted else "tabpfn"
    pkg = "tabpfn-client" if hosted else "tabpfn"
    out = subprocess.run([str(p), "-c", f"import {mod}; import importlib.metadata as m;"
                          f"print(m.version('{pkg}'))"], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"FATAL {p} cannot import {mod}:\n{out.stderr.strip()}")
    return str(p), ("hosted " if hosted else "local ") + out.stdout.strip()


def lazyqsar_bin() -> str:
    """The `gradi-lazyqsar` interpreter. Checked before any expensive work, not at first use."""
    p = Path(os.environ.get("GRADI_LAZYQSAR_BIN", DEFAULT_LAZYQSAR_BIN))
    if not p.exists():
        sys.exit(f"FATAL no lazy-qsar interpreter at {p}. Create it -- and do NOT install lazy-qsar "
                 "into `gradi`: it pins numpy==2.1.3 / scikit-learn==1.6.1 against gradi's 2.4.6 / "
                 "1.9.0 under torch, which would break stages 01 and 04.\n"
                 "  conda create -y -n gradi-lazyqsar python=3.11\n"
                 "  git clone --branch v3.4.4 https://github.com/ersilia-os/lazy-qsar && "
                 "~/miniconda3/envs/gradi-lazyqsar/bin/pip install './lazy-qsar[fit]'")
    out = subprocess.run([str(p), "-c", "import importlib.metadata as m;"
                          "print(m.version('lazyqsar'))"], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"FATAL {p} cannot import lazyqsar:\n{out.stderr.strip()}")
    return str(p), out.stdout.strip()


def stage04():
    """Import scripts/degradability/predict.py as a module, for its estimator and splitters.

    Deliberate: reimplementing `_pipe`/`_splitters`/`_assert_folds_disjoint` here would let this
    comparison silently drift from the thing it is comparing against. `scripts/` is not a package,
    hence the explicit spec load.
    """
    spec = importlib.util.spec_from_file_location(
        "stage04", REPO_ROOT / "scripts" / "degradability/predict.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.VERBOSE = False              # its helpers print through its own `say`; we print our own
    for name in ("_pipe", "_splitters", "_assert_folds_disjoint"):
        if not hasattr(mod, name):
            sys.exit(f"FATAL scripts/degradability/predict.py has no {name} -- it was refactored; "
                     "this comparison must be re-aligned rather than guessed at")
    return mod


def build_features() -> pd.DataFrame:
    """One frame indexed by uniprot_ac carrying every candidate feature, for S. aureus.

    `log_length` and the ESM-C block come from stage 04's own recipe; the ProteomeLM block is
    stage 07's layer-8 output. Joined INNER, so every feature set is scored on exactly the same
    proteins -- otherwise a difference could just be a difference in who was included.
    """
    prot = P.load(SPECIES)[["uniprot_ac", "sequence"]].set_index("uniprot_ac")
    feats = D.sequence_features(prot["sequence"])

    accs_e, mat_e = E.load(SPECIES)
    esmc = pd.DataFrame(mat_e, index=pd.Index(accs_e, name="uniprot_ac"),
                        columns=[f"e{i}" for i in range(mat_e.shape[1])])

    accs_p, mat_p = M.load(SPECIES)
    plm = pd.DataFrame(mat_p, index=pd.Index(accs_p, name="uniprot_ac"),
                       columns=[f"m{i}" for i in range(mat_p.shape[1])])

    t5 = E.load_prott5_frame(SPECIES)

    md = M.metadata(SPECIES)
    t5md = E.prott5_metadata(SPECIES)
    say(f"  esm-c      {mat_e.shape} from stage 01")
    say(f"  proteomelm {mat_p.shape} from stage 07: {md['model'].split('/')[-1]} "
        f"layer {md['layer']}/{md['n_layers']}, group_embeds={md['group_embeds_mode']}")
    say(f"  prott5     {t5.shape} from 01_embeddings_prott5: {t5md['model'].split('/')[-1]}, "
        f"pooling={t5md['pooling']}, {t5md['n_chunked']} windowed (>{t5md['max_len']} aa)")

    out = feats.join(esmc, how="inner").join(plm, how="inner").join(t5, how="inner")
    say(f"  joined     {len(out):,} proteins carry all three (of {len(prot):,} in the proteome)")
    return out


def columns_for(name: str, frame: pd.DataFrame) -> list[str]:
    esmc = [c for c in frame.columns if c.startswith("e") and c[1:].isdigit()]
    plm = [c for c in frame.columns if c.startswith("m") and c[1:].isdigit()]
    t5 = [c for c in frame.columns if c.startswith("t") and c[1:].isdigit()]
    return {"length": ["log_length"], "esmc": esmc, "proteomelm": plm, "prott5": t5,
            "both": esmc + plm}[name]


def make_folds(s04, X, y, groups, activator: str, folds: int, seed: int):
    """Materialise the split list ONCE per seed, so every estimator is scored on the same partition.

    Two estimators each calling `cv.split()` would give the same answer today -- StratifiedGroupKFold
    with a fixed `random_state` is deterministic -- but nothing enforces that, and a paired
    comparison whose pairing rests on an unstated invariant is one refactor from being unpaired.
    Materialising it makes the guarantee explicit and lets the fold assignment be handed to a worker
    in another process.
    """
    splits_per_seed, fold_matrix = [], []
    for si in range(D.N_CV_SEEDS):
        cv = dict(s04._splitters(folds, seed + si))["clustered"]
        # The guarantee behind every number here: no homology cluster spans two folds.
        s04._assert_folds_disjoint(cv, X, y, groups, f"{activator}", seed + si)
        splits = list(cv.split(X, y, groups=groups))
        assign = np.full(len(y), -1, dtype=int)
        for k, (_, te) in enumerate(splits):
            assign[te] = k
        if (assign < 0).any():
            sys.exit(f"FATAL {activator} seed {seed + si}: {int((assign < 0).sum())} proteins are "
                     "in no test fold -- the out-of-fold score would not cover the labeled set")
        splits_per_seed.append(splits)
        fold_matrix.append(assign)
    return splits_per_seed, np.vstack(fold_matrix)


def forest_oof(s04, X, y, splits_per_seed, seed: int) -> list[np.ndarray]:
    """Stage 04's forest, one out-of-fold probability vector per seed."""
    from sklearn.model_selection import cross_val_predict

    return [cross_val_predict(s04._pipe(seed), X, y, cv=splits, method="predict_proba")[:, 1]
            for splits in splits_per_seed]


def worker_oof(est: str, binary: str, X, y, fold_matrix, tag: str, refresh: bool):
    """An out-of-process estimator on the same folds. Cached -- a full grid is hours of fitting.

    `est` selects the worker; both speak the same bundle protocol (X, y, folds in; oof out), which
    is what keeps every head on byte-identical partitions.
    """
    worker = {"lazyqsar": WORKER, "tabpfn": TABPFN_WORKER}[est]
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    bundle = WORK_DIR / f"bundle_{tag}.npz"
    result = WORK_DIR / f"oof_{est}_{tag}.npz" if est != "lazyqsar" else WORK_DIR / f"oof_{tag}.npz"

    if result.exists() and not refresh:
        z = np.load(result, allow_pickle=False)
        oof = z["oof"]
        if oof.shape == fold_matrix.shape:
            say(f"    cached {result.name} ({float(z['fit_seconds'].sum()) / 60:.1f} min of "
                "fitting)")
            return [oof[si] for si in range(oof.shape[0])], z
        say(f"    cached {result.name} has shape {oof.shape}, need {fold_matrix.shape} -- refitting")

    np.savez_compressed(bundle, X=X.astype("float32"), y=y.astype(int), folds=fold_matrix)
    say(f"    dispatching to gradi-{est}: {X.shape} x {fold_matrix.shape[0]} seeds "
        f"x {int(fold_matrix.max()) + 1} folds")
    cmd = [binary, str(worker), "--in", str(bundle), "--out", str(result), "--tag", tag]
    if est == "tabpfn" and os.environ.get("GRADI_TABPFN_HOSTED", "1") != "0":
        cmd.append("--hosted")      # local weights need the licence ACCEPTED, not just a token
    proc = subprocess.run(cmd, text=True)
    if proc.returncode != 0 or not result.exists():
        sys.exit(f"FATAL the {est} worker failed for {tag} (exit {proc.returncode}). "
                 "Its output is above; this is not something to work around silently.")
    z = np.load(result, allow_pickle=False)
    oof = z["oof"]
    return [oof[si] for si in range(oof.shape[0])], z


def score(y, probs: list[np.ndarray]) -> dict:
    """Mean +/- SD over seeds, plus the bootstrap CI of the seed-averaged probability.

    PER-SEED values are kept, not just the summary. They are what a PAIRED comparison needs: every
    arm ran on byte-identical folds, so most of the seed-to-seed variation is shared partition
    difficulty rather than head-specific noise. Comparing two arms' independent SDs discards the
    pairing and is badly under-powered when that shared variance is large -- measured on ONC212,
    where both arms carry SD ~0.0075 because seed 1 is simply an easy split for every head.
    """
    from sklearn.metrics import roc_auc_score

    aucs = [roc_auc_score(y, p) for p in probs]
    prs = [D.pr_auc(y, p) for p in probs]
    mean_prob = np.mean(probs, axis=0)
    st = D.auroc_ci(y, mean_prob)
    return {
        "roc_auc": round(float(np.mean(aucs)), 4),
        "roc_auc_sd": round(float(np.std(aucs)), 4),
        "roc_auc_lo": st["roc_auc_lo"], "roc_auc_hi": st["roc_auc_hi"],
        "pr_auc": round(float(np.mean(prs)), 4),
        "pr_auc_sd": round(float(np.std(prs)), 4),
        "roc_auc_seeds": ";".join(f"{a:.4f}" for a in aucs),
        "pr_auc_seeds": ";".join(f"{a:.4f}" for a in prs),
    }, mean_prob


def evaluate(s04, binary, activator, lab, feats, estimators, feature_sets,
             folds, seed, refresh) -> tuple[list[dict], dict]:
    """Every (estimator, feature set) on identical folds. Returns (rows, seed-averaged OOF)."""
    sub = lab[lab[activator].notna()].copy()
    sub = sub[sub["uniprot_ac"].isin(feats.index)]
    X_all = feats.loc[sub["uniprot_ac"]]
    y = sub[activator].to_numpy(dtype=float).astype(int)
    groups = sub["cluster"].to_numpy()
    say(f"  {activator:<8} {len(sub):>5} labeled, {pd.Series(groups).nunique():>5} clusters, "
        f"{int(y.sum()):>4} positives ({100 * y.mean():.1f}%)")

    # StratifiedGroupKFold splits on y and groups only, so a placeholder of the right length is
    # enough -- materialising the 2,304-column matrix here just to be ignored would be a big copy.
    splits_per_seed, fold_matrix = make_folds(
        s04, np.zeros((len(y), 1)), y, groups, activator, folds, seed)

    rows, oof = [], {}
    for est in estimators:
        for fs in feature_sets:
            cols = columns_for(fs, X_all)
            X = X_all[cols].to_numpy(dtype=float)
            say(f"    {est:<9} {fs:<11} {len(cols):>5} feats")
            audit = {}
            if est == "forest":
                probs = forest_oof(s04, X, y, splits_per_seed, seed)
            else:
                probs, z = worker_oof(est, binary[est], X, y, fold_matrix,
                                      f"{activator}_{fs}", refresh)
                audit = {"fit_minutes": round(float(z["fit_seconds"].sum()) / 60, 1)}
                if est == "lazyqsar":
                    audit["lazyqsar_version"] = str(z["lazyqsar_version"])
                    # lazy-qsar's own ungrouped estimate. LEAKY -- kept only to show the gap.
                    audit["internal_oof_auc_leaky"] = round(
                        float(np.nanmean(z["internal_oof_auc"])), 4)
                    audit["portfolio"] = str(pd.Series(z["portfolio"].ravel()).mode().iat[0])
                else:
                    audit["tabpfn_version"] = str(z["tabpfn_version"])
                    audit["model_version"] = str(z["model_version"])
            st, mean_prob = score(y, probs)
            oof[(est, fs)] = mean_prob
            rows.append({"activator": activator, "estimator": est, "feature_set": fs,
                         "n_features": len(cols), "n": len(sub), "n_pos": int(y.sum()),
                         "base_rate": round(float(y.mean()), 4),
                         "n_cv_seeds": D.N_CV_SEEDS, **st, **audit})
            say(f"      AUROC {st['roc_auc']:.4f} +/- {st['roc_auc_sd']:.4f}   "
                f"PR-AUC {st['pr_auc']:.4f} +/- {st['pr_auc_sd']:.4f}"
                + (f"   [lazy-qsar's own leaky oof_auc_ {audit['internal_oof_auc_leaky']:.4f}, "
                   f"portfolio {audit['portfolio']}]"
                   if "internal_oof_auc_leaky" in audit else ""))
    return rows, oof


def redundancy(oof: dict, activators, estimators) -> None:
    """Do ESM-C and ProteomeLM rank the same proteins, or agree only on average?

    This is the number that decides whether the two are worth carrying separately. Equal AUROCs are
    compatible with either answer: two models can score identically while disagreeing protein by
    protein (in which case concatenating them should have helped a lot, and combining their outputs
    would be real ensembling) or while producing near-identical rankings (in which case ProteomeLM's
    proteome context is re-deriving what ESM-C already encodes, and a second column buys nothing).

    Stage 04 raised exactly this trap for the two ACTIVATORS -- predictions at rho 0.82-0.84 while
    the labels agree at only 0.52, so "both activators agree" is closer to one opinion than two.
    Same test, different axis.
    """
    from scipy.stats import spearmanr

    rule()
    say("REDUNDANCY - do the arms rank the same proteins?")
    rule()
    import itertools

    for act in activators:
        o = oof[act]
        for est in estimators:
            arms = [f for f in ("esmc", "proteomelm", "prott5") if (est, f) in o]
            for a, b in itertools.combinations(arms, 2):
                r = spearmanr(o[(est, a)], o[(est, b)]).statistic
                say(f"  {act:<8} {est:<9} spearman({a}, {b}) = {r:.3f}")
        for f in ("esmc", "proteomelm", "prott5"):
            if ("forest", f) in o and ("lazyqsar", f) in o:
                r = spearmanr(o[("forest", f)], o[("lazyqsar", f)]).statistic
                say(f"  {act:<8} {'across':<9} spearman(forest, lazyqsar) on {f} = {r:.3f}")
    say("")
    say("  High rho with equal AUROC means the second arm is re-deriving signal the first already")
    say("  carries -- not adding an independent view. Read it that way before treating either as")
    say("  second evidence.")


def verdict(res: pd.DataFrame, activators, estimators, feature_sets) -> None:
    """Each arm against forest + esmc -- the FORMER shipped configuration.

    The baseline stays the forest even though TabPFN now ships, because this table's job is to show
    WHY the estimator changed. Re-basing it on TabPFN would make the deltas near zero and destroy
    the evidence for the switch.
    """
    rule()
    say("VERDICT - every arm against the FORMER shipped configuration (forest + esmc)")
    rule()
    for act in activators:
        r = res[res.activator == act]
        base = r[(r.estimator == "forest") & (r.feature_set == "esmc")]
        if base.empty:
            continue
        b = base.iloc[0]
        say(f"  {act}   baseline forest+esmc = AUROC {b.roc_auc:.4f}+/-{b.roc_auc_sd:.4f}  "
            f"PR {b.pr_auc:.4f}+/-{b.pr_auc_sd:.4f}")
        for est in estimators:
            for fs in feature_sets:
                row = r[(r.estimator == est) & (r.feature_set == fs)]
                if row.empty or (est == "forest" and fs == "esmc"):
                    continue
                a = row.iloc[0]

                # Judge BOTH metrics. AUROC averages over the whole curve, including a tail no
                # shortlist will ever reach; PR-AUC at a 13.7%/24.6% base rate weights the top,
                # which is what a prioritized list actually consumes. They can disagree -- measured:
                # lazy-qsar matched the forest on ADEP4 AUROC (+0.0018) while beating it on PR-AUC
                # (+0.0207) -- so reporting only AUROC would have called a real gain "no difference".
                def judge(val, sd, bval, bsd):
                    d = val - bval
                    noise = float(np.hypot(bsd, sd))     # not vs 0: vs the reseeding noise
                    return d, ("BETTER" if d > 2 * noise else
                               "WORSE" if d < -2 * noise else "no difference")

                # PAIRED: same folds per seed, so test the per-seed differences. This is the
                # powerful test and the one the design earns; the unpaired fallback below is only
                # for rows written before per-seed values were stored.
                def paired(col):
                    try:
                        x = np.array([float(v) for v in a[f"{col}_seeds"].split(";")])
                        z = np.array([float(v) for v in b[f"{col}_seeds"].split(";")])
                    except Exception:
                        return None
                    if len(x) != len(z) or len(x) < 2:
                        return None
                    dif = x - z
                    # SD of the MEAN difference. A consistent shift survives here even when both
                    # arms are individually noisy, because the shared difficulty cancels.
                    se = float(np.std(dif, ddof=1) / np.sqrt(len(dif)))
                    return float(dif.mean()), se, int((dif > 0).sum()), len(dif)

                out = []
                for col, lab_ in (("roc_auc", "AUROC"), ("pr_auc", "PR")):
                    pr_ = paired(col)
                    if pr_ is not None:
                        d, se, wins, n_ = pr_
                        call = ("BETTER" if d > 2 * se else
                                "WORSE" if d < -2 * se else "no difference")
                        out.append(f"{lab_} {a[col]:.4f} {d:>+7.4f}+/-{se:.4f} "
                                   f"[{wins}/{n_} seeds] {call:<14}")
                    else:
                        d, call = judge(a[col], a[f"{col}_sd"], b[col], b[f"{col}_sd"])
                        out.append(f"{lab_} {a[col]:.4f} {d:>+7.4f} {call:<14} (unpaired)")
                say(f"    {est:<9} {fs:<11} " + "  ".join(out))
        say("")
    say("  The delta is the mean PER-SEED difference on identical folds, +/- its standard error;")
    say("  'BETTER' means it exceeds 2 SE. This is a PAIRED test -- shared partition difficulty")
    say("  cancels, which matters because on ONC212 both arms carry SD ~0.008 from the splits")
    say("  alone. `[k/n seeds]` is how many seeds the arm won outright: read it alongside the")
    say("  delta, since 5/5 with a small mean is a consistent shift, not a fluke.")


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--activator", nargs="+", default=list(D.ACTIVATORS),
                    choices=list(D.ACTIVATORS))
    ap.add_argument("--estimator", nargs="+", default=list(DEFAULT_ESTIMATORS),
                    choices=list(ESTIMATORS))
    ap.add_argument("--features", nargs="+", default=list(DEFAULT_FEATURE_SETS),
                    choices=list(FEATURE_SETS))
    ap.add_argument("--folds", type=int, default=D.N_SPLITS)
    ap.add_argument("--seed", type=int, default=D.SEED)
    ap.add_argument("--refresh", action="store_true", help="refit lazy-qsar, ignoring its cache")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    # Keep the printed and written order canonical whatever order the flags came in.
    estimators = [e for e in ESTIMATORS if e in args.estimator]
    feature_sets = [f for f in FEATURE_SETS if f in args.features]

    binary, lq_version, tp_version = {}, "-", "-"
    if "lazyqsar" in estimators:
        binary["lazyqsar"], lq_version = lazyqsar_bin()
    if "tabpfn" in estimators:
        binary["tabpfn"], tp_version = tabpfn_bin()

    rule("=")
    say("STAGE 04 sibling - ProteomeLM vs ESM-C, lazy-qsar vs the stage-04 forest")
    rule("=")
    say(f"  species      {SPECIES}   (the only one with measured labels)")
    say(f"  feature sets {', '.join(feature_sets)}")
    say(f"  estimators   {', '.join(estimators)}")
    say(f"  cv           cluster-grouped, {args.folds} folds, {D.N_CV_SEEDS} seeds, PAIRED on "
        "materialised splits")
    say(f"  forest       RandomForest{D.RF_PARAMS}, imported from stage 04")
    say(f"  lazy-qsar    {lq_version} in gradi-lazyqsar, dispatched across a process boundary")
    if "tabpfn" in estimators:
        say(f"  tabpfn       {tp_version} (TabPFN-3.5) in gradi-tabpfn -- NON-COMMERCIAL weights")
    say("  stage 04 NOW ships tabpfn + esmc: ADEP4 0.8738 +/- 0.0029 | ONC212 0.7671 +/- 0.0060")
    say("  the forest it replaced:            ADEP4 0.8645 +/- 0.0021 | ONC212 0.7517 +/- 0.0075")
    say("  honest length baselines there (LOGISTIC): ADEP4 0.776 | ONC212 0.687 -- the `length`")
    say("  rows below use the same estimator as everything else, which reads ~0.09 low for a forest")
    rule("=")

    s04 = stage04()
    rule()
    say("FEATURES")
    rule()
    feats = build_features()
    lab = D.load_labels()

    rule()
    say("CROSS-VALIDATE")
    rule()
    rows, oof = [], {}
    for act in args.activator:
        r, o = evaluate(s04, binary, act, lab, feats, estimators, feature_sets,
                        args.folds, args.seed, args.refresh)
        rows += r
        oof[act] = o
    res = pd.DataFrame(rows)

    redundancy(oof, args.activator, estimators)
    verdict(res, args.activator, estimators, feature_sets)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "proteomelm_vs_esmc.tsv"
    # MERGE rather than overwrite. A partial run (`--features length`) must not silently replace a
    # full grid's table with one row -- which is exactly what the first version of this did.
    KEY = ["activator", "estimator", "feature_set"]
    if out.exists():
        old = pd.read_csv(out, sep="\t")
        fresh = set(map(tuple, res[KEY].values))
        old = old[~old[KEY].apply(tuple, axis=1).isin(fresh)]
        n_kept = len(old)
        res = pd.concat([old, res], ignore_index=True)
    else:
        n_kept = 0
    res = res.sort_values(KEY).reset_index(drop=True)
    res.to_csv(out, sep="\t", index=False)
    rule()
    say(f"  wrote {out.relative_to(REPO_ROOT)}  ({len(res)} rows"
        + (f"; {n_kept} carried over from the previous run)" if n_kept else ")"))
    rule("=")
    say("comparison complete.")
    rule("=")


if __name__ == "__main__":
    main()
