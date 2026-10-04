"""Degradability -- is this protein a substrate of activated partnerless ClpP?

Stage 04. A binary classifier trained on the two *S. aureus* activated-ClpP screens using the
stage-01 ESM-C embeddings, then applied to the rest of *S. aureus* and the whole of *E. coli* and
*K. pneumoniae*. Every protein ends with a probability on one comparable scale, plus the measured
call wherever a screen made one.

Why this axis is the hard one
----------------------------
v1 built it from degron motifs and disorder and pre-registered a gate -- *does any sequence feature
beat protein length alone?* The answer was no. Disorder added **+0.009** over length, a 13-feature
gradient-boosted model scored **0.762** against logistic length-only at **0.775**, and degron motifs
took **0%** of SHAP attribution. Its retrospective's conclusion was blunt:

    A degradability axis needs new *data*, not new features.

Two things make this a different attempt rather than a re-derivation:

1. **The screens are the label now, not a feature.** v1 fed Conlon+Jacques into its composite, so it
   could never validate against them -- `degradability_datasets.md` section 10.1 says exactly that:
   "already the FEATURE in `10c` -> circular". Training on them directly is the non-circular use.
2. **`legacy/scripts/10f_activator_esm_head.py` already specified this stage** -- train on S. aureus
   activator labels with ESM-C, apply to E. coli -- and never ran. Its outputs are absent and its
   input `saureus_esmc600m.npz` never existed; v1 was blocked for want of S. aureus embeddings, which
   stage 01 now has. This stage executes that design and inherits its methodology.

Binary classification, and why the obvious objection does not apply
-------------------------------------------------------------------
The two activators' binary calls overlap at only **Jaccard 0.32**, which looks like a reason not to
classify. It is not. Jaccard measures *set overlap* -- where each experiment drew its line -- and
ONC212 simply drew a looser one, calling 165 hits ADEP4 misses against 28 the other way. A classifier
emits a probability and is scored on *ranking*, and on ranking the two agree well: one activator's
continuous readout predicts the other's binary call at **AUROC 0.815 / 0.877** over 1,022 shared
proteins. That threshold disagreement never touches a probability.

What classification buys is an exact comparison. v1's 0.7747, 0.6799, 0.6916 and 0.7650 are all
logistic-classifier AUROCs, so "beat length" becomes like-for-like rather than approximate. What it
costs is magnitude inside each class: a -9.4 log2FC and a -1.1 both become 1.

**The ceiling is AUROC 0.877 (adep4) / 0.815 (onc212), not 1.0.**

Two activators, never merged
----------------------------
`adep4_*` from Conlon 2013, `onc212_*` from Jacques 2020, separate models and separate columns. v1's
`10c` max-pooled them and its own retrospective flagged that as inflating the score.

The label is asymmetric, on purpose
-----------------------------------
ADEP4 is `log2FC <= -1 AND padj < 0.05`, ONC212 is `log2FC <= -1` alone -- Jacques published no
p-values at all. Each is the best its own paper supports, both are v1's audited definitions
(`legacy/scripts/10e_activator_features.py:119-131`), and that is what keeps these numbers comparable
to v1's. Unmeasured is NaN, never a negative.

The identifier problem, and why the join is by sequence
------------------------------------------------------
Conlon is strain **COL** (`SACOL#####` / `YP_18xxxx`), Jacques is **C0673** (`ODV*`), and the v2
proteome is **NCTC 8325** (`SAOUHSC_*` / `YP_498xxx`). Measured: **0 of 1,943 labeled proteins join
by any identifier** -- not RefSeq, not locus tag, not Jacques's own Mu50 UniProt column. Modern
RefSeq re-tagged COL as `SACOL_RS#####` keeping no `old_locus_tag`, and UniProt demoted the COL and
Mu50 proteomes as redundant, so a UniProt bridge would silently drop ~45% of Conlon.

So we go to sequence, which is CLAUDE.md's standing rule. Within species it is nearly lossless:
**1,873 / 1,943 (96.4%)** at >=95% identity and >=80% coverage, median identity 100%. DIAMOND comes
from `gradi-ortho` (`GRADI_DIAMOND_BIN` overrides).

What the E. coli and K. pneumoniae rows are not
----------------------------------------------
There are **no** activated-ClpP measurements for either species, anywhere, so `<act>_hit` is empty
for every one of their rows and every probability is a **ranking hypothesis** from a model trained in
another phylum. `nn_similarity` records how far each prediction reached in ESM-C space and
`evidence/domain_bands.tsv` prices it -- banding out-of-fold AUROC by that distance and reweighting
to the distances each species actually lands in. Read that file's caveats before quoting it.

Output
------
    data/processed/degradability/degradability_<species>.tsv
        uniprot_ac   adep4_prob   onc212_prob   nn_similarity

`_prob` is the model's probability for EVERY protein -- out-of-fold where labeled, so it is one
comparable scale across all 13,020. **The measured calls are not a column here**: they were empty
for every E. coli and K. pneumoniae row, and they live in `evidence/labels_saureus.tsv` with their
continuous log2FCs. `src.degradability.measured()` re-attaches them on `uniprot_ac`.

    evidence/labels_saureus.tsv          measured call + log2FC on uniprot_ac + join evidence
    evidence/seqmap_audit.tsv            every join: source acc, target, pident, coverage, verdict
    evidence/cv_<activator>.tsv          block x split scheme x {roc_auc, pr_auc, base_rate}
    evidence/cross_activator.tsv         the non-circular test, both directions, both modes
    evidence/domain_bands.tsv            AUROC by similarity band + per-species reweighting
    evidence/oof_<activator>_saureus.tsv out-of-fold probability per labeled protein, per block
    evidence/model_<activator>.npz       the final TRAINING SET + its CV metrics
                                          (TabPFN learns in context: the training set
                                          IS the model, and it outlives any pickle)
    evidence/manifest.tsv

Run with the `gradi` env:
    python scripts/degradability/predict.py                          # both activators, three species
    python scripts/degradability/predict.py --limit 200              # smoke test -> scratch/smoke_*
    python scripts/degradability/predict.py --activator adep4        # one activator only
    python scripts/degradability/predict.py --refresh                # rebuild the DIAMOND map too
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D
from src import tabpfn as T  # noqa: E402
from src import embeddings as E  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "degradability"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
DEFAULT_SPECIES = ("kpneumoniae", "ecoli", "saureus")
DEFAULT_DIAMOND_DIR = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

# Proteins v1 pre-registered a prediction for, and the prediction held (HISTORY.md section 7.5).
# All three should come out `onc212_hit == 0`, which tests the threshold as well as the value.
SPOT_CHECKS = {
    "dnaK": "ONC212 1.4x down -- NOT a hit; survives on cleavage evidence instead",
    "acpP": "ONC212 unchanged (+0.03) -- the canonical ClpP substrate, and not a hit here",
    "gyrB": "essentially no activated-ClpP evidence at all",
}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------- the sequence join

def ensure_diamond() -> str:
    """Locate the DIAMOND binary, or exit with the fix."""
    override = os.environ.get("GRADI_DIAMOND_BIN")
    if shutil.which("diamond") is None:
        for d in ([Path(override)] if override else [DEFAULT_DIAMOND_DIR]):
            if (d / "diamond").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("diamond")
    if exe is None:
        sys.exit(
            "diamond not found.\n"
            f"  looked on PATH and in {override or DEFAULT_DIAMOND_DIR}\n"
            "  it lives in the `gradi-ortho` env (osx-64, no arm64 build); point GRADI_DIAMOND_BIN\n"
            "  at a directory containing it. The labeled proteins come from three DIFFERENT\n"
            "  S. aureus strains and join to the reference proteome by sequence only -- 0 of 1,943\n"
            "  match by any identifier -- so this stage cannot run without it."
        )
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return f"{exe}  ({(out.stdout or out.stderr).strip().splitlines()[0]})"


def map_labels_to_proteome(refresh: bool) -> pd.DataFrame:
    """DIAMOND the labeled COL/C0673 sequences onto the NCTC 8325 reference proteome.

    Returns the audit frame: one row per labeled protein with its best target, identity, coverage and
    a verdict. Cached, because it is deterministic and takes ~20 s.
    """
    audit_path = EVIDENCE_DIR / "seqmap_audit.tsv"
    if audit_path.exists() and not refresh:
        audit = pd.read_csv(audit_path, sep="\t")
        say(f"  [cache] {audit_path.name}  {int((audit['verdict'] == 'direct').sum()):,} direct")
        return audit

    say(f"  diamond  : {ensure_diamond()}")
    target = P.load("saureus")[["uniprot_ac", "sequence"]]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        with (tmp / "target.faa").open("w") as f:
            for ac, seq in zip(target["uniprot_ac"], target["sequence"]):
                f.write(f">{ac}\n{seq}\n")
        subprocess.run(["diamond", "makedb", "--in", str(tmp / "target.faa"),
                        "-d", str(tmp / "target"), "--quiet"], check=True)
        say(f"  target   : {len(target):,} NCTC 8325 proteins")
        say(f"  query    : {D.LABEL_FASTA.relative_to(REPO_ROOT)}")
        hits = tmp / "hits.tsv"
        subprocess.run(
            ["diamond", "blastp", "-q", str(D.LABEL_FASTA), "-d", str(tmp / "target"),
             "-o", str(hits), "--max-target-seqs", "1", "--outfmt", "6",
             "qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore",
             "--quiet", "--threads", str(max(1, (os.cpu_count() or 2) - 1))],
            check=True)
        h = pd.read_csv(hits, sep="\t",
                        names=["saureus_acc", "uniprot_ac", "pident", "alen",
                               "qlen", "slen", "bitscore"])

    h["coverage"] = h["alen"] / h[["qlen", "slen"]].max(axis=1)
    best = h.sort_values("bitscore", ascending=False).drop_duplicates("saureus_acc")
    direct = (best["pident"] >= D.MIN_PIDENT) & (best["coverage"] >= D.MIN_COVERAGE)
    best = best.assign(verdict=np.where(direct, "direct", "below_threshold"))

    labels = pd.read_csv(D.LABEL_TABLE)[["saureus_acc", "sacol", "gene"]]
    audit = labels.merge(best[["saureus_acc", "uniprot_ac", "pident", "coverage", "verdict"]],
                         on="saureus_acc", how="left")
    audit["verdict"] = audit["verdict"].fillna("no_hit")
    audit["pident"] = audit["pident"].round(2)
    audit["coverage"] = audit["coverage"].round(4)
    audit.to_csv(audit_path, sep="\t", index=False)
    return audit


def build_labels(audit: pd.DataFrame) -> pd.DataFrame:
    """Measured calls keyed on `uniprot_ac`, with their log2FCs and the cluster for CV grouping."""
    raw = pd.read_csv(D.LABEL_TABLE)
    clusters = pd.read_csv(D.LABEL_CLUSTERS, sep="\t", names=["cluster", "saureus_acc"])

    wanted = list(D.ACTIVATORS.values()) + list(D.ACTIVATOR_CONTINUOUS.values())
    keep = audit[audit["verdict"] == "direct"][["saureus_acc", "uniprot_ac", "pident", "coverage"]]
    lab = (keep.merge(raw[["saureus_acc", "sacol", "gene", "length"] + wanted],
                      on="saureus_acc", how="left")
               .merge(clusters, on="saureus_acc", how="left"))
    # `adep4`/`onc212` are the binary labels; `<act>_log2fc` keeps the value behind each call.
    lab = lab.rename(columns={v: k for k, v in D.ACTIVATORS.items()})
    lab = lab.rename(columns={v: f"{k}_log2fc" for k, v in D.ACTIVATOR_CONTINUOUS.items()})
    # Two labeled proteins landing on one accession would double-weight it; keep the better join.
    lab = lab.sort_values(["pident", "coverage"], ascending=False).drop_duplicates("uniprot_ac")
    # A protein with no cluster assignment is its own group, never a shared one.
    lab["cluster"] = lab["cluster"].fillna(lab["saureus_acc"])
    lab.to_csv(EVIDENCE_DIR / "labels_saureus.tsv", sep="\t", index=False)
    return lab


# ---------------- features

def feature_frame(species: str) -> pd.DataFrame:
    """Sequence features plus ESM-C embeddings for one species, indexed by `uniprot_ac`."""
    prot = P.load(species)[["uniprot_ac", "sequence"]].set_index("uniprot_ac")
    feats = D.sequence_features(prot["sequence"])
    accs, mat = E.load(species)
    emb = pd.DataFrame(mat, index=pd.Index(accs, name="uniprot_ac"),
                       columns=[f"e{i}" for i in range(mat.shape[1])])
    out = feats.join(emb, how="inner")
    if len(out) != len(prot):
        say(f"    [!] {species}: {len(prot) - len(out)} proteins have no embedding and are dropped")
    return out


# ---------------- cross-validation

def _pipe(seed: int = D.SEED):
    """One estimator for every feature block, so a difference is about features not capacity.

    `RandomForestClassifier` with `class_weight="balanced"`, which carries the 13.7% / 24.6% positive
    rates without resampling. No `StandardScaler`: trees split on thresholds and are scale-invariant,
    so scaling would only cost time.

    Chosen over the balanced logistic v1's `10f` specified, on measured grounds. Under identical
    cluster-grouped CV on the ESM-C block the two are a wash -- logistic 0.873 vs RF 0.857 on ADEP4,
    logistic 0.744 vs RF 0.750 on ONC212, with overlapping 95% intervals either way. What breaks the
    tie is robustness: logistic needs its regularisation to be right, and getting it wrong is not
    loud. v1's `10f` hard-coded `C=1.0`, which on 1,152 standardised dimensions at n~1,677 saturated
    **52% of out-of-fold probabilities at exactly 0.000 or 1.000** -- unrankable ties across half the
    proteome -- while also costing accuracy (0.813 against 0.873 at C=0.001). 10f never ran, so that
    was never caught. A forest has no such cliff and saturates 0% of predictions untuned.

    (`HistGradientBoostingClassifier` was also measured and is not used: no better than either
    linear or forest, and it saturated 25-67% of its probabilities.)
    """
    from sklearn.ensemble import RandomForestClassifier

    return RandomForestClassifier(n_jobs=-1, random_state=seed, **D.RF_PARAMS)


# ---------------------------------------------------------------- the estimator seam
#
# Every estimator use in this stage goes through `predict_fold`. Before TabPFN there were five
# separate call sites, each constructing its own sklearn object; TabPFN cannot live in `gradi`
# (torch 2.14 + mlx against gradi's 2.12 under ESM-C), so each had to become an out-of-process call.
# Routing them through one function is what keeps `--estimator forest` bit-identical to the version
# that produced the shipped numbers -- the seam is inert for the forest by construction.

# TabPFN is TRANSVERSAL -- the dispatch machinery, the content-addressed cache, the retry and the
# credit guard live in `src/tabpfn.py` so essentiality, ligands or any other axis can call the
# project's default learner without importlib-loading a degradability script. This stage keeps only
# the seam: one function, so `--estimator forest` stays bit-identical to the shipped numbers.
TABPFN_CACHE = T.CACHE_DIR


def tabpfn_bin() -> str:
    """Kept as a name this stage's callers already use; the check itself is in `src/tabpfn.py`."""
    return T.binary()


def predict_fold(X_train, y_train, X_test, *, estimator: str, seed: int = D.SEED,
                 tag: str = "") -> np.ndarray:
    """Fit on (X_train, y_train), return P(class 1) for X_test. The single estimator call site.

    `forest` is stage 04's shipped RandomForest, unchanged and in-process. `tabpfn` dispatches to
    `src/tabpfn.py`, which serves from a cache keyed on the DATA -- hosted TabPFN bills ~10,000
    credits per call flat, so a re-run that recomputes nothing must also spend nothing.
    """
    if estimator == "forest":
        return _pipe(seed).fit(X_train, y_train).predict_proba(X_test)[:, 1]
    return T.predict_fold(X_train, y_train, X_test, seed=seed, tag=tag)


def predict_fold_regression(X_train, y_train, X_test, *, seed: int = D.SEED,
                            tag: str = "") -> np.ndarray:
    """TabPFN REGRESSION on a continuous target. Same cache, different task in the key."""
    return T.predict_fold_regression(X_train, y_train, X_test, seed=seed, tag=tag)


def oof_predict(X, y, splits, *, estimator: str, seed: int = D.SEED, tag: str = "") -> np.ndarray:
    """Out-of-fold probabilities over a materialised split list.

    Replaces `cross_val_predict`, which can only drive an in-process sklearn estimator. For the
    forest the two are equivalent by construction -- fit on train, predict test -- which is what
    makes `--estimator forest` reproduce the shipped numbers exactly.
    """
    if estimator != "forest":
        return T.oof_predict(X, y, splits, seed=seed, tag=tag)
    prob = np.full(len(y), np.nan)
    for k, (tr, te) in enumerate(splits):
        prob[te] = predict_fold(X[tr], y[tr], X[te], estimator=estimator, seed=seed,
                                tag=f"{tag}/fold{k}")
    if np.isnan(prob).any():
        sys.exit(f"FATAL {tag}: {int(np.isnan(prob).sum())} proteins were in no test fold")
    return prob


def _assert_folds_disjoint(cv, X, y, groups, activator: str, seed: int) -> None:
    """No homology cluster may appear in two folds. The guarantee behind the out-of-fold score.

    Checked per seed rather than once: the seed-averaged probability is only leakage-free if EVERY
    split respected the groups, and `StratifiedGroupKFold` balances classes as well as groups, so it
    is worth confirming it did not quietly trade one for the other.
    """
    seen: dict = {}
    for k, (_, te) in enumerate(cv.split(X, y, groups=groups)):
        for cl in set(groups[te]):
            if cl in seen:
                sys.exit(f"FAILED: {activator} seed {seed}: cluster {cl!r} appears in folds "
                         f"{seen[cl]} and {k}. The out-of-fold score would be leaky.")
            seen[cl] = k


def _saturated(prob) -> float:
    """Fraction of probabilities pinned at 0 or 1 -- i.e. unrankable ties. Watch this."""
    prob = np.asarray(prob, dtype=float)
    return float(((prob <= 0.001) | (prob >= 0.999)).mean())


def _splitters(folds: int, seed: int):
    """Stratified-grouped and stratified-plain, in that order.

    `StratifiedGroupKFold` rather than plain `GroupKFold`: at a 13.6% positive rate an unstratified
    fold can end up with too few positives to score an AUROC at all. It respects the homology groups
    AND the class balance. The plain scheme runs alongside purely so the leakage gap is visible --
    v1's `07d` used an ungrouped `cv=5` for its essentiality head, and this is what that costs.
    """
    from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold

    return (("clustered", StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)),
            ("plain", StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)))


def cross_validate(activator: str, lab: pd.DataFrame, feats: pd.DataFrame,
                   folds: int, seed: int,
                   estimator: str = "forest") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Repeated cluster-grouped CV. Returns (cv table, seed-averaged out-of-fold probabilities)."""
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import cross_val_predict

    sub = lab[lab[activator].notna()].copy()
    sub = sub[sub["uniprot_ac"].isin(feats.index)]
    X_all = feats.loc[sub["uniprot_ac"]]
    y = sub[activator].to_numpy(dtype=float).astype(int)
    groups = sub["cluster"].to_numpy()
    base_rate = float(y.mean())
    say(f"  {activator:<8} {len(sub):>5} labeled, {pd.Series(groups).nunique():>5} clusters, "
        f"{int(y.sum()):>4} positives ({100 * base_rate:.1f}%)")

    rows, oof_store, fold_of = [], {}, None
    for block in D.FEATURE_BLOCKS:
        cols = D.block_columns(block, X_all)
        X = X_all[cols].to_numpy(dtype=float)
        rec = {"activator": activator, "block": block, "estimator": estimator,
               "n_features": len(cols), "n": len(sub), "n_pos": int(y.sum()),
               "base_rate": round(base_rate, 4), "n_cv_seeds": D.N_CV_SEEDS}
        for scheme in ("clustered", "plain"):
            # The CV is REPEATED over seeds and averaged. Not decoration: a single-seed estimate
            # carries ~+/-0.004 of arbitrariness here -- larger than any effect this stage has tested
            # and rejected -- so quoting one seed makes runs look different when nothing changed.
            aucs, prs, probs = [], [], []
            for si in range(D.N_CV_SEEDS):
                cv = dict(_splitters(folds, seed + si))[scheme]
                # StratifiedKFold ignores `groups` and warns about it; pass it only where used.
                g = groups if scheme == "clustered" else None
                # Materialise the splits, then drive them ourselves: `cross_val_predict` can only
                # run an in-process sklearn estimator, and TabPFN lives in another interpreter.
                # For the forest this is identical work in the same order.
                splits = list(cv.split(X, y, groups=g))
                prob = oof_predict(X, y, splits, estimator=estimator, seed=seed,
                                   tag=f"{activator}/{block}/{scheme}/s{seed + si}")
                aucs.append(roc_auc_score(y, prob))
                prs.append(D.pr_auc(y, prob))
                probs.append(prob)
                if scheme == "clustered":
                    _assert_folds_disjoint(cv, X, y, groups, activator, seed + si)
                    if si == 0:
                        # Seed 0's folds, so the per-fold ROC curves have one split to draw.
                        fold_of = np.empty(len(y), dtype=int)
                        for k, (_, te) in enumerate(cv.split(X, y, groups=groups)):
                            fold_of[te] = k
            # Seed-averaged out-of-fold probability. Every seed's value for a protein comes from a
            # model that excluded it, so averaging is leakage-free and strictly lower-variance.
            mean_prob = np.mean(probs, axis=0)
            st = D.auroc_ci(y, mean_prob)
            rec[f"roc_auc_{scheme}"] = round(float(np.mean(aucs)), 4)
            rec[f"roc_auc_{scheme}_sd"] = round(float(np.std(aucs)), 4)
            rec[f"roc_auc_{scheme}_lo"] = st["roc_auc_lo"]
            rec[f"roc_auc_{scheme}_hi"] = st["roc_auc_hi"]
            rec[f"pr_auc_{scheme}"] = round(float(np.mean(prs)), 4)
            rec[f"pr_auc_{scheme}_sd"] = round(float(np.std(prs)), 4)
            rec[f"saturated_{scheme}"] = round(_saturated(mean_prob), 4)
            if scheme == "clustered":
                oof_store[block] = mean_prob
        rec["leakage_gap"] = round(rec["roc_auc_plain"] - rec["roc_auc_clustered"], 4)
        if rec["roc_auc_clustered_sd"] == 0.0:
            sys.exit(f"FAILED: {activator}/{block} reports a zero seed SD over "
                     f"{D.N_CV_SEEDS} seeds -- the seed loop collapsed to one split.")
        rows.append(rec)

    cv = pd.DataFrame(rows)
    oof = pd.DataFrame({"uniprot_ac": sub["uniprot_ac"].to_numpy(),
                        "sacol": sub["sacol"].to_numpy(),
                        "gene": sub["gene"].to_numpy(),
                        "cluster": groups,
                        "fold": fold_of,
                        "hit": y,
                        "log2fc": sub[f"{activator}_log2fc"].round(4).to_numpy()})
    for block, prob in oof_store.items():
        oof[f"oof_{block.replace('+', '_')}"] = np.round(prob, 4)
    return cv, oof


def cross_activator(lab: pd.DataFrame, feats: pd.DataFrame,
                    estimator: str = "forest") -> pd.DataFrame:
    """Fit on one activator, score the other's labeled set. Both directions, two modes.

    The only non-circular test available: two chemically distinct activators of the same enzyme.
    `degradability_datasets.md` section 10.4 lists it as "available now, not done".
    """
    cols = D.block_columns("esmc", feats)
    rows = []
    for train, test in (("adep4", "onc212"), ("onc212", "adep4")):
        tr = lab[lab[train].notna() & lab["uniprot_ac"].isin(feats.index)]
        Xtr = feats.loc[tr["uniprot_ac"], cols].to_numpy(dtype=float)
        ytr = tr[train].to_numpy(dtype=float).astype(int)
        base = lab[lab[test].notna() & lab["uniprot_ac"].isin(feats.index)]

        # Two modes, because they answer different questions and only one is comparable to v1.
        #
        #   `shared`          v1's setup: score every protein the other activator labeled, including
        #                     ones this model trained on. Sequence-circular, but it isolates the
        #                     question v1 asked -- does susceptibility transfer between two
        #                     chemistries? -- and is the only apples-to-apples comparison available.
        #   `cluster_disjoint` also drops any test protein sharing a homology cluster with training,
        #                     so it tests chemistry AND sequence novelty at once. Stricter, but the
        #                     adep4 -> onc212 direction runs out of data: adep4 covers almost every
        #                     cluster, so almost nothing of onc212 is left.
        for mode in ("shared", "cluster_disjoint"):
            te = base if mode == "shared" else base[~base["cluster"].isin(set(tr["cluster"]))]
            if len(te) < 30:
                say(f"    {train} -> {test} [{mode}]: only {len(te)} proteins left, skipped")
                continue
            prob = predict_fold(Xtr, ytr,
                                feats.loc[te["uniprot_ac"], cols].to_numpy(dtype=float),
                                estimator=estimator, tag=f"cross/{train}->{test}/{mode}")
            y = te[test].to_numpy(dtype=float).astype(int)
            st = D.auroc_ci(y, prob)
            rows.append({"train_on": train, "evaluate_on": test, "mode": mode,
                         "n_train": len(tr), "n_test": len(te), **st,
                         "pr_auc": D.pr_auc(y, prob), "base_rate": round(float(y.mean()), 4),
                         "v1_auroc_handbuilt": (D.V1_CROSS_ACTIVATOR_AUROC.get(f"{train}_to_{test}")
                                                if mode == "shared" else None)})
    return pd.DataFrame(rows)


def cutoff_sensitivity(lab: pd.DataFrame, feats: pd.DataFrame, activators: list[str],
                       folds: int, seed: int, estimator: str = "forest") -> pd.DataFrame:
    """Is the conclusion an artifact of the -1.0 log2FC cutoff? Sweep it and see.

    **`model` must not be read as "a stricter cutoff is better".** AUROC is not comparable across
    different labels: tightening the cutoff shrinks the positive class and makes it more extreme,
    which inflates AUROC mechanically (ADEP4 runs 0.827 at -0.5 up to 0.943 at -3.0 on 46 positives).
    Optimising it would buy a spurious +0.08 and a label chosen for being easy rather than meaningful.

    The column that matters is `gap`. `cross_assay` -- the other screen's own readout as a predictor
    of these calls -- inflates in lockstep, so their difference is the invariant, and it holds at
    every cutoff. That is what makes the shipped result cutoff-independent.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    cols = D.block_columns("esmc", feats)
    ll = np.log10(lab.set_index("uniprot_ac")["length"].astype(float))
    rows = []
    for a in activators:
        other = "onc212" if a == "adep4" else "adep4"
        for cut in D.CUTOFF_SWEEP:
            sub = lab[lab[f"{a}_log2fc"].notna() & lab["uniprot_ac"].isin(feats.index)]
            y = (sub[f"{a}_log2fc"] <= cut).to_numpy().astype(int)
            if y.sum() < 25:
                say(f"    {a:<7} {cut:>5.1f}: only {y.sum()} positives, skipped")
                continue
            g = sub["cluster"].to_numpy()
            cv = dict(_splitters(folds, seed))["clustered"]
            X = feats.loc[sub["uniprot_ac"], cols].to_numpy(dtype=float)
            prob = oof_predict(X, y, list(cv.split(X, y, groups=g)), estimator=estimator,
                               seed=seed, tag=f"cutoff/{a}/{cut}")
            model = float(roc_auc_score(y, prob))
            # The other screen's measured readout, ranking THESE calls. -log2FC because depletion is
            # negative, so more-negative must rank as more-likely-a-hit.
            m2 = sub[sub[f"{other}_log2fc"].notna()]
            y2 = (m2[f"{a}_log2fc"] <= cut).to_numpy().astype(int)
            xa = (float(roc_auc_score(y2, -m2[f"{other}_log2fc"]))
                  if 10 <= y2.sum() < len(y2) else float("nan"))
            # Length-only under BALANCED LOGISTIC, not the forest: a forest on one feature bins it
            # into steps and scores ~0.09 lower, which would flatter the embedding for free.
            lp = make_pipeline(StandardScaler(),
                               LogisticRegression(max_iter=3000, class_weight="balanced"))
            pl = cross_val_predict(lp, ll.reindex(sub["uniprot_ac"]).to_numpy().reshape(-1, 1),
                                   y, cv=cv, groups=g, method="predict_proba")[:, 1]
            rows.append({"activator": a, "cutoff": cut, "n": len(y), "n_pos": int(y.sum()),
                         "base_rate": round(float(y.mean()), 4), "model": round(model, 4),
                         "cross_assay": round(xa, 4) if xa == xa else None,
                         "gap": round(model - xa, 4) if xa == xa else None,
                         "length_logistic": round(float(roc_auc_score(y, pl)), 4),
                         "esmc_over_length": round(model - float(roc_auc_score(y, pl)), 4)})
    return pd.DataFrame(rows)


def domain_bands(activator: str, oof: pd.DataFrame, lab: pd.DataFrame,
                 feats: pd.DataFrame, species_list: list[str]) -> pd.DataFrame:
    """Price the extrapolation: out-of-fold AUROC by distance to the training set, reweighted.

    Quoting the within-S. aureus AUROC for an E. coli prediction would overstate it, because the
    labeled set is dominated by proteins whose nearest training neighbour is close, while a real
    query is further away. So band by cosine and reweight by where each species' queries land.
    """
    cols = D.block_columns("esmc", feats)
    train_acc = [a for a in lab.loc[lab[activator].notna(), "uniprot_ac"] if a in feats.index]
    T = feats.loc[train_acc, cols].to_numpy(dtype=float)
    T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-9

    def nn_sim(species: str, exclude_self: bool) -> pd.Series:
        f = feats if species == "saureus" else feature_frame(species)
        Q = f[cols].to_numpy(dtype=float)
        Q = Q / (np.linalg.norm(Q, axis=1, keepdims=True) + 1e-9)
        sim = Q @ T.T
        if exclude_self:
            pos = {a: i for i, a in enumerate(train_acc)}
            for r, a in enumerate(f.index):
                if a in pos:
                    sim[r, pos[a]] = -2.0        # never let a protein be its own neighbour
        return pd.Series(sim.max(axis=1), index=f.index)

    sa_sim = nn_sim("saureus", exclude_self=True)
    o = oof.set_index("uniprot_ac")
    o = o.assign(nn=sa_sim.reindex(o.index))
    per_band, rows = {}, []
    for lo, hi in D.SIMILARITY_BANDS:
        g = o[(o["nn"] >= lo) & (o["nn"] < hi)]
        st = (D.auroc_ci(g["hit"], g["oof_esmc"]) if len(g) >= 20
              else {"roc_auc": np.nan, "roc_auc_lo": np.nan, "roc_auc_hi": np.nan,
                    "n": len(g), "n_pos": int(g["hit"].sum()) if len(g) else 0})
        per_band[(lo, hi)] = st["roc_auc"]
        rows.append({"activator": activator, "scope": "saureus_labeled_oof",
                     "band": f"[{lo:.2f},{hi:.2f})", "n": st["n"], "n_pos": st["n_pos"],
                     "roc_auc": st["roc_auc"], "roc_auc_lo": st["roc_auc_lo"],
                     "roc_auc_hi": st["roc_auc_hi"], "share": np.nan, "expected_roc_auc": np.nan})

    for sp in species_list:
        sim = nn_sim(sp, exclude_self=(sp == "saureus"))
        expected, total = 0.0, 0.0
        for (lo, hi), auc in per_band.items():
            share = float(((sim >= lo) & (sim < hi)).mean())
            rows.append({"activator": activator, "scope": sp, "band": f"[{lo:.2f},{hi:.2f})",
                         "n": int(((sim >= lo) & (sim < hi)).sum()), "n_pos": np.nan,
                         "roc_auc": auc, "roc_auc_lo": np.nan, "roc_auc_hi": np.nan,
                         "share": round(share, 4), "expected_roc_auc": np.nan})
            if np.isfinite(auc):
                expected += share * auc
                total += share
        rows.append({"activator": activator, "scope": sp, "band": "REWEIGHTED", "n": len(sim),
                     "n_pos": np.nan, "roc_auc": np.nan, "roc_auc_lo": np.nan,
                     "roc_auc_hi": np.nan, "share": round(total, 4),
                     "expected_roc_auc": round(expected / total, 4) if total else np.nan})
    return pd.DataFrame(rows)


# ---------------- predict

def fit_final(activator: str, lab: pd.DataFrame, feats: pd.DataFrame, cv: pd.DataFrame,
              estimator: str = "forest"):
    """Bank the final training set, after the honest CV number is already banked.

    **The persisted artifact changed when TabPFN arrived, and it had to.** The forest pickled a
    fitted `sklearn.Pipeline`; TabPFN has no such object -- it learns IN CONTEXT, so its "model" is
    exactly the training set plus the configuration, and hosted inference keeps no local state at
    all. Persisting the training set is therefore the honest artifact for both: it is what actually
    determines the predictions, it survives a library upgrade that would invalidate a pickle, and
    for the forest it reproduces the same model given the same seed.

    Returns `(X_train, y_train)` -- fed to `predict_fold` per species by the caller.
    """
    sub = lab[lab[activator].notna() & lab["uniprot_ac"].isin(feats.index)]
    cols = D.block_columns("esmc", feats)
    Xtr = feats.loc[sub["uniprot_ac"], cols].to_numpy(dtype=float)
    ytr = sub[activator].to_numpy(dtype=float).astype(int)
    row = cv[cv["block"] == "esmc"].iloc[0]
    np.savez_compressed(
        EVIDENCE_DIR / f"model_{activator}.npz",
        accessions=sub["uniprot_ac"].to_numpy().astype(object),
        X_train=Xtr.astype("float32"), y_train=ytr,
        estimator=estimator, activator=activator, block="esmc",
        n_train=len(sub), n_pos=int(sub[activator].sum()),
        roc_auc_clustered=float(row["roc_auc_clustered"]),
        roc_auc_plain=float(row["roc_auc_plain"]),
        pr_auc_clustered=float(row["pr_auc_clustered"]),
        label=D.ACTIVATOR_PAPER[activator])
    return Xtr, ytr


# ---------------- main

def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES),
                    choices=list(DEFAULT_SPECIES),
                    help="species to score (human has no labels and is out of scope)")
    ap.add_argument("--activator", nargs="+", default=list(D.ACTIVATORS),
                    choices=list(D.ACTIVATORS))
    ap.add_argument("--folds", type=int, default=D.N_SPLITS)
    ap.add_argument("--seed", type=int, default=D.SEED)
    ap.add_argument("--estimator", default=D.DEFAULT_ESTIMATOR, choices=("forest", "tabpfn"),
                    help="tabpfn (default, measured better) or forest (the former shipped model)")
    ap.add_argument("--limit", type=int, help="only the first N labeled proteins (smoke test)")
    ap.add_argument("--refresh", action="store_true", help="rebuild the DIAMOND sequence map")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    pre = "smoke_" if args.limit else ""

    rule("=")
    say("STAGE 04 - degradability from the S. aureus activator screens")
    rule("=")
    say(f"  labels   : {D.LABEL_TABLE.relative_to(REPO_ROOT)}  (v1's audited *_bin columns)")
    say("  features : stage-01 ESM-C 600M (1152-dim) + length/composition baselines")
    say(f"  estimator: {args.estimator}"
        + (f"  RandomForest{tuple(D.RF_PARAMS.values())}" if args.estimator == "forest"
           else f"  TabPFN-{D.TABPFN_MODEL} via gradi-tabpfn (non-commercial weights)"))
    say(f"  activators: {', '.join(args.activator)}  (separate models, never merged)")
    say(f"  species  : {', '.join(args.species)}")
    say(f"  cv       : StratifiedGroupKFold({args.folds}) on MMseqs2 30%-id clusters "
        f"+ StratifiedKFold, seed {args.seed}")
    say(f"  cross-assay ref: AUROC {D.CROSS_ASSAY_AUROC}   (what the OTHER screen achieves; not a bound)")
    if args.limit:
        say(f"  limit    : {args.limit} labeled proteins (SMOKE TEST)")
    if args.dry_run:
        say("\n  --dry-run: nothing fitted, nothing written.")
        for sp in args.species:
            say(f"    {sp:<14} -> degradability_{sp}.tsv")
        return
    say()

    rule()
    say("MAP - labeled proteins onto the reference proteome, by sequence")
    rule()
    audit = map_labels_to_proteome(args.refresh)
    n_direct = int((audit["verdict"] == "direct").sum())
    frac = n_direct / len(audit)
    say(f"  {len(audit):,} labeled proteins -> {n_direct:,} direct "
        f"({100 * frac:.1f}%) at >={D.MIN_PIDENT}% id / >={100 * D.MIN_COVERAGE:.0f}% cov")
    for v, n in audit["verdict"].value_counts().items():
        say(f"    {v:<16} {n:>5}")
    ok = audit[audit["verdict"] == "direct"]
    say(f"  identity: median {ok['pident'].median():.1f}%  5th pct {ok['pident'].quantile(.05):.1f}%")
    if frac < D.MIN_MAP_FRACTION:
        sys.exit(f"FAILED: only {100 * frac:.1f}% of labeled proteins mapped, floor is "
                 f"{100 * D.MIN_MAP_FRACTION:.0f}%. Every number downstream depends on this join.")
    say()

    rule()
    say("LABELS")
    rule()
    lab = build_labels(audit)
    if args.limit:
        lab = lab.head(args.limit)
    say(f"  {len(lab):,} proteins on uniprot_ac, {lab['cluster'].nunique():,} homology clusters")
    from sklearn.metrics import roc_auc_score
    prot = P.load("saureus")[["uniprot_ac", "sequence"]].set_index("uniprot_ac")
    ll = D.sequence_features(prot["sequence"])["log_length"].reindex(lab["uniprot_ac"]).to_numpy()
    bad = []
    for a in D.ACTIVATORS:
        s = lab[a]
        say(f"    {a:<8} n={s.notna().sum():>5}  positives {int(s.sum()):>4} "
            f"({100 * s.mean():.1f}%)  [{D.ACTIVATOR_PAPER[a]}]")
    for a in args.activator:
        # Direction guard. Short proteins deplete under activated ClpP -- v1 measured raw `length` at
        # AUROC 0.2247, i.e. inverted -- so ranking SHORTER-first must retrieve hits. A flip here
        # would invert the whole deliverable while every other number still looked plausible.
        y = lab[a].to_numpy(dtype=float)
        keep = np.isfinite(y) & np.isfinite(ll)
        auc = float(roc_auc_score(y[keep].astype(int), -ll[keep]))
        flag = "ok" if auc > 0.5 else "<- FAIL (inverted: label or join is wrong)"
        say(f"    direction  AUROC(hit, shorter-first) for {a} = {auc:.3f}   {flag}")
        if auc <= 0.5:
            bad.append(f"{a}: AUROC(hit, -log_length) = {auc:.3f}, must exceed 0.5")
    if bad:
        sys.exit("FAILED direction check:\n  " + "\n  ".join(bad))
    say()

    feats = feature_frame("saureus")
    rule()
    say("CROSS-VALIDATE - the baseline ladder, cluster-grouped vs plain")
    rule()
    cvs, oofs = {}, {}
    for a in args.activator:
        cv, oof = cross_validate(a, lab, feats, args.folds, args.seed, args.estimator)
        cv.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}cv_{a}.tsv", sep="\t", index=False)
        oof.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}oof_{a}_saureus.tsv", sep="\t", index=False)
        cvs[a], oofs[a] = cv, oof
        say(f"    {'block':<12} {'nfeat':>6} {'ROC-AUC':>16} {'95% CI':>16} "
            f"{'PR-AUC':>15} {'plain':>7} {'leak':>7} {'sat%':>6}")
        for _, r in cv.iterrows():
            say(f"    {r['block']:<12} {r['n_features']:>6} "
                f"{r['roc_auc_clustered']:>8.3f} +/-{r['roc_auc_clustered_sd']:<5.3f} "
                f"[{r['roc_auc_clustered_lo']:>6.3f},{r['roc_auc_clustered_hi']:>6.3f}] "
                f"{r['pr_auc_clustered']:>7.3f} +/-{r['pr_auc_clustered_sd']:<5.3f} "
                f"{r['roc_auc_plain']:>7.3f} {r['leakage_gap']:>7.3f} "
                f"{100 * r['saturated_clustered']:>6.1f}")
        esm = float(cv.loc[cv["block"] == "esmc", "roc_auc_clustered"].iloc[0])
        sd = float(cv.loc[cv["block"] == "esmc", "roc_auc_clustered_sd"].iloc[0])
        xa = D.CROSS_ASSAY_AUROC[a]
        say(f"    -> this model {esm:.3f} +/- {sd:.3f} over {D.N_CV_SEEDS} CV seeds")
        say(f"       v1 length-only {D.V1_LENGTH_ONLY_AUROC[a]}   "
            f"cross-assay reference {xa}   gap {esm - xa:+.3f}")
        say(f"    -> base rate {float(cv['base_rate'].iloc[0]):.3f}, so a PR-AUC must be read "
            "against that, not against 0.5")
        say()

    rule()
    say("CROSS-ACTIVATOR - the only non-circular test")
    rule()
    xa = cross_activator(lab, feats, args.estimator)
    if not xa.empty:
        xa.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}cross_activator.tsv", sep="\t", index=False)
        say(f"  {'train':<8} {'test':<8} {'mode':<17} {'n':>6} {'ROC-AUC':>8} {'95% CI':>16} "
            f"{'PR-AUC':>7} {'v1':>7}")
        for _, r in xa.iterrows():
            v1 = f"{r['v1_auroc_handbuilt']:.4f}" if pd.notna(r["v1_auroc_handbuilt"]) else "-"
            say(f"  {r['train_on']:<8} {r['evaluate_on']:<8} {r['mode']:<17} {r['n_test']:>6} "
                f"{r['roc_auc']:>8.3f} [{r['roc_auc_lo']:>6.3f},{r['roc_auc_hi']:>6.3f}] "
                f"{r['pr_auc']:>7.3f} {v1:>7}")
        say("\n  `shared` is v1's setup and the only comparable one; `cluster_disjoint` also")
        say("  requires sequence novelty, which the adep4 -> onc212 direction cannot supply.")
    say()

    rule()
    say("CUTOFF SENSITIVITY - is the result an artifact of the -1.0 log2FC threshold?")
    rule()
    sens = cutoff_sensitivity(lab, feats, args.activator, args.folds, args.seed, args.estimator)
    sens.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}cutoff_sensitivity.tsv", sep="\t", index=False)
    say(f"  {'act':<8} {'cutoff':>7} {'n_pos':>6} {'model':>7} {'cross-assay':>12} "
        f"{'gap':>7} {'length':>7} {'ESM-C over len':>15}")
    for _, r in sens.iterrows():
        xa = f"{r['cross_assay']:.3f}" if pd.notna(r["cross_assay"]) else "-"
        gp = f"{r['gap']:+.3f}" if pd.notna(r["gap"]) else "-"
        star = "  <- shipped" if r["cutoff"] == D.ABUNDANCE_LOG2_HIT else ""
        say(f"  {r['activator']:<8} {r['cutoff']:>7.1f} {int(r['n_pos']):>6} "
            f"{r['model']:>7.3f} {xa:>12} {gp:>7} {r['length_logistic']:>7.3f} "
            f"{r['esmc_over_length']:>+15.3f}{star}")
    say("\n  `model` RISES with strictness for a trivial reason -- fewer, more extreme positives are")
    say("  easier to separate -- so it is NOT a tuning signal. `cross_assay` rises with it, and the")
    say("  `gap` between them is the invariant. A stable gap is what makes the result cutoff-free.")
    # The sweep and the main CV must agree about the shipped label, or one of them is wrong.
    for a in args.activator:
        row = sens[(sens["activator"] == a) & (sens["cutoff"] == D.ABUNDANCE_LOG2_HIT)]
        if row.empty:
            continue
        main = float(cvs[a].loc[cvs[a]["block"] == "esmc", "roc_auc_clustered"].iloc[0])
        if abs(float(row["model"].iloc[0]) - main) > 0.03:
            sys.exit(f"FAILED: {a} cutoff sweep gives {float(row['model'].iloc[0]):.3f} at "
                     f"{D.ABUNDANCE_LOG2_HIT} but the main CV gives {main:.3f}. The two paths "
                     "disagree about the label.")
    say()

    rule()
    say("DOMAIN - how far each prediction reaches, and what that is worth")
    rule()
    bands = pd.concat([domain_bands(a, oofs[a], lab, feats, args.species)
                       for a in args.activator], ignore_index=True)
    bands.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}domain_bands.tsv", sep="\t", index=False)
    for a in args.activator:
        b = bands[(bands["activator"] == a) & (bands["scope"] == "saureus_labeled_oof")]
        say(f"  {a}: out-of-fold AUROC by distance to nearest other training protein")
        for _, r in b.iterrows():
            v = f"{r['roc_auc']:.3f}" if np.isfinite(r["roc_auc"]) else "  n/a"
            say(f"    {r['band']:<14} n={int(r['n']):>5}  pos={int(r['n_pos']):>4}  AUROC={v}")
        for sp in args.species:
            rw = bands[(bands["activator"] == a) & (bands["scope"] == sp)
                       & (bands["band"] == "REWEIGHTED")]
            if len(rw):
                say(f"    -> {sp:<14} expected AUROC {float(rw['expected_roc_auc'].iloc[0]):.3f} "
                    f"after reweighting to its own distance distribution")
        say()

    rule()
    say("PREDICT")
    rule()
    models = {a: fit_final(a, lab, feats, cvs[a], args.estimator) for a in args.activator}
    meas_hit = {a: dict(zip(lab["uniprot_ac"], lab[a])) for a in D.ACTIVATORS}
    # Out-of-fold probability for the labeled proteins, so `_prob` is one comparable scale.
    oof_prob = {a: dict(zip(oofs[a]["uniprot_ac"], oofs[a]["oof_esmc"])) for a in args.activator}
    rows = []
    audit_rows = []
    for sp in args.species:
        f = feats if sp == "saureus" else feature_frame(sp)
        cols = D.block_columns("esmc", f)
        out = pd.DataFrame({"uniprot_ac": P.load(sp)["uniprot_ac"]})
        for a in D.ACTIVATORS:
            if a not in models:
                out[f"{a}_hit"] = pd.array([pd.NA] * len(out), dtype="Int64")
                out[f"{a}_prob"], out[f"{a}_source"] = np.nan, ""
                continue
            Xtr, ytr = models[a]
            fresh = pd.Series(predict_fold(Xtr, ytr, f[cols].to_numpy(dtype=float),
                                           estimator=args.estimator, tag=f"predict/{a}/{sp}"),
                              index=f.index)
            hit = out["uniprot_ac"].map(meas_hit[a])
            oofp = out["uniprot_ac"].map(oof_prob[a])
            out[f"{a}_hit"] = pd.array(hit.astype("Float64").round(0), dtype="Int64")
            # Out-of-fold where labeled, the full-fit prediction where not.
            out[f"{a}_prob"] = np.where(oofp.notna(), oofp,
                                        out["uniprot_ac"].map(fresh)).round(4)
            out[f"{a}_source"] = np.where(hit.notna(), "measured", "predicted")

        train_acc = sorted({x for a in args.activator
                            for x in lab.loc[lab[a].notna(), "uniprot_ac"]} & set(feats.index))
        T = feats.loc[train_acc, cols].to_numpy(dtype=float)
        T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-9
        Q = f[cols].to_numpy(dtype=float)
        Q = Q / (np.linalg.norm(Q, axis=1, keepdims=True) + 1e-9)
        out["nn_similarity"] = out["uniprot_ac"].map(
            pd.Series((Q @ T.T).max(axis=1), index=f.index)).round(4)

        # `_hit` and `_source` stay as WORKING columns -- the manifest below is computed from
        # them -- but they are not written. Both are exactly reconstructible from
        # `evidence/labels_saureus.tsv` via `D.measured()`, and as columns they were empty for the
        # 10,131 Ec/Kp rows that have no measurement anywhere.

        # The two standard columns every axis ships. Both are derivations of what is already in
        # `out`, computed through `src/degradability.py` so the file and the loader cannot drift.
        out["degradability_consensus"] = D.consensus(out).round(6)
        out["degradability_evidence"] = D.evidence(out)

        path = (SCRATCH_DIR / f"smoke_degradability_{sp}.tsv") if args.limit \
            else (OUT_DIR / f"degradability_{sp}.tsv")
        out[D.OUT_COLUMNS].to_csv(path, sep="\t", index=False)

        # The ladder says only "measured", so the agree/disagree detail would otherwise be lost.
        # It is the one thing a reader might expect a 3 to mean, so keep it checkable per protein.
        n_meas_act = out["adep4_hit"].notna().astype(int) + out["onc212_hit"].notna().astype(int)
        audit_rows.append(pd.DataFrame({
            "species": sp,
            "uniprot_ac": out["uniprot_ac"],
            "adep4_measured": out["adep4_hit"],
            "onc212_measured": out["onc212_hit"],
            "n_measured": n_meas_act,
            # NA, not False, where fewer than two measurements exist: an activator a protein was
            # never tested against is missing data, not a dissenting source.
            "labels_agree": pd.array(
                np.where(n_meas_act == 2, out["adep4_hit"] == out["onc212_hit"], None),
                dtype="boolean"),
            "nn_similarity": out["nn_similarity"],
            "band_validated": out["nn_similarity"] >= D.VALIDATED_SIMILARITY,
            "degradability_consensus": out["degradability_consensus"],
            "degradability_evidence": out["degradability_evidence"],
        }))

        # Median over PREDICTED rows only: a training protein is its own nearest neighbour, so
        # including the measured ones would report a trivial 1.000 for S. aureus.
        pr_rows = out[(out["adep4_source"] == "predicted") | (out["onc212_source"] == "predicted")]
        med = float(pr_rows["nn_similarity"].median()) if len(pr_rows) else float("nan")
        n_meas = int(out["adep4_hit"].notna().sum() + out["onc212_hit"].notna().sum())
        lv = out["degradability_evidence"].value_counts()
        say(f"  {sp:<14} {len(out):>6} rows   measured cells {n_meas:>5}   "
            f"nn_similarity median (predicted) {med:.3f}   "
            f"evidence " + " ".join(f"L{k}={int(lv.get(k, 0)):,}" for k in (1, 2, 3)))
        rows.append({"species": sp, "n": len(out), "n_expected": len(P.load(sp)),
                     "measured_cells": n_meas,
                     "adep4_measured": int(out["adep4_hit"].notna().sum()),
                     "adep4_measured_pos": int((out["adep4_hit"] == 1).sum()),
                     "onc212_measured": int(out["onc212_hit"].notna().sum()),
                     "onc212_measured_pos": int((out["onc212_hit"] == 1).sum()),
                     "nn_similarity_median_predicted": round(med, 4) if med == med else None,
                     "activators": ",".join(args.activator),
                     "path": str(path.relative_to(REPO_ROOT)),
                     "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds")})
    pd.DataFrame(rows).to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}manifest.tsv", sep="\t", index=False)

    audit = pd.concat(audit_rows, ignore_index=True)
    audit.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}consensus_audit.tsv",
                 sep="\t", index=False)
    both = audit[audit["n_measured"] == 2]
    if len(both):
        say(f"  both activators measured {len(both):,}   labels agree "
            f"{int(both['labels_agree'].sum()):,} ({100 * both['labels_agree'].mean():.1f}%)"
            f"   -- NOT in the ladder, which reads level 3 as 'measured at all'")
    say()

    rule()
    say("SPOT CHECKS - v1's pre-registered predictions, which held")
    rule()
    sa_path = (SCRATCH_DIR / "smoke_degradability_saureus.tsv") if args.limit \
        else (OUT_DIR / "degradability_saureus.tsv")
    fails = []
    if sa_path.exists():
        names = P.load("saureus")[["uniprot_ac", "gene_name"]]
        sa = pd.read_csv(sa_path, sep="\t").merge(names, on="uniprot_ac")
        l2 = lab.set_index("uniprot_ac")
        for gene, note in SPOT_CHECKS.items():
            r = sa[sa["gene_name"].str.lower() == gene.lower()]
            if r.empty:
                say(f"  {gene:<6} not present under that gene name  (skipped)")
                continue
            r = r.iloc[0]
            fc = l2["onc212_log2fc"].get(r["uniprot_ac"], np.nan)
            hit = r["onc212_hit"]
            ok = (not pd.notna(hit)) or int(hit) == 0
            say(f"  {gene:<6} onc212 hit={('-' if pd.isna(hit) else int(hit))} "
                f"(log2FC {fc:>6.2f})  prob {r['onc212_prob']:.3f}  "
                f"adep4 hit={('-' if pd.isna(r['adep4_hit']) else int(r['adep4_hit']))}   "
                f"{'ok' if ok else '<- FAIL'}")
            say(f"  {'':<6} v1: {note}")
            if not ok:
                fails.append(f"{gene}: onc212_hit={int(hit)}, v1 says it is not a hit")
    say()
    if len(args.activator) == 2:
        from scipy.stats import spearmanr
        for sp in args.species:
            d = pd.read_csv((SCRATCH_DIR / f"smoke_degradability_{sp}.tsv") if args.limit
                            else (OUT_DIR / f"degradability_{sp}.tsv"), sep="\t")
            p = d[(d["adep4_source"] == "predicted") & (d["onc212_source"] == "predicted")]
            if len(p) > 50:
                r = float(spearmanr(p["adep4_prob"], p["onc212_prob"]).statistic)
                note = "  <- suspiciously high: both may be tracking length" if r > 0.90 else ""
                say(f"  predicted adep4 vs onc212 probability, {sp:<14} rho {r:+.3f}{note}")
    say()

    rule()
    say("OUTPUTS")
    rule()
    for r in rows:
        say(f"  {r['path']:<54} {(REPO_ROOT / r['path']).stat().st_size / 1e3:>8.1f} kB")
    say("  evidence/ + scratch/    labels_saureus.tsv, seqmap_audit.tsv, cv_<activator>.tsv,")
    say("                cross_activator.tsv, domain_bands.tsv, oof_<activator>_saureus.tsv,")
    say("                model_<activator>.npz, manifest.tsv")
    say()

    bad = [r for r in rows if r["n"] != r["n_expected"]]
    if bad:
        sys.exit("FAILED: " + ", ".join(
            f"{r['species']} has {r['n']} rows, expected {r['n_expected']}" for r in bad))
    # Label integrity: a `measured` source must carry a call, and every protein must have a
    # probability -- an empty `_prob` would silently drop a protein from any ranking.
    for sp in args.species:
        d = pd.read_csv((SCRATCH_DIR / f"smoke_degradability_{sp}.tsv") if args.limit
                        else (OUT_DIR / f"degradability_{sp}.tsv"), sep="\t")
        for a in args.activator:
            m = d[d[f"{a}_source"] == "measured"]
            if m[f"{a}_hit"].isna().any():
                sys.exit(f"FAILED: {sp}/{a} has 'measured' rows with no hit call.")
            if not d[f"{a}_hit"].dropna().isin([0, 1]).all():
                sys.exit(f"FAILED: {sp}/{a}_hit contains values outside {{0,1}}.")
            p = d[f"{a}_prob"]
            if p.isna().any() or not p.between(0, 1).all():
                sys.exit(f"FAILED: {sp}/{a}_prob has missing or out-of-range values.")
            # A saturated probability is an unrankable tie. At C=1.0 this hit 59% and made the
            # score useless for more than half the proteome; 20% is a generous ceiling.
            sat = _saturated(p)
            if sat > 0.20:
                sys.exit(f"FAILED: {100 * sat:.1f}% of {sp}/{a}_prob is pinned at 0 or 1, so the "
                         "score cannot rank those proteins. Regularisation is too weak.")
    if fails:
        sys.exit("FAILED spot checks:\n  " + "\n  ".join(fails))
    rule("=")
    say("stage 04 complete.")
    rule("=")


if __name__ == "__main__":
    main()
