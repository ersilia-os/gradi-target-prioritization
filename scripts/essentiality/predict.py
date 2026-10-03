"""One model per essentiality endpoint, and the cross-species test the axis has never had.

Each screen from `screens.py` is fitted INDEPENDENTLY. Nothing is merged: these are different
organisms, different assays and different base rates, and a merged label would be a number with no
referent.

THE POINT OF THIS SCRIPT is not another per-organism AUROC -- it is `--score-on`. Because ECL8 gives
us MEASURED Klebsiella labels, a model fitted on E. coli can be applied to K. pneumoniae embeddings
and scored against real Kp data. No label is ever copied between species: the cross-species step is
the model's, and this measures whether it works. Geptop's 0.59-0.81 came from two unrelated species
and never from our anchor.

TWO THINGS THAT WOULD SILENTLY INFLATE THE SCORE, both guarded:

  1. PARALOG LEAKAGE. Train and test are the same organism here, so the twins to separate are
     paralogs. Folds are grouped on the stage-05 orthogroup -- a LEAKAGE CONTROL, not a label
     transfer; no essentiality value crosses a species boundary anywhere in this script. Proteins
     with no orthogroup become their own group (stage 04's idiom). `leakage_gap` reports what the
     grouping costs, because an ungrouped CV is what v1's essentiality head did.
  2. PAYING TWICE. Hosted TabPFN bills ~10,000 credits per CALL, flat. `--dry-run` prices the whole
     run through `src.tabpfn.would_cost` before a single credit is spent, and a re-run that
     recomputes nothing spends nothing.

THE DEFAULT ESTIMATOR HERE IS THE FOREST, NOT TabPFN -- see DEFAULT_ESTIMATOR below for why, and
note that this is a deliberate exception to CLAUDE.md's standing rule, granted for this axis only
while the datasets are still being characterised. `--estimator tabpfn` switches, through a seam
narrow enough that nothing else about the evaluation changes.

Run with the `gradi` env (TabPFN, if selected, runs across a process boundary in `gradi-tabpfn`):
    python scripts/essentiality/predict.py
    python scripts/essentiality/predict.py --endpoint essential_ecoli_k12_knockout
    python scripts/essentiality/predict.py --estimator tabpfn --dry-run
    python scripts/essentiality/predict.py --endpoint essential_ecoli_k12_knockout \
        --score-on essential_kpneumoniae_ecl8_tradis
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import embeddings as E  # noqa: E402
from src import orthology as O  # noqa: E402
from src import proteomelm as PLM  # noqa: E402
from src import degradability as D  # noqa: E402
from src import matrices as M  # noqa: E402
from src import tabpfn as T  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from screens import SCREENS  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
TRAINING_DIR = OUT_DIR / "training_sets"

# THE DEFAULT ESTIMATOR ON THIS AXIS IS THE FOREST, on the project owner's instruction: TabPFN
# only once these datasets have earned it. CLAUDE.md's standing "always TabPFN-3.5" rule is the
# general case; this is an explicit exception for the essentiality endpoints while they are still
# being characterised, and it is not a licence to reach for sklearn elsewhere.
#
# It is also what makes the characterisation affordable. The forest is free, runs in-process, and
# needs no token; the same 10-endpoint 5-seed sweep under hosted TabPFN is 250 calls -- ~2.5M
# credits, 12.5% of the monthly quota -- which is not a price worth paying to find out which
# datasets are learnable at all.
DEFAULT_ESTIMATOR = "forest"
ESTIMATORS = ("forest", "tabpfn")

N_FOLDS = 5
N_SEEDS = 5          # the axis-wide convention; a single seed carries ~+/-0.004 of arbitrariness
# `proteomelm_orthodb` is the SAME model with the functional encoding the paper trained with (the
# mean ESM-C vector of each protein's OrthoDB group) instead of the released inference fallback
# (each protein's own vector). Two entries, never merged: the whole point is to compare them.
FEATURES = ("esmc", "prott5", "proteomelm", "proteomelm_orthodb")

# THE ENDPOINT LIST IS screens.py's REGISTRY, not a second copy of it. Every column, the proteome
# its rows live in and its one-line description come from there, so renaming a screen cannot leave
# this script pointing at a file that no longer exists.
ENDPOINT_SPECIES = {s.column: s.features for s in SCREENS}
ENDPOINT_COMMENT = {s.column: s.comment for s in SCREENS}

ANCHORS = ("kpneumoniae", "ecoli", "saureus")
STRAIN_DIR = REPO_ROOT / "data" / "processed" / "embeddings" / "scratch" / "strains"
CLUSTER_DIR = OUT_DIR / "scratch" / "paralog_clusters"

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def load_endpoint(name: str) -> pd.DataFrame:
    p = TRAINING_DIR / f"{name}.tsv"
    if not p.exists():
        sys.exit(f"FATAL {p} missing -- run scripts/essentiality/screens.py first")
    d = pd.read_csv(p, sep="\t")
    # The training sets key on `key`, which is NOT a UniProt accession for the five non-anchor
    # screens. The internal name stays `uniprot_ac` only because every downstream function here
    # already uses it as the row key; the value is whatever the screen's own proteome uses.
    return d.rename(columns={"key": "uniprot_ac"})


def features_for(species: str, accessions: list[str], kind: str) -> tuple[list[str], np.ndarray]:
    """Feature matrix for a set of proteins, from whichever proteome they actually live in.

    A SCREEN STRAIN IS NOT A REGISTRY SPECIES. Its embeddings are training features, not a
    deliverable matrix: they sit under `scratch/strains/`, carry no canonical row order, and are
    keyed on whatever identifier the strain's own FASTA uses (a `KPNRH_*` locus tag, a
    `lcl|...` DEG header). That is deliberate -- label, sequence and vector then share one
    namespace and nothing has to be joined across annotations.
    """
    if species not in ANCHORS:
        f = STRAIN_DIR / f"{'embeddings' if kind == 'esmc' else kind}_{species}.npz"
        if not f.exists():
            flag = "--strain " + species
            tool = "esmc.py" if kind == "esmc" else f"{kind}.py"
            sys.exit(f"FATAL {f} missing -- run scripts/embeddings/{tool} {flag}")
        z = np.load(f, allow_pickle=True)
        accs = [str(x) for x in z["accessions"]]
        idx = {x: i for i, x in enumerate(accs)}
        keep = [x for x in accessions if x in idx]
        return keep, z["embeddings"][[idx[x] for x in keep]]
    if kind == "esmc":
        return E.vectors_for(species, accessions)
    if kind == "prott5":
        accs, mat = E.load_prott5(species)
        idx = {a: i for i, a in enumerate(accs)}
        keep = [a for a in accessions if a in idx]
        return keep, mat[[idx[a] for a in keep]]
    if kind.startswith("proteomelm"):
        mode = "orthodb" if kind.endswith("_orthodb") else "self"
        # PLM.vectors_for asserts the file's recorded mode matches, so a `self` matrix can never be
        # served as `orthodb` -- the two are the same shape over the same accessions and nothing
        # about a matrix's appearance says which it is.
        return PLM.vectors_for(species, accessions, mode=mode)
    sys.exit(f"FATAL unknown feature set {kind!r}")


def groups_for(species: str, accessions: list[str]) -> np.ndarray:
    """Orthogroup per protein; an unassigned protein becomes its own group (stage 04's idiom).

    This is a LEAKAGE CONTROL. Paralogs are near-duplicates in embedding space, and a random split
    would put a protein's twin in the training fold and report a score the model did not earn.
    """
    if species not in ANCHORS:
        # No OrthoFinder run covers a screen strain, so its groups come from an MMseqs2 clustering
        # of its own proteome at 30% id / 80% coverage -- stage 04's threshold, reused so the two
        # axes group on a comparable definition of "too similar to split". Measured on ECL8: 997 of
        # 5,178 proteins collapse into 520 paralog families, largest 22. Without this, a fifth of
        # the proteome could have its twin in the training fold.
        cl = CLUSTER_DIR / f"{species}_id30_cov80.tsv"
        if not cl.exists():
            sys.exit(f"FATAL {cl} missing -- run\n"
                     f"  python scripts/essentiality/paralog_clusters.py --strain {species}")
        d = pd.read_csv(cl, sep="\t")
        m = dict(zip(d["member"].astype(str), d["cluster"].astype(str)))
        return np.array([m.get(a, "") or f"_self_{a}" for a in accessions])
    og = O.load_dense(species)[["uniprot_ac", "orthogroup"]]
    m = dict(zip(og["uniprot_ac"], og["orthogroup"].fillna("")))
    return np.array([m.get(a, "") or f"_self_{a}" for a in accessions])


def splits_for(y: np.ndarray, groups: np.ndarray, folds: int, seed: int, scheme: str):
    from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
    if scheme == "grouped":
        cv = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)
        sp = list(cv.split(np.zeros(len(y)), y, groups=groups))
        seen: dict = {}
        for k, (_, te) in enumerate(sp):
            for g in set(groups[te]):
                if g in seen:
                    sys.exit(f"FAILED seed {seed}: orthogroup {g!r} appears in folds {seen[g]} "
                             f"and {k}. The out-of-fold score would be leaky.")
                seen[g] = k
        return sp
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    return list(cv.split(np.zeros(len(y)), y))


def score(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    from sklearn.metrics import average_precision_score, roc_auc_score
    return float(roc_auc_score(y, p)), float(average_precision_score(y, p))


# ------------------------------------------------------------------------- the estimator seam
#
# ONE call site, so switching estimator cannot change anything else about the evaluation -- same
# folds, same grouping, same metric. That is what makes a later forest-vs-TabPFN comparison a
# comparison of estimators rather than of two scripts.

def predict_fold(X_tr, y_tr, X_te, *, estimator: str, seed: int, tag: str) -> np.ndarray:
    """Fit on (X_tr, y_tr), return P(essential) for X_te.

    The forest is stage 04's, imported rather than re-specified: `D.RF_PARAMS` is
    `n_estimators=500, max_features=0.1, min_samples_leaf=3, class_weight="balanced"`. Copying
    those numbers into a second file is how two axes silently drift apart, and `class_weight` is
    load-bearing here -- base rates run 0.048 to 0.194 across these ten datasets.
    """
    if estimator == "forest":
        from sklearn.ensemble import RandomForestClassifier
        m = RandomForestClassifier(n_jobs=-1, random_state=seed, **D.RF_PARAMS)
        m.fit(X_tr, y_tr)
        return m.predict_proba(X_te)[:, 1]
    return T.predict_fold(X_tr, y_tr, X_te, tag=tag)


def oof_predict(X, y, splits, *, estimator: str, seed: int, tag: str) -> np.ndarray:
    if estimator != "forest":
        return T.oof_predict(X, y, splits, seed=seed, tag=tag)
    prob = np.zeros(len(y), dtype=float)
    for k, (tr, te) in enumerate(splits):
        prob[te] = predict_fold(X[tr], y[tr], X[te], estimator=estimator, seed=seed,
                                tag=f"{tag}/f{k}")
    return prob


def run_endpoint(name: str, kind: str, folds: int, seeds: int, schemes: tuple[str, ...],
                 dry: bool, estimator: str = DEFAULT_ESTIMATOR) -> dict | None:
    species = ENDPOINT_SPECIES[name]
    d = load_endpoint(name)
    accs = d["uniprot_ac"].astype(str).tolist()
    found, X = features_for(species, accs, kind)
    sub = d.set_index("uniprot_ac").loc[found]
    y = sub["label"].to_numpy(dtype=int)
    groups = groups_for(species, found)

    say(f"  {name:42s} {kind:10s} X {X.shape[0]:5d} x {X.shape[1]:4d}   "
        f"pos {int(y.sum()):4d}  base {y.mean():.4f}  groups {len(set(groups)):5d}")
    if len(found) < len(accs):
        say(f"    {len(accs) - len(found)} of {len(accs)} proteins have no {kind} embedding "
            "-- dropped, never zero-filled")

    n_fits = len(schemes) * seeds * folds
    if estimator == "tabpfn":
        calls = [(X[tr], y[tr], X[te]) for scheme in schemes for seed in range(seeds)
                 for tr, te in splits_for(y, groups, folds, seed, scheme)]
        n_un, credits = T.would_cost(calls)
        say(f"    {len(calls)} calls, {n_un} uncached -> {credits:,} credits")
    else:
        say(f"    {n_fits} forest fits, in-process, no credits")
    if dry:
        return None

    out = {}
    for scheme in schemes:
        aucs, prs, oofs = [], [], []
        for seed in range(seeds):
            sp = splits_for(y, groups, folds, seed, scheme)
            p = oof_predict(X, y, sp, estimator=estimator, seed=seed,
                            tag=f"ess/{name}/{kind}/{scheme}/s{seed}")
            a, pr = score(y, p)
            aucs.append(a); prs.append(pr); oofs.append(p)
            say(f"      {scheme:8s} seed {seed}: AUROC {a:.4f}  PR {pr:.4f}")
        out[scheme] = (np.mean(aucs), np.std(aucs, ddof=1) if seeds > 1 else 0.0,
                       np.mean(prs), np.mean(oofs, axis=0))

    g = out["grouped"]
    pl = out.get("plain")
    rec = {"endpoint": name, "features": kind, "estimator": estimator,
           "n": len(y), "n_pos": int(y.sum()),
           "base_rate": round(float(y.mean()), 4), "n_features": X.shape[1],
           "n_folds": folds, "n_seeds": seeds,
           "roc_auc_grouped": round(g[0], 4), "roc_auc_grouped_sd": round(g[1], 4),
           "pr_auc_grouped": round(g[2], 4),
           "roc_auc_plain": round(pl[0], 4) if pl else pd.NA,
           "pr_auc_plain": round(pl[2], 4) if pl else pd.NA,
           "leakage_gap": round(pl[0] - g[0], 4) if pl else pd.NA}
    say(f"    GROUPED AUROC {g[0]:.4f} +/- {g[1]:.4f}   PR {g[2]:.4f}"
        + (f"   leakage_gap {rec['leakage_gap']:+.4f}" if pl else "   (plain not run)"))
    oof = pd.DataFrame({"uniprot_ac": found, "label": y, "oof_prob": g[3], "group": groups})
    oof.to_csv(EVIDENCE_DIR / f"oof_{name}_{kind}_{estimator}.tsv", sep="\t", index=False)
    return rec


def cross_species(train: str, test: str, kind: str, dry: bool,
                  estimator: str = DEFAULT_ESTIMATOR) -> dict | None:
    """Fit on one organism, predict another from ITS OWN embeddings, score on ITS OWN labels.

    THE point of this axis. No label is copied between organisms and no orthology is involved: the
    model sees organism A's proteins and their labels, then is shown organism B's embeddings and
    asked. Because ECL8 gives us MEASURED Klebsiella labels, `--train keio_ess --score-on bn373_ess`
    is the first honest answer this project has to "does an essentiality model reach our anchor".

    ONE TabPFN call, not a CV: `fit` is free and `predict` bills flat, so the whole test costs
    ~10,000 credits regardless of how many proteins are scored.
    """
    tr_sp, te_sp = ENDPOINT_SPECIES[train], ENDPOINT_SPECIES[test]
    dtr, dte = load_endpoint(train), load_endpoint(test)
    f_tr, X_tr = features_for(tr_sp, dtr["uniprot_ac"].astype(str).tolist(), kind)
    f_te, X_te = features_for(te_sp, dte["uniprot_ac"].astype(str).tolist(), kind)
    y_tr = dtr.set_index("uniprot_ac").loc[f_tr, "label"].to_numpy(dtype=int)
    y_te = dte.set_index("uniprot_ac").loc[f_te, "label"].to_numpy(dtype=int)

    say(f"  {train} ({tr_sp}) -> {test} ({te_sp})   features={kind}")
    say(f"    train {X_tr.shape[0]:5d} x {X_tr.shape[1]}  base {y_tr.mean():.4f}")
    say(f"    test  {X_te.shape[0]:5d} x {X_te.shape[1]}  base {y_te.mean():.4f}")
    if X_tr.shape[1] != X_te.shape[1]:
        sys.exit(f"FATAL feature dims differ ({X_tr.shape[1]} vs {X_te.shape[1]}) -- the two sides "
                 "must use the same embedding model or the comparison is meaningless")
    if estimator == "tabpfn":
        n_un, credits = T.would_cost([(X_tr, y_tr, X_te)])
        say(f"    {n_un} uncached call -> {credits:,} credits")
    if dry:
        return None

    p = predict_fold(X_tr, y_tr, X_te, estimator=estimator, seed=0,
                     tag=f"ess/cross/{train}->{test}/{kind}")
    auc, pr = score(y_te, p)
    lift = pr / y_te.mean()
    say(f"    AUROC {auc:.4f}   PR {pr:.4f}   base {y_te.mean():.4f}   lift {lift:.2f}x")
    return {"train": train, "test": test, "features": kind, "estimator": estimator,
            "n_train": int(X_tr.shape[0]), "n_test": int(X_te.shape[0]),
            "base_train": round(float(y_tr.mean()), 4), "base_test": round(float(y_te.mean()), 4),
            "roc_auc": round(auc, 4), "pr_auc": round(pr, 4), "pr_lift": round(lift, 2)}


def score_proteome(species: str, kind: str, estimator: str, seeds: int,
                   folds: int) -> pd.DataFrame:
    """Every screen's model applied to one ANCHOR proteome -> `screens_<species>.tsv`.

    THESE ARE PREDICTIONS IN EVERY COLUMN AND EVERY ROW, and the file is named and documented so
    that cannot be misread. Each screen was measured on a strain that is not the anchor -- measured,
    five of the nine key onto the anchor at exactly 0 of 4,930 / 4,809 / 4,981 / 4,981 / 5,433 -- so
    transferring the labels was rejected in favour of transferring the MODEL. One comparable 0-1
    scale across all nine columns, and no measurement is implied anywhere.

    OUT-OF-FOLD WHERE THE PROTEIN WAS IN THAT SCREEN'S TRAINING SET. The four b-number E. coli
    screens sit 100% on the E. coli anchor, so a plain fit-then-predict would be scoring its own
    training data and the column would be optimistic exactly where it overlaps and honest
    elsewhere -- two different quantities in one column. The existing
    `evidence/oof_<endpoint>_<features>_<estimator>.tsv` values are substituted in for those rows,
    which is what keeps the column on one scale.

    The output is COMPLETE and CANONICAL: one row per anchor protein, in proteome order.
    """
    import numpy as np

    accs = list(M.canonical(species))
    found, X_target = features_for(species, accs, kind)
    if len(found) != len(accs):
        sys.exit(f"FATAL {species}: {len(accs) - len(found)} proteins have no {kind} embedding, so "
                 "the column could not be complete. Fill them where the embedding is built.")
    out = pd.DataFrame({"uniprot_ac": found})
    audit = []

    for s in SCREENS:
        d = load_endpoint(s.column)
        tr_accs = d["uniprot_ac"].astype(str).tolist()
        f_tr, X_tr = features_for(s.features, tr_accs, kind)
        y_tr = d.set_index("uniprot_ac").loc[f_tr, "label"].to_numpy(dtype=int)
        if X_tr.shape[1] != X_target.shape[1]:
            sys.exit(f"FATAL {s.column}: feature dims differ ({X_tr.shape[1]} vs "
                     f"{X_target.shape[1]}) -- both sides must use the same embedding model")

        probs = np.zeros((seeds, len(found)), dtype=float)
        for seed in range(seeds):
            probs[seed] = predict_fold(X_tr, y_tr, X_target, estimator=estimator, seed=seed,
                                       tag=f"ess/proteome/{s.column}->{species}/{kind}/s{seed}")
        col = probs.mean(axis=0)

        # substitute out-of-fold values wherever this anchor protein WAS in the training set
        oof_p = EVIDENCE_DIR / f"oof_{s.column}_{kind}_{estimator}.tsv"
        n_oof = 0
        if oof_p.exists():
            o = pd.read_csv(oof_p, sep="\t")
            lut = dict(zip(o["uniprot_ac"].astype(str), o["oof_prob"]))
            hit = np.array([a in lut for a in found])
            if hit.any():
                col[hit] = [lut[a] for a, h in zip(found, hit) if h]
                n_oof = int(hit.sum())
        overlap = n_oof / len(found)
        say(f"  {s.column:42s} -> {species:12s} mean {col.mean():.4f}  "
            f"oof-substituted {n_oof:5d} ({overlap:.1%})")

        out[s.column] = col.round(4)
        cv = EVIDENCE_DIR / "cv_endpoints.tsv"
        own = {}
        if cv.exists():
            c = pd.read_csv(cv, sep="\t")
            c = c[(c.endpoint == s.column) & (c.features == kind) & (c.estimator == estimator)]
            if len(c):
                own = c.iloc[0].to_dict()
        audit.append({"column": s.column, "scored_proteome": species,
                      "trained_on": s.features, "organism": s.organism, "assay": s.assay,
                      "features": kind, "estimator": estimator, "n_seeds": seeds,
                      "n_train": int(X_tr.shape[0]), "train_base_rate": round(float(y_tr.mean()), 4),
                      "n_scored": len(found),
                      "own_organism_roc_auc": own.get("roc_auc_grouped"),
                      "own_organism_pr_auc": own.get("pr_auc_grouped"),
                      "oof_substituted": n_oof,
                      "train_target_overlap": round(overlap, 4),
                      "is_same_species": s.features == species
                      or s.features.startswith(species[:2]) and species in s.features})
    return out, pd.DataFrame(audit)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--endpoint", nargs="+", default=list(ENDPOINT_SPECIES),
                    choices=list(ENDPOINT_SPECIES))
    ap.add_argument("--features", nargs="+", default=["prott5"], choices=list(FEATURES))
    ap.add_argument("--estimator", default=DEFAULT_ESTIMATOR, choices=list(ESTIMATORS),
                    help="forest (default) runs in-process and free; tabpfn bills ~10,000 credits "
                         "per uncached call. The default is the project owner's instruction for "
                         "this axis, not a general preference -- see CLAUDE.md.")
    ap.add_argument("--folds", type=int, default=N_FOLDS)
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--schemes", nargs="+", default=["grouped"], choices=["grouped", "plain"],
                    help="`plain` measures the leakage gap and DOUBLES the cost -- run it on one "
                         "endpoint to calibrate, not on every combination")
    ap.add_argument("--score-on", nargs="+", choices=list(ENDPOINT_SPECIES),
                    help="fit on --endpoint and score on THESE endpoints' own labels and "
                         "embeddings -- the cross-species test; no orthology, no label transfer")
    ap.add_argument("--score-proteome", nargs="+", metavar="SPECIES",
                    choices=["ecoli", "kpneumoniae", "saureus"],
                    help="apply EVERY screen's model to these anchor proteomes and write "
                         "screens_<species>.tsv -- one predicted probability column per screen, "
                         "complete and in canonical row order. Predictions throughout; see the "
                         "transfer audit before reading any column as evidence.")
    ap.add_argument("--dry-run", action="store_true", help="price the run; spend nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("scripts/essentiality/predict.py -- TabPFN per endpoint, orthogroup-grouped CV")
    rule("=")
    say(f"  endpoints  {', '.join(a.endpoint)}")
    say(f"  features   {', '.join(a.features)}")
    say(f"  cv         StratifiedGroupKFold({a.folds}) on stage-05 orthogroups x {a.seeds} seeds")
    say(f"  schemes    {', '.join(a.schemes)}")
    say(f"  estimator  " + ("RandomForest(**src.degradability.RF_PARAMS), in-process"
                              if a.estimator == "forest"
                              else f"TabPFN-{T.MODEL} via gradi-tabpfn (non-commercial weights)"))
    say("  cost       none -- the forest is free" if a.estimator == "forest"
        else f"  cost       ~{T.CREDITS_PER_CALL:,} credits per uncached call, flat")
    rule()
    say("ENDPOINTS")
    rule()

    if a.score_proteome:
        for kind in a.features:
            for sp in a.score_proteome:
                rule()
                say(f"SCORE PROTEOME  {sp}  ({kind}, {a.estimator})")
                rule()
                tbl, aud = score_proteome(sp, kind, a.estimator, a.seeds, a.folds)
                if a.dry_run:
                    continue
                tbl = M.reindex(tbl, sp)
                outp = OUT_DIR / f"screens_{sp}.tsv"
                tbl.to_csv(outp, sep="\t", index=False)
                ap_path = EVIDENCE_DIR / "screens_transfer_audit.tsv"
                if ap_path.exists():
                    prev = pd.read_csv(ap_path, sep="\t")
                    k = ["column", "scored_proteome", "features", "estimator"]
                    aud = pd.concat([prev[~prev.set_index(k).index.isin(aud.set_index(k).index)],
                                     aud])
                aud.to_csv(ap_path, sep="\t", index=False)
                say(f"  wrote {outp.relative_to(REPO_ROOT)}  ({tbl.shape[0]} x {tbl.shape[1]})")
                say(f"  wrote {ap_path.relative_to(REPO_ROOT)}")
        rule("=")
        return

    if a.score_on:
        rows = []
        for tr in a.endpoint:
            for te in a.score_on:
                if tr == te:
                    continue
                for kind in a.features:
                    r = cross_species(tr, te, kind, a.dry_run, a.estimator)
                    if r:
                        rows.append(r)
                    say("")
        if rows:
            EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
            out = EVIDENCE_DIR / "crossspecies.tsv"
            df = pd.DataFrame(rows)
            if out.exists():
                prev = pd.read_csv(out, sep="\t")
                k = [c for c in ("train", "test", "features", "estimator") if c in prev.columns]
                df = pd.concat([prev[~prev.set_index(k).index.isin(df.set_index(k).index)], df])
            df.to_csv(out, sep="\t", index=False)
            rule()
            say(df.to_string(index=False))
            say(f"\n  wrote {out.relative_to(REPO_ROOT)}")
        return

    recs = []
    for name in a.endpoint:
        for kind in a.features:
            r = run_endpoint(name, kind, a.folds, a.seeds, tuple(a.schemes), a.dry_run,
                             a.estimator)
            if r:
                recs.append(r)
            say("")

    if a.dry_run:
        rule("=")
        say("dry run: nothing computed, nothing spent.")
        return
    if not recs:
        return
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    cv = pd.DataFrame(recs)
    cv["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    out = EVIDENCE_DIR / "cv_endpoints.tsv"
    if out.exists():
        # A results file written before `estimator` existed would KeyError here and throw away a
        # completed run at the very last step -- which it did once. Merge on the keys both frames
        # actually have, never on the keys this version happens to write.
        prev = pd.read_csv(out, sep="\t")
        keys = [k for k in ("endpoint", "features", "estimator") if k in prev.columns]
        cv = pd.concat([prev[~prev.set_index(keys).index.isin(cv.set_index(keys).index)], cv])
    cv.to_csv(out, sep="\t", index=False)
    rule()
    say("RESULTS")
    rule()
    say(cv[["endpoint", "features", "estimator", "n", "base_rate", "roc_auc_grouped",
            "pr_auc_grouped", "leakage_gap"]].to_string(index=False))
    say(f"\n  wrote {out.relative_to(REPO_ROOT)}")
    rule("=")


if __name__ == "__main__":
    main()
