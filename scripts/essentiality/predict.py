"""TabPFN per essentiality endpoint, and the cross-species test the axis has never had.

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

Run with the `gradi` env (TabPFN itself runs across a process boundary in `gradi-tabpfn`):
    python scripts/essentiality/predict.py --dry-run
    python scripts/essentiality/predict.py --endpoint keio_ess --features esmc
    python scripts/essentiality/predict.py --endpoint keio_ess --score-on ecl8_ess
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
from src import tabpfn as T  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"

N_FOLDS = 5
N_SEEDS = 5          # the axis-wide convention; a single seed carries ~+/-0.004 of arbitrariness
FEATURES = ("esmc", "prott5", "proteomelm")

# Which anchor proteome each endpoint's proteins live in -- i.e. whose embeddings to use.
ECL8 = "kpneumoniae__ecl8__GCA_000315385.1"
KPNIH1 = "kpneumoniae__kpnih1__GCA_000281535.2"
STRAIN_DIR = REPO_ROOT / "data" / "processed" / "embeddings" / "scratch" / "strains"

ENDPOINT_SPECIES = {"keio_ess": "ecoli", "goodall_ess": "ecoli", "bw25113_ess": "ecoli",
                    "conservation": "ecoli", "bn373_ess": ECL8, "kpnih1_ess": KPNIH1}

# `conservation` is a 3-class ordinal (non/mid/core), but the worker returns predict_proba[:, 1] --
# binary only. Rather than silently mis-read a 3-class probability, it is modelled as CORE vs REST.
# The raw `pct_essential` stays in the screen table, so the ordinal version costs no re-ingest.
BINARISE = {"conservation": lambda y: (y == 2).astype(int)}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def load_endpoint(name: str) -> pd.DataFrame:
    p = EVIDENCE_DIR / f"screen_{name}.tsv"
    if not p.exists():
        sys.exit(f"FATAL {p} missing -- run scripts/essentiality/screens.py first")
    d = pd.read_csv(p, sep="\t")
    if name in BINARISE:
        d["label"] = BINARISE[name](d["label"])
    return d


def features_for(species: str, accessions: list[str], kind: str) -> tuple[list[str], np.ndarray]:
    # A screen strain is not a registry species: its embeddings are training features under
    # scratch/, carry no canonical row order, and only ESM-C has been run for them.
    if species not in ("kpneumoniae", "ecoli", "saureus"):
        if kind != "esmc":
            sys.exit(f"FATAL only ESM-C exists for strain {species!r}; "
                     f"run scripts/embeddings/esmc.py --strain {species} for others")
        f = STRAIN_DIR / f"embeddings_{species}.npz"
        if not f.exists():
            sys.exit(f"FATAL {f} missing -- run scripts/embeddings/esmc.py --strain {species}")
        z = np.load(f, allow_pickle=True)
        accs = [str(a) for a in z["accessions"]]
        idx = {a: i for i, a in enumerate(accs)}
        keep = [a for a in accessions if a in idx]
        return keep, z["embeddings"][[idx[a] for a in keep]]
    if kind == "esmc":
        return E.vectors_for(species, accessions)
    if kind == "prott5":
        accs, mat = E.load_prott5(species)
        idx = {a: i for i, a in enumerate(accs)}
        keep = [a for a in accessions if a in idx]
        return keep, mat[[idx[a] for a in keep]]
    if kind == "proteomelm":
        return PLM.vectors_for(species, accessions)
    sys.exit(f"FATAL unknown feature set {kind!r}")


def groups_for(species: str, accessions: list[str]) -> np.ndarray:
    """Orthogroup per protein; an unassigned protein becomes its own group (stage 04's idiom).

    This is a LEAKAGE CONTROL. Paralogs are near-duplicates in embedding space, and a random split
    would put a protein's twin in the training fold and report a score the model did not earn.
    """
    if species not in ("kpneumoniae", "ecoli", "saureus"):
        # No OrthoFinder run covers a screen strain, so its groups come from an MMseqs2 clustering
        # of its own proteome at 30% id / 80% coverage -- stage 04's threshold, reused so the two
        # axes group on a comparable definition of "too similar to split". Measured on ECL8: 997 of
        # 5,178 proteins collapse into 520 paralog families, largest 22. Without this, a fifth of
        # the proteome could have its twin in the training fold.
        cl = (REPO_ROOT / "data" / "processed" / "essentiality" / "scratch" / "paralog_clusters"
              / f"{species}_id30_cov80.tsv")
        if not cl.exists():
            sys.exit(f"FATAL {cl} missing -- run\n"
                     f"  python scripts/essentiality/paralog_clusters.py --strain {species}")
        d = pd.read_csv(cl, sep="\t")
        m = dict(zip(d["member"].astype(str), d["cluster"].astype(str)))
        return np.array([m.get(a, "") or f"_self_{a}" for a in accessions])
    og = O.load(species)[["uniprot_ac", "orthogroup"]]
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


def run_endpoint(name: str, kind: str, folds: int, seeds: int, schemes: tuple[str, ...],
                 dry: bool) -> dict | None:
    species = ENDPOINT_SPECIES[name]
    d = load_endpoint(name)
    accs = d["uniprot_ac"].astype(str).tolist()
    found, X = features_for(species, accs, kind)
    sub = d.set_index("uniprot_ac").loc[found]
    y = sub["label"].to_numpy(dtype=int)
    groups = groups_for(species, found)

    say(f"  {name:14s} {kind:11s} X {X.shape[0]:5d} x {X.shape[1]:4d}   "
        f"pos {int(y.sum()):4d}  base {y.mean():.4f}  groups {len(set(groups)):5d}")
    if len(found) < len(accs):
        say(f"    {len(accs) - len(found)} of {len(accs)} proteins have no {kind} embedding "
            "-- dropped, never zero-filled")

    calls = []
    for scheme in schemes:
        for seed in range(seeds):
            for tr, te in splits_for(y, groups, folds, seed, scheme):
                calls.append((X[tr], y[tr], X[te]))
    n_un, credits = T.would_cost(calls)
    say(f"    {len(calls)} calls, {n_un} uncached -> {credits:,} credits")
    if dry:
        return None

    out = {}
    for scheme in schemes:
        aucs, prs, oofs = [], [], []
        for seed in range(seeds):
            sp = splits_for(y, groups, folds, seed, scheme)
            p = T.oof_predict(X, y, sp, seed=seed, tag=f"ess/{name}/{kind}/{scheme}/s{seed}")
            a, pr = score(y, p)
            aucs.append(a); prs.append(pr); oofs.append(p)
            say(f"      {scheme:8s} seed {seed}: AUROC {a:.4f}  PR {pr:.4f}")
        out[scheme] = (np.mean(aucs), np.std(aucs, ddof=1) if seeds > 1 else 0.0,
                       np.mean(prs), np.mean(oofs, axis=0))

    g = out["grouped"]
    pl = out.get("plain")
    rec = {"endpoint": name, "features": kind, "n": len(y), "n_pos": int(y.sum()),
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
    oof.to_csv(EVIDENCE_DIR / f"oof_{name}_{kind}.tsv", sep="\t", index=False)
    return rec


def cross_species(train: str, test: str, kind: str, dry: bool) -> dict | None:
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
    n_un, credits = T.would_cost([(X_tr, y_tr, X_te)])
    say(f"    {n_un} uncached call -> {credits:,} credits")
    if dry:
        return None

    p = T.predict_fold(X_tr, y_tr, X_te, tag=f"ess/cross/{train}->{test}/{kind}")
    auc, pr = score(y_te, p)
    lift = pr / y_te.mean()
    say(f"    AUROC {auc:.4f}   PR {pr:.4f}   base {y_te.mean():.4f}   lift {lift:.2f}x")
    return {"train": train, "test": test, "features": kind,
            "n_train": int(X_tr.shape[0]), "n_test": int(X_te.shape[0]),
            "base_train": round(float(y_tr.mean()), 4), "base_test": round(float(y_te.mean()), 4),
            "roc_auc": round(auc, 4), "pr_auc": round(pr, 4), "pr_lift": round(lift, 2)}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--endpoint", nargs="+", default=list(ENDPOINT_SPECIES),
                    choices=list(ENDPOINT_SPECIES))
    ap.add_argument("--features", nargs="+", default=["esmc"], choices=list(FEATURES))
    ap.add_argument("--folds", type=int, default=N_FOLDS)
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--schemes", nargs="+", default=["grouped"], choices=["grouped", "plain"],
                    help="`plain` measures the leakage gap and DOUBLES the cost -- run it on one "
                         "endpoint to calibrate, not on every combination")
    ap.add_argument("--score-on", nargs="+", choices=list(ENDPOINT_SPECIES),
                    help="fit on --endpoint and score on THESE endpoints' own labels and "
                         "embeddings -- the cross-species test; no orthology, no label transfer")
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
    say(f"  estimator  TabPFN-{T.MODEL} via gradi-tabpfn (non-commercial weights)")
    say(f"  cost       ~{T.CREDITS_PER_CALL:,} credits per uncached call, flat")
    rule()
    say("ENDPOINTS")
    rule()

    if a.score_on:
        rows = []
        for tr in a.endpoint:
            for te in a.score_on:
                if tr == te:
                    continue
                for kind in a.features:
                    r = cross_species(tr, te, kind, a.dry_run)
                    if r:
                        rows.append(r)
                    say("")
        if rows:
            EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
            out = EVIDENCE_DIR / "crossspecies.tsv"
            df = pd.DataFrame(rows)
            if out.exists():
                prev = pd.read_csv(out, sep="\t")
                k = ["train", "test", "features"]
                df = pd.concat([prev[~prev.set_index(k).index.isin(df.set_index(k).index)], df])
            df.to_csv(out, sep="\t", index=False)
            rule()
            say(df.to_string(index=False))
            say(f"\n  wrote {out.relative_to(REPO_ROOT)}")
        return

    recs = []
    for name in a.endpoint:
        for kind in a.features:
            r = run_endpoint(name, kind, a.folds, a.seeds, tuple(a.schemes), a.dry_run)
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
        prev = pd.read_csv(out, sep="\t")
        keys = ["endpoint", "features"]
        cv = pd.concat([prev[~prev.set_index(keys).index.isin(cv.set_index(keys).index)], cv])
    cv.to_csv(out, sep="\t", index=False)
    rule()
    say("RESULTS")
    rule()
    say(cv[["endpoint", "features", "n", "base_rate", "roc_auc_grouped",
            "pr_auc_grouped", "leakage_gap"]].to_string(index=False))
    say(f"\n  wrote {out.relative_to(REPO_ROOT)}")
    rule("=")


if __name__ == "__main__":
    main()
