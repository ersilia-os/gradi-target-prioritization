"""ESM-C predictor of activated-ClpP susceptibility: train on S. aureus, apply to E. coli (docs §3.3c).

WHY THIS EXISTS, AND WHAT IT CANNOT DO
--------------------------------------
`10e` pre-registered a gate — "does any sequence feature beat protein length alone?" — and the answer
was **no**: disorder adds +0.009 over length, a gradient-boosted model on all 13 hand-built features
scores 0.762 against logistic length-only at 0.775, and degron motifs attribute 0%. On that basis this
script was initially NOT built. It is built now at explicit request, and the case for it is real: the
`10e` cross-activator test showed abundance susceptibility **does** generalise between two chemically
distinct activators (AUROC 0.692 / 0.765), so there is learnable signal. The open question is whether a
1152-dim language-model embedding captures more of it than length + composition do.

Two hard limits, stated because the output must not be over-read:

1. **The E. coli predictions cannot be validated.** No E. coli activated-ClpP labels exist — that is the
   gap in `docs/degradability_datasets.md` §10. Every E. coli row is therefore a ranking hypothesis, and
   the output file carries a `validation` column saying so.
2. **ESM embeddings encode species and length.** A head trained on S. aureus partly learns
   "S. aureus-ness", and mean-pooled embeddings correlate with length — the very feature that dominates
   here. So the honest test is not "is the AUROC high" but "does it beat length and the simple features
   under homology-aware CV". That comparison is the point of this script.

DESIGN
------
* **Cluster-grouped CV** on the MMseqs2 30%-identity clusters from `10e` (1,665 clusters over 1,943
  sequences). Plain CV is reported alongside so the leakage gap is visible — `07d` uses plain `cv=5`
  for its essentiality head, and this quantifies what that costs.
* **Every baseline re-run under the identical splits**: length only, simple features, ESM only,
  ESM + simple. A model is only credited with what it adds over the cheaper one.
* **Both targets, never pooled** (abundance and cleavage), because they behave differently.
* **Cross-activator held-out** for the ESM head too, since that is the only non-circular test available.

Outputs
-------
  output/results/other/activator_esm_cv.csv          model x target AUROC, clustered and plain
  output/results/other/activator_esm_cross.csv       cross-activator generalisation for the ESM head
  output/results/ecoli/ec_activator_predicted.csv    E. coli ranking hypothesis (UNVALIDATED)

Run with the `gradi` env, after 10e and after embedding the S. aureus set:
    python scripts/01a_esmc_embeddings.py \\
        --fasta data/processed/other/degradability/activator/saureus_activator_proteins.faa \\
        --out   data/processed/other/degradability/activator/saureus_esmc600m.npz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

RES = D.REPO_ROOT / "output" / "results" / "other"
SA_NPZ = D.xbac_processed_dir("activator") / "saureus_esmc600m.npz"
RNG_SEED = 0
N_SPLITS = 5
SIMPLE = ["log_length", "disorder_mean", "frac_charged", "gravy", "frac_PEST"]
TARGETS = ["conlon_abundance_bin", "jacques_abundance_bin",
           "conlon_cleavage_bin", "jacques_cleavage_bin"]


def _pipe():
    # Logistic with L2 on standardised features. Deliberately the SAME estimator for ESM and for the
    # baselines, so a difference in AUROC is attributable to the features and not to model capacity.
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0))


def load_inputs() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    f = RES / "activator_features.csv"
    if not f.exists():
        raise SystemExit(f"missing {f} — run scripts/10e_activator_features.py first")
    feats = pd.read_csv(f)
    if not SA_NPZ.exists():
        raise SystemExit(f"missing {SA_NPZ} — embed the S. aureus set first (see this file's docstring)")
    z = np.load(SA_NPZ, allow_pickle=True)
    accs = [str(a) for a in z["accessions"]]
    emb = z["embeddings"].astype(np.float32)
    print(f"[in] features {feats.shape}, embeddings {emb.shape}")
    return feats, emb, accs


def assemble(feats: pd.DataFrame, emb: np.ndarray, accs: list[str]) -> pd.DataFrame:
    idx = {a: i for i, a in enumerate(accs)}
    feats = feats[feats["saureus_acc"].isin(idx)].copy()
    feats["_row"] = feats["saureus_acc"].map(idx)
    print(f"[join] {len(feats)} proteins have both an embedding and features")
    return feats


def evaluate(df: pd.DataFrame, emb: np.ndarray) -> pd.DataFrame:
    """AUROC per (model, target) under clustered and plain CV, on identical row subsets."""
    rows = []
    for tgt in TARGETS:
        sub = df.dropna(subset=[tgt] + SIMPLE)
        y = sub[tgt].astype(int).to_numpy()
        if y.sum() < 25 or len(set(y)) < 2:
            continue
        E = emb[sub["_row"].to_numpy()]
        S = sub[SIMPLE].to_numpy()
        blocks = {"length only": sub[["log_length"]].to_numpy(),
                  "simple (5 feat)": S,
                  "ESM-C only": E,
                  "ESM-C + simple": np.hstack([E, S])}
        groups = sub["cluster"].to_numpy() if "cluster" in sub.columns else None
        for name, X in blocks.items():
            r = {"target": tgt, "model": name, "n": len(sub), "n_pos": int(y.sum()),
                 "n_features": X.shape[1]}
            for tag, cv, kw in (("clustered", GroupKFold(n_splits=N_SPLITS),
                                 {"groups": groups} if groups is not None else {}),
                                ("plain", KFold(N_SPLITS, shuffle=True, random_state=RNG_SEED), {})):
                if tag == "clustered" and groups is None:
                    r["auroc_clustered"] = np.nan
                    continue
                p = cross_val_predict(_pipe(), X, y, cv=cv, **kw, method="predict_proba")[:, 1]
                r[f"auroc_{tag}"] = round(roc_auc_score(y, p), 4)
            r["leakage_gap"] = round(r.get("auroc_plain", np.nan) - r.get("auroc_clustered", np.nan), 4)
            rows.append(r)
            print(f"    {tgt:24} {name:16} clustered {r.get('auroc_clustered', float('nan')):.3f} "
                  f"plain {r.get('auroc_plain', float('nan')):.3f}", flush=True)
    return pd.DataFrame(rows)


def cross_activator(df: pd.DataFrame, emb: np.ndarray) -> pd.DataFrame:
    """Train on one activator's labels, evaluate against the other's — for the ESM head."""
    pairs = [("conlon_abundance_bin", "jacques_abundance_bin"),
             ("jacques_abundance_bin", "conlon_abundance_bin"),
             ("conlon_cleavage_bin", "jacques_cleavage_bin"),
             ("jacques_cleavage_bin", "conlon_cleavage_bin")]
    rows = []
    for tr_lab, te_lab in pairs:
        sub = df.dropna(subset=[tr_lab, te_lab] + SIMPLE)
        ytr = sub[tr_lab].astype(int).to_numpy(); yte = sub[te_lab].astype(int).to_numpy()
        if len(sub) < 60 or ytr.sum() < 15 or len(set(yte)) < 2:
            continue
        groups = sub["cluster"].to_numpy() if "cluster" in sub.columns else None
        cv, kw = ((GroupKFold(n_splits=N_SPLITS), {"groups": groups}) if groups is not None
                  else (KFold(N_SPLITS, shuffle=True, random_state=RNG_SEED), {}))
        for name, X in (("ESM-C + simple", np.hstack([emb[sub["_row"].to_numpy()],
                                                     sub[SIMPLE].to_numpy()])),
                        ("simple (5 feat)", sub[SIMPLE].to_numpy())):
            oof = cross_val_predict(_pipe(), X, ytr, cv=cv, **kw, method="predict_proba")[:, 1]
            rows.append({"train_on": tr_lab, "evaluate_on": te_lab, "model": name,
                         "n_shared": len(sub),
                         "auroc_same_experiment": round(roc_auc_score(ytr, oof), 4),
                         "auroc_other_experiment": round(roc_auc_score(yte, oof), 4)})
    return pd.DataFrame(rows)


def apply_to_ecoli(df: pd.DataFrame, emb: np.ndarray, target: str) -> pd.DataFrame:
    """Fit on all S. aureus, score every E. coli protein. UNVALIDATED by construction."""
    sub = df.dropna(subset=[target] + SIMPLE)
    y = sub[target].astype(int).to_numpy()
    X = np.hstack([emb[sub["_row"].to_numpy()], sub[SIMPLE].to_numpy()])
    pipe = _pipe(); pipe.fit(X, y)

    z = np.load(D.processed_dir("ecoli", "embeddings") / "ec_esmc600m_embeddings.npz",
                allow_pickle=True)
    ec_acc = [str(a) for a in z["accessions"]]
    ec_emb = z["embeddings"].astype(np.float32)
    # The E. coli side needs the SAME simple features in the SAME order. They are recomputed from the
    # E. coli proteome rather than transferred, so nothing depends on orthology here.
    ec_simple = ecoli_simple_features()
    ec_simple = ec_simple[ec_simple["uniprot_accession"].isin(ec_acc)].reset_index(drop=True)
    row = {a: i for i, a in enumerate(ec_acc)}
    Xe = np.hstack([ec_emb[[row[a] for a in ec_simple["uniprot_accession"]]],
                    ec_simple[SIMPLE].to_numpy()])
    score = pipe.predict_proba(Xe)[:, 1]
    out = pd.DataFrame({
        "uniprot_accession": ec_simple["uniprot_accession"],
        "activator_predicted_score": np.round(score, 4),
        "activator_predicted_pct": np.round(100 * pd.Series(score).rank(pct=True), 1),
        "trained_on": target,
        "validation": "none — cross-species, unlabelled",
    })
    return out.sort_values("activator_predicted_score", ascending=False).reset_index(drop=True)


def ecoli_simple_features(refresh: bool = False) -> pd.DataFrame:
    """The five simple features for the E. coli proteome, computed the SAME way as for S. aureus.

    Deliberately reuses 10e's own feature functions (imported, not copied) so the training and
    application sides cannot drift. In particular disorder must come from **metapredict**, not from the
    pLDDT-derived disorder in `10d`: the head is trained on metapredict values, and scoring it with a
    differently-scaled disorder column would silently corrupt every prediction.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_activator_features", Path(__file__).with_name("10e_activator_features.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    seqs: dict[str, str] = {}
    acc = None
    for line in D.proteome_fasta("ecoli").read_text().splitlines():
        if line.startswith(">"):
            acc = D.acc_from_header(line)
            seqs[acc] = ""
        elif acc:
            seqs[acc] += line.strip()
    seqs = {a: v for a, v in seqs.items() if v}
    print(f"  [ecoli] computing simple features for {len(seqs)} proteins "
          f"(metapredict, cached)…", flush=True)
    cache = D.degradability_processed_dir("ecoli") / "ec_metapredict_disorder.json"
    dis = mod.disorder_features(seqs, refresh=refresh, cache=cache, id_col="uniprot_accession")
    seq = mod.sequence_features(seqs, id_col="uniprot_accession")
    return dis.merge(seq, on="uniprot_accession", how="inner")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-ecoli", action="store_true",
                    help="evaluate on S. aureus only; do not emit the unvalidated E. coli ranking")
    args = ap.parse_args()

    feats, emb, accs = load_inputs()
    df = assemble(feats, emb, accs)

    print("\n[cv] AUROC by model x target (clustered vs plain CV)…", flush=True)
    cvres = evaluate(df, emb)
    report_cv(cvres)

    print("\n[cross] cross-activator generalisation…", flush=True)
    cx = cross_activator(df, emb)
    report_cross(cx)

    RES.mkdir(parents=True, exist_ok=True)
    if not args.skip_ecoli:
        # Trained on the target with the most support and the best cross-activator transfer.
        print("\n[ecoli] scoring the E. coli proteome (UNVALIDATED — no E. coli labels exist)…",
              flush=True)
        ec = apply_to_ecoli(df, emb, "conlon_abundance_bin")
        dest = D.results_dir("ecoli") / "ec_activator_predicted.csv"
        ec.to_csv(dest, index=False)
        print(f"[out] {dest.relative_to(D.REPO_ROOT)}  ({len(ec)} rows) — "
              f"carries validation='none'")
    cvres.to_csv(RES / "activator_esm_cv.csv", index=False)
    cx.to_csv(RES / "activator_esm_cross.csv", index=False)
    print(f"\n[out] {(RES / 'activator_esm_cv.csv').relative_to(D.REPO_ROOT)}")
    print(f"[out] {(RES / 'activator_esm_cross.csv').relative_to(D.REPO_ROOT)}")


def report_cv(cv: pd.DataFrame) -> None:
    if cv.empty:
        return
    print(f"\n{'='*100}\nESM-C vs BASELINES — identical rows, identical estimator, "
          f"{N_SPLITS}-fold CV\n")
    for tgt, sub in cv.groupby("target", sort=False):
        base = sub.loc[sub["model"] == "length only", "auroc_clustered"]
        base = float(base.iloc[0]) if len(base) else np.nan
        print(f"--- {tgt}  (n={sub['n'].iloc[0]}, positives={sub['n_pos'].iloc[0]})")
        print(f"    {'model':18} {'feats':>6} {'clustered':>10} {'plain':>8} {'leak':>7} {'vs length':>10}")
        for _, r in sub.iterrows():
            delta = r["auroc_clustered"] - base
            print(f"    {r['model']:18} {r['n_features']:>6} {r['auroc_clustered']:>10.3f} "
                  f"{r['auroc_plain']:>8.3f} {r['leakage_gap']:>+7.3f} {delta:>+10.3f}")
        print()
    esm = cv[cv["model"] == "ESM-C + simple"].set_index("target")["auroc_clustered"]
    smp = cv[cv["model"] == "simple (5 feat)"].set_index("target")["auroc_clustered"]
    gain = (esm - smp).dropna()
    print(f"{'='*100}\nESM-C + simple minus simple alone: "
          + ", ".join(f"{t.replace('_bin','')} {v:+.3f}" for t, v in gain.items()))
    if gain.max() < 0.02:
        print("VERDICT: the 1152-dim embedding adds nothing meaningful over five hand-built features.")
    else:
        print("VERDICT: the embedding adds signal over the hand-built features on at least one target.")
    leak = cv["leakage_gap"].dropna()
    if len(leak):
        print(f"\nHomology leakage (plain CV minus clustered CV): mean {leak.mean():+.3f}, "
              f"max {leak.max():+.3f} — this is what `07d`'s plain cv=5 would over-report.")


def report_cross(cx: pd.DataFrame) -> None:
    if cx.empty:
        print("  not enough shared labelled proteins")
        return
    print(f"\n{'='*100}\nCROSS-ACTIVATOR GENERALISATION — the only non-circular test\n")
    print(f"  {'train on':24} {'model':16} {'same expt':>10} {'OTHER expt':>11}")
    for _, r in cx.iterrows():
        print(f"  {r['train_on']:24} {r['model']:16} {r['auroc_same_experiment']:>10.3f} "
              f"{r['auroc_other_experiment']:>11.3f}")


if __name__ == "__main__":
    main()
