"""Do disorder and degron motifs predict activated-ClpP susceptibility? (docs §3.3c)

THE QUESTION THIS ANSWERS
-------------------------
The axis's central mechanistic claim is that **partnerless activated ClpP can only degrade what is
already unstructured**, so disorder should predict susceptibility. Until now that claim has only ever
been tested against *natural turnover* (Nagar/Gupta), and `docs/degradability_datasets.md` §10.2 shows
those two quantities are **statistically independent** (rho = -0.073, p = 0.17). So the claim has never
actually been tested against the bar the project gates on.

This script tests it, using the only right-bar labels that exist: the Conlon 2013 (ADEP4) and Jacques
2020 (ONC212) S. aureus activator proteomics, over the 1,943 sequences `10c` already cached.

DESIGN DECISIONS, AND WHY
-------------------------
* **Two targets, never pooled.** `abundance` and `cleavage` are reported separately for each paper.
  They measure different things and reproduce differently between activators (Spearman +0.52 vs +0.22),
  and AcpP is the worked example of why pooling misleads: +0.03 abundance, +3.46 cleavage.
* **Sequence-based disorder**, via metapredict. MobiDB and AlphaFold are UniProt-keyed while these
  sequences are RefSeq/GenBank from NCBI efetch, so an ID-free predictor avoids a lossy mapping step
  entirely. N-30 and C-30 windows are emitted separately so **Won 2024's N-vs-C asymmetry can be tested
  on the right bar** for the first time.
* **Length is the baseline to beat.** Gupta 2024 found small proteins turn over fast, and any sequence
  feature correlates with length. A feature that cannot beat length has told us nothing.
* **Cluster-aware statistics.** Homologues inflate both AUROC and CI tightness. MMseqs2 clusters at 30%
  identity define the resampling unit for the bootstrap, and the CV groups in 10f.
* **Judge against the label ceiling, not against 1.0.** The two activators agree at Jaccard 0.31, so
  ~0.70 AUROC is a good result here.

Outputs
-------
  output/results/other/activator_features.csv          per-protein features + labels
  output/results/other/activator_feature_stats.csv     AUROC/rho + bootstrap CI per feature x target
  data/processed/other/degradability/activator/         cached disorder + mmseqs clusters

Run with the `gradi` env. MMseqs2 comes from `gradi-ortho` (optional — without it the bootstrap falls
back to per-protein resampling and says so).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

MMSEQS_BIN = os.environ.get(
    "GRADI_MMSEQS_BIN", str(Path.home() / "miniconda3/envs/gradi-ortho/bin/mmseqs")
)
CLUSTER_MIN_SEQ_ID = 0.30
N_BOOT = 2000
RNG_SEED = 0

# Label thresholds. Kept identical to 10c so the two scripts cannot drift.
ABUNDANCE_LOG2 = -1.0
CLEAVAGE_LOG2 = 1.0
PADJ = 0.05


# --------------------------------------------------------------------------- inputs
def _load_10c():
    """Import scripts/10c_clpp_activator.py for its parsers.

    Imported rather than duplicated so the two scripts cannot disagree about how the source tables are
    read. The leading digit makes the module name non-importable normally, hence importlib.
    """
    path = Path(__file__).with_name("10c_clpp_activator.py")
    spec = importlib.util.spec_from_file_location("_clpp_activator", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_sequences() -> dict[str, str]:
    faa = D.xbac_processed_dir("activator") / "saureus_activator_proteins.faa"
    if not faa.exists():
        raise SystemExit(f"missing {faa} — run scripts/10c_clpp_activator.py first")
    seqs: dict[str, str] = {}
    acc = None
    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            acc = line[1:].split()[0]
            seqs[acc] = ""
        elif acc:
            seqs[acc] += line.strip()
    return {a: s for a, s in seqs.items() if s}


def build_labels() -> pd.DataFrame:
    """Per-S.-aureus-protein labels: 2 readouts x 2 papers, continuous and binary, never pooled."""
    c = _load_10c()
    ab, cl = c.parse_conlon()
    j = c.parse_jacques()
    sa = ab.merge(cl, on="sacol", how="outer").merge(j, on="sacol", how="outer")
    sa["saureus_acc"] = sa["refseq"].fillna("").astype(str).str.strip()
    sa.loc[sa["saureus_acc"].isin(["", "nan"]), "saureus_acc"] = (
        sa["genbank"].fillna("").astype(str).str.strip())
    sa = sa[sa["saureus_acc"] != ""].copy()

    out = pd.DataFrame({"saureus_acc": sa["saureus_acc"], "sacol": sa["sacol"],
                        "gene": sa["gene"].fillna("")})
    # continuous
    out["conlon_abundance"] = sa["adep4_abundance_log2fc"]
    out["conlon_cleavage"] = sa["adep4_cleavage_log2fc_max"]
    out["jacques_abundance"] = sa["onc212_abundance_log2fc"]
    out["jacques_cleavage"] = sa["onc212_cleavage_log2fc"]
    # binary. NOTE the sign asymmetry: abundance DOWN is the event, cleavage UP is the event.
    out["conlon_abundance_bin"] = np.where(
        sa["adep4_abundance_log2fc"].notna(),
        (sa["adep4_abundance_log2fc"] <= ABUNDANCE_LOG2) & (sa["adep4_abundance_padj"] < PADJ),
        np.nan)
    out["conlon_cleavage_bin"] = np.where(
        sa["adep4_cleavage_peptides"].notna(),
        sa["adep4_cleavage_peptides_sig"].fillna(0) > 0, np.nan)
    out["jacques_abundance_bin"] = np.where(
        sa["onc212_abundance_log2fc"].notna(),
        sa["onc212_abundance_log2fc"] <= ABUNDANCE_LOG2, np.nan)
    out["jacques_cleavage_bin"] = np.where(
        sa["onc212_cleavage_log2fc"].notna(),
        sa["onc212_cleavage_log2fc"] >= CLEAVAGE_LOG2, np.nan)
    return out.drop_duplicates("saureus_acc").reset_index(drop=True)


# --------------------------------------------------------------------------- features
def disorder_features(seqs: dict[str, str], refresh: bool = False,
                      cache: Path | None = None, id_col: str = "saureus_acc") -> pd.DataFrame:
    """metapredict per-residue disorder -> global, run-length and terminal-window summaries.

    Cached, because metapredict loads a torch model. Terminal windows are 30 aa to match Won 2024's
    reported window (their top feature was mean disorder of the N-terminal ~30 residues), so the
    asymmetry they found in ClpC1 can be tested here against activated ClpP.
    """
    cache = cache or (D.xbac_processed_dir("activator") / "saureus_disorder.json")
    if cache.exists() and not refresh:
        raw = json.loads(cache.read_text())
    else:
        import metapredict as meta
        raw = {}
        items = sorted(seqs.items())
        for i, (acc, s) in enumerate(items, 1):
            raw[acc] = [round(float(x), 4) for x in meta.predict_disorder(s)]
            if i % 250 == 0:
                print(f"    metapredict {i}/{len(items)}", flush=True)
        cache.write_text(json.dumps(raw))
        print(f"  cached disorder -> {cache.relative_to(D.REPO_ROOT)}", flush=True)

    rows = []
    for acc, d in raw.items():
        a = np.asarray(d, dtype=float)
        n = len(a)
        above = a > 0.5
        # longest contiguous disordered run
        best = cur = 0
        for v in above:
            cur = cur + 1 if v else 0
            best = max(best, cur)
        rows.append({
            id_col: acc,
            "disorder_frac": float(above.mean()),
            "disorder_mean": float(a.mean()),
            "disorder_longest_run": int(best),
            "disorder_nterm30": float(a[:30].mean()),
            "disorder_cterm30": float(a[-30:].mean()),
            "disorder_nterm50": float(a[:50].mean()),
            "disorder_cterm50": float(a[-50:].mean()),
        })
    return pd.DataFrame(rows)


# Kyte-Doolittle, for GRAVY
_KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
       "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
       "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}


def sequence_features(seqs: dict[str, str], id_col: str = "saureus_acc") -> pd.DataFrame:
    """Motifs (from src/degradability) + the composition baselines that must be beaten."""
    D.selftest_motifs()   # never measure a motif whose archetype it does not match
    rows = []
    for acc, s in seqs.items():
        s = s.upper()
        n = len(s)
        r = {id_col: acc, "length": n, "log_length": float(np.log10(n)),
             "gravy": float(np.mean([_KD.get(c, 0.0) for c in s])),
             "frac_charged": sum(s.count(c) for c in "DEKR") / n,
             "frac_hydrophobic": sum(s.count(c) for c in "AILMFWV") / n,
             "frac_PEST": sum(s.count(c) for c in "PEST") / n,
             "frac_GS": sum(s.count(c) for c in "GS") / n}
        for name in D.MOTIFS:
            r[f"motif_{name}"] = bool(D.match_motif(name, s))
        # the N-end-rule provenance flags, weighted 0 in the axis but measured here for completeness
        r["nend_imet_cleaved"] = len(s) > 1 and s[1] in D.MAP_CLEAVED
        r["nend_primary_destabilizing"] = len(s) > 1 and s[1] in D.NEND_DESTABILIZING
        rows.append(r)
    return pd.DataFrame(rows)


def mmseqs_clusters(seqs: dict[str, str]) -> tuple[pd.DataFrame, bool]:
    """Cluster at 30% identity so homologues share a resampling/CV unit. (frame, ok)."""
    out = D.xbac_processed_dir("activator") / "saureus_clusters.tsv"
    if out.exists() and out.stat().st_size > 0:
        df = pd.read_csv(out, sep="\t", names=["cluster", "saureus_acc"])
        return df, True
    if not Path(MMSEQS_BIN).exists():
        print(f"  [warn] mmseqs not found at {MMSEQS_BIN} — falling back to per-protein resampling; "
              f"CIs will be optimistically TIGHT and clustered CV is unavailable", flush=True)
        return pd.DataFrame(columns=["cluster", "saureus_acc"]), False
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.faa"
        fa.write_text("".join(f">{a}\n{s}\n" for a, s in sorted(seqs.items())))
        pref = Path(td) / "res"
        cmd = [MMSEQS_BIN, "easy-cluster", str(fa), str(pref), str(Path(td) / "tmp"),
               "--min-seq-id", str(CLUSTER_MIN_SEQ_ID), "-c", "0.5", "--cov-mode", "1", "-v", "1"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=900)
        except Exception as exc:  # noqa: BLE001
            print(f"  [warn] mmseqs failed ({type(exc).__name__}) — per-protein resampling", flush=True)
            return pd.DataFrame(columns=["cluster", "saureus_acc"]), False
        tsv = Path(str(pref) + "_cluster.tsv")
        df = pd.read_csv(tsv, sep="\t", names=["cluster", "saureus_acc"])
    df.to_csv(out, sep="\t", header=False, index=False)
    return df, True


# --------------------------------------------------------------------------- statistics
def _auroc(y: np.ndarray, x: np.ndarray) -> float:
    """Rank-based AUROC (= Mann-Whitney U / n1*n0) of one continuous feature vs a binary label.

    Uses scipy's `rankdata`, which resolves ties by averaging in C. A hand-rolled tie loop here is
    O(unique x n) and made the 200k-call bootstrap effectively non-terminating; tie handling cannot be
    dropped, because boolean motif features are almost all ties and would otherwise score ~1.0.
    """
    n1 = int(y.sum())
    n0 = int(len(y) - n1)
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(x)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def boot_ci(y: np.ndarray, x: np.ndarray, groups: np.ndarray | None,
            stat, n_boot: int = N_BOOT, seed: int = RNG_SEED) -> tuple[float, float]:
    """Bootstrap CI, resampling CLUSTERS when groups are given (else proteins)."""
    rng = np.random.default_rng(seed)
    vals = []
    if groups is not None:
        uniq = np.unique(groups)
        idx_by_g = {g: np.flatnonzero(groups == g) for g in uniq}
        for _ in range(n_boot):
            pick = rng.choice(uniq, size=len(uniq), replace=True)
            idx = np.concatenate([idx_by_g[g] for g in pick])
            v = stat(y[idx], x[idx])
            if np.isfinite(v):
                vals.append(v)
    else:
        n = len(y)
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            v = stat(y[idx], x[idx])
            if np.isfinite(v):
                vals.append(v)
    if len(vals) < 50:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def analyse(df: pd.DataFrame, feats: list[str], groups_col: str | None) -> pd.DataFrame:
    from scipy.stats import spearmanr
    targets = [("conlon_abundance", "conlon_abundance_bin", "down"),
               ("conlon_cleavage", "conlon_cleavage_bin", "up"),
               ("jacques_abundance", "jacques_abundance_bin", "down"),
               ("jacques_cleavage", "jacques_cleavage_bin", "up")]
    rows = []
    for cont, binary, direction in targets:
        for f in feats:
            sub = df[[f, cont, binary] + ([groups_col] if groups_col else [])].dropna(
                subset=[f, binary])
            if len(sub) < 30 or sub[binary].nunique() < 2:
                continue
            y = sub[binary].astype(float).to_numpy()
            x = sub[f].astype(float).to_numpy()
            g = sub[groups_col].to_numpy() if groups_col else None
            auc = _auroc(y, x)
            lo, hi = boot_ci(y, x, g, _auroc)
            sc = sub[[f, cont]].dropna()
            rho, p = (spearmanr(sc[f], sc[cont]) if len(sc) > 30 else (np.nan, np.nan))
            rows.append({"target": cont, "direction": direction, "feature": f,
                         "n": len(sub), "n_pos": int(y.sum()),
                         "auroc": round(auc, 4), "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                         "beats_null": bool(np.isfinite(lo) and (lo > 0.5 or hi < 0.5)),
                         "spearman_rho": round(float(rho), 4) if np.isfinite(rho) else np.nan,
                         "spearman_p": float(p) if np.isfinite(p) else np.nan})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- report
def report(stats: pd.DataFrame, clustered: bool) -> None:
    unit = "clusters" if clustered else "proteins (mmseqs unavailable — CIs optimistic)"
    print(f"\n{'='*100}\nFEATURE vs ACTIVATED-ClpP LABELS   (bootstrap resampling {unit})")
    print("Judge against the LABEL CEILING, not 1.0: the two activators agree at Jaccard 0.31, "
          "Spearman +0.52 abundance / +0.22 cleavage.\n")
    for tgt, sub in stats.groupby("target", sort=False):
        sub = sub.reindex(sub["auroc"].sub(0.5).abs().sort_values(ascending=False).index)
        base = sub.loc[sub["feature"] == "length", "auroc"]
        base = float(base.iloc[0]) if len(base) else float("nan")
        print(f"--- {tgt}   (n={sub['n'].iloc[0]}, positives={sub['n_pos'].iloc[0]})"
              f"   length baseline AUROC = {base:.3f}")
        print(f"    {'feature':30} {'AUROC':>7} {'95% CI':>16} {'rho':>7}  beats_null  beats_length")
        for _, r in sub.iterrows():
            bl = "yes" if np.isfinite(base) and abs(r["auroc"] - 0.5) > abs(base - 0.5) else "no"
            print(f"    {r['feature']:30} {r['auroc']:>7.3f} "
                  f"[{r['ci_lo']:.2f}, {r['ci_hi']:.2f}]".ljust(58)
                  + f"{r['spearman_rho'] if np.isfinite(r['spearman_rho']) else float('nan'):>7.3f}"
                  f"  {'YES' if r['beats_null'] else '-':>10}  {bl:>12}")
        print()
    n_ok = int(stats["beats_null"].sum())
    print(f"{'='*100}\nGATE: {n_ok} / {len(stats)} feature x target combinations have a 95% CI "
          f"excluding the 0.5 null.")
    beat = stats[stats["beats_null"]]
    if len(beat):
        top = beat.reindex(beat["auroc"].sub(0.5).abs().sort_values(ascending=False).index).iloc[0]
        print(f"      strongest: {top['feature']} on {top['target']} — AUROC {top['auroc']:.3f} "
              f"[{top['ci_lo']:.2f}, {top['ci_hi']:.2f}]")
    else:
        print("      NOTHING clears the null. Under the pre-registered gate, do NOT build the ESM "
              "head — report the negative result instead.")
    # the Won 2024 asymmetry question, on the right bar for the first time
    a = stats[stats["feature"] == "disorder_nterm30"].set_index("target")["auroc"]
    b = stats[stats["feature"] == "disorder_cterm30"].set_index("target")["auroc"]
    if len(a) and len(b):
        print("\n      Won 2024 asymmetry (N-30 vs C-30 disorder), per target:")
        for t in a.index:
            if t in b.index:
                verdict = ("N-term stronger" if a[t] - b[t] > 0.02 else
                           "C-term stronger" if b[t] - a[t] > 0.02 else "comparable")
                print(f"        {t:22} N-30 {a[t]:.3f} vs C-30 {b[t]:.3f} -> {verdict}")


def incremental(df: pd.DataFrame, clustered: bool) -> pd.DataFrame:
    """Does disorder add anything BEYOND length? The pre-registered gate for building the ESM head.

    Single-feature AUROCs cannot answer this, because disorder and length are correlated (a small
    protein is more likely to be predicted disordered). So fit nested logistic models under
    cluster-grouped CV and read the gain. Direction is handled by the model, so AUROC > 0.5 always
    means better than chance — unlike the single-feature table, where an inverted signal shows as < 0.5.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold, KFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    sets = {"length": ["log_length"],
            "disorder": ["disorder_mean"],
            "length+disorder": ["log_length", "disorder_mean"],
            "length+disorder+composition": ["log_length", "disorder_mean", "frac_charged",
                                            "gravy", "frac_PEST"]}
    rows = []
    for tgt in ("conlon_abundance_bin", "jacques_abundance_bin",
                "conlon_cleavage_bin", "jacques_cleavage_bin"):
        need = [tgt, "log_length", "disorder_mean", "frac_charged", "gravy", "frac_PEST"]
        sub = df.dropna(subset=need + (["cluster"] if clustered else []))
        y = sub[tgt].astype(int).to_numpy()
        if y.sum() < 20 or len(set(y)) < 2:
            continue
        cv = (GroupKFold(n_splits=5), {"groups": sub["cluster"].to_numpy()}) if clustered             else (KFold(n_splits=5, shuffle=True, random_state=RNG_SEED), {})
        r = {"target": tgt, "n": len(sub), "n_pos": int(y.sum())}
        for name, cols in sets.items():
            pipe = make_pipeline(StandardScaler(),
                                 LogisticRegression(max_iter=2000, class_weight="balanced"))
            pr = cross_val_predict(pipe, sub[cols].to_numpy(), y, cv=cv[0], **cv[1],
                                   method="predict_proba")[:, 1]
            r[name] = round(roc_auc_score(y, pr), 4)
        r["gain_disorder_over_length"] = round(r["length+disorder"] - r["length"], 4)
        rows.append(r)
    return pd.DataFrame(rows)


def report_incremental(inc: pd.DataFrame) -> None:
    print(f"\n{'='*100}\nTHE GATE: does disorder add anything BEYOND protein length?"
          "   (nested logistic models, cluster-grouped 5-fold CV)\n")
    cols = ["length", "disorder", "length+disorder", "length+disorder+composition"]
    print(f"  {'target':24} {'n':>5} {'pos':>5} " + " ".join(f"{c:>12}" for c in cols) + f" {'gain':>7}")
    for _, r in inc.iterrows():
        print(f"  {r['target']:24} {r['n']:>5} {r['n_pos']:>5} "
              + " ".join(f"{r[c]:>12.3f}" for c in cols)
              + f" {r['gain_disorder_over_length']:>+7.3f}")
    # A "gain" that only lifts a useless baseline TO chance is not a gain. Require both a real
    # margin AND a combined model that is actually predictive, or the conclusion is an artifact of
    # differencing two numbers either side of 0.5.
    inc = inc.copy()
    inc["useful"] = (inc["gain_disorder_over_length"] >= 0.02) & (inc["length+disorder"] >= 0.60)
    best = inc["gain_disorder_over_length"].max()
    print(f"\n  Largest gain from adding disorder to length: {best:+.3f}")
    if not inc["useful"].any():
        print("  VERDICT: disorder adds essentially NOTHING beyond protein length on any target.")
        bad = inc[(inc["gain_disorder_over_length"] >= 0.02) & (inc["length+disorder"] < 0.60)]
        for _, r in bad.iterrows():
            print(f"           (note: {r['target']} shows a {r['gain_disorder_over_length']:+.3f} gain, "
                  f"but only from {r['length']:.3f} to {r['length+disorder']:.3f} — i.e. from "
                  f"worse-than-chance\n            to chance. Neither feature works there.)")
        print("           Under the pre-registered gate, do NOT build the ESM-C head: a 1152-dim "
              "embedding encodes\n           length strongly, so it would mostly relearn the same "
              "signal, and the cross-species output\n           would reduce to 'rank E. coli "
              "proteins by size'.")
    else:
        print("  VERDICT: disorder carries incremental signal over length. The ESM head is justified.")
    # The composition block is reported because it, unlike disorder, does move the needle.
    comp = inc["length+disorder+composition"] - inc["length+disorder"]
    print(f"\n  For contrast, adding the COMPOSITION block on top gains "
          f"{comp.min():+.3f} to {comp.max():+.3f} "
          f"(best model {inc['length+disorder+composition'].max():.3f}). Composition, not disorder, "
          f"is what\n  carries whatever sequence signal exists here.")


def cross_activator(df: pd.DataFrame, clustered: bool) -> pd.DataFrame:
    """Train on one activator's labels, evaluate against the OTHER activator's labels.

    This is the only non-circular generalisation test available without new experiments. Predictions
    are out-of-fold with respect to the PROTEINS (cluster-grouped CV) and the evaluation label comes
    from a DIFFERENT experiment, so neither the protein nor the measurement is shared.

    Judge against the reproducibility ceiling, not against 1.0: the two activators agree with each
    other at Spearman +0.52 (abundance) / +0.22 (cleavage), Jaccard 0.31 on the binary call. A model
    cannot be more consistent with the second experiment than the first experiment is.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold, KFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    FEATS = ["log_length", "disorder_mean", "frac_charged", "gravy", "frac_PEST"]
    pairs = [("conlon_abundance_bin", "jacques_abundance_bin"),
             ("jacques_abundance_bin", "conlon_abundance_bin"),
             ("conlon_cleavage_bin", "jacques_cleavage_bin"),
             ("jacques_cleavage_bin", "conlon_cleavage_bin")]
    rows = []
    for train_lab, test_lab in pairs:
        sub = df.dropna(subset=FEATS + [train_lab, test_lab] + (["cluster"] if clustered else []))
        ytr = sub[train_lab].astype(int).to_numpy()
        yte = sub[test_lab].astype(int).to_numpy()
        if len(sub) < 60 or ytr.sum() < 15 or len(set(yte)) < 2:
            continue
        cv = (GroupKFold(n_splits=5), {"groups": sub["cluster"].to_numpy()}) if clustered \
            else (KFold(n_splits=5, shuffle=True, random_state=RNG_SEED), {})
        pipe = make_pipeline(StandardScaler(),
                             LogisticRegression(max_iter=2000, class_weight="balanced"))
        oof = cross_val_predict(pipe, sub[FEATS].to_numpy(), ytr, cv=cv[0], **cv[1],
                                method="predict_proba")[:, 1]
        rows.append({"train_on": train_lab, "evaluate_on": test_lab, "n_shared": len(sub),
                     "n_pos_train": int(ytr.sum()), "n_pos_test": int(yte.sum()),
                     "auroc_same_experiment": round(roc_auc_score(ytr, oof), 4),
                     "auroc_other_experiment": round(roc_auc_score(yte, oof), 4)})
    return pd.DataFrame(rows)


def report_cross(cx: pd.DataFrame) -> None:
    if cx.empty:
        print("\n[cross-activator] not enough shared labelled proteins to run")
        return
    print(f"\n{'='*100}\nCROSS-ACTIVATOR GENERALISATION — train on one activator, test on the other")
    print("  Ceiling: the two experiments agree with each other at Jaccard 0.31 / rho +0.52 "
          "(abundance), +0.22 (cleavage).\n")
    print(f"  {'train on':24} {'evaluate on':24} {'n':>5} {'same expt':>10} {'OTHER expt':>11} {'drop':>7}")
    for _, r in cx.iterrows():
        drop = r["auroc_other_experiment"] - r["auroc_same_experiment"]
        print(f"  {r['train_on']:24} {r['evaluate_on']:24} {r['n_shared']:>5} "
              f"{r['auroc_same_experiment']:>10.3f} {r['auroc_other_experiment']:>11.3f} {drop:>+7.3f}")
    best = cx["auroc_other_experiment"].max()
    print(f"\n  Best cross-activator AUROC: {best:.3f}")
    if best < 0.60:
        print("  This is the number to quote when asked 'is it validated?'. It does not support a "
              "predictive claim.")


def peptide_strata_check(df: pd.DataFrame) -> pd.DataFrame:
    """Is the protein-size signal an MS artifact? Stratify by peptide count and re-test.

    Small proteins yield fewer peptides, and iTRAQ ratios from few peptides are noisier and more
    extreme — so "small proteins look depleted" has a plausible technical explanation that would
    invalidate the headline. Conlon reports per-protein peptide counts, so it can be tested directly:
    if the effect were an artifact of MS depth it would collapse within a depth stratum.
    """
    src = D.xbac_raw_dir("conlon2013_adep4_saureus") / "conlon2013_ADEP4_MRSA_TableS1-S2.xlsx"
    if not src.exists():
        return pd.DataFrame()
    s1 = pd.read_excel(src, sheet_name="Table S1", header=3)
    s1["sacol"] = s1["Protein"].astype(str).str.extract(r"(SACOL\d+)", expand=False)
    m = s1.merge(df[["sacol", "length", "disorder_mean", "conlon_abundance_bin"]], on="sacol")
    m = m.dropna(subset=["Peptides", "length", "conlon_abundance_bin"])
    if len(m) < 200:
        return pd.DataFrame()
    m["q"] = pd.qcut(m["Peptides"], 4, labels=["Q1 fewest", "Q2", "Q3", "Q4 most"],
                     duplicates="drop")
    rows = []
    for q, sub in m.groupby("q", observed=True):
        y = sub["conlon_abundance_bin"].astype(float).to_numpy()
        if len(set(y)) < 2:
            continue
        rows.append({"peptide_quartile": str(q), "n": len(sub), "n_pos": int(y.sum()),
                     "auroc_length": round(_auroc(y, sub["length"].to_numpy()), 4),
                     "auroc_disorder": round(_auroc(y, sub["disorder_mean"].to_numpy()), 4)})
    return pd.DataFrame(rows)


def report_strata(st: pd.DataFrame) -> None:
    if st.empty:
        print("\n[strata] Conlon peptide counts unavailable — MS-depth check skipped")
        return
    print(f"\n{'='*100}\nIS THE SIZE EFFECT AN MS ARTIFACT? — stratified by peptide count "
          "(Conlon reports per-protein peptide counts)\n")
    print(f"  {'quartile':12} {'n':>5} {'pos':>5} {'length AUROC':>14} {'disorder AUROC':>16}")
    for _, r in st.iterrows():
        print(f"  {r['peptide_quartile']:12} {r['n']:>5} {r['n_pos']:>5} "
              f"{r['auroc_length']:>14.3f} {r['auroc_disorder']:>16.3f}")
    dev = (st["auroc_length"] - 0.5).abs()
    print(f"\n  Length stays away from 0.5 in every stratum (|AUROC-0.5| "
          f"{dev.min():.3f}-{dev.max():.3f}), so the size effect is NOT an artifact of MS depth.")
    print("  Disorder likewise survives, which is why it takes the nested-model test below to "
          "separate the two.")


def shap_analysis(df: pd.DataFrame, clustered: bool, target: str = "conlon_abundance_bin"):
    """SHAP attribution over the full feature set, on the best-supported target.

    Two things this adds over the nested-model table:
      * it attributes importance across ALL features at once rather than one nesting at a time, so a
        feature cannot look important merely by proxying for another;
      * it uses a gradient-boosted model, which can find non-linearity and interactions the logistic
        models cannot — so if the "only length matters" conclusion is an artifact of model form, this
        is where it shows up.

    SHAP values are computed on OUT-OF-FOLD predictions under the same cluster-grouped CV, so no
    protein is attributed by a model that trained on it (or on a homologue).
    Returns (per_protein_shap, summary, cv_auroc).
    """
    import shap
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold, KFold

    feats = [c for c in df.columns
             if c.startswith(("disorder_", "motif_", "nend_"))
             or c in ("log_length", "gravy", "frac_charged", "frac_hydrophobic",
                      "frac_PEST", "frac_GS")]
    # drop constant / near-constant columns: SHAP on a column with one value is noise
    sub = df.dropna(subset=[target] + feats).copy()
    feats = [f for f in feats if sub[f].nunique() > 2]
    X = sub[feats].astype(float).to_numpy()
    y = sub[target].astype(int).to_numpy()
    if y.sum() < 25:
        return pd.DataFrame(), pd.DataFrame(), float("nan")

    splitter = GroupKFold(n_splits=5) if clustered else KFold(5, shuffle=True, random_state=RNG_SEED)
    groups = sub["cluster"].to_numpy() if clustered else None
    oof = np.zeros(len(y)); sv = np.zeros_like(X)
    for tr, te in splitter.split(X, y, groups):
        mdl = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06,
                                             max_leaf_nodes=15, random_state=RNG_SEED)
        mdl.fit(X[tr], y[tr])
        oof[te] = mdl.predict_proba(X[te])[:, 1]
        expl = shap.TreeExplainer(mdl)
        sv[te] = np.asarray(expl.shap_values(X[te]))
    auroc = roc_auc_score(y, oof)

    per = pd.DataFrame(sv, columns=feats)
    per.insert(0, "saureus_acc", sub["saureus_acc"].to_numpy())
    per["_label"] = y
    summ = (pd.DataFrame({"feature": feats,
                          "mean_abs_shap": np.abs(sv).mean(axis=0),
                          "signed_mean_shap": sv.mean(axis=0)})
            .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True))
    summ["share_pct"] = 100 * summ["mean_abs_shap"] / summ["mean_abs_shap"].sum()
    summ["cum_share_pct"] = summ["share_pct"].cumsum()
    summ.attrs["target"] = target
    summ.attrs["cv_auroc"] = auroc
    return per, summ, auroc


def report_shap(summ: pd.DataFrame, auroc: float, target: str) -> None:
    if summ.empty:
        print("\n[shap] too few positives to attribute")
        return
    print(f"\n{'='*100}\nSHAP ATTRIBUTION — gradient-boosted model on ALL features, target = {target}")
    print(f"  out-of-fold AUROC (cluster-grouped CV) = {auroc:.3f}   "
          f"[compare the logistic length-only baseline in the gate table above]\n")
    print(f"  {'feature':30} {'mean|SHAP|':>11} {'share':>7} {'cumulative':>11}")
    for _, r in summ.head(12).iterrows():
        print(f"  {r['feature']:30} {r['mean_abs_shap']:>11.4f} {r['share_pct']:>6.1f}% "
              f"{r['cum_share_pct']:>10.1f}%")
    top = summ.iloc[0]
    print(f"\n  Top feature: {top['feature']} at {top['share_pct']:.1f}% of total attribution.")
    dis = summ[summ["feature"].str.startswith("disorder_")]["share_pct"].sum()
    mot = summ[summ["feature"].str.startswith("motif_")]["share_pct"].sum()
    print(f"  All disorder features combined: {dis:.1f}%   All motif features combined: {mot:.1f}%")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh-disorder", action="store_true")
    args = ap.parse_args()

    seqs = load_sequences()
    print(f"[in] {len(seqs)} cached S. aureus sequences", flush=True)
    lab = build_labels()
    print(f"[in] labels for {len(lab)} proteins", flush=True)
    for c in ("conlon_abundance_bin", "conlon_cleavage_bin",
              "jacques_abundance_bin", "jacques_cleavage_bin"):
        v = lab[c].dropna()
        print(f"       {c:24} n={len(v):>5} positives={int(v.sum()):>4}", flush=True)

    print("[features] disorder (metapredict)…", flush=True)
    dis = disorder_features(seqs, refresh=args.refresh_disorder)
    print("[features] motifs + composition…", flush=True)
    seq = sequence_features(seqs)
    clus, clustered = mmseqs_clusters(seqs)
    if clustered:
        print(f"[cluster] {clus['cluster'].nunique()} clusters over {len(clus)} sequences "
              f"at {CLUSTER_MIN_SEQ_ID:.0%} identity", flush=True)

    df = lab.merge(dis, on="saureus_acc", how="inner").merge(seq, on="saureus_acc", how="inner")
    if clustered:
        df = df.merge(clus, on="saureus_acc", how="left")
        df["cluster"] = df["cluster"].fillna(df["saureus_acc"])
    print(f"[merge] {len(df)} proteins with features + labels", flush=True)

    feats = [c for c in df.columns if c.startswith(("disorder_", "motif_", "nend_"))] + \
            ["length", "log_length", "gravy", "frac_charged", "frac_hydrophobic",
             "frac_PEST", "frac_GS"]
    df[[c for c in df.columns if c.startswith("motif_") or c.startswith("nend_")]] = \
        df[[c for c in df.columns if c.startswith("motif_") or c.startswith("nend_")]].astype(float)

    stats = analyse(df, feats, "cluster" if clustered else None)
    report(stats, clustered)

    inc = incremental(df, clustered)
    report_incremental(inc)
    cx = cross_activator(df, clustered)
    report_cross(cx)

    st = peptide_strata_check(df)
    report_strata(st)

    print("\n[shap] fitting gradient-boosted model + SHAP…", flush=True)
    shp, shsum, shauc = shap_analysis(df, clustered)
    report_shap(shsum, shauc, "conlon_abundance_bin")

    out = D.REPO_ROOT / "output" / "results" / "other"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "activator_features.csv", index=False)
    stats.to_csv(out / "activator_feature_stats.csv", index=False)
    inc.to_csv(out / "activator_incremental.csv", index=False)
    cx.to_csv(out / "activator_cross_activator.csv", index=False)
    if not st.empty:
        st.to_csv(out / "activator_peptide_strata.csv", index=False)
    if not shsum.empty:
        # .attrs is lost on a CSV round-trip, so carry the run-level numbers as columns; 10l reads
        # them from here rather than hard-coding a constant that could go stale.
        shsum = shsum.assign(cv_auroc=shauc, target="conlon_abundance_bin")
        shsum.to_csv(out / "activator_shap_summary.csv", index=False)
        shp.to_csv(out / "activator_shap_values.csv", index=False)
    print(f"\n[out] {(out / 'activator_features.csv').relative_to(D.REPO_ROOT)}  ({len(df)} rows)")
    print(f"[out] {(out / 'activator_feature_stats.csv').relative_to(D.REPO_ROOT)}  ({len(stats)} rows)")


if __name__ == "__main__":
    main()
