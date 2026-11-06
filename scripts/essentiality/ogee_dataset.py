"""OGEE corpus -> (X, y, groups), and LEAVE-SPECIES-OUT cross-validation.

Turns `ogee_proteomes.py`'s labelled proteins into a matrix and scores it. The product is an
`ogee_ess` opinion that sits BESIDE `geptop_ess` and `deg_ess` -- a third, independent view, never
a replacement: Geptop is a reciprocal-best-hit ortholog score over 37 reference genomes, this is a
learned sequence model over 26 measured proteomes. Different machinery, different failure modes.

WHY LEAVE-SPECIES-OUT AND NOT THE PAPER'S SPLIT. Malbranke 2026 clusters proteins across all
genomes with MMseqs2 at 40% identity and assigns whole clusters to train/val/test, then holds out
S. cerevisiae and four E. coli strains once. Leave-species-out instead holds out EVERY taxon in
turn, and nothing is held out permanently -- the project owner's instruction. Two reasons it is the
better fit here:

  * it answers the question this axis actually asks. K. pneumoniae is absent from OGEE entirely, so
    the number we need is "given N measured bacteria, how well can we call the N+1th?" -- which is
    exactly what a held-out taxon measures, and exactly what `predict.py --score-on` measures for
    the published screens.
  * it gives 26 estimates instead of 1, so the SPREAD is visible. With base rates running 10.2x
    across taxa (0.048 to 0.493) a single held-out organism would be an accident of which one.

WHAT IT DOES AND DOES NOT CONTROL. No protein from the held-out taxon is in training. Its
ORTHOLOGS in the other 25 taxa are, and that is a real information channel -- `dnaA` appears ~26
times with ~26 labels. So these numbers answer "can we recognise a known family in a new organism",
not "can we predict a family never seen". That is the honest reading and it is the operationally
relevant one for Kp; the stricter question needs the paper's 40% clustering, which `--scheme
cluster` is left open for. Quote which one you ran.

Run with the `gradi` env:
    python scripts/essentiality/ogee_dataset.py --write-fasta      # then run esmc.py --strain
    python scripts/essentiality/ogee_dataset.py --embed-plan
    python scripts/essentiality/ogee_dataset.py                    # build + leave-species-out CV
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
from src import degradability as D  # noqa: E402
from src import matrices as M  # noqa: E402  -- RF_PARAMS, so both axes share one forest

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
FASTA_DIR = REPO_ROOT / "data" / "source" / "ncbi" / "ogee_proteomes"
STRAIN_EMB = REPO_ROOT / "data" / "processed" / "embeddings" / "scratch" / "strains"
TRAINING_DIR = OUT_DIR / "training_sets"
CORPUS = TRAINING_DIR / "ogee_corpus.tsv"

# ONE SEED, deliberately, and it is not the usual 5. The axis-wide 5-seed convention exists to
# average out FOLD-ASSIGNMENT noise -- but here the folds are fixed: each fold IS one taxon, so
# there is no partition randomness to average, only the forest's own, and per-seed SDs across this
# axis have run <= 0.004 throughout. Meanwhile the result being asked for is the SPREAD ACROSS 26
# TAXA, which one seed measures as well as five.
#
# It is also the difference between a usable run and an overnight one. Each fit is a 500-tree
# forest on ~76,000 x 1,152 (against ~4,000 x 1,024 for the per-screen CVs), measured at ~4 min;
# 26 folds x 5 seeds = 130 fits is 8-9 hours, 26 fits is ~1.7. If one taxon later looks
# surprising, reseeding that single fold costs 4 minutes.
N_SEEDS = 1
MIN_TEST_POS = 20      # a held-out taxon with fewer positives cannot support an AUPR worth quoting

VERBOSE = True


def say(m: str = "") -> None:
    if VERBOSE:
        print(m, flush=True)


def rule(c: str = "-", w: int = 108) -> None:
    say(c * w)


def load_corpus() -> pd.DataFrame:
    if not CORPUS.exists():
        sys.exit(f"FATAL {CORPUS} missing -- run scripts/essentiality/ogee_proteomes.py first")
    d = pd.read_csv(CORPUS, sep="\t")
    d["taxaID"] = d["taxaID"].astype(int)
    return d


def write_fastas(d: pd.DataFrame) -> None:
    """One FASTA per taxon, keyed on `key` (the assembly's protein id), for `esmc.py --strain`.

    ONE RECORD PER PROTEIN, not per label row: the corpus is already collapsed on
    (taxon, key), which is what fixed E. coli's two-namespace double count. Writing label
    rows would re-introduce it.
    """
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    say(f"  {'taxid':>8} {'organism':44s} {'proteins':>9}  file")
    for t, g in d.groupby("taxaID"):
        g = g.drop_duplicates("key")
        p = FASTA_DIR / f"{t}.faa"
        with p.open("w") as fh:
            for pid, seq in zip(g.key, g.sequence):
                fh.write(f">{pid}\n{seq}\n")
        say(f"  {t:>8} {str(g.organism.iloc[0])[:44]:44s} {len(g):>9,}  "
            f"{p.relative_to(REPO_ROOT)}")


def load_embeddings(taxa: list[int]) -> tuple[dict[int, dict[str, np.ndarray]], list[int]]:
    """key -> vector, per taxon. Missing taxa are RETURNED, not silently skipped."""
    have, missing = {}, []
    for t in taxa:
        f = STRAIN_EMB / f"embeddings_{t}.npz"
        if not f.exists():
            missing.append(t)
            continue
        z = np.load(f, allow_pickle=True)
        accs = [str(x) for x in z["accessions"]]
        mat = z["embeddings"]
        have[t] = {a: mat[i] for i, a in enumerate(accs)}
    return have, missing


def build(d: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    taxa = sorted(d.taxaID.unique())
    emb, missing = load_embeddings(taxa)
    if missing:
        say(f"  {len(missing)} taxa have no ESM-C yet: {missing[:8]}"
            f"{' ...' if len(missing) > 8 else ''}")
        say("  run:  python scripts/embeddings/esmc.py --strain " + " ".join(map(str, missing)))
        if not emb:
            sys.exit("FATAL no taxon has embeddings; nothing to build")
    rows, X = [], []
    for t in taxa:
        if t not in emb:
            continue
        lut = emb[t]
        g = d[d.taxaID == t]
        keep = g[g.key.astype(str).isin(lut)]
        # DROPPED, never zero-filled: a zero row is a real position in embedding space, and
        # inventing one puts a fabricated point in the training set.
        if len(keep) < len(g):
            say(f"    taxid {t}: {len(g) - len(keep)} of {len(g)} proteins have no vector, dropped")
        for pid, y, org in zip(keep.key.astype(str), keep.label, keep.organism):
            X.append(lut[pid])
            rows.append({"taxaID": t, "organism": org, "key": pid, "label": int(y)})
    meta = pd.DataFrame(rows)
    return (np.vstack(X).astype(np.float32), meta.label.to_numpy(dtype=int),
            meta.taxaID.to_numpy(dtype=int), meta)


def sd_txt(v: list[float]) -> str:
    """`+/-x` over several seeds, `(1 seed)` over one.

    NOT `+/-0.0000`: the SD of a single value is undefined, and numpy returns nan with a
    `Degrees of freedom <= 0` warning. Printing either would read as a measured zero variance.
    """
    import numpy as _np
    return f"+/-{_np.std(v, ddof=1):.4f}" if len(v) > 1 else "(1 seed)"


def forest(seed: int):
    from sklearn.ensemble import RandomForestClassifier
    # src.degradability.RF_PARAMS, imported not re-specified: n_estimators=500, max_features=0.1,
    # min_samples_leaf=3, class_weight="balanced". The class weight is load-bearing -- base rates
    # run 0.048 to 0.493 across these taxa.
    return RandomForestClassifier(n_jobs=-1, random_state=seed, **D.RF_PARAMS)


# Which OGEE taxon IS this anchor. `None` = the organism is absent from OGEE entirely.
ANCHOR_TAXID = {"ecoli": 83333, "saureus": 93061, "kpneumoniae": None}


def fit_corpus(X, y, seeds: int):
    """ONE model on the whole corpus, reused for every anchor.

    All 26 taxa, every time -- the project owner's instruction. Fitting per species would produce
    identical models and cost ~6 min each.
    """
    models = []
    for seed in range(seeds):
        m = forest(seed)
        m.fit(X, y)
        models.append(m)
    return models


def score_anchor(species: str, models, d: pd.DataFrame) -> pd.DataFrame:
    """OGEE essentiality for one anchor proteome -> `ogee_<species>.tsv`. Three columns.

        ogee_ess        the full-corpus model's probability. 0-1, never null.
        ogee_evidence   the MEASURED OGEE label for this exact protein: 1 essential, 0
                        non-essential, empty where OGEE has none.

    `ogee_evidence` is per-protein, not a per-table constant, and it is the honest way to carry
    the caveat that one model for all three anchors creates. E. coli K-12 and S. aureus NCTC 8325
    ARE OGEE taxa and ARE our anchors -- 100.0% and 97.5% of their corpus proteins are literally
    the same SEQUENCE as an anchor protein -- so for those rows the model was fitted on this
    protein with this label, and `ogee_ess` is closer to recall than to prediction. Where
    `ogee_evidence` is empty the value is a genuine out-of-corpus prediction. K. pneumoniae is
    absent from OGEE entirely, so its column is empty throughout.

    Read it exactly as `deg_ess` is read: **a null is UNMEASURED, never non-essential.**

    The labels are mapped onto the anchor BY EXACT SEQUENCE, the house join, and the rate is
    reported rather than assumed. Predictions run on the anchor's own ESM-C matrix -- same model,
    same 1,152 dims as the corpus -- so the output is complete by construction.
    """
    import numpy as np
    from src import embeddings as E
    from src import proteomes as P

    accs, X_anchor = E.vectors_for(species, list(M.canonical(species)))
    if len(accs) != len(M.canonical(species)):
        sys.exit(f"FATAL {species}: {len(M.canonical(species)) - len(accs)} proteins have no "
                 "ESM-C embedding, so the column could not be complete")
    col = np.mean([m.predict_proba(X_anchor)[:, 1] for m in models], axis=0)

    tx = ANCHOR_TAXID[species]
    label = pd.Series(pd.NA, index=range(len(accs)), dtype="Int64")
    if tx is not None:
        prot = P.load(species)[["uniprot_ac", "sequence"]]
        seq = dict(zip(prot.uniprot_ac, prot.sequence.astype(str)))
        lut = dict(zip(d[d.taxaID == tx].sequence.astype(str), d[d.taxaID == tx].label))
        label = pd.Series([lut.get(seq.get(a), pd.NA) for a in accs], dtype="Int64")
        n = int(label.notna().sum())
        say(f"    {n:,}/{len(accs):,} ({n / len(accs):.1%}) carry a measured OGEE label "
            f"(exact sequence); {int(label.sum() or 0):,} essential")
    else:
        say(f"    0/{len(accs):,} carry a measured OGEE label -- absent from OGEE entirely, so "
            "every value is an out-of-corpus prediction")
    say(f"    ogee_ess mean {col.mean():.4f}   top-decile cut {np.quantile(col, 0.9):.4f}")
    return pd.DataFrame({"uniprot_ac": accs, "ogee_ess": col.round(4),
                         "ogee_evidence": label})


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write-fasta", action="store_true", help="emit per-taxon FASTA and exit")
    ap.add_argument("--embed-plan", action="store_true", help="print the esmc command and exit")
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--score-anchor", nargs="+", metavar="SPECIES",
                    choices=["ecoli", "kpneumoniae", "saureus"],
                    help="fit on the corpus and write ogee_<species>.tsv. The species' own taxon "
                         "is held out where it is in the corpus (ecoli, saureus).")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    d = load_corpus()
    rule("=")
    say("ogee_dataset.py -- (X, y, groups) and leave-species-out CV for the ogee_ess corpus")
    rule("=")
    say(f"  corpus  {len(d):,} labelled proteins, {d.taxaID.nunique()} taxa, "
        f"{int(d.label.sum()):,} essential (base {d.label.mean():.3f})")
    say(f"  cv      leave-species-out x {a.seeds} seeds, RandomForest on ESM-C")
    say("  note    orthologs of the held-out taxon ARE in training; see the module docstring")
    rule()

    if a.write_fasta:
        say("PER-TAXON FASTA")
        rule()
        write_fastas(d)
        rule("=")
        say("next: python scripts/embeddings/esmc.py --strain "
            + " ".join(str(t) for t in sorted(d.taxaID.unique())))
        return
    if a.embed_plan:
        n = d.drop_duplicates(["taxaID", "key"]).shape[0]
        res = d.drop_duplicates(["taxaID", "key"]).sequence.str.len().sum()
        say(f"  {n:,} proteins, {res:,} residues -> ~{res / 1880 / 3600:.1f} h at 1,880 aa/s")
        return

    say("BUILD")
    rule()
    X, y, groups, meta = build(d)
    if a.score_anchor:
        rule()
        say("SCORE ANCHOR PROTEOMES")
        rule()
        say(f"  fitting {a.seeds} forest(s) on all {len(set(groups))} taxa, "
            f"{len(y):,} proteins ...")
        models = fit_corpus(X, y, a.seeds)
        for sp in a.score_anchor:
            say(f"  {sp}")
            t = score_anchor(sp, models, d)
            t = M.reindex(t, sp)
            o = OUT_DIR / f"ogee_{sp}.tsv"
            t.to_csv(o, sep="\t", index=False)
            say(f"    wrote {o.relative_to(REPO_ROOT)}  ({t.shape[0]} x {t.shape[1]})")
        rule("=")
        return
    say(f"  X {X.shape}   y {y.shape}   {len(set(groups))} groups   base {y.mean():.4f}")
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(SCRATCH_DIR / "ogee_dataset.npz", X=X, y=y, groups=groups,
                        key=meta.key.to_numpy(dtype=object))
    meta.to_csv(EVIDENCE_DIR / "ogee_dataset_rows.tsv", sep="\t", index=False)
    say(f"  wrote {(SCRATCH_DIR / 'ogee_dataset.npz').relative_to(REPO_ROOT)} and "
        f"{(EVIDENCE_DIR / 'ogee_dataset_rows.tsv').name}")

    from sklearn.metrics import average_precision_score, roc_auc_score
    rule()
    say("LEAVE-SPECIES-OUT")
    rule()
    say(f"  {'taxid':>8} {'organism':36s} {'n_test':>7} {'pos':>6} {'base':>6} "
        f"{'AUROC':>14} {'AUPR':>7} {'lift':>5}")
    recs, oof = [], np.full(len(y), np.nan)
    names = dict(zip(meta.taxaID, meta.organism))
    for t in sorted(set(groups)):
        te = groups == t
        tr = ~te
        if int(y[te].sum()) < MIN_TEST_POS or len(set(y[tr])) < 2:
            say(f"  {t:>8} {str(names[t])[:36]:36s} {te.sum():>7,} {int(y[te].sum()):>6} "
                f"  SKIPPED (<{MIN_TEST_POS} positives)")
            continue
        aucs, prs, ps = [], [], []
        for seed in range(a.seeds):
            m = forest(seed)
            m.fit(X[tr], y[tr])
            p = m.predict_proba(X[te])[:, 1]
            aucs.append(roc_auc_score(y[te], p))
            prs.append(average_precision_score(y[te], p))
            ps.append(p)
        oof[te] = np.mean(ps, axis=0)
        base = float(y[te].mean())
        say(f"  {t:>8} {str(names[t])[:36]:36s} {te.sum():>7,} {int(y[te].sum()):>6} "
            f"{base:>6.3f} {np.mean(aucs):>8.4f} {sd_txt(aucs):>10s} "
            f"{np.mean(prs):>7.4f} {np.mean(prs) / base:>5.2f}")
        recs.append({"taxaID": t, "organism": names[t], "n_test": int(te.sum()),
                     "n_pos": int(y[te].sum()), "base_rate": round(base, 4),
                     "roc_auc": round(float(np.mean(aucs)), 4),
                     "roc_auc_sd": (round(float(np.std(aucs, ddof=1)), 4)
                                    if len(aucs) > 1 else pd.NA),
                     "pr_auc": round(float(np.mean(prs)), 4),
                     "pr_lift": round(float(np.mean(prs)) / base, 2),
                     "n_train": int(tr.sum()), "n_seeds": a.seeds})
    if not recs:
        sys.exit("FATAL no taxon had enough positives to score")
    r = pd.DataFrame(recs)
    r["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    r.to_csv(EVIDENCE_DIR / "ogee_leave_species_out.tsv", sep="\t", index=False)

    rule()
    say("SUMMARY  -- the SPREAD is the result, not the mean")
    rule()
    say(f"  taxa scored {len(r)} of {d.taxaID.nunique()}")
    say(f"  AUROC  mean {r.roc_auc.mean():.4f}  median {r.roc_auc.median():.4f}  "
        f"min {r.roc_auc.min():.4f} ({r.loc[r.roc_auc.idxmin(), 'organism'][:30]})  "
        f"max {r.roc_auc.max():.4f}")
    say(f"  AUPR   mean {r.pr_auc.mean():.4f}  median {r.pr_auc.median():.4f}   "
        f"lift over base: median {r.pr_lift.median():.2f}x")
    pooled = np.isfinite(oof)
    say(f"  pooled over all held-out taxa: AUROC {roc_auc_score(y[pooled], oof[pooled]):.4f}  "
        f"AUPR {average_precision_score(y[pooled], oof[pooled]):.4f}")
    say("  READ THE PER-TAXON COLUMN, NOT THE POOLED NUMBER: pooling mixes organisms whose base "
        "rates differ 10x, which flatters the metric.")
    say(f"  wrote {(EVIDENCE_DIR / 'ogee_leave_species_out.tsv').relative_to(REPO_ROOT)}")
    rule("=")


if __name__ == "__main__":
    main()
