"""Load the stage-04 degradability calls, and the pieces the stage shares.

Stage 04 asks one question of every bacterial protein: **is it a substrate of activated partnerless
ClpP?** -- the recruiting handle for BacPROTAC-style degradation. It is a binary classifier, trained
on the two *S. aureus* activator screens and applied to everything else:

    data/processed/degradability/degradability_{kpneumoniae,ecoli,saureus}.tsv
        uniprot_ac   adep4_prob   onc212_prob   nn_similarity

One number per activator, on one scale
--------------------------------------
`<act>_prob` is the model's probability for **every** protein, and for a labeled one it is the
**out-of-fold** value from models that did not train on it, averaged over `N_CV_SEEDS` repeats of
the cluster-grouped split. That is why it can be ranked across all 13,020 proteins at once. Pinning
a measured hit to 1.0 would have put all 233 of them above every predicted protein however
confident, so a shortlist would have filled with whatever happened to be assayed rather than with
the best candidates.

**The measured calls are NOT in this table** (owner's call, 2026-10-03). They were two columns that
were empty for 10,131 of 13,020 proteins -- every E. coli and K. pneumoniae row -- and the `_source`
columns beside them read `predicted` for all but 1,871. Nothing is lost: the measurement lives in
`evidence/labels_saureus.tsv` with its continuous log2FCs and its cluster assignment, which is
strictly more than the bit was, and **`measured()` / `with_measured()` below re-attach it** keyed on
`uniprot_ac`. `hits()` still lets a measurement win over the model, exactly as before -- it reads
the labels itself rather than a column.

Two activators, never merged
----------------------------
`adep4_*` from Conlon 2013, `onc212_*` from Jacques 2020. They are kept apart because they disagree
about where the line sits: their binary calls overlap at **Jaccard 0.32**, ONC212 calling 165 hits
ADEP4 misses against only 28 the other way -- it is the more promiscuous compound at 30 uM. They
nonetheless *rank* proteins similarly (AUROC 0.82-0.88 either way), which is why a probability-based
classifier is not damaged by that threshold disagreement. v1's `10c` max-pooled the two and its own
retrospective called that inflation.

The ceiling, and it is not 1.0
------------------------------
How well a *perfect* model of one activator predicts the other's calls, measured on 1,022 shared
proteins: **AUROC 0.877** for ADEP4, **0.815** for ONC212. Judge nothing here against 1.0.

The label is asymmetric, on purpose
-----------------------------------
ADEP4 is `log2FC <= -1 AND padj < 0.05`; ONC212 is `log2FC <= -1` alone, because Jacques published no
p-values at all. Each is the best its own paper supports, and both are v1's audited definitions
(`legacy/scripts/10e_activator_features.py:119-131`), which is what makes these numbers directly
comparable to v1's. Do not "fix" it into symmetry.

What these numbers are not
--------------------------
*E. coli* and *K. pneumoniae* have **no** activated-ClpP measurements, anywhere -- every row for them
is a **ranking hypothesis** from a model trained in another phylum, and and no screen measured
any of them. `nn_similarity` records how far the extrapolation reached; `evidence/domain_bands.tsv`
prices it. Read that file's caveats before quoting its numbers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DEGRADABILITY_DIR = REPO_ROOT / "data" / "processed" / "degradability"
EVIDENCE_DIR = DEGRADABILITY_DIR / "evidence"
SCRATCH_DIR = DEGRADABILITY_DIR / "scratch"
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# ---------------------------------------------------------------- the label source

# v1's audited merge of both papers' supplementary tables. Read this rather than the .xlsx: the parse
# is non-obvious (Conlon's header is on row 4 of S1 and row 6 of S2, `Sample_ A` carries a literal
# space, Jacques writes "NA" as a *string*) and it passed an internal-consistency audit -- recomputing
# Conlon's `Average` from its raw intensities agreed to max |delta| = 0.0010.
LABEL_TABLE = REPO_ROOT / "output" / "results" / "other" / "activator_features.csv"
LABEL_FASTA = (REPO_ROOT / "data" / "processed" / "legacy" / "v1" / "other_degradability"
               / "activator" / "saureus_activator_proteins.faa")
LABEL_CLUSTERS = LABEL_FASTA.parent / "saureus_clusters.tsv"

# activator -> the BINARY label column in LABEL_TABLE. NaN there means unmeasured, never negative.
ACTIVATORS: dict[str, str] = {
    "adep4": "conlon_abundance_bin",     # log2FC <= -1 AND padj < 0.05   (Conlon 2013, strain COL)
    "onc212": "jacques_abundance_bin",   # log2FC <= -1 only              (Jacques 2020, C0673)
}

# The continuous log2FCs behind those calls. Kept for the audit trail and the spot checks -- a
# measured value of -0.51 explains a `_hit` of 0 in a way the bit alone cannot.
ACTIVATOR_CONTINUOUS: dict[str, str] = {
    "adep4": "conlon_abundance",
    "onc212": "jacques_abundance",
}

ACTIVATOR_PAPER: dict[str, str] = {
    "adep4": "Conlon 2013 (ADEP4, MRSA, iTRAQ) Table S1; log2FC <= -1 and padj < 0.05",
    "onc212": "Jacques 2020 (ONC212) Table S3, 24 h; log2FC <= -1 (no p-values published)",
}

# v1's thresholds, for the label-integrity audit only -- the calls themselves are read from
# LABEL_TABLE rather than recomputed.
ABUNDANCE_LOG2_HIT = -1.0
PADJ_HIT = 0.05

# What a SECOND, INDEPENDENT EXPERIMENT achieves at the same task: the AUROC of one activator's
# continuous readout used to predict the OTHER's binary call, over the 1,022 proteins both screens
# measured. The natural yardstick for a predictor here -- far more informative than 1.0.
#
# Deliberately NOT called a "ceiling". It is not a bound: the model exceeds it at three of the five
# cutoffs in `evidence/cutoff_sensitivity.tsv`. And v1's `degradability_datasets.md` section 10.3
# already uses "reproducibility ceiling" for a DIFFERENT quantity -- label-vs-label agreement,
# rho 0.52 / Jaccard 0.32 -- so reusing the word here would conflate two measurements.
CROSS_ASSAY_AUROC = {"adep4": 0.877, "onc212": 0.815}

# Effect-size cutoffs for the sensitivity sweep. The shipped label uses -1.0 and that is not up for
# optimisation: **AUROC is not comparable across different labels.** A stricter cutoff shrinks the
# positive class and makes it more extreme, which inflates AUROC mechanically -- ADEP4 runs 0.827 at
# -0.5 up to 0.943 at -3.0 -- so tuning on it would buy a spurious +0.08 and a label chosen for being
# easy rather than meaningful.
#
# The sweep earns its keep a different way: CROSS_ASSAY_AUROC inflates in lockstep, so the *gap*
# between them is the invariant, and it is stable at every cutoff. That makes the conclusion
# cutoff-independent. See `load_cutoff_sensitivity()`.
CUTOFF_SWEEP: tuple[float, ...] = (-0.5, -1.0, -1.5, -2.0, -3.0)

# v1's benchmarks on these exact labels, for the run log to beat.
V1_LENGTH_ONLY_AUROC = {"adep4": 0.7747, "onc212": 0.6799}
V1_GBM_13_FEATURE_AUROC = {"adep4": 0.762}
V1_CROSS_ACTIVATOR_AUROC = {"adep4_to_onc212": 0.6916, "onc212_to_adep4": 0.7650}

# ---------------------------------------------------------------- sequence mapping

# Conlon is strain COL (SACOL / YP_18xxxx), Jacques is C0673 (ODV*), and the v2 proteome is
# NCTC 8325 (SAOUHSC_* / YP_498xxx). Measured: **0 of 1,943 labeled proteins join by ANY
# identifier** -- not RefSeq, not locus tag, not Jacques's own Mu50 UniProt column. Modern RefSeq
# re-tagged COL as SACOL_RS##### keeping no `old_locus_tag`, and UniProt demoted COL and Mu50 as
# redundant. So the join is by sequence, per CLAUDE.md, and within species it is nearly lossless:
# 1,873/1,943 (96.4%) at these cuts, median identity 100%.
MIN_PIDENT = 95.0        # CLAUDE.md's ">=95% identity = direct" rule
MIN_COVERAGE = 0.80
MIN_MAP_FRACTION = 0.90  # below this the stage exits; a silent loss here corrupts everything

# ---------------------------------------------------------------- model

# ---------------------------------------------------------------- the estimator
#
# **TabPFN-3.5 is the shipped estimator**; the RandomForest below is retained so the former shipped
# numbers stay regenerable with `--estimator forest`. The swap is evidence-led, not preference:
# measured over a 3-estimator x 3-embedding x 2-activator grid on byte-identical cluster-grouped
# folds (`output/results/degradability/proteomelm_vs_esmc.tsv`), TabPFN beat the forest on
# **PR-AUC in 6 of 6 arms, each winning all 5 seeds** (+0.041/+0.014/+0.056 on ADEP4 across
# esmc/proteomelm/prott5, +0.033/+0.028/+0.021 on ONC212). On AUROC the gain is smaller and only
# clears with ESM-C (+0.009 / +0.015).
#
# **The gain is concentrated at the TOP of the ranking** -- roughly 4x larger on PR-AUC than AUROC
# -- which is exactly what a prioritized shortlist consumes and what this stage's Fisher enrichments
# are computed on. Report both metrics; AUROC alone called 4 of those 6 arms "no difference".
#
# Use the PAIRED test when comparing arms: every head runs identical folds, so shared partition
# difficulty cancels. The unpaired 2-SD bar is badly under-powered here (on ONC212 both arms carry
# SD ~0.008 from the splits alone) and mislabels real effects as noise.
#
# lazy-qsar 3.4.4 was tested in the same grid and REJECTED: its one real gain (ADEP4 PR +0.021,
# 5/5 seeds) does not transfer to ONC212 (+0.005, 3/5), and it is the least seed-stable of the three.
DEFAULT_ESTIMATOR = "tabpfn"
TABPFN_MODEL = "v3.5"          # the default checkpoint in tabpfn 9.0.0 / tabpfn-client 0.6.0

N_SPLITS = 5
SEED = 0
N_BOOTSTRAP = 1000

# The cluster-grouped CV is repeated over this many seeds and the results averaged. Not decoration:
# a SINGLE-seed estimate carries about +/-0.004 of pure arbitrariness here, which is LARGER than
# either of the two effects that were tested and rejected (localization features +0.002, a 10-point
# hyperparameter sweep +0.003). Seed 0 alone reported 0.857 on ADEP4 against a 5-seed mean of 0.862.
N_CV_SEEDS = 5

# The forest. `class_weight="balanced"` carries the 13.7% / 24.6% positive rates without resampling.
#
# `max_features=0.1` and `min_samples_leaf=3` came out marginally ahead of the defaults (sqrt, 1) in
# a 10-point sweep and held their sign across 5 CV reseedings -- but by only **+0.002 / +0.003**, at
# or below seed noise. Adopted because it is free, NOT because tuning mattered: the whole sweep
# spanned 0.854-0.862 on ADEP4, less than one confidence-interval half-width, and its winner differed
# between the two activators. Do not spend more time here; the headroom is in better labels.
RF_PARAMS: dict = {
    "n_estimators": 500,
    "max_features": 0.1,
    "min_samples_leaf": 3,
    "class_weight": "balanced",
}

# One feature block: the ESM-C embedding. **Adding more features was tested and does not help** --
# localization (11 features from stage 03) reaches 0.748 ROC-AUC on its own but adds only +0.002 to
# ESM-C on ADEP4 and -0.005 on ONC212; COG function (24) adds a further +0.001. The embedding already
# encodes localization, which is unsurprising: DeepLocPro is ESM-2-based and TMbed is ProtT5-based.
# The full table is in docs/degradability.md. Localization is still worth joining in for
# *filtering* downstream -- see the compartment hit rates there -- just not as a feature.
#
# One feature block: the ESM-C embedding. The baseline ladder that used to run here (length,
# composition, esmc+length) is retired -- it was there to answer "does the embedding earn its keep
# over protein length alone", and it did, measured under this exact CV: **+0.081 on ADEP4
# (0.857 vs 0.776) and +0.062 on ONC212 (0.750 vs 0.687)**, against v1's length-only 0.7747 / 0.6799.
# Those numbers are recorded in docs/degradability.md rather than recomputed every run.
#
# One trap worth keeping from that exercise: a RANDOM FOREST is the wrong estimator for a one-feature
# baseline. On log_length alone it scored 0.682 against logistic's 0.776 on identical splits, because
# 500 bootstrapped trees bin a single variable into steps instead of using it monotonically. Quoting
# the forest's length number would have credited the embedding with the baseline's handicap.
FEATURE_BLOCKS = ("esmc",)

# Kyte-Doolittle hydropathy, for the composition block.
KYTE_DOOLITTLE: dict[str, float] = dict(zip(
    "ARNDCQEGHILKMFPSTWYV",
    [1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5,
     3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2],
))
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CHARGED = set("DEKRH")
HYDROPHOBIC = set("AILMFWVC")

# Bands for pricing the cross-species extrapolation. Same idea as the control that turned a rejected
# k-NN's 86.8% headline into an honest 68%: band performance by donor distance, then reweight by the
# distances that actually occur.
SIMILARITY_BANDS: tuple[tuple[float, float], ...] = (
    (0.00, 0.80), (0.80, 0.90), (0.90, 0.95), (0.95, 1.01),
)

OUT_COLUMNS = ["uniprot_ac", "adep4_prob", "onc212_prob", "nn_similarity"]


def sequence_features(sequences: pd.Series) -> pd.DataFrame:
    """Length and composition features from raw sequences.

    One function, used for the training species and the prediction species alike. `10f` learned this
    the hard way: scoring E. coli with a differently-scaled feature column silently corrupts every
    prediction, so the two sides must never have separate implementations.
    """
    seqs = sequences.fillna("").astype(str)
    n = seqs.str.len().clip(lower=1)
    out = pd.DataFrame(index=sequences.index)
    out["log_length"] = np.log10(n)
    for aa in AMINO_ACIDS:
        out[f"frac_{aa}"] = seqs.str.count(aa) / n
    out["gravy"] = [
        sum(KYTE_DOOLITTLE.get(c, 0.0) for c in s) / max(len(s), 1) for s in seqs
    ]
    out["frac_charged"] = [sum(c in CHARGED for c in s) / max(len(s), 1) for s in seqs]
    out["frac_hydrophobic"] = [sum(c in HYDROPHOBIC for c in s) / max(len(s), 1) for s in seqs]
    return out


def block_columns(block: str, feature_frame: pd.DataFrame) -> list[str]:
    """Which columns of the assembled feature frame a named block uses."""
    comp = [c for c in feature_frame.columns if c.startswith("frac_") or c == "gravy"]
    esmc = [c for c in feature_frame.columns if c.startswith("e") and c[1:].isdigit()]
    if block == "length":
        return ["log_length"]
    if block == "composition":
        return comp
    if block == "esmc":
        return esmc
    if block == "esmc+length":
        return esmc + ["log_length"]
    raise ValueError(f"unknown feature block {block!r}; expected one of {FEATURE_BLOCKS}")


def auroc_ci(y_true, y_score, n_boot: int = N_BOOTSTRAP, seed: int = SEED) -> dict:
    """ROC-AUC with a bootstrap 95% CI, stratified so class balance is held fixed.

    The CI is not optional. v1's `disorder_auroc()` made the same demand in its docstring and for the
    same reason: on this data a point estimate alone cannot tell a result from a rumour -- most of the
    features v1 tried had intervals covering 0.5.
    """
    from sklearn.metrics import roc_auc_score

    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    ok = np.isfinite(y_true) & np.isfinite(y_score)
    y_true, y_score = y_true[ok].astype(int), y_score[ok]
    pos, neg = np.flatnonzero(y_true == 1), np.flatnonzero(y_true == 0)
    if len(pos) < 5 or len(neg) < 5:
        return {"roc_auc": float("nan"), "roc_auc_lo": float("nan"),
                "roc_auc_hi": float("nan"), "n": int(len(y_true)), "n_pos": int(len(pos))}
    auc = float(roc_auc_score(y_true, y_score))
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(n_boot):
        # Resample within each class so the base rate cannot drift between replicates.
        idx = np.concatenate([rng.choice(pos, len(pos), replace=True),
                              rng.choice(neg, len(neg), replace=True)])
        boot.append(roc_auc_score(y_true[idx], y_score[idx]))
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return {"roc_auc": round(auc, 4), "roc_auc_lo": round(float(lo), 4),
            "roc_auc_hi": round(float(hi), 4), "n": int(len(y_true)), "n_pos": int(len(pos))}


def pr_auc(y_true, y_score) -> float:
    """Average precision. Necessary alongside ROC-AUC at these base rates.

    A ROC-AUC of 0.87 at a 13.6% positive rate can still mean poor precision in the top ranks, and
    the top ranks are where a shortlist actually lives. Compare it against the base rate, which is
    what a random ranker would score.
    """
    from sklearn.metrics import average_precision_score

    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    ok = np.isfinite(y_true) & np.isfinite(y_score)
    y_true, y_score = y_true[ok].astype(int), y_score[ok]
    if y_true.sum() < 5 or y_true.sum() == len(y_true):
        return float("nan")
    return round(float(average_precision_score(y_true, y_score)), 4)


# ---------------------------------------------------------------- loaders

def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")


def _check_activator(activator: str) -> None:
    if activator not in ACTIVATORS:
        raise ValueError(f"unknown activator {activator!r}; expected one of {tuple(ACTIVATORS)}")


def _path(directory: Path, name: str) -> Path:
    path = directory / name
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/degradability/predict.py first")
    return path


def load(species: str) -> pd.DataFrame:
    """One species' table: measured call, model probability and source per activator."""
    _check(species)
    return pd.read_csv(_path(DEGRADABILITY_DIR, f"degradability_{species}.tsv"), sep="\t")


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column added back."""
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_labels() -> pd.DataFrame:
    """The measured *S. aureus* calls and their log2FCs on `uniprot_ac`, with the join evidence."""
    return pd.read_csv(_path(EVIDENCE_DIR, "labels_saureus.tsv"), sep="\t")


def load_seqmap_audit() -> pd.DataFrame:
    """Every labeled protein's sequence join: source accession, target, pident, coverage, verdict."""
    return pd.read_csv(_path(EVIDENCE_DIR, "seqmap_audit.tsv"), sep="\t")


def load_cv(activator: str) -> pd.DataFrame:
    """Cross-validation table: feature block x split scheme x {roc_auc, pr_auc, base_rate}."""
    _check_activator(activator)
    return pd.read_csv(_path(EVIDENCE_DIR, f"cv_{activator}.tsv"), sep="\t")


def load_cross_activator() -> pd.DataFrame:
    """The non-circular test: fit on one activator, score the other's labeled set, both ways."""
    return pd.read_csv(_path(EVIDENCE_DIR, "cross_activator.tsv"), sep="\t")


def load_domain_bands() -> pd.DataFrame:
    """Out-of-fold AUROC banded by `nn_similarity`, and the per-species reweighted expectation."""
    return pd.read_csv(_path(EVIDENCE_DIR, "domain_bands.tsv"), sep="\t")


def load_oof(activator: str) -> pd.DataFrame:
    """Out-of-fold probabilities for every labeled protein, one column per feature block."""
    _check_activator(activator)
    return pd.read_csv(_path(EVIDENCE_DIR, f"oof_{activator}_saureus.tsv"), sep="\t")


def load_cutoff_sensitivity() -> pd.DataFrame:
    """Model vs cross-assay reference vs length, across log2FC cutoffs.

    The point of this table is the `gap` column, not the `model` column. `model` rises with cutoff
    strictness for a trivial reason -- a smaller, more extreme positive class is easier to separate --
    so it must never be read as "a stricter cutoff is better". `cross_assay` rises with it, and the
    gap between them stays put, which is what makes the shipped conclusion independent of the cutoff.
    """
    return pd.read_csv(_path(EVIDENCE_DIR, "cutoff_sensitivity.tsv"), sep="\t")


def manifest() -> pd.DataFrame:
    """One row per species: counts, coverage, provenance."""
    return pd.read_csv(_path(EVIDENCE_DIR, "manifest.tsv"), sep="\t")


# A balanced forest at a 13.7% / 24.6% base rate stays conservative: `adep4_prob` tops out around
# 0.61 on K. pneumoniae, so **0.5 is not a meaningful cut** -- it selects only 20 of 5,728 Kp
# proteins. These are the out-of-fold thresholds that reproduce each activator's own training base
# rate, which is a defensible default if you must threshold at all. Ranking is usually better.
# The probability cut that reproduces each activator's training base rate on the LABELED set --
# an empirical quantile of the out-of-fold distribution, so it is a property of the ESTIMATOR'S
# CALIBRATION and must be re-derived whenever the estimator changes. Measured on TabPFN-3.5's OOF;
# the forest's were 0.394 / 0.426 and are wrong for TabPFN, whose probabilities span a much wider
# range (max 0.87-0.95 on Kp against the forest's ~0.61, so `>= 0.5` now selects 220 Kp proteins
# where it selected 20).
#
# **It reproduces the base rate on the LABELED set, not on a proteome.** Applied whole-proteome the
# same cut selects 14-21% (adep4) and 31-35% (onc212), because unlabeled proteins do not share the
# labeled set's score distribution. That is expected, not a defect -- but do not read the cut as
# "this many proteins are substrates".
BASE_RATE_THRESHOLD = {"adep4": 0.328, "onc212": 0.313}


def hits(df: pd.DataFrame, activator: str, threshold: float | None = None) -> pd.DataFrame:
    """Predicted-or-measured positives for one activator.

    A measured call wins over the model where one exists; elsewhere the probability is thresholded.
    `threshold` defaults to `BASE_RATE_THRESHOLD[activator]` -- the value that reproduces the
    training base rate -- **not** 0.5, which at these base rates selects almost nothing.

    Thresholding is a decision you are making, not one the data supplies. Prefer ranking on `_prob`
    and taking a top-N when you can; the model was validated on ranking (AUROC), not on any cut.
    """
    _check_activator(activator)
    if threshold is None:
        threshold = BASE_RATE_THRESHOLD[activator]
    # The measured calls are no longer a column, so read them from the labels. Behaviour is
    # unchanged: a measurement still wins over the model wherever one exists. For E. coli and
    # K. pneumoniae the map is empty, so every row is thresholded -- which is the truth about them.
    hit = df["uniprot_ac"].map(measured(activator))
    prob = df[f"{activator}_prob"]
    return df[np.where(hit.notna(), hit == 1, prob >= threshold)]


def measured(activator: str) -> pd.Series:
    """The MEASURED 1/0 calls for one activator, indexed by `uniprot_ac`. 1,677 adep4 / 1,045
    onc212, all *S. aureus* -- there are no activated-ClpP measurements for E. coli or
    K. pneumoniae anywhere, so mapping this onto either returns all-NA, which is correct.

    This replaces the `<act>_hit` column the deliverable used to carry. The source of truth is
    `evidence/labels_saureus.tsv`, which also holds the continuous log2FC behind each call and the
    cluster used for grouping -- a measured -0.51 explains a 0 in a way the bit never could.
    """
    _check_activator(activator)
    # The labels file names its call columns exactly `adep4` / `onc212`.
    lab = load_labels().set_index("uniprot_ac")[activator]
    return lab[lab.notna()].round(0).astype("Int64")


def with_measured(df: pd.DataFrame) -> pd.DataFrame:
    """`df` plus an `<act>_hit` column per activator, for the plots and audits that need the
    measured call beside the probability. A copy -- the deliverable on disk stays four columns."""
    out = df.copy()
    for activator in ACTIVATORS:
        out[f"{activator}_hit"] = out["uniprot_ac"].map(measured(activator))
    return out
