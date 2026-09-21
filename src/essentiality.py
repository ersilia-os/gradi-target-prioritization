"""Load the stage-07 essentiality artifacts.

Stage 07 has two halves. The **label corpus** (`essentiality/labels.py`,
`essentiality/deg_proteomes.py`) is measured essentiality from DEG: 173,048 labeled proteins over 38
species. The **predictors** score proteomes that have no measurement of their own — which is the
point, because *K. pneumoniae* HS11286, the anchor, has none.

Read these through the helpers rather than the raw TSVs: two of the columns are easy to misuse, and
both traps are enforced in the docstrings below.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ESSENTIALITY_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = ESSENTIALITY_DIR / "evidence"
SCRATCH_DIR = ESSENTIALITY_DIR / "scratch"
# Geptop is a PROKARYOTE method -- 37 prokaryote references, no eukaryotes -- so human is absent by
# construction, not by omission.
GEPTOP_SPECIES = ("kpneumoniae", "ecoli", "saureus")
GEPTOP_CUTOFF = 0.24          # the authors' default
# Our anchors that ARE among Geptop's own references, so their scores are circular.
GEPTOP_IN_REFERENCE_SET = ("ecoli", "saureus")


def _path(name: str) -> Path:
    path = ESSENTIALITY_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run the stage-07 scripts first")
    return path


def load_labels() -> pd.DataFrame:
    """The DEG label corpus: 173,048 rows, `essential` 0/1, with per-dataset provenance.

    `species` appears multiple times where DEG holds multiple screens (S. aureus 7x, P. aeruginosa
    4x). They are NOT consolidated: the base rate spans 15x across species and swings 117 -> 336 ->
    551 within P. aeruginosa alone, so union-vs-intersection is a choice to make explicitly.
    """
    return pd.read_csv(_path("evidence/labeled_proteins.tsv"), sep="\t")


def load_geptop(species: str) -> pd.DataFrame:
    """Geptop 2.0 scores for one proteome, keyed on `uniprot_ac`.

    **Two traps, both recorded as columns rather than left to the reader.**

    1. **A score of 0.0 has TWO meanings — read `geptop_evidence`, never the score alone.** On Kp,
       3,799 of 5,728 proteins (66.3%) score exactly zero, but they split:

         `essential_orthologs`       1,929 (33.7%) — essential RBH found, score > 0
         `orthologs_none_essential`  3,350 (58.5%) — has orthologs (median 7 references), none
                                     essential. **A CONFIDENT NON-ESSENTIAL CALL — evidence, not
                                     absence.**
         `no_orthologs`                449 ( 7.8%) — the method cannot see this protein at all.
                                     This, and only this, is missing information.

       So real coverage is **92.2%**, not 33.7%. An earlier version of this loader defined
       `geptop_informative` as `n_essential_rbh > 0`, which conflated the middle group with the last
       and understated coverage by 58 points. It is now `n_rbh > 0`.

       The zero-scoring rows are **tied**, so a zero means *unranked*, not "low rank" — never read a
       percentile off the zero block.
    2. **`geptop_score` is PROTEOME-RELATIVE.** Min-max runs within each species, so 0.6 in Kp and
       0.6 in Sa are not the same quantity. Use `geptop_score_raw` if you need something comparable
       across species — the same trap as stage 01's per-species t-SNE coordinates.

    And `geptop_in_reference_set` flags *E. coli* and *S. aureus*, whose own proteomes are two of
    Geptop's 37 references: their CV self-distance is 0, the method's `weight = 100` branch fires,
    and the score largely reads off the answer. Use the DEG measurement for those two.
    """
    if species not in GEPTOP_SPECIES:
        raise ValueError(f"{species!r} has no Geptop score; Geptop is prokaryote-only "
                         f"(expected one of {GEPTOP_SPECIES})")
    return pd.read_csv(_path(f"geptop_{species}.tsv"), sep="\t")


def load_geptop_all(species: tuple[str, ...] = GEPTOP_SPECIES) -> pd.DataFrame:
    frames = []
    for sp in species:
        df = load_geptop(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def geptop_reference_audit(species: str) -> pd.DataFrame:
    """Per-reference CV distance, weight, and RBH counts — the evidence behind each score.

    Read this before trusting a score: it shows which of the 37 references actually contributed,
    and the `weight` column makes the circularity visible (a self-reference weighs 100 against
    ~2.0–3.3 for everything else).
    """
    # Moved out of scratch/ in 2026-09: a scratch purge is meant to be free, and this table is
    # cited. See scripts/essentiality/geptop.py.
    return pd.read_csv(_path(f"evidence/geptop_reference_audit_{species}.tsv"), sep="\t")


def geptop_informative_only(species: str) -> pd.DataFrame:
    """Proteins Geptop had ANY orthology evidence for (`n_rbh > 0`), positive or negative.

    This is 92.2% of Kp, not 33.7% — it keeps the confident-negative group, which is a real call.
    Use `geptop_scored_only` if you want just the proteins with a non-zero score.
    """
    df = load_geptop(species)
    return df[df["geptop_informative"] == 1].copy()


def geptop_scored_only(species: str) -> pd.DataFrame:
    """Only proteins with a NON-ZERO score, i.e. at least one essential ortholog.

    Use this when a ranking must not be diluted by the tied zero block — but remember the excluded
    rows are mostly confident negatives, not unknowns.
    """
    df = load_geptop(species)
    return df[df["geptop_evidence"] == "essential_orthologs"].copy()


# ---------------------------------------------------------------- the merged table (the deliverable)

MERGE_SPECIES = ("kpneumoniae", "ecoli", "saureus")


def load(species: str) -> pd.DataFrame:
    """`essentiality_<species>.tsv` — the stage-07 deliverable, two headline columns per protein.

        geptop_ess   continuous 0–1, ALWAYS present. Geptop 2.0's orthology+phylogeny score.
        deg_ess      0 / 0.5 / 1 / NaN. Fraction of DEG screens ON THIS EXACT STRAIN that called the
                     protein essential. **NaN means unmeasured — which is all 5,728 of K. pneumoniae.**

    **`deg_ess` is three-valued only incidentally.** Each species happens to have exactly two
    screens today, so the fraction can only be 0, 0.5 or 1; a third screen would make it quarters.
    Read it as "fraction of screens agreeing", never as a fixed scale. The 0.5 bucket matters: on
    E. coli **490 of 695 essential calls rest on one screen only**, against 205 where both agree.

    Three things to hold when combining them:

    1. **`geptop_ess` is not a probability at the extremes.** Measured against out-of-set labels
       (*R. solanacearum*): scores of 0.70–1.00 are only ~70% essential, and a score of exactly 0 is
       still 3.1% essential. Monotonic and roughly calibrated mid-range (0.30–0.50 → 49.6%), so
       ranking is sound — but 1.0 is not certainty and 0.0 is not "no".
    2. **`geptop_ess` is proteome-relative** (min-max within each species). Use `geptop_score_raw`
       from `load_geptop()` for anything cross-species.
    3. **For E. coli and S. aureus `geptop_ess` is circular** — both are Geptop references,
       contributing 53.2% and 58.3% of their own score. For those two, `deg_ess` is the column that
       means something. Geptop exists for K. pneumoniae, which has no measurement at all.

    `essentiality` / `essentiality_source` merge the two (measurement where present, prediction
    otherwise) as a convenience. It **mixes units** — measured calls are pinned to 1.0/0.0 while
    predictions are continuous — so every measured essential outranks every prediction, and the
    column is NOT comparable across species. Prefer the two columns and your own rule.
    """
    if species not in MERGE_SPECIES:
        raise ValueError(f"{species!r} has no merged essentiality table "
                         f"(expected one of {MERGE_SPECIES})")
    return pd.read_csv(_path(f"essentiality_{species}.tsv"), sep="\t")


def load_all(species: tuple[str, ...] = MERGE_SPECIES) -> pd.DataFrame:
    """Every species in one frame, with a `species` column.

    Convenient, but **do not sort the result by `essentiality`**: E. coli and S. aureus are ~95%
    measured while K. pneumoniae is 100% predicted, so a pooled ranking favours the measured species
    for a reason that is not biology. Rank within a species.
    """
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------- the training sets, and the register

TRAINING_DIR = ESSENTIALITY_DIR / "training_sets"


def list_training_sets() -> list[str]:
    """Every clean training set on this axis, by name. One directory, no filename prefix to know."""
    if not TRAINING_DIR.is_dir():
        raise FileNotFoundError(
            f"{TRAINING_DIR} -- run scripts/essentiality/screens.py")
    return sorted(p.stem for p in TRAINING_DIR.glob("*.tsv"))


def load_training_set(name: str) -> pd.DataFrame:
    """One training set: `key`, `label`, `source_id`, `features_from` (+ screen-specific extras).

    **`key` IS NOT A UNIPROT ACCESSION** for the five screens measured on non-anchor strains. It is
    whatever identifier that strain's own proteome uses -- `CCN31837.1` (ECL8 EMBL protein id),
    `WP_038431262.1` (KPPR1 RefSeq), `KPNRH_00001` (a locus tag), `lcl|HG941718.1_prot_...` (a DEG
    FASTA header). `features_from` names the proteome the key indexes, and that proteome's
    embedding matrix is the only place the key means anything.

    Measured, so the trap is not theoretical: those five key onto the relevant anchor proteome at
    **0 of 4,930 / 4,809 / 4,981 / 4,981 / 5,433**. Joining one to `proteome_<species>.tsv` on
    accession returns an empty frame, not an error.
    """
    if name not in list_training_sets():
        raise FileNotFoundError(
            f"no training set {name!r}. Available: {', '.join(list_training_sets())}")
    return pd.read_csv(TRAINING_DIR / f"{name}.tsv", sep="\t")


def registry() -> pd.DataFrame:
    """Every essentiality dataset ever found, and whether it became a training set.

    `disposition` is the column to read: `training_set` / `held_back` / `refuted` /
    `corpus_only`. Most rows are NOT training sets, and `reason` carries why -- a positives-only
    gene list, a condition-dependent screen, an undeposited assembly, a DESeq run that did not
    converge. Generated by `scripts/essentiality/registry.py`, which exits non-zero unless it
    reconciles with `training_sets/` in both directions.
    """
    return pd.read_csv(_path("dataset_registry.tsv"), sep="\t")


def load_deg(species: str) -> pd.DataFrame:
    """Measured DEG essentiality for one anchor proteome: `deg_ess` + its provenance columns.

    `deg_ess` is the FRACTION of DEG screens on the exact anchor strain calling the protein
    essential, so with two screens it takes 0, 0.5 or 1 -- and **null means unmeasured**, never
    non-essential. K. pneumoniae is null throughout: DEG indexes no Klebsiella.
    """
    if species not in MERGE_SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {MERGE_SPECIES}")
    return pd.read_csv(_path(f"deg_{species}.tsv"), sep="\t")


def load_ogee(species: str) -> pd.DataFrame:
    """OGEE-trained essentiality PREDICTION for one anchor proteome.

    A third opinion beside `geptop_ess` and `deg_ess`, from a different machine: a learned model
    over 26 measured prokaryotic proteomes, against Geptop's reciprocal-best-hit ortholog score
    over 37 reference genomes.

    **Read `ogee_evidence`.** E. coli and S. aureus ARE in the OGEE corpus, so their values come
    from the leave-species-out fold that held them out (`leave_species_out`); K. pneumoniae is
    absent from OGEE entirely, so its values come from a full-corpus fit (`full_corpus`) and are a
    pure cross-species extrapolation.
    """
    if species not in MERGE_SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {MERGE_SPECIES}")
    return pd.read_csv(_path(f"ogee_{species}.tsv"), sep="\t")


def load_screens(species: str) -> pd.DataFrame:
    """One PREDICTED probability per published screen, for every protein of an anchor proteome.

    **These are predictions, not measurements, in every column and every row.** Each screen was
    measured on a strain that is not the anchor, so the column is a model trained on that strain
    and applied here. Even where the anchor protein happens to be in the screen's own training set
    the value is out-of-fold, so the whole column is on one comparable scale.

    `evidence/screens_transfer_audit.tsv` prices each column: the training strain, its measured
    own-organism grouped AUROC/AUPR, and the relevant cross-species figure. Read it before
    treating any column as evidence -- measured transfer is Ec->Kp 0.86 and Kp->Ec 0.93 mean
    AUROC, and it is better when the ASSAY matches.
    """
    if species not in MERGE_SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {MERGE_SPECIES}")
    return pd.read_csv(_path(f"screens_{species}.tsv"), sep="\t")
