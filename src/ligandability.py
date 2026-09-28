"""Load the stage-06 ligandability tables, and the pieces the stage and its plots share.

Stage 06 asks one question of every protein: **can a small molecule bind it?** Part 1 answers it
from measured bioactivity in ChEMBL 37. The anchor proteome is dark TrEMBL, so its accessions are
almost absent from ChEMBL; the axis therefore rests on sequence mapping plus homology transfer, and
every number here is "ligands measured against something that looks like this protein", never
"ligands measured against this protein" unless `direct_hit` is set.

On disk::

    data/processed/ligands/
        chembl_<species>.tsv        the deliverable, one row per protein, keyed on uniprot_ac
        evidence/ + scratch/
            chembl_targets.tsv      component_id x tid bridge: organism, superkingdom, target_type
            chembl_ligands.tsv      tid x parent_molregno: the best pchembl measured on that pair
            chembl_hits.tsv         the DIAMOND join: uniprot_ac x component_id, pident/qcov/scov
            chembl_targets.faa      the DIAMOND target database (one record per component_id)
            scaffolds.tsv           parent_molregno -> canonical smiles -> Murcko generic scaffold
            cutoff_sensitivity.tsv  every bucket count at pChEMBL >= 5 / 6 / 7
            control.tsv             v1's pre-registered numbers against v2's
            manifest.tsv

Why three evidence matrices and one dense table
------------------------------------------------
The expensive steps are the SQL pass over a 30 GB database and the DIAMOND search. Both are cached
as their own tables, and the buckets are computed *from* them. Re-banding the identity cuts, moving
the pChEMBL line or adding a bucket is therefore a recompute of `aggregate()`, not a re-run. The
dense per-protein table is where a **zero is readable**: a sparse hit list cannot distinguish "no
ligand" from "not searched".

Why the buckets are restricted to true Bacteria
------------------------------------------------
v1's first pass put no identity floor on a "non-human" bucket and reported **424** potent Kp
proteins, with rat `P97697` winning the bacterial slot for several of them. Restricting to true
Bacteria (`organism_class.l1`) plus a 40% identity floor corrected it to **175**. That restriction
is kept, and so that it is never silently applied the literal all-organism count ships beside it as
`allorg_n_compounds` / `allorg_n_scaffolds`, and human ships as its own liability block.

Counts are unions, not best-hit counts
---------------------------------------
v1 picked one ChEMBL target per bucket by `max(n_potent, best_pchembl)` and reported *its* counts,
which both under-counts a protein resembling several liganded targets and biases toward whichever
target happened to be screened hardest. Here a bucket's count is the union of distinct parent
compounds over every target in that bucket's pool. Provenance for the single closest target is kept
alongside as `best_target` / `best_pident` / `best_organism`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LIGAND_DIR = REPO_ROOT / "data" / "processed" / "ligands"
EVIDENCE_DIR = LIGAND_DIR / "evidence"
SCRATCH_DIR = LIGAND_DIR / "scratch"
# Human is deliberately NOT a query species here -- it is the reference pool the selectivity
# columns are measured against. Contrast src/orthology.py, where SPECIES is all four.
SPECIES = ("kpneumoniae", "ecoli", "saureus")

CHEMBL_VERSION = "37"

# ------- the chosen configuration

# Same-species matching is by ORGANISM NAME PREFIX, not tax_id. Measured: ChEMBL files strains
# under their own taxids, so tax_id 562 alone finds 65 E. coli single-protein targets while
# "Escherichia coli%" finds 225 -- K-12 lives at 83333. A two-word binomial prefix is used so that
# `Klebsiella aerogenes` is not swept in by a genus match.
SPECIES_ORGANISM = {
    "kpneumoniae": "Klebsiella pneumoniae",
    "ecoli": "Escherichia coli",
    "saureus": "Staphylococcus aureus",
}

# B and F both. F ("functional") carries most bacterial enzyme-inhibition data -- a gyrase
# supercoiling IC50 is F, not B -- so B-only would discard the majority of the antibacterial
# literature. The split is reported every run.
ASSAY_TYPES = ("B", "F")

# ChEMBL's confidence_score is target-type specific, and one gate cannot serve both tracks:
#   9 direct single protein   8 homologous single protein
#   7 direct complex subunits 6 homologous complex subunits
#   5/4 family, "multiple targets may be assigned"
# Measured on ChEMBL 37: >=8 removes 0 of 3,271,336 single-protein rows (it is subsumed by
# requiring pchembl_value at all), but >=8 returns EXACTLY ZERO protein complexes -- which would
# have shipped empty complex_* columns looking like a real biological zero.
MIN_CONFIDENCE_SINGLE = 8
MIN_CONFIDENCE_COMPLEX = 6

# PROTEIN FAMILY and PROTEIN COMPLEX GROUP sit at confidence 4-5, "multiple homologous protein
# targets may be assigned" -- the ligand's actual target is ambiguous even within the group, so
# they are measured and excluded rather than silently absent. See docs/ligands.md.
SINGLE_TYPES = ("SINGLE PROTEIN",)
COMPLEX_TYPES = ("PROTEIN COMPLEX",)
EXCLUDED_TYPES = ("PROTEIN COMPLEX GROUP", "PROTEIN FAMILY")

# pChEMBL 6 = 1 uM, the field-standard "genuine ligand" line and v1's choice, so the 175 Kp / 155 Ec
# controls stay comparable. 5 and 7 ship in cutoff_sensitivity.tsv. The cutoff is NOT tuned:
# counts are not comparable across cutoffs, so the sweep is reported, never optimised.
PCHEMBL_CUTOFFS = (5.0, 6.0, 7.0)
PCHEMBL_HEADLINE = 6.0

# Homology bands. 95/40 are CLAUDE.md house rules inherited from v1. 60 is this stage's own
# choice, set from stage 05's measured identities: Kp<->Ec RBH median is 86% and Kp<->human 37%,
# so 60 sits above the cross-kingdom noise floor and below the cross-genus ortholog median.
# It is the one arbitrary number in the design -- see docs/ligands.md.
DIRECT_PIDENT = 95.0
CLOSE_PIDENT = 60.0
REMOTE_PIDENT = 40.0

# Both coverages, not just the query's. v1 used --query-cover 50 alone, which admits a hit that
# covers half of our protein but a tenth of a much longer ChEMBL target.
MIN_QCOV = 50.0
MIN_SCOV = 50.0

# The distance ladder, tightest first. `direct` (>=95% identity) is "essentially this protein,
# possibly another strain"; `species` is the taxonomic cut; `close` and `remote` are transfer.
BUCKETS = ("direct", "species", "close", "remote")

# ------- v1's numbers, pre-registered

V1_KP_POTENT = 175  # Kp proteins with a potent bacterial ligand, after v1's rat-trap correction
V1_EC_POTENT = 155
V1_KP_POTENT_UNCORRECTED = 424  # the bug: no identity floor, "non-human" instead of Bacteria
V1_SPOT_KP_BLA = ("A0A0H3H184", "Q93LQ9", 100.0)  # HS11286 -> K. pneumoniae beta-lactamase

# Wide on purpose: v2 changes the bucket definition, adds quality gates and unions over a pool
# rather than taking one best hit, so equality would be the wrong test.
TOLERANCE_FRACTION = 0.50

DELIVERABLE_COLUMNS = (
    "uniprot_ac",
    "direct_hit",
    "best_target",
    "best_pident",
    "best_organism",
    "n_compounds_tested",
    "direct_n_compounds",
    "direct_n_scaffolds",
    "direct_best_pchembl",
    "species_n_compounds",
    "species_n_scaffolds",
    "species_best_pchembl",
    "close_n_compounds",
    "close_n_scaffolds",
    "close_best_pchembl",
    "remote_n_compounds",
    "remote_n_scaffolds",
    "remote_best_pchembl",
    "human_n_compounds",
    "human_n_scaffolds",
    "human_best_pchembl",
    "human_best_pident",
    "allorg_n_compounds",
    "allorg_n_scaffolds",
    "complex_n_compounds",
    "complex_best_target",
)

_BOOL = ("direct_hit",)
_INT = (
    "n_compounds_tested",
    "direct_n_compounds",
    "direct_n_scaffolds",
    "species_n_compounds",
    "species_n_scaffolds",
    "close_n_compounds",
    "close_n_scaffolds",
    "remote_n_compounds",
    "remote_n_scaffolds",
    "human_n_compounds",
    "human_n_scaffolds",
    "allorg_n_compounds",
    "allorg_n_scaffolds",
    "complex_n_compounds",
)
_FLOAT = (
    "best_pident",
    "direct_best_pchembl",
    "species_best_pchembl",
    "close_best_pchembl",
    "remote_best_pchembl",
    "human_best_pchembl",
    "human_best_pident",
)


# ------- loaders


def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")


def _path(directory: Path, name: str) -> Path:
    path = directory / name
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/ligands/chembl.py first")
    return path


# Identifier columns that MUST be integers. A string `component_id` here would make the DIAMOND
# join in the stage silently return zero rows rather than raising -- the worst kind of failure.
_ID_INT = ("component_id", "tid", "parent_molregno", "molregno")

# Bare column names in the evidence tables. They carry no `_`-suffix for the rules below to catch,
# and a string here compares fine against another string but raises against a float -- or, worse,
# joins to nothing. Named explicitly rather than pattern-matched.
_BARE_FLOAT = ("pchembl", "pident", "ppos", "qcov", "scov", "evalue", "bitscore")
_BARE_INT = ("alnlen", "qlen", "slen")


def _read(path: Path) -> pd.DataFrame:
    """Read a stage-06 TSV as strings, then coerce by column name."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in df.columns:
        if col in _BOOL:
            df[col] = df[col].map({"True": True, "False": False, "": False}).astype(bool)
        elif col in _ID_INT:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        elif (col in _INT or col in _BARE_INT or col.startswith("n_")
              or col.endswith(("_n_compounds", "_n_scaffolds"))):
            df[col] = pd.to_numeric(df[col].replace("", "0"), errors="coerce").astype("Int64")
        elif col in _FLOAT or col in _BARE_FLOAT or col.endswith(
            ("_pchembl", "_pident", "_qcov", "_scov", "_bitscore")
        ):
            df[col] = pd.to_numeric(df[col].replace("", None), errors="coerce")
    return df


def load(species: str) -> pd.DataFrame:
    """One species' ligandability evidence, one row per protein, keyed on uniprot_ac."""
    _check(species)
    return _read(_path(LIGAND_DIR, f"chembl_{species}.tsv"))


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column added back."""
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_targets() -> pd.DataFrame:
    """The ChEMBL side: component_id x tid, with organism, superkingdom and target_type."""
    return _read(_path(SCRATCH_DIR, "chembl_targets.tsv"))


def load_ligands() -> pd.DataFrame:
    """One row per (tid, parent_molregno): the best pChEMBL measured on that pair."""
    return _read(_path(SCRATCH_DIR, "chembl_ligands.tsv"))


def load_hits() -> pd.DataFrame:
    """The DIAMOND join: uniprot_ac x component_id with pident, qcov, scov, bitscore."""
    return _read(_path(EVIDENCE_DIR, "chembl_hits.tsv"))


def load_scaffolds() -> pd.DataFrame:
    """parent_molregno -> canonical_smiles -> Bemis-Murcko generic scaffold.

    An acyclic compound has no Murcko scaffold and gets the empty string. That is one shared
    bucket, not a missing value: it is never dropped and never counted as zero distinct scaffolds.
    """
    return _read(_path(EVIDENCE_DIR, "scaffolds.tsv"))


def load_cutoff_sensitivity() -> pd.DataFrame:
    """Bucket counts at every pChEMBL cutoff. Reported, never optimised."""
    return _read(_path(EVIDENCE_DIR, "cutoff_sensitivity.tsv"))


def control() -> pd.DataFrame:
    """v1's pre-registered numbers against v2's."""
    return _read(_path(EVIDENCE_DIR, "control.tsv"))


def manifest() -> pd.DataFrame:
    """One row per species: counts, coverage, provenance."""
    return _read(_path(EVIDENCE_DIR, "manifest.tsv"))


# ------- using the output


def evidence_level(df: pd.DataFrame) -> pd.Series:
    """The closest bucket in which a protein has a potent ligand.

    One of `direct`, `species`, `close`, `remote`, `none` -- ordered by how much of a leap the
    evidence requires. `direct` means a ChEMBL target at >=95% identity, i.e. the protein itself
    (possibly another strain); everything else is transfer from a homolog and should be read as a
    hypothesis, not a measurement.
    """
    out = pd.Series("none", index=df.index, dtype=object)
    for bucket in reversed(BUCKETS):  # remote first, so tighter buckets overwrite
        out[df[f"{bucket}_n_compounds"].fillna(0) > 0] = bucket
    return out


def selectivity_risk(df: pd.DataFrame) -> pd.Series:
    """Human ligand evidence weighed against bacterial -- higher means a worse selectivity story.

    `log2((human + 1) / (bacterial + 1))` on distinct scaffolds. Positive means the protein's human
    relatives are better liganded than its bacterial ones. This is a liability flag for triage, not
    a prediction of cross-reactivity: it says nothing about whether the binding sites are shared.
    """
    human = df["human_n_scaffolds"].fillna(0).astype(float)
    bact = df["remote_n_scaffolds"].fillna(0).astype(float)
    return pd.Series(np.log2((human + 1.0) / (bact + 1.0)), index=df.index)


def _coerce_precedents(d: pd.DataFrame) -> pd.DataFrame:
    """Force the numeric precedent columns to nullable numeric types.

    Defensive, and it earned its place: a round trip through TSV turns a column built with `pd.NA`
    into object dtype holding strings, so `best_pchembl_bacteria >= 6` raised TypeError and sorting
    ordered '9.02' above '10.1'. `src/precedents.py` now types them on the way out; this makes any
    file already on disk load correctly too.
    """
    for c in ("n_ligands_exact", "n_ligands_bacteria", "n_ligands_human",
              "n_targets_bacteria", "n_targets_human", "n_ligands_bacteria_complex",
              "n_compounds_assayed_bacteria", "n_compounds_potent_bacteria"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce").astype("Int64")
    for c in ("best_pident_bacteria", "best_pchembl_bacteria",
              "best_pident_human", "best_pchembl_human", "hit_rate_bacteria"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce").astype("Float64")
    for c in ("exact_route", "exact_target", "precedent_evidence"):
        if c in d.columns:
            d[c] = d[c].astype("string")
    return d


def load_effort(refresh: bool = False) -> pd.DataFrame:
    """Per-ChEMBL-target screening effort, from `scripts/ligands/effort.py`.

    One row per bacterial target: `n_compounds_assayed` (NO pChEMBL filter -- the denominator the
    axis never had), `n_compounds_assayed_any_relation`, `n_compounds_measurable`,
    `n_compounds_potent`, `n_compounds_reported_inactive`, `n_docs`, `n_assays`, `first_doc_year`,
    `last_doc_year`, `max_pchembl`, `hit_rate`.

    **`n_compounds_tested` in `chembl_<species>.tsv` is NOT this**, despite what an earlier version
    of `docs/ligands.md` said: it is computed from the already-pChEMBL-filtered ligand table, so it
    separates weakly-potent from never-measured rather than screened from unscreened.

    `hit_rate` is null, never 0, where nothing was assayed -- a zero would claim a measurement
    nobody made. Needs the 30.5 GB dump to regenerate; the extract itself is 869 rows.
    """
    path = SCRATCH_DIR / "chembl_effort.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing -- written by scripts/ligands/effort.py, which needs the ChEMBL dump "
            "restored (see data/raw/other/chembl/SOURCE.md).")
    return _read(path)


def load_precedents(species: str) -> pd.DataFrame:
    """Ligand precedent per protein: exact / bacterial / human counts.

    Written by `scripts/ligands/precedents.py --species`. Three counts and their provenance:

        n_ligands_exact      an EXACT match -- UniProt accession or identical sequence
        n_ligands_bacteria   UNIQUE molecules over every bacterial target passing the floors
        n_ligands_human      the same over human targets

    **`n_ligands_bacteria` counts MOLECULES, not target-compound pairs** -- distinct
    `parent_molregno` over the union, so a compound tested against three homologs counts once.

    **`n_ligands_human` is a LIABILITY, never add it to the bacterial count.** Measured on E. coli
    `clpP`: 136 bacterial against 210 human at 56.3% identity.

    **`n_ligands_bacteria_complex` is separate and not inside (b)** -- DNA gyrase is a
    `PROTEIN COMPLEX` in ChEMBL, and E. coli `gyrB` carries 666 single-protein ligands against
    1,412 complex ones. v1 dropped that track and made GyrA/GyrB look unliganded.

    These are POTENCY-MEASURABLE ligands only, so the axis does not say "has an antibiotic" --
    MIC and %-inhibition are absent by construction.
    """
    _check(species)
    return _coerce_precedents(_read(_path(LIGAND_DIR, f"precedents_{species}.tsv")))


def load_precedents_full(species: str) -> pd.DataFrame:
    """The same, with the provenance columns the deliverable omits.

    `precedents_<sp>.tsv` is deliberately four columns plus the key. This is the 13-column version:
    which route the exact match came through, how many targets each count unions over, best
    identity per side, and -- the one worth knowing about -- `n_ligands_bacteria_complex`.

    **The complex track is NOT in `n_ligands_bacteria`.** E. coli `gyrB` carries 666 single-protein
    ligands against 1,412 complex ones, and Kp `A0A0H3H0Y6` (gyrA) carries 1,410 complex
    against 131 single.
    Read the deliverable alone and DNA gyrase looks unliganded, which is the v1 error. Come here
    before concluding a target has no chemistry.
    """
    _check(species)
    return _coerce_precedents(_read(_path(EVIDENCE_DIR, f"precedents_full_{species}.tsv")))
