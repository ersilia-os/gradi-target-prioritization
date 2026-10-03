"""Load the structural ligandability axis -- can a small molecule bind this protein's fold?

    data/processed/pockets/pockets_<species>.tsv        THE DELIVERABLE (complete, canonical)
        uniprot_ac · p2rank_score · fpocket_score · p2rank_n_pockets · holo_identity ·
        pdb_n_structures · pdb_coverage · af_plddt

TWO KINDS OF COLUMN, NEVER MERGED
----------------------------------
* **Predicted, on the AlphaFold model**: `p2rank_score` (P2Rank 2.5.1 calibrated probability of
  the best pocket), `fpocket_score` (fpocket 4.0 druggability of the best pocket) and
  `p2rank_n_pockets` (P2Rank pockets at probability >= 0.5). Only pockets whose lining residues
  average pLDDT >= 70 count -- the single place model confidence enters.
* **Measured, in the PDB**: `holo_identity`, the % identity to the closest bacterial chain seen with
  a drug-like ligand bound in the aligned binding site (BioLiP, filtered by PLINDER's artefact
  list, PDBe cofactors, nucleotides and ECMDB metabolites -- see `scripts/pockets/holo.py`).
  And `pdb_n_structures` / `pdb_coverage`: how many PDB entries have a chain that IS this protein
  (>= 95% identity, ligand or not) and what fraction of its sequence they cover
  (`scripts/pockets/pdb_coverage.py`). The 95% rule counts near-identical proteins of other
  species, so a conserved enterobacterial protein inherits E. coli's structures.

A ZERO AND AN NA ARE DIFFERENT CLAIMS
---------------------------------------
With an AlphaFold model, 0 means the tools looked and found no admitted pocket. Without one the
pocket columns and `af_plddt` are NA. `holo_identity`, `pdb_n_structures` and `pdb_coverage` are
never NA: every protein was searched, and 0 means no drug-like holo structure in its bacterial family -- which on
these proteomes is ~88% of proteins, a fact about the PDB, not about the proteins.

HOW MUCH TO TRUST THE PREDICTED COLUMNS -- measured, see `scripts/pockets/merge.py`
-----------------------------------------------------------------------------------------
Predicted and measured columns share no input, so their agreement is a fair test. P2Rank separates
proteins with any drug-like holo evidence from those with none at AUROC ~0.70 on Kp; fpocket's
druggability manages ~0.58. **Prefer `p2rank_score`.** fpocket is kept because the owner asked for
both and it is a different method (geometry vs a trained model), not because it earns its place on
this test.

    from src import pockets as K
    df = K.load("kpneumoniae")
    K.druggability(df)        # convenience 0-1, derived on the fly, NOT stored
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
SPECIES = ("kpneumoniae", "ecoli", "saureus")
COLUMNS = ["uniprot_ac", "p2rank_score", "fpocket_score", "p2rank_n_pockets", "holo_identity",
           "pdb_n_structures", "pdb_coverage", "af_plddt"]


def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES} -- human is "
                         "not in this axis")


def _read(path: Path, hint: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run {hint} first")
    return pd.read_csv(path, sep="\t", low_memory=False)


# ---------------------------------------------------------------- the deliverable

def load(species: str) -> pd.DataFrame:
    """One species, 8 columns, canonical row order (asserted).

    **No `evidence` column** (owner's call, 2026-10-03). It was a function of two columns that ARE
    here, verified exactly reconstructible on all three species before removal:

        af_plddt is NA        <=>  no model, so the pocket columns could not be computed
        holo_identity > 0     <=>  a drug-like bacterial co-crystal was found for the family

    which is the whole of what `pdb+af` / `af_only` / `pdb_only` / `none` said. **An NA in the
    pocket columns is "could not look", not "looked and found nothing"** — a protein WITH a model
    and no admitted pocket gets 0. Never `fillna(0)`; that is the v1 mistake
    (`legacy/HISTORY.md:199`).

    The one thing the labels added and these columns cannot is the model's provenance — AlphaFold
    DB vs ESMFold for the proteins AFDB does not cover. That is `model_source` in
    `evidence/alphafold_<species>.tsv`, via `load_alphafold()`.
    """
    from src import matrices as M

    _check(species)
    d = _read(TASK_DIR / f"pockets_{species}.tsv", "scripts/pockets/merge.py")
    M.assert_canonical(d["uniprot_ac"], species, "structure table")
    d["p2rank_n_pockets"] = d["p2rank_n_pockets"].astype("Int64")
    return d[COLUMNS]


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    return pd.concat([load(s).assign(species=s) for s in species], ignore_index=True)


def druggability(df: pd.DataFrame) -> pd.Series:
    """Convenience 0-1: mean of the within-table percentile ranks of `p2rank_score` and
    `holo_identity`. **Derived on the fly and deliberately not stored** (the studiedness lesson:
    a stored blend becomes the number people quote and stops being checked).

    Why these two: they are the best-supported column of each kind -- P2Rank agrees with the PDB
    better than fpocket does (AUROC ~0.70 vs ~0.58). Percentile ranks rather than raw values,
    because the two are in different units and a raw average would weight by scale.
    Where there is no model, only `holo_identity` contributes. Call it on ONE species at a time:
    percentiles are relative to the table passed in.
    """
    p2 = df["p2rank_score"].rank(pct=True)
    ho = df["holo_identity"].rank(pct=True)
    return pd.concat([p2, ho], axis=1).mean(axis=1, skipna=True).rename("druggability")


# ---------------------------------------------------------------- evidence

def load_alphafold(species: str) -> pd.DataFrame:
    """Per protein: model status, length, mean pLDDT, ordered/disordered fractions."""
    _check(species)
    return _read(EVIDENCE_DIR / f"alphafold_{species}.tsv", "scripts/pockets/structures.py")


def load_pockets(species: str, admitted_only: bool = False) -> pd.DataFrame:
    """LONG: one row per pocket per tool, with lining residues and pocket pLDDT."""
    _check(species)
    d = _read(EVIDENCE_DIR / f"pocket_list_{species}.tsv", "scripts/pockets/predict.py")
    return d[d["admitted"]] if admitted_only else d


def load_pdb(species: str) -> pd.DataFrame:
    """Per protein: PDB entries / chains that are this protein, coverage, best identity, PDB ids."""
    _check(species)
    return _read(EVIDENCE_DIR / f"pdb_{species}.tsv", "scripts/pockets/pdb_coverage.py")


def load_holo(species: str) -> pd.DataFrame:
    """Per protein: best bacterial / any-organism holo chain, ligands, the QED-filtered variant."""
    _check(species)
    return _read(EVIDENCE_DIR / f"holo_{species}.tsv", "scripts/pockets/holo.py")


def load_ligand_classes() -> pd.DataFrame:
    """One row per BioLiP ligand code with every drug-likeness criterion as its own flag."""
    return _read(EVIDENCE_DIR / "ligand_classes.tsv", "scripts/pockets/holo.py")


def load_alphafill_comparison() -> pd.DataFrame:
    """AlphaFill drug-like transplants vs the BioLiP route, Kp and Ec (a measurement)."""
    return _read(EVIDENCE_DIR / "alphafill_comparison.tsv", "scripts/pockets/alphafill_check.py")


def druglike(codes) -> np.ndarray:
    """Boolean per chemical-component code under the shipped definition (unknown -> False)."""
    cl = load_ligand_classes().set_index("ligand")["druglike"]
    return np.array([bool(cl.get(str(c).upper(), False)) for c in codes])
