"""Load the structural ligandability axis -- can a small molecule bind this protein's fold?

    data/processed/pockets/pockets_<species>.tsv        THE DELIVERABLE (complete, canonical)
        uniprot_ac · p2rank_score · fpocket_score · n_ligands_pdb · n_ligands_alphafill ·
        n_pdb_structures · af_plddt

THREE KINDS OF COLUMN, NEVER MERGED
-------------------------------------
* **Predicted** on the AlphaFold model: `p2rank_score` (P2Rank 2.5.1 calibrated probability of the
  best pocket) and `fpocket_score` (fpocket 4.0 druggability of the best pocket). Only pockets
  whose lining residues average pLDDT >= 70 count -- the single place model confidence enters.
* **Measured** in the PDB: `n_ligands_pdb`, drug-like ligands seen bound to **this protein's own**
  structures, and `n_pdb_structures`, how many PDB entries are this protein at all (ligand or not;
  partial structures count).
* **Modelled** by a third party: `n_ligands_alphafill`, drug-like ligands AlphaFill transplanted
  onto this protein's AlphaFold model by structural superposition, from donors at ~30% median
  identity.

**NOT the same `n_ligands` as the ligands axis.** These two count distinct Bemis-Murcko
SCAFFOLDS physically seen in a structure (PDB) or superposed onto the model (AlphaFill), with no
potency involved, and they are DISJOINT. `ligands_<species>.tsv`'s `n_ligands_*` count distinct
MOLECULES with a measured pChEMBL >= 6 over a taxonomic pool, and those NEST. Structural evidence
vs assay evidence. **Never add the two ligand counts.** A co-crystal of this protein and a transplant from a remote
homolog are not the same evidence, and they differ by an order of magnitude in reach (tens of
proteins against ~1,500). There is deliberately no combined column.

**Both counts are NON-REDUNDANT**: distinct Bemis-Murcko generic scaffolds, the ChEMBL axis's own
definition (`scripts/ligands/chembl.py:_scaffold_chunk`, imported rather than restated). Counting
raw chemical-component codes overstates by 25-34%; those raw counts ship as `n_codes_*` in
`evidence/ligand_counts_<species>.tsv`.

A ZERO AND AN NA ARE DIFFERENT CLAIMS
---------------------------------------
With an AlphaFold model, a pocket score of 0 means the tools looked and found no admitted pocket.
Without one the pocket scores and `af_plddt` are NA. The ligand counts and `n_pdb_structures` are
never NA: every protein was searched, and 0 means nothing was found -- which for `n_ligands_pdb`
is ~95% of each proteome, a fact about what crystallographers have done, not about the proteins.
**Never `fillna(0)` the pocket scores**; that is the v1 mistake (`legacy/HISTORY.md:199`).

HOW MUCH TO TRUST THE PREDICTED COLUMNS -- measured, see `scripts/pockets/merge.py`
-----------------------------------------------------------------------------------------
Predicted and measured columns share no input, so their agreement is a fair test -- **but only
once protein LENGTH is controlled**, because length alone separates ligand-bearing proteins from
the rest at AUROC 0.65-0.67 (big proteins are crystallised more often and offer more surface).
Within length deciles, against proteins that really do carry a drug-like ligand in their own PDB
structures: **P2Rank 0.49 (Kp) / 0.56 (Ec) / 0.62 (Sa); fpocket 0.44 / 0.52 / 0.48.**
So the pocket scores are a **weak and inconsistent** guide to where ligands are actually found --
on K. pneumoniae, the anchor, P2Rank adds nothing over protein size. Treat them as a soft prior
and prefer the measured columns wherever they are non-zero.

**There is no `druggability()` helper.** One shipped briefly, averaging percentile ranks of
`p2rank_score` and the old `holo_identity`. It was removed on 2026-10-03 with that column: the
three kinds of evidence here are a weak prior, a sparse measurement and a third-party model, and
no defensible weighting of them exists. Combine at prioritisation time, across axes, where the
weighting is an explicit choice someone owns.

    from src import pockets as K
    df = K.load("kpneumoniae")
    K.load_ligands_pdb("kpneumoniae")          # which ligands, on which PDB chains
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
SPECIES = ("kpneumoniae", "ecoli", "saureus")
COLUMNS = ["uniprot_ac", "p2rank_score", "fpocket_score", "n_ligands_pdb", "n_ligands_alphafill",
           "n_pdb_structures", "af_plddt"]


def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES} -- human is "
                         "not in this axis")


def _read(path: Path, hint: str) -> pd.DataFrame:
    """**`NA` is sodium's PDB chemical-component code.** pandas reads it as a missing value by
    default, which silently nulled 12% of the transplant rows once; empty fields are still NaN."""
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run {hint} first")
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""], low_memory=False)


# ---------------------------------------------------------------- the deliverable

def load(species: str) -> pd.DataFrame:
    """One species, 7 columns, canonical row order (asserted).

    **No `evidence` column** (owner's call, 2026-10-03): it was a function of columns that ARE
    here -- `af_plddt` is NA exactly when no model exists, and `n_ligands_pdb > 0` exactly when a
    drug-like ligand was seen on this protein. The one thing the labels added and these columns
    cannot is the model's provenance, AlphaFold DB vs ESMFold; that is `model_source` in
    `evidence/alphafold_<species>.tsv`, via `load_alphafold()`.
    """
    from src import matrices as M

    _check(species)
    d = _read(TASK_DIR / f"pockets_{species}.tsv", "scripts/pockets/merge.py")
    M.assert_canonical(d["uniprot_ac"], species, "pockets table")
    return d[COLUMNS]


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    return pd.concat([load(s).assign(species=s) for s in species], ignore_index=True)


# ---------------------------------------------------------------- evidence

def load_alphafold(species: str) -> pd.DataFrame:
    """Per protein: model status, `model_source` (AlphaFold DB vs ESMFold), length, mean pLDDT."""
    _check(species)
    return _read(EVIDENCE_DIR / f"alphafold_{species}.tsv", "scripts/pockets/structures.py")


def load_pockets(species: str, admitted_only: bool = False) -> pd.DataFrame:
    """LONG: one row per pocket per tool, with lining residues and pocket pLDDT."""
    _check(species)
    d = _read(EVIDENCE_DIR / f"pocket_list_{species}.tsv", "scripts/pockets/predict.py")
    return d[d["admitted"]] if admitted_only else d


def load_pdb(species: str) -> pd.DataFrame:
    """Per protein: PDB entries that are this protein, coverage, best identity, PDB ids."""
    _check(species)
    return _read(EVIDENCE_DIR / f"pdb_{species}.tsv", "scripts/pockets/pdb_coverage.py")


def load_pdb_chains(species: str) -> pd.DataFrame:
    """LONG: every PDB chain that IS this protein (>= 95% identity), with its identity."""
    _check(species)
    return _read(EVIDENCE_DIR / f"pdb_chains_{species}.tsv", "scripts/pockets/pdb_coverage.py")


def load_ligands_pdb(species: str) -> pd.DataFrame:
    """LONG: drug-like ligand x this protein's own PDB chain (code, scaffold, pdb, chain)."""
    _check(species)
    return _read(EVIDENCE_DIR / f"ligands_pdb_{species}.tsv", "scripts/pockets/holo.py")


def load_ligands_alphafill(species: str) -> pd.DataFrame:
    """LONG: drug-like AlphaFill transplant (code, scaffold, donor pdb, donor identity,
    `local_rmsd`, clash count). **Unfiltered on RMSD** -- filter here for a stricter set."""
    _check(species)
    return _read(EVIDENCE_DIR / f"ligands_alphafill_{species}.tsv", "scripts/pockets/holo.py")


def load_transplants(species: str) -> pd.DataFrame:
    """LONG: EVERY AlphaFill transplant, before any drug-likeness filter."""
    _check(species)
    return _read(EVIDENCE_DIR / f"transplants_{species}.tsv", "scripts/pockets/alphafill.py")


def load_ligand_counts(species: str) -> pd.DataFrame:
    """Per protein: both scaffold counts, plus the raw code counts they collapse from."""
    _check(species)
    return _read(EVIDENCE_DIR / f"ligand_counts_{species}.tsv", "scripts/pockets/holo.py")


def load_ligand_classes() -> pd.DataFrame:
    """One row per ligand code with every drug-likeness criterion as its own flag, its
    Bemis-Murcko generic scaffold, and `in_biolip` (false for codes only AlphaFill uses)."""
    return _read(EVIDENCE_DIR / "ligand_classes.tsv", "scripts/pockets/holo.py")


def druglike(codes) -> list[bool]:
    """Per chemical-component code, under the shipped definition (unknown -> False)."""
    cl = load_ligand_classes().set_index("ligand")["druglike"]
    return [bool(cl.get(str(c).upper(), False)) for c in codes]
