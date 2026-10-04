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
           "n_pdb_structures", "af_plddt", "pockets_consensus", "pockets_evidence"]

# The columns the consensus ranks on -- the four that measure LIGANDABILITY.
#
# `af_plddt` is EXCLUDED: it is model quality, and including it would claim a confidently-modelled
# protein is more druggable. `n_pdb_structures` is EXCLUDED as the DENOMINATOR -- exactly what
# `n_assayed` is to `n_ligands` in the ligands axis: how often anyone crystallised this protein,
# not whether anything binds. In the consensus it would add a fame component that
# `studiedness/confounds.py` warns against double-counting; it does real work in the ladder.
# Same rule as degradability excluding `nn_similarity` (reach) and ligands excluding `n_assayed_*`.
CONSENSUS_COLUMNS = ("p2rank_score", "fpocket_score", "n_ligands_pdb", "n_ligands_alphafill")

# Spearman rho between `pockets_consensus` and protein length, measured 2026-10-04:
# Kp +0.663 · Ec +0.638 · Sa +0.654. The stage FAILS above this, which would mean the consensus had
# become more length-driven than `p2rank_score` itself (rho 0.72). It is a drift guard, NOT an
# endorsement -- see `consensus()`.
MAX_LENGTH_RHO = 0.75


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
    """One species, 9 columns, canonical row order (asserted).

    **The old four-label `evidence` column is gone** (owner's call, 2026-10-03), replaced since
    2026-10-04 by the standard `pockets_evidence` (1-3), which is a different quantity -- see
    `evidence()`. **The reconstructibility that justified removing it was checked against
    `holo_identity > 0`**, not against `n_ligands_pdb > 0`: `holo_identity` was itself deleted hours
    later the same day. Restating the claim in terms of the surviving column, as earlier versions of
    this docstring and of `docs/pockets.md` did, makes it a tautology rather than a verification.
    The one thing those labels added and these columns cannot is the model's provenance, AlphaFold
    DB vs ESMFold; that is `model_source` in `evidence/alphafold_<species>.tsv`, via
    `load_alphafold()`.
    """
    from src import matrices as M

    _check(species)
    d = _read(TASK_DIR / f"pockets_{species}.tsv", "scripts/pockets/merge.py")
    M.assert_canonical(d["uniprot_ac"], species, "pockets table")
    return d[COLUMNS]


def consensus(df: pd.DataFrame) -> pd.Series:
    """`pockets_consensus`: structural ligandability, 0-1, within-species.

    Mean percentile rank over `CONSENSUS_COLUMNS`; see that constant for what is left out and why.
    Plain `consensus.percentile_consensus` -- no zero-floor, unlike ligands: P2Rank reads exactly 0
    on 22.2% (Kp) / 19.4% (Ec) / 26.1% (Sa) of a proteome, large but nothing like the 97% that
    forced that deviation.

    **THIS COLUMN IS SUBSTANTIALLY A RANKING BY PROTEIN LENGTH. READ THIS BEFORE USING IT.**

        rho(p2rank_score, length)        Kp +0.723  Ec +0.693  Sa +0.736
        rho(THIS COLUMN, length)         Kp +0.663  Ec +0.638  Sa +0.654
        AUROC vs a measured PDB ligand -- p2rank_score   0.615 / 0.659 / 0.702
        AUROC vs a measured PDB ligand -- LENGTH ALONE   0.673 / 0.652 / 0.657
        p2rank_score WITHIN length deciles               0.494 / 0.561 / 0.619

    **On the anchor, protein length beats the pocket score, and within length deciles P2Rank sits
    at chance.** `p2rank_score` (rho 0.72) is nearly as length-confounded as `p2rank_n_pockets`
    (rho 0.84), which this axis DELETED for exactly that.

    **IT IS `druggability()` RETURNING, AND PRETENDING OTHERWISE WOULD BE DISHONEST.** That helper
    averaged percentile ranks of `p2rank_score` and the old `holo_identity` -- mechanically this
    function -- and was removed on 2026-10-03 because *"no defensible weighting exists across a
    weak prior, a sparse measurement and a third party's model"*, a decision named in CLAUDE.md's
    NO COMPOSITE SCORE, EVER section.

    **It ships on the project owner's instruction, 2026-10-04** -- *"do consensus based on all the
    columns. it may not be perfect, but it is something"* -- given AFTER the table above was
    measured and put to them. So this is a decision taken with the cost in view, not an oversight.

    The one thing that distinguishes the two: **`druggability()` was an intra-axis VERDICT, while
    `<axis>_consensus` exists for CROSS-AXIS COMPARABILITY** -- a reader stacking ten deliverables
    reads one scale without knowing any axis's internals. That justifies the column existing. It
    does not make the confound go away, which is why the confound is quoted here.

    **NEVER VALIDATE THIS AGAINST `n_ligands_pdb`**: that column is one of its inputs, so the 0.958
    AUROC such a check returns is circular and means nothing. The honest external comparison is the
    length row above.
    """
    from src import consensus as consensus_mod
    return consensus_mod.percentile_consensus(df, list(CONSENSUS_COLUMNS))


def evidence(df: pd.DataFrame) -> pd.Series:
    """`pockets_evidence`, the 1-3 ladder, as a nullable integer.

        3  MEASURED  -- this protein has its own PDB structure
        2  MODELLED  -- no structure, but AlphaFill transplanted a drug-like ligand
        1  PREDICTED -- pocket scores on a model, nothing else

    Counts: Kp 3,862 / 1,297 / 569 · Ec 2,013 / 497 / 1,893 · Sa 1,811 / 483 / 595.

    **"No model" is NOT a level -- only TWO proteins project-wide lack one** (Kp `irp1` 3,163 aa,
    Sa `ebh` 9,535 aa; both above AlphaFold DB's 2,700-aa cut and left unfolded on the owner's
    decision at ~18 h CPU each). A level built on two proteins would be degenerate.

    **AND THE LADDER DELIBERATELY DOES NOT DEPEND ON HAVING A MODEL** -- an experiment does not stop
    counting because AlphaFold declined the sequence. Sa `ebh` has **no model and 2 PDB structures**,
    so it is level 3, correctly: fragments of it have been crystallised even though the full-length
    fold was never predicted. Kp `irp1` has neither and is level 1. An assertion that model-less
    proteins must be level 1 was written during development and was WRONG; `ebh` is the case that
    caught it.

    Their `pockets_consensus` rests on the two COUNT columns alone, since
    `percentile_consensus` skips NaN row-wise -- two inputs instead of four, for two proteins.

    **It grades PROVENANCE, not outcome** -- the pattern `ligands_evidence` set, which is what keeps
    it independent of the consensus. Measured: within L2/L3 the two correlate NEGATIVELY (-0.37 Kp
    / -0.35 Ec / -0.43 Sa), because level 3 holds crystallised proteins that often have no drug-like
    ligand at all, while level 2 has a transplanted one by definition.

    **It is not the removed `evidence` column returning.** That one was four provenance labels
    verified reconstructible from `af_plddt` NA plus `holo_identity > 0` -- a column deleted hours
    later the same day. This ladder is built from `n_pdb_structures` and `n_ligands_alphafill`, and
    its purpose is cross-axis comparability.

    **`fillna(0)` here touches only the COUNT columns**, which have no NA by assertion upstream.
    Never `fillna(0)` a pocket score: an NA there means "could not look", and filling it was the v1
    mistake.
    """
    need = ("n_pdb_structures", "n_ligands_alphafill")
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise KeyError(f"pockets_evidence needs {missing}; present: {list(df.columns)}")
    num = lambda c: pd.to_numeric(df[c], errors="coerce").fillna(0)
    level = pd.Series(1, index=df.index, dtype="int64")
    level[num("n_ligands_alphafill") > 0] = 2
    level[num("n_pdb_structures") > 0] = 3
    return level


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
