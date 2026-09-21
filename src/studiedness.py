"""Load the studiedness axis -- how much is already known about each protein.

The axis is wanted in both directions: an uncharacterised target is a risk, but it is also the
novelty the GraDi collaboration is looking for. `novelty()` is the complement of the family score.

    data/processed/studiedness/studiedness_<species>.tsv        THE DELIVERABLE
        uniprot_ac · studiedness_own · studiedness_family · evidence

THE TWO COLUMNS ANSWER DIFFERENT QUESTIONS AND MUST NOT BE COLLAPSED
---------------------------------------------------------------------
`studiedness_own` is what is known about **this accession**. `studiedness_family` is what is known
about the **best-characterised homolog** reachable in SwissProt. Keeping both is what makes "dark
in K. pneumoniae, famous in E. coli" readable off one row -- and that gap is the normal case here,
not an edge case.

**`studiedness_own` is near-flat on Kp and Sa by design. That is the measurement, not a defect.**
Measured on UniProt 2026_03: K. pneumoniae HS11286 is 0.1% reviewed with 5,710 of 5,728 proteins
carrying exactly one PubMed id (the genome paper), and S. aureus NCTC 8325 has none for 2,532 of
2,889. Do not rank Kp or Sa on `studiedness_own`; rank on `studiedness_family`.

BOTH COLUMNS ARE ON ONE FIXED GLOBAL SCALE
-------------------------------------------
Each is `0.6 * min(1, log1p(n_pubs) / log1p(P_REF)) + 0.4 * (annotation_score - 1) / 4`, where
`P_REF` is a **single constant** over all 575,748 reviewed UniProt entries, recorded in
`evidence/scale.tsv`. It is deliberately NOT a per-proteome percentile: that is the documented
`geptop_score` trap, where a proteome-relative score means 0.6 in Kp is not 0.6 in Sa. Here the
three species, and the two columns, are directly comparable. Raw counts ship in `evidence/` so the
weights can be changed without re-running anything.

**`P_REF` is the 99th percentile (204 papers), not the 95th, and that was measured.** The SwissProt
publication distribution is extremely skewed -- P50 1, P75 3, P90 12, P95 38, P99 204, max 20,402 --
so anchoring at P95 clipped the top of the scale and left E. coli with a median `studiedness_own`
of exactly 1.000: unrankable ties over half a proteome, the same failure CLAUDE.md records for v1's
under-regularised logistic in stage 04. At P99 only ~1% of SwissProt saturates. A
`saturated > 20%` guard exits `transfer.py` if that ever stops being true.

`n_pubs` is the union of UniProt's `lit_pubmed_id` and NCBI's gene2pubmed -- measured, because
gene2pubmed resolves a well-studied organism ~4x more finely (E. coli 193 distinct values against
UniProt's 48) and gives 184 S. aureus proteins their first paper.

READ `evidence` BEFORE TREATING A LOW SCORE AS NOVELTY
-------------------------------------------------------
Nested identity bands, the shape the `ligands/` axis already uses, plus two tiers that both
score 0 and mean different things:

    swissprot_direct    >= 95% identity -- effectively the same protein
    swissprot_close     >= 60% identity
    swissprot_homolog   >= 40% identity -- the transfer band
    below_floor         a SwissProt hit exists, but under 40% identity or 50% coverage
    no_hit              NOTHING in SwissProt resembles this protein at all

The last two both score **0.0**, and the axis deliberately refuses to invent a number for either
-- but they are very different claims, so they are separate tiers. `no_hit` is the strongest
novelty signal the axis produces: nothing among 575,748 curated entries looks like this protein.
`below_floor` means it does have a distant relative whose literature is simply too far away to
carry. Without the split, a third of K. pneumoniae would be one undifferentiated tie at the
bottom of the ranking.

Neither is a missing measurement and neither may be imputed. The hard ceiling is DIAMOND's own
hit rate: *any* SwissProt hit exists for only **81.6% of Kp and 73.0% of Sa** (99.4% of Ec).

THE 40% FLOOR IS A MEASURED TRADE-OFF, AND 25% WAS TRIED AND REJECTED
----------------------------------------------------------------------
**Decoys do not decide it.** Composition-preserving shuffles of our own 13,020 sequences match
SwissProt at **0.0% at every floor down to 20%** (`evidence/decoy_calibration.tsv`), so the floor
is not defending against spurious homology -- DIAMOND's e-value already does that.

**The held-out control decides it, and it is a monotonic trade-off** between coverage and
fidelity (`evidence/floor_sensitivity.tsv`):

    floor    Kp cov   Sa cov   held-out control rho
     25%      77.6%    69.1%         0.4248
     40%      66.0%    56.0%         0.5338     <-- shipped
     60%      56.6%    40.5%         0.5487

25% buys ~12 points of K. pneumoniae coverage for 0.11 of control. 40% is CLAUDE.md's standing
annotation-transfer threshold and nothing measured here justifies departing from it. The whole
sweep ships, so the floor can be re-chosen without re-running DIAMOND.

DO NOT MAXIMISE THE CONTROL -- IT HAS A CEILING WELL BELOW 1
--------------------------------------------------------------
The control correlates the transferred family score against E. coli's OWN literature, and those
are deliberately different quantities: a protein with 3 papers of its own whose human homolog has
300 *should* score low on `own` and high on `family`. That gap is the entire point of the axis,
so perfect agreement is not the target and rho cannot reach 1.

Worse, rho rises monotonically with the floor (0.59 at 95%) for a nearly circular reason: the
control excludes *Escherichia* only, so at a high floor the surviving donors are largely
Salmonella and Shigella near-duplicates whose publication counts track E. coli's because they are
effectively the same proteins. **Tuning the floor to maximise rho would drive it to 95%, where
`family` collapses onto `own` and the axis stops doing anything.** It is a pass/fail check that
transfer carries real signal on held-out data, not a score to optimise.

`n_pubs` is the union of UniProt's `lit_pubmed_id` and NCBI's gene2pubmed -- measured, because
gene2pubmed resolves a well-studied organism ~4x more finely (E. coli 193 distinct values against
UniProt's 48) and gives 184 S. aureus proteins their first paper.

READ `evidence` BEFORE TREATING A LOW SCORE AS NOVELTY
-------------------------------------------------------
Nested identity bands, the shape the `ligands/` axis already uses, plus two tiers that both
score 0 and mean different things:

    swissprot_direct    >= 95% identity -- effectively the same protein
    swissprot_close     >= 60% identity
    swissprot_homolog   >= 40% identity -- the transfer band
    below_floor         a SwissProt hit exists, but under 40% identity or 50% coverage
    no_hit              NOTHING in SwissProt resembles this protein at all

The last two both score **0.0**, and the axis deliberately refuses to invent a number for either
-- but they are very different claims, so they are separate tiers. `no_hit` is the strongest
novelty signal the axis produces: nothing among 575,748 curated entries looks like this protein.
`below_floor` means it does have a distant relative whose literature is simply too far away to
carry. Without the split, a third of K. pneumoniae would be one undifferentiated tie at the
bottom of the ranking.

Neither is a missing measurement and neither may be imputed. The hard ceiling is DIAMOND's own
hit rate: *any* SwissProt hit exists for only **81.6% of Kp and 73.0% of Sa** (99.4% of Ec).

THE 25% FLOOR IS A MEASURED TRADE-OFF, NOT A CALIBRATION
---------------------------------------------------------
CLAUDE.md's >= 40% rule is for annotation **transfer** -- putting a GO term or an EC number on a
protein. This axis asks something weaker: *is this protein's family studied at all?* That is
homology detection, the task `orthology/orthodb.py` calibrated at 25%/50%.

**Decoys do not decide it.** Composition-preserving shuffles of our own 13,020 sequences match
SwissProt at **0.0% at every floor down to 20%** (`evidence/decoy_calibration.tsv`), so the floor
is not defending against spurious homology -- DIAMOND's e-value already does that.

**The held-out control decides it, and it is a genuine trade-off** -- coverage against fidelity,
measured in `evidence/floor_sensitivity.tsv` under band-first donor selection. A looser floor
reaches more of the proteome and carries literature further than it safely travels. 25% is the
chosen point because it sits at the coverage ceiling (DIAMOND finds *any* SwissProt hit for only
81.6% of Kp and 73.0% of Sa, so there is almost nothing left to gain below it) while the control
stays above its floor. **Read `evidence` and `donor_pident`**: a `swissprot_remote` score rests on
a 25-60% identity homolog and is a weaker claim than a `swissprot_direct` one. The sweep is
shipped so the floor can be re-chosen without re-running DIAMOND.

Unknome knownness (`load_unknome`) is a second opinion from a different quantity -- a weighted GO
term count over the protein's PANTHER family, not a literature count. It covers only Kp 71.2% /
Ec 78.2% / Sa 65.1%, capped by the PANTHER xref, which is why it is evidence and not a deliverable
column.

Two caveats that belong with any number quoted from here: a literature count measures how popular
the **organism** is as much as the protein, and studiedness is not druggability -- the `ligands/`
axis measures that separately.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
STUDIEDNESS_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = STUDIEDNESS_DIR / "evidence"
SCRATCH_DIR = STUDIEDNESS_DIR / "scratch"
# The three bacteria. Human is out of scope for a bacterial target-prioritization axis.
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# The blend. Literature carries more weight than curation depth because it is the quantity the
# axis is named for; annotation score is kept because it separates entries that share a pub count.
# Both components are in [0, 1], so the score is too.
W_LITERATURE = 0.6
W_ANNOTATION = 0.4

# The quantile of the SwissProt publication distribution that anchors the log scale. 99, not 95 --
# measured; see the module docstring. A run whose scores saturate above this fraction is refused.
SCALE_QUANTILE = 99.0
MAX_SATURATED = 0.20

# Nested identity bands for the SwissProt transfer, and the floor below which nothing transfers.
# 25% is decoy-calibrated (see the module docstring); 95% is CLAUDE.md's "direct" band; 60% is the
# ligands axis's "close". Coverage floor 50% on BOTH query and subject, so a single domain cannot
# claim a whole protein -- measured to cost only ~2 points of coverage.
IDENTITY_DIRECT = 95.0
IDENTITY_CLOSE = 60.0
IDENTITY_FLOOR = 40.0
COVERAGE_FLOOR = 50.0


# ---------------------------------------------------------------- the scale (shared by scripts)

def literature_score(n_pubs, p_ref: float):
    """Publication count on [0, 1], log-compressed against a FIXED global reference count.

    Log because the distribution spans four orders of magnitude (1 to 20,402) and the difference
    between 1 and 5 papers matters far more than between 200 and 204. Clipped at 1 so the handful
    of entries above `p_ref` do not compress everyone else -- but `p_ref` is chosen high enough
    (P99) that clipping is rare; see the module docstring on why P95 was wrong.
    """
    if p_ref <= 0:
        raise ValueError("p_ref must be positive; read it from evidence/scale.tsv")
    return np.minimum(1.0, np.log1p(np.asarray(n_pubs, dtype=float)) / np.log1p(p_ref))


def annotation_component(annotation_score):
    """UniProt's 1-5 annotation score on [0, 1]. 1 is the floor, not zero -- every entry has one."""
    a = np.clip(np.asarray(annotation_score, dtype=float), 1.0, 5.0)
    return (a - 1.0) / 4.0


def score(n_pubs, annotation_score, p_ref: float):
    """The studiedness blend. Used for BOTH `_own` and `_family` so the two are comparable."""
    return (W_LITERATURE * literature_score(n_pubs, p_ref)
            + W_ANNOTATION * annotation_component(annotation_score))


def scale() -> dict:
    """The fixed global constants the scores were built on."""
    path = EVIDENCE_DIR / "scale.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/studiedness/transfer.py first")
    row = pd.read_csv(path, sep="\t").iloc[0]
    return row.to_dict()


# ---------------------------------------------------------------- the deliverable

def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")


def _read(path: Path, hint: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run {hint} first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _numeric(df: pd.DataFrame, cols: tuple[str, ...]) -> pd.DataFrame:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load(species: str) -> pd.DataFrame:
    """The deliverable: one row per protein, canonical order, no nulls in either score.

    `studiedness_own` is near-flat on Kp and Sa -- see the module docstring. Rank on
    `studiedness_family`, and read `evidence` before calling a low score novelty.
    """
    _check(species)
    df = _read(STUDIEDNESS_DIR / f"studiedness_{species}.tsv",
               "scripts/studiedness/merge.py")
    return _numeric(df, ("studiedness_own", "studiedness_family"))


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All three species stacked, with `species` inserted first."""
    frames = []
    for sp in species:
        d = load(sp)
        d.insert(0, "species", sp)
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def novelty(species: str) -> pd.DataFrame:
    """`uniprot_ac` and `novelty` = 1 - studiedness_family, highest first.

    The axis read in the direction the collaboration usually wants it. `evidence == 'no_homolog'`
    rows sit at exactly 1.0 and are the strongest claim the axis makes -- nothing among 575,748
    curated entries resembles them.
    """
    d = load(species)
    out = d[["uniprot_ac", "evidence"]].copy()
    out["novelty"] = 1.0 - d["studiedness_family"]
    # Ties at 1.0 are broken toward the STRONGER claim: `no_hit` (nothing in SwissProt resembles
    # this protein) outranks `below_floor` (a distant relative exists, just too far to transfer
    # from). Both score 0 on studiedness_family, so without this the top of the list would be
    # ordered by accession, which means nothing.
    rank = {"no_hit": 0, "below_floor": 1}
    out["_tier"] = out["evidence"].map(rank).fillna(2)
    return (out.sort_values(["novelty", "_tier"], ascending=[False, True], kind="mergesort")
               .drop(columns="_tier").reset_index(drop=True))


# ---------------------------------------------------------------- the evidence

def load_own(species: str) -> pd.DataFrame:
    """Per-protein UniProt signals: annotation score, protein existence, both pub counts."""
    _check(species)
    df = _read(EVIDENCE_DIR / f"own_{species}.tsv", "scripts/studiedness/transfer.py")
    return _numeric(df, ("annotation_score", "n_pubs_uniprot", "n_pubs_gene2pubmed", "n_pubs"))


def load_transfer(species: str) -> pd.DataFrame:
    """The chosen SwissProt donor per protein, with its identity, organism and counts."""
    _check(species)
    df = _read(EVIDENCE_DIR / f"transfer_{species}.tsv", "scripts/studiedness/transfer.py")
    return _numeric(df, ("donor_pident", "donor_qcov", "donor_n_pubs",
                         "donor_n_pubs_uniprot", "donor_n_pubs_gene2pubmed",
                         "donor_annotation_score", "n_candidates",
                         "nearest_pident", "nearest_n_pubs",
                         "studiedness_own", "studiedness_family",
                         "studiedness_family_bacteria", "studiedness_family_any"))


def load_unknome(species: str) -> pd.DataFrame:
    """Unknome family knownness. A GO-term count, NOT a literature count -- see the docstring.

    Coverage is capped by the PANTHER xref (Kp 71.2% / Ec 78.2% / Sa 65.1%). A protein with no
    family has an EMPTY `unknome_knownness`, never a zero: a zero would claim the family is
    unstudied, which is a much stronger statement than "no PANTHER family".
    """
    _check(species)
    df = _read(EVIDENCE_DIR / f"unknome_{species}.tsv", "scripts/studiedness/unknome.py")
    return _numeric(df, ("unknome_knownness", "unknome_num_species"))


def load_gene2pubmed(species: str) -> pd.DataFrame:
    """NCBI literature counts beside UniProt's, per anchor protein."""
    _check(species)
    df = _read(EVIDENCE_DIR / f"gene2pubmed_{species}.tsv", "scripts/studiedness/gene2pubmed.py")
    return _numeric(df, ("n_pubs_uniprot", "n_pubs_gene2pubmed", "n_pubs"))


def load_route_comparison() -> pd.DataFrame:
    """SwissProt homology against the free four-species ortholog table. Nothing is merged."""
    return _read(EVIDENCE_DIR / "route_comparison.tsv", "scripts/studiedness/transfer.py")


def load_donor_scope_comparison() -> pd.DataFrame:
    """Bacteria-only donors against unrestricted donors, per species.

    Both scopes are computed on every run, so this is a measurement rather than an argument.
    `studiedness_family_bacteria` and `studiedness_family_any` both ship in `load_transfer()`;
    `donor_scope` there says which one became the deliverable.
    """
    return _read(EVIDENCE_DIR / "donor_scope_comparison.tsv", "scripts/studiedness/transfer.py")


def load_floor_sensitivity() -> pd.DataFrame:
    """Coverage and control score across identity floors -- the one arbitrary number, swept."""
    return _read(EVIDENCE_DIR / "floor_sensitivity.tsv", "scripts/studiedness/transfer.py")


def control() -> pd.DataFrame:
    """The E. coli held-out transfer control: predicted family score vs its own measured score.

    E. coli is the only anchor with real measured literature, so it is the only place the transfer
    mechanism can be tested. Every E. coli donor is removed from SwissProt and the family score
    recomputed; the correlation against `studiedness_own` is a genuine held-out result.
    """
    return _read(EVIDENCE_DIR / "control_ecoli_heldout.tsv", "scripts/studiedness/transfer.py")


def manifest() -> pd.DataFrame:
    """One row per species: coverage, evidence-tier counts, medians."""
    return _read(EVIDENCE_DIR / "manifest.tsv", "scripts/studiedness/merge.py")
