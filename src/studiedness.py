"""Load the studiedness axis -- how much is already known about each protein.

The axis is wanted in both directions: an uncharacterised target is a risk, but it is also the
novelty the GraDi collaboration is looking for. `novelty()` reads the same data the other way.

    data/processed/studiedness/studiedness_<species>.tsv        THE DELIVERABLE
        uniprot_ac · n_papers_uniprot_own · n_papers_uniprot_prokaryotic
                   · n_papers_pubtator_prokaryotic

THREE counts: SOURCE (uniprot = curated references a curator read; pubtator = PubTator3 text
mining) x SCOPE (own = this accession; prokaryotic = its best-studied prokaryotic SwissProt
homolog -- NOT literally this one accession: `_own` is the ligands axis's `exact` rule,
the union over accession + identical sequence + same species at >=95%, because a protein does
not stop being itself between strains) -- minus the pubtator/own cell, dropped 2026-10-03 because text mining does not rescue a
dark anchor: it was median 0 on Kp with only 2.2% of proteins above zero, and median 0 on Sa.
**Never summed, never max()-ed, never blended.** `evidence` is NOT here -- it ships in
`load_transfer()`.

THE NUMBER IS A PAPER COUNT. THAT IS THE WHOLE DEFINITION.
-----------------------------------------------------------
`n_papers_uniprot_prokaryotic` is the number of **curated references on the best-studied prokaryotic SwissProt
homolog** -- papers a UniProt curator actually read and used to annotate that protein. Nothing is
scaled, weighted or blended. A 5 is five papers. A 0 is zero papers, not "unknown".

**`_prokaryotic` here is NOT the ligands axis's `_bacterial`, and harmonizing them would be
wrong.** A donor qualifies if its lineage lacks Eukaryota -- Bacteria, Archaea **and phages** --
because of the 93 proteins a strict bacterial rule stranded with no donor, 70 lost theirs to a
virus. `ligands_<species>.tsv` restricts to true Bacteria, for its own measured reason (an
unrestricted pool gave 424 potent Kp proteins against a true 175). Different populations,
different names, both deliberate.

`n_papers_uniprot_own` is the same count on THIS accession. It is near-constant on Kp (the genome paper)
and Sa **by design** -- it is the measurement of darkness, and the gap between the two columns is
the axis's entire product. Do not rank Kp or Sa on `n_papers_uniprot_own`; rank on `n_papers_uniprot_prokaryotic`.

WHY NOT A 0-1 SCORE: THE BLEND THAT SHIPPED FIRST WAS UNINTERPRETABLE
-----------------------------------------------------------------------
An earlier version shipped `0.6*min(1, log1p(n)/log1p(204)) + 0.4*(annotation_score-1)/4`. It was
removed on 2026-09-22 after three measurements, and **must not be reintroduced**:

1. **The same value meant different things.** At 0.5: 13 papers if the donor's annotation score
   was 3, **83** if it was 1, **1** if it was 5. A protein with zero papers scored 0.4 when its
   donor was annotation-5.
2. **The halves double-counted.** UniProt's annotation score is largely a function of how much is
   known, so it correlated with the paper count at **r = 0.64-0.70**.
3. **The weights did the opposite of what they claimed.** Literature was weighted 0.6 "because it
   is the quantity the axis is named for", but the annotation component's spread was nearly
   double (sd 0.33-0.36 against 0.17-0.20), so annotation swung the score MORE.

Dropping it moved the ranking by spearman 0.86-0.91, so this was a real change, not a relabelling.
`scaled()` below returns a 0-1 version on the fly for anyone combining this axis with the others;
it is deliberately **not stored**, so there is exactly one source of truth on disk.

THE COUNT IS SWISSPROT-CURATED ONLY, AND THERE IS NO CEILING
--------------------------------------------------------------
NCBI gene2pubmed was measured and is NOT used for the shipped number (see
`load_gene2pubmed()`), though it ships beside it as evidence.

**There is no cap.** An earlier note in this project claimed SwissProt reference counts "saturate
at 58"; that was wrong and is corrected here. Verified against UniProt: the TSV export is not
truncated (human TP53 returns 225 PubMed ids), counts reach **225 overall and 119 among
prokaryotic entries**, and the 58 is simply *E. coli* GroEL (`P0A6F5`) happening to be the
most-curated donor selected in all three species.

The dynamic range is small because bacterial proteins are: median donor 4-6 papers, 37-48 distinct
values. **The ties sit in the poorly-studied bulk, not at the top** -- measured, the three highest
values are held by 1-2 proteins each while the largest single tie is 460 Kp proteins at 4 papers
(8.0% of the proteome). A shortlist reads the top, where the ranking is near-unique.

DONORS ARE PROKARYOTIC: BACTERIA, ARCHAEA AND PHAGES -- NEVER EUKARYOTES
-------------------------------------------------------------------------
A donor qualifies if its UniProt taxonomic lineage lacks `Eukaryota (domain)`. Without the
restriction, five of K. pneumoniae's ten highest-scoring proteins took human donors (HSPD1, HADHA,
CTPS1, LONP1, AFG3L2), so the score claimed "well studied because its human mitochondrial homolog
is" -- true, and the wrong question for an antibacterial target.

**Phages are deliberately IN, and "restrict to Bacteria" is the obvious rule that is wrong.** Of
the 93 proteins a strict bacterial rule stranded with no donor, **70 lost theirs to a virus**
(Escherichia phage lambda, P1). A prophage protein whose best-characterised relative is a lambda
protein is not novel, and scoring it 0 fails in the one direction this axis must not fail.

**Scope barely moves the numbers** -- the three scores correlate at rho 0.983-0.990 and only 254
of 5,728 Kp proteins change, because the log scale saturates above P_REF. It is kept because it
makes the number mean the right thing, not because it changes it. See `docs/studiedness.md` §2b.

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

# Reference count used by `scaled()` as the 0-1 anchor. The 99th percentile of the SwissProt
# publication distribution; nothing on disk depends on it.
SCALE_REFERENCE = 204.0

# No single NON-ZERO paper count may hold more than this share of the scored proteins, or the
# ranking has collapsed into ties. Zero is excluded: a third of Kp genuinely has no in-scope
# homolog, and that group is already split into `no_hit` and `below_floor`.
MAX_TIE_FRACTION = 0.25

# PubTator3 symbol-count inflation over the curated count, above which `pubtator_ambiguous` is
# set. **THERE IS NO NATURAL CUT HERE -- the ratio is a smooth heavy tail** (p50 0.5 · p75 1.9 ·
# p90 9.2 · p95 39.8 · p99 621), so this is a CONSERVATISM CHOICE, not an accuracy threshold --
# the same status `ligands/transfer_calibration.py` records for its identity bands. 50 sits at
# about p95.5. The RATIO ITSELF SHIPS (`n_papers_pubtator_prokaryotic_ratio`), so anyone who wants a
# different cut sets one without re-running anything.
PUBTATOR_RATIO_FLAG = 50.0

# Nested identity bands for the SwissProt transfer, and the floor below which nothing transfers.
# 25% is decoy-calibrated (see the module docstring); 95% is CLAUDE.md's "direct" band; 60% is the
# ligands axis's "close". Coverage floor 50% on BOTH query and subject, so a single domain cannot
# claim a whole protein -- measured to cost only ~2 points of coverage.
IDENTITY_DIRECT = 95.0
IDENTITY_CLOSE = 60.0
IDENTITY_FLOOR = 40.0
COVERAGE_FLOOR = 50.0


# ---------------------------------------------------------------- derived, never stored

def scaled(n_papers, ref: float = SCALE_REFERENCE):
    """Paper counts compressed to [0, 1], for combining this axis with the 0-1 ones.

    `log1p(n)/log1p(ref)` clipped at 1. Log because the distribution spans four orders of
    magnitude and the step from 1 to 5 papers matters far more than 200 to 204.

    **Derived on the fly and deliberately not shipped.** The deliverable carries the integer so
    there is one source of truth; anyone who needs a 0-1 column calls this and records the `ref`
    they used. It is NOT the old composite score -- that blended in UniProt's annotation score and
    was removed for being uninterpretable (see the module docstring).
    """
    if ref <= 0:
        raise ValueError("ref must be positive")
    return np.minimum(1.0, np.log1p(np.asarray(n_papers, dtype=float)) / np.log1p(ref))


def definition() -> dict:
    """What the shipped numbers actually mean: count source, donor scope, identity floors."""
    path = EVIDENCE_DIR / "definition.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/studiedness/transfer.py first")
    return pd.read_csv(path, sep="\t").iloc[0].to_dict()


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


def _counts(df: pd.DataFrame, cols: tuple[str, ...]) -> pd.DataFrame:
    """Paper counts as nullable integers -- a count of 4.0 invites being read as a score."""
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


# The three counts the consensus ranks. All three are PAPER COUNTS -- homogeneous units, which is
# the whole reason this is not the composite removed on 2026-09-22 (see `consensus()`).
CONSENSUS_COLUMNS = ("n_papers_uniprot_own", "n_papers_uniprot_prokaryotic",
                     "n_papers_pubtator_prokaryotic")

# The five transfer tiers, collapsed onto the 1-3 ladder. `swissprot_direct` is a donor at >=95%
# identity, i.e. effectively this protein's own literature.
EVIDENCE_DIRECT = ("swissprot_direct",)
EVIDENCE_HOMOLOG = ("swissprot_close", "swissprot_homolog")


def consensus(df: pd.DataFrame) -> pd.Series:
    """`studiedness_consensus`: how studied this protein is, 0-1, WITHIN-SPECIES.

    Mean percentile rank over `CONSENSUS_COLUMNS`, with proteins carrying **no literature at all**
    pinned to exactly 0 -- the zero-floor `src/ligandability.py` introduced, so a 0 means "nothing
    anywhere" rather than an arbitrary mid-scale rank. **It fires only on S. aureus** (1,258
    proteins): on Kp and Ec every protein has at least its genome paper in `n_papers_uniprot_own`,
    so nothing qualifies.

    **THIS IS NOT THE 0-1 COMPOSITE REMOVED ON 2026-09-22, and the difference is the whole point.**
    That one was `0.6 x log-scaled paper count + 0.4 x (annotation_score - 1)/4` -- a literature
    count blended with **UniProt's annotation-quality rating**, two different quantities on two
    scales with hand-invented weights. All three of its recorded faults are faults of that
    heterogeneity: the same value meant 13 papers at annotation 3 and 83 at annotation 1; the two
    halves double-counted at r 0.64-0.70; and the weights did the opposite of what the code claimed,
    because the annotation term's spread was nearly double. **Here every input is a count of papers**
    -- no annotation score, no hand-chosen weights, and percentile ranking removes the unequal-spread
    problem outright. `scaled()` stays unstored and is a different object again: one column,
    log-compressed against a caller-chosen `ref`, magnitude-preserving and cross-species comparable,
    where this is rank-only, reference-free and **within-species**.

    **BE HONEST ABOUT WHAT IT ADDS: little.** It tracks the axis's designated ranking column
    `n_papers_uniprot_prokaryotic` at rho **+0.974 Kp / +0.938 Ec / +0.963 Sa**. Anyone working
    inside this axis should still rank on that column; this exists so a reader stacking ten
    deliverables has one scale. The same candour `degradability_consensus` carries at rho 0.97.

    **`NEVER a max() ACROSS THEM` is respected.** That rule exists because a `max()` or a sum
    silently switches *which definition* a protein's number came from, row by row. An equal-weight
    rank mean uses all three for every protein, so no row changes definition.

    **The big tie block is the data, not this column.** Largest tie 34.3% (Kp) / 3.2% (Ec) / 43.5%
    (Sa) -- above the axis's own `MAX_TIE_FRACTION` of 0.25 on two species, **and so is the column
    it is built to summarise**: `n_papers_uniprot_prokaryotic` ties 34.4% / 14.2% / 44.2% and
    `n_papers_uniprot_own` 95.1% / 9.7% / 63.7%. A third of K. pneumoniae simply has no curated
    donor at 40% identity. This column is marginally *better* tied than the incumbent, not worse.

    **`n_papers_pubtator_prokaryotic` carries real blanks** (Kp 2,149 · Ec 190 · Sa 1,356) meaning
    "no in-scope donor had a GeneID" -- **not zero**. `percentile_consensus` skips them row-wise, so
    those proteins are ranked on the other two. Never fill them.
    """
    from src import consensus as consensus_mod
    cols = list(CONSENSUS_COLUMNS)
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"studiedness_consensus needs {missing}; present: {list(df.columns)}")
    has = pd.concat([pd.to_numeric(df[c], errors="coerce").fillna(0) for c in cols],
                    axis=1).gt(0).any(axis=1)
    out = pd.Series(0.0, index=df.index, dtype=float)
    if has.any():
        out.loc[has] = consensus_mod.percentile_consensus(df.loc[has], cols).to_numpy()
    return out


def evidence(species: str) -> pd.Series:
    """`studiedness_evidence`, the 1-3 ladder, indexed like `load_transfer(species)`.

        3  swissprot_direct                      a donor at >=95% -- this protein's own literature
        2  swissprot_close | swissprot_homolog   a real donor at >=40%
        1  below_floor | no_hit                  nothing in SwissProt resembles it

    Counts: Kp 1,961 / 3,403 / 364 · Ec 29 / 519 / 3,855 · Sa 1,276 / 1,281 / 332.

    **It takes a SPECIES, not a frame**, because the tiers live in
    `evidence/transfer_<species>.tsv` rather than in the deliverable -- unlike every other axis's
    `evidence()`. That file is 1:1 with the deliverable on `uniprot_ac`, verified, so the join
    cannot fan out.

    **THIS IS THE COLUMN'S REAL JUSTIFICATION AND IT IS STRONGER HERE THAN ON ANY OTHER AXIS.** The
    deliverable dropped its `evidence` column on 2026-10-03, which left **a 0 in
    `n_papers_uniprot_prokaryotic` ambiguous in the shipped table**: `no_hit` -- the strongest
    novelty claim this axis makes -- and `below_floor` both read 0, and the reader was told to join
    `load_transfer()` before treating a 0 as novelty. **Level 1 recovers exactly that class**: of
    K. pneumoniae's 1,969 proteins reading 0, **1,961 are level 1**; only 8 are "a donor exists and
    has no papers". So the pair puts the distinction back in the deliverable in the project-wide
    form.

    **It still merges `no_hit` with `below_floor`**, the stronger and weaker novelty claims, so the
    five-tier column stays the finer instrument -- kept in `evidence/consensus_audit.tsv` and in
    `load_transfer()`.

    **NOTHING HERE IS AN EXPERIMENT**, so a 3 is not experimental corroboration and a 1 means "no
    donor to look at" rather than "not yet measured" -- the pattern `function` and `orthology` set.

    **DO NOT EXPECT THE CONSENSUS TO RISE WITH THE LEVEL.** Median consensus by level is Kp 0.284 /
    0.603 / 0.741 but **Ec 0.028 / 0.555 / 0.501** and **Sa 0.275 / 0.698 / 0.646** -- level 2 above
    level 3 on two species. That is correct and is the proof the two columns are complementary: on
    E. coli `swissprot_direct` covers 3,855 proteins including obscure ones, while `close`/`homolog`
    are the conserved families whose donors are famous. **Evidence grades how directly the count was
    measured, not how large it is.**
    """
    _check(species)
    tiers = load_transfer(species).set_index("uniprot_ac")["evidence"]
    acc = load(species)["uniprot_ac"]
    t = acc.map(tiers)
    if t.isna().any():
        raise ValueError(
            f"{species}: {int(t.isna().sum())} proteins missing from "
            f"evidence/transfer_{species}.tsv -- run scripts/studiedness/transfer.py first")
    level = pd.Series(1, index=acc.index, dtype="Int64")
    level[t.isin(EVIDENCE_HOMOLOG)] = 2
    level[t.isin(EVIDENCE_DIRECT)] = 3
    return level


def load(species: str) -> pd.DataFrame:
    """The deliverable: one row per protein, canonical order, no nulls in either count.

    THREE literature counts, three definitions, NEVER summed or `max()`-ed together:

        n_papers_uniprot_own            curated references on THIS accession. Near-flat on Kp and
                                        Sa by design — it is the measurement of darkness.
        n_papers_uniprot_prokaryotic    curated references on the best-studied prokaryotic
                                        SwissProt homolog. **THE SHIPPED RANKING — rank on this.**
        n_papers_pubtator_prokaryotic   PubTator3 TEXT-MINED papers on a prokaryotic homolog,
                                        keyed on NCBI GeneID (NOT a gene symbol -- the symbol
                                        routes were measured and rejected). ITS OWN DONOR,
                                        chosen by PubTator count, which differs from the curated
                                        column's donor on 15.5% of Kp proteins.

    **The pubtator/own cell was measured and dropped** (2026-10-03). Text mining does not rescue a
    dark anchor: `n_papers_pubtator_own` was median 0 on Kp with only **2.2%** of proteins above
    zero and median 0 on Sa, against E. coli's median 8 and 91.5%. Where curation is silent,
    PubTator is silent too — the darkness is real, not an artifact of which corpus was searched.

    **The two homolog columns are near-redundant** (Spearman **0.786** on Kp), while
    `n_papers_uniprot_own` is nearly independent of both (0.19 and 0.09). So this is one strong
    axis — homolog literature, counted two ways — plus the own-protein darkness measure.

    **`n_papers_pubtator_prokaryotic` is NOT complete**: filled for 62.5% of Kp, 95.7% of Ec, 53.1%
    of Sa. An empty cell means the donor carries no gene symbol to look up, which is not a zero —
    see the note below.

    **`evidence` IS NOT IN THIS TABLE** (project owner, 2026-10-03) -- it ships in
    `load_transfer(species)`. A 0 in `n_papers_uniprot_prokaryotic` is therefore ambiguous here: it may be
    `no_hit` (nothing in 575,748 curated entries resembles the protein, the strongest novelty
    claim the axis makes) or `below_floor` (a curated relative exists but below 40% identity).
    Join `load_transfer()` on `uniprot_ac` before reading a 0 as novelty.

    **Why the third column exists**: it measurably transfers better -- held-out E. coli control,
    identical folds, **0.4054 against 0.3428** on the common subset, and it reaches 97.9% of
    scored Kp where a GeneID-keyed count reaches 1.8%.

    **Why it is NOT the ranking**: that control is E. coli-only, and E. coli gene symbols are
    precisely the ones that entered human nomenclature, so it is weakest where the failure mode
    lives. The count is also SPECIES-AGNOSTIC -- it cannot be donor-scoped the way the curated
    count is (see `docs/studiedness.md` 2b), so a residual eukaryotic contribution survives.

    **Two kinds of blank in `n_papers_pubtator_prokaryotic`, and they are different claims**: `0`
    means no donor at all (what `n_papers_uniprot_prokaryotic` also says, with `evidence` giving the kind),
    while **empty means a donor exists but carries no gene symbol to look up** -- not a measured
    zero. Kp 80 · Ec 10 · Sa 117 scored proteins.

    **`pubtator_ambiguous`** marks a lookup likely measuring a DIFFERENT protein: a human gene
    symbol (`crp` returns
    345,630 papers for human C-reactive protein, not the cAMP receptor protein. 3.8% of donors;
    they score WORSE than average, so this flags an artifact rather than explaining the column.
    """
    _check(species)
    df = _read(STUDIEDNESS_DIR / f"studiedness_{species}.tsv",
               "scripts/studiedness/merge.py")
    cols = [c for c in ("n_papers_uniprot_own", "n_papers_uniprot_prokaryotic",
                        "n_papers_pubtator_prokaryotic")
            if c in df.columns]
    df = _numeric(df, tuple(cols))
    for c in ("n_papers_uniprot_own", "n_papers_uniprot_prokaryotic"):
        if c in df.columns:
            df[c] = df[c].astype(int)
    # NULLABLE INTEGER, not float. The column carries blanks (no in-scope donor had a GeneID),
    # which forces float and prints a paper count as `14.0`. Int64 keeps the blanks AND the
    # integer semantics -- a count of papers is never fractional.
    if "n_papers_pubtator_prokaryotic" in df.columns:
        df["n_papers_pubtator_prokaryotic"] = df["n_papers_pubtator_prokaryotic"].astype("Int64")
    # The two standard columns. `_read()` returns EVERY column as `str` and this function coerces a
    # hand-kept list, with no prefix rule to piggyback on -- so an unregistered column stays text,
    # `>= 0.5` raises TypeError and a sort puts '0.9' above '0.12'. Fifth axis to need this.
    if "studiedness_consensus" in df.columns:
        df["studiedness_consensus"] = pd.to_numeric(df["studiedness_consensus"], errors="coerce")
    if "studiedness_evidence" in df.columns:
        df["studiedness_evidence"] = pd.to_numeric(
            df["studiedness_evidence"], errors="coerce").astype("Int64")
    return df


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All three species stacked, with `species` inserted first."""
    frames = []
    for sp in species:
        d = load(sp)
        d.insert(0, "species", sp)
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def novelty(species: str) -> pd.DataFrame:
    """The axis read the way the collaboration usually wants it: least-studied first.

    **Returns an ORDERING, not a novelty number.** The old version returned `1 - score`, which
    only worked while the axis shipped a bounded composite; inventing a bounded "novelty" from an
    unbounded count would be exactly the kind of opaque derived number this axis just removed.
    The sort order is the answer.

    Ties at 0 papers are broken toward the stronger claim: `no_hit` (nothing in SwissProt
    resembles this protein) ranks above `below_floor` (a distant relative exists, too far to
    transfer from).
    """
    # `evidence` left the deliverable on 2026-10-03, so the tier is pulled from the transfer
    # table. Without it the ties at 0 are unordered and the strongest novelty claim the axis
    # makes -- `no_hit` -- would be scattered through a third of the proteome.
    d = load(species)
    tiers = load_transfer(species)[["uniprot_ac", "evidence"]]
    out = d[["uniprot_ac", "n_papers_uniprot_prokaryotic"]].merge(
        tiers, on="uniprot_ac", how="left")
    rank = {"no_hit": 0, "below_floor": 1}
    out["_tier"] = out["evidence"].map(rank).fillna(2)
    return (out.sort_values(["n_papers_uniprot_prokaryotic", "_tier"], ascending=[True, True],
                            kind="mergesort")
               .drop(columns="_tier").reset_index(drop=True))


# ---------------------------------------------------------------- the evidence

def load_own(species: str) -> pd.DataFrame:
    """Per-protein UniProt signals: annotation score, protein existence, both pub counts."""
    _check(species)
    df = _read(EVIDENCE_DIR / f"own_{species}.tsv", "scripts/studiedness/transfer.py")
    return _counts(_numeric(df, ("annotation_score",)),
                   ("n_pubs_uniprot", "n_pubs_gene2pubmed", "n_papers_uniprot_own"))


def load_transfer(species: str) -> pd.DataFrame:
    """The chosen SwissProt donor per protein, with its identity, organism and counts."""
    _check(species)
    df = _read(EVIDENCE_DIR / f"transfer_{species}.tsv", "scripts/studiedness/transfer.py")
    df = _numeric(df, ("donor_pident", "donor_qcov", "donor_annotation_score",
                       "n_candidates", "nearest_pident"))
    return _counts(df, ("donor_n_pubs_uniprot", "donor_n_pubs_gene2pubmed", "nearest_n_pubs",
                        "n_papers_uniprot_own", "n_papers_uniprot_prokaryotic",
                        "n_papers_uniprot_prokaryotic", "n_papers_uniprot_bacteria",
                        "n_papers_uniprot_any"))


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
    """NCBI literature counts beside UniProt's, per anchor protein.

    **A MEASURED ALTERNATIVE, not the shipped number.** gene2pubmed is larger for 87-94% of
    donors (median 2.8x) and gives 184 S. aureus proteins their first paper, but the deliverable
    counts UniProt's curated references only -- papers a curator read, one consistent definition.
    Kept and documented rather than dropped, the `interpro2go` precedent.
    """
    _check(species)
    df = _read(EVIDENCE_DIR / f"gene2pubmed_{species}.tsv", "scripts/studiedness/gene2pubmed.py")
    return _counts(df, ("n_pubs_uniprot", "n_pubs_gene2pubmed", "n_pubs"))


def load_route_comparison() -> pd.DataFrame:
    """SwissProt homology against the free four-species ortholog table. Nothing is merged."""
    return _read(EVIDENCE_DIR / "route_comparison.tsv", "scripts/studiedness/transfer.py")


def load_donor_scope_comparison() -> pd.DataFrame:
    """Bacteria-only donors against unrestricted donors, per species.

    All THREE scopes are computed on every run, so this is a measurement rather than an argument.
    `n_papers_uniprot_prokaryotic`, `_bacteria` and `_any` all ship in `load_transfer()`, and
    `donor_scope` there names the one that became the deliverable. Switching is a column swap.
    """
    return _read(EVIDENCE_DIR / "donor_scope_comparison.tsv", "scripts/studiedness/transfer.py")


def load_floor_sensitivity() -> pd.DataFrame:
    """Coverage and control score across identity floors -- the one arbitrary number, swept."""
    return _read(EVIDENCE_DIR / "floor_sensitivity.tsv", "scripts/studiedness/transfer.py")


def load_pubtator(species: str) -> pd.DataFrame:
    """PubTator3 text-mined counts per anchor protein, by BOTH routes. Evidence, not the axis.

    Written by `scripts/studiedness/pubtator.py`. Columns:

        n_pubs_pubtator_geneid              keyed on NCBI GeneID, from the 756 MB bulk file
        n_pubs_pubtator_symbol_anyspecies   keyed on gene SYMBOL, from the PubTator3 API

    **THE TWO ARE NOT THE SAME DATASET FOR BACTERIA.** PubTator3 normalises bacterial gene
    mentions to species-agnostic symbol concepts rather than to strain GeneIDs, so the GeneID
    route is nearly blind outside E. coli while the symbol route carries real volume.

    **`_symbol_anyspecies` IS NAMED FOR ITS LIMITATION.** It pools every organism with a gene of
    that name -- `@GENE_CLPP` sums E. coli, S. aureus, human mitochondrial and plant CLPP -- so
    it **cannot be donor-scoped** and must never be merged into a scoped count. Using it without
    that caveat reintroduces "well studied because its human homolog is", which
    `docs/studiedness.md` 2b rejects for donor scope.

    **Three definitions mean three columns, never a `max()` across them** -- mixing definitions
    per protein is what killed the 0-1 composite on 2026-09-22.
    """
    _check(species)
    return _read(EVIDENCE_DIR / f"pubtator_{species}.tsv", "scripts/studiedness/pubtator.py")


def load_pubtator_route_comparison() -> pd.DataFrame:
    """What each PubTator3 route reaches per species: coverage, distinct values, largest tie."""
    return _read(EVIDENCE_DIR / "pubtator_route_comparison.tsv",
                 "scripts/studiedness/pubtator.py")


def load_confounds() -> pd.DataFrame:
    """How studiedness relates to every other axis -- a measurement, not a decision.

    Written by `scripts/studiedness/confounds.py`. One row per
    (species, axis, endpoint, studiedness_column), carrying spearman `rho`, the same rho over the
    SCORED proteins alone (`scored_rho`), and for binary endpoints AUROC, PR-AUC, the base rate,
    and length's AUROC/PR-AUC beside them as the baseline.

    **READ `rho` AND `scored_rho` TOGETHER.** A third of K. pneumoniae scores 0 in two tiers, and
    a relationship that collapses toward zero in `scored_rho` was a correlation with "did DIAMOND
    find a donor", not with studiedness. Degradability is exactly that case: Kp `adep4_prob` reads
    rho -0.152 overall and +0.008 scored-only.

    **AUROC BELOW 0.5 IS A DIRECTION, NOT A FAILURE** -- it is studiedness predicting the endpoint
    as named, so `goslim_unannotated` reading 0.153 means studiedness predicts being *annotated*
    at 0.847. Left uninverted so the direction of each relationship stays visible.

    **An AUROC without length's beside it means nothing** -- length is the generic confound on
    this project and scores 0.66-0.73 on ligand precedent by itself.
    """
    return _read(EVIDENCE_DIR / "confounds.tsv", "scripts/studiedness/confounds.py")


def control() -> pd.DataFrame:
    """The E. coli held-out transfer control: predicted family score vs its own measured score.

    E. coli is the only anchor with real measured literature, so it is the only place the transfer
    mechanism can be tested. Every E. coli donor is removed from SwissProt and the family score
    recomputed; the correlation against `n_papers_uniprot_own` is a genuine held-out result.
    """
    return _read(EVIDENCE_DIR / "control_ecoli_heldout.tsv", "scripts/studiedness/transfer.py")


def manifest() -> pd.DataFrame:
    """One row per species: coverage, evidence-tier counts, medians."""
    return _read(EVIDENCE_DIR / "manifest.tsv", "scripts/studiedness/merge.py")
