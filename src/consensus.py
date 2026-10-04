"""The two standard columns an axis ships: `<axis>_consensus` and `<axis>_evidence`.

Each deliverable answers one question in its own columns, on its own scale. These two make the
axes comparable without a reader knowing any axis's internals:

    <axis>_consensus   float 0-1   higher = more of whatever the axis measures
    <axis>_evidence    int 1-3     how well corroborated that value is

**THE RULE IS: EVIDENCE ALWAYS, CONSENSUS WHERE THE AXIS HAS A MAGNITUDE.** A consensus needs a
0-1 "how much of this property" to express; function and proteomes have none -- *how much function
does a protein have* and *how much identity* are not quantities -- so those two ship the evidence
column alone. A column is omitted, never shipped all-null: a null here means a real absence, and
an all-null column invites someone to fill it.

**EMBEDDINGS SHIPS NEITHER, by decision (project owner, 2026-10-04).** `projection_<sp>.tsv` is
two t-SNE coordinates: there is no magnitude to rank and no evidence to grade, because the
coordinates are a deterministic reduction of an embedding rather than a claim about the protein.
Forcing the pair on would produce two columns that say the same thing for all 13,020 rows. The
gap is deliberate -- do not "complete" it.

**PREFIXED WITH THE AXIS, and that is not decoration.** The tables share a row set and row order
so any two stack with `pd.concat(axis=1)`. A bare `consensus` would give ten identically-named
columns in a merged frame; pandas keeps duplicates silently and name-based selection returns a
DataFrame instead of a Series, so the error surfaces far from its cause.

THE CONSENSUS IS A MEAN OF PERCENTILE RANKS, NOT OF VALUES
-----------------------------------------------------------
An axis's inputs are usually not comparable as numbers -- essentiality averages a
proteome-relative min-max score with two probabilities on different calibrations. Ranking first
makes the mean scale-free, and it matches the project's standing rule: **rank within a species,
never across**. `<axis>_consensus` is therefore a within-species quantity and comparing it between
species is a mistake the number will not stop you making.

**Ties take average ranks.** A block of equal values shares one percentile and contributes no
ordering within itself. That is the honest behaviour and it matters here: `geptop_ess` is 66% ties
at exactly 0 on K. pneumoniae, so across two-thirds of the anchor it adds a constant and the other
inputs carry the ranking.

THE EVIDENCE LADDER IS COUNT + CONCORDANCE
--------------------------------------------
    3  two or more independent experimental sources, agreeing with each other AND with
       the consensus
    2  exactly one source, OR two or more that conflict -- with each other or with the consensus
    1  no experimental measurement at all: prediction only

**Level 0 does not exist.** Every protein in every axis has a prediction, so "no evidence
whatever" is unreachable and a 0 would mean the source count is broken.

**`<axis>_evidence` IS NOT PURELY EXPERIMENTAL -- do not read level 2 as "the experiment was
weak".** A protein measured twice, unanimously, whose consensus contradicts the measurements lands
at 2, not 3. That is deliberate (project owner, 2026-10-04): the two columns are meant to be read
together, so a high consensus beside a 3 means "corroborated, and the models agree". The cost is
real and worth stating -- a model error can lower the stated evidence for a sound experiment --
which is why each axis also writes a per-protein source breakdown to its `evidence/` directory,
where what was actually measured stays checkable.

CONCORDANCE USES A BASE-RATE CUT
----------------------------------
A consensus percentile is a ranking, not a probability of the positive class, so "does the
prediction agree with this measurement" needs a cut. The cut is the measured base rate: a positive
call agrees when the protein sits in the species' top `base_rate` fraction. This is the pattern
`src.degradability.hits()` already uses, for the same reason -- 0.5 is the wrong threshold when
positives are 10% of the proteome. **Derive the rate per run and record it; never hard-code it.**
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# The ladder. Named so a caller never writes a bare integer into a deliverable.
EVIDENCE_CORROBORATED = 3    # >=2 sources, unanimous, concordant with the consensus
EVIDENCE_PARTIAL = 2         # exactly one source, or several that conflict
EVIDENCE_PREDICTED = 1       # no experimental measurement
EVIDENCE_LEVELS = (EVIDENCE_PREDICTED, EVIDENCE_PARTIAL, EVIDENCE_CORROBORATED)


def percentile_consensus(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Mean within-species percentile rank over `columns`, as a float in [0, 1].

    Every column is ranked independently with `pct=True` and average tie handling, then averaged
    with equal weight. A column that is entirely null contributes nothing rather than poisoning
    the mean; a column that is entirely constant contributes a constant, which is correct -- it
    carries no ordering.

    Raises if a named column is missing, because a consensus silently computed over two of three
    sources is a different quantity wearing the same name.
    """
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"consensus inputs absent: {missing}. Present: {list(df.columns)}")
    ranks = pd.DataFrame(index=df.index)
    for c in columns:
        s = pd.to_numeric(df[c], errors="coerce")
        if s.notna().sum() == 0:
            continue
        ranks[c] = s.rank(pct=True, method="average", na_option="keep")
    if ranks.empty:
        raise ValueError(f"every consensus input was null or non-numeric: {columns}")
    return ranks.mean(axis=1, skipna=True).astype(float)


def base_rate_cut(consensus: pd.Series, base_rate: float) -> float:
    """The consensus value above which a protein counts as a predicted positive.

    `base_rate` is the measured fraction of positives. The cut is its empirical quantile, so the
    predicted-positive set is the same size as the measured one -- which is what makes
    "prediction agrees with measurement" a fair question rather than a thresholding artifact.
    """
    if not 0 < base_rate < 1:
        raise ValueError(f"base_rate must be in (0, 1), got {base_rate!r}")
    return float(consensus.quantile(1.0 - base_rate))


def evidence_level(n_sources: pd.Series, unanimous: pd.Series,
                   concordant: pd.Series) -> pd.Series:
    """The 1-3 ladder from its three ingredients, as a nullable integer.

    All three arguments are per-protein and aligned. `unanimous` and `concordant` are only
    consulted where `n_sources >= 2` and `n_sources >= 1` respectively, so their value where no
    source exists is irrelevant -- callers need not invent one.
    """
    n = pd.to_numeric(n_sources, errors="coerce").fillna(0).astype(int)
    out = np.where(
        n == 0, EVIDENCE_PREDICTED,
        np.where((n >= 2) & unanimous.fillna(False) & concordant.fillna(False),
                 EVIDENCE_CORROBORATED, EVIDENCE_PARTIAL))
    return pd.Series(out, index=n_sources.index, dtype="Int64")
