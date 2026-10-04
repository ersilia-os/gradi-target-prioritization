"""The shortlist predicates, defined once, because they are a scientific claim and not a style.

Every other file in `plotting/` copy-pastes its own style block -- that is the house convention in
`scripts/plots/` and it is deliberate. This module is the one exception: the funnel in
`shortlist.py` and the highlighting in `novelty.py` must mean the SAME thing by "a candidate", and
two copies of a filter cascade is how they stop meaning the same thing.

**There is no composite score here, and there must not be one.** Three axes each REMOVED one:
essentiality dropped `essentiality`/`essentiality_source` (it mixed a measured 1.0 with a predicted
1.0, and on Kp was a verbatim copy of `screens_ess_mean`); studiedness removed its 0-1 composite on
2026-09-22 (the two halves double-counted, r 0.64-0.70); pockets ships no `druggability()` because
"no defensible weighting exists across a weak prior, a sparse measurement and a third party's
model". A weighted score reintroduced in a plotting script would be exactly the thing those three
removals were for, and it would be invisible -- a ranked list looks equally plausible whatever
weights produced it.

So the shortlist is a CASCADE OF STATED PREDICATES. Each one is a single comparison with a written
rationale, applied in order, and every step reports its own count. The reader can disagree with any
rule and see what it cost: `cascade()` returns the steps, and `leave_one_out()` prices each filter
by dropping it.

Two of the rules are less obvious than they look:

- **`not_membrane`, NOT "cytoplasmic only".** The measured S. aureus hit rates are cytoplasm
  0.180/0.275, membrane 0.031/0.105 and **extracellular 0.049/0.346** -- a secreted protein transits
  the cytoplasm unfolded and IS reachable by activated ClpP, while a membrane protein is inserted
  co-translationally and never presents a soluble chain. Filtering to the cytoplasm would discard
  the compartment with the HIGHEST measured ONC212 hit rate.
- **`essential` uses a PROTEOME-WIDE quantile, not a fixed value.** `screens_ess_mean` is comparable
  within a species and never across one (`src/essentiality.py` says so), so a hard cut would mean a
  different stringency in each organism.

The thresholds are read from the stages themselves (`degradability.BASE_RATE_THRESHOLD`) rather than
restated, so an estimator change moves them here too -- those cuts are empirical quantiles of the
OOF distribution and `CLAUDE.md` requires them to be re-derived whenever the estimator changes.

Run with the `gradi` env. Imported, never executed as a worker.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import pandas as pd

# The proteome-wide percentile that defines "essential" for the shortlist. Not a measured constant
# -- a stringency choice, stated here so it is arguable.
ESSENTIAL_PERCENTILE = 90.0

# Which activator's probability gates the cascade. The two columns correlate at rho 0.89 while the
# labels agree at only rho 0.52, so they are closer to one opinion than two; gating on both would
# double-count. ADEP4 is the better-measured arm (AUROC 0.8738 / PR 0.6103 vs 0.7671 / 0.5803).
ACTIVATOR = "adep4"


@dataclass(frozen=True)
class Filter:
    """One predicate in the cascade, with the reason it is there."""

    name: str
    column: str
    rationale: str
    predicate: Callable[[pd.DataFrame], pd.Series]
    rule: str  # human-readable, echoed verbatim into the run log and onto the figure

    def apply(self, df: pd.DataFrame) -> pd.Series:
        return self.predicate(df).fillna(False).astype(bool)


@dataclass(frozen=True)
class Step:
    """The outcome of one filter: what it was, and what it cost."""

    name: str
    rule: str
    rationale: str
    n_before: int
    n_after: int

    @property
    def n_removed(self) -> int:
        return self.n_before - self.n_after


def _essential_cut(df: pd.DataFrame) -> float:
    """The proteome-wide essentiality cut. Computed on the FULL frame passed in, which must be the
    whole proteome -- computing it on an already-filtered frame would silently re-centre it."""
    return float(df["screens_ess_mean"].quantile(ESSENTIAL_PERCENTILE / 100.0))


def build_filters(df: pd.DataFrame) -> list[Filter]:
    """The cascade, in order. `df` must be the WHOLE proteome: the essentiality cut is a quantile
    over it, so passing a subset moves the threshold without saying so."""
    cut = _essential_cut(df)
    prob_col = f"{ACTIVATOR}_prob"
    from src import degradability as D

    deg_cut = D.BASE_RATE_THRESHOLD[ACTIVATOR]

    return [
        Filter(
            name="not_membrane",
            column="localization",
            rule='localization != "cytoplasmic_membrane"',
            rationale=(
                "Membrane proteins are inserted co-translationally and never present a soluble "
                "cytoplasmic chain: measured hit rate 0.031/0.105 against cytoplasm's 0.180/0.275. "
                "Secreted proteins are KEPT -- they transit the cytoplasm unfolded and show the "
                "highest measured ONC212 rate (0.346)."
            ),
            predicate=lambda d: d["localization"] != "cytoplasmic_membrane",
        ),
        Filter(
            name="selective",
            column="has_human_ortholog",
            rule="has_human_ortholog == False",
            rationale=(
                "No human ortholog by the union of OrthoFinder and RBH. The error direction "
                "matters: under-detecting human homology makes a target look MORE selective than "
                "it is, so the union (not either method alone) is the conservative choice."
            ),
            predicate=lambda d: ~d["has_human_ortholog"].astype(bool),
        ),
        Filter(
            name="essential",
            column="screens_ess_mean",
            rule=f"screens_ess_mean >= {cut:.4f}  (proteome-wide {ESSENTIAL_PERCENTILE:.0f}th pct)",
            rationale=(
                "Top decile of the screens-trained transfer models, which reach AUROC 0.89-0.96 on "
                "the three measured Kp screens. A quantile, not a fixed value, because the column "
                "is comparable within a species and never across one. Note these are NOT nine "
                "independent votes -- all nine read the same ProtT5 embedding."
            ),
            predicate=lambda d: d["screens_ess_mean"] >= cut,
        ),
        Filter(
            name="degradable",
            column=prob_col,
            rule=f"{prob_col} >= {deg_cut}  (base-rate cut)",
            rationale=(
                f"The re-derived base-rate cut for {ACTIVATOR.upper()}, not 0.5: probabilities top "
                "out at 0.868 on Kp, so 0.5 is the wrong threshold and it MOVED with the "
                "estimator. Kp rows are ranking hypotheses extrapolated from S. aureus labels, "
                "never measurements."
            ),
            predicate=lambda d: d[prob_col] >= deg_cut,
        ),
    ]


def cascade(df: pd.DataFrame) -> tuple[pd.DataFrame, list[Step]]:
    """Apply the filters in order. Returns the survivors and one `Step` per filter.

    `df` must be the whole proteome, one row per protein, already joined across the axes the
    predicates name (see `shortlist.py:load_joined`)."""
    steps: list[Step] = []
    cur = df
    for f in build_filters(df):
        before = len(cur)
        cur = cur[f.apply(cur)]
        steps.append(Step(f.name, f.rule, f.rationale, before, len(cur)))
    return cur.copy(), steps


def leave_one_out(df: pd.DataFrame) -> dict[str, int]:
    """What each filter costs: the survivor count with that one rule DROPPED.

    This is the honest way to show a cascade. A funnel alone makes the order look load-bearing when
    it is not -- the filters are conjunctive, so the final count is order-independent, and what a
    reader actually wants to know is which rule is doing the work."""
    out: dict[str, int] = {}
    filters = build_filters(df)
    for skip in filters:
        cur = df
        for f in filters:
            if f.name == skip.name:
                continue
            cur = cur[f.apply(cur)]
        out[skip.name] = len(cur)
    return out


# -- the data contract -------------------------------------------------------------------------
# The predicates above name columns from five different axes. The join that produces them lives
# here rather than in each figure script, because a filter and the frame it is applied to are one
# thing: a script that built the frame differently would run the same rules against different data
# and never say so.

#: axis module -> the columns the cascade (or the shortlist table) needs from it.
JOIN_COLUMNS: dict[str, tuple[str, ...]] = {
    "proteomes": ("gene_name", "protein_name"),
    "localization": ("localization", "cytoplasmic_fraction"),
    "orthology": ("has_human_ortholog", "bacterial_panel_orthologs"),
    "essentiality": ("geptop_ess", "proteomelm_ess", "screens_ess_mean"),
    "degradability": ("adep4_prob", "onc212_prob", "nn_similarity"),
    "studiedness": ("n_papers_uniprot_prokaryotic",),
    "ligands": ("n_ligands_bacterial", "n_assayed_bacterial"),
    "pockets": ("p2rank_score", "n_pdb_structures", "af_plddt"),
}


def load_joined(species: str = "kpneumoniae") -> pd.DataFrame:
    """One row per protein, in canonical order, carrying every column the cascade reads.

    Joins are INNER on `uniprot_ac` and the result is asserted to keep the full proteome: every
    axis ships one complete row per protein in the same order, so a short join means an axis has
    regressed, and silently dropping proteins is exactly the failure that invariant exists to
    catch. Loaded through the `src/` loaders only -- nothing is read by path and nothing is
    recomputed."""
    from src import (
        degradability,
        essentiality,
        ligandability,
        localization,
        orthology,
        pockets,
        proteomes,
        studiedness,
    )

    modules = {
        "proteomes": proteomes,
        "localization": localization,
        "orthology": orthology,
        "essentiality": essentiality,
        "degradability": degradability,
        "studiedness": studiedness,
        "ligands": ligandability,
        "pockets": pockets,
    }

    base = proteomes.load(species)
    n_proteome = len(base)
    df = base[["uniprot_ac", *JOIN_COLUMNS["proteomes"]]]

    for axis, cols in JOIN_COLUMNS.items():
        if axis == "proteomes":
            continue
        part = modules[axis].load(species)[["uniprot_ac", *cols]]
        df = df.merge(part, on="uniprot_ac", how="inner")
        if len(df) != n_proteome:
            raise AssertionError(
                f"{axis} join dropped rows for {species}: {len(df)} != {n_proteome}. "
                "Every axis must ship one complete row per protein in canonical order; "
                "run `python -m src.matrices` to find which axis regressed."
            )

    return df.reset_index(drop=True)
