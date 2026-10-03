"""Load the stage-03 localization predictions, and the vocabulary they share.

Stage 03 assigns every bacterial protein a subcellular compartment, **from sequence alone**. Two
predictors run, and each keeps its own file -- they answer different questions and are not merged:

    data/processed/localization/deeplocpro_{kpneumoniae,ecoli,saureus}.tsv
        uniprot_ac  localization  confidence          <- the compartment call

    data/processed/localization/tmbed_{kpneumoniae,ecoli,saureus}.tsv
        uniprot_ac  cytoplasmic_fraction  has_signal_peptide     <- the topology score

One row per protein in both, keyed on `uniprot_ac`, every protein present. `localization` is never
empty: DeepLocPro emits a softmax over six classes and always returns a call, so coverage is 100%
by construction and there is no `unknown` class. v1 shipped UniProt-curated labels only and left
59.6% of Kp / 48.6% of Ec at `unknown`, which silently dropped them from every shortlist.

Human is absent by construction: DeepLocPro is prokaryote-only.

The two files answer different questions
----------------------------------------
`localization` is a **compartment**, one of six labels. `cytoplasmic_fraction` is a **number in
[0, 1]** -- the share of residues TMbed places on the cytoplasmic side. They agree in the aggregate
(fraction falls monotonically as the compartment moves outward) but neither derives from the other,
and TMbed is not a localization classifier. Join them with `load()` when you want both.

No accessibility score is computed here. Turning a compartment plus a fraction into a single
"can Clp reach this" verdict is a modelling decision that belongs to whichever stage consumes it.

The Gram-positive remap
-----------------------
DeepLocPro takes a Gram group, and for `positive` it does NOT merely mask the two Gram-negative-only
classes -- `DeepLocPro/utils.py:remap_probabilities` **adds their probability into Extracellular**.
So *S. aureus* has four reachable classes, and its `extracellular` count absorbs whatever the model
wanted to call periplasmic. `evidence/deeplocpro_probabilities_<species>.tsv` holds the raw
un-remapped six-vector, which is the only place that effect is measurable. Read it before comparing
*S. aureus*'s `extracellular` share against either Gram-negative.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
LOCALIZATION_DIR = REPO_ROOT / "data" / "processed" / "localization"
EVIDENCE_DIR = LOCALIZATION_DIR / "evidence"
SCRATCH_DIR = LOCALIZATION_DIR / "scratch"
# The three bacteria. Human is not predictable here -- see the module docstring.
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# ---------------------------------------------------------------- vocabulary

# The canonical class set, ordered **inward -> outward** through the cell envelope, which is also
# the order Clp reachability falls in. Six labels and no more: unlike v1 there is no `unknown`
# (DeepLocPro always calls) and no side-less `membrane` (that was an artifact of parsing UniProt
# free text, which this stage does not do).
LOC_CLASSES: tuple[str, ...] = (
    "cytoplasm",
    "cytoplasmic_membrane",
    "periplasm",
    "outer_membrane",
    "cell_wall_surface",
    "extracellular",
)

# v1's NPG palette, carried over verbatim so figure colours cannot drift between the two versions.
LOC_CLASS_COLOR: dict[str, str] = {
    "cytoplasm": "#00A087",
    "cytoplasmic_membrane": "#F39B7F",
    "periplasm": "#E64B35",
    "outer_membrane": "#8491B4",
    "cell_wall_surface": "#7E6148",
    "extracellular": "#3C5488",
}

LOC_CLASS_ABBREV: dict[str, str] = {
    "cytoplasm": "Cyt",
    "cytoplasmic_membrane": "CM",
    "periplasm": "Peri",
    "outer_membrane": "OM",
    "cell_wall_surface": "CW",
    "extracellular": "Ext",
}

# DeepLocPro 1.0's own class order, verbatim from `DeepLocPro/deeplocpro.py:36`. The probability
# vector comes back in THIS order; do not reorder it without changing DLP_CANON with it.
DLP_LABELS: tuple[str, ...] = (
    "Cell wall & surface",
    "Extracellular",
    "Cytoplasmic",
    "Cytoplasmic Membrane",
    "Outer Membrane",
    "Periplasmic",
)

DEEPLOCPRO_MAP: dict[str, str] = {
    "cytoplasmic": "cytoplasm",
    "cytoplasmic membrane": "cytoplasmic_membrane",
    "periplasmic": "periplasm",
    "outer membrane": "outer_membrane",
    "cell wall & surface": "cell_wall_surface",
    "extracellular": "extracellular",
}

# Our labels, in DeepLocPro's vector order.
DLP_CANON: tuple[str, ...] = tuple(DEEPLOCPRO_MAP[lbl.lower()] for lbl in DLP_LABELS)

# Which Gram group each anchor belongs to. Kp and Ec are Gram-negative; Sa is Gram-positive and so
# has neither a periplasm nor an outer membrane.
GRAM_GROUP: dict[str, str] = {
    "kpneumoniae": "negative",
    "ecoli": "negative",
    "saureus": "positive",
}

# The classes a Gram-positive organism can actually occupy.
GRAM_POSITIVE_CLASSES: tuple[str, ...] = (
    "cytoplasm", "cytoplasmic_membrane", "cell_wall_surface", "extracellular",
)

# A Gram-negative outer-membrane barrel starts at 8 strands; one or two stray predicted strands are
# noise. Kept as a constant rather than a column -- the boolean is trivial to re-derive from
# `n_tm_strand` in evidence/tmbed_topology_<species>.tsv, and no consumer needs it yet.
BETA_BARREL_MIN_STRANDS = 8

# ESM-2 attention is quadratic in length. Localization signal is overwhelmingly N-terminal (signal
# peptides, TM topology), so the C-terminus is what gets cut. v1 measured only 2 Kp / 1 Ec proteins
# affected; the count is printed every run regardless.
DLP_MAX_LENGTH = 2000

# TMbed's per-residue alphabet at `--out-format 4` (`tmbed/tmbed.py:252`). Format 4 rather than v1's
# format 1: format 1 collapses `h` into `H` and `b` into `B`, discarding segment orientation, for no
# saving at all.
TMBED_ALPHABET: dict[str, str] = {
    "H": "transmembrane alpha-helix, inside -> outside",
    "h": "transmembrane alpha-helix, outside -> inside",
    "B": "transmembrane beta-strand, inside -> outside",
    "b": "transmembrane beta-strand, outside -> inside",
    "S": "signal peptide",
    "i": "not in a membrane, cytoplasmic side",
    "o": "not in a membrane, outside (periplasm / extracellular)",
}


def remap_gram_positive(probs: list[float]) -> list[float]:
    """Fold periplasm and outer-membrane probability into extracellular, DeepLocPro's own rule.

    Reproduces `DeepLocPro/utils.py:remap_probabilities` on our vector, which is in `DLP_LABELS`
    order. This is the whole of what `--group positive` does, and it is arithmetic on the output --
    not a different forward pass -- which is why stage 03 runs the model once, stores the raw
    vector, and applies this only when picking a label.
    """
    out = list(probs)
    out[1] = out[1] + out[4] + out[5]     # Extracellular += Outer Membrane + Periplasmic
    out[4] = 0.0
    out[5] = 0.0
    return out


def label_from_probs(probs: list[float], gram_group: str) -> tuple[str, float, float]:
    """(label, confidence, margin) from a raw DeepLocPro vector, applying the Gram remap first."""
    p = remap_gram_positive(probs) if gram_group in ("positive", "archaea") else list(probs)
    ranked = sorted(range(len(p)), key=lambda i: -p[i])
    return DLP_CANON[ranked[0]], p[ranked[0]], p[ranked[0]] - p[ranked[1]]


# ---------------------------------------------------------------- TMbed parsing

def _n_segments(labels: str, chars: str) -> int:
    """Count contiguous runs, counted **per letter** and summed.

    Per letter, not per character class: TMbed's Viterbi decoding gives a segment one orientation,
    so `HHHhhh` is two segments (one each way), which merged counting would report as one.
    """
    return sum(len([s for s in re.split(f"[^{c}]+", labels) if s]) for c in chars)


def tmbed_features(labels: str) -> dict:
    """Per-protein features from one TMbed label string.

    `cytoplasmic_fraction` and `has_signal_peptide` are the deliverable; the rest goes to
    `evidence/tmbed_topology_<species>.tsv` because the run is CPU-only hours and cannot be
    partially redone, so discarding a free number would be the expensive choice.
    """
    n = len(labels) or 1
    inside_runs = [s for s in re.split("[^i]+", labels) if s]
    return {
        "cytoplasmic_fraction": round(labels.count("i") / n, 4),
        "has_signal_peptide": "S" in labels,
        "n_tm_helix": _n_segments(labels, "Hh"),
        "n_tm_strand": _n_segments(labels, "Bb"),
        "cyto_longest_segment": max((len(r) for r in inside_runs), default=0),
        "tmbed_length": len(labels),
    }


def parse_tmbed(path: Path) -> dict[str, dict]:
    """Parse a TMbed 3-line prediction file into {accession: features}.

    The file is `>header` / sequence / labels, repeating. Headers are bare accessions because the
    query FASTA is built that way -- see `write_query_fasta` in scripts/localization/predict.py.
    """
    out: dict[str, dict] = {}
    lines = [ln.rstrip("\n") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    for i in range(0, len(lines) - 2, 3):
        header, seq, labels = lines[i], lines[i + 1], lines[i + 2]
        if not header.startswith(">"):
            continue
        if len(labels) != len(seq):
            raise ValueError(
                f"{path.name}: {header[1:]} has {len(seq)} residues but {len(labels)} labels -- "
                "the 3-line block is misaligned, do not trust this shard"
            )
        out[header[1:].strip()] = tmbed_features(labels)
    return out


def read_tmbed_labels(path: Path) -> dict[str, str]:
    """Parse a TMbed 3-line prediction file into {accession: raw label string}."""
    out: dict[str, str] = {}
    lines = [ln.rstrip("\n") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    for i in range(0, len(lines) - 2, 3):
        if lines[i].startswith(">"):
            out[lines[i][1:].strip()] = lines[i + 2]
    return out


# ---------------------------------------------------------------- loaders

def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")


def _path(directory: Path, name: str) -> Path:
    path = directory / name
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/localization/predict.py first")
    return path


def load_deeplocpro(species: str) -> pd.DataFrame:
    """One species' compartment calls: uniprot_ac, localization, confidence.

    **Evidence, not the deliverable** -- `load()` is the table this axis ships. Kept as its own
    file because the two predictors run and resume independently, so a `--only` run must still
    have somewhere to write.
    """
    _check(species)
    return pd.read_csv(_path(EVIDENCE_DIR, f"deeplocpro_{species}.tsv"), sep="\t",
                       dtype={"uniprot_ac": str, "localization": str, "confidence": float})


def load_tmbed(species: str) -> pd.DataFrame:
    """One species' topology score: uniprot_ac, cytoplasmic_fraction, has_signal_peptide.

    **Evidence, not the deliverable** -- same contract as `load_deeplocpro`.
    """
    _check(species)
    return pd.read_csv(_path(EVIDENCE_DIR, f"tmbed_{species}.tsv"), sep="\t",
                       dtype={"uniprot_ac": str, "cytoplasmic_fraction": float,
                              "has_signal_peptide": bool})


def load(species: str) -> pd.DataFrame:
    """`localization_<species>.tsv` -- the axis's single deliverable, complete and canonical.

        uniprot_ac  localization  cytoplasmic_fraction

    **`confidence` is NOT here** (owner's call, 2026-10-03) — byte-identical in
    `evidence/deeplocpro_<species>.tsv`, via `load_deeplocpro()`. **Know what goes with that.**
    DeepLocPro always returns a call and has no `unknown` class, so `localization` reads equally
    authoritative for every protein; `confidence` was the only column saying otherwise. **12–15% of
    calls sit below 0.7 and 2–4% below 0.5** — the winning class holding less than half the
    probability mass. Join it back before trusting a single label, and use
    `load_probabilities()` for the full un-remapped six-vector.

    **`has_signal_peptide` is NOT here either** (owner's call, 2026-10-03) — it ships byte-identically in
    `evidence/tmbed_<species>.tsv`, via `load_tmbed()`. What it said that `cytoplasmic_fraction`
    alone cannot is **why** a fraction is near zero: exported rather than membrane-buried. That
    distinction is mechanistically live — a secreted protein transits the cytoplasm unfolded and IS
    reachable by activated ClpP, where a membrane protein never does and is protected — so join it
    back before using localization as a degradability filter.

    **No `evidence` column, unlike the other axes.** Both predictors cover 100% of every proteome
    BY CONSTRUCTION, so the column was constant across all 13,020 proteins and said nothing;
    `merge.py` still exits non-zero if either track is silent. Do not re-add it as a constant.

    **Two predictors side by side, never reduced to one call.** DeepLocPro answers *which
    compartment*, TMbed *how much of the chain faces the cytoplasm*; neither derives from the
    other, and where they disagree the disagreement is the information -- TMbed corroborates
    `extracellular`, DeepLocPro's weakest class, from outside that model. **Prefer
    `cytoplasmic_fraction` over `localization` where a choice is forced.**

    **100% coverage is a property of the method, not evidence.** DeepLocPro always returns a call
    and has no `unknown` class, so a confident-looking label is not a measurement.

    **The Gram-positive trap.** On *S. aureus* the model runs in `positive` mode, which does not
    merely mask `periplasm` and `outer_membrane` -- it ADDS their probability mass into
    `extracellular`. Sa therefore has four reachable classes, and the absence of those two labels
    is STRUCTURAL, not missing data. `load_probabilities()` is the only place that moved mass can
    be measured.

    **E. coli K-12 is almost certainly in DeepLocPro's training set**, so any E. coli agreement
    number is a sanity check rather than validation; Kp and Sa are the honest test sets.

    No accessibility score is derived here -- that is a modelling decision for whichever stage
    consumes this, and v1's `clp_accessibility` ladder was consumed by nothing.
    """
    _check(species)
    return pd.read_csv(_path(LOCALIZATION_DIR, f"localization_{species}.tsv"), sep="\t")


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column added back."""
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_probabilities(species: str) -> pd.DataFrame:
    """The raw un-remapped six-class vector, plus margin and the truncation flag.

    For *S. aureus* this is the only place the Gram-positive remap's effect can be measured:
    `p_periplasm + p_outer_membrane` is exactly the mass that was moved into `extracellular`.
    """
    _check(species)
    return pd.read_csv(_path(EVIDENCE_DIR, f"deeplocpro_probabilities_{species}.tsv"), sep="\t")


def load_topology(species: str) -> pd.DataFrame:
    """TMbed's other features: n_tm_helix, n_tm_strand, cyto_longest_segment, tmbed_length."""
    _check(species)
    return pd.read_csv(_path(EVIDENCE_DIR, f"tmbed_topology_{species}.tsv"), sep="\t")


def load_labels(species: str) -> pd.DataFrame:
    """The raw per-residue TMbed label string, one row per protein.

    Kept so any future topology feature is derivable with zero recompute -- TMbed is CPU-only hours.
    """
    _check(species)
    return pd.read_csv(_path(EVIDENCE_DIR, f"tmbed_labels_{species}.tsv"), sep="\t", dtype=str)


def load_counts(species: str) -> pd.DataFrame:
    """Per-compartment counts for one species: all six classes, zeros included."""
    _check(species)
    return pd.read_csv(_path(EVIDENCE_DIR, f"deeplocpro_counts_{species}.tsv"), sep="\t")


def manifest(predictor: str = "deeplocpro") -> pd.DataFrame:
    """One row per species for the named predictor: counts, runtime, provenance."""
    if predictor not in ("deeplocpro", "tmbed"):
        raise ValueError(f"unknown predictor {predictor!r}; expected 'deeplocpro' or 'tmbed'")
    return pd.read_csv(_path(EVIDENCE_DIR, f"{predictor}_manifest.tsv"), sep="\t")
