"""Load the 2D projections of the stage-01 ESM-C embeddings.

`scripts/embeddings/projection.py` reduces each species' 1152-dim ESM-C matrix to two
coordinates, one row per protein:

    data/processed/embeddings/projection_{kpneumoniae,ecoli,saureus}.tsv
        uniprot_ac  tsne_x  tsne_y

Coverage is 100% by construction -- t-SNE returns a coordinate for every row it is given, and the
stage exits non-zero if a single protein goes missing.

The recipe, and why it is not up for re-derivation
--------------------------------------------------
openTSNE multiscale (perplexities 50/500), cosine, on PCA-50 of the z-scored embeddings, with
`dof=0.8`. v1 chose this over UMAP and PaCMAP in a 33-config sweep across four axes; the verdict is
tabulated in `docs/embeddings.md` (Part 2) and the original in
`legacy/docs/projection_exploration_log.md`. The sweep code was deleted, so those tables are the
record. `dof` is the load-bearing knob, not perplexity: 0.5-0.6 gives compact clusters flung too far
apart, 1.0 keeps them close but smears them, 0.8 is the balance.

Use `--method umap|pacmap` on the stage for a spot comparison; it writes to `evidence/ + scratch/` and never
touches the deliverable.

Coordinates are RELATIVE, and per-species
-----------------------------------------
The three species are three independent embeddings with no shared frame. `tsne_x = 40` in
*K. pneumoniae* has nothing to do with `tsne_x = 40` in *E. coli*: not the same axis, not the same
units, not even the same orientation. So:

- Plot one species per panel. `load_all()` adds a `species` column for convenience, but stacking the
  three into one scatter is meaningless, and shared `xlim`/`ylim` across panels is misleading.
- Never compute a distance, a neighbourhood or a cluster across species on these columns. Go back to
  `src/embeddings.py` and work in the 1152-dim space, which *is* shared.
- t-SNE distances are not metric even *within* a species: cluster membership and local neighbourhoods
  are meaningful, absolute distances and cluster sizes are not.

A cross-species-comparable frame would need something like v1's `01c` ESM Atlas coordinates (a single
global UMAP of the protein universe, keyed on sequence hash), which this stage does not produce.

No colouring is stored here
---------------------------
Just the coordinates. Colour by joining whatever you want to show -- every stage keys on
`uniprot_ac`, so it is a one-liner:

    proj = projections.load("kpneumoniae")
    cog  = function.load_cog("kpneumoniae")
    df   = proj.merge(cog[["uniprot_ac", "cog_category"]], on="uniprot_ac", how="left")

v1 persisted a `family` column (KMeans k=60 in a cosine UMAP-15d space) that its own shipped figure
then ignored in favour of a density fade. It is not reproduced: an arbitrary cluster id that no
downstream consumer used is not worth a column.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PROJECTION_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
EVIDENCE_DIR = PROJECTION_DIR / "evidence"
SCRATCH_DIR = PROJECTION_DIR / "scratch"
# The three bacteria. Human has no stage-01 embeddings, so there is nothing to project.
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# ---------------------------------------------------------------- the chosen configuration

METHOD = "opentsne"
PCA_COMPONENTS = 50  # PCA-50 measured indistinguishable from full 1152-d, and much faster
PERPLEXITIES = (50, 500)  # multiscale: local + global structure together
METRIC = "cosine"  # note: openTSNE's NN backends do NOT support metric="correlation"
DOF = 0.8  # the compact <-> detached knob; 0.7-0.8 only

# Optimisation schedule, both phases, verbatim from v1's `01b_esmc_projections.py`.
EARLY_ITER, EARLY_EXAGGERATION, EARLY_MOMENTUM = 250, 12, 0.5
LATE_ITER, LATE_MOMENTUM = 500, 0.8

# Alternatives, kept runnable via `--method` but never the deliverable. openTSNE is canonical.
COMPARISON_METHODS = ("umap", "pacmap")

# Trustworthiness floor for the stage's control. Calibrated, not guessed -- measured on *S. aureus*:
# random 2D coordinates score 0.4996, PCA-2 alone scores 0.7821, and the chosen recipe scores 0.9752.
# 0.90 sits clear of both degenerate baselines and well under what a healthy run achieves, so it
# fires on a genuinely broken configuration without tripping on run-to-run variation.
MIN_TRUSTWORTHINESS = 0.90
TRUSTWORTHINESS_K = 10

COORD_COLS = ("tsne_x", "tsne_y")
_DTYPE = {"uniprot_ac": str, "tsne_x": float, "tsne_y": float}


def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")


def _path(directory: Path, name: str) -> Path:
    path = directory / name
    if not path.exists():
        raise FileNotFoundError(
            f"{path} -- run scripts/embeddings/projection.py first"
        )
    return path


def load(species: str) -> pd.DataFrame:
    """One species' map: `uniprot_ac`, `tsne_x`, `tsne_y`. One row per protein."""
    _check(species)
    return pd.read_csv(
        _path(PROJECTION_DIR, f"projection_{species}.tsv"), sep="\t", dtype=_DTYPE
    )


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All three maps stacked, with a `species` column.

    For per-species faceting and joins only. The coordinate frames are independent, so a single
    scatter of this frame -- or a shared axis range across facets -- says nothing real. See the
    module docstring.
    """
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_method(species: str, method: str) -> pd.DataFrame:
    """A comparison projection from `--method umap|pacmap`, out of `evidence/ + scratch/`.

    Columns are `uniprot_ac`, `<method>_x`, `<method>_y` -- deliberately not `tsne_*`, so a
    comparison map can never be mistaken for the deliverable in a merged frame.
    """
    _check(species)
    if method not in COMPARISON_METHODS:
        raise ValueError(
            f"unknown comparison method {method!r}; expected one of {COMPARISON_METHODS}"
        )
    return pd.read_csv(
        _path(EVIDENCE_DIR, f"projection_{method}_{species}.tsv"),
        sep="\t",
        dtype={"uniprot_ac": str, f"{method}_x": float, f"{method}_y": float},
    )


def coords_for(species: str, accessions: list[str]) -> tuple[list[str], np.ndarray]:
    """Coordinates for a given accession list, in that order. Returns `(found, xy)`.

    Accessions with no coordinate are dropped from `found` rather than filled with zeros or NaN --
    mirroring `src.embeddings.vectors_for`, and for the same reason: **(0, 0) is a real position on
    the map**, near the dense centre, so a zero-filled row would silently plant a protein in the
    middle of the proteome.
    """
    df = load(species).set_index("uniprot_ac")
    found = [a for a in accessions if a in df.index]
    if not found:
        return [], np.zeros((0, 2), dtype=float)
    return found, df.loc[found, list(COORD_COLS)].to_numpy(dtype=float)


def manifest() -> pd.DataFrame:
    """The projection run manifest: per species n, method, params, KL, trustworthiness, seconds."""
    return pd.read_csv(_path(EVIDENCE_DIR, "projection_manifest.tsv"), sep="\t")
