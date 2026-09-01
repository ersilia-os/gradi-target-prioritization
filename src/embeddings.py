"""Load the stage-01 ESM-C embeddings.

Stage 01 writes one NPZ per species under `data/processed/01_embeddings/`. Use these helpers rather
than `np.load` directly: `accessions` is an **object array**, so a raw load needs
`allow_pickle=True`, and v1 open-coded the load-plus-row-index dance in two separate scripts.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
EMBEDDING_DIR = REPO_ROOT / "data" / "processed" / "01_embeddings"
ACCESSORY_DIR = EMBEDDING_DIR / "accessory"

MODEL_ID = "esmc_600m"
EMBED_DIM = 1152


def _path(species: str) -> Path:
    path = EMBEDDING_DIR / f"embeddings_{species}.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/01_embeddings.py first")
    return path


def load(species: str) -> tuple[np.ndarray, np.ndarray]:
    """`(accessions, matrix)` — a `(n,)` array of UniProt ACs and an `(n, 1152)` float32 matrix.

    The two are row-aligned: `matrix[i]` is the embedding of `accessions[i]`.
    """
    z = np.load(_path(species), allow_pickle=True)
    return z["accessions"].astype(str), z["embeddings"].astype(np.float32)


def load_frame(species: str) -> pd.DataFrame:
    """The embedding matrix as a DataFrame indexed by `uniprot_ac`, columns `e0..e1151`."""
    accs, mat = load(species)
    return pd.DataFrame(mat, index=pd.Index(accs, name="uniprot_ac"),
                        columns=[f"e{i}" for i in range(mat.shape[1])])


def load_lookup(species: str) -> dict[str, np.ndarray]:
    """`{accession: vector}` — for code that needs a handful of proteins, not the whole matrix."""
    accs, mat = load(species)
    return {a: mat[i] for i, a in enumerate(accs)}


def vectors_for(species: str, accessions: list[str]) -> tuple[list[str], np.ndarray]:
    """Rows for a given accession list, in that order. Returns `(found, matrix)`.

    Accessions with no embedding are dropped from `found` rather than filled with zeros or NaN — a
    silent zero row would look like a real position in embedding space.
    """
    accs, mat = load(species)
    index = {a: i for i, a in enumerate(accs)}
    found = [a for a in accessions if a in index]
    if not found:
        return [], np.zeros((0, mat.shape[1]), dtype=np.float32)
    return found, mat[[index[a] for a in found]]


def metadata(species: str) -> dict:
    """What produced this file: model, pooling, dim."""
    z = np.load(_path(species), allow_pickle=True)
    return {"model": str(z["model"]), "pooling": str(z["pooling"]), "dim": int(z["dim"]),
            "n": int(z["embeddings"].shape[0])}


def manifest() -> pd.DataFrame:
    """The stage-01 run manifest: per species n, dim, residues, device, seconds, sha256."""
    path = ACCESSORY_DIR / "manifest.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/01_embeddings.py first")
    return pd.read_csv(path, sep="\t")
