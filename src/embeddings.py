"""Load the stage-01 ESM-C embeddings.

Stage 01 writes one NPZ per species under `data/processed/embeddings/`. Use these helpers rather
than `np.load` directly: `accessions` is an **object array**, so a raw load needs
`allow_pickle=True`, and v1 open-coded the load-plus-row-index dance in two separate scripts.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
EMBEDDING_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
EVIDENCE_DIR = EMBEDDING_DIR / "evidence"
SCRATCH_DIR = EMBEDDING_DIR / "scratch"
MODEL_ID = "esmc_600m"
EMBED_DIM = 1152


def _path(species: str) -> Path:
    path = EMBEDDING_DIR / f"embeddings_{species}.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/embeddings/esmc.py first")
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
    path = EVIDENCE_DIR / "esmc_manifest.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/embeddings/esmc.py first")
    return pd.read_csv(path, sep="\t")


# --------------------------------------------------------------- ProtT5 (a stage-01 sibling)
#
# A SECOND language model over the same proteomes, written by `scripts/embeddings/prott5.py`.
# Kept beside the ESM-C helpers rather than merged into them: the two are different models with
# different widths (1152 vs 1024) and different pooling, and silently returning one where the caller
# expected the other is exactly the kind of mix-up the explicit names prevent.
#
# Unlike ESM-C, these are VALIDATED AGAINST AN EXTERNAL REFERENCE: UniProt publishes its own ProtT5
# vectors for E. coli (one of only 7 proteomes), and ours reproduce them at median cosine 1.000000 /
# worst 0.999994 over a 300-protein length-weighted sample. S. aureus and K. pneumoniae are NOT
# published by UniProt, which is why they are computed locally and why that control matters.

PROTT5_MODEL_ID = "Rostlab/prot_t5_xl_half_uniref50-enc"
PROTT5_DIM = 1024


def _prott5_path(species: str) -> Path:
    path = EMBEDDING_DIR / f"prott5_{species}.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/embeddings/prott5.py first")
    return path


def load_prott5(species: str) -> tuple[np.ndarray, np.ndarray]:
    """`(accessions, matrix)` — a `(n,)` array of UniProt ACs and an `(n, 1024)` float32 matrix."""
    z = np.load(_prott5_path(species), allow_pickle=True)
    return z["accessions"].astype(str), z["embeddings"].astype(np.float32)


def load_prott5_frame(species: str) -> pd.DataFrame:
    """The ProtT5 matrix as a DataFrame indexed by `uniprot_ac`, columns `t0..t1023`.

    `t` rather than `e`, so a frame carrying both models cannot collide.
    """
    accs, mat = load_prott5(species)
    return pd.DataFrame(mat, index=pd.Index(accs, name="uniprot_ac"),
                        columns=[f"t{i}" for i in range(mat.shape[1])])


def prott5_metadata(species: str) -> dict:
    """What produced this file, including which proteins were WINDOWED rather than single-pass.

    `n_chunked` is not cosmetic: a windowed protein's vector is the length-weighted mean of window
    means, so each window saw only its own context. One S. aureus protein (Q2FYJ6, 9,535 aa) is in
    that class. Check this before treating a long protein's embedding as equivalent to the rest.
    """
    z = np.load(_prott5_path(species), allow_pickle=True)
    chunked = z["chunked"] if "chunked" in z else np.zeros(z["embeddings"].shape[0], dtype=bool)
    return {"model": str(z["model"]), "pooling": str(z["pooling"]), "dim": int(z["dim"]),
            "n": int(z["embeddings"].shape[0]), "n_chunked": int(np.asarray(chunked).sum()),
            "max_len": int(z["max_len"]) if "max_len" in z else None}


def prott5_control() -> pd.DataFrame:
    """Per-protein cosine of our E. coli ProtT5 against UniProt's published vectors."""
    path = EVIDENCE_DIR / "prott5_control.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/embeddings/prott5.py first")
    return pd.read_csv(path, sep="\t")
