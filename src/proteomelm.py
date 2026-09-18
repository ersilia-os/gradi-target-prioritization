"""Load the ProteomeLM contextualised embeddings written by scripts/embeddings/proteomelm.py.

This mirrors `src/embeddings.py` deliberately -- same function names, same return shapes -- so code
that reads one can read the other by swapping the import. It exists for the same reason: the
`accessions` object array needs `allow_pickle=True`, and open-coding the load-plus-row-index dance
at every call site is how v1 accumulated three subtly different versions of it.

What these embeddings ARE, so they are not confused with stage 01's:

  * stage 01 (`src/embeddings.py`) is **per-protein ESM-C 600M**, 1152-d, mean over residues with
    BOS/EOS stripped. Each protein is independent of every other.
  * these are **ProteomeLM-L layer 8**, z-scored genome-wide, computed from ESM-C vectors pooled
    ProteomeLM's way (BOS/EOS INCLUDED -- the two conventions agree at cosine 0.999972, measured).
    They are **context-dependent**: the same protein embedded inside its full proteome versus a half
    proteome differs at cosine ~0.965. So a row here is only meaningful together with the proteome
    it was computed in, and rows from different runs are not interchangeable.

`group_embeds_mode` records the functional encoding used. It is `self` -- each protein's own ESM-C
vector -- which is the authors' own released inference default (`use_odb=False`), not the OrthoDB
hierarchy vectors used during training. Check it before comparing against anything.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTEOMELM_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
EVIDENCE_DIR = PROTEOMELM_DIR / "evidence"
SCRATCH_DIR = PROTEOMELM_DIR / "scratch"
MODEL_REPO = "Bitbol-Lab/ProteomeLM-L"
LAYER = 8
EMBED_DIM = 1152
GROUP_EMBEDS_MODE = "self"
LABELS = ("kpneumoniae", "ecoli", "saureus")


def _path(label: str) -> Path:
    path = PROTEOMELM_DIR / f"proteomelm_{label}.npz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} -- run scripts/embeddings/proteomelm.py first")
    return path


def load(label: str) -> tuple[np.ndarray, np.ndarray]:
    """(accessions (n,), embeddings (n, dim) float32), row-aligned."""
    z = np.load(_path(label), allow_pickle=True)
    return z["accessions"].astype(str), z["embeddings"].astype(np.float32)


def load_frame(label: str) -> pd.DataFrame:
    """Indexed by uniprot_ac, columns p0..p<dim-1>."""
    accs, mat = load(label)
    return pd.DataFrame(mat, index=pd.Index(accs, name="uniprot_ac"),
                        columns=[f"p{i}" for i in range(mat.shape[1])])


def load_lookup(label: str) -> dict[str, np.ndarray]:
    accs, mat = load(label)
    return {a: mat[i] for i, a in enumerate(accs)}


def vectors_for(label: str, accessions: list[str]) -> tuple[list[str], np.ndarray]:
    """Rows in the REQUESTED order. Missing accessions are DROPPED, never zero-filled -- a zero row
    is a real position in embedding space, and silently inventing one is how a downstream model ends
    up trained on fabricated points."""
    lut = load_lookup(label)
    found = [a for a in accessions if a in lut]
    if not found:
        return [], np.zeros((0, EMBED_DIM), dtype=np.float32)
    return found, np.vstack([lut[a] for a in found]).astype(np.float32)


def metadata(label: str) -> dict:
    z = np.load(_path(label), allow_pickle=True)
    return {
        "model": str(z["model"]),
        "layer": int(z["layer"]),
        "n_layers": int(z["n_layers"]),
        "dim": int(z["dim"]),
        "group_embeds_mode": str(z["group_embeds_mode"]),
        "esmc_pooling": str(z["esmc_pooling"]),
        "proteome_id": str(z["proteome_id"]),
        "n": int(z["embeddings"].shape[0]),
    }


def manifest() -> pd.DataFrame:
    path = EVIDENCE_DIR / "proteomelm_manifest.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/embeddings/proteomelm.py first")
    return pd.read_csv(path, sep="\t")
