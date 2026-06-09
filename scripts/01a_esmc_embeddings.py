"""Compute per-protein ESM-C 600M embeddings for the HS11286 proteome.

ESM-C (Cambrian, EvolutionaryScale) is the current-generation protein language model
for *representations*; the 600M variant emits 1152-dimensional embeddings and rivals
ESM2-3B quality. This is the sidecar embedding called for in docs/01_task_agnostic.md
§1.5 (originally specced as ESM2-650M; ESM-C 600M is the suggested upgrade).

For each protein we take the mean over residue embeddings (excluding the BOS/EOS special
tokens) to get one 1152-d vector. Output is a single NPZ keyed by UniProt accession.

Embeddings are a reusable transformed derivative of the proteome, so they are written to
data/processed/ (eosvc-tracked, not Git), keyed by the canonical UniProt accession.

Usage:
    python scripts/01_esmc_embeddings.py                 # all proteins, auto device
    python scripts/01_esmc_embeddings.py --limit 4       # smoke test (first 4)
    python scripts/01_esmc_embeddings.py --device cpu     # force a device
"""

from __future__ import annotations

# Let any op unsupported on Apple MPS fall back to CPU instead of crashing.
# Must be set before torch is imported.
import os

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

MODEL_ID = "esmc_600m"
EMBED_DIM = 1152
POOLING = "mean"

REPO_ROOT = Path(__file__).resolve().parents[1]
# organism -> (proteome file stem, output prefix)
ORGANISMS = {
    "kpneumoniae": ("UP000007841_HS11286", "kp"),
    "ecoli": ("UP000000625_EcoliK12", "ec"),
}


def pick_device(requested: str) -> str:
    if requested != "auto":
        return requested
    if torch.cuda.is_available():
        return "cuda"
    if (
        getattr(torch.backends, "mps", None) is not None
        and torch.backends.mps.is_available()
    ):
        return "mps"
    return "cpu"


def load_client(device: str):
    """Load ESM-C 600M; give an actionable hint if the gated download fails."""
    try:
        from esm.models.esmc import ESMC
    except ImportError as exc:
        sys.exit(
            f"Could not import the ESM-C model ({exc}).\n"
            "Install dependencies with:\n"
            "    python -m pip install -r requirements.txt\n"
            "(installs torch + esm + httpx)."
        )
    try:
        client = ESMC.from_pretrained(MODEL_ID).to(device).eval()
    except Exception as exc:  # noqa: BLE001 - surface a friendly hint, then re-raise
        msg = str(exc).lower()
        if any(
            k in msg for k in ("401", "403", "gated", "token", "authenticate", "login")
        ):
            sys.exit(
                f"Failed to download ESM-C weights ({MODEL_ID}).\n"
                "The HuggingFace repo may be gated. Accept the license on the model page, then:\n"
                "    huggingface-cli login      # or: export HF_TOKEN=...\n"
                f"Original error: {exc}"
            )
        raise
    return client


def embed_sequence(client, seq: str) -> np.ndarray:
    """Mean-pooled 1152-d embedding for one sequence (BOS/EOS excluded)."""
    from esm.sdk.api import ESMProtein, LogitsConfig

    with torch.no_grad():
        protein_tensor = client.encode(ESMProtein(sequence=seq))
        out = client.logits(
            protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
        )
        emb = out.embeddings[0]  # (L+2, 1152): BOS, residues..., EOS
        vec = emb[1:-1].mean(dim=0)  # mean over residues only
    return vec.float().cpu().numpy()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    ap.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    ap.add_argument("--limit", type=int, default=0, help="embed only first N (0 = all)")
    ap.add_argument("--out", type=Path, default=None, help="override output NPZ path")
    args = ap.parse_args()

    stem, prefix = ORGANISMS[args.organism]
    tsv_path = REPO_ROOT / "data" / "raw" / args.organism / "proteome" / f"{stem}.tsv"
    out_path = args.out or (
        REPO_ROOT
        / "data"
        / "processed"
        / args.organism
        / "embeddings"
        / f"{prefix}_esmc600m_embeddings.npz"
    )

    device = pick_device(args.device)
    print(f"Device: {device}")

    df = pd.read_csv(tsv_path, sep="\t")
    if args.limit:
        df = df.head(args.limit)
    print(f"Proteins to embed: {len(df)} (from {tsv_path.relative_to(REPO_ROOT)})")

    print(f"Loading model {MODEL_ID} ...")
    client = load_client(device)

    accessions: list[str] = []
    vectors: list[np.ndarray] = []
    skipped: list[str] = []
    for entry, seq in tqdm(
        zip(df["Entry"], df["Sequence"]), total=len(df), desc="embedding"
    ):
        try:
            vectors.append(embed_sequence(client, seq))
            accessions.append(entry)
        except Exception as exc:  # noqa: BLE001 - skip a bad protein, keep the run alive
            skipped.append(entry)
            print(f"  skipped {entry}: {exc}")

    embeddings = np.vstack(vectors).astype(np.float32)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_path,
        accessions=np.array(accessions, dtype=object),
        embeddings=embeddings,
        model=np.array(MODEL_ID),
        pooling=np.array(POOLING),
        dim=np.array(EMBED_DIM),
    )
    print(
        f"Wrote {embeddings.shape[0]} x {embeddings.shape[1]} embeddings "
        f"to {out_path.relative_to(REPO_ROOT)}"
    )
    if skipped:
        print(f"Skipped {len(skipped)} proteins: {', '.join(skipped)}")


if __name__ == "__main__":
    main()
