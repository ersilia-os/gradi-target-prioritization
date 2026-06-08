"""Project the ESM-C proteome embeddings to a single publication-quality 2D map.

openTSNE (multiscale, cosine) on PCA-50 of the z-scored ESM-C 600M embeddings, with
`dof=0.8` -- the tail weight that keeps families distinct but not flung apart. Points are
colored by KMeans families derived in a cosine UMAP space, so colors form coherent islands.

This is the chosen configuration after a broad sweep; see
docs/projection_exploration_log.md for what else was tried and why this won.

Output (PNG only): output/plots/kp_esmc600m_projection.png
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import colorcet as cc
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openTSNE import TSNEEmbedding, affinity, initialization
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

REPO_ROOT = Path(__file__).resolve().parents[1]
EMB_PATH = REPO_ROOT / "data" / "processed" / "embeddings" / "kp_esmc600m_embeddings.npz"
OUT_PATH = REPO_ROOT / "output" / "plots" / "kp_esmc600m_projection.png"
PALETTE = cc.glasbey_dark


def project(pca, seed):
    """openTSNE multiscale (50/500), cosine, dof=0.8 -- the chosen configuration."""
    aff = affinity.Multiscale(pca, perplexities=[50, 500], metric="cosine",
                              n_jobs=-1, random_state=seed)
    init = initialization.pca(pca, random_state=seed)
    emb = TSNEEmbedding(init, aff, dof=0.8, n_jobs=-1, random_state=seed)
    emb = emb.optimize(n_iter=250, exaggeration=12, momentum=0.5)
    emb = emb.optimize(n_iter=500, momentum=0.8)
    return np.asarray(emb)


def coherent_colors(z, seed, k):
    """KMeans(k) families in a cosine UMAP-15d space -> coherent per-point colors."""
    import umap
    space = umap.UMAP(n_components=15, n_neighbors=15, min_dist=0.0,
                      metric="cosine", random_state=seed).fit_transform(z)
    labels = KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(space)
    return np.array([matplotlib.colors.to_rgba(PALETTE[int(l) % len(PALETTE)]) for l in labels])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=60, help="KMeans families for coloring")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    data = np.load(EMB_PATH, allow_pickle=True)
    emb = data["embeddings"].astype(np.float32)
    print(f"Loaded {emb.shape[0]} x {emb.shape[1]} embeddings")

    z = StandardScaler().fit_transform(emb)
    pca = PCA(n_components=min(50, *emb.shape), random_state=args.seed).fit_transform(z)

    colors = coherent_colors(z, args.seed, args.k)
    print("Projecting (openTSNE multiscale, cosine, dof=0.8) ...")
    xy = project(pca, args.seed)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.scatter(xy[:, 0], xy[:, 1], s=7, c=colors, linewidths=0, alpha=0.8)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
