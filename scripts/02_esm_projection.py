"""Comprehensive openTSNE evaluation for the ESM-C proteome embeddings (PNG only).

openTSNE is the canonical projection method for this project. This script evaluates it
thoroughly across three axes, holding the coloring fixed so panels are comparable:

  A. INPUT / PCA-prior : full z-scored 1152-d, PCA-50/100/200, and L2-normalized full
  B. DISTANCE metric   : cosine, euclidean, correlation, manhattan
  C. t-SNE STRUCTURE   : single perplexity vs multiscale; tail weight `dof`; exaggeration

Coloring: KMeans(k) in a cosine UMAP space (UMAP used ONLY to derive stable colors), the
same labels on every panel -> coherent colored families, no gray noise.

Outputs (PNG only): output/plots/projection_sweep/<config>.png (300 dpi, full-res) and a
large many-panel contact sheet output/plots/projection_sweep_contactsheet.png to explore.
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import traceback
from pathlib import Path

import colorcet as cc
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, normalize

REPO_ROOT = Path(__file__).resolve().parents[1]
EMB_PATH = REPO_ROOT / "data" / "processed" / "embeddings" / "kp_esmc600m_embeddings.npz"
OUT_DIR = REPO_ROOT / "output" / "plots" / "projection_sweep"
SHEET = REPO_ROOT / "output" / "plots" / "projection_sweep_contactsheet.png"
PALETTE = cc.glasbey_dark


# ----------------------------- openTSNE -----------------------------
def otsne(X, perplexity=50, metric="cosine", exaggeration=None, dof=1.0, seed=0):
    from openTSNE import TSNE
    return np.asarray(
        TSNE(perplexity=perplexity, metric=metric, initialization="pca", dof=dof,
             exaggeration=exaggeration, n_jobs=-1, random_state=seed, verbose=False).fit(X)
    )


def otsne_multiscale(X, perplexities=(50, 500), metric="cosine", exaggeration=None, dof=1.0, seed=0):
    from openTSNE import TSNEEmbedding, affinity, initialization
    aff = affinity.Multiscale(X, perplexities=list(perplexities), metric=metric,
                              n_jobs=-1, random_state=seed)
    init = initialization.pca(X, random_state=seed)
    emb = TSNEEmbedding(init, aff, dof=dof, n_jobs=-1, random_state=seed)
    emb = emb.optimize(n_iter=250, exaggeration=12, momentum=0.5)
    emb = emb.optimize(n_iter=500, exaggeration=exaggeration, momentum=0.8)
    return np.asarray(emb)


# ----------------------------- coloring -----------------------------
def coherent_labels(z, seed, k):
    """Stable spatially-coherent families: KMeans(k) in a cosine UMAP-15d space."""
    import umap
    space = umap.UMAP(n_components=15, n_neighbors=15, min_dist=0.0,
                      metric="cosine", random_state=seed).fit_transform(z)
    return KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(space)


# ----------------------------- plotting -----------------------------
def style_ax(ax, xy, colors, title, point_size):
    ax.scatter(xy[:, 0], xy[:, 1], s=point_size, c=colors, linewidths=0, alpha=0.8)
    ax.set_title(title, fontsize=8)
    ax.set_aspect("equal")
    ax.axis("off")


def save_individual(name, xy, colors):
    fig, ax = plt.subplots(figsize=(11, 11))
    style_ax(ax, xy, colors, "", point_size=6)
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_configs(reps, seed):
    """Comprehensive openTSNE matrix across input/PCA, distance metric, and t-SNE structure."""
    cfg = {}
    # ---- Block A: representation x distance metric (single perplexity 50, dof 0.8) ----
    metrics = ["cosine", "euclidean", "manhattan"]  # robust across openTSNE NN backends
    for rep_name, X in reps.items():
        for m in metrics:
            cfg[f"A_{rep_name}_{m}"] = (lambda X=X, m=m: otsne(X, 50, m, dof=0.8, seed=seed))

    # ---- Block B: t-SNE structure on pca50 + cosine ----
    p50 = reps["pca50"]
    for perp in (15, 30, 50, 100):
        cfg[f"B_p{perp}"] = (lambda perp=perp: otsne(p50, perp, "cosine", dof=1.0, seed=seed))
    cfg["B_ms30_200"] = lambda: otsne_multiscale(p50, (30, 200), "cosine", seed=seed)
    cfg["B_ms50_500"] = lambda: otsne_multiscale(p50, (50, 500), "cosine", seed=seed)
    cfg["B_ms30_300_1000"] = lambda: otsne_multiscale(p50, (30, 300, 1000), "cosine", seed=seed)
    for d in (0.5, 0.7, 0.9):
        cfg[f"B_ms50_500_dof{d}"] = (lambda d=d: otsne_multiscale(p50, (50, 500), "cosine", dof=d, seed=seed))
    for e in (1.5, 2.0):
        cfg[f"B_ms50_500_exag{e}"] = (lambda e=e: otsne_multiscale(p50, (50, 500), "cosine", exaggeration=e, dof=1.0, seed=seed))

    # ---- Block C: compact-but-cohesive sweet spot (pca50 + cosine) ----
    cfg["C_p50_dof0.7_exag1.5"] = lambda: otsne(p50, 50, "cosine", exaggeration=1.5, dof=0.7, seed=seed)
    cfg["C_p30_dof0.7_exag2"] = lambda: otsne(p50, 30, "cosine", exaggeration=2.0, dof=0.7, seed=seed)
    cfg["C_ms50_500_dof0.7_exag1.5"] = lambda: otsne_multiscale(p50, (50, 500), "cosine", exaggeration=1.5, dof=0.7, seed=seed)
    cfg["C_ms50_500_dof0.8"] = lambda: otsne_multiscale(p50, (50, 500), "cosine", dof=0.8, seed=seed)
    cfg["C_ms100_800_dof0.8"] = lambda: otsne_multiscale(p50, (100, 800), "cosine", dof=0.8, seed=seed)
    cfg["C_ms50_500_dof0.8_exag1.5"] = lambda: otsne_multiscale(p50, (50, 500), "cosine", exaggeration=1.5, dof=0.8, seed=seed)
    return cfg


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=60, help="KMeans families for coloring")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", type=str, default="")
    args = ap.parse_args()

    data = np.load(EMB_PATH, allow_pickle=True)
    emb = data["embeddings"].astype(np.float32)
    if args.limit:
        emb = emb[: args.limit]
    print(f"Loaded {emb.shape[0]} x {emb.shape[1]} embeddings")

    z = StandardScaler().fit_transform(emb)
    pca200 = PCA(n_components=min(200, *emb.shape), random_state=args.seed).fit_transform(z)
    reps = {
        "full": z,
        "pca50": pca200[:, :50],
        "pca100": pca200[:, :100],
        "pca200": pca200,
        "l2full": normalize(z),
    }

    labels = coherent_labels(z, args.seed, args.k)
    print(f"KMeans coloring: {len(set(labels))} families")
    colors = np.array([matplotlib.colors.to_rgba(PALETTE[int(l) % len(PALETTE)]) for l in labels])

    configs = build_configs(reps, args.seed)
    if args.only:
        wanted = {n.strip() for n in args.only.split(",")}
        configs = {n: f for n, f in configs.items() if n in wanted}
    print(f"Total configs: {len(configs)}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results: dict[str, np.ndarray] = {}
    for name, fn in configs.items():
        try:
            print(f"  running {name} ...", flush=True)
            xy = np.asarray(fn(), dtype=float)
            results[name] = xy
            save_individual(name, xy, colors)
        except Exception:
            print(f"  FAILED {name}:\n{traceback.format_exc()}")

    names = sorted(results)  # groups A_/B_/C_ together
    ncol = 6
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.2 * nrow), squeeze=False)
    for ax, name in zip(axes.ravel(), names):
        style_ax(ax, results[name], colors, name, point_size=1.5)
    for ax in axes.ravel()[len(names):]:
        ax.axis("off")
    fig.suptitle(
        f"openTSNE comprehensive sweep — ESM-C 600M proteome (n={emb.shape[0]}, {len(set(labels))} families)\n"
        f"A: input/PCA x metric   B: t-SNE structure   C: compact-cohesive sweet spot",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(SHEET, dpi=140, facecolor="white")
    plt.close(fig)
    print(f"\nWrote {len(results)} PNGs to {OUT_DIR.relative_to(REPO_ROOT)}/ and contact sheet {SHEET.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
