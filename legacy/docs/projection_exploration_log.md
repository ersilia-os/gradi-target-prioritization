# ESM-C embedding projection — exploration log

Record of how the 2D projection of the HS11286 ESM-C 600M embeddings
(`data/processed/embeddings/kp_esmc600m_embeddings.npz`, 5,728 × 1152) was chosen.
The production script `scripts/02_esm_projection.py` now produces **only the final figure**;
this file documents the exploration behind that choice so it isn't lost.

## Goal

A beautiful, publication-quality 2D map that reveals the many protein families in the
proteome — clusters that are **distinct but not flung apart**, with an even point distribution.

## Final choice

**openTSNE multiscale, cosine, on PCA-50 of z-scored embeddings, `dof=0.8`**, colored by
KMeans(k=60) families computed in a cosine UMAP space. (Sweep name: `C_ms50_500_dof0.8`.)

- Method: openTSNE (user's preferred / canonical method).
- Input: StandardScaler → PCA-50. (PCA-50 ≈ full 1152-d — see below.)
- Affinities: `Multiscale(perplexities=[50, 500])`, `metric="cosine"`.
- Tail weight: `dof=0.8` (the cohesion knob — see below).
- Output: `output/plots/kp_esmc600m_projection.png`.

## What was tried, and what we learned

### Methods
- **openTSNE** — chosen. Multiscale affinities (local+global) + cosine + tunable `dof`.
- **UMAP** — good for cohesive continents but families blend; kept only as a comparison.
- **PaCMAP** — pretty continental gradient but split into 2 masses with filament tails; not used.

### Input / PCA-prior (axis A)
- Compared full z-scored 1152-d, PCA-50/100/200, and L2-normalized full.
- **PCA-50 ≈ full-dim**: the map is essentially identical, so PCA-50 is the correct fast
  default. (Validated the earlier doubt about pre-PCA — it is not the problem.)

### Distance metric (axis B)
- Compared cosine, euclidean, manhattan (correlation unsupported by openTSNE NN backends).
- **cosine** is clearly cleanest for these transformer embeddings; euclidean/manhattan diffuse.

### t-SNE structure (axis C) — the compact ↔ detached tension
- `dof` (tail weight) is the key knob:
  - **low `dof`** (e.g. 0.5–0.6, heavy tails) → tight, compact clusters but flung far apart
    ("satellites too detached").
  - **`dof`→1.0** → clusters stay close but diffuse/spread ("too spread out").
  - **`dof≈0.7–0.8`** → the balance: compact-ish families that stay connected. **0.8 chosen.**
- **Multiscale** perplexities `[50, 500]` give local + global structure together.
- Very large global perplexity (e.g. 200/1500) → too spread out. Avoid.
- **Exaggeration** > ~1.5 combined with heavy tails → ugly filament artifacts. Avoid for the
  final; mild exaggeration only tightens, at `dof≈1.0`.

### Coloring (made the maps readable)
Evolved through several attempts:
1. KMeans on PCA-50 → salt-and-pepper (global clusters don't match the map's local structure).
2. HDBSCAN on PCA-50 → only ~2 clusters (structure is non-linear).
3. HDBSCAN on a cosine UMAP space → coherent but ~30% gray noise washed out the core.
4. **KMeans(k≈60) in a cosine UMAP-15d space** → coherent colored families, every point
   colored, no gray. **This is the final coloring.** (UMAP is used here ONLY to derive
   stable colors, not for the projection itself.)

## Artifacts from the exploration (eosvc-tracked, gitignored)
- Final figure: `output/plots/kp_esmc600m_projection.png`
- (Historical) full 33-config contact sheet + per-config PNGs were produced under
  `output/plots/projection_sweep*`; regenerate by checking out an earlier version of
  `scripts/02_esm_projection.py` if a full sweep is ever needed again.

## Reproduce the final figure
```bash
conda activate gradi
python scripts/02_esm_projection.py        # writes output/plots/kp_esmc600m_projection.png
```
