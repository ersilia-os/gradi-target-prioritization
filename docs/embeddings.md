# Stage 01 — ESM-C 600M embeddings

`scripts/embeddings/esmc.py` · helpers `src/embeddings.py` · outputs `data/processed/embeddings/`

One mean-pooled 1152-dimensional vector per protein, for the three bacteria — and, in Part 2, the 2D
map of that space.

| species | proteins | residues | max length |
|---|---|---|---|
| *K. pneumoniae* HS11286 | 5,728 | 1,646,261 | 3,163 |
| *E. coli* K-12 | 4,403 | 1,354,442 | 2,339 |
| *S. aureus* NCTC 8325 | 2,889 | 797,397 | **9,535** |
| **total** | **13,020** | **3,798,100** | |

# Part 1 — the embeddings

## The model, and why not a bigger one

`esmc_600m` (`ESMC_600M_202412`), 1152-dim, mean-pooled over residues with BOS/EOS stripped.

**This is the largest ESM-C with downloadable weights.** Verified in the installed `esm` 3.2.1:
`pretrained.py` registers only `ESMC_300M` and `ESMC_600M`, while `esmc-6b` appears solely in
`utils/constants/models.py` as a name check for the hosted **Forge API** — which requires a token,
sends every sequence to a third party, has rate limits and cost, and cannot be reproduced offline.
600M is also what v1 used, so the existing ESM Atlas coordinates are on the same footing.

Weights come from the local HuggingFace cache
(`~/.cache/huggingface/hub/models--EvolutionaryScale--esmc-600m-2024-12`), already present from v1.

## Scope: no human, no ProtT5 (for now)

**Human is out.** It would have added 20,416 proteins and 11.4M residues — 75% of all residues — plus
every pathological length: titin 34,350 aa, MUC16 14,507, MUC3B 13,477, three proteins over 10k, 36
over 5k, 465 over 2k. Nothing in the code excludes it; `--species human` adds it back, but expect the
long tail to dominate the runtime and to need a chunking or truncation policy that this stage does not
currently have.

**ProtT5 is deferred, and the reason is not obvious.** UniProt *does* publish precomputed ProtT5-XL-U50
per-protein embeddings in HDF5, at
`ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/`. But only for
**8 proteomes**:

| | |
|---|---|
| *E. coli* K-12 `UP000000625_83333/per-protein.h5` | 10.1 MB ✅ |
| human `UP000005640_9606/per-protein.h5` | 47.4 MB ✅ |
| **Kp HS11286** | **404** ❌ |
| **Sa NCTC 8325** | **404** ❌ |
| all Swiss-Prot `uniprot_sprot/per-protein.h5` | 1.3 GB |

The others are mouse, rat, *C. elegans*, Arabidopsis and SARS-CoV-2. The Swiss-Prot bulk file is no
help for the anchor: HS11286 has **7 reviewed entries**, so it would fill 7 of 5,728 proteins.

So "ProtT5 as available in UniProt" cannot cover the anchor organism. Covering Kp and *S. aureus*
means running ProtT5-XL-U50 locally (~2.25 GB, and v1's TMbed experience says it is CPU-only on Apple
Silicon because TMbed gates GPU on `torch.cuda.is_available()`). That is a separate decision.

**If it is ever taken up, there is a free validation opportunity**: download UniProt's h5 for E. coli,
run our own ProtT5 on the same sequences, and require cosine ≈ 1. Same model, same input — if it does
not reproduce, our setup is wrong, and we would know before trusting Kp numbers.

Also note UniProt excludes sequences over 12k residues from its own embeddings, "due to limitation of
GPU memory".

## Measured throughput (MPS, torch 2.12)

Taken before any code was written, one sequence at a time:

```
   100 aa  0.30 s      257 aa  0.15 s      500 aa  0.26 s     994 aa  0.59 s
  2035 aa  1.98 s     3163 aa  3.24 s     9535 aa 22.0  s
```

≈1,700 aa/s in the typical 250–1,000 aa band, falling to ~975 aa/s by 3k as O(L²) attention bites.
**A full run is 42–63 minutes.**

**No length cap or chunking is needed.** Only 9 of the 13,020 proteins exceed 2,000 aa, and the
longest — `Q2FYJ6`, a 9,535 aa *S. aureus* ECM-binding protein — embeds in 22 s without trouble.

## Sharding and resume — the one thing changed from v1

v1's `01a_esmc_embeddings.py` accumulated every vector in a Python list and wrote a single
`np.savez_compressed` at the end, so **a killed run lost everything**. For a 45-minute job that is the
wrong trade, so this stage copies the pattern v1's own `09d_tmbed.py` got right:

- proteins sorted by accession for determinism, cut into shards of 250 (`--shard-size`);
- each shard written to `accessory/shards/embeddings_<species>_<shard_size>_<NNNN>.npz` via **atomic
  `tmp.rename()`**, so a partial file is never cached;
- **the shard size is part of the filename** — changing `--shard-size` must not silently reuse
  mismatched shards, a trap v1 documented;
- a restart skips completed shards. Verified: removing one shard of three and re-running recomputes
  only that shard and reproduces a **byte-identical** matrix.

## Output

```
data/processed/embeddings/
  embeddings_kpneumoniae.npz   embeddings_ecoli.npz   embeddings_saureus.npz
  accessory/
    manifest.tsv          species, device, n, dim, residues, seconds, sha256
    shards/               the resumable per-shard cache
```

Each NPZ:

```
accessions   object   (n,)        row-aligned to `embeddings`
embeddings   float32  (n, 1152)   mean-pooled, BOS/EOS stripped
model        scalar   "esmc_600m"
pooling      scalar   "mean"
dim          scalar   1152
```

**Load through `src/embeddings.py`**, not by hand — `accessions` is an object array, so a raw
`np.load` needs `allow_pickle=True`, and v1 open-coded the load-plus-row-index dance separately in
`07d` and `10f`:

| helper | returns |
|---|---|
| `load(species)` | `(accessions, matrix)`, row-aligned |
| `load_frame(species)` | DataFrame indexed by `uniprot_ac`, columns `e0..e1151` |
| `load_lookup(species)` | `{accession: vector}` for a handful of proteins |
| `vectors_for(species, accs)` | rows for a given list, in order; missing accessions are **dropped, not zero-filled** |
| `metadata(species)` / `manifest()` | what produced the file |

`vectors_for` drops rather than zero-fills deliberately: a zero row is a real position in embedding
space and would look like a legitimate protein to anything downstream.

## Running it

```bash
python scripts/embeddings/esmc.py                                  # the three bacteria
python scripts/embeddings/esmc.py --limit 20 --species saureus     # smoke test, seconds
python scripts/embeddings/esmc.py --species kpneumoniae --refresh  # discard shards, recompute
```

Flags: `--species` · `--shard-size 250` · `--device auto|cuda|mps|cpu` · `--limit N` · `--refresh` ·
`-q`.

**Env: `gradi`, not `gradi-loc`.** `gradi-loc` carries `fair-esm`, which claims the same top-level
`esm` package name as the EvolutionaryScale `esm` this needs. The script says so on ImportError rather
than failing obscurely.

`PYTORCH_ENABLE_MPS_FALLBACK=1` is set **before torch is imported** — order-sensitive, and silent if
you get it wrong.

### Guards

- **A failed protein is named, never swallowed.** v1 collected failures into a `skipped` list and
  carried on quietly; here the run prints them and **exits non-zero** if any species is incomplete,
  because a missing row would corrupt every downstream axis.
- **Row count is asserted against stage 00** — the summary prints `n` beside `expected` and flags a
  mismatch.

## Run log

### 2026-09-01 — first full run

All three species, **exact row counts, 0 skipped**. MPS, torch 2.12, `esmc_600m`.

| species | n | dim | residues | wall clock |
|---|---|---|---|---|
| *K. pneumoniae* | 5,728 | 1152 | 1,646,261 | 14.5 min |
| *E. coli* | 4,403 | 1152 | 1,354,442 | ~12 min |
| *S. aureus* | 2,889 | 1152 | 797,397 | **21.4 min** |

Verified: shapes `(n, 1152)` float32, all values finite, **row-aligned** to the stage-00 tables (not
merely the same set), no zero rows.

`ecoli` has 4,362 distinct rows out of 4,403. That is correct, not a failure: the proteome contains
41 duplicate sequences, which must produce identical embeddings.

Biology sanity check on `clpP`:

```
Kp vs Ec   0.995     same protein, close relatives
Kp vs Sa   0.972     same protein, distant
clpP vs gyrB (Kp)   0.920     different proteins, same organism
```

Ordering is right. Absolute cosines run high across the board, which is normal for mean-pooled PLM
embeddings — only the relative ordering is meaningful.

### One measurement that did not survive contact with the real run

Planning timed the longest protein, `Q2FYJ6` (9,535 aa), at **22 s** in isolation and concluded no
length cap was needed. In the actual run its shard took **886.8 s** against ~37 s for every other
shard — that one protein cost roughly 850 extra seconds, **14 of S. aureus's 21.4 minutes**.

The isolated timing was on a freshly-loaded model with clean MPS memory; in a warm process with
accumulated allocations the same sequence is ~40× slower, almost certainly falling back to CPU under
memory pressure (`PYTORCH_ENABLE_MPS_FALLBACK=1` makes that silent by design).

Consequences worth carrying: **time a long sequence in a warm process, not a fresh one**, and if human
is ever added — 465 proteins over 2,000 aa, titin at 34,350 — this effect will dominate the run
entirely, so a chunking or truncation policy becomes mandatory rather than optional.

### The shard cache is tied to the filename

Renaming outputs from `<species>.npz` to `embeddings_<species>.npz` also renamed the shard pattern,
which **orphaned all 41 cached shards** and started a full recompute. That is the "shard size is in
the filename" guard working in reverse: it prevents silent reuse of mismatched shards, and any naming
change invalidates the cache. The fix was to rename the shard files to match; a re-run then reported
`23 cached / 18 cached, 0 to compute`, confirming the cache was honoured.

# Part 2 — 2D projections

`scripts/embeddings/projection.py` · figures `scripts/plots/projection.py` ·
helpers `src/projections.py` · outputs `data/processed/embeddings/projection_<species>.tsv`

Two coordinates per protein, one map per species, so the proteome can be looked at and so later
stages can paint their own scores onto a fixed frame.

| species | proteins | trustworthiness (k=10) | KL | time |
|---|---|---|---|---|
| *K. pneumoniae* | 5,728 | **0.9790** | 1.274 | 66 s |
| *E. coli* | 4,403 | **0.9735** | 1.195 | 45 s |
| *S. aureus* | 2,889 | **0.9752** | 0.923 | 29 s |
| **total** | **13,020** | | | **2.3 min** |

## The recipe is inherited, not re-derived

**openTSNE multiscale (perplexities 50/500), cosine, on PCA-50 of the z-scored embeddings, with
`dof=0.8`.** v1 chose this in a 33-config sweep and the verdict is reproduced here because the sweep
code was deleted (`scripts/02_esm_projection.py`, commit `a9a7939`) — this table and
`legacy/docs/projection_exploration_log.md` are the whole record. **Do not re-run the sweep.**

| axis | tested | verdict |
|---|---|---|
| input | full 1152-d, PCA-50/100/200, L2-normalised | **PCA-50 ≈ full-dim** — the map is essentially identical, so PCA-50 is the correct fast default. This tested and dismissed an earlier worry that pre-PCA was distorting the result |
| metric | cosine, euclidean, manhattan | **cosine** clearly cleanest for transformer embeddings; euclidean and manhattan diffuse. `correlation` is **unsupported** by openTSNE's NN backends |
| `dof` | 0.5 → 1.0 | **0.8.** 0.5–0.6 gives compact clusters flung too far apart ("satellites too detached"); →1.0 keeps them close but diffuse. 0.7–0.8 is the balance |
| perplexity | single 15/30/50/100; multiscale (30,200)/(50,500)/(30,300,1000) | **multiscale [50, 500]** — local and global structure together. Very large global perplexity (200/1500) spreads the map out |
| exaggeration | 1.5, 2.0 | **off in the final phase.** Above ~1.5 with heavy tails it produces filament artifacts |
| method | openTSNE, UMAP, PaCMAP | **openTSNE.** UMAP: continents cohere but families blend. PaCMAP: pretty continental gradient, but splits into 2 masses with filament tails |

`dof` is the load-bearing knob, not perplexity — worth knowing before tuning anything else.

`--method umap|pacmap` re-runs the comparison into `accessory/` and never touches the deliverable.
It exists because v1's comparison code was lost; the columns are named `<method>_x`/`<method>_y` so a
comparison map cannot be mistaken for the real one in a merged frame.

## Two things from v1 that were deliberately not carried over

1. **v1's shipped figure never used the `family` column it wrote.** `01b` computed KMeans(k=60) in a
   cosine UMAP-15d space, persisted it, then coloured the PNG by `gaussian_kde` density instead. The
   exploration log's "colored by KMeans families" describes the sweep, not the script that shipped.
   No `family` column here: an arbitrary cluster id that no downstream consumer used is not worth
   persisting, and real labels (COG, localization) are one `merge` away.
2. **v1 fed t-SNE the PCA-50 but its colour-clustering the full 1,152-d matrix** — slower, and two
   different spaces answering one question. Here `prepare()` produces one matrix that the projection
   *and* the control both use.

## The control: trustworthiness

A t-SNE that diverged, or was handed the wrong metric, still returns a plausible-looking array of
finite numbers. So the guard is not a shape check but
`sklearn.manifold.trustworthiness(pca, xy, n_neighbors=10, metric="cosine")` — the fraction of each
protein's 10 nearest neighbours in the input space that are still among its neighbours on the map.

The floor is **calibrated, not guessed**. Measured on *S. aureus*:

| configuration | trustworthiness (k=10, cosine) |
|---|---|
| random 2D coordinates | 0.4996 |
| PCA-2 (first two components only) | 0.7821 |
| **the chosen recipe** | **0.9752** |
| floor (`MIN_TRUSTWORTHINESS`) | 0.90 |

0.90 sits clear of both degenerate baselines and well below what a healthy run achieves, so it fires
on a genuinely broken configuration without tripping on run-to-run variation. All three species
landed at 0.973–0.979.

KL divergence is recorded but **not** enforced — there is no meaningful absolute floor for it; it is
useful only as a drift signal between runs.

## Coordinates are relative, and per-species

The three species are three independent embeddings with **no shared frame**. `tsne_x = 40` in
*K. pneumoniae* has nothing to do with `tsne_x = 40` in *E. coli*: not the same axis, not the same
units, not even the same orientation. Consequences:

- **One species per panel.** `load_all()` adds a `species` column for joining and faceting, but a
  single scatter of all three — or a shared axis range across facets — implies a comparability that
  does not exist.
- **Never compute a distance, neighbourhood or cluster across species on these columns.** Go back to
  `src/embeddings.py` and work in the 1152-dim space, which *is* shared.
- Even *within* a species, t-SNE distances are not metric: cluster membership and local
  neighbourhoods are meaningful, absolute distances and cluster sizes are not.

A cross-species-comparable frame would need something like v1's `01c` ESM Atlas coordinates — a
single global UMAP of the protein universe keyed on sequence hash, fetched from the Biohub API. That
is a network-fetching stage of its own and is not produced here.

## Output

```
data/processed/embeddings/
    projection_kpneumoniae.tsv      280 kB
    projection_ecoli.tsv            198 kB
    projection_saureus.tsv          130 kB
    accessory/
        projection_manifest.tsv     params, KL, trustworthiness, timings
        projection_<method>_<sp>.tsv   only with --method umap|pacmap
        smoke_*                        only with --limit
```

```
uniprot_ac      str     UniProt accession, the canonical key
tsne_x          float   first  t-SNE coordinate (relative, per-species)
tsne_y          float   second t-SNE coordinate (relative, per-species)
```

Three columns and nothing else. **No colour column** — join what you want to show on `uniprot_ac`;
every stage keys on it, so it is a one-liner. Coverage is 100% by construction: t-SNE returns a
coordinate for every row it is given, and the stage exits non-zero if one goes missing.

| helper | returns |
|---|---|
| `load(species)` | `uniprot_ac`, `tsne_x`, `tsne_y` |
| `load_all(species=SPECIES)` | all three stacked, with a `species` column |
| `load_method(species, method)` | a `--method umap\|pacmap` comparison map from `accessory/` |
| `coords_for(species, accessions)` | `(found, xy)` in the given order; missing **dropped**, not zero-filled |
| `manifest()` | the run manifest |

`coords_for` drops rather than zero-fills for the same reason `embeddings.vectors_for` does:
**(0, 0) is a real position on the map**, near the dense centre, so a zero-filled row would silently
plant a protein in the middle of the proteome.

## Figures

`scripts/plots/projection.py` writes `output/plots/embeddings/`:
`projection_<species>.png` (one large map each) and `projection_all.png` (a 3-panel sheet).

Uncoloured, by request. Each map is a single stylia hue faded by 2D point density, so dense family
cores read saturated and sparse outliers pale — an aesthetic for an overplotted scatter, not an
encoding of any variable, which is why there is no legend. Kp `cobalt` and Ec `crimson` are v1's
hues, carried over so figures do not drift between versions; Sa gets `turquoise`.

Because nothing is coloured, the plotter is a thin consumer of the TSV: a COG-category or
degradability-scored version is a `merge` plus a `c=`, with no change to the stage.

Everything goes through stylia — `create_figure` / `label` / `save_figure`, `slide` format and
`ersilia` style, with no hand-set marker or font sizes. The first version of this plotter used raw
`plt.subplots` / `fig.savefig` / `ax.set_title` and a hard-coded font size; that is against house
style, and was corrected when stage 03's figures were written.

## Running it

```bash
~/miniconda3/envs/gradi/bin/python scripts/embeddings/projection.py
~/miniconda3/envs/gradi/bin/python scripts/embeddings/projection.py --species saureus
~/miniconda3/envs/gradi/bin/python scripts/embeddings/projection.py --method umap
~/miniconda3/envs/gradi/bin/python scripts/plots/projection.py
```

Flags: `--species` · `--method {opentsne,umap,pacmap}` · `--pca 50` · `--perplexities 50 500` ·
`--dof 0.8` · `--seed 0` · `--limit N` (smoke test → `accessory/smoke_*`) · `--refresh` ·
`--dry-run` · `-q/--quiet`.

Everything needed is already in `gradi` (openTSNE 1.0.4, umap-learn 0.5.12, pacmap 0.9.1,
scikit-learn 1.9.0, stylia 1.0.1). Nothing to install, no foreign env, no process boundary.

### Guards

- **Row count** — the output must have exactly as many rows as the species has embeddings. A
  projection that silently drops proteins would quietly corrupt every downstream landscape. Exits
  non-zero.
- **Finite coordinates** — no `NaN`/`inf`. Exits non-zero.
- **Accession order** — the output accessions must match the embeddings row-for-row, not merely as a
  set. Exits non-zero.
- **Trustworthiness ≥ 0.90**, per species. Exits non-zero, and the message quotes the random-2D and
  PCA-2 baselines so the number can be judged.
- Comparison (`--method`) and smoke (`--limit`) runs are **not** guarded — they are exploratory and
  never become the deliverable.

## Run log

### 2026-09-02 — first full run

All three species, `--seed 0`, default recipe. **2.3 min total** on CPU (no GPU path needed):

```
kpneumoniae     5728 x 1152  ->  PCA-50   trustworthiness 0.9790   KL 1.274   66s
ecoli           4403 x 1152  ->  PCA-50   trustworthiness 0.9735   KL 1.195   45s
saureus         2889 x 1152  ->  PCA-50   trustworthiness 0.9752   KL 0.923   29s
```

Verified: row counts, finite coordinates and accession order all match the embeddings for all three.

The maps meet v1's stated goal — "clusters distinct but not flung apart, with an even point
distribution". Each species resolves into one large continent with visible internal substructure plus
three to five detached satellites; no single blob, no filament spray. The three share that
architecture, which is what one would expect of three bacterial proteomes seen through the same model.

`prepare()` reproduced the planning measurement on *S. aureus* to four decimal places
(0.9752), which is the reassuring version of the `n_jobs=-1` caveat below: structure is stable even
though floats are not guaranteed bit-identical.

### v1's method verdict, now with a number attached

v1 chose openTSNE over UMAP and PaCMAP on looks — "families blend", "splits into 2 masses with
filament tails". Running all three through `--method` on *S. aureus* puts a measurement on that
judgement for the first time:

| method | trustworthiness (k=10, cosine) | time |
|---|---|---|
| **openTSNE** (multiscale 50/500, `dof=0.8`) | **0.9752** | 29 s |
| UMAP (`n_neighbors=15, min_dist=0.1`) | 0.9599 | 13 s |
| PaCMAP (defaults) | 0.9355 | 2 s |

The ranking is exactly the one v1 reached by eye, which is a mild but real vindication of a decision
that was, until now, only defensible by pointing at a picture. All three clear the 0.90 floor, so
this is a preference among good options, not a rescue from a bad one — and the differences are small
enough that they would not, on their own, have settled the choice. Speed runs the other way; the
canonical map is the slowest of the three and that is a fine trade at ~1 min per species.

Both comparison maps are kept under `accessory/projection_{umap,pacmap}_saureus.tsv` as the evidence.

### The COG join, and what 100% actually means

Spot-checking the `uniprot_ac` key against stage 02: the merge matches **100% of rows** for all three
species, and **79.05%** of Kp rows carry a non-empty `cog_category` — reproducing stage 02's
documented 79.1% Kp coverage exactly.

Both numbers are correct and they are different questions. `cog_<species>.tsv` carries a row for
every protein and an **empty string** for the unclassified, so *joining* is total while *annotation*
is not. Worth knowing before reading any coverage percentage off a merged frame: check for the empty
string, not for `NaN`.

### A cached re-run destroyed the manifest, once

First version of the caching path: a species whose table already exists is served from cache and its
manifest row is written from what *this* run knows — which, having fitted nothing, is no
trustworthiness, no KL and `seconds=0`. So simply running the stage twice replaced the real
measurements with blanks. The tables were fine; only the record of how good they were was gone, which
is the sort of loss that is invisible until someone needs the number.

`merge_manifest()` now carries a cached species' previous row forward and only lets a genuinely
recomputed species overwrite its own row. Two consequences worth keeping:

- the **summary prints the stored numbers** on a cached run, tagged `(cached)`, rather than blanks;
- the **trustworthiness floor is enforced against the stored score**, so a cached re-run still fails
  if the map on disk is bad — the guard does not quietly pass because this run measured nothing.

Rebuilding with `--refresh` reproduced all three scores exactly (0.9790 / 0.9735 / 0.9752, identical
KL), which is the reassuring counterpart to the `n_jobs=-1` caveat below.

### `n_jobs=-1` means the floats are not reproducible, only the structure

The seed threads into `PCA`, `affinity.Multiscale`, `initialization.pca` and `TSNEEmbedding`, so a
re-run on the same machine reproduces the map. But openTSNE's parallel reduction is not
order-deterministic, so **coordinates are not bit-identical across machines** even at a fixed seed.

This is why the control is trustworthiness rather than a checksum of the coordinates: a hash would
fail spuriously on another machine while telling you nothing about whether the map is any good.

### openTSNE silently clamps too-large perplexities

The smoke test (`--limit 300`) printed `Perplexity value 500 is too high. Using perplexity 99.67
instead` — openTSNE clamps to roughly `(n-1)/3` rather than failing. Harmless here, and it means the
default `[50, 500]` is safe at *S. aureus*'s n=2,889 (clamp point 962). But it does mean **a small
`--limit` is not running the production recipe**, so a smoke test validates plumbing, not the map.

### The exploration log's own "Reproduce" block no longer runs

`legacy/docs/projection_exploration_log.md` ends with a command referencing
`scripts/02_esm_projection.py` and `data/processed/embeddings/kp_esmc600m_embeddings.npz`. All of
those paths moved or were deleted. This section of this file supersedes it; the log is worth reading
for the *reasoning*, not for the commands.
