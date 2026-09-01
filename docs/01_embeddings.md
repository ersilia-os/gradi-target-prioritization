# Stage 01 — ESM-C 600M embeddings

`scripts/01_embeddings.py` · helpers `src/embeddings.py` · outputs `data/processed/01_embeddings/`

One mean-pooled 1152-dimensional vector per protein, for the three bacteria.

| species | proteins | residues | max length |
|---|---|---|---|
| *K. pneumoniae* HS11286 | 5,728 | 1,646,261 | 3,163 |
| *E. coli* K-12 | 4,403 | 1,354,442 | 2,339 |
| *S. aureus* NCTC 8325 | 2,889 | 797,397 | **9,535** |
| **total** | **13,020** | **3,798,100** | |

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
data/processed/01_embeddings/
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
python scripts/01_embeddings.py                                  # the three bacteria
python scripts/01_embeddings.py --limit 20 --species saureus     # smoke test, seconds
python scripts/01_embeddings.py --species kpneumoniae --refresh  # discard shards, recompute
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
