# Embeddings — ESM-C, ProtT5, ProteomeLM, and the 2D map

`scripts/embeddings/` · helpers `src/embeddings.py`, `src/proteomelm.py`, `src/projections.py` ·
outputs `data/processed/embeddings/`

Three project-wide embedding matrices over the same three bacterial proteomes, plus the 2D map of
the first one:

| part | what | script | file | dim |
|---|---|---|---|---|
| **1** | ESM-C 600M, per-protein | `embeddings/esmc.py` | `embeddings_<sp>.npz` | 1152 |
| **2** | ProtT5-XL-U50, per-protein | `embeddings/prott5.py` | `prott5_<sp>.npz` | 1024 |
| **3** | ProteomeLM-L layer 8, **contextualised** | `embeddings/proteomelm.py` | `proteomelm_<sp>[_<mode>].npz` | 1152 |
| **4** | t-SNE of Part 1 | `embeddings/projection.py` | `projection_<sp>.tsv` | 2 |

All four are complete and in canonical row order (`src/matrices.py`). Parts 1–3 are three different
answers to the same question, and *External models* in CLAUDE.md is the rule that keeps them
separable: none of them is substituted into a tool that was fitted on another.

| species | proteins | residues | max length |
|---|---|---|---|
| *K. pneumoniae* HS11286 | 5,728 | 1,646,261 | 3,163 |
| *E. coli* K-12 | 4,403 | 1,354,442 | 2,339 |
| *S. aureus* NCTC 8325 | 2,889 | 797,397 | **9,535** |
| **total** | **13,020** | **3,798,100** | |

# Part 1 — ESM-C 600M

`scripts/embeddings/esmc.py` · helpers `src/embeddings.py` · outputs
`data/processed/embeddings/embeddings_<species>.npz`

## The model, and why not a bigger one

`esmc_600m` (`ESMC_600M_202412`), 1152-dim, mean-pooled over residues with BOS/EOS stripped.

**That stripping is a convention, not a fact about the model, and Part 3 does it the other way.**
ProteomeLM pools over non-pad tokens, so BOS/EOS are *inside* its mean; the two conventions agree
at cosine **0.999970 median / 0.998291 worst** over the same proteins. Negligible in size, still a
different representation — which is why `embeddings/proteomelm.py` computes its own ESM-C rather
than reading these files. A tool's own convention wins; see Part 3 and *External models*.

**This is the largest ESM-C with downloadable weights.** Verified in the installed `esm` 3.2.1:
`pretrained.py` registers only `ESMC_300M` and `ESMC_600M`, while `esmc-6b` appears solely in
`utils/constants/models.py` as a name check for the hosted **Forge API** — which requires a token,
sends every sequence to a third party, has rate limits and cost, and cannot be reproduced offline.
600M is also what v1 used, so the existing ESM Atlas coordinates are on the same footing.

Weights come from the local HuggingFace cache
(`~/.cache/huggingface/hub/models--EvolutionaryScale--esmc-600m-2024-12`), already present from v1.

## Scope: no human

**Human is out.** It would have added 20,416 proteins and 11.4M residues — 75% of all residues — plus
every pathological length: titin 34,350 aa, MUC16 14,507, MUC3B 13,477, three proteins over 10k, 36
over 5k, 465 over 2k. Nothing in the code excludes it; `--species human` adds it back, but expect the
long tail to dominate the runtime and to need a chunking or truncation policy that this stage does not
currently have — see the `Q2FYJ6` measurement in the run log for how badly one long sequence behaves
in a warm process.

**ProtT5 is no longer deferred — it has its own part (Part 2).** The reason it needed one is still
worth keeping, because it is the reason a local run was unavoidable. UniProt *does* publish
precomputed ProtT5-XL-U50 per-protein embeddings in HDF5, at
`ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/` — but only for
**7 organism proteomes, plus an all-Swiss-Prot file**:

| | |
|---|---|
| *E. coli* K-12 `UP000000625_83333/per-protein.h5` | 10.1 MB ✅ |
| human `UP000005640_9606/per-protein.h5` | 47.4 MB ✅ |
| **Kp HS11286** | **404** ❌ |
| **Sa NCTC 8325** | **404** ❌ |
| all Swiss-Prot `uniprot_sprot/per-protein.h5` (not a proteome) | 1.3 GB |

The remaining five organisms are mouse, rat, *C. elegans*, Arabidopsis and SARS-CoV-2. The
Swiss-Prot bulk file is no help for the anchor: HS11286 has **7 reviewed entries**, so it would fill
7 of 5,728 proteins.

So "ProtT5 as available in UniProt" cannot cover the anchor organism — but *E. coli* MG1655 being on
that list is also the free validation opportunity, and Part 2 takes it: same model, same sequences,
so cosine must be ≈ 1 or our setup is wrong.

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
- each shard written to `scratch/shards_esmc/embeddings_<species>_<shard_size>_<NNNN>.npz` via
  **atomic `tmp.rename()`**, so a partial file is never cached;
- **the shard size is part of the filename** — changing `--shard-size` must not silently reuse
  mismatched shards, a trap v1 documented;
- a restart skips completed shards. Verified: removing one shard of three and re-running recomputes
  only that shard and reproduces a **byte-identical** matrix.

## Output

```
data/processed/embeddings/
  embeddings_kpneumoniae.npz   embeddings_ecoli.npz   embeddings_saureus.npz
  evidence/
    esmc_manifest.tsv     species, device, n, n_expected, dim, residues, seconds, skipped, sha256
  scratch/
    shards_esmc/          the resumable per-shard cache (safe to delete; a re-run rebuilds it)
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

# Part 2 — ProtT5-XL-U50

`scripts/embeddings/prott5.py` · worker `scripts/embeddings/workers/prott5.py` (env `gradi-loc`) ·
helpers `src/embeddings.py` · outputs `data/processed/embeddings/prott5_<species>.npz`

The second project-wide matrix: **1024-d float32, one vector per protein**, written beside the
ESM-C files over the same three proteomes. A sibling of Part 1, not a stage of its own.

| species | n | dim | windowed | file |
|---|---|---|---|---|
| *K. pneumoniae* | 5,728 | 1024 | 0 | 13.9 MB |
| *E. coli* | 4,403 | 1024 | 0 | 10.6 MB |
| *S. aureus* | 2,889 | 1024 | **1** | 7.0 MB |

## The control is the point, not an extra

*E. coli* MG1655 is one of the seven proteomes UniProt publishes ProtT5 vectors for (see *Scope*
above) **and it is our exact anchor**. Same model over the same sequences must therefore reproduce
the same vectors — a real, falsifiable check on a local inference pipeline, which an embedding
otherwise never gets. **Nothing is written if it fails**, because a plausible-looking vector that
differs from the published one is worse than none.

Scored on a 300-protein sample **weighted to the longest** (`--control-n`, half the sample is the
longest proteins in the proteome, half random), because that is where a pooling, truncation or
windowing bug actually shows; agreement on a 250-residue protein is nearly free.

| | |
|---|---|
| sample | 300 E. coli proteins, **18–2,339 aa** (single-pass agreement is established only to the longest protein checked) |
| median cosine | **1.000000** |
| worst cosine | **0.999994** |
| above the 0.99 floor | **100%** |

Per-protein cosines: `evidence/prott5_control.tsv` (300 rows), loadable with
`src.embeddings.prott5_control()`. The floor `COSINE_FLOOR = 0.99` is calibrated in this script's
own run, not inherited — same model, same sequences, so agreement should be essentially exact.

## It settled how residues are pooled

UniProt's README does not say. ProtT5's tokenizer appends an **EOS and no BOS**, so "mean over
residues" is ambiguous until you state whether EOS is in the mean. Both were embedded and scored
against the reference rather than argued about:

| pooling | median cosine vs UniProt | worst |
|---|---|---|
| **`mean_no_eos`** (the default) | **1.000000** | **0.999994** |
| `mean_with_eos` | — | 0.998837 |

`--compare-pooling` re-runs both and overwrites `evidence/prott5_pooling.tsv`, which is the record;
a single-pooling run leaves only the chosen row in it, so read the `pooling` column before quoting
the file. Residues are space-separated and **U/Z/O/B are mapped to `X`** — ProtT5's vocabulary has
no token for them, and selenocysteine is real here (*E. coli* `fdnG`, `fdhF`), so skipping that
would silently route rare residues to `<unk>`.

## Nothing is downloaded but a ~10 MB reference h5

The weights are already on disk: **TMbed ships `Rostlab/prot_t5_xl_half_uniref50-enc`** (1024-d,
fp16, encoder-only, 2.4 GB) inside its own package in the `gradi-loc` env. So inference runs
**across a process boundary** into that env — `scripts/embeddings/workers/prott5.py`, with
`GRADI_LOC_BIN` overriding the interpreter — and the only fetch is UniProt's reference h5, whose
byte count is checked against `Content-Length` rather than trusted from an HTTP 200.

The boundary is mandatory for a second reason: **ProtT5 needs `transformers==4.44.2`**. 5.x routes
`T5Tokenizer` through the tiktoken converter and dies with a spurious *"`tiktoken` is required to
read a `tiktoken` file"*, which reads like a missing dependency and is not one. `gradi` runs 4.48.1
and must keep doing so.

## This is the half, encoder-only build — which is why it must not be fed to someone else's tool

`prot_t5_xl_uniref50` (full) and `prot_t5_xl_half_uniref50-enc` are **different weights**, even
though both are called "ProtT5-XL-U50". Ours is the half, encoder-only one.

These vectors are valid for **our own** models — ProtT5 competed fairly in the degradability head
comparison, where it was the *best* features under TabPFN and the *worst* under a forest, which is
exactly the kind of interaction a shared, honestly-generated feature set is supposed to expose.
They are **never** substituted into an external predictor fitted on something else; SAFPred
(`docs/function.md`) is the worked rejection, and *External models* in CLAUDE.md is the rule.

## Sharding, and the one windowed protein

Shards of 250 into `scratch/shards_prott5_<species>/shard_250_<NNNNN>.npz`, resumable, shard size
in the filename so a changed `--shard-size` cannot reuse mismatched shards — the Part 1 pattern,
for the same reason (the encoder on MPS has stalled in this repo before).

Sequences above `--max-len 4096` are embedded in **non-overlapping windows** and pooled as the
length-weighted mean of window means, which *is* the mean over all residues; the windows differ
from a single pass only in that each sees its own context. Every such protein is **flagged in the
output** (`chunked`), never silently approximated. Measured: exactly **one** protein in the three
proteomes is windowed — *S. aureus* `Q2FYJ6`, 9,535 aa, the same protein that dominated Part 1's
runtime. Kp's longest is 3,163 aa and E. coli's 2,339 aa, so both are 0.

## Output

```
data/processed/embeddings/
  prott5_kpneumoniae.npz   prott5_ecoli.npz   prott5_saureus.npz
  evidence/
    prott5_control.tsv                  per-protein cosine vs UniProt -- the evidence
    prott5_pooling.tsv                  the pooling comparison
    prott5_ecoli_control[_mean_with_eos].npz   the control subsets themselves
  scratch/
    shards_prott5_<label>/              resume cache
    prott5_input_<stem>.tsv             what was handed to the worker
    strains/prott5_<label>.npz          screen-strain matrices (see below)
```

```
accessions   object   (n,)        row-aligned to `embeddings`
embeddings   float32  (n, 1024)   mean over residues, EOS excluded
model        scalar   "Rostlab/prot_t5_xl_half_uniref50-enc"
pooling      scalar   "mean_no_eos"
dim          scalar   1024
chunked      bool     (n,)        True where the protein was windowed
max_len      scalar   4096
```

| helper (`src/embeddings.py`) | returns |
|---|---|
| `load_prott5(species)` | `(accessions, matrix)`, row-aligned |
| `load_prott5_frame(species)` | DataFrame indexed by `uniprot_ac`, columns `t0..t1023` |
| `prott5_metadata(species)` | model, pooling, dim, how many were windowed |
| `prott5_control()` | the per-protein cosines against UniProt |

## `--strain`: training features, not deliverable matrices

`--strain LABEL` embeds a **screen strain** instead of a registry species, writing
`data/processed/embeddings/scratch/strains/prott5_<label>.npz`. Built so far:

```
prott5_kpneumoniae__ecl8__GCA_000315385.1.npz     prott5_DEG1048.npz
prott5_kpneumoniae__kpnih1__GCA_000281535.2.npz   prott5_DEG1056.npz
prott5_kpneumoniae__kppr1__GCF_000742755.1.npz    prott5_rh201207_bruchmann.npz
```

These are **not** project matrices and are deliberately not treated as such: a screen strain has no
UniProt proteome and no stage-00 table, so there is **no canonical row order**, and the key is
whatever identifier the strain's own FASTA or TSV uses (RefSeq/GenBank protein accession, locus
tag, or a DEG FASTA header). That is the point — label and feature vector then share one namespace
and the essentiality screens need no join at all. The proteome is resolved by label across three
locations (registry tier-D `.faa`, `kp_strains` `.tsv`, DEG `.faa`) rather than making every caller
know where it lives. The UniProt control is skipped for a strain run, since it exists only for
E. coli.

## Running it

```bash
~/miniconda3/envs/gradi/bin/python scripts/embeddings/prott5.py --species ecoli saureus kpneumoniae
~/miniconda3/envs/gradi/bin/python scripts/embeddings/prott5.py --species ecoli --compare-pooling
~/miniconda3/envs/gradi/bin/python scripts/embeddings/prott5.py --strain DEG1048
```

Flags: `--species` · `--strain LABEL ...` · `--pooling {mean_no_eos,mean_with_eos}` ·
`--compare-pooling` · `--control-n 300` · `--device {mps,cpu}` · `--refresh` · `-q`.

**Run it from `gradi`**, which dispatches to `gradi-loc`. Never activate `gradi-loc` to run the
stage, and never install ProtT5's `transformers` pin into `gradi`.

# Part 3 — ProteomeLM-L, contextualised

`scripts/embeddings/proteomelm.py` · gate `scripts/embeddings/orthodb_group_check.py` ·
helpers `src/proteomelm.py` · outputs `data/processed/embeddings/proteomelm_<species>[_<mode>].npz`

The third project-wide matrix, and **the only one where a protein's vector depends on the rest of
its proteome.** ProteomeLM (Malbranke, Zalaffi & Bitbol, PNAS 2026, `10.1073/pnas.2524201123`;
papers in `docs/papers/`) is a set-transformer over a whole proteome: one token per protein, each
token being that protein's mean-pooled ESM-C vector. The recipe:

1. sequences from stage 00 (`src.proteomes.load`), sorted by accession;
2. **ESM-C 600M, mean over NON-PAD tokens** — BOS/EOS *inside* the mean, ProteomeLM's own
   convention, not Part 1's;
3. one ProteomeLM-**L** forward over the entire proteome, `output_attentions=False` so the N×N
   attention is never materialised (N = 5,728 for Kp);
4. `hidden_states[8]`, **z-scored with the genome-wide mean/SD**, as the paper specifies for
   ProteomeLM-Ess's input.

**Layer 8 of 18 is the paper's own best configuration, not a guess**: *"the best performing version
of ProteomeLM-Ess is the one trained on the embeddings of layer 8 of ProteomeLM-L, yielding an AUC
of 0.93."* Intermediate layers beat the last one, consistently with their unsupervised PPI result.
A vector from the wrong layer is a different representation entirely and nothing about its shape
says so.

## It computes its own ESM-C, and must

This script does **not** import `src/embeddings.py` and does **not** read Part 1's `.npz`. The
pooling conventions differ — `ESMC_POOLING = "mean_with_bos_eos"` here, BOS/EOS stripped there — at
**cosine 0.999970 median / 0.998291 worst**. Small, and still a different convention; *External
models* says the tool's own convention wins.

The consequence is that Part 1 becomes a **control rather than a dependency**: `crosscheck_stage01`
reports the cosine against it every run and never acts on it, and agreement validates both
pipelines precisely because neither feeds the other. A missing Part 1 file is not a failure.

## NOT shardable

The ProteomeLM forward runs over the whole proteome at once, so **a shard boundary changes the
values** — the proteome *is* the context. Only the ESM-C step shards and resumes
(`scratch/shards_proteomelm/`, accession-keyed with a sequence sha256 so a changed sequence is
recomputed rather than silently served). The forward itself is 2–3 s per proteome, so there is
nothing to save.

Same reason `--limit` is **not a free smoke test**: truncating the proteome changes the values of
the proteins that remain. Smoke output goes to `scratch/smoke_*` and is not comparable to a full
run, and the script says so before it starts.

## `--group-embeds {self,orthodb}` — the load-bearing knob

ProteomeLM takes a per-protein **functional encoding** (`group_embeds`) which during TRAINING is
the **mean ESM-C embedding of the protein's OrthoDB orthologous group**; per the paper's SI §6 that
is the mechanism by which it beats its own `ProteomeLM-Discrete` ablation. It is an **additive
second branch** — `embedding_main(x) + embedding_encoder(group)` — so it changes every hidden
state, including layer 8.

| mode | what it passes | provenance |
|---|---|---|
| **`self`** (default) | `group_embeds=None`, so each protein is its own functional encoding | the authors' **released inference** default — `prepare_ppi(..., use_odb=False)`, annotated `# TODO: use odb on the fly` — **not** what the model was trained with. Every run before 2026-09-21 used it |
| **`orthodb`** | the real thing, from the authors' `group_vectors_*.pkl` joined to our own `orthodb_<species>.tsv` | what the model was TRAINED with. The authors' own route is UniProt's `xref_orthodb`, which is **0.0% on Kp**, so our DIAMOND-derived table substitutes for it |

`self` passes `None` rather than a self-copy deliberately, so the baseline stays **byte-identical**
to every run made before the flag existed.

**Mapped fraction, measured at `--min-group-size 50`: ecoli 88.2% · kpneumoniae 81.7% ·
saureus 76.9%.** Unmapped proteins **fall back to their own ESM-C vector**, which is the authors'
training dataloader's behaviour. **Always report that fraction** — where it is low the two modes
converge on the same input by construction, so a null result would say nothing about the encoding.
The run warns below 30%, and `group_mapped_frac` is stored in both the npz and the manifest.

## Run the release gate before trusting `orthodb` mode

OrthoDB group ids are *"not stable and re-used between releases"*. Ours are **`odb12v2`**; the
authors' pickles carry whatever release they trained on. **Had the id spaces differed, every lookup
would miss, every protein would fall back to self, and the run would silently reproduce `self`
while looking like a completed experiment** — the worst available outcome, so it is checked for its
own sake before any compute is spent.

`scripts/embeddings/orthodb_group_check.py` measures it. **PASSED**, at `min_group_size=50`
against 1,737,393 group vectors (`evidence/orthodb_group_vector_overlap.tsv`):

| species | `orthodb_og_domain` ids present | proteins covered | `orthodb_og_narrow` proteins covered |
|---|---|---|---|
| ecoli | **88.8%** | **85.5%** | 18.8% |
| kpneumoniae | **79.6%** | 75.0% | 54.1% |
| saureus | 82.5% | 75.8% | 56.2% |

That table also decided something unguessable: **use `orthodb_og_domain`, not `orthodb_og_narrow`**
— 85.5% vs 18.8% per-protein on E. coli, because narrow groups are clade-specific
(Enterobacteriaceae, *Klebsiella*) and mostly absent from a size-thresholded table. Both are still
tried per protein, **domain first**, mirroring the authors' `;`-separated multi-OG semantics where
the first group present in the table wins — and the union beats either column alone, which is how
Kp reaches 81.7% against 75.0% from domain alone.

**The four `group_vectors_*.pkl` files are DISJOINT SIZE BANDS, not nested supersets.** The suffix
is a group-size threshold and the authors' loader **merges** every file at or above
`min_group_size`, so a *lower* number loads *more* groups. Escalating `_200` → `_50` moved Kp
per-protein coverage 66.5% → 75.0% and E. coli 79.6% → 85.5%. `_0` is **18.1 GB and unpickles
whole (~18 GB RAM)** while our three proteomes use at most ~13,000 distinct groups — write a
filtering loader before reaching for it. Provenance, byte counts and md5s:
`data/source/proteomelm/SOURCE.md`.

## Four traps, all hit at least once

1. **`scripts/embeddings/proteomelm.py` SHADOWS the installed `proteomelm` package.** Python puts a
   script's own directory on `sys.path`, so a bare `import proteomelm` from here resolves to *this
   file* — verified — and any submodule import dies with `'proteomelm' is not a package`, which
   reads like a broken install rather than a name collision. `installed_proteomelm()` drops the
   directory from `sys.path` explicitly instead of relying on it to lose the race.
2. **The output filename must carry the mode.** It originally did not, so an `orthodb` run would
   have silently overwritten the `self` matrices — which are the comparison's baseline and are read
   by the stage-04 head comparison. `self` keeps the bare name (nothing downstream moves); every
   other mode gets a suffix.
3. **The two controls must carry the group tensor.** `check_permutation` permutes it alongside the
   proteins and `check_context` slices it alongside them. Otherwise permuting re-pairs every
   protein with a *different* protein's group vector, which is a real change to the input, and the
   control fails for the wrong reason.
4. **Group vectors are `bfloat16`** on disk while the model runs `.float()` — cast, or the forward
   dies.

## Two controls, and the one that justifies the stage

Both exit non-zero on failure.

**Permutation invariance.** There are no positional embeddings, so permuting the input and
un-permuting the output must change nothing: measured **max|diff| ~1e-05** against a
`PERM_TOLERANCE` of 1e-3. (`max_position_embeddings: 512` in the HF config is an inert inherited
DistilBert field, *not* a proteome-size cap.)

**Context sensitivity — the reason the stage exists.** The same proteins embedded inside the full
proteome versus inside half of it, median cosine, which must stay **below 0.999**; at ~1.0
ProteomeLM has collapsed to a per-protein encoder and there is nothing here that Part 1 does not
already give. Recorded per run in the manifests:

| species | `self` | `orthodb` |
|---|---|---|
| *E. coli* | 0.9783 | 0.9702 |
| *K. pneumoniae* | 0.9947 | 0.9857 |
| *S. aureus* | 0.9927 | **0.9924** |

Two things to read off it. **The margin is thin on the Gram-positive**: *S. aureus* sits ~0.007
from the floor in both modes and is the one to watch — a smaller proteome is a weaker context.
And the mode matters less than it might: the group vector is a per-protein input that dilutes the
contextual signal, so it is not a priori obvious which way this moves, and measured half-proteome
to half-proteome the two modes land within 0.009 of each other on every species.

A third number is in circulation and is **a different measurement, not a contradiction**: the
script's docstring records **0.9676** for E. coli, from 2,000 proteins embedded inside the full
proteome versus alone. The shipped control uses `len(esmc) // 2` (2,201 for E. coli), hence 0.9783.
Quote the manifest value with the species and the mode, or the docstring value with "2,000
proteins" — never one as if it were the other.

## Output

```
data/processed/embeddings/
  proteomelm_kpneumoniae.npz          proteomelm_kpneumoniae_orthodb.npz
  proteomelm_ecoli.npz                proteomelm_ecoli_orthodb.npz
  proteomelm_saureus.npz              proteomelm_saureus_orthodb.npz
  evidence/
    proteomelm_manifest.tsv           `self`    -- n, layer, perm_maxdiff, context_cosine, sha256
    proteomelm_manifest_orthodb.tsv   `orthodb` -- the same, plus min_group_size, group_mapped_frac
    orthodb_group_vector_overlap.tsv  the release gate
  scratch/
    shards_proteomelm/                the ESM-C resume cache (ProteomeLM itself does not shard)
    smoke_*                           only with --limit
    strains/proteomelm_<label>.npz    screen strains, `self` only
```

```
accessions         object   (n,)        canonical row order (src/matrices.py)
embeddings         float32  (n, 1152)   layer 8, z-scored genome-wide
model              scalar   "Bitbol-Lab/ProteomeLM-L"
layer / n_layers   scalar   8 / 18
group_embeds_mode  scalar   "self" | "orthodb"
min_group_size     scalar   50, or -1 under `self`
group_mapped_frac  scalar   the fraction that got a real group vector
esmc_pooling       scalar   "mean_with_bos_eos"
proteome_id        scalar   e.g. "UP000007841"
```

**Load through `src/proteomelm.py`** — `load`, `load_frame`, `load_lookup`, `vectors_for`,
`metadata`, `manifest`, each taking `mode="self"|"orthodb"`. The module mirrors `src/embeddings.py`
function-for-function, so code that reads one reads the other by swapping the import.

**The loaders ASSERT the file's recorded mode matches the one requested.** The two modes are the
same shape over the same accessions, so nothing about a matrix's appearance says which it is, and
fitting on one while scoring on the other would return plausible, well-formed, wrong numbers. The
filename carries the mode and `_read()` confirms the file agrees with its own name; it raises
rather than returning a representation the caller did not ask for. `vectors_for` drops missing
accessions rather than zero-filling, for the Part 1 reason.

`--strain` works here too, writing `scratch/strains/proteomelm_<label>.npz`, and **refuses
`--group-embeds orthodb`**: a screen strain has no OrthoDB table, so every protein would fall back
to its own ESM-C vector and the run would be `self` under a name saying otherwise.

## Two traps in the sequences, inherited and asserted

**Do not fetch bacteria with `reviewed:true`.** The paper's own `download_proteome` defaults to it;
for Kp HS11286 that returns a handful of entries instead of 5,728. Human is the opposite — its
unfiltered proteome is 147,506 TrEMBL-bloated entries — so the asymmetry is asserted, not assumed.

**An earlier version downloaded the proteomes from UniProt itself** and it was dropped as measured
pointless: the download was byte-identical to stage 00 (4,403/4,403 E. coli sequences, no extras
either way) while adding a real failure mode — UniProt's stream endpoint is chunked, has no
`Content-Length` to verify against, and dropped mid-transfer on the 5,728-protein Kp fetch
(`http.client.IncompleteRead`). Reading the local table cannot fail that way.

## Running it

```bash
~/miniconda3/envs/gradi/bin/python scripts/embeddings/orthodb_group_check.py --min-group-size 50
~/miniconda3/envs/gradi/bin/python scripts/embeddings/proteomelm.py
~/miniconda3/envs/gradi/bin/python scripts/embeddings/proteomelm.py --group-embeds orthodb
```

Flags: `--label {kpneumoniae,ecoli,saureus,human}` · `--size {XS,S,M,L}` (default **L**) ·
`--layer 8` · `--strain LABEL ...` · `--group-embeds {self,orthodb}` ·
`--min-group-size {0,10,50,200}` · `--shard-size 250` · `--device {auto,cuda,mps,cpu}` ·
`--limit N` · `--refresh` · `--dry-run` · `-q`.

**Env `gradi`, not `gradi-loc`** — its `fair-esm` claims the same top-level `esm` name as the
EvolutionaryScale package the ESM-C step needs. `PYTORCH_ENABLE_MPS_FALLBACK=1` is set before torch
is imported, as in Part 1.

# Part 4 — 2D projections

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

`--method umap|pacmap` re-runs the comparison into `evidence/` and never touches the deliverable.
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
    evidence/
        projection_manifest.tsv        params, KL, trustworthiness, timings
        projection_<method>_<sp>.tsv   only with --method umap|pacmap
        projection_<method>_manifest.tsv
    scratch/
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
| `load_method(species, method)` | a `--method umap\|pacmap` comparison map from `evidence/` |
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
`--dof 0.8` · `--seed 0` · `--limit N` (smoke test → `scratch/smoke_*`) · `--refresh` ·
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

Both comparison maps are kept under `evidence/projection_{umap,pacmap}_saureus.tsv` as the evidence.

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
