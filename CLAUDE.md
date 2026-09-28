# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Maintaining this file

**A pipeline step is not done until this file describes it.** In the same commit as the code, add or
update: the stage entry (what it does, what it writes, which conda env it needs), any new trap worth
warning about, and any new environment requirement. This file is the map; a stale map is worse than
none.

Do not let it grow into an inventory of everything ever tried — that is what `legacy/HISTORY.md` is
for. Keep it to what someone needs to work here *today*.

**Every path written here is a claim, and claims get checked.** When you add or edit an entry, the
paths in it must be paths that exist — a wrong path in the map is worse than no map, because it is
believed. Two mistakes already caught this way: a blanket rename put `smoke_*` files under
`evidence/` when smoke output is `scratch/` by definition, and a path rewrite produced
`scripts/localization/localization/workers/`. Both read fine; both were wrong. So:

```bash
P=~/miniconda3/envs/gradi/bin/python

$P tools/check_claude_md.py   # EVERY path this file names must exist. Non-zero exit if not.

# every script still imports and parses its CLI -- run from ELSEWHERE, which is what proves
# `parents[2]` and sys.path are right. Run plot scripts ONE AT A TIME (see the stylia trap).
for f in scripts/*/*.py scripts/*/*/*.py; do (cd /tmp && $P "$OLDPWD/$f" --help >/dev/null) || echo "BROKEN $f"; done

$P -c "from src import proteomes,embeddings,function,localization,degradability,orthology,\
ligandability,essentiality,projections,proteomelm,matrices,tabpfn; print('ok')"

$P -m src.matrices            # every matrix complete AND in canonical order
```

`check_claude_md.py` knows about three kinds of legitimate exception, so read its header before
adding one: dumps documented as deleted after use (ChEMBL, eggNOG, OrthoDB), `<placeholder>` and
`*` patterns (globbed, must match at least one real file), and a short list of paths quoted
**because they are wrong** — the cautionary examples above. If you add such an example, add it to
`NEGATIVE_EXAMPLES` rather than weakening the check.

**Four house rules for edits to this file**, all of which have been violated at least once:

1. **Name scripts by their new path** — `function/cog.py`, never `02_function_cog.py`. Same for
   docs: `docs/function.md`, not `docs/02_function.md`. The numbers are gone; see *Directory
   contract*.
2. **Name the tier, not `accessory/`** — that directory no longer exists. A file is in the task
   root, in `evidence/`, or in `scratch/`.
3. **Say which conda env a stage needs**, and whether it crosses a process boundary. The env splits
   are mandatory, not cosmetic, and each one exists because installing the dependency into `gradi`
   broke something specific.
4. **Quote measured numbers with what produced them.** "0.8738" alone is not maintainable; "ADEP4
   0.8738 under TabPFN-3.5, 5 seeds, cluster-grouped CV" survives an estimator change because the
   next person can tell whether it still applies.

**When a standing instruction from the project owner arrives, it becomes a `##` section**, not a
line inside a stage entry — *Supervised ML*, *COMPLETE matrices*, *External models* and *Directory
contract* all began that way. Stage entries describe one stage; standing rules bind every stage,
including ones not written yet, and burying one in a stage entry is how it gets missed.

## House style

**Scripts print generously.** Every script opens with a banner naming itself, its inputs and its
outputs; reports per-step progress with running counts; and closes with a summary/coverage table
carrying percentages. Nothing important happens silently. Offer `-q/--quiet` for non-interactive use.

Two corollaries, both learned the hard way in v1:

- **Verify counts against an authoritative source before trusting a payload.** An HTTP 200 is not
  evidence of data — v1 was repeatedly bitten by short or empty responses (a 61-byte 200, an empty
  202, a 200 from the wrong article). Assert, and fail loudly.
- **Never resolve ambiguity silently.** Emit the chosen value, its source, and the alternatives, plus
  an audit row. If a rule decided between candidates, record which rule.

**Figures use stylia** (`ersilia-os/stylia`), and it has a trap: `stylia/__init__.py` calls
`shutil.rmtree(matplotlib.get_cachedir())` **at import time**. Two plot scripts running
concurrently therefore delete each other's font cache and die with a `FileNotFoundError` from deep
inside `shutil`, naming `~/.matplotlib` rather than anything of ours. **Run plot scripts one at a
time.** The signature that identifies it: a *different* subset fails on each run, and every script
passes when run alone.

## Supervised ML: always TabPFN-3.5

**Whenever this project needs supervised machine learning — classification OR regression — use
TabPFN-3.5.** Not RandomForest, not a hand-tuned sklearn model, not an AutoML wrapper. This is a
standing instruction from the project owner and applies to every stage, existing and future.

**ONE STANDING EXCEPTION, granted 2026-09-21: the ESSENTIALITY ENDPOINTS run a RandomForest until
the project owner is convinced by the datasets.** Their words: *"for now, use the random forest
option (tabpfn only when i am convinced)."* So `scripts/essentiality/predict.py` defaults to
`--estimator forest` and `--estimator tabpfn` is opt-in. This is **scoped to that axis** — it is
not a licence to reach for sklearn anywhere else, and it lapses when the owner says so.

Two reasons it is also the right call on the merits. (1) The question being asked of these ten
datasets is *which of them is learnable at all*, which a free in-process estimator answers as well
as a paid one. (2) A 10-endpoint × 5-seed × 5-fold sweep is **250 hosted calls ≈ 2.5M credits,
12.5% of the monthly quota**, which is a lot to spend before knowing whether a dataset is worth
modelling. The forest is free and needs no token.

**Keep the comparison possible.** `predict.py` routes both estimators through ONE `predict_fold` /
`oof_predict` seam on identical folds, grouping and metrics, and the forest is
`src.degradability.RF_PARAMS` imported rather than re-specified — so a later forest-vs-TabPFN
result is a comparison of estimators, not of two scripts. Write `estimator` into every results row
(`evidence/cv_endpoints.tsv` does) or the numbers become unattributable.

- **Call it through `src/tabpfn.py`** — `predict_fold`, `predict_fold_regression`, `oof_predict`.
  It is transversal, so the dispatch, the content-addressed cache, the retry and the credit guard
  live in `src/` where any axis can reach them; nothing should shell out to the worker directly.
- Package `tabpfn` 9.0.0 (local weights) or `tabpfn-client` 0.6.0 (hosted); model `v3.5` is the
  default checkpoint of both. Env **`gradi-tabpfn`** — **never install it into `gradi`**, it pulls
  torch 2.14 against gradi's 2.12 and would break stage 01's ESM-C. Run it across a process
  boundary: `scripts/workers/tabpfn_cv.py`, `GRADI_TABPFN_BIN` overrides.
- **`GRADI_TABPFN_CACHE_ONLY=1` refuses to spend.** Any refactor that could disturb what the cache
  key hashes must be verified under it: if the stage still reproduces its numbers, the cache
  survived. This is the ONLY correct check — "verifying" the cache with a novel input always misses,
  which looks identical to a broken cache and proves nothing. Both mistakes have been made; the
  key change cost **1,540,000 credits**. Cache lives at `data/processed/tabpfn/cache`, shared
  across axes, filenames ARE the keys (so the directory can move safely).
- Stage 04 under that flag reproduces **ADEP4 0.8738 / PR 0.6103, ONC212 0.7671 / PR 0.5803**.
- **No PCA needed for our sizes, but know the real limits.** The installed `tabpfn` 9.0.0 declares
  `MAX_NUMBER_OF_SAMPLES = 10_000` / `50_000` and `MAX_NUMBER_OF_FEATURES = 500` / `2_000`
  depending on the inference config (`tabpfn/inference_config.py`). So 1,152-d embeddings go in
  whole and every estimator is compared on identical features — but a corpus of ~170k rows does NOT
  fit one call, which is why per-endpoint modelling is the right shape for a multi-screen axis
  rather than a compromise. **An earlier version of this file claimed "20,000 features and
  1,000,000 rows"; that was wrong** — measure against the installed package before sizing a corpus.

**The evidence, so nobody re-litigates it.** Measured on stage 04 over a 3-estimator × 3-embedding ×
2-activator grid with byte-identical cluster-grouped folds
(`output/results/degradability/proteomelm_vs_esmc.tsv`): TabPFN beat the hand-set
RandomForest on **PR-AUC in 6 of 6 arms, each winning all 5 seeds** (+0.041 / +0.014 / +0.056 on
ADEP4; +0.033 / +0.028 / +0.021 on ONC212). **lazy-qsar 3.4.4 was tested in the same grid and
rejected** — its one real gain (ADEP4 PR +0.021) does not transfer to ONC212, and it is the least
seed-stable of the three. This overturned stage 04's earlier "the estimator does not matter"
conclusion, which was true only for the models it had tried.

**Three rules that come with it:**

1. **Report PR-AUC alongside AUROC, never AUROC alone.** The gain is ~4× larger on PR-AUC —
   concentrated at the top of the ranking, which is what a prioritized shortlist consumes. AUROC
   alone called 4 of those 6 arms "no difference".
2. **Compare arms with a PAIRED test** — per-seed differences ± standard error on identical folds,
   plus `[k/5 seeds]`. Two independent SDs is badly under-powered when partition difficulty is
   shared, and mislabels real effects as noise.
3. **Estimator and features interact — do not pick embeddings with a convenient model.** ProtT5 is
   the *worst* features under a forest (AUROC −0.015, 0/5 seeds) and the *best* under TabPFN
   (+0.011, 5/5).

**Two open caveats**, neither resolved: the weights are **non-commercial licensed** (fine for
methods work; must be settled against how GraDi's outputs are used before shipping), and local
weights need the licence **accepted** at ux.priorlabs.ai — a valid `TABPFN_TOKEN` alone is not
enough, so only hosted works until someone clicks accept. Hosted bills ~10,000 credits per *call*
(flat, independent of data size; `fit` is free) against a 20M monthly quota, and uploads features
and labels to a third party — hence the content-addressed cache in the worker.

**The regressor is measured too, and also wins.** `scripts/degradability/regressor.py` puts
`TabPFNRegressor` against stage 04's ESM-C ridge on continuous log2FC, same folds, 5 seeds, paired:
ADEP4 rho +0.0365 / AUROC +0.0035, ONC212 rho +0.0263 / AUROC +0.0161 — **4 of 4 comparisons, each
winning all 5 seeds**. The ridge baseline reproduces the documented 0.515/0.874 and 0.363/0.753
first, which is what validates the harness before any credits are spent.

**An open lead, deliberately not a finding:** TabPFN *regression* ranks slightly better than TabPFN
*classification* on the same binary labels (0.8812 vs 0.8738 on ADEP4; 0.7720 vs 0.7671 on ONC212),
which would mean the continuous log2FC carries signal the `<= -1` cutoff throws away. Both gaps sit
inside the combined seed noise, so this needs a paired test between the two framings before anyone
acts on it.

## Every axis ends in COMPLETE matrices

**One row per protein, always.** A protein with no annotation is an **all-zero row carrying an
`evidence` label**, never a missing row. Standing instruction from the project owner; it applies to
every axis, existing and future.

**"Matrix" means a real feature matrix** — rows = proteins, one column per vocabulary term, binary
0/1 (or graded where the source is graded), plus one `evidence` column saying where the row came
from. A long-form table with `;`-packed term lists is the SOURCE, not the deliverable.

**AND IN THE SAME ORDER.** Every matrix in this project has the same rows in the same order:
the order `data/processed/proteomes/proteome_<species>.tsv` is written in. Standing instruction,
and it is not cosmetic — once it holds, any two axes can be placed side by side with `np.hstack` or
`pd.concat(axis=1)` and no join at all. Skipping a join that *would* have worked is not the point;
the point is the failure mode. Forgetting to write one does not raise, it silently misaligns a
column. And an embedding matrix has no accession column to join on in the first place: its rows are
positional, so alignment there is *only* ever by order.

Use `src/matrices.py`: `reindex()` / `reindex_arrays()` on the write side, `assert_canonical()` on
the read side, and `python -m src.matrices` to audit every matrix the project ships. Completeness
is checked in the same breath, because the two rules are one rule — a matrix missing a protein
cannot be in canonical order either. `reindex()` deliberately REFUSES to fill a missing protein: a
NaN row invented there is a protein the axis never measured, dressed as one it measured as unknown.
Fill it where the rows are built, and say what the fill means.

Four more rules that come with it:

1. **Full vocabulary as columns, not just observed terms**, so every species file has the same shape
   and the three stack. Terms an organism structurally cannot have stay as **structural zeros —
   that is information**: 8 of 97 GO-slim terms are eukaryote/plant concepts, COG `Y` is nuclear
   structure, and *S. aureus* has no periplasm or outer membrane because it is Gram-positive.
2. **Keep MULTI-LABEL.** The single "chosen" term discards a lot — measured: 34–52% of GO-annotated
   proteins carry more than one slim term, 12.5% of COG-classified proteins more than one letter.
3. **A zero means "not annotated", not "absent"** — for 26.3% of Kp that means nothing is known.
   Say so in the loader docstring; downstream must never read it as a measured negative.
4. **Verify with a ROUND-TRIP, not a shape check.** Reconstructing the term lists from the matrix
   must reproduce the source columns exactly. Shape checks pass on wrong matrices.

Status, from `python -m src.matrices`: **36/36 canonical.** (It was 31/33 while
`prott5_{kpneumoniae,ecoli}.npz` did not exist, then 33/33, and is now 36 because the audit also
covers `proteomelm_<species>_orthodb.npz` — a SECOND representation of the same proteins, same
shape, same accessions. An unaudited matrix is exactly the silent misalignment this rule exists to
catch, so a new representation goes into the audit list, not beside it.) Functional annotation is
done (`function/matrix.py`). Localization has two complete,
canonical tables but `deeplocpro_` is still single-label categorical and needs one-hot over the
6-class union before it is a matrix in the sense above.

## External models: never reuse an embedding you have not proven identical

**Before feeding any embedding we already hold into an EXTERNAL model, prove it is numerically
identical to what that model expects — same model TYPE and same model VERSION.** If it cannot be
proven, generate the embedding with the tool's own code path, however long that takes. Standing
instruction from the project owner; it applies to every stage.

**The bar is numerical identity, not similarity.** Cosine 0.99 is a *fail*. The check that passes
looks like `scripts/embeddings/prott5.py`, which reproduces UniProt's published ProtT5 vectors at
**median cosine 1.000000, worst 0.999994** over a 300-protein length-weighted sample, and exactly
1.000000 on a 24-protein probe. Anything short of that is an unproven assumption.

**Check all four, not just the first:**

1. **Model family** — ESM-C is not ESM-2 is not ESM-1b is not ProtT5.
2. **Model version and variant** — `prot_t5_xl_uniref50` (full) vs `prot_t5_xl_half_uniref50-enc`
   (half precision, encoder-only) are *different weights*, even though both are "ProtT5-XL-U50".
3. **Everything between the model and the vector** — pooling AND every other post-processing step.
   Each of these alone makes the numbers different:
   - **Pooling method** — mean over residues vs max vs a CLS/BOS token.
   - **Which tokens are pooled** — measured for ESM-C: ProteomeLM masks only pad tokens (BOS/EOS
     *included*) while stage 01 strips them; cosine 0.999970 median, 0.998291 worst. Small, and
     still a different convention. ProtT5 has EOS and no BOS, so "mean over residues" is ambiguous
     until you say whether EOS is in — our `mean_no_eos` beat `mean_with_eos` against UniProt
     (worst case 0.999994 vs 0.998837), which is how the choice was *measured* rather than assumed.
   - **Which LAYER** — stage 07 uses ProteomeLM-L **layer 8 of 18**, not the last. A vector from
     the wrong layer is a different representation entirely, and nothing about its shape says so.
   - **Normalisation** — z-scoring, L2, or none. Stage 07 z-scores genome-wide; a model expecting
     raw activations would be silently mis-fed.
4. **Residue encoding and truncation** — whether U/Z/O/B map to X, whether sequences are windowed,
   and at what length (ours windows above 4,096 aa and flags the one protein affected).

**Why this is a hard rule.** The failure is silent and unfalsifiable: a predictor handed a subtly
different representation from the one it was fitted on still returns plausible, well-formed numbers,
and nothing in its output says they are wrong. The saving is only compute.

**Worked example — SAFPred (stage 02).** Its shipped demo embeddings are **1280-d = ESM-1b**, not
ProtT5, so the authors' own reference data cannot even validate a ProtT5 substitution; and it calls
`bio_embeddings.ProtTransT5XLU50Embedder` on the *full* `prot_t5_xl_uniref50` while our ProtT5 comes
from TMbed's bundled *half, encoder-only* build. Verdict: **SAFPred generates its own embeddings.**
Our ProtT5 stays valid for OUR OWN models — it competed fairly in the stage-04 head comparison —
but is not fed to someone else's fitted tool.

## Status

Active **target-prioritization** analysis for the **GraDi** collaboration (BacPROTAC-style targeted
protein degradation in Gram-negatives; Prof. Erick Strauss, Stellenbosch). Ersilia's deliverable is a
prioritized list of proteins of interest. Ligand identification is out of scope.

**This is v2, a deliberate restart.** v1 is complete and frozen under `legacy/` — see the *Legacy*
section below before assuming anything about prior work.

**The consortium's own proteins of interest live in `src/interest.py`.** There was no
machine-readable list anywhere in the repo — the targets are stated in prose in two legacy documents
(`legacy/docs/03_degradability.md:420`, the ~30 envelope targets: LPS biosynthesis, LptA–G, BamA/D,
Lol/Lnt/LspA, Sec/YidC/LepB, FtsH; and `legacy/HISTORY.md:611`, the v5 proposal's GyrA/GyrB
"cytosolic", conflicting with the kick-off's "targets in the periplasm", never settled). That module
is the first time they are written down as data, and it is **an expansion of prose, flagged as a
draft**: every judgement is recorded in its `EXPANSION_NOTES`, and `PANEL` is a plain dict so
swapping in a curated file changes nothing downstream. Load through `annotate`, `coverage`,
`missing`. **Match is by gene symbol, so read `coverage()` before concluding a target is absent** —
Kp names 63.4% of its proteome and Sa 44.6%, and Sa's 11/43 match is mostly the *structural* absence
of LPS/Lpt/Bam in a Gram-positive rather than a naming gap.

## Species

Four, **one reference proteome each**. That proteome is the unit of analysis; nothing else is.

| species | proteome | strain | n | locus tags |
|---|---|---|---|---|
| *K. pneumoniae* | `UP000007841` | HS11286 | 5,728 | `KPHS_*` |
| *E. coli* | `UP000000625` | K-12 MG1655 | 4,403 | b-numbers |
| *S. aureus* | `UP000008816` | NCTC 8325 | 2,889 | `SAOUHSC_*` |
| *H. sapiens* | `UP000005640` | — | 20,416 | — |

**HS11286** is the anchor: the only *K. pneumoniae* proteome UniProt flags "Reference and
representative". Rejected as anchor — **ATCC 43816 / KPPR1** (no curated UniProt reference proteome)
and **MGH 78578** (historical reference; reviewed coverage is ~0.5% species-wide and propagates by
orthology anyway).

***S. aureus* is new in v2.** It is the organism the activated-ClpP proteomics (Conlon 2013 ADEP4,
Jacques 2020 ONC212) is native to — v1 could only reach that data by cross-phylum DIAMOND RBH at
median 42% identity — and the only one of the three bacteria with ClpC/McsB, the machine every
published BacPROTAC actually targets. Whether it becomes a full target organism or stays a reference
is deferred.

**Human must be fetched `reviewed:true`** — the unfiltered proteome is 147,506 TrEMBL-bloated entries.

## Identifier convention

**UniProt accession is the canonical key** in every dataset, `data/` and `output/` alike. Gene
symbols, locus tags, b-numbers, RefSeq and GeneID are retained as provenance columns — never as the
primary key. `KPHS_*` is a locus tag, not an accession.

**But join on `locus_tag`, not on gene name.** `locus_tag` is at 100% coverage while `gene_name` is at
18.4% on Kp, and essentially every published bacterial dataset (Tn-seq, TraDIS, CRISPRi, proteomics)
keys on locus tags. v1's single most expensive recurring bug was a gene-symbol join that silently
lost ~27% of known essentials.

**Map by sequence, not by accession, when reaching an external database.** HS11286 is a dark TrEMBL
proteome whose accessions rarely match those used by ChEMBL, BindingDB or PDB-SIFTS. v1's convention,
which worked: DIAMOND blastp, **≥95% identity = "direct"**, **≥40% floor for transfer**, and restrict
the bacterial bucket to *true Bacteria* rather than "non-human".

## Identifier mapping: spend what it takes, and track the rate

**Standing instruction from the project owner.** Identifier mapping is not plumbing to get past --
it is where this project's data is won or lost, and the losses are silent. A decimated join looks
exactly like a small dataset.

**Measured, on the essentiality screens:**

| screen | v1's route | v2's route | gain |
|---|---|---|---|
| Ramage 2017 KPNIH1 | gene symbol → 212/424 | GenBank `locus_tag` → **424/424** | **2.0×** |
| BN373 / ECL8 TraDIS | never attempted | GenBank `locus_tag` → 4,930/5,048 (97.7%) | new |
| Goodall 2018 | gene symbol | compendium symbol→b-number + symbol fallback → 97.6% | +5.3 pp |

**Four rules:**

1. **Try every route, MEASURE each, then pick.** Goodall: compendium route 96.4%, our own symbol
   index 93.2%, **union 98.6%** — so the union earned its keep. Never assume one route is enough.
2. **Prefer GenBank (`GCA_`) over RefSeq (`GCF_`) for any strain whose screen predates the current
   annotation.** PGAP re-annotation silently drops submitter locus tags — **263 of RefSeq's 281
   ECL8 misses were tags absent from the annotation entirely**, not proteinless genes, and GenBank
   was a strict superset both times (RefSeq added zero). This generalises the `saureus__col` trap
   already recorded under *Registry tiers*.
3. **NCBI GFF splits what you need across two features**: `old_locus_tag` sits on the **gene**,
   `protein_id` on the **CDS**, linked by `ID`/`Parent`. Walk both or the map comes out empty.
4. **Write the rate to an audit table every run, with a floor that exits non-zero.** See
   `data/processed/essentiality/evidence/screen_join_audit.tsv`.

**And know when to stop — then write the dead end down.** Eichelberger 2024's ECL8 screen could not
be joined at all: its `ecl8_#####` tags annotate a different assembly (475 of 5,165 numeric
suffixes shared, **0 of 5,074 coordinates matched**, symbols reached 3.7%). That is recorded in
`scripts/essentiality/screens.py` under `UNJOINABLE` so nobody re-derives it. It was superseded by
`bn373_ess`, the same organism and same assay keyed on deposited tags.

## Two-track persistence: Git vs eosvc

Enforced by `.gitignore`:

- **Git**: `src/`, `scripts/`, `docs/`, `assets/`, `legacy/` (code and docs), `LICENSE`, `README.md`,
  `install.sh`, `requirements.txt`, `access.json`.
- **eosvc (DVC + S3), not Git**: `data/` and `output/`. Use
  `~/miniconda3/bin/eosvc {view,upload,download,delete} --path <relative path>`.
- `tmp/` is local-only scratch, tracked by neither.
- The only generated artifacts in Git are `legacy/app/data/*.json`, the webapp payloads (an explicit
  `.gitignore` exception).

`eosvc upload` PUTs **one object at a time (~2 files/sec)**, so a directory of many small files is
hours, not minutes. Archive such trees before pushing — this is why `data/processed/legacy/v1/` holds
tar.gz files. And `eosvc view --path <file>` is unreliable for a single file; read the parent
directory listing's status column instead.

## Directory contract

**One folder per task, no stage numbers.** This replaced 32 flat `scripts/NN_*.py` files in
September 2026. The numbers encoded an order that was partly fiction — function, degradability and
ligands are largely independent — and they had drifted: the ProteomeLM *embeddings* script was
numbered `07_` under essentiality purely because that is when it was built, and wrote to a
`07_proteomelm/` directory matching no stage at all. **Execution order is deliberately left
unexpressed** until it is actually known; nothing in the layout implies one.

### The whole map

```
scripts/
  proteomes/      download.py                                     -> proteome_<sp>.tsv
  embeddings/     esmc.py  prott5.py  proteomelm.py  projection.py -> <model>_<sp>.npz
                  workers/prott5.py
  function/       cog.py  eggnog.py  goslim.py  matrix.py  deepgo.py
                                                                  -> goslim_matrix_/cog_matrix_<sp>.tsv
  localization/   predict.py  workers/deeplocpro.py               -> deeplocpro_/tmbed_<sp>.tsv
  orthology/      orthofinder.py  orthodb.py                      -> orthologs/neighbors/orthodb_<sp>
  degradability/  predict.py  enrichment.py  regressor.py
                  head_comparison.py  workers/lazyqsar_cv.py      -> degradability_<sp>.tsv
  essentiality/   labels.py  deg_proteomes.py  geptop.py  merge.py -> essentiality_<sp>.tsv
  ligands/        chembl.py  bindingdb.py                         -> chembl_<sp>.tsv
  studiedness/    fetch.py  gene2pubmed.py  unknome.py
                  transfer.py  merge.py                            -> studiedness_<sp>.tsv
  pockets/  interactome/                    README.md only -- real axes, no code yet
  plots/          10 scripts, ALL figures
  workers/        tabpfn_cv.py             transversal; every axis may call it

src/              flat. one module per task + matrices.py, tabpfn.py, interest.py, proteomelm.py
docs/             one .md per task, named for the task (docs/function.md, not docs/02_*.md)
tools/            one-shot migration scripts, kept for the record

data/source/<provider>/       uniprot deg eggnog orthodb cdd go sprofgo geptop ncbi unknome
data/processed/<task>/        THE DELIVERABLES -- nothing else at this level
data/processed/<task>/evidence/
data/processed/<task>/scratch/
data/processed/tabpfn/cache/  shared content-addressed predict cache
data/raw/                     the frozen v1 archive -- do NOT extend
output/{results,plots}/<task>/
```

### Where does a new file go?

| the thing | goes | the test |
|---|---|---|
| a stage script | `scripts/<task>/` | which axis does it advance? |
| a figure | `scripts/plots/` | **always** — never inside a task |
| an entry point for another conda env | `scripts/<task>/workers/`, or `scripts/workers/` if >1 task uses it | is it *imported*? then it is not a worker |
| a loader | `src/<task>.py` | stages never read another stage's files by path |
| something fetched from outside | `data/source/<provider>/` | **who published it**, not who wants it |
| the axis's answer | `data/processed/<task>/` | would you hand this to a collaborator? |
| a control, audit, CV table, manifest, vocabulary | `.../evidence/` | would you **cite or check** it? |
| a cache, shard dir, raw tool dump, smoke output | `.../scratch/` | would you **delete** it to reclaim space? |

**`data/source/` is by PROVIDER, not by task.** A source is not owned by one axis — the UniProt
proteomes feed every one of them, and filing them under `proteomes/` said otherwise. Where a file
came from is stable; which task wants it is not.

**Only deliverables sit at a task root.** `function/` holds exactly the two matrices it is
supposed to; the `cog_`/`goslim_`/`eggnog_` per-protein tables they are built *from* are evidence.

**`scratch/` is safe to purge and pointless to upload to eosvc** — it is ~18 GB of caches against
~90 small files in `evidence/`. That separation is the practical payoff of the split.

Two names that look like they belong to one tier and do not: **geptop's `cv_*.pkl` is a COMPOSITION
VECTOR cache** (scratch), not cross-validation; and a `{pre}`-prefixed file is evidence on a real
run and scratch under `--limit`, because `pre` is the smoke prefix.

Naming: `<artifact>_<species>.<ext>` — the artifact type leads, so a directory sorts by kind rather
than by organism.

### Four rules the layout depends on

1. **`REPO_ROOT = Path(__file__).resolve().parents[2]`** in every v2 script — they are one level
   deeper than they used to be. **`src/` keeps `parents[1]`** (it did not move) and **`legacy/`
   keeps `parents[1]`** (its `data`/`output` symlinks depend on it). These differences are
   deliberate; do not "fix" them.
2. **Plots live in their own top-level folder**, because figures are often comparative *across*
   categories and would otherwise have no home.
3. **Workers are filed by *how* they run, not by who calls them.** A worker is an entry point under
   a *different* conda interpreter, never imported. They do not go in `src/`: that is the package
   `gradi` imports from, and importing a worker would drag its foreign dependencies into the
   `gradi` process — the whole thing the env split prevents.
4. **A tool's own embeddings belong to the tool.** SPROF-GO's ProtT5 representations live under
   `function/scratch/`, TMbed's under `localization/`. `embeddings/` is only for the three
   first-class, project-wide matrices (ESM-C, ProtT5, ProteomeLM) — see *External models* for why a
   tool's embeddings must never be swapped for ours.

### Two documented exceptions — leave them alone

- **`data/raw/` still exists** and is the frozen v1 archive (`ecoli/`, `kpneumoniae/`, `human/`,
  `other/`, `legacy/`), referenced directly by v1 scripts through the `legacy/data -> data/`
  symlink and indexed by `data/raw/PROVENANCE.md`. v1 used an organism-first layout; do not extend
  it. It could not be renamed to `data/source/` without breaking the frozen pipeline.
- **`ligands/chembl.py` and `ligands/bindingdb.py` read their bulk dumps from
  `data/raw/other/{chembl,bindingdb}/`**, not `data/source/`. Those 11.9 GB are shared with v1
  scripts using the same path, so moving them for tidiness would break legacy. Both are deleted
  after each run and re-fetched from their `SOURCE.md` anyway.

## Pipeline (`scripts/`)

One folder per task (see *Directory contract*). Run with the `gradi` env.

- **`proteomes/download.py`** — fetch the four reference proteomes and build **one identified table
  per species**. Registry-driven from `src/proteome_registry.tsv` (41 rows, 4 tiers, every row carrying
  a required `why`).

  Its real job is fixing the naming gap: gene names cover only **18.4% of Kp** and 28.2% of
  *S. aureus* (vs 100% for E. coli), and names are what literature and databases key on. Filled in
  three labelled tiers, all **within the species** — `anchor` → `species_exact` (identical sequence)
  → `species_uniref90` (same UniRef90 cluster, via UniProt ID-mapping, no DIAMOND needed). Measured:
  **Kp 18.4% → 63.4%**, Sa 28.2% → 44.6%. Cross-species naming is deliberately deferred to a stage
  that has real orthology and can label it as the inference it is.

  `gene_name` is always **exactly one preferred name**, chosen by a total order (reviewed donor →
  real name over a `y###` placeholder → most attested → shortest, then alphabetical), so a tie can
  never leave it empty. Alternatives are kept in preference order and every fill is audited.

  Outputs **four tables and nothing else at the top level** —
  `data/processed/proteomes/proteome_<species>.tsv`, **9 columns**, keyed on `uniprot_ac`: identity and
  nothing else. Supporting detail lives in `evidence/`: `locus_tags_<species>.tsv`,
  `annotation_<species>.tsv` (the wide xref layer), `name_audit.tsv`, `registry.tsv`, `manifest.tsv`
  and the UniRef90 cache.

  **Load through `src/proteomes.py`** — `load`, `load_all`, `load_locus_tags`, `with_locus_tags`,
  `load_annotation`, `id_bridge` — rather than adding stacked artifacts. The first run wrote a
  `proteins.parquet` and an `id_bridge.tsv` that were verified pure derivations of the four tables:
  27 MB for zero new information. **Locus tags are the join key for published bacterial data**
  (Tn-seq, TraDIS, CRISPRi, proteomics), so reach for `with_locus_tags()` before joining anything
  external.

  CLI: `--tier A,B,C,D` · `--only LABEL` · `--refresh` · `--dry-run` · `-q`.
  Figures: `scripts/plots/proteomes.py` (stylia) → `gene_names.png`: coverage by fill tier, how
  the 363 contested fills were decided, and **identifier coverage — the figure behind the
  join-on-`locus_tag` rule**. Two corrections it surfaced: `locus_tag` is **98.4% on *S. aureus***,
  not 100% (so a locus-tag join is near-lossless, not lossless), and human has **no locus tags at
  all** — the rule is a bacterial one; use `uniprot_ac` for human.
  Details, traps and the run log: `docs/proteomes.md`.

- **`embeddings/esmc.py`** — ESM-C 600M mean-pooled embeddings, one 1152-dim vector per protein, for the
  **three bacteria** (13,020 proteins, 3.8M residues, ~45 min on MPS at ~1,700 aa/s). Writes
  `data/processed/embeddings/embeddings_<species>.npz` (`accessions` + `embeddings` float32, BOS/EOS stripped)
  plus a resumable shard cache in `scratch/shards_esmc/`. **Load through `src/embeddings.py`** — the
  `accessions` object array needs `allow_pickle=True`, and `vectors_for()` drops missing accessions
  rather than zero-filling, since a zero row is a real position in embedding space.
  Details: `docs/embeddings.md`.

  **`esmc_600m` is the largest ESM-C with local weights** — `esmc-6b` is Forge-API-only (token, third
  party, not reproducible offline). Weights are already in the HF cache from v1.

  **Two deliberate exclusions.** *Human* is out: it is 75% of all residues and holds every
  pathological length (titin 34,350 aa). `--species human` adds it back but would need a
  chunking/truncation policy this stage does not have. *ProtT5* is **no longer deferred** — it has
  its own sibling script, `embeddings/prott5.py`, described below. The reason it needed one: UniProt
  publishes precomputed ProtT5 h5 for **only 8 proteomes** — E. coli and human are there, **Kp
  HS11286 and Sa NCTC 8325 both 404**, and the all-Swiss-Prot file would cover 7 of 5,728 Kp
  proteins. Covering the anchor meant running ProtT5-XL-U50 locally, validated against UniProt's
  E. coli h5.

  Unlike v1's `01a`, this **shards and resumes** — v1 held everything in memory and wrote once at the
  end, so a killed run lost the lot. Shard size is in the cache filename so changing it cannot reuse
  mismatched shards.

- **`embeddings/prott5.py`** — **ProtT5-XL-U50 per-protein embeddings, and the control that
  makes them trustworthy.** A sibling of `esmc.py`, writing `prott5_<species>.npz` (1024-d, float32)
  beside it. The second of the three project-wide embedding matrices.

  **The control is the point, not an extra.** E. coli MG1655 is one of the 7 proteomes UniProt
  publishes precomputed ProtT5 vectors for, and it is our exact anchor — so the same model over the
  same sequences must reproduce the same vectors. Measured on a 300-protein sample **weighted to the
  longest** (18–2,339 aa, where a pooling or windowing bug actually shows): **median cosine
  1.000000, worst 0.999994, 100% above 0.99.** Nothing is written if that fails. This is the
  worked example behind *External models* above — the bar is numerical identity, not similarity.

  **It settled the one thing UniProt's README leaves unstated: how residues are pooled.** ProtT5 has
  an EOS token and no BOS, so "mean over residues" is ambiguous until you say whether EOS is in.
  Scored both ways against the reference: **`mean_no_eos` worst-case 0.999994 vs `mean_with_eos`
  0.998837**, so the choice is measured rather than assumed (`--compare-pooling` re-runs it;
  `evidence/prott5_pooling.tsv` is the record).

  **Nothing is downloaded but a 10 MB reference h5.** The weights are already on disk: TMbed ships
  `Rostlab/prot_t5_xl_half_uniref50-enc` (1024-d, fp16, encoder-only, 2.4 GB) inside its package in
  `gradi-loc`. Inference therefore runs **across a process boundary** into that env
  (`embeddings/workers/prott5.py`, `GRADI_LOC_BIN` overrides) — ProtT5 needs
  `transformers==4.44.2`, because 5.x routes `T5Tokenizer` through the tiktoken converter and dies
  with a spurious tiktoken error.

  **This is the half, encoder-only build — which is why it must not be fed to someone else's tool.**
  `prot_t5_xl_uniref50` (full) and `prot_t5_xl_half_uniref50-enc` are different weights. Our vectors
  are valid for OUR models (ProtT5 competed fairly in the stage-04 head comparison, where it was the
  *best* features under TabPFN and the *worst* under a forest) but are not substituted into an
  external predictor — see *External models*, and SAFPred as the worked rejection.

  Shards of 250 into `scratch/shards_prott5_<species>/` and resumes. Sequences above 4,096 aa are
  windowed, and the affected protein is named rather than silently truncated.
  CLI: `--species` · `--pooling {mean_no_eos,mean_with_eos}` · `--compare-pooling` · `--control-n` ·
  `--device {mps,cpu}` · `--refresh` · `-q`. Details: `docs/embeddings.md`.

- **`embeddings/proteomelm.py`** — **contextualised embeddings: the third project-wide matrix, and
  the only one where a protein's vector depends on the rest of its proteome.** ProteomeLM-L
  (Malbranke, Zalaffi & Bitbol, PNAS 2026; papers in `docs/papers/`) takes a whole proteome's ESM-C
  vectors in one forward pass and returns **layer 8 of 18**, z-scored genome-wide. Writes
  `proteomelm_<species>.npz` for the three bacteria.

  **Layer 8 of L is not a guess — it is the paper's own best configuration** for essentiality:
  *"the best performing version of ProteomeLM-Ess is the one trained on the embeddings of layer 8
  of ProteomeLM-L, yielding an AUC of 0.93."* Intermediate layers beat the last one, consistently
  with their unsupervised PPI result.

  **It computes its OWN ESM-C and must**: ProteomeLM pools over non-pad tokens, so **BOS/EOS are
  INSIDE the mean** (`ESMC_POOLING = "mean_with_bos_eos"`), while `embeddings/esmc.py` strips them.
  Measured difference cosine 0.999970 median / 0.998291 worst — negligible in size, still a
  different convention, and *External models* says the tool's own convention wins. Stage 01's npz
  is cross-checked and **reported, never depended on**.

  **NOT SHARDABLE.** A shard boundary changes the values, because the whole point is that the
  proteome is the context. The ESM-C *inputs* shard and resume; the ProteomeLM forward is one pass
  (2–3 s per proteome).

  **`--group-embeds {self,orthodb}` — the functional encoding, and this is the load-bearing knob.**
  ProteomeLM takes a per-protein `group_embeds` which during TRAINING is the **mean ESM-C embedding
  of the protein's OrthoDB orthologous group**; per the paper's SI §6 that is the mechanism by which
  it beats its own `ProteomeLM-Discrete` ablation. It is an **additive second branch**
  (`embedding_main(x) + embedding_encoder(group)`), so it changes every hidden state.
  - **`self`** (default, and what every run before 2026-09-21 used) passes `group_embeds=None`, so
    each protein is its own functional encoding. That is the authors' **released inference**
    default — `prepare_ppi(..., use_odb=False)`, annotated `# TODO: use odb on the fly` — not what
    the model was trained with.
  - **`orthodb`** uses the real thing, from the authors' `group_vectors_*.pkl` tables joined to our
    own `orthodb_<species>.tsv`. Their own route is UniProt's `xref_orthodb`, which is **0.0% on
    Kp**, so our DIAMOND-derived table substitutes for it.
  - Mapped fractions, measured: **ecoli 88.2% · kpneumoniae 81.7% · saureus 76.9%** at
    `--min-group-size 50`. Unmapped proteins fall back to their own ESM-C vector, which is the
    authors' training dataloader's behaviour. **Always report that fraction**: where it is low the
    two modes converge on the same input and a null result says nothing about the encoding.
  - Provenance, byte counts, md5s and the release check: `data/source/proteomelm/SOURCE.md`.

  **RUN `embeddings/orthodb_group_check.py` BEFORE TRUSTING `orthodb` MODE.** OrthoDB group ids are
  *"not stable and re-used between releases"*; ours are `odb12v2` and the authors' pickles carry
  whatever they trained on. Had the id spaces differed, every lookup would miss, every protein
  would fall back to self, and the run would **silently reproduce `self`** while looking like a
  completed experiment. Measured and **PASSED**: E. coli 88.8% of group ids present, Kp 79.6%.
  That check also decided something unguessable — **use `orthodb_og_domain`, not
  `orthodb_og_narrow`** (85.5% vs 18.8% per-protein on E. coli; narrow groups are clade-specific
  and mostly absent from a size-thresholded table). Both are tried, domain first, mirroring the
  authors' `;`-separated multi-OG semantics. Evidence:
  `data/processed/embeddings/evidence/orthodb_group_vector_overlap.tsv`.

  **Four traps, all hit at least once:**
  1. **`scripts/embeddings/proteomelm.py` SHADOWS the installed `proteomelm` package.** Python puts
     a script's directory on `sys.path`, so a bare `import proteomelm` resolves to *this file* and
     any submodule import dies with `'proteomelm' is not a package` — which reads like a broken
     install. `installed_proteomelm()` drops the directory explicitly.
  2. **The output filename must carry the mode.** It originally did not, so an `orthodb` run would
     have silently overwritten the `self` matrices that are the comparison's baseline and are read
     by the stage-04 head comparison. `self` keeps the bare name; other modes get a suffix.
  3. **The two controls must carry the group tensor.** `check_permutation` permutes it alongside
     and `check_context` slices it alongside; otherwise permuting the proteins re-pairs every
     protein with a *different* protein's group vector and the control fails for the wrong reason.
  4. **Group vectors are `bfloat16`** and the model runs `.float()` — cast, or the forward dies.

  **The four `group_vectors_*.pkl` files are DISJOINT SIZE BANDS, not nested supersets.** The
  suffix is a group-size threshold and the authors' loader MERGES every file at or above
  `min_group_size`, so a *lower* number loads *more* groups. `_0` is 18.1 GB and unpickles whole
  (~18 GB RAM) while our proteomes use ~13,000 distinct groups — write a filtering loader before
  reaching for it.

  **Two controls that are the reason to trust the matrix**, both exiting non-zero on failure:
  permutation invariance (there are no positional embeddings; measured max|diff| ~1e-05) and
  **context sensitivity** — the full-proteome vs half-proteome cosine, which must stay *below*
  0.999 or ProteomeLM has collapsed to a per-protein encoder and the stage has no reason to exist.
  Measured under `self`: E. coli 0.9676. **Under `orthodb` it RISES** (E. coli 0.9702,
  Kp 0.9857, **Sa 0.9924**) — expected, since the group vector is a per-protein input that dilutes
  the contextual signal, but *S. aureus* is close enough to the floor to watch.

  **Load through `src/proteomelm.py`** — `load`, `load_frame`, `load_lookup`, `vectors_for`,
  `metadata`, `manifest`, each taking `mode="self"|"orthodb"`. **The loaders ASSERT the file's
  recorded mode matches the one requested**: the two modes are the same shape over the same
  accessions, so nothing about a matrix says which it is, and fitting on one while scoring on the
  other would return plausible wrong numbers. `vectors_for` drops missing accessions rather than
  zero-filling.

  CLI: `--label` · `--size {XS,S,M,L}` · `--layer` · `--group-embeds {self,orthodb}` ·
  `--min-group-size {0,10,50,200}` · `--shard-size` · `--device` · `--limit` · `--refresh` ·
  `--dry-run` · `-q`. Details: `docs/embeddings.md`.

- **`embeddings/projection.py`** — **the 2D map of that space**, a sibling of the stage rather
  than a stage of its own (the `function/cog.py` / `function/eggnog.py` precedent), writing
  `data/processed/embeddings/projection_<species>.tsv`: `uniprot_ac · tsne_x · tsne_y`, three
  columns, 100% coverage, **no colour column** — join COG/localization/degradability on `uniprot_ac`.
  2.3 min for all three species; nothing to install.

  **The recipe is inherited from v1, not re-derived: openTSNE multiscale (perplexities 50/500),
  cosine, on PCA-50 of the z-scored embeddings, `dof=0.8`.** It won a 33-config sweep over input
  representation, metric, t-SNE structure and method (beating UMAP, whose families blend, and PaCMAP,
  which splits into two masses with filament tails). **Do not re-run that sweep** — the sweep code
  was deleted, so the verdict table in `docs/embeddings.md` Part 2 is the record. `dof` is the
  load-bearing knob, not perplexity: 0.5–0.6 flings clusters apart, 1.0 smears them. `--method
  umap|pacmap` re-runs a comparison into `evidence/` and never touches the deliverable.

  **The control is `trustworthiness` (k=10, cosine) with a floor of 0.90**, because a diverged t-SNE
  still returns finite, plausible-looking numbers. The floor is calibrated on measured baselines:
  random 2D scores **0.4996**, PCA-2 alone **0.7821**, this recipe **0.9752**. Measured: Kp 0.9790 ·
  Ec 0.9735 · Sa 0.9752.

  **Coordinates are relative and per-species** — three independent embeddings, no shared frame. One
  species per panel; never compute a distance or cluster across species on these columns (use the
  1152-dim space, which *is* shared). `n_jobs=-1` means the floats are not bit-reproducible across
  machines, only the structure — which is why the guard is trustworthiness and not a checksum.

  **Load through `src/projections.py`** — `load`, `load_all`, `load_method`, `coords_for`,
  `manifest`. `coords_for` drops missing accessions rather than zero-filling: **(0, 0) is a real
  position on the map**, in the dense centre.

  Two v1 mistakes not carried over: its shipped figure never used the `family` column it persisted
  (it coloured by density instead), and it fed t-SNE the PCA-50 while clustering the full 1,152-d
  matrix. Neither is reproduced.

  CLI: `--species` · `--method` · `--pca` · `--perplexities` · `--dof` · `--seed` · `--limit` ·
  `--refresh` · `--dry-run` · `-q`. Figures: `scripts/plots/projection.py` (stylia,
  uncoloured, density-faded). Details, traps and the run log: `docs/embeddings.md` Part 2.

- **`function/cog.py`** — **broad functional categories, from sequence.** COG is the first
  of two classifications this stage produces, which is why the stage is `02_function` while the
  artifact is `cog_<species>.tsv`; a further scheme gets `<scheme>_<species>.tsv` beside it and
  `<scheme>_*` files in `evidence/`.

  **The lookup route is not available**: UniProt's `eggnog` xref is **0.00% (n=0) on Kp HS11286**
  (92.98% Ec, 79.23% Sa) — the dark-TrEMBL problem landing on the anchor organism. So it is computed:
  **COGclassifier 2.0** runs RPS-BLAST against the CDD COG profile DB and maps best hit → CDD ID →
  COG ID → one of **26 COG2024 category letters** in 4 groups. ~208 MB of downloads, **9 min** for
  all 13,020 proteins. eggNOG-mapper would have been ~20 GB unpacked for the same letter.

  Measured, first full run: **Kp 79.1% classified / 72.8% informative · Ec 84.5% / 78.0% ·
  Sa 73.1% / 65.7%**. *Informative* excludes `R` (general function prediction only) and `S`
  (function unknown) — classified, but not an answer. For scale, v1's keyword `functional_class`
  heuristic put 75.4% of Kp outside `other`, without being a controlled vocabulary at all.

  **Coverage is not 100% and must not be forced to be.** NCBI's own curators cover **81.6%** of
  E. coli K-12; this stage covers 84.5%, and the union of the two is 84.9% — there is no headroom.
  The unclassified are short (median ~90 aa vs ~295), uncharacterised, and for Kp/Sa ~60% carry no
  Pfam and no InterPro either. Raising `--evalue` buys noise: on shuffled-sequence decoys the hit
  rate goes **0.4% → 24% → 84.3%** at e-value 1e-2 → 1 → 10, while real coverage goes 73% → 84% →
  98%. The default sits at the knee. Details: `docs/function.md`.

  **`rpsblast` is borrowed from the `gradi-prokka` env** by prepending its `bin` to `PATH`
  (`GRADI_RPSBLAST_BIN` overrides) — the `GRADI_DIAMOND_BIN` pattern. Do **not** `conda install
  blast` into `gradi`: there is no osx-arm64 build, so it would drag the whole env to osx-64 and
  take ESM-C down with it.

  **Human is excluded by construction**, not by choice: COG2024 is 2,296 genomes of bacteria and
  archaea with **zero eukaryotes**, so there is no honest letter to give it. `--species human` is not
  offered.

  **The E. coli control**: MG1655 is itself a COG2024 reference genome, so NCBI publishes curated
  assignments for our exact anchor proteome. Agreement **97.79% on category** (94.91% on COG id)
  over 3,575 comparable proteins; the script exits non-zero below a 90% floor, and below any failed
  spot check (`clpP/clpX/clpA/clpC/clpB` → `O`, `rpsA/rplB` → `J`).

  **Load through `src/function.py`** — `load_cog`, `load_cog_all`, `load_cog_counts`, `informative`.
  `cog_category` is one letter; `cog_category_all` keeps the full string, because 9–16% of classified
  proteins carry a multi-letter COG and taking the first is a real choice.

  CLI: `--species` · `--evalue` · `--threads` · `--limit` · `--refresh` · `--no-control` ·
  `--dry-run` · `-q`. Figures: `scripts/plots/function.py` (stylia).
  Details, traps and the run log: `docs/function.md`.

- **`function/eggnog.py`** + **`function/goslim.py`** — **the second scheme: GO slim.** COG
  answers *which orthologous group*; this answers *what does it do*, in **`goslim_prokaryote`**
  (97 terms, three aspects, one chosen term per aspect plus the full multi-label set). Two tiers,
  both real annotation: `curated` (UniProt's own GO) → `eggnog` (**eggNOG-mapper v2**).

  **Only orthology can add anything here, and that was measured twice.** The proteins COG misses
  carry no Pfam and no InterPro either (~60% Kp/Sa), so every domain-based route re-derives what
  stage 00 already has: **InterPro2GO as a tier filled 5 proteins out of 13,020**, and querying
  InterPro's *live* API gained a GO term for **0 of 15** GO-less Kp proteins — UniProt's electronic
  GO is already InterPro2GO-derived. InterProScan and `Pfam+pfam2go` would hit the same wall;
  InterProScan is `linux-64` only in any case. `interpro2go` is still scored every run as a
  **REJECTED ALTERNATIVE**; do not re-add it.

  **Coverage is NOT 100%, deliberately.** An earlier version reached 100% via a k-NN transfer in
  ESM-C space; it was removed on instruction to use well-established tools only. That was right
  twice: a hand-rolled k-NN is not citable, and its honest accuracy was far below its headline —
  leave-one-out on E. coli gave MF 86.8%, but **reweighted to the donor distances that actually
  occurred it was MF 68% / BP 53% / CC 66%**. A protein no established tool can annotate now gets an
  **empty row**, printed under `NOT ANNOTATED`. **Do not reintroduce it.**

  **`gradi-emapper` (osx-64, Rosetta)** — no arm64 build; installing into `gradi` would take ESM-C
  down with it. Two traps: bioconda puts `diamond`/`mmseqs` in the env `bin/` while emapper looks in
  `site-packages/eggnogmapper/bin/` (symlink them), and **`download_eggnog_data.py` is broken** —
  it fetches from the dead `eggnogdb.embl.de` and prints `Finished.` with exit 0 having downloaded
  nothing. Use `eggnog5.embl.de`, resume with `curl -C -` (it drops long transfers), and verify
  against `Content-Length`; note curl's own `--retry` restarts a `-C -` transfer from zero.

  **Measured, and the honest verdict: the gain is small.** Kp **69.7 → 73.7%**, Ec 87.0 → 88.7%,
  Sa 65.0 → 65.6% — **322 proteins out of 13,020**. The bottleneck is not orthology: eggNOG places
  ~73% of the unannotated proteins in a group, but only **2–18% of those groups carry any GO**.
  Ceilings estimated beforehand from eggNOG *xref* coverage (96%/82%) were far too optimistic
  because they assumed groups carry GO. **Test the OG→GO yield before paying for a database this
  size.** Quality when it does answer is the best of anything tried: MF 80.8 / BP 82.4 / CC 81.1%
  on the E. coli control, vs 68/53/66 for the rejected ESM-C tier.

  Kept anyway because the same run yields what nothing else does: **KEGG KO 64.1% on Kp** (new —
  UniProt carries only KEGG *gene* ids), preferred names 67.0% (above stage 00's 63.4%
  `gene_name`, so a candidate for closing that gap), EC 28.6%.

  **The three-way COG cross-check vindicates Part 1's tool**: against NCBI's curated COG2024 on
  E. coli, **COGclassifier agrees 97.8%** while emapper's `COG_category` agrees only 63.3% (the two
  tools agree with each other 62.7%). emapper's letters come from eggNOG's own OG→category mapping,
  not NCBI's — a different vocabulary. Do not swap COGclassifier out for it.

  The database unpacks to **50.6 GB** (`eggnog.db` alone is 41 GB — not the ~21 GB first estimated).
  Room was made by deleting `data/raw/other/chembl/chembl_37/` (28 GB), `PROVENANCE.md` Rule 1,
  archive verified with `tar -tzvf` first. **The database was then deleted after the run** — every
  derived artifact is saved and only re-running needs it back; recovery procedure and byte counts
  are in `data/source/eggnog/SOURCE.md` with `fetch_eggnog.sh`. It is public and
  re-derivable — **never upload it to eosvc.**

  **`eggnog_<species>.tsv`** is a first-class artifact from the same run: `kegg_ko`, `kegg_pathway`,
  `kegg_module`, `brite`, `ec`, `cazy`, `pfams`, `preferred_name`, `description`, and its own
  `cog_category` — an independent second opinion on the COG letter.

  **Load through `src/function.py`** — `load(species)` gives both schemes joined; also `load_cog`,
  `load_goslim`, `load_eggnog`, `load_goslim_terms`, `confident`, `informative`.

  Needs `goatools`. Details, traps and the run log: `docs/function.md`.

- **`function/matrix.py`** — **the stage-02 deliverable: two COMPLETE binary matrices per
  species.** `goslim_matrix_<species>.tsv` (n × 99: `uniprot_ac` + 97 GO-slim terms + `evidence`)
  and `cog_matrix_<species>.tsv` (n × 28: + 26 COG letters + `evidence`). **Recomputes nothing** —
  it reshapes `cog_` and `goslim_`; seconds to run, no database.

  Built from the `*_all` columns so **multi-label is preserved**, which is the point: 34–52% of
  GO-annotated proteins carry >1 slim term (max 7) and 12.4–12.6% of COG-classified proteins carry
  >1 letter (max 4). Reconciled exactly against the independently-computed
  `evidence/cog_counts_<species>.tsv`, which counts only the CHOSEN letter — *chosen + extra ==
  matrix total* to the unit (Kp 4,528 + 591 = 5,119).

  **Always-zero columns are kept and are information**: goslim Kp 16 / Ec 8 / Sa 19, cog 1
  everywhere (`Y`, nuclear structure). *S. aureus* has the most because it is Gram-positive.
  **A zero means NOT ANNOTATED** — 1,506 Kp proteins (26.3%) are all-zero because nothing is known;
  `evidence` (`curated`|`eggnog`|`none`, `cogclassifier`|`none`) separates that from a real negative.

  **`eggnog_<species>.tsv` is demoted, not deleted** — the goslim `eggnog` tier derives from its
  `gos` column, so deleting it makes the stage non-regenerable without the 50.6 GB database. Leaving
  the deliverable set: `pfams` 85.4%, `kegg_ko`/`brite` 64.1%, `description` 86.8%, `preferred_name`
  67.0%, `ec`, `kegg_pathway`, `kegg_module`, `cazy`, and emapper's independent COG call covering
  **480 Kp proteins COGclassifier misses**.

  **The COG vocabulary is now VENDORED** at `data/source/cdd/cog_func_category.tsv`: it was
  read from COGclassifier's installed package resource, and a matrix whose schema depends on a pip
  install is not reproducible. Column order follows that file (grouped, not alphabetical); it has no
  trailing newline, so `wc -l` says 25 for 26 rows.

  Verified by **round-trip** — rebuilding the `;`-joined lists from the matrix reproduces
  `goslim_*_all` and `cog_category_all` exactly for all 13,020 proteins.
  **Load through `src/function.py`** — `load_goslim_matrix`, `load_cog_matrix`,
  `load_goslim_matrix_all`, `load_cog_matrix_all`, `matrix_manifest`.
  CLI: `--species` · `--dry-run` · `-q`. Details: `docs/function.md`.

- **`degradability/predict.py`** — **is this protein a substrate of activated partnerless ClpP?** A binary
  **TabPFN-3.5** classifier on stage-01 ESM-C embeddings, trained on the two *S. aureus* activator
  screens and applied to the rest of Sa and all of Ec/Kp. (`--estimator forest` regenerates the
  former RandomForest numbers; see the retrofit note below.) Output
  `degradability_<species>.tsv`: `<act>_hit` (the measured 1/0, empty where unmeasured) ·
  `<act>_prob` (the model's probability for EVERY protein, **out-of-fold** where labeled, so it is
  one comparable scale across all 13,020) · `<act>_source` · `nn_similarity`.

  **This overturns v1's headline negative result.** v1's gate — *does any sequence feature beat
  protein length alone?* — answered no (disorder +0.009; a 13-feature GBM 0.762 vs logistic
  length-only 0.775; motifs 0% of SHAP) and concluded *"needs new data, not new features"*. Measured
  under cluster-grouped CV repeated over 5 seeds: **ADEP4 0.8738 ± 0.0029 (PR-AUC 0.6103)** against a
  properly-estimated length baseline of 0.776; **ONC212 0.7671 ± 0.0060 (PR-AUC 0.5803)** against
  0.687. The forest it replaced scored 0.8645 / 0.7517 (PR 0.5632 / 0.5481). On the cross-activator test — the only non-circular one, which v1 flagged "available now,
  not done" — **0.745 vs v1's 0.6916** and **0.848 vs 0.7650**.

  **The TabPFN retrofit, and the four traps it introduced.** The estimator lives in `gradi-tabpfn`
  and is reached through a **process boundary** (`predict_fold` in the stage → `scripts/workers/
  tabpfn_cv.py`); `GRADI_TABPFN_BIN` overrides, `GRADI_TABPFN_HOSTED=0` asks for local weights.
  1. **The hosted API is a third-party uptime dependency.** A GCS `500 InternalError` — server-side,
     nothing wrong with our data — killed a 40-minute run at call 64 of ~160. Every call now retries
     5× with 5/10/20/40 s backoff. A failed attempt is not billed.
  2. **Credits are metered but flat per CALL** — measured 10,000 whether the test set is 350 or
     2,889 rows, and `fit` is free. A full stage-04 run is ~160 calls ≈ 1.6M of a 20M monthly quota.
     The **content-addressed cache** (`data/processed/tabpfn/cache/`, keyed on sha256 of X_train/y_train/
     X_test) is therefore load-bearing: a re-run that recomputes nothing spends nothing, and an
     interrupted run resumes for free. **But the key format is silently fragile** — any change to
     what `_key` hashes invalidates every entry with no error, just universal misses. Measured:
     adding one field to the hashed string cost **1.54M credits (7.7% of the quota)** for numbers
     that reproduced identically. To check whether a call is cached, compute the key and test for
     the file; never "verify the cache" with a made-up input, which always misses and proves
     nothing.
  3. **Local weights need the licence ACCEPTED**, not just a token — verified: token valid,
     `accepted: False`, download refused. Until someone accepts at ux.priorlabs.ai only hosted
     works, which uploads features and labels to a third party.
  4. **`model_<activator>.npz` replaced `model_<activator>.joblib`.** TabPFN learns *in context*, so
     there is no fitted object to pickle — the artifact is the **training set + config**, which is
     what actually determines the predictions and survives a library upgrade.

  **The weights are non-commercial licensed.** Fine for methods work; unresolved for GraDi's
  deliverable, and it must be settled before these outputs ship externally.

  **Report mean ± SD over `N_CV_SEEDS = 5`, never one seed.** A single-seed estimate carries ±0.004
  of pure arbitrariness — larger than any of the three effects tested and rejected below — and seed 0
  alone read 0.857 against a 5-seed mean of 0.865. `_prob` is the **seed-averaged** out-of-fold value.

  **`CROSS_ASSAY_AUROC` = 0.877 / 0.815 is the yardstick, not 1.0** — what one activator's own
  measured readout achieves at predicting the other's calls. **Deliberately not called a "ceiling"**:
  the model exceeds it at three of five cutoffs, and v1's §10.3 uses that word for a different
  quantity (label-vs-label agreement, ρ 0.52 / Jaccard 0.32).

  **Three things were tested and none helps** — the useful result, because it puts the headroom in the
  labels not the method. (a) **Localization/COG features add +0.002 / −0.005**: the embedding already
  encodes them (DeepLocPro is ESM-2-based). (b) **Hyperparameters do not matter** — a 10-point sweep
  spans 0.854–0.862 and its winner does not survive reseeding. (c) **The log2FC cutoff must NOT be
  tuned**: AUROC is not comparable across labels, so tightening it inflates the score mechanically
  (ADEP4 0.827 → 0.943 from −0.5 to −3.0). `evidence/cutoff_sensitivity.tsv` sweeps it anyway
  because `CROSS_ASSAY_AUROC` rises in lockstep — the **gap** is the invariant and it is stable at
  every cutoff, which makes the result cutoff-independent.

  **A mechanistic finding for downstream filtering.** Hit rate by compartment: cytoplasm 0.180/0.275,
  **membrane 0.031/0.105**, cell wall 0.000, but **extracellular 0.049/0.346**. Membrane proteins are
  protected (co-translational insertion, never a soluble cytoplasmic chain); secreted ones are **not**
  (they transit the cytoplasm unfolded). So **"cytoplasmic only" is the wrong filter** — "not
  membrane" is closer. The screens disagree on extracellular, so that class is unsettled.
  Figures (stylia, both to `output/plots/degradability/`):
  `scripts/plots/degradability.py` — per-fold ROC + the cross-assay reference as its own curve, i.e.
  *how good is the model*. `scripts/plots/degradability_predictions.py` — *what did it predict, and
  what is it worth*: `predictions.png` (score distributions per species/activator + Sa split by
  measured status), `agreement.png` (the two activators' predictions at rho 0.83–0.84 while their
  labels agree at 0.52 — one opinion, not two), `extrapolation.png` (`nn_similarity`, band
  composition, and the measured per-band AUROC). **Read `extrapolation.png` panel B before quoting
  any expected AUROC: a fifth of Kp sits in bands with no estimate at all**, so the reweighted figure
  covers only 79.3% Kp / 87.5% Ec / 92.1% Sa.
  `scripts/plots/degradability_top.py` — **the named top of the ranking**: `top_<activator>.png`
  (the 20 highest-ranked proteins per species, labelled `gene_name` → eggNOG `preferred_name` →
  accession, coloured by COG group, `*` = a measured Sa hit) and `top_composition.png`. Three things
  it establishes: the top 100 is **94–97% cytoplasmic** (vs ~60% of the proteome) and 36–44%
  information-storage (vs ~16%), which is the right direction for an axis about Clp reach and was not
  built in; **9/13 measured proteins in Sa's ADEP4 top 20 are confirmed hits, 14/16 for ONC212**; and
  **the two activators' top-100 lists share only 24/100 on Kp** (32 Ec, 40 Sa) despite rho 0.835 —
  *global rank correlation is not shortlist agreement*, so do not treat either list as the answer.
  ONC212's top is almost pure ribosome, ADEP4's is broader (H-NS, cold-shock, PPIases, chaperones).
  CLI: `--top` · `--composition-top` · `--species` · `-q`.
  `scripts/degradability/enrichment.py` — **Fisher-exact COG enrichment of the full-proteome
  predictions**, `cog_fisher.png`. **Hits = the top 10% of each proteome** (573 Kp / 440 Ec /
  289 Sa), *not* a fixed count: the proteomes differ 2× in size, so an absolute N would make "the
  top" 1.7% of Kp against 3.5% of Sa and the odds ratios incomparable across species — and at a
  fixed top-100 most small categories had nothing in scope and returned OR 0.00, q 1.0, which looks
  like a finding and is not one. Two-sided, BH **within each species×activator×test-family block**
  (pooling COG with the panel made a category's q depend on how many panel families sat beside it).
  194 tests, 67 significant, raw 2×2 counts in
  `output/results/degradability/enrichment_fisher.tsv`.

  **Judge enrichment against the LABELS, not against the previous model — this was got wrong once.**
  The forest's shipped text claimed five categories "move the same way in all six
  species×activator combinations" (J up 3.1–10.2, O up 1.8–2.8; M/E/G down) and called that the
  axis's best sanity check. **Two of those were forest artifacts.** Measured directly on the
  labels, with no model involved: under ONC212, **COG O (chaperones) is NOT enriched among hits**
  — 17.1% against a 24.6% base rate, OR 0.63, p 0.42 — yet the forest predicted a 1.8–2.8
  enrichment. TabPFN predicts 0.86–1.03 there, i.e. it tracks the measurement. Likewise J under
  ADEP4: label OR **1.81**, TabPFN 1.47–2.59, forest 3.1–10.2 (over-enriched).

  So the criterion is **does the model's enrichment match the labels' enrichment**, and by it
  TabPFN is the better model. Using the incumbent's behaviour as the acceptance test silently
  assumes the incumbent was right; it wasn't.

  Still true and worth keeping: function and localization are not model features, so any agreement
  is recovered unaided. Measured under TabPFN — J up (1.47–5.93), M down (0.00–0.12), E down
  (0.09–0.20), all matching the labels' direction. **One caveat in the other direction: TabPFN
  OVER-depletes G** (predicts 0.04–0.30 where the labels show 0.74 / 0.65, both non-significant).
  `unclassified` is strongly enriched (OR 5.2–7.6) — tested and **not** explained by protein length
  (TabPFN tracks length at ρ −0.62 vs the forest's −0.61, a 0.02 difference), so it is unexplained
  rather than dismissed. Note both models track length far more strongly than the labels do
  (−0.61 vs −0.33): a standing caveat on this axis, forest and TabPFN alike. The activators differ
  exactly as their named top-20 lists do: ONC212 narrow (J and little else), ADEP4 broad (J + K
  transcription + X mobilome + S + V + D), with **X flipping sign** between them. **B chromatin is
  the largest OR (64) and the weakest evidence — 8 members, 7 in the top 10%: quote the count, not
  the ratio.** `unclassified` is tested rather than dropped and is *not* consistent (enriched on Kp,
  depleted on Sa/ONC212). COG category names are read from COGclassifier's bundled
  `resources/cog_func_category.tsv`, never hard-coded.

  `loc_fisher.png` — **the same test on the stage-03 compartments**, and the axis's premise confirmed
  from outside the model: `cytoplasmic_fraction == 1` enriched at **OR 11.5–13.1**, `cytoplasm` at
  7.9–10.3, while signal peptides (0.01–0.11), TM helices (0.02–0.07), cytoplasmic membrane and
  periplasm are all strongly depleted, and **not one of the 67 Kp / 66 Ec β-barrels reaches the top
  10%**. Note **TMbed's residue count out-predicts DeepLocPro's class label** (11.5–13.1 vs
  7.9–10.3) — the two share no machinery, so that is a fair comparison and it backs
  `docs/localization.md`'s advice to prefer the fraction per protein. ADEP4 depletes `OM`/`Ext`
  while ONC212 sits at chance there; since `extracellular` is stage 03's least trustworthy class,
  read that as possibly a fact about the label.

  `interest_panel.png` is secondary: the consortium panel runs OR 0.21–0.82 with every q > 0.2 — a
  trend toward depletion with no power to call it, agreeing in direction with **M** being the most
  depleted category. Highest members are the cytoplasm-facing ones (`secA` 96th percentile, `secB`
  94–95th, `yajC` 89–93rd, `lptB` 79–81st, GyrA/GyrB 60–70th on Sa).
  CLI: `--top-pct 10` · `--top-n` (absolute override) · `--multi` (count every COG letter, not just
  the first) · `--species` · `-q`.

  **Why classification is fine here, despite Jaccard 0.32.** The two activators' binary calls overlap
  at only 0.32, but that is a *threshold* disagreement (ONC212 at 30 µM is looser: 165 hits ADEP4
  misses vs 28 the other way). On *ranking* they agree at AUROC 0.82–0.88, and a probability-emitting
  classifier is scored on ranking. A regression variant was built first and lands in the same
  territory (0.874 / 0.753), so the conclusion does not hinge on the framing — see the run log.

  **Two estimator traps, both measured.** (1) **A forest is wrong for a one-feature baseline** — RF on
  `log_length` scores 0.682 against logistic's 0.776, because 500 trees bin one variable into steps;
  quoting it would have inflated the ESM-C delta to +0.175. (2) **v1's `10f` hard-coded
  `LogisticRegression(C=1.0)`, which is badly under-regularised here** — it pins **52% of
  probabilities at exactly 0.000/1.000** (unrankable ties over half the proteome) and costs accuracy
  (0.813 vs 0.873 at C=0.001). 10f never ran, so neither was ever caught. A **`saturated > 20%`
  guard** now exits the stage.

  **The labels join to NOTHING by identifier.** Conlon is strain COL (`SACOL`/`YP_18xxxx`), Jacques is
  C0673 (`ODV*`), the proteome is NCTC 8325: **0 of 1,943** match by RefSeq, locus tag or Jacques's own
  Mu50 UniProt column. So it is a **DIAMOND sequence join** — within species, nearly lossless:
  **1,873/1,943 (96.4%)** at ≥95% id, median identity 100%. Exits below a 90% floor.

  **The label is v1's audited `*_bin`, and asymmetric on purpose**: ADEP4 = `log2FC ≤ −1 AND padj <
  0.05`, ONC212 = `log2FC ≤ −1` alone (Jacques published no p-values). Unmeasured is NaN, never a
  negative. **Cleavage is deliberately unused** — it reproduces at ρ 0.22 between activators, and
  within one chemistry abundance vs cleavage agree at ρ 0.06. **Two activators, never merged.**

  **The two probability columns are NOT independent evidence.** Their predicted values correlate at
  **ρ 0.89** across all three species (the forest gave 0.83–0.86, so TabPFN made them *less*
  independent, not more) while the underlying labels agree at only ρ 0.52 / Jaccard 0.32 — both models read the same embedding, so the learnable signal is shared. A shortlist
  built on "both activators agree" is closer to one opinion than two. (Flag fires above 0.90.)

  **0.5 is the wrong threshold, and the right one MOVED with the estimator.** TabPFN's probabilities
  span 0.87–0.95 on Kp (the balanced forest topped out near 0.61), so `≥0.5` now selects 220 Kp
  proteins where it selected 20. `src.degradability.hits()` defaults to the **re-derived**
  base-rate-reproducing cuts (**0.328 / 0.313**, replacing the forest's 0.394 / 0.426) — these are
  empirical quantiles of the OOF distribution and are a property of the estimator's calibration, so
  re-derive them on any estimator change. They reproduce the base rate on the LABELED set, not on a
  proteome: applied whole-proteome the same cut takes 14–21% (adep4) / 31–35% (onc212). Better
  still, rank on `_prob` — the model was validated on ranking, not a cut.

  ***E. coli* and *K. pneumoniae* rows are ranking hypotheses, not measurements** — no activated-ClpP
  labels exist for either, so `_hit` is empty throughout. `nn_similarity` +
  `evidence/domain_bands.tsv` price the extrapolation (**Kp 0.814 / Ec 0.817** for ADEP4 under
  TabPFN, up from the forest's 0.752 / 0.754). Read
  cautiously: only two bands populate, and the premise that ESM-C cosine measures transferability is
  unvalidated and untestable without Gram-negative labels.

  **Load through `src/degradability.py`** — `load`, `load_all`, `hits`, `load_labels`,
  `load_seqmap_audit`, `load_cv`, `load_cross_activator`, `load_domain_bands`, `load_oof`,
  `manifest`. Name collides with the frozen `legacy/src/degradability.py`; v2 never imports that one.

  CLI: `--species` · `--activator {adep4,onc212}` · `--folds` · `--seed` · `--limit` · `--refresh` ·
  `--dry-run` · `-q`. ~6 min end to end (5 CV seeds + the cutoff sweep).
  Details, traps and the run log: `docs/degradability.md`.

- **`orthology/orthofinder.py`** — **who is an ortholog of whom, and how similar they are.** Two matrices,
  because v1 shipped only the first and its retrospective calls the result **trap 1**: the kp→ec
  table's `pident`/`coverage`/`bitscore` were *entirely empty*, so "any axis that thresholds transfer
  on percent identity silently drops everything". All four proteomes in one joint run — **33,436
  proteins, 5,064 orthogroups, 10 min**.

  Outputs `data/processed/orthology/`: **`orthologs.tsv`** (sparse, 29,844 pairs called by either
  method — `is_ortholog_orthofinder` and `is_rbh` **side by side, never merged**, plus
  `same_orthogroup`), **`neighbors.tsv`** (sparse, 146,723 rows — the **top-5 nearest neighbours in
  each target species** with `pident`, `ppos`, both coverages and `bitscore_norm`; self-hits
  dropped, so the within-species block is a protein's nearest *paralogs*), and
  **`orthology_<species>.tsv`** (DENSE, one row per protein).

  **A 0 here is a MEASURED 0.** OrthoFinder runs **de novo on our own FASTAs**, so every protein is
  either assigned an orthogroup or named in `Orthogroups_UnassignedGenes.tsv` — and the stage exits
  unless that accounts for each proteome exactly (measured: 33,436/33,436). **This is why no lookup
  database is used**: stage 02's eggNOG OGs are free and already on disk but cover only 91.5% Kp /
  96.3% Ec / 89.1% Sa and no human, so a zero there is an unknown in disguise; OrthoDB is worse on
  this dark anchor. A *sparse* matrix cannot express a zero at all — that is what the dense
  per-protein table is for.

  **OrthoFinder's recall depends on panel size; RBH's does not — know this before quoting a number.**
  v1's much-cited **55.5%** (3,179 Kp proteins with an *E. coli* ortholog) came from a **25-species**
  run. At four species OrthoFinder gives 2,568 (44.8%) — *not a defect*. **RBH is the like-for-like
  comparison and it passes: 3,074 vs v1's 3,003 (+2.4%), median identity 85.7% vs 86.0%.** A
  2-species smoke test made the effect unmistakable: OrthoFinder emitted just 631 ortholog groups.
  Adding tier C is the lever if recall matters more than run time (OrthoFinder 3's
  `--assign/--core` adds species without recomputing; the run dir is kept for that).

  **`--very-sensitive`, not v1's defaults, and the error direction is the reason**: under-detecting
  human homology makes a target look *more selective than it is*. Measured payoff — Kp proteins with
  any human hit **732 (12.8%) → 1,430 (25.0%)**. Identity and coverage are columns, never filters.

  **The spot check found something that belongs in front of the collaboration**: Kp `clpP` → **human
  mitochondrial CLPP at 56.3% identity / 79.2% positives**, an ortholog by *both* methods (and
  `clpX` → human CLPX at 41.1%). The degradation handle itself has a close human ortholog — and
  stage 04's ONC212 is an imipridone whose characterised human target *is* ClpP. `ftsZ`, by
  contrast, has **no human hit at all** — a measured zero.

  **Load through `src/orthology.py`** — `load_orthologs`, `load_neighbors`, `load`, `load_all`,
  `orthologs_of`, `neighbors_of`, `load_orthogroups`, `load_disagreement`, `control`, `manifest`.
  `SPECIES` here is all **four** — human is in scope, unlike `src/function.py`.

  **Two traps.** **OrthoFinder exits 0 when its dependency check fails** — no `diamond` on `PATH`
  means an ERROR block, an empty `Results_` dir and a *success* return code; only the
  results-directory check catches it, so `ensure_diamond()` must run *before* OrthoFinder. And the
  **launcher cannot be called directly**: its `#!/usr/bin/env python3` shebang resolves to the
  unrelated `ersilia` env and dies on `ete4` — which **is installed**; name `<env>/bin/python`
  explicitly rather than chasing the wrong bug.

  CLI: `--species` · `--top-k 5` · `--evalue` · `--sensitivity` · `--threads` · `--refresh` ·
  `--dry-run` · `-q`. Figures: `scripts/plots/orthology.py` (stylia).
  Details, traps and the run log: `docs/orthology.md`.

- **`orthology/orthodb.py`** — **absolute orthologous groups, from OrthoDB v12.2.** OrthoFinder's
  orthogroups are *de novo* and therefore **panel-dependent** — recall rises with the number of
  species in the run. For a grouping that does not move when the query changes, the groups must be
  defined elsewhere: **17,551 bacterial + 5,952 eukaryotic species, 13.0M groups over 990 levels.**
  Writes `orthodb_<species>.tsv` for **all four proteomes**, alongside OrthoFinder's, never instead.

  **v12 REPLACES v11; the two id spaces must never be joined.** OrthoDB's own README says an "OG
  unique id (**not stable and re-used between releases**)", so relabelling a v11 id as v12 would be
  silently wrong. Every row carries `orthodb_version`. The v11 artifacts were deleted, not kept
  alongside.

  Measured at the default `--reps 20`: **Kp 74.6 → 92.7% · Ec 88.8 → 95.5% · Sa 88.3 → 91.8% ·
  human → 96.0%**, and 100% *verdict* coverage is still guaranteed (33,436/33,436 — 8,865
  `assigned_by_uniprot`, 22,903 `assigned_by_sequence`, 1,668 `no_group`).

  **Three changes bought that, each measured.** (1) **Search all of OrthoDB, not one assembly** — v11
  DIAMONDed each species against its own OrthoDB organism, which for Kp is a single 4,975-protein
  strain. (2) **Take groups at every level, not just the domain** — plenty of families have a group at
  Enterobacteriaceae and none at Bacteria; restricting to `at2` left **510 Kp proteins** with a clean
  hit and no group. (3) **But transfer only at levels our own lineage passes through** — a Kp protein
  matching a *Bacillus* gene whose only group is Bacillus-level has learned nothing. Lineages come
  from `level2species` column 4 and live in `src/orthology.ORTHODB_LINEAGE`.

  **The identity floor is decoy-calibrated, not inherited.** 33,436 composition-preserving shuffles
  of our own sequences, searched against the same database, produced **0 hits at 40%/50% and just 2
  at 25%/50%**, while the 40% floor cost ~10 points of coverage on Kp and **~22 on S. aureus** — and
  E. coli correctness is **flat** across floors (73.3% at 40/50 vs 72.6% at 20/50). So the floor is
  **25% identity / 50% coverage**. CLAUDE.md's ≥40% rule is for annotation **transfer**, a stricter
  task — do not conflate them.

  **Filter on `orthodb_confidence`, NOT on identity** — the winning group's share of total bitscore
  at the domain level. On the E. coli control **>0.9 is 98.0% correct and <0.5 only 41.4%**, while
  95–100% identity alone is just 79.0%. But read it **with `orthodb_n_candidate_ogs` and
  `orthodb_match_pident`**: a single candidate takes 100% of the share trivially, and while
  `conf>0.9 & n_cand==1` still scores ~95%, **`n_cand==1` at 25–45% identity drops to 80.0%** — a
  lone candidate found at low identity is the one shape to distrust. Verified independently three
  ways (UniProt xref, group-name-vs-protein-name, and OrthoDB's own `/v12/blast`); see
  `docs/orthology.md`.

  **The honest limit, and why the obvious accuracy number is wrong.** Sequence-tier assignment scores
  **75.4% exact-id agreement** with OrthoDB's own answer on E. coli — but **OrthoDB maintains parallel
  Bacteria-level groups for the same family**, verified against its own API: `5287828at2`
  "Chromosomal replication control, initiator DnaA" and `9807019at2` "chromosomal replication
  initiation protein A" are both real, both Bacteria-level, both dnaA. **57.1% of the apparent errors
  are a parallel group with the same function name**, so corrected agreement is **89.8%**. Do not
  quote the exact-id figure as the error rate, and do not compare two proteins by group id alone —
  check the name.

  **`--reps` is the accuracy knob, and 20 is the measured default.** Going 5 → 20 moved exact-id
  accuracy **64.7 → 75.4%**, name-corrected **84.9 → 89.8%**, and nearly doubled the
  high-confidence calls (664 → 1,178) while raising *their* accuracy **95.5 → 98.0%**. Coverage
  barely moved (Kp 91.1 → 92.7%): the gain is trust, not reach. Costs ~31 min of DIAMOND instead of
  ~13, and an 11.9M-sequence database instead of 5.5M.

  **K. pneumoniae HS11286 is absent from OrthoDB entirely** (taxid 1125630 is not in `species.tab`;
  the only K. pneumoniae organism is `72407_0`), so **all** 5,216 of its assignments are sequence-tier
  and only 1,923 of them clear the 0.9 confidence bar. Low confidence does not mean
  wrong — Kp's `clpP` is correctly `9802800at2`, agreeing with all three species *and* with v11, at
  confidence 0.219 — it means parallel groups competed.

  **Accession joins into OrthoDB still do not work** — the fourth time this project has hit that wall.
  UniProt's `xref_orthodb`: Ec 96.4%, Sa 99.4%, **Kp 0.0%** (and *no* Klebsiella proteome carries it).
  OrthoDB's own mapped UniProt column is worse: Kp 0.0%, Ec 73.6%, Sa 93.2%, **human 17.0%**.

  **OrthoDB has NO root level spanning domains** — a bacterial group is `<n>at2`, human's `<n>at2759`,
  and they are **not comparable**. Cross-domain similarity is `neighbors.tsv`, not this table.

  The small tables are kept (135 MB); `genes`, `OG2genes` (streamed **twice**) and the 37 GB
  `aa_fasta` are **stream-filtered and never stored** — ~46 GB at a measured 16.6 MB/s, so a cold run
  is ~50 min of transfer plus ~31 min of DIAMOND at `--reps 20`.
  **Load through `src/orthology.py`** — `load_orthodb`, `load_orthodb_all`, `load_orthodb_long`.
  CLI: `--species` · `--reps` · `--sensitivity` · `--max-targets` · `--threads` · `--refresh` ·
  `--dry-run` · `-q`. Details and the run log: `docs/orthology.md`.

- **`ligands/chembl.py`** — **can a small molecule bind this protein?** The first
  ligandability track, from measured bioactivity in **ChEMBL 37**. Maps the three bacterial
  proteomes onto ChEMBL's target sequences with DIAMOND and counts the **non-redundant** ligands
  reachable at four nested distances. ~2 min with the dump; **a cached re-run is 0.3 min and needs
  no database at all.**

  Measured, first full run — proteins with a potent ligand (pChEMBL ≥ 6): **Kp 113 · Ec 96 · Sa 78**
  at `remote`, of which **21 / 65 / 32** have one on a ≥95%-identical target. **That is ~2% of each
  proteome, and it must not be forced upward** — it is a fact about how little of the bacterial
  proteome anyone has screened. v1's numbers (175 Kp / 155 Ec) reproduce within tolerance and in the
  right direction: v2 is *lower* because it adds `standard_relation='='`, `data_validity_comment IS
  NULL`, `potential_duplicate=0` and a subject-coverage floor v1 lacked.

  **Homology transfer is the whole game.** ChEMBL holds **21** Kp, **225** Ec and **80** Sa
  single-protein targets against proteomes of 5,728 / 4,403 / 2,889. v1's validation case reproduces
  exactly: HS11286 `A0A0H3H184` → ChEMBL `Q93LQ9` (Kp β-lactamase) at **100% identity**, invisible
  to accession matching.

  **Same-species matching is by organism NAME, not `tax_id`** — ChEMBL files strains under their own
  taxids, so `tax_id=562` finds **65** E. coli targets while `organism LIKE 'Escherichia coli%'`
  finds **225** (K-12 lives at 83333). Two-word binomial, so `Klebsiella aerogenes` is not swept in.

  **One `confidence_score` gate cannot serve both target types**, and this is load-bearing. `≥8`
  removes **0 of 3,271,336** single-protein rows (subsumed by requiring `pchembl_value` — keep it,
  but do not call it quality work), and returns **exactly zero protein complexes**, which would
  have shipped empty `complex_*` columns reading as a real biological zero. **DNA gyrase is a
  `PROTEIN COMPLEX` in ChEMBL** (Ec 713 compounds, Sa 491; topoisomerase IV 200; Mtb ClpP1P2 15) —
  GyrA/GyrB are in `src/interest.py`, and v1's SINGLE-PROTEIN-only rule made them look unliganded.
  Complexes run at `≥6`; `PROTEIN FAMILY`/`PROTEIN COMPLEX GROUP` (confidence 4–5, target ambiguous
  within the group) are **measured and excluded**, not silently absent.

  **Non-redundant means scaffolds.** v1 counted raw `molregno`, so salts double-counted. Here
  compounds collapse to the `molecule_hierarchy` parent and diversity is distinct **Bemis–Murcko
  generic scaffolds** — median **2.0 compounds per scaffold**, and Ec `folA` is **443 compounds but
  64 scaffolds**, `lpxC` 636 → 150, `ampC` 12,671 → 5,829. Acyclic compounds share one empty-string
  bucket; never dropped, never zero. Counts are **unions over the bucket's target pool**, not v1's
  single best hit.

  **Buckets are nested and restricted to true Bacteria**: `direct` ≥95% ⊆ `close` ≥60% ⊆ `remote`
  ≥40%, plus `species` (same organism name) which `close` is a superset of by construction. The
  restriction is v1's documented fix — an unrestricted "non-human" bucket with no identity floor
  gave **424** potent Kp proteins, rat `P97697` winning bacterial slots at ~30% id, against a true
  175. Nothing is resolved silently: `allorg_n_compounds`/`_n_scaffolds` ship the literal
  all-organism number, and human is its own liability block, **never merged into a bucket**. The
  **60% band is the one arbitrary number and turns out not to matter** — the identity distribution
  of liganded hits is bimodal (a mass at 40–45%, a spike at 95–100%, a trough between).

  **The selectivity finding belongs in front of the collaboration**: `clpP` carries **106 human
  compounds against 61 bacterial**. Stage 05 measures Kp clpP → human mitochondrial CLPP at 56.3%,
  and stage 04's ONC212 is an imipridone whose characterised human target *is* ClpP — so this axis
  sees the same liability from an independent direction. `ftsZ` is the control: **0 human**, real
  bacterial evidence, matching stage 05's measured zero.

  **A free cross-check that comes out clean**: `human_best_pident` (DIAMOND vs ChEMBL sequences)
  against stage 05's `neighbors_of(ac, "human")` (DIAMOND vs the human proteome) — different
  databases, same quantity, **r 0.979–0.988, median |diff| 0.0 pp**.

  **`ligands/precedents.py` + `src/precedents.py` — LIGAND PRECEDENT FOR ANY SEQUENCE.** A query
  tool, not a proteome stage: give it a sequence and it returns three counts in about a second.

  | | |
  |---|---|
  | `n_ligands_exact` | ligands on an EXACT match — UniProt accession, or identical sequence |
  | `n_ligands_bacteria` | **UNIQUE** ligands across bacterial targets passing the identity + both coverage floors |
  | `n_ligands_human` | the same over human targets — **a LIABILITY, never added to the bacterial count** |

  **It needs no database.** Three cached extracts totalling 82 MB —
  `scratch/chembl_targets.faa` (8,469 sequences, headers = `component_id`),
  `chembl_targets.tsv` (9,347 targets with `organism` + `superkingdom`) and
  `chembl_ligands.tsv` (2,591,526 rows) — against `chembl.py`'s 30.5 GB dump. The join is
  `sequence → DIAMOND → component_id → tid → parent_molregno`. Composition, measured: **1,305,242
  distinct compounds**, 110,019 on Bacteria targets and 1,101,363 on human.

  **"UNIQUE" is the whole point and it is not a per-target sum.** `parent_molregno` is ChEMBL's
  `molecule_hierarchy` parent, so counting DISTINCT values over the union of every homologous
  target counts MOLECULES — one compound tested against three homologs counts once. Summing
  per-target counts inflates it, worse the wider the band.

  **The complex track is reported separately and is NOT inside the bacterial count.** Measured:
  E. coli `gyrB` carries **666 single-protein ligands against 1,412 complex** ones, and Kp
  `A0A0H3H0Y6` (gyrA) carries **1,410 complex against 131 single**. DNA gyrase is a `PROTEIN COMPLEX` in
  ChEMBL, and v1 dropped that track and made GyrA/GyrB look unliganded.

  **Validated against the existing axis at matched semantics: 99.9% agreement on all three
  species, 0 chembl-only** (`evidence/precedent_control.tsv`, written every `--species` run). The
  two disagree at *default* settings for two deliberate reasons — the tool excludes the complex
  track, and counts every potency-measurable ligand rather than only pChEMBL ≥ 6 (66 of 67
  default-only Kp proteins sit below 6). Put both back and they reconcile.

  **A case measured and settled, so nobody re-investigates it.** Kp `A0A0H3GWM6` is **99.2%
  identical** to E. coli `P0ADG7` (32 compounds, pChEMBL 8.82) and correctly gets nothing: the hit
  exists (component 1947 = CHEMBL3630) at `qcov 100.0 / scov 26.6` — a **130-aa fragment** against
  a **488-aa** IMP dehydrogenase. `MIN_QCOV`/`MIN_SCOV = 50` rejected it. **Identity alone cannot
  transfer a ligand count**, and both floors are imported from `src/ligandability.py` rather than
  restated.

  Spot checks: `folA` 517 exact / 595 bacterial / 0 human · `clpP` 30 / 136 / **210 human** at 56.3%
  identity, the selectivity liability this axis already documents.
  CLI: `--sequence` · `--fasta` · `--accession` · `--species` (batch → `precedents_<sp>.tsv`,
  complete and canonical) · `--min-identity 40` · `--min-pchembl` · `-q`. DIAMOND from
  `gradi-ortho` via `GRADI_DIAMOND_BIN`. **Load through `src/ligandability.py`** —
  `load_precedents`.

  **`ligands/bindingdb.py` is a MEASUREMENT, not a deliverable.** v1 shipped BindingDB as a
  co-equal track (93 potent Kp proteins) **without ever measuring the overlap**. Measured here —
  cheaply, because BindingDB carries the target chain sequence *and* the ligand InChIKey inline:
  **33 proteins across all three species gain their first potent ligand, 11.5% over ChEMBL's 287**.
  Real but modest; ChEMBL 37 already ingests 2.68M BindingDB patent activities. **Not promoted** —
  if it ever is, columns go *beside* the ChEMBL ones, never merged.

  **`pchembl_value` only exists for `=` relations on IC50/EC50/Ki/Kd/Potency in nM**, so MIC and
  %-inhibition are absent by construction: **this axis does not say "has an antibiotic"**, and it is
  why the ribosome — whose drugs are measured as MIC — is largely missing from the complex track,
  which matters given stage 04's top-100 is ribosome-heavy.

  Needs **rdkit** (PyPI `macosx_11_0_arm64` wheel — must not flip `gradi` to osx-64, which would
  take ESM-C down with it). The 30.5 GB dump is deleted after each run:
  `data/raw/other/chembl/SOURCE.md` has the recovery procedure (see the ChEMBL exception below). **Never upload it to eosvc.**

  **The archive's integrity is VERIFIED, and the version is now ASSERTED rather than declared.**
  `chembl_37_sqlite.tar.gz` matches EBI's published sha256
  (`33c2037405…`, `releases/chembl_37/checksums.txt`) and its `Content-Length` byte-for-byte — so
  the download is the genuine complete release, not a truncated transfer. That had never been
  checked: the first run recorded byte counts only. Two gaps it exposed, both closed. (1)
  **`CHEMBL_VERSION = "37"` is a bare literal while `find_db` globs `chembl_*.db`**, so a different
  release extracted beside this one would have been consumed silently and labelled 37 in every
  manifest — `assert_version()` now reads the dump's own `version` table (filename as fallback) and
  exits non-zero on a mismatch. (2) **`SOURCE.md` pointed at the FTP `latest/` path**, which moves
  with every release; the pinned `releases/chembl_37/` URL replaces it.

  **The cached extracts were cross-checked against live ChEMBL**, an independent route into the
  same release — 4 targets under this stage's own predicate: `CHEMBL1293248` 24,681 activities and
  `CHEMBL2390811` 8 activities / 8 compounds both agree **exactly**, and the only two deltas
  (`CHEMBL5465386` +18, `CHEMBL2026` +26) are the `potential_duplicate` rows the extract drops on
  purpose. So the extracts are faithful, not merely assumed so.

  **Load through `src/ligandability.py`** — `load`, `load_all`, `load_targets`, `load_ligands`,
  `load_hits`, `load_scaffolds`, `load_cutoff_sensitivity`, `control`, `manifest`, plus
  `evidence_level` (the tightest bucket with a potent ligand) and `selectivity_risk`.

  CLI: `--species` · `--pchembl` · `--threads` · `--limit` · `--refresh` · `--dry-run` · `-q`.
  Figures: `scripts/plots/ligands.py` (stylia).
  Details, traps and the run log: `docs/ligands.md`.

- **`essentiality/labels.py`** + **`essentiality/deg_proteomes.py`** — **the essentiality training
  corpus.** Stage 07 predicts essentiality with **ProteomeLM-Ess** (Malbranke, Zalaffi & Bitbol, PNAS
  2026, `10.1073/pnas.2524201123`; papers in `docs/papers/`). These two scripts build its labels; the
  predictor itself is the next part.

  **OGEE v3's SERVER is gone, but the bulk file is RECOVERABLE — corrected 2026-09-18.** An earlier
  version of this file said the dataset was lost. The server diagnosis below is right and still
  worth keeping; the conclusion drawn from it was wrong, because nobody checked the Internet
  Archive. Two independent mirrors return byte-identical content (md5
  `b42f4a3358484b490e8bffe8c90edd00`, 1,151,931 B), both verified:

      https://web.archive.org/web/20250620231412id_/https://v3.ogee.info/static/files/gene_essentiality.txt.gz
      https://raw.githubusercontent.com/ThomasBeder/OGEE_ID_conversion/HEAD/Essential_gene_inormation.tar.gz

  The `id_` suffix on the Wayback URL is **mandatory** — without it you get the HTML wrapper, not
  the file. Staged at `data/source/ogee/`. CC BY 3.0.

  **Measured contents**: 255,162 rows x 8 columns (`dataset · taxaID · locus · gene · score ·
  essentiality · pmid · Ref_db`), 89 taxa, 124 datasets, `essentiality` in {E, NE, C, ND, ...}.
  Against the DEG-derived corpus (173,048 proteins, 20,194 essential, 38 species) that is roughly
  **2x the species and 2.2x the positives**, and it ships the negative class explicitly instead of
  reconstructing it from NCBI proteomes. **E. coli 16,979 rows / 1,523 essential; S. aureus 5,612 /
  671 — and *K. pneumoniae* ZERO**, the same structural gap DEG has. So this widens the corpus but
  does NOT close the anchor's hole; only the Tn-seq/TraDIS screens in
  `scripts/essentiality/screens.py` do that.

  **The lesson worth keeping: a dead server is not a lost dataset.** Check the Wayback Machine with
  `id_`, and check for a GitHub mirror, before recording anything as unobtainable.

  The original (correct) server diagnosis: `v3.ogee.info`
  completes TCP then aborts the TLS handshake. Verified against LibreSSL 3.3.6, OpenSSL 3.5.6 **and
  real Chrome** (`ERR_SSL_PROTOCOL_ERROR`), with SNI, without SNI, pinned to TLS 1.2, ALPN disabled,
  forced http/1.1, by direct IP; plain HTTP 308-redirects into the same endpoint. Server-side, so no
  client-side workaround and no browser trick helps. The OGEE paper names no FTP, mirror or deposit;
  Database Commons lists none; there is no GitHub copy (human subsets only) and no OGEE v4. v1 hit
  this in May (`legacy/HISTORY.md` trap 29). Maintainer for a manual ask: Weihua Chen, HUST.
  **So the corpus is DEG, OGEE's own prokaryotic upstream.** Reproduce the method exactly; substitute
  the corpus; keep the label source swappable.

  **`essentiality/labels.py`** fetches DEG's three bulk files and builds the label table with
  **every dataset and every column preserved** — 66 datasets, 42 species, 26,619 essential genes,
  each with **its protein sequence** (`DEG10.aa`, 1:1 by DEG gene id). The sequences are the point:
  the join downstream is **by sequence**, because DEG's identifier columns are patchy — gene symbol
  63.1%, GI 54.5%, COG 48.4%, **UniProt AC 45.7%**, locus tag 42.9%.

  **Exclusions are FLAGGED, never dropped** (`retained` says what training may use): 4 non-genome-wide
  methods (antisense RNA, MATT, insertion-duplication, transposon-hybridisation — "absent from the
  list" cannot mean non-essential there) and 11 condition-specific screens (tobramycin, murine
  pneumonia, bile, cholesterol, kanamycin — these are why *P. aeruginosa* PAO1 swings **117 → 336 →
  551** essential genes across its three datasets). **51 datasets retained over 38 species.**
  Consolidation to one label set per species is deliberately deferred — *S. aureus* appears 7×,
  *P. aeruginosa* 4×, *Salmonella* 4×, *E. coli* 4× — so the union-vs-intersection rule is chosen
  against measured counts rather than guessed.

  **Four traps in DEG, all measured.** (1) **326 of the 26,619 "sequences" are the literal string
  `Not available now.`** — the join to `DEG10.aa` succeeds, so a coverage check reads 100% while the
  payload is junk; only validating the amino-acid alphabet catches it, and *E. coli* O157:H7 is
  **14.1%** placeholder. A further 52 carry an **internal `*` stop** (frameshifted translations,
  mostly *N. gonorrhoeae*); a *trailing* stop is stripped, not rejected. Short proteins are NOT
  filtered — `rpmJ` is a real 37-aa essential gene. (2) **`DEG1058` (*S. suis*) is an in-vivo pig
  infection screen** whose condition string reads "Columbia blood base agar"; its paper is titled
  *"…**conditionally** essential genes for S. suis infection in pigs"* and every gene carries a note
  `Recovered from blood/cerebrospinal fluid/meninges`. Caught by reading DEG's **per-gene notes**, not
  its condition field — so the filter is data-driven. (3) **DEG's index disagrees with its own shipped
  rows** on 2 of 66 datasets (*S. oneidensis* 403 vs 402; *R. palustris* 522 vs **552**); printed
  every run. (4) Two columns are undocumented: col 10 packs `locus_tag:X;gi:N`, col 13 is the per-gene
  note that exposed trap 2.

  **`essentiality/deg_proteomes.py`** supplies the **negative class**, which DEG does not ship (38 of
  66 datasets report a non-essential *count* and no genes). It fetches each retained dataset's
  complete proteome from the RefSeq/EMBL replicon accessions DEG records — NCBI
  `efetch db=nuccore rettype=fasta_cds_aa`, which also yields locus tags — then labels it: essential =
  matched a DEG positive, non-essential = the rest. ProteomeLM needs the whole proteome anyway, since
  it contextualises every protein against every other.

  **Exact sequence matching is NOT enough, and the per-dataset rate is the control.** DEG's vintage
  annotation has drifted from the current replicon: exact match recovers only **59–95%** (median
  84.8%), and the unmatched are same-length point-substituted, not a start-codon or prefix artefact —
  measured, and `locus_tag` is **0%** on the older datasets so there is no identifier fallback. With
  **DIAMOND ≥95% identity** (the house threshold, borrowed from `gradi-ortho` via `GRADI_DIAMOND_BIN`)
  the median rises to **99.5%**, recovering +1,320 positives at median 100% identity. Without it,
  hundreds of genuinely essential genes would be silently labeled non-essential.

  **Do not scan for accessions with a `\b`-anchored regex.** DEG's replicon field mixes `,`, `;`,
  `, ` and bare spaces, and a scanning regex silently drops **`NZ_`-prefixed RefSeq accessions**
  because `_` is a word character so the boundary never fires. That cost *E. coli* O157:H7 — the
  largest single dataset, 1,071 positives — on the first run. Split on delimiters, then full-match
  each token.

  Measured: **49 of 51 datasets, 173,048 labeled proteins, 20,194 essential (11.67%), 38 species** —
  comparable in scale to the paper's 213,608 labels over 83 genomes, and bacteria-only. Three known
  gaps, all flagged in `proteome_join.tsv` rather than hidden: `DEG1003` *V. cholerae* joins at only
  **74.0%** (real annotation drift, below the 80% floor), `DEG1037` *S. pyogenes* MGAS5448 records no
  replicon at all, and `DEG1053` *B. cenocepacia* K56-2 records a WGS **master** accession for which
  `efetch` returns zero CDS.

  **Spot checks pass**: `dnaA`, `gyrB`, `rpoB`, `ftsZ`, `murC`, `mraY`, `rplF` all essential in
  *E. coli* MG1655; `lacZ` and `araB` not. **Watch the base rate before training** — it spans 15×
  across species, from *M. genitalium* **72.7%** essential (a minimal genome) to *M. avium* **4.7%**.

  Outputs `data/processed/essentiality/evidence/`: `deg_datasets.tsv` (66 × 22),
  `deg_genes.tsv` (26,619 × 28), `proteome_join.tsv` (51 × 14, the join control),
  `labeled_proteins.tsv` (173,048 × 11). Raw + provenance: `data/source/deg/SOURCE.md`
  and `proteomes/<deg_dataset_id>.faa`.

  **Useful for any gated supplement**: `https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/supplementaryFiles`
  returns a ZIP of the publisher's own files and **bypasses the Cloudflare 403s** on ASM, PNAS, Nature
  and Cell — better than v1's authenticated-Chrome route for open-access papers.

  CLI (both): `--refresh` · `--dry-run` · `-q`; part 2 adds `--limit N` (smoke test, writes only
  `scratch/smoke_*`). ~1 min and ~15 min respectively; both cache and re-run cheaply.

- **`essentiality/geptop.py`** — **can the organism live without this protein?** Geptop 2.0
  (Wen 2019, the latest release), **ported faithfully to Python 3**: for each protein, sum `1/d`
  over the 37 reference prokaryotes where it has a **reciprocal best hit** to a gene DEG lists
  essential, with `d` the **k=6 composition-vector (CV-tree) distance**. Writes
  `geptop_<species>.tsv` for the three bacteria.

  **The upstream cannot be run as shipped.** `github.com/RiversDong/geptop` is a 25 MB `.rar`
  containing Python-2 code, and it has three genuine defects beyond the language: an unwaited
  `os.popen(makeblastdb)` racing `blastp`, broken `except (Exception, err)` handlers, and a used
  `import pp` (Parallel Python, py2-only, dead). Those three are fixed and **nothing else** — the
  arithmetic is transcribed, including two quirks preserved on purpose (an unobserved 6-mer scores
  exactly `−1`; `Distance()` accumulates `Q` with the same product as `O`). v1's version substituted
  median RBH % identity for the CV distance; that was a real deviation and is not repeated.

  **The phylogeny term validates against taxonomy unaided** — same strain 0.064, same species 0.092,
  same family 0.300, same phylum 0.486, cross-phylum 0.496–0.499 — which is the best evidence the
  port is faithful. Distances **saturate near 0.5**, so `1/d` only spans ~2.0–3.3 across unrelated
  references.

  **Accuracy is 0.59–0.84 and the paper quotes the mean — do not quote 0.84.** Validated on DEG
  species that are not Geptop references, excluded at **genus** level (stricter than the paper's own
  leave-one-out, which scores organisms that ARE in its reference set): *R. solanacearum* **0.810**
  (0.841 on scored proteins) reproduces the published figure; *H. influenzae* manages only **0.587**.
  Its label set is the likely cause — 36.9% essential is an extreme outlier and it has MORE
  essential RBHs yet worse AUROC. **The Kp prediction inherits that variance.**

  **A score of 0 means two different things — read `geptop_evidence`, never the score alone.** On
  Kp: `essential_orthologs` 1,929 (33.7%), `orthologs_none_essential` **3,350 (58.5%) — a confident
  NON-essential call, evidence not absence**, `no_orthologs` 449 (7.8%) — the only true gap. So
  coverage is **92.2%**, not 33.7%; an earlier version of this stage defined `informative` as
  `n_essential_rbh > 0` and understated it by 58 points. The zero rows are **tied**: a zero means
  *unranked*, not "low rank".

  **Two traps.** `geptop_score` is **proteome-relative** (min-max within each species, so 0.6 in Kp
  ≠ 0.6 in Sa — use `geptop_score_raw`), and **two of our anchors ARE references**: `DataSet4` is
  *E. coli* MG1655, `DataSet13` is *S. aureus* NCTC 8325, and the self-contribution is **measured, not estimated: 53.2% of E. coli's total
  reference weight and 58.3% of S. aureus's comes from the organism ITSELF** (d 0.0113 / 0.0083,
  weight 88.7 / 121.0), against **3.6%** for Kp. `d` is not exactly 0 — our UniProt proteomes differ
  slightly from Geptop's RefSeq builds — so `if d == 0: weight = 100` never fires, yet `1/d` exceeds
  it anyway. Leave-one-out is
  deliberately **not** applied — DEG supplies a measured label for those two, and **Kp is absent
  from all 37** (weights 2.00–2.86, no self-weight). Every row carries `geptop_in_reference_set`.

  Spot checks (Kp): `ftsZ` **1.000, the highest-scoring protein of 5,728**; `gyrB` 0.971, `rpoB`
  0.943, `dnaA` 0.845, `secA` 0.729 all >98th percentile; `lacZ`/`araB`/`fadB` exactly 0. `clpP`
  0.231 (92.9th) sits **below** the 0.24 cutoff — correct, and worth knowing since ClpP is the
  degradation handle.

  **`reference_audit.tsv` is EVIDENCE and lives in `evidence/`, not in the run directory.** It
  used to be written under `scratch/geptop/<species>/`, where a routine `scratch/` purge deleted
  all three copies and cost a 90-minute re-run — the deliverables survived, but the only record of
  which references contributed did not. `scratch/` is documented as safe to delete, so anything
  that is cited cannot live there; the directory contract's own test settles it ("would you cite or
  check it" → `evidence/`). Now `evidence/geptop_reference_audit_<species>.tsv`.

  **`blastp`/`makeblastdb` are borrowed from `gradi-prokka`** (`GRADI_BLAST_BIN` overrides) — the
  `GRADI_RPSBLAST_BIN` pattern; do not install blast into `gradi`. **Human is excluded by
  construction**: all 37 references are prokaryotes, so there is no honest score, exactly as for
  COG. ~22 s per blastp run × 74 per proteome ≈ 30 min/species; composition vectors ~3 min once for
  all 37. BLAST output and CVs are both cached, so a re-run is a re-parse.
  **Load through `src/essentiality.py`** — `load_geptop`, `load_geptop_all`,
  `geptop_reference_audit`, `geptop_informative_only`, `geptop_scored_only`, `load_labels`.
  CLI: `--species` · `--validate DEG_SPECIES ...` · `--threads` · `--cutoff` · `--cv-jobs` ·
  `--refresh` · `--dry-run` · `-q`. Details: `docs/essentiality.md`.

- **`essentiality/ogee.py`** + **`ogee_proteomes.py`** + **`ogee_dataset.py`** — **a third
  opinion on essentiality, learned from 26 measured prokaryotic proteomes.** OGEE v3 supplies
  labels and no sequences; these three scripts characterise the corpus, attach sequences, and
  train. Output `ogee_<species>.tsv`: `ogee_ess` (0–1, never null) + `ogee_evidence` (the MEASURED
  OGEE label for that exact protein — 1, 0, or empty).

  **HALF OF OGEE IS UNUSABLE, and the reason is one we had already recorded from the other side.**
  Of 87 taxa with decided E/NE calls, **40 have ZERO negatives** — 17,743 positives-only entries —
  and **25 of those 40 come from one paper, PMID 29769716, Price 2018, the Fitness Browser
  RB-TnSeq collection.** RB-TnSeq cannot see essential genes by construction (no insertions survive,
  so they are absent rather than extreme), which makes its essential list positives-only by the
  nature of the assay. That inflates OGEE's overall base rate to **0.310** against 0.117 for our
  DEG corpus. Applying the project's own both-classes + prokaryote filters leaves **26 taxa /
  78,893 proteins / base 0.173**.

  **Leave-species-out, not the paper's single holdout** — the project owner's instruction, and the
  better fit: K. pneumoniae is absent from OGEE entirely, so the question is "given N measured
  bacteria, how well can we call the N+1th?", which is exactly what a held-out taxon measures. It
  also gives 26 estimates instead of 1, and **the spread is the result**: AUROC **0.529–0.940**,
  mean 0.782, median 0.805 (`evidence/ogee_leave_species_out.tsv`).

  **THE HEADLINE FINDING: `corr(base_rate, AUROC) = −0.643`.** Screens that call many genes
  essential are much harder to predict, and two same-organism pairs make it unarguable —
  *P. aeruginosa* PAO1 (base 0.081) scores **0.837** while PA14 (base 0.302) scores **0.634**;
  *Salmonella* Typhi CT18 (0.105) scores 0.782 while Typhimurium SL1344 (0.404) scores **0.529**,
  barely above chance. A screen calling 40% of genes essential is measuring fitness defect, not
  essentiality. **This is what prices the Kp column**: our three measured Kp screens run base
  0.075–0.106, squarely in the band where held-out taxa score **0.78–0.93**.

  **ONE model, fit on all 26 taxa, used for all three anchors** (owner's instruction). The caveat
  that creates is carried as DATA, per protein, in `ogee_evidence`: E. coli K-12 and S. aureus
  NCTC 8325 ARE OGEE taxa AND our anchors — **100.0% and 97.5% of their corpus proteins are
  literally the same SEQUENCE as an anchor protein** — so where `ogee_evidence` is non-null the
  model was fitted on that protein with that label and `ogee_ess` is closer to recall than
  prediction. Measured coverage: **Ec 4,193/4,403 (95.2%), Sa 2,815/2,889 (97.4%), Kp 0/5,728.**
  The effect is visible in the output: same model, top-decile cut **0.861 on E. coli against
  0.471 on Kp**. So `ogee_ess` is comparable WITHIN a species, never across.

  **Three traps.** (1) `locus` is 100% populated but heterogeneous, and the namespace is a property
  of the (taxon, dataset) BLOCK — pure at median 1.000, and 87 of 89 taxa use one throughout.
  **E. coli K-12 is one of the two exceptions and carries TWO DISJOINT namespaces** (PEC gene
  symbols, and b-numbers, zero shared ids), so a naive `(taxid, locus)` dedup counts every gene
  twice — 9,496 "genes" for a 4,403-gene organism. Dedup only AFTER resolving to a protein.
  (2) **A missing underscore cost three taxa entirely** — OGEE writes `HI0001`/`HP0001`/`MPN001`
  where the assembly writes `HI_0001`/`HP_0001`/`MPN_001`, scoring 0.000 until a
  punctuation-insensitive fallback was added (H. influenzae → 0.960, H. pylori → 0.994).
  (3) **GCA and GCF are not interchangeable here** — Synechococcus keys on RefSeq `SYNPCC7942_RS*`
  tags that only the GCF annotation carries, scoring 0.226 on GenBank and **0.969** on RefSeq. So
  the assembly is chosen **by measured join rate** across candidates, not by rule; per-candidate
  rates are in `evidence/ogee_proteome_join.tsv`.

  **`--seeds 1` is deliberate here**, against the axis-wide 5. The folds are FIXED (each fold is one
  taxon), so there is no partition randomness to average — only the forest's own, and per-seed SDs
  on this axis run ≤0.004. It is also 26 fits instead of 130: measured ~4 min per fit on
  76,000 × 1,152, i.e. 1.7 h against 8–9. The SD prints `(1 seed)` and stores `NA`, never
  `+/-0.0000`.

  **Load through `src/essentiality.py`** — `load_ogee`. CLI: `ogee.py --rule {any,all,majority}` ·
  `ogee_proteomes.py --candidates N` · `ogee_dataset.py --write-fasta | --embed-plan |
  --score-anchor SPECIES`. Provenance: `data/source/ogee/SOURCE.md`.

- **`essentiality/merge.py`** — **the headline summary**: `essentiality_<species>.tsv`, one column
  per evidence source. It also writes **`deg_<species>.tsv`**, which carries `deg_ess` and its four
  provenance columns; those used to be inlined into the headline and were split out so every
  evidence source has the same shape (its own file + one summary column). **`geptop_ess`** (0–1 continuous, never null — Geptop's prediction)
  and **`deg_ess`** (0 / 0.5 / 1 / **null** — the fraction of DEG screens on our exact strain calling
  it essential; null means unmeasured). Coverage: Kp `deg_ess` **null for all 5,728** (*Klebsiella*
  is absent from DEG), Ec 96.6% measured, Sa 92.7%.

  **Two screens of the same strain disagree, and the column records it rather than hiding it.**
  E. coli has DEG1018 (genetic footprinting) + DEG1019 (Keio knockout); *S. aureus* has DEG1017
  (TMDH) + DEG1061 (Tn-seq). On E. coli **490 of 695 essential calls rest on ONE screen**, only 205
  on both — so `--rule any` gives 695 and `--rule all` gives 205, a **3.4× spread from one choice**.
  `deg_ess = 0.5` is exactly that disagreement. Its three values are **incidental** (two screens
  today); a third screen makes it quarters.

  **`geptop_ess` and `deg_ess` are comparable for RANKING, not as VALUES.** Measured on out-of-set
  labels, `geptop_ess` is monotonic and near-calibrated mid-range (0.30–0.50 → 49.6% essential) but
  a predicted **1.0 means ~70%** and a predicted **0.0 means 3.1%**, where a measured 1.0 means it
  *is* essential. The convenience column `essentiality` merges them and therefore **mixes units** —
  measurements pinned to 1.0/0.0 outrank every prediction. **Rank within a species, never across**:
  Ec/Sa are ~95% measured while Kp is 100% predicted.

  Joined **by sequence, not accession** (99.4% Ec / 96.9% Sa exact match). **Load through
  `src/essentiality.py`** — `load`, `load_all`. CLI: `--species` · `--rule {any,all}` · `-q`.

  **This axis ships TWO tables at the task root, deliberately** — `essentiality_<species>.tsv` and
  `geptop_<species>.tsv`. The second is the raw Geptop output (`geptop_score`, `geptop_score_raw`,
  `geptop_essential`, `geptop_n_rbh`, `geptop_n_essential_rbh`, + evidence/informative/reference
  flags); it reads like evidence but the project owner chose it as a deliverable. Do not demote it
  to `evidence/` in a tidy-up.
  Details: `docs/essentiality.md`.

- **The essentiality axis has FOUR tiers, and each answers one question.** Restructured
  2026-09-21 on the project owner's instruction, because the old layout put the clean training sets
  among 40-odd audit tables and had no machine-readable record of which datasets were used.

  ```
  data/{raw,source}/**/SOURCE.md          34 dirs: what each dataset IS and whether v2 uses it
  data/processed/essentiality/
    dataset_registry.tsv                  47 rows: every dataset found + its DISPOSITION
    training_sets/                        10 files: the clean ML-ready sets, uniform schema
    <source>_<species>.tsv                one table per EVIDENCE SOURCE, per species
    essentiality_<species>.tsv            the headline summary, one column per source
    evidence/ · scratch/                  audits · caches
  ```

  **`dataset_registry.tsv` is the answer to "which of these did we train on, and why not the
  rest?"** `disposition` ∈ `training_set` (10) · `held_back` (30) · `refuted` (6) ·
  `corpus_only` (1), each with its `reason`. **GENERATED by `essentiality/registry.py`** from
  `screens.py`'s own `SCREENS`/`HELD_BACK`/`UNJOINABLE`, `screen_join_audit.tsv`,
  `deg_datasets.tsv` and `ogee_taxa.tsv` — a hand-maintained second copy is how a register stops
  agreeing with the code. It **exits non-zero unless it reconciles with `training_sets/` in both
  directions**: every declared training set must exist as a file, and every file must have a row.

  **Source directories are DECLARED, never inferred** — `Screen.source_dir` in `screens.py`. An
  earlier version guessed by substring and matched `essential_ecoli_bw25113_tradis_goodall` to
  `data/source/go`, because "go" is inside "goodall". A false provenance link reads exactly like a
  correct one, which is worse than a blank.

  **`training_sets/*.tsv` share one schema: `key · label · source_id · features_from`.** The key
  column used to be called `uniprot_ac` and that was **actively misleading**: for five of the nine
  screens it holds `CCN31837.1` (ECL8 EMBL), `WP_038431262.1` (KPPR1 RefSeq), `KPNRH_00001` (a
  locus tag) or `lcl|HG941718.1_prot_...` (a DEG FASTA header). Measured, those five key onto their
  species' anchor proteome at **0 of 4,930 / 4,809 / 4,981 / 4,981 / 5,433** — a join on accession
  returns an empty frame, not an error. `features_from` names the proteome whose embedding matrix
  the key indexes.

  **One file per evidence source, per species**, each complete and canonical:
  `geptop_<sp>.tsv` (prediction) · `deg_<sp>.tsv` (measured, null where unmeasured) ·
  `ogee_<sp>.tsv` (prediction) · `screens_<sp>.tsv` (9 predicted columns). The headline
  `essentiality_<sp>.tsv` carries **one column per source plus the merge** — the relationship
  `geptop_<sp>.tsv` always had with it, now applied to every source.

  **`screens_<sp>.tsv` IS PREDICTIONS THROUGHOUT.** Transferring the measured labels was
  considered and rejected on the owner's instruction: exact-sequence transfer recovers only
  **13.7–71.6%** of rows and **29.8–61.4%** of positives, so it would silently mislabel ~1,700
  measured essentials as non-essential. The model is transferred instead, giving one comparable
  0–1 scale. Where an anchor protein WAS in a screen's training set (the four b-number E. coli
  screens overlap the E. coli anchor 100%) the **out-of-fold** value is substituted, so the column
  is not part in-sample and part honest. `evidence/screens_transfer_audit.tsv` prices every
  column — training strain, own-organism grouped AUROC/AUPR, overlap — and must be read before any
  `_prob` is treated as evidence.

  **Load through `src/essentiality.py`**: `list_training_sets`, `load_training_set`, `registry`,
  `load_deg`, `load_ogee`, `load_screens`, plus the existing geptop and merge loaders.

- **`essentiality/screens.py`** + **`essentiality/predict.py`** + **`essentiality/summary.py`** —
  **published screens as independent endpoints, with a control that says the labels are not
  inverted.**

  `geptop_ess` is a prediction at 100% coverage and `deg_ess` was **0% on K. pneumoniae** — DEG
  indexes no *Klebsiella*, and neither does OGEE. The anchor's column was entirely predicted.
  These scripts fix that from screens v1 had already fetched and then lost to a gene-symbol join.

  **TEN training sets, one per SOURCE**, conditions aggregated within a source (the project
  owner's rule — never one endpoint per condition or per DEG dataset). Column names are the output
  matrix's names: lowercase, `_`-separated, self-explanatory, no `_ess` suffix.

  | column | organism | assay | n | pos | base | join | ribo |
  |---|---|---|---|---|---|---|---|
  | `essential_kpneumoniae_ecl8_tradis` | **Kp ECL8** | TraDIS/DESeq | 4,930 | 523 | 0.106 | 97.7% | 0.942 |
  | `essential_kpneumoniae_rh201207_tradis` | **Kp RH201207** | TraDIS | 4,981 | 471 | 0.095 | 92.9% | 0.938 |
  | `essential_kpneumoniae_atcc43816_tradis` | **Kp ATCC 43816** | TraDIS | 4,809 | 363 | 0.075 | 92.6% | 0.906 |
  | `essential_ecoli_k12_knockout` | Ec K-12 | **arrayed knockout** | 4,190 | 286 | 0.068 | 97.1% | 0.774 |
  | `essential_ecoli_mg1655_footprinting` | Ec MG1655 | footprinting | 4,253 | 604 | 0.142 | 98.5% | 0.538 |
  | `essential_ecoli_bw25113_tradis_goodall` | Ec BW25113 | TraDIS | 4,056 | 354 | 0.087 | 97.5% | 0.961 |
  | `essential_ecoli_bw25113_tnseq_choe` | Ec BW25113 | Tn-seq, LB only | 4,272 | 440 | 0.103 | 95.0% | 0.906 |
  | `essential_ecoli_st131_tradis` | **Ec ST131 EC958** | TraDIS | 4,981 | 300 | 0.060 | 100% | 0.741 |
  | `essential_ecoli_o157h7_tnseq` | **Ec O157:H7** | Tn-seq | 5,433 | 1,055 | 0.194 | 100% | 1.000 |
  | `core_essential_gammaproteobacteria` | clade | conservation | 4,256 | 205 | 0.048 | 99.0% | 0.906 |

  **Three filters define the set, each on instruction**: *both classes required* (a positives-only
  gene list is not a training set — drops Ramage 2017 and Paczosa 2020); *no condition-dependent
  data for now* (drops Choe's M9 arm, Rome 2026, Short 2020's serum screens, Bruchmann's
  `2hpi`/`6hpi`, every in-vivo mouse screen); *no duplicates* (DEG1019 is Keio again). Everything
  excluded stays parseable and is listed in `screens.py`'s `HELD_BACK` — **decisions, not gaps.**

  **Nothing is merged.** Base rates span **0.048–0.194**, and that spread is method and strain, not
  noise. The sharpest case is one organism and one dataset under two analyses: Goodall's own
  BW25113 calls give 354 positives and a DESeq reanalysis of the same reads gives 258.

  **`summary.py` IS THE CORRECTNESS CHECK, and it has already earned its place.** `num_positives`
  looks identical whether a label set is right or inverted, and an inverted column trains a
  confident, well-formed, exactly wrong model. Every column is scored against the **ribosome**
  (must be essential) and the **textbook dispensables** (`lacZ`, `araB`, `fadB`, flagellar,
  fimbrial — must not be). It writes `output/results/essentiality/screen_summary.tsv`.

  It caught a real one: **the compendium's `BW25113.out.DESeq.tsv` did not converge**, and the
  eleventh column was retired because of it. `padj` is NaN for **2,310 of 4,256 genes (54.3%)**
  against 7.0% for ECL8 through the same pipeline, and of the 189 genes with **zero insertion
  sites**, **108 are called `Unchanged`** (ECL8: 91 of 91 correctly `Reduced`). DESeq cannot
  compute a statistic for a gene with no insertions, so it falls through to non-essential — the
  label is inverted for exactly the genes that matter most. Ribosome recall 0.113 against
  Goodall's 0.961 on the same data. Nothing was lost; Goodall's own calls ship.

  **CALIBRATE A BIOLOGICAL CONTROL ON A MEASUREMENT, NOT ON 1.00 — this was got wrong once.** A
  first pass set the bar at 0.50 separation and failed Gerdes; it would have failed the gold
  standard too. **Keio, an arrayed knockout collection, reaches only 0.774**, and its misses
  (`rplA`, `rplI`, `rplK`, `rplY`, `rpmE/F/G/I`, `rpsF/O/T/U`) are genuinely dispensable in
  *E. coli*. So ~0.8 is the ceiling; the bars are set to catch breakage, not noise (recall ≥ 0.45,
  separation ≥ 0.40).

  **Features: ProtT5, measured not assumed.** On the Keio endpoint, 5 seeds, byte-identical folds,
  paired: ProtT5 − ESM-C **+0.0603 PR** [5/5 seeds]; ProteomeLM − ESM-C +0.0273 [5/5]; **ProtT5 −
  ProteomeLM +0.0330 PR [5/5] but −0.0002 AUROC [2/5, a tie]**. AUROC alone would have called the
  winner a coin flip — the same failure stage 04 records.

  **Every screen strain now has its own ProtT5 matrix**, in
  `data/processed/embeddings/scratch/strains/` (`prott5_<label>.npz`) — ECL8, KPPR1, KPNIH1,
  RH201207, DEG1048, DEG1056, via `prott5.py --strain`. A screen strain's embeddings are training
  features, not a deliverable matrix: no canonical row order, keyed on whatever identifier the
  strain's own FASTA uses.

  **The two DEG-derived columns key on `fasta_id`** — verbatim the header of
  `data/source/ncbi/deg_proteomes/<id>.faa` — so label and feature vector share one namespace and
  there is no join at all. Measured against the E. coli anchor by exact sequence, EC958 shares
  13.7% and O157:H7 23.0%: genuinely different proteomes. DEG1018 (Gerdes) is the opposite at
  **99.4%**, so it maps onto the anchor and needs no embedding of its own.

  **Grouped CV, and the gap is published.** Anchors group on stage-05 orthogroups; screen strains
  have no OrthoFinder run, so `essentiality/paralog_clusters.py` builds MMseqs2 clusters at 30% id
  / 80% coverage (stage 04's threshold), for all five strains. Measured on ECL8: 997 of 5,178
  proteins collapse into 520 families. **`ecl8` ungrouped 0.8743 → grouped 0.8705, leakage gap
  +0.0038 AUROC / +0.0096 PR.**

  **MMSEQS2 REWRITES FASTA HEADERS CONTAINING `|`.** It parses them as NCBI db-style fields, so
  DEG's `>lcl|HG941718.1_prot_CDN80371.1_1` comes back out of `_cluster.tsv` with the `lcl|`
  silently gone — and that id is exactly the key the screen table and the embedding share. The
  clustering then mapped to **0 of 4,981** proteins, which is indistinguishable from a proteome
  with no paralogs. `paralog_clusters.py` now writes index surrogates (`>0`, `>1`, …) and
  substitutes the real ids back, so no header convention can break it.

  **THE HEADLINE: an E. coli-trained model reaches K. pneumoniae, with no orthology anywhere.**
  `--score-on` fits on one organism and scores on another's OWN labels and OWN embeddings; no
  essentiality value crosses a species boundary. Measured under TabPFN
  (`evidence/crossspecies.tsv`):

      Ec -> Kp   mean AUROC 0.8597   best goodall -> ecl8   0.8934, PR 0.7359
      Kp -> Ec   mean AUROC 0.9338   best ecl8 -> keio      0.9732, PR 0.7431

  Geptop's measured transfer, for comparison, was 0.59–0.81 on two unrelated species and never on
  our anchor. **Transfer is better when the ASSAY matches**: Goodall (TraDIS) reaches Kp at 0.89
  while Keio (knockout) reaches 0.82–0.84, consistently in both directions.

  **Quote BOTH the grouped CV and the cross-species number, never one alone.** They answer
  different questions and the cross-species one is higher: grouped CV holds out whole paralog
  FAMILIES ("can we predict a family never seen?"), while cross-species trains on every family and
  tests on an organism where most genes have a homolog ("can we recognise a known family in a new
  organism?"). The second is not label leakage — nothing is copied — but gene-family overlap is a
  real information channel. The first is the honest generalisation measure; the second is the
  operationally relevant one for Kp.

  **The estimator is the FOREST by default on this axis** — see the exception under *Supervised
  ML* above. `--estimator tabpfn` switches through a seam narrow enough that nothing else about
  the evaluation changes. Hosted TabPFN bills ~10,000 credits per CALL, flat, and `fit` is free —
  so a whole cross-species evaluation is one call while a 5-seed 5-fold CV is 25. `--dry-run`
  prices any TabPFN run through `src.tabpfn.would_cost` before spending.

  **Per-dataset provenance lives beside the data**: a `SOURCE.md` in each of the 28 dataset
  directories under `data/raw/{ecoli,kpneumoniae,other}/essentiality/` and `data/source/{deg,ogee}/`,
  recording what the dataset is, whether v2 uses it, and — for the majority that are NOT used —
  exactly why. `docs/essentiality_screens.md` is the cross-cutting catalogue.

  CLI (`predict.py`): `--endpoint` · `--features {esmc,prott5,proteomelm}` ·
  `--estimator {forest,tabpfn}` · `--score-on` · `--seeds` · `--folds` ·
  `--schemes {grouped,plain}` · `--dry-run` · `-q`.
  CLI (`summary.py`): `--features`. CLI (`screens.py`): `--screen` · `--list` · `--dry-run` · `-q`.

- **`essentiality/fetch_screens.py`** — **downloads published screens, and writes down what a human
  must fetch by hand.** A screen we cannot download is a task, not a dead end; the real output is
  `evidence/screen_fetch_status.tsv` plus a `PLACEHOLDER.md` naming what to click.

  Ladder: direct URL → `curl -L` → Europe PMC `supplementaryFiles` → record manual instructions.
  **Publishers and repositories want OPPOSITE headers** — publishers 403 a bare client, while
  **figshare returns HTTP 202 forever to a full Chrome User-Agent and 200 to a short one**. Both are
  tried. Full route table and the recovered datasets: `docs/essentiality_screens.md`.

  **Three payloads that returned HTTP 200 and were not data**, all caught by content inspection:
  a Nature article page (407 KB of HTML), a PMC `/bin/` URL (a reCAPTCHA page), a figshare 202
  (0 bytes). The house rule earns its keep; assert on content, never on status.

- **`studiedness/fetch.py`** + **`gene2pubmed.py`** + **`unknome.py`** + **`transfer.py`** +
  **`merge.py`** — **how much is already known about this protein?** Wanted in both directions: an
  uncharacterised target is a risk, but it is also the novelty the collaboration is looking for
  (`src.studiedness.novelty()` reads it the other way). Deliverable
  `studiedness_<species>.tsv`: **4 columns** — `uniprot_ac · n_papers_own · n_papers_family ·
  evidence`, complete and canonical for the three bacteria.

  **THE NUMBER IS A PAPER COUNT — curated SwissProt references, nothing scaled or blended.**
  `n_papers_family` is the number of curated PubMed references on the best-studied prokaryotic
  SwissProt homolog. A 5 is five papers; a 0 is zero papers. **A 0-1 composite shipped first and
  was REMOVED on 2026-09-22 — do not reintroduce it.** It read
  `0.6*min(1, log1p(n)/log1p(204)) + 0.4*(annotation_score-1)/4`, and three measurements killed
  it: (1) the same value meant **13 papers at annotation 3, 83 at annotation 1, 1 at annotation
  5**, and a protein with zero papers scored 0.4; (2) the halves **double-counted** — UniProt's
  annotation score correlates with the paper count at **r 0.64-0.70**; (3) the weights did the
  opposite of what the code claimed — literature was weighted 0.6 "because it is the quantity the
  axis is named for" but the annotation component's spread was nearly **double** (sd 0.33-0.36 vs
  0.17-0.20). Dropping it moved the ranking by spearman 0.86-0.91.
  `src.studiedness.scaled()` derives a 0-1 column on the fly for combining with the other axes;
  it is deliberately **not stored**, so there is one source of truth on disk.

  **There is NO ceiling on the count, and an earlier claim in this file that SwissProt "saturates
  at 58" was wrong.** Verified against UniProt: the TSV is not truncated (human TP53 returns 225
  PubMed ids), counts reach **225 overall / 119 among prokaryotic entries**, and the 58 is just
  *E. coli* GroEL (`P0A6F5`) happening to be the most-curated donor picked in all three species.
  The range is small because bacterial proteins are: **37-48 distinct values, median donor 4-6**.
  Ties sit in the poorly-studied bulk, not at the top (largest non-zero tie 12.4-18.9% of scored
  proteins, guarded at 25%), so a shortlist reading the top gets a near-unique ranking.

  **gene2pubmed is a MEASURED ALTERNATIVE, not the count.** It is larger for 87-94% of donors
  (median 2.8x) and would lift Sa coverage 12.4% -> 18.7%, but `max(curated, gene-linked)`
  switches definition per row — Kp `rpoB` took NCBI's 350 while ~11% of donors took SwissProt's.
  One consistent definition beat the bigger number. Both components ship in `evidence/`.

  **The anchors' own literature is a dead column, and that is the whole shape of the axis.**
  Measured on UniProt 2026_03: **5,710 of 5,728 Kp proteins carry exactly one PubMed id** (the
  genome paper — 7 distinct values over the proteome) and 2,532 of 2,889 Sa proteins carry none,
  against E. coli's median 5 / 48 distinct. This confirms v1 (`legacy/HISTORY.md:60`) and sharpens
  it. `n_papers_own` still ships **because it IS the measurement of darkness** — keeping it
  beside `n_papers_family` is what makes "dark in *Klebsiella*, famous in *E. coli*" readable
  off one row. **Do not rank Kp or Sa on `_own`; rank on `_family`.**

  **The number is transferred by DIAMOND from the most-cited SwissProt homolog** (575,748 reviewed
  entries, ~14,000 species) — the house "map by sequence, not by accession" rule. The free
  four-species ortholog table was measured against it and **reaches an Ec-or-human ortholog for
  only Kp 68.5% / Sa 46.6%**: *S. aureus* is Gram-positive, so under half of it has an E. coli
  ortholog at all. Both routes ship in `evidence/route_comparison.tsv`; nothing is merged. **Read
  that table honestly — the free route reaches slightly MORE Kp (68.5% vs 66.0%)**, because it
  applies no identity or coverage floor, but it wins nothing elsewhere (Ec 22.3% vs 99.3%) and
  yields only the existence of an ortholog, not a literature count.

  **`--max-target-seqs` was the worst bug in the axis and the lesson generalises.** DIAMOND
  returns the top k by **bitscore** = the CLOSEST relatives, but this stage wants the **best-cited**
  homolog, which for a conserved protein is a distant model-organism entry. At k=50 Sa `groEL`'s
  50 candidates were all Staphylococci and Bacilli, its donor was *B. subtilis* GroEL with **9
  papers**, and E. coli GroEL (246) never appeared — **systematically understating studiedness
  precisely for the most conserved, most studied families**. Fixed with k=500 plus a **deep second
  pass at k=5,000 over only the 1,403 queries (12.6%) still capped**; after it, 0 remain capped.
  Held-out control at the 40% floor: **0.4666 → 0.5338 → 0.5436**. `n_candidates` ships per
  protein and the run prints the capped fraction.

  **Five evidence tiers, and TWO of them score 0 meaning different things.** `swissprot_direct`
  (≥95%) ⊃ `swissprot_close` (≥60%) ⊃ `swissprot_homolog` (≥40%), then **`below_floor`** (a hit
  exists but under the floor) and **`no_hit`** (nothing in SwissProt at all). Both score 0.0 and
  the axis refuses to invent a number for either, but `no_hit` is the strongest novelty claim it
  makes while `below_floor` merely has a too-distant relative — without the split a third of Kp is
  one undifferentiated tie. **A zero is an answer, not a gap; never impute it.** The hard ceiling
  is DIAMOND's own hit rate: any SwissProt hit exists for only **81.6% of Kp, 73.0% of Sa**.

  **The floor is 40% and the trade-off is monotonic, not a calibration.** Decoys settle nothing —
  composition-preserving shuffles of our own 13,020 sequences match at **0.0% at every floor down
  to 20%**, so DIAMOND's e-value already handles spurious homology. The held-out control decides
  it: 25% → Kp 77.6% / rho 0.4537; **40% → Kp 66.0% / rho 0.5436**; 60% → Kp 56.6% / rho 0.5490.
  40% is CLAUDE.md's standing transfer threshold and sits at the knee. **Preferring the closest
  band over the most-cited donor was tried and is WORSE** (0.4426 vs 0.4666).

  **THE CONTROL HAS A CEILING WELL BELOW 1 — DO NOT MAXIMISE IT.** It correlates the transferred
  family score against E. coli's OWN literature, which are deliberately different quantities: a
  protein with 3 papers whose human homolog has 300 *should* score low on `_own` and high on
  `_family`. And rho rises with the floor (0.5922 at 95%) nearly circularly, because the exclusion
  removes *Escherichia* only and the survivors at high identity are Salmonella/Shigella
  near-duplicates. Maximising it drives the floor to 95%, where `_family` collapses onto `_own`
  and the axis stops doing anything. Measured: **spearman 0.5411 / pearson 0.5351 over 2,945
  proteins**, floor 0.45, exits non-zero below.

  **`family − own` is the axis working, and the contrast is the clearest number here: +3 papers
  (median) on Kp** (dark anchor, known family) **against +0 on E. coli** (its own literature
  already is its family's). On Kp the median protein has **1 paper of its own and 4 on its
  family**; Kp `groEL` reads **1 and 58**. Coverage with a donor: Kp 65.8% · Ec 99.3% · Sa 55.8%.

  **DONOR SCOPE: not eukaryotic, and "restrict to Bacteria" is the obvious rule that is WRONG.**
  `--donor-scope prokaryotic` (default) admits any donor whose UniProt `lineage` lacks
  `Eukaryota (domain)` — Bacteria, Archaea **and phages**. The restriction exists because under
  `any`, five of Kp's ten highest-scoring proteins took human donors (HSPD1 934 papers, HADHA,
  CTPS1, LONP1, AFG3L2), making the claim "well studied because its human mitochondrial homolog
  is" — the ligands axis's recorded mistake in another guise (an unrestricted non-human bucket
  gave 424 apparent potent Kp proteins against a true 175). **But strict Bacteria strands
  prophage proteins: of the 93 proteins it left with no donor, 70 (75%) lost theirs to a VIRUS**
  — *Escherichia* phage lambda and P1 — and a Kp prophage protein whose best relative is a lambda
  protein is not novel. Admitting phages recovers 62 of Kp's 77 and 12 of Sa's 16.
  **All three scopes are computed every run** into `evidence/donor_scope_comparison.tsv`, with
  `n_papers_family_{prokaryotic,bacteria,any}` all in the transfer table — switching the
  deliverable is a column swap, not a re-run.
  **The honest headline is that scope barely matters numerically**: the three counts correlate at
  **rho 0.978–0.991** and only 215 of 5,728 Kp proteins change at all. It is kept because it makes
  the number mean the right thing, not because it moves it.
  The control cannot arbitrate — **prokaryotic 0.328 sp / 0.350 pe · bacteria 0.322 / 0.344 ·
  any 0.339 / 0.348**: spearman spans 0.006 and the two metrics *disagree* on the winner, which
  is the signal that it is noise. `prokaryotic` takes it on pearson.
  **Tiers are scope-consistent**: `below_floor` means an IN-SCOPE hit below the floor, so a
  protein whose only curated relative is eukaryotic reads `no_hit`, not a distant prokaryotic one.

  **Unknome is EVIDENCE, kept and not promoted.** Joined via the `panther` xref already on disk to
  the **cluster** table — and the trap is severe: the per-protein download reads
  **`knownness = 0.000` for all 1,882 S. aureus entries** including `clpP`, `rpoB`, `ftsZ`, while
  Sa `clpP`'s cluster `UKP00027` (shared with Ec `clpP`) reads 10.6. Joining it would have
  silently zeroed a proteome. Coverage is capped by the xref (Kp 71.2% · Ec 78.2% · Sa 65.1%) and
  Unknome has **no Klebsiella**. It counts **GO terms, not papers**, so it is a real second opinion
  (spearman vs `_family` 0.35–0.52) — but it is silent exactly where the literature route is also
  silent, so it cannot become the axis. **REJECTED and not to be re-added: training ESM-C on
  Unknome knownness** — knownness is a cluster property so the effective N is 15,589 not 1.9M, it
  approximates a lookup DIAMOND does exactly and auditably, mean-reversion would hide the novel
  proteins the axis exists to surface, and a fourth ESM-C column correlates with the other axes for
  reasons unrelated to studiedness. Reasons in full: `docs/studiedness.md` §5.

  **A finding for the collaboration: the consortium's own panel is NOT novel** — `src/interest.py`
  sits at the **80th percentile (median) on Kp, 74th on Ec, 69th on Sa**. Read
  `interest.coverage()` first; the panel matches by gene symbol.

  **Load through `src/studiedness.py`** — `load`, `load_all`, `novelty`, `load_own`,
  `load_transfer`, `load_unknome`, `load_gene2pubmed`, `load_route_comparison`,
  `load_floor_sensitivity`, `load_donor_scope_comparison`, `control`, `definition`, `scaled`,
  `manifest`.

  Run in order: `fetch.py` → `gene2pubmed.py` → `unknome.py` → `transfer.py` → `merge.py`. The
  `gradi` env throughout; DIAMOND borrowed from `gradi-ortho` via `GRADI_DIAMOND_BIN`. ~20 min
  cold. CLI: `--species` · `--refresh` · `--dry-run` · `-q`, plus `transfer.py`'s `--limit`
  (smoke, writes only `scratch/`), `--donor-scope {prokaryotic,bacteria,any}`, `--threads`,
  `--max-targets`, `--deep-targets`, `--sensitivity`, `--no-control`.
  Figures: `scripts/plots/studiedness.py` (stylia, one at a time) → `studiedness.png`,
  `control.png`, `interest_panel.png`.
  **Deletable after a run** and never to be uploaded to eosvc: `uniprot_sprot.fasta.gz` (89.5 MB),
  `gene2pubmed.gz` (273.8 MB), and `data/processed/studiedness/scratch/` (805 MB of DIAMOND hits).
  Details, traps and the run log: `docs/studiedness.md`.

Registry tiers: **A** the 4 anchors · **B** same-species name-donor pools · **C** a 26-species
comparator panel for orthology (v1's curated panel, with every species now pinned to an explicit
proteome id) · **D** evidence-bridge strains from **NCBI**, because **UniProt serves zero entries for
any proteome it does not flag "Reference proteome"** — verified for ECL8, KPNIH1, KPPR1, NJST258,
BW25113, COL, Mu50. C and D are fetched only on request.

Two tier-D traps worth knowing: **S. aureus COL must come from GenBank, not RefSeq** (`GCF_` retains
60 original `SACOL####` tags, `GCA_` retains 2,711), and **ECL8's GFF carries
`old_locus_tag=BN373_…,KPNEcl8_…`**, which should make v1's `gradi-prokka` re-annotation unnecessary.

## Legacy

**`legacy/` holds the complete v1 pipeline, frozen. Do not extend it.** Start at
**`legacy/HISTORY.md`** — the retrospective: what was built and what was not, the methodological
decisions and why, ~70 documented traps, which data sources were obtained and how, and which
artifacts look like data but are not.

Three things there must not be trusted (`legacy/HISTORY.md` §7): the **Kp legacy degradability TSV**
(four verified defects, including an inverted N-end rule supplying 680 of 698 `medium` calls), the
**E. coli degradability values in the webapp** (a deterministic MD5 mock), and
**`data/raw/legacy/clp_substrates/`** (45- and 35-row hand-curated substitutes, not the papers'
tables).

The archived scripts still run. They resolve paths via `parents[1]`, so `legacy/data` and
`legacy/output` are **symlinks** to the real trees one level up — that is what keeps them working
without editing 78 files. Without those symlinks a script exits 0 having written nothing. Run from
inside `legacy/`:

```bash
cd legacy && python scripts/09h_localization_plots.py --organism ecoli
```

The v1 webapp is still deployed from `legacy/app/` by `.github/workflows/pages.yml`.

## Setup

Conda environments (see `install.sh`):

- **`gradi`** (Python 3.11) — the main env; `bash install.sh` runs `pip install -r requirements.txt`.
  **Stage 00 needs nothing beyond `requests`, `pandas`, `tenacity` and `pyarrow`.** Stage 06 adds
  **`rdkit`** (Bemis–Murcko scaffolds) — install it from **PyPI**, whose `macosx_11_0_arm64` wheel
  keeps the env native; a conda install that flips `gradi` to osx-64 would take ESM-C down with it,
  the same trap as `blast` and `eggnog-mapper`.
- **`gradi-ortho`** (osx-64 bioconda, Rosetta on Apple Silicon) — OrthoFinder + DIAMOND, which have no
  arm64 build. Its DIAMOND binary was borrowed by many v1 scripts via `GRADI_DIAMOND_BIN`.
- **`gradi-prokka`** (osx-64, Rosetta) — BLAST+ suite. Stage 02 borrows its **`rpsblast`** via
  `PATH` (`GRADI_RPSBLAST_BIN` overrides). Installing `blast` into `gradi` would force the whole env
  to osx-64 — there is no arm64 build.
- **`gradi-pockets`** (osx-64, Rosetta) — `fpocket` + `openjdk=17` for P2Rank.
- **`gradi-pymol`** — `pymol-open-source`, for ray-traced structure cartoons.
- **`gradi-tabpfn`** (Python 3.11) — **TabPFN-3.5, the project's default supervised learner**
  (see the section at the top of this file). `tabpfn==9.0.0` + `tabpfn-client==0.6.0`. **This split
  is mandatory**: it pulls torch 2.14 + mlx against `gradi`'s 2.12, and installing it there would
  break stage 01's ESM-C — the same trap as `blast`, `eggnog-mapper` and `rdkit`. Stage 04 runs in
  `gradi` and reaches across a process boundary to `scripts/workers/tabpfn_cv.py`; never activate
  this env to run a stage. Needs `TABPFN_TOKEN` in the environment.
- **`gradi-lazyqsar`** (Python 3.11) — lazy-qsar 3.4.4, for the stage-04 head comparison only.
  **Rejected as an estimator** (see the top of this file); kept so the comparison is reproducible.
  Pins `numpy==2.1.3` / `scikit-learn==1.6.1`, so the same mandatory split applies.
- **`gradi-loc`** (Python 3.11) — DeepLocPro + TMbed, for **stage 03**. **This split is mandatory,
  not cosmetic**: DeepLocPro needs `fair-esm`, which claims the same top-level `esm` package as the
  EvolutionaryScale `esm` used for ESM-C embeddings. Installing it into `gradi` silently breaks
  stage 01. Stage 03 runs in `gradi` and reaches across a process boundary (`GRADI_LOC_BIN`
  overrides the location) — never activate this env to run the stage. Pins `setuptools<81`
  (DeepLocPro imports `pkg_resources`) and `transformers==4.44.2` (TMbed's ProtT5 `T5Tokenizer`
  dies on transformers 5.x with a spurious tiktoken error).

Do **not** use the machine's default `python3` — it resolves to an unrelated `ersilia` env. Use
`~/miniconda3/envs/gradi/bin/python`.

No build, lint or test commands are configured yet. Document them here when added.
