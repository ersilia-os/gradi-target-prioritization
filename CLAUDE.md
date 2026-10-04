# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

**This file is the MAP, not the record.** It holds what someone needs to work here *today*: what
each stage does, what it writes, which conda env it needs, and the traps that are silent. The
measured tables, the rejected alternatives and the run logs live in `docs/<task>.md`, one per axis;
everything ever tried in v1 lives in `legacy/HISTORY.md`. **When this file and a doc disagree, the
doc is the record** — fix this file.

## Maintaining this file

**A pipeline step is not done until this file describes it.** In the same commit as the code, add or
update: the stage entry (what it does, what it writes, which conda env), any new silent trap, any
new environment requirement.

**Every path written here is a claim, and claims get checked.**

```bash
P=~/miniconda3/envs/gradi/bin/python

$P tools/check_claude_md.py   # EVERY path this file names must exist. Non-zero exit if not.

# every script still imports and parses its CLI -- run from ELSEWHERE, which is what proves
# `parents[2]` and sys.path are right. Run plot scripts ONE AT A TIME (see the stylia trap).
for f in scripts/*/*.py scripts/*/*/*.py; do (cd /tmp && $P "$OLDPWD/$f" --help >/dev/null) || echo "BROKEN $f"; done

$P -c "from src import proteomes,embeddings,function,localization,degradability,orthology,\
ligandability,essentiality,projections,proteomelm,matrices,tabpfn,pockets; print('ok')"

$P -m src.matrices            # every matrix complete AND in canonical order
```

Two wrong paths already shipped and both read fine: a blanket rename put `smoke_*` files under
`evidence/` when smoke output is `scratch/` by definition, and a path rewrite produced
`scripts/localization/localization/workers/`. `check_claude_md.py` knows three legitimate
exceptions — read its header before adding one: dumps documented as deleted after use,
`<placeholder>`/`*` patterns, and `NEGATIVE_EXAMPLES`, paths quoted **because they are wrong**. Add
a cautionary example there rather than weakening the check.

**Four house rules for edits here**, all violated at least once:

1. **Name scripts by their new path** — `function/cog.py`, never `02_function_cog.py`; `docs/function.md`, not `docs/02_function.md`.
2. **Name the tier, not `accessory/`** — that directory no longer exists. A file is in the task root, in `evidence/`, or in `scratch/`.
3. **Say which conda env a stage needs**, and whether it crosses a process boundary. The env splits are mandatory; each exists because installing the dependency into `gradi` broke something specific.
4. **Quote measured numbers with what produced them.** "0.8738" alone is not maintainable; "ADEP4 0.8738 under TabPFN-3.5, 5 seeds, cluster-grouped CV" survives an estimator change.

**A standing instruction from the project owner becomes a `##` section**, not a line inside a stage
entry. Stage entries describe one stage; standing rules bind every stage, including ones not written
yet, and burying one in a stage entry is how it gets missed.

## House style

**Scripts print generously**: a banner naming the script, its inputs and outputs; per-step progress
with running counts; a closing coverage table with percentages. Nothing important happens silently.
Offer `-q/--quiet`.

- **Verify counts against an authoritative source before trusting a payload.** An HTTP 200 is not evidence of data — v1 was repeatedly bitten by short or empty responses (a 61-byte 200, an empty 202, a 200 from the wrong article). Assert on content, and fail loudly.
- **Never resolve ambiguity silently.** Emit the chosen value, its source and the alternatives, plus an audit row. If a rule decided between candidates, record which rule.

**Figures use stylia** (`ersilia-os/stylia`), which calls `shutil.rmtree(matplotlib.get_cachedir())`
**at import time** — so two plot scripts running concurrently delete each other's font cache and die
with a `FileNotFoundError` naming `~/.matplotlib`. **Run plot scripts one at a time.** The signature:
a *different* subset fails each run, and every script passes alone.

## Supervised ML: always TabPFN-3.5

**Whenever this project needs supervised machine learning — classification OR regression — use
TabPFN-3.5.** Not RandomForest, not a hand-tuned sklearn model, not an AutoML wrapper. Standing
instruction from the project owner; applies to every stage, existing and future.

**ONE STANDING EXCEPTION, granted 2026-09-21: the ESSENTIALITY ENDPOINTS run a RandomForest until
the project owner is convinced by the datasets** — *"for now, use the random forest option (tabpfn
only when i am convinced)."* **Scoped to that axis only**, and it lapses when the owner says so. It
is also right on the merits: the question there is *which of ten datasets is learnable at all*, and
a 10-endpoint × 5-seed × 5-fold sweep is **250 hosted calls ≈ 2.5M credits, 12.5% of the monthly
quota**. **Keep the comparison possible** — `predict.py` routes both estimators through ONE
`predict_fold`/`oof_predict` seam on identical folds, grouping and metrics, the forest is
`src.degradability.RF_PARAMS` imported rather than re-specified, and `estimator` is written into
every results row.

- **Call it through `src/tabpfn.py`** — dispatch, cache, retry and credit guard live in `src/` so any axis can reach them; nothing shells out to the worker directly.
- Package `tabpfn` 9.0.0 (local) or `tabpfn-client` 0.6.0 (hosted). Env **`gradi-tabpfn`** — **never install it into `gradi`**, it pulls torch 2.14 against gradi's 2.12 and would break stage 01's ESM-C. Cross a process boundary: `scripts/workers/tabpfn_cv.py`, `GRADI_TABPFN_BIN` overrides.
- **`GRADI_TABPFN_CACHE_ONLY=1` refuses to spend**, and is the ONLY correct way to verify the cache after a refactor: if the stage still reproduces its numbers, the cache survived. "Verifying" with a novel input always misses, which looks identical to a broken cache. Both mistakes have been made; the key change cost **1,540,000 credits**. Cache at `data/processed/tabpfn/cache`, shared across axes, filenames ARE the keys.
- Under that flag, degradability reproduces **ADEP4 0.8738 / PR 0.6103, ONC212 0.7671 / PR 0.5803**.
- **No PCA needed at our sizes, but know the real limits**: the installed `tabpfn` 9.0.0 declares `MAX_NUMBER_OF_SAMPLES` 10,000/50,000 and `MAX_NUMBER_OF_FEATURES` 500/2,000 depending on inference config. So 1,152-d embeddings go in whole, but ~170k rows do NOT fit one call — which is why per-endpoint modelling is the right shape. **An earlier version of this file claimed "20,000 features and 1,000,000 rows"; that was wrong** — measure against the installed package before sizing a corpus.

**The evidence, so nobody re-litigates it.** On a 3-estimator × 3-embedding × 2-activator grid with
byte-identical cluster-grouped folds (`output/results/degradability/proteomelm_vs_esmc.tsv`), TabPFN
beat the hand-set RandomForest on **PR-AUC in 6 of 6 arms, each winning all 5 seeds**; **lazy-qsar
3.4.4 was tested in the same grid and rejected**. The regressor wins too: against the ESM-C ridge on
continuous log2FC, same folds, paired — **4 of 4 comparisons, each winning all 5 seeds**.

**Three rules that come with it:**

1. **Report PR-AUC alongside AUROC, never AUROC alone.** The gain is ~4× larger on PR-AUC — concentrated at the top of the ranking, which is what a shortlist consumes. AUROC alone called 4 of those 6 arms "no difference".
2. **Compare arms with a PAIRED test** — per-seed differences ± standard error on identical folds, plus `[k/5 seeds]`. Two independent SDs is badly under-powered when partition difficulty is shared.
3. **Estimator and features interact — do not pick embeddings with a convenient model.** ProtT5 is the *worst* features under a forest and the *best* under TabPFN.

**Two open caveats**, neither resolved: the weights are **non-commercial licensed** (fine for methods
work; must be settled before GraDi's outputs ship), and local weights need the licence **accepted**
at ux.priorlabs.ai — a valid `TABPFN_TOKEN` alone is not enough, so only hosted works. Hosted bills
~10,000 credits per *call* (flat, independent of data size; `fit` is free) against a 20M monthly
quota, and uploads features and labels to a third party.

**An open lead, deliberately not a finding:** TabPFN *regression* ranks slightly better than TabPFN
*classification* on the same binary labels (0.8812 vs 0.8738 on ADEP4), which would mean the
continuous log2FC carries signal the `<= -1` cutoff throws away. Both gaps sit inside the combined
seed noise; needs a paired test between the two framings before anyone acts on it.

## Every axis ends in COMPLETE matrices

**One row per protein, always.** A protein with no annotation is an **all-zero row**, never a
missing row. Standing instruction from the project owner.

**An `evidence` column is the default, NOT an invariant** — narrowed 2026-10-03, when the owner
had studiedness drop it. Where an axis omits it, the tiers must still ship in that axis's
`evidence/` tree and the loader docstring must say what a 0 can no longer distinguish:
`studiedness_<sp>.tsv` has no `evidence`, so a `n_papers_uniprot_prokaryotic == 0` is `no_hit` or
`below_floor` and `load_transfer()` is the only way to tell — 1,961 Kp proteins.

**The invariant is one row per protein, complete, in canonical order.** The SHAPE a vocabulary axis
ships in is a per-axis choice, not a project rule — narrowed 2026-10-03, when function moved to
`;`-packed term columns. What does not change: the row set, the row order, and that a protein with
nothing known is a present row, not a missing one.

**When an axis ships packed, the term matrix stays beside it in `evidence/`**, because a packed
list **cannot express a structural zero**: it cannot tell a term that is merely unannotated from one
the organism cannot reach. 8 GO-slim terms are eukaryote/plant concepts, *S. aureus* has 19
unreachable because it is Gram-positive, COG `Y` is nuclear structure — all of which are **kept
columns** in the matrix and simply absent in the packed column. The two forms must be **provably
interchangeable**, asserted in both directions, not merely both present.

**AND IN THE SAME ORDER.** Every matrix has the same rows in the same order: the order
`data/processed/proteomes/proteome_<species>.tsv` is written in. Not cosmetic — once it holds, any
two axes stack with `np.hstack` or `pd.concat(axis=1)` and no join at all. The point is the failure
mode: forgetting to write one does not raise, it silently misaligns a column. And an embedding
matrix has no accession column to join on in the first place — its rows are positional, so alignment
there is *only* ever by order.

Use `src/matrices.py`: `reindex()` / `reindex_arrays()` on the write side, `assert_canonical()` on
the read side, `python -m src.matrices` to audit everything the project ships. Completeness is
checked in the same breath, because the two rules are one rule. **`reindex()` deliberately REFUSES
to fill a missing protein**: a NaN row invented there is a protein the axis never measured, dressed
as one it measured as unknown. Fill it where the rows are built, and say what the fill means.

1. **Full vocabulary as columns, not just observed terms**, so every species file has the same shape and the three stack. Terms an organism structurally cannot have stay as **structural zeros — that is information**: 8 of 97 GO-slim terms are eukaryote/plant concepts, COG `Y` is nuclear structure, and *S. aureus* has no periplasm or outer membrane because it is Gram-positive.
2. **Keep MULTI-LABEL.** A single "chosen" term discards a lot: 34–52% of GO-annotated proteins carry more than one slim term, 12.5% of COG-classified proteins more than one letter.
3. **A zero means "not annotated", not "absent"** — for 26.3% of Kp that means nothing is known. Say so in the loader docstring; downstream must never read it as a measured negative.
4. **Verify with a ROUND-TRIP, not a shape check.** Reconstructing the term lists from the matrix must reproduce the source columns exactly. Shape checks pass on wrong matrices.

Status, from `python -m src.matrices`: **60/60 canonical** (2026-10-03, after function's packed `function/function` joined its two matrices in the audit). A new representation of the same
proteins (e.g.
`proteomelm_<species>_orthodb.npz`) goes **into the audit list, not beside it** — an unaudited matrix
is exactly the silent misalignment this rule exists to catch. Localization now ships ONE complete canonical table,
`localization_<sp>.tsv`; its `localization` column is still single-label categorical and would need
one-hot over the 6-class union before it is a matrix in this sense.

## External models: never reuse an embedding you have not proven identical

**Before feeding any embedding we already hold into an EXTERNAL model, prove it is numerically
identical to what that model expects — same model TYPE and same model VERSION.** If it cannot be
proven, generate the embedding with the tool's own code path, however long that takes. Standing
instruction from the project owner.

**The bar is numerical identity, not similarity. Cosine 0.99 is a *fail*.** The check that passes
looks like `scripts/embeddings/prott5.py`, reproducing UniProt's published ProtT5 vectors at **median
cosine 1.000000, worst 0.999994**.

**Check all four, not just the first:**

1. **Model family** — ESM-C is not ESM-2 is not ESM-1b is not ProtT5.
2. **Model version and variant** — `prot_t5_xl_uniref50` (full) vs `prot_t5_xl_half_uniref50-enc` (half precision, encoder-only) are *different weights*, though both are "ProtT5-XL-U50".
3. **Everything between the model and the vector.** Each of these alone changes the numbers: **pooling method**; **which tokens are pooled** (ProteomeLM includes BOS/EOS, `embeddings/esmc.py` strips them — cosine 0.999970 median, small and still a different convention); **which LAYER** (ProteomeLM-L layer 8 of 18, and nothing about a vector's shape says which layer it came from); **normalisation**.
4. **Residue encoding and truncation** — whether U/Z/O/B map to X, whether sequences are windowed and at what length.

**Why this is a hard rule.** The failure is silent and unfalsifiable: a predictor handed a subtly
different representation from the one it was fitted on still returns plausible, well-formed numbers,
and nothing in its output says they are wrong. The saving is only compute.

**Worked example — SAFPred.** Its shipped demo embeddings are **1280-d = ESM-1b**, not ProtT5, so
the authors' own reference data cannot validate a ProtT5 substitution; and it calls the *full*
`prot_t5_xl_uniref50` while our ProtT5 is TMbed's *half, encoder-only* build. Verdict: **SAFPred
generates its own embeddings.** Ours stay valid for OUR OWN models but are not fed to someone else's
fitted tool.

## Status

Active **target-prioritization** analysis for the **GraDi** collaboration (BacPROTAC-style targeted
protein degradation in Gram-negatives; Prof. Erick Strauss, Stellenbosch). Ersilia's deliverable is a
prioritized list of proteins of interest. Ligand identification is out of scope.

**This is v2, a deliberate restart.** v1 is complete and frozen under `legacy/` — see *Legacy* below
before assuming anything about prior work.

**The consortium's own proteins of interest live in `src/interest.py`.** There was no
machine-readable list anywhere in the repo; the targets are stated in prose in two legacy documents,
and the kick-off's "targets in the periplasm" conflicts with the v5 proposal's cytosolic GyrA/GyrB,
never settled. That module is **an expansion of prose, flagged as a draft**: every judgement is in
its `EXPANSION_NOTES`, and `PANEL` is a plain dict so swapping in a curated file changes nothing
downstream. **Match is by gene symbol, so read `coverage()` before concluding a target is absent** —
Sa's 11/43 match is mostly the *structural* absence of LPS/Lpt/Bam in a Gram-positive, not a naming
gap.

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
and **MGH 78578** (historical reference; reviewed coverage ~0.5% species-wide).

***S. aureus* is new in v2.** It is the organism the activated-ClpP proteomics (Conlon 2013 ADEP4,
Jacques 2020 ONC212) is native to — v1 could only reach that data by cross-phylum DIAMOND RBH at
median 42% identity — and the only one of the three bacteria with ClpC/McsB, the machine every
published BacPROTAC actually targets. Whether it becomes a full target organism is deferred.

**Human must be fetched `reviewed:true`** — the unfiltered proteome is 147,506 TrEMBL-bloated entries.

## Identifier convention

**UniProt accession is the canonical key** in every dataset. Gene symbols, locus tags, b-numbers,
RefSeq and GeneID are provenance columns — never the primary key. `KPHS_*` is a locus tag, not an
accession.

**But join on `locus_tag`, not on gene name.** `locus_tag` is at 100% coverage on Kp while
`gene_name` is at 18.4%, and essentially every published bacterial dataset (Tn-seq, TraDIS, CRISPRi,
proteomics) keys on locus tags. v1's single most expensive recurring bug was a gene-symbol join that
silently lost ~27% of known essentials. Two corrections from the stage-00 figure: `locus_tag` is
**98.4% on *S. aureus***, not 100%, and **human has no locus tags at all** — the rule is bacterial;
use `uniprot_ac` for human.

**Map by sequence, not by accession, when reaching an external database.** HS11286 is a dark TrEMBL
proteome whose accessions rarely match ChEMBL, BindingDB or PDB-SIFTS. v1's convention, which
worked: DIAMOND blastp, **≥95% identity = "direct"**, **≥40% floor for transfer**, and restrict the
bacterial bucket to *true Bacteria* rather than "non-human".

## Identifier mapping: spend what it takes, and track the rate

**Standing instruction from the project owner.** Identifier mapping is where this project's data is
won or lost, and the losses are silent: a decimated join looks exactly like a small dataset.

| screen | v1's route | v2's route | gain |
|---|---|---|---|
| Ramage 2017 KPNIH1 | gene symbol → 212/424 | GenBank `locus_tag` → **424/424** | **2.0×** |
| BN373 / ECL8 TraDIS | never attempted | GenBank `locus_tag` → 4,930/5,048 (97.7%) | new |
| Goodall 2018 | gene symbol | compendium symbol→b-number + symbol fallback → 97.6% | +5.3 pp |

1. **Try every route, MEASURE each, then pick.** Goodall: compendium route 96.4%, our own symbol index 93.2%, **union 98.6%** — the union earned its keep. Never assume one route is enough.
2. **Prefer GenBank (`GCA_`) over RefSeq (`GCF_`) for any strain whose screen predates the current annotation.** PGAP re-annotation silently drops submitter locus tags — **263 of RefSeq's 281 ECL8 misses were tags absent from the annotation entirely**, not proteinless genes, and GenBank was a strict superset both times.
3. **NCBI GFF splits what you need across two features**: `old_locus_tag` on the **gene**, `protein_id` on the **CDS**, linked by `ID`/`Parent`. Walk both or the map comes out empty.
4. **Write the rate to an audit table every run, with a floor that exits non-zero.** See `data/processed/essentiality/evidence/screen_join_audit.tsv`.

**And know when to stop — then write the dead end down.** Eichelberger 2024's ECL8 screen could not
be joined at all (0 of 5,074 coordinates matched). Recorded in `scripts/essentiality/screens.py`
under `UNJOINABLE` so nobody re-derives it.

## Two-track persistence: Git vs eosvc

Enforced by `.gitignore`:

- **Git**: `src/`, `scripts/`, `docs/`, `assets/`, `legacy/` (code and docs), `LICENSE`, `README.md`, `install.sh`, `requirements.txt`, `access.json`.
- **eosvc (DVC + S3), not Git**: `data/` and `output/`. Use `~/miniconda3/bin/eosvc {view,upload,download,delete} --path <relative path>`.
- `tmp/` is local-only scratch, tracked by neither. The only generated artifacts in Git are `legacy/app/data/*.json`, the webapp payloads.

`eosvc upload` PUTs **one object at a time (~2 files/sec)**, so a directory of many small files is
hours, not minutes — archive such trees first. And `eosvc view --path <file>` is unreliable for a
single file; read the parent directory listing's status column instead.

**Never upload to eosvc**: the ChEMBL, eggNOG and OrthoDB dumps, SwissProt/gene2pubmed downloads, and
any `scratch/` tree. All are public and re-derivable, and `scratch/` is ~18 GB.

## Directory contract

**One folder per task, no stage numbers.** This replaced 32 flat `scripts/NN_*.py` files in
September 2026. The numbers encoded an order that was partly fiction — function, degradability and
ligands are largely independent — and they had drifted. **Execution order is deliberately left
unexpressed** until it is actually known; nothing in the layout implies one.

```
scripts/
  proteomes/      download.py                                     -> proteome_<sp>.tsv
  embeddings/     esmc.py  prott5.py  proteomelm.py  projection.py -> <model>_<sp>.npz
                  orthodb_group_check.py  workers/prott5.py
  function/       cog.py  eggnog.py  goslim.py  matrix.py  deepgo.py
                                                                  -> function_<sp>.tsv
  localization/   predict.py  merge.py  workers/deeplocpro.py     -> localization_<sp>.tsv
  orthology/      orthofinder.py  orthodb.py                      -> orthologs/neighbors/orthodb_<sp>
  degradability/  predict.py  enrichment.py  regressor.py
                  head_comparison.py  workers/lazyqsar_cv.py      -> degradability_<sp>.tsv
  essentiality/   labels.py  deg_proteomes.py  geptop.py  ogee.py
                  screens.py  predict.py  merge.py  registry.py   -> essentiality_<sp>.tsv
  ligands/        chembl.py  bindingdb.py  effort.py
                  ligands.py  transfer_calibration.py             -> ligands_<sp>.tsv
  studiedness/    fetch.py  gene2pubmed.py  unknome.py
                  transfer.py  merge.py                           -> studiedness_<sp>.tsv
  pockets/        structures.py  esmfold.py  predict.py
                  pdb_coverage.py  alphafill.py  holo.py  merge.py -> pockets_<sp>.tsv
  interactome/                              README.md only -- a real axis, no code yet
  plots/          10 scripts, ALL figures
  workers/        tabpfn_cv.py             transversal; every axis may call it

plotting/         TOP-LEVEL, cross-axis, presentation-only. parents[1], NOT parents[2].
                  filters.py = the shortlist predicates, defined once
                  palette.py = npg palette + article style, no ersilia branding
                                                            -> output/plots/presentation/

src/              flat. one module per task + matrices.py, tabpfn.py, interest.py, proteomelm.py
docs/             one .md per task, named for the task (docs/function.md, not docs/02_function.md)
tools/            one-shot migration scripts, kept for the record

data/source/<provider>/       uniprot deg eggnog orthodb cdd go sprofgo geptop ncbi unknome
                              alphafold alphafill biolip wwpdb sifts pdbe ecmdb plinder
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

**`data/source/` is by PROVIDER, not by task** — the UniProt proteomes feed every axis. Where a file
came from is stable; which task wants it is not. **Only deliverables sit at a task root.**

**`scratch/` is safe to purge and pointless to upload** — ~18 GB of caches against ~90 small files in
`evidence/`. But anything you would **cite** cannot live there: a routine purge once deleted the only
record of which Geptop references contributed and cost a 90-minute re-run.

Two names that look like they belong to one tier and do not: **geptop's `cv_*.pkl` is a COMPOSITION
VECTOR cache** (scratch), not cross-validation; and a `{pre}`-prefixed file is evidence on a real run
and scratch under `--limit`, because `pre` is the smoke prefix.

Naming: `<artifact>_<species>.<ext>` — the artifact type leads, so a directory sorts by kind.

### Four rules the layout depends on

1. **`REPO_ROOT = Path(__file__).resolve().parents[2]`** in every v2 script — they are one level deeper than they used to be. **`src/` keeps `parents[1]`** and **`legacy/` keeps `parents[1]`** (its `data`/`output` symlinks depend on it). Deliberate; do not "fix" them.
2. **Plots live in their own top-level folder**, because figures are often comparative *across* categories and would otherwise have no home.
3. **Workers are filed by *how* they run, not by who calls them.** A worker is an entry point under a *different* conda interpreter, never imported. They do not go in `src/`: importing one would drag its foreign dependencies into the `gradi` process, the whole thing the env split prevents.
4. **A tool's own embeddings belong to the tool.** SPROF-GO's ProtT5 representations live under `function/scratch/`, TMbed's under `localization/`. `embeddings/` is only for the three first-class, project-wide matrices — see *External models*.

### Two documented exceptions — leave them alone

- **`data/raw/` still exists** and is the frozen v1 archive, referenced directly by v1 scripts through the `legacy/data -> data/` symlink and indexed by `data/raw/PROVENANCE.md`. v1 used an organism-first layout; do not extend it.
- **`ligands/chembl.py` and `ligands/bindingdb.py` read their bulk dumps from `data/raw/other/{chembl,bindingdb}/`**, not `data/source/`. Those 11.9 GB are shared with v1 scripts using the same path.

## Pipeline (`scripts/`)

One folder per task. Run everything in the `gradi` env unless an entry says otherwise; where a stage
needs another env it crosses a **process boundary** and never asks you to activate it.

**Each entry keeps only what you need before opening anything: what it writes, which env, and the
traps that fail silently. The measured tables, the rejected alternatives, the controls, the CLI and
the run log are in `docs/<task>.md`** — named at the end of each entry, and that doc is the record.

- **`proteomes/download.py`** → `data/processed/proteomes/proteome_<species>.tsv`, four tables,
  **5 columns**, keyed on `uniprot_ac`: `gene_name` · `protein_name` · `sequence` ·
  **`proteome_evidence`** (owner's call, 2026-10-03 — it was 9).

  **`proteome_evidence` (1–3) REPLACED `is_reviewed`** (2026-10-04), which was nearly degenerate
  per species: Ec and human are **100% reviewed**, Kp is **7 of 5,728**, so on three of four
  proteomes the boolean separated nothing. **3** = the entry carries its OWN identity
  (SwissProt-reviewed, **or** a gene symbol on the anchor entry itself) AND a specific protein name
  · **2** = one of those · **1** = neither. Kp 2,205/2,469/1,054 · Ec 0/648/3,755 · Sa
  1,510/485/894 · human 0/490/19,926. **Ec and human have NO level 1 and that is correct** — every
  entry is curator-read. **The two FILLED naming tiers deliberately do NOT count as own identity**:
  that is this axis's own inference (Kp 18.4% → 63.4%), not the protein's record. Generic name =
  `uncharacterized|hypothetical|unknown function|DUF\d+` — **DUF is Domain of Unknown Function**,
  while `Oxidoreductase` and `N-acetyltransferase domain-containing protein` are real classes.
  `is_reviewed` survives byte-identically in `evidence/proteome_full_<species>.tsv`. **No
  `proteome_consensus`** — identity is not a magnitude. Rule in `src.proteomes.identity_evidence`,
  shared by the writer so there is ONE definition. **`sequence` STAYS in the deliverable deliberately**: the
  standing rule is *map by sequence, not by accession*, so the column every external join needs
  belongs in the table everything loads. The four provenance columns — `gene_name_source`,
  `gene_synonyms`, `refseq`, `geneid` — are in **`evidence/proteome_full_<species>.tsv`** via
  **`load_full()`**, and **`geneid` is NOT idle there**: `gene2pubmed.py` and `pubtator.py` key NCBI
  literature counts on it, and `id_bridge()` reads it. Audits and xrefs are `evidence/`, the
  UniRef90 cache `scratch/`. Registry-driven from `src/proteome_registry.tsv` (41 rows, 4 tiers,
  every row carrying a required `why`).

  Its real job is the naming gap — **Kp 18.4% → 63.4%**, Sa 28.2% → 44.6% — filled in three labelled
  tiers, all **within the species**. Cross-species naming is deferred to a stage with real orthology
  that can label it as the inference it is. **Load through `src/proteomes.py`** rather than adding
  stacked artifacts: a first run wrote a `proteins.parquet` and an `id_bridge.tsv` that were pure
  derivations, 27 MB for zero new information.

  Registry tiers: **A** the 4 anchors · **B** same-species name-donor pools · **C** a 26-species
  comparator panel · **D** evidence-bridge strains from **NCBI**, because **UniProt serves zero
  entries for any proteome it does not flag "Reference proteome"**. C and D are fetched on request.
  Details: `docs/proteomes.md`.

- **`embeddings/esmc.py`** → `embeddings_<species>.npz`, ESM-C 600M mean-pooled, 1152-d, BOS/EOS
  stripped, three bacteria, with a resumable shard cache in `scratch/shards_esmc/` whose shard size
  is in the filename so a change cannot reuse mismatched shards. **Human is deliberately out**: 75%
  of all residues and every pathological length, and adding it needs a chunking policy this stage
  does not have. **Load through `src/embeddings.py`** — `accessions` needs `allow_pickle=True`, and
  `vectors_for()` drops missing accessions rather than zero-filling, since a zero row is a real
  position in embedding space.

- **`embeddings/prott5.py`** → `prott5_<species>.npz`, 1024-d, the second project-wide matrix.
  Inference crosses into **`gradi-loc`** (`embeddings/workers/prott5.py`, `GRADI_LOC_BIN`), because
  TMbed already ships the weights there.

  **The control is the point, not an extra**: E. coli MG1655 is both our anchor and one of the
  proteomes UniProt publishes precomputed ProtT5 vectors for, so the same model over the same
  sequences must reproduce them — **median cosine 1.000000, worst 0.999994** over a 300-protein
  length-weighted sample, and nothing is written if that fails. It also settled a convention UniProt
  leaves unstated: ProtT5 has EOS and no BOS, and **`mean_no_eos` beat `mean_with_eos` 0.999994 vs
  0.998837** — measured, not assumed.

  **This is the half, encoder-only build** (`prot_t5_xl_half_uniref50-enc`, *different weights* from
  the full `prot_t5_xl_uniref50`): valid for OUR models, never substituted into an external
  predictor. `--strain` writes per-screen-strain matrices to `scratch/strains/` — training features,
  not deliverable matrices, with no canonical row order.

- **`embeddings/proteomelm.py`** → `proteomelm_<species>.npz`, the third project-wide matrix and the
  only one where a protein's vector depends on the rest of its proteome: ProteomeLM-L over a whole
  proteome's ESM-C vectors in one pass, **layer 8 of 18** (the paper's own best configuration),
  z-scored genome-wide.

  **It computes its OWN ESM-C and must** — ProteomeLM pools over non-pad tokens, so BOS/EOS are
  *inside* the mean while `esmc.py` strips them; stage 01's npz is cross-checked and **reported,
  never depended on**. **NOT SHARDABLE**: a shard boundary changes the values, because the proteome
  *is* the context.

  **`--group-embeds {self,orthodb}` is the load-bearing knob.** In training, `group_embeds` is the
  mean ESM-C embedding of the protein's OrthoDB group, added as a second branch, so it changes every
  hidden state. `self` (default) passes `None` — the authors' released *inference* default, not what
  the model was trained with. `orthodb` uses the real thing via our own `orthodb_<species>.tsv`,
  because UniProt's `xref_orthodb` is **0.0% on Kp**. Mapped: ecoli 88.2% · kpneumoniae 81.7% ·
  saureus 76.9%; the rest fall back to their own vector, so **always report that fraction** — where
  it is low the two modes converge and a null result says nothing.

  **Run `embeddings/orthodb_group_check.py` before trusting `orthodb` mode.** OrthoDB group ids are
  not stable between releases; had the id spaces differed, every lookup would miss and the run would
  **silently reproduce `self`** while looking like a completed experiment. Measured and passed.

  Four traps, each hit at least once: the script **shadows the installed `proteomelm` package**, so a
  bare `import proteomelm` resolves to the script and dies with `'proteomelm' is not a package`,
  which reads like a broken install; **the output filename must carry the mode**, or an `orthodb` run
  overwrites the `self` baseline; **both controls must carry the group tensor**, or permuting
  re-pairs each protein with another's group vector and the control fails for the wrong reason;
  group vectors are **`bfloat16`** against a `.float()` forward. Also: **the four
  `group_vectors_*.pkl` are DISJOINT SIZE BANDS, not nested supersets** — a *lower* `min_group_size`
  loads *more*, and `_0` is 18.1 GB unpickled whole.

  Two controls exit non-zero: permutation invariance, and **context sensitivity** — the full- vs
  half-proteome cosine must stay *below* 0.999, or ProteomeLM has collapsed to a per-protein encoder
  and the stage has no reason to exist (Kp 0.9947 and Sa 0.9927 sit close enough to watch).
  **Load through `src/proteomelm.py`, which ASSERTS the file's recorded mode matches the one
  requested**: the two modes are the same shape over the same accessions, so nothing about a matrix
  says which it is, and fitting on one while scoring on the other returns plausible wrong numbers.
  Details: `docs/embeddings.md`. Provenance: `data/source/proteomelm/SOURCE.md`.

- **`embeddings/projection.py`** → `projection_<species>.tsv`: `uniprot_ac · tsne_x · tsne_y`, 100%
  coverage, **no colour column** — join anything else on `uniprot_ac`. **This axis ships NEITHER
  standard column** (owner's call, 2026-10-04): a t-SNE coordinate has no magnitude to rank and no
  evidence to grade — it is a deterministic reduction of an embedding, not a claim about the
  protein — so `<axis>_consensus`/`<axis>_evidence` would be constant for all 13,020 rows. **The
  gap is deliberate; do not "complete" it.**

  **The recipe is inherited from v1, not re-derived** (openTSNE multiscale, cosine, PCA-50, `dof=0.8`,
  winner of a 33-config sweep). **Do not re-run that sweep** — the code was deleted, so the verdict
  table in `docs/embeddings.md` is the record. `dof` is the load-bearing knob, not perplexity. The
  control is **`trustworthiness` (k=10), floor 0.90**, because a diverged t-SNE still returns
  plausible numbers.

  **Coordinates are per-species with no shared frame** — never compute a distance across species on
  these columns; use the 1152-d space. `n_jobs=-1` means the floats are not bit-reproducible, only
  the structure, which is why the guard is trustworthiness and not a checksum. `coords_for` drops
  missing accessions rather than zero-filling: **(0, 0) is a real position**, in the dense centre.

- **`localization/predict.py`** + **`localization/merge.py`** → **`localization_<species>.tsv`**,
  **4 columns**, **from sequence alone**: `localization` (DeepLocPro 1.0),
  `cytoplasmic_fraction` (TMbed) and **`localization_evidence`** (1–3).

  **`localization_evidence`**: **3** = the two predictors CONCUR, DeepLocPro is confident (≥0.7)
  **and** a curated GO cellular-component term agrees · **2** = one of those · **1** = neither, or
  a GO term **contradicts** the call. Kp 1,009/3,200/1,519 · Ec 624/2,098/1,681 · Sa 589/1,606/694.
  **LEVEL 3 IS NOT "EXPERIMENTALLY LOCALIZED"** — nothing in this axis is an experiment. **The GO
  term is the only signal that is not a sequence model** (both predictors read the sequence), which
  is why it is worth the cross-axis read: it is present for 29.6–48.4% of proteins and agrees with
  DeepLocPro **90.8–94.8%** where present. `GO_CC_TO_COMPARTMENT` in `src/localization.py` maps the
  six goslim CC terms 1:1 onto the six classes. The weak class sinks as it should: `extracellular`
  reaches 3 for **0.4% on Kp** against cytoplasm's 15.1%.

  **`merge.py` now READS `function_<species>.tsv`** — the one cross-axis dependency in this stage,
  and it **exits non-zero if absent** rather than quietly computing a two-signal ladder under the
  same column name. Run `function/matrix.py` first.

  **Still no `localization_consensus`, and the reason is stronger than "no magnitude"**: collapsing
  a compartment and a fraction into one 0–1 number IS an accessibility score, which this axis
  deliberately does not compute — v1's `clp_accessibility` ladder was consumed by nothing. **`confidence` ships in `evidence/deeplocpro_<species>.tsv`**
  (owner's call, 2026-10-03), byte-identical — but DeepLocPro **always** returns a call, so the
  label now reads equally authoritative for every protein and nothing in the table says which calls
  are weak: **12–15% sit below 0.7 confidence, 2–4% below 0.5.** Join it back before trusting one
  label. **`has_signal_peptide` ships in `evidence/tmbed_<species>.tsv`**
  (owner's call, 2026-10-03), byte-identical. It is the one column that says **WHY** a fraction is
  near zero — exported rather than membrane-buried — which is exactly the degradability mechanism
  (a secreted protein transits the cytoplasm unfolded and IS reachable; a membrane protein never
  does), so `degradability/enrichment.py` joins it back explicitly. Env **`gradi-loc`** across a
  process boundary (`localization/workers/deeplocpro.py`, `GRADI_LOC_BIN`); the two tracks run
  **sequentially, never concurrently**. DeepLocPro always returns a call, so **100% coverage is a
  property of the method** and there is no `unknown` class.

  **ONE table, two predictors side by side — and `merge.py` arbitrates NOTHING.** They answer
  different questions and neither derives from the other, so where they disagree the disagreement
  is the information; TMbed corroborating `extracellular` from outside DeepLocPro is the only
  cross-check this axis has. **The two per-predictor tables live in `evidence/`**, not because they
  are secondary but because the tracks run and resume independently — TMbed alone is CPU-only hours
  — so `--only tmbed` must have somewhere to write. `merge.py` recomputes nothing, seconds, `gradi`.

  **NO `evidence` COLUMN HERE, and that is deliberate** (owner's call, 2026-10-03): both predictors
  cover 100% of every proteome **by construction**, so it was constant across all 13,020 proteins
  and said nothing. The axis-wide rule assumes an axis that can fail to annotate a protein; this one
  cannot. The check survives — `merge.py` **exits non-zero if either track is silent for a single
  protein**, which is a broken run, not a sparse one. **Do not re-add it as a constant.**

  **The Gram-positive trap — read before comparing *S. aureus*.** `positive` mode does not merely
  mask the two Gram-negative-only classes, it **adds their probability mass into `Extracellular`**.
  So Sa has four reachable classes, its `extracellular` count absorbs whatever the model wanted to
  call periplasmic, and a zero in the OM/periplasm columns is **structural, not missing data**.

  **E. coli K-12 is almost certainly in DeepLocPro's training set**, so any E. coli agreement number
  is a sanity check, not validation — Kp and Sa are the honest test sets. **`extracellular` is the
  weakest class**, which TMbed corroborates from outside the model; **prefer TMbed's fraction over
  DeepLocPro's label** where a choice exists. **Deliberately not computed here:** any accessibility
  score — that is a modelling decision for whichever stage consumes it, and v1's `clp_accessibility`
  ladder was consumed by nothing. Details: `docs/localization.md`.

- **`function/cog.py`** → `cog_<species>.tsv`, one of **26 COG2024 category letters** per protein via
  COGclassifier 2.0 over the CDD COG profile DB. **The lookup route is not available**: UniProt's
  `eggnog` xref is **0.00% (n=0) on Kp HS11286**.

  **Coverage is not 100% and must not be forced to be** (Kp 79.1% classified / 72.8% informative, Ec
  84.5% / 78.0%, Sa 73.1% / 65.7%; *informative* excludes `R` and `S`). NCBI's own curators reach
  81.6% of E. coli K-12, so there is no headroom, and **raising `--evalue` buys noise** — on shuffled
  decoys the hit rate goes 0.4% → 24% → 84.3% at 1e-2 → 1 → 10. The default sits at the knee.

  **`rpsblast` is borrowed from `gradi-prokka`** via `PATH` (`GRADI_RPSBLAST_BIN`); do not
  `conda install blast` into `gradi` — no osx-arm64 build, so it flips the env and takes ESM-C with
  it. **Human is excluded by construction**: COG2024 has zero eukaryotes. `cog_category` is one
  letter and `cog_category_all` the full string, because 9–16% carry a multi-letter COG and taking
  the first is a real choice.

- **`function/eggnog.py`** + **`function/goslim.py`** → GO slim in `goslim_prokaryote` (97 terms),
  two tiers of real annotation: `curated` → `eggnog` (**eggNOG-mapper v2**, env **`gradi-emapper`**,
  osx-64/Rosetta). Also writes `eggnog_<species>.tsv`.

  **Coverage is NOT 100%, deliberately.** An earlier version reached 100% by k-NN transfer in ESM-C
  space; it was removed on instruction to use well-established tools only, and that was right twice —
  it is not citable, and its honest reweighted accuracy was **MF 68% / BP 53% / CC 66%**. A protein
  no established tool can annotate gets an empty row. **Do not reintroduce it.**

  **`interpro2go` is a scored REJECTED ALTERNATIVE; do not re-add it** — as a tier it filled **5
  proteins out of 13,020**, because UniProt's electronic GO is already InterPro2GO-derived.

  **The gain is small** — Kp 69.7 → 73.7%, 322 proteins out of 13,020 — and the bottleneck is not
  orthology: eggNOG places ~73% of unannotated proteins in a group but only **2–18% of those groups
  carry any GO**. **Test the OG→GO yield before paying for a database this size.** Kept because the
  same run yields what nothing else does (KEGG KO 64.1% on Kp, `pfams` 85.4%, `description` 86.8%).

  **Do not swap COGclassifier for emapper's `COG_category`**: against NCBI's curated COG2024 on
  E. coli, COGclassifier agrees 97.8% and emapper 63.3% — a different vocabulary. The database
  unpacks to **50.6 GB** and **was deleted after the run** (recovery:
  `data/source/eggnog/SOURCE.md`). **`eggnog_<species>.tsv` is demoted, not deleted** — the goslim
  `eggnog` tier derives from its `gos` column, so deleting it makes the stage non-regenerable. Needs
  `goatools`.

- **`function/matrix.py`** → **`function_<species>.tsv`**, **4 columns**: `uniprot_ac` ·
  `cog_categories` · `goslim_terms` · **`function_evidence`**, both term columns `;`-joined in vocabulary order (owner's call,
  2026-10-03 — it was two wide matrices). **Recomputes nothing** — it reshapes; seconds, no database.

  **The matrices still ship, in `evidence/`** (n × 99 and n × 28) **and are still audited**, because
  they are what carries the **structural zeros the packed form cannot express** — 8 GO terms are
  eukaryote/plant, Sa has 19 unreachable as a Gram-positive, COG `Y` is nuclear structure. The stage
  asserts **BOTH directions**: the matrices round-trip to the long-form source, AND the packed
  columns re-expand to the matrices exactly. Use `load_goslim_matrix()` / `load_cog_matrix()` for a
  feature matrix or to tell "impossible" from "unknown".

  **`function_evidence` (1–3), and NO `function_consensus`** — the first axis where only half the
  standard pair applies. *How much function does a protein have* is not a quantity; the nearest
  candidate, annotation richness, is a STUDIEDNESS measure. **The convention is evidence always,
  consensus where the axis has a magnitude.** **3** = both schemes annotate it, the GO is
  UniProt-**curated** and the COG is **informative** (a letter outside `R`/`S`) · **2** = annotated
  but not corroborated · **1** = neither scheme. Kp 1,019/1,012/3,697 · Ec 328/776/3,299 · Sa
  672/535/1,682.

  **Level 2 is NOT "badly annotated"** — a protein with excellent curated GO but no COG hit caps
  there, and COG tops out near **81.6%** by NCBI's own curators, so a missing COG is usually the
  method's ceiling. **It is not a fame measure either**: 210 E. coli proteins named
  "Uncharacterized" sit at 3 (32.4% of that group, against Kp's 3.3%), because UniProt leaves them
  unnamed while curating their class — "Uncharacterized MFS-type transporter" carries GO:0005215.
  **GO EVIDENCE CODES WOULD BE BETTER AND THIS REPO HAS NONE**: no GAF anywhere, `go_id` is a bare
  list — obtaining them is a per-proteome GOA download, not a reshape. Don't go looking.

  **An EMPTY list means NOT ANNOTATED**, never "ruled out" — 1,506 Kp proteins (26.3%) carry no GO
  term because nothing is known about them.

  **NO EVIDENCE COLUMN, for either scheme.** COG never had one to carry: `cogclassifier` vs `none`
  is 1:1 with non-empty vs empty on all three species — measured, not assumed. GO-slim's was real
  but separated only **322 proteins of 13,020** (Kp 231 · Ec 75 · Sa 16) and survives
  byte-identically as `evidence` in `evidence/goslim_matrix_<sp>.tsv` and `goslim_source` in
  `evidence/goslim_<sp>.tsv`. **So the shipped table does NOT say whether a GO term is
  UniProt-curated or inferred from an eggNOG orthogroup** — read one of those two before treating a
  term as curated. Built from the `*_all` columns so **multi-label is preserved** (max 11 slim terms
  on Kp). **The COG vocabulary is VENDORED** at
  `data/source/cdd/cog_func_category.tsv`: a schema that depends on a pip install is not
  reproducible. Details: `docs/function.md`.

- **`degradability/predict.py`** → `degradability_<species>.tsv`, **6 columns**: `uniprot_ac` ·
  `adep4_prob` · `onc212_prob` · `nn_similarity` · `degradability_consensus` ·
  `degradability_evidence`. Is this protein a substrate of activated
  partnerless ClpP? A binary **TabPFN-3.5** classifier on ESM-C embeddings, trained on the two
  *S. aureus* activator screens and applied to the rest of Sa and all of Ec/Kp. `<act>_prob` is the
  **seed-averaged out-of-fold** probability for every protein, so it is one comparable scale across
  all 13,020.

  **The standard pair, added 2026-10-04** — the first axis with both a magnitude and an experiment,
  so unlike function and proteomes it ships both. **`degradability_consensus` ADDS VERY LITTLE and
  must not be sold as corroboration**: it tracks either activator alone at **rho 0.970–0.973**,
  because the two probabilities are themselves at rho 0.884. `nn_similarity` is deliberately **not**
  an input — it measures reach, not degradability, and belongs to the evidence column.

  `degradability_evidence`: **3** measured, in an activator screen · **2** predicted inside a band
  with a validated AUROC (`nn_similarity >= 0.90`) · **1** predicted beyond any validated band.
  Counts Kp **1,175 / 4,553 / 0** · Ec **541 / 3,862 / 0** · Sa **221 / 959 / 1,709**.
  **Kp and Ec cannot exceed 2 and that needed no special-casing** (owner's constraint): the screens
  are Sa-only, so the cap falls out of the ladder — asserted anyway, so a future non-Sa label set
  widens it loudly. **The 0.90 cut is READ OFF `evidence/domain_bands.tsv`, not chosen** — the two
  bands below it carry a NaN `roc_auc`, so level 1 means *never validated this far out*;
  `VALIDATED_SIMILARITY` derives from `SIMILARITY_BANDS[2][0]` so the two cannot drift. **HERE A 3
  IS ONE MEASUREMENT, deviating from the generic ladder in `src/consensus.py`** ("≥2 sources,
  unanimous, concordant"), because the two activators are not two sources: an activator a protein
  was never tested against is missing data, not dissent. The 823-agree / 190-disagree split among
  the 1,013 measured by both lives in `evidence/consensus_audit.tsv`.

  **"In the labels file" is NOT "measured" — the gap is 162 proteins.** `labels_saureus.tsv` holds
  1,871 rows of which **162 are NaN for BOTH activators**, sequence-mapped but never called, so
  level 3 reads 1,709. A membership test against the file would promote them to "measured" on the
  strength of a successful join alone.

  **THE MEASURED CALLS ARE NOT COLUMNS** (owner's call, 2026-10-03). `<act>_hit`/`<act>_source` were
  dropped: `_hit` was empty for **10,131 of 13,020** proteins and `_source` read `predicted` for all
  but 1,871. **Nothing was lost, and it was CHECKED, not assumed** — both are exactly reconstructible
  from `evidence/labels_saureus.tsv`, verified per species and per activator before the tables were
  rewritten. That file is strictly richer: it carries the continuous log2FC (a measured −0.51
  explains a `0` in a way the bit cannot) and the grouping cluster. **`src.degradability.measured()`
  and `with_measured()` put the calls back**, and **`hits()` is unchanged in behaviour** — a
  measurement still beats the model, it just reads the labels instead of a column.

  **This overturns v1's "needs new data, not new features"**: cluster-grouped CV over 5 seeds gives
  **ADEP4 0.8738 ± 0.0029 (PR 0.6103)** against a properly-estimated length baseline of 0.776, and
  **ONC212 0.7671 ± 0.0060 (PR 0.5803)** against 0.687. **Report mean ± SD over 5 seeds, never one** —
  a single-seed estimate carries ±0.004 of pure arbitrariness, larger than any effect this stage
  tested and rejected. **`CROSS_ASSAY_AUROC` = 0.877 / 0.815 is the yardstick, not 1.0**, and
  deliberately **not** called a ceiling: the model exceeds it at three of five cutoffs.

  **The log2FC cutoff must NOT be tuned** — AUROC is not comparable across labels, so tightening it
  inflates the score mechanically. The sweep ships anyway because `CROSS_ASSAY_AUROC` rises in
  lockstep: the **gap** is the invariant.

  **The labels join to NOTHING by identifier** (**0 of 1,943** by RefSeq, locus tag or Jacques's own
  UniProt column), so it is a **DIAMOND sequence join** at 96.4%, exiting below a 90% floor. The
  label is v1's audited `*_bin`, **asymmetric on purpose**. Unmeasured is NaN, never a negative.
  **Cleavage is deliberately unused. Two activators, never merged.** The artifact is
  `model_<activator>.npz`, **not `.joblib`** — TabPFN learns in context, so there is nothing to
  pickle; what determines the predictions is the training set + config.

  **The two probability columns are NOT independent evidence**: they correlate at **ρ 0.89** while
  the labels agree at ρ 0.52 / Jaccard 0.32 — both read the same embedding, so "both activators
  agree" is closer to one opinion than two.

  **0.5 is the wrong threshold and the right one MOVED with the estimator.** Probabilities top out at
  0.868 on Kp, so `≥0.5` now selects 220 Kp proteins where it selected 20. `src.degradability.hits()`
  defaults to the re-derived base-rate cuts (**0.328 / 0.313**, replacing the forest's 0.394 / 0.426)
  — empirical quantiles of the OOF distribution, so **re-derive them on any estimator change**. They
  reproduce the base rate on the LABELED set, not a proteome: whole-proteome the same cut takes
  14.4–20.9% (adep4) / 31.4–34.6% (onc212). Better still, rank on `_prob`.

  ***E. coli* and *K. pneumoniae* rows are ranking hypotheses, not measurements.** `nn_similarity` and
  `evidence/domain_bands.tsv` price the extrapolation (ADEP4 Kp 0.814 / Ec 0.817), but **a fifth of
  Kp sits in bands with no estimate at all**, and the premise that ESM-C cosine measures
  transferability is unvalidated and untestable without Gram-negative labels.

  **A mechanistic finding for downstream filtering**: hit rate is cytoplasm 0.180/0.275 but
  **membrane 0.031/0.105** and **extracellular 0.049/0.346**. Membrane proteins are protected
  (co-translational insertion, never a soluble cytoplasmic chain); secreted ones are not (they
  transit the cytoplasm unfolded). So **"cytoplasmic only" is the wrong filter** — "not membrane" is
  closer.

  **`degradability/enrichment.py`** tests COG and localization enrichment of the full-proteome
  predictions, with **hits = the top 10% of each proteome**, not a fixed count: the proteomes differ
  2× in size, so an absolute N makes the odds ratios incomparable. **Judge enrichment against the
  LABELS, not against the previous model — this was got wrong once**: two of the forest's five
  "consistent" categories were artifacts, and using the incumbent's behaviour as the acceptance test
  silently assumes the incumbent was right. **The clearest case is `extracellular`, which the forest
  depleted and TabPFN ENRICHES** (OR 2.08–3.80, significant on all three species under ONC212) —
  tracking the measured 0.346 hit rate above cytoplasm's 0.275, i.e. the same "secreted proteins
  transit the cytoplasm unfolded" mechanism. Two standing caveats: **B chromatin is the largest OR
  (64) and the weakest evidence — 8 members, quote the count, not the ratio**, and both estimators
  track length far more strongly than the labels do (−0.61 vs −0.33).

  Figures (stylia, one at a time) in `scripts/plots/`: `degradability.py`,
  `degradability_predictions.py` (**read `extrapolation.png` panel B before quoting any expected
  AUROC**) and `degradability_top.py`, which shows that **the two activators' top-100 lists share
  only 24/100 on Kp** despite ρ 0.89 — *global rank correlation is not shortlist agreement*.

  **Load through `src/degradability.py`**; the name collides with the frozen
  `legacy/src/degradability.py`, which v2 never imports. Details: `docs/degradability.md`.

- **`orthology/orthofinder.py`** → **`orthology_<species>.tsv`**, **4 columns**: `uniprot_ac` ·
  `has_human_ortholog` · `bacterial_panel_orthologs` · **`orthology_evidence`** (owner's call,
  2026-10-03 — it was the 29-column dense table). Also **`orthologs.tsv`** (sparse; `is_ortholog_orthofinder` and `is_rbh` **side by
  side, never merged**), **`neighbors.tsv`** (sparse, top-5 nearest neighbours per target species
  with identity and both coverages; self-hits dropped, so the within-species block is a protein's
  nearest *paralogs*) and **`evidence/orthology_<species>.tsv`**, the DENSE table, via
  `load_dense()`.

  **`orthology_evidence` (1–3), and NO `orthology_consensus`** (2026-10-04). **The missing
  consensus is not just "no magnitude"**: the two columns point in OPPOSITE prioritization
  directions — you want **no** human ortholog and **broad** bacterial conservation — so a mean is a
  weighting, not a summary, and this file forbids one outright. **3** = placed in a grouping AND
  both OrthoFinder and RBH found a bacterial ortholog AND the two do not conflict on the human
  call · **2** = placed but not corroborated · **1** = in NEITHER an OrthoFinder orthogroup nor an
  OrthoDB group, so both columns beside it are *could not look*. Kp **291/2,639/2,798** · Ec
  **110/1,537/2,756** · Sa **216/1,875/798**.

  **IT REQUIRES A POSITIVE FINDING, NOT ONLY A RELIABLE MEASUREMENT** (owner's call), so ***S.
  aureus* reaches 3 for only 27.6%** — the lone Gram-positive among the anchors, the same effect
  behind its 0.179 median. **A low level on Sa is partly its biology, not only our uncertainty.**
  **Nothing in this axis is an experiment**: a 3 is not experimental corroboration and a 1 means
  *could not look*, the pattern `function_evidence` set.

  **Measured, not assumed**: `in_orthogroup` **uniquely blocks 0 proteins at level 3** (implied by
  the OrthoFinder term; load-bearing only for level 1) and `in_orthodb` blocks 10/3/0 — the
  discrimination is the **RBH** term (301/189/218) and the **human conflict** (312/301/126). **RBH
  reaches only the other two ANCHORS, not the 28-species panel**, so the test is narrower than the
  column it grades. Of Kp's 951 human calls only **501 are found by both methods**; Ec **`tufA`**
  is the case to remember — conserved in all 28 and still capped at 2, because OrthoFinder calls a
  human ortholog and RBH does not (EF-Tu vs mitochondrial TUFM). **OrthoDB CANNOT be a third
  opinion on the human call and this was VERIFIED**: `<n>at2` and `<n>at2759` share **zero** ids.
  Conditions per protein in `evidence/evidence_audit.tsv`.

  **`orthology_evidence` matches NONE of `_read()`'s prefixes**, so it is registered in the Int64
  branch and in `DELIVERABLE_DTYPES` by name — the third column to hit that trap after
  `n_bacterial_orthologs` and `bacterial_panel_orthologs`, and the one check that catches it is
  that `.mean()` returns a float rather than a concatenated string.

  **`bacterial_panel_orthologs` IS A FRACTION (0–1), NOT A COUNT** — the name reads like a count,
  so check the scale. **It is OVER 28, NOT 26** — the 26 tier-C comparator proteomes plus the
  three bacterial anchors, minus this protein's own species. It counts **SPECIES**, never proteins,
  so a paralog pair does not inflate it. Median **0.536 Kp · 0.571 Ec · 0.179 Sa** — the last is
  *S. aureus*'s Gram-positive isolation against a mostly Gram-negative panel, not a defect.
  **Human has NO deliverable** (the panel is bacterial) and `load()` refuses it by name. All
  four proteomes in one joint run. v1 shipped only the first, with its identity and coverage columns
  *entirely empty*, so any axis thresholding transfer on identity silently dropped everything.

  **A 0 here is a MEASURED 0**: OrthoFinder runs **de novo on our own FASTAs**, so every protein is
  either in an orthogroup or named in `Orthogroups_UnassignedGenes.tsv`, and the stage exits unless
  that accounts for each proteome exactly. A *sparse* matrix cannot express a zero at all — that is
  what the dense table is for.

  **`--panel full` runs the 26-species tier-C comparator panel alongside the anchors** — 30
  proteomes, **140,396 proteins, 47 min**. The comparators inform the orthogroups and get no tables
  and no pairwise DIAMOND (16 searches, not 900); they exist so the conservation count means
  something. Default is `anchors`, the prior behaviour.

  **The two columns the panel is for**, in `orthology_<species>.tsv`:
  **`has_human_ortholog`** (union of both methods) and **`n_bacterial_orthologs`** — how many of
  the 28 bacterial proteomes share this protein's orthogroup, counted over SPECIES not proteins,
  with `bacterial_panel_size` beside it so 12-of-28 cannot be misread as 12-of-3. Validated against
  biology: `ftsZ`/`gyrB`/`dnaA`/`secA` **28**, `rpoB` 27, `clpP` 16, `lacZ` 5–8. **Sa's median is 5
  against Kp/Ec's 15–16** — correct, it is the lone Gram-positive in a mostly Gram-negative panel,
  not a defect. It is **orthogroup-based, not the union**: RBH to 26 comparators would need 156
  more DIAMOND searches (~5 h) for a second opinion on the same question.

  **OrthoFinder's recall depends on panel size; RBH's does not — know this before quoting a number.**
  v1's much-cited 55.5% came from a 25-species run; at four species OrthoFinder gives 44.8%, *not a
  defect*, and RBH is the like-for-like comparison that passes (3,074 vs v1's 3,003).
  **`--very-sensitive`, and the error direction is the reason**: under-detecting human homology makes
  a target look *more selective than it is*. Identity and coverage are columns, never filters.

  **WIDENING THE PANEL MOVES BOTH DIRECTIONS AT ONCE, and the human one goes DOWN.** Measured in
  `evidence/panel_expansion.tsv`, written every run:

  | | in_orthogroup | has_human_ortholog |
  |---|---|---|
  | Kp | 4,322 → **5,289** (+967) | 1,184 → **951** (−233) |
  | Ec | 3,646 → **4,247** (+601) | 980 → **838** (−142) |
  | Sa | 1,754 → **2,262** (+508) | 700 → **624** (−76) |

  **An earlier version of this plan asserted a human drop would be a regression. That was wrong**,
  and the reason is worth keeping: "recall rises with panel size" is a **bacteria↔bacteria**
  statistic, and the panel adds 26 bacteria and **zero eukaryotes**. More bacterial resolution lets
  OrthoFinder separate orthologs from out-paralogs on the bacterial side, so marginal
  bacteria↔human calls are withdrawn. The evidence that it is refinement, not breakage: **at ≥60%
  identity 13/13 Kp, 3/3 Sa and 11/12 Ec human relationships survive**, at ≥50% it is 56/60, 23/23,
  57/59 — the losses sit at 40–50%, where ortholog-vs-paralog is genuinely ambiguous (median
  identity of the kept 33.5% vs 28.8% for the dropped).

  **This run is also the clearest argument for keeping both methods.** Kp `clpP` → human CLPP
  (56.3%) is documented as an ortholog by *both*; OrthoFinder now calls it **`of=0`** and only RBH
  still finds it. At ≥50% identity RBH rescues 9 Kp proteins OrthoFinder misses. The union holds
  where the single methods would not.

  **Two traps.** **OrthoFinder exits 0 when its dependency check fails** — no `diamond` on `PATH`
  gives an ERROR block, an empty `Results_` dir and a *success* return code, so `ensure_diamond()`
  must run *before* it. And the **launcher cannot be called directly**: its `#!/usr/bin/env python3`
  resolves to the unrelated `ersilia` env and dies on `ete4`, which **is installed** — name
  `<env>/bin/python` explicitly rather than chasing the wrong bug. `SPECIES` here is all **four**.

- **`orthology/orthodb.py`** → `orthodb_<species>.tsv` for all four proteomes, alongside
  OrthoFinder's and never instead: OrthoFinder's orthogroups are *de novo* and panel-dependent, so a
  grouping that does not move when the query changes must be defined elsewhere. Coverage at the
  default `--reps 20`: Kp 92.7% · Ec 95.5% · Sa 91.8% · human 96.0%, with 100% *verdict* coverage.

  **v12 REPLACES v11; the two id spaces must never be joined** — an OG id is *"not stable and re-used
  between releases"*. Every row carries `orthodb_version`.

  **Filter on `orthodb_confidence`, NOT on identity** — on the E. coli control >0.9 is 98.0% correct
  and <0.5 only 41.4%, while 95–100% identity alone is just 79.0%. But read it **with
  `orthodb_n_candidate_ogs`**: a single candidate takes 100% of the share trivially, and `n_cand==1`
  at 25–45% identity drops to 80.0% — a lone candidate found at low identity is the one shape to
  distrust. The floor is **25% identity / 50% coverage**, decoy-calibrated; the project's ≥40% rule
  is for annotation **transfer**, a stricter task, so do not conflate them.

  **Do not compare two proteins by group id alone — check the name.** OrthoDB maintains parallel
  Bacteria-level groups for the same family, and **57.1% of the apparent sequence-tier errors are a
  parallel group with the same function name**: exact-id agreement 75.4%, corrected 89.8%. Do not
  quote the exact-id figure as the error rate.

  **K. pneumoniae HS11286 is absent from OrthoDB entirely**, so all 5,312 of its assignments are
  sequence-tier and only 1,923 clear the 0.9 bar — low confidence is not wrong, it means parallel
  groups competed. **Accession joins into OrthoDB do not work** (`xref_orthodb` is 0.0% on Kp), and
  **OrthoDB has NO root level spanning domains**: `<n>at2` and `<n>at2759` are **not comparable**, so
  cross-domain similarity is `neighbors.tsv`. Details: `docs/orthology.md`.

- **`ligands/chembl.py`** → **`evidence/chembl_<species>.tsv`**: measured bioactivity from
  **ChEMBL 37**, mapping
  the three bacterial proteomes onto ChEMBL's target sequences with DIAMOND and counting
  **non-redundant** ligands at four nested distances. A cached re-run needs no database.

  Proteins with a potent ligand: Kp 113 · Ec 96 · Sa 78. **That is ~2% of each proteome and must not
  be forced upward** — it is a fact about how little of the bacterial proteome anyone has screened.
  **Homology transfer is the whole game**: ChEMBL holds 21 Kp single-protein targets against 5,728.

  **Same-species matching is by organism NAME, not `tax_id`** — ChEMBL files strains under their own
  taxids, so `tax_id=562` finds 65 E. coli targets where the name finds **225**.

  **One `confidence_score` gate cannot serve both target types.** `≥8` returns **exactly zero protein
  complexes**, which would ship empty `complex_*` columns reading as a real biological zero — and
  **DNA gyrase is a `PROTEIN COMPLEX` in ChEMBL**, with GyrA/GyrB in `src/interest.py`. Complexes run
  at `≥6`; `PROTEIN FAMILY`/`PROTEIN COMPLEX GROUP` are **measured and excluded**, not silently
  absent.

  **Buckets are nested and restricted to true Bacteria** (`direct` ≥95% ⊆ `close` ≥60% ⊆ `remote`
  ≥40%). The restriction is v1's documented fix: an unrestricted "non-human" bucket with no identity
  floor gave **424** potent Kp proteins against a true 175. **Human is its own liability block, never
  merged into a bucket** — `clpP` carries 106 human compounds against 61 bacterial.

  **`pchembl_value` only exists for `=` relations on IC50/EC50/Ki/Kd/Potency in nM**, so MIC and
  %-inhibition are absent by construction: **this axis does not say "has an antibiotic"**, which is
  why the ribosome is largely missing — and degradability's top-100 is ribosome-heavy.

  **The version is ASSERTED, not declared, and the assertion was wrong first.** `version` IS NOT ONE
  ROW — it holds 11, `ChEMBL_37` is not first, and `LIKE 'ChEMBL_%'` also matches
  `ChEMBL_Structure_Pipeline`. The first `assert_version()` **rejected the correct database**: it had
  passed a synthetic one-row fixture, i.e. it tested the assumption rather than the schema.

  Needs **rdkit from PyPI** (a conda install flips `gradi` to osx-64 and takes ESM-C with it). The
  30.5 GB dump is deleted after each run; recovery in `data/raw/other/chembl/SOURCE.md`.
  **Load through `src/ligandability.py`** — **`load()` is the DELIVERABLE** (`ligands_<sp>.tsv`),
  `load_chembl()` the evidence table, `load_full()` the 19-column provenance view. Swapped
  2026-10-03: `load()` used to return the ChEMBL table, which stopped making sense once that was
  demoted. **`load_ligands()` is NOT any of these** — it is the raw `scratch/chembl_ligands.tsv`
  extract, 2.6M rows. Details: `docs/ligands.md`.

- **`ligands/effort.py`** → `scratch/chembl_effort.tsv`, **the DENOMINATOR and the axis's only real
  negatives**: `chembl.py` and `ligands.py` both require `pchembl_value IS NOT NULL`, so a
  compound somebody assayed that did NOT work is invisible, and "nobody screened this" collapses into
  "people screened it and nothing worked". Of **453 bacterial targets with ≥10 compounds assayed, 148
  (32.7%) never reached pChEMBL 6**, and 208 targets were invisible entirely.

  **TWO denominators ship and neither is merged.** An `IC50 > 100 µM` is the clearest statement in
  the database that a compound does not bind, so an effort count must include it — but a *ratio*
  needs its denominator from the numerator's own population. `hit_rate` is null, never 0, where
  nothing was assayed. It imports `assert_version` and `_activity_where` from `chembl.py` rather than
  restating them and **exits non-zero if either clause it relaxes is gone** — a denominator over a
  different population from the numerator is wrong in a way no shape check could see.

- **`ligands/validate_api.py`** — **the counts, checked against the LIVE ChEMBL API.** Everything
  else in this axis descends from three cached extracts, so a wrong SQL would be agreed with by
  every downstream check: the `chembl.py` control shares the extracts and the assertions share the
  code. This asks a different machine over HTTP. **88/88 comparisons over 49 proteins match
  exactly**, all three species, 1 to 12,438 compounds. Round 1 checks `n_ligands_own` against the
  resolved exact target (30 cases); round 2 checks `n_ligands_bacterial` and
  `n_measured_bacterial` as **UNIONS over every bacterial homolog** (29 proteins × 2; Kp `KPC-2`
  unions 52 targets) — that is where a double-count or a dropped `tid` would show.

  **It also checks the union against the per-target SUM**, which it must stay below: KPC-2 276 vs
  365, `bla` 269 vs 327, `folA` 443 vs 529 — and **the API agrees with the union, never the sum**,
  which is the distinct-molecule claim verified from outside.

  Two things that make the comparison legitimate: the activity endpoint does **not expose
  `confidence_score`**, which is fine HERE and nowhere else because `>= 8` removes 0 of 3,271,336
  single-protein rows, and only SINGLE PROTEIN tids are compared on both sides. A target too large
  to page is **skipped and named, never truncated** — Ec `ampC` (12,438 potent over 16 targets)
  is the one, and its exact count was verified in round 1 anyway.

  **It found a real defect**: `exact_target` round-tripped as `""` rather than NA, so
  `exact_target.notna()` was True for all 13,020 proteins and useless as a filter — the run picked
  up the whole proteome as test cases. Counts were unaffected; `_coerce_precedents` now maps empty
  to `pd.NA` and the column agrees with `exact_route != "none"` at 12 / 182 / 69.
  Evidence: `evidence/precedent_api_validation.tsv`. CLI: `--per-species 8` · `--round exact union`
  · `-q`. Needs network, no dump. ~15 min.

- **`ligands/transfer_calibration.py`** — **the bands 95/60/40 CANNOT be calibrated, and that is the
  result.** Over 1,582 ChEMBL target pairs where both ligand sets are known: compound-set overlap
  does not transfer at any identity (median Jaccard ~0.00 in every band — two near-identical targets
  are one enzyme screened twice against different libraries, so **Jaccard measures campaign
  coincidence**); `P(potent | neighbour potent)` is **FLAT** from 25% to 100%, because **62% of
  bacterial ChEMBL targets already carry a potent compound** — a protein enters ChEMBL when somebody
  believed it was druggable; and neither species nor RBH adds anything beyond identity, **so no
  orthology criterion is added**.

  **So the floor controls COVERAGE, not transfer reliability — document it as a conservatism choice,
  never as an accuracy threshold.** Calibrating it honestly needs proteins nobody chose to screen,
  which ChEMBL by construction does not contain.

- **`ligands/ligands.py`** + **`src/precedents.py`** — **THE AXIS DELIVERABLE**,
  `ligands_<species>.tsv`, **10 columns**; also a query tool for ANY sequence, needing no database
  (three cached extracts, 82 MB, ~1 s).

  **`chembl_<species>.tsv` was DEMOTED to `evidence/` on 2026-10-03** (owner's call), because
  measured against `ligands_<species>.tsv` it is largely a duplicate: `n_ligands_bacterial` ==
  `remote_n_compounds` and `n_ligands_human` == `human_n_compounds` **exactly** — same 113 Kp
  proteins, 5,286 vs 5,286 compounds, ρ 1.0. **What ONLY the chembl table has, so read it from
  `evidence/` when you need it**: `*_n_scaffolds` (Kp's 5,286 potent compounds are **1,593 Murcko
  scaffolds**, 3.3 per scaffold — the difference between 50 starting points and 50 analogues of one
  series), the identity bands broken out (direct 21 / close 79 / remote 113 proteins), the match
  provenance (`best_target`/`best_pident`/`best_organism`) and `allorg_*`. **What ONLY the
  deliverable has is the denominator**, `n_assayed_*` — the axis's only real negatives.

  **The standard pair, added 2026-10-04, divides the axis's own two questions.**
  **`ligands_consensus`** is the OUTCOME: mean within-species percentile over `n_ligands_own` ·
  `n_ligands_bacterial` · `best_pactivity_bacterial`. **`n_ligands_human` is EXCLUDED** — a
  liability pointing the other way — and so is `n_assayed_*`, which is effort, not ligandability.
  **THE ZERO BLOCK IS PINNED TO 0, a deliberate deviation from `percentile_consensus()`**: ranking
  runs only within the ~3% with any evidence (Kp 180 · Ec 160 · Sa 119), because under average-rank
  ties the other 97% would read **0.491**, mid-scale, and "nobody looked" would appear moderately
  ligandable. For scale, `essentiality_consensus`'s largest tie block is **8**.
  **`percentile_consensus()` itself is unchanged.**

  **`ligands_evidence` grades PROVENANCE, NOT OUTCOME** — **3** somebody assayed THIS protein ·
  **2** only a bacterial homolog was assayed · **1** nothing in ChEMBL, so the 0 beside it is an
  open question, never a measured negative. Kp **5,453/264/11** · Ec **4,151/79/173** · Sa
  **2,726/94/69**. The split keeps the two columns independent: overall ρ 0.79–0.85 but **only
  0.11–0.24 within the evidence-bearing subset**, and **62 Ec proteins sit at evidence 3 with
  consensus 0** — assayed directly, nothing potent, which no single column can say.
  **KP'S 11 IS THE FINDING**: `bla`, `KPC-2`, `blaSHV-11`, `blaCTX-M-14`, `ybtE`, `rfbD`, `rpsR`,
  `atsA`, `dxs`, `uppS`, `acpP` — **the only K. pneumoniae proteins anyone has screened directly
  are the resistance enzymes.** **It is NOT the dropped four-tier `precedent_evidence`**, which was
  an identity-band vocabulary; every `<axis>_evidence` is derivable from its own axis, and the
  point is cross-axis comparability. Both columns match **no rule in `_read()`** and are registered
  by name in `_coerce_precedents()` plus a new `DELIVERABLE_DTYPES` assertion — the trap this axis
  first paid for as `best_pchembl`.

  **A RENAME IN GIT DOES NOT REACH A FILE IN eosvc, and this axis proved it.** Commit `e881e47`
  renamed six deliverable columns on 2026-10-03; `data/` is gitignored, so
  `evidence/precedents_full_<sp>.tsv` kept `n_ligands`, `n_assayed`, `best_pactivity_bacteria`,
  `n_targets_bacteria`, `best_pident_bacteria`, `n_measured` — and **`ligands/validate_api.py` was
  silently broken for a day**, reading `n_ligands_own` off a file that had no such column. The
  2026-10-04 regeneration repaired it. **After renaming a column, regenerate every file that
  carries it**, and prefer re-running a stage to rewriting a deliverable in place.

  **TWO QUESTIONS, NOT ONE, and a single count conflates them**: `n_ligands_*` is POTENT (pChEMBL ≥ 6)
  and `n_assayed_*` is "has anyone looked", each over three scopes — this protein, the bacterial pool,
  and human (**a LIABILITY, never summed into the bacterial count**). **A 0 against 158 assayed
  compounds is a measured discouragement; a 0 against 0 is an open question.** `n_assayed_*` is **NA,
  never 0**, when the effort extract is missing, because a 0 would claim nobody ever assayed the
  protein — the opposite evidence from "we do not know". (The old four-tier `precedent_evidence`
  column was dropped as redundant once the denominator shipped: `no_homolog` is just
  `n_targets_bacterial == 0`.)

  **"UNIQUE" is the whole point and it is not a per-target sum** — counting DISTINCT
  `parent_molregno` over the union of homologous targets counts MOLECULES, where summing per-target
  counts inflates, worse the wider the band. **The complex track is reported separately and is NOT
  inside the bacterial count**: Kp gyrA carries 1,410 complex ligands against 131 single, and v1
  dropped that track and made GyrA/GyrB look unliganded.

  **`exact` IS SPECIES-LEVEL, not byte-level — a protein does not stop being itself between
  strains**; it is the union of accession, identical sequence and same species at ≥95%, and 95 is a
  statement about protein IDENTITY, **not** an accuracy threshold. Kp `pyrH` is the case to remember
  — 158 compounds assayed against a 98.3%-identical target and not one potent, previously
  indistinguishable from "nobody looked".

  **Identity alone cannot transfer a ligand count**: Kp `A0A0H3GWM6` is 99.2% identical to E. coli
  `P0ADG7` and correctly gets nothing, because the hit is a 130-aa fragment against a 488-aa protein
  and `MIN_QCOV`/`MIN_SCOV = 50` rejected it.

- **`ligands/bindingdb.py` is a MEASUREMENT, not a deliverable.** v1 shipped BindingDB as a co-equal
  track **without ever measuring the overlap**; measured, it adds 33 proteins across all three
  species, 11.5% over ChEMBL's 287. **Not promoted**; if it ever is, columns go *beside* the ChEMBL
  ones, never merged.

- **`pockets/structures.py`** → **`esmfold.py`** → **`predict.py`** → **`pdb_coverage.py`** →
  **`alphafill.py`** → **`holo.py`** → **`merge.py`** — **structural ligandability**: can a small
  molecule bind this fold? Deliverable `pockets_<species>.tsv`, **7 columns**, complete and
  canonical for the three bacteria: `p2rank_score` · `fpocket_score` (predicted, on AlphaFold v6
  models) · `n_ligands_pdb` (**measured**: drug-like ligands in this protein's OWN PDB structures)
  · `n_ligands_alphafill` (**modelled**: drug-like ligands AlphaFill transplanted onto its model)
  · `n_pdb_structures` (PDB entries that ARE this protein, ligand or not; partial structures
  count; coverage is in `evidence/pdb_<sp>.tsv`) · `af_plddt`.

  **NEVER SUM THE TWO LIGAND COUNTS** — a co-crystal of this protein (88 Kp / 308 Ec / 90 Sa
  proteins) and a transplant from a ~30%-identity homolog (1,533 / 1,196 / 704) are different
  evidence. Both are **non-redundant = distinct Bemis-Murcko generic scaffolds**, reusing
  `scripts/ligands/chembl.py:_scaffold_chunk` rather than a second copy; raw code counts overstate
  by 25–34% and ship as `n_codes_*`.

  **No `evidence` column** (owner's call, 2026-10-03): it was a function of two shipped columns —
  `af_plddt` is NA exactly when no model exists, `n_ligands_pdb > 0` exactly when a drug-like
  ligand was seen on this protein — verified exactly reconstructible on all three species before
  removal. **An NA in the pocket columns is "could not look", NOT "looked and found nothing"**: a
  protein WITH a model and no admitted pocket gets 0, so never `fillna(0)` — that is the v1
  mistake. The one thing the labels added is the model's provenance (AlphaFold DB vs ESMFold),
  which is `model_source` in `evidence/alphafold_<species>.tsv`. All `gradi`; fpocket/P2Rank from `gradi-pockets`, DIAMOND from `gradi-ortho`. ~40 min
  cold. **Load through `src/pockets.py`.** Details: `docs/pockets.md`.

  **Confidence enters ONCE**: a pocket counts only if its residues average pLDDT ≥ 70. v1 applied
  pLDDT twice. Neither tool reads it (P2Rank's `alphafold` config drops B-factor), so the filter is
  the only use. **Reproduces v1 exactly**: any P2Rank pocket on 4,542 Kp / 3,589 Ec proteins.
  AlphaFold models are used only if their sequence equals the proteome's (0 mismatches).

  **34 proteins have no AlphaFold DB model, and no other accession has one either** (checked via
  UniParc identical-sequence groups): AFDB skips < 16 aa (21 Ec micro-peptides), selenocysteine
  (fdhF/fdnG/fdoG), pseudogenes written with X, and > 2,700 aa. `esmfold.py` folds the 32 up to
  2,700 aa with ESMFold v1 (HF port, in `gradi`, CPU; U→C, X kept — **ESMFold writes no atoms for
  an X but keeps the numbering, so read models by residue number**). The two giants (Kp irp1, Sa
  ebh) stay NA by the **owner's decision** (~18 h of CPU). `model_source` records the predictor.

  **PDB coverage is by SEQUENCE**: DIAMOND vs every `pdb_seqres` protein chain, ≥ 95% identity,
  ≥ 50% of the chain aligned → **Kp 569 (9.9%)** proteins with a structure, against **30 (0.5%)**
  by v1's SIFTS accession route; Ec 1,893 (43.0%), Sa 595 (20.6%). ≥ 95% counts near-identical
  proteins of other species (Kp rpoB inherits E. coli's 411 RNAP entries), and coverage is of
  SEQRES, not resolved residues.

  **"Drug-like" is built from published sources, not a denylist** — BioLiP, PLINDER's artefact
  list, PDBe cofactor classes, a nucleotide SMARTS, and **ECMDB metabolites (owner's choice,
  2026-10-03)**; each is a flag in `evidence/ligand_classes.tsv`. **QED ≥ 0.2 and Ro3 "fragment"
  were measured and REJECTED: both delete antibiotics** (novobiocin 0.184, rifampicin 0.109, the
  aminoglycosides; fosfomycin, D-cycloserine). The vocabulary covers AlphaFill's codes too, with
  `in_biolip` as a flag rather than a requirement. Two traps: an `[R1]` ring SMARTS misses
  cyclic-di-GMP (use `[R]`), and BioLiP alone leaks detergents/cryoprotectants (hence PLINDER).
  A third, worth remembering project-wide: **`NA` is sodium's chemical-component code**, and
  pandas reads it as a missing value — it silently nulled 12% of the transplant rows once.

  **A FUSION CONSTRUCT CARRIES ITS PARTNER'S LIGANDS.** A chain counts as this protein at ≥ 95%
  identity over ≥ 50% **of the chain**, which also admits a fusion where this protein is the
  larger half — Ec `malE` (MBP, a crystallisation chaperone) reads 41 ligands that are really the
  fusion partners' (its chains sit at median `chain_coverage` 0.62 against fabI's 0.989). Only
  3.5–5.0% of rows come from chains below 0.9, and a 0.9 floor fixes malE (41 → 2) but strips
  Sa `gyrA` (16 → 4), whose fluoroquinolone complexes are genuine GyrB–GyrA fusions — so the floor
  is **measured and NOT applied**; `chain_coverage` ships per chain instead.

  **DO NOT TRANSFER LIGANDS BY SEQUENCE — that is reinventing AlphaFill, and worse.** Until
  2026-10-03 this axis shipped `holo_identity` (DIAMOND vs BioLiP holo chains, a binding-site span
  test, the best hit's % identity). Measured on the same proteomes with the same drug-likeness
  rule: **AlphaFill reaches 1,533 Kp / 1,196 Ec proteins, that route reached 184 / 172.** It was
  also bimodal rather than continuous, and never checked the site was conserved. Deleted with it:
  `alphafill_check.py`, whose "AlphaFill adds little" verdict came from imposing a 40% identity
  floor **AlphaFill does not use** (its donors sit at ~30% median identity by design). AlphaFill
  transplants are **not filtered on `local_rmsd`** — it publishes the metric, not a threshold —
  and the ligand is **`analogue_id`, not `compound_id`** (ANP→ATP, ACO→CoA; 2.7% of transplants).

  **Control for LENGTH before quoting any pocket-vs-PDB agreement** — length alone predicts a
  measured ligand at AUROC 0.65–0.67 (big proteins are crystallised more and have more surface).
  Within length deciles, against `n_ligands_pdb > 0`: **P2Rank 0.494 (Kp) / 0.561 (Ec) / 0.619
  (Sa), fpocket 0.435 / 0.523 / 0.477** — so **on the anchor the pocket scores add nothing over
  protein size**, and they are a soft prior at best. `merge.py` reports `*_len`; quote that, never
  the raw AUROC. **No pocket count** (removed 2026-10-03): the P2Rank authors set no probability
  cutoff, and without one the count is size (ρ 0.84 with length, no signal within length); a
  0–0.7 sweep never beat `p2rank_score`. **No `druggability()` helper either** — no defensible
  weighting exists across a weak prior, a sparse measurement and a third party's model.

- **`essentiality/labels.py`** + **`essentiality/deg_proteomes.py`** — the training corpus, from
  **DEG**: 49 of 51 datasets, 173,048 labeled proteins, 20,194 essential (11.67%), 38 species, **each
  gene with its protein sequence**, because the join downstream is **by sequence** (DEG's UniProt AC
  is 45.7%, locus tag 42.9%).

  **Exclusions are FLAGGED, never dropped**: 4 non-genome-wide methods, where "absent from the list"
  cannot mean non-essential, and 11 condition-specific screens — which is why *P. aeruginosa* PAO1
  swings **117 → 336 → 551** essential genes across its three datasets.

  **Exact sequence matching is NOT enough, and the per-dataset rate is the control**: DEG's vintage
  annotation has drifted, so exact match recovers only 59–95% and `locus_tag` is 0% on the older
  datasets, while **DIAMOND ≥95% lifts the median to 99.5%**. Without it, hundreds of genuinely
  essential genes are silently labeled non-essential.

  **Two traps worth carrying here.** **326 of the 26,619 "sequences" are the literal string `Not
  available now.`** — the join succeeds, so a coverage check reads 100% while the payload is junk;
  only validating the amino-acid alphabet catches it. And **do not scan for accessions with a
  `\b`-anchored regex**: it silently drops `NZ_`-prefixed accessions because `_` is a word character,
  which cost E. coli O157:H7, the largest single dataset, on the first run.

  **Watch the base rate before training** — it spans 15× across species. **OGEE v3's server is gone
  but the bulk file is RECOVERABLE** from two byte-identical mirrors (the Wayback URL's `id_` suffix
  is mandatory): **a dead server is not a lost dataset.** It widens the corpus but does NOT close the
  anchor's hole — **K. pneumoniae is absent from OGEE as from DEG**. For any gated supplement, Europe
  PMC's `supplementaryFiles` endpoint **bypasses the Cloudflare 403s** on ASM, PNAS, Nature and Cell.

- **`essentiality/geptop.py`** → `geptop_<species>.tsv`. Geptop 2.0 **ported faithfully to Python 3**
  (three genuine defects fixed and **nothing else**; v1 substituted median RBH identity for the
  composition-vector distance, a real deviation not repeated). **Accuracy is 0.59–0.81 over all
  proteins (0.61–0.84 over scored ones), and the paper quotes a mean of 0.84 — do not quote 0.84 as
  ours.** *H. influenzae* manages only 0.587, and **the Kp prediction inherits that variance.**

  **A score of 0 means two different things — read `geptop_evidence`, never the score alone.** On Kp,
  58.5% is `orthologs_none_essential`, **a confident NON-essential call, evidence not absence**, and
  only 7.8% is `no_orthologs`: coverage is **92.2%**, not the 33.7% an earlier version implied. Zero
  rows are **tied** — unranked, not low-ranked.

  **Two traps.** `geptop_score` is **proteome-relative**, so use `geptop_score_raw` across species;
  and **two anchors ARE references**, contributing 53.2% of E. coli's own reference weight and 58.3%
  of S. aureus's against 3.6% for Kp. Leave-one-out is deliberately **not** applied: DEG supplies a
  measured label for those two, and Kp is absent from all 37. **`blastp` borrowed from
  `gradi-prokka`** (`GRADI_BLAST_BIN`); **human excluded by construction**.

- **`essentiality/ogee.py`** + **`ogee_proteomes.py`** + **`ogee_dataset.py`** → `ogee_<species>.tsv`:
  `ogee_ess` (0–1, never null) + `ogee_evidence` (the MEASURED OGEE label, or empty).

  **HALF OF OGEE IS UNUSABLE**: of 87 taxa with decided calls, **40 have ZERO negatives**, 25 of them
  from the Fitness Browser RB-TnSeq collection, which cannot see essential genes by construction. The
  project's filters leave 26 taxa / 78,893 proteins / base 0.173.

  **Leave-species-out, not a single holdout** (owner's instruction, and the better fit since Kp is
  absent from OGEE): 26 estimates instead of 1, and **the spread is the result** — AUROC 0.529–0.940.
  **THE HEADLINE: `corr(base_rate, AUROC) = −0.643`** — screens calling many genes essential are much
  harder to predict, and *P. aeruginosa* PAO1 (base 0.081) scores 0.837 while PA14 (0.302) scores
  0.634. A screen calling 40% of genes essential is measuring fitness defect, not essentiality.
  **This prices the Kp column**: our Kp screens run base 0.075–0.106, where held-out taxa score
  0.78–0.93.

  **ONE model for all three anchors** (owner's instruction), with the caveat carried as DATA per
  protein: E. coli K-12 and Sa NCTC 8325 are OGEE taxa *and* our anchors, so where `ogee_evidence` is
  non-null `ogee_ess` is closer to recall than prediction. Same model, top-decile cut **0.861 on
  E. coli against 0.471 on Kp** — so **`ogee_ess` is comparable WITHIN a species, never across**.
  **`--seeds 1` is deliberate here** against the axis-wide 5: the folds are FIXED, one taxon each, so
  there is no partition randomness to average.

- **`essentiality/proteomelm_ess.py`** → `proteomelm_ess_<species>.tsv`: the **paper's own
  essentiality head**, deferred for months as "weights unreleased" and released to us by Cyril
  Malbranke on 2026-10-01. `proteomelm_ess` is p(essential), 0–1, never null; `proteomelm_ess_rank`
  is within-proteome; `proteomelm_ess_evidence` is the column to read first. **Apache-2.0** — a far
  easier licence than TabPFN's, so this one can ship.

  **THE SCORE MEANS SOMETHING DIFFERENT IN EACH ANCHOR**, from the authors' own `genomes.tsv`
  (staged at `data/source/proteomelm/ess_genomes.tsv`): **ecoli `held_out`** (their Fig. 5B genome,
  290 E / 3,969 NE), **saureus `in_training`** (taxid 93061 — our exact 2,889-protein proteome is in
  their cross-validation set, so that column is closer to RECALL), **kpneumoniae `unseen_species`**
  (no *K. pneumoniae* in their 89 genomes; the only *Klebsiella* is *K. michiganensis*,
  positives-only). Comparable WITHIN a species, never across — the rule `ogee_ess` already carries.

  **It must NOT be fed our own `proteomelm_<species>.npz`.** Same backbone, same layer 8, but ours
  are **z-scored genome-wide** and **window** sequences above 4,096 aa, while the head takes **raw**
  `hidden_states[8]` with sequences **truncated** at 4,096. The worker recomputes ESM-C and the
  backbone pass end to end. (The difference that turned out NOT to exist: their `group_embeds=x` is
  our `self` mode, because the model does `if group_embeds is None: group_embeds =
  inputs_embeds.clone()` — checked rather than assumed.)

  **Env `gradi-plm-ess`**, across a process boundary (`essentiality/workers/proteomelm_ess.py`,
  `GRADI_PLM_ESS_BIN`). The head ships only in the authors' git build, which pulls **torch 2.14**
  against `gradi`'s 2.12 — the collision that breaks stage 01's ESM-C. ~7–12 min per proteome on MPS.

  **CLASS 0 IS ESSENTIAL** (`id2label: {0: essential}`), so an off-by-one yields a confident,
  well-formed, exactly inverted column. The stage runs a ribosome-vs-dispensables polarity control
  and **exits non-zero** below 0.80: measured **Kp 0.9370 · Ec 0.9717 · Sa 0.9918**.

  **Measured on labels the authors never saw** (`evidence/proteomelm_ess_validation.tsv`). **We
  reproduce their headline**: they report 0.952 held-out on E. coli, we measure **0.9726** against
  Keio on our exact anchor — which is what validates the whole chain. Kp, scoring each screen
  **strain itself** so no identifier mapping is involved (100% key overlap on all three):
  **ATCC 43816 0.9473 / AUPR 0.748 · ECL8 0.8323 / 0.643 · RH201207 0.8201 / 0.597**.

  **But our own pipeline still wins on Kp**, and that is the result to keep: assay-matched
  Goodall→Kp reaches **0.9597 / 0.8845 / 0.8903** on the identical three endpoints, beating this
  head on **3 of 3, on both AUROC and AUPR**. ProteomeLM-Ess beats only the assay-MISmatched Keio
  transfer, and only on RH201207 (0.820 vs 0.810) — the same "transfer is better when the ASSAY
  matches" finding again. So it ships as an independent fifth opinion, not as a replacement.

  **It is genuinely independent**, which is why it is worth a column: ρ 0.44 with `geptop_ess`, 0.37
  with `ogee_ess`, 0.32 with `screens_ess_mean` on Kp, and top-500 shortlists overlap only ~335/500.
  Contrast degradability's two activator columns at ρ 0.89, which are one opinion wearing two hats.

  **The base-rate effect reproduces here, from outside this project**: across eight screens in two
  species the AUROC falls near-monotonically as the screen's base rate rises (Ec 0.048 → 0.980 down
  to 0.103 → 0.645; Kp 0.076 → 0.947 down to 0.106 → 0.832). That is independent confirmation of the
  `corr(base_rate, AUROC) = −0.643` this axis measured on OGEE, from a different model on different
  labels. **Their training set has the same flaw we documented**: 37 of their 82 cross-validation
  genomes have ZERO negatives — the positives-only RB-TnSeq artifact.

  **A finding for the collaboration**: the consortium's own panel sits at the **97.1st percentile
  (median)** of this score on Kp — `lpxL` 99.9, `lptG` 99.8, GyrA/GyrB 99.4–99.8, `lnt`, `lolC`,
  `secA`, `yidC`, `lptD`, `lptB` all above 98.5. Beside studiedness (80th percentile — not novel)
  and degradability (OR 0.21–0.82, trending depleted), the panel is **the right biology and the
  wrong chemistry for a degrader**. Details: `docs/essentiality.md`.

- **`essentiality/screens.py`** + **`predict.py`** + **`summary.py`** — published screens as
  independent endpoints, **TEN training sets, one per SOURCE**, conditions aggregated within a source
  (owner's rule). Three filters define the set, each on instruction: *both classes required*, *no
  condition-dependent data for now*, *no duplicates*; everything excluded stays parseable in
  `HELD_BACK` — **decisions, not gaps**. **Nothing is merged**: base rates span 0.048–0.194, and that
  spread is method and strain, not noise.

  **`summary.py` IS THE CORRECTNESS CHECK.** `num_positives` looks identical whether a label set is
  right or inverted, and an inverted column trains a confident, well-formed, exactly wrong model, so
  every column is scored against the **ribosome** (must be essential) and the **textbook
  dispensables** (must not be). It caught a real one: a compendium DESeq table **did not converge**,
  and of the 189 genes with **zero insertion sites, 108 were called `Unchanged`** — DESeq cannot
  compute a statistic for a gene with no insertions and falls through to non-essential, inverting the
  label for exactly the genes that matter most.

  **CALIBRATE A BIOLOGICAL CONTROL ON A MEASUREMENT, NOT ON 1.00 — got wrong once.** A first pass set
  the bar at 0.50 separation and would have failed the gold standard: **Keio reaches only 0.774**,
  and its misses are genuinely dispensable in *E. coli*.

  **TEN training sets, NINE screen endpoints.** `ogee_corpus` is the tenth and feeds `ogee_ess`
  instead, so `screens_ess_mean` is the mean of **9** model probabilities, not 10 — and they are **not 9
  independent votes**: all nine read the same ProtT5 embedding, so a high value is one correlated
  opinion, not a consensus count.

  **Features: ProtT5, measured not assumed** — on Keio, paired over identical folds, ProtT5 −
  ProteomeLM is **+0.0330 PR [5/5 seeds] but −0.0002 AUROC [2/5, a tie]**: AUROC alone would have
  called the winner a coin flip.

  **MMSEQS2 REWRITES FASTA HEADERS CONTAINING `|`** — DEG's `>lcl|…` comes back with `lcl|` silently
  gone, which is exactly the key the screen table and the embedding share. The clustering then mapped
  to **0 of 4,981** proteins, indistinguishable from a proteome with no paralogs;
  `paralog_clusters.py` now writes index surrogates and substitutes the real ids back.

  **THE HEADLINE: an E. coli-trained model reaches K. pneumoniae with no orthology anywhere** —
  `--score-on` fits on one organism and scores on another's OWN labels and OWN embeddings, giving
  **Ec → Kp mean AUROC 0.8597** and **Kp → Ec 0.9338**, better when the **ASSAY matches**.
  **Quote BOTH the grouped CV and the cross-species number, never one alone**: grouped CV holds out
  whole paralog FAMILIES, while cross-species trains on every family and tests on an organism where
  most genes have a homolog. The second is not label leakage — nothing is copied — but gene-family
  overlap is a real information channel. The first is the honest generalisation measure; the second
  is the operationally relevant one for Kp. Catalogue: `docs/essentiality_screens.md`.

- **`essentiality/strain_homologs.py`** → **`merge.py`** + **`registry.py`** →
  `essentiality_<species>.tsv`, **6 columns**: `uniprot_ac` + three predictors (`geptop_ess`,
  `proteomelm_ess`, `screens_ess_mean`) + **`essentiality_consensus`** (0–1) and
  **`essentiality_evidence`** (1–3) — plus `deg_<species>.tsv` for the measurement. Joined **by
  sequence, not accession**.

  **THE TWO STANDARD COLUMNS EVERY AXIS WILL SHIP** — convention in **`src/consensus.py`**, which
  the other nine axes reuse. `<axis>_consensus` is the mean **within-species percentile rank** of
  that axis's predictors (ranked first because the inputs are not comparable as values; **never
  compare it across species**). `<axis>_evidence` is **count + concordance**: **3** = ≥2
  independent experimental sources, unanimous, AND agreeing with the consensus · **2** = one
  source, or several that conflict · **1** = no measurement. **Level 0 cannot occur** — every
  protein has a prediction.

  **`<axis>_evidence` IS NOT PURELY EXPERIMENTAL** (owner's call, 2026-10-04): a protein measured
  twice, unanimously, whose consensus contradicts it lands at **2**. Do not read 2 as "the
  experiment was weak". Concordance uses a **base-rate cut, never 0.5** — essentials are 11–17% of
  a proteome — recorded per run in `evidence/consensus_audit.tsv`.

  **`strain_homologs.py` is why Kp has any evidence at all.** DEG has **no Klebsiella** and exact
  sequence matching against ALL of DEG reaches **29 of 5,728 Kp proteins (0.5%)**. A DIAMOND join
  at **≥95% identity / ≥50% subject coverage** to the three screened Kp strains reaches **4,725
  (82.5%), 4,390 in 2+ strains**. **It is a FLAG, NOT A LABEL TRANSFER** — the measured call never
  enters the consensus, so the rejected label-transfer stays rejected. **Kp reaches level 3 for
  3,973 proteins and NOT ONE is measured on HS11286**; `evidence/strain_homologs_<species>.tsv`
  records which strains covered each, so that stays checkable.

  **`essentiality`/`essentiality_source` were dropped** (owner's call, 2026-10-03), so **which
  column to rank on is now an explicit choice**. Two reasons, both measured: the merge **MIXED
  UNITS** — a measured call pinned to 1.0/0.0 against a continuous prediction, so every measured
  essential outranked every prediction by construction — and on **Kp, the anchor, it was a verbatim
  copy of `screens_ess_mean` for all 5,728 rows**. Rebuild it with
  `np.where(deg_essential_any.notna(), deg_essential_any, screens_ess_mean)`.

  **`geptop_evidence` and `geptop_in_reference_set` moved to `geptop_<species>.tsv`**, where they
  were byte-identical duplicates. **This matters**: a `geptop_ess` of 0 still has TWO meanings —
  `orthologs_none_essential` (58.5% of Kp, a confident NON-essential call) vs `no_orthologs` (7.8%)
  — so read `geptop_evidence` from the per-source file before reading a zero. **Every dropped column
  was verified byte-identical or exactly derivable before the tables were rewritten.**

  **`ogee_ess` is NOT in the summary** (owner's call, 2026-10-03); `ogee_<species>.tsv` and the OGEE
  scripts are untouched. **The obvious reason is WRONG and was measured**: it does not duplicate
  `proteomelm_ess` — those two run **rho 0.33–0.64** sharing only **296–384 of their top 500**, i.e.
  genuinely different opinions. It is redundant with **`screens_ess_mean`** (rho **0.488 / 0.550 /
  0.365** against ProteomeLM's 0.321 / 0.316 / 0.127), which both drives the merge and validates
  better — **AUROC 0.89–0.96** on the measured Kp screens against OGEE's leave-species-out
  **0.529–0.940** and a Kp top-decile cut of 0.471 where Ec reads 0.861.

  **`deg_ess` and `essentiality_rule` are NOT in the summary either** (same call). `deg_ess`
  was a **byte-identical duplicate** of the column in `deg_<species>.tsv` — checked per species
  before the rewrite — and that file is RICHER, carrying `deg_essential_any` and `deg_essential_all`
  side by side where the summary held only whichever `--rule` picked. **Nothing about `essentiality`
  changed**: the measurement is still inside it wherever `essentiality_source == "measured"` (Ec
  4,253 · Sa 2,678). `essentiality_rule` carried **nothing** — a 1:1 function of
  `essentiality_source`; the rule is a property of the RUN and is in
  `evidence/essentiality_merge_manifest.tsv` beside BOTH counts, which is what keeps the 3.4× spread
  visible.

  **Two screens of the same strain disagree, and the column records it rather than hiding it**: on
  E. coli 490 of 695 essential calls rest on ONE screen, so `--rule any` gives 695 and `--rule all`
  gives 205 — a **3.4× spread from one choice**. **`geptop_ess` and `deg_ess` are comparable for
  RANKING, not as VALUES** — a predicted 1.0 means ~70% where a measured 1.0 means it *is* essential,
  so the convenience column `essentiality` **mixes units**. **Rank within a species, never across.**

  **`essentiality`'s FALLBACK is `screens_ess_mean`, not `geptop_ess` — changed 2026-10-03.** Measured
  on K. pneumoniae, the only anchor where every predictor is honest: the screens-trained transfer
  models reach **AUROC 0.89–0.96** on the three measured Kp screens against Geptop's validated
  0.59–0.81, and `geptop_ess` left **3,799 of 5,728 Kp proteins (66.3%) tied at exactly 0** — i.e.
  two-thirds of the anchor unranked in the headline column, on an axis consumed by ranking. After
  the change, Kp ties at 0 are **0**. `essentiality_source` now reads `predicted_screens_ess_mean`.
  **Do NOT re-derive this from an E. coli comparison**: `geptop_ess` is **53.2% self-derived on
  E. coli and 58.3% on S. aureus** (both are Geptop reference genomes), so it looks like the best
  predictor there — AUPR 0.977 against Keio — and that number is circular. Only the fallback
  changed; `geptop_ess` keeps its column and its file.

  **This axis ships TWO tables at the task root, deliberately** — the second, `geptop_<species>.tsv`,
  reads like evidence but the project owner chose it as a deliverable; do not demote it in a tidy-up.

  **The axis has FOUR tiers**, restructured 2026-09-21 on the owner's instruction:

  ```
  data/{raw,source}/**/SOURCE.md          34 dirs: what each dataset IS and whether v2 uses it
  data/processed/essentiality/
    dataset_registry.tsv                  47 rows: every dataset found + its DISPOSITION
    training_sets/                        10 files: the clean ML-ready sets, uniform schema
    <source>_<species>.tsv                one table per EVIDENCE SOURCE, per species
    essentiality_<species>.tsv            the headline summary, one column per source
    evidence/ · scratch/                  audits · caches
  ```

  **`dataset_registry.tsv` is GENERATED by `registry.py`** from `screens.py`'s own declarations — a
  hand-maintained second copy is how a register stops agreeing with the code — and it **exits
  non-zero unless it reconciles with `training_sets/` in both directions**. **Source directories are
  DECLARED, never inferred**: an earlier version guessed by substring and matched a Goodall screen to
  `data/source/go`, because "go" is inside "goodall", and a false provenance link reads exactly like
  a correct one.

  **`training_sets/*.tsv` share one schema: `key · label · source_id · features_from`.** The key
  column used to be `uniprot_ac` and that was **actively misleading** — for five of the nine screens
  it holds a RefSeq id, an EMBL id, a locus tag or a DEG FASTA header, which key onto their species'
  anchor proteome at **0 of ~5,000**: a join on accession returns an empty frame, not an error.

  **`screens_<sp>.tsv` IS PREDICTIONS THROUGHOUT.** Transferring the measured labels was considered
  and rejected on the owner's instruction — exact-sequence transfer recovers only 13.7–71.6% of rows,
  so it would silently mislabel ~1,700 measured essentials as non-essential. Where an anchor protein
  WAS in a screen's training set the **out-of-fold** value is substituted, so the column is not part
  in-sample and part honest; `evidence/screens_transfer_audit.tsv` must be read before any `_prob` is
  treated as evidence.

  **`essentiality/fetch_screens.py`** downloads published screens and **writes down what a human must
  fetch by hand** — the real output is `evidence/screen_fetch_status.tsv` plus a `PLACEHOLDER.md`
  naming what to click. **Publishers and repositories want OPPOSITE headers**: publishers 403 a bare
  client, while **figshare returns HTTP 202 forever to a full Chrome User-Agent and 200 to a short
  one**. **Three payloads returned HTTP 200 and were not data**, all caught by content inspection.
  Details: `docs/essentiality.md`.

- **`studiedness/fetch.py`** + **`gene2pubmed.py`** + **`pubtator.py`** + **`unknome.py`** +
  **`transfer.py`** + **`merge.py`** → `studiedness_<species>.tsv`: **three counts** — `n_papers_uniprot_own` ·
  **`n_papers_uniprot_prokaryotic`** (the ranking) · `n_papers_pubtator_prokaryotic`.
  **The names say SOURCE and DONOR SCOPE and claim nothing more — "family" was dropped
  2026-10-03 because neither `_prokaryotic` column aggregates a family**: each reads ONE donor,
  and not the same one (the two donors agree on 84.5% of Kp, differ on 15.5%). **`_prokaryotic`
  names the DONOR POOL, not a species search** — for 99.8% of Kp the donor is another organism,
  usually E. coli K-12.

  **REJECTED, built and removed the same day: `n_papers_pubtator_own`** — 98% empty on Kp, 85% on
  Sa; **a zero-filled column reads as noise, not as evidence.** The free-text rebuild
  (`<gene_name> AND @SPECIES_<taxid>`) was measured at **22/24 non-zero on a Kp sample against
  2/24** and rejected anyway: co-occurrence not identity, caps at gene-name coverage (Kp 63.4%),
  ~8,000 API calls, and it would imply the same definition as an exact curated count.
  **`n_papers_pubtator_prokaryotic` STAYS** — it is 41% empty against the curated column's 34%,
  **84% of that emptiness is shared** (proteins with no donor at all), and it is the **best
  predictor on the held-out control, 0.3722 against 0.3398**. Do not confuse the two.

  **THREE COUNTS, THREE DEFINITIONS, NEVER a `max()` ACROSS THEM** — that is what killed the 0-1
  composite on 2026-09-22. **Rank on `n_papers_uniprot_prokaryotic`.**

  **`evidence` IS NOT IN THE DELIVERABLE** (owner, 2026-10-03) — it stays in
  `evidence/transfer_<sp>.tsv`. So **a 0 in `n_papers_uniprot_prokaryotic` is ambiguous in this table**:
  `no_hit` (nothing among 575,748 curated entries resembles it — the strongest novelty claim the
  axis makes) and `below_floor` both read 0. On Kp that is 1,961 proteins. Join `load_transfer()`
  before reading a 0 as novelty.

  **`n_papers_uniprot_own` USES THE LIGANDS AXIS'S `exact` RULE** (2026-10-03): the UNION of PMIDs over
  **accession + identical sequence + same species ≥95%**, because *a protein does not stop being
  itself between strains* and the two axes must mean the same thing by "this protein".
  **Union, not max** — the analogue of pooling distinct molecules. **S. aureus proteins with any
  literature went 357 → 1,049** (NCTC 8325 is the anchor while most Sa curation sits under
  Newman/USA300/Mu50/N315); Ec median 5 → 8; Kp distinct values 7 → 32. Held-out control
  **0.3280 → 0.3398**. **`own` can now EXCEED `family`** (Ec median 8 vs 6) because `own` unions
  every strain entry while `family` reads ONE donor — the old "+3 on Kp, +0 on Ec" framing is
  withdrawn. **Still do not rank Kp or Sa on `_own`; rank on `n_papers_uniprot_prokaryotic`.**

  **`pubtator.py` → `n_papers_pubtator_prokaryotic`, keyed on NCBI GeneID.** A GeneID names one gene in
  one organism, so **the species is already in the key** and there is no symbol ambiguity: Kp
  `crp` reads **376** where the gene-SYMBOL route read **345,630** (human C-reactive protein).
  No network calls — the 756 MB `gene2pubtator3.gz` bulk file.

  **THE SELECTION RULE MUST MATCH THE VALUE READ, AND IT IS WORTH HALF THE SIGNAL.** Reading
  PubTator off the donor chosen for most CURATED papers scores **−0.0132** on the held-out
  control; `pubtator_donors()` picks the donor with most PUBTATOR papers and scores **0.3722**,
  beating the curated column's **0.3398** on identical folds. A mismatched rule looks exactly like
  "this source is weak" — the largest measurement error made on this axis, made twice before being
  caught. `pubtator_donor_*` records which donor was used.

  **The gene-SYMBOL routes were built, measured and REJECTED** (all-species, species-scoped with
  strain→species rollup, and free-text): `@GENE_<SYMBOL>` is missing for most bacterial genes
  (0 for `fnbA`, which has 1,158 papers; 31.5% of donor symbols return 0, 446 of them on donors
  with ≥5 curated papers), and free text matches abbreviations and strain names (`tam` → the
  Translocation Assembly Module, `nimR` → `AH1-NIMR`). **`NOT @SPECIES_9606` is a trap** — it
  strips 76-85% of virulence/clinical literature, exactly our biology.

  **A blank is NOT a zero**: blank = no in-scope donor carried a GeneID (unreachable); 0 =
  PubTator annotated 36M abstracts and never linked that gene. Coverage of scored proteins ~95%.

  **Three data traps.** `gene2pubtator3.gz` is **PMID-sorted**, so its head looks human-centric —
  that is how an earlier untested rejection happened; **the bulk file and the `search/` API
  disagree** (`ftsZ` 13 vs 2), the bulk being complete; and the download **truncated at 88 MB of
  756 MB while `curl` exited 0**.

  **THE NUMBER IS A PAPER COUNT** — curated SwissProt references, nothing scaled or blended. **A
  0-1 composite shipped first and was REMOVED on 2026-09-22; do not reintroduce it**: the same
  value meant 13 papers at annotation 3 and 83 at annotation 1, the two halves **double-counted**
  (r 0.64–0.70), and the weights did the opposite of what the code claimed. `scaled()` derives a
  0-1 column on the fly and is deliberately **not stored**.

  **`--max-target-seqs` was the worst bug in the axis and the lesson generalises.** DIAMOND
  returns the top k by **bitscore** = the CLOSEST relatives, but this stage wants the
  **best-cited** homolog, which for a conserved protein is a distant model-organism entry. At
  k=50 Sa `groEL`'s donor was *B. subtilis* GroEL with 9 papers and E. coli GroEL (246) never
  appeared — **systematically understating studiedness precisely for the most conserved families.**

  **Five evidence tiers, and TWO of them score 0 meaning different things**: `swissprot_direct`
  (≥95%) ⊃ `swissprot_close` (≥60%) ⊃ `swissprot_homolog` (≥40%), then **`below_floor`** (an
  in-scope hit under the floor) and **`no_hit`** (nothing in SwissProt at all). `no_hit` is the
  strongest novelty claim the axis makes. **A zero is an answer, not a gap; never impute it.**

  **The floor is 40%, and the control's own optimum agrees with it** — rho peaks there and falls
  away on both sides. **THE CONTROL IS A FLOOR-CHECK, NOT AN OBJECTIVE**: it correlates
  deliberately different quantities, so it has a ceiling well below 1. Measured: **spearman
  0.3398 / pearson 0.3164 over 2,945 of 4,403 proteins**, floor 0.25, exiting non-zero below.

  **DONOR SCOPE: not eukaryotic, and "restrict to Bacteria" is the obvious rule that is WRONG.**
  `prokaryotic` (default) admits any donor whose lineage lacks Eukaryota — Bacteria, Archaea
  **and phages**. Under `any`, five of Kp's ten highest-scoring proteins took human donors; but
  **strict Bacteria strands prophage proteins**, 70 of the 93 it leaves donorless having lost
  theirs to a virus.

  **Unknome is EVIDENCE, kept and not promoted**, and its trap is severe: the per-protein download
  reads **`knownness = 0.000` for all 1,882 S. aureus entries** including `clpP`, `rpoB` and
  `ftsZ`, while the **cluster** table reads 10.6 for the same protein — joining the per-protein
  file would have silently zeroed a proteome. **REJECTED and not to be re-added: training ESM-C on
  Unknome knownness.**

  **`studiedness/confounds.py`** → `evidence/confounds.tsv`: the axis measured against every
  other axis, fitting nothing. **Read `rho` beside `scored_rho`** (scored tiers only) — every
  degradability correlation **collapses to ~0 under it** (−0.15 → +0.008), so it was the
  zero-tier block, not the probabilities. **AUROC below 0.5 is a DIRECTION, not a failure.**
  **THE UNPRICED CONFOUND IS ESSENTIALITY** — `geptop_ess` rho **0.38–0.50**, and alone among the
  axes it **survives stratification** (0.35–0.39), so stacking the two double-counts.
  **Degradability and localization are clean.**

  **`studiedness/orthodb_transfer.py` is a REJECTED ALTERNATIVE, kept so the comparison is
  reproducible.** It attacks the axis's hard ceiling — 1,961 Kp proteins (34%) and 1,276 Sa (44%)
  have NO curated donor at 40% identity — via OrthoDB orthogroups, which 79% of them do have.
  **The transferred count carries NO information**: held-out E. coli with the anchor's genus struck
  out of the MEMBER pool, reached-only spearman **−0.1068** (narrow −0.0852, domain −0.0059), while
  a bare reached-vs-not flag scores **0.5058** — so the apparent 0.4140 over all rows is presence/
  absence, not literature. **Filling those cells would be worse than leaving them blank.** Reach
  was capped anyway at 22.9% because **only 34% of SwissProt is in OrthoDB**.
  **THE TRAP, worth remembering beyond this axis: the first control read 0.8351 and was
  self-correlation** — the sequence control excludes Escherichia from the DONOR pool, but a group
  contains the query's OWN entry, so **98.3% of "donors" were Escherichia**. A leak is invisible in
  a correlation and was caught only by listing donor organisms; `validate()` now exits non-zero if
  any same-genus donor survives. **Keep `scratch/orthodb_{gene2ac,og_members}.tsv`** — they cost
  two streamed passes over ~9 GB and make a re-run seconds.

  **A finding for the collaboration: the consortium's own panel is NOT novel** — `src/interest.py`
  sits at the 80th percentile (median) on Kp. Details: `docs/studiedness.md`.

## Presentation figures live in `plotting/`, not `scripts/plots/`

**Two folders, two jobs, and they must not be merged.** `scripts/plots/` holds the PER-AXIS
DIAGNOSTICS — written to defend a stage, cited from `docs/<task>.md`. **`plotting/` is top-level,
CROSS-AXIS and presentation-only**: it joins degradability to localization to essentiality to
orthology to answer the only question the GraDi consortium asks, *which K. pneumoniae proteins
should we pursue*. Outputs go to `output/plots/presentation/`, never to `output/plots/<task>/`.
The record is `docs/presentation.md`.

**`REPO_ROOT = Path(__file__).resolve().parents[1]` in `plotting/`, NOT `parents[2]`.** Rule 1 of
*Four rules the layout depends on* says `parents[2]` because a stage script sits two levels
deep under `scripts/<task>/`; a `plotting/` script is ONE level deep. Copying the `parents[2]` line resolves `REPO_ROOT` to the repo's
PARENT, and the resulting ImportError names `src`, not the path — so it reads as a broken conda env.
The check that proves it is running each script `--help` from ELSEWHERE.

**NO COMPOSITE SCORE, EVER — a plotting script is exactly where one would sneak back in.** Three
axes each REMOVED one: essentiality dropped `essentiality`/`essentiality_source` (it MIXED UNITS, a
measured 1.0 against a predicted 1.0, and on Kp was a verbatim copy of `screens_ess_mean`);
studiedness removed its 0-1 composite on 2026-09-22 (the two halves double-counted, r 0.64–0.70);
pockets ships no `druggability()` because *"no defensible weighting exists"*. A weighted rank is
invisible once drawn — a ranked list looks equally plausible whatever weights produced it.

**So the shortlist is a CASCADE OF STATED PREDICATES**, defined ONCE in `plotting/filters.py` and
imported, never restated: `build_filters()` · `cascade()` · `leave_one_out()` · `load_joined()`.
Measured on Kp, 2026-10-04:
**5,728 → 4,294 (not membrane) → 3,497 (no human ortholog) → 235 (essentiality top decile) →
59 (ADEP4 ≥ 0.328)**.

**`leave_one_out()` is not decoration — the filters are CONJUNCTIVE, so the funnel overstates the
early rules.** Without `essential` the cascade gives 792; without `not_membrane` it gives **63
against 59**. The membrane rule costs almost nothing because the degradability model already
learned the mechanism, which is why `not_membrane` is a statement about mechanism and not a filter
that earns its keep numerically.

**`plotting/` USES NONE OF STYLIA'S ERSILIA BRANDING — npg PALETTE, `article` STYLE, BLACK LABELS.**
Standing instruction from the project owner, 2026-10-04, given in two parts.

**The palette**: colours come from `stylia.CategoricalPalette("npg")` via **`plotting/palette.py`**,
never from `stylia.NamedColors()` (the ersilia plum/orange/mint set). Roles are SEMANTIC —
`PRIMARY`/`SECONDARY`/`TERTIARY`/`ACCENT`/`MUTED`/`INK`/`BACKDROP`/`SEQUENTIAL` — because `NPG[5]`
at a call site says nothing about intent and changes meaning silently if the palette is reordered.
`SEQUENTIAL` exists so heatmaps are npg too: `BuPu`/`viridis` are not, and a deck that is npg
everywhere except its heatmaps is the drift this module prevents.

**The style**: **`stylia.set_style("article")`, NOT `("ersilia")`** — the ersilia style paints every
text element AND every spine plum (`#50285A`), and the owner asked for black labels. `set_format`
and `set_style` are different knobs from the palette; only `article` gives black. `PAL.INK` is
likewise plain black, for in-plot annotations.

**NO `abc=` PANEL LETTERS.** `stylia.label()` takes an `abc` argument that stamps A/B/C on each
panel; `plotting/` does not use it (owner's call). Figures are referred to by what the panel shows,
never by a letter — a letter in a docstring goes stale the moment a panel is reordered. Note the
per-axis figures under `scripts/plots/` DO still use `abc=`; this rule is `plotting/` only.

**`plotting/palette.py` is the SECOND shared module, and centralising it is a deliberate departure
from `scripts/plots/`**, where style constants are copy-pasted by design. The reason the exception
holds: a palette stops being per-script taste once the owner sets it as a rule, and **a deck
rendered half in one palette and half in another is a defect no single script can see**.
`plotting/filters.py` is the other shared module, for the same class of reason — it carries a
claim, not a preference.

**`src/localization.py:LOC_CLASS_COLOR` IS ALREADY npg** (#E64B35, #00A087, #3C5488, #F39B7F,
#8491B4, #7E6148 are ggsci `pal_npg`) and is **imported, never redefined** — restating those six in
`plotting/` is how the deck drifts away from the localization axis's own figures.

**Capitalise axis labels** (`"Proteins remaining"`, not `"proteins remaining"`) — a standing
convention, and one these figures got wrong on the first pass.

**Four stylia traps, every one silent, all four hit while building this axis:**

1. **`create_figure(width=, height=)` takes FRACTIONS OF THE FORMAT SIZE, not inches.** `width=13`
   asks for thirteen slide-widths and yields a **4-gigapixel, 16 MB PNG** — and the script exits 0.
2. **`label(..., xlabel=None)` writes the literal placeholder `"X-axis / Units"`** onto the figure.
   Pass `""` for no label.
3. **`label(..., xlabel=None)` writes the placeholder, and `set_style` gates `NamedColors()`.**
   Without `set_style("ersilia")` first, `NamedColors()` returns `ArticleColors` and `NC.plum`
   raises `AttributeError`, reading as a broken install — moot here, since `plotting/` takes its
   colours from `palette.py` and must not call `NamedColors()` at all.
4. **Run plot scripts ONE AT A TIME** (the known `rmtree` trap). A tight loop over `plotting/*.py`
   reports a *different* spurious subset as broken each run.

Each script prints its own numbers and its own caveats after the figures, and **no figure asserts a
value its own script did not compute** — that is how these ten were validated against the measured
tables already in this file.

## Legacy

**`legacy/` holds the complete v1 pipeline, frozen. Do not extend it.** Start at
**`legacy/HISTORY.md`** — the retrospective: what was built and what was not, the methodological
decisions and why, ~70 documented traps, which data sources were obtained and how, and which
artifacts look like data but are not.

Three things there must not be trusted (`legacy/HISTORY.md` §7): the **Kp legacy degradability TSV**
(four verified defects, including an inverted N-end rule supplying 680 of 698 `medium` calls), the
**E. coli degradability values in the webapp** (a deterministic MD5 mock), and
**`data/raw/legacy/clp_substrates/`** (45- and 35-row hand-curated substitutes, not the papers' tables).

The archived scripts still run. They resolve paths via `parents[1]`, so `legacy/data` and
`legacy/output` are **symlinks** to the real trees one level up — that is what keeps them working
without editing 78 files. Without those symlinks a script exits 0 having written nothing. Run from
inside `legacy/`:

```bash
cd legacy && python scripts/09h_localization_plots.py --organism ecoli
```

The v1 webapp is still deployed from `legacy/app/` by `.github/workflows/pages.yml`.

## Setup

Conda environments (see `install.sh`). **Do not use the machine's default `python3`** — it resolves to
an unrelated `ersilia` env. Use `~/miniconda3/envs/gradi/bin/python`.

| env | what it holds | why it is separate |
|---|---|---|
| **`gradi`** (3.11) | the main env; `bash install.sh` runs `pip install -r requirements.txt`. Adds **`rdkit`** from PyPI for ligands | — |
| **`gradi-ortho`** (osx-64, Rosetta) | OrthoFinder + DIAMOND | no arm64 build. Its DIAMOND is borrowed project-wide via `GRADI_DIAMOND_BIN` |
| **`gradi-prokka`** (osx-64, Rosetta) | BLAST+ suite | no arm64 build; `rpsblast`/`blastp` borrowed via `GRADI_RPSBLAST_BIN` / `GRADI_BLAST_BIN` |
| **`gradi-emapper`** (osx-64, Rosetta) | eggNOG-mapper v2 | no arm64 build |
| **`gradi-loc`** (3.11) | DeepLocPro + TMbed | **`fair-esm` claims the same top-level `esm` package as EvolutionaryScale's ESM-C.** Pins `setuptools<81` and `transformers==4.44.2` |
| **`gradi-tabpfn`** (3.11) | `tabpfn==9.0.0` + `tabpfn-client==0.6.0` | torch 2.14 against gradi's 2.12. Needs `TABPFN_TOKEN` |
| **`gradi-lazyqsar`** (3.11) | lazy-qsar 3.4.4, for the degradability head comparison only | pins `numpy==2.1.3` / `scikit-learn==1.6.1`. **Rejected as an estimator**; kept so the comparison is reproducible |
| **`gradi-plm-ess`** (3.11) | ProteomeLM-Ess, the authors' **git** build (`proteomelm @ git+https://github.com/Bitbol-Lab/ProteomeLM` plus `httpx`, which `esm` needs and the install misses) | pulls **torch 2.14** against gradi's 2.12. The model code is numerically equivalent to gradi's installed build for our use — the split is about pip, not the model |
| **`gradi-pockets`** (osx-64, Rosetta) | `fpocket` 4.0 + `openjdk=17` for P2Rank 2.5.1 (tarball in `tmp/tools/p2rank_2.5.1`) | no arm64 build. `pockets/predict.py` runs in `gradi` and calls it across a process boundary (`FPOCKET_BIN`, `P2RANK_DIR`, `POCKETS_JAVA_HOME`) |
| **`gradi-pymol`** | `pymol-open-source` | ray-traced structure cartoons |

**Every split above is mandatory, not cosmetic.** The recurring failure is the same one each time:
installing a package with no osx-arm64 build (`blast`, `eggnog-mapper`) flips `gradi` to osx-64, or a
package claims a name or a torch version `gradi` already uses — and either takes stage 01's ESM-C down
with it. A stage that needs a foreign env **reaches across a process boundary to a `workers/` script**;
never activate one of these to run a stage.

No build, lint or test commands are configured yet. Document them here when added.
