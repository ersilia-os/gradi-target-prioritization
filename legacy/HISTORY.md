# GraDi target prioritization — v1 history and learnings

This is the retrospective for the **first** version of the target-prioritization pipeline, archived
here in `legacy/` in September 2026. Everything under `legacy/` is **frozen**: read it, cite it, copy
ideas out of it, but do not extend it.

Read this document before touching anything in `legacy/`. Sections 5 and 7 are the ones that will
save you time — §5 is the trap list, §7 is the list of artifacts that look like data but are not.

---

## 1. What this was

Ersilia's target-selection contribution to **Gr-ADI** — *"Exploring BacPROTACs as a new paradigm for
antibacterial discovery"* (Prof. Erick Strauss, Stellenbosch). The deliverable was a prioritized list
of proteins of interest for *K. pneumoniae* and *E. coli*, scored along five axes. Ligand
identification was explicitly out of scope.

**Anchor organisms.** *K. pneumoniae* HS11286 (`UP000007841`, 5,728 proteins) as the anchor;
*E. coli* K-12 MG1655 (`UP000000625`, 4,403) as a co-equal second organism and the curation hub;
human (`UP000005640`, reviewed canonical ~20,416) for selectivity.

**Persistence was split in two, deliberately.** Git tracked `src/ scripts/ docs/ assets/ app/` plus
the usual repo furniture. `data/` and `output/` were tracked by **eosvc** (DVC + S3) and gitignored.
`tmp/` was local-only scratch. The only generated artifacts in Git were `app/data/{kp,ec}.json`, the
webapp payloads.

**Directory contract.** `data/{raw,processed}` and `output/{results,plots}`, each bucketed
**organism-first** (`kpneumoniae/`, `ecoli/`, `human/`, `other/`, `legacy/`). v2 replaces this with a
stage-numbered convention (`data/processed/00_proteomes/`), so do not carry the organism-first layout
forward.

**Chronology.** Diagramming and spec work in May 2026 → a v1 annotation memo (2026-05-13) →
task-agnostic annotation `01`–`05` (June) → ligandability `06*` (late June) → essentiality `07*`
(July, including a 2026-07-13 pass that made *E. coli* first-class) → the webapp `08*` (late July) →
localization `09*` (early August) → degradability `10*` (4–6 August, left unfinished).

Consortium meetings: #1 2026-05-14 kick-off · #2 06-12 task-agnostic · #3 06-26 ligandability ·
#4 07-14 essentiality · #5 07-24 first draft of the prioritization browser.

---

## 2. The five axes — built vs not

### Axis 1 · Task-agnostic annotation (`docs/01_task_agnostic.md`)

Per-protein evidence consumed by every other axis. Ten specced tracks; **seven shipped**.

| Track | Status | Result |
|---|---|---|
| 1.0 reference proteome | done | `00a` |
| 1.1a PANTHER family | done | Kp 4,081/5,728 (71.2%), Ec 3,445/4,403 (78.2%) |
| 1.1b InterPro domains | done | Kp 4,973 (86.8%), Ec 4,204 (95.5%) |
| 1.2a PDB coverage | done | **Kp 30/5,728 = 0.5%**; Ec 1,779 = 40.4% |
| 1.2b AlphaFold pLDDT/PAE | done | Kp 5,727/5,728; median mean-pLDDT ≈ 91; ~45% multidomain |
| 1.3a BV-BRC PLFam/PGFam | **never built** | — |
| 1.3b within-Kp pan-genome class | **never built** | the memo calls conservation "a placeholder" |
| 1.3c cross-species broad-spectrum | **replaced** | done via 03a/03c OrthoFinder, not BV-BRC/OrthoDB |
| 1.3d selectivity vs human | done | 03a DIAMOND human tier + 03c |
| 1.4 bibliometric popularity | done, **re-specced** | HS11286 is bibliometrically dark (UniProt pub counts flat at 1 = the genome paper), so studiedness is **transferred from the best-characterised ortholog**; Europe PMC dropped |
| 1.5 embeddings | **upgraded** | ESM-2 650M → **ESM-C 600M** (1152-d, mean-pooled) |
| 1.5b ESM Atlas coords | added | 5,727/5,728 located |

Never built from the suggestions list: Foldseek/ProstT5, PPanGGOLiN, eggNOG-mapper v2, MobiDB-lite at
this axis, SaProt, canSAR.

### Axis 2 · Ligandability (`docs/02_ligandability.md`) — **complete, both organisms**

`ligandability_score = 0.45·binding + 0.30·structural + 0.25·pocket`, with **evidence-driven tiers,
not score thresholds**.

Tracks: 2.1a ortholog expansion (reuses the 03a table) · 2.1b ChEMBL 37 + BindingDB · 2.2a PDB
co-crystals · 2.2b AlphaFill · 2.3b fpocket 4.0 + P2Rank 2.5.1, pLDDT-weighted · 2.4 disorder filter.
**2.3a AF2Bind was deferred** — `06f` emits a NaN placeholder so the merge schema stays uniform.

Results: Kp 3,335 tractable / 884 partial / 1,509 intractable; Ec 2,787 / 619 / 997. Hard evidence
429 Kp / 1,240 Ec. Prime shortlist 1,949 Kp (306 with hard evidence).

### Axis 3 · Degradability (`docs/03_degradability.md`) — **most documented, least complete**

Re-specced twice: 2026-08-05 (E. coli first) and 2026-08-06 (the ClpP protease decision, §4).

Built: 3.0 machinery census · 3.1 initiation regions · 3.2a–d degron rules (partial) · 3.3 measured
turnover/attribution · 3.3c activated-ClpP proteomics.
**Not built**: 3.2e–i envelope rules, 3.4 biophysics, 3.5 abundance, 3.6 assembly state,
3.7 compartment routing, **3.8 the composite merge**, 3.9 `tpd_advantage`.

**There is no `*_degradability.csv` on disk for either organism.** The two-bar composite exists only
as prose. The webapp's `comp_degradability` is not this axis — see §7.

### Axis 4 · Essentiality (`docs/04_essentiality.md`) — **complete, both organisms**

`0.40·experimental + 0.20·transfer + 0.40·predictor`, **renormalised over available tracks, never
zero-filled**.

Tracks: 4.1a/b/c Kp Tn-seq + CRISPRi (`07b`, `07l`) · 4.2a *E. coli* transfer via the
Enterobacteriaceae-TraDIS compendium (`07c`) · 4.3a **ProteomeLM-Ess** as primary (`07d`; the `-Ess`
head is unreleased upstream, so a logistic head was trained locally on curated EcoGene labels —
5-fold CV **AUROC 0.809**) · 4.3b Geptop 2.0 reimplemented with DIAMOND (`07e`) · 4.3d FBA on
iYL1228/iML1515 (`07f`). **4.3c DeeplyEssential deferred** (`07g` NaN placeholder — no released
weights, py2/TF1.6, no licence).

Plus a prediction-free "publication view" (`07l`/`07m`) and *E. coli* first-class parity
(`07n`/`07o`, ingesting Keio, Goodall TraDIS, four CRISPRi screens, and RB-TnSeq's 280-condition
matrix).

Tiers — Kp: 401 essential / 835 likely / 4,492 non; Ec: 402 / 135 / 3,866. **154 Kp "prime" targets**
(essential ∧ tractable ∧ broad-selective).

Not built: 4.4 graded CRISPRi-vulnerability refit, 4.5 synthetic lethality, OGEE/DEG direct lookup,
Bacformer.

### Axis 5 · Expression and localization (`docs/05_expression_and_localization.md`)

**§5.1 localization: complete at 100% coverage, both organisms** (was 40.4% Kp / 51.4% Ec).
**§5.2 expression: never started.**

Tracks: 5.1a UniProt with an ECO experimental split (`09a`) · 5.1b **PSORTb 3.0 taken precomputed
from PSORTdb** (`09b`) · 5.1c **DeepLocPro** as primary predictor (`09c`) · 5.1d TMbed topology
(`09d`, 133 β-barrels) · 5.1e SignalP 6.0 — **not installed (licensed)**, so `09e` runs a calibrated
lipobox fallback (P 0.74 / R 0.82 / F1 0.78) · 5.1f STEPdb 2.0 + ortholog transfer (`09f`), **which is
what takes Kp from 1 experimentally-evidenced localization to ~1,680**.

Compartments (Kp): cytoplasm 3,470 · inner membrane 1,425 · extracellular 389 · periplasm 234 ·
outer membrane 193 · cell-wall surface 17. Zero null `clp_accessibility`.

---

## 3. Pipeline stage-by-stage

Canonical outputs are `output/results/<organism>/<prefix>_*.csv`, prefixes `kp` / `ec`.

| Stage | Scripts | Ran? | Canonical output |
|---|---|---|---|
| **00** | `00a_fetch_proteomes`, `00b_proteome_descriptors` | both | `data/raw/<org>/proteome/*.{fasta,tsv}` |
| **01** | `01a_esmc_embeddings`, `01b`, `01c_esmatlas_coords`, `01d` | all | `<prefix>_esmc600m_embeddings.npz`, `_projection.csv`, `_esmatlas_coords.csv` |
| **02** | `02a_interpro`, `02b_panther`, `02c` | all | `data/processed/<org>/families/*.csv` |
| **03** | `03a_orthology_general` (**`gradi-ortho`**), `03b`, `03c_orthology_focused`, `03d` | all | `data/processed/other/orthology/<prefix>_orthologs_long.tsv`, `three_way_orthogroups.tsv` |
| **04** | `04a_alphafold_structures`, `04b`, `04c_pdb_coverage`, `04d` | all | `<prefix>_alphafold_structure.csv`, `<prefix>_pdb_coverage.csv` |
| **05** | `05a_popularity`, `05b` | both | `<prefix>_popularity.csv` |
| **06** | `06a`–`06o` (15) | all but `06f` | **`<prefix>_ligandability.csv` + `_shortlist.csv`** |
| **07** | `07a`–`07o` (15) | all but `07g`; `07j` retired into `07k` | **`<prefix>_essentiality.csv` + `_shortlist.csv`** |
| **08** | `08a_webapp_export`, `08b_validate_export` | both | `app/data/{kp,ec}.json` — the only generated files in Git |
| **09** | `09a`–`09l` (12) | all; `09e` in fallback | **`<prefix>_localization.csv` + `_shortlist.csv`** |
| **10** | `10a`–`10o` (14 files) | **partial** | `<prefix>_deg_degrons.csv`, `_clpp_activator.csv`, `_disorder.csv`, `_deg_measured.csv`. **No merge.** |

**Two deliberate NaN placeholders**, `06f_af2bind.py` and `07g_deeplyessential.py`. They emit a NaN
column so the merge schema stays uniform whether or not the track exists. This pattern is worth
keeping.

**Stage-10 numbering is broken** — the single worst thing about resuming that axis:

- **Two files share `10e`**: `10e_measured_turnover.py` and `10e_activator_features.py`.
- **Two files share `10l`**: `10l_activator_plots.py` and `10l_activator_feature_plots.py`.
- Gaps at `10g`, `10h`, `10i`.
- **Three mutually inconsistent numbering schemes exist**: the filesystem; `docs/03_degradability.md`
  §"Planned outputs" (which names `10e_ecoli_turnover`, `10f_biophysics`, `10g_assembly_state`,
  `10h_compartment_channels`, `10i_degradability_merge` — none of which exist under those names); and
  `docs/degradability_report.md` §7 (an older scheme again: `10b_terminal_disorder`,
  `10c_degron_motifs`, `10d_ecoli_turnover`).
- `10k`/`10o` plots are **E. coli only**; `--organism kpneumoniae` runs but was never checked.
- The `10b` CSVs on disk are **stale** — from 4 Aug, predating the 6 Aug Flynn N-M1/N-M3
  implementation. Re-running `10b` picks them up with no code change.

---

## 4. Decisions and their rationale

**Anchor strain = HS11286.** `UP000007841`, `GCF_000240185.1`, `KPHS_*`, 5,728 proteins (5,316
chromosomal + 412 plasmid). **The only Kp proteome UniProt flags "Reference and representative".**
Rejected: **ATCC 43816 / KPPR1** (the originally suggested strain — popular in mouse in-vivo
essentiality work but has no curated UniProt reference proteome) and **MGH 78578** (the historical
reference holding most reviewed entries, but reviewed coverage is ~0.5% species-wide and propagates
by orthology anyway).

**UniProt accession is the canonical key, always** — every dataset, `data/` and `output/`. Locus tags,
b-numbers, JW ids and RefSeq are retained as provenance columns only. `KPHS_*` is a locus tag, not an
accession.

**Mapping rule**: resolve all orthology onto HS11286 where an ortholog exists; carry evidence from
other strains (KPPR1, KPNIH1, ECL8, BW25113) as an annotation on the HS11286 protein with the source
strain preserved. Species-level resources (ChEMBL taxid 573) are strain-agnostic.

**Sequence mapping, not accession matching — the single most consequential technical decision.**
HS11286 is a dark TrEMBL proteome whose accessions rarely equal those used by ChEMBL, BindingDB and
PDB-SIFTS. Exact-accession matching silently missed *direct* Kp data including the clinically central
SHV/OXA-48/CTX-M/NDM/AmpC β-lactamases. The fix: DIAMOND blastp against external target sequences,
with **≥95% identity = "direct"**. Validated on `A0A0H3H184` → ChEMBL `Q93LQ9` at 100% identity, 104
potent compounds plus PDB 5eec/6d15 co-crystals, all previously invisible.
- **Refinement (2026-06-25)**: the bacterial bucket must be *true Bacteria*, not "non-human", plus a
  **40% identity floor** for transfer. Without it, heavily-screened rat/mouse/eel/parasite targets at
  ~30% identity inflated counts (the first pass gave 424 Kp ChEMBL-potent proteins with rat `P97697`
  as a "best hit"). Corrected: 175 Kp / 93 BindingDB.

**Two orthology methods for two distances.** OrthoFinder (phylogenetic orthogroups) within bacteria;
**DIAMOND similarity for bacteria→human**, because orthology inference is unreliable across that
distance and the subtractive-genomics standard is homology-flagging.

**Renormalise, don't zero-fill.** A sub-score entirely absent for a protein is *dropped* and the
remaining weights renormalised — a missing measurement lowers confidence, it must not push the score
toward zero. This is `07h`'s idiom and became the house rule; `docs/03_degradability.md` explicitly
says to follow `07h`, **not** `06g`'s zero-fill.

**Evidence-driven tiers, not score thresholds** (`06g`, `07h`). E.g. `tractable` if there is any ≤1 µM
activity or an own/≥95%-identity co-crystal, regardless of composite score.

**ProteomeLM: train our own head.** The published `-Ess` head is unreleased. Labels are the **curated
EcoGene set, not OGEE**, which side-steps the OGEE-parroting concern. Reuses the 01a ESM-C embeddings,
so the forward pass is 2–3 s per proteome.

**Geptop reimplemented with DIAMOND** — the upstream py2/NCBI-BLAST `.rar` is unbuildable headless.
Reference weight = median RBH %identity as a data-derived proxy for Geptop's composition-vector
distance (a documented simplification), cutoff 0.24.

**Compendium over gated primary.** Goodall 2018's ASM/PMC supplement is unfetchable headless, but its
data lives cleaner and richer in the **Enterobacteriaceae-TraDIS compendium** (the same group's own
GitHub release), which additionally upgrades §4.2a to a graded 12-genome cross-species consensus.

**Localization: evidence before prediction.** Precedence is `uniprot_experimental` →
`ortholog_transfer`/`stepdb` → `stepdb_curated` → `uniprot_curated` → `deeplocpro` → `psortb` →
`unknown`. TMbed/SignalP may **overrule a predicted call only — never an experimental or curated
one** — in the two cases they settle definitively (β-barrel ⇒ outer membrane; Sec/SPII lipoprotein
sorted by the Lol "+2 rule"). Consequence: **`psortb` never wins**, by design; its role is
corroboration via `predictor_agreement`, not adjudication.

**DeepLocPro replaces "run PSORTb ourselves"; PSORTb is kept but taken precomputed.** Two independent
stacks (SVM/motif/BLAST vs a protein language model) means agreement is genuine corroboration rather
than two views of one embedding. DeepLocPro beats PSORTb 3.0 on the post-2010 Gram-negative benchmark
(accuracy 0.74 vs 0.34, macro-F1 0.75 vs 0.35, MCC 0.69 vs 0.30).

**`clp_accessibility` ladder**: cytoplasm 1.0 · inner membrane with ≥30% cytoplasm-facing residues
0.6, else 0.2 (0.4 when topology is missing) · periplasm 0.2 · outer membrane / β-barrel /
extracellular / cell surface 0.0.

**Categorical transfer needed its own rule.** `transfer_ecoli_to_kp` reduces with max/mean, which is
meaningless for a class label; `transfer_categorical_ecoli_to_kp` uses **donor consensus and abstains
on ties** (only 81 of 3,179 anchors have >1 E. coli ortholog).

### The protease decision (2026-08-06) — the axis-defining call

From the funded Gr-ADI research vision (`Draft-v5(Final)`): the handle is **ClpP engaged directly by a
small-molecule activator, with NO unfoldase partner**. WP2 is literally named *WP2_ClpPELs* (ClpP-
Engaging Ligands); ClpC is explicitly ruled out for Gram-negatives. Consequences:

1. **Activated ClpP cannot unfold anything.** The activator occupies the ClpX/ClpA docking cleft,
   opening the pore, but nothing pulls. The ~5/~20/~37 aa initiation-region rule the whole spec was
   built on is an **unfoldase rule**. Under this modality the substrate must **already be
   unstructured** — so global disorder, `two_domain_architecture` and conformational stability become
   first-class, and the nine-channel protease survey is demoted to a fallback register.
2. **Two bars, never averaged.** The WP1 criterion changed between drafts (v3 "partnerless activated
   ClpP" → v5 "substrate of a ClpP protease *complex*") while the validation assay stayed partnerless.
   So the axis emits `clpP_complex_substrate` (selection, well populated) and `partnerless_clpP`
   (validation, sparse) plus `bar_disagreement` — "naming the project's likeliest wet-lab
   disappointments in advance IS the deliverable".
3. **ClpC/McsB are absent from both organisms**, so every published (ClpC1-based) BacPROTAC describes
   a machine they do not have.
4. **Signal-peptide status flips sign.** A Sec precursor is translocation-competent and therefore
   unfolded, so it marks a real (time-limited) **pre-export window** for a cytoplasmic protease.

**E. coli is PRIMARY for degradability — an inversion of every other axis.** Every experimental
dataset in the field is E. coli; Kp values are transferred, capped at **3,179/5,728 = 55.5%** by
orthogroup coverage. Building Ec first also buys an internal validation set. The machinery is 76–99%
identical between the two (re-derived by pairwise alignment against a 20–32% random-pair baseline).

**Cleavage weighted above abundance (0.65/0.35)** in `10c`: a 24 h abundance drop conflates
degradation with growth arrest, regulon change and resynthesis; a rise in endogenous-protease-generated
peptides is a direct proteolysis product. The case that settles it: **AcpP under ONC212 = +0.03
abundance (nothing) but +3.46 cleavage** — an abundance-only score would have called the project's own
second target a non-substrate.

**Engageability vs depletability pull in opposite directions** on natural half-life. Engageability
enters the score; **`resynthesis_burden` is a separate column, shortlist filter only**. GroEL is the
worked negative control (both termini face the barrel interior, hugely abundant, CLIPPERs could only
deplete it ~40%).

**`compartment_handle` is a router, not a filter** — periplasmic proteins are scored on the DegP
channel, not zeroed. Zero/NaN is reserved for `handle == none`.

**Filing convention for manually-supplied papers**: **data** → `data/raw/<organism|other>/degradability/<key>/`
with a `SOURCE.md`; **article PDFs** → the shared pool `data/raw/other/degradability/literature/`. Both
under `data/`, so paywalled binaries never enter a public Git history. `SOURCE.md` is invisible to the
fetch ladder's `_has_data()` so it cannot be mistaken for the dataset.

---

## 5. Traps and pitfalls

The complete collection. Grouped by kind; the infrastructural ones bite first.

### 5.1 Infrastructural

1. **The 03a kp→ec ortholog table has entirely empty `pident` / `coverage` / `bitscore`.** It is
   OrthoFinder orthogroup membership only. Any axis that thresholds transfer on percent identity
   **silently drops everything**. `pident` exists *only* for the human/DIAMOND tier.
2. **`fair-esm` collides with EvolutionaryScale `esm`** — both claim the top-level `esm` package.
   Installing DeepLocPro into `gradi` silently breaks `01a`. Hence the mandatory `gradi-loc` split.
3. **Do not use the machine's default `python3`** — it resolves to an unrelated `ersilia` env.
4. **PSORTb under Docker hangs ~19 h under Rosetta on Apple Silicon** and never emitted a single
   prediction. Route around it via PSORTdb's precomputed per-genome tables.
5. **DeepLocPro's own `-d mps` flag is broken** — `EnsembleModel.embed_batch()` gates device placement
   on `torch.cuda.is_available()`, so weights move to MPS while tokens stay on CPU. `09c` drives
   `EnsembleModel` directly instead of the CLI.
6. **TMbed runs on CPU** for the same reason. `09d` shards and caches instead, and **the shard cache
   key includes the shard size** so changing `--shard-size` cannot silently reuse mismatched shards.
7. **Load-bearing pins**: `setuptools<81` (DeepLocPro imports `pkg_resources`), `transformers==4.44.2`
   (TMbed's ProtT5 `T5Tokenizer` dies on transformers 5.x with a spurious tiktoken error),
   `xlrd>=2.0.1` (legacy BIFF `.xls` supplements openpyxl cannot open — and pandas' error points at
   the *missing package*, not the format).
8. **Runtime extrapolation trap**: TMbed measured ~17 min/shard while DeepLocPro competed for CPU;
   with the machine to itself it is ~4 min/shard. Do not extrapolate an ETA from a contended run.
9. **DeepLocPro truncates sequences at 2,000 residues** (ESM-2 cost is quadratic). Affected proteins
   carry `dlp_truncated`.
10. **DeepLocPro's GitHub LICENSE is CC BY-NC-SA 4.0** (non-commercial) while the paper states
    CC BY 4.0 — needs confirming with partners before any commercial deliverable.

### 5.2 Data acquisition

11. **MobiDB's bulk endpoint answers `405` to HEAD but `200` to GET.** Do not conclude it is broken
    from a HEAD probe.
12. **ASM journals return 403 to non-browser clients; PMC fronts binaries with a reCAPTCHA.** The
    working route was an authenticated Chrome session via chrome-devtools MCP `evaluate_script`
    same-origin fetch. Used for Jana 2023, Goodall 2018, Hawkins 2020, Rousset 2021.
13. **Nichols 2011 is the one un-fetched E. coli set** — its PMC copy sits behind an image reCAPTCHA,
    not solved on principle. Redundant with the RB-TnSeq 280-condition matrix.
14. **"The journal is paywalled" ≠ "the data is unreachable".** Conlon 2013, Cappelletti 2021 and
    To 2021 are all gated at the article but open at the data. Always probe Europe PMC
    `supplementaryFiles` and `media.springernature.com` first.
15. **Two Springer hosts behave differently**: `media.springernature.com` serves Conlon's supplement
    free; `static-content.springer.com` 403s to bare curl.
16. **figshare `ndownloader` answers HTTP 202 with a zero-length body** while preparing a file — and
    `raise_for_status()` does **not** raise on 202. Fixed in `10a`.
17. **`.xlsx` and `.zip` share the `PK` magic bytes.** The zip test ran first, so workbooks were
    opened as archives, found no data member, and were reported **gated despite HTTP 200**. Now
    discriminated by a top-level `[Content_Types].xml`.
18. **Nested archives were silently dropped** — MDPI ships the supplement as an inner `…-s001.zip`;
    `.zip` isn't in `DATA_EXTS`, so the outer pass reported "archive held no data members".
    `_extract_data_members()` now recurses one depth-capped level.
19. **PDF-only / text-only payloads were rejected** by the accept branch (MobiDB TSV, To 2021 SI PDF
    both bounced). Fixed with `_wants_raw()`/`RAW_EXTS`, opt-in by manifest.
20. **A 200 from the wrong PMCID is indistinguishable from a 200 from the right one.**
    `to2021_refoldability` pointed at `PMC8382223` — an unrelated *Anal Chem* paper about bioprinting
    hydrogels — and fetched its SI. Caught by a title sweep over cached PMC XML. **Verify PMCIDs by
    article title.**
21. **Article PDFs are not data.** Niwa 2012 and Calloni 2012 PDFs contain 0 b-numbers in extracted
    text — all numbers are in separate supplements. Verify by counting identifiers in the extracted
    text before treating a PDF as ingestible.
22. **MobiDB returns a header-only 61-byte HTTP 200 for a wrong proteome id.** `_fetch_mobidb()`
    therefore requires the header to start `acc` **and** ≥1,000 rows.
23. **PSORTdb POST needs a body** or it returns `411 Invalid Request` (`--data "id="` suffices). Its
    search endpoint is slow (>2 min for a genome-name query) — cache the assembly ids (`446671`
    HS11286, `449203` MG1655), never look them up in a loop.
24. **Do not use the PSORTdb bulk download** — 1.2 GB gzipped, **18.7 GB uncompressed as one TSV**,
    all Gram-negatives concatenated with no genome column. The per-genome route returns the same rows
    in seconds.
25. **STEPdb: pick the right file.** Only *"E coli K-12 strain MG1655 basic proteome.csv"* carries
    per-protein subcellular classes. The similarly-named *"Protein location and pseudogenes.csv"* is a
    pseudogene table, `@`-delimited, with padded junk columns.
26. **STEPdb is semicolon-delimited with semicolons inside free-text cells**, so a minority of rows
    shift columns; `classify_stepdb()` returns `unknown` for unrecognised patterns so shifted rows drop
    rather than mis-assign. It has **no evidence flag** — a non-empty `Annotation References` is the
    bench-support proxy (2,114/3,897).
27. **SignalP 6.0's direct `sw_request` CGI URL returns `Error 24`** — go through the service page and
    accept the licence.
28. **RegulonDB now fails** (both `network_tf_gene.txt` and `network_sigma_gene.txt`, 0/2 files) — the
    site changed again. Annotation-only and weight 0, so it gates nothing.
29. **Dead or hostile resources**: **PortEco's domain is now an online-gambling spam site — never cite
    or link it**; EcoGene returns HTTP 000; QSbio/anti-QSalign serve an empty Drupal shell; the
    original eSOL `tanpaku.org` host is dead (use the NBDC mirror); **OGEE v3's TLS certificate is
    misconfigured** and its SPA has no working API; **DEG does not index *K. pneumoniae* at all**;
    PHDB is `http://`-only with no bulk export; bioRxiv blocks scripted fetch (500).
30. **eSOL is b-number-keyed with no UniProt column**, so `solubility_esol` is **0% transferable to
    Kp** — not the ~31% arithmetic estimate.

### 5.3 Scientific and analytical

31. **The peripheral-membrane trap in UniProt CC text.** `"Cytoplasm. Cell inner membrane; Peripheral
    membrane protein."` (the DnaK archetype) was classified `inner_membrane` because the keyword
    cascade tested "inner membrane" first and **CC order carries no priority**. That understates
    exactly what the axis measures. 43 E. coli proteins carry peripheral+cytoplasmic-side; **40 were
    mis-binned**. Caught by the known-marker spot check — the argument for keeping that check.
32. **β-barrel threshold is ≥8 predicted strands.** TolC scores 4 and is correctly *not* flagged (it
    is a trimer contributing 4 per monomer).
33. **231 Kp "inner-membrane" proteins have `cyto_residue_fraction = 1.0` and no TM helix** —
    peripheral, membrane-associated, scored 0.6 deliberately. Easy to misread as a bug.
34. **DeepLocPro has training-set leakage on E. coli** (trained on UniProt 2023_03 + PSORTdb 4.0). Its
    84.1% agreement with UniProt-curated Ec labels is a **sanity check, not independent validation**;
    the meaningful number is DeepLocPro vs PSORTb at 92.6%.
35. **DeepLocPro is least decisive on extracellular proteins** (0.63 vs 0.82–0.93 elsewhere, leaking
    0.17 to cytoplasm) — and extracellular is where Kp has *no* experimental evidence. The axis's
    weakest corner.
36. **The lipobox parameters were tuned on E. coli and applied unchanged to Kp.** Tightening hydropathy
    to ≥2.0 collapses recall to 0.53 — lipoprotein h-regions are shorter and less hydrophobic than
    classic Sec/SPI.
37. **The naive N-end rule was implemented backwards.** `nterm_destabilizing` flagged residue 2 ∈
    {L,F,Y,W}, but MAP excises the initiator Met **only when residue 2 is small**, so a bulky residue 2
    means Met is *retained* and the mature N-terminus is Met (stabilizing). All 684 flagged Kp proteins
    have a bulky residue 2. It supplied **680 of 698 `medium`-tier calls**.
38. **The "ssrA-like" regex cannot match ssrA.** `[YAFWLIVM]?[ALV]A[ALV]A$` requires an alternating
    `X-A-X-A` tail; ssrA ends `…YALAA`, whose last four are `ALAA`. Fires on 5/5,728.
39. **`NM1_RX = ^M?[AILVMFW]{2,}` fires on 1,760/5,728 = 30.7%** — no discriminative power. And that
    regex is **not** Flynn's N-M1 (`T-X-K-[ILV]`, anchored 1–4 residues in), so the stated reason for
    dropping NM1/NM3 does not survive contact with the paper.
40. **`cterm_cm1_strict` is a strict subset of `cterm_cm1_broad`** — summing the two double-counts all
    11 proteins. Only **42 of 4,403** proteins carry any weighted motif.
41. **`degron_score` is exposure, not motif** — ρ(exposure, degron_score) = **0.999**. The composite is
    0.45·motif + 0.55·exposure but the motif term is non-zero for only 42 proteins.
42. **`cterm_cm1_broad` carries weight 0.40 on an odds ratio whose 95% CI (0.76–11.99) covers 1**,
    resting on 2 fast-turnover proteins out of 14 labelled hits — and it supplies all 42
    motif-contributing proteins. Meanwhile `nend_imet_cleaved` is the *only* feature whose CI excludes
    1 (1.18–3.13) and is correctly weighted 0 (41% prevalence = a prior, not a motif).
43. **Not one of the 11 ssrA-like proteins has an exposed C-terminus** (best `ppx` 0.49 against a 0.5
    bar). The two with the highest `degron_score` get it from the *opposite* terminus.
44. **The activated-ClpP evidence is *S. aureus*** (Conlon/Jacques), transferred by cross-phylum
    DIAMOND RBH at **median 42% identity**, reaching only 608 Kp / 609 Ec (10.6% / 13.8%). No native
    Gram-negative dataset exists — the Gr-ADI SoW flags the gap itself. `10o` labels every panel that
    uses it.
45. **The two bars are statistically independent** — `activator_evidence` vs Nagar half-life:
    **ρ = −0.073, p = 0.17, n = 352**, no monotonic trend by class (fast-vs-stable Mann–Whitney
    p = 0.38). So the 3,262-protein turnover set **cannot validate `partnerless_clpP`**, and `10d`'s
    best AUROC 0.587 is evidence about a different thing.
46. **And the direction may be actively wrong**: a good degrader target is *stable natively but
    degradable when forced*, so optimising against natural half-life could select the wrong proteins.
47. **The reproducibility ceiling is low.** ADEP4 vs ONC212 agree ρ +0.52 on abundance but only
    **+0.22 on cleavage**, **Jaccard 0.31** on the binary ≥2×-down call. Judge nothing against 1.0.
    Also `10c` **max-pools across the two activators**, which over two poorly-agreeing measurements
    inflates the score — flagged as needing revision. The 0.65/0.35 weighting has never been fitted,
    and cleavage is the *less* reproducible readout.
48. **The two half-life datasets barely agree** — Gupta vs Nagar **ρ = 0.11** across 866 shared
    proteins. Do not treat "measured half-life" as one quantity.
49. **The growth-correction cliff: 1,148 → 907 → 54.** After correcting for dilution by growth, only
    **54** E. coli proteins retain a positive proteolytic half-life. Most are diluted faster than they
    are proteolysed.
50. **ADEP4 abundance vs ADEP4 cleavage agree at ρ = 0.06** (n = 288, only 16 satisfy both) — the
    disagreement is between *readouts*, not chemistries.
51. **Nothing computed predicts activated-ClpP cleavage.** Every feature tried (motif score,
    degron_score, exposure, disorder fraction, longest internal run, initiation length, length) has a
    bootstrap AUROC interval covering 0.5, most point estimates *below* it. `10e`'s pre-registered
    gate: disorder adds **+0.009 over protein length alone**; a gradient-boosted model on 13 features
    (0.762) **loses** to logistic length-only (0.775); motifs attribute 0% of SHAP.
52. **SHAP splits credit among correlated features** — disorder collects ~45% of attribution yet adds
    nothing in nested models. The nested models measure what is *added*.
53. **~40% of E. coli cytoplasmic proteolysis is attributable to none of ClpP, Lon or HslV**
    (Gupta 2024).
54. **Gupta's KO panel has NO ΔftsH, and it is 13 conditions not 14** — measured from the fetched
    TableS1 (3,263 rows). Four documents had carried the error. FtsH attribution is unavailable from
    Gupta.
55. **Our Gupta protease attribution is wider than the paper's** (ClpP 189 / Lon 130 / HslV 29 vs their
    64 / 14 / 1) because they applied replicate-level significance testing the table does not expose —
    prefer the continuous `gupta_log2_*` columns.
56. **Two annotation traps for ClpC**: UniProt `A0ABY6X745` is auto-labelled `clpC` in
    *K. quasivariicola* but is a ClpB/ClpK misannotation (949 aa, lacks IPR001943); and *Appl Microbiol*
    2026;6:63 "clpC family genes in Kp heat survival" is almost certainly detecting *clpK*/*clpB* —
    **must not be cited as evidence of ClpC in Kp**.
57. **Kp ClpP is annotated at 194 aa vs E. coli's 207** — the HS11286 entry may be N-terminally
    truncated, which matters if the model is used structurally.
58. **`clpA` and `sspB` gene symbols are absent from HS11286** and would be missed by a symbol grep —
    hence track 3.0's insistence on a sequence-based census.
59. **A gene-symbol join loses ~27% of essentials.** HS11286 lacks canonical symbols for *lpxC, folA,
    dnaN, dnaB/C/E, fabB/D/G/I, accB/C* and ~80 more. This is the single most expensive recurring cost
    in v1 — v2's stage 00 exists largely to fix it.
60. **AlphaFold pLDDT misreads obligate complex subunits** — ordered in situ, disordered alone.
61. **`10b` uses `TERM_WINDOWS = (10, 20)`; the spec asks 15/30/50** — those exist only in `10d`'s
    `_disorder.csv`.
62. **Won 2024's N-vs-C asymmetry already failed to transfer** (`10d`: N-30 0.553 vs C-30 0.548) — a
    warning against importing it.
63. **LAMP-D / Stevens-Cullinane 2025 must be excluded from any label set** — it degrades NDM-1 in a
    Gram-negative but by **photo-oxidative backbone scission, explicitly without host proteolytic
    machinery**. A negative control, not a positive example.
64. **ProHL is a direct competitor** (deep-learning bacterial half-life classifier in E. coli) —
    **read its corrigendum before benchmarking**.
65. **`10o` panel 1 is a drawn schematic** — the only non-measured panel across either degron slide.
66. **Kp coverage estimates in `degradability_datasets.md` §8 are arithmetic, not measured — and both
    of the two checked were wrong** (`halflife_min` 48.0% not ~31%; `solubility_esol` 0% not 31%),
    erring in *both* directions.
67. **Schweke 2024 is *E. coli* O157:H7, not K-12** — intersection with `UP000000625` is **exactly
    zero**. Use Seq2Symm instead.
68. **STRING's E. coli taxon is 511145** (not 83333, not 562), and join via UniProt `xref_string`,
    **not** the aliases file (which maps one node to three accessions).
69. **Documentation drift**: `README.md` links `docs/pipeline.md`, which does not exist.
    `docs/mermaid_style.md` used the deleted `03_annotate_clp_degradability.py` as its worked example.
    The TraDIS compendium is described as "12 genomes" in two docs and "13 TIS libraries incl.
    *Enterobacter*" in two others — **the decoded reality is 12, with no Enterobacter**.
70. **The WebSearch budget was exhausted mid-session** on the degradability work; the remainder ran on
    Europe PMC REST + PubMed E-utilities + curl, which are strictly better for date-bounded sweeps —
    but they index abstracts/full text, so a purely web-surfaced resource (a lab's unpublished
    supplementary site, a bare GitHub data repo) could still have been missed.

---

## 6. Data sources — obtained vs never obtained

### Obtained and in use

**Proteomes** — UniProt REST stream: `UP000007841` (HS11286, 5,728), `UP000000625` (E. coli K-12,
4,403, 100% Swiss-Prot), `UP000005640` (human, **reviewed canonical ~20,416**, deliberately not the
~147k unfiltered set).

**Annotation** — InterPro 108.0 (`ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list`) ·
PANTHER 19.0 HMM_classifications (`data.pantherdb.org/ftp/hmm_classifications/current_release/`) ·
AlphaFold DB v6 API · PDBe SIFTS `best_structures` + `ligand_monomers` · Biohub ESM Atlas
(`POST https://biohub.ai/esm/protein/api/v1alpha1/umap/coords/by-hash`, ≤10 hashes/request, keyed by
`md5(sequence)`).

**Ligandability** — ChEMBL 37 SQLite (29 May 2026) · BindingDB All TSV 2025-04 · `pdb_seqres` full
(1,069,489 chains) · AlphaFill API 2.1.1 · fpocket 4.0 · P2Rank 2.5.1 (standalone Java tarball under
`tmp/tools/`).

**Essentiality** — Enterobacteriaceae-TraDIS compendium (GitHub `Gardner-BinfLab`, fetched
2026-07-13; `giant-tab_final.tsv`, PMID 39207104) · Eichelberger/Short 2024 ECL8 (eLife CDN — the only
fully scriptable publisher in the set) · Mike & Bachman 2023 + Bachman 2015/2025 (Europe PMC supp
zips) · **Jana/Zhu 2023 full supplementary via authenticated Chrome** (which also contains,
re-tabulated inside it, the Ramage 2017 KPNIH1 424-gene set and the Bachman 2015 in-vivo lists) ·
ProteomeLM-M (git, Apache-2.0) · `Geptop_v2.0.rar` (extracted with `unar`) · iYL1228 + iML1515 BiGG
models · **Keio via PEC `PECData.dat`** · Goodall 2018 / Hawkins 2020 / Rousset 2021 **via Chrome** ·
Rousset 2018 / Wang 2018 / Cui 2018 via Europe PMC · **RB-TnSeq Fitness Browser** direct from
`morgannprice.org/FEBA/Keio/` (3,790 genes × 3,511 conditions, 247 MB, reduced to a min-fitness summary
plus a 280-condition matrix; the raw file is not retained).

**The 12 TraDIS compendium genomes** (names only in `07l_publication_essentiality.py:129-143`, never in
the docs): Kp **ECL8** · Kp **RH201207** · Ec **BW25113** · Ec **ST131 EC958** · Ec **UPEC ST131
NCTC13441** · *C. rodentium* **ICC168** · *S.* Typhi **Ty2** · *S.* Tm **A130** · *S.* Tm **D23580** ·
*S.* Tm **SL3261** · *S.* Tm **SL1344** · *S.* Enteritidis **P125109**. Decoding: TraDIS log-ratio
≤ −0.5 = essential (F1 = 0.75 vs EcoGene). Coverage 3,170 Kp proteins, 249 core-essential (≥80% of 12).

**Localization** — UniProt REST (`cc_subcellular_location,ft_signal,ft_transmem,ft_lipid`) · PSORTdb
per-genome POST (`db.psort.org/search/results/download?assembly=446671`) · UniProt async ID-mapping
(RefSeq→UniProt, chunked at 5,000) · STEPdb 2.0 CSV · ESM-2 650M + ProtT5-XL-U50 weights (~4.75 GB,
auto-downloaded into `gradi-loc`).

**Degradability — the 2026-08-06 fetch run landed 172 MB across 13 datasets** via the `10a` manifest:
Gupta 2024 (Springer CDN, TableS1 3,263×42, header at row 3) · Cappelletti 2021 LiP-MS (Europe PMC,
60.6 MB) · Mateus 2018 TPP · Mateus 2020 (Supp Data 5, the 121-strain Keio panel) · MacKrell 2026 ·
Schmidt 2016 (30 sheets) · Györkei 2022 · Niwa 2022 (nested zip) · **eSOL**
(`dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip`, 4,132×28 — and it **also contains Niwa 2012's
chaperone experiment**, so Niwa 2012 needs no separate fetch) · Nagar 2021 (sd002 = 1,150 × 190, the
188-feature matrix) · Calloni 2012 (manual) · **MobiDB bulk for both proteomes**
(`mobidb.org/api/download?proteome=<UPID>&format=tsv`; Kp 85,234 rows / 5,728 accessions, Ec 73,479 /
4,404). Plus **Conlon 2013** (free at `media.springernature.com/…MOESM97_ESM.xlsx`) and **Jacques 2020**
(GSA figshare `10.25386/genetics.11873841`, article 11873841, file 21767271). Plus, supplied manually
as PDFs in `data/raw/other/degradability/literature/`: **Flynn 2003, Neher 2006, Ziemski 2021 Table S1,
Niwa 2012, Calloni 2012**.

### Never obtained

- **Nichols 2011** chemical-genomics S-score matrix — PMC image reCAPTCHA, not solved on principle.
- **Cain 2017** NJST258 Tn-seq — not fetched, low priority (the compendium lacks NJST258).
- **Feng 2013** (ACS) — 403, no PMC record. Skipped (ClpC-family).
- **Meltome Atlas per-protein T_m file** — the browser returns 200 but the *specific* file was never
  confirmed. **Do not cite until checked.**
- **PHDB bulk export** — no route exists.
- **Bienvenut 2015** E. coli mature N-terminome — Wiley paywalled, PRIDE holds RAW only.
- **CPB-ChaFRADIC C-terminome** (Chen 2020) — ACS paywall.
- **To 2021 / To 2022 machine-readable tables** — ACS-gated; only PDFs. To 2022's Europe PMC zip is
  **figures only** (200 OK, 380 kB, 12 jpg/gif, no data).
- **Humbard 2013** N-degradome, **Westphal 2012 / Arends 2016** FtsH traps — no Europe PMC supp zip.
- **Klimecka 2021** — paper open, but no consolidated table (12×4 numbers would need hand-curation
  from figures).
- **Bhat 2013** *Caulobacter* ClpP substrates — PMC3681837 returns `idIsNotOpenAccess`.
- **OGEE v3 / DEG** — no working programmatic endpoint; both repackage the primary screens anyway.
- **BV-BRC PATtyFam pan-genome counts** — the per-PLFam aggregation over ~30k Kp genomes was deferred
  and never done.
- **Illenseher 2025** Kp-native iBAQ abundance (`PXD052921`) and **Muselius 2020** Kp Δ*lon* proteome
  (`PXD015623`, *the only Kp-native proteolysis dataset in existence*, 26 accumulators) — routes
  verified, never ingested.
- **The Gr-ADI WP1 assay results** — the one dataset that would actually validate the degradability
  axis. Does not exist yet.
- **Erick Strauss's e-mail of 14 May 2026, "regarding proteolysis rates"** — flagged in
  `degradability_report.md` §12.8 as **unread**.

---

## 7. Do NOT trust these

1. **The Kp `comp_degradability` in the webapp is a frozen legacy file**,
   `data/processed/legacy/klebsiella_pneumoniae_clp_degradability.tsv`, written by a since-deleted
   regex script (`scripts/03_annotate_clp_degradability.py`, removed in `8327de3`). **Four verified
   defects**: an inverted N-end rule supplying 680/698 `medium` calls; an ssrA regex that cannot match
   ssrA (fires on 5/5,728); a CM2 with no signal; NM1 at 30.7% prevalence. The AlphaFold modulation the
   spec promised was never implemented. Net signal: 5 ssrA-like, 21 trap-flagged, 10 in `high`.
2. **E. coli degradability in the webapp is a deterministic MD5 mock** —
   `08a_webapp_export.py:degradability_frame()`, `h = _hash01(a, "deg")`. Flagged
   `PROVISIONAL = { degradability: ["ec"] }` in `app/config.js:37`. Its tier thresholds (0.5/0.15)
   differ from real Kp's (0.50/0.25), giving a purely artefactual 30-fold discrepancy in `high` rate
   (5.1% ec vs 0.17% kp). **The mock is honest — it is just not information.**
3. **`data/raw/legacy/clp_substrates/` are hand-curated substitutes, not the papers' tables.**
   `flynn2003_*.tsv` is **45 rows**, explicitly not a verbatim copy of Flynn 2003 Tables 1–2;
   `nagar2021_*.tsv` is **35 rows** against the real 1,149. Their column names
   (`ecoli_clp_trapped`, `ecoli_halflife_class`) imply otherwise. **Must stop being described as
   Flynn's census.** Superseded and unused by the `10e` measured layer.
4. **`degradability_report.md` §12 supersedes parts of §1–§7** — read §12 first. Its earlier "the
   consortium is deliberately protease-agnostic" conclusion is explicitly **retracted** (it came from
   Strauss's *other* grant, the TELL/ME proposal).
5. **The Gr-ADI proposal's own stated evidence for its targets does not fully hold**, measured from the
   actual tables: under ONC212, DnaK is 1.4× not ≥2×, and **AcpP is unchanged (+0.03)**. Both survive
   on cleavage evidence instead. **GyrB has essentially no activated-ClpP evidence at all**
   (`partnerless_clpP` = 0.0015; it maps cleanly at 44.8% identity, so this is absence of signal, not
   absence of mapping) — and GyrA/GyrB are the declared second targets. This was a **pre-registered**
   prediction that held.
6. **The degradability axis has no usable validation set for its own gate.** Measured, not argued.
7. **The 09e SignalP track is a fallback**, recorded as `signalp_source=lipobox`; it finds lipoproteins
   only and does not type Tat/SPI, Tat/SPII or Sec/SPIII.
8. **`08b_validate_export.py` requires `comp_degradability` and `degradability_tier`** while
   `comp_degradability` is **absent from the payload `components` list** in both `app/data/kp.json` and
   `ec.json` — an internal inconsistency.
9. **`clp_accessibility` is computed but never consumed by the degradability axis.**
10. **Disorder exists in the repo only as whole-protein pLDDT fractions, consumed solely as a
    *ligandability penalty*** in `06g`, despite the spec citing disordered termini as the dominant
    recognition signal.
11. **Stale artifacts**: `data/raw/<org>/localization/<prefix>_localization.tsv` (superseded) and
    `data/raw/<org>/localization/psortb/` (empty — the Docker run never produced output).
12. Two consortium notes conflict and were never settled: the kick-off suggested "going after targets
    in the periplasm", but ClpP is cytoplasmic, and the v5 proposal points the other way (GyrA/GyrB,
    "cytosolic"). The reconciliation offered was the pre-export window. Separately, **FtsH is listed as
    a consortium *target*** while also being one of only three experimentally essential proteases and a
    candidate handle.

---

## 8. Environments

| Env | Why it exists | Used by |
|---|---|---|
| **`gradi`** (py3.11) | the main env; `pip install -r requirements.txt` | everything except the four below |
| **`gradi-ortho`** (osx-64 bioconda, **Rosetta**) | OrthoFinder + DIAMOND have no arm64 build | `03a`, `03c` — and **its DIAMOND binary is borrowed by `06a/06b/06c`, `07b/07e/07f`, `09f`, `10c`** via `GRADI_DIAMOND_BIN` |
| **`gradi-pockets`** (osx-64, Rosetta) | `fpocket` + `openjdk=17` (JRE for P2Rank) | `06e` only (the script itself runs in `gradi`) |
| **`gradi-pymol`** | `pymol-open-source` for ray-traced cartoons | `06n` only, via subprocess |
| **`gradi-loc`** (py3.11) | **mandatory split, not cosmetic** — DeepLocPro needs `fair-esm`, which claims the same top-level `esm` package as EvolutionaryScale `esm` used by `01a`. Pins `setuptools<81`, `transformers==4.44.2`, `tokenizers<0.20` | `09c`, `09d` only |
| **`gradi-prokka`** (osx-64, Rosetta) | optional — re-annotate the ECL8 genome, whose `ecl8_*` locus tags come from a non-deposited Prokka annotation | `07b` only |

Other facts: **`unar`** (`brew install unar`) is needed for `Geptop_v2.0.rar`. **`poppler`** was
missing and would unblock PDF table extraction — one line, low risk. **ProteomeLM is not on PyPI**
(`pip install "git+https://github.com/Bitbol-Lab/ProteomeLM.git"`, Apache-2.0). **No build, lint or
test commands were ever configured** — a `.ruff_cache/` exists but no config.

**Everything is cached and resumable** — DIAMOND hits, ProteomeLM ctx embeddings and head, the FBA KO
table, DeepLocPro per-protein cache, TMbed shard cache, PDBe/AlphaFill responses, activator sequences.
A killed run resumes on re-invocation.

---

## 9. Index of `legacy/docs/`

**Specs** — `01_task_agnostic.md` (245) · `02_ligandability.md` (131) · `03_degradability.md` (488) ·
`04_essentiality.md` (143) · `05_expression_and_localization.md` (97).

**Run logs** — `ligandability_log.md` · `essentiality_log.md` · `localization_log.md` ·
`degradability_degron_log.md` · `degradability_measured_log.md` · `projection_exploration_log.md` ·
`{,ec_}interpro_annotation_log.md` · `{,ec_}panther_annotation_log.md` (the unprefixed pair is Kp).

**Reports and references** — `essentiality_report.md` (359) · `essentiality_references.md` (206) ·
`degradability_report.md` (937 — **§12 supersedes parts of §1–§7**) ·
`degradability_references.md` (985) · `degradability_datasets.md` (352 — the dataset/tool inventory and
column schemas; start here to build anything degradability-related) ·
`degradability_downloads.md` (489 — download routes) · `degradability_overview.pdf`.

**Conventions** — `mermaid_style.md` (146).

---

## 10. What to carry forward

Three things, in order of value.

1. **The `07h` scoring contract.** Renormalise over available tracks rather than zero-filling; return
   `np.nan` from a sub-score that has no evidence; keep weight constants as auditable module-level
   names; emit `_sources` provenance columns beside every score; drive tiers off evidence, not
   thresholds. This is the repo's best-earned pattern and every later axis was told to copy it.

2. **Sequence-first identifier resolution.** DIAMOND ≥95% identity = "direct", ≥40% floor for
   transfer. This is what made the dark HS11286 proteome workable at all, and the empty-`pident`
   orthology table (trap 1) is the counterpart that silently defeats naive transfer. v2's stage 00
   exists largely to front-load this: a gene-symbol join loses ~27% of essentials (trap 59).

3. **The degradability negative result — the most valuable thing in the repo.** Partnerless activated
   ClpP has no usable sequence degron (only 42/4,403 proteins carry a weighted motif;
   ρ(exposure, degron_score) = 0.999, so the score measures exposure not motif; 0 of 11 ssrA-like
   proteins has an exposed C-terminus). The only right-bar data is cross-phylum *S. aureus* at ~12%
   coverage and median 42% identity. The two bars are statistically independent (ρ = −0.073, p = 0.17),
   so the turnover set cannot validate the assay. And nothing computed predicts cleavage — disorder
   adds +0.009 over protein length alone, and a 13-feature GBM loses to logistic length-only.

   **Preserving that chain of measurement is what stops the next attempt from rebuilding a
   degron-based score.** If v2 wants a degradability axis, it needs new *data*, not new features.
