# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status

Active **target-prioritization** analysis for the **GraDi** collaboration: a per-protein
annotation/prioritization pipeline anchored on *K. pneumoniae* HS11286, with *E. coli* K-12 and
human as comparison organisms. `scripts/` holds a numbered pipeline (fetch proteomes →
language-model embeddings & 2D maps → family/structure annotation → cross-species orthology); the
five prioritization axes are specified in `docs/01_…`–`docs/05_…`. `data/` and `output/` are
organized **organism-first** (see *Directory contract*), and `requirements.txt` / `install.sh`
are populated (see *Setup*).

## Two-track persistence: Git vs. eosvc

This repo deliberately splits what is tracked where, and the split is enforced by `.gitignore`:

- **Tracked in Git**: `src/`, `scripts/`, `notebooks/`, `assets/`, `LICENSE`, `README.md`, `install.sh`, `requirements.txt`, `access.json`.
- **Tracked by [eosvc](https://github.com/ersilia-os/eosvc) (DVC + S3), NOT Git**: `data/` and `output/`. Both directories are listed in `.gitignore`.
- `tmp/` is local-only scratch space (gitignored).

Consequence: when adding datasets or generated results, write them under `data/` or `output/` so they go through eosvc, not Git. Do not commit data files into Git directly. Empty directories are preserved with `.gitkeep` so the structure survives an empty checkout.

`access.json` declares the visibility of `data` and `output` for eosvc (`"public"` or otherwise). Update it if access requirements change.

## Identifier convention

Always use **UniProt identifiers** (UniProt accessions, e.g. `P12345`) as the canonical identifier for proteins and genes in **every** dataset — both `data/` inputs and `output/` results. When a source provides only other identifiers (gene symbols, Ensembl/Entrez/RefSeq IDs, etc.), map them to UniProt accessions and use the UniProt accession as the primary key. Original identifiers may be retained as additional columns for provenance, but they must not replace the UniProt accession.

## Anchor strain

The single *K. pneumoniae* anchor strain is **HS11286** (UniProt proteome `UP000007841`; NCBI `GCF_000240185.1`; locus-tag prefix `KPHS_*`; 5,728 proteins). It is the **only** *K. pneumoniae* proteome UniProt flags as a "Reference and representative proteome", so it has the most complete annotation and the cleanest cross-references. The FASTA lives at `data/raw/kpneumoniae/proteome/UP000007841_HS11286.fasta` (fetch with `scripts/00a_fetch_proteomes.py`, which also fetches the E. coli K-12 `UP000000625` and human `UP000005640` reference proteomes used for orthology/selectivity).

Why not the alternatives:
- **ATCC 43816 / KPPR1** (the originally suggested strain) is popular only in mouse *in-vivo* essentiality work, not annotation — it has no curated UniProt reference proteome. Avoid it as the anchor.
- **MGH 78578** is the *historical* reference and holds most of the species' reviewed (Swiss-Prot) entries, but reviewed coverage is not a meaningful strain criterion here: across *all* K. pneumoniae strains there are only ~960 reviewed vs ~187k unreviewed entries (~0.5%), and reviewed annotation propagates across strains via orthology anyway.

**Mapping rule:** resolve all orthology and identifier mappings **onto HS11286** whenever a corresponding ortholog exists. Evidence and data from other strains (e.g. KPPR1, KPNIH1, ECL8 essentiality; *E. coli* K-12 orthology) should be carried as annotations on the HS11286 protein (keyed by its UniProt accession), preserving the source strain/ID as provenance columns rather than as a separate primary key. Note `KPHS_*` is a locus tag, not a UniProt accession — keep both, but the UniProt accession is canonical (see *Identifier convention*).

Resources keyed at the **species** level (e.g. ChEMBL, taxid 573) are strain-agnostic; the anchor-strain choice does not affect them.

## Directory contract

The two eosvc-tracked trees keep their template roles, but are organized **organism-first**:

- `data/raw/` — original, untouched inputs.  `data/processed/` — cleaned/transformed derivatives.
- `output/results/` — numerical results, logs, text.  `output/plots/` — figures.

Inside each of those four, content is bucketed **by organism, then by data type**:
- `kpneumoniae/`, `ecoli/`, `human/` — the focal organisms. E.g. `data/raw/kpneumoniae/proteome/`,
  `data/raw/kpneumoniae/{interpro,panther}/`, `data/processed/kpneumoniae/{embeddings,families,alphafold}/`,
  `output/{results,plots}/kpneumoniae/`. (Per-organism subfolders are created as each track is run,
  so not every organism has every subfolder yet.)
- `other/` — cross-/multi-species inputs that aren't a single focal organism; e.g. the orthology
  panel and its OrthoFinder run live in `data/{raw,processed}/other/orthology/`.
- `legacy/` — parked pre-reorganization material; exempt from the organism scheme — don't extend it.

Code/doc buckets: `scripts/` (numbered pipeline, below), `notebooks/` (exploration), `src/`
(reusable modules), `assets/` (static resources), `docs/` (the five prioritization-axis specs
`01_…`–`05_…` plus per-step logs).

### Pipeline (`scripts/`)

Numbered so the stage is obvious; the `0Nx` letter groups variants of a stage. Run with the
`gradi` env unless noted. Per-organism scripts take `--organism {kpneumoniae,ecoli}` (a few
also `human`); they default to `kpneumoniae`.

- `00a_fetch_proteomes.py` — fetch the kp/ecoli/human reference proteomes from UniProt (FASTA + TSV).
- `00b_proteome_descriptors.py` — 2×2 descriptor overview figure for the three proteomes.
- `01a_esmc_embeddings.py` — ESM-C 600M per-protein embeddings (NPZ).
- `01b_esmc_projections.py` — 2D openTSNE map of the embeddings (PNG; single faded stylia hue per organism).
- `01c_esmatlas_coords.py` — absolute Biohub ESM Atlas coordinates.
- `01d_esm_projections.py` — 2×2 plot comparing relative (openTSNE) vs absolute (Atlas) projections.
- `02a_interpro_annotation.py` / `02b_panther_annotation.py` — InterPro / PANTHER family annotation.
- `02c_family_plots.py` — family-overview barplots (kp + ecoli).
- `03a_orthology_general.py` — cross-species ortholog "synonym" table for any focal anchor (`--organism`;
  OrthoFinder + DIAMOND; needs the **`gradi-ortho`** env).
- `03b_general_orthology_plots.py` — per-anchor 2×2 overview of the 03a orthology mapping
  (`--organism`; writes `output/plots/03b_general_orthology_{kp,ec}.png`).
- `03c_orthology_focused.py` — focused 3-way orthology of kp × ecoli × human (OrthoFinder). Writes
  TABLES only (orthogroup-membership/Venn data, per-protein selectivity categories, broad-spectrum+
  human-selective shortlist, RBH %identity); plotting is done downstream from these tables.
- `03d_orthology_focused_plots.py` — one stylia slide (npg) from the 03c tables: UpSet of orthogroup
  membership, selectivity-category bars, RBH %identity distributions, proteome composition.
- `04a_alphafold_structures.py` — AlphaFold model availability / pLDDT / domain summary.
- `04b_alphafold_plots.py` — plots of the AlphaFold structural annotation.
- `04c_pdb_coverage.py` — experimental PDB structure coverage per protein (PDBe SIFTS).
- `04d_pdb_plots.py` — plots of the PDB structural-coverage annotation.
- `05a_popularity.py` / `05b_popularity_plots.py` — "studiedness" score transferred via orthology.
- `06a_chembl_bioactivity.py` / `06b_bindingdb_bioactivity.py` — ligandability §2.1b: # molecules
  tested and # potent (≤1 µM, pChEMBL/pAff ≥ 6) per protein, direct + ortholog-expanded, from
  local ChEMBL SQLite / BindingDB TSV dumps (needs the dumps under `data/raw/other/`).
- `06c_pdb_cocrystals.py` — §2.2a: drug-like PDB co-crystal ligands (direct + ortholog).
- `06d_alphafill_ligands.py` — §2.2b: AlphaFill transplanted-ligand evidence (alphafill.eu API).
- `06e_pockets.py` — §2.3b: fpocket + P2Rank pocket detection on the AlphaFold models, pLDDT-weighted
  (needs the **`gradi-pockets`** env + P2Rank; the script runs in `gradi`).
- `06f_af2bind.py` — §2.3a: AF2Bind binding-site prediction (scaffold; **deferred**, NaN placeholder).
- `06g_ligandability_merge.py` — §2.4 + composite: disorder filter, per-track sub-scores, and the
  final `ligandability_score` + `ligandability_tier` (`output/results/<org>/<prefix>_ligandability.csv`).
- `06h_ligandability_plots.py` — composite ligandability slide (`output/plots/06h_ligandability.png`).
- `06i_chembl_plots.py` / `06j_bindingdb_plots.py` / `06k_pdb_cocrystal_plots.py` /
  `06l_alphafill_plots.py` / `06m_pocket_plots.py` — per-resource slides (one **2×3 6-panel** figure
  per `--organism`, every panel single-organism — no content shared between the kp and ec slides).
  Stylia slide, NPG palette. Outputs `output/plots/06{i,j,k,l}_*_{kp,ec}.png` and
  `06m_pocket_{kp,ec}.png`. The structural slides (06k/06l) use the **strict drug-like** ligand tier
  (see `src/ligandability.py`: cofactors/nucleotides, amino acids, sugars, lipids/detergents,
  buffers/cryo/solvents excluded; the broad "any bound ligand" counts are retained as
  `*_n_ligand_any` / `*_has_ligand` columns). `06m` is the AlphaFold-structure druggability /
  binding-site poster (fpocket + P2Rank pockets, pLDDT-weighted `pocket_consensus_score` from 06e).
  Shared ligandability helpers live in `src/ligandability.py`.
- `06n_structure_snapshots.py` — ray-traced AlphaFold cartoon snapshots of the top druggable targets
  (coloured by per-residue pLDDT; top P2Rank pocket as green sticks/surface), 6 per organism →
  `output/plots/06n_structures_{kp,ec}.png`. Runs in `gradi` (target selection + montage) and shells
  out to the **`gradi-pymol`** env (PyMOL) for rendering via `scripts/_06n_pymol_render.py`.
- `06o_ligandability_landscape.py` — capstone synthesis (`output/plots/06o_landscape_{kp,ec}.png`):
  ligandability projected onto the ESM-C protein-universe map + the prioritization to the **prime**
  shortlist (broad-spectrum + human-selective + tractable), the evidence basis of the prime set, and
  a "neglected & druggable" view crossing ligandability with bibliometric studiedness (05a popularity).
- `07a–07k` — **essentiality** axis (docs §4). Emits a graded `essentiality_score` [0–1] +
  `essentiality_tier` per protein (`output/results/<org>/<prefix>_essentiality.csv` + `_shortlist.csv`).
  `07a` robust fetcher (Enterobacteriaceae-TraDIS compendium + open supp tables; ladder publisher-CDN →
  Europe-PMC-supp-zip → NCBI-OA-tarball → placeholder). Tracks: `07b` Kp Tn-seq/CRISPRi (ECL8 + KPPR1,
  gene-symbol mapped), `07c` E. coli EcoGene-essential transfer + graded broad-spectrum %essential,
  `07d` **ProteomeLM-Ess** (primary; backbone over the 01a ESM-C embeddings + a logistic head we train
  on E. coli labels, since the `-Ess` head is unreleased), `07e` **Geptop 2.0** reimplemented with
  DIAMOND, `07f` **FBA** single-gene deletion (iYL1228 kp / iML1515 ec, `cobra`), `07g` DeeplyEssential
  (deferred placeholder). `07h` merges into the graded composite (missing tracks renormalised, not
  zero-filled). `07i/07j/07k` stylia NPG slides (predictors · summary · cross-axis landscape).
  **`07l/07m` = the publication (experimental-only, prediction-free) view**: `07l` consolidates the Kp
  Mobile-CRISPRi-seq library + in-vivo (Jana 2023), KPNIH1/ECL8 Tn-seq, and the 12-genome
  Enterobacteriaceae-TraDIS cross-species essential matrix into `<prefix>_ess_publications.csv`; `07m`
  plots the dedicated CRISPRi/experimental slide. Jana 2023's ASM-gated tables were fetched via an
  authenticated Chrome session (chrome-devtools MCP `evaluate_script` same-origin fetch).
  **`07n/07o` give E. coli first-class parity** (E. coli is a target organism, not just a transfer
  reference): `07n_ecoli_experimental.py` ingests the major E. coli screens — Keio KO (PEC), Goodall
  TraDIS, Rousset 2018/2021 & Wang 2018 & Cui 2018 & Hawkins 2020 CRISPRi, and RB-TnSeq/Fitness Browser
  (280-condition antibiotic/stress matrix) — into `ec_ess_experimental.csv`; `07o_condition_plots.py`
  is the condition/stress slide (E. coli antibiotic sensitivity; Kp host-niche urine/serum/in-vivo).
  These feed E. coli's `evidence_experimental` (07h, 0.40 axis) and also transfer onto Kp via ortholog
  (07c). Gated E. coli sets (Goodall ASM, Hawkins Cell, Rousset 2021 Springer) were fetched via Chrome;
  Nichols 2011 (PMC reCAPTCHA) is the one un-fetched set. Shared helpers in `src/essentiality.py`
  (reuses `src/ligandability.py`; adds `jw_to_uniprot`, `gene_aliases_to_uniprot`, `transfer_ecoli_to_kp`).
  Run log: `docs/essentiality_log.md`.
  Needs `cobra`/`scikit-learn`/`openpyxl` + ProteomeLM from git (see `install.sh`); DIAMOND from
  `gradi-ortho`; `unar` for the Geptop `.rar`; optional `gradi-prokka` env.
- `08a_webapp_export.py` / `08b_validate_export.py` — build and validate the `app/data/{kp,ec}.json`
  payloads consumed by the webapp (the only generated artifacts tracked in Git, not eosvc).
- `09a`–`09h` — **localization** axis (docs §5.1). Emits `localization` + a graded
  `clp_accessibility` [0–1] per protein (`output/results/<org>/<prefix>_localization.csv` +
  `_shortlist.csv`), at **100% coverage** for both organisms (was 40% Kp / 51% Ec).
  Tracks: `09a` UniProt (curated + ECO experimental split + lipid anchors), `09b` **PSORTb 3.0 taken
  precomputed from PSORTdb** (per-genome download + RefSeq→UniProt id-mapping; we do *not* run
  PSORTb — the Docker image hangs for ~19 h under Rosetta), `09c` **DeepLocPro** (primary predictor;
  ESM-2 650M, beats PSORTb 3.0 on the post-2010 Gram-neg benchmark), `09d` **TMbed** topology
  (β-barrels + cytoplasm-facing residue fraction), `09e` SignalP 6.0 / lipobox Sec/SPII typing +
  Lol "+2 rule", `09f` **STEPdb 2.0** for Ec and its ortholog transfer onto Kp (this is what takes Kp
  from *1* experimentally-evidenced localization to ~1,680). `09g` merges them under an
  evidence-before-prediction precedence with topology overrides. **Four slides per organism:**
  `09h` axis summary · `09i` **compartment atlas** (small multiples over the ESM-C map — compartments
  occupy visibly distinct territory in sequence space) · `09j` evidence & predictor behaviour ·
  `09k` topology & envelope architecture. Compartment order/palette are shared via
  `LOC_CLASS_ORDER`/`LOC_CLASS_COLOR` in `src/localization.py` so they cannot drift between figures.
  Shared helpers in `src/localization.py`. `09c`/`09d` need the **`gradi-loc`** env (see *Setup*).
  Run log, calibration evidence and caveats: `docs/localization_log.md`. Access routes, assembly ids
  and identifier bridges: `docs/localization_downloads.md`. Two traps worth knowing before touching
  any axis: the 03a kp→ec ortholog table has **empty `pident`/`coverage`/`bitscore`** (OrthoFinder
  orthogroups only), so identity thresholds on transfer silently drop everything; and `fair-esm`
  collides with the EvolutionaryScale `esm` used by `01a`, hence the separate env.
- `10a`–`10c` — **degradability** axis (docs §3), partially built; see the dedicated section below for the
  protease decision and the two-bar design. `10a_fetch_degradability.py` robust fetcher (9-dataset manifest;
  ladder publisher-CDN → Europe-PMC → NCBI-OA → placeholder, plus `--stage` to ingest browser-obtained files);
  `10b_degrons.py` degron motifs + terminal pLDDT exposure/initiation-region length at full proteome coverage;
  `10c_clpp_activator.py` the **activated-ClpP (partnerless)** evidence track — Conlon 2013 (ADEP4) + Jacques
  2020 (ONC212) *S. aureus* proteomics, mapped on by NCBI-efetch → DIAMOND RBH (needs DIAMOND from
  `gradi-ortho`), writing `output/results/<org>/<prefix>_clpp_activator.csv`. Shared helpers in
  `src/degradability.py` (reuses `src/ligandability.py` and `src/essentiality.py`).
  **Degron slides (E. coli only so far):** `10k_degron_plots.py` — do sequence degrons work at all
  (prevalence, odds ratios *with intervals*, motif × accessibility, what drives `degron_score`, the
  ssrA-like set). `10o_degron_relevance.py` — do they work for *partnerless activated ClpP*, which is
  the question that decides the axis; answer: no, and the mechanism is why. Headlines: only **42 of
  4,403** proteins carry a weighted motif, ρ(exposure, `degron_score`) = **0.999** so the score is
  exposure not motif, **0 of the 11** ssrA-like proteins has an exposed C-terminus, and only **8** of
  the 609 activator-evidence proteins carry a weighted motif. Note `cterm_cm1_broad` holds weight
  0.40 on an OR whose 95% CI (0.76–11.99) **covers 1**. Run log: `docs/degradability_degron_log.md`.
  ⚠ The activated-ClpP evidence is ***S. aureus*** data (Conlon/Jacques) transferred by RBH at median
  42% identity — no native Gram-negative dataset exists — so `10o` labels every panel using it.
- `10e_measured_turnover.py` — the **measured** turnover/attribution layer: Nagar 2021 half-lives
  (materialises `D.load_nagar()`, which had only ever been used in memory), Gupta 2024 (13 conditions
  + the ΔclpP/Δlon/ΔhslV/triple/ΔsmpB panel) and Niwa 2022 (Lon vs ClpXP vs HslUV) →
  `<prefix>_deg_measured.csv`. Ec 74.7% / Kp 48.3% (transferred). **MacKrell 2026 deliberately not
  ingested** — no published half-life column, only raw timecourses; `sd04` is ML predictions.
  Slides `10l` (activator) · `10m` (turnover + attribution) · `10n` (cross-dataset landscape).
  Run log and findings: `docs/degradability_measured_log.md`. Headlines worth knowing: the two
  half-life datasets agree at only **ρ = 0.11**, ADEP4 abundance vs cleavage at **ρ = 0.06**, and
  growth correction leaves just **54** E. coli proteins with a positive proteolytic half-life.

### Degradability (docs §3) — partially built, stage `10*`

**The protease is settled: the recruiting handle is `ClpP`, engaged directly by a small-molecule activator,
with NO unfoldase partner.** The funded Gr-ADI research vision names it (WP2 is called *WP2_ClpPELs* —
ClpP-Engaging Ligands) and explicitly rules out ClpC for Gram-negatives. Consequences that govern the whole
axis:

- Activated ClpP **cannot unfold anything** — the activator occupies the ClpX/ClpA docking cleft, opening the
  pore, but nothing pulls. So the initiation-region rule (~5/~20/~37 aa) is an **unfoldase** rule and applies to
  ClpXP/ClpAP/Lon/HslUV/FtsH only. Under the Gr-ADI modality the substrate must **already be unstructured**:
  disorder, low stability, nascent chain, or a natively dynamic assembly. Disorder and
  `two_domain_architecture` are therefore first-class features, not minor ones.
- The proposal's WP1 criterion (c) **changed between drafts** (v3 "partnerless activated ClpP" → v5 "substrate
  of a ClpP protease *complex*") while the validation assay stayed partnerless. So the axis emits **two bars
  that are never averaged**: `clpP_complex_substrate` (selection, well populated) and `partnerless_clpP`
  (validation, sparse), plus `bar_disagreement`.
- ClpC/McsB are **absent** from both organisms, so every published (ClpC1-based) BacPROTAC describes a machine
  they do not have.

Built and run: `src/degradability.py` (helpers; corrected degron motifs with archetype self-tests;
growth-corrected turnover; per-paper trap weights with a ClpC cap), `10a_fetch_degradability.py` (9-dataset
manifest + fetch ladder + `--stage` for browser-obtained files), `10b_degrons.py` (motifs + terminal pLDDT
exposure, full proteome coverage), `10c_clpp_activator.py` (**track 3.3c** — the ADEP4/ONC212 activated-ClpP
evidence, transferred onto Kp/Ec by NCBI-efetch → DIAMOND RBH; 608/5,728 Kp and 609/4,403 Ec, ~11–14% coverage,
which is the honest ceiling). Not built: global disorder/architecture, turnover, biophysics, assembly state, the
compartment router, the merge, plots.

**The webapp still serves the OLD values — do not trust `comp_degradability` yet:**

- **Kp**: a frozen legacy file, `data/processed/legacy/klebsiella_pneumoniae_clp_degradability.tsv`, written by
  a since-deleted regex script (`scripts/03_annotate_clp_degradability.py`, removed in `8327de3`; the `03x` slot
  is now orthology). Its "experimental" inputs under `data/raw/legacy/clp_substrates/` are **hand-curated 45-
  and 35-protein substitutes**, not the papers' supplementary tables (see their `SOURCE.md`). Net signal: 5
  proteins with an ssrA-like C-terminus, 21 with trap evidence, 10 in the `high` tier. It has **four verified
  defects**, including an **inverted N-end rule** that supplies 680 of its 698 `medium` calls.
- **E. coli**: a **deterministic MD5 mock** (`08a_webapp_export.py:degradability_frame`), flagged
  `PROVISIONAL` in `app/config.js:37`. Not data.

Do not extend either; retire both when the merge lands.

Docs, in the order they are useful:
- **`docs/degradability_datasets.md`** — start here to build anything. The full dataset inventory (size, key,
  access, orthology-transfer verdict, effort), the computational-tool table, the PDF-extraction strategy, and the
  **column schemas for Ec and Kp separately** plus a 12-column minimum viable set. Headline: **one MobiDB bulk
  `curl` supplies `disorder_fraction` *and* the Pfam/Gene3D boundaries for `two_domain_architecture` at 100% of
  both proteomes** (⚠ its endpoint answers `405` to HEAD — use GET).
- `docs/03_degradability.md` — the spec; see the protease-decision block at the top and the 2026-08-06 composite
  revision (two bars, `bar_disagreement`, track 3.3c).
- `docs/degradability_report.md` — audit, Gr-ADI project context, the four-target check (**§12 supersedes parts
  of §1–§7**).
- `docs/degradability_references.md` — citations (§10 = project documents + the activator layer).
- `docs/degradability_downloads.md` — download routes (§A″ activator, §A‴ in-hand literature + route corrections).
- `data/raw/other/degradability/literature/` — **PDFs of three previously-gated papers**: Flynn 2003 and Neher 2006 (the E. coli ClpXP
  trap sets, ~60 + ~100 proteins) and Ziemski 2021 Table S1. Flynn's real N-M1 (`T-X-K-[ILV]`, 1–4
  residues in) and N-M3 consensuses are now **implemented** in `src/degradability.py` with archetype
  self-tests, at **weight 0 pending measurement** — see `degradability_datasets.md` §3.1. They are
  absent from the current `*_deg_degrons.csv` vintage; re-run `10b` to materialise them.

## Setup

Conda environments (commands documented in `install.sh`):

- **`gradi`** (Python 3.11) — the main env for everything except orthology:
  `conda create -y -n gradi python=3.11 && conda activate gradi && bash install.sh`
  (`install.sh` runs `pip install -r requirements.txt`: pandas, pyarrow, requests, biopython,
  torch, esm, openTSNE, umap-learn, matplotlib, colorcet, networkx, …).
- **`gradi-ortho`** (osx-64 bioconda; runs under Rosetta on Apple Silicon) — OrthoFinder + DIAMOND
  for the orthology scripts (`03a_orthology_general.py`, `03c_orthology_focused.py`) only.
  Created via micromamba; see the command block in `install.sh`.
- **`gradi-pockets`** (osx-64 bioconda; Rosetta on Apple Silicon) — `fpocket` + `openjdk=17` (JRE for
  P2Rank) for `06e_pockets.py` only. P2Rank itself is a standalone Java tarball under `tmp/tools/`.
  See `install.sh`.
- **`gradi-pymol`** (conda-forge `pymol-open-source`) — ray-traced AlphaFold cartoons for
  `06n_structure_snapshots.py` only (invoked from `gradi` via subprocess). See `install.sh`.
- **`gradi-loc`** (Python 3.11) — DeepLocPro + TMbed for the localization predictors
  (`09c_deeplocpro.py`, `09d_tmbed.py`) only. **This env must stay separate from `gradi`**:
  DeepLocPro needs `fair-esm`, which claims the same top-level `esm` package as the
  EvolutionaryScale `esm` that `01a_esmc_embeddings.py` relies on. Also pins `setuptools<81`
  (DeepLocPro imports `pkg_resources`) and `transformers==4.44.2` (TMbed's ProtT5 tokenizer breaks
  on transformers 5.x). See `install.sh`.

Do NOT use the machine's default `python3` (it resolves to an unrelated `ersilia` env). No build,
lint, or test commands are configured yet — document them here when added.
