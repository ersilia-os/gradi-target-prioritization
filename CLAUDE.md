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

Do NOT use the machine's default `python3` (it resolves to an unrelated `ersilia` env). No build,
lint, or test commands are configured yet — document them here when added.
