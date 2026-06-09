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
- `04a_alphafold_structures.py` — AlphaFold model availability / pLDDT / domain summary.
- `04b_alphafold_plots.py` — plots of the AlphaFold structural annotation.
- `04c_pdb_coverage.py` — experimental PDB structure coverage per protein (PDBe SIFTS).
- `04d_pdb_plots.py` — plots of the PDB structural-coverage annotation.

## Setup

Two conda environments (commands documented in `install.sh`):

- **`gradi`** (Python 3.11) — the main env for everything except orthology:
  `conda create -y -n gradi python=3.11 && conda activate gradi && bash install.sh`
  (`install.sh` runs `pip install -r requirements.txt`: pandas, pyarrow, requests, biopython,
  torch, esm, openTSNE, umap-learn, matplotlib, colorcet, networkx, …).
- **`gradi-ortho`** (osx-64 bioconda; runs under Rosetta on Apple Silicon) — OrthoFinder + DIAMOND
  for the orthology scripts (`03a_orthology_general.py`, `03c_orthology_focused.py`) only.
  Created via micromamba; see the command block in `install.sh`.

Do NOT use the machine's default `python3` (it resolves to an unrelated `ersilia` env). No build,
lint, or test commands are configured yet — document them here when added.
