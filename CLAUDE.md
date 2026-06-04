# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status

This repo is currently a freshly-instantiated copy of the Ersilia analysis template (`eos-analysis-template`). There is no project code, no dependencies declared in `requirements.txt`, and no logic in `install.sh` — these are empty placeholders awaiting the actual project work. `scripts/`, `notebooks/`, and `src/` (not yet created) are expected destinations for new code. The repo name implies the eventual analysis will be a **target prioritization** project for the **GraDi** collaboration.

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

The single *K. pneumoniae* anchor strain is **HS11286** (UniProt proteome `UP000007841`; NCBI `GCF_000240185.1`; locus-tag prefix `KPHS_*`; 5,728 proteins). It is the **only** *K. pneumoniae* proteome UniProt flags as a "Reference and representative proteome", so it has the most complete annotation and the cleanest cross-references. The FASTA lives at `data/raw/proteome/UP000007841_HS11286.fasta` (fetch with `scripts/fetch_proteome_hs11286.py`).

Why not the alternatives:
- **ATCC 43816 / KPPR1** (the originally suggested strain) is popular only in mouse *in-vivo* essentiality work, not annotation — it has no curated UniProt reference proteome. Avoid it as the anchor.
- **MGH 78578** is the *historical* reference and holds most of the species' reviewed (Swiss-Prot) entries, but reviewed coverage is not a meaningful strain criterion here: across *all* K. pneumoniae strains there are only ~960 reviewed vs ~187k unreviewed entries (~0.5%), and reviewed annotation propagates across strains via orthology anyway.

**Mapping rule:** resolve all orthology and identifier mappings **onto HS11286** whenever a corresponding ortholog exists. Evidence and data from other strains (e.g. KPPR1, KPNIH1, ECL8 essentiality; *E. coli* K-12 orthology) should be carried as annotations on the HS11286 protein (keyed by its UniProt accession), preserving the source strain/ID as provenance columns rather than as a separate primary key. Note `KPHS_*` is a locus tag, not a UniProt accession — keep both, but the UniProt accession is canonical (see *Identifier convention*).

Resources keyed at the **species** level (e.g. ChEMBL, taxid 573) are strain-agnostic; the anchor-strain choice does not affect them.

## Directory contract

The template prescribes specific roles — keep new files in the right bucket:

- `data/raw/` — original, untouched inputs. `data/processed/` — cleaned/transformed derivatives.
- `output/results/` — numerical results, logs, text outputs. `output/plots/` — figures.
- `scripts/` — standalone preprocessing/automation scripts.
- `notebooks/` — Jupyter exploration and prototyping.
- `src/` — reusable modules (directory not yet created; create when needed).
- `assets/` — static resources (images, figures used in docs).

## Setup

`install.sh` and `requirements.txt` are intentionally empty; populate them as dependencies are introduced. No build, lint, or test commands are configured yet — when adding them, document them here.
