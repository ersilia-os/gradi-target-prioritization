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
