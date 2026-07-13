#!/usr/bin/env bash
# Install dependencies for K. pneumoniae target annotation (v1).
# Assumes a Python >= 3.10 environment is already active (conda or venv).
# This project uses a dedicated conda env named 'gradi' (Python 3.11):
#     conda create -y -n gradi python=3.11 && conda activate gradi && bash install.sh
#
# Orthology (scripts/04_orthology.py) needs OrthoFinder + DIAMOND, installed in a
# SEPARATE bioconda env 'gradi-ortho' (osx-64; runs under Rosetta on Apple Silicon):
#     CONDA_SUBDIR=osx-64 conda create -y -n gradi-ortho -c bioconda -c conda-forge \
#         orthofinder diamond pandas pyarrow requests biopython tenacity python=3.11
#     conda activate gradi-ortho && python scripts/04_orthology.py
#
# Ligandability pocket detection (scripts/06e_pockets.py) needs fpocket + a JRE for P2Rank,
# in a SEPARATE bioconda env 'gradi-pockets' (osx-64; Rosetta on Apple Silicon):
#     CONDA_SUBDIR=osx-64 conda create -y -n gradi-pockets -c conda-forge -c bioconda \
#         fpocket openjdk=17 python=3.11
# P2Rank is a standalone Java tool (download the release tarball; gitignored under tmp/):
#     mkdir -p tmp/tools && cd tmp/tools \
#       && curl -sL -o p2rank.tar.gz https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz \
#       && tar xzf p2rank.tar.gz
# 06e finds these via env vars (defaults shown for this machine): FPOCKET_BIN
# (=~/miniconda3/envs/gradi-pockets/bin/fpocket), P2RANK_DIR (=tmp/tools/p2rank_2.5.1),
# POCKETS_JAVA_HOME (=~/miniconda3/envs/gradi-pockets/lib/jvm). The 06e script itself runs in `gradi`.
#
# Structure snapshots (scripts/06n_structure_snapshots.py) ray-trace AlphaFold cartoons with PyMOL,
# in a SEPARATE env 'gradi-pymol':
#     conda create -y -n gradi-pymol -c conda-forge pymol-open-source
# 06n runs in `gradi` (target selection + montage) and shells out to `gradi-pymol` for rendering
# (scripts/_06n_pymol_render.py).
#
# Ligandability bioactivity (scripts/06a, 06b) needs bulk dumps (eosvc/gitignored, NOT Git):
#     ChEMBL SQLite -> data/raw/other/chembl/   (ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_NN_sqlite.tar.gz)
#     BindingDB TSV -> data/raw/other/bindingdb/ (bindingdb.org/rwd/bind/downloads/BindingDB_All_YYYYMM_tsv.zip)
#
# Essentiality predictors (scripts/07*) run in `gradi`, but with two extra pieces:
#   - ProteomeLM backbone (07d, track 4.3a) is NOT on PyPI — install from git (Apache-2.0):
#         pip install "git+https://github.com/Bitbol-Lab/ProteomeLM.git"
#     (cobra / scikit-learn / openpyxl are in requirements.txt). 07d reuses the ESM-C 600M
#     embeddings from 01a as ProteomeLM's input; the -Ess head is trained locally on E. coli labels.
#   - ECL8 essentiality (07b, track 4.1a) needs a re-annotation of the ECL8 genome (the paper's
#     `ecl8_*` locus tags are a non-deposited Prokka annotation), in a SEPARATE bioconda env
#     'gradi-prokka' (osx-64; Rosetta on Apple Silicon):
#         CONDA_SUBDIR=osx-64 conda create -y -n gradi-prokka -c bioconda -c conda-forge prokka
#     07b shells out to it if present; otherwise it falls back to a gene-symbol bridge.
#   - Geptop (07e, track 4.3b) reuses DIAMOND + BLAST from `gradi-ortho`; its reference sets come
#     from the Geptop_v2.0.rar, extracted with `unar` (brew install unar).
#   - Geptop/FBA/strain-mapping (07b/07e/07f) all reuse the DIAMOND binary from `gradi-ortho`.
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
