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
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
