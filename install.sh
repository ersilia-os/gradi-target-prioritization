#!/usr/bin/env bash
# Install dependencies for K. pneumoniae target annotation (v1).
# Assumes a Python >= 3.10 environment is already active (conda or venv).
# This project uses a dedicated conda env named 'gradi' (Python 3.11):
#     conda create -y -n gradi python=3.11 && conda activate gradi && bash install.sh
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
