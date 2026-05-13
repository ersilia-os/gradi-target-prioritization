#!/usr/bin/env bash
# Install dependencies for K. pneumoniae target annotation (v1).
# Assumes a Python >= 3.10 environment is already active (conda or venv).
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
