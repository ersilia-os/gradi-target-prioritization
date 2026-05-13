#!/usr/bin/env python3
"""Build output/results/kp_target_annotation_v1.{csv,parquet}.

Usage (from repo root):
    python scripts/build_annotation.py

Inputs (already staged or downloadable per data/raw/.../SOURCE.md):
    data/raw/klebsiella_pneumoniae_proteome.tsv     (UniProt HS11286)
    data/raw/escherichia_coli_proteome.tsv          (UniProt MG1655)
    data/raw/bvbrc/hs11286_features.tsv             (BV-BRC PATtyFam)
    data/raw/essentiality/literature/eichelberger2024_ECL8/*.xlsx
    data/raw/essentiality/literature/zhu2023_kp_crispri/curated_highlights.tsv
    data/raw/essentiality/literature/bachman2015_KPPR1/*.xlsx    (manual stage)
    data/raw/essentiality/literature/ramage2017_KPNIH1/*.xlsx    (manual stage)
    data/raw/essentiality/literature/goodall2018_Ec_BW25113/*.xlsx (manual stage)
"""

from pathlib import Path
import sys

# Allow running from any cwd
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.assemble import main

if __name__ == "__main__":
    main()
