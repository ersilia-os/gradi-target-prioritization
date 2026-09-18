#!/usr/bin/env bash
# Install dependencies for K. pneumoniae target annotation (v1).
# Assumes a Python >= 3.10 environment is already active (conda or venv).
# This project uses a dedicated conda env named 'gradi' (Python 3.11):
#     conda create -y -n gradi python=3.11 && conda activate gradi && bash install.sh
#
# Orthology (scripts/orthology/orthofinder.py) needs OrthoFinder + DIAMOND, installed in a
# SEPARATE bioconda env 'gradi-ortho' (osx-64; runs under Rosetta on Apple Silicon):
#     CONDA_SUBDIR=osx-64 conda create -y -n gradi-ortho -c bioconda -c conda-forge \
#         orthofinder diamond pandas pyarrow requests biopython tenacity python=3.11
#     python scripts/orthology/orthofinder.py        # runs in `gradi`, shells out to gradi-ortho
#   Do NOT call the orthofinder launcher directly: its `#!/usr/bin/env python3` shebang
#   resolves to an unrelated env. The stage names <env>/bin/python explicitly.
#
# Ligandability pocket detection (legacy/scripts/06e_pockets.py (v1 only -- see scripts/pockets/README.md)) needs fpocket + a JRE for P2Rank,
# in a SEPARATE bioconda env 'gradi-pockets' (osx-64; Rosetta on Apple Silicon):
#     CONDA_SUBDIR=osx-64 conda create -y -n gradi-pockets -c conda-forge -c bioconda \
#         fpocket openjdk=17 python=3.11
# P2Rank is a standalone Java tool (download the release tarball; gitignored under tmp/):
#     mkdir -p tmp/tools && cd tmp/tools \
#       && curl -sL -o p2rank.tar.gz https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz \
#       && tar xzf p2rank.tar.gz
# that v1 script finds these via env vars (defaults shown for this machine): FPOCKET_BIN
# (=~/miniconda3/envs/gradi-pockets/bin/fpocket), P2RANK_DIR (=tmp/tools/p2rank_2.5.1),
# POCKETS_JAVA_HOME (=~/miniconda3/envs/gradi-pockets/lib/jvm). It runs from inside `legacy/`.
#
# Structure snapshots (legacy/scripts/06n_structure_snapshots.py) ray-trace AlphaFold cartoons with PyMOL,
# in a SEPARATE env 'gradi-pymol':
#     conda create -y -n gradi-pymol -c conda-forge pymol-open-source
# it runs from inside `legacy/` (target selection + montage) and shells out to `gradi-pymol` for rendering
# (legacy/scripts/_06n_pymol_render.py).
#
# Stage 03 localization (scripts/localization/predict.py) needs a SEPARATE env 'gradi-loc' for both its
# predictors. This split is MANDATORY, not cosmetic: DeepLocPro depends on `fair-esm`, which installs
# the same top-level `esm` package as the EvolutionaryScale `esm` that `gradi` uses for the stage-01
# ESM-C embeddings. Installing DeepLocPro into `gradi` silently breaks stage 01.
#     conda create -y -n gradi-loc python=3.11
#     conda activate gradi-loc
#     pip install "torch>=2.0" pandas numpy h5py sentencepiece tqdm typer biopython "setuptools<81"
#     pip install "transformers==4.44.2" "tokenizers<0.20"
#     pip install git+https://github.com/BernhoferM/TMbed.git
#     pip install git+https://github.com/Jaimomar99/deeplocpro.git
# Two pins that are load-bearing:
#   - `setuptools<81` — DeepLocPro imports `pkg_resources`, removed from newer setuptools.
#   - `transformers==4.44.2` — TMbed loads ProtT5 through `T5Tokenizer`; transformers 5.x routes
#     that through the tiktoken converter and dies with "`tiktoken` is required to read a
#     `tiktoken` file".
# DeepLocPro downloads ESM-2 650M and TMbed downloads ProtT5-XL-U50 (~2.25 GB) on first run.
#
# Do NOT activate this env to run stage 03. The stage runs in `gradi` and reaches across a process
# boundary: `scripts/localization/workers/deeplocpro.py` under the gradi-loc interpreter (stdlib + torch only), plus
# the `tmbed` console script, whose shebang already points there. GRADI_LOC_BIN overrides where the
# stage looks. DeepLocPro's own `-d mps` flag is broken (`embed_batch()` gates device placement on
# `torch.cuda.is_available()`), so the stage drives `EnsembleModel` itself — measured ~5.6 prot/s on
# MPS. TMbed gates GPU on CUDA too, so on macOS it runs ProtT5 on CPU: expect hours, sharded and
# resumable. Run the two tracks SEQUENTIALLY — contending for the CPU inflated v1's TMbed shard
# time from ~4 min to ~17 min.
#
# Ligandability bioactivity (scripts/ligands/chembl.py, bindingdb.py) needs bulk dumps (eosvc/gitignored, NOT Git):
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
#
# TabPFN-3.5 is the project's DEFAULT supervised learner (see CLAUDE.md, "Supervised ML"). It runs
# in a SEPARATE env 'gradi-tabpfn'. This split is MANDATORY: tabpfn pulls torch 2.14 + mlx against
# `gradi`'s torch 2.12, and installing it there would break stage 01 (ESM-C) -- the same trap as
# blast, eggnog-mapper and rdkit.
#
#     conda create -y -n gradi-tabpfn python=3.11
#     ~/miniconda3/envs/gradi-tabpfn/bin/pip install "tabpfn==9.0.0" "tabpfn-client==0.6.0"
#
# Needs TABPFN_TOKEN in the environment (https://ux.priorlabs.ai -> account). Two things to know:
#   - LOCAL weights additionally need the licence ACCEPTED on the Licenses tab. A valid token is NOT
#     enough -- verified: token valid, `accepted: False`, download refused. Until then only the
#     hosted route works (`GRADI_TABPFN_HOSTED=1`, the default), which UPLOADS features and labels
#     to Prior Labs and bills ~10,000 credits per call against a 20M monthly quota.
#   - The weights are NON-COMMERCIAL licensed. Fine for methods work; resolve before shipping.
#
# Stage 04 runs in `gradi` and shells out to scripts/workers/tabpfn_cv.py under this interpreter
# (GRADI_TABPFN_BIN overrides). Do NOT activate this env to run a stage.
#
# The lazy-qsar estimator comparison (scripts/degradability/head_comparison.py, --estimator lazyqsar)
# needs a SEPARATE env 'gradi-lazyqsar'. This split is MANDATORY, not cosmetic: lazy-qsar pins
# numpy==2.1.3 and scikit-learn==1.6.1, while `gradi` runs numpy 2.4.6 / scikit-learn 1.9.0 under
# torch 2.12. Installing it into `gradi` downgrades both underneath torch and breaks stage 01
# (ESM-C) and stage 04 (the forest) together -- the same trap as blast, eggnog-mapper and rdkit.
#
#     conda create -y -n gradi-lazyqsar python=3.11
#     git clone --branch v3.4.4 https://github.com/ersilia-os/lazy-qsar
#     ~/miniconda3/envs/gradi-lazyqsar/bin/pip install './lazy-qsar[fit]'
#
# The `[fit]` extra is required -- the base install is inference-only (onnxruntime, no sklearn), so
# `LazyClassifier.fit` raises on a bare install. Pin the TAG, not `main`: the class selects its own
# lr/xgb/rf/svc portfolio per fit, so an unpinned upgrade can silently change which models produced
# a published number.
#
# Do NOT activate this env to run the comparison. It runs in `gradi` and reaches across a process
# boundary to scripts/degradability/workers/lazyqsar_cv.py; GRADI_LAZYQSAR_BIN overrides where it looks.
# The worker takes its FOLDS FROM THE CALLER and never uses lazy-qsar's own `oof_auc_`, which comes
# from ungrouped internal splits that put sequence homologs on both sides -- measured +0.014
# optimistic on a one-feature matrix, and the reason the folds are computed in `gradi` by stage 04's
# own splitters and passed in.
#
# Stage 02 GO slim (scripts/function/eggnog.py) needs eggNOG-mapper in a SEPARATE env
# 'gradi-emapper' (osx-64; Rosetta on Apple Silicon -- there is no osx-arm64 build, and installing
# it into `gradi` would drag that whole env to osx-64 and take ESM-C down with it):
#
#     CONDA_SUBDIR=osx-64 conda create -y -n gradi-emapper -c conda-forge -c bioconda \
#         eggnog-mapper=2.1.15
#
#   Two traps, both hit on the first run:
#
#   1. bioconda installs diamond/mmseqs into the env's bin/, but emapper looks for them inside
#      site-packages/eggnogmapper/bin/. Link them:
#        E=$CONDA_PREFIX; B=$E/lib/python3.11/site-packages/eggnogmapper/bin; mkdir -p $B
#        for t in diamond mmseqs hmmsearch hmmscan phmmer; do ln -sf $E/bin/$t $B/$t; done
#
#   2. `download_eggnog_data.py` DOES NOT WORK. It fetches from eggnogdb.embl.de, which no longer
#      resolves, and then prints "Finished." with exit status 0 having downloaded nothing. The live
#      host is eggnog5.embl.de. Fetch the two files yourself, RESUMABLY (the server drops long
#      transfers), and verify each against its Content-Length before decompressing:
#        http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz            -> eggnog.db
#        http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz -> eggnog_proteins.dmnd
#        http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz      -> eggnog.taxa.db
#      into data/raw/function/eggnog/ (~21 GB unpacked). Public and re-derivable: do NOT eosvc it.
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
