# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Maintaining this file

**A pipeline step is not done until this file describes it.** In the same commit as the code, add or
update: the stage entry (what it does, what it writes, which conda env it needs), any new trap worth
warning about, and any new environment requirement. This file is the map; a stale map is worse than
none.

Do not let it grow into an inventory of everything ever tried — that is what `legacy/HISTORY.md` is
for. Keep it to what someone needs to work here *today*.

## House style

**Scripts print generously.** Every script opens with a banner naming itself, its inputs and its
outputs; reports per-step progress with running counts; and closes with a summary/coverage table
carrying percentages. Nothing important happens silently. Offer `-q/--quiet` for non-interactive use.

Two corollaries, both learned the hard way in v1:

- **Verify counts against an authoritative source before trusting a payload.** An HTTP 200 is not
  evidence of data — v1 was repeatedly bitten by short or empty responses (a 61-byte 200, an empty
  202, a 200 from the wrong article). Assert, and fail loudly.
- **Never resolve ambiguity silently.** Emit the chosen value, its source, and the alternatives, plus
  an audit row. If a rule decided between candidates, record which rule.

## Status

Active **target-prioritization** analysis for the **GraDi** collaboration (BacPROTAC-style targeted
protein degradation in Gram-negatives; Prof. Erick Strauss, Stellenbosch). Ersilia's deliverable is a
prioritized list of proteins of interest. Ligand identification is out of scope.

**This is v2, a deliberate restart.** v1 is complete and frozen under `legacy/` — see the *Legacy*
section below before assuming anything about prior work.

## Species

Four, **one reference proteome each**. That proteome is the unit of analysis; nothing else is.

| species | proteome | strain | n | locus tags |
|---|---|---|---|---|
| *K. pneumoniae* | `UP000007841` | HS11286 | 5,728 | `KPHS_*` |
| *E. coli* | `UP000000625` | K-12 MG1655 | 4,403 | b-numbers |
| *S. aureus* | `UP000008816` | NCTC 8325 | 2,889 | `SAOUHSC_*` |
| *H. sapiens* | `UP000005640` | — | 20,416 | — |

**HS11286** is the anchor: the only *K. pneumoniae* proteome UniProt flags "Reference and
representative". Rejected as anchor — **ATCC 43816 / KPPR1** (no curated UniProt reference proteome)
and **MGH 78578** (historical reference; reviewed coverage is ~0.5% species-wide and propagates by
orthology anyway).

***S. aureus* is new in v2.** It is the organism the activated-ClpP proteomics (Conlon 2013 ADEP4,
Jacques 2020 ONC212) is native to — v1 could only reach that data by cross-phylum DIAMOND RBH at
median 42% identity — and the only one of the three bacteria with ClpC/McsB, the machine every
published BacPROTAC actually targets. Whether it becomes a full target organism or stays a reference
is deferred.

**Human must be fetched `reviewed:true`** — the unfiltered proteome is 147,506 TrEMBL-bloated entries.

## Identifier convention

**UniProt accession is the canonical key** in every dataset, `data/` and `output/` alike. Gene
symbols, locus tags, b-numbers, RefSeq and GeneID are retained as provenance columns — never as the
primary key. `KPHS_*` is a locus tag, not an accession.

**But join on `locus_tag`, not on gene name.** `locus_tag` is at 100% coverage while `gene_name` is at
18.4% on Kp, and essentially every published bacterial dataset (Tn-seq, TraDIS, CRISPRi, proteomics)
keys on locus tags. v1's single most expensive recurring bug was a gene-symbol join that silently
lost ~27% of known essentials.

**Map by sequence, not by accession, when reaching an external database.** HS11286 is a dark TrEMBL
proteome whose accessions rarely match those used by ChEMBL, BindingDB or PDB-SIFTS. v1's convention,
which worked: DIAMOND blastp, **≥95% identity = "direct"**, **≥40% floor for transfer**, and restrict
the bacterial bucket to *true Bacteria* rather than "non-human".

## Two-track persistence: Git vs eosvc

Enforced by `.gitignore`:

- **Git**: `src/`, `scripts/`, `docs/`, `assets/`, `legacy/` (code and docs), `LICENSE`, `README.md`,
  `install.sh`, `requirements.txt`, `access.json`.
- **eosvc (DVC + S3), not Git**: `data/` and `output/`. Use
  `~/miniconda3/bin/eosvc {view,upload,download,delete} --path <relative path>`.
- `tmp/` is local-only scratch, tracked by neither.
- The only generated artifacts in Git are `legacy/app/data/*.json`, the webapp payloads (an explicit
  `.gitignore` exception).

`eosvc upload` PUTs **one object at a time (~2 files/sec)**, so a directory of many small files is
hours, not minutes. Archive such trees before pushing — this is why `data/processed/legacy/v1/` holds
tar.gz files. And `eosvc view --path <file>` is unreliable for a single file; read the parent
directory listing's status column instead.

## Directory contract

`data/{raw,processed}` and `output/{results,plots}`, bucketed **by stage number**:

```
data/raw/00_proteomes/{uniprot,ncbi}/     as fetched, + <label>.SOURCE.md
data/processed/00_proteomes/              the stage output
data/raw/PROVENANCE.md                    what every raw dataset is, and if you could get it again
```

v1 used an organism-first layout (`data/processed/kpneumoniae/...`); **do not extend it.** Its
remaining contents are v1 source data, indexed by `data/raw/PROVENANCE.md`.

## Pipeline (`scripts/`)

Numbered so the stage is obvious. Run with the `gradi` env.

- **`00_download_proteomes.py`** — fetch the four reference proteomes and build **one identified table
  per species**. Registry-driven from `src/proteome_registry.tsv` (41 rows, 4 tiers, every row carrying
  a required `why`).

  Its real job is fixing the naming gap: gene names cover only **18.4% of Kp** and 28.2% of
  *S. aureus* (vs 100% for E. coli), and names are what literature and databases key on. Filled in
  three labelled tiers, all **within the species** — `anchor` → `species_exact` (identical sequence)
  → `species_uniref90` (same UniRef90 cluster, via UniProt ID-mapping, no DIAMOND needed). Measured:
  **Kp 18.4% → 63.4%**, Sa 28.2% → 44.6%. Cross-species naming is deliberately deferred to a stage
  that has real orthology and can label it as the inference it is.

  `gene_name` is always **exactly one preferred name**, chosen by a total order (reviewed donor →
  real name over a `y###` placeholder → most attested → shortest, then alphabetical), so a tie can
  never leave it empty. Alternatives are kept in preference order and every fill is audited.

  Outputs `data/processed/00_proteomes/`: `<species>_identity.tsv` (the deliverable),
  `<species>_annotation.tsv`, `proteins.parquet`, `id_bridge.tsv`, `name_audit.tsv`, `manifest.tsv`,
  `registry.tsv`. CLI: `--tier A,B,C,D` · `--only LABEL` · `--refresh` · `--dry-run` · `-q`.
  Details, traps and the run log: `docs/00_proteomes.md`.

Registry tiers: **A** the 4 anchors · **B** same-species name-donor pools · **C** a 26-species
comparator panel for orthology (v1's curated panel, with every species now pinned to an explicit
proteome id) · **D** evidence-bridge strains from **NCBI**, because **UniProt serves zero entries for
any proteome it does not flag "Reference proteome"** — verified for ECL8, KPNIH1, KPPR1, NJST258,
BW25113, COL, Mu50. C and D are fetched only on request.

Two tier-D traps worth knowing: **S. aureus COL must come from GenBank, not RefSeq** (`GCF_` retains
60 original `SACOL####` tags, `GCA_` retains 2,711), and **ECL8's GFF carries
`old_locus_tag=BN373_…,KPNEcl8_…`**, which should make v1's `gradi-prokka` re-annotation unnecessary.

## Legacy

**`legacy/` holds the complete v1 pipeline, frozen. Do not extend it.** Start at
**`legacy/HISTORY.md`** — the retrospective: what was built and what was not, the methodological
decisions and why, ~70 documented traps, which data sources were obtained and how, and which
artifacts look like data but are not.

Three things there must not be trusted (`legacy/HISTORY.md` §7): the **Kp legacy degradability TSV**
(four verified defects, including an inverted N-end rule supplying 680 of 698 `medium` calls), the
**E. coli degradability values in the webapp** (a deterministic MD5 mock), and
**`data/raw/legacy/clp_substrates/`** (45- and 35-row hand-curated substitutes, not the papers'
tables).

The archived scripts still run. They resolve paths via `parents[1]`, so `legacy/data` and
`legacy/output` are **symlinks** to the real trees one level up — that is what keeps them working
without editing 78 files. Without those symlinks a script exits 0 having written nothing. Run from
inside `legacy/`:

```bash
cd legacy && python scripts/09h_localization_plots.py --organism ecoli
```

The v1 webapp is still deployed from `legacy/app/` by `.github/workflows/pages.yml`.

## Setup

Conda environments (see `install.sh`):

- **`gradi`** (Python 3.11) — the main env; `bash install.sh` runs `pip install -r requirements.txt`.
  **Stage 00 needs nothing beyond `requests`, `pandas`, `tenacity` and `pyarrow`.**
- **`gradi-ortho`** (osx-64 bioconda, Rosetta on Apple Silicon) — OrthoFinder + DIAMOND, which have no
  arm64 build. Its DIAMOND binary was borrowed by many v1 scripts via `GRADI_DIAMOND_BIN`.
- **`gradi-pockets`** (osx-64, Rosetta) — `fpocket` + `openjdk=17` for P2Rank.
- **`gradi-pymol`** — `pymol-open-source`, for ray-traced structure cartoons.
- **`gradi-loc`** (Python 3.11) — DeepLocPro + TMbed. **This split is mandatory, not cosmetic**:
  DeepLocPro needs `fair-esm`, which claims the same top-level `esm` package as the EvolutionaryScale
  `esm` used for ESM-C embeddings. Installing it into `gradi` silently breaks embedding generation.
  Pins `setuptools<81` and `transformers==4.44.2`.

Do **not** use the machine's default `python3` — it resolves to an unrelated `ersilia` env. Use
`~/miniconda3/envs/gradi/bin/python`.

No build, lint or test commands are configured yet. Document them here when added.
