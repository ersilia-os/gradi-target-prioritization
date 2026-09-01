# Ligandability assessment log (docs §2)

Implemented `scripts/06a–06h` + `src/ligandability.py`. Per-protein ligandability for
K. pneumoniae HS11286 (5,728) and E. coli K-12 (4,403), keyed by UniProt accession.

## Data versions
- **ChEMBL 37** SQLite dump (29 May 2026) — `data/raw/other/chembl/` (11,481 protein targets with B/F activity).
- **BindingDB** All TSV, release 2025-04 — `data/raw/other/bindingdb/`.
- **pdb_seqres** (RCSB derived_data, full) — `data/raw/other/pdb/` (1,069,489 protein chains).
- **AlphaFold** models v6 (from 04a) + **AlphaFill** API (alphafill_version 2.1.1).
- Pocket tools: **fpocket 4.0**, **P2Rank 2.5.1** (`-c alphafold`).

## Key design decision — sequence mapping, not accession
HS11286 is a dark TrEMBL proteome; its accessions ≠ the reviewed/other-strain accessions used by
ChEMBL / BindingDB / PDB-SIFTS. Exact-accession matching missed *direct* K. pneumoniae data
(notably the SHV/OXA-48/CTX-M/NDM/AmpC β-lactamases). Fixed by mapping proteins to external target
**sequences** with DIAMOND (≥95% id = "direct"). Validated: HS11286 `A0A0H3H184` → ChEMBL Q93LQ9
(*K. pneumoniae* β-lactamase, 100% id, 104 potent compounds) + PDB 5eec/6d15… inhibitor co-crystals
— all previously invisible to accession matching.

**Refinement (2026-06-25):** the bacterial bucket must be *true Bacteria*, not "non-human", and
transfer needs an identity floor. Without this, heavily-screened rat/mouse/eel/parasite targets at
~30% id inflated counts (a first pass gave 424 kp ChEMBL-potent proteins, with rat P97697 as the
"best" hit for several proteins). Fixed: restrict to Bacteria via ChEMBL `organism_class`
(BindingDB bridged by bacterial genus names) + a 40% identity floor for transfer. Corrected kp
ChEMBL-potent = 175, BindingDB = 93 — top targets now all real (gyrA/parE/gyrB/recA/def from
E. coli, P. aeruginosa, M. tuberculosis, S. aureus).

## Coverage (kp / ec)
- ChEMBL ≤1 µM: 175 / 155 proteins.   BindingDB ≤1 µM: 93 / ~.
- PDB drug-like co-crystal (any): 2,525 / 2,157.   AlphaFill drug-like transplant: 2,993 / 2,460.
- Druggable pocket (`evidence_pocket ≥ 0.5`): 2,032 / 1,706.
- Tiers — kp: 3,335 tractable / 884 partial / 1,509 intractable; ec: 2,787 / 619 / 997.
- Hard evidence (≤1 µM activity or own/≥95%-id co-crystal): 429 kp / 1,240 ec.
- Prime (broad-spectrum, human-selective) kp: 3,015 → 1,948 tractable, **223 with hard evidence**.

## Spot-checks (kp)
Top by score: pyrD (DHODH), lpxH, leuS, map, KPHS_07570 — all real antibacterial targets.
gyrA 0.77 / gyrB 0.91 / β-lactamase 0.84 — all tractable.
Prime shortlist (`kp_ligandability_shortlist.csv`, 1,949 targets, 306 with hard evidence): top =
lpxH, aroA, murC, mtnN, coaA, menB, dxr (with fosmidomycin/FOM via AlphaFill) — textbook
broad-spectrum, human-selective, ligandable antibacterial targets.

## Reproduce
```
# dumps (see install.sh) -> data/raw/other/{chembl,bindingdb,pdb}; envs: gradi, gradi-ortho, gradi-pockets
conda run -n gradi python scripts/06a_chembl_bioactivity.py   --organism kpneumoniae --refresh
conda run -n gradi python scripts/06b_bindingdb_bioactivity.py --organism kpneumoniae --refresh
conda run -n gradi python scripts/06c_pdb_cocrystals.py        --organism kpneumoniae
conda run -n gradi python scripts/06d_alphafill_ligands.py     --organism kpneumoniae
conda run -n gradi python scripts/06e_pockets.py               --organism kpneumoniae
conda run -n gradi python scripts/06f_af2bind.py               --organism kpneumoniae   # deferred placeholder
conda run -n gradi python scripts/06g_ligandability_merge.py   --organism kpneumoniae
conda run -n gradi python scripts/06h_ligandability_plots.py
# repeat 06a-06g with --organism ecoli; outputs in output/results/<org>/<prefix>_ligandability.csv
```

## Deferred
AF2Bind (2.3a) — scaffolded, needs GPU. PocketMiner/CryptoBank, FTMap, PASSer, canSAR/DoGSiteScorer,
BioLiP2 — see suggestions in `02_ligandability.md`.
