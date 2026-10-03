# Structural ligandability (`pockets/`)

**Can a small molecule bind this protein's fold?** The bioactivity track (`docs/ligands.md`) asks
whether somebody *measured* a potent compound; this asks whether the *structure* says one could
bind. Two answers that share no input and are never merged:

- **Predicted** — pockets on the AlphaFold model, by fpocket 4.0 and P2Rank 2.5.1.
- **Measured** — the closest bacterial PDB chain seen with a drug-like ligand in the same site.

Deliverable `data/processed/pockets/pockets_<species>.tsv`, complete and canonical for the three
bacteria (`python -m src.matrices`):

| column | kind | meaning |
|---|---|---|
| `uniprot_ac` | key | |
| `p2rank_score` | predicted, 0–1 | P2Rank (`-c alphafold`) calibrated probability of the best admitted pocket |
| `fpocket_score` | predicted, 0–1 | fpocket druggability score of the best admitted pocket |
| `p2rank_n_pockets` | int | admitted P2Rank pockets at probability ≥ 0.5 |
| `holo_identity` | measured, 0–100 | % identity to the closest bacterial chain with a drug-like ligand in the aligned site; 0 = none |
| `pdb_n_structures` | measured, int | PDB entries with a chain that IS this protein (≥ 95% identity, ≥ 50% of the chain aligned), ligand or not; 0 = none |
| `pdb_coverage` | measured, 0–1 | fraction of the protein's residues covered by those chains (SEQRES) |
| `af_plddt` | confidence | mean AlphaFold pLDDT; **NA means no model** |

### There is no `evidence` column

Dropped on the project owner's instruction, **2026-10-03**. It held
`pdb+af` · `af_only` · `pdb_only` · `none`, and every one of those is a function of two columns
that are still here — verified exactly reconstructible on all three species before removal:

| old label | = |
|---|---|
| `pdb+af` | `af_plddt.notna()` and `holo_identity > 0` |
| `af_only` | `af_plddt.notna()` and `holo_identity == 0` |
| `pdb_only` | `af_plddt.isna()` and `holo_identity > 0` |
| `none` | `af_plddt.isna()` and `holo_identity == 0` |

**The one thing it added that these columns cannot say** is where the model came from — AlphaFold
DB or ESMFold, for the proteins AFDB does not cover. That is `model_source` in
`evidence/alphafold_<species>.tsv`, via `load_alphafold()`. An ESMFold row's `af_plddt` is
ESMFold's pLDDT on the same 0–100 scale, so the column stays comparable either way.

Load through `src/pockets.py` (`load`, `load_all`, `druggability`, `load_alphafold`,
`load_pockets`, `load_pdb`, `load_holo`, `load_ligand_classes`, `load_alphafill_comparison`,
`druglike`).

## Pipeline

```
structures.py      AlphaFold DB v6, validated by sequence   -> evidence/alphafold_<sp>.tsv
esmfold.py         ESMFold for proteins AFDB lacks          -> scratch/esmfold/, evidence/esmfold_<sp>.tsv
structures.py      (again, to pick the ESMFold models up)
predict.py         fpocket + P2Rank (gradi-pockets)         -> evidence/pocket_list_<sp>.tsv
holo.py            ligand classes + DIAMOND vs BioLiP       -> evidence/{ligand_classes,holo_<sp>}.tsv
pdb_coverage.py    DIAMOND vs every PDB protein chain       -> evidence/pdb_<sp>.tsv
alphafill_check.py measurement only                         -> evidence/alphafill_comparison.tsv
merge.py           the deliverable + checks                 -> pockets_<sp>.tsv
```

All run with `gradi`. fpocket/P2Rank live in `gradi-pockets` (osx-64, Rosetta; `FPOCKET_BIN`,
`P2RANK_DIR`, `POCKETS_JAVA_HOME` override) and DIAMOND in `gradi-ortho` (`GRADI_DIAMOND_BIN`).
Cold run ~40 min (pockets ~35, holo ~3); every step caches.

## Structures

AlphaFold DB is still v6 (checked 2026-10-03), so v1's cached Kp/Ec models were seeded into
`data/source/alphafold/` and re-validated; Sa was downloaded. **A model is used only if its sequence
equals the proteome's exactly** — 0 mismatches in all three species.

| species | models | no model | median pLDDT | % proteins ≥ 70 |
|---|---|---|---|---|
| Kp | 5,727 / 5,728 | 1 | 91.1 | 93.3 |
| Ec | 4,371 / 4,403 | 32 | 91.3 | 96.6 |
| Sa | 2,888 / 2,889 | 1 | 91.4 | 93.6 |

### The 34 proteins AlphaFold DB does not model, and ESMFold

All 34 are AFDB 404s, and **no other accession rescues any of them**: UniParc groups identical
sequences across accessions, and none of the 34 has an AFDB model under ANY accession sharing its
exact sequence, even where hundreds do (pheM 1,062, fdhF 126). They are AFDB's exclusion rules:

| why | proteins |
|---|---|
| shorter than 16 aa | 21 E. coli micro-peptides (8–15 aa: leader peptides, small proteins) |
| selenocysteine (U) | E. coli fdhF, fdnG, fdoG |
| in-frame stop written as X (pseudogenes) | E. coli mdtQ, efeU, ybfI, yhdW, ycgI, ybfG, ypjI |
| newer than AFDB v6 | E. coli xtpA (2025 entry) |
| longer than 2,700 aa | Kp irp1 (3,163 aa), Sa ebh (9,535 aa) |

E. coli has the most because Swiss-Prot curates its leader peptides and pseudogene fragments, which
the Kp and Sa proteomes mostly do not list.

`esmfold.py` folds the 32 up to 2,700 aa with **ESMFold v1** (Meta's `facebook/esmfold_v1` via the
Hugging Face port, already in `gradi`; CPU). Selenocysteine → cysteine for prediction (standard),
X kept as ESMFold's unknown token, both recorded. pLDDT is rescaled to 0–100 so the ≥ 70 pocket
admission means the same thing on both predictors. **The two giants are not folded — owner's
decision, 2026-10-03**: measured cost scales ~L^2.3 (85 s at 221 aa, 20 min at 715 aa, 71 min at
1,015 aa), so ~18 h even in 1,400-aa fragments. Their pocket columns stay NA. MPS was tested and
is no faster (103 s vs 85 s, identical structure: CA RMSD 0.001 Å).

`model_source` in `evidence/alphafold_<sp>.tsv` says which predictor each model came from. ESMFold
is less accurate than AlphaFold2 on average; treat those 32 rows' pocket columns accordingly.

## Pockets

**Confidence enters once.** A pocket is *admitted* when its lining residues average pLDDT ≥ 70
(AlphaFold's "confident" band). Neither tool reads pLDDT — P2Rank's `alphafold` config drops
B-factor as a feature and fpocket ignores it — so this filter is the only use. v1 used pLDDT twice,
inside a consensus score *and* as a whole-protein disorder penalty.

| species | tool | pockets | admitted | proteins with an admitted pocket |
|---|---|---|---|---|
| Kp | fpocket | 110,828 | 91.4% | 5,512 (96.2%) |
| Kp | P2Rank | 21,755 | 97.1% | 4,456 (77.8%) |
| Ec | fpocket | 91,124 | 91.7% | 4,283 (98.0%) |
| Ec | P2Rank | 18,298 | 96.6% | 3,543 (81.1%) |
| Sa | fpocket | 53,974 | 89.1% | 2,784 (96.4%) |
| Sa | P2Rank | 10,050 | 95.5% | 2,137 (74.0%) |

**Reproduces v1 exactly**: proteins with any P2Rank pocket, 4,542 Kp / 3,589 Ec — v1's numbers to
the unit (`legacy/docs/ligandability_log.md`).

**The count is P2Rank's, not fpocket's**: fpocket finds a pocket on ~97% of proteins, so its count
carries no signal. Raw counts and every pocket's residues stay in `evidence/pocket_list_<sp>.tsv`.

**Zero is not missing, and the NA is now the only thing that says so.** With a model and no
admitted pocket the columns are **0** — the tools looked and found nothing. Without a model they
are **NA** — nobody could look. Those are different claims, and `fillna(0)` collapses them, which
is the v1 mistake (`legacy/HISTORY.md:199`). Affected rows are few (Kp 1, Ec 30, Sa 1) but they
are exactly the proteins with no structural information at all.

## What counts as a drug-like ligand

v1 used a hand-built denylist of ~300 codes that grew whenever an artefact surfaced
(`legacy/src/ligandability.py:197-297`). It is replaced by established sources, each answering one
question, each a separate flag in `evidence/ligand_classes.tsv`:

| criterion | source | codes left | sites left |
|---|---|---|---|
| biologically relevant | BioLiP (Yang 2013) | 42,616 | 989,058 |
| not peptide / RNA / DNA | BioLiP classes | 42,613 | 714,573 |
| parseable structure | wwPDB CCD SMILES | 42,364 | 695,821 |
| organic (has carbon) | RDKit | 42,073 | 450,557 |
| not a polymer unit / saccharide | CCD `_chem_comp.type` | 40,623 | 385,737 |
| not a cofactor | PDBe cofactor classes (ex-EBI CoFactor) | 40,287 | 328,937 |
| not a nucleotide | 5'-phosphorylated nucleoside SMARTS (sc-PDB's category) | 39,186 | 264,501 |
| not an artefact | PLINDER badlist (265 codes, Apr 2024) | 39,082 | 259,656 |
| not an endogenous metabolite | ECMDB (3,760 E. coli metabolites), InChIKey skeleton | 38,595 | 239,888 |

**Two established filters were measured and rejected, because both delete antibiotics:**

- **QED ≥ 0.2** (PLINDER's default) fails clorobiocin 0.086, novobiocin 0.184, rifampicin 0.109,
  paromomycin 0.114, kanamycin 0.167. Natural products are large and polar; QED penalises exactly
  that. Kept as `holo_identity_qed` in `evidence/holo_<sp>.tsv`.
- **PLINDER's Ro3 "fragment"** (MW < 300, cLogP < 3, ≤ 3 HBD/HBA) drops fosfomycin and
  D-cycloserine.

**The metabolite filter was the project owner's choice (2026-10-03)** among four measured options
(proteins with bacterial holo evidence, Kp / Ec / Sa, before the PLINDER and ECMDB steps):
keep metabolites 855 / 793 / 395 · drop Ro3 fragments 691 / 653 / 323 · QED 762 / 700 / 356 ·
drop metabolites → **608 / 559 / 307** in the shipped table. ECMDB was checked to contain none of 13
common antibiotics, so it removes pyruvate, glycerol, adenine, 2-oxoglutarate without touching
small xenobiotics.

Two traps hit while building it:

1. **`[R1]` in a ring SMARTS misses fused macrocycles.** Cyclic-di-GMP's ribose atoms are also in
   the macrocycle, so a nucleotide pattern written with `[R1]` let C2E through — it was the single
   most frequent "drug-like" winner on Kp. Use `[R]`.
2. **BioLiP's curation leaks artefacts.** C8E, LDA (detergents), GOL, PEG, TRS were "best ligand"
   for dozens of proteins until PLINDER's list was added on top.

**Residual, documented rather than patched**: a few lipids and detergents PLINDER does not list
(cardiolipin CDN, LMNG `LMN`, PIPES `PIN`) still win for a handful of proteins. Growing a local
denylist is the v1 failure mode; the fix, if it matters, is a better published list.

## Holo evidence

DIAMOND (`--sensitive`, `-k 2000`) of each proteome against the 51,169 distinct sequences of the
128,554 BioLiP chains carrying a drug-like ligand. A hit counts when identity ≥ 40%
(`src.ligandability.REMOTE_PIDENT`), the PDB chain is ≥ 50% aligned (`MIN_SCOV`), and **every
binding-site residue lies inside the alignment**. **There is no query-coverage floor**: co-crystals
are often domain constructs — E. coli GyrB/clorobiocin (1kzn) is a 186-residue domain of an
804-residue protein, 23% query coverage, which `MIN_QCOV = 50` would have rejected. BioLiP gives
the site residues, so the site is tested directly. **Bacterial** = every SIFTS taxid of the chain
descends from Bacteria (NCBI taxdump).

| species | holo > 0 | own co-crystal (≥ 95) | any organism > 0 | median identity when > 0 |
|---|---|---|---|---|
| Kp | 608 (10.6%) | 80 | 668 | 74.6 |
| Ec | 559 (12.7%) | 295 | 611 | 96.7 |
| Sa | 307 (10.6%) | 86 | 347 | 57.0 |

Spot checks (E. coli, all `holo_identity` 100): folA/methotrexate, gyrB, rpoB/rifampicin, fabI,
murA/fosfomycin, lpxC, ampC, acrB/erythromycin, def. On Kp, the β-lactamase reaches 100% and
D-cycloserine (alr), fosfomycin (dxr) and rifampicin (rpoB) all carry.

**The ribosome is only partly recovered** — 9–11 ribosomal proteins per species (tetracycline,
spectinomycin, …) out of ~55. Ribosome antibiotics mostly contact rRNA, not protein, so a protein
table under-represents them by construction. The gap `docs/ligands.md` flags is narrowed, not
closed.

**Not reconciled with v1's 2,525 Kp "drug-like co-crystals."** v1 had no site-coverage test, no
taxonomy restriction on its sequence routes, and kept metabolites and sugars; the closest v2 number
(`any organism`, metabolites kept) is ~900. The difference has not been decomposed.

## Experimental structures (`pdb_n_structures`, `pdb_coverage`)

How many PDB entries exist for this protein, ligand or not, and how much of it they cover —
requested by the project owner 2026-10-03. DIAMOND of each proteome against all 1,103,516 protein
chains of the wwPDB's `pdb_seqres` (187,359 distinct sequences). A chain IS this protein at ≥ 95%
identity (`src.ligandability.DIRECT_PIDENT`) with ≥ 50% of the PDB chain aligned (`MIN_SCOV`, so a
co-crystallised peptide or a fusion partner does not count). No query-coverage floor: a one-domain
structure counts, and `pdb_coverage` says how much it covers.

| species | proteins with a structure | median entries when > 0 | coverage ≥ 0.9 |
|---|---|---|---|
| Kp | **569 (9.9%)** | 3 | 460 |
| Ec | 1,893 (43.0%) | 3 | 1,695 |
| Sa | 595 (20.6%) | 2 | 477 |

**By sequence, not accession, is the whole gain on Kp**: v1's SIFTS accession route found 30
(0.5%). E. coli agrees with v1 (1,779, 40.4%). Spot checks: Ec rpoB 411 entries, folA 177, ampC
133, lacZ 89; Sa clpP 44, ftsZ 50; Sa gyrB 2 entries at coverage 0.36 (domain constructs).

Two things to read into the numbers. **≥ 95% counts near-identical proteins of other species** —
Kp rpoB inherits E. coli's 411 RNA-polymerase entries, the intended "same protein" semantics.
**Coverage is of the deposited sequence (SEQRES)**, so loops absent from the density still count;
observed-residue coverage would need per-chain SIFTS residue mappings.

## AlphaFill: measured and left out

The plan was PDB ∪ AlphaFill. `alphafill_check.py` measured it on v1's cached responses with the
same ligand classes and floor: AlphaFill reaches **37 Kp / 21 Ec** proteins the BioLiP route misses
even with no taxonomy restriction (out of 668 / 611), and the ligands behind the gain are mostly
guanidine (GAI) and a mercury soak (EMT), plus a few real paromomycin contacts. It also gives no
donor taxonomy without a second lookup. Not in the deliverable; Sa was not fetched.

## How much to trust the pocket scores

Predicted and measured columns share no input, so agreement is a fair test (`merge.py`): AUROC of
the pocket score for separating proteins with drug-like holo evidence from those without.

| species | P2Rank, own (≥ 95) | fpocket, own | P2Rank, family (> 0) | fpocket, family |
|---|---|---|---|---|
| Kp | 0.642 | 0.538 | 0.699 | 0.575 |
| Ec | 0.669 | 0.588 | 0.667 | 0.572 |
| Sa | 0.712 | 0.526 | 0.663 | 0.527 |

**P2Rank carries real signal; fpocket's druggability barely does.** fpocket is kept because the
owner asked for both scores and it is a different method, not because it earns its place on this
test. Note the test is biased in fpocket's disfavour in one respect — the best fpocket pocket over
~20 per protein saturates (median ~0.5), compressing its range.

`src.pockets.druggability()` derives a 0–1 convenience score on the fly — the mean within-species
percentile rank of `p2rank_score` and `holo_identity` — and is deliberately not stored.

## Open leads

- **AF2BIND** (Gazizov … Ovchinnikov, Nat Methods 2026): one AF2 pass per protein with 20 bait
  residues, per-residue binding probability. Needs ColabDesign/JAX in its own env; time a 20-protein
  probe before committing, then test it against `holo_identity` exactly as P2Rank was tested above.
- **Decompose the v1 gap** (2,525 vs ~900 on Kp) criterion by criterion.
- **Human** is out of scope here; pocket similarity to human homologs is a selectivity question for
  its own stage.
