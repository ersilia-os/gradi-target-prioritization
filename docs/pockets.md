# Structural ligandability (`pockets/`)

**Can a small molecule bind this protein's fold?** The bioactivity track (`docs/ligands.md`) asks
whether somebody *measured* a potent compound; this asks whether the *structure* says one could
bind. Three kinds of answer that share no input and are never merged:

- **Predicted** — pockets on the AlphaFold model, by fpocket 4.0 and P2Rank 2.5.1.
- **Measured** — drug-like ligands in this protein's *own* PDB structures, and how many
  structures it has at all.
- **Modelled** — drug-like ligands AlphaFill transplants onto its AlphaFold model from remote
  homologs.

Deliverable `data/processed/pockets/pockets_<species>.tsv`, complete and canonical for the three
bacteria (`python -m src.matrices`):

| column | kind | meaning |
|---|---|---|
| `uniprot_ac` | key | |
| `p2rank_score` | predicted, 0–1 | P2Rank (`-c alphafold`) calibrated probability of the best admitted pocket |
| `fpocket_score` | predicted, 0–1 | fpocket druggability score of the best admitted pocket |
| `n_ligands_pdb` | measured, int | non-redundant drug-like ligands bound in this protein's **own** PDB structures; 0 = none |
| `n_ligands_alphafill` | modelled, int | non-redundant drug-like ligands AlphaFill transplanted onto its model; 0 = none |
| `pdb_n_structures` | measured, int | PDB entries with a chain that IS this protein (≥ 95% identity, ≥ 50% of the chain aligned), ligand or not; 0 = none |
| `af_plddt` | confidence | mean AlphaFold pLDDT; **NA means no model** |

### There is no `evidence` column

Dropped on the project owner's instruction, **2026-10-03**. It held
`pdb+af` · `af_only` · `pdb_only` · `none`, and every one of those is a function of two columns
that are still here — verified exactly reconstructible on all three species before removal:

| old label | = |
|---|---|
| `pdb+af` | `af_plddt.notna()` and `n_ligands_pdb > 0` |
| `af_only` | `af_plddt.notna()` and `n_ligands_pdb == 0` |
| `pdb_only` | `af_plddt.isna()` and `n_ligands_pdb > 0` |
| `none` | `af_plddt.isna()` and `n_ligands_pdb == 0` |

**The one thing it added that these columns cannot say** is where the model came from — AlphaFold
DB or ESMFold, for the proteins AFDB does not cover. That is `model_source` in
`evidence/alphafold_<species>.tsv`, via `load_alphafold()`. An ESMFold row's `af_plddt` is
ESMFold's pLDDT on the same 0–100 scale, so the column stays comparable either way.

**The two ligand counts are never summed**, and both are *non-redundant*: distinct Bemis–Murcko
generic scaffolds, reusing the ChEMBL axis's own `_scaffold_chunk` rather than a second copy.
Counting raw PDB ligand codes instead overstates by 25–34%; those raw counts ship as `n_codes_*`
in `evidence/ligand_counts_<species>.tsv`.

Load through `src/pockets.py` (`load`, `load_all`, `load_alphafold`, `load_pockets`, `load_pdb`,
`load_pdb_chains`, `load_ligands_pdb`, `load_ligands_alphafill`, `load_transplants`,
`load_ligand_counts`, `load_ligand_classes`, `druglike`). There is **no `druggability()` helper**:
one shipped briefly and was removed with `holo_identity`, because no defensible weighting exists
across a weak prior, a sparse measurement and a third party's model. Combine at prioritisation
time, across axes, where the weighting is an explicit choice someone owns.

## Pipeline

```
structures.py      AlphaFold DB v6, validated by sequence   -> evidence/alphafold_<sp>.tsv
esmfold.py         ESMFold for proteins AFDB lacks          -> scratch/esmfold/, evidence/esmfold_<sp>.tsv
structures.py      (again, to pick the ESMFold models up)
predict.py         fpocket + P2Rank (gradi-pockets)         -> evidence/pocket_list_<sp>.tsv
pdb_coverage.py    DIAMOND vs every PDB protein chain       -> evidence/pdb_<sp>.tsv, pdb_chains_<sp>.tsv
alphafill.py       fetch + flatten transplants              -> evidence/transplants_<sp>.tsv
holo.py            drug-likeness vocabulary + both counts   -> evidence/{ligand_classes,ligands_*,ligand_counts}_<sp>.tsv
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

The vocabulary covers **BioLiP's codes plus the 121 AlphaFill transplants that BioLiP never
lists** (22 of them drug-like). `in_biolip` is a flag, not a requirement — for those 121, BioLiP's
own curation cannot apply and only the remaining criteria do.

| criterion | source | codes left | sites left |
|---|---|---|---|
| all codes | BioLiP (Yang 2013) + AlphaFill | 42,737 | 989,058 |
| not peptide / RNA / DNA | BioLiP classes | 42,734 | 714,573 |
| parseable structure | wwPDB CCD SMILES | 42,485 | 695,821 |
| organic (has carbon) | RDKit | 42,185 | 450,557 |
| not a polymer unit / saccharide | CCD `_chem_comp.type` | 40,733 | 385,737 |
| not a cofactor | PDBe cofactor classes (ex-EBI CoFactor) | 40,396 | 328,937 |
| not a nucleotide | 5'-phosphorylated nucleoside SMARTS (sc-PDB's category) | 39,292 | 264,501 |
| not an artefact | PLINDER badlist (265 codes, Apr 2024) | 39,113 | 259,656 |
| not an endogenous metabolite | ECMDB (3,760 E. coli metabolites), InChIKey skeleton | 38,617 | 239,888 |

**Two established filters were measured and rejected, because both delete antibiotics:**

- **QED ≥ 0.2** (PLINDER's default) fails clorobiocin 0.086, novobiocin 0.184, rifampicin 0.109,
  paromomycin 0.114, kanamycin 0.167. Natural products are large and polar; QED penalises exactly
  that. Kept as the `druglike_qed` flag in `evidence/ligand_classes.tsv`.
- **PLINDER's Ro3 "fragment"** (MW < 300, cLogP < 3, ≤ 3 HBD/HBA) drops fosfomycin and
  D-cycloserine.

**The metabolite filter was the project owner's choice (2026-10-03)** among four options measured
at the time on the since-removed sequence-transfer route (proteins reached, Kp / Ec / Sa): keep
metabolites 855 / 793 / 395 · drop Ro3 fragments 691 / 653 / 323 · QED 762 / 700 / 356 · drop
metabolites 608 / 559 / 307. Those counts are historical — the route they were measured on is
gone — but the ranking they establish is what the choice rests on. ECMDB was checked to contain none of 13
common antibiotics, so it removes pyruvate, glycerol, adenine, 2-oxoglutarate without touching
small xenobiotics.

Two traps hit while building it:

1. **`[R1]` in a ring SMARTS misses fused macrocycles.** Cyclic-di-GMP's ribose atoms are also in
   the macrocycle, so a nucleotide pattern written with `[R1]` let C2E through — it was the single
   most frequent "drug-like" winner on Kp. Use `[R]`.
2. **BioLiP's curation leaks artefacts.** C8E, LDA (detergents), GOL, PEG, TRS were "best ligand"
   for dozens of proteins until PLINDER's list was added on top.

**Residual, documented rather than patched**: a few lipids, detergents and tiny molecules that
PLINDER does not list still pass — LMNG (`AV0`), octadecane (`8K6`), guanidine (`GAI`), cardiolipin
(`CDN`), PIPES (`PIN`). Growing a local denylist is the v1 failure mode; the fix, if it matters, is
a better published list.

## Ligands bound: the two counts

### `n_ligands_pdb` — this protein's own structures (MEASURED)

`pdb_coverage.py` already matched our proteins to PDB **chains** at ≥ 95% identity. Those chains
*are* this protein, so BioLiP's ligand rows attach straight to them by `(pdb, chain)`: **no
alignment of our own, no transfer, and no binding-site test to make.** Keep the drug-like ones,
collapse to scaffolds, count.

| species | proteins with ≥ 1 | ligand codes | scaffolds | redundancy collapse |
|---|---|---|---|---|
| Kp | 88 (1.5%) | 544 | 367 | 32.5% |
| Ec | 308 (7.0%) | 1,289 | 855 | 33.7% |
| Sa | 90 (3.1%) | 345 | 230 | 33.3% |

**A zero here is overwhelmingly "nobody has solved it with a drug bound", not "undruggable"** —
~95% of each proteome. Compare `pdb_n_structures`: Ec has 1,893 proteins with a structure but only
308 with a drug-like ligand in one.

Named spot checks, enforced by `merge.py` (a broken join still yields plausible counts, so the
check has to be specific molecules): Ec `folA` carries **methotrexate and trimethoprim** among 23;
Sa `gyrB` carries **novobiocin**; Sa `fabI` carries **triclosan** among 24; Kp `clpP` carries the
ClpP activator `KHS`; Ec `ftsZ` carries **none**, which is correct — it has no drug co-crystal of
its own, while Sa `ftsZ` has 9.

### `n_ligands_alphafill` — transplanted onto the model (MODELLED)

AlphaFill (Hekkelman et al., *Nat Methods* 2023) superposes homologous PDB structures onto an
AlphaFold model and transplants their ligands. Donors are at **≥ 25% identity over ≥ 85 aligned
residues** — in our data a median of ~30% — so this is genuinely remote transfer, done
structurally rather than by sequence.

| species | proteins with ≥ 1 | ligand codes | scaffolds | redundancy collapse |
|---|---|---|---|---|
| Kp | 1,533 (26.8%) | 4,729 | 3,461 | 26.8% |
| Ec | 1,196 (27.2%) | 3,336 | 2,488 | 25.4% |
| Sa | 704 (24.4%) | 2,106 | 1,532 | 27.3% |

**Nothing is filtered on quality** (owner's call). AlphaFill publishes `local_rmsd` — the backbone
fit within 6 Å, median 0.84 Å here, p75 2.0 Å — as a metric but **no threshold**, so inventing one
would be ours, not theirs. Every transplant ships with its RMSD, donor identity and clash count in
`evidence/ligands_alphafill_<sp>.tsv`; filter there.

**`analogue_id`, not `compound_id`, is the ligand.** Where the donor holds a close analogue,
AlphaFill records the donor's own compound as `compound_id` and what it stands in for as
`analogue_id` — ANP→ATP, ACO→CoA, AGS→ATP, TGG→glutathione, 2.7% of transplants. The analogue is
what AlphaFill asserts the protein binds, so it is what gets classified; conveniently it also
resolves these cases to the parent cofactor, which the drug-likeness filter then correctly removes.

### Why we do not transfer ligands ourselves

Until 2026-10-03 this axis shipped `holo_identity`: DIAMOND against BioLiP's holo chains, a
binding-site span test, and the best hit's % identity as the score. It was removed because it was
a **hand-rolled AlphaFill**, and a worse one. Measured on the same proteomes with the same
drug-likeness rule: **AlphaFill 1,533 Kp / 1,196 Ec proteins; the sequence route 184 / 172.**

Three further faults, recorded so the idea is not revived:

1. **The scale was bimodal, not continuous** — nonzero values were either ≥ 95% (the protein's own
   structure) or 40–60% (a distant relative), 32–57% of them below 60%.
2. **Identity is a weak axis here anyway** — `ligands/transfer_calibration.py` measured
   P(potent | neighbour potent) as flat from 25% to 100% and concluded the bands control coverage,
   not reliability.
3. **It never checked the binding site was conserved** — only that the site's residue *span* fell
   inside the alignment. Realigning properly showed sites *are* better conserved than whole
   proteins (median site identity 67–80% at 40–50% global identity), so the column also *understated*
   its own evidence while keeping the 16–21% of low-identity cases where the pocket genuinely differs.

The same deletion removed `alphafill_check.py`, whose conclusion that "AlphaFill adds little" came
from imposing a 40% identity floor **AlphaFill does not use**.

**The ribosome stays only partly covered** — ribosome antibiotics mostly contact rRNA, not protein,
so a protein table under-represents them by construction. The gap `docs/ligands.md` flags is
narrowed, not closed.

**v1's "2,525 Kp drug-like co-crystals" is not comparable** and was never reconciled: v1 had no
site test, no taxonomy restriction, and kept metabolites, cofactors and sugars.

## Experimental structures (`pdb_n_structures`)

How many PDB entries exist for this protein, ligand or not — requested by the project owner
2026-10-03. A `pdb_coverage` column (how much of the protein those entries cover) shipped briefly
and was dropped from the deliverable the same day on the owner's instruction; it stays per protein
in `evidence/pdb_<sp>.tsv`. DIAMOND of each proteome against all 1,103,516 protein
chains of the wwPDB's `pdb_seqres` (187,359 distinct sequences). A chain IS this protein at ≥ 95%
identity (`src.ligandability.DIRECT_PIDENT`) with ≥ 50% of the PDB chain aligned (`MIN_SCOV`, so a
co-crystallised peptide or a fusion partner does not count). No query-coverage floor: a one-domain
structure counts.

**Partial structures are found and counted.** Measured: 312 / 903 / 259 entry matches (Kp / Ec /
Sa), 1–6% of all, cover less than half the protein. Chains aligned over fewer than 30 residues are
15–62 matches per species, and only **3 proteins per species** are matched through them alone —
genuinely tiny proteins with real structures (Ec P0AD89, the 24-aa TnaC leader peptide, solved in
the ribosome, 7oj0; Sa P0C7Y1, a 21-aa phenol-soluble modulin). So no extra length floor is needed.

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

## How much to trust the pocket scores

Predicted and measured columns share no input, so agreement is a fair test (`merge.py`): AUROC of
the pocket score for separating proteins that **have a drug-like ligand in their own PDB
structures** from those that do not. **The raw number is confounded by protein length**: length
ALONE scores 0.65–0.67, because big proteins are both crystallised more often and offer more
surface for pockets. The column to quote is the AUROC *within length deciles*.

| species | length alone | P2Rank raw | **P2Rank within length** | fpocket raw | **fpocket within length** |
|---|---|---|---|---|---|
| Kp | 0.673 | 0.615 | **0.494** | 0.511 | **0.435** |
| Ec | 0.652 | 0.659 | **0.561** | 0.582 | **0.523** |
| Sa | 0.657 | 0.702 | **0.619** | 0.534 | **0.477** |

**The pocket scores are a weak and inconsistent guide to where drug-like ligands are actually
found.** On *K. pneumoniae* — the anchor — P2Rank adds **nothing** over protein size (0.494), and
fpocket is below chance. Only *S. aureus* shows a clear signal (0.619). Treat the pocket columns as
a soft prior and prefer the measured columns wherever they are non-zero.

These numbers are stricter than the ones this section carried before 2026-10-03 (P2Rank 0.54–0.58).
That is the label changing, not the model: the old label was "a drug-like ligand somewhere in the
family at ≥ 40% identity", the new one is "a drug-like ligand on this protein itself".

### There is no pocket count — removed 2026-10-03, after measuring it

`p2rank_n_pockets` counted admitted pockets at probability ≥ 0.5. **The P2Rank authors recommend no
cutoff**: PrankWeb lists every pocket and the papers evaluate by rank (top-n), never by probability.
The probability is calibrated on HOLO4K as P(x) = Tₓ / (Tₓ + Fₓ), so 0.5 is only the point where true
and false pockets balance. (The ">0.8 / 0.5–0.8 / <0.5" bands seen online are PocketDock's, a
third-party tool.) Swept on all three species, AUROC within length deciles:

| cutoff | none | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|
| Kp | 0.485 | 0.512 | 0.548 | 0.565 | 0.574 | 0.570 | 0.551 |
| Ec | 0.477 | 0.500 | 0.535 | 0.546 | 0.561 | 0.564 | 0.540 |
| Sa | 0.432 | 0.453 | 0.497 | 0.512 | 0.522 | 0.546 | 0.550 |
| ρ with length | 0.84 | 0.77 | 0.72 | 0.68 | 0.64 | 0.59 | 0.51 |

Without a cutoff the count is protein size and nothing else; the sum of probabilities (a cutoff-free
"expected number of sites") is no better (0.50–0.56). Signal appears only from 0.3 to 0.5 and never
exceeds `p2rank_score`'s. So the column was removed rather than re-tuned. Every pocket is still in
`evidence/pocket_list_<sp>.tsv`.

There is no combined score. See the note under the column table: `druggability()` was removed
with `holo_identity`, because no defensible weighting exists across a weak prior, a sparse
measurement and a third party's model.

## Open leads

- **AF2BIND** (Gazizov … Ovchinnikov, Nat Methods 2026): one AF2 pass per protein with 20 bait
  residues, per-residue binding probability. Needs ColabDesign/JAX in its own env; time a 20-protein
  probe before committing, then test it against `n_ligands_pdb` within length deciles, exactly as
  P2Rank is tested above. The bar it has to clear is low — P2Rank manages 0.494 on Kp.
- **Why is Kp so much worse than Sa** (P2Rank 0.494 vs 0.619 within length)? With only 88 Kp
  proteins carrying a measured ligand it may be sample size, or it may be that Kp's ligand-bearing
  proteins are the conserved enzymes whose structures came from E. coli.
- **Decompose the v1 gap** (v1's 2,525 Kp "drug-like co-crystals" against this axis's 88 own-PDB
  proteins) criterion by criterion. The definitions differ in at least four ways.
- **Human** is out of scope here; pocket similarity to human homologs is a selectivity question for
  its own stage.
