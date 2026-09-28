# Stage 06 — ligandability, part 1: measured bioactivity

`scripts/ligands/chembl.py` · `scripts/ligands/bindingdb.py` ·
`scripts/plots/ligands.py` · loader `src/ligandability.py`

Asks one question of every protein: **can a small molecule bind it?** Answered from ChEMBL 37 by
mapping the three bacterial proteomes onto ChEMBL's target sequences with DIAMOND, then counting
the **non-redundant** ligands reachable at four levels of evolutionary distance.

| species | n | direct (≥95% id) | same species | close (≥60%) | remote (≥40%) | human | complex |
|---|---|---|---|---|---|---|---|
| *K. pneumoniae* | 5,728 | 21 | 5 | 79 | **113** | 95 | 16 |
| *E. coli* | 4,403 | 65 | 66 | 73 | **96** | 86 | 15 |
| *S. aureus* | 2,889 | 32 | 32 | 36 | **78** | 56 | 10 |

Counts are proteins with ≥1 potent ligand (pChEMBL ≥ 6 = 1 µM). Buckets are **nested**:
`remote ⊇ close ⊇ species`, `remote ⊇ direct`.

![coverage](../output/plots/ligands/chembl_coverage.png)

**The headline is that this axis is sparse, and that is the finding.** 2.0% of Kp, 2.2% of Ec and
2.7% of Sa carry any potent-ligand evidence at all. It must not be forced upward — see *Coverage is
low, and correctly so*.

## Homology transfer is the whole game

ChEMBL 37 holds **21** *K. pneumoniae*, **225** *E. coli* and **80** *S. aureus* single-protein
targets, against proteomes of 5,728 / 4,403 / 2,889. Exact lookup is not a strategy, and the anchor
proteome is dark TrEMBL besides. So the axis rests on CLAUDE.md's house rule — **map by sequence,
never by accession** — with the resemblance reported next to every count.

It works: v1's validation case reproduces exactly. HS11286 `A0A0H3H184` maps to ChEMBL `Q93LQ9`
(*K. pneumoniae* β-lactamase) at **100% identity**, and that protein is invisible to accession
matching.

The payoff is visible in the `direct` column above. Kp has **21** proteins with a potent ligand on a
≥95%-identical target but only **5** whose liganded target is labelled *K. pneumoniae* — so 16 Kp
proteins have an essentially identical liganded target filed under some other organism name.

## Same-species matching is by organism NAME, not tax_id

ChEMBL files strains under their own taxids. Measured:

| query | single-protein targets |
|---|---|
| `tax_id = 562` | 65 |
| `organism LIKE 'Escherichia coli%'` | **225** (K-12 lives at 83333) |

A taxid match would have undercounted the best-represented of the three species roughly threefold.
The prefix is the two-word binomial, so `Klebsiella aerogenes` is not swept in by a genus match.

## One confidence gate cannot serve both target types

`confidence_score` is target-type specific — 9/8 direct/homologous **single protein**, 7/6
direct/homologous **complex subunits**, 5/4 family-level. Two measured consequences:

- **`≥ 8` removes 0 of 3,271,336 single-protein rows.** It is entirely subsumed by requiring a
  `pchembl_value` at all. It is kept as an explicit guard, but it is *not* quality work and must not
  be described as such.
- **`≥ 8` returns exactly zero protein complexes.** One gate for both tracks would have shipped
  empty `complex_*` columns reading as a real biological zero.

That second point is load-bearing, because **DNA gyrase is a `PROTEIN COMPLEX` in ChEMBL**:

| ChEMBL target | organism | compounds |
|---|---|---|
| `CHEMBL2094139` DNA gyrase | *E. coli* | 713 |
| `CHEMBL3038482` DNA gyrase | *S. aureus* | 491 |
| `CHEMBL3038508` Topoisomerase IV | *S. aureus* | 200 |
| `CHEMBL4662931` ClpP1P2 | *M. tuberculosis* | 15 |

GyrA/GyrB are in the consortium's own interest panel. Under a SINGLE-PROTEIN-only rule — v1's rule —
E. coli gyrB shows 295 compounds; with the complex track it also reaches 666 more, and gyrA goes
from 79 to 79 + 659. `ClpP1P2` matters for the same reason the degradation handle does.

`PROTEIN FAMILY` and `PROTEIN COMPLEX GROUP` sit at confidence 4–5, *"multiple homologous protein
targets may be assigned"* — the ligand's real target is ambiguous even within the group. They are
**measured and excluded**, not silently absent.

## Non-redundant means scaffolds

v1 counted `COUNT(DISTINCT molregno)`, so a compound and its hydrochloride salt counted twice. Here
compounds collapse to their `molecule_hierarchy` parent, and the diversity number is distinct
**Bemis–Murcko generic scaffolds**.

![redundancy](../output/plots/ligands/chembl_redundancy.png)

It changes the picture materially — a median of **2.0 compounds per scaffold**, and far more at the
top of the ranking:

| protein | compounds | scaffolds | ratio |
|---|---|---|---|
| *E. coli* `ampC` β-lactamase | 12,671 | 5,829 | 2.2 |
| *E. coli* `lpxC` | 636 | 150 | 4.2 |
| *E. coli* `folA` | 443 | **64** | 6.9 |

folA looks like 443 ligands and is really 64 chemical ideas. An acyclic compound has no Murcko
scaffold and joins one shared empty-string bucket — never dropped, never counted as zero.

## Counts are unions over a pool, not one best hit

v1 picked one target per bucket by `max(n_potent, best_pchembl)` and reported *its* counts, which
under-counts a protein resembling several liganded targets and biases toward whichever target was
screened hardest. Here a bucket's count is the union of distinct parent compounds over every target
in the pool. The closest target is still reported as `best_target` / `best_pident` / `best_organism`.

## The buckets are restricted to true Bacteria — and the alternative ships anyway

This is the fix for v1's worst bug in this axis. Its first pass put **no identity floor** on a
"non-human" bucket and reported **424** potent Kp proteins, with rat `P97697` winning the bacterial
slot for several of them; restricting to true Bacteria (`organism_class.l1`) plus a 40% floor
corrected it to **175**.

The user's own framing for this stage was three buckets over "any species". The buckets are built
over Bacteria regardless, because the alternative is a known, documented bug — but nothing is
resolved silently: the literal all-organism count ships as `allorg_n_compounds` /
`allorg_n_scaffolds`, and human ships as its own liability block. Read the difference rather than
trusting the restriction.

## The 60% band is the one arbitrary number, and it turns out not to matter

`close` was set at 60% from stage 05's measured identities (Kp↔Ec RBH median **86%**, Kp↔human
**37%**), to sit above the cross-kingdom noise floor and below the cross-genus ortholog median.

Panel C of `chembl_coverage.png` shows why the choice is low-stakes: the identity distribution of
liganded hits is **strongly bimodal** — a large mass at 40–45% and a spike at 95–100%, with a trough
between roughly 55 and 95. Very little evidence lives near 60, so moving the band moves few proteins.

`close` is also a **superset of `species` by construction**. The bands alone do not nest — a
same-species paralog at 45% identity is in `species` but below the 60% cut — and a non-monotonic
ladder is a footgun for anything downstream that subtracts one bucket from another. The consequence,
stated plainly: **a same-species target between 40 and 60% identity is counted as `close`.**

## Coverage is low, and correctly so

Roughly 2% of each proteome has potent-ligand evidence. That is a fact about how little of the
bacterial proteome anyone has ever screened, not a mapping failure, and there is no honest way to
raise it. The v1 controls agree within tolerance and in the expected direction — v2 is *lower*,
because it adds `standard_relation = '='`, `data_validity_comment IS NULL`, `potential_duplicate = 0`
and a subject-coverage floor v1 did not have:

| quantity | v1 | v2 |
|---|---|---|
| Kp proteins with a potent bacterial ligand | 175 | 113 |
| Ec proteins with a potent bacterial ligand | 155 | 96 |

`n_compounds_tested` exists so that **screened-and-clean can be told from never-screened** — a
protein with 200 compounds tested and none potent is a different object from one nobody ever tried.

## The selectivity finding

![selectivity](../output/plots/ligands/chembl_selectivity.png)

Good antibacterial targets sit on the x-axis: `lpxC`, `folA`, `gyrA`/`gyrB`, `ampC` all have
substantial bacterial evidence and **zero** human evidence.

**`clpP` does not.** It carries **106 human compounds against 61 bacterial**. Stage 05 measures Kp
clpP → human mitochondrial CLPP at 56.3% identity, and stage 04's ONC212 is an imipridone whose
characterised human target *is* ClpP. The degradation handle the whole BacPROTAC strategy depends on
has a druggable human ortholog, and this axis sees it from a completely independent direction. It
belongs in front of the collaboration.

`ftsZ` is the control in the other direction: **0 human compounds**, matching stage 05's measured
zero human hit, with real bacterial evidence.

## Cross-check: two independent routes to human identity

`human_best_pident` is measured here by DIAMOND against **ChEMBL target sequences**; stage 05's
`neighbors_of(ac, "human")` measures it against the **human reference proteome**. Different
databases, different runs, same quantity:

| species | n comparable | r | median &#124;diff&#124; |
|---|---|---|---|
| *K. pneumoniae* | 117 | 0.988 | 0.0 pp |
| *E. coli* | 110 | 0.979 | 0.0 pp |
| *S. aureus* | 67 | 0.980 | 0.0 pp |

This is free validation and it is as clean as it could be.

## BindingDB: measured, modest, not promoted

`scripts/ligands/bindingdb.py` exists to answer *does BindingDB add anything?* before
anything is built on it. v1 shipped BindingDB as a co-equal track reporting 93 potent Kp proteins
against ChEMBL's 175 — **without ever measuring the overlap**, so that 93 may have been almost
entirely redundant.

The measurement is cheap because BindingDB carries the **target chain sequence inline** (column 39)
and the **ligand InChIKey** (column 4), so both halves of the join are already in the file.

| species | ChEMBL | BindingDB | overlap | **gains first ligand** | BindingDB ligands | InChIKeys not in ChEMBL |
|---|---|---|---|---|---|---|
| *K. pneumoniae* | 113 | 87 | 72 | **15** | 2,480 | 600 (24%) |
| *E. coli* | 96 | 78 | 64 | **14** | 2,424 | 700 (29%) |
| *S. aureus* | 78 | 56 | 52 | **4** | 1,742 | 434 (25%) |

**33 proteins across all three species gain their first potent ligand — 11.5% over ChEMBL's 287.**
Real but modest, and much of the overlap is structural: ChEMBL 37 already ingests BindingDB patent
bioactivity (13,835 assays, 2.68M activities).

The verdict is deliberately left as a measurement rather than folded in. Promoting it means adding
`bindingdb_*` columns beside the ChEMBL ones — **never merged**, on the stage-04 principle that two
sources with 25% overlap are not two independent opinions.

The number is robust: widening the bacterial genus map from 78 genera (derived from hit targets
only) to 124 (from the full `organism_class` table) moved bacterial target sequences 423 → 425 and
left the verdict unchanged at 33 / 11.5%.

## Outputs

`data/processed/ligands/chembl_<species>.tsv`, one row per protein, keyed on `uniprot_ac`:

| column | meaning |
|---|---|
| `direct_hit` | a ChEMBL target at ≥95% identity exists (whether or not it has ligands) |
| `best_target`, `best_pident`, `best_organism` | the closest bacterial single-protein target |
| `n_compounds_tested` | parent compounds with any pChEMBL over the `remote` pool |
| `direct_*`, `species_*`, `close_*`, `remote_*` | `n_compounds`, `n_scaffolds`, `best_pchembl` per bucket at pChEMBL ≥ 6 |
| `human_n_compounds/_n_scaffolds/_best_pchembl/_best_pident` | selectivity liability — never merged into a bucket |
| `allorg_n_compounds`, `allorg_n_scaffolds` | the literal all-organism count, so the Bacteria restriction is visible |
| `complex_n_compounds`, `complex_best_target` | component of a liganded complex |

`accessory/`: `chembl_targets.tsv` (component↔target bridge) · `chembl_ligands.tsv` (2.59M
target×compound pairs) · `chembl_hits.tsv` (the raw DIAMOND join, unfiltered by band) ·
`chembl_targets.faa` · `scaffolds.tsv` (1.3M parent compounds → SMILES, InChIKey, molecule ChEMBL
id, generic scaffold) · `organism_class.tsv` · `cutoff_sensitivity.tsv` · `control.tsv` ·
`manifest.tsv` · `bindingdb_*`.

| helper | returns |
|---|---|
| `load(species)` / `load_all()` | the deliverable |
| `load_targets` / `load_ligands` / `load_hits` / `load_scaffolds` | the accessory matrices |
| `load_cutoff_sensitivity` / `control` / `manifest` | the run record |
| `evidence_level(df)` | `direct` / `species` / `close` / `remote` / `none` — the tightest bucket with a potent ligand |
| `selectivity_risk(df)` | `log2((human+1)/(bacterial+1))` on scaffolds; a triage flag, not a cross-reactivity prediction |

## Cutoff sensitivity

Reported, never optimised — stage 04's lesson, that counts are not comparable across cutoffs:

| pChEMBL | Kp remote | Ec remote | Sa remote |
|---|---|---|---|
| ≥ 5 (10 µM) | 146 | 129 | 107 |
| **≥ 6 (1 µM)** | **113** | **96** | **78** |
| ≥ 7 (100 nM) | 82 | 69 | 60 |

6 is the field-standard "genuine ligand" line and v1's choice, so the controls stay comparable.
Requiring `pchembl_value IS NOT NULL` already performs most quality filtering for free.

## Running it

```bash
python scripts/ligands/chembl.py --dry-run
python scripts/ligands/chembl.py --limit 200 --species ecoli   # -> accessory/smoke_*
python scripts/ligands/chembl.py                               # ~2 min with the dump
python scripts/ligands/bindingdb.py                            # ~8 min, first pass
python scripts/plots/ligands.py
```

The 30.5 GB dump is needed only for the first run; **a cached re-run takes 0.3 min and needs no
database at all**. Restore it with `tar -xzf data/raw/other/chembl/chembl_37_sqlite.tar.gz -C
data/raw/other/chembl` — see `data/raw/other/chembl/SOURCE.md`.

## Traps

- **The rat trap.** An unrestricted "non-human" bucket with no identity floor gave 424 potent Kp
  proteins against a true 175, rat `P97697` winning bacterial slots at ~30% identity. Never widen
  the bucket without both the Bacteria restriction and the identity floor.
- **`confidence_score ≥ 8` returns zero protein complexes**, by definition of ChEMBL's scale. Using
  one gate for both target types silently empties the complex columns.
- **`tax_id = 562` finds a third of E. coli's targets.** Match the organism name.
- **`chembl_37_blast.fa.gz` is target-redundant** — `>CHEMBL1907607_O09028`, one record per target,
  so a shared sequence repeats. Build the FASTA from `component_sequences`.
- **BindingDB's TSV is ragged**: measured 50 / 62 / 74 / 86 columns in the first 200k rows (50 base
  plus 12 per extra target chain). `pd.read_csv` with a fixed header breaks; stream it with `csv`
  and raise `csv.field_size_limit`.
- **BindingDB has no taxonomy id.** v1 bridged the superkingdom via ChEMBL genus names but *warned
  and returned an empty set* when the file was missing, silently emptying its bacterial bucket. Here
  the map comes from `accessory/organism_class.tsv` and a missing or empty map exits non-zero; the
  count of unclassifiable organism strings (207) is printed every run.
- **Do not use Python's `hash()` for a cache key.** String hashing is randomised per process, so a
  builtin-hash sequence id differs between the run that writes a cache and the run that reads it.
  `hashlib.sha1` throughout.
- **`pchembl_value` only exists** for `=` relations on IC50/EC50/XC50/AC50/Ki/Kd/Potency in nM. MIC
  and %-inhibition are absent by construction, so whole-cell antibacterial screening does not
  contribute — **this axis does not say "has an antibiotic"**. It is also why the ribosome, whose
  drugs are overwhelmingly measured as MIC, is largely absent from the complex track.
- **rdkit must stay arm64.** The PyPI `macosx_11_0_arm64` wheel is fine; a conda install that flips
  `gradi` to osx-64 would take ESM-C down with it.

## Run log

### 2026-09-03 — first full run, three species

Extraction: 2,520,321 single-protein (target, compound) pairs over **8,167 targets** and 1,267,120
parent compounds; 71,205 complex pairs over **469 targets**. 8,469 distinct target sequences in the
DIAMOND database. Assay split **61.5% B / 38.5% F** — confirming that a binding-only filter would
discard over a third of the evidence, which for bacterial enzymology is where most of it lives.

Scaffolds: 1,303,734 of 1,305,242 parent compounds have a structure (99.9%); **187,119 distinct
generic scaffolds**, 5,338 acyclic. ~90 s on 9 cores.

Coverage guards: the subject-coverage floor v1 lacked drops 10.7% of Sa's ≥40%-identity hits after
the query-coverage filter has already run. Kept — a hit covering half our protein and a tenth of a
much longer ChEMBL target is not evidence about our protein.

**Four bugs the assertions caught, all of which would have shipped silently.**

1. **The buckets were not nested.** `species` (same species, ≥40%) is not inside `close` (≥60%), so a
   same-species paralog at 45% broke the ladder. Fixed by making `close` a superset by construction.
2. **`component_id` was read as a string**, which would have made the DIAMOND merge return **zero
   rows** rather than raising. The loader now coerces identifier columns explicitly, and bare-named
   numeric columns (`pchembl`, `qcov`, `scov`, `pident`) alongside them.
3. **`human_best_pident` was never computed** — it fell through to the NaN fallback, and the stage-05
   cross-check silently compared 0 proteins and reported "not comparable" rather than failing.
4. **The scaffold cache could never be reused.** Its check required the cache to cover every
   requested molregno, but ~1,500 parent compounds have no structure in ChEMBL at all, so the subset
   test always failed and every run re-derived 1.3M scaffolds — and demanded the 30 GB dump back.

A fifth, in the plots: seeding `np.random.default_rng(0)` once per axis gave x and y **identical**
jitter offsets, drawing spurious diagonal streaks through the selectivity cloud.

**`evidence_level` was wrong before the `direct` bucket existed.** It inferred "direct" from
`direct_hit & species_n_compounds > 0`, which called only 5 Kp proteins direct against 30 with a
≥95% hit. Adding a real `direct` bucket (potent compounds on targets at ≥95% identity) gives 21 —
and exposes that 16 Kp proteins have an essentially identical liganded target filed under another
organism's name.

The database was deleted after the run (30.5 GB reclaimed); `data/raw/other/chembl/SOURCE.md`
records the recovery procedure and byte counts. It is public and re-derivable — **never upload it
to eosvc.**

## Open leads

- **Promote BindingDB**, if 33 proteins is judged worth it. Columns beside, never merged.
- **The newer BindingDB is already on disk** — `data/raw/legacy/bindingdb/BindingDB_All_202605_tsv.zip`
  (2026-05) is a year newer than the 2025-04 in use, and `PROVENANCE.md` flags it as the upgrade
  candidate. Re-running the gain measurement against it is ~8 min.
- **Structure-based ligandability** — PDB co-crystals, AlphaFill transplants, fpocket/P2Rank pockets.
  v1 built all of them (`legacy/scripts/06c`–`06g`) and they reach far more proteins than
  bioactivity does (v1: 2,525 Kp with a drug-like co-crystal against 175 with a potent ligand). That
  is the obvious part 2, and the composite score belongs after it, not before.
- **The ribosome is missing** and the reason is structural, not incidental: its drugs are measured
  as MIC, which carries no pChEMBL. Given stage 04's top-100 is ribosome-heavy, a targeted route to
  ribosome-binding evidence (PDB co-crystals will find it) matters more here than elsewhere.
- **`allorg` versus `remote` is unexploited.** The gap between them is a measured statement about how
  much apparent ligandability comes from non-bacterial homologs, which is exactly the quantity v1's
  rat trap got wrong. Worth a figure.

---

# Part 3 — ligand precedent for an arbitrary sequence

`scripts/ligands/precedents.py` (CLI) and `src/precedents.py` (library) answer a different question
from Part 1. Part 1 asks *what does our proteome have*, needs the 30.5 GB ChEMBL dump, and produces
a per-species table. This asks *what about **this** sequence* — any sequence, in about a second,
from 82 MB of cached extracts.

```
(a) n_ligands_exact      ligands on an EXACT match: UniProt accession, or identical sequence
(b) n_ligands_bacteria   UNIQUE ligands across BACTERIAL targets, by identity
(c) n_ligands_human      ligands on HUMAN orthologs
```

## What it runs on

| file | size | holds |
|---|---|---|
| `scratch/chembl_targets.faa` | 5.2 MB | 8,469 sequences, headers = `component_id` — the DIAMOND subject DB |
| `scratch/chembl_targets.tsv` | 1.0 MB | 9,347 targets with `organism` and `superkingdom` |
| `scratch/chembl_ligands.tsv` | 76 MB | 2,591,526 rows, `tid → parent_molregno` |

`sequence → DIAMOND → component_id → tid → parent_molregno`. Measured composition: **1,305,242
distinct compounds**, of which **110,019** sit on Bacteria targets and **1,101,363** on human.
`superkingdom` gives the (b)/(c) split directly — 703 bacterial targets against 5,546 human.

## Three things the counts mean, and one they do not

**(b) counts MOLECULES, not target-compound pairs.** `parent_molregno` is ChEMBL's
`molecule_hierarchy` parent, so salts are already collapsed; taking the DISTINCT set over the union
of every passing target means a compound tested against three homologs counts **once**. Summing
per-target counts would inflate it, and the wider the identity band the worse it gets — which is
precisely the regime this tool is for.

**(c) is a liability, not a precedent.** A ligand-bearing human ortholog says the fold is druggable
*and* that hitting it may be dangerous. E. coli `clpP` is the case: **136 bacterial against 210
human** at 56.3% identity, the human mitochondrial CLPP ortholog. The CLI prints a warning when (c)
exceeds (b). Never sum them.

**The complex track is separate and is not inside (b).** E. coli `gyrB`: **666 single-protein
ligands, 1,412 complex**. Kp `A0A0H3H0Y6` (gyrA): **1,410 complex against 131 single.** DNA gyrase is a
`PROTEIN COMPLEX` in ChEMBL and v1's single-protein-only rule made GyrA/GyrB look unliganded.

**It does not say "has an antibiotic".** `pchembl` is 100% populated in the extract, so every count
is of potency-measurable ligands — `=` relations on IC50/EC50/Ki/Kd/Potency in nM. MIC and
%-inhibition are absent by construction, so a ribosomal protein looking empty here is a fact about
assay type, not about biology.

## Identity alone cannot transfer a ligand count

The case that proves it, recorded so it is not re-investigated. Kp `A0A0H3GWM6` is **99.2%
identical** to E. coli `P0ADG7`, which carries 32 compounds at pChEMBL 8.82 — and correctly gets
**nothing**. The raw hits show the match *was* found (component 1947 = CHEMBL3630, 99.2% identity)
with `qcov 100.0` but **`scov 26.6`**: the Kp entry is a **130-aa fragment** against a **488-aa**
IMP dehydrogenase. `MIN_QCOV = MIN_SCOV = 50.0` rejected it, exactly as the comment above those
constants anticipates.

Both floors are **imported** from `src/ligandability.py`, never restated, so the tool cannot drift
from the axis it sits beside.

## The control

`evidence/precedent_control.tsv`, written on every `--species` run: agreement with
`chembl_<sp>.tsv` **at matched semantics** (pChEMBL ≥ 6, complex track included).

| species | both positive | precedents-only | chembl-only | agreement |
|---|---|---|---|---|
| kpneumoniae | 113 | 4 | **0** | **99.9%** |
| ecoli | 96 | 4 | **0** | **99.9%** |
| saureus | 78 | 4 | **0** | **99.9%** |

At *default* settings the two differ, and both reasons are design rather than defect: the tool
excludes the complex track, and counts every measurable ligand rather than only pChEMBL ≥ 6 — 66 of
67 default-only Kp proteins sit below 6. Restore both and they reconcile, which is what the control
records.

## Batch output

`--species` runs a whole anchor proteome through the same code path and writes
`precedents_<species>.tsv`, complete and in canonical row order. Measured:

| species | any exact | any bacterial | any human |
|---|---|---|---|
| kpneumoniae | 6 | 179 | 114 |
| ecoli | 112 | 159 | 106 |
| saureus | 27 | 117 | 66 |

Kp's exact count is low by construction — ChEMBL holds only 21 *K. pneumoniae* targets, and the
anchor is a dark TrEMBL proteome whose accessions are largely absent. That is the same fact that
makes homology transfer the whole game in Part 1.

Load through `src/ligandability.py` — `load_precedents(species)`.

## The bug an adversarial audit found, and why the first control missed it

The first version of this tool shipped with a defect worth recording, because the shape of it
recurs.

**`component_id → tid` is ONE-TO-MANY and was collapsed with a dict.** `chembl_targets.tsv` has
9,347 rows over **8,469 distinct `component_id`** — **490 components map to more than one `tid`**,
up to 15. `dict(zip(component_id, tid))` keeps the last and drops the rest silently. Measured cost:
**21 proteins understated by 3,324 ligands.**

The worst case is the one the design was explicitly built to prevent. Component 166 (E. coli gyrA)
carries `tid 53` (SINGLE PROTEIN, 117 ligands) **and** `tid 104721` (PROTEIN COMPLEX). The complex
row came last, so the single-protein ligands vanished and **gyrA read `n_ligands_bacteria = 0`** —
v1's "GyrA/GyrB look unliganded" error, re-entering through a different door. `chembl.py:481` does
the same join correctly with a one-to-many merge.

| protein | | shipped | correct |
|---|---|---|---|
| E. coli gyrA | | **0** | **131** |
| E. coli gyrB | | 371 | **666** |
| E. coli folA | bacterial | 183 | **595** |
| E. coli folA | best pChEMBL | 9.02 | **10.92** |

**The control could not see it, by construction.** It compared `>0` against `>0`, so it reported
99.93% agreement while the underlying counts disagreed on 18 of 309 ligand-bearing Kp proteins by
1,019 ligands. A presence check cannot detect a magnitude error. The control now compares totals
(`n_count_mismatch`, `total_ligand_delta`, `max_abs_delta`).

**And the control's "matched semantics" were not matched.** It added this tool's complex track to
one side, while `chembl.py:452` defines the remote pool as `single & bact` — single-track only.
That inflated the apparent delta to +1,784 and hid the true residual. With both sides on the single
track and the collapse fixed, the two implementations agree **exactly**: 0 count mismatches, 0
total delta, 0 either-way, on all three species.

Three smaller defects from the same audit, all fixed: identical sequences shared by several
components were likewise collapsed (10 exact rows understated — Kp `bla` read 24 against a true
241); `best_pchembl_bacteria` was computed across both tracks while the count was single-only,
putting "0 ligands, best pChEMBL 9.68" on 14 rows; and an id containing a space was silently
truncated by DIAMOND to a zero count, now refused.

## One honest caveat about (c)

The 40% identity floor is inherited from a rule calibrated for **bacterial** transfer, and it is
marginal across kingdoms. Measured on Kp: bacterial hits sit at a median **84.1%** identity, human
hits at **45.3%**, with **53 of 114 in the 40–45% band** and none above 80%. So (c) is
systematically weaker evidence than (b). `clpP` at 56.3% is a real ortholog; a 40.1% human hit is
not the same claim. Filter on `best_pident_human` in `precedents_full_<sp>.tsv` before relying on it.
