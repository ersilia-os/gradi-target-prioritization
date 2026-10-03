# Stage 06 — ligandability, part 1: measured bioactivity

`scripts/ligands/chembl.py` · `effort.py` · `transfer_calibration.py` · `bindingdb.py` ·
`precedents.py` · `scripts/plots/ligands.py` · loaders `src/ligandability.py`, `src/precedents.py`

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

## The 95/60/40 bands CANNOT be calibrated, and that is the result

`close` was set at 60% from stage 05's measured identities (Kp↔Ec RBH median **86%**, Kp↔human
**37%**), to sit above the cross-kingdom noise floor and below the cross-genus ortholog median.
Panel C of `chembl_coverage.png` argued the choice was low-stakes on distributional grounds: the
identity of liganded hits is **strongly bimodal** — a mass at 40–45%, a spike at 95–100%, a trough
between — so very little evidence lives near 60 and moving the band moves few proteins.

**`scripts/ligands/transfer_calibration.py` supersedes that argument by measuring the thing the
bands are supposed to encode.** It is self-contained inside ChEMBL, which is the only place the
measurement is possible: for a **pair** of bacterial targets we know *both* ligand sets, so "does
ligand evidence travel at this identity?" has an answer. **1,582 pairs over 687 sequences and 131
species.** Three negatives, all of them useful.

### 1. Compound-set overlap does not transfer at any identity

**Median Jaccard is ~0.00 in every band, 95–100% included.** Among the pairs that share anything at
all it is **0.018–0.036**, with no trend across the ladder.

Two near-identical ChEMBL targets are **one enzyme screened twice against different libraries** —
they share the protein, not the chemistry. So Jaccard measures *campaign coincidence*, not
transferability, and cannot calibrate a transfer rule. That negative is the reason the conditional
below is the statistic.

### 2. `P(potent | neighbour potent)` is FLAT from 25% to 100% identity

Pairs where both targets carry ≥5 assayed compounds, by the identity of the pair:

| identity band | n_pairs | P(potent \| neighbour potent) | lift over base |
|---|---|---|---|
| 25–30 | 60 | 0.866 | 1.40× |
| 30–40 | 178 | 0.929 | 1.50× |
| 40–50 | 106 | 0.976 | 1.58× |
| 50–60 | 38 | 0.844 | 1.36× |
| 60–70 | 23 | 0.952 | 1.54× |
| 70–80 | 10 | 0.941 | 1.52× |
| 80–90 | **2** | 0.667 | 1.08× |
| 90–95 | **3** | 1.000 | 1.62× |
| 95–100 | 36 | 0.957 | 1.55× |

Base rate **0.619**. **No decay, and a lift of only ~1.36–1.58× across the whole range.** The cause
is selection, and it is structural: **62% of bacterial ChEMBL targets already carry a potent
compound**, because a protein enters ChEMBL when somebody believed it was druggable.

**So the identity floor controls COVERAGE, not transfer reliability. Document it as a conservatism
choice; never as an accuracy threshold.** Calibrating it honestly would need proteins nobody chose
to screen, which ChEMBL by construction does not contain.

**Read `n_pairs` before quoting a band.** 80–90 rests on **2 pairs** and 90–95 on **3** — they are
noise, and the population is bimodal for the same reason the identity histogram is: a bacterial
target's nearest ChEMBL relative is either the same enzyme in another strain or a different family,
rarely anything between.

### 3. Neither species nor RBH adds anything beyond identity

At 95–100%, same-species scores **0.955** against cross-species **0.960** — indistinguishable. RBH
is *below* the unrestricted conditional in five of the nine bands (0.44 against 0.93 at 30–40%, on
10 pairs). **So no orthology criterion is added to the bucket rule**, and the useful part is
knowing the simple rule was not leaving anything on the table before anyone builds the complicated
one.

Evidence `evidence/transfer_calibration.tsv` (band × population × split), pairs in
`scratch/transfer_pairs.tsv`. CLI: `--min-compounds 5` · `--threads` · `--dry-run` · `-q`. DIAMOND
from `gradi-ortho`. ~2 min.

### The bands still nest, and that part is not arbitrary

`close` is a **superset of `species` by construction**. The bands alone do not nest — a
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

`n_compounds_tested` was meant to make **screened-and-clean tellable from never-screened** — a
protein with 200 compounds tested and none potent is a different object from one nobody ever tried.

**CORRECTED 2026-09-28: it does not do that, and an earlier version of this page claimed it did.**
`chembl.py:594` computes it from the ligand table, which is already filtered to
`pchembl_value IS NOT NULL`. So it separates *weakly potent* from *never measured* — a compound
that was assayed and produced no measurable potency was never in the table to be counted. The
distinction the sentence promises needs a denominator the axis did not have.

`ligands/effort.py` supplies it. Measured on ChEMBL 37, bacterial targets:

| | |
|---|---|
| bacterial SINGLE PROTEIN targets | 1,073 |
| …with B/F activity at confidence ≥ 8, any relation | 987 |
| …restricted to `standard_relation '='` | **869** |
| …with any potency-measurable compound — *all the axis could see* | **661** |
| …with a compound at pChEMBL ≥ 6 | 430 |

**208 bacterial targets have nothing potency-measurable at all**, so they were invisible to this
axis entirely, and **226,038 compound-target pairs were assayed with no measurable potency**.

The useful part is that this recovers the axis's only real negatives: of **453 bacterial targets
with ≥ 10 compounds assayed, 148 (32.7%) never reached pChEMBL 6** (at ≥ 5 compounds, 251 of 603 —
41.6%). A protein somebody tried and failed to drug is a measured discouragement; a protein nobody
has opened is an open question. Before this they were the same zero.

### TWO denominators ship, and neither is merged

A denominator is only a denominator against a stated population, so `scratch/chembl_effort.tsv`
(869 bacterial targets × 17 columns) carries both:

| column | population | why it exists |
|---|---|---|
| `n_compounds_assayed` | **population-identical to `chembl.py`'s predicate** — B/F, confidence ≥ 8, `standard_relation = '='` | this is what `hit_rate` divides by |
| `n_compounds_assayed_any_relation` | the same, with the `standard_relation = '='` clause dropped | an inequality is still a measurement |
| `n_compounds_reported_inactive` | the explicit `>` rows | the non-binders, counted |

**An `IC50 > 100 µM` is the clearest statement in the database that a compound does not bind**, so
an *effort* count must include it. Dropping the `=` clause alone recovers **118 targets**, **6,682
`>` rows over 4,834 compounds** and **32,093 null-relation rows**. But a *ratio* needs its
denominator drawn from the numerator's own population — hence both columns, side by side, neither
merged into the other. `hit_rate` is **null, never 0**, where nothing was assayed.

**How the two stay population-identical.** `effort.py` does not restate `chembl.py`'s predicate; it
**imports `assert_version` and `_activity_where` from it by explicit spec load**, and **exits
non-zero if either clause it relaxes is no longer there**. That guard is not decoration: a
denominator computed over a different target population from the numerator is wrong in a way no
shape check could ever see — the columns line up, the ratio is finite, and it means nothing.

CLI: `--refresh` · `--dry-run` · `-q`. ~2 min, and it needs the 30.5 GB dump restored.

## REJECTED — a "ChEMBL precedence model" on embeddings

Proposed 2026-09-28 (*"train on chembl, X = proteins, y = num_ligands, and then learn a model"*)
and **not built**, by the project owner's decision once the numbers were in. Recorded with them so
it is not re-proposed — the `interpro2go` and Unknome precedents: a rejected alternative is written
down, not silently dropped.

1. **`y = num_ligands` is screening effort, near-definitionally.** Within bacterial ChEMBL,
   `rho(n_ligands, n_activities) = +0.985` and the median is **1.00 activity per compound**. A
   model of y is a model of how many analogues someone published. The repo's own numbers already
   said this from the other side: `ampC` is 12,671 compounds over 5,829 scaffolds, `folA` 443 over
   64 — a lead-optimisation series moves y by hundreds and the number of chemotypes by ~1.
2. **There are no zeros to learn from.** Every ChEMBL target has `y >= 1` by construction, so the
   model never sees a protein that resisted ligand discovery and cannot output a calibrated zero.
   Applied to a proteome that is ~96% never-screened, it would be extrapolating onto a support it
   has never seen. This is the OGEE positives-only trap (`CLAUDE.md`: 40 taxa with zero negatives
   because RB-TnSeq cannot see essential genes) in a different database.
3. **It approximates a lookup we already ship exactly.** `src/precedents.py` answers "does a
   similar protein have ligands" in ~1 s and *names* the donor target, organism and identity, at
   100% agreement with `chembl.py`. This is Unknome §5 reason (2), and it applies harder here
   because the lookup is already a canonical deliverable rather than a possibility.
4. **The honest baseline is not chance, and it is a lookup too.** Measured on our three proteomes
   against "has any measurable bacterial ligand" (now `n_measured_bacterial > 0`; base rate 3.1–4.1%): protein length alone scores AUROC
   0.66–0.71, **studiedness `n_papers_family` scores 0.83–0.86 / AP 0.15–0.17**. A model would have
   to beat a citation count that names its own donors. (`best_pident_bacteria` scores 0.999 and is
   **circular** — it IS the label's definition at the 40% floor. It is not a baseline.)
5. **Mean reversion points the wrong way for this consortium.** GraDi wants novel targets. A model
   trained on ChEMBL rewards the already-prosecuted families — gyrase, DHFR, PBPs, FabI, LpxC — and
   hands a genuinely unexplored envelope protein a middling score that reads as "moderately
   promising" rather than "nobody knows".
6. **Sample size, under the grouping this repo requires.** Pooled by OrthoFinder orthogroup, the
   positives are **200 distinct liganded families against 6,189 unliganded**, and 656 bacterial
   ChEMBL sequences collapse similarly (E. coli and M. tuberculosis alone are ~35% of them).

**What was built instead**, because the same investigation produced it: the effort denominator
(`ligands/effort.py`) and the four evidence tiers above. Both are exact, auditable and name their
sources. If the question is ever revisited, the version worth measuring is **"given a family was
screened, did it ever yield a potent ligand"** — that one has real negatives (of 453 bacterial
targets with ≥10 compounds assayed, 148 never reached pChEMBL 6) and is conditioned on effort by
construction, rather than being a popularity model wearing a ligandability label.

## Precedent: species-level exact, and what a zero means

**The exact count (`n_ligands` / `n_measured`) is species-level, not byte-level.** It was "identical sequence, or an accession
match", so a single substitution in another isolate demoted the same enzyme to a homolog and its
ligands left the exact count. It is now the UNION of three routes — accession, identical sequence,
and **same species at ≥ 95% identity** (two-word binomial, the same rule `chembl.py` uses for its
`species` bucket) — and `exact_route` records which fired.

Measured: E. coli **103 → 106** proteins with exact evidence, S. aureus **27 → 51**, Kp 5 → 6.
The case that makes it concrete:

| gene | before | after | identity |
|---|---|---|---|
| `def` (peptide deformylase) | **0** | **196** | 99.2% |
| `thyA` | 52 | 234 | 100% |
| `nfsA` | **0** | 14 | 98.3% |
| `lacZ` | 3 | 81 | 100% |

A peptide deformylase inhibitor programme was invisible to the exact count because ChEMBL's entry
is a different *E. coli* strain. **95% here is a claim about protein identity, not about how far
evidence travels** — `transfer_calibration.py` measured the potency conditional as flat from 25% to
100% (see *The 95/60/40 bands CANNOT be calibrated*), so no identity number on this axis is an
accuracy threshold.

Two spot checks, read from `evidence/precedents_full_ecoli.tsv` as shipped (measurable-potency
counts, exact / bacterial / human): **`folA` 517 / 595 / 0** — a heavily prosecuted antibacterial
target with no human liability at all — and **`clpP` 30 / 136 / 210 human at 56.3% identity**,
which is the degradation handle's selectivity problem arriving from a third independent direction
(stage 05's orthology and Part 1's selectivity figure are the other two).

### The table: potent counts and assayed counts

The deliverable asks two questions per scope, not one — **"did anyone find a sub-micromolar
binder"** and **"did anyone look"**:

| column | meaning |
|---|---|
| `n_ligands` | potent (pChEMBL ≥ 6) on **this protein** |
| `n_ligands_bacterial` | potent over the bacterial pool |
| `n_ligands_human` | potent over human targets — liability |
| `n_assayed` | compounds **assayed** against this protein, any outcome |
| `n_assayed_bacterial` | assayed over the bacterial pool |
| `n_assayed_human` | assayed over human targets |
| `best_pactivity_bacteria` | max pChEMBL over the bacterial pool |

Measured, proteins with at least one:

| | Kp | Ec | Sa |
|---|---|---|---|
| potent, this protein | 5 | 64 | 32 |
| **potent, bacterial** | **113** | **96** | **78** |
| any measurable potency (`n_measured_bacterial`, evidence table) | 180 | 160 | 119 |
| assayed, bacterial | 275 | 252 | 163 |

**The potent figures are the same 113 / 96 / 78 `chembl.py` reports**, so the two tables now agree
on what "has a ligand" means. They did not before: the counts were *any measurable potency*, which
includes a pChEMBL of 4.2 — a weak binder nobody would call a ligand.

**Every count is DISTINCT MOLECULES** — distinct `parent_molregno` over the **union** of the pool's
targets, never a sum of per-target counts. E. coli `folA` hits 6 bacterial targets and reads 443
where the per-target sum is 529; the 86 compounds tested against several homologs are counted once.
`gyrB`: union 295, sum 351.

**`*_bacterial` INCLUDES the exact match.** folA reads 388 potent on itself and 443 bacterial, and
the 443 *contains* the 388 — the same nesting `chembl.py` uses for direct ⊆ close ⊆ remote.
**Never sum the two.** Guaranteed per row:
`n_ligands ≤ n_ligands_bacterial ≤ n_measured_bacterial ≤ n_assayed_bacterial`.

**`best_pactivity_bacteria` is ChEMBL's `pchembl_value`**, renamed because the axis speaks of
activity rather than of one database's column name. Do not look for a `pactivity` field in ChEMBL.

**`n_assayed*` is NA, never 0, when the effort extract is absent** — it needs the 30.5 GB dump,
while the rest runs off cached extracts. A 0 would claim nobody ever assayed the protein.

### Reading a zero

A 0 in `n_ligands_bacterial` means one of three things, and `n_assayed_bacterial` says which:

| | |
|---|---|
| assayed > 0 | **a measured discouragement** — Kp `pyrH`: 158 compounds assayed against a 98.3%-identical target, none potent. Also Ec `polA` 36, `mrcB` 34, `phoA` 16 |
| assayed = 0, `n_targets_bacteria` > 0 | a homolog exists, nobody has opened it |
| `n_targets_bacteria` = 0 | nothing in ChEMBL within the floors — ~95% of each proteome |

A four-way `precedent_evidence` category used to encode this. **It was dropped**: "158 assayed, 0
potent" says strictly more than a label, and every category is recoverable from the counts plus
`n_targets_bacteria` in the evidence table.

**Those assayed numbers could not exist without `effort.py`'s 326 recovered sequences.**
`chembl.py` builds `chembl_targets.faa` after the pChEMBL filter, so a target whose compounds were
all assayed and none measurable has no sequence at all and DIAMOND cannot reach it. They carry no
ligand rows by construction, which is why the potent counts reproduce 113 / 96 / 78 exactly after
the change.

## Validated against the live ChEMBL API

Everything else in this axis descends from three cached extracts that `chembl.py` wrote from the
dump. If that SQL were subtly wrong, every downstream check would agree with it — the control
shares the extracts, the assertions share the code. `ligands/validate_api.py` asks a different
machine the same questions over HTTP.

**88/88 comparisons over 49 proteins match exactly**, all three species, spanning 1 to 12,438
compounds.

| round | what it checks | cases | result |
|---|---|---|---|
| exact | `n_ligands` against the resolved exact target | 30 | 30/30 |
| union | `n_ligands_bacterial` over every bacterial homolog | 29 | 29/29 |
| union | `n_measured_bacterial`, same pools | 29 | 29/29 |

Round 2 is the harder one: Kp `KPC-2` unions **52 targets**, `ctx-m-14` 49, `blaSHV-11` 39. A
dropped `tid` or a double-count would show there and nowhere else.

**The union is also checked against the per-target sum**, which it must stay below:

| | targets | union | per-target sum |
|---|---|---|---|
| Kp `KPC-2` | 51 | **276** | 365 |
| Kp `bla` | 36 | **269** | 327 |
| Ec `folA` | 6 | **443** | 529 |

**The API agrees with the union, never with the sum** — the distinct-molecule claim, verified from
outside the code that makes it.

Two things that make the comparison legitimate. The activity endpoint does **not expose
`confidence_score`**, which is acceptable *here and nowhere else*: `>= 8` removes 0 of 3,271,336
single-protein rows, so the gate is a no-op, and only SINGLE PROTEIN tids are compared on both
sides. And a target too large to page honestly is **skipped and named rather than truncated** —
E. coli `ampC` (12,438 potent compounds over 16 targets) is the only one, and its single-target
exact count was verified in round 1 regardless.

**The run found a defect worth recording.** Selecting test cases by `exact_target.notna()` picked
up 12,767 of 13,020 proteins. `pd.NA` written to TSV returns as `""`, and under `string` dtype an
empty string is a valid non-null value — so that column could not be used to ask "does this
protein have an exact match". The counts were never affected. `_coerce_precedents` now maps empty
to `pd.NA`, and the column agrees with `exact_route != "none"` at 12 / 182 / 69.

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

`data/processed/ligands/evidence/chembl_<species>.tsv`, one row per protein, keyed on
`uniprot_ac`. **Demoted from the task root on 2026-10-03** (owner's call): measured against
`ligands_<species>.tsv`, which is now the axis deliverable, two of its columns are exact
duplicates — `remote_n_compounds` == `n_ligands_bacterial` and `human_n_compounds` ==
`n_ligands_human`, same 113 Kp proteins, 5,286 vs 5,286 compounds, ρ 1.0.

**It is kept, and read from `evidence/`, for what precedents cannot express:**

| only here | why it matters |
|---|---|
| `*_n_scaffolds` | Kp's 5,286 potent compounds are **1,593 Murcko scaffolds** — 3.3 per scaffold. 50 analogues of one series is not 50 starting points, and precedents has no scaffold column at all |
| the bands split out | direct ≥95% (21 Kp proteins) · close ≥60% (79) · remote ≥40% (113); precedents collapses to exact + bacterial pool |
| `best_target` · `best_pident` · `best_organism` | which ChEMBL target the match actually came from |
| `allorg_*` | the unrestricted count, kept as the comparison that justified restricting to true Bacteria (424 vs a true 175) |

**What precedents has and this does not** is the denominator, `n_assayed*` — the axis's only real
negatives. Columns:

| column | meaning |
|---|---|
| `direct_hit` | a ChEMBL target at ≥95% identity exists (whether or not it has ligands) |
| `best_target`, `best_pident`, `best_organism` | the closest bacterial single-protein target |
| `n_compounds_tested` | parent compounds with any pChEMBL over the `remote` pool — **not** a screening denominator, see above |
| `direct_*`, `species_*`, `close_*`, `remote_*` | `n_compounds`, `n_scaffolds`, `best_pchembl` per bucket at pChEMBL ≥ 6 |
| `human_n_compounds/_n_scaffolds/_best_pchembl/_best_pident` | selectivity liability — never merged into a bucket |
| `allorg_n_compounds`, `allorg_n_scaffolds` | the literal all-organism count, so the Bacteria restriction is visible |
| `complex_n_compounds`, `complex_best_target` | component of a liganded complex |

**`evidence/`** — cite-or-check: `chembl_hits.tsv` (the raw DIAMOND join, unfiltered by band) ·
`scaffolds.tsv` (1.3M parent compounds → SMILES, InChIKey, molecule ChEMBL id, generic scaffold) ·
`organism_class.tsv` · `cutoff_sensitivity.tsv` · `control.tsv` · `precedent_control.tsv` ·
`precedents_full_<species>.tsv` · `transfer_calibration.tsv` · `chembl_effort_funnel.tsv` ·
`chembl_effort_relations.tsv` · `bindingdb_gain.tsv` · `manifest.tsv`.

**`scratch/`** — regenerable caches, safe to purge: `chembl_targets.tsv` (component↔target bridge)
· `chembl_ligands.tsv` (2.59M target×compound pairs) · `chembl_targets.faa` · `chembl_assayed.tsv`
· `chembl_effort.tsv` and `chembl_effort_targets.{tsv,faa}` · `transfer_pairs.tsv` ·
`bindingdb_ligands.tsv` · `bindingdb_targets.faa` · `hits/` · `smoke_*`.

The three `scratch/` extracts the precedents tool runs on total 82 MB and are what let it answer
without the dump — see Part 3.

| helper | returns |
|---|---|
| `load(species)` / `load_all()` | the deliverable |
| `load_targets` / `load_ligands` / `load_hits` / `load_scaffolds` | the supporting matrices |
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
python scripts/ligands/chembl.py --limit 200 --species ecoli   # -> scratch/smoke_*
python scripts/ligands/chembl.py                               # ~2 min with the dump
python scripts/ligands/effort.py                               # ~2 min, needs the dump
python scripts/ligands/transfer_calibration.py                 # ~2 min, DIAMOND from gradi-ortho
python scripts/ligands/bindingdb.py                            # ~8 min, first pass
python scripts/ligands/ligands.py --species kpneumoniae     # ~1 s per sequence, no dump
python scripts/plots/ligands.py
```

Full CLIs:

| script | flags |
|---|---|
| `chembl.py` | `--species` · `--pchembl` · `--threads` · `--limit` · `--refresh` · `--dry-run` · `-q` |
| `effort.py` | `--refresh` · `--dry-run` · `-q` |
| `transfer_calibration.py` | `--min-compounds 5` · `--threads` · `--dry-run` · `-q` |
| `precedents.py` | `--sequence` · `--fasta` · `--accession` · `--organism NAME` · `--species` (batch → `ligands_<sp>.tsv`, complete and canonical) · `--min-identity 40` · `--min-pchembl` · `-q` |

DIAMOND comes from `gradi-ortho` via `GRADI_DIAMOND_BIN` in all of them.

The 30.5 GB dump is needed only for the first run; **a cached re-run takes 0.3 min and needs no
database at all**. Restore it with `tar -xzf data/raw/other/chembl/chembl_37_sqlite.tar.gz -C
data/raw/other/chembl` — see `data/raw/other/chembl/SOURCE.md`. Extraction takes **30 s**, not the
~4 min first recorded.

### The archive is verified, and the version is now ASSERTED rather than declared

`chembl_37_sqlite.tar.gz` matches **EBI's published sha256** (`33c2037405…`, from
`releases/chembl_37/checksums.txt`) and that URL's `Content-Length` byte-for-byte, so what is on
disk is the genuine complete release and not a truncated or resumed transfer. **That had never been
checked** — the first run recorded byte counts only, and a byte count does not detect corruption.
Two gaps the check exposed, both closed:

1. **`CHEMBL_VERSION = "37"` was a bare literal while `find_db` globs `chembl_*.db`.** A different
   release extracted beside this one would have been consumed silently and labelled 37 in every
   manifest. **`assert_version()`** now reads the dump's own `version` table (filename as fallback)
   and exits non-zero on a mismatch.
2. **`SOURCE.md` pointed at the FTP `latest/` path**, which moves with every release. The pinned
   `releases/chembl_37/` URL replaces it.

**The cached extracts reproduce byte-for-byte from the restored dump.** `extract_chembl(refresh=
True)` into a temp directory returned all three sha256-identical: `chembl_targets.tsv` 9,347 × 9,
`chembl_ligands.tsv` 2,591,526 × 6, `chembl_targets.faa` 8,469 records. The dump self-identifies as
`ChEMBL_37` dated 2026-05-01 and its counts match the release notes. Independently, *before* the
restore, 4 targets were checked against live ChEMBL: `CHEMBL1293248` 24,681 activities and
`CHEMBL2390811` 8/8 agree **exactly**; the two deltas (`CHEMBL5465386` +18, `CHEMBL2026` +26) are
the `potential_duplicate` rows the extract drops on purpose.

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
  the map comes from `evidence/organism_class.tsv` and a missing or empty map exits non-zero; the
  count of unclassifiable organism strings (207) is printed every run.
- **ChEMBL's `version` table IS NOT ONE ROW — it holds 11, and `ChEMBL_37` is not first.** A
  `fetchone()` returns `Bioassay Ontology 2.0`, and `LIKE 'ChEMBL_%'` does not disambiguate either,
  because `ChEMBL_Structure_Pipeline 1.2.0` matches it. Only `ChEMBL_<digits>` exactly is the
  release. **The first `assert_version()` got this wrong and rejected the correct database** — it
  had passed a synthetic one-row fixture, i.e. it tested the assumption rather than the schema,
  which is the house rule about asserting on content wearing a different hat.
  Two upstream versions worth knowing, both read from that same table: **Swiss-Prot 2025_03**
  supplies `component_sequences` (the sequences `precedents.py` searches), and **RDKit 2022.09.4**
  did the salt stripping behind `molecule_hierarchy` — which is what the whole `parent_molregno`
  collapse rests on.
- **A denominator must come from the numerator's population.** `effort.py` relaxes exactly one
  clause of `chembl.py`'s predicate and imports the rest rather than restating it, exiting non-zero
  if either clause it relaxes has moved. Two counts drawn from different target populations line up
  perfectly in a spreadsheet and mean nothing.
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
- **Structure-based ligandability is now its own axis** — `docs/pockets.md`,
  `pockets_<species>.tsv`. Under a stricter, established drug-like definition and a binding-site
  coverage test it gives 608 Kp proteins with bacterial holo evidence (v1's 2,525 is not
  reconciled; see that doc). AlphaFill was measured there and left out.
- **The ribosome is missing** and the reason is structural, not incidental: its drugs are measured
  as MIC, which carries no pChEMBL. Given stage 04's top-100 is ribosome-heavy, a targeted route to
  ribosome-binding evidence matters more here than elsewhere. PDB co-crystals recover it only
  partly (9–11 ribosomal proteins per species in `docs/pockets.md`): the drugs mostly contact rRNA.
- **`allorg` versus `remote` is unexploited.** The gap between them is a measured statement about how
  much apparent ligandability comes from non-bacterial homologs, which is exactly the quantity v1's
  rat trap got wrong. Worth a figure.

---

# Part 3 — ligand precedent for an arbitrary sequence

`scripts/ligands/ligands.py` (CLI) and `src/precedents.py` (library) answer a different question
from Part 1. Part 1 asks *what does our proteome have*, needs the 30.5 GB ChEMBL dump, and produces
a per-species table. This asks *what about **this** sequence* — any sequence, in about a second,
from 82 MB of cached extracts.

```
n_ligands              POTENT (pChEMBL >= 6) on THIS protein -- accession, identical
                       sequence, or same species at >= 95%
n_ligands_bacterial    potent over the bacterial pool   (INCLUDES this protein)
n_ligands_human        potent over human targets        (LIABILITY, never summed in)
n_assayed              compounds ASSAYED against this protein, whatever the outcome
n_assayed_bacterial    assayed over the bacterial pool
n_assayed_human        assayed over human targets
best_pactivity_bacteria   max pChEMBL over the bacterial pool
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

**They count MOLECULES, not target-compound pairs.** `parent_molregno` is ChEMBL's
`molecule_hierarchy` parent, so salts are already collapsed; taking the DISTINCT set over the union
of every passing target means a compound tested against three homologs counts **once**. Measured:
E. coli `folA` reads 443 over 6 targets where the per-target sum is 529. Summing would inflate it,
and the wider the identity band the worse it gets — precisely the regime this tool is for. This
holds for the assayed side too, which is why `effort.py` emits `(target, compound)` pairs rather
than per-target totals.

**`n_ligands_human` is a liability, not a precedent.** A ligand-bearing human ortholog says the
fold is druggable *and* that hitting it may be dangerous. E. coli `clpP` is the case: **136
bacterial measurable against 210 human** at 56.3% identity, the human mitochondrial CLPP ortholog.
The CLI warns when human exceeds bacterial. Never sum them.

**The complex track is separate and is not inside `n_ligands_bacterial`.** E. coli `gyrB`: **666
single-protein measurable ligands, 1,412 complex**. Kp `A0A0H3H0Y6` (gyrA): **1,410 complex against
131 single.** DNA gyrase is a `PROTEIN COMPLEX` in ChEMBL and v1's single-protein-only rule made
GyrA/GyrB look unliganded. `n_ligands_bacterial_complex` is in the evidence table.

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
`evidence/chembl_<sp>.tsv` **at matched semantics** (pChEMBL ≥ 6, complex track included).

| species | both positive | precedents-only | chembl-only | agreement | count mismatches | total ligand delta |
|---|---|---|---|---|---|---|
| kpneumoniae | 113 | **0** | **0** | **100%** | 0 | 0 |
| ecoli | 96 | **0** | **0** | **100%** | 0 | 0 |
| saureus | 78 | **0** | **0** | **100%** | 0 | 0 |

*(These replace a "precedents-only 4 / agreement 99.9%" table that this page carried after the
agreement had already become exact — the adversarial-audit section below records the fix, and the
shipped `precedent_control.tsv` has read `agreement 1.0`, `n_count_mismatch 0`,
`total_ligand_delta 0` and `max_abs_delta 0` on all three species since. `matched_semantics` in
that file states the comparison: `pchembl>=6, SINGLE track both sides (chembl.py:452)`.)*

At *default* settings the two differ, and both reasons are design rather than defect: the tool
excludes the complex track, and counts every measurable ligand rather than only pChEMBL ≥ 6 — 66 of
67 default-only Kp proteins sit below 6. Restore both and they reconcile, which is what the control
records.

## Batch output

`--species` runs a whole anchor proteome through the same code path and writes
`ligands_<species>.tsv`, complete and in canonical row order. Measured:

| species | any exact | any bacterial | any human |
|---|---|---|---|
| kpneumoniae | 6 | 179 | 114 |
| ecoli | 112 | 159 | 106 |
| saureus | 27 | 117 | 66 |

Kp's exact count is low by construction — ChEMBL holds only 21 *K. pneumoniae* targets, and the
anchor is a dark TrEMBL proteome whose accessions are largely absent. That is the same fact that
makes homology transfer the whole game in Part 1.

Load through `src/ligandability.py` — **`load(species)`** for the deliverable, `load_all()` for
all three stacked, `load_full(species)` for the 19-column provenance view, and
`load_chembl(species)` / `load_chembl_all()` for the evidence table. Swapped on 2026-10-03 so
`load()` means the deliverable here as it does on every other axis; it used to return the ChEMBL
table. Note `load_ligands()` is unrelated — it is the raw `scratch/chembl_ligands.tsv` extract.

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
