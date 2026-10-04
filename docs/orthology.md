# Stage 05 — orthology

`scripts/orthology/orthofinder.py` · figures `scripts/plots/orthology.py` · loader `src/orthology.py` ·
outputs `data/processed/orthology/`

Two questions about the same protein pairs, answered by two methods and never merged: **are they
orthologs**, and **how similar are they**. All four proteomes in one run — 33,436 proteins.

| | proteins | in an orthogroup | with paralogs | →Kp | →Ec | →Sa | →Hs |
|---|---|---|---|---|---|---|---|
| *K. pneumoniae* | 5,728 | 4,322 | 2,175 | — | 3,784 | 1,716 | 1,184 |
| *E. coli* | 4,403 | 3,646 | 1,373 | 3,385 | — | 1,412 | 980 |
| *S. aureus* | 2,889 | 1,754 | 898 | 1,322 | 1,251 | — | 700 |
| *H. sapiens* | 20,416 | 16,282 | 15,791 | 1,555 | 1,499 | 1,298 | — |

**5,064 orthogroups**, ~410 of them spanning all four species. 10 min end to end.

## The deliverable is four columns

`orthology_<species>.tsv`, the three bacteria only, on the project owner's instruction
(2026-10-03):

| column | meaning |
|---|---|
| `has_human_ortholog` | the selectivity liability |
| `bacterial_panel_orthologs` | **a fraction, 0–1, not a count** — the share of the bacterial panel sharing this protein's orthogroup |
| `orthology_evidence` | 1–3, how well corroborated the two columns beside it are |

**The denominator is 28, not 26.** The 26 tier-C comparator proteomes plus the three bacterial
anchors, minus this protein's own species. It counts **species**, never proteins, so a paralog pair
does not inflate it. `bacterial_panel_size` ships in the dense table so the number stays
interpretable if the panel ever changes — 12 of 28 and 12 of 3 are different claims.

| | median `prop` | `has_human_ortholog` |
|---|---|---|
| Kp | 0.536 | 16.6% |
| Ec | 0.571 | 19.0% |
| Sa | **0.179** | 21.6% |

*S. aureus*'s low median is its Gram-positive isolation against a panel that is mostly
Gram-negative — biology, not a defect.

**A 0 is MEASURED.** OrthoFinder runs *de novo* on our own FASTAs and the stage exits unless every
protein is accounted for, so "no bacterial ortholog" is a finding, not a lookup miss.

### `orthology_evidence`, and why there is no `orthology_consensus`

Added **2026-10-04**, the axis's share of the project-wide pair in `src/consensus.py`.

**There is no consensus column, and "no magnitude" is only half the reason.** The two columns point
in *opposite* prioritization directions — you want **no** human ortholog and **broad** bacterial
conservation — so a mean of them is not a summary, it is a weighting, and CLAUDE.md forbids one
outright (*"NO COMPOSITE SCORE, EVER"*; three axes have already removed one). `has_human_ortholog`
is a boolean besides. The convention set by function and proteomes holds: **evidence always,
consensus where the axis has a single magnitude.**

```
3  placed, AND both OrthoFinder and RBH independently found a bacterial ortholog,
   AND the two methods do not conflict on the human call
2  placed by at least one grouping, but not corroborated
1  placed by NEITHER grouping -- nothing could be looked up, so both columns beside it
   are "could not look", not "looked and found nothing"
```

| | L1 | L2 | L3 |
|---|---|---|---|
| Kp | 291 | 2,639 | 2,798 |
| Ec | 110 | 1,537 | 2,756 |
| Sa | 216 | 1,875 | **798** |

"Placed" means an OrthoFinder orthogroup **or** an OrthoDB group (`orthodb_verdict` of
`assigned_by_sequence`/`assigned_by_uniprot`). OrthoDB is external and panel-independent, which is
exactly what OrthoFinder's *de novo* grouping is not.

**LEVEL 3 REQUIRES A POSITIVE FINDING, NOT ONLY A RELIABLE MEASUREMENT** (project owner,
2026-10-04), chosen over a reliability-only ladder with the cost in view: ***S. aureus* reaches 3
for only 798 proteins (27.6%)** because it is the lone Gram-positive among the anchors and both its
comparators are Gram-negative — the same structural effect behind its 0.179 median above. **So a
low level on Sa is partly its biology, not only our uncertainty.** That is the one way to misread
this column.

**Two things measured rather than assumed.** `in_orthogroup` **uniquely blocks 0 proteins at level
3 on all three species** — it is implied by the OrthoFinder bacterial term, and is load-bearing
only for level 1. `in_orthodb` is nearly redundant too (uniquely blocking 10 Kp / 3 Ec / 0 Sa). The
discrimination comes from the **RBH** term (uniquely blocking 301 / 189 / 218) and the **human
conflict** term (312 / 301 / 126); the two groupings' real job is separating level 1. And **RBH
corroboration reaches only the other two anchors, not the 28-species panel** — the pairwise DIAMOND
searches ran on anchors alone — so the test is *narrower* than the column it grades.

**The human-conflict term is what earns its keep.** `has_human_ortholog` is the union of the two
methods, and of Kp's 951 human calls only **501 are found by both** — 270 by OrthoFinder alone, 180
by RBH alone. The worked example is E. coli **`tufA`**: conserved in all 28 panel species
(`bacterial_panel_orthologs` 1.0) and still capped at **2**, because OrthoFinder calls a human
ortholog and RBH does not — EF-Tu against mitochondrial TUFM, exactly the 40–50% band where
ortholog-vs-paralog is genuinely ambiguous. It is demoted not for being poorly characterised but
because the *decision-critical* claim is contested, which is the behaviour wanted.

**OrthoDB CANNOT be a third opinion on the human call — verified, not assumed.** Bacterial groups
are `<n>at2` and human `<n>at2759`; the two id sets share **zero** members, and OrthoDB has no root
level spanning domains, so a bacterial protein is never searched against eukaryotic groups. The
axis has exactly **two** opinions on human orthology. `best_identity_human` is not a third either:
RBH is derived from the same DIAMOND search.

**Nothing in this axis is an experiment**, so a 3 is not experimental corroboration and a 1 means
"could not look" rather than "not yet measured" — the pattern `function_evidence` set.

Every condition is kept per protein in **`evidence/evidence_audit.tsv`**, because the ladder
collapses five booleans into one integer and a bare 2 never says which half failed.

**Polarity control**: median `bacterial_panel_orthologs` rises with the level — Kp 0.00 / 0.32 /
0.64, Ec 0.00 / 0.36 / 0.64, Sa 0.00 / 0.11 / 0.86 — and **level 1's median is exactly 0.00**,
which is what makes "could not look" legible.

**Human has no deliverable** — the panel columns are bacterial — and `load()` refuses it by name
rather than returning something misleading. Use `load_dense("human")`.

**Everything else moved to `evidence/orthology_<species>.tsv`**, 29 columns, via `load_dense()`:
`orthogroup`, `in_orthogroup`, `orthogroup_size`, `n_paralogs`, `searched`, the five
per-target-species columns for each of the four proteomes, `n_bacterial_orthologs` and
`bacterial_panel_size`. Three consumers read it there: `essentiality/predict.py` (orthogroup for
paralog grouping), `studiedness/transfer.py` and `plots/orthology.py`.

## The 26-species panel, and the two columns it is for

`--panel full` adds tier C — 26 curated bacterial comparators, pinned in
`src/proteome_registry.tsv` — to the OrthoFinder run. **30 proteomes, 140,396 proteins, 47 min.**
The comparators inform the orthogroups and nothing else: no tables of their own, and pairwise
DIAMOND stays on the four anchors (16 searches, not 900).

Two columns in `orthology_<species>.tsv`:

| column | meaning |
|---|---|
| `has_human_ortholog` | union of OrthoFinder and RBH, unchanged definition |
| `n_bacterial_orthologs` | how many of the **28** bacterial proteomes share this protein's orthogroup |
| `bacterial_panel_size` | 28 — so 12-of-28 cannot be misread as 12-of-3 |

`n_bacterial_orthologs` counts **species, not proteins**: a paralog pair does not make a protein
more conserved. It is **orthogroup-based rather than the union** used by `has_human_ortholog` —
RBH against 26 comparators would cost 156 further DIAMOND searches (~5 h) for a second opinion on
the same question.

**Validated against biology**, which is the check that matters for a conservation column:

| | ftsZ | gyrB | dnaA | secA | rpoB | clpP | lacZ |
|---|---|---|---|---|---|---|---|
| Kp | 28 | 28 | 28 | 28 | 27 | 16 | 5 |
| Ec | 28 | 28 | 28 | 28 | 27 | 16 | 8 |
| Sa | 28 | 28 | 28 | — | 27 | 7 | — |

Core essential genes sit at the panel maximum; a metabolic gene like `lacZ` does not. **S. aureus's
median is 5 against Kp/Ec's 15–16** — correct, not a defect: it is the only Gram-positive in a
mostly Gram-negative panel.

### Widening the panel moves two things at once, and human goes DOWN

Written to `evidence/panel_expansion.tsv` every run, because a shipped column shifting silently is
the thing to avoid:

| | in_orthogroup | has_human_ortholog |
|---|---|---|
| Kp | 4,322 → **5,289** (+967) | 1,184 → **951** (−233) |
| Ec | 3,646 → **4,247** (+601) | 980 → **838** (−142) |
| Sa | 1,754 → **2,262** (+508) | 700 → **624** (−76) |
| human | 16,282 → 16,507 (+225) | — |

The coverage gain is the point: **+2,301 anchor proteins gained an orthogroup**, the unassigned set
falling 7,432 → 5,131.

**The human drop was predicted to be a regression, and that prediction was wrong.** "Recall rises
with panel size" is a **bacteria↔bacteria** statistic (v1's 55.5% Kp↔Ec at 25 species vs 44.8% at
4); this panel adds 26 bacteria and **no eukaryotes**. Better resolution on the bacterial side lets
OrthoFinder separate orthologs from out-paralogs, so marginal bacteria↔human calls are withdrawn.

Refinement rather than breakage, measured:

| human identity | Kp survive | Ec | Sa |
|---|---|---|---|
| ≥ 60% | **13/13** | 11/12 | **3/3** |
| ≥ 50% | 56/60 | 57/59 | **23/23** |
| ≥ 40% | 195/224 | 184/206 | 118/128 |

The losses concentrate at 40–50%, where ortholog-vs-paralog is genuinely ambiguous — median
identity 33.5% among the kept against 28.8% among the dropped.

### The run that justifies keeping both methods

Kp `clpP` → human mitochondrial CLPP at 56.3% identity is documented here as an ortholog by *both*
methods. After the expansion OrthoFinder calls it **`of=0`**, and only RBH still finds it — so
`has_human_ortholog` stays `True` solely because the axis takes the union. At ≥50% identity RBH
rescues **9 Kp proteins** OrthoFinder misses (9 Ec, 4 Sa). Had either method been adopted alone,
this panel change would have silently dropped real human liabilities.

## Coverage: there is no "not in the database" failure here

The axis cannot lose a protein to a lookup miss, because OrthoFinder runs **de novo on our own
FASTAs**. Verified: `orthology_<sp>.tsv` has exactly one row per proteome protein with `searched`
true — 5,728 / 4,403 / 2,889 / 20,416, **100%**. Compare the lookup routes on the same proteins:
OrthoDB reaches 92.7% of Kp, eggNOG 91.5%.

A protein with **no orthogroup is a measured singleton, not a gap**. Before the expansion, of Kp's
1,406 such proteins only **214 (15.2%)** had any cross-species DIAMOND neighbour at all; the other
85% genuinely had no detectable homolog among three bacteria. Those ~15% are most of what the wider
panel recovered.

**A trap in the loaders, found twice.** `orthodb_og_domain.notna()` used to report 5,728 of 5,728
on Kp — but 416 of those are the empty string, so the real coverage is 5,312 (92.7%). `pd.NA`
written to TSV returns as `""`, and under `string` dtype that is a valid non-null value. The same
defect was fixed in `src/ligandability._coerce_precedents` for `exact_target`, so it is a property
of this repo's TSV round-trips. Both loaders now map empty to `pd.NA`. Relatedly,
`src/orthology._read` coerced numerics from a hand-kept column list, so `n_bacterial_orthologs`
came back as **text** and `<= panel_size` raised `TypeError`; the rule is now the `n_` prefix, so a
new count column cannot be forgotten.

## Why this stage exists in this shape

v1 built orthology twice and its retrospective records the outcome as **trap 1**, its headline
orthology bug:

> The 03a kp→ec ortholog table has entirely empty `pident` / `coverage` / `bitscore`. It is
> OrthoFinder orthogroup membership only. Any axis that thresholds transfer on percent identity
> **silently drops everything**.

So this stage ships a discrete matrix *and* a continuous one. No consumer ever has to threshold on a
column nobody populated.

## A zero here is a measured zero

The requirement: when `is_ortholog` is 0, that must mean *we looked and they are not orthologs*,
never *this accession was not in the database*.

**It rules out accession LOOKUP, not databases.** The distinction is the whole thing, and an earlier
draft of this document got it wrong. A lookup returning nothing because the accession is not
cross-referenced is an *unknown*. A search returning nothing because the sequence matched nothing is
a *measured* negative. On this anchor the two diverge violently:

| route | Kp coverage |
|---|---|
| UniProt's eggNOG xref (lookup) | **0.00%** |
| OrthoDB REST by accession (lookup) | **0** — even for reviewed accessions |
| v1's OrthoDB gene-symbol pivot | 18.3% — i.e. Kp's 18.4% native gene-name coverage |
| eggNOG-mapper, run over the sequences (search) | **91.5%** |
| DIAMOND vs OrthoDB, one assembly (search, v11) | 74.6% |
| DIAMOND vs **all** of OrthoDB (search, v12, `--reps 20`) | **92.7%** |

So the rule is *search, never look up* — CLAUDE.md's standing convention — and both group sources
here obey it: OrthoFinder de novo, and OrthoDB by sequence (see *Absolute groups* below).

### stage 02's eggNOG groups are a measured source too — the retraction, stated once

**eggNOG-mapper's zeros are measured, not unknown, and anything saying otherwise is the retracted
version.** An earlier draft of this document ruled eggNOG out for "9–11% unknowns" and that
reasoning was wrong: stage 02 *ran* emapper over every protein, so a protein with no orthologous
group is a searched negative exactly like an OrthoFinder unassigned gene. The unknown problem
belongs to the *xref lookup* route — UniProt's eggNOG xref is 0.00% on Kp while emapper reached
91.5% over the same proteins. Conflating search with lookup wrongly eliminated the one source that
was already absolute *and* already computed. (Recorded again in the 2026-09-03 run log, correction
1; this is the considered position and it stands.)

Per-species, from stage 02's own run log — groups assigned by emapper:

| | Kp HS11286 | Ec K-12 | Sa NCTC 8325 | human |
|---|---|---|---|---|
| eggNOG OG (search) | 91.5% | 96.3% | 89.1% | **not run** |

**So the reason to run OrthoFinder de novo is not that eggNOG's zeros are soft.** Two real ones
survive. (1) **Coverage of all four species.** Stage 02 is the three bacteria only — COG2024 has no
eukaryotes and human was excluded there by construction — so eggNOG answers nothing about the human
column, which is the whole off-target question this axis exists to price. (2) **A measured zero for
every protein, by construction rather than by coverage.** OrthoFinder either assigns an input
protein an orthogroup or names it in `Orthogroups_UnassignedGenes.tsv`, and the stage exits
non-zero unless those two account for each proteome exactly; emapper's 89–96% leaves 4–11% that is
*searched* but not *reconciled against the input set*, which is a weaker guarantee even though it
is the same kind of evidence. And neither source resolves pairwise orthologs or paralogs at all —
only OrthoFinder's `Orthologues/` does.

**OrthoFinder run de novo on our own four FASTAs satisfies the requirement by construction.** Every
input protein is either assigned an orthogroup or named in `Orthogroups_UnassignedGenes.tsv`; both
are measured outcomes. Measured, first run:

```
kpneumoniae      5728 /   5728 placed by OrthoFinder
ecoli            4403 /   4403 placed by OrthoFinder
saureus          2889 /   2889 placed by OrthoFinder
human           20416 /  20416 placed by OrthoFinder
```

The stage **exits non-zero** unless that is exact, and every accession is a query in all four DIAMOND
searches (`searched` records it). A *sparse* matrix cannot express a zero at all — absence of a row
*is* the zero — which is why the dense per-protein table exists.

## Method

### Discrete — OrthoFinder 3.1.5, de novo, all four species jointly

A **full run, not `-og`**: `-og` stops after orthogroups and never writes `Orthologues/`, which is the
only place pairwise orthologs appear. FASTAs are written from the `sequence` column of
`proteome_<species>.tsv` rather than from `data/raw/`, so they cannot drift from the table everything
else joins on.

`Orthologues/` files are **many-to-many** — a cell lists several genes on each side. That is
co-orthology, and it is precisely what RBH structurally cannot represent; every listed combination
becomes a pair.

### Continuous — DIAMOND 2.2.1, 16 ordered species pairs

```
diamond blastp --very-sensitive --evalue 1e-3 --max-target-seqs {5, or 6 on the self-pair}
  --outfmt 6 qseqid sseqid pident ppos length qlen slen qcovhsp scovhsp evalue bitscore
```

**`--very-sensitive`, against v1's plain defaults, and the reason is the direction of the error.**
Under-detecting human homology makes a target look *more selective than it is* — which advances a
candidate that in fact has a human ortholog. It paid off exactly as predicted: **Kp proteins with any
human hit went from v1's 732 (12.8%) to 1,430 (25.0%)**. v1 was calling ~700 Kp proteins human-free
that have a detectable human homolog.

Identity and coverage are **columns, never filters**, so a missing row means "no detectable homology
at e-value 1e-3" — a far stronger claim than v1's silently-thresholded table.

RBH is derived from the same searches (rank 1 is the best hit, so no extra pass is needed), with no
identity or coverage floor, matching v1's `03c` so the control is like-for-like.

## Outputs

```
data/processed/orthology/
  orthologs.tsv            29,844 rows   sparse: pairs called by EITHER method
  neighbors.tsv           146,723 rows   sparse: top-5 per (protein, target species)
  orthology_<species>.tsv   4 files      DENSE: one row per protein -- where a 0 is readable
  orthodb_<species>.tsv     4 files      DENSE: absolute groups (see below)
  evidence/   orthogroups.tsv, method_disagreement.tsv, control.tsv, manifest.tsv,
              orthodb_manifest.tsv, orthodb_transfer_audit.tsv, orthodb_decoy_hits.tsv
  scratch/    orthofinder/ (the native run), hits/, fasta/, and the OrthoDB slices --
              purgeable, and the reason a re-run is cheap rather than free
```

`orthologs.tsv` — `query_ac · query_species · target_ac · target_species ·
is_ortholog_orthofinder · is_rbh · same_orthogroup · orthogroup`

`neighbors.tsv` — the above plus `rank · pident · ppos · alnlen · qcov · scov · evalue · bitscore ·
bitscore_norm`. Self-hits are dropped, so the within-species block is a protein's nearest
**paralogs**. `bitscore_norm` divides by the query's self-bitscore, making it comparable across
proteins of very different length.

`orthology_<species>.tsv` — `uniprot_ac · orthogroup · in_orthogroup · orthogroup_size ·
n_paralogs · searched`, then **three counts per target species** (`n_orthologs_<sp>` either,
`n_orthologs_of_<sp>` OrthoFinder, `n_orthologs_rbh_<sp>` RBH) plus `best_identity_<sp>`,
`best_bitscore_<sp>`, and `has_human_ortholog`.

Three counts and not one because **the two methods respond to different things** — see the panel-size
effect below. A single merged count would hide which of the two moved.

| helper | returns |
|---|---|
| `load_orthologs()` | the discrete matrix |
| `load_neighbors()` | the continuous matrix |
| `load(species)` / `load_all()` | the dense per-protein table(s) |
| `orthologs_of(ac)` / `neighbors_of(ac)` | one protein's rows |
| `load_orthogroups()` / `load_disagreement()` / `control()` / `manifest()` | the `evidence/` tables |

**`SPECIES` in `src/orthology.py` is all FOUR** — `("kpneumoniae", "ecoli", "saureus", "human")`.
Human is in scope on this axis, unlike `src/function.py`, where COG2024 has no eukaryotes to give
it an honest letter. Anything iterating `SPECIES` across axes must not assume the two agree.

## What the two methods disagree about

**35.2% of ortholog pairs are called by both**; 15,944 by OrthoFinder alone, 3,398 by RBH alone.
v1 computed both on the Kp–Ec pair (3,179 vs 3,003) and **never compared them** — a free validation
set left unused. Disagreements are written to `evidence/method_disagreement.tsv` (19,342 rows).

Agreement falls with evolutionary distance — Ec–Sa 50%, Kp–Ec 45%, and only 21–24% for the human
pairs. That independently vindicates v1's methodological note, *"orthology inference is unreliable
across that distance"*, while still giving both readings for human rather than picking one.

### OrthoFinder's recall depends on panel size; RBH's does not

The most important thing to know before quoting a number from this stage. Measured here:

| Kp proteins with an *E. coli* ortholog | v1 | v2 |
|---|---|---|
| OrthoFinder | 3,179 (55.5%) — **25-species run** | 2,568 (44.8%) — 4-species run |
| DIAMOND RBH | 3,003 | **3,074 (53.7%)** |
| either | — | 3,784 (66.1%) |

**RBH reproduces v1 within 2.4%**, and RBH median identity lands at 85.7% against v1's 86.0%. That is
the like-for-like comparison and it passes.

**OrthoFinder does not, and should not.** v1's widely-quoted 55.5% came from a *25-species* run;
orthogroup inference needs phylogenetic signal, and a small panel starves it. A 2-species smoke test
made this unmistakable: OrthoFinder emitted just **631** ortholog groups for Kp×Ec alone, against
2,568 Kp proteins in the 4-species run. **So the 55.5% figure is not a target this stage can hit at
four species, and falling short of it is not a defect.** If OrthoFinder recall matters more than run
time, add the tier-C panel — that is the lever, and it is cheaper than it looks: **OrthoFinder 3's
incremental mode (`--assign` / `--core`) adds species against an existing core without recomputing
it**, which is exactly why `scratch/orthofinder/Results_Sep03/` is kept rather than purged with the
rest of `scratch/`.

## Absolute groups — OrthoDB v12.2

`scripts/orthology/orthodb.py` · outputs `orthodb_<species>.tsv`

OrthoFinder's orthogroups are **panel-dependent** (above), which makes them the wrong instrument for
"how broadly conserved is this target". This adds the other kind: groups defined once, elsewhere,
over a fixed set of genomes — **13.0M groups over 990 levels, 17,551 bacterial and 5,952 eukaryotic
species**. They sit *alongside* OrthoFinder's, never instead: OrthoFinder remains the only source of
pairwise ortholog/paralog resolution.

**v12 replaces v11 outright; the two id spaces must never be joined.** OrthoDB's own README states
that an "OG unique id (**not stable and re-used between releases**)", so relabelling a v11 id as v12
would be silently wrong. Every row carries `orthodb_version`, and the v11 artifacts were deleted
rather than kept alongside. Downloads are at
`https://data.orthodb.org/v12/download/odb_data_dump/` with `odb12v2_*` naming — note that
`/v12/download/odb12v0_*` 404s, which makes the path easy to miss.

### Coverage

| species | v12 assigned | v11 | |
|---|---|---|---|
| *K. pneumoniae* | **92.7%** (5,312 / 5,728) | 74.6% | **+18.1** |
| *E. coli* | **95.5%** (4,204 / 4,403) | 88.8% | +6.7 |
| *S. aureus* | **91.8%** (2,652 / 2,889) | 88.3% | +3.5 |
| *H. sapiens* | **96.0%** (19,600 / 20,416) | — | — |

**100% *assigned* is not achievable from any absolute database, and is not the goal.** Forcing it
means loosening thresholds until noise fills the gap — the failure stage 02 measured for COG, where
the shuffled-decoy hit rate went 0.4% → 84.3% as the e-value relaxed.

**100% *verdict* coverage is guaranteed.** Every accession leaves with one of a closed set that never
contains "unknown" — 33,436 / 33,436, and the stage exits otherwise:

| verdict | meaning | n |
|---|---|---|
| `assigned_by_uniprot` | the protein **is** an OrthoDB gene; its own groups, no inference | 8,865 |
| `assigned_by_sequence` | groups transferred from homologs by bitscore vote | 22,903 |
| `no_group` | searched, no group — **a measured zero** | 1,668 |

### What lifted coverage: three changes, each measured

1. **Search all of OrthoDB, not one assembly.** v11 DIAMONDed each proteome against *its own*
   OrthoDB organism — for Kp a single 4,975-protein strain, which is why it capped at 74.6%. v12
   searches representatives of every domain-level group (`--reps 20`, 11,858,928 sequences).
2. **Take groups at every level, not just the domain.** Plenty of families have a group at
   Enterobacteriaceae and none at Bacteria. Measured: **510 Kp proteins** had a clean hit and still
   returned `no_group` while only `at2` groups were considered. Worth +9.4 points on its own.
3. **But transfer only at levels our own lineage passes through.** A Kp protein matching a *Bacillus*
   gene whose only group is Bacillus-level has learned nothing about Kp; (2) alone inflates coverage
   with meaningless assignments. Lineages come from `level2species` column 4 and are pinned in
   `src/orthology.ORTHODB_LINEAGE`, e.g. Kp `{2, 1224, 1236, 91347, 543, 570, 72407}`.

### The identity floor is decoy-calibrated, not inherited

33,436 composition-preserving shuffles of our own sequences were searched against the same database.
At 40%/50% they produced **0 hits**; at 25%/50%, **2**; at 20%/20%, 5. Meanwhile the 40% floor cost
~10 points of coverage on Kp and **~22 on *S. aureus***, and E. coli correctness is **flat** across
floors (73.3% at 40/50 vs 72.6% at 20/50). So the floor is **25% identity / 50% coverage**.

CLAUDE.md's ≥40% rule is for annotation **transfer**, a stricter task than group assignment. Do not
conflate them.

### Accession joins into OrthoDB do not work — now measured a fourth time

| route | Kp | Ec | Sa | human |
|---|---|---|---|---|
| UniProt's `xref_orthodb` | **0.0%** | 96.4% | 99.4% | — |
| OrthoDB's own mapped `uniprot_id` | **0.0%** | 73.6% | 93.2% | **17.0%** |
| **DIAMOND vs all of OrthoDB** (`--reps 20`) | **92.7%** | 95.5% | 91.8% | 96.0% |

**HS11286 is absent from OrthoDB entirely** — taxid 1125630 is not in `species.tab`, and the only
*K. pneumoniae* organism is `72407_0`. No *Klebsiella* proteome in UniProt carries an OrthoDB xref at
all (MGH 78578 is also 0.0%). So every species is assigned by sequence, the standing rule.

### Filter on `orthodb_confidence`, never on identity

`orthodb_confidence` is the winning group's share of total bitscore at the domain level. On the
E. coli control it is the only signal that separates:

| vote margin | n | correct | | best-hit identity | n | correct |
|---|---|---|---|---|---|---|
| **>0.9** | 1,178 | **98.0%** | | 95–100% | 214 | 79.0% |
| 0.75–0.9 | 265 | 90.9% | | 60–80% | 735 | 72.8% |
| 0.5–0.75 | 772 | 71.2% | | 35–45% | 415 | 45.1% |
| <0.5 | 1,397 | 41.4% | | 25–35% | 179 | 38.0% |

Median confidence: Ec **1.00**, Sa **1.00**, human 0.76, **Kp 0.53** — Kp has no OrthoDB
representation, so all **5,312** of its assignments are sequence-tier (at `--reps 20`; the shipped
`orthodb_kpneumoniae.tsv` is 5,312 `assigned_by_sequence` + 416 `no_group` = 5,728).

### The accuracy number, and why the obvious one is wrong

Sequence-tier assignment scores **75.4% exact-id agreement** with OrthoDB's own answer on E. coli.
That figure is an **underestimate**, because **OrthoDB maintains parallel Bacteria-level groups for
the same family**. Verified against its own v12 API:

- `5287828at2` "Chromosomal replication control, initiator DnaA" **and** `9807019at2` "chromosomal
  replication initiation protein A" — both real, both Bacteria-level, both *dnaA*
- `5287163at2` "DNA gyrase subunit B" **and** `9802808at2` "DNA gyrase subunit b"

**57.1% of the apparent errors are a parallel group carrying the same function name**, so corrected
agreement is **89.8%**. Two consequences: never compare two proteins by group id alone — check the
name — and low confidence does not imply wrong. Kp's `clpP` is correctly `9802800at2`, agreeing with
all three species *and* with v11, at confidence 0.219.

**This depresses the shared-core statistic.** Bacteria-level groups shared by all three bacteria:
**384** (Kp 4,231 · Ec 3,222 · Sa 2,088; union 6,816), against v11's 629 — not because conservation
changed, but because parallel groups split Kp/Ec from Sa. Use group *names*, or `neighbors.tsv`, for
core-genome questions.

### Verified against three sources that did not feed the assignment

1. **UniProt's own `xref_orthodb`** (never used by the script — it reads OrthoDB's `genes.tab`
   instead). Any-shared-group agreement, name-corrected for parallel groups: *E. coli* **94.8%**,
   *S. aureus* **99.4%**. **Human scores 0%** and that is UniProt's problem, not ours: its human
   OrthoDB ids come from a different release. Checked by hand — `A0A1B0GTW7` is CIROP, alternative
   name *leishmanolysin*, and our group is "leishmanolysin"; `A0JNW5` is BLTP3B, formerly
   *UHRF1BP1L*, and our group is "UHRF1-binding protein 1". Both correct; UniProt's xref named
   unrelated proteins.
2. **Group name vs the protein's own UniProt name** — independent of OrthoDB ids entirely. Human
   **84.0%**. The bacteria score 58–62%, but that is the *metric's* ceiling, not our error: for
   Ec/Sa the `assigned_by_uniprot` tier **is** OrthoDB's own answer with no inference from us, and it
   scores only **56–57%**, while our sequence tier scores **62–67%** — at or above it.
3. **OrthoDB's own server-side search** (`/v12/blast`, RapSearch2 — different algorithm, different
   machine), on 12 random Kp sequence-tier calls. Six matched on a function word; of the rest, three
   were ours being vaguer but correct (`pbpG` → "penicillin-binding protein", `betA` →
   "oxidoreductase", `bamE` → "membrane protein") and **three were genuinely wrong** (`cysC` assigned
   "translation elongation factor", `yehS` "Alpha-D-phosphohexomutase", `ycaP` "Beta-ketoacyl-").
   A ~25% error rate on Kp's hardest tier. **That sample was taken at `--reps 5`**; the check was
   repeated after the switch to 20 (below).

**The weak corner, and it is not where I first guessed.** Two of those three Kp errors carried
`orthodb_confidence = 1.0`, which looked like the metric failing — a single candidate group takes
100% of the bitscore share trivially. Stratifying on E. coli shows it does **not** fail:
`conf>0.9 & n_cand==1` is **95.0%** accurate and `conf>0.9 & n_cand>=2` is **98.0%**. What *is* weak
is **one candidate at low identity**: `n_cand==1` scores **80.0%** at 25–45% identity against 97–99%
above 45%. So read `orthodb_confidence` together with `orthodb_n_candidate_ogs` **and**
`orthodb_match_pident`; a lone candidate found at 30% identity is the one shape to distrust.

### Why the rest have no group, and which zero is stronger

| | `no_sequence_match` | `orthodb_singleton` |
|---|---|---|
| *K. pneumoniae* | 416 | 0 |
| *E. coli* | 140 | 59 |
| *S. aureus* | 72 | 165 |
| *H. sapiens* | 731 | 85 |

`orthodb_singleton` is the **stronger** zero: the protein *is* in OrthoDB and OrthoDB places it in no
group at any level. `no_sequence_match` is weaker — nothing in OrthoDB matched. Kp shows 0 singletons
because it has no OrthoDB gene to be a singleton *of*. The unassigned are short: *S. aureus* median
**52 aa** (82% under 100 aa), *E. coli* 92 aa — small hypothetical ORFs.

### OrthoDB has no root level

Its top levels are Bacteria (2), Archaea (2157), Eukaryota (2759) and Viruses (10239) — **none
spanning domains**. So **nothing in this table answers "does this target have a human ortholog"**; a
bacterial group is `<n>at2` and human's `<n>at2759`, and they are not comparable. That question is
answered twice already, by `neighbors.tsv` and by OrthoFinder's joint run.

### Columns

`uniprot_ac · orthodb_version · orthodb_og_domain · orthodb_domain_level · orthodb_og_narrow ·
orthodb_narrow_level · orthodb_n_levels · orthodb_og_name · orthodb_verdict ·
orthodb_no_group_reason · orthodb_match_pident · orthodb_confidence · orthodb_n_candidate_ogs ·
orthodb_gene_id`

Two levels are pivoted into columns — the **domain** and the **narrowest** assigned level.
`scratch/orthodb_groups_long.tsv` keeps all of them: **210,707 rows over 22 levels and 153,669
distinct orthogroups**. Load with `load_orthodb`, `load_orthodb_all`, `load_orthodb_long`;
`orthodb_confidence` and `orthodb_match_pident` come back as floats, so the documented filter works.

### Cost

The four small tables are kept (135 MB). `genes.tab` (4.5 GB), `OG2genes.tab` (4.5 GB, streamed
**twice** — once for representatives, once for all levels of the genes actually matched) and
`aa_fasta` (37 GB) are **stream-filtered and never stored**. Measured throughput 16.6 MB/s, so a cold
run is ~50 min of transfer, ~31 min of DIAMOND `--very-sensitive`, and 34 s of `makedb`.

### Running it

```bash
~/miniconda3/envs/gradi/bin/python scripts/orthology/orthodb.py --dry-run
~/miniconda3/envs/gradi/bin/python scripts/orthology/orthodb.py
~/miniconda3/envs/gradi/bin/python scripts/orthology/orthodb.py --reps 5    # faster, less accurate
```

Flags: `--species` · `--reps 20` · `--sensitivity very-sensitive` · `--max-targets 25` · `--threads` ·
`--refresh` · `--dry-run` · `-q`. Every slice is cached, so a re-run that only changes the assignment
logic is ~2 min. DIAMOND comes from `gradi-ortho`; `GRADI_DIAMOND_BIN` overrides.

### Traps specific to this script

- **The v12 download path is not where you would guess.** `/v12/download/odb12v0_*` 404s; the real
  location is `/v12/download/odb_data_dump/` with `odb12v2_*` naming.
- **`defaultdict` membership tests create keys.** `kept[og] < reps and og in dom` silently
  materialised an entry for all 13.0M groups before the domain check ran. Order the conditions the
  other way.
- **Restricting to domain-level groups looks like missing data.** It cost 510 Kp proteins that had
  perfectly good hits — the symptom is `no_group` on proteins whose DIAMOND hit passed every floor.
- **An OG id encodes its own level** (`9781621at2` → 2), so ranking specificity needs no lookup over
  the 13.0M-row `OGs.tab` — only `levels.tab` for the species counts.
- **`orthodb_confidence` must load as float.** `src/orthology._read` reads everything as `str` by
  default, so the documented `> 0.9` filter raised `TypeError` until the numeric columns were
  declared.


## Running it

```bash
~/miniconda3/envs/gradi/bin/python scripts/orthology/orthofinder.py
~/miniconda3/envs/gradi/bin/python scripts/orthology/orthofinder.py --species kpneumoniae ecoli --top-k 3
~/miniconda3/envs/gradi/bin/python scripts/orthology/orthofinder.py --refresh
~/miniconda3/envs/gradi/bin/python scripts/plots/orthology.py
```

Flags: `--species` · `--top-k 5` · `--evalue 1e-3` · `--sensitivity` · `--threads` · `--refresh` ·
`--dry-run` · `-q`. Without `--refresh` the OrthoFinder run directory and every DIAMOND hit file are
reused; OrthoFinder is the expensive step.

The stage runs in **`gradi`** and shells out to **`gradi-ortho`** (osx-64, Rosetta) for both tools.
`GRADI_DIAMOND_BIN` (a *directory*) and `GRADI_ORTHOFINDER_ENV` (an *env root*) override.

### Guards

- **Coverage** — assigned ∪ unassigned must equal each proteome exactly. This is the guard that makes
  a zero a measured zero; without it the matrix is uninterpretable. Exits non-zero.
- **No self-hits** in `neighbors.tsv`, and `rank` never exceeds `--top-k`. Exits non-zero.
- **RBH symmetry** — `is_rbh(a,b) == is_rbh(b,a)`. Exits non-zero.
- **No within-species OrthoFinder call** — orthology is a between-species relation, so a
  within-species `True` would be a category error. Measured 0. Exits non-zero.
- **v1 reproduction** — written to `evidence/control.tsv` every run, reported not enforced, because
  the OrthoFinder half is legitimately panel-dependent (above).

## Traps

- **OrthoFinder exits 0 when its dependency check fails.** If `diamond` is not on `PATH` when
  OrthoFinder starts, it prints an ERROR block, writes an empty `Results_<date>/` and **returns
  success**. The return code is not evidence; the results-directory check is what catches it. This
  bit during development: `ensure_diamond()` was called *after* the OrthoFinder step, so the PATH
  injection came too late. It now runs first.
- **The OrthoFinder launcher cannot be called directly.** Its shebang is `#!/usr/bin/env python3`,
  which on this machine resolves to the unrelated `ersilia` env, so `orthofinder -h` dies with
  `ModuleNotFoundError: No module named 'ete4'`. **`ete4` 4.3.0 *is* installed in `gradi-ortho`** —
  there is nothing to install, and an install attempt would be chasing the wrong bug. Name the
  interpreter: `<env>/bin/python <env>/bin/orthofinder`.
- **OrthoFinder treats every file in `-f <dir>` as a proteome.** A stale FASTA from an earlier
  `--species` subset silently joins the run. The stage deletes non-selected FASTAs before starting.
- **`Orthologues/` is many-to-many.** Reading one gene per cell would drop the co-orthologs, which
  are the whole reason to prefer OrthoFinder over RBH.
- **`is_ortholog_orthofinder` is meaningless within a species** (it is always `False` there, by
  construction). Within-species relatedness is `same_orthogroup` / `n_paralogs`.
- **The dense `n_orthologs_*` counts are unions of the two methods.** For a method-specific number
  use `n_orthologs_of_*` or `n_orthologs_rbh_*`; for the per-pair detail, `orthologs.tsv`.

## Run log

### 2026-09-03 — first full run

OrthoFinder 3.1.5 + DIAMOND 2.2.1, 4 species, 33,436 proteins, `--very-sensitive -e 1e-3`,
top-5 per target species, `-t 9`. **10.0 min** total (OrthoFinder ~7 min, 16 DIAMOND searches ~70 s,
of which human×human alone is 38 s). All guards passed.

Coverage exact in all four species — 33,436 / 33,436 placed. 5,064 orthogroups; 146,723 neighbour
rows (well under the 669k ceiling, because many proteins have fewer than five hits in a distant
species); 29,844 ortholog pairs.

**The `--very-sensitive` change did what it was chosen to do**: Kp proteins with any human hit
**732 → 1,430 (12.8% → 25.0%)**. RBH median identity to human came out slightly *lower* than v1
(34.1 vs 37.0 for Kp), which is the same effect seen from the other side — the extra pairs found are
the marginal ones, and they pull the median down.

### The spot check, and one result worth acting on

Kp → *E. coli*, rank 1, all RBH and all OrthoFinder orthologs: `clpP` **99.0%**, `clpX` 98.8%,
`rpsA` 98.9%, `ftsZ` 98.7%, `gyrB` 95.0% — consistent with v1's measurement that the ClpP machinery
is 76–99% identical between the two.

Then the human column:

| Kp protein | closest human protein | identity | positives | ortholog? |
|---|---|---|---|---|
| `clpP` | **CLPP** (Q16740), mitochondrial ClpP | **56.3%** | 79.2% | both methods |
| `clpX` | **CLPX** (O76031), mitochondrial ClpX | 41.1% | 54.7% | both methods |
| `gyrB` | TOP2A (P11388) | 24.7% | 43.3% | both methods |
| `rpsA` | PDCD11 (Q14690) | 25.3% | 45.1% | RBH only |
| `ftsZ` | **no hit** | — | — | **measured zero** |

**The degradation handle itself has a close human ortholog.** Human mitochondrial ClpP is a genuine
ortholog of bacterial ClpP at 56.3% identity / 79.2% positives, called by both methods — and stage
04's own second activator, ONC212, is an imipridone whose characterised human target *is* ClpP. This
is not a subtle inference; it falls out of the first spot check. It bears on the whole
BacPROTAC-style strategy and belongs in front of the collaboration, not in a table.

`ftsZ` is the counter-example and shows the requirement working: **no human hit at all**, and because
every accession was searched, that is a measured zero rather than a gap.

### The panel-size effect was found by the smoke test

Running `--species kpneumoniae ecoli` first, as a two-species smoke test, is what surfaced it:
OrthoFinder returned only **631** ortholog groups and the agreement with RBH was 14%. It looked like
a parsing bug. It was not — `Orthologues/kpneumoniae__v__ecoli.tsv` genuinely had 631 rows.
Orthogroup inference needs species to have signal to work with. Worth remembering before reading any
OrthoFinder count as an absolute: **it is a function of the panel, not only of the biology.**

### 2026-09-03 — OrthoDB added, because OrthoFinder's groups are panel-dependent  **[v11 — SUPERSEDED, see 2026-09-15]**

`scripts/orthology/orthodb.py`. The panel-size finding above is a limitation of *de novo*
inference, not of OrthoFinder: separating an ortholog from an in-paralog means reconciling gene trees
against a species tree, and a small panel starves that of signal. Every de novo method has it. If the
grouping must not move when the query changes, the groups have to be defined elsewhere.

**0.9 min** once the tables are cached; the streamed 5.1 GB filter dominates a cold run.

**Three corrections I made to my own earlier reasoning, all worth keeping:**

1. **eggNOG's zeros are measured, not unknown.** An earlier draft ruled eggNOG out for "9–11%
   unknowns". That was wrong: stage 02 *ran* emapper over every protein, so a protein with no OG is a
   searched negative. The unknown problem belongs to the *xref lookup* route — which is why UniProt's
   eggNOG xref is 0.00% on Kp while emapper reached 91.5%. Conflating search with lookup wrongly
   eliminated the one source that was already absolute *and* already computed.
2. **OrthoDB has no root level.** I recommended it partly on "one absolute space across all four
   species". It has 1,003 levels topping out at Bacteria / Archaea / Eukaryota — **none spanning
   domains** — so it cannot compare a bacterium to human at all. The decision still stood, because
   human selectivity is already answered by `neighbors.tsv` and OrthoFinder, and OrthoDB's real
   strength is bacterial breadth over 17,551 species.
3. **v1's OrthoDB artifact does not show native HS11286 membership.** It has a self-row for all 5,728
   Kp proteins, which looks like full coverage, but only **1,050 (18.3%)** carry a group — almost
   exactly Kp's 18.4% native gene-name coverage, because the file was built by pivoting through gene
   symbols. The self-rows are placeholders.

**The accession join is a lottery, and that is why every species goes by sequence.** OrthoDB's own
`uniprot_id` matched only 62.1% of *E. coli* and 54.8% of human — it builds from RefSeq assemblies
whose UniProt cross-references point at other substrains and at unreviewed entries our reviewed-only
human proteome excludes. *S. aureus* fared well (95.4%) only because OrthoDB happens to carry our
exact strain. Switching to DIAMOND against each species' own OrthoDB assembly, joined on
`protein_id` (100% populated), lifted **Ec 60.6% → 88.8%** and **human 54.2% → 66.5%** at median
identity 100.0%.

**It reproduces OrthoDB's own answers, checked two ways.** `clpP` came out `9802800at2` in all three
bacteria — the same id the public REST API returns for *E. coli* `P0A6G7`. And Kp `A0A0H3GPY7`
(*fadB*) got `5389341at2`, exactly what v1's gene-symbol route assigned it among the 18.3% it
reached.

**629 Bacteria-level orthogroups are shared by all three bacteria** — and that number does not move
if we change which species we run, which was the whole point.

> **Superseded.** On v12 the same intersection gives **384**, not because conservation changed but
> because OrthoDB carries parallel groups for one family and they split Kp/Ec from Sa. The
> panel-independence claim still holds; the *id-intersection method* does not. See the 2026-09-15
> entry.

One thing not done: the residual `no_group` is dominated by `no_sequence_match` (Kp 1,393,
human 6,729), i.e. the OrthoDB assembly has no counterpart, not OrthoDB having no group. For Kp that
is structural — the donor strain carries 4,975 genes against our 5,728. Searching all of OrthoDB
(the 22 GB `all_og_fasta`) would recover some of it; whether that is worth 22 GB is untested.

### 2026-09-15 — rebuilt on OrthoDB v12.2, replacing v11

v11 was replaced rather than migrated: OrthoDB's README says an "OG unique id (**not stable and
re-used between releases**)", so no id can be carried across. All four proteomes were reassigned.

**Coverage: Kp 74.6 → 91.1% · Ec 88.8 → 95.1% · Sa 88.3 → 91.3% · human → 95.5%.**

The first attempt was **worse than v11 on the anchor** — 71.6% against 74.6% — which the "must beat
the baseline" guard caught. Diagnosing it produced the three changes that matter. Kp's losses split
cleanly: 92.8% of its proteins got *some* DIAMOND hit, 80.5% survived the floors, 71.6% were
assigned. So one third of the loss was the floor and one third was the domain-level restriction,
and only the last third was genuinely unmatched.

Fixing the level restriction took Kp 71.6 → 81.0%, and adding the lineage rule plus the recalibrated
floor took it to 91.1%. The lineage rule matters on its own terms, not just for coverage: without it,
step 2 was handing Kp proteins groups defined at *Bacillus* or *Streptococcus* level, which say
nothing about Kp.

**The decoy control is what killed the 40% floor.** It had been inherited from CLAUDE.md's
annotation-transfer rule without anyone asking what it bought. 33,436 shuffled sequences produced 0
hits at 40%/50% and 2 at 25%/50% — so the floor was buying nothing and costing ~22 points on
*S. aureus*. This is the stage-02 decoy pattern reused, and it should be reused again anywhere a
threshold is inherited rather than measured.

**Three things I got wrong, worth recording.**

1. **I reported 64.7% accuracy before checking what a disagreement meant.** OrthoDB maintains
   *parallel* Bacteria-level groups for the same family — `5287828at2` and `9807019at2` are both
   real, both Bacteria-level, both *dnaA* — and **57.1% of the apparent errors are exactly that**.
   True agreement was **84.9%** (at `--reps 5`). The lesson is that a disagreement rate against a reference is
   meaningless until you look at the disagreements.
2. **My first evaluation was biased by my own caching.** Groups had only been fetched for
   *best-by-identity* hit genes, so comparing alternative rules scored them on incomplete data — it
   showed "best by bitscore" resolving *fewer* proteins than "best by identity", which is impossible.
   Fixed by fetching groups for every hit gene.
3. **A `\b`-anchored regex silently dropped `NZ_`-prefixed accessions** in the sibling essentiality
   work the same day, because `_` is a word character. Different script, same class of bug: a
   scanning regex that fails quietly on a subset of well-formed input.

**`--reps 20` was then tested and became the default** — see the entry below. And the shared-by-all-three-bacteria
count falls to **384** from v11's 629 — not because conservation changed, but because parallel groups
split Kp/Ec from Sa. **Group-id intersection across species is not a safe way to find core genes in
v12**; use names, or `neighbors.tsv`.

### 2026-09-15 — `--reps 20` measured, and made the default

The open question from the entry above was whether indexing more members per group would help. It
does, and the gain is in **trust rather than reach**:

| | `--reps 5` | `--reps 20` |
|---|---|---|
| representative sequences | 5,477,084 | **11,858,928** |
| Kp coverage | 91.1% | 92.7% |
| exact-id accuracy (E. coli, n≈3,100) | 64.7% | **75.4%** |
| name-corrected accuracy | 84.9% | **89.8%** |
| high-confidence calls (conf>0.9) | 664 | **1,178** |
| accuracy *of those* | 95.5% | **98.0%** |
| DIAMOND search | 13 min | 31 min |

Coverage moved 1.6 points; accuracy moved 10.7. That is the expected shape — with 5 members per
group the true group often has no close representative, so the vote is decided among neighbours. The
high-confidence count on Kp went 1,198 → **1,923** and on human 7,951 → **10,763**.

**Checked against OrthoDB's own API, level-matched** (`/v12/blast`, RapSearch2 — a different
algorithm on their server, searching *all* 147M genes rather than our subset). On a cached 49-protein
sample: *E. coli* **100%**, *S. aureus* **100%**, human 83%, **Kp 71%**; overall **86%** against 84%
at `--reps 5`. The comparison has to be **level-matched** — a protein sits in a nested set of groups,
one per level, so comparing whole sets penalises resolving a different level and reads far worse than
it is (72% instead of 86% on the same data).

**We are deliberately not reproducing the API exactly.** It searches all of OrthoDB with RapSearch2
and takes its best hit's LCA cluster; we search ~8% of it with DIAMOND, vote over 25 hits, and
restrict to our own lineage. The search space is the remaining gap and the reason Kp sits at 71% —
`--reps 50` would narrow it further at roughly linear cost in search time.
