# Stage 07 — essentiality

`scripts/essentiality/labels.py` · `scripts/essentiality/deg_proteomes.py` ·
`scripts/embeddings/proteomelm.py` · `scripts/essentiality/geptop.py` ·
loader `src/essentiality.py`

One question per protein: **can the organism live without it?** The anchor makes this hard —
*K. pneumoniae* HS11286 has **no essentiality measurement of its own**, so for the species that
matters most every number is a prediction.

| part | what it gives | state |
|---|---|---|
| label corpus (DEG) | 173,048 labeled proteins, 38 species, 20,194 essential | built |
| ProteomeLM embeddings | proteome-contextualised vectors, layer 8 of ProteomeLM-L | built |
| **Geptop 2.0** | **orthology + phylogeny score for all three bacteria** | **built** |
| ProteomeLM-Ess head | the paper's own essentiality head | deferred (weights unreleased) |

---

## Geptop 2.0 — orthology and phylogeny

### What it computes

    score_raw(p) = SUM over 37 reference prokaryotes r of  1/d(query, r)
                   for each r where p has a RECIPROCAL BEST HIT to a gene DEG lists essential
    score(p)     = (score_raw(p) - min) / (max - min)     <-- min-max, WITHIN this proteome

`d` is the **k=6 composition-vector (CV-tree) distance**. v1 substituted median RBH % identity for
it; that was a real deviation and is not repeated.

### It is a port, and here is exactly what changed

`data/source/geptop/geptop2.py` is kept verbatim as the reference but **cannot be
run** — Python 2 only, plus three genuine defects. The arithmetic is transcribed unchanged; only
these were fixed:

1. `os.popen(makeblastdb)` followed immediately by `os.popen(blastp)` with **no wait** — a real
   race. Now `subprocess.run(..., check=True)`, sequenced.
2. `except (Exception, err)` / `except (Exception, IOError)` — handlers that swallow failures and
   print a class as a value.
3. `import pp` (Parallel Python, py2-only, dead) genuinely used at its line 152 → `concurrent.futures`.

Also: each BLAST database is built **once**, where the original rebuilt the same 37 references 74
times per proteome.

**Two quirks were deliberately PRESERVED** because "fixing" them would change every score: a 6-mer
reachable from an observed 5-mer but never itself observed gets exactly `−1` rather than a computed
value; and `Distance()` accumulates `Q` with the same product as `O` in its shared branch, so `Q`
is not the usual sum of squares.

### The phylogeny term validates against known taxonomy

Unprompted, from transcribed arithmetic — the strongest evidence the port is faithful:

| CV distance | pair |
|---|---|
| 0.064 | *S. aureus* N315 vs NCTC 8325 (same species, different strain) |
| 0.092 | *S.* Typhi vs Typhimurium (same species) |
| 0.300 | *E. coli* vs *Salmonella* (same family) |
| 0.486 | *E. coli* vs *P. aeruginosa* (same phylum) |
| 0.496–0.499 | *E. coli* vs *B. subtilis* / *S. aureus* / *M. tuberculosis* / *M. genitalium* |

Distances **saturate near 0.5** cross-phylum, so `1/d` only spans ~2.0–3.3 across unrelated
references. That matters for the next point.

### Accuracy: 0.59–0.84, and the paper quotes the mean

Validated on DEG species that are **not Geptop references, excluded at GENUS level** — stricter
than the paper's own leave-one-out, which evaluates organisms that *are* in the reference set.

| species | n | base rate | AUROC | scored-subset AUROC |
|---|---|---|---|---|
| *R. solanacearum* (DEG1057) | 5,121 | 8.7% | **0.810** | **0.841** |
| *H. influenzae* (DEG1005) | 1,721 | 36.9% | 0.587 | 0.611 |
| paper's claim | — | — | 0.84 (mean) | — |

**Read this as 0.59–0.84, not 0.84.** Ralstonia reproduces the published figure; *H. influenzae*
does not, and its label set is the likely reason — 36.9% essential is an extreme outlier (most
bacteria are 5–15%), and it has MORE essential RBHs (66.2%) yet *worse* AUROC, which is the
signature of a screen calling too much essential rather than of a failing predictor. **The Kp
prediction inherits that variance: quoting 0.84 for Kp would be quoting a best case.**

### Spot checks (K. pneumoniae)

| gene | score | percentile | essential RBHs |
|---|---|---|---|
| `ftsZ` | **1.000** | 100.0 | 33 |
| `gyrB` | 0.971 | 99.9 | 32 |
| `rpoB` | 0.943 | 99.9 | 31 |
| `dnaA` | 0.845 | 99.0 | 28 |
| `secA` | 0.729 | 98.1 | 24 |
| `clpP` | 0.231 | 92.9 | 8 |
| `lacZ`, `araB`, `fadB` | 0.000 | — | 0 |

`ftsZ` is the single highest-scoring protein of 5,728. `clpP` sits high but **below** the 0.24
cutoff — correct, ClpP is dispensable in most bacteria, and it is the degradation handle the whole
collaboration depends on.

### Coverage — a zero means two different things

**This was got wrong once and the wrong number was reported.** `geptop_score == 0` splits:

| `geptop_evidence` | Kp | meaning |
|---|---|---|
| `essential_orthologs` | 1,929 (33.7%) | essential RBH found → score > 0 |
| `orthologs_none_essential` | 3,350 (58.5%) | has orthologs (median 7 refs), none essential — **a confident NON-essential call, evidence not absence** |
| `no_orthologs` | 449 (7.8%) | the method cannot see this protein — the only true gap |

So coverage is **92.2%**, not 33.7%. An earlier version defined `geptop_informative` as
`n_essential_rbh > 0`, conflating the middle group with the last and understating coverage by 58
points. It is now `n_rbh > 0`. **The zero-scoring rows are TIED**, so a zero means *unranked*, never
"low rank" — do not read a percentile off the zero block.

### Two traps that will bite

1. **`geptop_score` is PROTEOME-RELATIVE.** Min-max runs within each species, so 0.6 in Kp and 0.6
   in Sa are not the same quantity. Use `geptop_score_raw` for anything cross-species — the same
   trap as stage 01's per-species t-SNE coordinates.
2. **Two of our anchors ARE Geptop references, and the effect is measured, not estimated.**
   `DataSet4.faa` is *E. coli* K-12 MG1655 and `DataSet13.faa` is *S. aureus* NCTC 8325 — our exact
   stage-00 proteomes, with 296 and 350 of their own genes in DEG2. Measured share of total
   reference weight coming from the organism ITSELF:

   | species | nearest reference | d | weight | share of all weight |
   |---|---|---|---|---|
   | kpneumoniae | DataSet4 (E. coli) | 0.3497 | 2.86 | **3.6%** — clean |
   | ecoli | DataSet4 (**itself**) | 0.0113 | 88.71 | **53.2%** |
   | saureus | DataSet13 (**itself**) | 0.0083 | 120.98 | **58.3%** |

   So over half of each of those scores is self-derived. Note `d` is **not** exactly zero — our
   UniProt proteomes differ slightly from Geptop's RefSeq builds of the same strains — so the
   authors' `if d == 0: weight = 100` branch never fires, yet `1/d` reaches 89 and 121 anyway,
   *worse* than the hard-coded constant in Sa's case. **Leave-one-out is deliberately NOT applied** — for those
   two DEG supplies a measured label anyway, and Kp, the actual anchor, is absent from all 37
   references (its weights run 2.00–2.86, no self-weight). Every row carries
   `geptop_in_reference_set`; accuracy is measured on out-of-set species instead.

### Practicalities

`blastp`/`makeblastdb` are borrowed from the **`gradi-prokka`** env (`GRADI_BLAST_BIN` overrides) —
the `GRADI_RPSBLAST_BIN` pattern. Do **not** install blast into `gradi`: no osx-arm64 build, it
would drag the env to osx-64 and take ESM-C down with it.

**Human is excluded by construction**, not by choice: the 37 references are all prokaryotes with no
eukaryote among them, so there is no honest score to give it — the same reason stage 02 does not
offer `--species human` for COG.

Measured: ~22 s per blastp run, 74 runs per proteome → ~30 min per species. Composition vectors are
~3 min for all 37 references, cached. Both BLAST output and CVs are cached, so a re-run is a
re-parse.

**Load through `src/essentiality.py`** — `load_geptop`, `load_geptop_all`,
`geptop_reference_audit` (per-reference distance, weight and RBH counts — the evidence behind each
score), `geptop_informative_only`, `geptop_scored_only`, `load_labels`.

CLI: `--species` · `--validate DEG_SPECIES ...` · `--threads` · `--cutoff` · `--cv-jobs` ·
`--refresh` · `--dry-run` · `-q`.

---

## The layout — four tiers, restructured 2026-09-21

The axis used to keep its clean training sets among forty-odd audit tables and had no
machine-readable record of which datasets were actually used. Four tiers now, each answering one
question.

```
data/{raw,source}/**/SOURCE.md        34 dirs — what each dataset IS, and whether v2 uses it
data/processed/essentiality/
  dataset_registry.tsv                every dataset found + its DISPOSITION          (47 rows)
  training_sets/                      the clean ML-ready sets, one schema            (10 files)
  geptop_<sp>.tsv                     evidence source: orthology prediction
  deg_<sp>.tsv                        evidence source: DEG measurement
  ogee_<sp>.tsv                       evidence source: OGEE-trained prediction
  screens_<sp>.tsv                    evidence source: 9 published-screen predictions
  essentiality_<sp>.tsv               the headline — one column per source + the merge
  evidence/                           audits, controls, CV tables, manifests
  scratch/                            caches. PURGEABLE, so nothing cited may live here
```

### `dataset_registry.tsv` — which datasets became training sets, and why not the rest

| disposition | n | meaning |
|---|---|---|
| `training_set` | 10 | the 9 published screens + the OGEE corpus |
| `held_back` | 30 | downloaded and parseable, deliberately unused — positives-only, condition-dependent, CRISPRi, duplicates |
| `refuted` | 6 | measured and killed; the `reason` column carries the evidence |
| `corpus_only` | 1 | DEG — 66 datasets feeding one column, not one column each |

**Generated, never hand-written.** `registry.py` assembles it from `screens.py`'s own
`SCREENS`/`HELD_BACK`/`UNJOINABLE`, `screen_join_audit.tsv`, `deg_datasets.tsv` and
`ogee_taxa.tsv`. A second hand-maintained copy is how a register stops agreeing with the code it
describes. It **reconciles both directions** and exits non-zero on a mismatch: every declared
training set must exist as a file in `training_sets/`, and every file there must have a row.

**Source directories are declared (`Screen.source_dir`), never inferred.** The first version
guessed by substring and matched `essential_ecoli_bw25113_tradis_goodall` to `data/source/go` —
"go" is inside "goodall". A false provenance link reads exactly like a correct one.

### `training_sets/` — one schema, and a key that is not an accession

`key · label · source_id · features_from`, plus any screen-specific extras.

**`key` is NOT a UniProt accession for five of the nine screens.** It is whatever identifier the
screened strain's own proteome uses:

| screen | key looks like | keys onto its anchor proteome |
|---|---|---|
| the 4 b-number E. coli screens | `P69924` | 100% |
| `kpneumoniae_ecl8_tradis` | `CCN31837.1` | **0 of 4,930** |
| `kpneumoniae_atcc43816_tradis` | `WP_038431262.1` | **0 of 4,809** |
| `kpneumoniae_rh201207_tradis` | `KPNRH_00001` | **0 of 4,981** |
| `ecoli_st131_tradis` | `lcl|HG941718.1_prot_...` | **0 of 4,981** |
| `ecoli_o157h7_tnseq` | `lcl|NZ_CP008957.1_prot_...` | **0 of 5,433** |

That is deliberate, not a defect: label, sequence and embedding then share one namespace and
nothing is joined across annotations. But a join to `proteome_<species>.tsv` on accession returns
an **empty frame, not an error**, which is why the column was renamed off `uniprot_ac` and
`features_from` now names the proteome the key indexes.

### `screens_<sp>.tsv` — predictions in every column and every row

Each of the nine screens was measured on a strain that is not the anchor, so **the model is
transferred, not the label**. Transferring labels was considered and rejected: exact-sequence
transfer onto the anchor recovers only **13.7–71.6%** of rows and **29.8–61.4%** of positives, so
it would silently mislabel roughly 1,700 measured essential genes as non-essential — the failure
mode `deg_proteomes.py` already documents for DEG.

Where an anchor protein *was* in a screen's training set — the four b-number E. coli screens
overlap the E. coli anchor 100% — the **out-of-fold** value is substituted. Without that the column
would be in-sample where it overlaps and honest where it does not: two different quantities under
one name.

**Read `evidence/screens_transfer_audit.tsv` before treating a column as evidence.** It gives, per
column, the training strain, the measured own-organism grouped AUROC/AUPR, and the train/target
overlap. For scale, measured cross-species transfer is **Ec→Kp 0.8597** and **Kp→Ec 0.9338** mean
AUROC, and it is better when the *assay* matches: Goodall TraDIS reaches Kp at 0.89 while Keio
knockout reaches 0.82–0.84.

### One trap the restructure created, and the rule it produced

Purging `scratch/geptop/` to reclaim 5.7 GB also deleted the three `reference_audit.tsv` files —
which `src.essentiality.geptop_reference_audit()` reads and which are the only record of which of
the 37 references contributed. The deliverables survived; the evidence did not, and it cost a
90-minute re-run.

The files were **misfiled**. `scratch/` is documented as safe to delete, so anything cited cannot
live there, and the directory contract's own test settles it: *would you cite or check it* →
`evidence/`. They are now `evidence/geptop_reference_audit_<species>.tsv`.

## The headline — `essentiality_<species>.tsv`

### The merged column, and the units it mixes


Two headline columns per protein, written by `scripts/essentiality/merge.py`:

| column | domain | meaning |
|---|---|---|
| **`geptop_ess`** | 0–1 continuous, **never null** | Geptop 2.0 orthology+phylogeny score |
| **`deg_ess`** | 0 / 0.5 / 1 / **null** | fraction of DEG screens on this exact strain calling it essential; **null = unmeasured** |

| species | `geptop_ess` | `deg_ess` |
|---|---|---|
| kpneumoniae | 5,728 / 5,728 | **0 / 5,728 — no measurement exists** |
| ecoli | 4,403 / 4,403 | 4,253 (96.6%) |
| saureus | 2,889 / 2,889 | 2,678 (92.7%) |

Supporting columns keep the evidence: `deg_n_datasets`, `deg_n_essential`, `deg_essential_any`,
`deg_essential_all`, `geptop_evidence`, `geptop_in_reference_set`. `essentiality` /
`essentiality_source` merge the two as a convenience.

### Where `deg_ess` comes from

DEG covers **our exact anchor strains** for two species and not the third — which is the whole
reason Geptop exists:

| species | datasets | essential |
|---|---|---|
| ecoli | DEG1018 (genetic footprinting) + DEG1019 (single-gene knockout) | 614 / 296 |
| saureus | DEG1017 (TMDH) + DEG1061 (Tn-seq) | 345 / 284 |
| kpneumoniae | **none — *Klebsiella* is absent from DEG entirely** | — |

Joined **by sequence, not accession** (the house rule): 99.4% of E. coli and 96.9% of *S. aureus*
DEG rows match a proteome sequence exactly.

### Two screens of the same strain disagree, and `deg_ess` records that

| `deg_ess` | meaning | E. coli | *S. aureus* |
|---|---|---|---|
| 0.0 | neither screen | 3,558 | 2,315 |
| **0.5** | **one screen only — they disagree** | **490** | 116 |
| 1.0 | both screens agree | 205 | 247 |

**On E. coli, 490 of 695 essential calls (70%) rest on a single screen**; only 205 are agreed by
both. *S. aureus* is much better (247 of 363 agreed). The two E. coli screens are genetic
footprinting (2003) and Keio single-gene knockout (2006) and they genuinely disagree — so a merged
column that flattens 0.5 to 1 asserts far more than the data supports. **`--rule any` gives 695
essential, `--rule all` gives 205: a 3.4× spread from one choice.** Neither is "the truth".

**`deg_ess` being three-valued is INCIDENTAL** — each species happens to have exactly two screens
today. A third would make it quarters. Read it as a fraction, not a fixed scale.

### Are `geptop_ess` and `deg_ess` comparable?

**For ranking, yes. As values, no.** Measured on the out-of-set validation species:

| `geptop_ess` | n | observed % essential |
|---|---|---|
| exactly 0 | 3,343 | **3.1%** |
| 0.10–0.20 | 163 | 20.2% |
| 0.30–0.50 | 119 | 49.6% |
| 0.50–0.70 | 103 | 59.2% |
| 0.70–1.00 | 128 | **69.5%** |
| *measured essential* | — | **100%** by definition |

Monotonic, and nearly a probability in mid-range — but a predicted **1.0 means ~70%**, not 100%,
and a predicted **0.0 means 3.1%**, not 0%. A measured 1.0 means it *is* essential. So the merged
`essentiality` column is deliberately a **hybrid**: pinning measurements to 1.0/0.0 asserts more
confidence than any prediction carries. Correct in principle (measurement beats prediction), but it
means every measured essential outranks every prediction.

### Three rules for using the table

1. **Rank within a species, never across.** Ec and Sa are ~95% measured, Kp is 100% predicted; a
   pooled sort favours the measured species for a reason that is not biology.
2. **For E. coli and *S. aureus*, `geptop_ess` is circular** (53.2% / 58.3% self-derived). Use
   `deg_ess` there. Geptop's job is **K. pneumoniae**, which has no alternative.
3. **A `geptop_ess` of 0 is a tie, not a rank** — 58.5% of Kp zeros are confident negatives and
   7.8% are "no information"; `geptop_evidence` separates them.

**Load through `src/essentiality.py`** — `load(species)`, `load_all()`.
CLI: `--species` · `--rule {any,all}` · `-q`.


### A bug worth recording: the composition-vector cache collided

`_cv_cached()` originally keyed on `fasta.stem`. Every query proteome is written to
`<dir>/query.faa`, so **all species and all validation runs collided on one `cv_query.pkl`**.
Whichever ran first wrote it; everything afterwards silently read that organism's composition vector
as its own and used the wrong phylogeny weights throughout. No error, no warning.

**How it was caught:** three different species reported the same nearest reference at the same
distance, and *S. aureus* was reported at 0.354 from *Salmonella* — which is *E. coli*'s distance to
*Salmonella*, not a cross-phylum one.

**Why it did little damage in the end**, which is also why it was easy to miss: CV distances
saturate near 0.5 cross-phylum, so for out-of-set queries the weights are nearly uniform (2.00–2.86)
whichever vector is used — the validation numbers moved by <0.005. It was severe only for the
in-set species, where the correct weight is 89–121 instead of 2.86. Kp was unaffected entirely
(it ran first; re-run agreed to max |diff| 0.00e+00).

The key is now a hash of the resolved path, and the pre-check in `main()` shares one
`_cv_cache_path()` function with `_cv_cached()` — they had already drifted once, with the pre-check
reporting "37 cached" while the real key recomputed all of them.
