# Stage 07 — essentiality

`scripts/essentiality/labels.py` · `deg_proteomes.py` · `geptop.py` · `ogee.py` ·
`ogee_proteomes.py` · `ogee_dataset.py` · `merge.py` · `registry.py` ·
`scripts/embeddings/proteomelm.py` · loader `src/essentiality.py`

The published screens — the ten training sets, what was held back and what was refuted — are in
**`docs/essentiality_screens.md`**.

One question per protein: **can the organism live without it?** The anchor makes this hard —
*K. pneumoniae* HS11286 has **no essentiality measurement of its own**, so for the species that
matters most every number is a prediction.

| part | what it gives | state |
|---|---|---|
| label corpus (DEG) | 173,048 labeled proteins, 38 species, 20,194 essential | built |
| ProteomeLM embeddings | proteome-contextualised vectors, layer 8 of ProteomeLM-L | built |
| **Geptop 2.0** | **orthology + phylogeny score for all three bacteria** | **built** |
| **OGEE-trained model** | **a second learned opinion, 26 measured prokaryotic proteomes** | **built** |
| ProteomeLM-Ess head | the paper's own essentiality head | deferred (weights unreleased) |

---

## The DEG label corpus — `labels.py` + `deg_proteomes.py`

Stage 07 predicts essentiality with **ProteomeLM-Ess** (Malbranke, Zalaffi & Bitbol, PNAS 2026,
doi `10.1073/pnas.2524201123`; papers in `docs/papers/`). These two scripts build its labels —
positives from DEG, negatives from the screened organisms' own proteomes.

### `labels.py` — every dataset, every column, and the sequence

Fetches DEG's three bulk files and builds the label table with **nothing dropped**: 66 datasets,
42 species, **26,619 essential genes**, each carrying **its protein sequence** (`DEG10.aa`, 1:1 by
DEG gene id).

**The sequences are the point.** The join downstream is **BY SEQUENCE**, because DEG's identifier
columns are patchy:

| identifier | coverage |
|---|---|
| gene symbol | 63.1% |
| GI | 54.5% |
| COG | 48.4% |
| UniProt AC | 45.7% |
| locus tag | 42.9% |

**Exclusions are FLAGGED, never dropped** — the `retained` column says what training may use:

- **4 non-genome-wide methods** — antisense RNA, MATT, insertion-duplication,
  transposon-hybridisation. "Absent from the list" cannot mean non-essential when the method never
  interrogated the whole genome.
- **11 condition-specific screens** — tobramycin, murine pneumonia, bile, cholesterol, kanamycin.
  This is why *P. aeruginosa* PAO1 swings **117 → 336 → 551** essential genes across its three
  datasets.

**51 datasets retained over 38 species.**

**Consolidation to one label set per species is deliberately deferred.** *S. aureus* appears 7×,
*P. aeruginosa* 4×, *Salmonella* 4×, *E. coli* 4× — so the union-vs-intersection rule gets chosen
against measured counts rather than guessed in advance. (The same choice, made downstream for the
anchor strains, is the 3.4× `--rule any`/`--rule all` spread recorded under *the headline* below.)

### Four DEG traps

1. **326 of the 26,619 "sequences" are the literal string `Not available now.`** The join to
   `DEG10.aa` *succeeds*, so a coverage check reads 100% while the payload is junk. Only validating
   the amino-acid alphabet catches it — and *E. coli* O157:H7 is **14.1% placeholder**. A further
   **52 carry an internal `*` stop** (frameshifted translations, mostly *N. gonorrhoeae*); a
   *trailing* stop is stripped, not rejected. **Short proteins are NOT filtered** — `rpmJ` is a
   real 37-aa essential gene.
2. **`DEG1058` (*S. suis*) is an in-vivo pig infection screen** whose condition string reads
   "Columbia blood base agar". Its paper is titled *"…**conditionally** essential genes for
   S. suis infection in pigs"* and every gene carries a note `Recovered from
   blood/cerebrospinal fluid/meninges`. Caught by reading DEG's **per-gene notes**, not its
   condition field — so the filter is data-driven, not a hand-maintained blocklist.
3. **DEG's index disagrees with its own shipped rows** on 2 of 66 datasets — *S. oneidensis* 403 vs
   402, *R. palustris* 522 vs **552**. Printed every run.
4. **Two columns are undocumented**: col 10 packs `locus_tag:X;gi:N`, and col 13 is the per-gene
   note that exposed trap 2.

### `deg_proteomes.py` — the negative class, which DEG does not ship

**38 of 66 datasets report a non-essential *count* and no genes.** So each retained dataset's
complete proteome is fetched from the RefSeq/EMBL replicon accessions DEG records — NCBI
`efetch db=nuccore rettype=fasta_cds_aa`, which also yields locus tags — and then labeled:
essential = matched a DEG positive, non-essential = the rest. ProteomeLM needs the whole proteome
anyway, since it contextualises every protein against every other.

**Exact sequence matching is NOT enough, and the per-dataset rate is the control.** DEG's vintage
annotation has drifted from the current replicon: exact match recovers only **59–95%**
(median **84.8%**), and the unmatched are **same-length point-substituted**, not a start-codon or
prefix artefact — measured, and `locus_tag` is **0%** on the older datasets so there is no
identifier fallback. With **DIAMOND ≥95% identity** (borrowed from `gradi-ortho` via
`GRADI_DIAMOND_BIN`) the median rises to **99.5%**, recovering **+1,320 positives at median 100%
identity**. Without it, hundreds of genuinely essential genes are silently labeled non-essential.

**Do not scan for accessions with a `\b`-anchored regex.** DEG's replicon field mixes `,`, `;`,
`, ` and bare spaces, and a scanning regex silently drops **`NZ_`-prefixed RefSeq accessions**:
`_` is a word character, so the boundary never fires. That cost *E. coli* O157:H7 — the largest
single dataset, 1,071 positives — on the first run. **Split on delimiters, then full-match each
token.**

### What came out

**49 of 51 datasets · 173,048 labeled proteins · 20,194 essential (11.67%) · 38 species** —
comparable in scale to the paper's 213,608 labels over 83 genomes, and bacteria-only.

Three known gaps, **flagged in `proteome_join.tsv` rather than hidden**:

| dataset | organism | what happened |
|---|---|---|
| `DEG1003` | *V. cholerae* N16961 | joins at only **74.0%** — real annotation drift, below the 80% floor |
| `DEG1037` | *S. pyogenes* MGAS5448 | records **no replicon at all** |
| `DEG1053` | *B. cenocepacia* K56-2 | records a WGS **master** accession; `efetch` returns zero CDS |

**Spot checks pass** on *E. coli* MG1655: `dnaA`, `gyrB`, `rpoB`, `ftsZ`, `murC`, `mraY`, `rplF`
all essential; `lacZ` and `araB` not.

**Watch the base rate before training** — it spans **15× across species**, from *M. genitalium*
**72.7%** essential (a minimal genome) to *M. avium* **4.7%**. A model fit across the corpus is fit
across that spread; see the OGEE section below, where the same quantity turns out to predict how
hard a held-out organism is.

### Practicalities

Outputs, all in `data/processed/essentiality/evidence/`:

| file | shape | what it is |
|---|---|---|
| `deg_datasets.tsv` | 66 × 22 | every dataset, with `retained` and `exclusion_reason` |
| `deg_genes.tsv` | 26,619 × 28 | every essential gene, with its sequence |
| `proteome_join.tsv` | 51 × 14 | **the join control** — exact/DIAMOND/total rate per dataset |
| `labeled_proteins.tsv` | 173,048 × 11 | the corpus itself |

Raw and provenance: `data/source/deg/SOURCE.md`, proteomes under
`data/source/ncbi/deg_proteomes/<deg_dataset_id>.faa`.

CLI (both): `--refresh` · `--dry-run` · `-q`; `deg_proteomes.py` adds `--limit N` (smoke test,
writes only `scratch/smoke_*`). ~1 min and ~15 min respectively; both cache and re-run cheaply.

---

## Geptop 2.0 — orthology and phylogeny

**Geptop 2.0 is Wen 2019**, the latest release. The upstream is `github.com/RiversDong/geptop` — a
**25 MB `.rar` of Python-2 code**, which is the whole reason this is a port rather than a
dependency.

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

### Accuracy: 0.59–0.81 over all proteins, and the paper quotes the mean

Validated on DEG species that are **not Geptop references, excluded at GENUS level** — stricter
than the paper's own leave-one-out, which evaluates organisms that *are* in the reference set.

| species | n | base rate | AUROC | scored-subset AUROC |
|---|---|---|---|---|
| *R. solanacearum* (DEG1057) | 5,121 | 8.7% | **0.810** | **0.841** |
| *H. influenzae* (DEG1005) | 1,721 | 36.9% | 0.587 | 0.611 |
| paper's claim | — | — | 0.84 (mean) | — |

**Read this as 0.59–0.81 over all proteins (0.61–0.84 over scored proteins), not 0.84.** Those are
two different metrics and `geptop_validation.tsv` keeps them in two columns — do **not** quote
"0.59–0.84" as one range, which silently takes the floor from one and the ceiling from the other.
Ralstonia reproduces the published figure; *H. influenzae*
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

## The OGEE-trained model — `ogee.py` + `ogee_proteomes.py` + `ogee_dataset.py`

A **third opinion** on essentiality, learned from **26 measured prokaryotic proteomes**. OGEE v3
supplies labels and no sequences; these three scripts characterise the corpus, attach sequences,
and train. Recovery of the dataset itself from a dead server is in
`docs/essentiality_screens.md` §5.

Output `ogee_<species>.tsv`, two columns beside the accession:

| column | domain | meaning |
|---|---|---|
| `ogee_ess` | 0–1, **never null** | the model's prediction |
| `ogee_evidence` | 1 / 0 / **empty** | the **MEASURED** OGEE label for that exact protein |

### Half of OGEE is unusable, and the reason is one we had already recorded from the other side

Of **87 taxa with decided E/NE calls, 40 have ZERO negatives** — **17,743 positives-only
entries** — and **25 of those 40 come from one paper**, PMID **29769716** (Price 2018), the
Fitness Browser RB-TnSeq collection. **RB-TnSeq cannot see essential genes by construction**: no
insertions survive, so an essential gene is *absent* from the table rather than carrying an
extreme value, which makes its essential list positives-only by the nature of the assay. That is
the same finding `docs/essentiality_screens.md` §5 records about `feba.db`, arriving here from the
label side instead of the data side.

It inflates OGEE's overall base rate to **0.310**, against **0.117** for the DEG corpus. Applying
the project's own **both-classes + prokaryote** filters leaves **26 taxa / 78,893 proteins /
base 0.173**.

### Leave-species-out, not the paper's single holdout

The project owner's instruction, and the better fit: **K. pneumoniae is absent from OGEE
entirely**, so the question this axis actually asks is *"given N measured bacteria, how well can
we call the N+1th?"* — which is exactly what a held-out taxon measures. It also gives **26
estimates instead of 1**, and the spread is the result:

**AUROC 0.529–0.940, mean 0.782, median 0.805** (`evidence/ogee_leave_species_out.tsv`).

### The headline finding: `corr(base_rate, AUROC) = −0.643`

**Screens that call many genes essential are much harder to predict.** Two same-organism pairs
make it unarguable — same species, same genome, different screen:

| taxon | base rate | AUROC |
|---|---|---|
| *P. aeruginosa* PAO1 | 0.081 | **0.837** |
| *P. aeruginosa* UCBPP-PA14 | 0.302 | **0.634** |
| *Salmonella* Typhi CT18 | 0.105 | **0.782** |
| *Salmonella* Typhimurium SL1344 | 0.404 | **0.529** — barely above chance |

A screen calling 40% of genes essential is measuring **fitness defect, not essentiality**.

**This is what prices the Kp column.** Our three measured Kp screens run base **0.075–0.106**,
squarely in the band where held-out taxa score **0.78–0.93**.

### One model for all three anchors, and the caveat carried as data

**ONE model, fit on all 26 taxa, used for all three anchors** — the project owner's instruction.
The caveat that creates is carried **per protein**, in `ogee_evidence`: *E. coli* K-12 and
*S. aureus* NCTC 8325 **are OGEE taxa AND our anchors**. Measured, **100.0%** and **97.5%** of
their corpus proteins are literally the **same SEQUENCE** as an anchor protein. So wherever
`ogee_evidence` is non-null, the model was fitted on that protein with that label, and `ogee_ess`
is closer to **recall than prediction**.

| species | proteins with a measured OGEE label |
|---|---|
| ecoli | 4,193 / 4,403 (95.2%) |
| saureus | 2,815 / 2,889 (97.4%) |
| **kpneumoniae** | **0 / 5,728** |

It is visible in the output: same model, top-decile cut **0.861 on E. coli against 0.471 on Kp**.
**So `ogee_ess` is comparable WITHIN a species, never across** — the same rule `geptop_score`
carries for the same kind of reason.

### Three traps

1. **`locus` is 100% populated but heterogeneous**, and the namespace is a property of the
   **(taxon, dataset) block** — pure at median 1.000, and 87 of 89 taxa use one throughout.
   **E. coli K-12 is one of the two exceptions and carries TWO DISJOINT namespaces** (PEC gene
   symbols, and b-numbers, **zero shared ids**), so a naive `(taxid, locus)` dedup counts every
   gene twice — **9,496 "genes" for a 4,403-gene organism**. **Dedup only AFTER resolving to a
   protein.**
2. **A missing underscore cost three taxa entirely.** OGEE writes `HI0001`/`HP0001`/`MPN001`
   where the assembly writes `HI_0001`/`HP_0001`/`MPN_001`; they scored **0.000** until a
   punctuation-insensitive fallback was added (*H. influenzae* → 0.960, *H. pylori* → 0.994).
3. **GCA and GCF are not interchangeable here.** *Synechococcus* keys on RefSeq
   `SYNPCC7942_RS*` tags that only the GCF annotation carries: **0.226 on GenBank, 0.969 on
   RefSeq**. So the assembly is chosen **by measured join rate across candidates, not by rule** —
   the inverse of the GenBank-over-RefSeq preference the screens need, which is why neither is a
   default. Per-candidate rates: `evidence/ogee_proteome_join.tsv`.

### `--seeds 1` is deliberate here

Against the axis-wide 5. **The folds are FIXED** — each fold is one taxon — so there is no
partition randomness to average, only the forest's own, and per-seed SDs on this axis run
**≤0.004**. It is also **26 fits instead of 130**: measured ~4 min per fit on 76,000 × 1,152, i.e.
**1.7 h against 8–9**. The SD prints `(1 seed)` and stores `NA`, never a misleading `+/-0.0000`.

### Practicalities

**Load through `src/essentiality.py`** — `load_ogee`.

CLI: `ogee.py --rule {any,all,majority}` · `ogee_proteomes.py --candidates N` ·
`ogee_dataset.py --write-fasta | --embed-plan | --score-anchor SPECIES`.
Provenance: `data/source/ogee/SOURCE.md`. Evidence: `evidence/ogee_leave_species_out.tsv`,
`ogee_proteome_join.tsv`, `ogee_taxa.tsv`.

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

**`merge.py` also writes `deg_<species>.tsv`**, carrying `deg_ess` and its four provenance
columns. Those used to be inlined into the headline and were **split out so every evidence source
has the same shape** — its own file, plus one summary column in `essentiality_<species>.tsv`.
That is the relationship `geptop_<species>.tsv` always had with the headline, now applied
uniformly to DEG, OGEE and the screens.

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

### Two tables at the task root, deliberately

This axis ships **both `essentiality_<species>.tsv` and `geptop_<species>.tsv`** at the task root.
The second is the raw Geptop output — `geptop_score`, `geptop_score_raw`, `geptop_essential`,
`geptop_n_rbh`, `geptop_n_essential_rbh`, plus the evidence / informative / reference flags. It
*reads* like evidence, and the directory contract's test would file it there, but **the project
owner chose it as a deliverable. Do not demote it to `evidence/` in a tidy-up.**

**Load through `src/essentiality.py`** — `load(species)`, `load_all()`, plus the per-source and
registry loaders: `list_training_sets`, `load_training_set`, `registry`, `load_deg`, `load_ogee`,
`load_screens`, and the Geptop set (`load_geptop`, `load_geptop_all`, `geptop_reference_audit`,
`geptop_informative_only`, `geptop_scored_only`, `load_labels`).
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
