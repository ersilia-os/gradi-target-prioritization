# Stage 02 — functional annotation

`scripts/function/cog.py` · `scripts/function/eggnog.py` · `scripts/function/goslim.py` ·
`scripts/plots/function.py` ·
loader `src/function.py`

Broad functional classes per protein, assigned **from sequence**. Two schemes, deliberately
answering different questions with different coverage:

| scheme | question | tool | vocabulary |
|---|---|---|---|
| **COG** (`cog_<species>.tsv`) | which orthologous group | COGclassifier / RPS-BLAST vs NCBI CDD | 26 COG2024 letters, 4 groups |
| **GO slim** (`goslim_<species>.tsv`) | what does it do | UniProt curation + **eggNOG-mapper v2** | 97 `goslim_prokaryote` terms, 3 aspects |

Both come from established tools, and **neither reaches 100%**. That is deliberate. Roughly a
quarter of each proteome is short, uncharacterised and carries no Pfam or InterPro either — nothing
annotates it, and an empty row is the honest answer. Coverage figures are in each part's run log.

`eggnog_<species>.tsv` sits beside them: everything else the same eggNOG-mapper run produced —
KEGG KO and pathways, EC, BRITE, CAZy, PFAMs, preferred names, and an independent second opinion on
the COG letter.

A third scheme would get `<scheme>_<species>.tsv` beside them in `evidence/`, its audits and
vocabularies in `evidence/` too, and its raw tool dumps in `scratch/`.

---

# Part 1 — COG functional categories

## Why COG, and why computed rather than looked up

The 26 COG functional categories, in four super-groups, are the standard broad functional vocabulary
in bacterial genomics: one letter per protein. That is the point. InterPro gives thousands of
domains and PANTHER thousands of families — informative, but not *broad*.

**The lookup route fails on the anchor organism.** Measured on the stage-00 annotation layer:

| eggNOG xref coverage | Kp HS11286 | Ec K-12 | Sa NCTC 8325 | human |
|---|---|---|---|---|
| `data/processed/proteomes/evidence/annotation_<species>.tsv` | **0.00% (n = 0)** | 92.98% | 79.23% | 91.48% |

Zero on *K. pneumoniae*. HS11286 is a dark TrEMBL proteome, so every accession-based route to COG
dies on exactly the species that matters most. A sequence-based one does not — the same finding that
forces DIAMOND-by-sequence everywhere else in this project.

**What v1 had here, and why it is not inherited.** Checked, not assumed:

| v1 artifact | what it was | why not reused |
|---|---|---|
| `02a_interpro_annotation.py` | InterPro/Pfam via UniProt xref, Kp 86.8% | thousands of domains, not a broad class |
| `02b_panther_annotation.py` | PANTHER 19.0 via UniProt xref, Kp 71.2% | worse coverage than InterPro, and PANTHER *is* an InterPro member DB; 88 MB file downloaded twice |
| `functional_class` in `08a_webapp_export.py:93-128` | hand-written substring keyword map, 11 ordered first-match-wins rules | `other` = 24.6% Kp / 17.0% Ec; rule order load-bearing and undocumented; materialised only inside the webapp JSON |
| eggNOG-mapper v2 | proposed at `legacy/docs/01_task_agnostic.md:240` | **never built** (`HISTORY.md:64`) |

v1 never computed COG at all. This stage is that gap closed, with a much lighter tool than the one
originally proposed.

## The tool

**COGclassifier 2.0.0** (MIT, `pip install cogclassifier`). Protein FASTA in, category letter out:

```
RPS-BLAST vs the CDD COG profile DB  ->  best hit  ->  CDD ID  ->  COG ID  ->  letter  ->  group
```

on **COG2024** — 5,050 COGs, 26 categories, 4 groups. `cog_definition.tsv` (cog-24.def.tab) and
`cog_func_category.tsv` (cog-24.fun.tab) ship inside the package; only the profile DB is downloaded.

| | |
|---|---|
| downloads | `Cog_LE.tar.gz` 200 MB + `cddid.tbl.gz` 7.8 MB, cached in `data/source/cdd/` |
| eggNOG-mapper would have been | ~11.4 GB compressed and **50.6 GB unpacked** (measured in Part 2; the ~20 GB estimated here beforehand was wrong), for the same letter |
| runtime | **9 min wall clock** for all 13,020 proteins, 9 threads, Rosetta |

We use the **library API only**, never its CLI and never its Altair charts — this repo plots with
stylia. And we drive `RpsBlast` directly rather than `CogClassifier.run()`, because the latter writes
its hit table to a tempdir; keeping it is what makes the run resumable.

### rpsblast is borrowed from another env

`rpsblast` is not in `gradi`. It comes from **`gradi-prokka`** (`rpsblast 2.17.0+`, Mach-O x86_64,
runs under Rosetta) by prepending that env's `bin` to `PATH` — the same pattern v1 used for
`GRADI_DIAMOND_BIN`. Override with `GRADI_RPSBLAST_BIN=<dir>`.

**Do not conda-install blast into `gradi`.** There is no osx-arm64 build, so it would drag the whole
environment to osx-64 and take ESM-C down with it.

## Scope: the three bacteria, and human cannot be added

COG2024 is **2,296 genomes of bacteria and archaea, zero eukaryotes**. Human would hit only the
conserved core and read as misleadingly sparse, so it is excluded by construction, not by choice.
`--species human` is not even offered. Host/off-target comparison is orthology — a later stage.

## Outputs

```
data/processed/function/
    function_<species>.tsv                 THE DELIVERABLE — 3 columns, keyed on uniprot_ac
    evidence/
        goslim_matrix_<species>.tsv        the same content as a matrix   (n, 99)
        cog_matrix_<species>.tsv           the same content as a matrix   (n, 28)
        cog_<species>.tsv                  the per-protein table — 8 columns, keyed on uniprot_ac
        cog_counts_<species>.tsv           per-letter counts, all 26 categories, zeros included
        cog_control_ecoli.tsv              the ground-truth comparison, row by row
        cog_crosscheck.tsv                 COGclassifier vs emapper vs NCBI (Part 2)
        cog_manifest.tsv
    scratch/
        cog_rpsblast_<species>.tsv         raw outfmt-6 hits (the resumable cache)
        .cog_rpsblast_<species>.key        md5(accession set | evalue) — the cache key
data/source/cdd/                  Cog_LE/, cddid.tbl.gz, cog_func_category.tsv, SOURCE.md
output/plots/function/cog_categories.png
```

| column | notes |
|---|---|
| `uniprot_ac` | join key |
| `cog_id` | e.g. `COG0542` |
| `cog_category` | **exactly one letter** — the broad class. Use this one. |
| `cog_category_all` | the whole string. COG4862 is `KTN`. |
| `cog_group` | one of the four super-groups |
| `cog_name` | e.g. `ATP-dependent Clp protease ClpA` |
| `cog_evalue`, `cog_identity` | RPS-BLAST evidence, so a stricter cutoff can be applied downstream **without re-running** |

Every protein has a row. One with no COG keeps its row with empty strings — absence of a category is
a fact about the protein, not a missing record.

Load through `src/function.py` (`load_cog`, `load_cog_all`, `load_cog_counts`, `informative`), which
reads with `keep_default_na=False` like `src/proteomes.py`, so missing values are `""` and not NaN.

## Run log — 2026-09-01, first full run

COGclassifier 2.0.0 · COG2024 · rpsblast 2.17.0+ · e-value 1e-2 · 9 threads · 9 min wall clock.

| species | n | classified | informative | R + S | multi-letter | no COG |
|---|---|---|---|---|---|---|
| Kp HS11286 | 5,728 | 4,528 (79.1%) | 4,170 (72.8%) | 358 | 570 | 1,200 |
| Ec K-12 | 4,403 | 3,719 (84.5%) | 3,435 (78.0%) | 284 | 461 | 684 |
| Sa NCTC 8325 | 2,889 | 2,112 (73.1%) | 1,898 (65.7%) | 214 | 240 | 777 |

*informative* excludes `R` (general function prediction only) and `S` (function unknown). Those are
classified, but they are not an answer, and reporting them inside a coverage figure would flatter it.

For scale: v1's keyword `functional_class` put 75.4% of Kp outside `other`. COG reaches a comparable
79.1% while being an actual controlled vocabulary rather than an ordered list of substrings.

### The E. coli control — 97.79% agreement

*E. coli* K-12 MG1655 (`GCF_000005845.2`) is itself one of the 2,296 COG2024 reference genomes, so
NCBI publishes **curated COG assignments for our exact anchor proteome**. This is a real ground
truth, not a second prediction. Same spirit as stage 00's E. coli naming control.

| | |
|---|---|
| NCBI curated assignments for this genome | 3,594 proteins |
| comparable (both sides assigned) | **3,575** |
| COG **id** agreement | 94.91% |
| COG **category** agreement | **97.79%** (floor 90%, enforced — the script exits non-zero below it) |

The 79 letter disagreements are the expected best-hit-vs-curated difference on multi-domain
proteins: `b0238` (hypoxanthine phosphoribosyltransferase) best-hits a PRTase profile (`COG2236`,
`H`) where NCBI curates the whole protein as `COG0503`/`F`. Category agreement being *higher* than
id agreement is the signal that the chain is sound — near-miss COGs usually share a category.

The row-by-row comparison is in `evidence/cog_control_ecoli.tsv` with `cog_id_agrees` and
`cog_category_agrees` flags.

### Why coverage is not 100%, and why it should not be

Asked and answered with numbers, 2026-09-01. Four separate checks:

**Nothing is lost in the pipeline.** Proteins with an RPS-BLAST hit but no final category —
i.e. lost in the CDD → COG ID → definition chain, which COGclassifier drops with a `logger.debug` —
number **2 (Kp), 0 (Ec), 3 (Sa)**. The gap is "no hit", not "hit discarded".

**COG does not span a proteome.** On E. coli K-12, against NCBI's own curated assignments:

| | coverage |
|---|---|
| NCBI curated COG2024 | 81.6% |
| this stage | **84.5%** |
| either | 84.9% |

We are *above* the curated reference, and the union adds 0.4pp. There is nothing to recover.

**The unclassified are a coherent class, not a random sample:**

| | no COG | classified |
|---|---|---|
| median length | 86–98 aa | 288–300 aa |
| under 100 aa | 51–56% | ~6% |
| "Uncharacterized" / DUF | 42–64% | 2–10% |
| no Pfam either | 25% Ec, **61–62%** Kp/Sa | 1–3% |
| no InterPro either | 24% Ec, **59–60%** Kp/Sa | 0.3–1.5% |

Short, uncharacterised, largely accessory-genome. For Kp and Sa roughly 60% of them are invisible to
Pfam and InterPro too — nothing annotates these, not just COG.

**Loosening the e-value buys noise.** Every S. aureus sequence was shuffled (same length, same
amino-acid composition, no real structure) and re-run as a decoy set:

| e-value | real | shuffled decoys |
|---|---|---|
| 1e-10 | 68.7% | 0.0% |
| 1e-5 | 70.9% | 0.0% |
| **1e-2 (default)** | **73.2%** | **0.4%** |
| 1e-1 | 76.0% | 3.6% |
| 1 | 84.3% | 24.0% |
| 10 | 97.5% | **84.3%** |

At e-value 10 you can report ~98% coverage — and so can nonsense. **Do not raise `--evalue` to make
the coverage table look better.** The default sits at the knee.

Beware a second trap here: measuring agreement against NCBI at looser thresholds looks reassuring
(97.95% at 1e-10 falling only to 97.30% at e-value 10) because that comparison is restricted to the
~3,550 proteins NCBI *did* assign, a set that barely grows. All the extra coverage lands in the
bucket NCBI's curators declined to assign, where there is no ground truth at all.

**The one genuine headroom** is E. coli-shaped: only 25% of its unclassified proteins lack a Pfam
domain, so ~75% of them carry a family COG simply has no group for. That is a scope gap, and a
second classification in this stage could close it. For Kp and Sa the same figure is ~61%, so there
the gap is real dark matter and a second scheme would add much less.

### Spot checks, enforced every run

`clpP`, `clpX`, `clpA`/`clpC`, `clpB` → `O`; `rpsA`, `rplB` → `J`. All 18 passed. These are not
matters of opinion: if a Clp protease subunit is not `O`, the CDD → COG → letter chain is broken
whatever the coverage table says, and the script exits non-zero.

### What the figure confirms

`output/plots/function/cog_categories.png` — Kp and Ec track closely across all 26 categories
(both Enterobacteriaceae, as they must), while Sa diverges. Two free sanity signals fall out:
*S. aureus* has essentially **zero `N` (cell motility)**, which is correct — it is non-motile — and
`Y` (nuclear structure) and `Z` (cytoskeleton) are empty for all three, as they must be for bacteria.

## Traps

1. **`cog_category` is a choice, not a formality.** 9–16% of classified proteins carry a multi-letter
   COG. COG orders those letters by importance and we take the first, but `cog_category_all` keeps
   the string and the count is printed every run. Never resolve it silently.
2. **`R` and `S` are classified but uninformative.** Anything reporting "COG coverage" must say which
   of the two numbers it means. `src.function.informative()` drops them.
3. **The hit cache is keyed, not counted.** A hit table cannot tell you what it covers, because a
   query with no hit leaves no row — so a 50-protein smoke-test cache would happily serve a
   2,889-protein run. The key is `md5(sorted accessions | evalue)` in
   `.cog_rpsblast_<species>.key`; a mismatch prints `[rerun]`. This was a real bug, caught before
   the first full run.
4. **`cog_definition.tsv` is ragged** — 4 to 10 tab-separated fields per line (5,043 of 5,050 rows at
   7). `pd.read_csv` rejects it. Parse it with the package's own `CogDefinitionRecord`.
5. **The CDD profile DB is volume-split.** `Cog_LE/` unpacks to `Cog.00.*`, `Cog.01.*` behind a
   `Cog.pal` alias — there is no `Cog.pn`. `-db Cog_LE/Cog` resolves it.
6. **The control joins on locus tag, not RefSeq.** `cog-24.cog.csv` column 0 is the locus tag — for
   MG1655 the b-number, which stage 00 stores at 100% coverage. Column 2 is a RefSeq accession, but
   for this genome NCBI uses `NP_` while UniProt cross-references mostly `WP_`, so joining on it
   returns **zero** overlapping rows. "Join on `locus_tag`, not gene name" applies to RefSeq too.
7. **The query FASTA is built from the stage-00 table, not `data/source/uniprot/proteomes/*.fasta`.** Those
   carry `sp|A5A616|MGTS_ECOLI` headers, so every hit would need parsing back to an accession.
   Building it here makes `QUERY_ID` *be* the join key and guarantees the query set is exactly the
   rows in the deliverable.
8. **`cog-24.cog.csv` is 637 MB** and is never downloaded — the control streams it through `grep`
   and keeps ~3.6k lines. It takes ~2 min. `--no-control` skips it.
9. **`--limit` never touches the real deliverable.** A smoke test writes everything under
   `scratch/smoke_*` — table, hit cache, counts, manifest — and skips the control, which would be
   meaningless on a truncated proteome. Without this, `--limit 50` silently replaced a 2,889-row
   `cog_saureus.tsv` with 50 rows. Caught before the first commit, not after.

---
---

# Part 2 — GO slim (`goslim_prokaryote`)

`scripts/function/eggnog.py` (the tool) · `scripts/function/goslim.py` (the vocabulary) ·
97 terms, purpose-built for bacteria.

## Why eggNOG-mapper, and why nothing else would do

COG's gap is a property of COG, not a bug (Part 1). Crucially, the proteins it misses are short,
uncharacterised, and for Kp/Sa **~60% carry no Pfam and no InterPro either**.

So a second scheme only helps if it uses a **different kind of evidence**. That the domain route
does not was measured twice, independently:

| test | result |
|---|---|
| InterPro2GO as a GO tier | filled **5 proteins out of 13,020** |
| Kp proteins with InterPro but no GO, having an InterPro entry that carries a GO mapping | **13 of 825 (1.6%)** |
| InterPro **live API** vs UniProt's xref, on GO-less Kp proteins | **0 of 15** would gain a GO term |

The cause is that UniProt's electronic GO is *already* InterPro2GO-derived. InterProScan,
`Pfam-A + pfam2go` and the InterPro API therefore all re-derive what stage 00 already has. InterPro's
live release does carry entries UniProt's xrefs lag behind — but those new entries have no GO
mappings, so the headroom is zero.

**Orthology is the only route left**, and that means **eggNOG-mapper v2** — the field standard for
bacterial functional annotation, and the tool v1 identified at `legacy/docs/01_task_agnostic.md:240`
and never built (`HISTORY.md:64`).

Other candidates, checked not assumed:

| candidate | verdict |
|---|---|
| **eggNOG-mapper 2.1.15** | bioconda **`osx-64`** build → runs under Rosetta. **Chosen.** |
| InterProScan 5 | bioconda **`linux-64` only**. Not installable here. |
| KofamScan | 1.5 GB and a fine tool, but emits **KEGG KO, not GO** — a different vocabulary. |
| `deepgoplus` 1.0.2 | pins **Python <3.8** (TensorFlow era). Dead on arm64 / py3.11. |
| `deepgo2` (DeepGO-SE) | pins **`dgl==1.1.2+cu117`** — CUDA-only, no arm64 build. |

### The rejected ESM-C tier — recorded so it is not reinvented

An earlier version of this stage reached **100% coverage** by transferring GO from the nearest
curated neighbour in ESM-C embedding space (stage 01). It was removed on the explicit instruction to
use well-established tools only, and that was the right call twice over:

- A hand-rolled k-NN is not a citable method, has no external validation, and no users but us.
- **Its honest accuracy was well below its headline.** Leave-one-out on E. coli gave MF 86.8% — but
  that donor bank is dominated by proteins with near-identical homologs, while a real query is
  unannotated *precisely because it is unusual*. Reweighted to the donor distances that actually
  occurred, it was **MF 68% · BP 53% · CC 66%**.

**Coverage is therefore no longer 100%, deliberately.** A protein no established tool can annotate
gets an empty row instead of a guess. `function/goslim.py` prints those under **NOT ANNOTATED**
and says so.

## Two tiers

Every row carries the tier that produced it in `goslim_source`:

| tier | what | Kp | Ec | Sa |
|---|---|---|---|---|
| `curated` | UniProt's own GO, mapped to slim | 69.7% | 87.0% | 65.0% |
| `eggnog` | eggNOG-mapper v2, orthology-based | **73.7%** | **88.7%** | **65.6%** |

`curated` sits below the raw GO coverage (72.9 / 90.3 / 67.4%) because a protein can carry GO that
maps to **no slim term at all** — 186 such in Kp. That is not a fill and is not counted as one.

## One primary term per aspect

GO is multi-label across three aspects, so each gets one chosen term plus the full set — the shape
COG has. The choice is an explicit total order, so a tie can never leave it empty:

```
deepest in the GO DAG (most specific)  ->  rarest in this proteome (most informative)  ->  lowest GO id
```

| column | notes |
|---|---|
| `goslim_mf` / `_bp` / `_cc` | the chosen GO id per aspect |
| `goslim_mf_name` / … | its label, e.g. `peptidase activity` |
| `goslim_mf_all` / … | the complete multi-label set, `;`-joined |
| `goslim_source` | `curated` \| `eggnog` \| empty |
| `eggnog_og`, `eggnog_evalue`, `eggnog_score` | evidence for `eggnog` rows |

A protein can have MF and no BP. Absence of an aspect is a fact, not a missing value.

`eggnog_<species>.tsv` is a first-class artifact from the same run, carrying everything else
eggNOG-mapper produced — `kegg_ko`, `kegg_pathway`, `kegg_module`, `brite`, `ec`, `cazy`, `pfams`,
`preferred_name`, `description`, and its own `cog_category` (an independent second opinion on
Part 1's letter).

**Load through `src/function.py`.** `load(species)` returns both schemes joined — the usual entry
point; `load_cog`, `load_goslim`, `load_eggnog` read one table each, `load_goslim_terms` the
97-term vocabulary, and `confident` / `informative` are the two filters described above
(`confident` keeps curated rows plus eggNOG rows below an e-value cut; `informative` drops COG `R`
and `S`). Everything reads with `keep_default_na=False`, so a missing value is `""` and not NaN.

## The control

Same trick as Part 1, and scored **identically to the method it replaced** so the numbers are
directly comparable: E. coli has curated GO for ~90% of its proteome, so for every protein with both
curated GO and an eggNOG-mapper annotation, compare the two per aspect. n = 3,445 proteins.

| method | MF | BP | CC |
|---|---|---|---|
| **eggNOG-mapper** | **80.8%** (J 0.714) | **82.4%** (J 0.593) | **81.1%** (J 0.750) |
| ESM-C k-NN — rejected, reweighted | 68% | 53% | 66% |
| InterPro2GO — rejected | 77.4% (J 0.681) | 60.5% (J 0.445) | 13.1% (J 0.121) |

**When eggNOG-mapper answers, it answers better than anything else tried.** It just rarely answers
where UniProt has not already — see the run log.

## Environment and its traps

`gradi-emapper` (osx-64, Rosetta). eggnog-mapper has no osx-arm64 build, and installing it into
`gradi` would drag that env to osx-64 and take ESM-C with it — the same reasoning as `blast`.

**`goslim.py` needs `goatools`** in `gradi` — it is what supplies `GODag` and `mapslim`, i.e. the
whole slim mapping. Nothing else in this stage depends on it.

1. **bioconda puts `diamond`/`mmseqs` in the env's `bin/`**, but emapper looks for them inside
   `site-packages/eggnogmapper/bin/`. Symlink them; `install.sh` records the loop.
2. **`download_eggnog_data.py` does not work.** It fetches from `eggnogdb.embl.de`, which no longer
   resolves — and then prints `Finished.` and exits **0** having downloaded nothing. The live host is
   `eggnog5.embl.de`.
3. **`eggnog5.embl.de` drops long transfers.** The first attempt died at 16% with
   `curl: (18) transfer closed with 5689882877 bytes remaining`. Fetch resumably (`curl -C -`) and
   loop until the byte count equals `Content-Length`. An exit code is not evidence of data. **And do
   not reach for curl's own `--retry`**: it restarts a `-C -` transfer from zero, so the loop has to
   be yours, re-invoking `curl -C -` each time.
4. The files unpack to **50.6 GB**, of which `eggnog.db` alone is **41.4 GB** and
   `eggnog_proteins.dmnd` 9.3 GB (`data/source/eggnog/SOURCE.md` holds the byte counts). **~21 GB was
   the planning estimate and it was wrong by more than half** — measure the *unpacked* size against
   free disk before committing to a database, not the download size. Public and re-derivable —
   **do not upload it to eosvc.**
5. Making room for it meant deleting `data/raw/other/chembl/chembl_37/` (28 GB), which
   `PROVENANCE.md` Rule 1 documents as re-derivable; `chembl_37_sqlite.tar.gz` was verified with
   `tar -tzvf` (exit 0, correct member, 30,480,314,368 bytes) **before** the directory was removed.

## Run log — 2026-09-01/02

emapper-2.1.15 · eggNOG DB 5.0.2 · diamond 2.0.15 · diamond mode, 9 cpu · **15 min** for all three
proteomes (5.5 / 4.8 / 4.8 min). The 50.6 GB database took ~4 h to fetch.

### eggNOG-mapper's own output

| species | orthologous group | GO | KEGG ko | EC | preferred name |
|---|---|---|---|---|---|
| Kp HS11286 | 91.5% | 46.4% | 64.1% | 28.6% | 67.0% |
| Ec K-12 | 96.3% | 78.7% | 74.5% | 30.7% | 92.9% |
| Sa NCTC 8325 | 89.1% | 20.1% | 55.1% | 27.9% | 56.4% |

Kp reaching 91.5% confirms the prediction that its 0% eggNOG *xref* was an artefact of being dark
TrEMBL, not a real absence — eggNOG-mapper computes the group from sequence and finds one for nearly
every protein.

### GO-slim coverage: the gain is small, and here is exactly why

| species | curated | + eggNOG | gain |
|---|---|---|---|
| Kp | 69.7% | **73.7%** | +4.0 pp (231 proteins) |
| Ec | 87.0% | **88.7%** | +1.7 pp (75) |
| Sa | 65.0% | **65.6%** | +0.6 pp (16) |

**322 proteins out of 13,020.** The bottleneck is not orthology assignment:

| | lacked GO | got an eggNOG group | **that group carried GO** | survived the slim |
|---|---|---|---|---|
| Kp | 1,737 | 1,261 (73%) | **249 (14%)** | 231 |
| Ec | 571 | 453 (79%) | **101 (18%)** | 75 |
| Sa | 1,010 | 703 (70%) | **24 (2%)** | 16 |

eggNOG places ~three-quarters of these proteins in groups perfectly well. **The groups themselves
have no GO.** This is the same fact that killed InterPro2GO, a third time: these proteins are
genuinely unannotated in the source databases and every tool inherits that. Ceiling estimates made
beforehand from eggNOG *xref* coverage (96% Ec / 82% Sa) were far too optimistic because they
assumed groups carry GO. **Test the OG→GO yield before paying for a database this size.**

### Three-way COG cross-check

`evidence/cog_crosscheck.tsv` — the same run emits `COG_category`, so Part 1's tool can be checked
against a real competitor on E. coli, with NCBI's curated COG2024 as ground truth:

| | agreement vs NCBI | coverage |
|---|---|---|
| **COGclassifier** (Part 1) | **97.8%** | 84.5% |
| emapper `COG_category` | 63.3% | 94.1% |

The two tools agree with each other only 62.7%. emapper's letters come from eggNOG's own OG→category
mapping, not from NCBI's COG2024 assignment — a legitimately different thing, but **not** the
vocabulary Part 1 is committed to. COGclassifier stays.

### Was it worth it?

For GO alone: 50.6 GB and 28 GB of deleted ChEMBL bought **+4.0 pp on the anchor organism**. Poor
value, and cheaply foreseeable — the OG→GO yield could have been tested first.

The tier is kept anyway because the additions are real, they score better than any alternative
tried, and the same run produced things nothing else provides: **KEGG KO at 64.1% on Kp** (entirely
new — UniProt carries only KEGG *gene* ids), preferred names at 67.0% (above stage 00's 63.4%
`gene_name` coverage, so a candidate source for closing that gap), and the COG cross-check above.

**The database was deleted after the run** (see `data/source/eggnog/SOURCE.md`, which keeps
the byte counts, the recovery procedure and `fetch_eggnog.sh`). Every derived artifact is saved;
only re-running needs the 50.6 GB back.

## Traps

1. **Coverage is not 100% and must not be made to look like it.** An empty row means no established
   tool could annotate that protein.
2. **`eggnog` rows are a tool's inference, not curation.** `src.function.confident()` keeps curated
   rows plus eggNOG rows below an e-value cut.
3. **`GODag` keys include `alt_id`s.** Counting the slim by iterating the DAG gives 153 terms, not
   97. Dedupe on `node.id`.
4. **Curated coverage < raw GO coverage** — GO that maps to no slim term is not a fill (186 in Kp).
5. **Obsolete GO ids**: UniProt xrefs lag the GO release; `mapslim` raises on retired ids, caught and
   treated as unmappable.
6. **`interpro2go` is downloaded but is not a tier** — it scores a rejected alternative every run
   under REJECTED ALTERNATIVE. Do not re-add it to the chain.
7. **Do not reintroduce the ESM-C transfer.** See above.

## Running it

```bash
python scripts/function/eggnog.py                 # the tool -- eggnog_<species>.tsv
python scripts/function/goslim.py                 # the vocabulary -- goslim_<species>.tsv
python scripts/function/goslim.py --no-eggnog     # curated tier only
```

CLI (goslim): `--species` · `--no-eggnog` · `--no-control` · `--refresh` · `--dry-run` · `-q`.
CLI (eggnog): `--species` · `--cpu` · `--limit` · `--refresh` · `--dry-run` · `-q`.


# Running the COG stage

```bash
python scripts/function/cog.py                                # three bacteria + control, ~9 min cold
python scripts/function/cog.py --species saureus --limit 50   # smoke test -> scratch/smoke_*
python scripts/function/cog.py --refresh                      # re-download and re-run RPS-BLAST
python scripts/function/cog.py --dry-run
python scripts/plots/function.py
```

CLI: `--species` · `--evalue` · `--threads` · `--limit` · `--refresh` · `--no-control` ·
`--dry-run` · `-q/--quiet`.

The script exits non-zero if a spot check fails or if E. coli category agreement drops below 90%.

---

## The deliverable: two complete matrices

`scripts/function/matrix.py` reduces this stage to **one packed table per species, with both
matrices beside it**. It **recomputes nothing** — it reshapes `cog_<species>.tsv` and
`goslim_<species>.tsv`, and runs in seconds.

    function_<species>.tsv          uniprot_ac · cog_categories · goslim_terms · function_evidence
    evidence/goslim_matrix_<sp>.tsv uniprot_ac + 97 GO-slim term columns + evidence   (n, 99)
    evidence/cog_matrix_<sp>.tsv    uniprot_ac + 26 COG letter columns   + evidence   (n, 28)

### `function_evidence`, and why there is no `function_consensus`

The axis's half of the standard pair (`src/consensus.py`). There is no consensus column because
"how much function does a protein have" is not a quantity — the nearest candidate, annotation
richness, measures how well *studied* it is. **Evidence always, consensus where the axis has a
magnitude.**

| level | rule | Kp | Ec | Sa |
|---|---|---|---|---|
| 3 | both schemes, GO **curated**, COG **informative** (a letter outside `R`/`S`) | 3,697 | 3,299 | 1,682 |
| 2 | at least one scheme, but a quality test fails | 1,012 | 776 | 535 |
| 1 | neither scheme annotates it | 1,019 | 328 | 672 |

**Level 3 means corroborated**, matching `essentiality_evidence`. The two schemes are genuinely
independent — COGclassifier runs rpsblast against CDD profiles; GO-slim comes from UniProt
curation or eggNOG. emapper's own `COG_category` is *not* a third opinion: 63.3% agreement with
NCBI's curated COG2024 against COGclassifier's 97.8%.

**Level 2 is "annotated, not corroborated" — never "badly annotated".** A protein with excellent
curated GO but no COG hit caps at 2, and COG coverage is bounded near 81.6% by NCBI's own curators
on *E. coli*, so a missing COG is usually the method's ceiling.

**It is not a fame measure.** 210 *E. coli* proteins named "Uncharacterized" sit at level 3 —
32.4% of that group, against 3.3% on Kp — because UniProt leaves them unnamed while curating their
class: "Uncharacterized MFS-type transporter YhhS" carries `GO:0005215`. The column measures
annotation support, not how much anyone has written about the protein.

**The signal we would rather have does not exist here.** GO evidence codes — experimental
`EXP`/`IDA`/`IMP` against electronic `IEA` — are absent: `goslim_<species>.tsv` carries only
`goslim_source`, `data/source/go/` holds the OBO files and `interpro2go` but no GAF, and UniProt's
`go_id` xref is a bare list. Obtaining them means a per-proteome GOA download.

The ladder also keeps the **curated-vs-eggNOG** distinction alive after `goslim_evidence` was
dropped from the shipped table for separating only 322 of 13,020 proteins.

### Why both forms ship

The packed table became the deliverable on the project owner's instruction, **2026-10-03**. Both
term columns are `;`-joined in vocabulary order, and an **empty string means no term** — a complete
row, not a missing one.

**The matrices stay because a packed list cannot express a structural zero.** It cannot distinguish
a term that is merely unannotated from one the organism cannot reach:

| | GO-slim columns always zero | why |
|---|---|---|
| Kp | 16 / 97 | 8 eukaryote/plant concepts + 8 not seen here |
| Ec | 8 / 97 | the eukaryote/plant 8 |
| **Sa** | **19 / 97** | the 8, plus 11 Gram-negative envelope terms it structurally lacks |

COG `Y` (nuclear structure) is always zero in all three, for the same kind of reason. Those are
**kept columns** in the matrix and simply absent from the packed column, so anything that needs to
tell *impossible* from *unknown* — or that wants an `hstack`-able feature matrix — reads
`load_goslim_matrix()` / `load_cog_matrix()`.

**The two forms are provably interchangeable, asserted in both directions**: the matrices
round-trip to the long-form source, and the packed columns re-expand to the matrices exactly. A
shape check would pass on a wrong table; these do not.

**Neither scheme ships an evidence column** (owner's call, 2026-10-03).

COG never had one to carry: `cogclassifier` vs `none` is 1:1 with non-empty vs empty on all three
species — measured, not assumed — so it only restated the term column.

GO-slim's was a real distinction, but a thin one: `curated` vs `eggnog` separates **322 proteins
out of 13,020** — Kp 231, Ec 75, Sa 16 — and it survives **byte-identically in two files that stay
on disk and stay audited**, as `evidence` in `evidence/goslim_matrix_<species>.tsv` and as
`goslim_source` in `evidence/goslim_<species>.tsv`.

**The cost, stated plainly:** the shipped table no longer says whether a protein's GO terms came
from UniProt curation or from an eggNOG orthogroup. That matters for the eggNOG tier specifically,
which is orthology-inferred rather than curated — read one of those two files before treating a
term as curated.

| species | goslim annotated | all-zero | cog annotated | all-zero |
|---|---|---|---|---|
| kpneumoniae | 4,222 (73.7%) | 1,506 | 4,528 (79.1%) | 1,200 |
| ecoli | 3,907 (88.7%) | 496 | 3,719 (84.5%) | 684 |
| saureus | 1,895 (65.6%) | 994 | 2,112 (73.1%) | 777 |

**Complete by construction**: one row per protein, always. An unannotated protein is an all-zero row
with `evidence == "none"`, never a missing row.

### Why multi-label, and what it recovers

Built from the `*_all` columns, not the single chosen term. That is the main gain, and it is large:

| | proteins with >1 term | max |
|---|---|---|
| GO slim MF | 34–41% of MF-annotated | 7 |
| GO slim BP | 32–52% of BP-annotated | 7 |
| COG letters | 12.4–12.6% of classified | 4 |

Reconciled exactly against `evidence/cog_counts_<species>.tsv`, which counts only the **chosen**
letter: *chosen + extra-from-multi-label == matrix total*, to the unit — Kp 4,528 + 591 = 5,119,
Ec 3,719 + 477 = 4,196, Sa 2,112 + 255 = 2,367. The gap between those two numbers is precisely the
information the old headline column discarded.

### Columns that are always zero are information, not padding

All 97 slim terms and all 26 COG letters appear in every species file, so the three stack without
alignment. Always-zero columns per species: goslim **Kp 16 · Ec 8 · Sa 19**; cog **1 everywhere**.

- **8 slim terms never occur anywhere** — thylakoid, photosynthesis, extracellular matrix, protein
  tag activity, histone binding, extracellular structure organization, nutrient reservoir, cell
  adhesion mediator. All eukaryote/plant concepts.
- **COG `Y` (nuclear structure) is always zero** — no bacterium has one.
- ***S. aureus* has 19 always-zero slim columns**, more than the Gram-negatives, because it is
  Gram-positive. That is biology showing through the schema.

Dropping these would give the species files different shapes and make them un-stackable.

### A zero means NOT ANNOTATED, not absent

For **26.3% of K. pneumoniae** an all-zero row means nothing is known, not that the function was
ruled out. `evidence` (`curated` | `eggnog` | `none`; `cogclassifier` | `none`) is what separates
the two. Downstream must not read a zero as a measured negative.

### What the two-matrix deliverable drops

`eggnog_<species>.tsv` **stays on disk but is no longer a deliverable**. It cannot be deleted — the
goslim `eggnog` tier is derived from its `gos` column, so removing it makes this stage
non-regenerable without re-running the 50.6 GB database. Leaving the deliverable set: `pfams`
(85.4% Kp), `kegg_ko`/`brite` (64.1%), `description` (86.8%), `preferred_name` (67.0%), `ec`
(28.6%), `kegg_pathway`, `kegg_module`, `cazy`, and emapper's independent COG call, which covers
**480 Kp proteins COGclassifier misses**. Recorded so the loss is deliberate and findable.

### The COG vocabulary is vendored, on purpose

The 26 letters were read at runtime from COGclassifier's **installed package resource**. A matrix
whose column schema is defined by a pip install is not reproducible, so it is now vendored at
`data/source/cdd/cog_func_category.tsv` with a `SOURCE.md`. Column order follows that file
(grouped: information storage → cellular processes → metabolism → poorly characterised), not
alphabetical. **The file has no trailing newline, so `wc -l` reports 25 for 26 rows.**

### Verification

The real check is a **round-trip**, not a shape check: reconstructing the `;`-joined term lists from
the matrix reproduces `goslim_*_all` and `cog_category_all` exactly, for every protein. Passing for
all 13,020. Shape checks pass on a wrong matrix; this does not.

**Load through `src/function.py`** — `load(species)` for the packed deliverable, `load_long` for
the joined long-form source, `load_goslim_matrix` / `load_cog_matrix` (+ their `_all` variants) for
the matrices, `matrix_manifest`.

CLI: `--species` · `--dry-run` · `-q`. There is no `--refresh` and no `--limit`: the script reads
two finished tables and reshapes them in seconds, so a re-run is the refresh.
