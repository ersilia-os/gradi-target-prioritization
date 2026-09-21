# Essentiality screens: what we have, what is new, what needs a human

Focus is ***E. coli* and *K. pneumoniae***. *S. aureus* findings are recorded at the end but were
deprioritised on instruction.

Every claim here is measured from the file on disk or verified against a live URL. Where a number
came from a paper rather than the data, it says so.

---

## 1. In use as training endpoints (6)

One endpoint per **source**, conditions aggregated within a source. Join rates are re-measured on
every run and written to `data/processed/essentiality/evidence/screen_join_audit.tsv`.

| endpoint | organism / strain | assay | rows | pos | base | join |
|---|---|---|---|---|---|---|
| `keio_ess` | *E. coli* K-12 (PEC/Keio) | arrayed knockout | 4,190 | 286 | 0.068 | 97.1% b-number |
| `goodall_ess` | *E. coli* BW25113 | TraDIS | 4,056 | 354 | 0.087 | 97.6% symbol→b-number |
| `bn373_ess` | ***K. pneumoniae* ECL8** | TraDIS/DESeq | 4,930 | 523 | 0.106 | 97.7% GenBank tag |
| `kpnih1_ess` | ***K. pneumoniae* KPNIH1** | Tn-seq | 5,485 | 412 | 0.075 | **424/424 listed** |
| `conservation` | Gammaproteobacteria | clade aggregate | 4,256 | core 205 / mid 530 / non 3,521 | — | 99.0% |
| `bw25113_ess` *(control)* | *E. coli* BW25113 | TraDIS/DESeq | 4,199 | 258 | 0.061 | 98.7% |

**Why each matters**

- **`keio_ess`** — the only **non-transposon** assay we hold. Every other screen infers essentiality
  from missing insertions; Keio actually tried to build each knockout and recorded which could not
  be made. Different failure modes, so it is genuinely independent evidence rather than a replicate.
- **`goodall_ess`** — the standard modern E. coli TraDIS reference, and the screen whose insertion
  index defines what "essential" means in most later papers.
- **`bn373_ess`** — **our only measured *Klebsiella* essentiality.** Before this the anchor organism
  had a 100%-predicted column. It is the compendium's uniform DESeq reprocessing of ECL8 TraDIS.
- **`kpnih1_ess`** — a **second, independent Kp genetic background** (ST258, the carbapenem-resistant
  lineage) against ECL8's K2-ST375. Two backgrounds let us ask whether an essential call is
  lineage-specific rather than assuming it.
- **`conservation`** — how conserved essentiality is across 12 Gammaproteobacteria genomes, as
  core / mid / non. The only endpoint whose label is a property of the gene *family* rather than one
  organism, which is why a per-protein embedding is straightforwardly the right feature for it.
- **`bw25113_ess`** — demoted to control: same compendium as `bn373_ess`, so it is the cleanest
  possible partner for a cross-species test (identical analysis pipeline on both sides).

**A finding that shapes how these should be read.** `goodall_ess` and `bw25113_ess` are the *same
strain and same assay*, yet agree at Jaccard **0.554** — while `keio_ess`, a completely *different*
assay, agrees with `goodall_ess` at **0.730**. The analysis choice moves the essential call more
than the assay does. The same pattern appears in ECL8 (authors' call 373 positives vs the
compendium's DESeq 562 on the same experiment). **Do not treat two reprocessings of one screen as
two opinions.**

---

## 2. Downloaded today, not yet endpoints

### *E. coli*

**Nichols 2011, Cell 144:143-156 — "Phenotypic landscape of a bacterial cell"**
`data/raw/ecoli/essentiality/nichols2011_chemgen/1-s2.0-S0092867410013747-mmc2.xls`
Verified: **3,980 genes × 325 columns** of continuous S-scores, sheet `TableS2-FinalData`, rows
keyed `ECK####-GENENAME`.
*Why it matters*: the widest condition-resolved fitness matrix for E. coli in existence, and the one
resource that can answer **"essential under which stress"** rather than "essential". v1 recorded it
as unobtainable — PMC serves it behind a reCAPTCHA. The Elsevier CDN is unprotected; the PII
(`S0092867410013747`) is not derivable from the DOI, which is why it stayed lost for a year.

**Choe 2025, iScience — CRISPRi across 13 conditions**
`data/raw/ecoli/essentiality/choe2025_ecoli/mmc3.xlsx` — **4,199 rows × 20**, continuous `ER`
(enrichment ratio), carries `b number`, `KEIO*` and `COG†` columns. **Two-row header.**
*Why it matters*: condition-resolved essentiality with a native b-number column (99.1% join), from a
different technology (CRISPRi) than our transposon screens.

**Choe 2023, mSystems — Tn-seq in LB *and* M9**
`data/raw/ecoli/essentiality/choe2023_ecoli/msystems.00896-22-s0002.xlsx` — **4,499 rows × 20**,
`IPKM` per medium.
*Why it matters*: two media is the minimum needed to separate "essential" from "essential in rich
medium" — the confound sitting under every single-condition screen we hold.

### Cross-organism

**Fitness Browser, February 2024 release** (figshare `10.6084/m9.figshare.25236931`, **CC BY 4.0**)
`data/raw/other/essentiality/fitness_browser_2024/` — `feba.db.gz` (2.3 GB → 7.4 GB sqlite),
`aaseqs.gz` (**221,030 protein sequences**), plus per-organism strain-fitness tables.
Schema: `GeneFitness(orgId, locusId, expName, fit, t)`; `Keio` = *E. coli* BW25113, 168 experiments,
3,789 genes.
*Why it matters*: the only large corpus that ships **graded fitness with sequences attached** — no
identifier join at all, directly embeddable, and larger than the whole DEG corpus (173,048).

> **Three caveats, all load-bearing.** (1) **RB-TnSeq structurally cannot see essential genes**: a
> gene with no surviving insertions has no fitness value, so it is *absent* rather than extreme.
> This is a FITNESS resource; an essentiality endpoint built from it would measure the wrong thing.
> (2) **There is no *K. pneumoniae* in it.** The single Klebsiella, `Koxy`, is
> *K. michiganensis* M5al per the database's own Organism table — a comparator, never a Kp label.
> (3) **I downloaded the wrong file first.** `feba.db` (2.3 GB → 7.4 GB) contains
> `GeneFitness(orgId, locusId, expName, fit, t)` and **no essentiality signal at all**. The actual
> essentiality calls are a 2.2 MB file on a different, un-gated host:
> `https://genomics.lbl.gov/supplemental/bigfit/essential_proteins.tab` — 13,869 calls across 32
> organisms with `nPosCentral` and insertion `dens`, **positives only** (negatives = genome minus
> the list). Staged alongside. Keep `feba.db` only if we want conditional fitness; for essentiality
> the small file is the whole payload.
>
> A newer release also exists (July 2026: 62 organisms, 9,704 experiments, figshare
> `10.6084/m9.figshare.32865896`) — ours is the February 2024 one.

**OGEE v3 — recovered, and it overturns a documented project belief**
`data/source/ogee/gene_essentiality.txt.gz` — md5 `b42f4a3358484b490e8bffe8c90edd00`.
CLAUDE.md recorded OGEE v3 as permanently lost: the server aborts its TLS handshake, verified
against three TLS stacks *and* real Chrome, with no mirror found. That diagnosis was correct and the
conclusion was still wrong — **nobody checked the Internet Archive.** Two independent mirrors return
byte-identical content. The `id_` suffix on the Wayback URL is mandatory.

Measured: **255,162 rows, 89 taxa, 124 datasets**; E. coli 16,979 rows / 1,523 essential,
S. aureus 5,612 / 671, **K. pneumoniae ZERO**. Against the DEG corpus that is ~2× the species and
~2.2× the positives, with the negative class shipped explicitly rather than reconstructed.
CC BY 3.0.

*Why it matters*: it roughly doubles the cross-species corpus behind `conservation`. *Why it does
not solve the main problem*: Kp is absent here exactly as it is from DEG. **Only the Tn-seq/TraDIS
screens close the anchor's gap.**

*The lesson*: a dead server is not a lost dataset. Check Wayback with `id_`, and check for a GitHub
mirror, before recording anything unobtainable.

---

## 3. Downloaded and verified — *K. pneumoniae*, seven new strains

All fetched via Europe PMC or a publisher CDN, row counts verified by reading the file.

| screen | strains | rows | readout | locus tags → assembly |
|---|---|---|---|---|
| **Bruchmann 2021** NAR | RH201207 | 5,390 | 3-state **+** logFC/q | `KPNRH_` → FR997879 |
| | ATCC 43816 | 5,217 | 3-state + logFC/q | `VK055_RS` → NZ_CP009208.1 |
| **Short 2020** IAI | B5055 | 5,309 | 3-state + contrasts | `BN49_` → FO834906 |
| | NTUH-K2044 | 5,413 | " | `KP1_` → AP006725 |
| | ATCC 43816 | 5,191 | " | `VK055_` → CP009208 (**GenBank, not _RS**) |
| | RH201207 | 5,726 | " | `RH201207_` → LT216436 |
| **Rome 2026** AAC | NJST258_2 (ST258) | 5,125 | TRANSIT ES/GD/NE/GA + log2FC | `KPNJ2_RS` → GCF_000597905.1 |
| **Gray 2024** eLife | 4-strain bridge | 385/642/608/476 | binary | see caveats |

**Why this is the important block.** DEG contains **zero** *Klebsiella*. Before today the anchor had
one measured screen; it now has **seven genetic backgrounds**, including three independent ST258
(carbapenem-resistant lineage) reads. That is enough to ask whether an essential call is
lineage-specific rather than assuming it.

**Bruchmann is the pick of them** — plain CSV, three-state call *and* continuous logFC/q, a built-in
cross-strain locus map, and COG letters, on two backgrounds we had nothing for.

### Two assembly traps, both verified against the deposited records

- **The same strain appears under two incompatible assemblies.** RH201207 is `RH201207_*` on
  `LT216436` in Jana/Short, and `KPNRH_*` on `FR997879` in Bruchmann/Gray — **different assemblies
  whose chromosomes differ by 901 bp. Do not merge them.**
- **ATCC 43816 appears as both `VK055_` (GenBank CP009208) and `VK055_RS` (RefSeq
  NZ_CP009208.1).** Not interchangeable — the same trap that cost hours on ECL8 today.

### Dropped after checking

**Jung 2019 (MH258, 5,473 genes)** — its ENA project `PRJEB31265` **exists but has released zero
records**, checked across `assembly`, `wgs_set`, `sequence`, `read_run` and `analysis`. Its
`gene_####` IDs are PATRIC output with no public counterpart and its coordinates are against an
unavailable assembly. Joinable only via free-text description. **Not worth it.**

## 4. Needs a human

**`paczosa2020_kppr1`** — *K. pneumoniae* ATCC 43816, neutropenic vs WT mouse, 166/194 genes
depleted. Readout is a **normalized competitive index** (`WT-mean-nCI`, `PMN-mean-nCI`), not logFC,
so it needs converting before mixing with other screens.

Confirmed unreachable by automation: Europe PMC reports not-OA, curl 403s, **and an in-page
`fetch()` from an already-loaded ASM tab also returns 403** — Cloudflare blocks XHR separately from
navigation. No SRA/BioProject exists; the reads were never deposited.

> **Navigate a browser directly to these two URLs** (do not try to script them):
> - `https://journals.asm.org/doi/suppl/10.1128/iai.00034-20/suppl_file/iai.00034-20-sd001.xlsx` (657 KB, Tables S1–S7)
> - `https://journals.asm.org/doi/suppl/10.1128/iai.00034-20/suppl_file/iai.00034-20-s0002.pdf` (642 KB)
>
> Save both to `data/raw/kpneumoniae/essentiality/paczosa2020_kppr1/`.

**`mazzuoli2025_saureus`** (bioRxiv) — 4 strains incl. NCTC8325-4, 307-gene core essentialome.
bioRxiv exposes only 3 PDF figures and `aureobrowse.veeninglab.com` is a Shiny shell with no static
export. Re-check after journal publication, or email the lab.

---

## 5. Measured and rejected — do not re-derive

- **`ecl8_ess` (Eichelberger's own call) is unjoinable.** Their `ecl8_#####` tags annotate an
  assembly nobody deposited: **475 of 5,165** numeric suffixes shared with `BN373_`, **0 of 5,074**
  coordinates matched, gene symbols reached **3.7%** (the deposited ECL8 annotation carries symbols
  for only 144 CDS). Superseded by `bn373_ess` — same organism, same assay, 97.7% join.
- **`cain2017_njst258` does not exist.** No NJST258 TraDIS screen was ever published; NJST258_1/2
  are assemblies only. The real Cain 2017 paper is RH201207 (`srep42483`) and its sole supplement is
  a 6-page PDF with no per-gene table.
- **Rosconi 2022** is *S. pneumoniae*, not *S. aureus*.
- **Bae 2004** is an arrayed *C. elegans* killing screen — binary virulence, not fitness.

---

## 6. Download routes that work, and three that lie

| route | verdict |
|---|---|
| `ebi.ac.uk/europepmc/webservices/rest/<PMCID>/supplementaryFiles` | **the default** — publisher's own ZIP, one curl, no browser |
| `ars.els-cdn.com/content/image/1-s2.0-<PII>-mmcN.xls` | works; CDN unprotected while cell.com is gated |
| `static-content.springer.com/esm/...` | works |
| `journals.plos.org/.../article/file?id=...&type=supplementary` | works |
| `journals.asm.org/doi/suppl/...` | **403 Cloudflare** — use Europe PMC, or a real browser |
| `pmc.ncbi.nlm.nih.gov/articles/instance/<n>/bin/<file>` | **serves a reCAPTCHA page as HTTP 200**, 21 KB of HTML |
| figshare `ndownloader.figshare.com/files/<id>` | works **only with a short or absent User-Agent** — a full Chrome UA gets HTTP 202 forever |

**Three traps, all of which cost real time today:**

1. **An HTTP 200 is not evidence of data.** A Nature article URL returned 407 KB of HTML and was
   recorded as a successful download until the content was checked. The fetcher now rejects anything
   that parses as markup, regardless of size.
2. **Publishers and repositories want opposite headers.** Publishers 403 a bare client; figshare
   202s a browser-like one. Try both.
3. **Check that a DOI resolves to the paper you think it does.** `10.1038/srep44322` was carried in
   a note as "Cain 2017 K. pneumoniae"; it is an unrelated neurovascular paper.

---

## 7. *S. aureus* — recorded, deprioritised

Not pursued on instruction, but the research found the gap is closeable and the best item is
specific: **Santiago 2018 Nat Chem Biol**, a **2,870 gene × 91 drug-condition continuous matrix** on
an NCTC 8325-derived strain whose `Gene` column is literally `SAOUHSC_00001` (96.4% join to our
proteome). Verified download (2,009,555 bytes):
`https://static-content.springer.com/esm/art%3A10.1038%2Fs41589-018-0041-4/MediaObjects/41589_2018_41_MOESM3_ESM.zip`
Nature's own article page is auth-gated and the `.xlsx` variants 403 — only this `.zip` works.

Also identified: **Liu 2024 mSystems** (CRISPRi essentialome on NCTC8325-4, `PMC11265419`) and
**Mårli 2025** (same library, BHI vs milk). Checked for redundancy rather than assumed: Spearman
ρ 0.39–0.41 between the two papers' log2FC, Jaccard 0.58 on essential calls — **semi-independent
replicates of one library, not two opinions.**
