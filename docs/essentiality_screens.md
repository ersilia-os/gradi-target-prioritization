# Essentiality screens: what is used, what is held back, and what must not be re-derived

Focus is ***E. coli* and *K. pneumoniae***. *S. aureus* findings are recorded at the end but were
deprioritised on instruction.

Every claim here is measured from the file on disk or verified against a live URL. Where a number
came from a paper rather than from the data, it says so.

**Most of the expensive knowledge in this document is NEGATIVE** — a dead end, a wrong DOI, a file
that is not what its size suggests, a label that is inverted. That is exactly what gets
rediscovered at cost, so sections 2, 5 and 6 are the ones to read before starting any new search.

Per-directory provenance lives beside the data, in a `SOURCE.md` in each dataset folder.

---

## 1. In use — ten training sets

One column per **source**, conditions aggregated within a source. Nothing is merged: these are
different strains, different assays and different analyses, and the base rate alone spans
**0.048–0.194** across the set. That spread is information about method and strain, not noise to
average away.

Join rates are re-measured every run into
`data/processed/essentiality/evidence/screen_join_audit.tsv`; the full table with comments is
`output/results/essentiality/screen_summary.tsv`.

| column | organism / strain | assay | n | pos | base | join | ribo |
|---|---|---|---|---|---|---|---|
| `essential_kpneumoniae_ecl8_tradis` | **Kp ECL8** (K2-ST375) | TraDIS/DESeq | 4,930 | 523 | 0.106 | 97.7% `BN373_` tag | 0.942 |
| `essential_kpneumoniae_rh201207_tradis` | **Kp RH201207** | TraDIS | 4,981 | 471 | 0.095 | 92.9% locus_tag | 0.938 |
| `essential_kpneumoniae_atcc43816_tradis` | **Kp ATCC 43816** | TraDIS | 4,809 | 363 | 0.075 | 92.6% locus_tag | 0.906 |
| `essential_ecoli_k12_knockout` | Ec K-12 (PEC/Keio) | **arrayed knockout** | 4,190 | 286 | 0.068 | 97.1% b-number | 0.774 |
| `essential_ecoli_mg1655_footprinting` | Ec MG1655 (DEG1018) | genetic footprinting | 4,253 | 604 | 0.142 | 98.5% exact sequence | 0.538 |
| `essential_ecoli_bw25113_tradis_goodall` | Ec BW25113 | TraDIS | 4,056 | 354 | 0.087 | 97.5% symbol→b-number | 0.961 |
| `essential_ecoli_bw25113_tnseq_choe` | Ec BW25113 | Tn-seq (LB only) | 4,272 | 440 | 0.103 | 95.0% b-number | 0.906 |
| `essential_ecoli_st131_tradis` | **Ec ST131 EC958** | TraDIS | 4,981 | 300 | 0.060 | 100% DEG `fasta_id` | 0.741 |
| `essential_ecoli_o157h7_tnseq` | **Ec O157:H7** | Tn-seq | 5,433 | 1,055 | 0.194 | 100% DEG `fasta_id` | 1.000 |
| `core_essential_gammaproteobacteria` | Gammaproteobacteria | clade conservation | 4,256 | 205 | 0.048 | 99.0% b-number | 0.906 |

`ribo` is the label-polarity control — see §1.3.

### 1.1 The three filters that define this set

Each was applied on instruction from the project owner.

1. **Both classes required.** A positives-only gene list is not a training set. Drops Ramage 2017
   (KPNIH1, 424 genes) and Paczosa 2020 (310 hits).
2. **No condition-dependent data, for now.** Drops Choe's M9 arm, Rome 2026 (iron-depleted),
   Short 2020's four serum screens, Bruchmann's `2hpi`/`6hpi` in-host columns, and every in-vivo
   mouse screen.
3. **No duplicate experiments.** DEG1019 is the Keio collection again under DEG's curation; PEC's
   is kept.

### 1.2 What the spread says, and why nothing is merged

Base rates run **0.060 to 0.194** — a 3.2× spread. The two extremes are the oldest assay (Gerdes
2003 genetic footprinting, 0.142) and a different pathotype (O157:H7, 0.194), so the variation is
method and strain, not measurement error.

The sharpest illustration is **one organism, one dataset, two analyses**: Goodall's own BW25113
calls give 354 positives, and the compendium's DESeq reanalysis of *the same reads* gives 258. A
merged E. coli essentiality label would silently pick one of those.

### 1.3 The polarity control — `scripts/essentiality/summary.py`

`num_positives` looks identical whether a label set is right or inverted, and an inverted column
trains a confident, well-formed, exactly wrong model. So every column is checked against biology
that cannot be in dispute: **ribosomal proteins must be essential** (`rps*`/`rpl*`/`rpm*`) and the
**textbook dispensables must not** (`lacZ`, `araB`, `fadB`, flagellar, fimbrial, sugar catabolism).

**Both bars are calibrated on measurements, not on 1.00 — this was got wrong once.** An initial
0.50 separation bar failed Gerdes and would have failed the gold standard too:

| screen | ribosome recall | reading |
|---|---|---|
| Goodall 2018 TraDIS | 0.961 | the cleanest E. coli screen here |
| **Keio arrayed knockout** | **0.774** | **the gold standard, and the realistic ceiling** |
| Gerdes 2003 footprinting | 0.538 | noisy, correctly polarised (8× over dispensables) |
| Ghomi BW25113 DESeq | **0.113** | **broken — retired, see §5** |

Keio's misses (`rplA`, `rplI`, `rplK`, `rplY`, `rpmE/F/G/I`, `rpsF/O/T/U`) are genuinely
dispensable in *E. coli*, so ~0.8 is the biological ceiling. The bars are set to catch **breakage,
not noise**: recall ≥ 0.45, separation ≥ 0.40.

### 1.4 Identifier notes specific to these ten

- **GenBank vs RefSeq is per-paper, not a default.** ECL8 and KPNIH1 need GenBank (PGAP dropped the
  submitter tags: 263 of RefSeq's 281 ECL8 misses were tags absent entirely). Paczosa keys on
  `VK055_RS*`, which **is** RefSeq, and ATCC 43816 uses the `GCF_` assembly for that reason.
- **NCBI GFF splits what you need across two features**: `old_locus_tag` sits on the **gene**,
  `protein_id` on the **CDS**, linked by `ID`/`Parent`. Walk both or the map comes out empty.
- **The same strain can need different records per paper.** RH201207 is `RH201207_*` on `LT216436`
  (Short 2020) and `KPNRH_*` on `FR997879`/`GCF_905477585.1` (Bruchmann 2021) — two assemblies of
  one strain, chromosomes 901 bp apart. Staged separately; never merge them.
- **The two DEG-derived columns key on `fasta_id`**, verbatim the header of
  `data/source/ncbi/deg_proteomes/<id>.faa`, so label and embedding share one namespace and there
  is no join at all. Measured against the E. coli anchor by exact sequence, EC958 shares only
  13.7% and O157:H7 23.0% — genuinely different proteomes, embedded in their own namespace.
  DEG1018 (Gerdes) is the opposite case at **99.4%**, so it maps onto the anchor and needs no
  embedding of its own.

---

## 2. Held back — downloaded, parseable, deliberately unused

These are **decisions, not gaps**. Each is on disk with a `SOURCE.md`.

| dataset | why held back |
|---|---|
| `ramage2017_kpnih1` | positives-only — 424 listed genes, no measured negatives |
| `paczosa2020_kppr1` | positives-only (310 hits) **and** condition-dependent (in-vivo, neutropenic) |
| `short2020_serum` | condition-dependent — serum resistance, 4 Kp strains |
| `rome2026_njst258` | condition-dependent — iron-depleted medium |
| `choe2023_m9` | condition-dependent — minimal medium; its extra 131 essentials are auxotrophies |
| `bruchmann_2hpi_6hpi` | condition-dependent — in-host timepoints from the same TraDIS sheets |
| `bachman2015/2023/2025_KPPR1` | condition-dependent — murine lung in-vivo fitness |
| `deg1019_keio` | duplicate — the Keio collection again under DEG's curation |
| CRISPRi: `cui2018`, `wang2018`, `rousset2018`, `rousset2021`, `hawkins2020`, `jana2023` | different quantity — CRISPRi measures **knockdown** fitness, with polar operon effects and guide-efficiency-dependent range. If ever used, they should form one endpoint per organism, not one per paper. |

The identifier work on Ramage is worth keeping even though the dataset is not: joining on GenBank
locus tags took it from v1's **212/424** by gene symbol to **424/424**.

---

## 3. Downloaded — wider corpora, not yet used

| corpus | size | why it is here |
|---|---|---|
| **OGEE v3** | 255,162 rows, 89 taxa, 124 datasets | **recovered from a dead server** — see §5. ~2× DEG's species and ~2.2× its positives, and it ships the negative class explicitly. **But *K. pneumoniae* is ZERO**, so it does not close the anchor's gap. |
| **DEG** | 49 sets, 173,048 rows, 38 species | the corpus behind `labels.py`/`deg_proteomes.py`; two of its sets are promoted to columns above |
| **Nichols 2011** | 3,980 × 324 continuous S-scores | the richest conditional E. coli resource. Its rows are Keio **deletion strains**, so by construction it covers only non-essential genes and cannot supply an essentiality label even in principle. |
| **Choe 2025** | 4,199 × 13 conditions | second-best conditional E. coli source |
| **Fitness Browser**, Feb 2024 | 221,030 sequences | conditional phenotypes only — see §5 for why it holds no essentiality signal |
| **iML1515 / iYL1228** | 2 metabolic models | FBA essentiality is a **prediction**; circular as a label, but legitimate as an independent comparator to score against |

---

## 4. Identified, not downloaded — needs a human

**`paczosa2020_kppr1` — RESOLVED, then held back for a different reason.** The files were
fetched by hand and live in `data/raw/kpneumoniae/essentiality/paczosa2020_kppr1/`. It is not used,
but the reason is no longer access: it is **positives-only** (310 published hits, no measured
negatives) **and** condition-dependent (in-vivo, neutropenic vs WT mouse). Two independent
disqualifications under the current filters.

The access finding is still worth keeping, because it is the hardest case met so far: Europe PMC
reports not-OA, curl 403s, **and an in-page `fetch()` from an already-loaded ASM tab also returns
403** — Cloudflare blocks XHR separately from navigation. No SRA/BioProject exists; the reads were
never deposited. A human navigating a browser directly to the URL is the only route.

**`mazzuoli2025_saureus`** (bioRxiv) — 4 strains incl. NCTC8325-4, 307-gene core essentialome.
bioRxiv exposes only 3 PDF figures and `aureobrowse.veeninglab.com` is a Shiny shell with no static
export. Re-check after journal publication, or email the lab.

---

## 5. Measured and refuted — do not re-derive

- **`ecl8_ess` (Eichelberger's own call) is unjoinable.** Their `ecl8_#####` tags annotate an
  assembly nobody deposited: **475 of 5,165** numeric suffixes shared with `BN373_`, **0 of 5,074**
  coordinates matched, gene symbols reached **3.7%** (the deposited ECL8 annotation carries symbols
  for only 144 CDS). Superseded by `bn373_ess` — same organism, same assay, 97.7% join.
- **OGEE v3's server is dead, but the file is NOT** — and CLAUDE.md recorded it as permanently
  lost on the strength of the server diagnosis alone. The diagnosis was right (`v3.ogee.info`
  aborts the TLS handshake server-side; verified against LibreSSL, OpenSSL and real Chrome, with
  and without SNI, pinned to TLS 1.2, by direct IP). The conclusion was wrong: **nobody checked
  the Internet Archive.** Two mirrors return byte-identical content, md5
  `b42f4a3358484b490e8bffe8c90edd00`, 1,151,931 bytes — the Wayback URL (**the `id_` suffix is
  mandatory**, or you get the HTML wrapper) and a GitHub mirror. *A dead server is not a lost
  dataset.*
- **`cain2017_njst258` does not exist.** No NJST258 TraDIS screen was ever published; NJST258_1/2
  are assemblies only. The real Cain 2017 paper is RH201207 (`srep42483`) and its sole supplement is
  a 6-page PDF with no per-gene table.
- **The compendium's `BW25113.out.DESeq.tsv` did not converge.** It looked like a free second
  opinion on Goodall's data under a uniform pipeline. It is not: `padj` is **NaN for 2,310 of
  4,256 genes (54.3%)** — against **7.0%** for ECL8 through the same pipeline — and of the 189
  genes with **zero insertion sites**, **108 are called `Unchanged`** (ECL8: 91 of 91 correctly
  `Reduced`). A gene with no insertions is the most essential thing a transposon library can show;
  DESeq cannot compute a statistic for it, so it falls through to non-essential. **The label is
  inverted for exactly the genes that matter most.** Caught by the ribosome control at 0.113
  against Goodall's 0.961 on the same data. Retired; the parser is kept in `screens.py` with the
  evidence. ECL8's column is unaffected and ships.
- **RB-TnSeq cannot see essential genes, by construction.** No insertions survive in an essential
  gene, so it is **absent** from the fitness table rather than carrying an extreme value.
  `feba.db` holds no essentiality signal. Two related corrections: the 2.3 GB `feba.db.gz` is the
  wrong file (the real one is the 2.2 MB `essential_proteins.tab`, a 1,000× difference), and
  **`Koxy` is *Klebsiella michiganensis* M5al, not *K. pneumoniae***.
- **BV-BRC "essentiality" is FBA prediction, not measurement** — circular as an ML label.
- **Jung 2019** — ENA project `PRJEB31265` exists but has released zero records.
- **MMseqs2 rewrites FASTA headers containing `|`.** It parses them as NCBI db-style fields, so
  DEG's `>lcl|HG941718.1_prot_CDN80371.1_1` comes back out of `_cluster.tsv` with the `lcl|`
  silently gone. That id is the key the screen table and the embedding share, so the clustering
  mapped to **0 of 4,981** proteins — indistinguishable from a proteome with no paralogs.
  `paralog_clusters.py` now writes index surrogates and substitutes the real ids back.
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
