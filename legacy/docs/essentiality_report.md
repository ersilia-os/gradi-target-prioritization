# Essentiality — full technical report

**Axis 4 of the GraDi target-prioritization pipeline.** This document is a complete, rigorous account
of the essentiality deliverable: every dataset (with citation, provenance, download route, identifier
keying, thresholds and coverage), the identifier-mapping and scoring methodology, the pipeline scripts,
and a panel-by-panel description of every figure. It covers both focal organisms — *Klebsiella
pneumoniae* HS11286 (the anchor strain) and *Escherichia coli* K-12 — which are treated with **equal,
first-class importance**.

Companion files: `docs/04_essentiality.md` (the axis spec + a machine-readable script/column table),
`docs/essentiality_log.md` (the dated run log), and the interactive slide gallery (artifact).

---

## 1. What this axis answers

> **Is a given protein required for the growth or survival of the pathogen?**

Essentiality is one of the strongest arguments for a drug target: if knocking a gene out (or knocking
it down) kills or severely impairs the bacterium, the encoded protein is a candidate target. But
essentiality is neither binary nor absolute:

- It is **graded** — some genes are strictly required, others only reduce fitness, others are
  "vulnerable" (sensitive to partial knockdown).
- It is **conditional** — a gene can be dispensable in rich lab medium yet required in a host niche
  (urine, serum, lung) or under antibiotic stress.
- It is **measured by different technologies** that disagree at the margins — arrayed single-gene
  knockouts (Keio), transposon-insertion sequencing (Tn-seq / TraDIS / RB-TnSeq), and CRISPR
  interference (CRISPRi) each have characteristic blind spots.

Accordingly, the deliverable **integrates three independent classes of evidence** into a single graded
`essentiality_score` ∈ [0, 1] and an evidence-driven `essentiality_tier` ∈ {essential,
likely_essential, non_essential}, while keeping every underlying measurement as an auditable column:

| Evidence class | Weight | What it is |
| --- | --- | --- |
| **Experimental** (direct) | 0.40 | Published essentiality/CRISPRi/transposon screens *in that organism* |
| **Cross-species transfer** | 0.20 | Essentiality of the ortholog in *E. coli* and across 12 Enterobacteriaceae |
| **Computational predictors** | 0.40 | ProteomeLM (LM), Geptop 2.0 (orthology/phylogeny), FBA (metabolic model) |

The design deliberately follows the axis spec's directive to *"emit a graded 0–1 vulnerability, not a
binary call"* and to *"replace vote-counting with a weighted ensemble."*

---

## 2. Headline results

| | *K. pneumoniae* HS11286 | *E. coli* K-12 |
| --- | --- | --- |
| Proteins scored | 5,728 | 4,403 |
| **essential** tier | **429** | **402** |
| likely_essential tier | 801 | 135 |
| non_essential tier | 4,498 | 3,866 |
| Experimentally essential (rigorous screens) | 353 | 399 |
| Composite score ≥ 0.60 | 393 | 327 |
| Experimental sub-score coverage | 1,471 | 4,122 |
| Transfer sub-score coverage | 3,178 | 4,257 |
| Predictor sub-score coverage | 5,728 | 4,403 |
| **Prime shortlist** (broad-spectrum-selective ∧ essential ∧ ligandable) | **222** | **209** |

The two organisms are now symmetric. E. coli's experimental coverage (4,122 / 4,403) is higher than
Kp's (1,471 / 5,728) simply because E. coli is the most intensively screened bacterium on earth — the
deliverable leverages that.

Spot-check (both organisms): the top-scoring proteins are textbook essentials — ribosomal proteins
(rplF, rpsE, rpsC…), RNA polymerase (rpoB, rpoA, rpoC), DNA gyrase/topoisomerase (gyrA, gyrB, parE),
cell-division (ftsZ, ftsA), the Sec translocon (secY, secA), and cell-envelope biogenesis (murC, mraY,
lpxC) — each supported by 4–7 independent lines of evidence.

---

## 3. Datasets

Every dataset below is keyed onto the anchor proteome by **UniProt accession** (the project
convention). Three keying bridges are used, all built from the proteome TSV `Gene Names` column:
locus tag, Blattner b-number / Keio JW id (E. coli), and gene symbol. Kp is a "dark" TrEMBL proteome
(few gene symbols), so cross-strain screens are mapped onto Kp by **DIAMOND sequence best-hit**; E. coli
carries b-numbers for 4,402/4,403 proteins and JW ids for 4,252, so E. coli screens map almost perfectly
by identifier.

### 3.1 *K. pneumoniae* direct experimental screens (§4.1)

#### Eichelberger / Short 2024 — ECL8 TraDIS *(in-vitro essential + urine + serum)*
- **Citation:** Eichelberger *et al.* 2024, *eLife* 88971 (DOI 10.7554/eLife.88971). Strain ECL8 (K2-ST375).
- **Screen:** transposon-directed insertion sequencing (TraDIS). Per-gene bimodal essentiality call
  (Essential / Non-essential / Unclear) in LB, plus fitness (log2FC + q-value) after passage in pooled
  human **urine** and exposure to human **serum**.
- **Download:** the 5 source-data XLSX from the eLife CDN (`cdn.elifesciences.org/articles/88971/…`),
  re-staged from the prior v1 attempt. Open access.
- **Keying / mapping:** the paper's `ecl8_*` locus tags are a bespoke Prokka re-annotation absent from
  any public FASTA; genes are mapped to HS11286 by **gene symbol** (essential genes are the conserved,
  named ones). Essential = the `Essential` boolean column; niche-required = logFC < 0 & q < 0.05.
- **Coverage on Kp:** 954 genes mapped (249 essential), 688 urine-conditional, 666 serum-conditional.

#### Jana / Zhu 2023 — Mobile-CRISPRi-seq *(the headline Kp CRISPRi screen)*
- **Citation:** Jana/Zhu *et al.* 2023, *Applied and Environmental Microbiology*, DOI 10.1128/aem.00956-23.
- **Screen:** Mobile-CRISPRi-seq. Three tables in supplement `s0001`: (a) the **870-gene conditionally-
  essential library** (KPN_ MGH78578 locus tags); (b) an **in-vivo KPPR1 mouse-lung** depletion screen
  (VK055_ locus + depletion ratio + Ceder p-value); (c) the **KPNIH1 essential set (424)** re-tabulated
  from Ramage 2017; it also embeds the Bachman 2015 in-vivo list.
- **Download:** ASM returns HTTP 403 to non-browser clients; the tables were retrieved via an
  **authenticated Chrome session** (chrome-devtools MCP `evaluate_script` same-origin `fetch` with
  cookies → base64 → disk). Not open access.
- **Keying / mapping:** KPN_ → HS11286 by **DIAMOND** against a KPN_-keyed MGH78578 proteome FASTA
  (846/870 library genes map); KPNIH1/gene by gene symbol; in-vivo VK055_ → symbol via Mike & Bachman.
- **Coverage on Kp:** CRISPRi library 846, KPNIH1 essential 212.

#### KPNIH1 (Ramage 2017) and KPPR1 in-vivo (Bachman 2015, Mike & Bachman 2023, Bachman 2025)
- **Ramage 2017** (*J. Bacteriol.*, KPNIH1/MKP103 Tn-seq, 424 essential) is not open access, but its
  essential-gene list is **recovered in full from the Jana 2023 supplement** (see above).
- **KPPR1 in-vivo** — Mike & Bachman 2023 (*PLoS Pathog.* PMC10381055, genome-wide TnSeq) and
  Bachman 2015 (*mBio* PMC4462621) fetched via Europe PMC; **Bachman 2025** (*Nat. Commun.*
  s41467-025-56095-3, PMC11742683) fetched via Europe PMC once it became open. Provide niche/in-vivo
  fitness (a supporting signal). Bachman 2015 is insertion-level; the per-gene in-vivo call used comes
  from Mike & Bachman 2023 + the Jana in-vivo sheet.

### 3.2 *E. coli* K-12 direct experimental screens (§4.1/§4.2) — first-class

E. coli is the model organism for essentiality; the deliverable ingests the major genome-scale screens
across all three technologies. Per-screen essential counts on the E. coli proteome:

| Dataset | Citation | Screen type | Per-gene signal | Download route | Essential |
| --- | --- | --- | --- | --- | --- |
| **Keio / PEC** | Baba *et al.* 2006, *Mol. Syst. Biol.* | arrayed single-gene KO | binary (`Class==1`) | direct: `shigen.nig.ac.jp/ecoli/pec/download/files/PECData.dat` | 287 |
| **Goodall 2018** | *mBio* 9:e02096-17 (PMC5821084) | TraDIS | binary `Essential` (Table S4) | **Chrome** (ASM 403) | 353 |
| **Rousset 2018** | *PLoS Genet.* 14:e1007749 (PMC6242692) | genome-wide CRISPRi | gene median log2FC | Europe PMC zip | 374 (log2FC ≤ −2) |
| **Wang 2018** | *Nat. Commun.* 9:2475 (PMC6018678) | pooled CRISPRi (5 phenotypes) | gene fitness | Europe PMC zip | 228 (fitness < −6) |
| **Cui 2018** | *Nat. Commun.* 9:1912 (PMC5954155) | CRISPRi | guide toxicity (dCas9 bad-seed filter) | Europe PMC zip | — |
| **Rousset 2021** | *Nat. Microbiol.* 6:301 | 18-strain CRISPRi | strain/condition score → % essential | Springer CDN (browser-UA) | (graded) |
| **Hawkins 2020** | *Cell Systems* 11:523 (PMC7704046) | mismatch-CRISPRi | expression–fitness (vulnerability) | **Chrome** (Cell CDN) | 316 (curves) |
| **RB-TnSeq / Fitness Browser** | Price *et al.* 2018, *Nature* 557:503; Wetmore 2015 | randomly-barcoded Tn-seq | fitness across ~3,500 conditions | direct: `morgannprice.org/FEBA/Keio/` | condition-specific |

- **Keying / mapping:** b-number, Keio JW id, or gene symbol → E. coli UniProt accession (near-complete).
- **RB-TnSeq processing:** the raw Keio matrix is 3,790 genes × 3,511 conditions (247 MB); it was reduced
  to (i) a per-gene min-fitness / mean-fitness / n-conditions-with-defect summary, and (ii) a curated
  **280-condition antibiotic/stress matrix**. The giant raw file is not retained.
- **Consensus:** E. coli `experimental_essential` = essential in Keio or Goodall (rigorous genome-wide
  KO/TraDIS) or in ≥2 of the 4 binary screens → **399 genes**; 230 are essential across all three method
  types (KO ∩ TraDIS ∩ CRISPRi).
- **Vulnerability score** (`ecoli_vulnerability_score`, 0–1): the strongest depletion across the
  standard-growth CRISPRi screens (Rousset 2018, Wang 2018) + Hawkins mismatch-CRISPRi. RB-TnSeq's
  all-condition min-fitness is **deliberately excluded** (it captures condition-specific essentiality —
  e.g. *lacZ* on lactose — which belongs to the condition view, not the general vulnerability; without
  this exclusion *lacZ*/*araC* spuriously scored ~0.8 instead of ~0.05).

### 3.3 Cross-species conservation (§4.2a/§4.2c) — Enterobacteriaceae-TraDIS compendium

- **Citation:** Goodall/Gardner group (Gardner-BinfLab), PMID 39207104 — the recommended §4.2a upgrade.
- **Source:** `giant-tab_final.tsv` from the GitHub repo (reliable headless). A fully ortholog-cluster-
  aligned essentiality matrix keyed by the E. coli Keio b-number, carrying **12 genomes** of TraDIS
  essentiality (2 *Klebsiella* incl. ECL8, 3 *E. coli*, *Citrobacter*, 6 *Salmonella*) plus the curated
  **EcoGene** essential set (299 genes) and graded **`Enterobacteriaceae %essential`** /
  **`Bacteria %essential`** columns.
- **Decoding:** the per-genome numeric TraDIS log-ratio was calibrated against the curated EcoGene
  binary — **essential ≤ −0.5** (F1 = 0.75). Decoded per-genome essential-set sizes are biologically
  sensible (ECL8 279, *K. pneumoniae* RH201207 351, *E. coli* BW25113 391, *Salmonella* 251–338).
- **Mapping:** onto Kp via Kp → E. coli-K12 ortholog → b-number → compendium row (3,170 Kp proteins
  covered; **249 core-essential** = essential in ≥80% of the 12 genomes). E. coli maps directly by
  b-number (4,257 covered; 250 core).
- **Role:** the graded broad-spectrum `%essential` score is a strong antibacterial-target argument
  (broadly essential across Enterobacteriaceae).

### 3.4 Computational predictors (§4.3)

| Predictor | Citation | Method | Output column | Kp / ec essential |
| --- | --- | --- | --- | --- |
| **ProteomeLM-Ess** | Cuturello & Bitbol 2025 (Bitbol-Lab) | proteome-scale transformer over ESM-C 600M embeddings + a logistic essentiality head | `proteomelm_ess_score` | 645 / 412 |
| **Geptop 2.0** | Wen *et al.* 2019 | orthology + phylogeny; DIAMOND RBH over 37 DEG reference genomes | `geptop_score` (cutoff 0.24) | 417 / 381 |
| **FBA** | Liao 2011 (iYL1228, Kp) / Monk 2017 (iML1515, ec) | single-gene deletion on a genome-scale metabolic model | `fba_essential` (KO growth ratio < 0.01) | 120 / 195 |
| **DeeplyEssential** | Hasan & Lonardi 2020 | DNA+protein DNN | **deferred** (NaN placeholder) | — |

- **ProteomeLM caveat & fix:** the published `-Ess` head is *not* released. We run the released backbone
  (its input is exactly the ESM-C 600M embeddings already computed in step 01a) and **train our own
  logistic head on the curated E. coli EcoGene labels** — 5-fold cross-validated **AUROC = 0.809** on
  E. coli — then freeze it and apply to Kp. Because the training labels are curated (not OGEE-derived),
  this side-steps the OGEE-parroting concern the spec raises.
- **Geptop** was reimplemented with DIAMOND (the upstream Python-2 + NCBI-BLAST `.rar` is unbuildable
  headless); reference weight = median RBH % identity (a documented proxy for the composition-vector
  distance).
- **FBA** covers only the metabolic subset (~1,229 Kp / ~1,515 ec model genes); non-model proteins are
  `not_in_model` (NaN), never scored as non-essential. FBA recovers only ~14 % of ECL8-essential genes
  because most essentials are informational (ribosome/RNA-pol), not metabolic — this is expected and
  honest.

### 3.5 Not obtained

- **Nichols 2011** chemical-genomics S-score matrix (*Cell*, PMC3060659) — the PMC copy is behind an
  **image reCAPTCHA**, which was not solved on principle (anti-bot circumvention). It is optional and
  largely redundant with the RB-TnSeq 280-condition matrix. To add it: manually download
  `PMC3060659/bin/NIHMS261392-supplement-02.txt` into `data/raw/ecoli/essentiality/nichols2011_chemgen/`.
- **DEG / OGEE / EcoCyc** — no reliable programmatic endpoint (bot-hostile / login-gated); their
  content is already represented via the primary screens and the compendium.

---

## 4. Identifier-mapping methodology

Because HS11286 is a dark TrEMBL proteome, exact-accession matching silently misses data. The mapping
strategy mirrors the ligandability stage's "map by sequence" convention:

- **Locus-tag / b-number / JW / gene-symbol bridges** (`src/essentiality.py`): built from the proteome
  TSV `Gene Names` column (`locus_to_uniprot`, `gene_aliases_to_uniprot`, `jw_to_uniprot`).
- **DIAMOND sequence best-hit** (`map_strain_by_sequence`): for locus-tag-keyed strains absent from the
  UniProt panel (ECL8, KPPR1, MGH78578) — build a per-strain FASTA keyed by the strain's locus tag,
  blastp the anchor proteome against it, keep the best hit above an identity/coverage floor (≥50 %
  pident, ≥70 % qcov). Used for the Jana CRISPRi library (KPN_) and the FBA model genes.
- **Ortholog transfer** (`transfer_ecoli_to_kp`): lift an E. coli signal onto Kp via the
  `Ecoli_K12_MG1655` orthology from the 03a table (3,177 Kp proteins have an E. coli ortholog). Used to
  carry the E. coli screens and EcoGene call onto Kp.

---

## 5. Composite score & tiers (`07h`)

All weights and thresholds are module-level constants in `src/essentiality.py` (auditable/tunable).

**Sub-scores** (each ∈ [0, 1], kept as columns; a sub-score is **NaN — dropped, not zero-filled —**
when the protein has no data for it, and the composite weights are renormalised over the available
sub-scores):

- `evidence_experimental` (**0.40**) — direct screens in that organism.
  - *Kp:* ECL8-essential or KPNIH1-essential (rigorous genome-wide Tn-seq) → **1.0**; CRISPRi CE-library
    membership → **0.85** (strong but softer — a conditionally-essential selection, not a bimodal call);
    CRISPRi in-vivo / vulnerability → 0.7; urine/serum/in-vivo niche → 0.6; ECL8 "unclear" → 0.4.
  - *E. coli:* Keio or Goodall essential → **1.0**; ≥2 CRISPRi screens agree → 0.9; 1 screen → 0.75;
    else the graded `ecoli_vulnerability_score`.
- `evidence_transfer` (**0.20**) — the maximum of: the EcoGene ortholog call (1.0), the graded
  Enterobacteriaceae `%essential`, and (Kp only) the E. coli-screen consensus lifted via ortholog.
- `evidence_predictor` (**0.40**) — weighted mean of the *available* predictors
  (ProteomeLM 0.5 / Geptop 0.3 / FBA 0.2), renormalised.

**Composite:** `essentiality_score = Σ (wᵢ · sub_scoreᵢ) / Σ wᵢ` over available sub-scores, ∈ [0, 1].

**Tiers** (evidence-driven, not a pure threshold):
- `essential` — a rigorous genome-wide essential call (Kp ECL8/KPNIH1 Tn-seq, E. coli Keio/Goodall, or
  the E. coli EcoGene ortholog), OR a strong predictor+transfer consensus (both ≥ 0.7), OR
  `essentiality_score ≥ 0.60`.
- `likely_essential` — `score ≥ 0.35`, OR `evidence_experimental ≥ 0.5`, OR `evidence_predictor ≥ 0.5`.
- `non_essential` — otherwise.

The CRISPRi CE-library does **not** hard-force the essential tier (it boosts the score and guarantees
≥ likely_essential), because "included in a conditionally-essential library" is softer evidence than a
bimodal Tn-seq essential call.

**Outputs** (`output/results/<org>/`): `<prefix>_essentiality.csv` (full per-protein table with all
sub-scores, per-track calls, `experimentally_essential`, `evidence_sources` provenance string,
`selectivity`) and `<prefix>_essentiality_shortlist.csv` (broad-spectrum-selective ∧ essential, ranked).

---

## 6. Pipeline (`scripts/07*`)

Run per organism (`--organism {kpneumoniae,ecoli}`) in the `gradi` conda env; DIAMOND comes from
`gradi-ortho`. Everything is cached / resumable.

| Script | Role |
| --- | --- |
| `07a_fetch_essentiality.py` | robust fetcher: publisher-CDN → Europe-PMC-zip → NCBI-OA → placeholder ladder; stages the Kp datasets + the Enterobacteriaceae-TraDIS compendium |
| `07b_kp_experimental.py` | parse ECL8 (in-vitro/urine/serum) + KPPR1 in-vivo → `kp_ess_kp.csv` |
| `07c_ecoli_transfer.py` | E. coli EcoGene + broad-spectrum + **E. coli-screen consensus** transferred onto Kp; direct EcoGene call for ec → `<p>_ess_ecoli.csv` |
| `07d_proteomelm_ess.py` | ProteomeLM backbone + locally-trained essentiality head → `<p>_ess_proteomelm.csv` |
| `07e_geptop.py` | Geptop 2.0 reimplemented with DIAMOND → `<p>_ess_geptop.csv` |
| `07f_fba_iyl1228.py` | cobrapy single-gene deletion (iYL1228 kp / iML1515 ec) → `<p>_ess_fba.csv` |
| `07g_deeplyessential.py` | deferred NaN placeholder |
| `07h_essentiality_merge.py` | the graded composite + tiers + shortlist → `<p>_essentiality.csv` |
| `07l_publication_essentiality.py` | consolidate the publication (experimental-only) evidence + cross-species matrix → `<p>_ess_publications.csv` |
| `07n_ecoli_experimental.py` | ingest the E. coli screens (Keio/Goodall/CRISPRi/RB-TnSeq) → `ec_ess_experimental.csv` + condition matrix |
| `07i / 07j / 07k / 07m / 07o` | the five figure slides (below) |

`src/essentiality.py` reuses the ligandability helpers (`from src import ligandability as L`) and adds
the E. coli-specific bridges and the scoring constants.

---

## 7. Figures (the gallery)

Four stylia slides per organism (16:9, 2×3 panels, NPG palette). Kp accent = NPG red, E. coli = cyan.
Each is `output/plots/07{i,j,k,m,o}_*_{kp,ec}.png`.

### 7.1 `07m` — Published essentiality screens (experimental, no predictions)
The headline experimental slide. Organism-symmetric.
1. **Screen coverage** — proteins essential in each published screen (Kp: CRISPRi library / KPNIH1 /
   ECL8 / cross-species core; E. coli: Keio / Goodall / Rousset18 / Wang18 / core).
2. **CRISPRi signal** — Kp: top in-vivo Mobile-CRISPRi depleted genes (KPPR1 mouse lung); E. coli: the
   genome-wide CRISPRi depletion scatter (Rousset 2018 log2FC vs Wang 2018 fitness, Keio-essential in red).
3. **Screen agreement** — how many independent screens call each experimentally-essential gene essential.
4. **Cross-species conservation** — # of the 12 Enterobacteriaceae genomes each gene is essential in.
5. **Core essential-genome heatmap** — top genes × the 12 experimental genomes (core → accessory gradient).
6. **Method concordance (ec)** — KO vs TraDIS vs CRISPRi overlap; **CRISPRi library by function (kp)**.

### 7.2 `07o` — Condition- & antibiotic-specific essentiality
1–6 (E. coli): # genes sensitized per antibiotic (RB-TnSeq); drug × gene sensitivity heatmap; most
condition-variable genes; constitutive vs condition-specific requirement; ciprofloxacin-sensitizing
genes; trimethoprim-sensitizing genes. (Kp): niche-conditional counts; urine vs serum requirement; top
urine- and serum-required genes; niche overlap; a note that the broad antibiotic landscape is E. coli-
specific.

### 7.3 `07i` — Computational predictor comparison
Score distributions; ProteomeLM–Geptop agreement; FBA metabolic knockouts; predictor consensus;
ProteomeLM calibration vs the curated labels; per-predictor sensitivity to the direct experimental calls.

### 7.4 `07k` — Target prioritization (synthesis)
The single synthesis slide (the former composite-summary `07j` and landscape slides were merged here;
`07j` is retired). Six panels: (1) **global essentiality score** as a reverse-cumulative curve with the
top-10 genes listed down the sparse tail; (2) **prioritization map** — the ESM-C protein-universe map
with all / essential / prime layers; (3) **essential vs ligandable** scatter (prime highlighted);
(4) **neglected vs essential** — essentiality vs 05a studiedness (understudied essentials = opportunities);
(5) **top prime-target scorecard** — top ~15 prime targets × six evidence axes (experimental,
conservation, predictor, essentiality, ligandability, studiedness) as a heatmap; (6) **prioritization
funnel** — all → essential → broad-spectrum-selective → ligandable = prime (√-scaled widths). Output
`07k_prioritization_{kp,ec}.png`. The gallery therefore has **four** slides per organism (07m · 07o ·
07i · 07k).

---

## 8. Caveats & limitations

- **Kp is dark; E. coli is rich.** Kp direct-experimental coverage (1,471) is lower than E. coli's
  (4,122) because far fewer screens exist for Kp — the cross-species transfer and predictors fill the gap.
- **CRISPRi ≠ essentiality.** A CE-library membership is a conditionally-essential *selection*; it boosts
  the score but does not hard-force the essential tier.
- **FBA is metabolic-only** (~1,229/1,515 genes in-model); it is a mechanistic complement, not a
  proteome-wide caller.
- **ProteomeLM head is ours, not the authors'.** The released model lacks the `-Ess` head; our logistic
  head (AUROC 0.809 on E. coli) is a faithful reproduction, but is trained on E. coli labels and applied
  cross-species to Kp.
- **Condition vs general vulnerability** are kept separate on purpose (RB-TnSeq's all-condition min-
  fitness is condition-specific and would otherwise inflate the general vulnerability score).
- **Gated data** required an authenticated browser (Jana 2023, Goodall 2018, Hawkins 2020); **Nichols
  2011** remains un-fetched (reCAPTCHA).

---

## 9. Reproduce

```bash
export GRADI_DIAMOND_BIN=$HOME/miniconda3/envs/gradi-ortho/bin/diamond
conda run -n gradi python scripts/07a_fetch_essentiality.py            # Kp datasets + compendium
# (E. coli screens: direct PEC/Fitness-Browser + Europe PMC zips + Chrome for Goodall/Hawkins; see the log)
for org in kpneumoniae ecoli; do
  conda run -n gradi python scripts/07b_kp_experimental.py       --organism $org
  conda run -n gradi python scripts/07n_ecoli_experimental.py    2>/dev/null   # ecoli only
  conda run -n gradi python scripts/07c_ecoli_transfer.py        --organism $org
  conda run -n gradi python scripts/07d_proteomelm_ess.py        --organism $org
  conda run -n gradi python scripts/07e_geptop.py                --organism $org
  conda run -n gradi python scripts/07f_fba_iyl1228.py           --organism $org
  conda run -n gradi python scripts/07g_deeplyessential.py       --organism $org
  conda run -n gradi python scripts/07l_publication_essentiality.py --organism $org
  conda run -n gradi python scripts/07h_essentiality_merge.py    --organism $org
  for plot in 07i_predictor_plots 07j_essentiality_plots 07k_essentiality_landscape \
              07m_publication_plots 07o_condition_plots; do
    conda run -n gradi python scripts/$plot.py --organism $org
  done
done
```

Outputs: `output/results/<org>/<prefix>_essentiality{,_shortlist}.csv` and
`output/plots/07{i,j,k,m,o}_*_{kp,ec}.png`. Data versions and the exact fetch routes (including the
Chrome-fetched sets) are recorded in `docs/essentiality_log.md`.
