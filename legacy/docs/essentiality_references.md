# Essentiality axis — references

Consolidated list of every published screen, computational predictor, metabolic model, tool,
and reference database used (or explicitly considered) in the **essentiality** axis of the GraDi
target-prioritization pipeline (docs §4, `scripts/07*`). References are split by organism —
**K. pneumoniae** and **E. coli** — followed by a **Shared methods & tools** section for the
organism-agnostic predictors, cross-species compendia, and infrastructure applied to both, and a
**watch-list** of documented-but-not-yet-ingested resources.

Notes:
- The K. pneumoniae anchor strain is **HS11286** (UniProt proteome `UP000007841`, 5,728 proteins);
  all Kp screens (ECL8, KPPR1, KPNIH1, MGH 78578) are mapped onto it by gene symbol or DIAMOND.
- **"Jana 2023" and "Zhu 2023" are the same paper** (*AEM*, `10.1128/aem.00956-23`); the code/docs
  use both names (also "Zhu/Wang 2023" in the curated-highlights file).
- Each entry notes the script / doc section that consumes it, and flags items that were
  **not obtained**, **deferred**, or are **watch-list** rather than used.

---

## 1. K. pneumoniae

### Direct experimental screens (§4.1)

- **Eichelberger / Short et al. 2024** — ECL8 (K2-ST375) TraDIS (Tn-seq); >554k unique insertions;
  in-vitro essential + urine + serum niche fitness. The deepest Kp essentiality call set.
  *eLife* 88971 · DOI [10.7554/eLife.88971](https://doi.org/10.7554/eLife.88971) · PMC11349299.
  Data via eLife CDN (`elife-88971-fig{1,4,6}-data1-v1.xlsx`). Used in `07a`/`07b`, §4.1a.

- **Jana / Zhu et al. 2023** — Kp Mobile-CRISPRi-seq; ~870 conditionally-essential genes (KPPR1S),
  in-vivo KPPR1 depletion screen (conditions: trimethoprim 25% MIC, polymyxin B, ceftriaxone,
  mouse lung). Supplement `s0001` also re-tabulates the Ramage 2017 KPNIH1 424-gene set and the
  Bachman 2015 KPPR1 in-vivo lists. *Appl Environ Microbiol* · DOI
  [10.1128/aem.00956-23](https://doi.org/10.1128/aem.00956-23) · PMC10617577. Fetched via an
  authenticated Chrome session (ASM 403s non-browser clients). Used in `07a`/`07l`/`07m`, §4.1c.

- **Ramage et al. 2017** — KPNIH1 / MKP103 (ST258, CR-Kp) 424-gene essential set (LB + chloramphenicol;
  Tn-seq + arrayed-library consensus). *J Bacteriol* · DOI
  [10.1128/JB.00352-17](https://doi.org/10.1128/JB.00352-17) · PMC5637181. Paper not open-access;
  the 424-gene set is recovered in full from the Jana 2023 `s0001` supplement. §4.1a.

- **Bachman et al. 2015** — KPPR1 (ATCC 43816) in-vivo lung Tn-seq (InSeq), C57BL/6 mouse pneumonia,
  24 h post-inoculation. *mBio* · DOI [10.1128/mBio.00775-15](https://doi.org/10.1128/mBio.00775-15) ·
  PMC4462621. Staged (insertion-level; not separately aggregated). §4.1b.

- **Mike & Bachman 2023** — KPPR1 tissue-specific in-vivo Tn-seq across blood, spleen, liver, lung
  (per-gene defect calls: log2FC < 0 & p < 0.05). Primary in-vivo source; also supplies the
  VK055_→gene-symbol bridge. *PLoS Pathog* · PMC10381055
  ([article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10381055/)). Used in `07a`/`07b`/`07l`, §4.1b.

- **Bachman et al. 2025** — KPPR1 bacteremic-dissemination Tn-seq. *Nat Commun* · DOI
  [10.1038/s41467-025-56095-3](https://doi.org/10.1038/s41467-025-56095-3) · PMC11742683
  (open access; fetched via Europe PMC). Staged, not aggregated this round. §4.1b.

- **Paczosa et al. 2020** — KPPR1 lung fitness, WT vs neutropenic hosts. *Infect Immun* ·
  PMID [31988174](https://pubmed.ncbi.nlm.nih.gov/31988174/). Listed resource; not in the fetch
  manifest. §4.1b.

- **Cain et al. 2017** — NJST258 (ST258) "secondary resistome" under colistin / imipenem /
  ciprofloxacin. *Sci Rep* srep42483 · PMC5309761
  ([article](https://www.nature.com/articles/srep42483)). **Not fetched** (compendium lacks
  NJST258; low priority). §4.1a.

### Metabolic model

- **Liao et al. 2011 (iYL1228)** — genome-scale metabolic reconstruction of *K. pneumoniae*
  MGH 78578 (1,229 genes, `KPN_` locus tags; MGH78578 UniProt proteome `UP000000265`).
  PMID [21478289](https://pubmed.ncbi.nlm.nih.gov/21478289/) · BiGG model
  [`iYL1228.json`](http://bigg.ucsd.edu/static/models/iYL1228.json). cobrapy single-gene deletion,
  KO/WT biomass ratio < 0.01 → in-silico essential; `KPN_` genes mapped onto HS11286 by DIAMOND.
  `07f`, §4.3d.

---

## 2. E. coli

### Direct experimental screens (`07n` / `07o`)

- **Baba et al. 2006 (Keio) via PEC** — arrayed single-gene KO essential set (`Class==1` in
  PECData.dat; 287 essential). Also the curated EcoGene 299-gene set (see below). *Mol Syst Biol* ·
  DOI [10.1038/msb4100050](https://doi.org/10.1038/msb4100050) · PMC1681482. Direct download:
  `shigen.nig.ac.jp/ecoli/pec/download/files/PECData.dat`. §4.2a.

- **Goodall et al. 2018** — BW25113 TraDIS (Keio ∩ PEC ∩ Goodall consensus; 353 essential).
  *mBio* 9:e02096-17 · DOI [10.1128/mBio.02096-17](https://doi.org/10.1128/mbio.02096-17) ·
  PMC5821084. Table S4. Fetched via Chrome (ASM 403). `07n`; transferred to Kp in `07c`.

- **Rousset et al. 2018** — genome-wide CRISPRi depletion (gene median log2FC; 374 essential at
  log2FC ≤ −2). *PLoS Genet* 14:e1007749 · PMC6242692. Europe PMC supp zip. `07n`.

- **Wang et al. 2018** — pooled CRISPRi, 5 phenotypes (gene fitness; 228 essential at fitness < −6).
  *Nat Commun* 9:2475 · PMC6018678. Europe PMC supp zip. `07n`.

- **Cui et al. 2018** — CRISPRi; dCas9 guide-toxicity ("bad-seed") filter reference. *Nat Commun*
  9:1912 · PMC5954155. Europe PMC supp zip. `07n`.

- **Rousset et al. 2021** — CRISPRi across 18 E. coli strains (fraction-essential; per-strain
  essential at score ≤ −3). *Nat Microbiol* 6:301. Springer CDN MOESM (browser-UA). `07n`.

- **Hawkins et al. 2020** — mismatch-CRISPRi vulnerability curves (min relative fitness; 316 curves).
  *Cell Systems* 11:523 · PMC7704046. Fetched via Chrome (Cell CDN). `07n`.

- **Price et al. 2018 / Wetmore et al. 2015** — RB-TnSeq / Fitness Browser (Keio E. coli),
  3,790 genes × ~3,500 conditions; per-gene min-fitness + the ~280-condition antibiotic/stress
  matrix (powers `07o`; deliberately excluded from the general vulnerability score). *Nature* 557:503;
  method Wetmore 2015. Direct: `morgannprice.org/FEBA/Keio/`. `07n`/`07o`.

- **EcoGene essential call (BW25113, 299 genes)** — curated binary essential set delivered via the
  Enterobacteriaceae-TraDIS compendium column `"EcoGene Essentiality: Escherichia coli BW25113"`;
  the strict §4.2a transfer set lifted onto Kp, and the ProteomeLM training labels. `07c`/`07d`.

- **Nichols et al. 2011** — chemical-genomics S-score matrix. *Cell* · PMC3060659.
  **Not obtained** (PMC copy behind an image reCAPTCHA; redundant with the RB-TnSeq 280-condition
  matrix). §4.

### Metabolic model

- **Monk et al. 2017 (iML1515)** — genome-scale reconstruction of *E. coli* K-12 MG1655
  (Blattner `b####` genes). BiGG model
  [`iML1515.json`](http://bigg.ucsd.edu/static/models/iML1515.json). cobrapy single-gene deletion
  (KO growth ratio < 0.01; 195 essential, metabolic subset). `07f`, §4.3d.

---

## 3. Shared methods, tools & cross-species resources

Organism-agnostic; applied to both E. coli and K. pneumoniae.

### Cross-species experimental compendium

- **Enterobacteriaceae-TraDIS compendium** (Goodall / Gardner group, Gardner-BinfLab) —
  ortholog-cluster-aligned essentiality matrix (`giant-tab_final.tsv`) across 12 genomes
  (2 *Klebsiella* incl. ECL8, 3 *E. coli* incl. BW25113 & NCTC13441 UPEC ST131, *Citrobacter*,
  *Salmonella*). Carries the curated EcoGene call + graded `Enterobacteriaceae %essential` /
  `Bacteria %essential`; the §4.2a/§4.2c backbone for broad-spectrum conservation-of-essentiality.
  PMID [39207104](https://pubmed.ncbi.nlm.nih.gov/39207104/) · GitHub
  [Gardner-BinfLab/Enterobacteriaceae-TraDIS](https://github.com/Gardner-BinfLab/Enterobacteriaceae-TraDIS)
  (raw base `raw.githubusercontent.com/Gardner-BinfLab/Enterobacteriaceae-TraDIS/master`).
  `07a`/`07c`/`07l`.

### Computational predictors

- **ProteomeLM** (Cuturello & Bitbol 2025, Bitbol-Lab; Apache-2.0) — proteome-scale transformer over
  the 01a ESM-C 600M embeddings. The upstream `-Ess` head is unreleased, so a logistic head is
  trained locally on the E. coli EcoGene labels (5-fold CV AUROC 0.809) then applied to both proteomes.
  Primary predictor (§4.3a). Code: [github.com/Bitbol-Lab/ProteomeLM](https://github.com/Bitbol-Lab/ProteomeLM) ·
  weights HuggingFace `Bitbol-Lab/ProteomeLM-M`. `07d`.

- **Geptop 2.0** (Wen et al. 2019) — orthology + phylogeny essentiality scoring; reimplemented with
  DIAMOND RBH over 37 DEG reference proteomes, cutoff 0.24.
  [cefg.uestc.cn/geptop](http://cefg.uestc.cn/geptop) · reference bundle `Geptop_v2.0.rar`
  (`datasets2/*.faa` + `DEG2` essential IDs, extracted with `unar`). §4.3b, `07e`.

- **DeeplyEssential** (Hasan & Lonardi 2020) — DNA + protein deep neural network trained on DEG.
  **Deferred** — NaN placeholder (no released weights; py2/TF1.6; license-less).
  [github.com/ucrbioinfo/DeeplyEssential](https://github.com/ucrbioinfo/DeeplyEssential). §4.3c, `07g`.

- **FBA / COBRApy** (`cobra >= 0.29`) — flux-balance single-gene deletion on the genome-scale
  metabolic models (iYL1228 Kp / iML1515 Ec) from BiGG (`bigg.ucsd.edu/static/models/`). §4.3d, `07f`.
  (Predictor consensus weights: ProteomeLM 0.5 · Geptop 0.3 · FBA 0.2.)

### Supporting tools

- **DIAMOND** — accelerated protein aligner; search engine for Geptop RBH and strain→HS11286 mapping.
  [github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond).
- **OrthoFinder** — orthogroup inference (the 07-series reuses the precomputed 03a orthology table).
  [github.com/davidemms/OrthoFinder](https://github.com/davidemms/OrthoFinder).
- **ESM-C 600M** (EvolutionaryScale, Cambrian) — per-protein embeddings feeding ProteomeLM.
  [evolutionaryscale.ai/blog/esm-cambrian](https://www.evolutionaryscale.ai/blog/esm-cambrian) ·
  Biohub ESM Atlas [biohub.ai/esm/protein/atlas](https://biohub.ai/esm/protein/atlas).

### Reference databases (checked; not directly ingested)

- **OGEE v3** (Gurumayum et al. 2021) — bacterial essential-gene DB; ProteomeLM's upstream training
  set (flagged for training-set-parroting concern). No working programmatic endpoint (TLS/SPA) →
  skipped for direct lookup. PMID [33084874](https://pubmed.ncbi.nlm.nih.gov/33084874/) ·
  [v3.ogee.info](https://v3.ogee.info) · [NAR paper](https://academic.oup.com/nar/article/49/D1/D998/5934414).
- **DEG 15** (Database of Essential Genes) — does **not** directly index *K. pneumoniae*; used only as
  Geptop's reference set (`DEG2`) and via E. coli. [tubic.org/deg](https://tubic.org/deg/).

### Data-fetch infrastructure (`07a`)

- Europe PMC supplementary-files REST — `ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles`
- NCBI Open Access service — `ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi`
- UniProt REST proteome stream — `rest.uniprot.org/uniprotkb/stream`
- BiGG model host — `bigg.ucsd.edu/static/models/{model}.json`
- Publisher CDNs (eLife / Springer / ASM Cell) — gated tables fetched via an authenticated Chrome
  session (chrome-devtools MCP same-origin fetch).

---

## 4. Watch-list / suggested (documented, not yet ingested)

Referenced in `docs/04_essentiality.md` as proposed enrichments, not part of the current results:

- **Bosch & Rock 2021** — graded CRISPRi vulnerability (*Mtb*); proposed refit of Jana 2023 CRISPRi
  titration into a per-gene Kp vulnerability score.
  [PMC8382161](https://pmc.ncbi.nlm.nih.gov/articles/PMC8382161/).
- **Lin et al. 2025** — hvKp Tn-seq in *Galleria*; to break the KPPR1 monoculture in §4.1b.
  [Front. Cell. Infect. Microbiol.](https://www.frontiersin.org/journals/cellular-and-infection-microbiology/articles/10.3389/fcimb.2025.1643224/full).
- **Insua et al. 2021** — MDR cKp in-vivo. [PMID 33512418](https://pubmed.ncbi.nlm.nih.gov/33512418/).
- **Liu / van Opijnen 2024** — synthetic lethality; informs a proposed Kp synthetic-lethality flag (§4.5).
  [Nat. Microbiol.](https://www.nature.com/articles/s41564-024-01759-x).
- **Christen et al. 2025** — InducTn-seq; watch-item, no Kp dataset yet.
  [PMID 40148565](https://pubmed.ncbi.nlm.nih.gov/40148565/).
- Cross-species vulnerability transfer refs (non-Kp, §4.2c enrichment): **Geisinger/Wang 2023**
  (*A. baumannii*), **Poulsen 2019** (*P. aeruginosa*), **Bosch 2021** (*Mtb*).
