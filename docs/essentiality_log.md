# Essentiality assessment log (docs §4)

Implemented `scripts/07a–07k` + `src/essentiality.py`. Per-protein essentiality for
K. pneumoniae HS11286 (5,728) and E. coli K-12 (4,403), keyed by UniProt accession. Emits a graded
[0,1] `essentiality_score` + `essentiality_tier` (essential / likely_essential / non_essential),
per the spec's "graded vulnerability, not a binary call" directive.

## Pipeline
`07a` fetch → `07b` Kp experimental · `07c` E. coli transfer → `07d` ProteomeLM · `07e` Geptop ·
`07f` FBA · `07g` DeeplyEssential(deferred) → `07h` merge → `07i` (predictors) · `07k` (6-panel target-
prioritization synthesis) · `07m` (publications) · `07o` (conditions) plots — 4 slides/organism; the old
`07j` composite-summary slide is retired (folded into `07k`). `src/essentiality.py`
reuses the ligandability helpers (`from src import ligandability as L`) and adds the locus-tag /
gene-symbol bridges, the DIAMOND-by-sequence strain mapper, and the scoring constants.

**Publication (experimental-only) sub-deliverable:** `07l_publication_essentiality.py` consolidates
*only measured* essentiality from published screens into `<prefix>_ess_publications.csv`, and
`07m_publication_plots.py` renders a dedicated prediction-free slide (CRISPRi library + in-vivo, KPNIH1/
ECL8 Tn-seq, cross-species core-essentialome heatmap, CRISPRi-library-by-function). This is the
headline for the "essentiality from publications" view; the Kp Mobile-CRISPRi-seq screen is central.

## Data versions / sources
- **Enterobacteriaceae-TraDIS compendium** (Goodall/Gardner, PMID 39207104), GitHub `Gardner-BinfLab`
  (fetched 2026-07-13) — the §4.2a/§4.2c backbone. `giant-tab_final.tsv` gives, per E. coli Keio
  b-number: the curated **EcoGene essential call (299 genes)** + graded **Enterobacteriaceae/Bacteria
  %essential** across E. coli, Klebsiella (incl. ECL8), Salmonella, Citrobacter, Enterobacter.
- **Eichelberger/Short 2024 ECL8** (eLife 88971) — §4.1a in-vitro essential (bimodal call) + urine/
  serum niche fitness. Re-staged from `data/raw/legacy/`.
- **Mike & Bachman 2023 KPPR1** (PLoS Pathog, PMC10381055, Europe PMC supp zip) — §4.1b in-vivo TnSeq.
- **Bachman 2015 KPPR1** (mBio, PMC4462621, Europe PMC supp zip) — §4.1b (insertion-level; not yet
  parsed into per-gene calls this round).
- **Jana/Zhu 2023 Mobile-CRISPRi-seq** (AEM aem.00956-23) — **full supplementary tables obtained**
  (2026-07-13) via an authenticated Chrome session (ASM 403s non-browser clients; `evaluate_script`
  fetch with session cookies). `s0001` = the 870-gene conditionally-essential CRISPRi library (KPN_),
  the in-vivo KPPR1 depletion screen (ratio + Ceder p), and — re-tabulated inside it — the **Ramage 2017
  KPNIH1 essential set (424)** and **Bachman 2015 KPPR1 in-vivo** gene lists. Staged under
  `data/raw/kpneumoniae/essentiality/jana2023_crispri/`.
- **Bachman 2025** (Nat Commun s41467-025-56095-3) — now open (PMC11742683); fetched via Europe PMC.
- **ProteomeLM-M** backbone (`Bitbol-Lab/ProteomeLM-M`, Apache-2.0) over the 01a ESM-C 600M embeddings.
- **Geptop 2.0** reference data (`Geptop_v2.0.rar`: 37 DEG reference proteomes + `DEG2` essential ids).
- **iYL1228** (Kp MGH78578) & **iML1515** (E. coli) BiGG metabolic models; `cobra` FBA.

## Key design decisions
- **Compendium over the gated Goodall Table S4.** Goodall 2018's ASM/PMC supp is unfetchable headless
  (ASM 403, no Europe PMC supp, OA tarball 404). Its data lives — cleaner and richer — in the
  Enterobacteriaceae-TraDIS compendium (GitHub raw), which *is* the Goodall/Gardner group's own
  release and upgrades §4.2a to a graded cross-species consensus. Fetcher ladder: publisher CDN →
  Europe PMC supp zip → NCBI OA tarball → placeholder (never blocks compute).
- **ProteomeLM `-Ess` head is unreleased**, so we train our own logistic head on the ProteomeLM
  contextualised embeddings using the curated E. coli EcoGene labels (from 07c), then freeze and apply
  to Kp. **5-fold CV AUROC = 0.809** on E. coli. Because the labels are curated (not OGEE-derived),
  this side-steps the OGEE-parroting concern. The expensive input (ESM-C 600M) already existed from
  01a, so the forward is ~2–3 s per proteome.
- **Geptop reimplemented with DIAMOND** (the py2/NCBI-BLAST `.rar` is unbuildable headless): best-hit
  ortholog per query in each of the 37 references, DEG-essential membership, weighted by the reference's
  median RBH %identity (a data-derived proxy for Geptop's composition-vector distance — documented
  simplification), min-max normalised, cutoff 0.24.
- **Cross-strain mapping by sequence.** ECL8/KPPR1/MGH78578 are locus-tag-keyed and absent from the
  UniProt orthology panel; FBA's `KPN_` genes are mapped onto HS11286 by DIAMOND (build a `KPN_`-keyed
  MGH78578 FASTA, blastp HS11286 against it). Kp screens (ECL8, KPPR1) are mapped by **gene symbol**
  (normalised, paralog suffix stripped) — adequate because essential genes are the named, conserved ones.
- **Graded composite (07h), missing tracks renormalised (not zero-filled).**
  `essentiality_score` = renormalised weighted sum over *available* sub-scores
  `0.40·experimental + 0.20·transfer + 0.40·predictor`, where `predictor` = weighted mean of available
  predictors (ProteomeLM 0.5 / Geptop 0.3 / FBA 0.2). Tier is evidence-driven: a direct Kp essential
  call, or strong predictor+transfer consensus, or score ≥ 0.60 → `essential`.

## Coverage (kp / ec)
- Tiers — kp: **383 essential / 368 likely / 4,977 non**; ec: 315 / 65 / 4,023.
- Experimentally essential (Kp ECL8 or E. coli EcoGene): 345 kp / 299 ec.
- Predictors — ProteomeLM essential: 645 kp / 412 ec; Geptop (≥0.24): 417 / 381; FBA-essential:
  120 / 195 (metabolic subset only: 1,237 kp / 1,515 ec proteins are in the model).
- Sub-score coverage kp: experimental 956, E. coli transfer 3,177, predictor 5,728.
- Predictor↔experiment: of 249 mapped Kp-ECL8-essential genes, ProteomeLM recovers 81%, Geptop 87%,
  FBA 14% (FBA covers only metabolism). ProteomeLM↔Geptop score r = 0.59.
- Cross-axis: **154 kp prime targets** (essential ∧ ligandable-tractable ∧ broad-spectrum-selective),
  143 ec. Shortlist `<prefix>_essentiality_shortlist.csv` = broad-selective ∧ essential (197 kp / 168 ec).

## Spot-checks (kp)
- Top by composite: rplF, ftsZ, gyrB, rplE, rpsG, rpoB, rpsE, secY, der, infB, rpoA — all textbook
  essential, each with ≥4 independent evidence sources.
- Top prime (essential × ligandable, broad-selective): rpsD, murC, ftsZ, rplF, mraY, murD, parC, gyrA,
  dxr — classic broad-spectrum, human-selective, druggable antibacterial targets (cell wall, gyrase/
  topo, ribosome, dxr=fosmidomycin).
- FBA essentials are all metabolic (murE, glmU, mraY, murC, murG, lpxK, fabA, asd, ilvD…), as expected.

## Reproduce
```
# gradi env (+ cobra, scikit-learn, openpyxl, ProteomeLM from git); DIAMOND from gradi-ortho
export GRADI_DIAMOND_BIN=$HOME/miniconda3/envs/gradi-ortho/bin/diamond
conda run -n gradi python scripts/07a_fetch_essentiality.py                                 # once (network)
for org in kpneumoniae ecoli; do
  conda run -n gradi python scripts/07b_kp_experimental.py   --organism $org
  conda run -n gradi python scripts/07c_ecoli_transfer.py    --organism $org
  conda run -n gradi python scripts/07d_proteomelm_ess.py    --organism $org   # trains head on ec first
  conda run -n gradi python scripts/07e_geptop.py            --organism $org
  conda run -n gradi python scripts/07f_fba_iyl1228.py       --organism $org
  conda run -n gradi python scripts/07g_deeplyessential.py   --organism $org   # deferred placeholder
  conda run -n gradi python scripts/07h_essentiality_merge.py --organism $org
  conda run -n gradi python scripts/07l_publication_essentiality.py --organism $org  # experimental-only
  conda run -n gradi python scripts/07i_predictor_plots.py   --organism $org
  conda run -n gradi python scripts/07j_essentiality_plots.py --organism $org
  conda run -n gradi python scripts/07k_essentiality_landscape.py --organism $org
  conda run -n gradi python scripts/07m_publication_plots.py  --organism $org        # publications slide
done
# outputs: output/results/<org>/<prefix>_essentiality{,_shortlist}.csv ; output/plots/07{i,j,k}_*.png
```
Everything is cached / resumable (DIAMOND hits, ProteomeLM ctx embeddings + head, FBA KO table,
downloaded files), so a killed run resumes on re-invocation.

## E. coli made first-class (2026-07-13)

E. coli was previously only a transfer/reference organism (all-NA experimental table; EcoGene demoted
to the 0.20 transfer axis; 4/6 publication panels were "Kp-specific" placeholders). It is now a
first-class experimental target with its own screens, symmetric with Kp.

**E. coli screens ingested (`07n_ecoli_experimental.py` → `ec_ess_experimental.csv`), keyed by
b-number / Keio-JW / gene symbol (near-complete mapping):**
- **Keio KO** (Baba 2006) via **PEC** `PECData.dat` (`Class==1`) — direct download. 287 essential.
- **Goodall 2018** BW25113 **TraDIS** (mBio Table S4) — fetched via **Chrome** (ASM 403). 353 essential.
- **Rousset 2018** genome-wide **CRISPRi** (PLoS Genet, PMC6242692) — Europe PMC zip; gene median log2FC.
- **Wang 2018** pooled **CRISPRi** (Nat Commun, PMC6018678) — Europe PMC zip; gene fitness (essential < −6).
- **Cui 2018** CRISPRi (PMC5954155) — Europe PMC zip (dCas9-toxicity filter reference).
- **Rousset 2021** 18-strain **CRISPRi** (Nat Microbiol) — Springer CDN MOESM (browser-UA); multistrain frac.
- **Hawkins 2020** mismatch-**CRISPRi** vulnerability curves (Cell Systems) — fetched via **Chrome** (Cell CDN).
- **RB-TnSeq / Fitness Browser** (Price 2018) — direct `morgannprice.org/FEBA/Keio/`; **3,790 genes ×
  3,511 conditions** (247 MB) processed to a per-gene min-fitness summary + a **280-condition
  antibiotic/stress matrix** (the giant raw file is not retained). Powers the 07o condition view.

**Integration:** `07c` now also lifts the E. coli screen consensus onto Kp via ortholog
(`ec_screens_essential_transfer`, 370 Kp proteins). `07l` populates real E. coli columns (first-class,
not compendium-only). `07h` feeds E. coli's own screens into `evidence_experimental` (the 0.40 axis) —
coverage rose from ~0 to **4,122/4,403** E. coli proteins. `07m` is now organism-symmetric (real
E. coli screen panels); `07o_condition_plots.py` adds a condition/stress slide (E. coli antibiotic
sensitivity from RB-TnSeq; Kp host-niche urine/serum/in-vivo). The general E. coli vulnerability score
deliberately excludes RB-TnSeq's all-condition min-fitness (that is condition-specific → 07o), so
lacZ/araC score ~0.05 not ~0.8.

**E. coli coverage:** 399 experimentally essential (Keio 287 / Goodall 353 / Rousset18 374 / Wang18 228;
346 in ≥2 screens, 230 in all 3 method types KO∩TraDIS∩CRISPRi). Composite tiers: **402 essential /
135 likely / 3,866 non**. RB-TnSeq: vancomycin 742 / gentamicin 388 / ciprofloxacin 279 genes sensitized.

**Only one dataset ungettable:** **Nichols 2011** chemical-genomics S-score matrix (Cell) — the PMC
copy sits behind an image reCAPTCHA (not solved on principle); optional/redundant with the RB-TnSeq
280-condition matrix. To add it manually: download `PMC3060659/bin/NIHMS261392-supplement-02.txt` in a
browser into `data/raw/ecoli/essentiality/nichols2011_chemgen/`.

## Publication-track coverage (kp)
- Experimentally essential (union of CRISPRi library ∪ KPNIH1 ∪ ECL8): **859** HS11286 proteins.
- CRISPRi 870-gene library mapped: 846/870. KPNIH1 essential (Ramage): 212 (by symbol). ECL8 essential: 249.
- Cross-species: 3,170 Kp proteins covered by the 12-genome compendium; **249 core-essential**
  (essential in ≥80% of genomes). Per-genome experimental essential-set sizes 251–391 (decoded).
- Screen agreement: 193 genes essential in all 3 direct Kp screens (CRISPRi ∩ KPNIH1 ∩ ECL8).

## Deferred / documented gaps
- **DeeplyEssential (4.3c)** — no released weights, license-less, py2/TF1.6; NaN placeholder (07g).
- **Bacformer** — optional extra predictor; Apple-Silicon flash-attn friction + genome-order prep; not run.
- **Bachman 2015 / 2025** staged; per-gene in-vivo comes from Mike&Bachman 2023 + the Jana in-vivo sheet
  (Bachman raw insertion tables not separately aggregated).
- **Ramage 2017** full paper is not OA, but its 424-gene KPNIH1 essential set is recovered from the Jana
  2023 supp. **Cain 2017** NJST258 not fetched (compendium lacks NJST258; low priority).
- **OGEE v3 / DEG** — no working programmatic endpoint; skipped (both repackage the primary screens).
- Spec extensions 4.4 (graded CRISPRi-vulnerability refit) / 4.5 (synthetic-lethality) — future work.

## CRISPRi/publication evidence folded into the composite (07h)
The 07l publication calls now feed `evidence_experimental` (07h reads `<prefix>_ess_publications.csv`):
CRISPRi CE-library membership → 0.85, KPNIH1 Tn-seq essential / ECL8-essential → 1.0, CRISPRi in-vivo /
vulnerability → 0.7, niche → 0.6 (max over signals). This raised the experimental sub-score coverage
from 956 → **1,471** kp proteins. The **hard** essential-tier override (`experimentally_essential`) stays
reserved for rigorous genome-wide calls (ECL8 / KPNIH1 Tn-seq, E. coli EcoGene) — CRISPRi CE-library
membership boosts the score and guarantees ≥ `likely_essential`, but does not blanket-force `essential`
(a conditionally-essential selection is softer than a bimodal Tn-seq call). Net kp tiers:
**401 essential / 835 likely_essential / 4,492 non_essential**; broad-selective∧essential shortlist 207.
`evidence_sources` now records `CRISPRi` / `KPNIH1` / `CRISPRi_invivo` provenance.
