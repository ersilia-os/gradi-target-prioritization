# Degradability — download checklist

Everything to fetch for the degradability axis, with the **exact file** to look for and a **route verified on
2026-08-04** (HTTP 200/206 unless marked otherwise). Companion to `docs/degradability_report.md` §4, which
explains *why* each dataset is wanted; this file is purely operational.

## The one trick that does most of the work

For every open-access paper below, a **single URL returns a zip of all supplementary files**:

```
https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/supplementaryFiles
```

Verified working for 15 of the papers here. It sidesteps two traps: publisher pages 403 scripted clients
(pnas.org, cell.com, jbc.org, wiley, ACS all returned 403), and direct
`pmc.ncbi.nlm.nih.gov/articles/instance/<id>/bin/<file>` URLs return an HTML interstitial rather than the file.
Where Europe PMC has no zip, the NCBI OA tarball is the fallback:

```
https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=<PMCID>   →  ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/…tar.gz
```

Legend: **✅** verified fetchable headlessly · **🔑** needs an authenticated browser session
(chrome-devtools `evaluate_script` same-origin fetch, as used for the gated essentiality tables) · **⚠** route
works but the exact file needs confirming.

---

## A. Priority — these six change the conclusions

| # | Paper | What you get | Exact file(s) | Route |
|---|---|---|---|---|
| A1 | **Nagar 2021** *mSystems* 6:e01296-20 · PMC7857536 | Pulsed-SILAC **half-lives for 1,149 E. coli proteins** + a ready-made **188-feature matrix** (incl. 128 PPI-network features) + stability classes. Keyed by UniProt. **The single biggest unlock** — CC-BY, fully automatable, and it replaces the hand-curated 35-row stand-in we currently use. | `mSystems.01296-20-sd001.xlsx` (raw time course, 1,602)<br>`mSystems.01296-20-sd002.xlsx` (**half-lives, 1,149**)<br>`mSystems.01296-20-sd003.xlsx` (stability classes) | ✅ Europe PMC zip (4.3 MB, 24 files) |
| A2 | **Gupta 2024** *Nat Commun* 15:5890 · PMC11246515 | `k_deg` + half-life for **~3,200 proteins × 13 conditions**, **plus protease attribution** from Δ*clpP*/Δ*lon*/Δ*hslV*/triple/Δ*smpB* (**no Δ*ftsH*** — see the 2026-08-06 correction), plus ~600 measured in-vivo N-termini. The axis's primary label set. | `41467_2024_49920_MOESM4_ESM.xlsx` (**Supp Data 1 — turnover**)<br>`…MOESM5_ESM.xlsx` (**Supp Data 2 — protease assignment**)<br>`…MOESM7_ESM.xlsx` (**Supp Data 4 — N-termini**)<br>`…MOESM6/8/9/10/11_ESM.xlsx` (rest) | ✅ Europe PMC zip (34.7 MB) or direct Springer `static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-49920-8/MediaObjects/<file>` |
| A3 | **Won 2024** *Nat Commun* 15:4065 · PMC11094019 | The **only validated quantitative bacterial degradability model**. Supp Data 1 = 348-protein feature matrix; Supp Data 2 = **72 targets with measured degradation constants** (the labels to test against); Supp Data 3 = predicted scores. Plus code on Zenodo. | `41467_2024_48506_MOESM4_ESM.xlsx`<br>`…MOESM5_ESM.xlsx`, `…MOESM6_ESM.xlsx`, `…MOESM8_ESM.xlsx` | ✅ Europe PMC zip (18.1 MB)<br>✅ Zenodo `10.5281/zenodo.11004395` (code + raw) |
| A4 | **DEtox** — Beardslee & Schmitz 2024 *eLife* · PMC10862746 | The **empirical C-terminal ClpXP degron PSSM** from ~100,000 measured pentapeptide tags, run in WT/Δ*clpX*/Δ*clpP*. Replaces the hand-written ssrA regex (which cannot even match ssrA). | `media-2.xlsx`, `media-3.xlsx` | ✅ Europe PMC zip (35.7 MB) |
| A5 | **Cragan 2025** *JBC* 301:108365 · PMC11986505 | The **Lon C-degron consensus** `x–[L/I]–[L/I/V]–H-COOH`, conserved across E. coli / Yersinia / Mycoplasma. A rule the old spec had nothing for. | `mmc1.xlsx` (+ `mmc2.docx`, `mmc3.pdf`) | ✅ Europe PMC zip (7.4 MB) |
| A6 | **Sen 2025 "N-FIVE"** bioRxiv 2025.05.22.655665v2 · PMC12258705 | The **N-degron P1–P5 stability model** from ~2.2 M variants — the thing that replaces (and explains the failure of) our inverted N-end-rule feature. **The data and model are on GitHub, not in the supplement.** | GitHub: `KunjapurLab/N-terminal-cluster-stability`<br>PMC zip has only `NIHPP…supplement-1.pdf` | ✅ GitHub<br>✅ Europe PMC zip (2.9 MB, PDF only) |

---

## A′. E. coli-first additions — found 2026-08-05, all routes tested

Everything in this section came out of the E. coli-focused pass and is **new since section A was written**.
Section A remains correct; this is what an E. coli-first build needs on top of it.

| # | Resource | What you get | Route |
|---|---|---|---|
| **A′1** | **Gupta 2024 Table S1 — read it again, more carefully** | The earlier entry undersold this. Table S1 is **half-lives for 3,263 proteins across 13 conditions** (measured from the fetched file 2026-08-06; earlier text said 3,262 x 14), and the conditions include **`clpP N-lim6`, `lon N-lim6`, `hslV N-lim6`, `Triple (ΔclpPΔlonΔhslV)` and `smpB N-lim6`**. Identifiers are `sp\|P00350\|6PGD_ECOLI` — **UniProt accessions, no mapping work at all**. A per-protein change in half-life on protease deletion is a **ground-truth label, not a feature**. This single file is the reason to build E. coli first. | ✅ `…/41467_2024_49920_MOESM4_ESM.xlsx` (1.45 MB) |
| **A′2** | **TRAINSPOTTER** — Van Damme et al., *NAR* 54:gkag587 (2026), PMC13245402, CC-BY | Deformylation-assisted N-terminomics in **K-12** (CAG12184): **1,082 N-termini → 806 initiation sites → 729 proteins**, with iMet-excised vs retained, alternative start codons, N-terminal extensions/truncations and Nt-acetylation. Table S1 columns are ingestion-ready and carry **both UniProt accession and b-number**, plus `NH2-M/unique peptide` — the observed mature N-terminus. Fixes the "annotated Met1 is often wrong" problem that the whole N-degron channel depends on. **Already downloaded to the session scratchpad.** | ✅ Europe PMC zip → `gkag587_supplemental_files.zip` → `Table S1.xlsx` (891 KB). The OUP CDN link 403s; PMC OA tarball 404s. Raw: PXD005901 |
| **A′3** | **E. coli PeptideAtlas build 585** — Jachmann et al., *J Proteome Res* 25:1027 (2026), PMC12887991 | `coordinate_mapping.txt` (99 MB): for **every peptide observed in any aggregated E. coli MS experiment**, its start/end position and preceding/following residue. Filter `start_pos ∈ {1,2}` → empirical mature N-termini; `end_pos == length` → native C-termini; non-tryptic preceding residue → candidate in-vivo cleavage sites. A **proteome-scale** termini prior for one download, complementing A′2's depth with breadth. Peptide-level evidence, so absence is weak evidence and tryptic bias applies. | ✅ `peptideatlas.org/builds/ecoli/202411/` — no login |
| **A′4** | **Li 2014** — *Cell* 157:624, PMID 24766808 | **Absolute protein synthesis rates**, molecules per cell per generation — the true "will it stay degraded" denominator, better than abundance. Table S1 = 4,095 genes (**3,041 unbracketed = high-confidence**); Table S4 = mRNA + translation efficiency (2,387 usable); Table S3 = 177 rows of complex stoichiometry. Identifier is a **bare lowercase gene name**, so budget a synonym crosswalk for ~360 rows. | ✅ **`ars.els-cdn.com/content/image/1-s2.0-S0092867414002323-mmc{1,3,4}.xlsx`** — this CDN pattern works for every Cell-press paper tested, unlike PMC `/bin/` (reCAPTCHA) and cell.com (403) |
| **A′5** | **Bienvenut 2015** — *Proteomics* 15:2503, PMID 26017780 | The canonical E. coli **mature N-terminome with quantitative modification fractions**: >1,000 N-termini, **56% initiator-Met removed, 10% Nt-acetylated, 5% N-formyl retained**, plus 140 signal-peptide-cleaved membrane-protein termini. Gives a *probability* the protein presents a given residue, not a boolean. | 🔑 Wiley supplementary is paywalled and **PRIDE holds RAW only** (PXD001979/002012/001983 = 166 RAW + Mascot `.dat`, no processed table). Institutional access or email the authors |
| **A′6** | **CPB-ChaFRADIC C-terminome** — Chen et al., *Anal Chem* 92 (2020), PMID 32441514, PXD018520 | The **only native E. coli C-terminome**: 604 canonical + 818 neo-C-termini. Our C-degron channel is currently calibrated entirely on *engineered* library tags — this is the real distribution those PSSMs should be scored against, and the 818 neo-termini are direct evidence of in-vivo endoproteolysis. | 🔑 ACS paywall for tables; PRIDE for raw |
| **A′7** | **Klimecka 2021 ssrA benchmark** — *Molecules* 26:5936, PMC8512704 | 12 eGFP–degron variants each measured four ways: **SspB K_D by MST** (0.08–0.33 µM), ClpX ATPase stimulation, **in vitro K_M/V_max ± SspB** (AANDENYALAA: 3.38 → 0.47 µM with SspB), and in vivo time courses in WT/Δ*clpX*/Δ*clpP*/Δ*sspB*. The only source that **separates adaptor binding from protease engagement** — exactly the decomposition a recruited-protease score needs. | ✅ paper open; ⚠ **no consolidated table** — ~12×4 numbers must be hand-curated from Fig. 4e/Fig. 7 |
| **A′8** | **Complex Portal · STRING · RCSB · PDBe · Seq2Symm** | See §F (F9–F13). Together these turn assembly state from "the weakest aspect" into a solved one for E. coli. | ✅ all |

### New modality and rule papers (2025–2026)

| Resource | What it adds | Route |
|---|---|---|
| **Han et al. 2025**, *Nat Commun* 16, `10.1038/s41467-025-66221-w`, PMID 41413196 | **A second ClpXP handle**: "bacNID" grafts a target-binding peptide *and* an **SspB-binding peptide** onto a gold nanoparticle, hijacking the SspB adaptor rather than ClpX directly. Includes resistance-evolution and an in vivo wound model. Peptide-on-nanoparticle, not a small molecule — permeability is solved by the carrier. | ✅ OA |
| **Coriano-Ortiz et al. 2026**, bioRxiv `10.64898/2026.06.17.732960` | **Lon's first working handle in E. coli.** Orthogonal *Mesoplasma florum* Lon with a LOV2 domain inserted at nearly every codon (726 variants screened); light-tunable, and it **degrades in stationary phase**, decoupling degradation from growth dilution. | ✅ preprint |
| **Iqbal, Keller & Ghanbarpour 2026**, *Cell Rep* 45:117231, PMID 41964960 | **FtsH degradability is gated, not intrinsic.** FtsH sits in an **HflK/HflC** cage; cryo-EM of a disulfide-locked closed state, and cells locked shut fail to recover from aminoglycosides. The FtsH channel needs an accessibility term alongside the degron. | ✅ OA |
| **Islam et al. 2026**, *Biophys J*, PMID 41520171 | For **ClpAP, unfolding is rate-limiting, not translocation** — 12.0 aa s⁻¹ (ClpA) vs 33.2 (ClpAP); translocation is 8–24× faster than unfolding either way. So the ClpAP term should be local mechanical stability next to the degron, not chain length. | 🔑 |
| **González et al. 2025**, *Nat Commun*, PMID 40993131 | In-cell NMR: **Prc cleaves membrane-bound NDM-1 first, and DegP then processes the fragments.** For β-lactamases the periplasmic route may be Prc→DegP, not DegP alone — so the DegP channel needs a Prc-accessibility term. | ✅ OA |
| **Roy, Nandakumar & Chaba 2025**, bioRxiv `10.1101/2025.11.03.686187` | **DegP is redox-gated**: a disulfide in its own protease domain acts as a sensor, and substrates carrying free cysteines can allosterically convert it. The DegP score needs a Cys/DSB term. | ✅ preprint |
| **Lee et al. 2026 "ProHL"**, *J Microbiol Biotechnol*, PMID 42046921 | **The direct competitor** — a deep-learning bacterial protein half-life classifier validated in E. coli. ⚠ **Read the corrigendum first** (`10.4014/jmb.2026.3602.c02`, PMID 42124426) before benchmarking against it. | ✅ OA |
| **Izert-Nowakowska et al. 2026**, bioRxiv `10.64898/2026.03.07.710301` | The **validation assay**, from the CLIPPERs group: pBAD eGFP–degron fusions with a plate end-point and a 96-well kinetic readout for high-throughput degron screening. Treat any degron set they report as a labelled test set. | ✅ preprint |
| **Stevens-Cullinane et al. 2025**, *JACS* 147, PMID 41308195 | ⚠ **Negative control, not a training example.** A Ru(II) photosensitiser degrades NDM-1 in a Gram-negative with 53-fold meropenem potentiation — but by **photo-oxidative backbone scission, explicitly without recruiting host proteolytic machinery**. Must not contaminate a "degradable by AAA+ protease" label. | 🔑 |
| **flDPnn3** (*J Mol Biol* 2026, PMID 41500381) · **DisProt 2026** (*NAR*, PMID 41249866) | Proteome-scale disorder prediction and curated disorder ground truth. Relevant because Won 2024's top feature used flDPnn specifically. | ✅ |

---

## B. The biophysics layer (E. coli, transferred onto Kp by ortholog)

| # | Paper | What you get | Exact file(s) | Route |
|---|---|---|---|---|
| B1 | **Cappelletti 2021** *Cell* 184:545 · PMC7836100 | **LiP-MS protease-accessibility density** for ~1,900 proteins across 8 carbon sources. Mechanistically the most apt single column — limited proteolysis measures exactly the accessible-flexible-region property. Strain BW25113, same as B2. | `mmc2.xlsx` … `mmc7.xlsx` / `mmc3.xls` (Table S3 is the E. coli per-protein/peptide table) | ✅ Europe PMC zip (34 MB, 32 files) ⚠ confirm which mmc is Table S3 |
| B2 | **Mateus 2018** *Mol Syst Biol* 14:e8242 · PMC6056769 | **Melting temperature (Tm) for 1,738 proteins**, intact cells. CC-BY. Residualise on localisation before use — Tm tracks compartment. | `MSB-14-e8242-s002.xlsx` (**Dataset EV1 = per-protein Tm**) ⚠ confirm against the file list; the zip also has s003/s008/s009/s011/s013 | ✅ Europe PMC zip (44.1 MB, 32 files) |
| B3 | **Mateus 2020** *Nature* 588:473 · PMC7612278 | Thermal-stability scores across **121 Keio deletion backgrounds** → the derived feature worth having is the **variance** (conformational plasticity), arguably better than a single Tm. | `41586_2020_3002_MOESM8_ESM.xlsx` (Supp Data 5, 17.5 MB) | ✅ direct Springer static-content |
| B4 | **eSOL / Niwa 2009** *PNAS* 106:4201 | **Chaperone-free % solubility for 3,198 proteins** — the best coverage of any dataset here, and keyed on **b-numbers** so joins are trivial. | `esol.zip` | ✅ `https://dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip`<br>(the original tanpaku.org site is dead) |
| B5 | **Niwa 2012** *PNAS* 109:8937 | **ΔSolubility with GroEL / trigger factor / DnaK** = per-protein chaperone dependence. | supplementary tables | 🔑 pnas.org 403s headlessly |
| B6 | **Györkei 2022** *Sci Rep* 12:6547 · PMC9023497 | **In-vivo solubility limit + 3-class aggregation-rate label** for 2,577 cytosolic proteins. | `41598_2022_10427_MOESM3…MOESM11_ESM.xlsx` | ✅ Europe PMC zip (3.0 MB)<br>✅ raw images zip at `group.szbk.u-szeged.hu/sysbiol/…/Proteome_wide_landscape_of_solubility_limits_in_a_bacterial.zip` |
| B7 | **To 2021** *JACS* 143:11435 | **(Non-)refoldability**: 396/1,198 (33%) non-refoldable, protein- **and domain-level** tables. Kinetic-stability proxy. | article + SI tables | ✅ open PDF at `par.nsf.gov/servlets/purl/10311595`<br>🔑 ACS SI tables (pubs.acs.org 403s); bioRxiv 2020.08.28.273110 is the other open copy |
| B8 | **To 2022** *PNAS* 119 | Chaperone-**assisted** refolding map — separates intrinsic from rescued refoldability. | supplementary datasets | 🔑 pnas.org 403s |
| B9 | **Schmidt 2016** *Nat Biotechnol* 34:104 · PMC4888949 | **Copies per cell × 22 conditions**, 2,359 proteins. The abundance term (Won 2024 measured r = −0.69 against degradation rate). | `NIHMS65833-supplement-Supplementary_tables.xlsx` (16 MB) | ✅ Europe PMC zip (24.1 MB) |
| B10 | **Calloni 2012** *Cell Rep* 1:251 | **DnaK interactome (~700 proteins)** + GroEL client classes I/II/III (30/80/42). Chaperone dependence as metastability. | supplementary tables | 🔑 cell.com 403s |
| B11 | **Jarzab 2020** Meltome Atlas *Nat Methods* 17:495 | Cross-species Tm — useful only as a prior to sanity-check ortholog transfer; for E. coli it duplicates B2. | `41592_2020_801_MOESM4_ESM.xlsx` (13.9 MB) | ✅ direct Springer ⚠ file identity inferred from size — open it before trusting |

---

## C. Kp-native anchors — the only *Klebsiella* data that exists

There is **no** Kp turnover, meltome or LiP-MS dataset. These two exist to check that the E. coli transfer
isn't systematically off.

| # | Paper | What you get | Exact file(s) | Route |
|---|---|---|---|---|
| C1 | **Illenseher 2025** *Front Microbiol* 16:1528869 · PMC12127431 | **Kp-native iBAQ abundance**, >2,800 proteins, ± heat and oxidative stress. Strain ATCC BAA-2146 Δ*wza*; ids are `Kpn2146_####` locus tags so needs DIAMOND mapping onto HS11286. | `Table_1.xlsx`, `Table_2.xlsx` | ✅ Europe PMC zip (8.1 MB) |
| C2 | **Muselius 2020** *Front Microbiol* 11:546 · PMC7194016 | **The only Kp proteolysis dataset in existence**: Δ*lon* vs WT (K52), 2,074 proteins quantified, **26 accumulate = candidate Kp Lon substrates**. | `Table_7.XLSX` (Δ*lon* vs WT) — also `Table_1`–`Table_6.XLSX` | ✅ Europe PMC zip (1.2 MB) |

---

## D. Modality and rules — read these, no data to ingest

| # | Paper | Why | Route |
|---|---|---|---|
| D1 | **Izert-Nowakowska 2025** CLIPPERs *EMBO Rep* 26:3994 · PMC12373786 | The only demonstration of degrading a native protein in a Gram-negative, and the source of the substrate-requirement checklist. Read the Appendix for the targets that **failed**. | ✅ Europe PMC zip (176 MB — large; the appendix PDFs are `44319_2025_510_MOESM1/15_ESM.pdf`) |
| D2 | **Chai-Danino 2026** *Nat Commun* 17:3067 | The FtsH rule: lipid-facing polar residues in TM helices, sufficient even for folded proteins. Opens the 1,320 inner-membrane Kp proteins. | ✅ Nature open access; preprint bioRxiv 2023.12.12.571171 |
| D3 | **Taylor 2026** DegP glues, bioRxiv `10.64898/2026.03.30.715243` | The periplasmic handle — the most drug-like Gram-negative result there is. **Read manually**, bioRxiv returned 500 to scripted fetch. | 🔑 manual |
| D4 | **Saunders 2020** *PNAS* 117:28005 + **Kenniston 2005** *PNAS* 102:1390 | The 5 / 20 / 37-residue initiation-region numbers. | 🔑 pnas.org 403s (both are old enough for PMC: PMC547888 for Kenniston) |
| D5 | **Nie 2025** *Chem Commun* 61:13149 | The cautionary tale — 1,669 Da, and the periplasmic target was relocalised to the cytoplasm to make it work. | 🔑 |
| D6 | Reviews: **Junker & Clausen 2026** *Trends Biochem Sci* 51:249 · **Petkov 2023** *Biochem J* 480:1719 (PMC10657178) · **Bazzacco 2025** *FEMS Microbiol Rev* · **Izert 2021** *Front Mol Biosci* 8:669762 (PMC8138137) | Orientation; Petkov is the one that states the Gram-positive ClpC/ClpE vs Gram-negative ClpX/ClpA split. | ✅ the two PMC ones; 🔑 the others |

---

## E. Substrate-trap sets (labels for validation)

| # | Source | What you get | Route |
|---|---|---|---|
| E1 | **Niwa 2022** *Molecules* 27:3772 · PMC9228906 | **Lon vs ClpXP vs HslUV attribution** over ~80 obligate GroE substrates — the best public "which protease eats which protein" comparison. | ✅ Europe PMC zip; data in `molecules-27-03772-s001.zip` |
| E2 | **Neher 2006** *Mol Cell* 22:193 | ClpXP trap ± DNA damage; 25% of the SOS proteome; 9 confirmed substrates. | 🔑 cell.com |
| E3 | **Westphal 2012** *JBC* 287:42962 · PMC3522291 + **Arends 2016** *Proteomics* 16:3161 | FtsH trap sets (15 and >50 candidates). | Europe PMC has **no** supp zip for PMC3522291 → use the article page; Arends is 🔑 Wiley |
| E4 | **Humbard 2013** *JBC* 288:28913 · PMC3789986 | **N-degradome**: >100 ClpS binders, Table 1 = measured neo-N-termini. The evidence that most N-degrons are generated by prior proteolysis. | Europe PMC has **no** supp zip; jbc.org PDF 403s → 🔑 or the PMC article page |
| E5 | **Flynn 2003** *Mol Cell* 11:671 | The five ClpX motif classes **and the trapped-substrate census**. Main-text tables, not a supplement — so an authenticated session should get them. Currently the one primary source we have never actually read. | 🔑 **highest-value single fetch** |
| E6 | Cross-bacterial traps: **Lunge 2020** (PMC7363115), **Bhat 2013** (PMC3681837), **Graham 2013** (PMC3807464), **Feng 2013** (ACS), **Ziemski 2021** (Wiley) | The old spec's §3.3b pool. **All gated** — the three PMC ids return `idIsNotOpenAccess`, the other two have no PMC record. Lowest priority: they're small label sets and the open E. coli turnover data supersedes them. | 🔑 all |

---

## F. Databases

**Bottom line: there is no database of bacterial protein degradability.** Nothing aggregates degrons,
half-lives and protease assignment for a bacterium the way DEG does for essentiality or ChEMBL does for
bioactivity. The axis has to be assembled from the papers in §A–§E. What *does* exist is a set of
single-property databases, two of which are genuinely worth using instead of recomputing.

### Worth using

| # | Database | What it gives | Covers Kp HS11286? | Route |
|---|---|---|---|---|
| **F1** | **MobiDB** | Per-accession **disorder** — a consensus plus ~13 individual predictors, including `prediction-disorder-alphafold`, with `start..end` regions and `content_fraction`. This is the disorder layer for track 3.1 without running anything ourselves, and it is citable. Verified against our own accessions: ClpX `A0A0H3GSR4` (len 424, 14 feature tracks), ClpP `A0A0H3GKH6` (13 tracks), *E. coli* ClpX `P0A6H1` (108 tracks — E. coli is far more richly annotated, incl. derived LIP/binding-mode features). | **Yes** — and **proteome-wide bulk download works**: `https://mobidb.org/api/download?proteome=UP000007841&format=tsv` returns TSV keyed by accession | ✅ |
| **F2** | **PHDB** — Protein Homeostasis Database (Ramakrishnan 2020, *Bioinformatics* 36:948) | **4,305 E. coli proteins × >100 proteostatic parameters** in one place: 13 whole-proteome chaperone studies, abundance, translation rate, solubility, net charge, secondary-structure content, contact order, TANGO aggregation propensity, IUPred disorder. If it has a bulk export this collapses several §B rows into one fetch. The paper describes a web UI only. | E. coli only → transfer by ortholog | ✅ site up (`phdb.switchlab.org`); **bulk export still unconfirmed** |
| F3 | **eSOL** | Chaperone-free % solubility, 3,198 E. coli proteins, b-number keyed. Same content as §B4 but browsable. | E. coli only | ✅ `togodb.biosciencedbc.jp/togodb/view/esol` (bulk zip in §B4) |
| F4 | **STEPdb 2.0** | E. coli localisation / topology / abundance / solubility / disorder, UniProt-keyed spreadsheets. Several columns are *predicted*, not experimental — label them as such. | E. coli only | ✅ `stepdb.eu` |
| F5 | **AlphaKnot 2.0** | Knot topology of AlphaFold-predicted models → the soft topological-resistance penalty in track 3.4. | Yes (AFDB-wide) | ✅ `alphaknot.cent.uw.edu.pl` |
| F6 | **PaxDb** (Kp = species 667127) | ppm abundance, integrated across studies. Secondary to §C1 for Kp and §B9 for E. coli. Site is JS-rendered → use the bulk release files. | Yes | ✅ `pax-db.org/downloads` |
| F7 | **MEROPS** | The peptidase reference: Lon = clan **S16**, ClpP = **S14**, with per-peptidase cleavage-specificity and substrate entries. Protease-centric, so it answers "what does Lon cut" rather than "is this protein degradable" — useful for the cleavage-site side of tracks 3.2b/3.4, not as a per-protein column. | protease entries only | ✅ `ebi.ac.uk/merops` (S16.001, S14.001 verified) |

### Checked and *not* usable for this axis

| Database | Why not |
|---|---|
| **TopFIND** | **RESOLVED, and it is empty for our purpose.** Live at `topfind.clip.msl.ubc.ca` (v4.1); the bulk dump really does download — `/assets/TopFIND_20211220.sql.gz`, 232 MB, no login. But 10 of 10 sampled E. coli entries returned evidence type *"electronic annotation"* with **zero literature PMIDs**, 0 cleavages and 0 substrates. Its E. coli layer is a re-projection of UniProt processing annotation plus MEROPS — and MEROPS is itself empty here (next row). Frozen at UniProt 2019_07. Gene-name search 500s; only accessions work. **Resolved and deliberately not ingested.** |
| **MEROPS** | **Quantified dead end.** The undocumented substrate file works — `ftp.ebi.ac.uk/pub/databases/merops/current_release/Substrate_search.txt`, 29 MB, 108,196 rows with P4–P4′ and a physiological/synthetic flag. Counted for our four proteases: **Lon (S16.001) = 2 real substrates** (SulA, CcdA), **ClpP (S14.*) = 0 physiological — the only entry is bovine insulin B-chain**, **FtsH (M41.001) = 1** (β-casein), **HslV (T1.006/7) = 0 rows entirely.** Ingest only as a generic cross-species P4–P4′ prior. |
| **PortEco / EcoGene / EcoliWiki** | **PortEco's domain is now an online-gambling spam site** — do not cite or link it. **EcoGene** returns HTTP 000, server gone. **EcoliWiki** is alive but holds no degradation, half-life or protease-substrate curation. |
| **EcoCyc** | Alive, and its release notes do advertise protein turnover — but **that dataset *is* Gupta 2024**, which we already have. Not exposed via the free `getxml` API; flat files need a BioCyc licence request. No net-new data. |
| **QSbio / anti-QSalign** | **Dead.** `qsbio.org` redirects to a Weizmann page that now serves an empty Drupal shell with no data links at all. The only surviving QSbio assignments are second-hand, in Schweke's `Q1qsbio` columns and Complex Portal's `qsproteome:` cross-references. |
| **Babu 2018** cell-envelope complexes | Triple dead end, all tested: not open access (`idIsNotOpenAccess`), the IntAct per-PMID export is empty (7 lines, 0 accessions), and a 48-way probe of Springer MOESM filenames returned zero hits. |
| **DEGRONOPEDIA** | Architecture is exactly right — degrons mapped to nearby disorder ("unfolding seeds") plus simulated proteolysis at newly generated termini. But the Methods list **11 supported organisms and *E. coli* is not among them**; all 46,743 degrons are ubiquitin-system motifs. No batch queries. Cite as the design precedent for the degron track, never as a data source. |
| **TPDdb** (*Nucleic Acids Res* 2026, 54:D1683) | A new curated database of **targeted protein degraders** — degradation activity, binding affinity, cytotoxicity per compound. Compound-centric and overwhelmingly human/oncology; no bacterial degraders of substance exist to curate. Also returned 403 to us. Watch it, don't depend on it. |
| **DEG / OGEE** | Essentiality, not degradability — already covered by the `07*` axis. |
| **PRIDE / ProteomeXchange** | Raw mass spec (see §G), not per-protein annotation. |

### Code repositories

**Gupta 2024**: `github.com/wuhrlab/ProteinTurnoverEcoli` — analysis notebooks only; the tidy tables are the
§A2 supplements. **N-FIVE**: `github.com/KunjapurLab/N-terminal-cluster-stability` — the model *and* the data.
**Mateus 2020**: `github.com/fstein/EcoliTPP`. **Won 2024**: Zenodo `10.5281/zenodo.11004395`.

| # | Also | What it gives | Route |
|---|---|---|---|
| F8 | **Schweke 2024** *Cell* 187:999 — homo-oligomer atlas | ⚠ **CORRECTION.** An earlier version of this file listed this as an *E. coli* resource. It is **E. coli O157:H7, not K-12**: the `org == "ec"` rows are 2,181 accessions whose **intersection with UP000000625 is exactly zero**, and only 820/2,181 match K-12 even by gene name. The table also lists *only* predicted homomers, so absence conflates "monomer" with "not predicted". Usable for K-12 only after an ortholog mapping step, yielding ~800–1,400 assignments. | supp downloads fine from the Elsevier CDN (`ars.els-cdn.com/content/image/1-s2.0-S009286742400059X-mmc3.xlsx`, 2.7 MB, `DATA_TABLE_S2`, 8,195 rows × 73 cols) |
| **F9** | **Seq2Symm** — Kshirsagar et al., *Nat Commun* 16:1969 (2025), `10.1038/s41467-025-57148-3`, PMC11868566 | **Use this instead of F8.** Sequence-only (protein-language-model) homo-oligomer symmetry classifier. The only route to **100% coverage of our exact 4,403 K-12 and 5,728 Kp sequences with zero strain mismatch and zero identifier mapping** — run locally, like `07d_proteomelm_ess.py` and `09c_deeplocpro.py` already are. | ✅ `github.com/microsoft/seq2symm` |
| **F10** | **RCSB assembly API** | **This is how assembly state gets solved.** The GraphQL `rcsb_struct_symmetry` field returns `oligomeric_state` literally as `"Homo 2-mer"` plus `stoichiometry` (`["A2"]`) and symmetry (`C2`). Search by taxonomy 83333 → **7,671 biological assemblies** for E. coli K-12; batch the GraphQL calls and join to UniProt. | ✅ `data.rcsb.org/graphql` + `search.rcsb.org` |
| **F11** | **PDBe `interface_residues`** | The mechanistically honest version of "is the terminus free". Returns interface residues **in UniProt numbering, derived from biological assemblies** (not the isolated monomer), *with the partner's identity* — so self-burial in a homomer (the GroEL case) is distinguishable from burial against a different subunit. ~4,403 calls, ~1,768 return data. | ✅ `ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/<ACC>` |
| **F12** | **Complex Portal** (EBI/IntAct) | 324 curated E. coli complexes, **781 unique participant accessions, with explicit stoichiometry** (`P23909(2)`) and a literal `Complex assembly` column (`Homodimer`, `Homotetramer`). Small, curated, gold-standard — the validation set for F10/F9. Also carries the only surviving QSbio ids. | ✅ `ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/83333.tsv` (554 KB) |
| **F13** | **STRING v12.0** | Per-protein interaction degree — Nagar 2021's single most predictive feature class. **Taxon is 511145** (not 83333, not 562). At `combined_score ≥ 700`: **4,006 nodes = 91% of the proteome**, median degree 10. Use `protein.physical.links` for the "buried in a complex" signal specifically. Join via UniProt `xref_string` (4,087/4,403), **not** via the aliases file — that maps one node to three accessions. | ✅ `stringdb-downloads.org/download/protein.physical.links.v12.0/511145.…txt.gz` |
| **F14** | **E. coli PeptideAtlas** build 585 — Jachmann et al., *J Proteome Res* 25:1027 (2026), PMC12887991 | The sleeper hit. `coordinate_mapping.txt` (99 MB) gives, for **every peptide observed in any aggregated E. coli MS experiment**, start/end position and the preceding/following residue. Filter `start_pos ∈ {1,2}` → empirical mature N-termini; `end_pos == length` → native C-termini; non-tryptic preceding residue → candidate in-vivo cleavage sites. A proteome-scale termini prior for one download. Evidence is peptide-level, so absence is weak evidence. | ✅ `peptideatlas.org/builds/ecoli/202411/` |

---

## G. Raw mass-spec (only if you want to reprocess)

Not needed for the axis — the processed tables above are sufficient. Listed for completeness:
**PXD042444** (Gupta 2024) · **PXD022112** (Nagar 2021) · **PXD062881** (MacKrell 2026) ·
**PXD022297** (Cappelletti LiP) · **PXD009495** (Mateus 2018 TPP) · **PXD016589** (Mateus 2020) ·
**PXD052921** (Illenseher, Kp) · **PXD015623** (Muselius, Kp Δ*lon*) · **PXD045730** (CLIPPERs) ·
**PXD011929** (Meltome, partial).

---

## H. Explicitly not worth downloading

- **Tsuboyama 2023** mega-scale ΔG (Zenodo `10.5281/zenodo.7992926`, ~1 GB). The 331 natural domains are
  40–72 aa, selected for cDNA-display compatibility, and not mapped to any proteome. Using it would mean
  training a ΔG predictor, which turns an experimental column into a predicted one. There is no proteome-scale
  experimental ΔG for any bacterium — Tm (B2) is the stand-in.
- **Wang 2016** *Salmonella* turnover (870 proteins, paywalled, no PMC) — too small to matter.
- **Piazza 2018** LiP-SMap — the openly available slice is only Data S8–S10; B1 supersedes it.

---

## Suggested order

**Revised 2026-08-05 for an E. coli-first build.** Fetch in this order; the first four are a day's work and
take the axis from unbuilt to defensible.

1. **A′1** (Gupta Table S1) — UniProt-keyed half-lives across 13 conditions *including the protease knockouts*.
   This is the label set. Nothing else should be built before it.
2. **A′2 + A′3** (TRAINSPOTTER, PeptideAtlas) — the real mature N-/C-termini. Without these the degron
   channels are scored on annotated termini that are often wrong.
3. **A4, A5, A6** (DEtox, Cragan, N-FIVE) — the three rule sets, all open, all small.
4. **F13 + F10/F11** (STRING, RCSB assemblies, PDBe interface residues) — the assembly-state and connectivity
   features, which turn out to be strong for E. coli rather than weak.
5. **A′4, B4** (Li synthesis rates, eSOL) — the stay-down denominator and the best-coverage biophysics column.
6. **A′7** (Klimecka benchmark) — calibration, once something needs calibrating.
7. Everything else as needed. **C1/C2** only when the Kp arm starts.

### Superseded ordering (Kp-first, kept for the record)

1. **A1, A2** — one Europe PMC call each. These two alone replace the fabricated inputs the current axis uses
   and give the ~3,200-protein label set the whole axis can be validated against.
2. **A3–A6** — the rule sets and the validated model. All open.
3. **B4, B9** — eSOL and copies-per-cell: trivial, high coverage, b-number keyed.
4. **C1, C2** — the two Kp-native anchors, so the transfer can be checked.
5. **B1, B2** — LiP and Tm, the two best biophysical columns.
6. **E5 (Flynn 2003)** via an authenticated session — the one primary source never actually read.
7. Everything else as needed.

---

## A″. Activated-ClpP (partnerless) evidence — added 2026-08-06

**These are the two most important downloads in this file** and neither was in it before. They are the only
datasets that measure the criterion the Gr-ADI proposal actually gates on ("activated ClpP … in the absence of
an unfoldase partner"), and the proposal picked its first two targets (DnaK, AcpP) off them.

Both are now automated in `scripts/10a_fetch_degradability.py` — `--only conlon2013_adep4_saureus` and
`--only jacques2020_onc212_saureus`. Nothing to do by hand. Recorded here because the working routes are
non-obvious and two of the three obvious ones look like paywalls but are not.

| Dataset | Route that works | Route that fails, and how it fails |
| --- | --- | --- |
| **Conlon 2013** *Nature* 503:365 (ADEP4) — Tables S1 (abundance, 1,712 proteins) + S2 (partially tryptic = **cleavage**, 631 proteins) | `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature12790/MediaObjects/41586_2013_BFnature12790_MOESM97_ESM.xlsx` — **free, 703 kB, no session needed** | Europe PMC → `"PMC4031760 is not open access one"`. `static-content.springer.com` (the *other* Springer host) → **403** to bare curl. PMC article page → **reCAPTCHA**. All three read as "paywalled"; the article is, the supplement is not. |
| **Jacques 2020** *Genetics* 214:1103 (ONC212) — Table S3 (1,620 proteins, abundance + **cleavage**) | GSA figshare `10.25386/genetics.11873841`, article id `11873841`, file id `21767271` | Europe PMC `supplementaryFiles` for PMC7153937 returns **only figure images** — a 200 OK with no data in it. figshare `ndownloader` answers **HTTP 202 with a zero-length body** while preparing the file, which `raise_for_status()` passes. |

Two ladder bugs these exposed, both now fixed in `10a`:
1. **HTTP 202 + empty body was banked as success** — 202 means "not ready, retry", and `raise_for_status()`
   does not raise on it. Now raises so tenacity backs off; an empty 200 is rejected the same way.
2. **`.xlsx` files were routed to the archive branch.** Office Open XML and ZIP share the `PK` magic bytes, and
   the zip test ran first — so a workbook was opened as an archive, found no member ending in a data extension
   (its real members are `xl/worksheets/*.xml`), and the dataset was reported **gated despite HTTP 200**. Now
   discriminated by the presence of a top-level `[Content_Types].xml`.

### Identifier bridge — why this needed sequence, not accessions

Both papers key on 2013-era **SACOL** locus tags (*S. aureus* COL). Do not waste time on these routes:

- Current RefSeq `GCF_000012045.1` has re-tagged everything `SACOL_RS#####` and keeps **no `old_locus_tag`** —
  the SACOL tags are simply not recoverable from the modern assembly.
- UniProt has **demoted the COL and Mu50 proteomes as redundant**: only 939 and 1,033 entries survive, so a
  UniProt-first bridge silently drops ~45% of Conlon.
- Gene-symbol matching fails on every hypothetical protein.

What works, and is what `10c_clpp_activator.py` does: in-paper accession (`YP_*` for Conlon, `ODV*` for Jacques)
→ **NCBI efetch** (still serves the superseded accessions, and the FASTA description even carries the SACOL tag)
→ **DIAMOND reciprocal best hit** against the Kp/Ec proteomes. 1,943 sequences cached once under
`data/processed/other/degradability/activator/`. Cross-phylum, so RBH only, with `pident`/`coverage` retained.

### Still worth fetching, in priority order

| Priority | Dataset | Why | Route |
| --- | --- | --- | --- |
| **1** | **Flynn 2003** *Mol Cell* 11:671, Tables 1–2 | The *E. coli* ClpXP selection bar is currently carried by a 45-row hand-curated stand-in. Table 2 also holds the five authoritative CM/NM consensuses. **Main-text tables, not a paywalled supplement.** | browser session → `--stage` |
| **2** | **Neher 2006** *Mol Cell* 22:193 | Second E. coli ClpXP trap set (DNA-damage condition); roughly doubles the trap positives. | browser session → `--stage` |
| **3** | **Meltome Atlas** (Jarzab 2020, *Nat Methods* 17:495) | The only proteome-scale conformational-stability signal, and stability is a first-class Regime B feature. **Verify E. coli coverage in the actual file before relying on it.** | open access |
| — | Feng 2013 (ACS), Ziemski 2021 (Wiley) | **Not worth further effort.** Both gated with no PMC record, both ClpC-family (a machine neither organism has), and Ziemski is the weakest evidence type in the pool (interaction screen). | — |

### The one experiment that would replace all of this

**Activator-treated proteomics in E. coli or K. pneumoniae** — ADEP4, ONC212 or a ClpPEL, treated vs untreated,
*with the partially/non-tryptic peptide analysis included* so the cleavage readout survives. Today the entire
partnerless bar rests on two *S. aureus* datasets transferred across a phylum boundary at 10.6% / 13.8%
coverage. One experiment in the right organism would replace that with direct evidence at full proteome
coverage. The Ersilia SoW already anticipates asking the consortium for exactly this.

---

## A‴. In hand as of 2026-08-06 — and route corrections

### `data/raw/other/degradability/literature/` — three previously-gated papers now on disk

| File | Is | Previously listed as |
|---|---|---|
| `PIIS1097276503000601.pdf` | **Flynn 2003** *Mol Cell* 11:671, "Proteomic Discovery of Cellular Substrates of the ClpXP Protease Reveals Five Classes of ClpX-Recognition Signals" (13 pp) | **priority 1** to fetch, gated |
| `PIIS1097276506001687.pdf` | **Neher 2006** *Mol Cell* 22:193, "Proteomic Profiling of ClpXP Substrates after DNA Damage…" (12 pp) | **priority 2** to fetch, gated |
| `febs15335-sup-0001-TableS1.pdf` (+ identical `.zip`) | **Ziemski 2021** *FEBS J* Table S1, Mtb ClpCP interaction screen (7 pp) | "not worth chasing" |

Flynn **60** + Neher **~100** proteins takes the E. coli ClpXP selection bar from a **45-row hand-curated
stand-in** to roughly **130–160 real proteins**. Extraction strategy, the recovered motif consensuses, and the
resulting correction to `src/degradability.py` are in **`docs/degradability_datasets.md` §3.1 and §5** — read
those before parsing. Short version: a stdlib `zlib` pass extracts the text with no install, but table structure
is lost, so **vocabulary-match against `L.load_genes('ecoli')` rather than parsing table positions**.
`brew install poppler` would make this materially easier and is not currently installed.

Ziemski remains **deprioritised** (ClpC-family, capped at `CLPC_ONLY_CAP = 0.45`, describes a machine neither
organism encodes) — now recorded as *in hand*, not *to fetch*.

### Route corrections found while probing (2026-08-06)

- **MobiDB bulk: `HEAD` returns `405`, `GET` returns `200`.** Do not conclude it is broken from a HEAD probe.
  Verified: `mobidb.org/api/download?proteome=UP000007841&format=tsv` → 6.31 MB / 85,234 rows / **5,728
  accessions**; `UP000000625` → 5.15 MB / 73,479 rows / **4,404 accessions**. This single file carries 8 disorder
  predictors + a consensus, AlphaFold pLDDT, **Pfam/Gene3D domain boundaries**, LIP regions, low-complexity, TM
  and signal-peptide tracks, each with `content_fraction`. **It is the cheapest high-value fetch in this axis.**
- **PHDB: `https://phdb.switchlab.org` fails; `http://` returns 200.** Alive but plain-HTTP only, web UI, and
  **no documented bulk export** — so treat the "4,305 proteins × >100 proteostatic parameters" shortcut as
  unavailable until someone confirms an export path.
- **Verified ✅** on this date: eSOL NBDC zip (`application/zip`); Gupta 2024 Zenodo `10895828` (one 51.5 MB zip)
  and its GitHub; Sen 2025 N-FIVE GitHub; STEPdb; `ecoliTPP.shiny.embl.de`; the Meltome browser; To 2021 via
  NSF-PAR (`application/pdf`); and all six PRIDE accessions `PXD042444` / `PXD022112` / `PXD009495` /
  `PXD016589` / `PXD022297` / `PXD062881` (titles match).
- **Still ⚠ unverified:** the Meltome Atlas *per-protein* supplementary file (identified by size, not content)
  and any PHDB bulk export.

---

## A⁗. 2026-08-06 (later) — two more papers in hand, one dataset turned out to be free, and the literature moved

### Location change

`docs/literature/` → **`data/raw/other/degradability/literature/`**. These are publisher PDFs (Cell Press,
Wiley, PNAS): ~9.7 MB of paywalled binaries that must **not** enter a public Git history. The new path is under
`data/`, which `.gitignore` line 3 excludes and **eosvc** tracks — consistent with the two-track rule in
`CLAUDE.md`. Verified with `git check-ignore` after the move. All doc references were repointed.

Contents (6 files):

| File | Paper |
|---|---|
| `PIIS1097276503000601.pdf` | Flynn 2003, *Mol Cell* 11:671 — ClpXP trap, 60 proteins + the five CM/NM consensuses |
| `PIIS1097276506001687.pdf` | Neher 2006, *Mol Cell* 22:193 — ClpXP trap ± DNA damage, ~100 proteins |
| `niwa-et-al-2012-…-cell-free-translation-system.pdf` | Niwa 2012, *PNAS* 109:8937 — chaperone effects |
| `1-s2.0-S2211124711000179-main.pdf` | Calloni 2012, *Cell Rep* 1:251 — DnaK interactome |
| `febs15335-sup-0001-TableS1.pdf` (+ `.zip`) | Ziemski 2021, *FEBS J* Table S1 — Mtb ClpCP screen |

### ⚠ Article PDFs are not data — check before planning around one

Both papers added in this batch are **article text only**. Measured by extracting their text and counting
identifiers: Niwa 2012 has **0 b-numbers / 8 gene-like symbols**; Calloni 2012 has **0 b-numbers / 7**. Neither
carries its per-protein table:

- Niwa 2012's 788-protein × 4-condition matrix is **Dataset S1 / Table S1**, a separate PNAS supplementary file.
- Calloni 2012's 674 DnaK interactors are **Table S1**, a separate Cell Press `mmc` file.

**Rule of thumb for this axis: assume the numbers are in a supplement, and verify by counting identifiers in the
extracted text before treating a PDF as ingestible.**

### ✅ …but Niwa 2012 needs no supplementary at all — eSOL already contains it

The paper says the data is "freely accessible" via the eSol database, and that is literally true: the
**NBDC eSOL CSV** (`dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip`, 197 KB, no login) carries
`Minus Sol (%)`, `TF Sol (%)`, `GroE Sol (%)` and `KJE Sol (%)` — chaperone-free, trigger-factor, GroEL/ES and
DnaK/DnaJ/GrpE solubility — for **exactly 788 proteins**, matching the paper's own count, plus matching yield
columns. 4,132 rows total, keyed by `JW_ID` / `B number` / gene name. So **Niwa 2009 (eSOL) and Niwa 2012
(chaperone dependence) are one download**, and the PNAS `DCSupplemental` route (which answers 403 anyway) is
unnecessary.

### Still genuinely blocked

| Paper | Status | Priority |
|---|---|---|
| ~~**Calloni 2012 Table S1**~~ | ✅ **RESOLVED — supplied manually 2026-08-06.** Filed at `data/raw/ecoli/degradability/calloni2012_dnak/` with a `SOURCE.md`. ScienceDirect `mmc*` still 403; there is no automated route, so treat it as manual-only | done |
| **Feng 2013** *J Proteome Res* 12:547 | ACS, **403**, no PMC record | **Skip** — ClpC-family, capped at 0.45, describes a machine neither organism encodes |
| **Meltome Atlas** per-protein T_m file | Browser 200 but the *specific* file was never confirmed | **Low** — Mateus 2018 covers E. coli; Meltome is only a cross-species sanity check. Needs a human to identify the right file |
| **PHDB** bulk export | `http://` only, no documented export | **Skip** — no route exists |

### Everything else in the minimum-viable-12 is automatable (verified 2026-08-06)

Probed and confirmed HTTP 200 with a real payload: **Gupta 2024** Table S1 (`media.springernature.com` MOESM4,
1.45 MB) · **Cappelletti 2021 LiP-MS** (Europe PMC `PMC7836100`, 34 MB — *in the OA subset despite being Cell
Press*) · **Mateus 2018** (`PMC6056769`, 44 MB) · **Mateus 2020** (Springer MOESM5/6/7) · **To 2021**
(`PMC8382223`, 1.34 MB — *OA despite being ACS*) · **To 2022** (`PMC9704704`) · **MacKrell 2026**
(`PMC12974527`, 18 MB) · **Schmidt 2016** (`PMC4888949`, 24 MB) · **Györkei 2022** (`PMC9023497`) ·
**Niwa 2022** (`PMC9228906`) · **eSOL** · **MobiDB bulk** (both proteomes).

**Lesson worth keeping:** "the journal is paywalled" does not imply "the data is unreachable". Three of the
above (Cappelletti, To 2021, and Conlon 2013 earlier) are gated at the article but open at the data. Always probe
Europe PMC `supplementaryFiles` and `media.springernature.com` before declaring something blocked.

### Filing convention for manually-supplied papers — settled 2026-08-06

Two distinct destinations, and the split matters because one tree is public-safe and the other is not:

| What | Where | Why |
|---|---|---|
| **Data** (supplementary tables, per-protein numbers) | `data/raw/<organism>/degradability/<key>/` for organism-native, `data/raw/other/degradability/<key>/` for cross-species | Matches the `10a` manifest's `organism` field, so the fetch ladder finds it as `cached` and never re-downloads. One directory per dataset, `<firstauthor><year>_<slug>`, plus a `SOURCE.md` |
| **Article PDFs** (the papers themselves) | `data/raw/other/degradability/literature/` | A shared pool — they are provenance, not data, and are not per-dataset. Gitignored, eosvc-tracked, so paywalled binaries stay out of the public Git history |

Worked example, this batch: Calloni 2012's **workbook** went to
`data/raw/ecoli/degradability/calloni2012_dnak/calloni2012_TableS1-S8_dnak_interactome.xls` (E. coli-native, so the
`ecoli` tree — alongside `nagar2021_halflife`, `flynn2003_clpxp_trap`, `neher2006_clpxp_sos`), while its
**article PDF** stayed in `literature/`. It was registered in the `10a` manifest as
`Dataset("calloni2012_dnak", "ecoli", …, open_access=False)`, and `--only calloni2012_dnak` now reports
`ok (cached)` without attempting a fetch. `SOURCE.md` is invisible to `_has_data()` (`.md` is not in `DATA_EXTS`),
so it cannot be mistaken for the dataset.

**New dependency:** `xlrd>=2.0.1`, pinned in `requirements.txt`. Calloni's workbook is legacy BIFF `.xls`, which
`openpyxl` cannot open — a real trap, because pandas' error message points at the missing package rather than the
format. Confirmed genuine BIFF from the magic bytes, and the file's OLE metadata (`author: tchen`,
`last saved by: Hayer-Hartl, Manajit`) independently confirms provenance.

---

## A⁵. Fetch run of 2026-08-06 — what landed, and three corrections

Ran via the `10a` manifest (not ad-hoc `curl`), so it is reproducible: `python scripts/10a_fetch_degradability.py`.
**172 MB landed across 13 datasets.** Row counts below were measured from the fetched files.

| Dataset | Route | Files | MB | Verified content |
|---|---|---|---|---|
| `gupta2024_turnover` | Springer CDN | 2 | 1.7 | `TableS1` **3,263 rows × 42 cols**; header at row **3**, short names at row 4; keys `sp\|A5A614\|YCIZ_ECOLI` |
| `cappelletti2021_lipms` | Europe PMC | 4 | 60.6 | `mmc2.xls` 15 sheets (+ mmc3, mmc4) |
| `mateus2018_tpp` | Europe PMC | 14 | 43.5 | `MSB-14-e8242-s002.xlsx` = **"Dataset EV1 E. coli"** (the T_m table) |
| `mateus2020_tpp121` | Springer CDN | 2 | 0.4 | Supp Data 5 (MOESM5); MOESM6 is 60 MB and not needed |
| `mackrell2026_boncat` | Europe PMC | 6 | 17.5 | `sd01/02/03.xlsx`; sd02 = 1,959 × 23 |
| `schmidt2016_abundance` | Europe PMC | 3 | 23.9 | `NIHMS65833-supplement-Supplementary_tables.xlsx`, **30 sheets** |
| `gyorkei2022_solubility` | Europe PMC | 12 | 2.8 | 11 MOESM xlsx + docx |
| `niwa2022_lon_clpxp` | Europe PMC (**nested**) | 6 | 5.0 | `Supplementary_Dataset_S1/S2/S3.xlsx`; S2 = 2,292 × 106 |
| `esol_solubility` | NBDC | 1 | 0.7 | `esol.csv` **4,132 × 28**; 3,173 `Solubility (%)`, **788** with all four chaperone columns |
| `nagar2021_halflife` | cached | 4 | 3.6 | `sd002.xlsx` = 1,150 × **190** (the 188-feature matrix) |
| `calloni2012_dnak` | manual | 2 | 1.4 | 8 sheets, Table S1 = 698 rows |
| `mobidb` (Kp) | GET | 1 | 6.3 | **85,234 rows / 5,728 accessions** |
| `mobidb` (Ec) | GET | 1 | 5.2 | **73,479 rows / 4,404 accessions** |

### ⚠ Correction 1 — Gupta's KO panel has no Δ*ftsH*, and it is 13 conditions not 14

Measured from the fetched `TableS1`: the knockout conditions are **Δ*clpP*, Δ*lon*, Δ*hslV*,
Δ*clpP*Δ*lon*Δ*hslV*, Δ*smpB*** — all N-lim chemostat at 6 h doubling. **`ftsH` appears nowhere in the sheet**
(both sheets searched). Earlier text in this file, `degradability_references.md` §3.1, `degradability_datasets.md`
§1 and `03_degradability.md` all listed Δ*ftsH*; all four are now corrected. Counted 13 `Average of half-lives`
columns and **3,263** data rows, settling the 3,262-vs-3,200 / 13-vs-14 discrepancy that had been carried in two
places. Consequence: **FtsH substrate attribution is not available from Gupta** — for that channel the only
sources are Westphal 2012 and Arends 2016, both Regime A and therefore low priority.

### ⚠ Correction 2 — `to2021_refoldability` had the WRONG PMCID, and To 2021 is not open access

The entry was first written with `PMC8382223`, and my earlier claim that "To 2021 is in the PMC OA subset despite
being ACS" rested on a 200 from that ID. **`PMC8382223` is a different paper** — *"Evaluation of the Higher Order
Structure of Biotherapeutics Embedded in Hydrogels for Bioprinting and Drug Release"*, Anal Chem. It fetched
`ac1c01850_si_001.pdf`, which was deleted.

The real record is **`PMC8650709`**, `isOpenAccess: N`. The machine-readable tables are ACS-gated: bioRxiv
`10.1101/2020.08.28.273110` (both v1 and v2) ships only `media-1.pdf`, and NSF-PAR `purl/10311595` is the
accepted-manuscript PDF. The manifest now fetches the **open bioRxiv SI PDF** and says plainly that the
per-protein table still needs a browser session or PDF text extraction.

**Generalisable lesson, now also a comment in `10a`: a 200 from the wrong PMCID is indistinguishable from a 200
from the right one.** Verify a PMCID by its article title. `_fetch_fulltext()` caches the article XML for
open-access datasets, which makes this checkable after the fact — a title sweep over those cached XMLs found this
error, and confirmed the other 10 were correct.

### ⚠ Correction 3 — To 2022's Europe PMC zip is figures only

`PMC9704704` returns HTTP 200 and a well-formed 380 kB zip containing **12 jpg/gif figure files and no data**.
So "the endpoint worked" and "we got the dataset" are different claims. Left as a placeholder; low priority
because eSOL already supplies quantitative chaperone dependence for 788 proteins.

### Two ladder bugs fixed by this run

1. **Nested archives were silently dropped.** MDPI ships the whole supplement as an inner `…-s001.zip`, and
   since `.zip` is not in `DATA_EXTS` the outer pass reported *"archive held no data members"* for a bundle that
   contained everything. `_extract_data_members()` now recurses **one level** (depth-capped so a zip bomb cannot
   spin). This is what unblocked Niwa 2022's three Supplementary Datasets.
2. **PDF-only and text-only payloads were rejected.** The accept branch handled spreadsheets and archives only,
   so MobiDB's TSV and To 2021's SI PDF both bounced. Added `_wants_raw()` / `RAW_EXTS` — accepted **only** when
   the manifest opts in by naming the extension in `filename`, with `_is_html_shell()` still guarding, so an
   error page still cannot be banked as data.

### Also

- **MobiDB is fetched by `_fetch_mobidb()`, not from `MANIFEST`.** Two redundant `Dataset` entries were removed
  in favour of that route, which validates harder than the generic ladder can: it requires the header to start
  `acc` **and** ≥1,000 rows, because a wrong proteome id returns a header-only 61-byte HTTP 200.
- **RegulonDB now fails** (`0/2 files`, both `network_tf_gene.txt` and `network_sigma_gene.txt`). Annotation-only
  and weighted 0, so it does not gate anything — but the site appears to have changed again.
- **New dependency:** `xlrd>=2.0.1` (pinned) for Calloni's legacy BIFF workbook.
