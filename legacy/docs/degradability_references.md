# Degradability axis — references

Consolidated list of every paper, dataset, tool and database used (or explicitly considered and rejected) in
the **degradability** axis of the GraDi target-prioritization pipeline (`docs/03_degradability.md`, planned
`scripts/10*`). Sections: **1** the degrader modality itself, **2** the recognition rules, **3** *E. coli*
proteome-scale datasets (the transfer layer), **4** *K. pneumoniae*-native resources, **5** reviews, **6**
tools & methods, **7** rejected / not obtained.

Notes:

- **Nothing in this list is ingested yet.** The axis has no pipeline; this file is the sourcing plan. See
  `docs/degradability_report.md` for per-dataset access routes and coverage.
- The Kp anchor strain is **HS11286** (UniProt `UP000007841`, 5,728 proteins). Because no *K. pneumoniae*
  turnover, meltome or LiP-MS dataset exists, §3 is *E. coli* K-12 and everything transfers onto HS11286 by
  ortholog — with a measured ceiling of **3,179 / 5,728 (55.5%)**.
- Every DOI, PMCID and PRIDE accession below was resolved programmatically (Crossref / NCBI E-utilities) and
  every download URL was probed (HTTP 200/206) on **2026-08-04**.
- Each entry notes the track (`3.x`) that consumes it.

---

## 1. The degrader modality

- **Morreale FE, Kleine S, Leodolter J, et al. 2022** — the original **BacPROTACs**: bivalent degraders that
  bind the ClpC N-terminal domain (cyclomarin-A-derived ligand) and recruit neo-substrates to ClpCP, which the
  binding event simultaneously activates. *Cell* 185:2338–2353.e18 ·
  DOI [10.1016/j.cell.2022.05.009](https://doi.org/10.1016/j.cell.2022.05.009) · PMC9240326.
  **Foundational but not transferable**: ClpC is Firmicute/Actinobacteria-restricted and absent from Kp — see
  §3.0 of the spec. Context for the whole axis.

- **Izert-Nowakowska MA, Klimecka MM, Antosiewicz A, et al. 2025** — **CLIPPERs** (Clp-Interacting Peptidic
  Protein Erasers): bait–linker–anchor peptides where the anchor is the **SspB C-terminal "XB" motif** binding
  the ClpX zinc-binding domain. Degraded endogenous, **untagged, degronless** GroEL in live *E. coli*
  (~17% at 1 h, ~40% at 6 h, 42 °C); abolished in Δ*clpX* and Δ*clpP*. Also the source of the
  substrate-requirement checklist (accessible terminus, monomeric/unassembled, unstructured, low copy number,
  low stability) inferred from the prototypes that failed. *EMBO Rep* 26(16):3994–4016 ·
  DOI [10.1038/s44319-025-00510-9](https://doi.org/10.1038/s44319-025-00510-9) · PMID 40562793 · PMC12373786 ·
  PRIDE **PXD045730** · EMDB EMD-19687 · PDB 8S32. **The paper that licenses this axis for a Gram-negative.**
  Tracks 3.0, 3.1, 3.5, 3.6, 3.7.

- **Taylor EK, Barbosa PS, Kadambi T, et al. 2026** — **DegP molecular glues**: tazobactam itself accelerates
  DegP-mediated degradation of TEM β-lactamases, enhanced by *linkerless* incorporation of dipeptide motifs
  enriched among DegP substrates; improved piperacillin synergy plus oral-dosing PK. A glue, not a chimera —
  small, no permeability penalty, and acting in the compartment the target occupies. **bioRxiv preprint** ·
  DOI [10.64898/2026.03.30.715243](https://doi.org/10.64898/2026.03.30.715243). Basis of the periplasmic
  channel. Track 3.7.

- **Nie Q, Wu JW, Zhang K, et al. 2025** — NacssrA-1: nacubactam covalently linked to the full ssrA degron
  `AANDENYALAA`, reported to degrade CTX-M-14 and resensitise to cefotaxime. *Chem Commun* 61:13149–13152 ·
  DOI [10.1039/D5CC02015H](https://doi.org/10.1039/d5cc02015h). **Cited as a cautionary counter-example, not
  as support**: the molecule is 1,669 Da and CTX-M β-lactamases are periplasmic while ClpXP is cytoplasmic —
  the target had been engineered to remain in the cytoplasm. Motivates compartment matching as a hard filter.

- **Davis JH, Baker TA, Sauer RT 2011** — small-molecule control of degradation using **split adaptors**: SspB
  divided into core-domain–FRB (recognises a DAS+4-tagged substrate) and XB-tail–FKBP12 (binds ClpX), assembled
  by rapamycin. The *E. coli* precedent for chemically addressing the ClpX/SspB axis, and the design the
  mycobacterial screen (§2) later reused. *ACS Chem Biol* 6(11):1205–13 ·
  DOI [10.1021/cb2001389](https://doi.org/10.1021/cb2001389) · PMID 21866931.

- **Brötz-Oesterhelt H, Vorbach A 2021** — ADEP reprogramming of ClpP: the drug binds the apical H-pockets,
  widens the axial pore and blocks Clp-ATPase docking, converting ClpP into an unregulated protease.
  *Front Mol Biosci* 8:690902 ·
  DOI [10.3389/fmolb.2021.690902](https://doi.org/10.3389/fmolb.2021.690902).
  **Not a modality for this axis** (it is proteome-wide, not targeted), but important *orthogonal* mechanistic
  evidence about which physical states are degradable — and the cautionary permeability datapoint (isolated
  *E. coli* ClpP is exquisitely ADEP-sensitive, whole-cell MIC >64 µg/mL, because ADEP is an AcrAB-TolC
  substrate).

- **Silber N, Pan S, Schäkermann S, et al. 2020** — ADEP-activated ClpP unfolds FtsZ specifically via its
  **conformationally flexible N-terminal domain**. *mBio* 11(3) ·
  DOI [10.1128/mBio.01006-20](https://doi.org/10.1128/mBio.01006-20) · PMID 32605984 · PMC7327170.
  Independent support for the N-terminal-flexibility argument in track 3.1.

---

## 2. Recognition rules (tracks 3.1, 3.2, 3.7)

### 2.1 The quantitative degradability model

- **Won HI, Zinga S, Kandror O, et al. 2024** — **the only validated quantitative bacterial degradability
  predictor.** 72 native *M. smegmatis* proteins screened for rapamycin-induced (FRB/FKBP) degradation by
  ClpC1P1P2; Lasso over 485 descriptors fitted on 54, held-out r > 0.6 (P = 0.0057) on 18. Top feature =
  **mean disorder propensity of the N-terminal ~30 residues** (flDPnn, threshold 0.2); the C-terminal 30 aa
  showed **no** correlation. Steady-state abundance anti-correlates with degradation rate at **r = −0.69,
  P = 6.6 × 10⁻⁹**. Windows tested: N-/C-terminal 15 and 30 aa plus full length. Most native proteins tested
  *were* degradable. *Nat Commun* 15:4065 ·
  DOI [10.1038/s41467-024-48506-8](https://doi.org/10.1038/s41467-024-48506-8) · PMID 38744895 · PMC11094019 ·
  code + data Zenodo [10.5281/zenodo.11004395](https://doi.org/10.5281/zenodo.11004395) ·
  Supplementary Data 1 (348-protein feature matrix), 2 (72 targets + degradation constants), 3 (predicted
  scores). Supplies the window scheme for 3.1 and the abundance term for 3.5.
  *Transfer caveat*: ClpC1 in an Actinobacterium — the N-vs-C asymmetry may not carry to ClpX, which reads both
  termini and whose best-characterised degron is C-terminal.

### 2.2 C-terminal degrons

- **Beardslee PC, Schmitz KR 2024** — **DEtox** (Degron Enrichment by Toxin): VapC bearing randomized
  C-terminal pentapeptides, ~100,000 unique tags of 3.2 × 10⁶ theoretical mapped by NGS, run in WT, Δ*clpX* and
  Δ*clpP*. **~1% of all random 5-mers are functional degrons.** Consensus (−5→−1)
  `[L/F/Y/W]-[L/A/R]-L-A-A`; terminal **Ala-Ala in 91%** of the top 100; Leu favoured at −3, aromatic at −5;
  depleted: polar@−5, bulky/β-branched@−4, Pro/His/Gly@−3. Best single tag `FKLVA`, 639-fold enriched.
  **ClpXP overwhelmingly dominates** short-C-degron recognition; Lon and FtsH contribute little. *eLife*
  reviewed preprint · DOI [10.7554/eLife.98528.1](https://doi.org/10.7554/eLife.98528.1) · PMC10862746 ·
  bioRxiv 2024.01.29.576913. **The empirical PSSM that replaces the hand-written ssrA regex.** Track 3.2a.

- **Cragan M, Puri N, Karzai AW 2025** — substrate recognition and cleavage-site preferences of **Lon**.
  Establishes a widely-distributed C-terminal Lon degron consensus **`x–[L/I]–[L/I/V]–H-COOH`**, invariantly
  ending in His; His@−1 strongly preferred, Asp/Arg/Lys@−1 strongly disfavoured (H→D drops degradation to
  ≤0.3 min⁻¹·hexamer⁻¹); affinities span K_D 0.24 → >50 µM. Cleavage-site P1 preference Phe ≫ Leu ≈ Ala,
  products 7–35 aa. **Conserved across *E. coli*, *Yersinia pestis* and *Mycoplasma pneumoniae*** — which is
  what licenses transfer to *Klebsiella*. *J Biol Chem* 301(4):108365 ·
  DOI [10.1016/j.jbc.2025.108365](https://doi.org/10.1016/j.jbc.2025.108365) · PMID 40023398 · PMC11986505.
  **A rule the previous spec lacked entirely.** Track 3.2b.

- **Gur E, Sauer RT 2008** — recognition of misfolded proteins by Lon: the positive determinant is a **cluster
  of aromatic side chains that is buried in the native fold**; the secondary determinant is the *absence* of
  small polar residues. Lon is therefore not a generic "any unfolded protein" protease. *Genes Dev*
  22(16):2267–77 · DOI [10.1101/gad.1670908](https://doi.org/10.1101/gad.1670908) · PMID 18708584.
  Track 3.2b / 3.4.

- **Gur E, Sauer RT 2009** — degrons program the speed and efficiency of Lon. The canonical **sul20C** degron
  (SulA C-terminal 20 aa, `ASSHATRQLSGLKIHSNLYH` — itself an instance of the 2025 consensus) is ~9× faster than
  the β20 degron at saturation; degradation is positively cooperative (Hill ~2). Supplies the loss-of-function
  controls (sul20C-YA, sul20C-HA). *PNAS* 106(44):18503–8 ·
  DOI [10.1073/pnas.0910392106](https://doi.org/10.1073/pnas.0910392106) · PMID 19841274.

- **Flynn JM, Neher SB, Kim YI, Sauer RT, Baker TA 2003** — the classic ClpP-trap census: >50 trapped *E. coli*
  substrates and the **five ClpX-recognition signal classes** (C-motif 1 `LAA-COOH` ssrA-type; C-motif 2
  `RRKKAI-COOH` MuA-type; N-motif 1 `polar–T/φ–φ–basic–φ`; N-motif 2 `NH₂-Met–basic–φ–φ–φ`; N-motif 3
  `φ–x–polar–x–polar–x–basic–polar`). 45% of ClpX-dependent trapped proteins had ssrA- or MuA-like C-termini.
  *Mol Cell* 11(3):671–83 ·
  DOI [10.1016/S1097-2765(03)00060-1](https://doi.org/10.1016/s1097-2765(03)00060-1) · PMID 12667450.
  **NOT OBTAINED — see §7.** The motif strings above are quoted from secondary sources (PMC2394798); the
  positional anchoring and the trapped-substrate table are behind the paywall. Track 3.2d, at low weight and
  explicitly flagged.

### 2.3 N-terminal degrons

- **Sen S, Corrado N, Tiso AR, Kho KK, Kunjapur AM 2025** — combinatorial mutagenesis of N-terminal sequences
  (~2.2 M variants, FACS + NGS, protein-stability-index readout) **expands the *E. coli* N-degron rule from
  position 1 to positions 1–5**. P1 destabilizing potency `F≈R < L < W≈K≈Y`. **Pro@P2 rescues (+1.05 PSI,
  effect size 0.71)** and **Gly@P2 rescues (+0.63)**; **Gln@P2 destabilizes as strongly as Leu@P1 (−0.51)**;
  clustered bulky FLWY over P2–P5 destabilizes (mean −0.56); **net negative charge over P2–P5 stabilizes** and
  at charge −4 rescues even FLYWRK P1 substrates (ClpS's pocket is flanked by Asp36/Glu41); multiple Gly/Ser in
  P2–P5 cumulatively stabilize. Novel non-canonical degrons found with Cys, Gln or His at P1 in permissive
  context. All P1–P5 effects vanish in **Δ*clpS***; the Arg/Lys effect specifically vanishes in **Δ*aat***.
  Model **"N-FIVE"** + code at
  [github.com/KunjapurLab/N-terminal-cluster-stability](https://github.com/KunjapurLab/N-terminal-cluster-stability).
  **bioRxiv preprint** 2025.05.22.655665v2 ·
  DOI [10.1101/2025.05.22.655665](https://doi.org/10.1101/2025.05.22.655665) · PMID 40661368 · PMC12258705.
  **Replaces `nterm_destabilizing` and explains why the P1-only rule fails.** Track 3.2c.

- **Humbard MA, Surkov S, De Donatis GM, Jenkins LM, Maurizi MR 2013** — **the N-degradome of *E. coli***.
  Immobilised wild-type ClpS affinity column on native lysate, specific elution with an N-degron peptide;
  >100 proteins differentially bound, 32 of 37 sequenced N-terminal peptides bore an N-degron, ~20–30%
  Aat-dependent. **Central finding for this axis: most N-degron-bearing proteins had been N-terminally
  truncated by endoproteases first** — the degron is generated post-translationally, so scanning annotated
  ORF termini systematically under-counts substrates. Table 1 lists the measured neo-N-termini.
  *J Biol Chem* 288(40):28913–28924 ·
  DOI [10.1074/jbc.M113.492108](https://doi.org/10.1074/jbc.M113.492108) · PMC3789986
  (**not in the PMC OA subset** — use the JBC PDF at
  `https://www.jbc.org/article/S0021-9258(20)48862-3/pdf`). Tracks 3.2c, 3.3; and the "latent degron" caveat.

- **Schuenemann VJ, Kralik SM, Albrecht R, et al. 2009** — structural basis of N-end rule substrate
  recognition by ClpS: the N-degron side chain inserts into a pre-formed deep hydrophobic cleft (~200 Å³,
  sized for Leu/Phe/Tyr, expandable for Trp). *EMBO Rep* 10(5):508–14 ·
  DOI [10.1038/embor.2009.62](https://doi.org/10.1038/embor.2009.62) · PMID 19373253.

- **Dougan DA, Truscott KN, Zeth K 2010** — review of the bacterial N-end rule pathway; confirms Leu/Phe/Trp/Tyr
  as the primary destabilizing residues, Arg/Lys as secondary (via Aat), and that bacteria have **no**
  ClpS-independent branch analogous to eukaryotic UBR. *Mol Microbiol* 76:545–558 ·
  DOI [10.1111/j.1365-2958.2010.07120.x](https://doi.org/10.1111/j.1365-2958.2010.07120.x).

- **Hirel PH, Schmitter MJ, Dessen P, et al. 1989** — the classical **Met-excision rule**: the extent of
  N-terminal Met removal decreases monotonically with the side-chain length of the penultimate residue.
  Operationally, Met is removed when residue 2 ∈ {Ala, Cys, Gly, Pro, Ser, Thr, Val}. *PNAS* 86(21):8247–51 ·
  DOI [10.1073/pnas.86.21.8247](https://doi.org/10.1073/pnas.86.21.8247) · PMID 2682640. Gate for track 3.2c.

- **Frottin F, Martinez A, Peynot P, et al. 2006** — refines the above against 862 *E. coli* proteins with
  experimentally known in-vivo N-termini: **Val and Thr at P1′ are cleaved much less efficiently** than
  Ala/Cys/Gly/Pro/Ser (some never), and P2′–P4′ identity can strongly slow the reaction — so this is not a pure
  position-2 lookup. *Mol Cell Proteomics* 5(12):2336–49 ·
  DOI [10.1074/mcp.M600225-MCP200](https://doi.org/10.1074/mcp.M600225-MCP200) · PMID 16963780. Track 3.2c.

### 2.4 Initiation regions, engagement mechanics, internal degrons

- **Fei X, Bell TA, Barkow SR, et al. 2020** — cryo-EM basis of ClpXP recognition and unfolding of ssrA-tagged
  substrates; the pore-1/pore-2/RKH loop architecture and the closed→open axial-channel transition.
  *eLife* 9:e61496 · DOI [10.7554/eLife.61496](https://doi.org/10.7554/eLife.61496) · PMID 33089779 · PMC7652416.

- **Saunders RA, Stinson BM, Baker TA, Sauer RT 2020** — **the initiation-region length numbers**: a "short"
  (~5-residue) degron forms a recognition complex contacting only two pore-1 loops in the **closed**-channel
  state, from which unfolding can proceed directly; a "long" (~20-residue) degron requires progression to the
  **open**-channel engaged complex gripping with five pore-1 loops. *PNAS* 117:28005–28013 ·
  DOI [10.1073/pnas.2010804117](https://doi.org/10.1073/pnas.2010804117). Track 3.1.

- **Kenniston JA, Baker TA, Sauer RT 2005** — partitioning between unfolding and release determines substrate
  selectivity and partial processing. Establishes the **~37-residue ClpX-pore-to-ClpP-active-site reach** (from
  the tail lengths of partial products) and that **longer unstructured tails promote retention during futile
  unfolding attempts** — i.e. tail length modulates commitment, not just binding. *PNAS* 102:1390–1395 ·
  DOI [10.1073/pnas.0409634102](https://doi.org/10.1073/pnas.0409634102) · PMC547888. Track 3.1.

- **Lyu Y, Bolstad I, Davis JH, Sauer RT, Ghanbarpour A 2026** — an open axial channel enhances degradation of
  **specific degron classes**: a pore-2-loop-deletion ClpX variant degrades N-motif-1/2/3 and C-motif-2
  substrates better than wild type, while **C-motif-1 (ssrA) is preferentially handled by the closed channel**.
  Implication for this axis: the N-motifs and C-motif-2 are the **natural** substrate route; ssrA is the
  engineered special case. *Protein Sci* 35(1):e70372 ·
  DOI [10.1002/pro.70372](https://doi.org/10.1002/pro.70372) · PMID 41432315. Tracks 3.2a, 3.2d.

- **Hoskins JR, Yanagihara K, Mizuuchi K, Wickner S 2002** — ClpAP and ClpXP degrade proteins with tags
  located in the **interior** of the primary sequence. *PNAS* 99:11037–11042 ·
  DOI [10.1073/pnas.172378899](https://doi.org/10.1073/pnas.172378899). Track 3.1 (internal degrons).

- **Camberg JL, Viola MG, Rea L, et al. 2014** — *E. coli* FtsZ has **two independent, additive** ClpXP sites:
  a C-terminal one (residues 366–383; critical Arg379/Lys380) and an **internal one in the disordered linker
  (residues 349–358, `QEQKPVAKVV`)**, ~30 residues upstream of the C-terminus. *PLoS One* 9(4):e94964 ·
  DOI [10.1371/journal.pone.0094964](https://doi.org/10.1371/journal.pone.0094964) · PMID 24722340 · PMC3983244.
  Track 3.1; also the ground-truth spot check for the internal-degron feature.

- **Baker TA, Sauer RT 2012** — ClpXP as an ATP-powered unfolding and degradation machine; the standard review
  of degron classes, adaptors and unfolding mechanics. *Biochim Biophys Acta* 1823(1):15–28 ·
  DOI [10.1016/j.bbamcr.2011.06.007](https://doi.org/10.1016/j.bbamcr.2011.06.007) · PMID 21736903.

### 2.5 The membrane channel (FtsH)

- **Chai-Danino M, Ravensary-Modin N, Vladimirov VI, et al. 2026** — **membrane-embedded polar residues target
  membrane proteins for degradation by FtsH.** Lipid-facing polar residues inside the transmembrane region —
  buried in the protein core when correctly folded — signal misfolding; recognition requires FtsH's own TM
  domain. Crucially, such residues can trigger degradation of a **folded** protein and **do not require** the
  extended cytosolic tail otherwise needed. *Nat Commun* 17:3067 ·
  DOI [10.1038/s41467-026-69829-8](https://doi.org/10.1038/s41467-026-69829-8) · PMID 41730867 ·
  preprint bioRxiv 2023.12.12.571171. **Opens the inner-membrane compartment (1,320 Kp proteins).** Track 3.7.

- **Chiba S, Akiyama Y, Ito K 2002** — membrane-protein degradation by FtsH can be initiated from **either
  end**, length-dependently and largely sequence-independently: **>~20 residues** of cytoplasmically exposed
  tail at the N-terminus, or **~10** at the C-terminus. *J Bacteriol* 184(17):4775–82 ·
  DOI [10.1128/JB.184.17.4775-4782.2002](https://doi.org/10.1128/JB.184.17.4775-4782.2002) · PMID 12169602.
  Track 3.7 thresholds.

---

## 3. *E. coli* K-12 proteome-scale datasets (the transfer layer)

### 3.1 Turnover, half-life, protease attribution (track 3.3)

- **Gupta M, Johnson ANT, Cruz ER, et al. 2024** — global protein turnover: heavy-isotope pulse + TMTproC over
  **~3,200 proteins × 13 conditions** (strain NCM3722), plus a **protease-knockout panel** (Δ*clpP*, Δ*lon*,
  Δ*hslV*, triple, Δ*smpB* — **not** Δ*ftsH*, see §11/§12 corrections) assigning substrates to six categories — predominantly ClpP (64), Lon (14),
  HslV (1, UhpA), additive (82), redundant (41), still degrading in the triple KO (~100). Findings that shape
  this axis: **disordered proteins have significantly shorter half-lives (p ≈ 1e-208)**; small (<10 kDa)
  proteins, transcriptional regulators and Fe–S proteins are short-lived; **rapidly degrading proteins showed
  no enrichment for destabilizing N-terminal residues**; **cytoplasmic proteins are selectively degraded while
  membrane and periplasmic proteins are largely stable**; and **~40% of active degradation persists in the
  triple knockout** — a major proteolysis pathway remains unidentified. Supplementary Data 4 gives ~600
  measured in-vivo N-termini. *Nat Commun* 15:5890 ·
  DOI [10.1038/s41467-024-49920-8](https://doi.org/10.1038/s41467-024-49920-8) · PMID 39003262 · PMC11246515 ·
  PRIDE **PXD042444** · code [github.com/wuhrlab/ProteinTurnoverEcoli](https://github.com/wuhrlab/ProteinTurnoverEcoli)
  · Zenodo [10.5281/zenodo.10895828](https://doi.org/10.5281/zenodo.10895828).
  **The primary label set for the axis.** Replaces the hand-curated half-life stand-in.

- **Nagar N, Ecker N, Loewenthal G, et al. 2021** — pulsed-SILAC half-lives for **1,149** *E. coli* proteins
  (of 1,602 quantified): 741 stable, 334 slow, 72 fast. Data Set S2 is a ready-made **188-feature matrix**
  (4 disorder features, N-end rule, known C-/N-terminal degrons, pI, MW, aromaticity, GRAVY, plus 128 node2vec
  PPI-network features). LASSO AUC 0.72 for fast-vs-rest. **Two results this axis must respect: the
  PPI-network features were the most informative, and the known degrons — N-end rule and ClpXP motifs — showed
  no significant association with stability.** *mSystems* 6(1):e01296-20 ·
  DOI [10.1128/mSystems.01296-20](https://doi.org/10.1128/mSystems.01296-20) · PMID 33531410 · PMC7857536 ·
  PRIDE **PXD022112**. Track 3.3, plus the §Caveats and the PPI suggestion.

- **MacKrell EJ, Lomenick B, Qiu Y, et al. 2026** — the newest global degradation dataset: BONCAT
  (Aha/Anl-NLL-MetRS) + TMT + FAIMS, time-resolved, in **both exponential (1,810 proteins) and stationary
  (1,339) phase**; 88 and 56 pronouncedly unstable. Validated unstable proteins include PdeH, **ClpS itself**,
  and all four A-type Fe–S carriers (IscA, ErpA, NfuA, SufA); mutagenesis of PdeH's **N-terminal extension
  abolished ClpXP recognition**. Supplementary S4 gives ML predictions for unmeasured proteins.
  *PNAS* 123 · DOI [10.1073/pnas.2515265123](https://doi.org/10.1073/pnas.2515265123) · PMC12974527 ·
  PRIDE **PXD062881**. Best independent validation label set; adds stationary-phase coverage the SILAC studies
  miss.

- **Neher SB, Villén J, Oakes EC, et al. 2006** — quantitative (SILAC) ClpXP-trap proteomics ± DNA damage:
  **25% of the SOS-response proteome are ClpXP substrates.** Nine confirmed: LexA fragments, UmuD/UmuD′, RecN,
  UvrA, DinD, DinI, Fur, SulA, YebG. *Mol Cell* 22(2):193–204 ·
  DOI [10.1016/j.molcel.2006.03.007](https://doi.org/10.1016/j.molcel.2006.03.007) · PMID 16630889.
  An **open** trap set replacing part of the hand-curated Flynn stand-in. Track 3.3.

- **Niwa T, Chadani Y, Taguchi H 2022** — shotgun proteomics under GroE depletion, comparing **Lon vs ClpXP vs
  HslUV**: Lon is the dominant protease degrading the ~80 obligate GroE substrates. The best public
  "which protease eats which *E. coli* protein" side-by-side. *Molecules* 27(12):3772 ·
  DOI [10.3390/molecules27123772](https://doi.org/10.3390/molecules27123772) · PMID 35744894 · PMC9228906.
  Tracks 3.3, 3.4.

- **Westphal K, Langklotz S, Thomanek N, Narberhaus F 2012** — FtsH-trap approach: 15 putative substrates,
  4 validated novel (IscS, DadA, FdoH, YfgM). *J Biol Chem* 287(51):42962–71 ·
  DOI [10.1074/jbc.M112.388470](https://doi.org/10.1074/jbc.M112.388470) · PMID 23091052 · PMC3522291.

- **Arends J, Thomanek N, Kuhlmann K, et al. 2016** — in-vivo trapping of FtsH substrates by label-free
  quantitative proteomics: **>50 putative substrates** across exponential and stationary phase.
  *Proteomics* 16:3161–3172 · DOI [10.1002/pmic.201600316](https://doi.org/10.1002/pmic.201600316) ·
  PMID 27766750. Track 3.3 (FtsH channel labels).

### 3.2 Thermal stability (track 3.4)

- **Mateus A, Bobonis J, Kurzawa N, et al. 2018** — **thermal proteome profiling in bacteria**; the canonical
  *E. coli* Tm reference. 1,831 proteins identified, **1,738 with a fitted apparent Tm** (strain BW25113,
  intact cells + lysate). Dataset EV1 is the per-protein Tm table. Note for use: **Tm follows a
  cell-surface → cytoplasm high-to-low gradient**, so Tm partly encodes compartment and should be residualised
  on localization. *Mol Syst Biol* 14(7):e8242 ·
  DOI [10.15252/msb.20188242](https://doi.org/10.15252/msb.20188242) · PMID 29980614 · PMC6056769 (CC-BY) ·
  PRIDE **PXD009495**.

- **Mateus A, Hevler J, Bobonis J, et al. 2020** — the functional proteome landscape of *E. coli*: TPP across
  **121 strains** (mostly Keio deletions); Supplementary Data 5 gives abundance + thermal-stability scores for
  **1,764** proteins. The derived feature worth extracting is **stability-score variance across the 121
  perturbations** — a conformational-plasticity measure arguably more informative than a single Tm.
  *Nature* 588:473–478 · DOI [10.1038/s41586-020-3002-5](https://doi.org/10.1038/s41586-020-3002-5) ·
  PMID 33299184 · PMC7612278 · PRIDE **PXD016589** · code
  [github.com/fstein/EcoliTPP](https://github.com/fstein/EcoliTPP) · browser `ecoliTPP.shiny.embl.de`.

- **Jarzab A, Kurzawa N, Hopf T, et al. 2020** — **Meltome Atlas**: ~48,000 proteins, 13 species including
  *E. coli* and five other prokaryotes. Value here is a cross-species Tm prior to sanity-check ortholog
  transfer; for *E. coli* itself it largely duplicates Mateus 2018. *Nat Methods* 17:495–503 ·
  DOI [10.1038/s41592-020-0801-4](https://doi.org/10.1038/s41592-020-0801-4) · PMID 32284610 ·
  PRIDE **PXD011929** (partial) · browser `meltomeatlas.proteomics.wzw.tum.de`.
  *Caveat*: the per-protein table's exact supplementary file was identified by size, not content — verify
  before committing (see `docs/degradability_report.md` §4).

### 3.3 Proteolytic accessibility (track 3.4)

- **Cappelletti V, Hauser T, Piazza I, et al. 2021** — dynamic 3D proteomes by **LiP-MS**: structural
  fingerprints for **~1,900** *E. coli* proteins (BW25113 — the same strain as Mateus 2018, so a clean join)
  across 8 carbon sources, at peptide resolution. Yields per-protein **protease-accessibility density**
  (LiP peptides normalised by tryptic-peptide count), half-tryptic fraction, and count of significantly
  changing peptides = local conformational plasticity. **The most mechanistically apt single column available
  for this axis** — limited proteolysis reports exactly the accessible, flexible-region property that AAA+
  engagement requires. *Cell* 184:545–559.e22 ·
  DOI [10.1016/j.cell.2020.12.021](https://doi.org/10.1016/j.cell.2020.12.021) · PMC7836100 ·
  PRIDE **PXD022297**.

### 3.4 Solubility, aggregation, refoldability, chaperone dependence (track 3.4)

- **Niwa T, Ying BW, Saito K, et al. 2009** — **eSOL**: chaperone-free solubility (% soluble) for **3,198**
  *E. coli* proteins individually synthesised in the PURE reconstituted translation system; bimodal
  distribution. Keyed on **b-numbers / JW IDs**, so joins are trivial. Best coverage of any dataset in this
  list. *PNAS* 106(11):4201–6 ·
  DOI [10.1073/pnas.0811922106](https://doi.org/10.1073/pnas.0811922106) · PMID 19251648.
  Download: `https://dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip` (NBDC LSDB Archive mirror — the
  original `tanpaku.org` / `tp-esol.genes.nig.ac.jp` site is inactive).

- **Niwa T, Kanamori T, Ueda T, Taguchi H 2012** — global analysis of chaperone effects in the same cell-free
  system: per-protein **ΔSolubility with GroEL / trigger factor / DnaK**, i.e. a chaperone-dependence measure.
  *PNAS* 109:8937–8942 · DOI [10.1073/pnas.1201380109](https://doi.org/10.1073/pnas.1201380109).

- **Györkei Á, Daruka L, Balogh D, et al. 2022** — proteome-wide **in-vivo solubility limits**: 2,577 cytosolic
  proteins screened (ASKA GFP library, K-12 AG1); per-protein solubility threshold plus a 3-class label
  {soluble, rapidly aggregating, slowly aggregating}. *Sci Rep* 12:6547 ·
  DOI [10.1038/s41598-022-10427-1](https://doi.org/10.1038/s41598-022-10427-1) · PMID 35449391 · PMC9023497.

- **To P, Whitehead B, Tarbox HE, Fried SD 2021** — **non-refoldability is pervasive**: of **1,198** *E. coli*
  proteins, **396 (33%) are non-refoldable** after 2 h; protein- *and* domain-level tables. A kinetic-stability
  / metastability proxy complementary to Tm. *J Am Chem Soc* 143(30):11435–11448 ·
  DOI [10.1021/jacs.1c03270](https://doi.org/10.1021/jacs.1c03270) · PMID 34308638.
  ACS-gated — use the NSF-PAR (`par.nsf.gov/servlets/purl/10311595`) or bioRxiv 2020.08.28.273110 copies.

- **To P, Xia Y, Lee SO, et al. 2022** — proteome-wide map of **chaperone-assisted** refolding in a
  cytosol-like milieu; the follow-up that separates intrinsic from chaperone-rescued refoldability.
  *PNAS* 119 · DOI [10.1073/pnas.2210536119](https://doi.org/10.1073/pnas.2210536119).

- **Calloni G, Chen T, Schermann SM, et al. 2012** — DnaK as the central hub of the *E. coli* chaperone
  network: **~700 DnaK interactors**, of which ~180 are aggregation-prone, and a GroEL client
  classification (**class I 30, class II 80, class III 42 obligate**). Chaperone dependence as an intrinsic
  metastability proxy. *Cell Rep* 1(3):251–64 ·
  DOI [10.1016/j.celrep.2011.12.007](https://doi.org/10.1016/j.celrep.2011.12.007) · PMID 22832197.

- **Ramakrishnan R, Houben B, Kreft Ł, Botzki A, Schymkowitz J, Rousseau F 2020** — **PHDB** (Protein
  Homeostasis Database): **4,305** *E. coli* proteins × >100 proteostatic parameters, aggregating 13
  whole-proteome chaperone studies plus abundance, translation rate, solubility, net charge, amino-acid
  composition, secondary-structure content, contact order, TANGO aggregation propensity and IUPred disorder.
  A potential one-stop shortcut for several track-3.4 columns **if** a bulk export exists — the paper
  describes a web interface only. *Bioinformatics* 36(3):948–949 ·
  DOI [10.1093/bioinformatics/btz628](https://doi.org/10.1093/bioinformatics/btz628) · `phdb.switchlab.org`.

### 3.5 Abundance (track 3.5)

- **Schmidt A, Kochanowski K, Vedelaar S, et al. 2016** — the quantitative and condition-dependent *E. coli*
  proteome: **2,359** proteins (~55% of ORFs, >95% of proteome mass) in **molecules/cell** across **22
  conditions**, biological triplicate. *Nat Biotechnol* 34(1):104–10 ·
  DOI [10.1038/nbt.3418](https://doi.org/10.1038/nbt.3418) · PMID 26641532 · PMC4888949 (free author
  manuscript carries the tables). Supplies the abundance term Won 2024 measured at r = −0.69.

### 3.6 Topology / localization covariates

- **STEPdb 2.0** — *E. coli* K-12 subcellular localization, signal peptides, TM topology, abundance,
  solubility, disorder and heat resistance as UniProt-keyed spreadsheets. `stepdb.eu` · PMC4256514.
  Useful covariate block; note several columns are *predicted* (TMHMM/SignalP/Phobius/IUPred), not
  experimental, and must be labelled as such.

---

## 4. *K. pneumoniae*-native resources

- **Illenseher M, Hentschker C, Gesell Salazar M, et al. 2025** — global quantitative proteome of a
  multi-resistant Kp strain: **>2,800** proteins in the cell pellet (54% of the annotated proteome) plus 2,300
  in the exoproteome, **iBAQ**, under control, heat stress (37→42 °C) and oxidative stress (0.3 mM H₂O₂).
  Strain **ATCC BAA-2146 Δ*wza*** (NDM-1, ST11, capsule-deficient); identifiers are NCBI locus tags
  `Kpn2146_####` (genome CP006659.2), so a DIAMOND sequence mapping onto HS11286 is required.
  *Front Microbiol* 16:1528869 ·
  DOI [10.3389/fmicb.2025.1528869](https://doi.org/10.3389/fmicb.2025.1528869) · PMC12127431 ·
  PRIDE **PXD052921**. Track 3.5 — the **Kp-native abundance anchor**.

- **Muselius B, Sukumaran A, Yeung J, Geddes-McAlister J 2020** — **the only Kp-native proteolysis dataset**:
  Δ*lon* vs WT vs complemented (K52 serotype), under iron limitation and repletion. 2,074 proteins quantified;
  59 significantly changed in Δ*lon* vs WT, of which **26 accumulate — candidate Kp Lon substrates** (e.g.
  BamA 8.91 log₂, LuxR 7.02, CmlA2 6.65). Only Δ*lon* was tested; no ClpP knockout.
  *Front Microbiol* 11:546 · DOI [10.3389/fmicb.2020.00546](https://doi.org/10.3389/fmicb.2020.00546) ·
  PMC7194016 · PRIDE **PXD015623**. Track 3.2b validation anchor. Companion secretome:
  *Microbiol Resour Announc* 2025 · DOI [10.1128/mra.00311-25](https://doi.org/10.1128/mra.00311-25).

- **Lin NM, Marino EC, Schlotmann JM, Rosen DA 2026** — differential contributions of ClpX and ClpP to
  pulmonary virulence in classical and hypervirulent Kp: **neither is essential in vitro** (Δ*clpX* and Δ*clpP*
  grow like WT over 18 h in ATCC 43816 and TOP52); in vivo Δ*clpX* is lung-attenuated in **both** pathotypes
  while Δ*clpP* largely retains virulence — so **ClpX has ClpP-independent roles** and the two are not
  interchangeable. ClpX loss reduces capsule in hvKp and increases piliation in both.
  *Infect Immun* 94(3) · DOI [10.1128/iai.00680-25](https://doi.org/10.1128/iai.00680-25).
  Track 3.0 — the resistance-liability framing.

- **Pan YJ, Lin TL, Hsu CR, et al. 2011** — Dictyostelium screen implicating Kp *clpX* in phagocytosis
  resistance and virulence; Δ*clpX* shows reduced capsule gene expression, and the phenotype is not explained
  by capsule loss alone. *Infect Immun* 79(3):997–1006 ·
  DOI [10.1128/IAI.00906-10](https://doi.org/10.1128/IAI.00906-10) · PMID 21173313 · PMC3067516.

- **Bojer MS, Struve C, Ingmer H, et al. 2010** — **ClpK**: a plasmid-encoded Clp ATPase conferring heat
  resistance and proposed to underlie nosocomial persistence in *Klebsiella*. *PLoS ONE* 5:e15467 ·
  DOI [10.1371/journal.pone.0015467](https://doi.org/10.1371/journal.pone.0015467).
  Present in HS11286 as `A0A0H3GSF4` (KPHS_23030, 884 aa). Explicitly **not** a degradation channel — a
  ClpB-class disaggregase with no ClpP partner and no characterised substrate motif. Review:
  "Caseinolytic proteins (Clp) in the genus *Klebsiella*: special focus on ClpK", PMC8746953.

- **UniProt HS11286 proteome** `UP000007841` — the accessions in the §3.0 machinery table were retrieved live
  from `rest.uniprot.org` and cross-checked against `data/processed/other/orthology/kp_orthologs_long.tsv`
  (the repo's own OrthoFinder output) for the two entries whose gene symbols are missing (`clpA`, `sspB`).

---

## 5. Reviews

- **Junker S, Clausen T 2026** — BacPROTAC-induced protein degradation as a new antibiotic concept, from the
  lab that invented the modality. Lists the remaining obstacles as **cell entry, substrate selectivity and PK**;
  notably does not resolve the Gram-negative handle problem. *Trends Biochem Sci* 51(3):249–260 ·
  DOI [10.1016/j.tibs.2025.12.012](https://doi.org/10.1016/j.tibs.2025.12.012) · PMID 41708383.

- **Bazzacco A, Mercorelli B, Loregian A 2025** — PROTACs and beyond for microbial pathogens; the broadest
  coverage. Identifies **only ClpXP and ClpC1P1P2** as demonstrated recruitable handles; Gram-negative
  coverage is thin (CLIPPERs and the CTX-M work only) and it does not treat ClpAP/ClpS/Lon/HslUV/FtsH as
  handles or discuss OM permeability. *FEMS Microbiol Rev* 49 ·
  DOI [10.1093/femsre/fuaf046](https://doi.org/10.1093/femsre/fuaf046).

- **Petkov R, Camp AH, Isaacson RL, Torpey JH 2023** — targeting bacterial degradation machinery. Explicit on
  the phylogenetic split — **Gram-positives use ClpC/ClpE, Gram-negatives use ClpX/ClpA** — and judges **Lon
  "unsuitable for antibiotic development at present"** and HslUV mechanistically under-characterised.
  *Biochem J* 480(21):1719–1731 · DOI [10.1042/BCJ20230191](https://doi.org/10.1042/BCJ20230191) · PMC10657178.
  The basis for not giving HslUV its own channel.

- **Jeong YC, Kim SH, Moon S, Kim H, Lee C 2026** — proteostasis-targeted antibacterial strategies; frames the
  four ubiquitous eubacterial proteolytic complexes (Lon, HslUV, ClpXP, FtsH) and states the criterion this
  axis operationalises: *"endogenous proteins with disordered termini are effective BacPROTAC targets."*
  *J Microbiol* 64(3):e2511007 · DOI [10.71150/jm.2511007](https://doi.org/10.71150/jm.2511007).

- **Izert MA, Klimecka MM, Górna MW 2021** — applications of bacterial degrons and degraders; still the best
  catalogue of degron tools and screens (ssrA, DAS+4, SspB systems).
  *Front Mol Biosci* 8:669762 · DOI [10.3389/fmolb.2021.669762](https://doi.org/10.3389/fmolb.2021.669762) ·
  PMC8138137.

- **Sarwan DN, Khedekar PB, Bhole RP, Chikhale RV 2025** — PROTAC technology against AMR. Useful for the
  Gram-negative liabilities: permeability of large bifunctional molecules, and resistance via ClpC2/ClpC3
  buffering or loss of the recruiter. *npj Antimicrob Resist* 3:94 ·
  DOI [10.1038/s44259-025-00136-w](https://doi.org/10.1038/s44259-025-00136-w) · PMID 41345312.

- **Sarathy JP, Aldrich CC, Go ML, Dick T 2023** — "PROTAC antibiotics: the time is now"; the event-driven vs
  occupancy-driven argument that underpins track 3.9. *Expert Opin Drug Discov* 18(4):363–370 ·
  DOI [10.1080/17460441.2023.2178413](https://doi.org/10.1080/17460441.2023.2178413) · PMID 37027333.

---

## 6. Tools & methods

- **Moreno J, Nielsen H, Winther O, Teufel F 2024** — **DeepLocPro**: prokaryote-specific subcellular
  localization (6 classes: cytoplasm, cytoplasmic membrane, periplasm, outer membrane, cell wall/surface,
  extracellular), trained on UniProt + PSORTdb, benchmarking above the PSORTb 3.0 ensemble, CPU-friendly.
  *Bioinformatics* 40(12):btae677 ·
  DOI [10.1093/bioinformatics/btae677](https://doi.org/10.1093/bioinformatics/btae677) · PMID 39540738 ·
  PMC11645106 · `services.healthtech.dtu.dk/services/DeepLocPro-1.0/` · `ku.biolib.com/deeplocpro`.
  Planned `09c` — needed to close the 3,414/5,728 localization gap that blocks track 3.7.

- **Emenecker RJ, Griffith D, Holehouse AS 2021** — **metapredict**: fast, pip-installable consensus disorder
  predictor. *Biophys J* 120:4312–4319 ·
  DOI [10.1016/j.bpj.2021.08.039](https://doi.org/10.1016/j.bpj.2021.08.039). Sequence-side disorder for
  track 3.1 (second opinion alongside pLDDT; also covers the one Kp protein with no AlphaFold model).

- **MobiDB** — per-accession disorder consensus plus ~13 individual predictors (including an AlphaFold-derived
  track), with regions and content fractions. **Verified to cover *K. pneumoniae* HS11286**, and a
  proteome-wide bulk download works:
  `https://mobidb.org/api/download?proteome=UP000007841&format=tsv`. Preferable to running a predictor
  ourselves for track 3.1's sequence-side disorder — it is a consensus, it is citable, and it needs no new
  environment. `mobidb.org` · already flagged as a wanted addition in `docs/01_task_agnostic.md`.

- **Erdős G, Dosztányi Z 2024** — **AIUPred**: energy-estimation + deep-learning disorder predictor with a
  standalone Python package. *Nucleic Acids Res* 52:W176–W181 ·
  DOI [10.1093/nar/gkae385](https://doi.org/10.1093/nar/gkae385) · `github.com/doszilab/AIUPred`.
  Alternative to metapredict; already flagged in `docs/02_ligandability.md` as the preferred disorder upgrade.

- **Hu G, Katuwawala A, Wang K, et al. 2021** — **flDPnn**: the disorder predictor Won 2024 actually used
  (threshold 0.2). *Nat Commun* 12 ·
  DOI [10.1038/s41467-021-24773-7](https://doi.org/10.1038/s41467-021-24773-7). Web server only; exact parity
  is not achievable locally, which is one reason to re-fit rather than reuse Won's coefficients.

- **Rubach P, Sikora M, Jarmolinska AI, et al. 2024** — **AlphaKnot 2.0**: knotting of AlphaFold-predicted
  models, with a database. *Nucleic Acids Res* 52:W187–W193 ·
  DOI [10.1093/nar/gkae443](https://doi.org/10.1093/nar/gkae443). Source of the topology feature in track 3.4.

- **San Martín Á, Rodriguez-Aliaga P, Molina JA, et al. 2017** — knots can impair degradation by ATP-dependent
  proteases: a knot tightened against a mechanically stable domain stalls translocation and releases partially
  degraded product; ATP cost rises 30–1000-fold. *PNAS* 114:9864–9869 ·
  DOI [10.1073/pnas.1705916114](https://doi.org/10.1073/pnas.1705916114) · PMC5604015.
  **Counter-evidence to weigh**: **Sivertsson EM, Jackson SE, Itzhaki LS 2019** — ClpXP *can* readily degrade
  3₁- and 5₂-knotted proteins. *Sci Rep* 9 ·
  DOI [10.1038/s41598-018-38173-3](https://doi.org/10.1038/s41598-018-38173-3) · PMC6382783.
  → implement as a **soft penalty, not a veto**.

- **Schweke H, Pacesa M, Levin T, et al. 2024** — atlas of protein homo-oligomerization across domains of
  life: AlphaFold2-based homomer prediction over four proteomes including *E. coli*; **~45% of a bacterial
  proteome forms homomers**. *Cell* 187:999–1010.e15 ·
  DOI [10.1016/j.cell.2024.01.022](https://doi.org/10.1016/j.cell.2024.01.022). Optional input to track 3.6 if
  its bulk table proves fetchable; the primary route is UniProt subunit annotation + PDB assembly
  stoichiometry.

- **Szulc NA, Stefaniak F, Piechota M, et al. 2024** — **DEGRONOPEDIA**: proteome-wide degron inspection,
  mapping degrons to disordered regions ("unfolding seeds") and simulating proteolysis to find degrons in
  newly formed termini. *Nucleic Acids Res* 52(W1):W221–W232 ·
  DOI [10.1093/nar/gkae238](https://doi.org/10.1093/nar/gkae238) · `degronopedia.com`.
  **Eukaryote/ubiquitin-oriented, so not directly usable** — cited because its "degron + nearby disorder +
  proteolysis-generated termini" architecture is exactly the design this axis adopts for bacteria.

- **Burgos R, Weber M, Martinez S, et al. 2020** — protein quality control and regulated proteolysis in
  *Mycoplasma pneumoniae*: conditional Lon and FtsH depletion yielding **62** and **34** candidate substrates
  plus protein half-lives. *Mol Syst Biol* 16 ·
  DOI [10.15252/msb.20209530](https://doi.org/10.15252/msb.20209530) · PMID 33320415 · PMC7737663.
  Cross-species trap set for track 3.3.

- **Bhat NH, Vass RH, Stoddard PR, et al. 2013** — ClpP substrates in *Caulobacter crescentus*.
  *Mol Microbiol* 88:1083–1092 · DOI [10.1111/mmi.12241](https://doi.org/10.1111/mmi.12241).

- **Feng J, Michalik S, Varming AN, et al. 2013** — ClpP-trap substrates in *Staphylococcus aureus* (~70
  proteins). *J Proteome Res* 12:547–558 · DOI [10.1021/pr300394r](https://doi.org/10.1021/pr300394r).

- **Lunge A, Gupta R, Choudhary E, et al. 2020** — *M. tuberculosis* ClpC1 regulates a subset of proteins
  **having intrinsically disordered termini** — the observation the previous version of this spec cited, now
  superseded quantitatively by Won 2024. *J Biol Chem* 295:9455–9473 ·
  DOI [10.1074/jbc.RA120.013456](https://doi.org/10.1074/jbc.RA120.013456).

---

## 7. Rejected, deferred, or not obtained

> **⚠ Superseded in part, 2026-08-06.** Three entries below are **no longer "not obtained"** — the PDFs are on
> disk in `data/raw/other/degradability/literature/`, supplied manually:
>
> | Paper | File | What it unblocks |
> | --- | --- | --- |
> | **Flynn 2003** | `PIIS1097276503000601.pdf` | 60-protein ClpXP trap census **and the real positional motif consensuses** (see below) |
> | **Neher 2006** | `PIIS1097276506001687.pdf` | ~100-protein SILAC ClpXP trap ± DNA damage, with a *clpX*⁻ control |
> | **Ziemski 2021** | `febs15335-sup-0001-TableS1.pdf` | Mtb ClpCP interaction screen — still deprioritised (ClpC-family, capped 0.45) |
>
> Flynn + Neher take the *E. coli* ClpXP selection bar from the 45-row hand-curated stand-in to roughly
> **130–160 real proteins**. Extraction strategy and the resulting code correction:
> `docs/degradability_datasets.md` §3.1 and §5.

- **Flynn et al. 2003 — the positional constraints are now available, and they contradict our own exclusion.**
  *(Historical note: cell.com and sciencedirect.com both returned HTTP 403, the paper has no PMC record, and the
  MIT DSpace copy blocked scripted fetch, so this was previously quotable only from secondary sources — Kim &
  Kim, PMC2394798 — for the motif strings but not their positions.)* With the paper in hand, the published
  consensuses are:
  **N-M1** strict `T₁-X₂-K₃-[ILV]₄`, **located 1–4 residues from the N terminus** (looser:
  `polar-T/φ-φ-basic-φ`), 18 of the 60 trapped proteins · **N-M2** 12 · **N-M3**
  `φ-X-polar-X-polar-X-basic-polar`, 10 · **C-M1** ssrA-like, critical region terminal `-Ala-Ala` enriched
  **7-fold** · **C-M2** MuA-like. ~90% of trapped proteins carry a candidate signal.
  **This invalidates the stated reason `src/degradability.py` gives for dropping NM1/NM3** — that reason cites
  `^M?[AILVMFW]{2,}` firing on 30.7% of the Kp proteome, but that regex is not NM1, and a four-residue anchored
  pattern cannot fire at that prevalence. Re-implement from the real consensus and re-measure before weighting.
  Also supplies a **positive-control peptide set**: ssrA `AANDENYALAA`, YdaM `KNDGRNRVLAA`, Crl `DFRDEPVKLTA`,
  LldD `ALAPMAKGNAA`, MuA `ILEQNRRKKAI`, YbaQ `ARREERAKKVA`, RibB `AYRQAHERKAS`, Gcp `RWPLAELPAA`, and Dps
  `FLWFIESNIE` (a negative for the `-AA` rule).
  Unchanged either way: the **hand-curated 45-protein table in `data/raw/legacy/clp_substrates/` must stop being
  described as Flynn's census.**

- **The whole cross-bacterial Clp-trap pool is gated.** Verified against the NCBI OA service: Lunge 2020
  (PMC7363115), Bhat 2013 (PMC3681837) and Graham 2013 (PMC3807464) all return `idIsNotOpenAccess`; Feng 2013
  (ACS) and Ziemski 2021 (Wiley) have no PMC record at all. **Ziemski 2021's Table S1 is now on disk** (`data/raw/other/degradability/literature/`), so the remaining un-obtained member of this pool is Feng 2013. The old §3.3b track therefore depended entirely on
  content no automated route can reach, which is a second reason (besides label count) to base the axis on the
  open *E. coli* turnover datasets instead. By contrast **Nagar 2021 is CC-BY and fully automatable** — the
  Europe PMC `supplementaryFiles` endpoint returns a zip directly.

- **Tsuboyama K, Dauparas J, Chen J, et al. 2023** — mega-scale folding stability (~776,000 curated ΔG values).
  *Nature* 620:434–444 · DOI [10.1038/s41586-023-06328-6](https://doi.org/10.1038/s41586-023-06328-6) ·
  Zenodo [10.5281/zenodo.7992926](https://doi.org/10.5281/zenodo.7992926).
  **Rejected for direct use**: the 331 natural domains are 40–72 aa, selected for cDNA-display compatibility,
  and not mapped to any proteome — using it would mean training a predictor, making the result a *predicted*
  rather than experimental feature. There is **no** proteome-scale experimental ΔG for any bacterium; Tm
  (§3.2) is the practical stand-in.

- **HslUV as a separate channel — dropped.** Substrate list is tiny (SulA, RcsA, RpoH, TraJ, RNase R, YbaB);
  recognition needs only "an aromatic π-system or a cation at the degron position" (PMC10743992), which is
  degenerate with the Lon aromatic-cluster feature. An *E. coli* proteome-microarray screen (PMID 27864322)
  yielded essentially one validated new substrate. Judged mechanistically under-characterised by Petkov 2023.

- **Piazza I, et al. 2018** — LiP-SMap on *E. coli* lysate; 7,345 predicted binding sites.
  *Cell* 172:358–372 · DOI [10.1016/j.cell.2017.12.006](https://doi.org/10.1016/j.cell.2017.12.006) ·
  PRIDE **PXD006543**. **Deferred** — the openly available slice (Mendeley Data
  [10.17632/nhsktkcs3d.1](https://doi.org/10.17632/nhsktkcs3d.1)) is only Data S8–S10; the peptide-level LiP
  tables are paywalled. Cappelletti 2021 (§3.3) supersedes it for this purpose.

- **Wang Z, Han QQ, Zhou MT, et al. 2016** — dynamic-SILAC protein turnover in *Salmonella* Typhimurium
  (870 proteins in DMEM, 311 during macrophage infection; median t½ 99 → 69 min).
  *J Basic Microbiol* 56:801–811 · DOI [10.1002/jobm.201500315](https://doi.org/10.1002/jobm.201500315) ·
  PMID 26773230. **Low priority** — paywalled, no PMC, small coverage; the only enteric-pathogen turnover set.

- **Nichols BP, et al. 2011** chemical-genomic profiles and other reCAPTCHA-gated PMC content — noted here only
  because the essentiality axis hit the same wall; not needed for this axis.

- **PaxDb** (*K. pneumoniae* species id 667127) — ppm abundance, integrated and clean, but the site is
  JS-rendered so it needs the bulk release files rather than page scraping. Secondary to Illenseher 2025 for
  Kp and to Schmidt 2016 for *E. coli*. `pax-db.org/downloads` · PaxDb v6.0, *Nucleic Acids Res* 54:D427 (2026),
  PMC12807614.

- **No *Klebsiella* ClpP/ClpX/Lon degradome exists**, and there is no Kp ClpP structure (UniProt
  `A0A0H3GKH6` has an AlphaFoldDB entry and no PDB cross-reference). No Salmonella or Pseudomonas degradome
  either. All recognition rules transfer from *E. coli* by orthology — defensible for Lon (conserved
  E. coli/Yersinia/Mycoplasma per Cragan 2025) but with an explicit caution for ClpX, whose specificity is
  demonstrably **not** universally conserved (*S. mutans* ClpXP recognises a tripeptide `LPF` that *E. coli*
  ClpXP does not; PMC5143411).

---

## 8. E. coli-first additions (2026-08-05 pass)

Added after the decision to build *E. coli* first. Everything here is new relative to §1–§7, and the four
envelope-protease rule sets were missed originally because they were deprioritised as *Klebsiella*-irrelevant.

### 8.1 The envelope protease rules — four channels the earlier spec had nothing for

- **RseP / YaeL — intramembrane metalloprotease (site-2 family).** Cleaves within the bilayer; requires
  **helix-destabilising residues in the substrate TM segment**, which also stabilise the substrate–RseP
  interaction, plus (for RseA) an exposed C-terminal hydrophobic residue generated by the prior DegS cut. The
  TM helix is unwound by strand addition to RseP's intramembrane β-sheet and clamped by a conserved Asn.
  Substrates include RseA, FecR and remnant signal peptides. **Essential in both organisms** (measured:
  Ec 0.88, Kp 0.82, both experimental). Key sources: *EMBO J* 2004 (RseP cleaves TM sequences);
  *eLife* 2015 (membrane-reentrant β-hairpin loop and selective cleavage); *PNAS* 2009 (C-terminal hydrophobic
  requirement); *Sci Adv* 2022 `10.1126/sciadv.abp9011` (mechanism); *JBC* (substrate recognition and binding).
  Track 3.2f.

- **GlpG — rhomboid intramembrane serine protease.** Recognition is a **sequence motif around the cleavage
  site**, and the authors state explicitly that it is *more strictly required than* TM-helix-destabilising
  residues, with a preference for a **small side chain at P1 and a negative charge at P1′** — and that
  genome-wide substrate prediction from the motif is feasible. Sources: *Nat Struct Mol Biol* / PMC2941825
  (recognition motif in rhomboid substrates); PMID 17501925 (sequence features required for GlpG cleavage);
  PMC4253528 (substrate–peptide complex structures, S1–S4 subsites). Track 3.2g.

- **DegS — the natural molecular glue.** A trimeric periplasmic serine protease with a PDZ domain per subunit.
  Peptides ending **`φ-x-Phe`** — the YxF motif of outer-membrane porins — bind the PDZ domain and activate
  proteolysis by **relieving inhibitory contacts** that otherwise capture loop L3. The terminal Phe is the key
  interaction; Phe/Tyr preferred at −2. Crucially the motif is **inaccessible in folded, membrane-embedded
  OMPs** and exposed only on misfolding, so the feature is conditional on foldedness. **Essential in both
  organisms** (Ec 0.88, Kp 0.60, both experimental). This is the natural precedent for the engineered DegP
  glues. Sources: *Cell* 2007 (allosteric activation of DegS); *Genes Dev* 2007 21:2659; *Mol Cell* 2008
  (OMP peptides modulate DegS by differential binding); PMC2764547 (relief-of-inhibition mechanism);
  *Nat Chem Biol* 2012 (shared energy landscape). Track 3.2h.

- **DegP and Prc / Tsp.** DegP cleaves **between paired hydrophobic residues** (PMC139609, PapA pilin) and
  switches between chaperone and protease. Prc/Tsp is C-terminal-specific and degrades **ssrA-tagged proteins
  that get exported** to the periplasm. Two 2025 refinements below. Track 3.2i.

- **Petkov et al. 2023** *Biochem J* 480:1719 (PMC10657178) — cited here for a **negative** reason: the
  field's main review on targeting bacterial degradation machinery **does not discuss FtsH, DegP, DegS, RseP or
  ClpS at all**. It is ClpCP-centric. That is the evidence that the envelope proteases are unexplored space
  rather than something the literature has already covered.

### 8.2 Real E. coli protein termini — the layer the N-degron channel depends on

- **Van Damme et al. 2026 — TRAINSPOTTER.** Deformylation-assisted N-terminomics in *E. coli* K-12 (CAG12184):
  **1,082 N-termini → 806 translation-initiation sites → 729 proteins**, resolving closely spaced start codons
  that ribosome profiling cannot. Table S1 carries **UniProt accession + b-number** plus the observed
  N-terminal peptide, iMet-excised vs retained, alternative/near-cognate starts, N-terminal
  extensions/truncations and Nt-acetylation. *Nucleic Acids Res* 54:gkag587 ·
  DOI [10.1093/nar/gkag587](https://doi.org/10.1093/nar/gkag587) · PMID 42258537 · PMC13245402 · CC-BY ·
  raw PXD005901. **The single best fix for "the annotated Met1 is often wrong."** Track 3.2c.

- **Bienvenut, Giglione & Meinnel 2015.** SILProNaQ positional proteomics ± actinonin: >1,000 unique N-termini
  with **quantitative modification fractions — 56% initiator-Met removed, 10% Nt-acetylated, 5% N-formyl
  retained** — plus 140 signal-peptide-cleaved mature termini. Gives a *probability* rather than a boolean.
  *Proteomics* 15:2503 · DOI [10.1002/pmic.201500027](https://doi.org/10.1002/pmic.201500027) · PMID 26017780.
  ⚠ Table is Wiley-paywalled and PRIDE holds **RAW only** (PXD001979 / PXD002012 / PXD001983 = 166 RAW files,
  no processed identification table). Institutional access or author contact required.

- **Chen et al. 2020 — CPB-ChaFRADIC.** The **only native E. coli C-terminome**: 604 canonical + 818
  neo-C-termini, described by the authors as the largest for E. coli. *Anal Chem* 92 ·
  DOI [10.1021/acs.analchem.0c00762](https://doi.org/10.1021/acs.analchem.0c00762) · PMID 32441514 ·
  PXD018520. Matters because our C-degron channels are calibrated entirely on *engineered* library tags.

- **Jachmann et al. 2026 — E. coli PeptideAtlas build 585.** 40 datasets, >73 M spectra, evidence for **4,755
  proteins** (1,376 with no prior protein-level support) and >10,000 PTM sites. Its `coordinate_mapping.txt`
  gives per-peptide start/end position and flanking residues proteome-wide, from which empirical N-/C-termini
  and non-tryptic cleavage sites can be derived. *J Proteome Res* 25:1027 ·
  DOI [10.1021/acs.jproteome.5c00902](https://doi.org/10.1021/acs.jproteome.5c00902) · PMID 41568995 ·
  PMC12887991. Also the natural **observability denominator** for a per-protein score.

- **Cartwright, Jha & Smith 2025** — structure and mechanism of the aminoacyl-tRNA-protein L/F- and
  R-transferases. Pins down the enzymology of **Aat**, which installs Leu/Phe on N-terminal Arg/Lys to *create*
  ClpS-recognised N-degrons. Consequence for track 3.2c: the rule is not "ends up with a bulky hydrophobic
  N-terminus" but "**is an Aat substrate**". *J Mol Biol* ·
  DOI [10.1016/j.jmb.2025.169210](https://doi.org/10.1016/j.jmb.2025.169210) · PMID 40381981.

### 8.3 Assembly state and synthesis flux — both stronger for E. coli than assumed

- **Kshirsagar et al. 2025 — Seq2Symm.** Sequence-only protein-language-model classifier of homo-oligomer
  symmetry. The only route to **100% coverage of our exact proteomes with no strain mismatch and no identifier
  mapping**; fits the repo's existing local-predictor pattern. *Nat Commun* 16:1969 ·
  DOI [10.1038/s41467-025-57148-3](https://doi.org/10.1038/s41467-025-57148-3) · PMID 40016259 · PMC11868566 ·
  code `github.com/microsoft/seq2symm`. **Preferred over Schweke 2024 (see §9).** Track 3.6.

- **RCSB / PDBe assembly APIs.** `rcsb_struct_symmetry` returns `oligomeric_state` literally as `"Homo 2-mer"`
  with `stoichiometry` and symmetry, over **7,671 E. coli K-12 biological assemblies**. PDBe's
  `uniprot/interface_residues/<ACC>` returns interface residues **in UniProt numbering, derived from biological
  assemblies, with partner identity** — so self-burial in a homomer (the GroEL case) is distinguishable from
  burial against a different subunit. Together these are the honest implementation of track 3.6.

- **Complex Portal** (EBI/IntAct) — 324 curated E. coli complexes, 781 participant accessions, **with explicit
  stoichiometry** and a `Complex assembly` column. *Nucleic Acids Res* 2022, PMC8689886. Gold-standard
  validation set for the above.

- **STRING v12.0**, taxon **511145**. Per-protein interaction degree: 4,006 nodes at `combined_score ≥ 700`
  = **91% of the proteome**. Directly addresses Nagar 2021's finding that PPI-network properties outrank every
  sequence motif for predicting half-life. Use the `physical.links` subnetwork for the assembly-burial signal
  specifically, and join via UniProt `xref_string` (4,087/4,403) rather than the ambiguous aliases file.

- **Li, Burkhardt, Gross & Weissman 2014.** Absolute protein **synthesis rates** (molecules per cell per
  generation) for 4,095 genes, 3,041 above threshold, plus mRNA and translation efficiency. The true
  resynthesis-flux denominator for track 3.5, better than abundance. *Cell* 157:624 ·
  DOI [10.1016/j.cell.2014.02.033](https://doi.org/10.1016/j.cell.2014.02.033) · PMID 24766808 · PMC4006352.
  Identifier is a bare gene name; ~360 rows need a synonym crosswalk. Supplements download from the Elsevier
  CDN (`ars.els-cdn.com/content/image/1-s2.0-{PII}-mmcN.xlsx`), which works where PMC `/bin/` and cell.com do
  not. Companion: **Irshad & Sharma 2023**, *Biophys Rep* 3:100131 (PMC10542608, OA) — translation
  **initiation rates** for 1,533 / 2,709 genes keyed on **b-numbers**, a cleaner identifier.

- **Onyenemezu et al. 2026** — proteome-wide protein stability under acid stress: pH50 for >90% of **1,675**
  E. coli proteins (range 2.28–6.33, median 5.11); ~9% acid-stable throughout. Notably **stability stratifies
  by compartment** — periplasmic proteins are enriched among the acid-stable, cytoplasmic 4.5–5.5,
  inner-membrane 5.5–6.0 — with no correlation to pI or MW. A biophysical lability axis orthogonal to turnover.
  *J Ind Microbiol Biotechnol* · DOI [10.1093/jimb/kuag016](https://doi.org/10.1093/jimb/kuag016) ·
  PMID 42384035. ⚠ not open access.

### 8.4 New modalities and rule changes (2025–2026)

- **Han et al. 2025 — bacNID.** A **second ClpXP handle**: a gold nanoparticle co-grafted with a
  target-binding peptide and an **SspB-binding peptide**, hijacking the SspB adaptor rather than binding ClpX
  directly. Includes resistance-evolution experiments and an in vivo wound-infection model. Peptide-on-carrier,
  not a small molecule. *Nat Commun* 16 ·
  DOI [10.1038/s41467-025-66221-w](https://doi.org/10.1038/s41467-025-66221-w) · PMID 41413196 · OA.

- **Coriano-Ortiz et al. 2026 — orthogonal light-gated Lon.** *Mesoplasma florum* Lon with a LOV2 domain
  inserted at nearly every codon (726 variants screened by FACS), orthogonal to native E. coli machinery,
  intensity-tunable — and it **degrades in stationary phase as well as exponential**, decoupling measured
  degradation from division-based dilution. Lon's first working handle in E. coli. bioRxiv ·
  DOI [10.64898/2026.06.17.732960](https://doi.org/10.64898/2026.06.17.732960).

- **Brown et al. 2026 — BacPROTACs outperform inhibitors.** First heterobifunctional BacPROTACs against an
  essential *M. tuberculosis* target (PptT), built by recycling existing inhibitors into degraders, with
  markedly better activity than the parents; also publishes a generalisable pipeline for measuring degradation
  in bacteria. Still a ClpC1 handle, so not transferable to Enterobacteriaceae — but the strongest
  "degradation beats inhibition" evidence in bacteria, and the assay cascade to validate any score. bioRxiv ·
  DOI [10.64898/2026.06.12.731830](https://doi.org/10.64898/2026.06.12.731830).

- **Iqbal, Keller & Ghanbarpour 2026.** FtsH forms an inner-membrane complex with the SPFH proteins
  **HflK/HflC**; a disulfide-crosslinked closed state was solved by cryo-EM, and cells locked closed — or with
  FtsH-binding-disrupting HflK/C mutants — **fail to recover from aminoglycosides**. Structures from
  tobramycin-treated cells show two openings that facilitate substrate entry under stress. **FtsH degradability
  is therefore gated, not intrinsic**, and track 3.2e needs an accessibility term alongside the degron.
  *Cell Rep* 45:117231 · DOI [10.1016/j.celrep.2026.117231](https://doi.org/10.1016/j.celrep.2026.117231) ·
  PMID 41964960 · OA.

- **Islam et al. 2026.** Direct rate constants on titin-I27 tandem substrates: unfolding + translocation at
  **12.0 aa s⁻¹ for ClpA alone vs 33.2 for ClpAP**, with translocation 8–24× faster than unfolding either way,
  and ClpAP unfolding ~3× faster than ClpA alone. **For ClpAP unfolding is rate-limiting, not translocation** —
  so the ClpAP term should be dominated by local mechanical stability next to the degron rather than chain
  length. *Biophys J* · DOI [10.1016/j.bpj.2026.01.007](https://doi.org/10.1016/j.bpj.2026.01.007) ·
  PMID 41520171.

- **González et al. 2025.** In-cell NMR of NDM-1 degradation under zinc starvation: **Prc cleaves
  membrane-bound NDM-1 at specific residues and DegP then processes the Prc-generated peptides** — a concerted
  two-protease periplasmic pathway at residue resolution in live cells. For β-lactamases the periplasmic route
  may be **Prc→DegP rather than DegP alone**, so track 3.2i needs a Prc-accessibility term. Intersects directly
  with the Taylor 2026 TEM glues. *Nat Commun* ·
  DOI [10.1038/s41467-025-62340-6](https://doi.org/10.1038/s41467-025-62340-6) · PMID 40993131 · OA.

- **Roy, Nandakumar & Chaba 2025.** DegP's protease-domain disulfide is a **redox sensor**: under long-chain
  fatty acids or alkaline pH it accumulates in the thiol, protease-active form, via two routes (impaired
  disulfide-bond machinery; thiol-containing substrates converting DegP allosterically). So DegP degradability
  is conditional on redox state **and on the substrate's own free-cysteine content**. bioRxiv ·
  DOI [10.1101/2025.11.03.686187](https://doi.org/10.1101/2025.11.03.686187).

- **Lee et al. 2026 — ProHL.** Deep-learning bacterial protein half-life classifier, validated in E. coli.
  **The direct competitor** to this axis's predictive component. *J Microbiol Biotechnol* ·
  DOI [10.4014/jmb.2601.01071](https://doi.org/10.4014/jmb.2601.01071) · PMID 42046921 · OA.
  ⚠ **Read the corrigendum before benchmarking**: DOI 10.4014/jmb.2026.3602.c02, PMID 42124426.

- **Izert-Nowakowska et al. 2026.** Degron evaluation in *E. coli* with a fluorescent reporter — pBAD
  eGFP–degron fusions, a plate end-point assay and a 96-well kinetic assay for high-throughput degron
  screening. From the CLIPPERs group, i.e. our closest comparator: ingest as **protocol**, and treat any degron
  set they publish as a labelled test set. bioRxiv ·
  DOI [10.64898/2026.03.07.710301](https://doi.org/10.64898/2026.03.07.710301).

- **Klimecka et al. 2021.** The ssrA calibration benchmark: 12 eGFP–degron variants, each measured by SspB
  K_D (MST, 0.08–0.33 µM), ClpX ATPase stimulation, in vitro K_M/V_max **± SspB**, and in vivo time courses in
  WT / Δ*clpX* / Δ*clpP* / Δ*sspB*. The only source that **separates adaptor binding from protease
  engagement**. *Molecules* 26:5936 · DOI [10.3390/molecules26195936](https://doi.org/10.3390/molecules26195936) ·
  PMC8512704. ⚠ no consolidated table — values must be hand-curated from Fig. 4e and Fig. 7.

- **Cronan & Kuzminov 2024.** Degron-controlled degradation in E. coli: DAS+4 variants and destabilised SspB
  alleles, with **15 native genes tagged and a graded outcome** — full inactivation for *recA*/*ruvB*/*dnaG*,
  partial for *polA*/*recB*/*gfp*/*priA*, none for *recC*/*rnhA*/*holD*/*lacZ*/*ligA*. A small **validation
  set** for the axis, plus a real mechanistic constraint: ClpXP capacity saturates. *ACS Synth Biol* 13:669 ·
  DOI [10.1021/acssynbio.3c00768](https://doi.org/10.1021/acssynbio.3c00768) · PMC10659297.

- **Stevens-Cullinane et al. 2025 — LAMP-D.** ⚠ **A negative control, not a training example.** A
  heterobifunctional Ru(II) photosensitiser degrades NDM-1 in a Gram-negative (>100-fold inhibition gain on
  irradiation, 53-fold meropenem potentiation) — but explicitly **without recruiting host proteolytic
  machinery**: this is photo-oxidative backbone scission. It must not contaminate a "degradable by AAA+
  protease" label. *JACS* 147 · DOI [10.1021/jacs.5c12405](https://doi.org/10.1021/jacs.5c12405) ·
  PMID 41308195.

- **Aldikacti et al. 2026.** Tn-seq across strains lacking key chaperones and proteases under multiple
  proteotoxic stresses; shows that major quality-control players **mask** correlations between transcriptomic
  response and gene fitness. A genetic-dependency layer, and an explicit warning against transcriptome-derived
  proxies for proteostasis. *Cell Rep* 45:116892 ·
  DOI [10.1016/j.celrep.2025.116892](https://doi.org/10.1016/j.celrep.2025.116892) · PMID 41575849 · OA.

- **Reitzel & Geddes-McAlister 2026.** Δ*lon* vs WT *K. pneumoniae* extracellular proteome — the only
  2025–26 proteome-scale protease-deletion dataset in Kp. *Microbiol Resour Announc* ·
  DOI [10.1128/mra.00311-25](https://doi.org/10.1128/mra.00311-25) · PMID 41778832 · OA.

- Tooling: **flDPnn3** (*J Mol Biol* 2026, PMID 41500381 — proteome-scale disorder, and the predictor family
  Won 2024's top feature used); **DisProt 2026** (*NAR*, PMID 41249866 — curated disorder ground truth);
  **PEGASUS** (*Protein Sci* 2025, PMID 40671366 — MD-derived flexibility from sequence); **zsasa** (bioRxiv
  2026 — proteome-scale SASA).

### 8.5 Verified negative findings

Worth recording because they bound the search and are defensible in a paper:

- **No bacterial degron database exists**, and none was released in 2025–26. DEGRONOPEDIA's Methods list
  **11 supported organisms and *E. coli* is not among them**; all 46,743 of its degrons are ubiquitin-system
  motifs.
- **No proteome-scale degron-tagging library exists for E. coli**, re-checked through 2026 including CRISPRi
  repurposing. The yeast (>5,600 ORFs, ~90% degraded) and human equivalents both depend on ubiquitin/E3
  machinery bacteria lack. DEtox — a *pentapeptide-tag* library, not an ORF library — remains the state of the
  art.
- **No proteome-scale E. coli ClpS substrate list exists.** Humbard 2013 remains the only one.
- **No new E. coli turnover atlas appeared in 2025–26** (queried four ways). Gupta 2024 and MacKrell 2026
  remain the only proteome-wide resources.
- **No chemical degrader in 2025–26 uses a ClpA/ClpS, FtsH or DegP recruitment handle in a Gram-negative with
  a true small molecule.** Taylor 2026 and Nie 2025 are the frontier.
- **No retractions, corrections or failed replications** among the references in §1–§7, checked by a dedicated
  publication-type sweep over 2025–2026.
- **MEROPS is quantifiably empty for our proteases**: Lon (S16.001) = 2 real substrates (SulA, CcdA);
  ClpP (S14.*) = **0 physiological, the only entry being bovine insulin B-chain**; FtsH (M41.001) = 1
  (β-casein); HslV (T1.006/7) = **0 rows**. This also explains why TopFIND's E. coli cleavage layer is empty,
  since TopFIND ingests MEROPS.
- **TopFIND is resolved but not worth ingesting.** Live at `topfind.clip.msl.ubc.ca` (v4.1), bulk dump at
  `/assets/TopFIND_20211220.sql.gz` (232 MB, no login) — but 10 of 10 sampled E. coli entries carry evidence
  type *"electronic annotation"* with **zero literature PMIDs**, 0 cleavages and 0 substrates. Frozen at
  UniProt 2019_07.
- **PortEco's domain is now an online-gambling spam site** — must never be cited or linked. **EcoGene** is
  gone (HTTP 000). **EcoliWiki** is alive but holds nothing relevant. **EcoCyc**'s only turnover dataset *is*
  Gupta 2024. **QSbio / anti-QSalign** now serve an empty page.

---

## 9. Corrections to earlier entries in this file

- **Schweke et al. 2024** (§6) was described as covering *E. coli*. It covers **E. coli O157:H7, not K-12**:
  the 2,181 `org == "ec"` accessions have an **intersection with UP000000625 of exactly zero**, and only
  820/2,181 match K-12 even by gene name. The table also contains only predicted homomers, so absence conflates
  "monomer" with "not predicted". **Superseded by Seq2Symm (§8.3) for our purposes.**
- **DEGRONOPEDIA** (§6) was described as "eukaryote/ubiquitin-oriented". Stronger than that: *E. coli* is not
  among its 11 supported organisms at all.
- The claim that **FtsH has degrader precedent** appeared in an earlier working note. It does not. Only ClpXP
  and DegP have demonstrated degrader chemistry; Lon has an engineered genetic handle (§8.4).
- The **Kp machinery accessions** in §4 were originally identified through the repo's OrthoFinder table. They
  have since been re-derived independently by pairwise sequence alignment: all 16 machinery pairs at
  **76–99% identity** against a **20–32% random-pair baseline** (ClpX 98.8%, ClpP 99.0%, FtsH 98.0%, Lon 98.4%,
  ClpA 96.8%). The assignments stand, but they no longer depend on the ortholog table being correct.

---

## 10. Project documents and the activated-ClpP layer (2026-08-06 pass)

Everything in §1–§9 is external literature. This section adds (a) the **Gr-ADI project documents**, which are
not literature and are the authority on what the axis is *for*, and (b) the two datasets that measure the
criterion the project actually gates on.

### 10.1 Project documents — the authority on scope

Not citable externally; recorded because they settle questions the literature cannot.

| Document | Where | Why it matters |
| --- | --- | --- |
| **Gr-ADI research vision, v5 Final** — "Exploring BacPROTACs as a new paradigm for antibacterial discovery" | Drive `SecondRound/ProposalSections/Gr-ADI_ResearchVision_TPD-vs-Gram-Negs_Draft-v5(Final)_MDF.docx` | **Names the protease: ClpP, activated, partnerless.** Rules out ClpC for Gram-negatives. WP1 criteria (a)–(c); WP2_ClpPELs; first targets DnaK + AcpP; second GyrA/GyrB; Objective 3 = carbapenemases pre-export. |
| Gr-ADI proposal, v3 draft | Drive `FirstRound/TPD-vs-Gram-Negs_Draft-v3.docx` | Superseded. Its criterion (c) is the *partnerless* wording; its first targets were **FabI** + DnaK. Useful only for the v3→v5 diff. |
| **Ersilia Scope of Work**, 26 Mar 2026 | Drive `SoW/260326_GrADI_SoW_Ersilia` | Ersilia leads the WP1 target-selection workflow. Names Conlon 2013 + Jacques 2020 as the *only* available criterion-(c) datasets and states the team "would be keen to obtain similar datasets for *E. coli* and *K. pneumoniae*". **Draft — its WP1.1 paragraph is truncated mid-sentence.** |
| Kick-off notes, 24 Apr 2026; team notes 27 Jul 2026 | Drive `GrADI_Notes` | "targets in the periplasm"; "degradability > vulnerability? Are there AA motifs…"; Edkins: "partly disordered and partly ordered (two domains)". |
| Munich meeting, Day 2, 16 Apr 2026 | Drive `Docs/2026-04-16_Gr-ADI Munich Meeting Presentation_Day 2.pdf` | Kim Lewis's ~30 cell-envelope targets, incl. **FtsH as a target**. |
| Strauss, RSWVF\R2\262048 (TELL/ME) | Drive `Docs/RSWVF_Strauss-ResearchProposal.pdf` | A **different** grant, covering the target-engaging-ligand end only. It never names a protease — which is what misled an earlier draft of our report into calling the workstream "protease-agnostic". Do not read it as the Gr-ADI scope. |

### 10.2 Track 3.3c — activated-ClpP (partnerless) proteomics

The only measurements of "degraded by activated ClpP in the absence of an unfoldase partner". Both *S. aureus*;
no equivalent exists for E. coli or K. pneumoniae.

- **Conlon BP, Nandy S, Ashraf S, et al.** *Activated ClpP kills persisters and eradicates a chronic biofilm
  infection.* **Nature** 503:365–370 (2013). DOI `10.1038/nature12790` · PMID 24226776 · PMC4031760.
  ADEP4 vs untreated MRSA, iTRAQ. **Table S1** = 1,712 proteins, fully tryptic → abundance (`Average` =
  log2(ADEP4/control), with an adjusted p-value). **Table S2** = 2,382 peptide rows over 631 proteins,
  *partially* tryptic (one terminus generated by an endogenous protease) → **direct cleavage products**. Keyed
  by SACOL locus tag + `YP_` RefSeq accessions.
  **Access:** the article is paywalled and Europe PMC answers `"Article with id PMC4031760 is not open access
  one"`; `static-content.springer.com` ESM paths answer **403** to bare curl; PMC now fronts the article with
  **reCAPTCHA**. But the workbook is **freely served** from
  `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature12790/MediaObjects/41586_2013_BFnature12790_MOESM97_ESM.xlsx`
  — note the host. Now fetched automatically by `10a` (`conlon2013_adep4_saureus`).
- **Jacques S, van der Sloot AM, Huard S, et al.** *Imipridone Anticancer Compounds Ectopically Activate the
  ClpP Protease and Represent a New Scaffold for Antibiotic Development.* **Genetics** 214:1103–1120 (2020).
  DOI `10.1534/genetics.120.303153` · PMID 32094149 · PMC7153937.
  30 µM ONC212 vs untreated, label-free LC-MS/MS-FAIMS, 10 min / 40 min / 24 h. **Table S3** = 1,620 *S. aureus*
  proteins with `24H_log2_fold-change` (abundance) and
  `10-40_minutes_non-tryptic_peptides_log2_fold-change` (cleavage). **Carries no p-values.** Keyed by GenBank
  `ODV*`, with `Mu50_homolog` (UniProt) and `COL_homolog` (SACOL) columns — the latter makes SACOL the join key
  between the two papers.
  **Access:** the Europe PMC `supplementaryFiles` endpoint for PMC7153937 returns **only figure images**. Data
  is on the GSA figshare deposit **`10.25386/genetics.11873841`** (article id 11873841); Table S3 is figshare
  file `21767271`. figshare's `ndownloader` answers HTTP **202 with an empty body** while preparing a file, so a
  naive `raise_for_status()` banks an empty download as success. Raw MS: PRIDE **PXD016119**.
  Table S2 (human NALM-6) is deliberately *not* ingested — it is the human-mitochondrial-ClpP off-target
  picture, not target evidence.

### 10.3 ClpP activator chemistry named by the proposal

Cited here because Ersilia's *other* WP2 deliverable (the ClpP/activator pocket assessment) rests on them, and
because the activator pocket is the ClpX/ClpA docking cleft.

- **Brötz-Oesterhelt H, et al.** *Dysregulation of bacterial proteolytic machinery by a new class of
  antibiotics.* **Nat Med** 11:1082–1087 (2005). DOI `10.1038/nm1306`. The ADEPs.
- **Li DHS, et al.** *Acyldepsipeptide Antibiotics Induce the Formation of a Structured Axial Channel in ClpP:
  A Model for the **ClpX/ClpA-Bound State** of ClpP.* **Chem Biol** 17:959–969 (2010).
  DOI `10.1016/j.chembiol.2010.07.008`. The title is the point: the activator site is the unfoldase interface.
- **Sass P, et al.** *Antibiotic acyldepsipeptides activate ClpP peptidase to degrade the cell division protein
  FtsZ.* **PNAS** 108:17474–17479 (2011). DOI `10.1073/pnas.1110385108`. The canonical partnerless substrate.
- **Leung E, et al.** *Activators of Cylindrical Proteases as Antimicrobials…* **Chem Biol** 18:1167–1178
  (2011). DOI `10.1016/j.chembiol.2011.07.023`. The ACP series; ACP6-12 activates *E. coli* ClpP, K_d ≈ 0.27 µM.
- **Wei B, et al.** *Anti-infective therapy using species-specific activators of Staphylococcus aureus ClpP.*
  **Nat Commun** 13:6909 (2022). DOI `10.1038/s41467-022-34753-0`. (R)-ZG197.
- **Lin F, et al.** *Structure-Based Design and Development of Phosphine Oxides as a Novel Chemotype for
  Antibiotics that Dysregulate Bacterial ClpP Proteases.* **J Med Chem** 67:15131–15147 (2024).
  DOI `10.1021/acs.jmedchem.4c00773`.
- **Zhang T, et al.** *Structure-Guided Development of ClpP Agonists with Potent Therapeutic Activities against
  Staphylococcus aureus Infection.* **J Med Chem** 68:1810–1823 (2025). DOI `10.1021/acs.jmedchem.4c02562`.
  The ZG/ZY series.
- **Brötz-Oesterhelt H & Vorbach A.** *Reprogramming of the Caseinolytic Protease by ADEP Antibiotics…*
  **Front Mol Biosci** 8:690902 (2021). DOI `10.3389/fmolb.2021.690902`. Review; the compressed→extended
  conformational switch.
- **Comajuncosa-Creus A, Jorba G, Barril X & Aloy P.** *Comprehensive detection and characterization of human
  druggable pockets through binding site descriptors.* **Nat Commun** 15:7917 (2024).
  DOI `10.1038/s41467-024-52146-3`. The method the proposal names for the ClpP pocket comparison.

### 10.4 Still to verify — do not cite until checked

- **Jarzab A, Kurzawa N, Hopf T, et al.** *Meltome atlas — thermal proteome stability across the tree of life.*
  **Nat Methods** 17:495–503 (2020). DOI `10.1038/s41592-020-0801-4`. Reported to include *E. coli*; this is the
  only candidate route to proteome-scale conformational stability, which Regime B needs. **E. coli coverage is
  unverified — check the actual supplementary file before relying on it.**
