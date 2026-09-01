# Degradability — dataset & tool inventory, and the column schema

Companion to `03_degradability.md` (spec), `degradability_report.md` (audit + Gr-ADI context) and
`degradability_references.md` (citations). **This file answers three operational questions the others do not:
which datasets exist, how to process each one, and which columns should end up in the output table.**

Written **2026-08-06**. Every `n` is quoted from `degradability_references.md` (which records it from the paper);
every download route was HEAD- or GET-probed on that date and is marked ✅ verified or ⚠ unverified. Coverage
figures marked **✓** were measured in this repo, not estimated.

Reference proteomes: **Kp HS11286 `UP000007841` = 5,728** · **Ec K-12 `UP000000625` = 4,403**.

---

## 0. The three facts that govern everything below

1. **The handle is partnerless activated ClpP** (see `03_degradability.md`, protease-decision block). So the
   feature layer that matters is *Regime B* — disorder, conformational stability, accessibility — not the
   initiation-region rule, which is a *Regime A* (unfoldase) rule.
2. **Two bars, never merged**: `clpP_complex_substrate` (the written v5 criterion, well populated) and
   `partnerless_clpP` (the assay that gates progression, sparse), plus `bar_disagreement`.
3. **Almost every experimental dataset in this field is *E. coli*.** Ec gets direct evidence on nearly every
   column; **Kp inherits a shadow capped at 3,179/5,728 = 55.5%**. That asymmetry is why §4 gives two schemas.

---

## 1. Turnover and protease attribution — the label layer

| Dataset | What it gives | n | Key | Access | Kp | Effort |
|---|---|---|---|---|---|---|
| **Gupta 2024** *Nat Commun* 15:5890 | **The primary label set.** Half-lives, **~3,200 proteins × 13 conditions**, plus a protease-KO panel (Δ*clpP*/Δ*lon*/Δ*hslV*/triple/Δ*smpB* (**no Δ*ftsH*** — see the 2026-08-06 correction)) assigning substrates to 6 categories (ClpP 64, Lon 14, HslV 1, additive 82, redundant 41, ~100 still degrading in the triple KO). Supp. Data 4 = ~600 measured in-vivo N-termini. | ~3,200 | UniProt | ✅ Zenodo `10895828` (single 51.5 MB zip) + GitHub `wuhrlab/ProteinTurnoverEcoli`; PRIDE `PXD042444` ✅ | transfer | 1 d |
| **Nagar 2021** *mSystems* 6:e01296-20 | Pulsed-SILAC half-lives (741 stable / 334 slow / 72 fast) **plus a ready-made 188-feature matrix** (4 disorder features, N-end rule, known degrons, pI, MW, GRAVY, 128 node2vec PPI features). | 1,149 | UniProt | ✅ already fetched (`10a`) | transfer | **done** |
| **MacKrell 2026** *PNAS* 123 | BONCAT+TMT, time-resolved. **The only source covering stationary phase.** 88 / 56 pronouncedly unstable. Supp. S4 = ML predictions for unmeasured proteins. | 1,810 exp / 1,339 stat | UniProt | ✅ PRIDE `PXD062881` | transfer | half d |
| **Flynn 2003** *Mol Cell* 11:671 | ClpXP substrate trap **+ the five CM/NM recognition consensuses**. ~90% of trapped proteins carry a candidate signal. Per-motif: N-M1 18 · N-M2 12 · N-M3 10. Plus validated C-terminal peptides (see §3.1). | 60 | gene symbol | ✅ **in hand** — `data/raw/other/degradability/literature/PIIS1097276503000601.pdf` | transfer | see §5 |
| **Neher 2006** *Mol Cell* 22:193 | ClpXP trap ± DNA damage, **SILAC-quantified**, with a *clpX*⁻ control (Table S1) to subtract false positives. Half the set moves >3× after damage; **25% of the SOS regulon** are substrates. | ~100 | gene symbol | ✅ **in hand** — `data/raw/other/degradability/literature/PIIS1097276506001687.pdf` | transfer | see §5 |
| Niwa 2022 *Molecules* 27:3772 | Best public **Lon vs ClpXP vs HslUV** side-by-side; Lon dominates the ~80 obligate GroE substrates. | ~80 | gene symbol | ✅ open, PMC9228906 | transfer | half d |
| Westphal 2012 · Arends 2016 | FtsH substrate traps — 15 putative / 4 validated, then **>50** across growth phases. | ~50 | gene symbol | ✅ / ⚠ Wiley | transfer | low value* |
| Burgos 2020 *Mol Syst Biol* 16 | *M. pneumoniae* Lon/FtsH conditional depletion — 62 + 34 candidates and half-lives. | 96 | gene symbol | ✅ open, PMC7737663 | caution | low value* |

\* **Low value under the ClpP decision** — FtsH and Lon are Regime A channels the project is not recruiting.
Keep as provenance and cross-check, not scored evidence.

---

## 2. Biophysics — the Regime B feature layer

The block the earlier write-up wrongly said we had nothing for. It exists, it is substantial, and it is E. coli.

| Dataset | What it gives | n | Key | Access | Effort |
|---|---|---|---|---|---|
| **Cappelletti 2021 LiP-MS** *Cell* 184:545 | **The most mechanistically apt single column available** — limited proteolysis measures exactly the accessible/flexible-region property engagement requires. Per-protein protease-accessibility density (LiP peptides ÷ tryptic count), half-tryptic fraction, and count of significantly changing peptides = local conformational plasticity. 8 carbon sources, peptide resolution. **Same strain (BW25113) as Mateus 2018 ⇒ clean join.** | ~1,900 | UniProt | ✅ PRIDE `PXD022297` | 1 d |
| **Mateus 2018 TPP** *Mol Syst Biol* 14:e8242 | The canonical Ec melting-temperature reference — **1,738 fitted apparent T_m** (Dataset EV1). **Must-handle caveat: T_m follows a cell-surface → cytoplasm high-to-low gradient, so it partly encodes compartment and must be residualised on localization before use.** | 1,738 | UniProt | ✅ CC-BY PMC6056769; PRIDE `PXD009495` | half d |
| **Mateus 2020** *Nature* 588:473 | TPP across **121 strains**. The derived feature worth extracting is **stability-score variance across the 121 perturbations** — a conformational-plasticity measure arguably more informative than a single T_m. | 1,764 | UniProt | ✅ PRIDE `PXD016589`; `ecoliTPP.shiny.embl.de` ✅; GitHub `fstein/EcoliTPP` | 1 d |
| **eSOL** Niwa 2009 *PNAS* 106:4201 **+ Niwa 2012 chaperone data** | Chaperone-free solubility (% soluble), PURE system, strongly bimodal — **and, verified 2026-08-06, the same CSV also carries Niwa 2012's chaperone experiment**: `Minus Sol (%)` / `TF Sol (%)` / `GroE Sol (%)` / `KJE Sol (%)` (+ matching yields) for **exactly 788 proteins**, matching the paper's own "788 proteins × 4 conditions". **So Niwa 2012 needs no separate supplementary fetch.** 4,132 rows total; `Solubility (%)` populated for 3,173. | 4,132 rows · 3,173 solubility · **788 chaperone** | `JW_ID` / `B number` / gene name | ✅ `dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip` (197 KB, verified) | 2 h |
| **To 2021** *JACS* 143:11435 | **396 of 1,198 (33%) non-refoldable** after 2 h — a kinetic-stability / metastability proxy complementary to T_m. Protein- *and* domain-level tables. | 1,198 | UniProt | ✅ ACS-gated, but NSF-PAR PDF `par.nsf.gov/servlets/purl/10311595` verified 200 | half d |
| To 2022 | Chaperone-*assisted* refolding in a cytosol-like milieu — separates intrinsic from chaperone-rescued refoldability. | ~1,200 | UniProt | ✅ Europe PMC `PMC9704704` (0.4 MB) | half d |
| **Calloni 2012** *Cell Rep* 1:251 | **DnaK interactome, 8 sheets.** Table S1 = the interactor set (**688 rows**; paper cites 674 high-confidence) and also carries **`GroEL class`**, **`Oligomeric state`**, `Sol (%)`, `Localization`, `TMD`, `SP`, `SCOP fold`, `PD/BG ratio`. Table S8 = 1,133 proteins aggregating in the mutant insoluble fraction. | 688 (+7 more sheets) | `EG` (EcoGene) + `Gene Name` — **no UniProt**, join via `gene_aliases_to_uniprot` | ✅ **in hand, supplied manually** — `data/raw/ecoli/degradability/calloni2012_dnak/` (ScienceDirect `mmc*` answers 403, no automated route). Legacy BIFF: needs **`xlrd>=2.0.1`**, not openpyxl; **header row differs per sheet** (S1 at 7, S8 at 8) | half d |
| Györkei 2022 *Sci Rep* 12:6547 | In-vivo solubility limits — per-protein solubility threshold + a 3-class label {soluble, rapidly aggregating, slowly aggregating}, ASKA GFP library. | 2,577 cytosolic | b-number | ✅ Europe PMC `PMC9023497` (3 MB) | half d |
| Meltome Atlas Jarzab 2020 | Cross-species T_m prior, 13 species. **For E. coli it largely duplicates Mateus 2018**, so it is a transfer sanity-check, not a primary source. | ~48,000 | mixed | ⚠ **per-protein file still unverified** (identified by size, not content). Browser ✅ 200 | skip for now |
| PHDB Ramakrishnan 2020 | Would be a one-stop shortcut: **4,305 Ec proteins × >100 proteostatic parameters** aggregating 13 whole-proteome chaperone studies + TANGO aggregation + IUPred disorder. | 4,305 | UniProt | ⚠ **`http://phdb.switchlab.org` only — `https://` fails.** Web UI; **no documented bulk export** | unknown |

**If you build only two things from this section:** LiP-MS (accessibility) and Mateus 2018 T_m residualised on
localization.

---

## 3. Recognition rules — implement, don't join

Different processing strategy: the work is writing and validating a matcher, not merging a table.

| Rule | Consensus / model | How measured | Bar | Strategy |
|---|---|---|---|---|
| **DEtox PSSM** Beardslee & Schmitz 2024 | `[L/F/Y/W]-[L/A/R]-L-A-A` at the C-terminus; terminal Ala-Ala in **91%** of the top 100; best tag `FKLVA` 639× enriched | VapC bearing **~100,000 random C-pentapeptides**, NGS in WT / Δ*clpX* / Δ*clpP*. **~1% of random 5-mers are functional degrons**; ClpXP overwhelmingly dominates | ClpXP (selection) | Implement as a **PSSM, not a regex**; score all proteins, don't threshold. Partly present as CM1 |
| **Lon C-degron** Cragan 2025 | `x–[L/I]–[L/I/V]–H-COOH`, invariantly His-terminal; H→D collapses degradation to ≤0.3 min⁻¹ | Biochemistry + K_D series 0.24 → >50 µM. **Conserved across E. coli / Y. pestis / M. pneumoniae** — licenses Kp transfer | Lon (Regime A) | Cheap regex, high specificity. Low priority under the ClpP decision |
| **N-FIVE** Sen 2025 | Extends the N-degron rule from P1 to **P1–P5**. P1 potency `F≈R < L < W≈K≈Y`; Pro@P2 **rescues** (+1.05); Gln@P2 destabilizes as much as Leu@P1; **net negative charge over P2–P5 stabilizes** | **~2.2 M variants**, FACS+NGS, protein-stability-index. All effects vanish in Δ*clpS*; the Arg/Lys effect vanishes in Δ*aat* | ClpAP+ClpS (selection) | **Use their model + code** (✅ GitHub `KunjapurLab/N-terminal-cluster-stability`), gated on Met-excision |
| **Won 2024** — the only validated one | Lasso over 485 descriptors. **Top feature = mean disorder of the N-terminal ~30 aa** (flDPnn, threshold 0.2); C-terminal 30 aa showed **no** correlation. Abundance anti-correlates at **r = −0.69** | 72 native *M. smegmatis* proteins, rapamycin-induced ClpC1P1P2 degradation; fitted on 54, held out r > 0.6 on 18 | ClpC1 — **transfer caution** | **Re-fit, don't reuse coefficients** (ClpC1/Actinobacterium; flDPnn is web-only so feature parity is unobtainable) |

> **Caveat spanning all four.** Every rule here describes how an *unfoldase* selects substrates.
> **Partnerless activated ClpP has no sequence degron at all** — the ADEP mechanism bypasses recognition. And
> two independent proteome-scale studies (Nagar 2021; Gupta 2024) found **no association between the classic
> motifs and measured half-life**; Gupta specifically reports rapidly degrading proteins show *no* enrichment for
> destabilizing N-terminal residues. Implement for completeness and provenance; do not expect them to carry the
> score. Only ~1% of either proteome carries a weighted motif (§5 of `degradability_report.md`).

### 3.1 Flynn's real consensuses — and a correction to our own code

`src/degradability.py` drops Flynn's N-motifs, arguing NM1 "fires on 30.7% of the Kp proteome" and citing
`^M?[AILVMFW]{2,}`. **That regex is not NM1.** Read from the paper (now in `data/raw/other/degradability/literature/`):

- **N-M1** strict consensus **`T₁-X₂-K₃-[ILV]₄`, located 1–4 residues from the N terminus**; looser form
  `polar-T/φ-φ-basic-φ`. 18 trapped proteins (σO, Dps + 16).
- **N-M3** `φ-X-polar-X-polar-X-basic-polar`. 10 proteins (Crl, DksA + 8).
- **N-M2** 12 proteins (DadA, IscS, OmpA + 9).
- **C-M1** ssrA-like; critical region terminal `-Ala-Ala`, enriched **7-fold**. **C-M2** MuA-like.

These are **far more specific** than the legacy pattern, so the stated reason for dropping NM1/NM3 does not
survive contact with the paper. They should be re-implemented from the real consensus and re-measured against the
Nagar/Gupta labels — the same archetype-assertion discipline `selftest_motifs()` already enforces for CM1.

**Ready-made positive-control set** (validated C-terminal peptides, Flynn Fig. 3/6): ssrA `AANDENYALAA` ·
YdaM `KNDGRNRVLAA` · Crl `DFRDEPVKLTA` · LldD `ALAPMAKGNAA` · MuA `ILEQNRRKKAI` · YbaQ `ARREERAKKVA` ·
RibB `AYRQAHERKAS` · Gcp `RWPLAELPAA` · Dps `FLWFIESNIE` (a negative for the -AA rule).

---

## 4. Computational tools

Formatted like `localization_log.md` §Sources — the tool, what it gives, and the quirk that will bite you.

| Tool | What it gives | Availability / env | Have it? | Verdict |
|---|---|---|---|---|
| **MobiDB bulk** | **Verified 2026-08-06, both proteomes, 100% of accessions.** Kp 6.31 MB / 85,234 rows / **5,728 accs**; Ec 5.15 MB / 73,479 rows / **4,404 accs**. One file contains: **8 disorder predictors** (dis465, disHL, glo, th_50, espX, iups, espN, iupl) + a `priority` consensus; **AlphaFold pLDDT**; AlphaFold-derived disorder; **Pfam / Gene3D / merged domain boundaries**; LIP (linear interacting peptide) regions; low-complexity; TM and signal-peptide tracks — each with `content_fraction` and `content_count`. | ✅ `mobidb.org/api/download?proteome=<UPID>&format=tsv`. **⚠ HEAD returns 405 — use GET.** No env, no install | 1 `curl` | **Do this first.** Supplies the global disorder fraction *and* the domain boundaries for `two_domain_architecture` — the two things previously called the highest-value remaining build |
| AlphaFold pLDDT (local) | Per-residue confidence; terminal exposure + initiation-region run length. Computed by `10b` | `gradi`; `data/processed/<org>/pockets/pdb/` | ✅ 5,727/5,728 Kp | Primary structural feature. **Caveat: misreads obligate complex subunits** — ordered in situ, disordered alone |
| Chainsaw · Merizo · SWORD2 | Domain segmentation from structure | pip/conda, new env | ✗ | **Probably unnecessary** — MobiDB domain calls cover **4,829 Kp / 4,053 Ec**. Only if those prove too coarse |
| metapredict · AIUPred | Local sequence-side disorder | pip; GitHub | ✗ | Optional; MobiDB already gives 8 predictors + consensus. Useful for the 1 Kp protein with no AlphaFold model |
| flDPnn | The predictor Won 2024 used (threshold 0.2) | **web server only** | ✗ | **Do not chase** — parity unobtainable; re-fit Won instead |
| TMbed · SignalP 6.0 · DeepLocPro · PSORTb | Topology, signal-peptide/lipobox typing, prokaryote localization. **Signal-peptide status is now a positive feature** (pre-export window), so this block is load-bearing | `gradi-loc` (must stay separate from `gradi` — `fair-esm` collides with EvolutionaryScale `esm`) | ✅ 100% both organisms | Reuse as-is |
| AlphaKnot 2.0 | Knotting of AlphaFold models; a knot tightened against a stable domain can stall translocation (ATP cost ↑30–1000×) | web + database | ✗ | **Soft penalty, never a veto** — Sivertsson 2019 showed ClpXP *can* degrade 3₁/5₂-knotted proteins. Low priority |
| RCSB `oligomeric_state` + PDBe interfaces | Assembly state; a homo-oligomer buries its termini (GroEL = worked negative control) | REST; pattern exists in the ligandability axis | partial | Do it for Ec (**1,779** with experimental stoichiometry); Kp has only **30**, so fall back to Schweke 2024 homomer prior or UniProt subunit text |
| `poppler` (`pdftotext`/`pdftoppm`) | PDF table extraction; also what the Read tool needs to render PDFs | `brew install poppler` | ✗ **missing** | One line, low risk. Add to `install.sh` — see §5 |

**Total new environments needed: zero. Total new downloads to bring the Regime B feature layer to full
coverage: one.**

---

## 5. Processing strategy for PDF-only sources

`data/raw/other/degradability/literature/` now holds Flynn 2003, Neher 2006 and Ziemski 2021 Table S1. Together Flynn (60) + Neher (~100)
take the Ec ClpXP selection bar from a **45-row hand-curated stand-in** to roughly **130–160 real proteins**.

1. Text extraction needs **nothing installed** — a stdlib `zlib` pass over `stream…endstream` objects works
   (~68 kB of text from Flynn). **But table structure is lost.**
2. So do **not** attempt positional table parsing. **Vocabulary-match instead**: extract candidate tokens and
   match against the Ec gene-symbol vocabulary from `L.load_genes('ecoli')` (or the proteome TSV), then resolve
   to UniProt accessions. **Report ambiguities; never guess.**
3. `brew install poppler` makes this materially easier and unblocks PDF rendering.
4. **Ziemski stays deprioritised** — ClpC-family, capped at `CLPC_ONLY_CAP = 0.45`, describes a machine neither
   organism has. Recorded as *in hand*, not *to fetch*.

---

## 6. Orthology transfer — what it buys, and what it costs

| Regime | Route | Measured ceiling |
|---|---|---|
| **direct** | Ec-native, or **recomputed** from Kp sequence/structure | full coverage |
| **transfer** | Ec → Kp by ortholog | **3,179 / 5,728 = 55.5%** ✓ |
| **caution** | cross-phylum, or describes an absent machine | **10.6% Kp / 13.8% Ec** ✓ (the S. aureus activator data) |

**Safe to transfer:** sequence-derived rules and structure-derived features — don't transfer them at all,
**recompute** on the Kp sequence/model for 100% coverage. Conserved biochemistry (the Lon degron, explicitly
conserved across three genera) transfers as a *rule* even though the data does not.

**Not safe:** measured per-protein values (half-life, T_m, solubility) are the ones capped at 55.5%, and the cap
is on *which proteins*, not on accuracy. Two known limitations: `transfer_ecoli_to_kp` reduces with max/mean,
which is **meaningless for a class label** — reuse `transfer_categorical_ecoli_to_kp` (donor consensus,
**abstains on ties**) as localization does; and the Kp↔Ec orthology is OrthoFinder orthogroup membership only,
carrying **no `pident`/`coverage`/`bitscore`**, so there is no identity to threshold on.

---

## 7. Column schema — *E. coli*

Conventions the schema enforces: (1) the two bars never merge; (2) missing ≠ zero — absent sub-scores are `NaN`
and remaining weights renormalise (`07h` idiom, not `06g` zero-fill), with a `*_confidence` = retained weight
mass; (3) evidence and prediction stay separate columns, and `insufficient_evidence` is a first-class tier
distinct from `low`; (4) provenance travels with the value (donor accession + method); (5) gates multiply, they
do not score.

| Column | Source | Type | Coverage | In score? |
|---|---|---|---|---|
| **Regime B — partnerless activated ClpP (gates progression)** | | | | |
| `disorder_fraction` | MobiDB `prediction-disorder-priority` | prediction | 100% ✓ | yes — primary |
| `two_domain_architecture` | MobiDB domain ∧ disorder | prediction | 2,000 ✓ (45%) | yes — primary |
| `cterm_init_region_len` · `nterm_*` · exposure | AlphaFold pLDDT (`10b`) | prediction | 4,326 ✓ (98%) | yes |
| `lip_accessibility` | Cappelletti 2021 | **evidence** | ~1,900 (43%) | yes |
| `tm_apparent` · `tm_residualised` | Mateus 2018 | **evidence** | 1,738 (39%) | yes |
| `stability_variance` | Mateus 2020 | **evidence** | 1,764 (40%) | tie-break |
| `non_refoldable` | To 2021 | **evidence** | 1,198 (27%) | tie-break |
| `solubility_esol` | eSOL | **evidence** | 3,198 (73%) | tie-break |
| `activator_evidence` · `_cleaved` | Conlon + Jacques (`10c`) | **evidence** | 609 ✓ (13.8%) | veto/confirm only |
| **Regime A — ClpP-complex substrate (the written criterion)** | | | | |
| `degron_cm1_pssm` | DEtox | rule | 100% | yes (low weight) |
| `degron_nm1` · `_nm3` | **Flynn 2003 real consensus** (§3.1) | rule | 100% | re-measure first |
| `ndegron_nfive` | Sen 2025 model + code | rule | 100% | yes (low weight) |
| `clpxp_trapped` | **Flynn 60 + Neher ~100** | **evidence** | ~140 (3%) | confirm only |
| `halflife_min` · `kdeg_growth_corrected` | Gupta 2024 + Nagar 2021 | **evidence** | ~3,200 (73%) | yes |
| `protease_attribution` | Gupta KO panel | **evidence** | ~200 (5%) | **label, not feature** |
| **Gates · modifiers · composite** | | | | |
| `localization` · `clp_accessibility` | `09g` | mixed | 100% ✓ | gate (×) |
| `compartment_handle` · `preexport_window` | `09e`/`09d` + the router | mixed | 100% ✓ | router |
| `assembly_state` | RCSB `oligomeric_state` | **evidence** | 1,779 (40%) | penalty |
| `resynthesis_burden` | Schmidt 2016, molecules/cell | **evidence** | 2,359 (54%) | **shortlist filter only** |
| `partnerless_clpP` · `clpP_complex_substrate` · `bar_disagreement` · `*_confidence` · `degradability_tier` | `10i` merge | composite | 100% | **the output** |

---

## 8. Column schema — *K. pneumoniae*

Same rows, so the two tables read as a diff.

| Column | Route to a Kp value | Coverage | Status |
|---|---|---|---|
| **Regime B — recomputed natively, so Kp is NOT disadvantaged here** | | | |
| `disorder_fraction` | **Recompute** — MobiDB Kp bulk | 5,728 ✓ 100% | direct |
| `two_domain_architecture` | **Recompute** — MobiDB domain ∧ disorder | 2,644 ✓ (46%) | direct |
| `cterm_init_region_len` · `nterm_*` | **Recompute** — Kp AlphaFold models | 5,727 ✓ 100% | direct |
| `lip_accessibility` | transfer from Ec LiP-MS | ~1,050 (18%) | transferred |
| `tm_apparent` | transfer from Mateus 2018 | ~965 (17%) | transferred |
| `stability_variance` · `non_refoldable` | transfer | ~980 · ~665 (17% · 12%) | transferred, thin |
| `solubility_esol` | transfer from eSOL | ~1,775 (31%) | transferred |
| `activator_evidence` · `_cleaved` | cross-phylum RBH (`10c`) | 608 ✓ (10.6%) | veto/confirm only |
| **Regime A — rules recompute; measurements do not** | | | |
| `degron_cm1_pssm` · `_nm1` · `_nm3` · `ndegron_nfive` | **Recompute on Kp sequence** | 5,728 ✓ 100% | direct |
| `clpxp_trapped` | transfer from Flynn + Neher | ~78 (1.4%) | **too thin to rank** |
| `halflife_min` · `kdeg_*` | transfer from Gupta + Nagar | ~1,775 (31%) | transferred |
| `protease_attribution` | transfer from Gupta KO panel | ~110 (1.9%) | label only |
| **Gates · modifiers · composite** | | | |
| `localization` · `clp_accessibility` · `compartment_handle` | **Native** `09g` (incl. STEPdb transferred at `09f`) | 5,728 ✓ 100% | direct |
| `assembly_state` | experimental stoichiometry — **Kp has only 30**; fall back to Schweke prior / UniProt subunit text | 30 (0.5%) | **near-empty** |
| `resynthesis_burden` | Illenseher 2025 Kp proteome (>2,800) — **relative, not molecules/cell** | ~2,800 (49%) | rescaled |
| composite columns | `10i`, same contract | 100% | the output |

**Transferred figures are the 55.5% ortholog ceiling applied to each source's own coverage — treat them as upper
bounds.** Kp is fully competitive on everything *recomputable* and thin on everything *measured*.
**Practical consequence: rank Kp on the recomputable columns; use the transferred ones only to confirm or veto.
Do not let a 31%-coverage half-life column drive a Kp ranking.**

---

## 9. Minimum viable set — the twelve columns to build first

Step numbers refer to the build order in `degradability_report.md` / the deck slide 25.

| # | Column | From | Step | Effort |
|---|---|---|---|---|
| 1 | `disorder_fraction` | MobiDB bulk — one `curl`, 100% both organisms | 1 | 2 h |
| 2 | `two_domain_architecture` | MobiDB domain ∧ disorder — **Edkins's criterion made computable** | 1 | 3 h |
| 3 | termini + exposure | `10b` | — | **done** |
| 4 | `activator_evidence` · `_cleaved` | `10c` | — | **done** |
| 5 | `localization` · `clp_accessibility` | `09g` | — | **done** |
| 6 | `compartment_handle` · `preexport_window` | `09e`/`09d` + the router | 4 | 1 d |
| 7 | `halflife_min` · `kdeg_growth_corrected` | Gupta 2024 (Zenodo) + Nagar 2021 | 6 | 1 d |
| 8 | `lip_accessibility` | Cappelletti, PRIDE `PXD022297` | 6 | 1 d |
| 9 | `tm_residualised` | Mateus 2018 EV1, residualised on localization | 6 | half d |
| 10 | `degron_cm1_pssm` + corrected `_nm1`/`_nm3` | DEtox + Flynn's real consensus (§3.1) | 5 | half d |
| 11 | `resynthesis_burden` | Schmidt 2016 (Ec) · Illenseher 2025 (Kp) | 6 | half d |
| 12 | the composite + `bar_disagreement` + tiers | `10i`, renormalise-don't-zero-fill | 2 | 1 d |

**Then validate before visualising** — correlate against Gupta/Nagar measured turnover across ~3,200 Ec proteins.
If that fails, the construction is wrong. Retire the Ec MD5 mock and the frozen Kp legacy file in the same commit
as the merge.

---

## 10. Validation: what we can and cannot test this against

Measured 2026-08-06. **Headline: the axis currently has no usable validation set for the bar that gates
progression, and that is now a measurement rather than a worry.**

### 10.1 The candidate label sets, and what each is a label *for*

| Label set | n | What it measures | Right bar? | Usable as validation? |
|---|---|---|---|---|
| **Gupta 2024** | 3,262 Ec | *natural* in-vivo half-life, 13 conditions, + which *endogenous* protease | ✗ | wrong quantity — §10.2 |
| **Nagar 2021** | 1,149 Ec (743 stable / 334 int / **72 fast**) | *natural* turnover, 3 classes | ✗ | this is what the `10d` AUROC used |
| **MacKrell 2026** | 1,810 exp / 1,339 stat | *natural* turnover, independent replicate | ✗ | same category |
| **Conlon + Jacques** | 1,949 Sa → **609 Ec / 608 Kp** | **activated partnerless ClpP** — right mechanism | ✓ | **already the FEATURE in `10c`** → circular |
| **Won 2024** | 72 *M. smegmatis* | **induced-proximity degradation** — right *modality*, measured rate constants | modality ✓, ClpC1/Actinobacterium ✗ | the only external set that exists |
| **Gr-ADI WP1 assay** | 0 | purified target + activated ClpP, no unfoldase | ✓✓ exactly right | **does not exist yet** |

### 10.2 The measurement that settles it — the two bars are statistically independent

609 E. coli proteins carry a transferred activator score; **352** also have a Nagar half-life. On that overlap:

| Activator readout | Spearman rho vs natural half-life | p | n |
|---|---|---|---|
| `activator_evidence` | **−0.073** | **0.17** | 352 |
| ADEP4 abundance log2FC | +0.134 | 0.014 | 336 |
| ONC212 cleavage log2FC | −0.065 | 0.27 | 292 |

And `activator_evidence` by Nagar class shows **no monotonic trend**: fast (n=14) **0.483** · intermediate
(n=125) **0.540** · stable (n=213) **0.433**. Fast vs stable: **Mann–Whitney p = 0.38**.

**Consequence.** `10d` reports the disorder features at best **AUROC 0.587 [0.49, 0.68]** against Nagar, with
only 4/22 features' CIs excluding 0.5 and 20 statistically indistinguishable from the best. Since the two
quantities are independent, that AUROC is **not weak evidence about the right thing — it is evidence about a
different thing.** Read it as characterising *natural instability*, and nothing more.

> **And the direction may be actively wrong.** A good degrader target is **stable natively but degradable when
> forced** — if it already turns over fast, the drug adds nothing. So optimising a degradability score against
> natural half-life could select the *wrong* proteins. This is the engageability-vs-depletability tension from
> §"Two questions hiding in one score" landing in the validation set, and it is the reason the big turnover
> datasets cannot simply be pressed into service as labels.

### 10.3 The reproducibility ceiling, and it is low

Before asking a model to predict activator susceptibility, ask how well the **two experiments agree with each
other**. On shared *S. aureus* proteins:

| Readout | Spearman rho | n |
|---|---|---|
| abundance (ADEP4 vs ONC212) | **+0.52** | 404 |
| cleavage (ADEP4 vs ONC212) | **+0.22** — weak | 241 |
| binary ≥2×-down call | overlap **37 / 121** → **Jaccard 0.31** (both 37 · ADEP4-only 8 · ONC212-only 76) | 404 |

Two chemically distinct activators of the *same enzyme* agree on under a third of their calls. **No predictor
should be judged against 1.0 — judge against this ceiling.** The ADEP4-only 8 vs ONC212-only 76 asymmetry also
says the two differ in stringency (ONC212 at 30 µM is far more promiscuous), so "activated-ClpP substrate" is
partly a **compound-specific** property, not a clean protein property.

Two consequences for the current implementation:
- **`10c` pools across activators by taking the max**, which over two poorly-agreeing measurements inflates the
  score. Needs revisiting — an agreement-weighted or intersection-based pool would be more honest.
- The **0.65 / 0.35 cleavage-over-abundance weighting** was argued from mechanism plus one example (AcpP) and has
  never been fitted. Note that cleavage is the *less* reproducible of the two readouts (rho 0.22 vs 0.52), which
  argues against weighting it higher on reproducibility grounds even though it is better on mechanistic ones.

### 10.4 Three tests that can actually be run

| # | Test | What it establishes | Cost | Status |
|---|---|---|---|---|
| 1 | **Held-out across activators** — fit/threshold on Conlon (ADEP4), test on Jacques (ONC212), then reverse | Whether the feature set generalises across two activators of the same enzyme. The only **non-circular** use of the right-bar data. Judge against the Jaccard-0.31 / rho-0.52 ceiling | half a day | **available now, not done** |
| 2 | **Won 2024 as an external set** — rank its 72 proteins by our features vs their measured degradation constants | Whether features transfer to a genuinely induced-proximity readout. Weak, but the only external test in existence — and Won's own N-vs-C asymmetry **already failed to transfer** (`10d`: N-30 0.553 vs C-30 0.548), which is a warning | 1 day | not done |
| 3 | **Ask Gr-ADI for the first WP1 assay results as a held-out set** | **The real answer.** Even 20–30 proteins with a binary call would outweigh all 3,262 turnover labels | a conversation | the SoW anticipates it |

### 10.5 What is honestly not validated

- **The composite** — `10i` does not exist, so the two-bar score has never been tested end to end.
- **The cleavage/abundance weighting** and **max-pooling across activators** (§10.3).
- **The Kp coverage figures in §8** are arithmetic (55.5% ceiling × source coverage), not measured. Two have now
  been measured by intersecting each source's accessions with the ortholog table, and **both were wrong**:
  `halflife_min` is **2,750 (48.0%)**, not ~31%; and `solubility_esol` is **0 — not transferable at all**, because
  eSOL is b-number-keyed with no UniProt column. The estimates erred in *both* directions, so the rest should be
  measured the same way before anyone relies on them.
- **Effort estimates** throughout this file are judgement, not measurement.

### 10.6 What is validated

To be fair to the other direction — these are real and they caught real errors:

- **Provenance**, byte-level: magic bytes (BIFF vs OOXML vs HTML), OLE metadata, and **article title from the
  cached PMC XML** — the last of which caught `to2021_refoldability` pointing at an unrelated *Anal Chem* paper.
- **Internal consistency**: recomputing Conlon's `Average` from its raw intensity columns (−9.65 vs −9.45
  reported) confirmed the sign convention; eSOL's 788 chaperone rows match the paper's own count; Gupta's
  condition list contradicted our docs and **caught the spurious Δ*ftsH*** and the 13-vs-14 discrepancy.
- **Cross-source agreement**: `10d`'s disorder fraction vs the independently computed ligandability one agrees to
  **max |Δ| = 0.0010**; the top activated-ClpP hits are trigger factor plus ribosomal proteins, independently
  matching Jacques's own reported ribosome enrichment.
- **One pre-registered prediction held**: GyrA/GyrB were predicted to satisfy the written criterion and fail the
  partnerless bar *before* `10c` was run; GyrB came out at 0.0015 with no cleavage evidence.
