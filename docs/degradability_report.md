# Degradability — full technical report

**Axis 3 of the GraDi target-prioritization pipeline.** Unlike the companion reports for ligandability and
essentiality, this document describes an axis that **has not been built**. It is therefore three things at
once: an *audit* of what the repository currently serves and why it is not usable (§2), the *conceptual
argument* for what "degradable" means in an Enterobacterium (§3), and the *sourcing and design plan* for the
real pipeline (§4–§7), followed by the validation checks the build must pass and the gaps that cannot be
closed (§8–§9).

Companion files: `docs/03_degradability.md` (the axis spec + track table) and
`docs/degradability_references.md` (every citation, with DOIs). A run log
(`docs/degradability_log.md`) will be created when there are runs to log.

All repository counts in this document were measured directly from the files on disk on **2026-08-04**, and
every external URL and identifier was resolved programmatically on the same date. **§12 (2026-08-06) supersedes
parts of §1–§7** — read it first.

---

## 13. Validation — the axis has no usable label set for its own gate (2026-08-06)

Measured, not argued. Full working in **`docs/degradability_datasets.md` §10**; the three results that matter:

1. **The two bars are statistically independent.** Over the 352 E. coli proteins with both a transferred
   activator score and a Nagar half-life, `activator_evidence` vs natural half-life is **rho = −0.073,
   p = 0.17**, with no monotonic trend across Nagar's classes (fast 0.483 · intermediate 0.540 · stable 0.433;
   fast-vs-stable Mann–Whitney **p = 0.38**). So the ~3,262-protein turnover label set **cannot validate the
   `partnerless_clpP` bar**, and `10d`'s best AUROC of 0.587 is evidence about a *different quantity*.
2. **It may be worse than useless.** A good degrader target is stable natively and degradable when forced, so
   optimising against natural half-life could select the wrong proteins outright.
3. **The reproducibility ceiling is low.** ADEP4 vs ONC212 agree at rho **+0.52** on abundance, only **+0.22**
   on cleavage, and **Jaccard 0.31** on the binary ≥2×-down call. No predictor should be judged against 1.0.
   This also undermines `10c`'s max-pooling across the two activators and the untested 0.65/0.35
   cleavage-over-abundance weighting.

The only non-circular test available today is **held-out across activators** (train on one, test on the other);
the real validation set is the **Gr-ADI WP1 assay**, which does not exist yet and is worth asking for.

---

## 12. Gr-ADI project context, and the corrections it forces (2026-08-06)

Everything above §12 was written from the literature outwards, without checking the project's own documents.
Reading the Gr-ADI Drive changes the framing and retracts one of this report's conclusions.

### 12.1 The protease is specified, and it is ClpP alone

Source: `SecondRound/ProposalSections/Gr-ADI_ResearchVision_TPD-vs-Gram-Negs_Draft-v5(Final)_MDF.docx`
(the funded research vision) plus `SoW/260326_GrADI_SoW_Ersilia` (Ersilia's scope of work, a draft — its
WP1.1 paragraph is truncated mid-sentence).

- *"…reduce target protein levels by **selectively engaging the bacterial ClpP protease**."*
- *"Small molecule BacPROTACs that **activate the ClpP protease** and accelerate the rate of degradation of
  selected targets through a proximity effect…"*
- WP2 is named **WP2_ClpPELs** — ClpP-Engaging Ligands.
- ClpC is explicitly ruled out: *"The dCymM-based ligand that engages the ClpC unfoldase … **only attaches to
  the Mtb ClpC**, and it is therefore not suitable for use in developing BacPROTACs that target Gram-negative
  bacteria. However, elucidation of the MOA of the **acyldepsipeptide (ADEP)** class … showed that it acts
  through dysregulation and **direct activation of ClpP**."*
- Chemistry: **ONC212** as the first ClpPEL scaffold; an **(R)-ZG197** derivative already tested against
  **purified Kp ClpP**; ACP6-12 (Ec ClpP, K_d ≈ 0.27 µM); the ZG/ZY series (Zhang 2025).
- Targets: first **DnaK** and **AcpP**; second **GyrA/GyrB** ("clinically validated *cytosolic* targets");
  Objective 3 degrades **carbapenemases** *"prior to their translocation to the periplasm"*.
- **FabI is superseded** — a v3 target (with fabimycin as TEL), replaced by AcpP in v5.

Ersilia's other named computational deliverable, for the record, is not degradability at all: the
**structure-based ClpP/activator pocket assessment** in WP2_ClpPELs, using the pocket-descriptor method of
Comajuncosa-Creus et al. (*Nat Commun* 15:7917, 2024). Note that the pocket in question **is** the ClpX/ClpA
docking cleft, so the two deliverables share a structural object.

### 12.2 ⚠ Correction to this report: the "protease-agnostic" reading was wrong

An earlier draft of this analysis concluded that the consortium's BacPROTAC workstream was *deliberately
protease-agnostic* and that the recruiting end was an open question for us to answer. **That was wrong.** It was
drawn from Strauss's *other* grant (RSWVF\R2\262048, the TELL/ME proposal), which builds the target-engaging-
ligand end only and indeed never names a protease, and whose supporting review (Birkholtz, Olivier, Welcome &
Strauss, *Curr Opin Chem Biol* 91:102655, 2026) is Mtb/ClpC1-centred. The Gr-ADI research vision is a different
document and it names the protease outright. The reason ClpC1 pervades the consortium's material is that ClpC1
is where its published expertise lies — not that ClpC1 is in scope.

### 12.3 The mechanistic consequence: two regimes, and the spec encodes the wrong one

Activated ClpP **has no unfoldase**. ADEP-class activators bind the hydrophobic clefts where the ClpX/ClpA
(L/I)GF loops dock, flipping ClpP compressed→extended so the axial pore opens and the catalytic triad orders —
but nothing pulls. Consequently:

- **Regime A (with unfoldase — ClpXP/ClpAP/Lon/HslUV/FtsH):** a loose end suffices; the machine unfolds. This is
  the ~5 / ~20 / ~37 aa initiation-region rule that §3 of this report is built on.
- **Regime B (partnerless activated ClpP — the Gr-ADI modality):** the substrate must **already be unstructured
  enough to diffuse into an open pore** — substantial intrinsic disorder, low conformational stability, a
  nascent chain, or a natively dynamic assembly.

So the central mechanistic rule in §3 describes the regime the project is *not* using. Three things follow, and
they are now first-class rather than incidental:

1. **Disorder is close to being the whole requirement**, not one feature of nine. The remark in the kick-off
   meeting that "a disordered region may be a good indicative" was correct, for this specific reason.
2. **Edkins's criterion — "partly disordered and partly ordered (two domains)" (27 Jul 2026) — is the required
   architecture**: disorder is ClpP's entry point, the folded domain is where the TEL binds. Implement as
   `two_domain_architecture` from the per-residue pLDDT profiles already cached for 100% of both proteomes.
3. **Conformational stability enters as a first-class feature and the repo has none.** The candidate source is
   the Meltome Atlas (Jarzab et al., *Nat Methods* 2020), which reportedly includes *E. coli* — **unverified;
   check the actual supplementary file before claiming it.**

Corroboration that we are measuring the right mechanism: Jacques 2020 reports the 200 most-depleted proteins
under ONC212 are enriched for ribosome-related functions (as was true for ADEP4), and our own `10c` transfer
puts **trigger factor** and a run of ribosomal proteins at the top — the nascent-chain signature.

### 12.4 The criterion moved between drafts, so there are two bars

| | WP1 criterion (c) |
| --- | --- |
| v3 draft | "show sensitivity to degradation by activated ClpP **in the absence of an unfoldase partner**" |
| **v5 Final** | "be a **known, preferred substrate of a ClpP protease complex**." |
| validation assay (unchanged) | "the activated ClpP of the target organism **only (i.e. in the absence of an unfoldase partner)**" … "**only targets that show degradation in this assay will progress**" |

"A ClpP protease *complex*" admits ClpXP and ClpAP, so the **selection** bar is well populated (essentially all
E. coli ClpP substrate knowledge is ClpXP/ClpAP knowledge) while the **validation** bar that decides progression
has almost no data. They disagree systematically, because a good ClpXP substrate is good *precisely because* a
motor unfolds it. The composite therefore emits `clpP_complex_substrate`, `partnerless_clpP` and
`bar_disagreement` and **never averages them** — see `docs/03_degradability.md` §Composite, revision 2026-08-06.

### 12.5 New track 3.3c, built: the two activator datasets

Both fetched via the `10a` manifest and transferred by `10c_clpp_activator.py`.

| Source | Readouts | n |
| --- | --- | --- |
| **Conlon 2013** *Nature* 503:365 (ADEP4, iTRAQ, MRSA) | S1 = abundance (`Average` = log2 ADEP4/ctrl + adj. p); S2 = **partially tryptic peptides** = direct cleavage | 1,712 proteins / 2,382 peptide rows (631 proteins) |
| **Jacques 2020** *Genetics* 214:1103 (30 µM ONC212, FAIMS) | `24H_log2_fold-change` = abundance; `10-40_minutes_non-tryptic_peptides_log2_fold-change` = cleavage. **No p-values.** | 1,620 proteins |

Union on SACOL: **1,949** *S. aureus* proteins; **404** depleted ≥2× by at least one activator; **953** with
cleavage evidence. Transferred: **608 / 5,728 Kp (10.6%)** and **609 / 4,403 Ec (13.8%)**, of which 330 / 334
carry cleavage evidence. That ceiling is set by the cross-phylum RBH step and is the honest number — this is
evidence for roughly one protein in eight, not a proteome-wide annotation.

**Cleavage is weighted above abundance (0.65 / 0.35)**, because a 24 h abundance drop conflates degradation with
growth arrest, regulon change and resynthesis whereas endogenous-protease-generated peptides are a direct
product of proteolysis.

Access notes worth preserving, since two of three obvious routes look like paywalls and are not:
- **Conlon**: article paywalled, Europe PMC answers "is not open access", `static-content.springer.com` answers
  403 to bare curl, and PMC now fronts PMC4031760 with **reCAPTCHA** — but the supplementary workbook is
  **free** on `media.springernature.com/original/springer-static/esm/...`. Now fetched automatically.
- **Jacques**: the Europe PMC `supplementaryFiles` endpoint for PMC7153937 returns **only figure images**; the
  data lives solely on the GSA figshare deposit `10.25386/genetics.11873841`. figshare's `ndownloader` answers
  HTTP **202** with an empty body while it prepares a file.

Identifier route (the obvious ones fail): 2013 **SACOL** locus tags → current RefSeq has re-tagged everything
`SACOL_RS*` with **no `old_locus_tag`**, and UniProt demoted the COL/Mu50 proteomes (939 / 1,033 entries left).
Working route: paper accession (`YP_*` / `ODV*`) → **NCBI efetch** sequence → **DIAMOND RBH** → Kp/Ec UniProt,
with `pident`/`coverage` kept on every row.

### 12.6 ⚠ The proposal's stated evidence for its own targets does not fully hold

v5 states that both proteomic analyses show DnaK "reduced by at least 2-fold", and that AcpP was "reduced by at
least 2-fold upon treatment of **either** ClpP activator". Measured from the actual tables:

| Target | ADEP4 abundance | pctile | ADEP4 cleavage | ONC212 abundance | ONC212 cleavage | `partnerless_clpP` |
| --- | --- | --- | --- | --- | --- | --- |
| **AcpP** (SACOL1247) | −2.49 (p=0.004) ✓ | 3.8% | absent from S2 | **+0.03 ✗** | +3.46 ✓ | **1.00** |
| **DnaK** (SACOL1637) | −1.41 (p=0.004) ✓ | 9.0% | **25 peptides, 15 sig.**, max +4.79 | **−0.51 ✗** (1.4×) | +1.78 ✓ | **0.90** |
| GyrA (SACOL0006) | +0.19 (p=0.24, ns) | 74.1% | absent from S2 | +0.70 | +2.00 ✓ | 0.65 |
| **GyrB** (SACOL0005) | −0.01 (p=0.96) | 48.9% | 1 peptide, **0 sig.**, −0.85 | absent from S3 | — | **0.0015** |

Two conclusions, and they point in opposite directions:

1. **The quantitative wording is wrong** — under ONC212, DnaK is 1.4× (not ≥2×) and AcpP is *unchanged*
   (`+0.03`). But **both have strong ONC212 cleavage evidence**, so the biological conclusion survives by a
   different and arguably better readout. The targets are fine; the stated justification is not.
2. **GyrB has essentially no activated-ClpP evidence at all** — zero abundance change (dead-centre of the
   distribution), its one cleavage peptide negative and non-significant, and absent from the ONC212 table. It
   maps successfully (44.8% id to Kp), so this is absence of signal, not absence of mapping. GyrA/GyrB are the
   declared *second* targets, which makes this worth raising before TEL synthesis effort is committed.

This was a **pre-registered** prediction: before running `10c` the expectation recorded in the plan was that
GyrA/GyrB would satisfy the written criterion while failing the partnerless assay. That is what happened.

### 12.7 Consortium constraints, and one tension to resolve

- Kick-off (24 Apr 2026): *"They suggested going after **targets in the periplasm**. If we can get into the
  cytoplasm, that would be better."* — sensible for permeability, but **ClpP is cytoplasmic** and cannot reach a
  folded periplasmic protein. The v5 proposal points the other way (GyrA/GyrB, "cytosolic"). This should be
  settled as a decision rather than left as two conflicting notes. The reconciliation is the pre-export window
  (§12.3, and the routing table in the spec).
- Kick-off: *"Perhaps **degradability > vulnerability**? Are there **AA motifs** that could signal
  degradation?"* — honest answer: **for partnerless ClpP, no.** There is no sequence degron for activated ClpP;
  the ADEP mechanism bypasses recognition entirely. Motifs exist for ClpXP, i.e. the bar we are not tested on,
  and at proteome scale they are near-empty anyway (§5, and only ~1% of either proteome carries a weighted
  motif).
- Munich Day 2 (16 Apr 2026): Kim Lewis's WP lists ~30 cell-envelope targets — LPS biosynthesis
  (LpxA/B/C/D/H/K, WaaA), LPS transport (MsbA, LptA–G), β-barrel assembly (BamA, BamD), lipoprotein
  trafficking (LspA, Lnt, LolA–E), secretion (SecA/E/Y, YidC, LepB, **FtsH**). Nearly all envelope, so nearly
  all reachable only through the pre-export window, if at all.
- **FtsH is listed as a consortium *target*** while also being one of only three experimentally essential
  proteases in both organisms and, in our earlier analysis, a candidate *handle*. Under the ClpP decision it is
  no longer a handle candidate — which simplifies matters — but one WP inhibiting a protease another WP might
  recruit is worth a sentence at a leadership meeting.
- *A. baumannii* is in the consortium strain panel. Not in our scope; it is also Gram-negative and therefore
  ClpC-less, so the handle analysis transfers unchanged when needed. (Our orthology panel already includes it,
  2,240 Kp orthologs.)
- **`clpP` is dispensable in both organisms** and ADEP resistance via `clpP` inactivation is documented, so the
  handle is itself a resistance liability. Human mitochondrial ClpP is a real off-target; the activator
  literature has mapped the Gram-positive / Gram-negative / human-mito selectivity determinants.

### 12.8 Loose end

Erick Strauss's e-mail of **14 May 2026, "regarding proteolysis rates"**, is unread by us. Proteolysis rates are
precisely the turnover layer this axis needs.

---

## 1. What this axis answers

> **If we could chemically recruit one of this organism's own proteases to a given protein, would that protein
> actually be destroyed?**

This is a *target-side* question, deliberately separate from the medicinal chemistry. It does not ask whether
a degrader molecule can be made or delivered — those are real and currently unsolved problems (§9) — but
whether the protein is physically amenable to being grabbed, unfolded and translocated by a AAA+ protease, or
(in the periplasm) glued into a protease-recognised conformation.

Three properties make this axis structurally different from the other four:

- **It is handle-dependent, and the handle depends on compartment.** There is no single "degradability". A
  cytoplasmic protein is a ClpXP question; an inner-membrane protein is an FtsH question; a periplasmic protein
  is a DegP question. The axis therefore emits a `compartment_handle` alongside the score.
- **Its evidence base is almost entirely *E. coli*.** No *K. pneumoniae* turnover, meltome or limited-proteolysis
  dataset exists. E. coli is therefore the **primary** organism for this axis and Kp the transfer organism —
  an inversion of every other axis in the repo.
- **It has never been validated in the target organism.** No *K. pneumoniae* protein has ever been chemically
  degraded. The nearest result is CLIPPER-mediated depletion of GroEL in *E. coli*, at roughly 40%. This axis
  is a mechanistic prior, and should be presented as one.

---

## 2. Audit: what the repository currently serves

### 2.1 There is no pipeline

`docs/03_degradability.md` specified nine tracks; none is implemented as written. The implementation that once
existed — `scripts/03_annotate_clp_degradability.py` and `src/degradability.py` — was deleted in commit
`8327de3` ("removed all junk scripts"). The `03x` script slot is now occupied by orthology
(`03a`–`03d`), so the axis has no home in the numbered pipeline; stage `10x` is free and is the proposed home.

Consequences visible elsewhere in the repo:

- No `output/results/<org>/<prefix>_degradability.csv` exists, unlike every other axis.
- `docs/mermaid_style.md` still uses the deleted script and its `clp_degradability_score` column as its worked
  style example.
- Degradability is **absent from `CLAUDE.md` entirely**, so an agent reading the project instructions has no
  way to learn that the axis is unimplemented.

### 2.2 What Kp values the webapp actually shows

`scripts/08a_webapp_export.py:degradability_frame()` reads a frozen legacy file,
`data/processed/legacy/klebsiella_pneumoniae_clp_degradability.tsv` (5,728 rows, 22 columns), produced by the
deleted regex script. That script scored:

```python
CM1_RX = re.compile(r"[YAFWLIVM]?[ALV]A[ALV]A$")        # ssrA-family, last 5 residues
CM2_RX = re.compile(r"[RK]{2,}[A-Z]{0,2}[AVILMG]{1,3}$") # MuA-family, last 8 residues
NEND_DESTABILIZING = set("LFYW")                         # residue at position 2
NM1_RX = re.compile(r"^M?[AILVMFW]{2,}")
NM2_RX = re.compile(r"^M?[KR][AILVMFW]")
NM3_RX = re.compile(r"^M?[TS][AILVMFW]")

degron_feature_score = min(1.0, 0.40*cterm_ssra_like + 0.25*nterm_destabilizing
                              + 0.15*cterm_mua_like + 0.05*(nm1+nm2+nm3))
score = degron_feature_score + (0.40 if ecoli_clp_trapped)
                             + (0.20 if halflife=='fast' else 0.10 if halflife=='slow')
tier  = 'high' if score >= 0.50 else 'medium' if score >= 0.25 else 'low'
```

The AlphaFold modulation the spec promised was **never implemented** — the score is pure sequence regex plus
gene-symbol transfer.

### 2.3 The "experimental" inputs are hand-curated substitutes

`data/raw/legacy/clp_substrates/SOURCE.md` is explicit about this, and it is worth restating because the
column names (`ecoli_clp_trapped`, `ecoli_halflife_class`) imply otherwise:

| File | Rows | What it actually is |
| --- | --- | --- |
| `flynn2003_ecoli_clp_substrates.tsv` | **45** | *Not* Flynn's supplementary table. A hand-curated set drawn from the paper's narrative plus the Sauer & Baker 2011 / Baker & Sauer 2012 reviews, written because "the supplementary tables are paywalled and PMC is currently rate-limiting bot access". |
| `nagar2021_ecoli_halflives.tsv` | **35** | *Not* Data Set S2 (1,149 rows). Previously-reported fast degraders named in the paper text plus Flynn overlap, plus a few stable housekeeping negatives. |

Access date recorded: 2026-05-13. Independently confirmed this session: Flynn 2003 remains unobtainable
(cell.com and sciencedirect.com return 403, not in PMC, DSpace blocks scripted fetch), **but Nagar's Data Set
S2 is freely retrievable** via the PMC OA service (§4).

### 2.4 The resulting signal is close to empty

Measured over all 5,728 rows of the legacy file:

| Feature | Proteins firing | Note |
| --- | --- | --- |
| `cterm_ssra_like` | **5** | 0.09% of the proteome |
| `cterm_mua_like` | 183 | |
| `nterm_destabilizing` | **684** (11.9%) | the only broadly-firing feature — and the one the literature contradicts (§5.1) |
| `ecoli_clp_trapped` | **21** | limited by the 45-row curated table |
| `ecoli_halflife_class` known | **15** | limited by the 35-row curated table |
| `clp_degradability_score` > 0 | 2,756 | almost entirely driven by `nterm_destabilizing` and `cterm_mua_like` |
| tier = `high` | **10** | |
| tier = `medium` | 698 | |
| tier = `low` | 5,020 | |

So ~90% of the non-zero signal in the current Kp axis traces to a single feature that has no proteome-scale
association with measured protein half-life. It is worse than that, though: three of the four regexes are
outright defective, and each defect was verified directly against the file on disk.

**Defect 1 — the N-end-rule feature is inverted.** `nterm_destabilizing` flags proteins whose residue 2 is
L, F, Y or W. But the N-end rule requires a destabilizing residue at **position 1 of the mature protein**, and
methionine aminopeptidase excises the initiator Met **only when residue 2 is small** (Ala, Cys, Gly, Pro, Ser,
Thr, Val — Hirel 1989, Frottin 2006). A bulky residue 2 therefore means the initiator Met is **retained**, so
the mature N-terminus is Met, which is *stabilizing*. Checked on all 5,728 rows: of the 684 flagged proteins,
**684 have a bulky residue 2** (L 419, F 157, Y 80, W 28) — i.e. every single one is a protein that does *not*
expose an N-degron. The feature selects exactly the complement of what it intends.

**Defect 2 — that inverted feature carries almost the whole tiering.** Of the 698 proteins in the `medium`
tier, **680 (97.4%)** are there because of `nterm_destabilizing`. So the Kp axis's principal output is not
merely unsupported by the literature; it is anti-correlated with the mechanism by construction.

**Defect 3 — the "ssrA-like" regex cannot match ssrA.** `CM1_RX = [YAFWLIVM]?[ALV]A[ALV]A$` requires an
alternating `X-A-X-A` tail. The *E. coli* ssrA tag ends `…YALAA`, whose final four residues are `ALAA`, which
does not alternate. Confirmed: the regex returns no match against `AANDENYALAA`, the very sequence it is
documented to encode. This is why it fires on only 5 proteins.

**Defect 4 — `NM1` has no discriminative power.** `NM1_RX = ^M?[AILVMFW]{2,}` fires on **1,760 / 5,728 =
30.7%** of the proteome. At that hit rate it cannot separate anything, which is presumably why it was given a
weight of 0.05.

### 2.5 E. coli degradability is a mock

`degradability_frame()` generates E. coli values from an MD5 hash of the accession — deterministic, so
reproducible, but synthetic:

```python
h = _hash01(a, "deg")
score = round(0.0 if h < 0.82 else (h - 0.82) / 0.18 * 0.7, 2)
cterm = _hash01(a, "ct") > 0.94; nterm = _hash01(a, "nt") > 0.88
trapped = _hash01(a, "tp") > 0.985
```

Two consequences:

- The mock uses tier thresholds `high ≥ 0.5 / medium ≥ 0.15`, while the real Kp tiering uses `0.50 / 0.25`.
  The two organisms' tiers are therefore **not comparable**: the mock calls **225 / 4,403 (5.1%)** EC proteins
  `high` against **10 / 5,728 (0.17%)** for real Kp — a 30-fold discrepancy that is purely an artefact.
- `app/app.js:207` hard-returns `0` evidence backing for degradability precisely because of the mock, so the
  axis contributes no confidence anywhere in the UI.

It is correctly flagged: `app/config.js:37` sets `PROVISIONAL = { degradability: ["ec"] }`, the methods panel
says so in prose, and the cells render hatched. The mock is honest — it is just not information.

### 2.6 Other current-state facts

- `comp_degradability` ships with `weight: 0, on: false` in `app/config.js:93`. Only the `◆ Degrader` preset
  (`WEIGHT_PRESETS.degrader`) gives it any weight. It is also **absent from the payload `components` list** in
  both `app/data/kp.json` and `app/data/ec.json`, while `scripts/08b_validate_export.py` nonetheless *requires*
  `comp_degradability` and `degradability_tier` to be present as columns.
- Disorder exists in the repo only as whole-protein pLDDT fractions (`af_frac_low_plddt`,
  `af_frac_very_low_plddt` from `04a`), and is consumed only as a **ligandability penalty** in `06g`
  (`disorder_frac ≥ 0.50 and not has_hard_evidence → tier = "intractable"`). Nothing wires disorder into
  degradability, despite the spec citing "disordered termini as the dominant recognition signal".
- `clp_accessibility` (from `09a`/`09b`) is computed but never used by the degradability axis, and is **blank
  for 3,414 / 5,728 (60%)** of Kp proteins.

---

## 3. The conceptual argument

### 3.1 ClpC-based BacPROTACs do not transfer to *Klebsiella*

Every published small-molecule BacPROTAC (Morreale 2022 and successors) recruits **ClpC1** through a
Cyclomarin-A-derived ligand, and the natural degron the system reads is **phosphoarginine**, written by the
kinase **McsB**. Both are Firmicutes/Actinobacteria proteins.

This was verified by enumerating the **entire ClpA/ClpB (Hsp100) family** per proteome — not by gene-symbol
matching, which is unreliable on a TrEMBL proteome:

| Proteome | Hsp100 members | ClpC? |
| --- | --- | --- |
| *K. pneumoniae* HS11286 | ClpA (759 aa), ClpB (823), ClpK (884) | **no** |
| *E. coli* K-12 | ClpA (758), ClpB (857) | **no** |
| *B. subtilis* (control) | ClpC (810), ClpE (699) | yes |
| *M. tuberculosis* (control) | ClpC1 `P9WPC9` (848) | yes |

Also absent from both Kp and Ec: `clpE`, `clpL`, `mcsA`/`mcsB`, `pafA`/`pup`/`mpa`/`dop` (the Actinobacterial
pupylation route), and any 20S proteasome α/β pair.

Two annotation traps worth documenting so they do not re-enter the dataset later:

1. UniProt entry `A0ABY6X745` is automatically labelled `clpC` in *K. quasivariicola*. Its Pfam/InterPro set
   (PF02861 Clp_N, PF10431, PF00004, PF07724, PF17871) is the generic ClpB/ClpC Hsp100 signature; at 949 aa it
   is longer than any authentic ClpC and it lacks IPR001943, which *B. subtilis* ClpC carries. It is a
   ClpB/ClpK misannotation.
2. *Appl Microbiol* 2026;6:63 screens "*clpC* family genes" in *K. pneumoniae* heat survival. Given (1) and the
   absence of ClpC from HS11286, this is almost certainly detecting *clpK*/*clpB*. It must not be cited as
   evidence of ClpC in Kp.

### 3.2 What Kp does have, and what is essential

Accessions retrieved live from `rest.uniprot.org`; the two entries with missing gene symbols were identified
through the repo's own OrthoFinder output (`data/processed/other/orthology/kp_orthologs_long.tsv`);
essentiality read from `output/results/kpneumoniae/kp_essentiality.csv`.

| Machine | Accession | Locus | Length | `essentiality_score` / tier | Experimental? |
| --- | --- | --- | --- | --- | --- |
| ClpP | `A0A0H3GKH6` | KPHS_11400 | 194 | 0.267 likely_essential | no (ECL8) |
| **ClpX** | `A0A0H3GSR4` | KPHS_11410 | 424 | 0.330 likely_essential | no (ECL8) |
| ClpA | `A0A0H3GQX4` | KPHS_17930 | 759 | — | — |
| ClpS | `A0A0H3GKY1` | KPHS_17920 | 105 | 0.006 non_essential | no |
| SspB | `A0A0H3GTU5` | KPHS_47670 | 164 | — | — |
| Lon | `A0A0H3GJ60` | KPHS_11420 | 741 | 0.344 likely_essential | no (ECL8, serum) |
| HslU | `A0A0H3GK74` | KPHS_00780 | 444 | 0.033 non_essential | no |
| HslV | `A0A0H3GLF6` | KPHS_00790 | 176 | 0.323 likely_essential | no |
| **FtsH** | `A0A0H3GTR0` | KPHS_47270 | 644 | **0.948 essential** | **yes (ECL8)** |
| DegP | `A0A0H3GJM8` | KPHS_09100 | 439 | — | — |
| ClpB | `A0A0H3GTV3` | KPHS_39850 | 823 | 0.291 likely_essential | no |
| ClpK | `A0A0H3GSF4` | KPHS_23030 | 884 | — | — |
| SmpB | `A0A0H3H1S8` | KPHS_40610 | 160 | 0.237 non_essential | no |

Notes: Kp ClpP is annotated at 194 aa against *E. coli*'s 207 — the HS11286 entry may be N-terminally
truncated, which matters if the model is ever used structurally. `clpP-clpX-lon` form a syntenic cluster as in
E. coli, and `clpS-clpA` are adjacent. The `clpA` and `sspB` symbols are **absent from HS11286's gene-name
column** and would be missed by a symbol grep — hence track 3.0's insistence on sequence-based census.

Additional periplasmic proteases present: DegQ `A0A0H3GZ27`, DegS `A0A0H3GY66`, Prc/Tsp `A0A0H3GQ33`. There is
also an uncharacterised "proteasome-type protease" `A0A0H3GUL5` (KPHS_33090, 248 aa, PF00227, Anbu/BPH-like,
distinct from HslV) that nobody has functionally described.

**Resistance framing.** ClpP, ClpX, Lon, HslUV and ClpS are all dispensable in Kp — confirmed in vitro by Lin
et al. 2026, where Δ*clpX* and Δ*clpP* grow like wild type. This cuts both ways. It is *desirable* that the
recruited machine is not itself a lethal target, because then the degrader's lethality comes from the
neo-substrate rather than from poisoning proteostasis. But it also means **loss-of-function of the recruited
protease is an available resistance route**. FtsH is the one exception (essential, experimentally, in ECL8) but
is membrane-anchored and therefore restricted in substrate scope. This should be stated in any target
narrative, not omitted.

### 3.3 The disorder hypothesis, sharpened

The proposition "a disordered region indicates degradability" is correct in direction and wrong in detail. The
useful feature is not global disorder fraction.

**(a) It must be terminal.** Won et al. 2024 screened **72 native *M. smegmatis* proteins** for
rapamycin-induced (FRB/FKBP) degradation by ClpC1P1P2 and fitted a Lasso model over 485 descriptors on 54 of
them, validating on 18 held out (r > 0.6, P = 0.0057). The top-ranked feature was the **mean disorder
propensity of the N-terminal 30 residues** (flDPnn, threshold 0.2), and the paper states it plainly: *"the more
disordered the N-terminal sequence is, the more efficiently TPD could act."* The C-terminal 30 residues showed
**no** correlation. Windows tested were N-/C-terminal 15 and 30 aa plus full length — the scheme track 3.1
adopts.

**(b) It is confirmed proteome-wide.** Gupta et al. 2024 measured turnover for ~3,200 *E. coli* proteins:
*"disordered proteins had significantly shorter half-lives than ordered proteins"*, **p ≈ 1e-208**. Also
short-lived: proteins under 10 kDa, transcriptional regulators, and Fe–S-cluster proteins. And a compartment
result that shapes the architecture: **cytoplasmic proteins are selectively degraded while membrane and
periplasmic proteins are largely stable** — so the periplasmic DegP route is not about natural turnover, it is
about glue-induced recognition of a non-native conformation.

**(c) There are hard, per-terminus length thresholds.** The unstructured segment is an *initiation region*
that the unfoldase must thread before it can pull:

| Threshold | Value | Source |
| --- | --- | --- |
| ClpXP recognition, closed axial channel (2 pore-1 loops) | **~5 residues** | Saunders 2020 *PNAS* |
| ClpXP engagement, open axial channel (5 pore-1 loops) | **~20 residues** | Saunders 2020; Fei 2020 *eLife* |
| ClpX pore → ClpP active-site reach | **~37 residues** | Kenniston 2005 *PNAS* (tail length of partial products) |
| Longer tails → higher commitment vs futile release | continuous | Kenniston 2005 |
| FtsH, cytosolic N-terminal tail | **>~20 residues** | Chiba 2002 *J Bacteriol* |
| FtsH, cytosolic C-terminal tail | **~10 residues** | Chiba 2002 |

Degradation rate also depends on **which** terminus is pulled, because unfolding force requirements differ by
geometry (*PNAS* 2017, directional pulling). So the feature must be computed per terminus.

**(d) Internal disordered loops count.** Internally placed ssrA tags are degraded by ClpAP and ClpXP
(Hoskins 2002), and *E. coli* FtsZ carries two independent, additive ClpXP sites — a C-terminal one (residues
366–383, critical Arg379/Lys380) and an **internal one in the disordered linker at residues 349–358,
`QEQKPVAKVV`** (Camberg 2014). MacKrell 2026 adds that mutating PdeH's N-terminal extension abolishes ClpXP
recognition.

Independent corroboration from a different modality: ADEP-activated ClpP degrades nascent chains with
susceptibility **inversely proportional to folding speed**, and destroys FtsZ specifically via its
conformationally flexible N-terminal domain (Silber 2020). Same physics, different trigger.

**Two caveats that must travel with this claim.** First, intrinsic disorder is not a universal accelerant of
proteolysis — disordered regions can be context-dependently protease-resistant; the argument here is specific
to ATP-dependent *processive* AAA+ proteases, where an accessible unstructured initiation region is a
mechanical requirement. Second, the Won 2024 N-vs-C asymmetry is a ClpC1 result in an Actinobacterium; ClpX
reads both termini and its best-characterised degron is C-terminal, so N- and C-terminal disorder should be
carried as separate features rather than importing the asymmetry.

### 3.4 Motifs alone do not work

This is the most important negative result for the axis design. Nagar et al. 2021 fitted half-lives for 1,149
*E. coli* proteins over 188 features and found the most informative were **protein-interaction-network
properties** (node connectivity, node2vec embeddings), followed by disorder propensity, mass and pI — with
**no significant association for the N-end rule or the ClpXP recognition motifs** applied to annotated ORF
termini. Gupta 2024 independently reports that *"rapidly degrading proteins showed no enrichment for
previously reported destabilizing N-terminal amino residues."*

Two reconciliations, both of which change the design:

1. **Real degrons are usually latent.** Humbard et al. 2013 pulled >100 ClpS-binding *E. coli* proteins off a
   native lysate and found that most had been **N-terminally truncated by endoproteases first** — the degron
   is generated post-translationally, so scanning the annotated Met1 terminus systematically under-counts
   substrates. LexA is the conformational analogue: its autocleavage fragment's *new* C-terminus
   (Val82-Ala83-Ala84-COOH) is a bona fide C-motif-1 degron that does not exist in the intact protein
   (Neher 2003).
2. **Accessibility, adaptor availability and complex burial dominate**, which is exactly what the
   PPI-connectivity result encodes.

DEtox reaches the same conclusion from the opposite direction: *"efficient proteolysis of proteins lacking
ssrA-like tags likely requires adaptors or multivalent interactions."*

**Design consequence:** every motif channel is scored as **motif × terminal accessibility × adaptor
presence**, never as motif presence alone; and PPI connectivity is carried as its own covariate (§10).

### 3.5 Three handles, three compartments

| Compartment | Kp proteins | Handle | Evidence |
| --- | --- | --- | --- |
| Cytoplasm | 513 called + most of the 3,414 unknown | **ClpXP** via the ClpX zinc-binding domain | **CLIPPERs** (Izert-Nowakowska 2025): anchor = the SspB "XB" motif (K_D ~65–70 nM for the ClpX ZBD, cf. PDB 2DS8); degraded endogenous, untagged, degronless GroEL in live *E. coli*; abolished in Δ*clpX*/Δ*clpP*. ssrA was deliberately rejected as an anchor because it makes the degrader itself a substrate. |
| Inner membrane | **1,320** (1,206 `inner_membrane` + 114 generic `membrane`) | **FtsH** | Classical cytosolic-tail thresholds (Chiba 2002) **plus** the 2026 result that **lipid-facing polar residues inside TM helices** target membrane proteins to FtsH — even when folded, and without the long tail (Chai-Danino 2026). |
| Periplasm / OM / secreted | **481** (113 + 84 + 268 + 16) | **DegP** molecular glues | Taylor et al. 2026: tazobactam itself accelerates DegP-mediated degradation of TEM β-lactamases, enhanced by *linkerless* incorporation of DegP-substrate dipeptide motifs; improved piperacillin synergy plus oral-dosing PK. A glue, not a chimera — small, no linker, no permeability penalty, acting in the target's own compartment. |

CLIPPER magnitudes, stated honestly: ~17% GroEL depletion at 1 h and ~40% at 6 h (42 °C), ~40% in 10 h
in vitro, 20–25% growth inhibition at 30 °C. It is a plasmid-encoded peptide tool, not a drug; the authors
name cytoplasmic peptide delivery as "the ultimate limitation."

**The cautionary counter-example.** Nie et al. 2025 conjugated nacubactam to the full ssrA degron
(`AANDENYALAA`) and reported degradation of CTX-M-14 with cefotaxime resensitisation. The molecule is
**1,669 Da**, and CTX-M β-lactamases are **periplasmic while ClpXP is cytoplasmic** — the target had been
engineered to remain in the cytoplasm. This is why `compartment_handle` is a hard routing decision in this
axis rather than a soft weight: a compartment mismatch is not a penalty, it is a category error.

### 3.6 Engageability versus depletability

Two distinct properties are hiding inside one score, and they pull in **opposite** directions on natural
half-life:

- **Engageability** — can the protease grab, unfold and translocate it? Favours an accessible unstructured
  terminal initiation region of sufficient length, low thermodynamic/kinetic stability, monomeric state, small
  size, and the right compartment. Correlates with a **short** natural half-life.
- **Depletability** — once degraded, does it stay gone long enough to matter? Favours **low synthesis flux**:
  low copy number, slow resynthesis. Correlates with a **long** natural half-life.

So "fast turnover ⇒ good target" is wrong as a single term. Engageability is built from structure, termini,
stability and protease attribution, and enters `degradability_score`. **`resynthesis_burden` is emitted as a
separate column and used only in shortlist filtering.**

The evidence for the depletability side is quantitative: Won 2024 found higher steady-state abundance
anti-correlates with induced degradation rate at **r = −0.69, P = 6.6 × 10⁻⁹**, and the CLIPPERs authors state
independently that *"targets with a high protein copy number may also be more resistant to the TPD approach."*

**GroEL is the worked negative control.** It fails on both counts: both its N- and C-termini face the interior
of the barrel (so only dissociated or newly-translated subunits can be engaged), and it is among the most
abundant proteins in the cell. It is precisely the protein CLIPPERs could only deplete ~40%. Any implementation
of this axis that scores GroEL highly is wrong, and that is a testable assertion (§8).

A further practical modifier: chemically induced degradation is **partial**. Prefer targets with a steep
dose–response — high CRISPRi transcriptional vulnerability, which `07b`/`07l` already provide — over targets
that are merely binary-essential.

### 3.7 Why degrade rather than inhibit

Targeted degradation is event-driven rather than occupancy-driven, removes non-catalytic and scaffolding
function, and does not require an active-site pocket. Its *unique* value is therefore concentrated on proteins
that are **essential but poorly ligandable** — a quadrant the repo already scores on two other axes. Track 3.9
emits this as `tpd_advantage = essentiality × (1 − ligandability) × degradability` so that the quadrant is a
sortable column rather than something a user must infer by combining two sliders.

---

## 4. Dataset inventory and access routes

**Access route note.** Direct PMC URLs of the form
`https://pmc.ncbi.nlm.nih.gov/articles/instance/<id>/bin/<file>.xlsx` return an HTML interstitial to scripted
clients. The reliable programmatic route for any open-access paper is the **PMC OA service**:

```
https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC6056769
  → https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/26/2d/PMC6056769.tar.gz   (all supplementary files)
```

This is already the third rung of `scripts/07a_fetch_essentiality.py`'s fetch ladder, so `10a` can reuse it
verbatim. Every route in the tables below was verified on 2026-08-04; ✅ means an HTTP 200/206 response was
received, and for PMC entries that the OA service returned a tarball path.

### 4.1 Rule sets (tracks 3.2, 3.7)

| Source | What it supplies | Route | Status |
| --- | --- | --- | --- |
| **DEtox** — Beardslee & Schmitz 2024, PMC10862746 | C-terminal ClpXP degron PSSM: `[L/F/Y/W]-[L/A/R]-L-A-A` over positions −5→−1; ~1% of random 5-mers are degrons; Ala-Ala at −2/−1 in 91% of top 100; depletions polar@−5, bulky/β-branched@−4, P/H/G@−3; best tag `FKLVA` 639×. Run in WT/Δ*clpX*/Δ*clpP* → attribution is essentially all ClpXP. | PMC OA tarball | ✅ |
| **Cragan 2025** — *JBC* 301(4):108365, PMC11986505 | Lon C-degron consensus `x–[L/I]–[L/I/V]–H-COOH`; His@−1 preferred, D/R/K@−1 disfavoured (H→D → ≤0.3 min⁻¹·hexamer⁻¹); K_D 0.24 → >50 µM; conserved E. coli / Yersinia / Mycoplasma. Also cleavage-site P1 preference Phe ≫ Leu ≈ Ala. | PMC OA tarball | ✅ |
| **Sen 2025 "N-FIVE"** — bioRxiv, PMC12258705 | N-degron P1–P5 stability model from ~2.2 M variants; effect sizes in §5.1; ΔclpS / Δaat dependency controls; model + code at `github.com/KunjapurLab/N-terminal-cluster-stability`. **Preprint — label as such.** | PMC OA tarball + GitHub | ✅ |
| **Chai-Danino 2026** — *Nat Commun* 17:3067 | FtsH rule: lipid-facing polar residues within TM helices, sufficient even for folded proteins without a long cytosolic tail. | Nature (open) / bioRxiv 2023.12.12.571171 | — |
| **Taylor 2026** — bioRxiv `10.64898/2026.03.30.715243` | DegP-substrate dipeptide motifs; the periplasmic glue concept. **Preprint; bioRxiv blocks scripted fetch** — read manually. | manual | ⚠ |
| **Won 2024** — *Nat Commun* 15:4065, PMC11094019 | Feature recipe + labels. Supp Data 1 = 348-protein feature matrix; Supp Data 2 = 72 targets with degradation constants `dG/dt`; Supp Data 3 = predicted degradability for 348 conserved essential proteins. Code/data Zenodo `10.5281/zenodo.11004395`. | `static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-48506-8/MediaObjects/41467_2024_48506_MOESM4..10_ESM.xlsx` | ✅ |
| **Flynn 2003** — *Mol Cell* 11:671, PMID 12667450 | The five ClpX motif classes and the trapped-substrate census. **Not reachable headlessly**: cell.com and sciencedirect.com both return 403 and there is no PMC record. The content is in the **article body**, not a separate supplement, so the route is an authenticated browser session against the article HTML — the same chrome-devtools `evaluate_script` same-origin fetch already used for the gated essentiality tables — not the OA services. | authenticated session | ⚠ |
| Cross-bacterial Clp-trap pool — Lunge 2020 (PMC7363115), Bhat 2013 (PMC3681837), Graham 2013 (PMC3807464) | The old §3.3b trap pool. **All three return `idIsNotOpenAccess`** from the NCBI OA service; Feng 2013 (ACS) and Ziemski 2021 (Wiley) have no PMC record at all. So the entire cross-bacterial trap channel needs an authenticated session, which is part of why it is demoted below the open E. coli turnover data. | authenticated session | ⚠ |
| Met-excision rules — Hirel 1989 (PMC298257), Frottin 2006 (PMID 16963780) | MAP specificity gate for the N-degron channel: Met removed when residue 2 ∈ {A,C,G,P,S,T,V}, with V and T much less efficiently cleaved and P2′–P4′ context able to slow the reaction. | open | ✅ |

### 4.2 E. coli proteome-scale datasets (tracks 3.3–3.5)

| Source | Column(s) | N | Strain | Identifier | Route | Status |
| --- | --- | --- | --- | --- | --- | --- |
| **Gupta 2024** *Nat Commun* 15:5890, PMC11246515, PXD042444 | `k_deg` (h⁻¹) + half-life × 13 conditions; protease attribution (ClpP 64 / Lon 14 / HslV 1 / additive 82 / redundant 41 / triple-KO-persistent ~100); Supp Data 4 = ~600 measured in-vivo N-termini | ~3,200 | NCM3722 | gene name + UniProt | `…/41467_2024_49920_MOESM4..11_ESM.xlsx` | ✅ |
| **Nagar 2021** *mSystems*, PMC7857536, PXD022112 | half-life + "actively degraded" label + **188-feature matrix** incl. 128 node2vec PPI features (Data Set S2, 1,149 rows) | 1,149 | K-12 JW2806-1 | UniProt + gene + ORF | PMC OA tarball, **or** the Europe PMC `supplementaryFiles` endpoint, which returns a zip directly (`sd001` raw time course, `sd002` half-lives, `sd003` stability classes) | ✅ **the single biggest unlock — CC-BY and fully automatable** |
| **MacKrell 2026** *PNAS* 123, PMC12974527, PXD062881 | half-life, **exponential and stationary** phase; 88 / 56 pronouncedly unstable; Supp S4 = ML predictions for unmeasured proteins | 1,810 / 1,339 | MG1655 | gene + UniProt + EcoCyc | PMC OA tarball | ✅ |
| **Cappelletti 2021** *Cell*, PMC7836100, PXD022297 | **LiP-MS accessibility density**, half-tryptic fraction, changing-peptide count across 8 carbon sources | ~1,900 | BW25113 | gene name | PMC OA tarball | ✅ |
| **Mateus 2018** *Mol Syst Biol*, PMC6056769, PXD009495 | **Tm (°C)** (Dataset EV1) — residualise on localization before use | 1,738 | BW25113 | UniProt + gene | PMC OA tarball | ✅ |
| **Mateus 2020** *Nature*, PMC7612278, PXD016589 | thermal-stability score across 121 Keio backgrounds → **variance = conformational plasticity** (Supp Data 5) | 1,764 | Keio panel | UniProt | `…/41586_2020_3002_MOESM8_ESM.xlsx` | ✅ |
| **Jarzab 2020** Meltome Atlas, PXD011929 | cross-species Tm prior for sanity-checking ortholog transfer | ~48k total | 6 prokaryotes | mixed | `…/41592_2020_801_MOESM4_ESM.xlsx` | ⚠ file identity inferred from size — verify contents |
| **Schmidt 2016** *Nat Biotechnol*, PMC4888949 | **copies/cell × 22 conditions** | 2,359 | BW25113 | gene + UniProt + b-number | PMC OA tarball | ✅ |
| **eSOL / Niwa 2009** *PNAS* | **% solubility** (chaperone-free PURE system), bimodal | **3,198** | K-12 | **b-number / JW ID** | `dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip` | ✅ (original tanpaku.org site is dead) |
| **Niwa 2012** *PNAS* | ΔSolubility with GroEL / TF / DnaK = chaperone dependence | ~800–3,000 | K-12 | b-number | open | ✅ |
| **Györkei 2022** *Sci Rep*, PMC9023497 | in-vivo solubility limit + 3-class aggregation-rate label | 2,577 | K-12 AG1 | UniProt + name | Springer supp + `group.szbk.u-szeged.hu/...zip` | ✅ |
| **To 2021** *JACS* + **To 2022** *PNAS* | (non-)refoldability, protein- and domain-level; 396/1,198 (33%) non-refoldable | 1,198 | K-12 | gene | ACS gated → NSF-PAR `par.nsf.gov/servlets/purl/10311595` or bioRxiv 2020.08.28.273110 | ⚠ |
| **Calloni 2012** *Cell Rep* | DnaK interactome (~700); GroEL client class I/II/III = 30/80/42 | ~700 | K-12 | gene | open | ✅ |
| **Niwa 2022** *Molecules*, PMC9228906 | Lon vs ClpXP vs HslUV attribution over ~80 obligate GroE substrates | ~80 | K-12 | gene | PMC OA | ✅ |
| **Neher 2006** *Mol Cell* | ClpXP trap ± DNA damage; 25% of the SOS proteome; 9 confirmed substrates | ~50 | K-12 | gene | open | ✅ |
| **Westphal 2012** + **Arends 2016** | FtsH trap sets (15 + >50 candidates) | >50 | K-12 | gene | PMC3522291 / Wiley | ✅ / ⚠ |
| **Humbard 2013** *JBC*, PMC3789986 | N-degradome: >100 ClpS binders; Table 1 measured neo-N-termini; ~20–30% Aat-dependent | >100 | K-12 | gene | **not in the PMC OA subset** → `jbc.org/article/S0021-9258(20)48862-3/pdf` | ⚠ |
| **PHDB** — Ramakrishnan 2020 | 13 chaperone studies + abundance + translation rate + solubility + TANGO + IUPred + contact order, in one place | **4,305** | K-12 | UniProt | `phdb.switchlab.org` — **web UI only per the paper; bulk export unverified** | ⚠ |
| **STEPdb 2.0** | localization/topology/abundance/solubility/disorder spreadsheets | proteome | K-12 | UniProt | `stepdb.eu` | ⚠ several columns are predicted, not experimental |

### 4.3 *K. pneumoniae*-native anchors

| Source | Column(s) | N | Strain | Identifier | Status |
| --- | --- | --- | --- | --- | --- |
| **Illenseher 2025** *Front Microbiol* 16:1528869, PMC12127431, PXD052921 | **iBAQ abundance** ± heat / oxidative stress (cell pellet >2,800; exoproteome 2,300) | >2,800 | ATCC BAA-2146 Δ*wza* | NCBI locus tags `Kpn2146_####` (CP006659.2) → needs DIAMOND mapping via `map_strain_by_sequence` | ✅ |
| **Muselius 2020** *Front Microbiol* 11:546, PMC7194016, PXD015623 | **the only Kp-native proteolysis data**: Δ*lon* vs WT; 59 significant of 2,074, of which **26 accumulate = candidate Lon substrates** (BamA 8.91 log₂, LuxR 7.02, CmlA2 6.65) | 2,074 | K52 serotype | gene + UniProt | ✅ |

---

## 5. Demoted and dropped, with the measurements

### 5.1 The naive N-end rule (`nterm_destabilizing`, weight 0.25)

Fires on **684 / 5,728 (11.9%)** Kp proteins — the only broadly-firing feature in the legacy score, and
therefore the source of most of its non-zero signal (**680 of the 698 `medium`-tier calls**). Two independent
proteome-scale datasets find no association with measured half-life (Gupta 2024, quoted verbatim in §3.4;
Nagar 2021). It is also ClpS-dependent, and ClpS is `essentiality_score` 0.006 in Kp — present but entirely
dispensable.

And, as §2.4 Defect 1 establishes, the implementation was **inverted**: all 684 flagged proteins have a bulky
residue 2, which means their initiator Met is retained and their mature N-terminus is Met. So the axis was not
just weighting a weak feature heavily — it was ranking the wrong proteins. Any replacement must gate on
Met-excision first and only then evaluate position 1.

Replaced by the P1–P5 model, which explains the failure: position 1 alone is not the determinant. From Sen
2025 (~2.2 M variants, protein-stability-index readout):

| Determinant | Effect |
| --- | --- |
| P1 destabilizing potency order | `F ≈ R < L < W ≈ K ≈ Y` |
| **Pro at P2** | **+1.05 PSI units** (P < 1e-99, effect size 0.71) — strongest single stabilizer, i.e. rescues a destabilizing P1 |
| **Gly at P2** | **+0.63 PSI** (P < 1e-99, ES 0.47) |
| **Gln at P2** | **−0.51 PSI** — as strong a degron determinant as Leu at P1 |
| Clustered bulky FLWY over P2–P5 | mean **−0.56 PSI** |
| Net negative charge over P2–P5 | stabilizes; at charge −4 it rescues even FLYWRK P1 substrates (ClpS's pocket is flanked by Asp36 and Glu41) |
| Multiple Gly/Ser in P2–P5 | cumulatively stabilizing |
| Non-canonical P1 | Cys, Gln, His can form degrons in permissive context |
| Dependency | all P1–P5 effects vanish in Δ*clpS*; the R/K effect vanishes in Δ*aat* |

So a protein with Leu at P1 and Pro at P2 is *stable*, and the legacy feature would have called it degradable.

### 5.2 The ssrA C-terminal regex (`cterm_ssra_like`, weight 0.40)

Fires on **5 / 5,728** Kp proteins — and per §2.4 Defect 3 the regex **cannot match the ssrA tag it is
documented to encode**, so that hit rate is not even a measurement of what it claims to measure. Separately,
the hit rate *should* be low on the biology: ssrA/C-motif-1 is the *engineered* degron, appended
post-translationally by tmRNA, not something genomes encode. Lyu et al. 2026 make the point
mechanistically: ClpXP's closed-channel conformation preferentially handles C-motif-1 (ssrA) while the open
channel favours N-motifs 1–3 and C-motif-2, i.e. **the N-motifs and C-motif-2 are the natural substrate route
and ssrA is the special case**. Replaced by the DEtox PSSM, which was derived from ~100,000 measured tags
rather than written by hand, and which quantifies *how much* each position contributes.

### 5.3 The hand-curated Flynn and Nagar stand-ins

45 and 35 rows respectively (§2.3), yielding trap flags for 21 Kp proteins and half-life classes for 15.
Replaced by Gupta 2024 (~3,200 proteins with `k_deg` *and* protease attribution), Nagar Data Set S2 (1,149
rows, freely retrievable), MacKrell 2026 (adds stationary phase), and open trap sets: Neher 2006 (ClpXP/SOS),
Niwa 2022 (Lon vs ClpXP vs HslUV), Westphal 2012 + Arends 2016 (FtsH), Humbard 2013 (N-degradome).

### 5.4 ClpK as a degron channel (old §3.4)

The previous spec treated ClpK as a *Klebsiella*-specific channel needing "a dedicated rule set / HMM".
ClpK **is** present in HS11286 (`A0A0H3GSF4`, KPHS_23030, 884 aa) — an earlier gene-symbol grep missed it — and
it does confer thermotolerance and nosocomial persistence (Bojer 2010). But it is a ClpB-class **disaggregase
with no ClpP partner** and no characterised substrate motif, so it cannot deliver anything to a peptidase.
Kept as a documentation note.

### 5.5 The ESM-2 classifier on pooled trap labels (old §3.5)

Would have been trained on ~45–150 positive labels. Gupta / Nagar / MacKrell provide ~3,200 *continuous*
labels, which turns this into a properly-powered regression, and Won 2024 supplies a validated feature recipe
to fit. Recommended reframing: **re-fit the Won 2024 feature set on E. coli `k_deg`**, then use the 72-protein
mycobacterial screen as an orthogonal held-out test.

### 5.6 HslUV as its own channel

Dropped. Substrate list is tiny (SulA, RcsA, RpoH, TraJ, RNase R, plus YbaB from a proteome-microarray screen
that yielded essentially one validated hit). There is no sequence consensus; the only computable determinant is
"an aromatic π-system or a cation at the degron position in a locally unstructured segment", which is
degenerate with the Lon aromatic-cluster feature. Petkov 2023 judges it mechanistically under-characterised.

### 5.7 OrthoDB cross-bacterial expansion (old §3.3a)

Redundant. `scripts/03a_orthology_general.py` and `03c_orthology_focused.py` already produce the ortholog
panel, and `src/essentiality.py:transfer_ecoli_to_kp` already implements the transfer.

---

## 6. Composite design

```
route  compartment_handle    cytoplasm → clpxp · inner_membrane → ftsh · periplasm/OM/secreted → degp
                             · unknown → NA (score NaN, not 0)
0.35   evidence_engagement   terminal + internal initiation regions, per terminus (3.1)
0.25   evidence_degron       handle-matched channel, each × terminal accessibility × adaptor presence (3.2, 3.7)
0.25   evidence_turnover     measured k_deg / instability / protease attribution / trap-set membership (3.3)
0.15   evidence_biophysics   LiP accessibility, Tm, refoldability, solubility, aggregation, chaperone
                             dependence, knot penalty, assembly state (3.4, 3.6)

degradability_score = Σ(w·e) / Σ(w over channels with a measurement)
degradability_tier  ∈ {high, medium, low}
resynthesis_burden    separate column — shortlist filter only, deliberately NOT in the score (§3.6)
tpd_advantage       = essentiality × (1 − ligandability) × degradability
degradability_sources ";"-joined provenance tokens, in 07h's style
```

Two design commitments worth stating explicitly:

1. **Missing tracks are renormalised, not zero-filled.** A protein with no LiP measurement should not be
   penalised as though it had been measured and found inaccessible. Follow
   `scripts/07h_essentiality_merge.py`'s idiom — numerator/denominator accumulators, sub-score builders
   returning `np.nan` rather than `0.0` when there is no measurement — and *not* `06g`'s zero-fill. The
   load-bearing contract is on the sub-score builders.
2. **`compartment_handle` is a router, not a multiplier.** Periplasmic proteins are scored on the DegP
   channel; they are not zeroed. Zero (or rather NaN) is reserved for `compartment_handle == none`. This is a
   change from the natural first instinct, and it is what the Taylor 2026 result buys.

Tier values stay `high` / `medium` / `low` because `app/config.js:355` already maps those abbreviations and
`app/config.js:598` already defines `stable` / `moderate` / `labile` display bands over `comp_degradability`.
No front-end churn is needed there.

---

## 7. Planned pipeline

Stage `10x` (`03x` is occupied by orthology). New `src/degradability.py` mirroring `src/essentiality.py`:
import `src/ligandability.py` and re-export the generic machinery (`ORGANISMS`, `run_diamond_blastp`,
`load_accessions`, `load_genes`, `results_dir`, `processed_dir`, `af_cif_path`, `transfer_ecoli_to_kp`,
`gene_aliases_to_uniprot`, `map_strain_by_sequence`) so each `10*` script needs one import. Standard CLI
pattern: `ArgumentParser(description=__doc__)`, `--organism` with `choices=list(D.ORGANISMS)` defaulting to
`kpneumoniae`, `_, prefix = D.ORGANISMS[org]`, `print(..., flush=True)`, and all-NA same-schema tables for
non-applicable tracks.

| Script | Purpose | Expected coverage |
| --- | --- | --- |
| `09c_deeplocpro.py` | **Prerequisite.** Fill the 3,414/5,728 blank localization calls with DeepLocPro (prokaryote-specific, ESM-based, CPU-friendly, benchmarks above the PSORTb 3.0 ensemble). PSORTb is parked because it crashes under Rosetta emulation. Without this the compartment router is unusable for 60% of the proteome. | → ~100% |
| `10a_fetch_degradability.py` | Reuse `07a`'s `Dataset` dataclass, the CDN → Europe-PMC → **NCBI-OA-tarball** ladder, magic-byte sniffing and non-blocking `PLACEHOLDER.txt` failure. Manifest = §4. | — |
| `10b_terminal_disorder.py` | **The core, and computable today with zero downloads.** Per-residue pLDDT from the AlphaFold CIFs already on disk (**5,727/5,728** kp, **4,371/4,403** ec) → `nterm{15,30,50}_plddt_mean`, `cterm{15,30,50}_plddt_mean`, leading/trailing runs below pLDDT 70 and 50, `internal_disorder_max_run`, and `initiation_region_len_{n,c}` bucketed against the 5 / 20 / 37 aa thresholds. Terminal burial via `Bio.PDB.SASA.ShrakeRupley` (biopython 1.87 present in `gradi`). `metapredict` for the sequence-side second opinion and the one protein with no model. | **~100%** |
| `10c_degron_motifs.py` | DEtox ClpXP PSSM · Lon C-degron · N-degron P1–P5 (with the MAP Met-excision gate and Aat handling) · Flynn N-motifs at low weight. Each gated on 3.1 accessibility and on adaptor presence. Plus the internal-degron scan. | 100% |
| `10d_ecoli_turnover.py` | Gupta / Nagar / MacKrell + trap sets. Direct for ec via `gene_aliases_to_uniprot`; kp via `transfer_ecoli_to_kp`. | ~3.2k ec, ≤55% kp |
| `10e_biophysics.py` | LiP, Tm (localization-residualised), stability-score variance, eSOL, aggregation class, refoldability, chaperone dependence, knot topology. | varies |
| `10f_abundance_burden.py` | Schmidt copies/cell; Kp iBAQ via DIAMOND. → `abundance_percentile`, `resynthesis_burden`. | ~50% |
| `10g_assembly_state.py` | UniProt subunit text via the ec ortholog + PDB assembly stoichiometry (extends `04c`; `kp_pdb_coverage.csv` currently carries only `pdb_n_chains`). | partial |
| `10h_compartment_channels.py` | ClpXP / FtsH (tail thresholds + TM lipid-facing polar residues) / DegP (substrate motifs, signal peptide, foldedness) → `compartment_handle`. | ~100% |
| `10i_degradability_merge.py` | The composite of §6 → `<prefix>_degradability.csv` + `_shortlist.csv`. | 100% |
| `10j/10k/10l_*_plots.py` | Stylia NPG slides, 2×3 per organism: degron/initiation channels · evidence summary · cross-axis landscape. | — |

Environment additions needed: `metapredict` and `peptides` in `requirements.txt`. Everything else is already
present in the `gradi` env (verified: biopython 1.87, pandas 2.3.3, numpy 2.4.6, scikit-learn 1.9.0,
matplotlib 3.10.9, stylia, torch 2.12.0, esm 3.2.1, openpyxl 3.1.5).

### Webapp changes the build will require

- `scripts/08a_webapp_export.py` — replace `degradability_frame()` (L251–318) with
  `_read_subset(rdir / f"{prefix}_degradability.csv", DEG_COLS)`; delete the `_hash01` mock; extend `DEG_COLS`
  with `compartment_handle`, `resynthesis_burden`, `tpd_advantage`; add `comp_degradability` to the payload
  `components` list, which currently omits it.
- `app/config.js` — remove `PROVISIONAL = { degradability: ["ec"] }` (L37) and `provisionalOrgs: ["ec"]`
  (L426); give `comp_degradability` a non-zero default weight and `on: true` (L93); add column definitions for
  the new headline columns. A `compartment_handle` legend-as-filter on the Map would follow the pattern
  established in commit `0618882`.
- `app/app.js` — `axisBacking()` (L207) hard-returns 0 for degradability because of the mock; make it
  evidence-based like the other axes. Bump `?v=N` on any edited asset.

---

## 8. Validation checks the build must pass

Recorded here so they are not lost between this documentation pass and the implementation.

1. **Coverage.** `10b` produces non-null rows for 5,727/5,728 kp and 4,371/4,403 ec proteins. Hand-check one
   IDP-tailed protein and one compact enzyme.
2. **Positive controls.** Gupta's fast set (RpoS, LpxC, DnaQ; t½ 0.7–1.2 h) must land in the top decile of
   `evidence_engagement` for E. coli. FtsZ must fire `internal_degron_hit` at the linker (residues 349–358).
   SulA must fire the Lon C-degron (its sul20C tail `…LSGLKIHSNLYH` matches the consensus).
3. **The negative control.** GroEL (`groL`) must score low on assembly state and high on `resynthesis_burden`
   (§3.6). If GroEL ranks well, the implementation has not encoded the CLIPPERs lesson.
4. **Compartment routing.** TEM / CTX-M β-lactamase orthologs must route to `degp`, never `clpxp` — the
   Nie 2025 error must be impossible by construction.
5. **The axis's own internal validation set.** `evidence_engagement` and `evidence_degron` must correlate
   significantly and positively with measured `k_deg` across the ~3,200 E. coli proteins. If they do not, the
   construction is wrong. This is the single most valuable reason to build E. coli first, and it is a check the
   legacy pipeline never had.
6. **Kp-native anchors.** The 26 Δ*lon* accumulators (PXD015623) should be enriched among high
   `cdegron_lon_score` Kp proteins. Kp iBAQ should correlate with Schmidt copies/cell across ortholog pairs.
7. **Cross-organism comparability.** `10i` tier distributions must be comparable between kp and ec. The mock's
   30-fold `high`-rate discrepancy (§2.5) is precisely the defect being fixed.
8. **Export.** `08a` then `08b` exits 0; `08b` already requires `comp_degradability` ∈ [0,1] and
   `degradability_tier`. Then open `app/` and exercise the `◆ Degrader` preset, confirming no hatched
   provisional cells remain.

---

## 9. Gaps that cannot be closed

- **No *K. pneumoniae* protein has ever been chemically degraded.** The axis is a prior. The nearest
  experimental result in any Gram-negative is ~40% depletion of GroEL in *E. coli* by a plasmid-expressed
  peptide.
- **No Kp turnover, meltome or LiP-MS dataset exists**, and none for *Salmonella* or *Pseudomonas* either
  (the one *Salmonella* dynamic-SILAC set covers 870 proteins and is paywalled). Every biophysical column is
  E. coli-derived and transferred, capped at **3,179/5,728 (55.5%)** ortholog coverage. The two Kp-native
  anchors in §4.3 exist only to check that the transfer is not systematically off.
- **No Kp degradome and no Kp ClpP structure.** UniProt `A0A0H3GKH6` has an AlphaFoldDB entry and no PDB
  cross-reference; an RCSB full-text search for "Klebsiella pneumoniae ClpP" returns only unrelated entries.
  Structural work would have to be modelled on E. coli ClpXP.
- **Rule transfer is defensible for Lon, less so for ClpX.** The Lon C-degron consensus is conserved across
  E. coli / Yersinia / Mycoplasma. But ClpX specificity is demonstrably *not* universally conserved —
  *S. mutans* ClpXP recognises a tripeptide `LPF` that E. coli ClpXP does not (PMC5143411). Flag the ClpX rules
  as an assumption.
- **Flynn 2003 is not reachable headlessly** (§4.1) — but it is not lost. Its motif definitions and
  trapped-substrate census are main-text content behind the Cell Press paywall, so an authenticated browser
  session should retrieve them; only the fully automated OA routes fail. The same is true of the whole
  cross-bacterial trap pool (Lunge 2020, Bhat 2013, Graham 2013 all return `idIsNotOpenAccess`; Feng 2013 and
  Ziemski 2021 have no PMC record). Until that is done, the C-terminal channels rest on DEtox and the Lon
  consensus, both of which are open — which is the better basis anyway.
- **~40% of E. coli cytoplasmic proteolysis is attributable to none of ClpP, Lon or HslV** (Gupta 2024). So
  protease attribution is informative but incomplete, and measured instability does not imply that a known
  handle is responsible.
- **No proteome-scale experimental ΔG exists for any bacterium.** Tm is the practical stand-in; Tsuboyama 2023
  is explicitly rejected for direct use (see `docs/degradability_references.md` §7).
- **Delivery, permeability and efflux are unsolved for this modality in Gram-negatives.** CLIPPERs are
  plasmid-expressed peptides. The one reported small-molecule Gram-negative "BacPROTAC" is 1,669 Da with an
  artificially relocalised target. ADEP is the cautionary calibration: isolated E. coli ClpP is exquisitely
  sensitive, yet whole-cell MIC is >64 µg/mL because ADEP is an AcrAB-TolC substrate — and AcrAB-TolC
  (RamA/RamR-regulated) is the dominant constitutive RND pump in *K. pneumoniae*. Any target narrative from
  this axis should carry an efflux/permeability caveat, and DegP's compartment is attractive partly *because*
  it sidesteps the inner membrane entirely.

---

## 10. Suggested extensions (not yet specified as tracks)

- **PPI-network connectivity as a first-class covariate.** It was the single most informative feature class in
  Nagar 2021, above every sequence motif, and it also encodes the "buried in a complex" intuition that track
  3.6 only approximates. Available from STRING via the E. coli ortholog.
- **Re-fit the Won 2024 Lasso on E. coli `k_deg`.** The published coefficients come from 54 mycobacterial
  proteins; the same feature recipe over ~3,200 continuous E. coli labels would be far better powered, and the
  72-protein screen then becomes an orthogonal held-out test set. This is the highest-value modelling step and
  it needs no new data beyond §4.
- **tmRNA-tagging propensity** as a 3.2 sibling: rare-codon density at the 3′ end, RNAfold MFE across the stop
  codon, internal anti-Shine–Dalgarno hits. Gupta 2024's Δ*smpB* condition supplies the labels.
- **σ32 / σS regulon membership** via the E. coli ortholog (RegulonDB) as a proteostasis-context covariate.
- **A ClpS-recruiting handle is unexplored.** ClpA and ClpS are both present in Kp and the ClpS N-degron pocket
  is structurally defined (Schuenemann 2009), but there is zero chemical precedent — the most interesting open
  medicinal-chemistry option, and one this axis could help motivate.
- **Latent-degron simulation.** Humbard 2013's finding that most N-degrons are generated by prior
  endoproteolysis suggests scoring the degrons that *would* exist at plausible internal cleavage sites, which
  is what DEGRONOPEDIA does for eukaryotes. Speculative, but it addresses the axis's largest known sensitivity
  gap.

---

## 11. E. coli-first build order, and corrections to this report (2026-08-05)

The axis is now to be built **E. coli first**. That decision triggered a dedicated E. coli literature and
resource pass, which changed enough that parts of this report above are now wrong. Corrections first.

### 11.1 Corrections to §§1–10

| Claim above | Correction |
| --- | --- |
| Assembly state is "the weakest aspect", ~0.5% coverage | True for **Kp** (30 PDB structures). **False for E. coli**, where it is a different regime entirely: 1,768 proteins with a PDB cross-reference, **7,671 biological assemblies** whose `oligomeric_state` RCSB returns literally as `"Homo 2-mer"`, 324 curated Complex Portal complexes **with explicit stoichiometry**, and — via PDBe `interface_residues` — actual N-/C-terminal burial **in the biological assembly** with partner identity, so self-burial in a homomer is separable from burial against another subunit. Plus Seq2Symm for 100% coverage. |
| Schweke 2024 is an E. coli homo-oligomer resource | It is **E. coli O157:H7**. Intersection with UP000000625 is **exactly zero**; 820/2,181 by gene name only. Use **Seq2Symm** instead. |
| The recognise gate has four rule channels | It has **nine**, in **two architectures**. The four envelope channels (RseP, GlpG, DegS, DegP/Prc) were missed because they are irrelevant to a ClpC-centric reading. See `docs/03_degradability.md` §3.2. |
| The initiation-region requirement gates the axis | It gates the **five AAA+ channels only**. The four envelope proteases are endoproteolytic — they cut rather than thread, and need an exposed cleavage site instead. Gating the periplasm on initiation regions would be a category error. |
| Kp machinery accessions were found "via ortholog" | Re-derived independently by pairwise alignment: 16/16 pairs at **76–99% identity** vs a **20–32% random-pair baseline**. The assignments stand but no longer depend on the ortholog table. |
| Only FtsH is undeletable | In **both** organisms, **three** proteases are experimentally essential: **FtsH, DegS and RseP** — all in the envelope, and none with any degrader precedent. |
| "Only three handles have degrader precedent (ClpXP, FtsH, DegP)" | **FtsH has none.** Two have demonstrated chemistry (ClpXP, DegP); Lon has an engineered genetic handle (mfLon); the other six are mechanistic priors. |
| No proteome-scale ΔG exists for any bacterium (§9) | Still true, but there is now a partial substitute: **proteome-wide acid-stability pH50 for 1,675 E. coli proteins**, which also stratifies by compartment. |

### 11.2 The one file that justifies E. coli-first

Gupta 2024's Table S1 is richer than §4.2 recorded. It is **3,262 proteins × 14 conditions**, keyed on
**UniProt accessions** (`sp|P00350|6PGD_ECOLI`), and the conditions include **`clpP N-lim6`, `lon N-lim6`,
`hslV N-lim6`, `Triple (ΔclpPΔlonΔhslV)` and `smpB N-lim6`**.

A per-protein change in half-life upon deleting a specific protease is not a feature of a degradability model —
it is a **ground-truth label for that protease's channel**. That makes E. coli the only organism in which the
axis can be trained and validated at all, and it should be the first thing fetched.

### 11.3 Build order

1. **Labels** — Gupta 2024 Table S1. Nothing else before it.
2. **Real termini** — TRAINSPOTTER (1,082 N-termini, 729 proteins, UniProt+b-number) and the PeptideAtlas
   `coordinate_mapping.txt` for proteome-scale empirical N-/C-termini. The N-degron and C-degron channels are
   currently scored on annotated termini that are frequently wrong.
3. **Rules** — DEtox, Cragan Lon, N-FIVE, plus the four envelope rules now specified.
4. **Structure/network features** — terminal pLDDT + SASA (free, on disk); STRING degree (taxon 511145, 91%
   coverage); RCSB `oligomeric_state` + PDBe interface residues.
5. **Flux** — Li 2014 synthesis rates; eSOL solubility.
6. **Calibration** — Klimecka 2021, which uniquely separates adaptor binding from protease engagement.

### 11.4 Two things to design around

**A competitor exists.** ProHL (*J Microbiol Biotechnol* 2026) is a deep-learning bacterial protein half-life
classifier validated in E. coli. Read its corrigendum before benchmarking. Our differentiator is not
half-life prediction — it is per-channel, compartment-routed *recruitability*, with an explicit
demonstrated-vs-prior distinction. Half-life is an input to that, not the output.

**One published result must be excluded from any label set.** LAMP-D (*JACS* 2025) degrades NDM-1 in a
Gram-negative and would look like a positive example, but it works by **photo-oxidative backbone scission,
explicitly without recruiting host proteolytic machinery**. It is a negative control for this axis.

### 11.5 Search-method note

The WebSearch budget for the session that produced this pass was exhausted partway through. The remainder ran
on the **Europe PMC REST** and **PubMed E-utilities** APIs plus direct `curl`, which are strictly better for
systematic date-bounded literature sweeps and remain available. Any future pass should start there rather than
with general web search. Caveat: those APIs index abstracts and full text, so a purely web-surfaced resource —
a lab's unpublished supplementary site, a bare GitHub data repo — could still be missed.
