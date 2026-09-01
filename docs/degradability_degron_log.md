# Degron / motif layer — run log (docs §3.2)

Covers the **computed sequence-degron** track: `scripts/10b_degrons.py` and the two slides built from
it, `10k_degron_plots.py` and `10o_degron_relevance.py`. Scope of this run: **E. coli K-12 only**.

`10b:29` promises to copy its validation into `docs/degradability_log.md`. That file was never
written, so until now the numbers survived only as comment blocks in `src/degradability.py`
(`:561-583`, `:514`). This log is that file, under a name that matches the other axis logs.

## What was measured

Source: `output/results/ecoli/ec_deg_degrons.csv` (4,403 proteins, 29 columns, written Aug 4).
Labels: Nagar 2021 pulsed-SILAC turnover classes, **direct** for E. coli — 1,148 labelled proteins,
baseline P(fast) = 0.063.

| motif | proteins | % of proteome | labelled hits | fast among them | OR | **95% CI** | weight |
|---|---:|---:|---:|---:|---:|---|---:|
| `cterm_cm1_strict` (φ-A-A$) | 11 | 0.25% | 4 | 2 | **15.24** | 2.60 – 89.44 | 1.00 |
| `cterm_cm1_broad` (-A-A$) | 42 | 0.95% | 14 | 2 | 3.02 | **0.76 – 11.99** | 0.40 |
| `cterm_cm2` (MuA) | 115 | 2.61% | 28 | 0 | 0.25 | 0.02 – 4.20 | 0.00 |
| `nterm_nm2` | 104 | 2.36% | 26 | 1 | 0.87 | 0.16 – 4.57 | 0.00 |
| `nend_primary_destabilizing` | 455 | 10.33% | 98 | 8 | 1.44 | 0.68 – 3.03 | 0.00 |
| `nend_imet_cleaved` | 1,793 | 40.72% | 543 | 45 | 1.92 | **1.18 – 3.13** | 0.00 |

Point estimates reproduce `data/processed/ecoli/degradability/ec_degron_validation.csv` exactly —
`10k` recomputes them from Nagar and warns on any disagreement; none fired.

### The intervals change the reading, and 10b never had them

`10b` reports point estimates only. Adding a Wald interval on log-OR reorders the picture twice:

1. **`cterm_cm1_broad` carries weight 0.40 on an interval that covers 1** (0.76 – 11.99). Its OR of
   3.02 rests on 2 fast-turnover proteins out of 14 labelled hits. The weight is not *wrong* — the
   direction is right and the mechanism is real — but it is not statistically supported, and it
   supplies 42 of the 42 proteins that get any motif contribution at all.
2. **`nend_imet_cleaved` is the only feature whose interval excludes 1** (1.18 – 3.13) and it is
   correctly weighted **0** — at 41% prevalence it is a prior, not a motif. Keeping it at 0 while
   `cm1_broad` sits at 0.40 is defensible only on mechanism, not on these numbers. Worth revisiting
   when `10e` measures the motifs against the activator labels.

`cterm_cm1_strict`'s interval excludes 1, but it rests on **4 labelled proteins**, 2 of them fast.

## The three findings that matter

**1. Only 42 of 4,403 proteins (1.0%) carry a weighted motif.** Not 53 — `cm1_strict` is a strict
subset of `cm1_broad` (every φ-AA ending also ends -AA), so summing the two double-counts all 11.

**2. `degron_score` is terminal exposure, not motifs.** ρ(`degron_exposure_score`, `degron_score`)
= **0.999** over 4,371 proteins with a model. The composite is `0.45·motif + 0.55·exposure`, but the
motif term is non-zero for 42 proteins, so in practice the score is the exposure.

**3. Not one of the 11 ssrA-like proteins has an exposed C-terminus.** Best is `ppx` at 0.49,
against a 0.5 bar; the proteome-wide rate is 30% of termini. And the two with the highest
`degron_score` — `ppx` (0.63) and `dgcC` (0.56) — get it from their **N**-terminus (0.74 and 0.98),
i.e. from the opposite end to the motif they carry. The full set:

```
ppx  ···PEIAA  C 0.49  N 0.74      yniA ···RLLAA  C 0.04  N 0.04
dgcC ···TEVAA  C 0.05  N 0.98      dgcM ···RVLAA  C 0.03  N 0.24
ytfR ···NAIAA  C 0.36  N 0.40      pxpA ···IVVAA  C 0.03  N 0.04
ulaC ···TNAAA  C 0.13  N 0.04      yhfA ···EVVAA  C 0.03  N 0.03
recN ···ELLAA  C 0.12  N 0.05      sodB ···KNLAA  C 0.03  N 0.03
hycD ···SLLAA  C 0.05  N 0.10
```

## …and why none of it decides the axis (`10o`)

The Gr-ADI handle is **partnerless activated ClpP**. A degron is read by an AAA+ unfoldase or its
adaptor — never by ClpP itself — and the activator's mechanism is to occupy the ClpX/ClpA docking
cleft. So every motif in the table above is read by the subunit this modality removes.

**Read the provenance before the numbers.** The only activated-ClpP proteomics that exists — Conlon
2013 (ADEP4) and Jacques 2020 (ONC212) — was measured in ***S. aureus***. There is no equivalent
dataset for E. coli or K. pneumoniae; the Gr-ADI SoW notes that gap itself. Every E. coli value below
is a DIAMOND-RBH ortholog transfer (`10c`), covering **609/4,403 (13.8%) at a median 42% identity**.
So these are S. aureus measurements attributed to E. coli orthologs, and `10o` labels every panel
that uses them accordingly rather than presenting them as native E. coli data.

Measured consequences, all from `10o`:

- Of the **609** E. coli proteins carrying transferred activated-ClpP evidence, only **8** carry a
  weighted degron motif, and only **3** of those are cleaved. The two evidence layers barely
  intersect, so the data cannot adjudicate degrons either way — which is why the mechanism argument,
  not a statistic, is the load-bearing one.
- Cleavage rate by motif status is flat: baseline 334/609 = 55%; `cm1_broad` 3/8, `cm2` 8/11,
  `nm2` 7/18, `nend_destab` 24/57, `imet` 184/311 — **every Wilson interval covers the baseline**.
  `cm1_strict` is omitted from that panel: only 2 of its 11 proteins have activator evidence.
- **Nothing computed predicts activated-ClpP cleavage.** AUROC with 1,000-fold stratified bootstrap:
  every feature tried — motif score, `degron_score`, terminal exposure, global disorder fraction,
  longest internal disordered run, initiation-region length, sequence length — has an interval
  covering 0.5, and most point estimates sit *below* it. Part of the ceiling is the readout itself:
  within these 609, ADEP4 abundance and ADEP4 cleavage agree at **ρ = 0.03**.
- **ClpC and McsB are absent from the proteome** (measured against the 00a proteome TSV, not
  asserted). ClpP, ClpX, ClpA, ClpS, SspB, Lon, HslU, FtsH, SmpB are all present. Every published
  BacPROTAC recruits ClpC1, so all of them describe a machine this organism does not have.

Bottom line for the composite: degron evidence belongs to the `clpP_complex_substrate` bar
(Regime A, selection) and contributes nothing to `partnerless_clpP` (Regime B, validation). It must
not be averaged into the latter — see the two-bar revision in `docs/03_degradability.md`.

## Slides

| script | output | content |
|---|---|---|
| `10k_degron_plots.py` | `10k_degrons_ec.png` | prevalence · OR-with-intervals · motif × accessibility · exposure map · what drives the score · the ssrA-like set |
| `10o_degron_relevance.py` | `10o_degron_relevance_ec.png` | who reads a degron · cleavage rate by motif · AUROC vs the activator · evidence overlap · machine inventory · degron_score ECDF, cleaved vs not |

Both default to `--organism ecoli`; `--organism kpneumoniae` runs but was **not** produced or checked
in this pass. For Kp the Nagar labels are ortholog-transferred, so panel 2 of `10k` is weaker there.

## Caveats

- **The CSVs are stale.** `ec_deg_degrons.csv` is from Aug 4 and predates the Aug-6 implementation of
  Flynn's real N-M1 (`T-X-K-[ILV]`, 1–4 residues in) and N-M3 consensuses in `src/degradability.py`.
  Neither column exists in the data, so neither slide plots them; both scripts derive the motif list
  from the file's own header and annotate the omission, so **re-running `10b` picks them up with no
  code change**. Prevalence is already known from the source comment: NM1 22/4,403 (0.50%),
  NM3 196/4,403 (4.45%). Both are weighted 0 pending measurement.
- `10b` uses `TERM_WINDOWS = (10, 20)`; the spec (docs §3.1) asks for 15/30/50. The 15/30/50 windows
  do exist, in `10d`'s `_disorder.csv`, which is what `10j` and `10o` panel 3 use.
- The activated-ClpP layer is *S. aureus* data transferred by RBH (see above), so it is biased
  toward conserved core proteins. That is the honest ceiling for this evidence, not a fixable gap —
  a native Gram-negative activated-ClpP proteomics dataset does not exist.
- `10o` panel 1 is a **drawn schematic**, the only non-measured panel across either slide. Everything
  else on both figures is computed from files in `output/results/ecoli/`.
