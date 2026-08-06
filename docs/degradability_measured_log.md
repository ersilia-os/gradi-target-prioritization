# Measured-degradability layer — run log (docs §3.3)

Implemented `scripts/10e` + the three slides `10l`/`10m`/`10n`. Covers only the **measured**
evidence — what was observed in a lab — as distinct from the computed layers (`10b` degrons,
`10d` disorder) that `10j`/`10k` already plot.

## What was ingested

`10e_measured_turnover.py` → `output/results/<org>/<prefix>_deg_measured.csv`.

| source | what it measures | Ec proteins |
|---|---|---|
| **Nagar 2021** | pulsed-SILAC half-lives, 3 classes | 1,148 |
| **Gupta 2024** | half-lives × 13 conditions **+ protease-KO panel** (ΔclpP/Δlon/ΔhslV/triple/ΔsmpB) | 3,259 |
| **Niwa 2022** | GroE depletion in WT / Δlon / ΔclpPX / ΔhslVU | 2,005 |

Nagar was already parsed by `D.load_nagar()` but had only ever been used in memory as a validation
label; `10e` materialises it for the first time. Coverage: **Ec 3,291/4,403 (74.7%)**,
**Kp 2,767/5,728 (48.3%)** — Kp entirely by ortholog transfer, under the 55.5% orthogroup ceiling.

**MacKrell 2026 was deliberately not ingested.** Its `sd02`/`sd03` carry raw abundance timecourses
with no published half-life column, and `sd04` is ML predictions. Fitting decay curves ourselves
would manufacture numbers the paper did not publish. Deferred, not forgotten.

### Two derivations are ours, not the papers'

- `gupta_log2_<protease>` = log2(KO half-life / WT half-life) at N-lim6; attribution at a
  **0.5 log2** (1.41×) stabilisation threshold. Gupta's published categories are not a column in
  Table S1. Our net is wider than theirs — **ClpP 189 / Lon 130 / HslV 29 / additive 74 /
  redundant 224**, against their ClpP 64 / Lon 14 / HslV 1 / additive 82 / redundant 41 — because
  they applied replicate-level significance testing the table does not expose. Prefer the continuous
  `gupta_log2_*` columns when the distinction matters.
- `niwa_rescue_<protease>` = Foldchange(KO) − Foldchange(WT). Niwa's fold-changes are GroE-on vs
  GroE-off *within* each background, so a blunted depletion in a protease KO implicates that
  protease.

## Slides

| script | output | content |
|---|---|---|
| `10l` | `10l_activator_{kp,ec}.png` | ADEP4/ONC212 — the only measurement of partnerless ClpP |
| `10m` | `10m_turnover_{kp,ec}.png` | half-lives and protease attribution |
| `10n` | `10n_landscape_{kp,ec}.png` | cross-dataset coverage, overlap, ESM-C map, AUROC accountability |

## What the plots actually show

Several results are uncomfortable and are plotted as-is rather than smoothed:

- **The two half-life datasets barely agree — Spearman ρ = 0.11** across 866 shared proteins
  (Gupta 2024 vs Nagar 2021). Different methods, different dynamic ranges; Gupta's values cluster in
  10–100 min while Nagar's span five orders of magnitude. Do not treat "measured half-life" as one
  quantity.
- **The two activator readouts are nearly orthogonal — ρ = 0.06** between ADEP4 abundance depletion
  and ADEP4 cleavage-peptide evidence (n = 288); only 16 proteins satisfy both. The two activators
  agree far better with each other (ADEP4 vs ONC212 abundance, **ρ = 0.51**, n = 414, 33 both
  depleted). So the disagreement is between *readouts*, not between chemistries.
- **ClpP-attributed proteins are shorter-lived, but only on one dataset.** Median Nagar half-life
  163 min vs 254 min for other measured proteins (Mann-Whitney p = 3.2e-5). On Gupta's *own*
  half-life the same proteins sit at 0.69 h against a 0.74 h proteome median — essentially no
  separation. The KO panel and the half-life column in the same paper do not tell the same story.
- **The growth-correction cliff: 1,148 → 907 → 54.** Of Nagar's measured half-lives, 241 are
  infinite (never measurably decayed), and after correcting for dilution by growth only **54**
  E. coli proteins retain a positive proteolytic half-life. Most E. coli proteins are diluted out
  faster than they are proteolysed, so their turnover carries almost no proteolytic signal.
- **ClpXP rescue is rare in the Niwa assay** — 9 proteins, against Lon 67 and HslUV 87.
- **Computed features predict measured turnover weakly**: of the twelve strongest disorder features,
  6/12 bootstrap intervals cross 0.5, and the best (`cterm_init_len_70`, `cterm_plddt_15`) reach
  only ~0.59. Consistent with what `10j` already reported.
- **The named targets reproduce the pre-registered values exactly** — AcpP 1.00, DnaK 0.90,
  GyrA 0.65, **GyrB 0.0015**. GyrB maps cleanly to Kp, so its near-zero score is absence of signal,
  not absence of mapping.

Totals: any measured evidence **Ec 2,815/4,403 (64%)**, **Kp 2,569/5,728 (45%)**; Clp-specific
evidence (Gupta ClpP ∪ Niwa ClpXP ∪ activator) **714 Ec / 700 Kp**.

## Caveats carried into every panel

- **No Kp protein has ever been measured.** Every Kp value is transferred — from E. coli by
  orthogroup, or from *S. aureus* by cross-phylum RBH at a median 42% identity. `10l` panel 5 and
  `10n` panel 6 exist so this is never off-screen.
- The activator layer reaches only **608 Kp / 609 Ec (11–14%)**.
- The two-bar composite (`clpP_complex_substrate`, `partnerless_clpP`, `bar_disagreement`) is
  **prose only** — `10i_degradability_merge.py` does not exist. These slides plot measured inputs
  and deliberately compute no composite.
- Flynn 2003 and Neher 2006 — the two weight-1.00 ClpXP trap sets — are still placeholders. Their
  PDFs are in `data/raw/other/degradability/literature/`; extracting them would be the single
  biggest improvement to the measured layer.
- The legacy `data/raw/legacy/clp_substrates/*.tsv` stand-ins (45 and 35 rows) are **not** used
  anywhere here; the real Nagar data supersedes them.
