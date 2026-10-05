# `plotting/` — the GraDi consortium deck

The record for the presentation figures. `CLAUDE.md` is the map; this file is what each figure
claims, the numbers it printed when it was built, and the caveats that must travel with it.

**Audience: the GraDi consortium** (Prof. Erick Strauss, Stellenbosch), not an internal methods
review. Kp-first; E. coli appears where the contrast earns its place; S. aureus appears only as the
source of the measured degradability labels.

**These are not the per-axis diagnostics.** Those live in `scripts/plots/` and are cited from
`docs/<task>.md`. `plotting/` is cross-axis and presentation-only, and the two must not be merged:
a diagnostic defends a stage, a deck figure answers "which protein should we pursue".

## The constraint: no composite score, anywhere

Three axes each REMOVED a composite, and reintroducing one in a plotting script would undo all
three invisibly — a ranked list looks equally plausible whatever weights produced it.

| axis | what was removed | why |
|---|---|---|
| essentiality | `essentiality` / `essentiality_source` | mixed units (measured 1.0 vs predicted 1.0); on Kp a verbatim copy of `screens_ess_mean` |
| studiedness | the 0-1 composite (2026-09-22) | the two halves double-counted, r 0.64–0.70 |
| pockets | `druggability()` | "no defensible weighting exists across a weak prior, a sparse measurement and a third party's model" |

So the shortlist is a **cascade of stated predicates**, defined once in `plotting/filters.py`, each
a single comparison with a written rationale and its own count.

## The cascade, as measured on 2026-10-04

| step | rule | n |
|---|---|---|
| start | Kp HS11286 proteome | 5,728 |
| not membrane | `localization != "cytoplasmic_membrane"` | 4,294 |
| selective | `has_human_ortholog == False` | 3,497 |
| essential | `screens_ess_mean >= 0.299` (proteome-wide 90th pct) | 235 |
| degradable | `adep4_prob >= 0.328` (`BASE_RATE_THRESHOLD`) | **59** |

**Leave-one-out is what matters**, because the filters are conjunctive and the funnel
overstates the early ones: without `not_membrane` 63 · without `selective` 104 · without `essential`
792 · without `degradable` 235. So **essentiality does the work, and `not_membrane` costs almost
nothing (59 → 63)** — the degradability model has already learned the membrane mechanism, which is
independent confirmation of the compartment panel of `degradability.png`.

Survivors are led by `hns`, `infA`, `greA`, `nusG`, `nusA`, `xseB`, `minE`; 49 of 59 carry a gene
symbol, median 12 papers, **0 with a potent ligand**.

## The figures

| script | output | the claim |
|---|---|---|
| `landscape.py` | `landscape.png` | the unit of analysis is a complete embedded proteome; membrane proteins form their own island, 21% is COG-unclassified |
| `coverage.py` | `axis_coverage.png` | what each axis can say, and where a 0 means *unknown* |
| `degradability.py` | `degradability.png` | 0.5 is the wrong cut; secreted proteins ARE reachable; the Sa→Kp extrapolation priced |
| `essentiality.py` | `essentiality.png` | three genuinely different opinions, no verdict |
| `interest_panel.py` | `interest_panel.png` | the consortium's panel is the right biology and the wrong chemistry |
| `novelty.py` | `novelty_vs_essentiality.png` | the unexplored-and-essential quadrant, with its confound priced |
| `selectivity.py` | `selectivity.png` | neither orthology method alone; the union is conservative |
| `ligandability.py` | `ligandability.png` | effort vs potency — a 0 is not always a 0 |
| `pockets.py` | `pockets.png` | structure is not the bottleneck; pocket scores add nothing over length on the anchor |
| `essentiality_agreement.py` | `essentiality_agreement.png` | two screens correlating, the 12-genome essentialome, and ProteomeLM-Ess ranked |
| `consensus_evidence.py` | `consensus_evidence.png` | what `<axis>_consensus` and `<axis>_evidence` are, worked on essentiality |
| `uniprot_annotation.py` | `uniprot_annotation.png` | E. coli is annotated and Kp is not — with human orthology as the control |
| `studiedness_essentiality.py` | `studiedness_essentiality.png` | why a zero differs by organism, what the literature actually talks about, and why "unexplored" and "essential" are not independent |
| `degradability_cv.py` | `degradability_cv.png` | does the model work? ROC and PR with the across-fold band, and where the scores land per proteome |
| `enrichment.py` | `enrichment.png` | what kind of protein ClpP reaches: COG function and compartment/topology |
| `projection_pair.py` | `projection_pair.png` | one map read twice: compartment, then the top 250 essential |
| `druggability_consensus.py` | `druggability_consensus.png` | several views per axis become one column — and the pocket one is largely protein length |
| `cavities.py` | `cavities.png` | what a top-ranked pocket actually looks like, on three proteins |
| `shortlist.py` | `shortlist.png`, `shortlist_kpneumoniae.tsv` | five rules, 59 proteins, and what each rule cost |

## Look: npg palette, `article` style, black labels, no panel letters

Owner's instruction, 2026-10-04, in two parts. The deck uses stylia but none of stylia's ersilia
branding.

**Palette** — `stylia.CategoricalPalette("npg")` through **`plotting/palette.py`**, never
`stylia.NamedColors()` (plum/orange/mint). Semantic names: `PRIMARY` · `SECONDARY` · `TERTIARY` ·
`ACCENT` · `MUTED` · `INK` · `BACKDROP` · `SEQUENTIAL`, because `NPG[5]` at a call site says nothing
about intent. `SEQUENTIAL` keeps heatmaps npg too — `BuPu` and `viridis` are not.

**Style** — `stylia.set_style("article")`, not `("ersilia")`. The ersilia style paints every text
element and spine plum (`#50285A`); article gives black. `PAL.INK` is plain black for in-plot text.

**No `abc=` panel letters.** `stylia.label()` can stamp A/B/C on panels; this deck does not. Panels
are referred to by what they show — a letter goes stale the moment a panel is reordered. The
per-axis figures under `scripts/plots/` still use `abc=`; this is a `plotting/` rule only.

`plotting/palette.py` is the second shared module, beside `filters.py`, and centralising it is a
deliberate departure from `scripts/plots/`'s copy-paste convention: a palette stops being
per-script taste once it is a rule, and a deck rendered half in one palette and half in another is
a defect no single script can see.

**`src/localization.py:LOC_CLASS_COLOR` is already npg** (#E64B35, #00A087, #3C5488, #F39B7F,
#8491B4, #7E6148 are ggsci `pal_npg`) and is imported rather than redefined, so the deck and the
localization axis's own plots cannot drift apart.

**Axis labels are capitalised** — a standing convention these figures got wrong on the first pass.

## Measured values these figures printed

Each reproduces a number `CLAUDE.md` already records, which is how the scripts were validated.

- **druggability_consensus**: ligand inputs agree with each other at ρ 0.23–0.58 but land at
  0.17–0.98 against `ligands_consensus`; pocket inputs 0.20–0.29 against each other, 0.16–0.83
  against `pockets_consensus`. ρ(`pockets_consensus`, length) **0.663 / 0.638 / 0.654**; AUROC vs
  a measured PDB ligand — length alone **0.673 / 0.652 / 0.657**, P2Rank *within length deciles*
  **0.489 / 0.565 / 0.614**. Reproduces `docs/pockets.md` to within the decile binning.
- **cavities**: irp2 consensus 0.981 / 2,035 aa / 1 PDB ligand · acrB 0.968 / 1,048 / 3 ·
  gyrB 0.891 / 805 / **14**. All three `pockets_evidence` 3.
- **coverage**: GO-slim 73.7%, COG 79.1%, geptop 92.2%, studiedness ambiguous 1,961 (`no_hit`
  1,058 + `below_floor` 903), ligands assayed 275 of 5,728.
- **degradability**: ADEP4 max 0.868, cut 0.328 selects 941 (16.4%), at 0.5 selects 220; ONC212
  33.1%; ρ(adep4, onc212) **0.895**. Measured hit rates cytoplasm 0.180/0.275, membrane
  0.031/0.105, extracellular 0.049/0.346.
- **essentiality**: Geptop 3,799 tied at 0 (66.3%); ρ screens/ProteomeLM 0.32, ProteomeLM/Geptop
  0.44; top-500 shared 333–369.
- **interest panel**: 40 of 43 symbols matched (missing `lpxC`, `kdsC`, `kdsD` — a naming gap, not
  an absence). Median percentile: essentiality 96.0, studiedness 80.3, ADEP4 **49.6**, ONC212
  **42.4**, conservation 78.0, pocket 58.5.
- **novelty**: 78 proteins are ≤5 papers AND essentiality top decile; 14 of them are on the
  shortlist. ρ(log papers, `screens_ess_mean`) = **0.122** — far weaker than `geptop_ess`'s
  documented 0.38–0.50, so ranking on `screens_ess_mean` buys a less circular quadrant.
- **selectivity**: human orthologs Kp 951 / Ec 838 / Sa 624 (union); **RBH-only 180 on Kp**.
- **ligands**: potent Kp 113 / Ec 96 / Sa 78; screened-but-nothing-potent Kp 64; never looked 5,453.
- **projection pair**: the top 250 essential have mean pairwise map distance **24.3 against 50.9**
  for random draws of the same size — they genuinely cluster, measured before being drawn. They are
  COG J (translation) 90 of 250, and 78.8% cytoplasmic against a 60.5% proteome background, with
  outer membrane enriched (8.4% vs 4.1%) and extracellular absent entirely.
- **essentiality agreement**: three levels. (1) **Two measured CRISPRi screens** (Rousset 2018 vs
  Wang 2018, E. coli, n=3,719) — **Pearson +0.908 but Spearman +0.393**, and the decomposition is
  the finding: ρ **+0.726 among Keio's 274 essentials** against **+0.252 among the 3,445
  dispensables**. The screens agree about what is essential, not about how dispensable the rest is;
  the essential tail anchors the Pearson. **Never quote the Pearson alone.** (2) **Cross-species
  essentialome**: of 3,170 Kp proteins with a call across 12 Enterobacteriaceae genomes, **447 are
  essential in ≥1** — **159 in all twelve, 65 in eleven, 89 in exactly one**. A solid core block
  over a sparse species-specific tail, which is why screens are never merged. (3) **ProteomeLM-Ess ranked**, top
  twelve named: `leuS` 0.9997, `metG`, `aspS`, `rpoB`, `murG`, `lpxL`*, `argS`, `rpoC`, `valS`,
  `lptG`*, `rpoD`, `gyrB`* — tRNA synthetases, RNA polymerase, peptidoglycan and gyrase, and **three consortium panel targets (`lpxL`, `lptG`, `gyrB`) in the top twelve** — reported in the run log, not highlighted on the panel. Median 0.1001; 204
  of 5,728 above 0.9. Geptop is not drawn but is measured: ρ **+0.439** against ProteomeLM-Ess,
  **+0.663 where Geptop > 0**, with **66.3% ties at exactly 0** holding the first number down. **Panels 1–2 read legacy v1 tables** (`output/results/`,
  2026-07) — the only route to two continuous screens on one key and to the 12-genome matrix; panel
  3 is v2.
- **consensus + evidence**: the convention explained on the one axis that has several predictors
  *and* several experiments. (1) The consensus is a **blend, not a copy** — on Kp the three
  predictors agree with each other at ρ 0.32–0.44 but each sits at **0.69 / 0.78 / 0.73** against
  the consensus. (2) **Same ladder, opposite provenance** — tiers are Kp **1,003/752/3,973** and Ec
  **110/929/3,364**, but Kp's 4,725 covered proteins are **100% proxy** (DEG has no Klebsiella; the
  calls come from ≥95% counterparts in ATCC 43816 / ECL8 / RH201207, and **not one tier-3 protein is
  measured on HS11286**) where Ec's are **0% proxy** — two DEG datasets on the anchor strain plus
  the OGEE label. (3) **The two columns are near-orthogonal and mildly anti-correlated**,
  ρ **−0.216 Kp / −0.349 Ec**: tier 3 is 65% below consensus 0.5 because concordance is tested with
  equality, so corroborated NON-essentials land there, while tier 2 holds the most high scorers
  (18.0% above 0.9). Concordance uses a **base-rate cut, never 0.5** (0.824 Kp / 0.738 Ec).
- **UniProt annotation**: `proteomes_evidence` tiers 1/2/3 are **Kp 2,205 / 2,469 / 1,054** against
  **Ec 0 / 648 / 3,755** — Ec is **85.3%** at the top tier where Kp is **18.4%**, and Kp has 2,205
  proteins at tier 1 where Ec has none. Same fact in other columns: SwissProt-reviewed 0.1% vs
  100%, gene symbol 63.4% vs 100%. **The control panel is what makes this a statement about
  curation** — five comparability measures, Kp vs Ec: in an orthogroup **92.3 / 96.5**, conserved in
  ≥50% of the 28-species panel **57.3 / 66.6**, core (≥90%) **8.6 / 10.6**, human ortholog
  **16.6 / 19.0**, has a paralog **36.2 / 32.9**. Every one within ~9 pp, at median identity to the
  human ortholog of 33.5% vs 33.0% — so the 4.6× annotation gap is about work done, not about what
  is there to find. (`bacterial_panel_orthologs` is a **fraction** over 28 proteomes, not a count.) `is_reviewed` is the blunter number and the worse
  column — degenerate per species, which is why `proteomes_evidence` replaced it.
- **studiedness vs essentiality**: ρ(`studiedness_consensus`, `essentiality_consensus`) is
  **+0.332 Kp / +0.440 Ec** over all proteins, **+0.433 / +0.448** at `studiedness_evidence >= 2`,
  and **+0.639 / +0.456** at evidence 3 — it **strengthens** as the literature evidence improves,
  so it is not an artifact of the zero block. Controlling for protein length changes nothing
  (+0.364 Kp, +0.437 Ec). Median essentiality across studiedness bins: Kp 0.41 → 0.82 (8 bins, the
  34% zero block is one tie), Ec 0.39 → 0.85 (10 bins). Cross-species: Kp median 4 papers with
  **34.4% at zero** against Ec's 6 and 0.7%; evidence tiers **1,961/3,403/364** vs **29/519/3,855**,
  i.e. Ec is 87.6% on its own curated literature where Kp is 6.4%. **The deck implication**: the
  "unexplored and essential" quadrant in `novelty_vs_essentiality.png` is genuinely thinner than it
  looks. Association only — people study essential genes, and easily-studied genes are easier to
  call essential. **The three panels do not rest on one source and the axis labels now say so**:
  the violins and the evidence tiers are UniProt only, while the relationship panel uses
  `studiedness_consensus`, a rank mean over `n_papers_uniprot_own`, `n_papers_uniprot_prokaryotic`
  and `n_papers_pubtator_prokaryotic`. The source does not change the SIGN but on Kp it changes the
  MAGNITUDE — ρ against `essentiality_consensus` is **consensus +0.332 · uniprot_prokaryotic +0.312
  · uniprot_own +0.126 · pubtator +0.417**, a 3.3× spread with **PubTator strongest**, matching the
  axis's own held-out control (0.3722 vs 0.3398). On Ec all three agree closely (0.381–0.417).
  `n_papers_pubtator_prokaryotic` carries **real blanks, not zeros** (Kp 2,149 = 37.5%; Ec 190),
  which `percentile_consensus()` skips row-wise — so the consensus rests on three columns for some
  proteins and two for others. The middle panel ranks every E. coli protein by PubTator mentions,
  log-log, naming the top 10 — **rpoB 2,541 · recA 633 · rpoS 439 · crp 376 · ftsZ 370 · hscA 353 ·
  dnaK 353 · dnaA 352 · sdiA 345 · lacI 330**, against a median of 10; the top 10 hold 6.0% of all
  mentions in the proteome. **It is E. coli and not Kp on purpose**: Kp's PubTator counts are
  borrowed from whichever donor the transfer found, usually the E. coli protein, so the two species'
  top values are literally identical and ranking Kp would draw E. coli's literature under a
  Klebsiella label.
- **degradability CV**: cluster-grouped, 5 folds x 5 seeds, TabPFN-3.5 on ESM-C. ADEP4
  **0.8738 ± 0.0029 AUROC / 0.6103 ± 0.0091 PR** on n=1,677 (base 0.137); ONC212
  **0.7671 ± 0.0060 / 0.5803 ± 0.0075** on n=1,045 (base 0.246). Length-only baseline 0.775 / 0.680;
  cross-assay yardstick 0.877 / 0.815 (a yardstick, **not** a ceiling). `leakage_gap` is +0.001 and
  −0.005 — grouping cost nothing on this label set, which is a property of the labels, not a licence
  to drop grouping. **The drawn curve is not the quoted number**: one OOF realisation gives 0.8761 /
  0.7706 against the 5-seed means, both inside one SD. The shaded band is the **min-to-max across
  the 5 CV folds, not a confidence interval** (per-fold AUROC 0.849–0.892 ADEP4, 0.742–0.799
  ONC212). The third panel is deployment, not validation: median score and % above the cut are
  Kp 0.080/16.4%, Ec 0.072/14.4%, Sa-seen 0.054/13.8%, Sa-unseen 0.098/30.9% under ADEP4 — only
  Sa-seen is validated. The length and cross-assay levels are drawn on the ROC panel as
  **iso-AUROC reference curves** (binormal, area exactly 0.775 / 0.877 etc.) and labelled on the
  curve rather than in the legend — an AUROC is an area, so a scalar cannot be drawn on an ROC axis
  any other way, and these are LEVELS, not the baseline's measured path.
- **enrichment**: hits are the top decile (Kp 573 of 5,728). Compartment is the clean story —
  `wholly cytoplasmic` is the strongest enrichment in the table, every membrane/export feature is
  depleted (TM helix, signal peptide, beta-barrel, inner membrane, periplasm), and **`extracellular`
  is ENRICHED** (Kp +1.9 log2 under ONC212, OR 3.80). Function: COG `M` cell wall/envelope is the
  strongest depletion (−6.2), `G`, `N`, `E`, `C`, `I` all strongly depleted; `J` translation
  enriched. Two rows need a health warning — see below.
- **pockets**: PDB 569 / 1,893 / 595; own drug-like ligands 88 / 308 / 90; AlphaFill 1,533 / 1,196 /
  704. Length-controlled AUROC P2Rank Kp **0.496** (raw 0.615), fpocket Kp **0.434** (raw 0.511).

## Measured but not drawn

Reproducible from the deliverables; recorded here so they are not lost.

- **degradability ↔ pockets = −0.571 Kp / −0.557 Ec / −0.561 Sa**, falling to ≈ −0.22 controlling
  for length but surviving on all three. The two chemistry-facing axes pull apart: a protein with a
  good pocket is one ClpP is less likely to degrade, so a shortlist demanding both is fighting
  itself.
- **Protein length is a latent axis** across the prioritisation: ρ −0.68 (degradability), +0.66
  (pockets), +0.32 (studiedness) on all three species, while essentiality (−0.04) and ligands
  (+0.12) are nearly clean.
- **degradability ↔ studiedness flips sign under length control** (Ec −0.053 → +0.176; Sa −0.101 →
  +0.210) — a suppression effect.
- **`ligands_consensus` is ~97% ties at exactly 0 on Kp** (the zero block is pinned by design), so
  any Spearman involving it is dominated by that tie block and is not an ordinary correlation.

## Caveats that must travel with the deck

1. **Kp degradability rows are ranking hypotheses**, extrapolated from S. aureus labels. The premise
   that ESM-C cosine measures transferability is unvalidated without Gram-negative labels.
2. **The two activator columns are one opinion**, ρ 0.89, while the labels agree at ρ 0.52.
3. **`screens_ess_mean` is nine models reading the SAME ProtT5 embedding** — one correlated
   opinion, not nine votes. Comparable within a species, never across.
4. **A 0 in studiedness is `no_hit` or `below_floor`** and the shipped table cannot tell them apart.
   1,961 Kp proteins. Join `load_transfer()`.
5. **Pocket scores add nothing over protein length on the anchor.** Quote the length-controlled
   number.
6. **The ligands axis does not say "has an antibiotic"** — `pchembl_value` excludes MIC, which is
   why the ribosome looks empty there while leading the degradability ranking.
7. **The two largest enrichment bars are the weakest evidence in that figure.** COG `B` has 8
   members, and both `B` (median 97 aa) and `unclassified` (88 aa) are short-protein categories
   against a classified median of 297 aa, while `spearman(length, adep4_prob) = -0.677`. Within
   length deciles `unclassified` falls from a 3.85x to a 2.64x rate ratio — a real residual, but a
   third of the bar is protein size. Measured 2026-10-04.
8. **`enrichment.png` draws only categories significant under at least one activator** (Kp: 19 of
   25; `A D O R V W` carry no evidence either way). That is a legibility choice, not a statistical
   one — `--all` keeps them and the run log names every omitted category with its counts. It is
   also what makes the figure slide-shaped: 25 rows forced an aspect of 1.62, taller than 16:9.
   Final shape is a 3.17 band, matching `degradability_cv.png` so the deck reads consistently.
9. **Categories under 5 members are dropped from `enrichment.png`, as a correctness fix.** The
   Haldane +0.5 correction in `log2_or_adj` outvotes the data on a group of one: COG `Z`
   (cytoskeleton, 1 member, 0 hits, raw OR 0) comes out at **+1.58, reading as ENRICHED**. No row
   with >= 5 members flips. `src.degradability.MIN_ENRICHMENT_GROUP` carries the threshold.
10. **`src/interest.py` is an expansion of prose, flagged as a draft**, matched by gene symbol
   against a proteome whose `gene_name` is 63.4% covered.

## Traps found while building this

Four in the stylia layer, all silent, all now in every script's docstring:

1. **`REPO_ROOT` is `parents[1]` here, not `parents[2]`.** `plotting/` is one level deep.
   `parents[2]` resolves to the repo's PARENT and the ImportError names `src`, not the path — it
   reads as a broken conda env.
2. **`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.**
   `width=13` asks for thirteen slide-widths: a 4-gigapixel, 16 MB PNG, and the script exits 0.
3. **`stylia.label(..., xlabel=None)` writes the literal placeholder `"X-axis / Units"`.** Pass `""`.
4. **`stylia.set_style("ersilia")` must precede `NamedColors()`.** Without it `NamedColors()` returns
   `ArticleColors`, whose palette is `amber/cobalt/crimson/...` — every `NC.plum` raises
   `AttributeError`.

And four more from `cavities.py`, where the structures are images rather than plots. **Each fix
caused the next**, and every one of them is invisible once correct:

5. **Per-panel `zoom` destroys comparability.** Filling each frame independently made an 805-aa
   protein (gyrB) render LARGER than a 1,048-aa one (acrB), contradicting the lengths in the
   titles. The worker now renders twice: once to learn the camera distance each molecule needs,
   once with the largest of them for every panel.
6. **A common camera distance puts the smallest molecule in the depth fog.** gyrB came out visibly
   pale beside irp2 purely because it is smaller — the exact misreading the common scale exists to
   prevent. `depth_cue`, `fog` and `ray_trace_fog` are all off.
7. **A union crop box leaves the narrow panels empty on one side**, and a per-panel box re-expands
   the small ones and undoes the common scale.
8. **`imshow` STRETCHES an image to fill its axes.** So cropping each panel to its own width is
   only safe when the axes carry `width_ratios` equal to those widths. Without that the narrow
   panels are silently stretched back and the common scale is lost with no visible symptom.

And the known one, which did bite: **run plot scripts ONE AT A TIME.** stylia `rmtree`s the
matplotlib cache dir at import, so a tight loop over the scripts produces a *different* subset of
spurious failures each run.

## Run

```bash
P=~/miniconda3/envs/gradi/bin/python

for f in plotting/*.py; do $P "$f" || break; done   # ONE AT A TIME

$P plotting/shortlist.py          # must print 5728 -> 4294 -> 3497 -> 235 -> 59
$P plotting/cavities.py           # reuses cached renders; --refresh re-runs pymol
```

**`cavities.py` crosses a process boundary** into `gradi-pymol` via `conda run`, calling
`scripts/pockets/workers/pymol_render.py`. It is the only figure here that is not pure matplotlib.
Renders cache in `data/processed/pockets/scratch/renders/`, so a re-run is seconds unless
`--refresh` is passed. **Never activate `gradi-pymol` to run it** — the worker exists so the stage
does not have to.

Outputs land in `output/plots/presentation/`. Nothing here recomputes anything: every figure reads
the stages' own tables through `src/` loaders.
