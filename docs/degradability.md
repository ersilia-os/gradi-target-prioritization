# Stage 04 — degradability

`scripts/degradability/predict.py` · `scripts/plots/degradability.py` · loader `src/degradability.py`

One question per protein: **is it a substrate of activated partnerless ClpP?** A binary
**TabPFN-3.5** classifier on stage-01 ESM-C embeddings, trained on the two *S. aureus* activator
screens and applied to the rest of *S. aureus* and all of *E. coli* and *K. pneumoniae*.

| activator | paper | labeled | positives | out-of-fold ROC-AUC | PR-AUC | cross-assay ref |
|---|---|---|---|---|---|---|
| `adep4` | Conlon 2013 (ADEP4) | 1,677 | 230 (13.7%) | **0.8738 ± 0.0029** | **0.6103** | 0.877 |
| `onc212` | Jacques 2020 (ONC212) | 1,045 | 257 (24.6%) | **0.7671 ± 0.0060** | **0.5803** | 0.815 |

**Superseded numbers.** Until 2026-09-17 this stage shipped a RandomForest(500) scoring
0.8645 / 0.7517 (PR 0.5632 / 0.5481). `--estimator forest` reproduces those exactly; the retrofit
is recorded in the run log. **Read PR-AUC alongside ROC-AUC** — TabPFN's gain is ~4× larger there,
i.e. concentrated at the top of the ranking, which is what a shortlist consumes.

ROC-AUC is the mean ± SD over **5 repeats** of the cluster-grouped split. A single-seed estimate
carries ±0.004 of pure arbitrariness — larger than any of the three effects tested below — so one
seed is not reported.

**This overturns v1's headline negative result.** v1 pre-registered the gate *"does any sequence
feature beat protein length alone?"* and answered no — disorder added +0.009 over length, a
13-feature GBM lost to logistic length-only (0.762 vs 0.775), degron motifs took 0% of SHAP. Its
retrospective concluded *"a degradability axis needs new **data**, not new features."* Measured under
identical cluster-grouped CV, the ESM-C forest beats a properly-estimated length baseline by
**+0.081** (ADEP4, 0.857 vs 0.776) and **+0.062** (ONC212, 0.750 vs 0.687).

And on ADEP4 it essentially matches **what a second independent experiment achieves**: 0.865 against
0.877, where 0.877 is ONC212's own measured readout used to predict ADEP4's calls.

**On terminology.** That 0.877 was previously called a "ceiling" here. It is not one — the model
exceeds it at three of the five cutoffs in the sensitivity table below — and v1's
`degradability_datasets.md` §10.3 already uses "reproducibility ceiling" for a *different* quantity
(label-vs-label agreement, ρ 0.52 / Jaccard 0.32). It is now the **cross-assay reference**:
`CROSS_ASSAY_AUROC`, what the other screen achieves at the same task.

![ROC curves](../output/plots/degradability/roc_curves.png)

## Binary classification, and why the obvious objection does not apply

The two activators' binary calls overlap at only **Jaccard 0.32**, which looks like a reason not to
classify. It is not, and the distinction matters:

| | value | what it says |
|---|---|---|
| Jaccard of the two binary calls | 0.323 | the experiments disagree on **where to draw the line** |
| ADEP4 readout → ONC212 calls, AUROC | 0.815 | but they agree well on **which proteins rank as susceptible** |
| ONC212 readout → ADEP4 calls, AUROC | 0.877 | |

ONC212 simply drew a looser line — 165 hits ADEP4 misses against 28 the other way, it being the more
promiscuous compound at 30 µM. That is a *threshold* disagreement, and a classifier that emits a
probability and is scored on ranking never touches it.

Classification also makes the v1 comparison exact: 0.7747, 0.6799, 0.6916 and 0.7650 are all
classifier AUROCs. What it costs is magnitude inside each class — a −9.4 log2FC and a −1.1 both
become `1`. A regression variant of this stage was built and measured first; see the run log.

**Judge nothing against 1.0.** The cross-assay reference is 0.877 / 0.815, and it is drawn as its own
curve on the figure rather than as a horizontal line, because an AUROC is an area and not a
true-positive rate.

## The label — verified against v1's code

`legacy/scripts/10e_activator_features.py:119-131`, with `ABUNDANCE_LOG2 = -1.0`, `PADJ = 0.05`:

```python
out["conlon_abundance_bin"]  = (adep4_abundance_log2fc <= -1.0) & (adep4_abundance_padj < 0.05)
out["jacques_abundance_bin"] =  onc212_abundance_log2fc <= -1.0        # Jacques published no p-values
```

Unmeasured stays `np.nan`, never 0 — a missing readout is not a negative. Verified arithmetically:
233 + 1,479 = 1,712, exactly the proteins carrying an ADEP4 value.

The definitions are **asymmetric** — ADEP4 significance-filtered, ONC212 effect-size only — because
that is the best each paper supports. Do not "fix" it into symmetry or the numbers stop being
comparable to v1's.

The **cleavage** readout is deliberately unused: it reproduces between activators at only ρ = 0.22,
and within a single chemistry abundance and cleavage agree at **ρ = 0.06**. Not two views of one
quantity.

## The identifier problem — 0 of 1,943 join

Conlon is strain **COL** (`SACOL#####` / `YP_18xxxx`), Jacques is **C0673** (`ODV*`), the v2 proteome
is **NCTC 8325** (`SAOUHSC_*` / `YP_498xxx`).

| route | result |
|---|---|
| `saureus_acc` (`YP_*`/`ODV*`) → proteome `refseq` | **0 / 1,943** |
| `sacol` → `locus_tag_all` (`SAOUHSC_*`) | **0 / 1,943** |
| Jacques `Mu50_homolog` UniProt → `uniprot_ac` | **0 / 534** |
| `gene` → `gene_name` | 559 / 691 symbols, but only 692 of 1,943 rows carry one |

The obvious repairs are dead ends (`legacy/scripts/10c_clpp_activator.py:47-56`): modern RefSeq
re-tagged COL as `SACOL_RS#####` keeping **no** `old_locus_tag`, and UniProt demoted the COL and Mu50
proteomes as redundant — 939 and 1,033 entries survive — so a UniProt bridge would silently drop ~45%
of Conlon.

So the join is **by sequence**, CLAUDE.md's standing rule, at its `≥95% identity = direct` threshold:

```
1,943 labeled proteins -> 1,873 direct (96.4%) at >=95.0% id / >=80% cov
  direct  1873   below_threshold  54   no_hit  16
identity: median 100.0%   5th pct 99.3%
```

Loosening to ≥50% would gain only 21 proteins, so the strict cut costs nothing. Unlike v1's
cross-phylum hop at median 42% identity, this one is free. **The stage exits below 90% mapped.**

## Method

```python
RandomForestClassifier(n_estimators=500, class_weight="balanced", n_jobs=-1, random_state=0)
```

No `StandardScaler` — trees split on thresholds and are scale-invariant. `class_weight="balanced"`
carries the 13.7% / 24.6% positive rates without resampling.

**Cross-validation: two 5-fold stratified schemes side by side.**

- **`StratifiedGroupKFold(5, shuffle=True, random_state=0)`** grouped on the MMseqs2 30%-identity
  clusters (1,461 clusters over 1,677 ADEP4 proteins). Stratified keeps the base rate in every fold;
  grouped stops homologues straddling the split. **This is the reported number.**
- **`StratifiedKFold(5, shuffle=True, random_state=0)`**, plain, with `leakage_gap` emitted. Measured
  at **−0.010 to +0.002** — negligible, because the cluster count is close to the protein count.
  v1's `07d` used an ungrouped `cv=5`; this is what that costs, and here it costs nothing.
- **ROC-AUC** with a stratified bootstrap 95% CI, plus **PR-AUC** and the base rate beside it.

**Cross-activator**, the only non-circular test — `degradability_datasets.md` §10.4 listed it as
*"available now, not done"*:

| train | test | mode | n | ROC-AUC | 95% CI | v1 hand-built |
|---|---|---|---|---|---|---|
| adep4 | onc212 | shared | 1,045 | **0.747** | [0.710, 0.783] | 0.6916 |
| onc212 | adep4 | shared | 1,677 | **0.843** | [0.818, 0.865] | 0.7650 |
| onc212 | adep4 | cluster_disjoint | 580 | 0.836 | [0.797, 0.871] | — |

`shared` is v1's setup and the comparable one. `cluster_disjoint` also demands sequence novelty; the
**adep4 → onc212 direction cannot supply it** (ADEP4 covers almost every cluster, leaving 29 ONC212
proteins) and is reported as skipped rather than quietly dropped.

## Two things the estimator search settled

Both are cheap traps that would have gone unnoticed.

**1. A forest is the wrong estimator for a one-feature baseline.** On `log_length` alone, RF scored
0.682 (ADEP4) against balanced logistic's 0.776 on identical rows and splits, and saturated 33% of
its probabilities — 500 bootstrapped trees bin a single variable into steps instead of using it
monotonically. Quoting RF's length number would have credited the embedding with **+0.175** when the
honest figure is **+0.081**. The baseline ladder is retired from the pipeline, but that is why the
delta above is quoted against the *logistic* baseline.

**2. Logistic at `C=1.0` saturates and loses accuracy.** v1's `10f` hard-coded it, and 10f never ran,
so it was never caught. On 1,152 standardised dimensions at n≈1,677 the sigmoid pins **52% of
out-of-fold probabilities at exactly 0.000 or 1.000** — unrankable ties across half the proteome —
while also costing real accuracy:

| C | ROC-AUC (ADEP4) | saturated |
|---|---|---|
| 1.0 | 0.813 | 52.3% |
| 0.01 | 0.855 | 1.4% |
| 0.001 | 0.873 | 0.0% |

A forest has no such cliff and saturates 0% untuned, which is what broke the tie: on ESM-C the two
estimators are otherwise a wash (logistic 0.873 vs RF 0.857 on ADEP4; 0.744 vs 0.750 on ONC212, CIs
overlapping both ways). `HistGradientBoostingClassifier` was also measured and rejected — no better,
and it saturated 25–67%. **A `saturated > 20%` guard now exits the stage**, so this cannot regress
silently.

## Three questions asked, three negatives

All three were measured under this stage's own cluster-grouped CV. None improves the model, and that
is the useful result: it locates the remaining headroom in the labels rather than the method.

### 1. Localization and function features add nothing

| block | features | ADEP4 | ONC212 |
|---|---|---|---|
| **esmc** (shipped) | 1152 | **0.857** [0.827, 0.882] | **0.750** [0.716, 0.784] |
| localization alone | 11 | 0.748 | 0.646 |
| function (COG) alone | 24 | 0.685 | 0.612 |
| esmc + localization | 1163 | 0.859 (**+0.002**) | 0.745 (**−0.005**) |
| esmc + loc + function | 1187 | 0.860 (+0.003) | 0.750 (+0.000) |

*(single-seed, seed 0, which is why `esmc` reads 0.857 here against the 0.865 headline.)*

Localization *alone* is a decent predictor, but the embedding already encodes it — DeepLocPro is
ESM-2-based and TMbed is ProtT5-based, so this is redundancy, not new information. v1's own lesson,
recurring: *"SHAP splits credit among correlated features … the nested models measure what is
added."* **The model stays ESM-C-only**, with no dependency on stages 02/03.

**But the compartment breakdown is worth more than the AUROC**, and it does not gate the way
intuition says:

| class | ADEP4 hit rate | ONC212 hit rate |
|---|---|---|
| cytoplasm | 0.180 | 0.275 |
| extracellular | 0.049 | **0.346** |
| **cytoplasmic_membrane** | **0.031** | **0.105** |
| cell_wall_surface | 0.000 (n=15) | 0.000 (n=7) |

**Membrane proteins are protected in both screens** — co-translational insertion means they never
present a soluble cytoplasmic chain for an open ClpP pore to take. **Secreted proteins are not
protected**, because they transit the cytoplasm as unfolded chains, which is exactly what activated
ClpP prefers. So if you filter on localization downstream, **"cytoplasmic only" is the wrong rule** —
"not membrane" is closer to what the data says. The two screens disagree on `extracellular`
(0.049 vs 0.346, n=82/52), so treat that class as unsettled. Join via `src/localization.py`.

### 2. Hyperparameters do not matter

A 10-point sweep over `max_features` ∈ {sqrt, log2, 0.02, 0.05, 0.10}, `min_samples_leaf` ∈ {1, 3}
and `criterion` ∈ {gini, entropy} spanned **0.854–0.862** (ADEP4) and **0.742–0.750** (ONC212) — less
than one CI half-width, with the winner differing between activators. Reseeding the CV showed it was
selection noise:

| config | ADEP4, 5 CV seeds | ONC212, 5 CV seeds |
|---|---|---|
| default (`sqrt`, msl=1) | 0.8621 ± 0.0036 | 0.7483 ± 0.0052 |
| sweep winner (`0.1`, msl=3) | 0.8645 ± 0.0021 | 0.7517 ± 0.0075 |

**+0.002 / +0.003, at or below seed noise.** `RF_PARAMS` adopts the winner because it is free and
consistent in sign — not because tuning mattered. Do not spend more time here.

### 3. The log2FC cutoff must not be tuned — but sweeping it proves the result is cutoff-free

**AUROC is not comparable across different labels.** Tightening the cutoff shrinks the positive class
and makes it more extreme, which inflates AUROC mechanically. Optimising on that would buy a spurious
+0.08 and a label chosen for being easy rather than meaningful.

`accessory/cutoff_sensitivity.tsv`, recomputed every run:

| activator | cutoff | n_pos | model | cross-assay | **gap** | length | ESM-C − length |
|---|---|---|---|---|---|---|---|
| adep4 | −0.5 | 379 | 0.827 | 0.804 | **+0.024** | 0.729 | +0.099 |
| adep4 | **−1.0** | 234 | 0.862 | 0.870 | **−0.008** | 0.780 | +0.082 |
| adep4 | −1.5 | 136 | 0.892 | 0.873 | **+0.020** | 0.843 | +0.049 |
| adep4 | −2.0 | 96 | 0.904 | 0.875 | **+0.028** | 0.843 | +0.060 |
| adep4 | −3.0 | 46 | 0.943 | 0.956 | **−0.014** | 0.871 | +0.072 |
| onc212 | −0.5 | 351 | 0.718 | 0.778 | −0.060 | 0.661 | +0.057 |
| onc212 | **−1.0** | 257 | 0.746 | 0.816 | −0.070 | 0.687 | +0.059 |
| onc212 | −1.5 | 191 | 0.764 | 0.846 | −0.082 | 0.717 | +0.048 |
| onc212 | −2.0 | 151 | 0.779 | 0.860 | −0.081 | 0.739 | +0.040 |
| onc212 | −3.0 | 96 | 0.796 | 0.873 | −0.077 | 0.786 | **+0.009** |

**`model` rises with strictness for a trivial reason and is not a tuning signal.** The column that
matters is `gap`: the cross-assay reference rises in lockstep, so their difference is the invariant —
ADEP4 within ±0.03 at **every** cutoff, ONC212 a consistent ~−0.07 below. That makes the shipped
conclusion **cutoff-independent**, which is far stronger than a single-cutoff claim.

One more thing the sweep shows: **ESM-C's margin over length shrinks as the cutoff tightens** —
ONC212 falls from +0.057 at −0.5 to **+0.009** at −3.0. The most dramatically depleted proteins are
largely just the smallest ones, so the embedding earns its keep in the moderate-effect middle rather
than at the extremes.

**The cutoff stays at −1.0**: it is v1's audited definition (keeping 0.7747 comparable), it is what
both papers effectively used, and it retains 234/257 positives where −3.0 leaves 46. A guard asserts
that the sweep's −1.0 row agrees with the main CV number, so the two code paths cannot drift apart
about what the label is.

## Outputs

```
data/processed/degradability/
    degradability_<species>.tsv               8 columns, keyed on uniprot_ac, every protein present
    accessory/
        labels_saureus.tsv                    measured call + log2FC on uniprot_ac + join evidence
        seqmap_audit.tsv                      every join: source acc, target, pident, cov, verdict
        cv_<activator>.tsv                    ROC-AUC, PR-AUC, base rate, leakage_gap, saturation
        cross_activator.tsv                   both directions, both modes
        cutoff_sensitivity.tsv                model vs cross-assay vs length, per cutoff
        domain_bands.tsv                      AUROC by similarity band + per-species reweighting
        oof_<activator>_saureus.tsv           out-of-fold probability + fold index per labeled protein
        model_<activator>.joblib              fitted forest + its CV metrics
        manifest.tsv
output/plots/degradability/roc_curves.png
```

| column | meaning |
|---|---|
| `<act>_hit` | the **measured** call, `1`/`0`, empty where that screen did not measure this protein |
| `<act>_prob` | the model's probability, for **every** protein — **out-of-fold** where labeled |
| `<act>_source` | `measured` \| `predicted` |
| `nn_similarity` | cosine to the nearest *S. aureus* training protein in ESM-C space |

**The three-column split is what keeps `_prob` honest.** It is one comparable scale across all 13,020
proteins, because a labeled protein carries its out-of-fold value rather than a 1.0 that would
outrank every uncertain prediction. Ground truth stays in `_hit`.

### Thresholding: 0.5 is the wrong cut

A balanced forest at a 13.7% base rate stays conservative — `adep4_prob` tops out at 0.614 on
*K. pneumoniae*, so `>= 0.5` selects **20 of 5,728**. The thresholds that reproduce each activator's
own training base rate are **0.328** (adep4) and **0.313** (onc212) — re-derived for TabPFN; the forest's were 0.394 / 0.426, which is what
`src.degradability.hits()` defaults to (320 Kp / 311 Ec / 359 Sa for adep4). Better still: rank on
`_prob` and take a top-N. The model was validated on ranking, not on any cut.

| helper | returns |
|---|---|
| `load(species)` / `load_all()` | the shipped table |
| `hits(df, activator, threshold=None)` | predicted-or-measured positives, base-rate threshold by default |
| `load_labels()` | measured calls + log2FCs on `uniprot_ac` |
| `load_seqmap_audit()` | every sequence join, including the 70 that failed |
| `load_cv(activator)` / `load_cross_activator()` / `load_domain_bands()` | the validation tables |
| `load_cutoff_sensitivity()` | the cutoff sweep — read the `gap` column, not `model` |
| `load_oof(activator)` | out-of-fold probabilities + fold index |
| `manifest()` | per-species counts and provenance |

## Pricing the cross-species leap

There are **no** activated-ClpP measurements for *E. coli* or *K. pneumoniae*, anywhere, so nothing
validates those rows. The stage bands out-of-fold AUROC by cosine distance to the nearest
*S. aureus* training protein (self excluded) and reweights by where each species' queries land — the
control that turned a rejected k-NN's 86.8% headline into an honest 68% in stage 02.

Reweighted expectation: **Kp 0.752 / Ec 0.754** (ADEP4), **Kp 0.659 / Ec 0.660** (ONC212). Median
`nn_similarity` on predicted rows is 0.928 (Kp), 0.936 (Ec).

**Read this in the direction that does not flatter it.** Only two bands populate — nothing lands
below cosine 0.90 — so the banding has almost no dynamic range to detect a decay with distance. And
the premise that ESM-C cosine measures transferability is **unvalidated and untestable** without
Gram-negative labels.

## Figures

Two scripts, answering different questions. Both write `output/plots/degradability/`.

| script | figure | question |
|---|---|---|
| `plots/degradability.py` | `roc_curves.png` | **how good is the model** — per-fold ROC on the one species with labels, against the cross-assay reference |
| `plots/degradability_predictions.py` | `predictions.png` | **what did it predict** — score distributions per species and activator, plus *S. aureus* split by measured status |
| | `agreement.png` | **are the two activators independent evidence?** No |
| | `extrapolation.png` | **what are the Kp/Ec numbers worth?** |
| `plots/degradability_top.py` | `top_adep4.png`, `top_onc212.png` | **which proteins are at the top**, named |
| | `top_composition.png` | what the top of the list is made of, and how much the two shortlists share |
| `degradability/enrichment.py` | `cog_fisher.png` | **Fisher-exact enrichment of COG categories** at the top of the ranking |
| | `loc_fisher.png` | **Fisher-exact enrichment of the localization compartments** and the TMbed topology features |
| | `interest_panel.png` | where the **GraDi consortium's proteins of interest** actually rank |
| | `enrichment_fisher.tsv` | every test, with raw 2x2 counts (`output/results/degradability/`) |

### `predictions.png`

Panels A/C are the predicted-probability distribution per species against the base-rate cut; B/D are
*S. aureus* split into measured hit / measured non-hit / unmeasured. B and D are the only panels in
the axis where a predicted score sits on the same axis as a real one, and the separation is visible:
ADEP4 measured hits centre near 0.45 against 0.20 for measured non-hits, while ONC212's two groups
(0.44 vs 0.32) sit much closer — which is the AUROC difference (0.865 vs 0.750) seen as a
distribution rather than a curve.

At the base-rate cut the selected fraction is **13.1% Kp · 14.5% Ec · 15.4% Sa** for ADEP4 and
**18.7 / 18.7 / 25.5%** for ONC212. Those are close to each activator's training base rate by
construction — that is what the cut is for — and are *not* evidence about how many real substrates
each proteome holds.

Note how far the distributions sit from 1.0: the ADEP4 median is 0.22–0.26 and almost nothing exceeds
0.7. This is the same fact as *0.5 is the wrong threshold*, seen directly.

### `agreement.png`

adep4 vs onc212 predicted probability, one panel per species, hexbinned. A tight near-linear ridge:
**rho +0.835 Kp · +0.833 Ec · +0.835 Sa** over all rows. The underlying *labels* agree at only
rho 0.52 / Jaccard 0.32, so the models agree far more than the experiments do — both read the same
ESM-C embedding, so the shared signal is the learnable part, not corroboration. **A shortlist built
on "both activators agree" is close to one opinion, not two.**

One number to keep straight: the stage itself prints rho over the **predicted rows only**, which for
*S. aureus* is a different population (n=1,212, rho **+0.866**) from the figure's all-rows
n=2,889 / +0.835. Kp and Ec have no measured rows, so both routes agree exactly there. The figure
puts `n` in each title for this reason.

### `extrapolation.png`

Panel A is `nn_similarity` per species with the band edges drawn: medians **0.928 Kp · 0.936 Ec ·
1.000 Sa** (Sa's is 1.000 because the measured proteins are their own nearest neighbours; the
manifest's 0.973 is the predicted-only median). Kp and Ec are genuinely further from the training set
than Sa, and Kp furthest.

Panel B is the honest part, and it is the one worth looking at before quoting any expected AUROC:
**a fifth of *K. pneumoniae* is unpriced.** Only two of four bands carry a measured AUROC, so the
reweighted figure covers **79.3% of Kp · 87.5% of Ec · 92.1% of Sa** (ADEP4) and says nothing at all
about the rest. The `expected_roc_auc` values — ADEP4 0.780 Kp / 0.785 Ec, ONC212 0.772 / 0.770 — are
conditional on being in a priced band.

Panel C is why even the priced part is shaky. ADEP4 rises with proximity as one would hope
(0.757 [0.695, 0.823] → 0.847 [0.812, 0.883]), but **ONC212 falls** (0.784 [0.699, 0.861] → 0.732
[0.694, 0.775]) — the opposite of the premise the reweighting rests on. The confidence intervals
overlap, so the non-monotonicity is within noise; that is the charitable reading, and it is also an
admission that with two bands and 118 proteins in the smaller one there is not enough here to
establish the trend either way. **The premise that ESM-C cosine measures transferability remains
unvalidated, and is untestable without Gram-negative labels.**

### `top_adep4.png` / `top_onc212.png`

The 20 highest-ranked proteins per species, labelled and coloured by COG group. Labels resolve
`gene_name` -> eggNOG `preferred_name` -> accession, which names 20/20 for Kp and Ec and 16/20 for
Sa: the top of the ranking is unusually well annotated because high-scoring proteins are small,
conserved and abundant. **Ranking, not thresholding** — the model was validated on AUROC.

**The two activators rank different biology, and it is obvious on sight.** ONC212's top is almost
pure ribosome — `rplX rplR rplC rpsQ rpsP rplA rplW rpsH`, one COG group filling the panel. ADEP4's
is broader: ribosomal proteins and translation factors (`rplK rpsA infA infB`) but also H-NS
chromatin proteins (`hns stpA`), cold-shock proteins (`cspB cspC cspI cspG`), FKBP-type PPIases
(`fklB slyD`), glutaredoxins (`grxA grxC`) and chaperones (`grpE groES`). That difference is the
mechanism behind the low shortlist overlap below, and it is the strongest argument for having kept
the two activators as separate models.

A `*` marks a protein the *S. aureus* screens actually measured as a hit. **This is the closest
thing to a validation of the top of the list, and it holds**: of Sa's top 20, ADEP4 has 13 measured
of which **9 are confirmed hits**, and ONC212 has 16 measured of which **14 are confirmed**. Kp and
Ec carry no `*` at all, because they have no measured labels — every bar there is a hypothesis.

### `top_composition.png`

Panels A and B put the top 100's composition against the whole proteome's (the black tick), so
enrichment reads as bar-vs-tick.

- **The top is overwhelmingly cytoplasmic: 94% Kp / 95% Ec / 97% Sa against ~60% of the proteome.**
  For an axis whose whole premise is that cytoplasmic Clp has to reach the substrate, this is the
  reassuring direction, and it was not built in — localization is not a feature of the model
  (`docs` §"Three questions asked, three negatives" measured that adding it changes nothing).
- **Information storage and processing is 36–44% of the top 100 against ~16% of the proteome** —
  translation and transcription machinery, which is what both screens found.

Panel C is the one that changes how the ranking should be used. Despite rho 0.835 between the two
activators' scores, their **top 100 lists share only 24/100 on Kp**, 32 on Ec, 40 on Sa — and at
top-10 the Kp lists share *nothing*. Chance is ~2/100, so the agreement is real and far above
random; it is simply nowhere near interchangeable. **Global rank correlation is not shortlist
agreement**, and a shortlist is what gets taken. Read alongside `agreement.png`: the two activators
are close to one opinion in aggregate and clearly two at the sharp end, and both statements matter.

### A localization mis-call the top list exposed

`rplK` — ribosomal protein L11, the **highest-ranked ADEP4 protein in both Gram-negatives** — is
called `outer_membrane` by DeepLocPro in Kp *and* Ec, while TMbed puts it at
`cytoplasmic_fraction = 1.0` with no signal peptide. A ribosomal protein is cytoplasmic; the
compartment call is wrong, reproducibly, in both species (Sa gets it right). Nothing downstream
breaks, but it is a concrete instance of the advice in `docs/localization.md`: the compartment
column is good in aggregate and fallible per protein, so check `cytoplasmic_fraction` before
trusting an individual call — especially for a protein you are about to put in a shortlist.

### `cog_fisher.png` — what the top of the ranking is made of, tested

**"The top" is the top 10% of each proteome** — 573 Kp / 440 Ec / 289 Sa. Not a fixed count: the
proteomes differ by 2x in size, so one absolute N would mean the top 1.7% of Kp against the top 3.5%
of Sa and the odds ratios would not be comparable between species. `--top-pct` changes the fraction,
`--top-n` forces an absolute count.

One 2x2 Fisher exact test per COG category over the **full proteome**: in-top-10% vs not,
in-category vs not. Two-sided, Benjamini-Hochberg **within each species x activator x test-family
block**. Category descriptions come from COGclassifier's own bundled
`resources/cog_func_category.tsv`, never a hard-coded list. All 194 tests with raw 2x2 counts in
`output/results/degradability/enrichment_fisher.tsv`; 67 are significant at FDR 0.05.

Odds ratios, `*` = FDR < 0.05:

| | category | Kp | Ec | Sa | Kp | Ec | Sa |
|---|---|---|---|---|---|---|---|
| | | **ADEP4** | | | **ONC212** | | |
| **J** | Translation, ribosomal structure | 4.24\* | 4.30\* | 3.08\* | 5.81\* | 7.64\* | 10.17\* |
| **O** | Protein turnover, chaperones | 2.79\* | 2.69\* | 2.16\* | 1.83\* | 2.06\* | 2.16\* |
| **B** | Chromatin structure and dynamics | 63.7\* | 64.1\* | 9.02 | 63.7\* | 27.4\* | 9.02 |
| **K** | Transcription | 1.94\* | 2.46\* | 2.77\* | 0.69 | 1.25 | 1.85\* |
| **X** | Mobilome: prophages, transposons | 2.06\* | 7.77\* | 2.24 | 0.43 | 0.90 | 0.00\* |
| **S** | Function unknown | 2.51\* | 2.00\* | 1.87 | 1.24 | 1.43 | 1.28 |
| **M** | Cell wall/membrane/envelope biogenesis | 0.20\* | 0.12\* | 0.06\* | 0.06\* | 0.12\* | 0.12\* |
| **E** | Amino acid transport and metabolism | 0.14\* | 0.11\* | 0.12\* | 0.23\* | 0.18\* | 0.12\* |
| **G** | Carbohydrate transport and metabolism | 0.17\* | 0.16\* | 0.30\* | 0.23\* | 0.18\* | 0.30\* |
| **C** | Energy production and conversion | 0.19\* | 0.15\* | 0.48 | 0.29\* | 0.28\* | 0.57 |
| **I** | Lipid transport and metabolism | 0.11\* | 0.14\* | 0.53 | 0.22\* | 0.28\* | 0.65 |
| **N** | Cell motility | 0.00\* | 0.08\* | 0.00 | 0.00\* | 0.08\* | 0.00 |

**A robust core, significant in all six species x activator combinations.** Up: **J** translation
and ribosome, **O** protein turnover and chaperones. Down: **M** envelope biogenesis, **E** amino
acid metabolism, **G** carbohydrate metabolism. Five categories, the same direction and roughly the
same magnitude in three organisms and two chemistries — the strongest internal evidence the axis has
that the model learned something real. **Neither function nor localization is a model feature** (the
negatives above measured that adding them changes nothing), so this is the ESM-C embedding recovering
the screens' biology unaided.

**B is the largest effect and the least trustworthy.** OR 64 on both Gram-negatives, with 7 of 8
chromatin proteins in the top 10% — H-NS, StpA, the HU subunits. The direction is almost certainly
real (those proteins are small, abundant, nucleoid-associated and appear by name at the top of the
ADEP4 lists), but **the category has 8 members**, so the odds ratio itself carries no useful
precision. Quote the count, not the ratio.

**The two activators differ, and it is the same difference the named lists show.** ONC212's profile
is narrow: J at OR 5.8–10.2 and almost nothing else moving. ADEP4's is broad — J plus **K**
transcription (significant in all three), **X** mobilome, **S** function unknown, **V** defense,
**D** cell division. **X flips sign between activators** (2.1–7.8 up for ADEP4, 0.0–0.9 for ONC212,
significantly *depleted* in Sa). This is `top_adep4.png` vs `top_onc212.png` — pure ribosome versus
ribosome plus chromatin, cold-shock and chaperones — restated as a statistic.

**`unclassified` is tested, not dropped**, and it does not behave consistently: enriched on Kp
(1.61\* ADEP4, 2.37\* ONC212) and on Ec ONC212 (1.70\*), but significantly *depleted* on Sa ONC212
(0.50\*). Dropping these proteins would have changed the background for every other test, so they
stay in; the inconsistency is a caution against reading the unclassified fraction as a single
phenomenon.

Two smaller decisions: `cog_category` is the **first letter** of a possibly multi-letter assignment
(9–16% of classified proteins carry more than one — `docs/function.md`), and `--multi` re-runs
counting a protein in every letter it carries. And the odds ratio is **reported raw but plotted with
a Haldane-Anscombe (+0.5) correction**, because a zero cell cannot be drawn on a log axis; the TSV
keeps the uncorrected value and the counts.

### `loc_fisher.png` — the top of the ranking is cytoplasmic, and Clp can reach it

The same 2x2 Fisher test, same top 10%, applied to the stage-03 compartments plus the TMbed features
that Clp reachability actually turns on. Ten groups per species: the six DeepLocPro compartments
(which partition the proteome exactly, since DeepLocPro always calls), plus `has_signal_peptide`,
`n_tm_helix > 0`, `n_tm_strand >= 8` (a beta-barrel) and `cytoplasmic_fraction == 1`.

| group | ADEP4 | ONC212 | reading |
|---|---|---|---|
| **wholly cytoplasmic** (`cytoplasmic_fraction == 1`) | **OR 11.5–13.1\*** | 11.5–13.1\* | the single strongest enrichment in the axis |
| **Cyt** (cytoplasm) | OR 7.9–10.3\* | 7.9–10.3\* | |
| **OM** (outer membrane) | 0.54\* Ec, 0.70 Kp | 0.99–1.02, q ≈ 1 | ADEP4 mildly avoids it; ONC212 is *exactly at chance* |
| **Ext** (extracellular) | 0.10–0.38\* | 0.62–0.83, n.s. | significant for ADEP4 only |
| **CM** (cytoplasmic membrane) | 0.08–0.09\* | 0.05–0.13\* | |
| **Peri** (periplasm) | 0.08–0.10\* | 0.07–0.10\* | |
| **signal peptide** | 0.04–0.11\* | 0.01–0.07\* | an export signal is close to disqualifying |
| **TM helix** | 0.05–0.07\* | 0.02–0.07\* | |
| **beta-barrel** | **0.00\*** | **0.00\*** | **not one of the 67 Kp / 66 Ec barrels is in the top 10%** |
| **CW** (cell wall surface) | 0.50, n.s. | 0.00, n.s. | Sa only, n=19 — no power |

**This is the axis's premise confirmed from the outside.** Cytoplasmic Clp should reach cytoplasmic
proteins and nothing else, and the ranking behaves exactly that way — without localization ever being
a model feature.

**TMbed's residue count beats DeepLocPro's class label at predicting a high score.**
`cytoplasmic_fraction == 1` gives OR 11.5–13.1 against 7.9–10.3 for the `cytoplasm` compartment
itself. The two predictors share no machinery, so this is a genuine comparison: *how much of the
protein sits in the cytoplasm* separates degradable from not slightly better than *which compartment
it was assigned to*. Consistent with `docs/localization.md`'s advice to prefer the fraction for
individual proteins.

**Zero beta-barrels in the top 10%**, on both Gram-negatives. A clean, mechanistically sensible
absolute: a 16–26-strand barrel folded into the outer membrane is not something a cytoplasmic
protease processes.

**Where the two activators differ is the outer membrane.** For ADEP4, `OM` and `Ext` are depleted;
for ONC212 they sit at chance (OR 0.99–1.02, q ≈ 1). Read cautiously rather than mechanistically:
`docs/localization.md` records that `extracellular` is the least trustworthy compartment (TMbed
disputes 58% of those calls), so a null there may be telling us about the label rather than the
chemistry.

Compartments absent from a species — `CW` in a Gram-negative, `Peri`/`OM` in a Gram-positive, and
beta-barrels in Sa — are skipped rather than tested, which is why those points are missing from the
figure. That is a structural zero, not a failed test.

### `interest_panel.png` — where the consortium panel sits

Secondary to the COG analysis above, and reported for completeness. The GraDi targets are recorded
only in prose (`legacy/docs/03_degradability.md:420`, the ~30 envelope targets; `legacy/HISTORY.md:611`,
the v5 GyrA/GyrB proposal); `src/interest.py` writes them down as data and is **explicitly a draft**
— see its `EXPANSION_NOTES`. Coverage by gene symbol: Ec 43/43, Kp 41/43, Sa 11/43, where Sa's gap is
mostly the *structural* absence of LPS/Lpt/Bam in a Gram-positive.

At the top 10% the panel is **consistently below expectation but nowhere near significant**: 3/44,
1/44, 1/43, 2/43, 1/12, 0/12 against 1.2–4.4 expected, odds ratios 0.21–0.82, every q > 0.2. So the
honest statement is a *trend* toward depletion with no power to call it — not the "OR 0.00, q = 1.0"
that a fixed top-100 produced by having almost no panel proteins in scope at all. It agrees in
direction with **M** (envelope biogenesis) being the most strongly depleted COG category, and with
the warning already at `legacy/docs/03_degradability.md:420` that these targets are almost all
envelope and reachable only through the pre-export window.

The members that do rank high are the cytoplasm-facing ones — `secA` 96th percentile, `secB` 94–95th,
`yajC` 89–93rd, `lptB` (the Lpt complex's cytoplasmic ABC ATPase) 79–81st, GyrA/GyrB 60–70th on Sa.

## Run log

### 2026-09-02 — Random Forest classifier, shipped

DIAMOND 2.2.1 · scikit-learn 1.9.0 · 5 folds · seed 0 · **30 s** wall clock, plots ~10 s.
ADEP4 0.857 [0.829, 0.882], PR-AUC 0.535 at a 0.137 base rate. ONC212 0.750 [0.714, 0.783], PR-AUC
0.539 at 0.246. Leakage gap +0.002 / −0.010. Saturation 0.0% on both.

**Spot checks — v1's pre-registered predictions, all held** (`HISTORY.md` §7.5). In binary form these
test the threshold as well as the value, which the continuous version could not:

| gene | ONC212 | ADEP4 | v1's prediction |
|---|---|---|---|
| `dnaK` | hit **0** (log2FC −0.51) | hit 1 | "1.4× down, NOT ≥2×" ✅ |
| `acpP` | hit **0** (log2FC +0.03) | hit 1 | "unchanged (+0.03)" ✅ exact |
| `gyrB` | unmeasured | hit 0 | "essentially no evidence" ✅ |

`acpP` matching to the second decimal is a hard check on the whole chain: xlsx parse → v1's merge →
DIAMOND sequence join → label table.

**The shared-confound diagnostic, which v1 had no way to run — and it is the least comfortable
number here.** The two activators' predicted probabilities correlate at **ρ = 0.824 (Kp) / 0.837 (Ec)
/ 0.844 (Sa)**, while their *labels* agree at only ρ 0.52 continuous / Jaccard 0.32. Two models
trained on moderately-agreeing targets are producing predictions that agree far more than the
experiments do.

Some of that is expected and benign: a model output is smoother than a noisy measurement, and both
models read the same 1,152-dim embedding, so whatever signal is learnable from it is shared. But it
means **`adep4_prob` and `onc212_prob` are not two independent lines of evidence** — treating
agreement between them as corroboration would be overreading. The measured AUROCs are unaffected
(they are scored against held-out labels), but a shortlist built on "both activators agree" is closer
to one opinion than two.

The flag fires above 0.90, so this does not stop the run. For comparison the regression variant
measured 0.72–0.80 on the same diagnostic, i.e. slightly *better* separated than the forest.

### 2026-09-02 — prediction figures, and what they made explicit

`scripts/plots/degradability_predictions.py` added: three figures covering all three species.
Nothing refitted — the tables and `accessory/domain_bands.tsv` are read as written.

The stage's own outputs had to be regenerated first: `degradability_<species>.tsv` were absent while
`accessory/` (models included) was intact, so a plain re-run rebuilt them deterministically in ~30 s
and reproduced every recorded number — n 5,728 / 4,403 / 2,889, `nn_similarity` medians, the spot
checks, and the band AUROCs all identical to the manifest.

Three things the figures turned from prose into numbers:

- **A fifth of *K. pneumoniae* is unpriced.** "Only two bands are populated" was already written
  down, but not what it costs: the reweighted expected AUROC covers **79.3% of Kp**, 87.5% of Ec,
  92.1% of Sa, and is silent on the rest. Any quote of "expected AUROC 0.78 on Kp" is conditional on
  a protein being in a priced band, and one in five is not.
- **The non-monotonicity is ONC212's alone, and it is within noise.** ADEP4 rises with proximity
  (0.757 → 0.847) as the premise predicts; ONC212 falls (0.784 → 0.732). Plotting the bootstrap CIs
  shows they overlap, so the reversal is not established — but neither is the trend, on 118 proteins
  in the smaller band.
- **The two activators' predictions agree more than their experiments do**, visibly: a single tight
  ridge at rho 0.835 against labels that agree at 0.52. The figure is the clearest statement of the
  "not independent evidence" caveat in the axis.

One bookkeeping trap found while writing it: the stage prints rho over **predicted rows only**, so
its *S. aureus* figure (n=1,212, +0.866) is not the same quantity as the all-rows +0.835 the plot
shows. Both are right; they are different populations, and Kp/Ec cannot show the discrepancy because
every row there is predicted. The figure carries `n` in each panel title so the two cannot be
conflated.

### 2026-09-02 — the named top of the list

`scripts/plots/degradability_top.py` added, on request to see the actual proteins rather than
distributions. `--top` and `--composition-top` set the two list lengths (20 and 100 by default).

**The ranking survives being read.** Both activators put small, abundant, cytoplasmic
translation-machinery proteins at the top, which is what the two screens measured and what the
length correlation predicts. Where truth exists it is largely right: **9 of the 13 measured proteins
in Sa's ADEP4 top 20 are confirmed hits, and 14 of 16 for ONC212.**

**The two activators rank visibly different biology.** ONC212's top 20 is almost pure ribosome;
ADEP4's mixes ribosome with H-NS chromatin proteins, cold-shock proteins, FKBP PPIases and
chaperones. This is the qualitative face of a number already recorded — and it retires any
temptation to merge them.

**The finding that changes how to use the output**: rho 0.835 between the activators' scores, but
their **top-100 lists share only 24/100 on Kp** (32 Ec, 40 Sa), and their top-10 on Kp share
nothing. Global rank correlation and shortlist agreement are different quantities, and the axis had
only measured the first. Both readings are needed: aggregate near-duplication argues against
treating "both activators agree" as corroboration, while the shortlist divergence argues against
picking either list as *the* answer.

**A cross-stage error surfaced.** `rplK` (ribosomal L11) is the top-ranked ADEP4 protein in both
Gram-negatives and DeepLocPro calls it `outer_membrane` in both, while TMbed reports
`cytoplasmic_fraction = 1.0` and no signal peptide. The compartment call is simply wrong, and
reproducibly so (Sa is correct). Recorded under Figures; it is a live example of why
`docs/localization.md` says to check the fraction before trusting an individual compartment call.

No plot-layout lesson worth generalising, except one that keeps recurring: in a panel whose bars
reach 94% there is no in-axes space for a legend, so legends here go **below** their own panel.
Placing them under each axes rather than under the figure keeps two neighbouring legends from
colliding.

### 2026-09-02 — localization enrichment

Same Fisher machinery, third test family: the six stage-03 compartments plus four TMbed-derived
features. `loc_fisher.png`.

**The axis's premise, confirmed from outside the model.** `cytoplasmic_fraction == 1` is enriched at
**OR 11.5–13.1** and the `cytoplasm` compartment at 7.9–10.3, while signal peptides (0.01–0.11), TM
helices (0.02–0.07), the cytoplasmic membrane (0.05–0.13) and the periplasm (0.07–0.10) are all
strongly depleted. **Not one of the 67 Kp / 66 Ec beta-barrels is in the top 10%.** Localization is
not a model feature, so the ranking arrived here on its own.

**A cross-predictor result worth keeping**: TMbed's residue count out-predicts DeepLocPro's class
label. `cytoplasmic_fraction == 1` (OR 11.5–13.1) beats the `cytoplasm` call itself (7.9–10.3) at
identifying high-scoring proteins. The two share no machinery, so this is a fair comparison, and it
supports `docs/localization.md`'s advice to prefer the fraction for individual proteins.

**One activator difference, read cautiously.** ADEP4 depletes `OM` and `Ext`; ONC212 sits at chance
there (OR 0.99–1.02, q ≈ 1). Since `extracellular` is the least trustworthy compartment in stage 03
— TMbed disputes 58% of those calls — this may be a fact about the label rather than the chemistry.

### 2026-09-02 — Fisher enrichment of the COG categories

`scripts/degradability/enrichment.py` added. Fisher exact per COG category over the full
proteome, **top 10% as hits** (573 Kp / 440 Ec / 289 Sa), two-sided, BH within each
species x activator x test-family block. 194 tests, 67 significant at FDR 0.05.

**Five categories move the same way in all six species x activator combinations.** Enriched: **J**
translation and ribosome (OR 3.1–10.2), **O** protein turnover and chaperones (1.8–2.8). Depleted:
**M** envelope biogenesis (0.06–0.20), **E** amino acid metabolism (0.11–0.23), **G** carbohydrate
metabolism (0.16–0.30). Three organisms, two chemistries, same direction and similar magnitude.
Since neither function nor localization is a model feature — measured earlier as adding nothing —
this is the embedding recovering the screens' biology unaided, and it is the best sanity check the
axis has.

**The two activators' profiles differ exactly as their named top-20 lists do.** ONC212 is narrow: J
at 5.8–10.2 and little else. ADEP4 is broad: J plus K transcription (significant in all three
species), X mobilome, S function unknown, V defense, D cell division. **X flips sign** between the
two (2.1–7.8 up for ADEP4; 0.0–0.9 for ONC212, significantly depleted in Sa).

**B (chromatin) is the biggest number and the weakest evidence.** OR 64 on both Gram-negatives with
7 of 8 members in the top 10%. The direction is credible — H-NS/StpA/HU appear by name at the top of
the ADEP4 lists — but with 8 members the ratio has no precision. Quote the count.

Two decisions worth recording:

- **Top 10% rather than a fixed count.** The proteomes differ 2x in size, so an absolute N defines
  "the top" as 1.7% of Kp against 3.5% of Sa and the odds ratios stop being comparable across
  species. This also gave the test real power: at a fixed top-100 most small categories had zero
  members in scope and returned OR 0.00 with q = 1.0, which looks like a finding and is not one.
- **BH within each test family, not across the block.** The COG categories and the consortium panel
  are different questions asked of the same ranking; pooling them made a COG category's q-value
  depend on how many panel families happened to be tested beside it.

`unclassified` is tested rather than dropped (removing it would change the background for every other
test) and does *not* behave consistently — enriched on Kp, depleted on Sa/ONC212. Do not read the
unclassified fraction as one phenomenon.

The consortium panel is reported alongside but is secondary: at the top 10% it is consistently below
expectation (OR 0.21–0.82) with every q > 0.2, i.e. a trend toward depletion with no power to call
it. `src/interest.py` had to be written to test it at all, since the targets existed only as prose;
it is a draft, with every expansion judgement recorded in its `EXPANSION_NOTES`.

### Superseded: the regression variant (2026-09-01)

This stage was first built as a **regression** on the continuous abundance log2FC, on the argument
that thresholding costs reproducibility (ρ 0.52 → Jaccard 0.31). **That argument was wrong** — it
measured set overlap, whereas a classifier is scored on ranking, where the two activators agree at
AUROC 0.82–0.88. Kept here because the numbers are informative:

| | length | ESM-C ridge | v1's label-agreement ρ |
|---|---|---|---|
| ADEP4 | ρ 0.319, AUROC 0.779 | ρ **0.515**, AUROC 0.874 | ρ 0.52 |
| ONC212 | ρ 0.253, AUROC 0.688 | ρ **0.363**, AUROC 0.753 | ρ 0.52 |

Cross-activator 0.777 / 0.864. So regression and classification land in the same territory
(ADEP4 0.874 vs 0.857; ONC212 0.753 vs 0.750), which is the reassuring outcome — the conclusion does
not hinge on the framing. Regression's one real advantage was retaining magnitude; classification's
is comparability with v1 and with the published literature, which is all AUROC.


### 2026-09-17 — retrofit to TabPFN-3.5

**What changed.** The estimator, and only the estimator. Features (ESM-C), labels, homology
clusters, folds, seeds and every control are unchanged.

**Why.** A 3-estimator × 3-embedding × 2-activator grid on byte-identical cluster-grouped folds
(`output/results/degradability/proteomelm_vs_esmc.tsv`, 18 arms). TabPFN beat the forest on
**PR-AUC in 6 of 6 arms, each winning all 5 seeds**: +0.041 / +0.014 / +0.056 on ADEP4 across
esmc/proteomelm/prott5, +0.033 / +0.028 / +0.021 on ONC212. On ROC-AUC the gain is smaller and
clears only with ESM-C. **lazy-qsar 3.4.4 was tested in the same grid and rejected** — its one real
gain (ADEP4 PR +0.021, 5/5 seeds) does not transfer to ONC212 (+0.005, 3/5).

**Use the PAIRED test.** Every head ran identical folds, so shared partition difficulty cancels.
The unpaired 2-SD bar called 1 of 4 arms "no difference" that the paired test resolves decisively;
on ONC212 both arms carry SD ~0.008 from the splits alone.

**Two corrections the paired test forced on earlier conclusions.**
- **ProteomeLM is not "equal" to ESM-C — it is worse.** Paired, `forest × ProteomeLM` loses to
  `forest × ESM-C` on ADEP4 by −0.0039 ROC-AUC and −0.0221 PR-AUC, **winning 0 of 5 seeds on both**.
  It never wins anywhere in the grid.
- **Estimator and features interact.** ProtT5 is the *worst* features under a forest (−0.0152,
  0/5 seeds) and the *best* under TabPFN (+0.0107, 5/5). Benchmarking embeddings with one
  convenient model would have discarded the best pairing. ESM-C stays regardless: TabPFN × ProtT5
  wins ADEP4 but TabPFN × ESM-C wins ONC212, a split decision, and ESM-C adds no new dependency.

**The acceptance test was wrong, and that is the most important finding here.** It was framed as
"do the five consistent COG directions survive?" — but that sentence described the *forest's*
behaviour, so using it to judge a new model assumed the incumbent was correct. Measured against the
**labels** instead: under ONC212, COG O (chaperones) is **not** enriched among measured hits (17.1%
vs a 24.6% base rate, OR 0.63, p 0.42), yet the forest predicted 1.8–2.8. TabPFN predicts
0.86–1.03 — it tracks the measurement. Same for J under ADEP4: label OR 1.81, TabPFN 1.47–2.59,
forest 3.1–10.2. **Judge enrichment against the labels, never against the previous model.**

**Caveats, stated rather than buried.** TabPFN over-depletes **G** (predicts 0.04–0.30 where the
labels show 0.74 / 0.65, both non-significant). `unclassified` is strongly enriched (OR 5.2–7.6) and
this is **unexplained** — tested against protein length and rejected (TabPFN ρ −0.62 vs forest
−0.61, a 0.02 difference). The two activators' predictions now correlate at **ρ 0.89** (forest
0.83–0.86), so they are *less* independent than before. **Forest and TabPFN agree on only 62–64% of
the Kp top 10%** — 208–219 of 573 proteins differ, so any shortlist built on the old outputs has
genuinely moved. Both models track length far more strongly than the labels do (−0.61 vs −0.33),
which is a standing caveat on the axis, not a property of the swap.

**Operational.** ~160 hosted calls ≈ 1.6M of a 20M monthly credit quota; cost is flat per *call*
(10,000 whether the test set is 350 or 2,889 rows) and `fit` is free. A GCS 500 killed the first
run at call 64, so calls now retry 5× with backoff, and the content-addressed cache makes any
interruption cost only the calls in flight. `model_<activator>.joblib` became
`model_<activator>.npz`: TabPFN learns in context, so the artifact is the training set + config.

## Traps

- **The whole label set joins to nothing by identifier.** Three strains, zero ID matches. Go through
  sequence; see the §"identifier problem" before trying anything else.
- **ONC212 has no p-values.** Its calls rest on effect size alone, which is part of why it is the
  weaker label (0.750 vs 0.857).
- **Jacques writes `"NA"` as a literal string**, not an empty cell; without `na_values=["NA"]` the
  column arrives as `object`. Conlon's header is on row 4 of Table S1 and row 6 of Table S2, and
  `Sample_ A` has a literal space after the underscore. This stage sidesteps all of it by reading
  v1's audited merge rather than the xlsx.
- **`legacy/HISTORY.md` §7's warning about 45-/35-row hand-curated substitutes does NOT apply here.**
  That concerns `data/raw/legacy/clp_substrates/` — *E. coli*, ClpXP, a different axis. The
  *S. aureus* supplements are the publishers' own binaries, `fetch_status.tsv` marks both
  `ok/cached`, row counts match the papers, and §10.6 records a passing internal-consistency audit
  (recomputing Conlon's `Average` from raw intensities, max |Δ| = 0.0010).
- **A good degrader target is stable natively but degradable when forced.** Do not validate this axis
  against natural half-life: v1 measured `activator_evidence` vs Nagar half-life at **ρ = −0.073,
  p = 0.17** and warned that *"optimising a degradability score against natural half-life could
  select the **wrong** proteins."*
- **ESM-C embeddings encode species and length**, so a model trained on *S. aureus* partly learns
  "*S. aureus*-ness". The length comparison rules out length specifically; it does **not** rule out
  species identity, which cannot be tested without Gram-negative labels.
- **Every *E. coli* and *K. pneumoniae* row is a ranking hypothesis**, not a measurement. v1's `10f`
  carried a literal `validation = "none — cross-species, unlabelled"` column for this reason.

## Open leads

- **Cleavage as its own axis.** 625 (Conlon) and 1,405 (Jacques) proteins have cleavage log2FCs on
  disk, unused. Between-activator reproducibility is only ρ = 0.22, so expect less — but it is the
  more mechanistically direct readout.
- **The Gr-ADI WP1 assay.** `degradability_datasets.md` §10.4: *"Even 20–30 proteins with a binary
  call would outweigh all 3,262 turnover labels."* The measurement that would let a Gram-negative
  prediction be checked at all.
- **Won 2024** — 72 *M. smegmatis* proteins with measured induced-proximity degradation rate
  constants: the only external set, and the right *modality*. Weak (ClpC1, Actinobacteria) but real.
- **Feng 2013** (*S. aureus* ClpXP/ClpCP) is a placeholder on disk, ACS-gated, never obtained.
