# Localization axis log (docs §5.1)

Implemented `scripts/09a–09h` + `src/localization.py`. Per-protein subcellular localization for
K. pneumoniae HS11286 (5,728) and E. coli K-12 (4,403), keyed by UniProt accession. Emits a
`localization` class plus a graded [0,1] `clp_accessibility` — reachability by the cytoplasmic Clp
machinery, the gating requirement for BacPROTAC-style degradation.

## What changed and why

The axis previously shipped **UniProt-curated labels only**, leaving `localization = unknown` and a
null `clp_accessibility` for **3,414/5,728 Kp (59.6%)** and **2,139/4,403 Ec (48.6%)** proteins,
which silently excluded all of them from the webapp's `degrader` preset. `09b` was written to fill
that tail with PSORTb 3.0 under Docker; it never produced a single prediction (the container hung
~19 h under Rosetta, commit `12ded33`).

Two decisions fixed it:

1. **DeepLocPro replaces "run PSORTb ourselves" as the predictor.** It is ESM-2 650M + an
   attention-pooled head over a 20-checkpoint ensemble, is pure PyTorch (so it runs natively on
   Apple Silicon — the Rosetta problem simply does not arise), always returns a call, and beats
   PSORTb 3.0 on the post-2010 Gram-negative benchmark: accuracy 0.74 vs 0.34, macro-F1 0.75 vs
   0.35, multiclass MCC 0.69 vs 0.30 (Moreno et al. 2024, *Bioinformatics* 40:btae677).
2. **PSORTb is kept, but taken precomputed.** The Brinkman lab publishes PSORTb 3.0 results for
   every RefSeq genome via PSORTdb, including both anchors. Same predictor, same version, zero local
   compute — and because its SCL-BLAST / motif / SVM modules share no machinery with a protein
   language model, agreement with DeepLocPro is genuine corroboration rather than two views of one
   embedding.

## Pipeline

`09a` UniProt → `09b` PSORTdb · `09c` DeepLocPro · `09d` TMbed · `09e` SignalP/lipobox ·
`09f` STEPdb + ortholog transfer → `09g` merge → `09h` slide.
`src/localization.py` reuses the ligandability and essentiality helpers and adds the class
vocabulary, the per-source translation tables, TMbed parsing, the Lol "+2 rule", the calibrated
lipobox, the categorical ortholog transfer, and the accessibility ladder.

## Sources

| Track | Source | Notes |
|---|---|---|
| 09a | UniProt REST (`cc_subcellular_location`, `ft_signal`, `ft_transmem`, `ft_lipid`) | Splits **experimental** (ECO:0000269 / ECO:0007744) from curated-by-similarity. Kp: **1** experimental vs Ec: **803** — the finding that motivates 09f. |
| 09b | PSORTdb per-genome download, `assembly=446671` (HS11286) / `449203` (MG1655) | POST needs a body or it 411s. Keyed by RefSeq (`ref|YP_…`); joined via UniProt ID-mapping, **100% of Kp / 93.4% of Ec** land inside the reference proteome. |
| 09c | DeepLocPro 1.0, `gradi-loc` env | ~2.3 proteins/s on MPS. |
| 09d | TMbed (Rostlab, Apache-2.0) | Per-residue `S/H/B/i/o`. Chosen over DeepTMHMM, which ships only via BioLib (Docker or cloud). |
| 09e | SignalP 6.0 (licensed) or lipobox fallback | See caveats. |
| 09f | STEPdb 2.0 (`stepdb.eu`, K-12 basic proteome CSV) + 03a orthology | 3,625 Ec proteins mapped, 1,880 with a cited reference. |

## Implementation notes worth keeping

- **`gradi-loc` must stay separate from `gradi`.** DeepLocPro needs `fair-esm`, which installs the
  same top-level `esm` package as the EvolutionaryScale `esm` that `01a_esmc_embeddings.py` uses.
- **DeepLocPro's own `-d mps` flag is broken.** `EnsembleModel.embed_batch()` gates device placement
  on `torch.cuda.is_available()`, so with `-d mps` the weights move to MPS while the tokens stay on
  CPU → *"Placeholder storage has not been allocated on MPS device"*. `09c` therefore drives
  `EnsembleModel` directly rather than shelling out to the CLI, which also buys per-protein caching
  and the full probability vector (used as `localization_confidence` and `dlp_margin`).
- **Pins:** `setuptools<81` (DeepLocPro imports `pkg_resources`), `transformers==4.44.2` (TMbed's
  `T5Tokenizer` load fails on transformers 5.x with a spurious tiktoken error).
- **TMbed runs on CPU here** — it gates GPU on `torch.cuda.is_available()`. Patching it to MPS would
  mean overriding a module-level device global plus `torch.use_deterministic_algorithms`, which is
  not cheaply verifiable, so it was left alone. `09d` shards and caches instead, and the shard cache
  key includes the shard size so changing `--shard-size` cannot silently reuse mismatched shards.
- **Categorical ortholog transfer needed its own rule.** `E.transfer_ecoli_to_kp` reduces with
  max/mean, meaningless for a class label. The kp→ec orthology in 03a is OrthoFinder orthogroup
  membership only — `pident`/`coverage`/`bitscore` are entirely empty — so there is no identity to
  threshold on either. `transfer_categorical_ecoli_to_kp` uses donor consensus and **abstains on
  ties** (only 81 of 3,179 anchors have >1 E. coli ortholog).
  **This limitation is not specific to localization** — any axis planning to threshold kp→ec
  transfer on percent identity will hit the same empty columns and should check before designing
  around them.
- **The peripheral-membrane trap in the UniProt CC text.** A protein annotated
  `"Cytoplasm. Cell inner membrane; Peripheral membrane protein."` (DnaK is the archetype) was being
  classified `inner_membrane`, because the keyword cascade tested "inner membrane" before
  "cytoplasm" and CC order carries no priority. That understates exactly what this axis measures: a
  peripherally-attached protein on the cytoplasmic face is as ClpXP-reachable as a soluble one, and
  it is the distinction `classify_stepdb()` already drew explicitly for STEPdb's `F1` class — so the
  UniProt parse was also *internally inconsistent* with the STEPdb parse. `classify_uniprot()` now
  resolves peripheral membrane proteins by the side they sit on before the compartment cascade
  runs. **43 E. coli proteins** carry `peripheral membrane protein` + `cytoplasmic side`; 40 were
  mis-binned. Caught by the known-marker spot check, which is the argument for keeping that check.
- **The β-barrel threshold is ≥8 predicted strands** (`tmbed_features`). Gram-negative OM barrels
  start at 8; one or two stray predicted strands are noise. Validated by the strand counts below —
  and by TolC, which sits at 4 because it is a trimer contributing 4 per monomer, and is correctly
  not flagged.
- **The topology overrides are a safety net, not a workhorse.** Across both organisms the β-barrel
  override fired **once** and the Lol "+2 rule" **5 times**, because the higher-precedence sources
  had already got those proteins right. `membrane_resolved` is the one that does real work (329 Kp /
  24 Ec) — that is UniProt's side-less "Membrane" label being resolved from topology. Worth knowing
  before anyone concludes the override machinery is redundant: its value is bounding the worst case,
  not its hit rate.
- **A DIAMOND fallback was deliberately not written** in 09b — see
  `docs/localization_downloads.md` §C. At 100% / 93.4% accession coverage it would have been dead
  code that never executed; there is a loud coverage assert and a documented remedy instead.

## Merge precedence (09g)

`uniprot_experimental` → `ortholog_transfer` / `stepdb` → `stepdb_curated` → `uniprot_curated` →
`deeplocpro` → `psortb` → `unknown`.

TMbed/SignalP may overrule a **predicted** call only (never an experimental or curated one), in the
two cases they settle definitively: a transmembrane β-barrel is an outer-membrane protein
(mis-calling one as cytoplasmic is this axis's worst-case error), and a Sec/SPII lipoprotein is
sorted by the Lol "+2 rule" — PSORTb's characteristic failure is binning OM lipoproteins as
periplasmic. The provisional `membrane` label (UniProt "Membrane", no side stated) is resolved from
topology regardless of source.

`clp_accessibility`: cytoplasm 1.0; inner membrane 0.6 if ≥30% of residues face the cytoplasm else
0.2 (0.4 when topology is missing); periplasm 0.2; outer membrane / β-barrel / extracellular /
cell surface 0.0; unknown → blank.

## Result

All of §5.1 is complete, run end-to-end and exported. Coverage went from **40.4% Kp / 51.4% Ec** to
**100% / 100%**, with zero null `clp_accessibility`.

| | Kp | Ec |
|---|---|---|
| localized | 5,728 (100%) | 4,403 (100%) |
| experimental | 1,681 | 2,031 |
| curated | 1,173 | 708 |
| predicted | 2,874 | 1,664 |
| `clp_accessibility` ≥ 0.5 (shortlist) | 4,174 | 3,139 |
| β-barrels (TMbed) | 67 | 66 |

| `clp_accessibility` | meaning | Kp | Ec |
|---|---|---|---|
| 1.0 | cytoplasm (incl. peripheral IM, cytoplasmic face) | 3,470 | 2,645 |
| 0.6 | inner membrane, ≥30% cytoplasm-facing | 704 | 494 |
| 0.2 | IM mostly buried, or periplasmic | 955 | 842 |
| 0.0 | OM / β-barrel / extracellular / cell surface | 599 | 422 |

| compartment | Kp | Ec |
|---|---|---|
| cytoplasm | 3,470 | 2,645 |
| inner_membrane | 1,425 | 1,099 |
| extracellular | 389 | 204 |
| periplasm | 234 | 240 |
| outer_membrane | 193 | 181 |
| cell_wall_surface | 17 | 34 |

Winning source: Kp `deeplocpro` 2,874 · `ortholog_transfer` 1,680 · `uniprot_curated` 1,173 ·
`uniprot_experimental` 1. Ec `deeplocpro` 1,664 · `stepdb` 1,229 · `uniprot_experimental` 802 ·
`uniprot_curated` 708. **`psortb` never wins** — DeepLocPro always has a call and outranks it — which
is by design: PSORTb's role here is corroboration through `predictor_agreement`, not adjudication.

Mean `predictor_agreement` 0.895 (Kp) / 0.919 (Ec).

### Runtime

- **09c DeepLocPro** ~2.3 proteins/s on MPS → ~75 min for both proteomes.
- **09d TMbed** ~4 min per 250-protein shard with the machine to itself, ≈2.5 h for both. An earlier
  ~17 min/shard measurement was taken while DeepLocPro was still competing for CPU — worth
  remembering before extrapolating an ETA from a contended run.
- Both are resumable, so neither needs babysitting.

### Validation

- **Known-marker spot check: 20/20 for E. coli** (OmpA/C/F, LamB, BtuB, LptD, FhuA, TolC, MalE,
  DsbA, OppA, MglB, SecY, AcrB, AtpB, FtsW, GroL, DnaK, RpoB, ClpP, TufA); 9/9 of the markers that
  exist in Kp.
- **β-barrel strand counts are textbook**: OmpA 8, OmpC/F 16, LamB 18, BtuB 22, FhuA 22, LptD 26.
  TolC scores 4 strands and is correctly *not* flagged — it is a trimer contributing 4 strands per
  monomer — and is still called outer-membrane by another source.
- **Topology tracks biology**: mean TM helices are 5.2 for inner-membrane proteins and 0.0 for
  cytoplasmic ones, and β-barrels peak sharply at the outer membrane.
- **DeepLocPro vs PSORTb**: 92.6% (Ec) / 92% (Kp) agreement between two fully independent predictors.
- **Invariant checked**: all 133 β-barrels across both organisms have `clp_accessibility == 0.0`.
  22 of them carry a non-`outer_membrane` label — all `uniprot_curated`, and all autotransporter /
  usher-family proteins whose *passenger* domain genuinely is extracellular while the β-barrel
  translocator sits in the OM. The score is unaffected because `clp_accessibility()` short-circuits
  on `is_beta_barrel` before consulting the label, and the override deliberately does not overrule a
  curated call.

## Caveats

- **DeepLocPro training-set leakage on E. coli.** It was trained on UniProt 2023_03 + PSORTdb 4.0, so
  K-12 proteins are almost certainly in its training set. Its 84.1% agreement with UniProt-curated
  E. coli labels is a **sanity check, not independent validation**. The meaningful number is
  DeepLocPro vs PSORTb — two fully independent predictors — at **92.6%**.
- **SignalP 6.0 is not installed.** It is licensed software requiring manual acceptance of the DTU
  academic terms (route in `docs/localization_downloads.md` §F); `09e` runs without it using a
  lipobox motif. The fallback finds lipoproteins only — it does not type Tat/SPI, Tat/SPII or
  Sec/SPIII. Set `SIGNALP6_BIN` and re-run `09e` + `09g` to upgrade the track.

  **How the lipobox was calibrated** — against UniProt's 99 E. coli lipid-anchor entries as truth, so
  nobody re-derives it. The motif alone is unusable; what makes it work is adding the other two
  defining features of a signal peptide (charged n-region, hydrophobic h-region):

  | rule | n called | precision | recall | F1 |
  |---|---|---|---|---|
  | broad motif `[LVIFMTA][ASTVILGMF][GAS]C`, Cys ≤40 | 256 | 0.36 | 0.94 | 0.52 |
  | Prosite PS51257 `[LVI][ASTVI][GAS]C`, Cys ≤40 | 163 | 0.52 | 0.85 | 0.64 |
  | + Cys in [14, 35] | 138 | 0.59 | 0.83 | 0.69 |
  | + h-region hydropathy ≥ 1.0 (Kyte-Doolittle, 12 aa) | 117 | 0.70 | 0.83 | 0.76 |
  | **+ K/R in the first 8 residues** ← shipped | **109** | **0.74** | **0.82** | **0.78** |

  Tightening hydropathy to ≥2.0 collapses recall to 0.53 — the h-region of a lipoprotein signal
  peptide is shorter and less hydrophobic than a classic Sec/SPI one, so do not push it further.
  Parameters live in `src/localization.py` (`LIPOBOX_*`). Note they were tuned on E. coli and applied
  unchanged to Kp, which is the usual mild overfitting risk; the filters are textbook signal-peptide
  architecture rather than arbitrary constants, which is the argument that it transfers.
- **STEPdb evidence tier is a proxy.** Its CSV mixes hand-curated and inferred assignments without an
  explicit flag, so a cited reference in `Annotation References` is used as the stand-in for bench
  support (2,114/3,897 rows). The file is semicolon-delimited with semicolons inside free-text cells,
  so a minority of rows shift; those fall through to `unknown` rather than to a wrong compartment.
- **DeepLocPro's GitHub LICENSE is CC BY-NC-SA 4.0** (non-commercial), while the paper states
  CC BY 4.0. Worth confirming against the GraDi collaboration's partners before this goes into a
  commercial deliverable.
- **Sequences are truncated at 2,000 residues** for DeepLocPro (ESM-2 cost is quadratic in length).
  Localization signal is overwhelmingly N-terminal; affected proteins carry `dlp_truncated`.

## Figures

Four slides per organism, all stylia "slide" format, NPG palette, every panel single-organism. The
compartment order and palette live in `src/localization.py` (`LOC_CLASS_ORDER`, `LOC_CLASS_COLOR`)
so colours cannot drift between figures.

| Script | Output | What it is for |
|---|---|---|
| `09h` | `09h_localization_{kp,ec}.png` | The axis summary: composition, coverage gained, evidence tier, DeepLocPro-vs-PSORTb concordance, topology, accessibility ladder. |
| `09i` | `09i_atlas_{kp,ec}.png` | **Compartment atlas** — small multiples over the ESM-C map, one panel per compartment. |
| `09j` | `09j_evidence_{kp,ec}.png` | Where the calls came from: source × compartment, corroboration depth, agreement, confidence, DeepLocPro decisiveness and what it conflates. |
| `09k` | `09k_topology_{kp,ec}.png` | Structural evidence: the 0.30 cut, β-barrel strand counts, TM helices per compartment, architecture composition, export signals, overrides. |

Two results worth reading off the figures:

- **Localization is legible in sequence space.** On the ESM-C map the compartments occupy visibly
  distinct territory — cytoplasm dominating one lobe, inner membrane a separate cluster, periplasm
  and outer membrane concentrated in their own regions. The axis is describing structure the protein
  language model already sees, not a label bolted on top.
- **DeepLocPro is least decisive on extracellular proteins.** Its within-model probability profile
  (09j panel 6) is diagonal at 0.90 cytoplasm / 0.93 inner membrane / 0.86 periplasm / 0.82 outer
  membrane but only **0.63 extracellular**, leaking 0.17 back to cytoplasm. Extracellular is also
  where Kp has *no* experimental evidence at all, so that class is the axis's weakest corner and
  should be treated with the most caution.

One thing 09k surfaces that is easy to misread as a bug: **231 Kp proteins labelled inner-membrane
have `cyto_residue_fraction = 1.0` and no TM helix at all** (193 of them also have zero UniProt
transmembrane segments). These are peripheral, membrane-associated proteins rather than
membrane-spanning ones. Scoring them 0.6 is intended — they sit on the cytoplasmic face and are
Clp-reachable — and the panel annotates them rather than hiding the spike.

## Artifacts and where they live

Code, docs and the webapp payload are in Git (commit `c7c7222`, branch `localization-axis`).
Everything under `data/` and `output/` goes through **eosvc, and has not been pushed yet**:

| | path | regenerable? |
|---|---|---|
| per-track + merged tables (20 CSVs) | `output/results/<org>/<prefix>_loc*.csv`, `<prefix>_localization*.csv` | yes, but cheap to keep — this is the deliverable |
| slides | `output/plots/09h_localization_{kp,ec}.png` | yes |
| PSORTdb tables + id-mapping cache | `data/raw/<org>/localization/psortdb/` | yes (network) |
| STEPdb CSV | `data/raw/ecoli/localization/stepdb/` | yes (network) |
| DeepLocPro per-protein cache | `data/processed/<org>/localization/deeplocpro/` | yes (~75 min compute) |
| TMbed shard cache | `data/processed/<org>/localization/tmbed/` | yes (~2.5 h compute) |

The two prediction caches are the bulk of the ~53 MB. Stale artifacts from the pre-rework track that
nothing reads any more are listed in `docs/localization_downloads.md` §H.

## Open leads (not done)

Carried over from the §5 spec, unchanged by this build:

- **PRED-TMBB2 β-barrel cross-check** — a second opinion on the call that matters most. TMbed is
  currently the only source of it.
- **Structure-derived surface exposure (new 5.3)** — DSSP/SASA on the AlphaFold models already
  fetched by `04a`. Would give a continuous exposure measure rather than a compartment proxy, and
  `src/degradability.py` already has an `exposure` helper.
- **The whole of §5.2** — expression/abundance (PaxDb, condition-specific Kp proteomics). Localization
  answers *can Clp reach it*; §5.2 answers *is it there at all*, and the axis is only half its name
  until that lands.
