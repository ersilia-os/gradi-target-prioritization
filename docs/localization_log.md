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

`clp_accessibility` distribution (Kp): 0.0 → 599, 0.2 → 955, 0.6 → 704, 1.0 → 3,470.

TMbed took ~4 min per 250-protein shard once it had the machine to itself (~2.5 h for both
proteomes), not the ~17 min/shard measured while DeepLocPro was still competing for CPU.

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
  academic terms; `09e` runs without it using a lipobox motif calibrated against UniProt's 99 E. coli
  lipid-anchor entries (**precision 0.74 / recall 0.82**; the bare Prosite motif alone is 0.36
  precision). The fallback finds lipoproteins only — it does not type Tat/SPI, Tat/SPII or Sec/SPIII.
  Set `SIGNALP6_BIN` and re-run `09e` + `09g` to upgrade the track.
- **STEPdb evidence tier is a proxy.** Its CSV mixes hand-curated and inferred assignments without an
  explicit flag, so a cited reference in `Annotation References` is used as the stand-in for bench
  support (2,114/3,897 rows). The file is semicolon-delimited with semicolons inside free-text cells,
  so a minority of rows shift; those fall through to `unknown` rather than to a wrong compartment.
- **DeepLocPro's GitHub LICENSE is CC BY-NC-SA 4.0** (non-commercial), while the paper states
  CC BY 4.0. Worth confirming against the GraDi collaboration's partners before this goes into a
  commercial deliverable.
- **Sequences are truncated at 2,000 residues** for DeepLocPro (ESM-2 cost is quadratic in length).
  Localization signal is overwhelmingly N-terminal; affected proteins carry `dlp_truncated`.
