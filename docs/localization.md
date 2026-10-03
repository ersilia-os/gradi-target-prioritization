# Stage 03 — localization

`scripts/localization/predict.py` → `scripts/localization/merge.py` · worker `scripts/localization/workers/deeplocpro.py` · loader `src/localization.py`

One subcellular compartment per protein, plus a topology score, assigned **from sequence alone** —
no accession lookup, no annotation transfer, no curated database. The same method for all three
species, so the label means the same thing in each.

**One table, two predictors side by side** — `localization_<species>.tsv`:

| column | from | what it is |
|---|---|---|
| `localization` | DeepLocPro 1.0 | the compartment call |
| `cytoplasmic_fraction` | TMbed | the topology score |

**`confidence` is not in the deliverable** (owner's call, 2026-10-03) — byte-identical in
`evidence/deeplocpro_<species>.tsv`. This one has teeth: DeepLocPro **always** returns a call and
has no `unknown` class, so every row looks equally decided, and `confidence` was the only thing
distinguishing them:

| | median | < 0.7 | < 0.5 |
|---|---|---|---|
| Kp | 0.975 | 841 (14.7%) | 246 (4.3%) |
| Ec | 0.982 | 524 (11.9%) | 108 (2.5%) |
| Sa | 0.964 | 421 (14.6%) | 55 (1.9%) |

Below 0.5 the winning class holds less than half the probability mass. Join it back before
trusting a single label; `load_probabilities()` gives the full un-remapped six-vector.

**`has_signal_peptide` is not in the deliverable either** (owner's call, 2026-10-03). It ships
byte-identically in `evidence/tmbed_<species>.tsv` via `load_tmbed()` — Kp 604 (10.5%), Ec 564
(12.8%), Sa 211 (7.3%).

**Join it back before using localization as a degradability filter.** It is the only column that
says *why* `cytoplasmic_fraction` is near zero — **exported** rather than **membrane-buried** — and
those have opposite consequences for activated ClpP: a secreted protein transits the cytoplasm
unfolded and is reachable (measured hit rate 0.049/0.346), a membrane protein is co-translationally
inserted and never presents a soluble chain (0.031/0.105). `degradability/enrichment.py` and
`plots/localization.py` both join it explicitly for that reason.

**No `evidence` column, unlike every other axis.** Both predictors cover 100% of every proteome by
construction — DeepLocPro always returns a call and TMbed labels every residue — so the column was
constant across all 13,020 proteins and carried no information. The check it would have encoded is
not lost: `merge.py` derives per-row coverage and **exits non-zero if either track is silent for a
single protein**, which is a broken run rather than a sparse one. Do not re-add it as a constant.

**The two are never reduced to one call.** They answer different questions and neither derives
from the other, so `merge.py` stacks the columns and arbitrates nothing — where they disagree,
that disagreement is the information, and TMbed corroborating `extracellular` from outside
DeepLocPro is the only cross-check this axis has. **Prefer `cytoplasmic_fraction` over
`localization` where a choice is forced.**

The per-predictor tables stay in `evidence/`, not because they are secondary but because the two
tracks run and resume independently — TMbed alone is CPU-only hours — so `predict.py --only tmbed`
must have somewhere to write. `merge.py` recomputes nothing and takes seconds.

## Why this axis exists, and why sequence only

The cytoplasmic Clp machinery can only degrade what it can physically reach, so localization is the
gate on BacPROTAC-style targeting. A protein with no assignment is a protein silently dropped from
every downstream shortlist — which is exactly what v1 did. Its axis shipped UniProt-curated labels
only, leaving **59.6% of Kp** and **48.6% of Ec** at `localization = unknown`, and those proteins
vanished from the webapp's `degrader` preset without a word.

v1 eventually reached 100% coverage with a six-track merge (`legacy/scripts/09a`–`09g`). Four of
those tracks are unusable here, and the reason each one goes is worth keeping:

| v1 track | why it goes |
|---|---|
| `09a` UniProt `cc_subcellular_location` | not sequence — and annotation depth is uneven, so the label would mean different things per species. UniProt has **1** experimental localization for Kp HS11286 against **803** for Ec K-12. |
| `09b` PSORTb via PSORTdb | keyed on RefSeq assembly ids, and in v1 it **never won a single call** — DeepLocPro always outranked it. |
| `09e` SignalP 6.0 | licensed, never installed. The shipped lipobox fallback was E.-coli-tuned and found lipoproteins only. |
| `09f` STEPdb 2.0 + ortholog transfer | *E. coli* only, and needs orthology this pipeline does not have yet. |

None of the four helps *S. aureus*, which is new in v2. What is left — DeepLocPro and TMbed — is what
was actually load-bearing, and is genuinely sequence-only.

**Deliberately not computed here:** any accessibility score. Turning a compartment plus a fraction
into one "can Clp reach this" verdict is a modelling decision, and it belongs to whichever stage
consumes it. v1 computed a `clp_accessibility` ladder that, per `legacy/HISTORY.md` §7 item 9,
**nothing ever consumed**.

## The six compartments

DeepLocPro 1.0 (Moreno et al. 2024, *Bioinformatics* 40:btae677) is ESM-2 650M with an
attention-pooled head, ensembled over 20 nested-CV checkpoints, emitting a softmax over six
prokaryotic compartments. It beats PSORTb 3.0 across the board on the post-2010 Gram-negative
benchmark — accuracy **0.74 vs 0.34**, macro-F1 0.75 vs 0.35, MCC 0.69 vs 0.30 — and, unlike PSORTb,
**always returns a call**. That is what makes 100% coverage a property of the method rather than
something to chase, and it is why there is no `unknown` class.

Ordered inward → outward through the cell envelope, which is also the order Clp reachability falls
in (`LOC_CLASSES` in `src/localization.py`):

| label | DeepLocPro string | abbrev | colour |
|---|---|---|---|
| `cytoplasm` | `Cytoplasmic` | Cyt | `#00A087` |
| `cytoplasmic_membrane` | `Cytoplasmic Membrane` | CM | `#F39B7F` |
| `periplasm` | `Periplasmic` | Peri | `#E64B35` |
| `outer_membrane` | `Outer Membrane` | OM | `#8491B4` |
| `cell_wall_surface` | `Cell wall & surface` | CW | `#7E6148` |
| `extracellular` | `Extracellular` | Ext | `#3C5488` |

Two changes from v1's vocabulary, which declared eight labels:

- **`inner_membrane` → `cytoplasmic_membrane`.** *S. aureus* has no inner membrane, it has one
  membrane. `cytoplasmic_membrane` is DeepLocPro's own term and is correct in either Gram type.
- **`membrane` and `unknown` are gone.** `membrane` (side unstated) was an artifact of parsing
  UniProt free text, which this stage does not do; `unknown` was a coverage gap that cannot exist
  when the predictor always calls.

The palette is v1's NPG one, carried over verbatim so figure colours cannot drift between versions.

### The Gram-positive trap — read this before comparing *S. aureus*

DeepLocPro takes a Gram group, and `positive` does **not** merely mask the two Gram-negative-only
classes. `DeepLocPro/utils.py:remap_probabilities` **adds their probability mass into
`Extracellular`**:

```python
def remap_probabilities(probs):
    '''Sum outer membrane and periplasm into extracellular'''
    probs_fixed[1] = probs_fixed[1] + probs_fixed[4] + probs_fixed[5]   # 1 = Extracellular
    probs_fixed[4] = 0   # Outer Membrane
    probs_fixed[5] = 0   # Periplasmic
```

So *S. aureus* has **four** reachable classes, and its `extracellular` count absorbs whatever the
model wanted to call periplasmic. A zero in the Kp/Ec-only columns for Sa is **structural, not
missing data**. v1 never met this — it had no Gram-positive organism.

The remap is arithmetic on the output vector, not a different forward pass. So the model runs
**once**, the raw un-remapped six-vector is persisted to
`evidence/deeplocpro_probabilities_<species>.tsv`, and the remap is applied only when choosing a
label. Costs nothing, and it keeps the inflation measurable: `p_periplasm + p_outer_membrane` is
exactly the mass that was moved. The script prints the mean every run.

## TMbed — what the score actually is

TMbed (Rostlab, Apache-2.0) is **not a second localization classifier**. It labels every residue with
one of seven states, from ProtT5-XL-U50 embeddings. Alphabet at `--out-format 4`
(`tmbed/tmbed.py:252`):

| label | meaning |
|---|---|
| `H` / `h` | transmembrane α-helix, inside→outside / outside→inside |
| `B` / `b` | transmembrane β-strand, inside→outside / outside→inside |
| `S` | signal peptide |
| `i` | not in a membrane, **cytoplasmic** side |
| `o` | not in a membrane, outside (periplasm / extracellular) |

One protein comes back as a string the same length as its sequence:

```
iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii        soluble cytoplasmic
SSSSSSSSSSSSooooooooooooooooooooooooooooooo       secreted: signal peptide, then all outside
iiiiiHHHHHHHHHHHHHHHHHHooooooHHHHHHHHHHHiiiii     2-TM protein, both termini cytoplasmic
SSSSSSSSSoBBBBBBBBBoooooobbbbbbbbbooooBBBBBBBo    outer-membrane β-barrel
```

**`cytoplasmic_fraction` = `count('i') / length`** — one number in [0, 1] answering literally *what
fraction of this protein sits in the cytoplasm*, which is what Clp reachability reduces to.

| value | reads as | measured mean, Kp |
|---|---|---|
| ~1.0 | wholly cytoplasmic; no membrane insertion, no export signal | `cytoplasm` 0.992 |
| ~0.4–0.5 | membrane protein with substantial cytoplasm-facing loops or domains | `cytoplasmic_membrane` 0.445 |
| ~0.2 | β-barrel — the periodic turns between strands count as cytoplasmic | `outer_membrane` 0.239 |
| ~0.0 | exported through Sec: signal peptide, then everything outside | `periplasm` 0.064 |

Note the third row. A barrel is **not** near zero, because its strands are joined by short
cytoplasmic turns, so it scores *above* a fully-exported soluble protein. Use `n_tm_strand` from
`evidence/tmbed_topology_<species>.tsv` to tell the two apart — the fraction alone will not.

**`has_signal_peptide` = `'S' in labels`** says *why* a fraction is near zero — exported rather than
membrane-buried — and makes one contradiction visible: a signal-peptide-bearing protein that
DeepLocPro called `cytoplasm` is almost certainly a mis-call.

**`--out-format 4`, not v1's format 1.** Format 1 collapses `h` into `H` and `b` into `B`, discarding
segment orientation, in the same 3-line file, for no saving at all. `i`/`o` counting is unaffected,
so `cytoplasmic_fraction` is identical — but helix and strand segments must then be counted **per
letter** (`HHHhhh` is two segments, one each way), which `_n_segments()` does.

### Why the other features go to `evidence/` rather than being dropped

TMbed is **CPU-only here** — it gates GPU on `torch.cuda.is_available()` (`tmbed/embed.py:34`) — so
the run is hours and **cannot be partially redone**. Extracting numbers from the label string costs
nothing, so discarding a free one would be the expensive choice.

- **`n_tm_helix`** disambiguates the score: `cytoplasmic_fraction = 1.0` with `n_tm_helix = 0` is a
  soluble cytoplasmic protein, `= 0.55` with `n_tm_helix = 6` is a polytopic membrane protein
  reachable by its loops. v1 found **231 Kp proteins** at fraction 1.0 with **zero** TM helices —
  peripheral, membrane-associated, and indistinguishable from a bug without this number.
- **`n_tm_strand`** is β-barrel evidence, and the one place TMbed can catch DeepLocPro's worst-case
  error (a barrel called cytoplasmic). It is the axis's best-validated signal, and the strand-count
  spot check reads from here — so this evidence column is load-bearing for the guards even though
  it is not in the deliverable. `BETA_BARREL_MIN_STRANDS = 8`.
- **`cyto_longest_segment`** — the longest contiguous cytoplasmic run. A threading handle needs
  contiguity, not just total mass.
- **The raw per-residue label string** — so *any* future topology feature is derivable with zero
  recompute.

## Two environments

Both predictors need **`gradi-loc`**, never `gradi`: DeepLocPro depends on `fair-esm`, which claims
the same top-level `esm` package as the EvolutionaryScale `esm` that stage 01 uses for ESM-C.
Installing it into `gradi` silently breaks embedding generation.

`gradi-loc` also carries much newer pandas and numpy than `gradi`, so the split is enforced by
process boundary, not by import:

- `scripts/localization/predict.py` runs in **`gradi`** and owns every table, the vocabulary and the guards.
- `scripts/localization/workers/deeplocpro.py` runs under the **`gradi-loc`** interpreter, imports **stdlib + torch
  only** (no pandas, no `src/`), and communicates one-way by file: FASTA in, one JSON per protein out.
- **TMbed needs no worker** — its `tmbed` console script already has `gradi-loc`'s interpreter in its
  shebang, so the stage shells the binary directly.

Override the env location with `GRADI_LOC_BIN`. This is a stronger separation than stage 02's, which
only borrows the `rpsblast` binary.

### Why we drive DeepLocPro's model instead of its CLI

1. **Its `-d mps` flag is broken.** `EnsembleModel.embed_batch()` gates device placement on
   `torch.cuda.is_available()` (`DeepLocPro/model.py:45`), so with `-d mps` the weights move to MPS
   while the tokens stay on CPU → *"Placeholder storage has not been allocated on MPS device"*.
   `predict_one()` is that method rewritten to put both on the same device.
2. **No resume.** The CLI writes one timestamped `results_<YYYYmmdd-HHMMSS>.csv` at the end; this is
   ~13k single-sequence forward passes and must survive an interruption.
3. It discards the probability vector, and writes a per-protein attention plot we do not want.

## Scope: the three bacteria, and human cannot be added

13,020 proteins. Human is out and cannot be otherwise — DeepLocPro is prokaryote-only, and a
eukaryotic call would need DeepLoc 2.x and a disjoint vocabulary. `--species human` is not offered.
This matches stages 01 and 02.

## Outputs

```
data/processed/localization/
    localization_<species>.tsv                THE DELIVERABLE — 5 columns, keyed on uniprot_ac
    evidence/
        deeplocpro_<species>.tsv                 the compartment call, per-predictor record
        tmbed_<species>.tsv                      the topology score, per-predictor record
        deeplocpro_probabilities_<species>.tsv   raw UN-REMAPPED 6-vector, margin, truncated
        deeplocpro_counts_<species>.tsv          per-compartment counts, all six, zeros included
        tmbed_topology_<species>.tsv             n_tm_helix, n_tm_strand, cyto_longest_segment, tmbed_length
        tmbed_labels_<species>.tsv               the raw per-residue label string
        {deeplocpro,tmbed}_manifest.tsv
    scratch/
        deeplocpro_cache/<species>_L2000/*.json  one JSON per protein (resumable)
        tmbed_shards/<species>_250_<NNNN>.pred   shard cache (resumable)
        .tmbed_<species>.key                     md5(accession set | shard size | out format)
        smoke_*                                  --limit output
```

**The split follows the directory contract's cite-or-delete test.** The six tables and the two
manifests are all cited or checked — the two per-predictor tables because `merge.py` is built from
them and the remap figure recomputes against them; `deeplocpro_probabilities_<species>.tsv` is what keeps the
Gram-positive remap measurable (the `gram_remap.png` figure recomputes the raw argmax from it),
`tmbed_topology_<species>.tsv` carries the strand count the β-barrel guard reads, and
`tmbed_labels_<species>.tsv` is the raw string every future topology feature derives from — so
they are `evidence/`. The two caches and the smoke output are the opposite: large, regenerable,
and deleted to reclaim space, so they are `scratch/`.

No raw directory: both models' weights come from the `gradi-loc` install and the HF/torch caches,
not from a fetch this stage performs.

| helper | returns |
|---|---|
| `load(species)` | **the deliverable** — both predictors, complete and canonical |
| `load_deeplocpro(species)` | the compartment call alone, from `evidence/` |
| `load_tmbed(species)` | the topology score alone, from `evidence/` |
| `load_all(species=SPECIES)` | all species stacked, with a `species` column |
| `load_probabilities(species)` | the raw un-remapped 6-vector, margin, truncated |
| `load_topology(species)` | `n_tm_helix`, `n_tm_strand`, `cyto_longest_segment`, `tmbed_length` |
| `load_labels(species)` | the raw per-residue label string |
| `load_counts(species)` | per-compartment counts |
| `manifest(predictor)` | `"deeplocpro"` or `"tmbed"` |

Vocabulary and parsing also live there: `LOC_CLASSES`, `LOC_CLASS_COLOR`, `LOC_CLASS_ABBREV`,
`DEEPLOCPRO_MAP`, `DLP_CANON`, `GRAM_GROUP`, `GRAM_POSITIVE_CLASSES`, `BETA_BARREL_MIN_STRANDS`,
`DLP_MAX_LENGTH`, `remap_gram_positive()`, `label_from_probs()`, `tmbed_features()`,
`parse_tmbed()`, `read_tmbed_labels()`.

## Figures

`scripts/plots/localization.py` writes `output/plots/localization/`. Three figures, each
answering one question the stage raises:

| figure | panels | what it shows |
|---|---|---|
| `composition.png` | A composition, B confidence | what the three proteomes are made of, and which compartment the model is least sure about |
| `cross_check.png` | A cytoplasmic fraction, B signal peptide, C topology | the stage's real validation — DeepLocPro and TMbed share no machinery, so their agreement is evidence |
| `gram_remap.png` | A before/after, B moved mass | how much of *S. aureus*'s `extracellular` is the remap rather than the model |

Compartment colours come from **`LOC_CLASS_COLOR`, not a stylia palette**: the class-to-colour
mapping is pinned in the vocabulary so it cannot drift between v1 and v2, or between figures here.
Everything is ordered `LOC_CLASSES` — inward → outward through the envelope, which is also the order
Clp reachability falls in.

Two deliberate choices worth knowing:

- **`cross_check.png` panel A is a violin, not a box.** `extracellular` is bimodal, so a box reports
  a median (0.96) that describes almost nothing in the class and hides the split entirely. The mean
  is marked separately because it is the number the run log quotes.
- **Panels pool the three species.** The compartment vocabulary is shared, and the structural
  absences (`cell_wall_surface` only in Sa, `periplasm`/`outer_membrane` only in the
  Gram-negatives) are visible as such. Per-species numbers are printed to stdout, and
  `composition.png` panel A keeps them separate.

`gram_remap.png` is the only figure that recomputes anything — the raw argmax labels, from the
un-remapped six-vector in `evidence/deeplocpro_probabilities_saureus.tsv`, which the stage persists
for exactly this purpose.

## Running it

```bash
python scripts/localization/predict.py                              # both predictors, 3 species
python scripts/localization/predict.py --species ecoli --limit 20   # smoke test -> scratch/smoke_*
python scripts/localization/predict.py --only deeplocpro            # skip TMbed's hours
python scripts/localization/predict.py --only tmbed --refresh       # rebuild every shard
python scripts/localization/predict.py --dry-run
```

CLI: `--species` · `--only {deeplocpro,tmbed}` · `--shard-size` · `--device` · `--threads` ·
`--limit` · `--refresh` · `--dry-run` · `-q/--quiet`.

Run with the **`gradi`** env — it shells out to `gradi-loc` itself. Do not activate `gradi-loc` and
run the stage from there; it imports pandas through `src/`.

### Guards

- **Row-count identity.** Every predictor must return exactly as many rows as the stage-00 table has
  sequences, per species, or the run exits non-zero. This is what makes "every protein has an
  assignment" a checked fact rather than an assumption.
- **Vocabulary.** Every label must be in `LOC_CLASSES`; a Gram-positive species must not receive
  `periplasm` or `outer_membrane` (if it does, the remap did not run); probability vectors must sum
  to 1.0 ± 1e-3.
- **Label/sequence length.** TMbed truncates nothing, so a label string whose length differs from
  the sequence means the wrong sequence was scored, or a 3-line block is misaligned. Both are hard
  failures.
- **β-barrel calls, enforced.** The strand *count* is printed against the textbook value, but what
  is enforced is the *call*: `ompA/ompC/ompF/lamB/btuB/fhuA/lptD` must be ≥ 8 strands and **TolC must
  not be**, because it is a trimer contributing 4 strands per monomer. A drift here means TMbed is
  misconfigured.
- **Marker panel, enforced as a rate.** v1's 20-marker E. coli set, printed per gene, gated at
  `MARKER_MIN_AGREEMENT = 0.80`. Not per gene: DeepLocPro's own benchmark accuracy is 0.74, so
  demanding 20/20 from a predictor is a spurious failure waiting to happen.
- **Truncation is reported, always.** Sequences over `DLP_MAX_LENGTH = 2000` are cut at the
  C-terminus before ESM-2 (attention is quadratic; localization signal is overwhelmingly
  N-terminal). `truncated` is an evidence column, but the count is printed every run — silently
  truncating is what the house style forbids.
- **Cache keys.** `DLP_MAX_LENGTH` is in the DeepLocPro cache **directory** name and the shard size
  is in each TMbed shard **file** name, so changing either cannot silently reuse mismatched results.
  The TMbed key sidecar additionally covers the accession set. `--limit` writes everything under
  `scratch/smoke_*`, so a smoke test can never clobber a full run.
- **Cross-check, reported not enforced.** Mean `cytoplasmic_fraction`, TM helices, strands and
  signal-peptide rate per DeepLocPro compartment. Two independent models, so the separation is the
  evidence they agree: helices concentrate in `cytoplasmic_membrane` (5.2 against 0.0 in
  `cytoplasm`), strands only in `outer_membrane`, signal peptides in `periplasm` and
  `outer_membrane`. A collapse in that separation means one of them is broken.

  **`cytoplasmic_fraction` is NOT monotone across the six compartments, and should not be read as
  if it were.** A β-barrel has periodic cytoplasmic *turns* between its strands, so
  `outer_membrane` (0.239 on Kp) sits well above a fully-exported soluble periplasmic protein
  (0.064). The fraction is monotone among the *soluble* compartments only. This is geometry, not a
  defect — do not "fix" it.

## Traps

- **`gradi-loc` must stay separate from `gradi`.** `fair-esm` vs EvolutionaryScale `esm`, same
  top-level package name. This is the single most expensive mistake available in this stage.
- **Load-bearing pins:** `setuptools<81` (DeepLocPro imports `pkg_resources`, removed in newer
  setuptools) and `transformers==4.44.2` (TMbed loads ProtT5 through `T5Tokenizer`; transformers 5.x
  routes that through the tiktoken converter and dies with *"`tiktoken` is required to read a
  `tiktoken` file"*).
- **Do not extrapolate a TMbed ETA from a contended run.** v1 measured ~17 min/shard while
  DeepLocPro was competing for the CPU, against ~4 min/shard with the machine to itself. This stage
  therefore runs the two tracks **sequentially**, never concurrently.
- **DeepLocPro was trained on UniProt 2023_03 + PSORTdb 4.0**, so *E. coli* K-12 is almost certainly
  in its training set. Any *E. coli* agreement number is a **sanity check, not independent
  validation**. Kp HS11286 and Sa NCTC 8325 are the honest test sets. v1 measured 84.1% agreement
  with UniProt-curated Ec labels and correctly refused to call it validation.
- **DeepLocPro is least decisive on `extracellular`, and TMbed says so independently.** v1 measured
  its within-model probability profile as diagonal at 0.90 cytoplasm / 0.93 inner membrane / 0.86
  periplasm / 0.82 outer membrane but only **0.63 extracellular**, leaking 0.17 back to cytoplasm.
  The cross-check corroborates it from outside the model: Kp's `extracellular` proteins average
  `cytoplasmic_fraction` **0.673** with only **21%** carrying a signal peptide — i.e. TMbed reads a
  large share of them as cytoplasmic, and a protein cannot be secreted with no export signal and
  every residue inside. Treat this class with the most caution; it is where the two predictors
  disagree most. With the Gram-positive remap folding two more classes into it, it is weakest of all
  on *S. aureus*.
- **`cell_wall_surface` is never predicted on a Gram-negative** — zero rows for Kp and Ec, in v1
  and again here. It is not a broken class: *S. aureus* gets **19** calls, `spa` (LPXTG-anchored
  protein A) among them. A zero in that column is a fact about Gram-negatives.
- **Klebsiella's OmpC/OmpF are OmpK36/OmpK35** and are not named `ompC`/`ompF` in HS11286, so those
  two markers are legitimately skipped for Kp. `phoE` stands in as a named Kp outer-membrane porin.
- **The same gene symbol can name a different protein in another species**, and a hand-written
  marker panel is exactly where that bites. The first full run failed its β-barrel check on
  `saureus/fhuA`: in *E. coli* `fhuA` is the 22-strand *"Ferrichrome outer membrane transporter"*,
  while in *S. aureus* it is *"Iron compound ABC transporter, ATP-binding protein"* — cytoplasmic,
  0 strands, and TMbed was right. `CLAUDE.md`'s *"join on locus_tag, not on gene name"* applies to
  spot-check panels too. The β-barrel panel is now **Gram-negative only**, and a Gram-positive gets
  the opposite invariant instead: **zero** β-barrels, because it has no outer membrane to hold one.
  That turns a check that could not apply into one that means something.
- **DeepLocPro's GitHub `LICENSE` is CC BY-NC-SA 4.0 (non-commercial)** while the paper states
  CC BY 4.0. Worth confirming with the GraDi partners before this reaches a commercial deliverable.

## Run log

### 2026-09-01 — first full run

DeepLocPro 1.0.0 · tmbed 1.0.2 · torch 2.13.0 · transformers 4.44.2 · MPS for DeepLocPro, CPU for
TMbed · **263 min wall clock** (DeepLocPro ~38 min at 5.5–6.0 proteins/s; TMbed ~225 min over 53
shards at a mean 4.2 min/shard). All 13,020 proteins, both predictors, **zero failures**.

| species | n | Cyt | CM | Peri | OM | CW | Ext |
|---|---|---|---|---|---|---|---|
| *K. pneumoniae* HS11286 | 5,728 | 3,463 | 1,434 | 339 | 232 | 0 | 260 |
| *E. coli* K-12 | 4,403 | 2,601 | 1,155 | 253 | 206 | 0 | 188 |
| *S. aureus* NCTC 8325 | 2,889 | 1,785 | 810 | — | — | 19 | 275 |

Mean confidence 0.891 / 0.909 / 0.888. Truncated at 2,000 aa: **2 Kp, 1 Ec, 6 Sa**.
β-barrels (≥ 8 strands): **67 Kp, 66 Ec, 0 Sa**. Signal peptides: 604 / 564 / 211.

### It reproduces v1 exactly

Every DeepLocPro count above is **identical** to v1's `09c`, and both β-barrel counts are identical
to v1's `09d` (67 + 66 = the 133 v1 reported). So are the truncation counts (2 Kp, 1 Ec) and the
mean TM helices in `cytoplasmic_membrane` (5.2). That is worth stating because the driver was
rewritten from scratch and TMbed's output format changed from 1 to 4 — same models, same proteomes,
same answers, through different code.

### Validation

- **Marker panel: 47/47 = 100%**, across all three species, against an 80% floor. Every outer-membrane
  porin, periplasmic binding protein, membrane transporter and cytoplasmic chaperone landed in the
  right compartment — including *S. aureus*'s `spa` → `cell_wall_surface` and `hly` → `extracellular`.
- **β-barrel strand counts are textbook**, on both Gram-negatives: OmpA 8, OmpC/OmpF 16, LamB 18,
  BtuB 22, FhuA 22, LptD 26 — and TolC at 4, correctly *not* flagged.
- **The two predictors separate cleanly**, which is the real cross-check since they share no
  machinery. Kp:

  | compartment | n | cyto_frac | TM helices | strands | signal peptide |
  |---|---|---|---|---|---|
  | `cytoplasm` | 3,463 | 0.992 | 0.0 | 0.0 | 0.1% |
  | `cytoplasmic_membrane` | 1,434 | 0.445 | 5.2 | 0.0 | 4.1% |
  | `periplasm` | 339 | 0.064 | 0.0 | 0.0 | 89.4% |
  | `outer_membrane` | 232 | 0.239 | 0.1 | 5.2 | 78.9% |
  | `extracellular` | 260 | 0.673 | 0.1 | 0.0 | 20.8% |

  Helices only in the membrane, strands only in the outer membrane, signal peptides only where
  something was exported. *S. aureus*'s 19 `cell_wall_surface` proteins come out at cyto_frac 0.041
  with **84.2%** carrying a signal peptide — textbook for Sec-exported LPXTG-anchored proteins.

### `extracellular` is the weakest class, and TMbed says so from outside the model

Kp's `extracellular` calls average `cytoplasmic_fraction` **0.673** with only **20.8%** carrying a
signal peptide. A protein cannot be secreted with no export signal and every residue on the
cytoplasmic side, so TMbed is disputing a large share of them. This corroborates, independently,
the weakness v1 measured *inside* the model — DeepLocPro's own probability profile is 0.82–0.93 on
every other class but **0.63** on extracellular, leaking 0.17 back to cytoplasm. Ec is somewhat
better (0.503 / 41.5%), Sa in between (0.550 / 34.9%).

**Treat `extracellular` as the least trustworthy compartment**, and prefer
`has_signal_peptide == False & cytoplasmic_fraction > 0.5` as a flag for "probably actually
cytoplasmic".

**How many that flag catches, measured** (`scripts/plots/localization.py`): **419 of 723
(58.0%)** `extracellular` calls across the three species — Kp **177/260 (68.1%)**, Ec 91/188
(48.4%), Sa 151/275 (54.9%). So on the anchor organism roughly **two thirds** of the class is
disputed by TMbed, which is a much stronger statement than "a large share" and makes
`extracellular` unusable as a filter on Kp without the flag applied.

The violin in `cross_check.png` panel A shows why a summary statistic misleads here: the
distribution is **bimodal**, with one mass at `cytoplasmic_fraction` ≈ 0 (genuinely exported) and
another at ≈ 1.0 (not exported at all). The class mean of 0.582 sits in the empty valley between
them and describes almost no protein in it — and the *median*, at 0.96, is worse still. Do not
reduce this class to one number.

### The Gram-positive remap, quantified

The remap moved a mean of **0.033** probability mass into `extracellular` — small on average, but
above 0.1 for **210** of 2,889 proteins and as high as 0.987 for one. Its effect on the labels:
**32 of the 275** `extracellular` calls belong to proteins whose *un-remapped* argmax was
`periplasm` or `outer_membrane`.

**Corrected 2026-09-02** (`scripts/plots/localization.py`, figure 3): an earlier version of this
section concluded "275 with the remap against 243 without — roughly 12% inflation". That
understated it. The un-remapped `extracellular` argmax count is **226**, so the inflation is
**226 → 275, +49 calls (+21.7%)**, and **49 labels change** in total, not 32.

The 32 was arrived at by subtracting only the flips whose raw argmax was `periplasm` or
`outer_membrane`. But the remap adds that mass *into* `extracellular`, so it can push
`extracellular` past **any** other winner, not just past the two classes it drained. The full
cross-tabulation of raw argmax → shipped label:

| raw argmax | → shipped `extracellular` |
|---|---|
| `outer_membrane` | 31 |
| `cytoplasm` | 13 |
| `cytoplasmic_membrane` | 3 |
| `cell_wall_surface` | 1 |
| `periplasm` | 1 |
| **total flipped** | **49** |

The 13 cytoplasm→extracellular flips are the ones the original framing could not see, and they are
the least comfortable of the set: a protein the model called cytoplasmic outright became
`extracellular` on the strength of Gram-negative-only probability mass. Reach for
`evidence/deeplocpro_probabilities_saureus.tsv` before trusting a borderline *S. aureus*
`extracellular` call.

Reassuringly, DeepLocPro largely already knows *S. aureus* is not Gram-negative: its raw argmax
distribution is 1,798 cytoplasm / 813 cytoplasmic_membrane / 226 extracellular / 31 outer_membrane /
20 cell_wall_surface / **1** periplasm. The remap is a safety net, not a workhorse.

### `cell_wall_surface` is a real class, on the right organism

Zero calls on both Gram-negatives — in v1 and again here. **19 on *S. aureus***, `spa` among them.
A zero in that column is a fact about Gram-negatives, not a broken class.

### One check failed, and it was the check that was wrong

The first run exited non-zero on `saureus/fhuA`: 0 strands where the panel expected ≥ 8. TMbed was
right and the panel was wrong — see the gene-symbol trap above. Fixed by restricting the β-barrel
panel to Gram-negatives and giving Gram-positives the *zero-barrels* invariant instead. The re-run
was **0.1 min**: every prediction came from cache, which is what the caches are for.

### Residual gaps worth knowing

- **No independent ground truth was used.** UniProt is deliberately not consulted, so the numbers
  above are internal consistency plus a marker panel, not a measured accuracy. The marker panel is
  20 well-known proteins, not a benchmark. If an accuracy figure is ever needed, UniProt's
  ECO-gated experimental localizations are the set to score against — 803 for E. coli, but **1** for
  Kp HS11286, and *E. coli* is almost certainly in DeepLocPro's training data.
- **`cytoplasmic_membrane` hides a real distinction.** Kp has proteins in it at
  `cytoplasmic_fraction = 1.0` with **zero** TM helices — peripheral, membrane-associated rather
  than membrane-spanning, and fully Clp-reachable. v1 counted 231 such proteins. Anything that
  treats `cytoplasmic_membrane` as uniformly half-accessible will be wrong about them; that is what
  `cytoplasmic_fraction` and `n_tm_helix` are for.
- **No figures yet.** `scripts/plots/localization.py` does not exist. v1's `09h`–`09l` are prior
  art but four of five panels read `clp_accessibility`, which this stage deliberately does not
  produce.

### 2026-09-02 — figures, and two numbers they corrected

`scripts/plots/localization.py` added. Nothing was re-predicted; the three figures read the
stage's own outputs, including the un-remapped probability vectors. Two claims in this document did
not survive being plotted:

- **The Gram-positive remap is bigger than reported.** 49 labels change, not 32, and the un-remapped
  `extracellular` count is 226, not 243 — so the inflation is +21.7%, not ~12%. The original figure
  counted only `periplasm`/`outer_membrane`-origin flips and missed 17 that came from elsewhere,
  including **13 from `cytoplasm`**. Corrected in place above, with the full cross-tabulation.
- **"TMbed disputes a large share of `extracellular`" is 58.0%** — 419 of 723 across the three
  species, and **68.1% on Kp**. Worth having as a number, because two thirds is a different
  proposition from "a large share" when deciding whether the class is usable as a filter.

Both were arithmetic-on-the-outputs findings, available since the first run and simply never
computed. Plotting a class forces you to look at its distribution, which is how the bimodality — and
therefore the uselessness of its mean *and* its median — became obvious.

Everything else reproduced exactly: per-compartment `cytoplasmic_fraction`, signal-peptide share, TM
helix and strand means all match the first-run table to three decimals, and the composition
percentages match `evidence/deeplocpro_counts_<species>.tsv`.
