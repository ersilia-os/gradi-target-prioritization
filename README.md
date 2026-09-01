# Gr-ADI: Target prioritization for *K. pneumoniae*, *E. coli* and *S. aureus*

> ⚠️ **Work in progress.** This repository is an active prototype. Scores, data and methods are
> provisional — do not treat any result here as final.

Computational target selection for the Gr-ADI project *"Exploring BacPROTACs as a new paradigm for
antibacterial discovery"*, led by **Prof. Erick Strauss** (Stellenbosch University).

This repository covers Ersilia's contribution to **target selection**: a workflow that prioritises
proteins of interest (PoI) as candidates for targeted protein degradation. Ligand identification
against the prioritized PoIs is a separate workstream and is **out of scope** here.

## Background

The project explores targeted protein degradation (TPD) as a new modality for antibacterial discovery.
Ersilia leads the computational target-selection workflow, applying integrative chemo- and
bioinformatic approaches to nominate degradable, druggable and biologically meaningful targets.

**Deliverable.** A prioritized list of PoIs.

## Status: v2, a deliberate restart

A first complete version of the prioritization logic was built between May and August 2026, covering
five scoring axes and a web browser for the results. It is **finished and frozen** under
[`legacy/`](legacy/).

v2 restarts the pipeline while keeping what v1 established. If you are picking this up, read
**[`legacy/HISTORY.md`](legacy/HISTORY.md)** first — it records what was built and what was not, the
methodological decisions and their rationale, roughly 70 documented traps, which data sources were
obtained and how, and which artifacts look like data but are not.

The three things v1 most wants passed forward:

1. **The scoring contract** — renormalise a composite over the tracks that actually have evidence
   rather than zero-filling a missing one, and drive tiers off evidence rather than score thresholds.
2. **Sequence-first identifier resolution** — *K. pneumoniae* HS11286 is a dark TrEMBL proteome, and
   matching external databases by accession silently misses even its clinically central
   β-lactamases. Map by sequence.
3. **A negative result worth not re-deriving** — partnerless activated ClpP has no usable sequence
   degron. Only 42 of 4,403 E. coli proteins carry a weighted motif, the composite tracks terminal
   exposure rather than any motif (ρ = 0.999), and nothing computed predicts cleavage: protein length
   alone beats a 13-feature model. A degradability axis needs new *data*, not new features.

### What v2 has so far

| stage | what |
|---|---|
| [`scripts/00_download_proteomes.py`](scripts/00_download_proteomes.py) | the four reference proteomes, and **one identified table per species** — see [`docs/00_proteomes.md`](docs/00_proteomes.md) |

Stage 00 exists mainly to fix the problem that cost v1 the most: the anchor proteomes are badly
under-named. Gene names cover only **18.4%** of *K. pneumoniae* HS11286 and 28.2% of *S. aureus*
NCTC 8325, against 100% for *E. coli* K-12 — and gene names are what literature and databases key on.
Filling them from the rest of each species (identical sequences, then UniRef90 clusters, never across
species) takes *K. pneumoniae* to **63.4%** and *S. aureus* to 44.6%.

Four species, one reference proteome each: *K. pneumoniae* HS11286 · *E. coli* K-12 MG1655 ·
***S. aureus* NCTC 8325** (new in v2 — the organism the activated-ClpP proteomics is native to, and
the only one of the three with the ClpC machinery every published BacPROTAC targets) · human, for
selectivity.

## Repository layout

```
scripts/   the numbered v2 pipeline
src/       shared helpers + proteome_registry.tsv
docs/      one document per stage
legacy/    the complete v1 pipeline, frozen — start at HISTORY.md
data/      inputs (eosvc/S3, not Git); see data/raw/PROVENANCE.md
output/    results and figures (eosvc/S3, not Git)
```

`data/` and `output/` are versioned with [eosvc](https://github.com/ersilia-os/eosvc) (DVC + S3), not
Git. [`data/raw/PROVENANCE.md`](data/raw/PROVENANCE.md) indexes every raw dataset and — importantly —
whether you could obtain it again, since a good deal of it came through authenticated browser
sessions or by hand from paywalled articles.

## The v1 target browser

The prioritization browser built for v1 is still live at
https://ersilia-os.github.io/gradi-target-prioritization/ and is deployed from `legacy/app/`. It
serves v1 numbers, including two degradability columns that `legacy/HISTORY.md` §7 documents as
untrustworthy.

## Meetings

- 26/05/14: [Meeting #1](https://docs.google.com/presentation/d/1ktqv42ylLPgQo6vBqlrP5tt2mJl2Mrk_qCTjA0cztTs/edit?usp=drivesdk). Kick-off; workflow diagrams.
- 26/06/12: [Meeting #2](https://docs.google.com/presentation/d/18RxzTKev5Cop0QIokVffumbeuKxct-54t6emjNE2n2A/edit?usp=sharing). Task-agnostic annotation of the proteomes.
- 26/06/26: [Meeting #3](https://docs.google.com/presentation/d/1_w6N2veARYSRlDvryVdt-O0DD93AiV-iSKjaDv00N68/edit?usp=sharing). Ligandability assessment.
- 26/07/14: [Meeting #4](https://docs.google.com/presentation/d/1gqcBd9pLYAknGwxpmVM3p7lFEUnMwRAYD2RQmBywWlE/edit?usp=sharing). Essentiality annotation and predictions.
- 26/07/24: [Meeting #5](https://docs.google.com/presentation/d/11yDNqMQUHVPZ8q0vn-_kKcUD9v02bo4CdBhp66HPn18/edit?usp=sharing). First draft of the target prioritization browser.

## About Ersilia

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling
sustainable research in the Global South.

![Ersilia Logo](assets/Ersilia_Brand.png)
