# Gr-ADI — Target prioritization for *K. pneumoniae* and *E. coli*

Computational target selection for the Gr-ADI project *"Exploring BacPROTACs as a new paradigm for antibacterial discovery"*, led by Prof. Erick Strauss (Stellenbosch University). This repository covers Ersilia's contribution to **WP1.1 — Target selection**: building a workflow that prioritises proteins of interest (PoI) in *Klebsiella pneumoniae* and *Escherichia coli* as candidates for targeted protein degradation (BacPROTACs).

## Background

Gr-ADI is a consortium-level effort to explore targeted protein degradation (TPD) as a new modality for antibacterial discovery against Gram-negative pathogens. Ersilia leads the computational target selection workflow, applying integrative chemo- and bioinformatic approaches to nominate degradable, druggable, and biologically meaningful targets in *K. pneumoniae* and *E. coli*. Ligand identification against the prioritized PoIs is handled in a separate workstream and is **out of scope** here.

**Deliverable.** A one-off, prioritized list of PoIs for *K. pneumoniae* and *E. coli*.

## Meetings

- 26/05/14: [Meeting #1](https://docs.google.com/presentation/d/1ktqv42ylLPgQo6vBqlrP5tt2mJl2Mrk_qCTjA0cztTs/edit?usp=drivesdk). Kick-off meeting.

## Prioritization criteria

Targets are scored and ranked along complementary axes:

- **Essentiality / vulnerability** — is the target required for fitness or survival?
- **Degradability** — structural and sequence features compatible with BacPROTAC-mediated degradation (e.g. engagement by ClpP/ClpC/ClpX machinery).
- **Ligandability** — presence of tractable pockets for small-molecule recruiters.
- **Novelty** — degree to which the target is unexplored relative to known antibacterial space.
- **Expression and localization** — proteomics-based evidence (initial datasets from *S. aureus* treated with ADEP4 / ONC212; ideally extended to *E. coli* and *K. pneumoniae* through the consortium).

## Getting started

```bash
git clone https://github.com/ersilia-os/gradi-target-prioritization
cd gradi-target-prioritization
bash install.sh
```

Code is tracked in Git; `data/` and `output/` are tracked via [eosvc](https://github.com/ersilia-os/eosvc) (DVC + S3).

## About Ersilia

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South.

![Ersilia Logo](assets/Ersilia_Brand.png)
