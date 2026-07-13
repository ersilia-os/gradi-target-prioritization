# Gr-ADI: Target prioritization for *K. pneumoniae* and *E. coli*

Computational target selection for the Gr-ADI project *"Exploring BacPROTACs as a new paradigm for antibacterial discovery"*, led by **Prof. Erick Strauss** (Stellenbosch University).

This repository covers Ersilia's contribution to **target selection**: building a workflow that prioritises proteins of interest (PoI) in *Klebsiella pneumoniae* and *Escherichia coli* as candidates for targeted protein degradation (BacPROTACs).

## Background

This project explores targeted protein degradation (TPD) as a new modality for antibacterial discovery against Gram-negative pathogens.

Ersilia leads the computational target selection workflow, applying integrative chemo- and bioinformatic approaches to nominate degradable, druggable, and biologically meaningful targets in *K. pneumoniae* and *E. coli*. Ligand identification against the prioritized PoIs is handled in a separate workstream and is **out of scope** here.

**Deliverable.** Prioritized list of PoIs for *K. pneumoniae* and *E. coli*.

## Meetings

- 26/05/14: [Meeting #1](https://docs.google.com/presentation/d/1ktqv42ylLPgQo6vBqlrP5tt2mJl2Mrk_qCTjA0cztTs/edit?usp=drivesdk). Kick-off meeting in which workflow diagrams are discussed.
- 26/06/12: [Meeting #2](https://docs.google.com/presentation/d/18RxzTKev5Cop0QIokVffumbeuKxct-54t6emjNE2n2A/edit?usp=sharing). Task-agnostic annotation of the *Klebsiella pneumoniae* and *Escherichia coli* proteomes.
- 26/06/26: [Meeting #3](https://docs.google.com/presentation/d/1_w6N2veARYSRlDvryVdt-O0DD93AiV-iSKjaDv00N68/edit?usp=sharing). Ligandability assessment.

## Prioritization criteria

Targets are scored along the axes documented in [`docs/`](docs/). See [`docs/pipeline.md`](docs/pipeline.md) for the index and diagram-style legend, and each section below for the full track breakdown, Mermaid diagram, and references.

1. **[Task-agnostic per-protein annotation](docs/01_task_agnostic.md)**. Proteome-wide evidence consumed by every downstream axis: PANTHER / InterPro family-and-domain classification, PDB + AlphaFold structural quality, BV-BRC conservation (within-Kp, cross-species, vs human), bibliometric novelty (UniProt + Europe PMC), and ESM-2 embeddings.
2. **[Ligandability](docs/02_ligandability.md)**. Can a small-molecule recruiter engage the target? Combines binding-affinity transfer from bacterial orthologs (OrthoDB + ChEMBL / BindingDB), structural ligand evidence (PDB co-crystals + AlphaFill) and pocket / binding-site prediction (AF2Bind + fpocket + P2Rank), with disorder as a negative signal.
3. **[Degradability](docs/03_degradability.md)**. Is the target a plausible Clp-protease substrate? Combines sequence-based degron motifs (modulated by AlphaFold pLDDT exposure), *E. coli* substrate-trap and half-life transfer (Flynn 2003 + Nagar 2021), cross-bacterial trap evidence (*Caulobacter*, *M. tuberculosis* ClpC1, *S. aureus*) and an ESM-2-based learned classifier.
4. **[Essentiality](docs/04_essentiality.md)**. Is the target required for fitness or survival? Combines direct *K. pneumoniae* Tn-seq (in vitro / in vivo) and CRISPRi with *E. coli* essentiality transfer and computational predictors (ProteomeLM-Ess, Geptop 2.0, DeeplyEssential, FBA on iYL1228).
5. **[Expression and localization](docs/05_expression_and_localization.md)**. Is the target expressed and accessibly localised? Combines UniProt / PSORTb / DeepLocPro subcellular calls with Kp proteomics abundance and cross-species expression transfer from *S. aureus* (ADEP4 / ONC212) and *E. coli*.

## About Ersilia

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South.

![Ersilia Logo](assets/Ersilia_Brand.png)
