# Ligandability assessment

Combines direct and ortholog-mediated ligand evidence (ChEMBL + BindingDB), structural ligand evidence (PDB co-crystals + AlphaFill) and pocket / binding-site prediction (AF2Bind, fpocket, P2Rank) to score whether a small-molecule recruiter could engage the target.

```mermaid
%%{init: {'theme':'base','themeVariables':{'primaryColor':'#FAD782','primaryBorderColor':'#50285A','primaryTextColor':'#50285A','lineColor':'#50285A','secondaryColor':'#8CC8FA','tertiaryColor':'#BEE6B4','clusterBkg':'#F0F0EE','clusterBorder':'#B0B0AE','titleColor':'#50285A','fontFamily':'Inter, system-ui, sans-serif'}}}%%
flowchart LR
    classDef source    fill:#AA96FA,stroke:#50285A,stroke-width:1.5px,color:#1F0F2E
    classDef dataset   fill:#FAD782,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef method    fill:#8CC8FA,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef embedding fill:#AA96FA,stroke:#50285A,stroke-width:1.5px,color:#1F0F2E
    classDef result    fill:#BEE6B4,stroke:#50285A,stroke-width:2px,color:#50285A
    classDef tagnostic fill:#DCA0DC,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef stub      fill:#FAA08C,stroke:#50285A,stroke-width:1.5px,stroke-dasharray:6 3,color:#50285A
    classDef planned   fill:#D2D2D0,stroke:#7A7A78,stroke-width:1px,stroke-dasharray:5 5,color:#5A5A58

    P["<i>Klebsiella pneumoniae</i> proteome"]:::source

    STR("Structures (PDB + AlphaFold)"):::tagnostic

    subgraph BIOACT [" Binding affinity measurements "]
        direction LR
        ORTH["<b>2.1a</b> · OrthoDB ortholog expansion"]:::method
        LIG["<b>2.1b</b> · ChEMBL + BindingDB<br/><sub>bioactivity</sub>"]:::method
    end

    subgraph STRLIG [" Structural ligand evidence "]
        direction LR
        COCR["<b>2.2a</b> · PDB co-crystals"]:::method
        AFILL["<b>2.2b</b> · AlphaFill ligand transfer"]:::method
    end

    subgraph POCKETS [" Pocket / binding-site prediction "]
        direction LR
        BIND["<b>2.3a</b> · AF2Bind<br/><sub>binding-site prediction</sub>"]:::method
        POCK["<b>2.3b</b> · fpocket / P2Rank<br/><sub>pocket detection</sub>"]:::method
    end

    P --> ORTH
    ORTH --> LIG
    P --> LIG

    P --> STR

    STR --> COCR
    STR --> AFILL

    STR --> BIND
    STR --> POCK

    STR --> DIS["<b>2.4</b> · Disorder filter"]:::method

    LIG   --> SCORE(["Ligandability annotation"]):::result
    COCR  --> SCORE
    AFILL --> SCORE
    BIND  --> SCORE
    POCK  --> SCORE
    DIS   --> SCORE
```

## Tracks

| ID | Title | Description | Resources |
| --- | --- | --- | --- |
| 2.1a | OrthoDB ortholog expansion | Fan each Kp protein into a bacterial-wide ortholog set so sparse direct activity data on HS11286 can be lifted from related bacteria. | OrthoDB |
| 2.1b | ChEMBL + BindingDB bioactivity | Ki / Kd / IC50 records pulled both directly for the Kp protein and (more importantly) across its bacterial ortholog set, drawn from both ChEMBL and BindingDB. | ChEMBL, BindingDB |
| 2.2a | PDB co-crystals | Scan PDB for protein-ligand co-crystal chains covering the Kp protein or its orthologs — direct empirical evidence of binding. | RCSB PDB |
| 2.2b | AlphaFill ligand transfer | Plausible ligands grafted onto the AlphaFold model from homologous PDB entries. | AlphaFill |
| 2.3a | AF2Bind binding-site prediction | Per-residue ligand-binding predictions from AlphaFold2's pair representation. | AF2Bind |
| 2.3b | fpocket / P2Rank pocket detection | Pocket detection on PDB / AlphaFold structures — volume, hydrophobicity, druggability scores. | fpocket, P2Rank |
| 2.4 | Disorder filter | Penalises proteins with a high fraction of low-pLDDT residues (proxy for intrinsic disorder). | AlphaFold DB |

## Key resources

| Resource | Description | Tracks |
| --- | --- | --- |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | Manually-curated bioactivity database (Ki / Kd / IC50, pChEMBL). | 2.1b |
| [BindingDB](https://www.bindingdb.org/) | Public protein–ligand binding-affinity database; complements ChEMBL with additional measured activities. | 2.1b |
| [AlphaFill](https://alphafill.eu/) | Companion to AlphaFold DB that transplants ligands and cofactors from homologous PDB entries onto AlphaFold models. | 2.2b |
| [AF2Bind](https://github.com/sokrypton/af2bind) | Predicts small-molecule binding residues using AlphaFold2's pair representation. | 2.3a |
| [fpocket](https://github.com/Discngine/fpocket) | Fast open-source pocket-detection algorithm based on Voronoi tessellation. | 2.3b |
| [P2Rank](https://github.com/rdk/p2rank) | Machine-learning predictor of ligand-binding sites from protein structure. | 2.3b |

## Suggestions

- **[PocketMiner](https://www.nature.com/articles/s41467-023-36699-3)** — cryptic / induced pocket prediction (GVP-GNN); new sub-track 2.3c. *Bacterial caveat:* eukaryote-leaning training — sanity-check on bacterial apo/holo pairs first.
- **[CryptoBank PLM](https://www.science.org/doi/10.1126/sciadv.ady6364)** ([cryptobankdb.com](https://cryptobankdb.com/)) — sequence-based cryptic prior; PLM head is organism-agnostic but training is eukaryote-dominated. Use as secondary signal.
- **[FTMap](https://pmc.ncbi.nlm.nih.gov/articles/PMC4762777/)** — computational solvent / hot-spot mapping for ternary-complex landing pads; new track 2.5. Organism-agnostic.
- **[PASSer](https://academic.oup.com/nar/article/51/W1/W427/7145694)** — allosteric site prediction; orthogonal BacPROTAC handle. *Bacterial caveat:* training is eukaryote-heavy — positive signal only, no penalty for absence.
- **[canSAR](https://academic.oup.com/nar/article/53/D1/D1287/7899530) + [DoGSiteScorer](https://www.zbh.uni-hamburg.de/en/forschung/amd/software/dogsitescorer.html)** — attach a calibrated druggability number to §2.3b pockets; canSAR 2024 covers AlphaFold across organisms (confirm Kp coverage on first use).
- **§2.2a PDB co-crystals → [BioLiP2](https://academic.oup.com/nar/article/52/D1/D404/7233921)** — pre-curated (biological-unit-aware, filters artefacts, weekly sync).
- **§2.4 disorder filter: keep pLDDT-fraction or swap to [AIUPred](https://academic.oup.com/nar/article/52/W1/W176/7673484)** for one calibrated number.
