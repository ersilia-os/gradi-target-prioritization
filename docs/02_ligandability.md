# Ligandability assessment

Part 2 of the GraDi target-prioritization pipeline. See
[`pipeline.md`](./pipeline.md) for the index and the diagram style legend.

Ligandability asks whether a small-molecule recruiter could engage the target —
the key tractability question for BacPROTAC discovery. We combine three
positive signals — known ligands transferred from bacterial orthologs (the
lift, since *K. pneumoniae* itself has almost no direct ChEMBL coverage),
structure-based pocket druggability, and a family/fold prior — with one
negative signal — predicted disorder.

```mermaid
%%{init: {'theme':'base','themeVariables':{'primaryColor':'#FAD782','primaryBorderColor':'#50285A','primaryTextColor':'#50285A','lineColor':'#50285A','secondaryColor':'#8CC8FA','tertiaryColor':'#BEE6B4','fontFamily':'Inter, system-ui, sans-serif'}}}%%
flowchart LR
    classDef source    fill:#AA96FA,stroke:#50285A,stroke-width:1.5px,color:#1F0F2E
    classDef dataset   fill:#FAD782,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef method    fill:#8CC8FA,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef result    fill:#BEE6B4,stroke:#50285A,stroke-width:2px,color:#50285A,font-weight:bold
    classDef tagnostic fill:#DCA0DC,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef stub      fill:#FAA08C,stroke:#50285A,stroke-width:1.5px,stroke-dasharray:6 3,color:#50285A
    classDef planned   fill:#D2D2D0,stroke:#7A7A78,stroke-width:1px,stroke-dasharray:5 5,color:#5A5A58

    P{{"<i>K. pneumoniae</i> protein<br/><sub>UniProt accession</sub>"}}:::source

    STR("Structures (PDB + AlphaFold)<br/><sub>from task-agnostic layer</sub>"):::tagnostic
    FAM("PANTHER + InterPro<br/><sub>from task-agnostic layer</sub>"):::tagnostic
    CONS("Cross-strain conservation<br/><sub>BV-BRC PATtyFams (from task-agnostic)</sub>"):::tagnostic

    P --> ORTH["OrthoDB ortholog expansion<br/><sub>bacterial-level groups; pivot via SwissProt</sub><br/><sub>scripts/02_fetch_orthologs.py</sub>"]:::method
    ORTH --> CHEM["ChEMBL bioactivity across orthologs<br/>Ki/Kd/IC50, pChEMBL ≥ 5 / ≥ 6<br/><sub>scripts/03_fetch_chembl_ligands.py</sub>"]:::method

    P --> STR
    P --> FAM
    P --> CONS

    STR --> POCK["Pocket detection<br/><sub>volume · hydrophobicity · druggability</sub><br/><sub>(planned)</sub>"]:::planned
    STR --> DIS["Disorder filter<br/><sub>AlphaFold pLDDT fractions</sub><br/><sub>(planned, negative signal)</sub>"]:::planned
    FAM --> PRIOR["Family / fold tractability prior<br/><sub>PANTHER + InterPro lookup</sub><br/><sub>(planned)</sub>"]:::planned

    CHEM  --> SCORE(["Composite ligandability score → tier<br/><sub>known ligands · pocket · family prior − disorder</sub>"]):::result
    POCK  --> SCORE
    PRIOR --> SCORE
    DIS  -.->|"penalty"| SCORE
    CONS -.->|"confidence modifier (planned)"| SCORE

    SCORE -.-> NEXT["→ final target ranking<br/><sub>(covered in subsequent diagram)</sub>"]:::planned
```

## Tracks

| Track | Input | Resource | Script | Output |
| --- | --- | --- | --- | --- |
| OrthoDB ortholog expansion | gene symbol | OrthoDB via UniProt (SwissProt pivot) | `scripts/02_fetch_orthologs.py` | `data/processed/klebsiella_pneumoniae_orthodb_orthologs.tsv` |
| ChEMBL bioactivity | UniProt + orthologs | ChEMBL REST (`/target`, `/activity`) | `scripts/03_fetch_chembl_ligands.py` | `data/processed/chembl_ligand_counts.tsv` |
| Pocket detection | AlphaFold / PDB structure | _tool TBD_ | _planned_ | _planned_ |
| Family / fold prior | PANTHER + InterPro IDs | curated tractability lookup | _planned_ | _planned_ |
| Disorder filter | AlphaFold pLDDT fractions | already in `*_structural_coverage.tsv` | _planned_ | _planned_ |

The ChEMBL track is the most consequential pillar in the implemented portion:
HS11286 is ~99% TrEMBL, so direct ChEMBL coverage is sparse — almost all the
ligand signal arrives through the OrthoDB expansion step that fans each Kp
protein into a bacterial-wide ortholog set.

---

**Prev:** [Task-agnostic per-protein annotation](./01_task_agnostic.md) ·
**Next:** [Degradability assessment](./03_degradability.md) ·
[Essentiality / vulnerability assessment](./04_essentiality.md)
