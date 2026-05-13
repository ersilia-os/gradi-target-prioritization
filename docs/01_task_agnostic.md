# Task-agnostic per-protein annotation

This layer produces general-purpose protein annotation not directly related to the PoI requirements.

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

    SRC["<b>1.0</b> · <i>K. pneumoniae</i> HS11286<br/>proteome (5,728 proteins)"]:::source
    SEQ["<i>sequence</i>"]:::dataset
    ID["<i>identifier</i>"]:::dataset
    SRC --> SEQ
    SRC --> ID

    subgraph FAMDOM [" Family &amp; domain annotation "]
        direction LR
        PAN["<b>1.1a</b> · PANTHER family / subfamily"]:::method
        INT["<b>1.1b</b> · InterPro domains"]:::method
    end
    ID --> PAN
    ID --> INT

    subgraph STRUC [" Structural annotation "]
        direction LR
        PDB["<b>1.2a</b> · PDB coverage"]:::method
        AF["<b>1.2b</b> · AlphaFold pLDDT"]:::method
    end
    ID --> PDB
    ID --> AF

    subgraph CONSERVE [" Conservation "]
        direction LR
        CONS_IDS["<b>1.3a</b> · BV-BRC protein families<br/><sub>PLFam (within-Kp)<br/>PGFam (cross-species)</sub>"]:::method
        CONS_KP["<b>1.3b</b> · Within-Kp pan-genome class"]:::method
        CONS_XS["<b>1.3c</b> · Cross-species broad-spectrum"]:::method
        CONS_SEL["<b>1.3d</b> · Selectivity vs human"]:::method
    end
    ID --> CONS_IDS
    ID --> CONS_KP
    SEQ --> CONS_KP
    ID --> CONS_XS
    SEQ --> CONS_XS
    ID --> CONS_SEL
    SEQ --> CONS_SEL

    ID --> POP["<b>1.4</b> · Bibliometric / popularity"]:::method

    SEQ --> ESM2(["<b>1.5</b> · ESM2 embeddings"]):::embedding

    PAN  --> T(["Task-agnostic annotation"]):::result
    INT  --> T
    PDB  --> T
    AF   --> T
    CONS_IDS --> T
    CONS_KP  --> T
    CONS_XS  --> T
    CONS_SEL --> T
    POP  --> T
```

## Tracks

| ID | Title | Description | Resources |
| --- | --- | --- | --- |
| 1.0 | Reference proteome | The *K. pneumoniae* HS11286 proteome (5,728 proteins) — anchor that every downstream track derives from. | UniProt |
| 1.1a | PANTHER family / subfamily | Functional protein-family classification from PANTHER HMMs. | UniProt, PANTHER |
| 1.1b | InterPro domains | Domain-composition annotation. | UniProt, InterPro |
| 1.2a | PDB coverage | Fraction of residues covered by experimentally-resolved PDB chains. | PDB, PDBe SIFTS |
| 1.2b | AlphaFold pLDDT | Predicted-structure confidence summarised across the protein (high / confident / low residue fractions). | AlphaFold DB |
| 1.3a | BV-BRC protein families | PLFam (within-Kp) and PGFam (global) cluster IDs from the BV-BRC PATtyFam pan-genome system. | BV-BRC |
| 1.3b | Within-Kp pan-genome class | Is the gene core / soft-core / shell / cloud across *K. pneumoniae* strains? | BV-BRC, OrthoFinder / DIAMOND |
| 1.3c | Cross-species broad-spectrum | Phyletic spread across bacterial pathogens (ESKAPE-E) — broad-spectrum signal. | BV-BRC, OrthoDB, BLAST |
| 1.3d | Selectivity vs human | Does a close human ortholog exist? Inverse signal — high similarity is a safety red flag. | UniProt (human), OrthoDB, eggNOG |
| 1.4 | Bibliometric / popularity | How well-studied the protein is, combining UniProt annotation depth and literature counts → tier: dark / studied / well_studied. | UniProt, Europe PMC |
| 1.5 | ESM2 embeddings | Standalone per-protein 1280-d language-model vector — kept separately, not joined into the task-agnostic annotation. | ESM2-650M |

## Key resources

| Resource | Description | Tracks |
| --- | --- | --- |
| [UniProt](https://www.uniprot.org/) | Universal protein sequence and annotation knowledgebase; source of the reference proteome, gene names, cross-references and curation depth. | 1.0, 1.1a, 1.1b, 1.2a, 1.2b, 1.3d, 1.4 |
| [PANTHER](https://www.pantherdb.org/) | Phylogeny-based protein family / subfamily classification built from HMMs. | 1.1a |
| [InterPro](https://www.ebi.ac.uk/interpro/) | Integrated database of protein domains, families and functional sites from member signature databases. | 1.1b |
| [RCSB PDB](https://www.rcsb.org/) | Archive of experimentally-determined 3D macromolecular structures. | 1.2a |
| [PDBe SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) | Residue-level mapping between PDB chains and UniProt sequences. | 1.2a |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/) | EMBL-EBI archive of AlphaFold2-predicted structures with per-residue pLDDT confidence. | 1.2b |
| [BV-BRC](https://www.bv-brc.org/) | Bacterial / viral bioinformatics resource (formerly PATRIC); supplies the PATtyFam protein families (PLFam, PGFam) and pan-genome context. | 1.3a, 1.3b, 1.3c |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | Tool for inferring orthogroups across a set of proteomes. | 1.3b |
| [DIAMOND](https://github.com/bbuchfink/diamond) | Fast accelerated protein sequence aligner, BLAST-compatible; used as the search engine for ortholog detection. | 1.3b, 1.3c |
| [OrthoDB](https://www.orthodb.org/) | Hierarchical catalog of orthologous gene groups across bacteria, eukaryotes and viruses. | 1.3c, 1.3d |
| [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) | Sequence similarity search service against NCBI reference databases. | 1.3c |
| [eggNOG](http://eggnog5.embl.de/) | Orthology resource with hierarchical functional annotation across taxa. | 1.3d |
| [Europe PMC](https://europepmc.org/) | Open literature database (PubMed + preprints + full text) used for publication / mention counts. | 1.4 |
| [ESM-2](https://github.com/facebookresearch/esm) | Meta AI protein language model; the 650M-parameter variant supplies the per-protein 1280-d embeddings. | 1.5 |
