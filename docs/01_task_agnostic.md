# Task-agnostic per-protein annotation

Part 1 of the GraDi target-prioritization pipeline. See
[`pipeline.md`](./pipeline.md) for the index and the diagram style legend.

This layer produces per-protein evidence that is independent of the downstream
prioritization axes. Each track below runs once per reference proteome and
writes a TSV under `data/processed/` keyed by UniProt accession (or, for
BV-BRC-anchored tracks, by locus tag). Nine tracks — structural annotation
(PDB coverage + AlphaFold pLDDT), the family/domain pair (PANTHER +
InterPro), conservation (BV-BRC protein families plus three planned flavors:
within-Kp pan-genome, cross-species broad-spectrum, and selectivity vs
human), and bibliometric popularity — are joined to form the task-agnostic
annotation table that all task-specific scorers consume. ESM2 embeddings are
kept as a separate per-protein artifact (a vector per protein), not as a
column in the joined table.

```mermaid
%%{init: {'theme':'base','themeVariables':{'primaryColor':'#FAD782','primaryBorderColor':'#50285A','primaryTextColor':'#50285A','lineColor':'#50285A','secondaryColor':'#8CC8FA','tertiaryColor':'#BEE6B4','clusterBkg':'#F0F0EE','clusterBorder':'#B0B0AE','titleColor':'#50285A','fontFamily':'Inter, system-ui, sans-serif'}}}%%
flowchart LR
    classDef source    fill:#AA96FA,stroke:#50285A,stroke-width:1.5px,color:#1F0F2E
    classDef dataset   fill:#FAD782,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef method    fill:#8CC8FA,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef embedding fill:#8CC8FA,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef result    fill:#BEE6B4,stroke:#50285A,stroke-width:2px,color:#50285A,font-weight:bold
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
        CONS_IDS["<b>1.3a</b> · BV-BRC protein families<br/><sub>PLFam (within-Kp) · PGFam (cross-species)</sub>"]:::method
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

    PAN  --> T(["Task-agnostic chunk"]):::result
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

| ID | Track | Input | Resource | Script | Output |
| --- | --- | --- | --- | --- | --- |
| **1.0** | *K. pneumoniae* HS11286 reference proteome | proteome ID | UniProt UP000007841 | `scripts/00_download_proteome.py` | `data/raw/<slug>_proteome.tsv` (accession · gene_names · sequence) |
| **1.1a** | PANTHER family / subfamily | UniProt xref | PANTHER HMM library | `scripts/01_annotate_panther.py` | `data/processed/<slug>_panther.tsv` |
| **1.1b** | InterPro domains | UniProt xref | InterPro (via UniProt xref_interpro) | _planned_ | _planned_ |
| **1.2a** | PDB coverage | accession | PDBe SIFTS bulk mapping | `scripts/02_structural_coverage.py` | `pdb_*` columns of `data/processed/<slug>_structural_coverage.tsv` |
| **1.2b** | AlphaFold pLDDT | accession | AlphaFold DB per-prediction API | `scripts/02_structural_coverage.py` | `afdb_*` columns of `data/processed/<slug>_structural_coverage.tsv` |
| **1.3a** | BV-BRC protein families (PLFam / PGFam) | locus_tag | BV-BRC `genome_feature` table — PATtyFam pan-genome clustering | `src/conservation.py` | `plfam_id` (within-Kp), `pgfam_id` (global), `has_plfam` |
| **1.3b** | Within-Kp pan-genome class | PLFam id *or* sequence | BV-BRC `kp_plfam_counts.tsv` (ID path) *or* OrthoFinder / DIAMOND across Kp reference panel (sequence path) | _planned_ | `kp_n_genomes`, `kp_conservation_class` (core / soft-core / shell / cloud) |
| **1.3c** | Cross-species broad-spectrum | PGFam id / OrthoDB group *or* sequence | BV-BRC PGFam counts + OrthoDB sizes (ID path) *or* BLAST against ESKAPE-E panel (sequence path) | _planned_ | `xs_n_species`, `xs_breadth_class` |
| **1.3d** | Selectivity vs human | sequence *or* accession | human proteome (UP000005640) BLAST/DIAMOND (sequence) *or* cross-kingdom OrthoDB / eggNOG (ID) | _planned_ | `human_ortholog_uniprot`, `human_identity_pct` *(inverse signal)* |
| **1.4** | Bibliometric / popularity | accession + gene_symbol | UniProt annotation depth + Europe PMC search | `scripts/02_annotate_popularity.py` | `popularity_tier`: dark / studied / well_studied — `data/processed/<slug>_popularity.tsv` |
| **1.5** | ESM2 embeddings *(standalone)* | sequence | ESM2-650M (1280-d) | _planned_ | _planned (standalone vector store, not joined)_ |

The reference proteome (UniProt **UP000007841**, *K. pneumoniae* HS11286,
5,728 proteins; columns: accession · gene_names · sequence) is produced by
`scripts/00_download_proteome.py` (UniProt stream API →
`data/raw/<slug>_proteome.tsv`). The **task-agnostic chunk** is the result of
joining the nine non-standalone tracks above by UniProt accession (with
`locus_tag` as the join key for BV-BRC-anchored tracks). Conservation is
listed here because it is per-protein and task-agnostic; downstream sections
(notably [essentiality](./04_essentiality.md)) treat it as a confidence
modifier rather than a primary signal. The within-Kp pan-genome class also
acts as a strain-coverage filter, and selectivity-vs-human is an inverse
signal carried into the final ranking as a safety axis.

---

**Next:** [Ligandability assessment](./02_ligandability.md) ·
[Degradability assessment](./03_degradability.md) ·
[Essentiality / vulnerability assessment](./04_essentiality.md)
