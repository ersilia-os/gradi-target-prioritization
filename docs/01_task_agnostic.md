# Task-agnostic per-protein annotation

Part 1 of the GraDi target-prioritization pipeline. See
[`pipeline.md`](./pipeline.md) for the index and the diagram style legend.

This layer produces per-protein evidence that is independent of the downstream
prioritization axes. Each track below runs once per reference proteome and
writes a TSV under `data/processed/` keyed by UniProt accession (or, for
BV-BRC conservation, by locus tag); the five outputs are joined to form the
task-agnostic annotation table that all task-specific scorers consume.

```mermaid
%%{init: {'theme':'base','themeVariables':{'primaryColor':'#FAD782','primaryBorderColor':'#50285A','primaryTextColor':'#50285A','lineColor':'#50285A','secondaryColor':'#8CC8FA','tertiaryColor':'#BEE6B4','fontFamily':'Inter, system-ui, sans-serif'}}}%%
flowchart TD
    classDef source    fill:#AA96FA,stroke:#50285A,stroke-width:1.5px,color:#1F0F2E
    classDef dataset   fill:#FAD782,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef method    fill:#8CC8FA,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef result    fill:#BEE6B4,stroke:#50285A,stroke-width:2px,color:#50285A,font-weight:bold
    classDef tagnostic fill:#DCA0DC,stroke:#50285A,stroke-width:1.5px,color:#50285A
    classDef stub      fill:#FAA08C,stroke:#50285A,stroke-width:1.5px,stroke-dasharray:6 3,color:#50285A
    classDef planned   fill:#D2D2D0,stroke:#7A7A78,stroke-width:1px,stroke-dasharray:5 5,color:#5A5A58

    SRC["UniProt reference proteome<br/>UP000007841 — <i>K. pneumoniae</i> HS11286<br/>5,728 proteins<br/><sub>columns: accession · gene_names · sequence</sub>"]:::source

    SRC --> ESM2["ESM2 embeddings<br/><sub>per-protein 1280-d vector</sub><br/><sub>(planned)</sub>"]:::planned
    SRC --> STR["Structural coverage<br/>PDB (PDBe SIFTS) + AlphaFold pLDDT<br/><sub>scripts/02_structural_coverage.py</sub>"]:::method
    SRC --> PAN["PANTHER family / subfamily<br/><sub>UniProt xref_panther</sub><br/><sub>scripts/01_annotate_panther.py</sub>"]:::method
    SRC --> INT["InterPro domains<br/><sub>UniProt xref_interpro / InterProScan</sub><br/><sub>(planned)</sub>"]:::planned
    SRC --> CONS["Cross-strain conservation<br/>BV-BRC PATtyFams (PLFam / PGFam)<br/><sub>data/raw/bvbrc/hs11286_features.tsv</sub><br/><sub>src/conservation.py</sub>"]:::method

    ESM2 --> T["Task-agnostic per-protein annotation table<br/><sub>joined by UniProt accession (locus_tag for BV-BRC)</sub>"]:::result
    STR  --> T
    PAN  --> T
    INT  --> T
    CONS --> T

    T -.-> NEXT["→ task-specific layer:<br/>ligandability · degradability · essentiality<br/><sub>(covered in the other sections)</sub>"]:::planned
```

## Tracks

| Track | Input | Resource | Script | Output |
| --- | --- | --- | --- | --- |
| ESM2 embeddings | sequence | ESM2-650M (1280-d) | _planned_ | _planned_ |
| Structural coverage | accession | PDBe SIFTS, AlphaFold DB | `scripts/02_structural_coverage.py` | `data/processed/<slug>_structural_coverage.tsv` |
| PANTHER family / subfamily | UniProt xref | PANTHER HMM library | `scripts/01_annotate_panther.py` | `data/processed/<slug>_panther.tsv` |
| InterPro domains | UniProt xref / sequence | InterPro / InterProScan | _planned_ | _planned_ |
| Cross-strain conservation | locus_tag | BV-BRC PATtyFams (PLFam / PGFam) | `src/conservation.py` | `plfam_id`, `pgfam_id`, `has_plfam` in the assembled table |

The reference proteome itself is produced by `scripts/00_download_proteome.py`
(UniProt stream API → `data/raw/<slug>_proteome.tsv`). BV-BRC conservation is
listed here because it is per-protein and task-agnostic; downstream sections
(notably [essentiality](./04_essentiality.md)) treat it as a confidence
modifier rather than a primary signal.

---

**Next:** [Ligandability assessment](./02_ligandability.md) ·
[Degradability assessment](./03_degradability.md) ·
[Essentiality / vulnerability assessment](./04_essentiality.md)
