# Pipeline overview

Diagrams describing the GraDi target-prioritization pipeline.

Subsequent sections (essentiality / ligandability / degradability scoring, and
final ranking) will be added as the corresponding workstreams come online.

**Legend.** Dashed boxes denote annotations that are planned but not yet
implemented in this repository.

---

## 1. Task-agnostic per-protein annotation

This layer produces per-protein evidence that is independent of the downstream
prioritization axes. Each track below runs once per reference proteome and
writes a TSV under `data/processed/` keyed by UniProt accession; the four
outputs are joined on accession to form the task-agnostic annotation table
that all task-specific scorers consume.

```mermaid
flowchart TD
    classDef planned stroke-dasharray: 5 5,stroke:#888,color:#555
    classDef sink fill:#eef,stroke:#446

    SRC["UniProt reference proteome<br/>UP000007841 — <i>K. pneumoniae</i> HS11286<br/>5,728 proteins<br/><sub>columns: accession · gene_names · sequence</sub>"]

    SRC --> ESM2["ESM2 embeddings<br/><sub>per-protein 1280-d vector</sub><br/><sub>(planned)</sub>"]:::planned
    SRC --> STR["Structural coverage<br/>PDB (PDBe SIFTS) + AlphaFold pLDDT<br/><sub>scripts/02_structural_coverage.py</sub>"]
    SRC --> PAN["PANTHER family / subfamily<br/><sub>UniProt xref_panther</sub><br/><sub>scripts/01_annotate_panther.py</sub>"]
    SRC --> INT["InterPro domains<br/><sub>UniProt xref_interpro / InterProScan</sub><br/><sub>(planned)</sub>"]:::planned

    ESM2 --> T["Task-agnostic per-protein annotation table<br/><sub>joined by UniProt accession</sub>"]:::sink
    STR  --> T
    PAN  --> T
    INT  --> T

    T -.-> NEXT["→ task-specific layer:<br/>essentiality · ligandability · degradability<br/><sub>(covered in subsequent diagrams)</sub>"]
```

### Tracks

| Track | Input | Resource | Script | Output |
| --- | --- | --- | --- | --- |
| ESM2 embeddings | sequence | ESM2-650M (1280-d) | _planned_ | _planned_ |
| Structural coverage | accession | PDBe SIFTS, AlphaFold DB | `scripts/02_structural_coverage.py` | `data/processed/<slug>_structural_coverage.tsv` |
| PANTHER family / subfamily | UniProt xref | PANTHER HMM library | `scripts/01_annotate_panther.py` | `data/processed/<slug>_panther.tsv` |
| InterPro domains | UniProt xref / sequence | InterPro / InterProScan | _planned_ | _planned_ |

The reference proteome itself is produced by `scripts/00_download_proteome.py`
(UniProt stream API → `data/raw/<slug>_proteome.tsv`).
