# Pipeline overview

Diagrams describing the GraDi target-prioritization pipeline.

Subsequent sections (essentiality, degradability scoring, and final ranking)
will be added as the corresponding workstreams come online.

**Legend.**
- Solid boxes denote tracks implemented in this repository.
- Dashed boxes denote planned tracks not yet implemented.
- Tan-filled boxes denote inputs pre-computed by the task-agnostic layer
  of Section 1.
- Light-grey boxes denote curated reference data files committed under
  `data/raw/` (e.g. Flynn 2003, Nagar 2021).

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

---

## 2. Ligandability assessment

Ligandability asks whether a small-molecule recruiter could engage the target —
the key tractability question for BacPROTAC discovery. We combine three
positive signals — known ligands transferred from bacterial orthologs (the
lift, since *K. pneumoniae* itself has almost no direct ChEMBL coverage),
structure-based pocket druggability, and a family/fold prior — with one
negative signal — predicted disorder.

```mermaid
flowchart TD
    classDef planned stroke-dasharray: 5 5,stroke:#888,color:#555
    classDef sink fill:#eef,stroke:#446
    classDef tagnostic fill:#fdf6e3,stroke:#b58900,color:#5b4500

    P["<i>K. pneumoniae</i> protein<br/><sub>UniProt accession</sub>"]

    STR["Structures (PDB + AlphaFold)<br/><sub>from task-agnostic layer</sub>"]:::tagnostic
    FAM["PANTHER + InterPro<br/><sub>from task-agnostic layer</sub>"]:::tagnostic

    P --> ORTH["OrthoDB ortholog expansion<br/><sub>bacterial-level groups; pivot via SwissProt</sub><br/><sub>scripts/02_fetch_orthologs.py</sub>"]
    ORTH --> CHEM["ChEMBL bioactivity across orthologs<br/>Ki/Kd/IC50, pChEMBL ≥ 5 / ≥ 6<br/><sub>scripts/03_fetch_chembl_ligands.py</sub>"]

    P --> STR
    P --> FAM

    STR --> POCK["Pocket detection<br/><sub>volume · hydrophobicity · druggability</sub><br/><sub>(planned)</sub>"]:::planned
    STR --> DIS["Disorder filter<br/><sub>AlphaFold pLDDT fractions</sub><br/><sub>(planned, negative signal)</sub>"]:::planned
    FAM --> PRIOR["Family / fold tractability prior<br/><sub>PANTHER + InterPro lookup</sub><br/><sub>(planned)</sub>"]:::planned

    CHEM  --> SCORE["Composite ligandability score → tier<br/><sub>known ligands · pocket · family prior − disorder</sub>"]:::sink
    POCK  --> SCORE
    PRIOR --> SCORE
    DIS  -.->|"penalty"| SCORE

    SCORE -.-> NEXT["→ final target ranking<br/><sub>(covered in subsequent diagram)</sub>"]
```

### Tracks

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

## 3. Degradability assessment (Clp proteases)

Degradability asks how susceptible the target is to bacterial Clp-protease
degradation — directly informative for BacPROTAC design since the recruiter
hands the substrate to ClpC/ClpX/ClpP for proteolysis. There is no
proteome-wide Clp-substrate measurement for *K. pneumoniae*, so the score
combines two evidence streams: rule-based degron motifs computed directly on
each Kp sequence, plus experimental evidence transferred from *E. coli* K-12
by gene-symbol orthology. Note that this layer does **not** consume the
task-agnostic outputs (structures, families, embeddings) — it works directly
on the raw sequence plus two curated reference TSVs.

```mermaid
flowchart TD
    classDef planned stroke-dasharray: 5 5,stroke:#888,color:#555
    classDef sink fill:#eef,stroke:#446
    classDef ref fill:#f5f5f5,stroke:#666,color:#333

    P["<i>K. pneumoniae</i> protein<br/><sub>UniProt accession · sequence · gene_symbol</sub>"]

    P --> CTERM["C-terminal degrons<br/><sub>CM1 ssrA-like (-LAA family), CM2 MuA-like</sub><br/><sub>scripts/03_annotate_clp_degradability.py</sub>"]
    P --> NTERM["N-terminal degrons<br/><sub>N-end rule (L/F/Y/W at pos 2), Flynn NM1/2/3</sub><br/><sub>scripts/03_annotate_clp_degradability.py</sub>"]

    P --> ORTH["Gene-symbol → <i>E. coli</i> K-12 ortholog<br/><sub>direct symbol match (not OrthoDB — Flynn/Nagar evidence is E. coli-specific)</sub><br/><sub>scripts/03_annotate_clp_degradability.py</sub>"]
    ORTH --> FLYNN["Flynn 2003 ClpXP/ClpAP trap census<br/><sub>data/raw/clp_substrates/flynn2003_ecoli_clp_substrates.tsv</sub>"]:::ref
    ORTH --> NAGAR["Nagar 2021 <i>E. coli</i> half-lives<br/><sub>data/raw/clp_substrates/nagar2021_ecoli_halflives.tsv</sub>"]:::ref

    P -.-> CLPK["ClpK paralog handling<br/><sub>Kp-specific heat-shock Clp; no E. coli ortholog</sub><br/><sub>(planned)</sub>"]:::planned
    P -.-> ESM2["ESM2-based degradability ML<br/><sub>per-protein embedding → learned classifier</sub><br/><sub>(planned; see Section 1)</sub>"]:::planned

    CTERM --> SCORE["Composite clp_degradability_score → tier<br/><sub>rule features (≤1.0) + 0.4 if trapped + 0.2/0.1 by t½ class</sub><br/><sub>src/degradability.py → src/assemble.py</sub>"]:::sink
    NTERM --> SCORE
    FLYNN --> SCORE
    NAGAR --> SCORE
    CLPK -.-> SCORE
    ESM2 -.-> SCORE

    SCORE -.-> NEXT["→ final target ranking<br/><sub>(covered in subsequent diagram)</sub>"]
```

### Tracks

| Track | Input | Resource | Script | Output |
| --- | --- | --- | --- | --- |
| C-terminal degrons | sequence | CM1 (ssrA-family) + CM2 (MuA-family) regex | `scripts/03_annotate_clp_degradability.py` | `cterm_ssra_like`, `cterm_mua_like` in `data/processed/klebsiella_pneumoniae_clp_degradability.tsv` |
| N-terminal degrons | sequence | N-end rule + Flynn 2003 NM1/2/3 | `scripts/03_annotate_clp_degradability.py` | `nterm_destabilizing`, `nterm_nm1/2/3` |
| Gene-symbol → *E. coli* ortholog | gene_symbol | local symbol match against `data/raw/escherichia_coli_proteome.tsv` | `scripts/03_annotate_clp_degradability.py` | `ecoli_ortholog_uniprot`, `ortholog_status` |
| Flynn 2003 trap census | E. coli gene_symbol | curated from Flynn 2003 + Sauer/Baker reviews | `scripts/03_annotate_clp_degradability.py` | `ecoli_clp_trapped`, `ecoli_clp_class` |
| Nagar 2021 half-lives | E. coli gene_symbol | curated from Nagar 2021 pulsed-SILAC | `scripts/03_annotate_clp_degradability.py` | `ecoli_halflife_class`, `ecoli_halflife_min` |
| ClpK paralog handling | sequence + Kp-specific HMM | Klebsiella ClpK literature | _planned_ | _planned_ |
| ESM2-based degradability ML | ESM2 embedding | re-train of Nagar 2021's 188-feature classifier | _planned_ | _planned_ |

Two architectural notes worth highlighting:

- **Why gene-symbol matching, not OrthoDB.** The two reference TSVs are
  *E. coli K-12-specific* (Flynn 2003 and Nagar 2021 both used MG1655). A
  bacterial-wide OrthoDB expansion (as used by the ligandability layer) would
  not add evidence — the labels only exist for E. coli. Simple symbol matching
  covers the ~18% of HS11286 entries that carry a canonical gene symbol; the
  rest fall through to the rule-based score alone.
- **Why this layer doesn't consume task-agnostic outputs.** Clp recognition
  is dominated by short linear motifs (terminal degrons) and biophysical
  flexibility. The current v1 captures the motif signal directly from
  sequence; structure-based / ESM2-based extensions are explicitly planned
  (dashed) — they would close the loop back to Section 1's task-agnostic
  layer once built.
