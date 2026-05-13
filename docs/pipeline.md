# Pipeline overview

Diagrams describing the GraDi target-prioritization pipeline.

A final-ranking section will be added once the corresponding workstream comes
online.

**Legend.** Diagrams use the canonical Mermaid style defined in
[`docs/mermaid_style.md`](./mermaid_style.md), anchored on the Ersilia brand
palette. Seven node classes:

- `:::source` (purple) — anchor proteome / starting entity.
- `:::dataset` (yellow) — curated reference data file (e.g. Flynn 2003, Nagar 2021).
- `:::method` (blue) — script / compute step.
- `:::result` (mint, bold) — output, score, or sink.
- `:::tagnostic` (pink) — pre-computed input from the task-agnostic layer of Section 1.
- `:::stub` (orange, dashed border) — parser exists in the codebase but the data file is not yet staged in `data/raw/`.
- `:::planned` (gray, dashed border) — not yet implemented.

---

## 1. Task-agnostic per-protein annotation

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

    T -.-> NEXT["→ task-specific layer:<br/>essentiality · ligandability · degradability<br/><sub>(covered in Sections 2–4)</sub>"]:::planned
```

### Tracks

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
(notably essentiality) treat it as a confidence modifier rather than a
primary signal.

---

## 2. Ligandability assessment

Ligandability asks whether a small-molecule recruiter could engage the target —
the key tractability question for BacPROTAC discovery. We combine three
positive signals — known ligands transferred from bacterial orthologs (the
lift, since *K. pneumoniae* itself has almost no direct ChEMBL coverage),
structure-based pocket druggability, and a family/fold prior — with one
negative signal — predicted disorder.

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

    P["<i>K. pneumoniae</i> protein<br/><sub>UniProt accession</sub>"]:::source

    STR["Structures (PDB + AlphaFold)<br/><sub>from task-agnostic layer</sub>"]:::tagnostic
    FAM["PANTHER + InterPro<br/><sub>from task-agnostic layer</sub>"]:::tagnostic

    P --> ORTH["OrthoDB ortholog expansion<br/><sub>bacterial-level groups; pivot via SwissProt</sub><br/><sub>scripts/02_fetch_orthologs.py</sub>"]:::method
    ORTH --> CHEM["ChEMBL bioactivity across orthologs<br/>Ki/Kd/IC50, pChEMBL ≥ 5 / ≥ 6<br/><sub>scripts/03_fetch_chembl_ligands.py</sub>"]:::method

    P --> STR
    P --> FAM

    STR --> POCK["Pocket detection<br/><sub>volume · hydrophobicity · druggability</sub><br/><sub>(planned)</sub>"]:::planned
    STR --> DIS["Disorder filter<br/><sub>AlphaFold pLDDT fractions</sub><br/><sub>(planned, negative signal)</sub>"]:::planned
    FAM --> PRIOR["Family / fold tractability prior<br/><sub>PANTHER + InterPro lookup</sub><br/><sub>(planned)</sub>"]:::planned

    CHEM  --> SCORE["Composite ligandability score → tier<br/><sub>known ligands · pocket · family prior − disorder</sub>"]:::result
    POCK  --> SCORE
    PRIOR --> SCORE
    DIS  -.->|"penalty"| SCORE

    SCORE -.-> NEXT["→ final target ranking<br/><sub>(covered in subsequent diagram)</sub>"]:::planned
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

---

## 4. Essentiality / vulnerability assessment

Essentiality asks "is this target required for fitness or survival?";
vulnerability sharpens it by asking "how much depletion is tolerated?" via
graded CRISPRi knockdown. The layer combines two evidence pathways — direct
*K. pneumoniae* measurements (Tn-seq, CRISPRi) and *E. coli* ortholog
inference — into a per-condition consensus call. It is the most data-source-
heavy axis: several papers contribute, with a mix of "data loaded today",
"parser written but file not staged" stubs, and several "explore more"
extensions worth flagging.

```mermaid
flowchart TD
    classDef planned stroke-dasharray: 5 5,stroke:#888,color:#555
    classDef stub fill:#fdf2dc,stroke:#a8772b,color:#5b4500,stroke-dasharray: 3 2
    classDef sink fill:#eef,stroke:#446
    classDef tagnostic fill:#fdf6e3,stroke:#b58900,color:#5b4500

    subgraph DIRECT [" Direct experimental evidence "]
        TNS["Kp Tn-seq screens<br/><sub>Eichelberger 2024 ECL8 ✓ loaded (3 flavors)</sub><br/><sub>Bachman 2015 KPPR1 · Ramage 2017 KPNIH1 — stubs</sub>"]
        CRI["Kp CRISPRi screens<br/><sub>Zhu 2023 Mobile-CRISPRi-seq highlights ✓ loaded (8 genes)</sub>"]
        ECO["E. coli TraDIS<br/><sub>Goodall 2018 BW25113 (stub — feeds Ec-inference path)</sub>"]:::stub
    end

    DEG["DEG 15 / OGEE v2<br/><sub>consolidated essential-gene databases</sub><br/><sub>(planned)</sub>"]:::planned
    FBA["FBA / iYL1228 metabolic model<br/><sub>in silico single-gene knockout</sub><br/><sub>(planned)</sub>"]:::planned
    ML["ML essentiality predictor<br/><sub>DeeplyEssential / ESM2-based</sub><br/><sub>(planned; ties to Section 1 ESM2)</sub>"]:::planned

    TNS --> PARSE["<i>src/essentiality.py</i> · per-paper parsers<br/><sub>emit long-form rows: (gene_symbol, call, score, source, condition, flavor)</sub>"]
    CRI --> PARSE
    ECO --> PARSE
    DEG -.-> PARSE
    FBA -.-> PARSE
    ML  -.-> PARSE

    PARSE --> FLAV["Per-flavor consensus by gene_symbol<br/><sub>5 flavors: in_vitro · in_vivo_lung · in_vivo_urine · in_vivo_serum · vulnerability</sub><br/><sub>CALL_PRIORITY: essential ▶ fitness_defect ▶ unclear ▶ fitness_advantage ▶ non_essential</sub><br/><sub><i>src/assemble.py</i> · _join_flavor_block</sub>"]

    PARSE -->|"E. coli essentials"| ECINF["Kp ↔ Ec gene-symbol match<br/><sub><i>src/assemble.py</i> · _join_ec_inferred</sub>"]
    ECINF -.-> ORTHO_UP["Upgrade: OrthoDB-based Ec transfer<br/><sub>reuse Section 2's OrthoDB pipeline</sub><br/><sub>(planned; lifts Ec coverage beyond ~18%)</sub>"]:::planned

    CONS["Cross-strain conservation<br/><sub>BV-BRC PATtyFams (from Section 1)</sub>"]:::tagnostic

    FLAV  --> OUT["ess_*_call · score · sources<br/><sub>5 flavor triples + ess_Ec_inferred</sub>"]:::sink
    ECINF --> OUT
    ORTHO_UP -.-> OUT
    CONS -.->|"confidence modifier (planned)"| OUT

    OUT -.-> NEXT["→ final target ranking<br/><sub>(covered in subsequent diagram)</sub>"]
```

### Tracks

| Track | Source(s) | Condition / flavor | Status | Output columns |
| --- | --- | --- | --- | --- |
| Kp Tn-seq (in vitro) | Eichelberger 2024 ECL8; Ramage 2017 KPNIH1 (stub) | `in_vitro_essential` | loaded + stub | `ess_in_vitro_*` |
| Kp Tn-seq (in vivo lung) | Bachman 2015 KPPR1 (stub) | `in_vivo_lung` | stub | `ess_in_vivo_lung_*` |
| Kp Tn-seq (in vivo urine) | Eichelberger 2024 | `in_vivo_urine` | loaded | `ess_in_vivo_urine_*` |
| Kp Tn-seq (in vivo serum) | Eichelberger 2024 | `in_vivo_serum` | loaded | `ess_in_vivo_serum_*` |
| Kp CRISPRi vulnerability | Zhu 2023 (highlights) | `vulnerability_crispri` | loaded (8 genes) | `ess_vulnerability_*` |
| E. coli ortholog inference | Goodall 2018 (stub) → `_join_ec_inferred` | `in_vitro_essential_Ec` | stub | `ess_Ec_inferred_call/via/sources` |
| DEG / OGEE consolidated DBs | DEG 15, OGEE v2 | cross-validation | _planned_ | _planned_ |
| FBA in silico knockout | iYL1228 (or newer) | computational | _planned_ | _planned_ |
| ML essentiality predictor | DeeplyEssential / ESM2-based | computational | _planned_ | _planned_ |
| OrthoDB-based Ec transfer | OrthoDB groups (Section 2) | upgrade of `_join_ec_inferred` | _planned_ | replaces `via=gene-symbol match` |

Two architectural notes:

- **Consensus priority is intentional.** `CALL_PRIORITY` (`essential > fitness_defect > unclear > fitness_advantage > non_essential`) means a single "essential" call wins over any number of weaker calls. This is deliberately optimistic: a target seen as essential under one condition by one source is flagged, even if other conditions disagree — the downstream final-ranking step is responsible for weighting per-condition specificity.
- **The Ec-inference path uses simple symbol matching today.** `_join_ec_inferred` does a lowercase gene-symbol intersection between Kp `kp_gene_symbol` and the E. coli essentiality table. Since only ~18% of HS11286 entries carry a canonical gene symbol, an OrthoDB-based upgrade (reusing the Section 2 pipeline) is queued as the highest-leverage extension here.
