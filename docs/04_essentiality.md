# Essentiality / vulnerability assessment

Part 4 of the GraDi target-prioritization pipeline. See
[`pipeline.md`](./pipeline.md) for the index and the diagram style legend.

Essentiality asks "is this target required for fitness or survival?";
vulnerability sharpens it by asking "how much depletion is tolerated?" via
graded CRISPRi knockdown. The layer combines two evidence pathways — direct
*K. pneumoniae* measurements (Tn-seq, CRISPRi) and *E. coli* ortholog
inference — into a per-condition consensus call. It is the most data-source-
heavy axis: several papers contribute, with a mix of "data loaded today",
"parser written but file not staged" stubs, and several "explore more"
extensions worth flagging.

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

    subgraph DIRECT [" Direct experimental evidence "]
        direction LR
        TNS[("Kp Tn-seq screens<br/><sub>Eichelberger 2024 ECL8<br/>✓ loaded (3 flavors)</sub><br/><sub>Bachman 2015 KPPR1<br/>Ramage 2017 KPNIH1<br/>— stubs</sub>")]:::dataset
        CRI[("Kp CRISPRi screens<br/><sub>Zhu 2023 Mobile-CRISPRi-seq<br/>highlights ✓ loaded (8 genes)</sub>")]:::dataset
        ECO[/"E. coli TraDIS<br/><sub>Goodall 2018 BW25113<br/>(stub — feeds Ec-inference path)</sub>"\]:::stub
    end

    DEG["DEG 15 / OGEE v2<br/><sub>consolidated essential-gene databases</sub><br/><sub>(planned)</sub>"]:::planned
    FBA["FBA / iYL1228 metabolic model<br/><sub>in silico single-gene knockout</sub><br/><sub>(planned)</sub>"]:::planned
    ML["ML essentiality predictor<br/><sub>DeeplyEssential / ESM2-based</sub><br/><sub>(planned; ties to the<br/>task-agnostic ESM2 track)</sub>"]:::planned

    TNS --> PARSE["<i>src/essentiality.py</i> · per-paper parsers<br/><sub>emit long-form rows:<br/>(gene_symbol, call, score,<br/>source, condition, flavor)</sub>"]:::method
    CRI --> PARSE
    ECO --> PARSE
    DEG -.-> PARSE
    FBA -.-> PARSE
    ML  -.-> PARSE

    PARSE --> FLAV["Per-flavor consensus by gene_symbol<br/><sub>5 flavors: in_vitro · in_vivo_lung<br/>in_vivo_urine · in_vivo_serum<br/>vulnerability</sub><br/><sub>CALL_PRIORITY:<br/>essential ▶ fitness_defect ▶ unclear<br/>▶ fitness_advantage ▶ non_essential</sub><br/><sub><i>src/assemble.py</i> · _join_flavor_block</sub>"]:::method

    PARSE -->|"E. coli essentials"| ECINF["Kp ↔ Ec gene-symbol match<br/><sub><i>src/assemble.py</i> · _join_ec_inferred</sub>"]:::method
    ECINF -.-> ORTHO_UP["Upgrade: OrthoDB-based Ec transfer<br/><sub>reuse the ligandability layer's<br/>OrthoDB pipeline</sub><br/><sub>(planned; lifts Ec coverage<br/>beyond ~18%)</sub>"]:::planned

    CONS("Cross-strain conservation<br/><sub>BV-BRC PATtyFams (from task-agnostic)</sub>"):::tagnostic

    FLAV  --> OUT(["ess_*_call · score · sources<br/><sub>5 flavor triples + ess_Ec_inferred</sub>"]):::result
    ECINF --> OUT
    ORTHO_UP -.-> OUT
    CONS -.->|"confidence modifier (planned)"| OUT

    OUT -.-> NEXT["→ final target ranking<br/><sub>(covered in subsequent diagram)</sub>"]:::planned
```

## Tracks

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
| OrthoDB-based Ec transfer | OrthoDB groups (from ligandability layer) | upgrade of `_join_ec_inferred` | _planned_ | replaces `via=gene-symbol match` |

Two architectural notes:

- **Consensus priority is intentional.** `CALL_PRIORITY`
  (`essential > fitness_defect > unclear > fitness_advantage > non_essential`)
  means a single "essential" call wins over any number of weaker calls. This
  is deliberately optimistic: a target seen as essential under one condition
  by one source is flagged, even if other conditions disagree — the
  downstream final-ranking step is responsible for weighting per-condition
  specificity.
- **The Ec-inference path uses simple symbol matching today.**
  `_join_ec_inferred` does a lowercase gene-symbol intersection between Kp
  `kp_gene_symbol` and the E. coli essentiality table. Since only ~18% of
  HS11286 entries carry a canonical gene symbol, an OrthoDB-based upgrade
  (reusing the [ligandability layer's](./02_ligandability.md) pipeline) is
  queued as the highest-leverage extension here.

---

**Prev:** [Degradability assessment](./03_degradability.md) ·
[Task-agnostic per-protein annotation](./01_task_agnostic.md) ·
[Ligandability assessment](./02_ligandability.md)
