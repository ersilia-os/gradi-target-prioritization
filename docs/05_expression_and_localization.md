# Expression and localization

Evaluates if the target is present in meaningful amounts and physically reachable by the cytoplasmic Clp machinery.

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

    ORTHO("Cross-species orthologs"):::tagnostic

    subgraph LOC [" Subcellular localization "]
        direction LR
        UPLOC["<b>5.1a</b> · UniProt subcellular location"]:::method
        PSORT["<b>5.1b</b> · PSORTb 3.0"]:::method
        DLP["<b>5.1c</b> · DeepLocPro"]:::method
    end

    subgraph EXPR [" Expression evidence "]
        direction LR
        PAXKP["<b>5.2a</b> · Kp proteomics abundance"]:::method
        SAU["<b>5.2b</b> · <i>S. aureus</i> ADEP4 / ONC212 transfer"]:::method
        ECEXP["<b>5.2c</b> · <i>E. coli</i> expression transfer"]:::method
    end

    P --> UPLOC
    P --> PSORT
    P --> DLP
    P --> PAXKP
    ORTHO --> SAU
    ORTHO --> ECEXP

    UPLOC --> T(["Expression &amp; localization annotation"]):::result
    PSORT --> T
    DLP   --> T
    PAXKP --> T
    SAU   --> T
    ECEXP --> T
```

## Tracks

| ID | Title | Description | Resources |
| --- | --- | --- | --- |
| 5.1a | UniProt subcellular location | Curated localization labels where available. | UniProt |
| 5.1b | PSORTb 3.0 | Rule-based prokaryotic localization predictor (Gram-aware). | PSORTb |
| 5.1c | DeepLocPro | ML-based prokaryotic localization predictor. | DeepLocPro |
| 5.2a | Kp proteomics abundance | Per-protein abundance from PaxDb / public Kp datasets. | PaxDb |
| 5.2b | *S. aureus* ADEP4 / ONC212 transfer | Clp-activator proteomics from *S. aureus*, transferred via ortholog match — empirical Clp-accessibility readout. | OrthoDB, published proteomics |
| 5.2c | *E. coli* expression transfer | *E. coli* abundance lifted onto Kp via OrthoDB. | PaxDb, OrthoDB |

## Key resources

| Resource | Description | Tracks |
| --- | --- | --- |
| [PSORTb](https://www.psort.org/psortb/) | Gold-standard rule-based prokaryotic subcellular-localization predictor. | 5.1b |
| [DeepLocPro](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/) | Deep-learning subcellular-localization predictor for prokaryotes. | 5.1c |
| [PaxDb](https://pax-db.org/) | Integrated absolute-abundance proteomics across organisms. | 5.2a, 5.2c |
