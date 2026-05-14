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

## Suggestions

_Audit findings from a 2026-05 literature review; not yet wired into the diagram or Tracks table. Two biggest gaps: (a) the structural-accessibility / TM-topology layer is absent, (b) condition-specific in-vivo Kp expression is PaxDb-baseline only. The §5 sink also conflates target-side Clp-accessibility with recruiter-side delivery — these should be separated architecturally._

### Add

- **[DeepTMHMM](https://www.biorxiv.org/content/10.1101/2022.04.08.487609v1)** as new **5.1d**. Resolves "buried in lipid bilayer vs cytoplasmic domain" — the missing axis for Clp accessibility. Covers α-helical TM, β-barrels and signal peptide in one pass. PSORTb / DeepLocPro give compartment but not how much of the protein faces the cytoplasm.
- **[SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) + [LipoP](https://services.healthtech.dtu.dk/services/LipoP-1.0/)** as new **5.1e**. Discriminates Sec/SPI vs Sec/SPII lipoprotein vs Tat — fixes PSORTb's recurring mis-bin of OM-lipoproteins as periplasmic. Resolves the four-way OM-lipo / IM-lipo / periplasmic-soluble / cytoplasmic split.
- **Structure-derived surface exposure (new track 5.3)** — [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) / SASA on AlphaFold + per-target pocket druggability number. AF2-derived RSA correlates ~0.815 with native; the 2025 ["pocketome universe"](https://www.cell.com/cell/fulltext/S0092-8674(22)00593-1) approach already runs this proteome-wide for 11 species. Currently absent.
- **Condition-specific / in-vivo Kp expression (new 5.2d)** — curated [PRIDE](https://www.ebi.ac.uk/pride/) Kp datasets (e.g. [PXD047744](https://www.ebi.ac.uk/pride/archive/projects/PXD047744)), colistin / meropenem / serum-stress proteomics, Tn-seq fitness in lung ([Bachman 2015 *mBio*](https://journals.asm.org/doi/10.1128/mbio.00775-15)) and urine + serum ([Short 2024 *eLife*](https://elifesciences.org/articles/88971)). PaxDb baseline is the wrong endpoint for an infection drug — niche-expressed targets are missed.
- **[PRED-TMBB2](https://academic.oup.com/bioinformatics/article/32/17/i665/2450774) β-barrel cross-check** (5.1d sibling) — orthogonal call for OM β-barrels (MCC 0.92). A β-barrel mis-call as cytoplasmic is the worst-case error.
- **Mtb cyclomarin / ecumicin + B. subtilis ADEP proteomics → §5.2b** — mechanistically closer to BacPROTAC than ADEP4 (recruits ClpC, doesn't just dysregulate ClpP). See the [BacPROTAC Gram-negative review (Front Chem 2024)](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2024.1358539/full).

### Upgrade

- **§5.1b PSORTb 3.0: keep, but cross-check the lipoprotein call with SignalP 6.0 + LipoP.** PSORTb's lipoprotein call is its weakest point.
- **Compose an explicit "Clp-accessibility score"** in the §5 sink rather than vote-counting: cytoplasmic = 1.0, IM with ≥30% cytoplasmic domain = 0.6, periplasmic-soluble = 0.2, IM-mostly = 0.2, OM-lipoprotein / β-barrel / extracellular = 0.0.
- **§5.2 normalisation: harmonise on [iBAQ](https://www.nature.com/articles/nature10098)** across §5.2a / 5.2b / 5.2c before ortholog lift; tag per-organism evidence quality (*S. aureus* and *E. coli* transfers have very different reliability).
- **Architectural — split the §5 sink** into two sinks: target-side Clp-accessibility vs compound-side recruiter delivery. The latter ([eNTRy rules](https://www.nature.com/articles/nature22308), OmpF / OmpC porins, AcrAB-TolC efflux, siderophore conjugation) is a compound property, not a target property.

### Skip

- [MULocDeep](http://mu-loc.org/), SCLpred-EMS, [TargetP-2.0](https://services.healthtech.dtu.dk/services/TargetP-2.0/) (eukaryotic-biased), [CELLO 2.5](http://cello.life.nctu.edu.tw/) / [BUSCA](http://busca.biocomp.unibo.it/) / [LocTree3](https://rostlab.org/services/loctree3/) (superseded), [Phobius](https://phobius.sbc.su.se/) / MEMSAT-SVM / TOPCONS / TMbed / DeepTMpred (DeepTMHMM covers it), BOMP / BetaTM / BOCTOPUS2 (PRED-TMBB2 better), [ProteomicsDB](https://www.proteomicsdb.org/) / MaxQB (human-centric — PaxDb + direct PRIDE pulls are right for bacteria), standalone transcriptomics (proteomics + Tn-seq more directly answers "is the protein there"; transcriptomics adds noise).
