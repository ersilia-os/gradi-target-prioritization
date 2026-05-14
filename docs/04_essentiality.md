# Essentiality

Determine if the target is required for fitness or survival in Klebsiella pneumoniae.

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

    subgraph DIRECT [" Direct experimental evidence "]
        direction LR
        TNS_IV["<b>4.1a</b> · Kp Tn-seq (in vitro)"]:::method
        TNS_VV["<b>4.1b</b> · Kp Tn-seq (in vivo)"]:::method
        CRI["<b>4.1c</b> · Kp CRISPRi"]:::method
    end

    ECINF["<b>4.2a</b> · <i>E. coli</i> essentiality transfer"]:::method

    subgraph PRED [" Computational predictions "]
        direction LR
        PLM["<b>4.3a</b> · ProteomeLM-Ess"]:::method
        GEP["<b>4.3b</b> · Geptop 2.0"]:::method
        DEEP["<b>4.3c</b> · DeeplyEssential"]:::method
        FBA["<b>4.3d</b> · FBA on iYL1228"]:::method
    end

    P --> TNS_IV
    P --> TNS_VV
    P --> CRI
    ORTHO --> ECINF
    P --> PLM
    P --> GEP
    P --> DEEP
    P --> FBA

    TNS_IV --> T(["Essentiality annotation"]):::result
    TNS_VV --> T
    CRI    --> T
    ECINF  --> T
    PLM    --> T
    GEP    --> T
    DEEP   --> T
    FBA    --> T
```

## Tracks

| ID | Title | Description | Resources |
| --- | --- | --- | --- |
| 4.1a | Kp Tn-seq (in vitro) | Transposon-insertion essentiality in laboratory media (LB) and under antibiotic stress; deepest available Kp essentiality calls. | Short 2024 (ECL8), Cain 2017 (NJST258 ST258) |
| 4.1b | Kp Tn-seq (in vivo) | Fitness in animal infection models — lung, urine, blood, spleen, liver — yielding niche-specific essentiality. | Bachman 2015, Paczosa 2020, Mike & Bachman 2023, Bachman 2025 (all KPPR1) |
| 4.1c | Kp CRISPRi | Mobile-CRISPRi-seq with graded knockdown across ~870 conditionally-essential genes in KPPR1S. | Jana 2023 |
| 4.2a | *E. coli* essentiality transfer | Keio + TraDIS three-way consensus (Keio ∩ PEC ∩ Goodall) lifted onto Kp via ortholog. | Goodall 2018, Keio (Baba 2006) |
| 4.3a | ProteomeLM-Ess | Whole-proteome transformer with supervised essentiality head; LM-based bacterial SOTA, trained on OGEE v3. | Bitbol Lab 2025 |
| 4.3b | Geptop 2.0 | Orthology + phylogeny-based essentiality scoring on a proteome FASTA. | Wen et al. 2019 |
| 4.3c | DeeplyEssential | DNA + protein deep neural network trained on DEG. | Hasan & Lonardi 2020 |
| 4.3d | FBA on iYL1228 | In silico single-gene knockout on a Kp genome-scale metabolic reconstruction. | Liao et al. 2011 |

## Key resources

| Resource | Description | Tracks |
| --- | --- | --- |
| [Short et al. 2024 (eLife)](https://doi.org/10.7554/eLife.88971.3) | ECL8 TraDIS, >554k unique insertions, LB / urine / serum. | 4.1a |
| [Cain et al. 2017 (Sci Rep)](https://www.nature.com/articles/srep42483) | NJST258 "secondary resistome" under colistin / imipenem / ciprofloxacin. | 4.1a |
| [Bachman et al. 2015 (mBio)](https://journals.asm.org/doi/10.1128/mbio.00775-15) | KPPR1 InSeq in C57BL/6 mouse pneumonia. | 4.1b |
| [Paczosa et al. 2020 (IAI)](https://pubmed.ncbi.nlm.nih.gov/31988174/) | KPPR1 lung fitness in WT vs neutropenic hosts. | 4.1b |
| [Mike & Bachman 2023 (PLoS Pathog)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10381055/) | KPPR1 tissue-specific fitness across blood, spleen, liver, lung. | 4.1b |
| [Bachman et al. 2025 (Nat Commun)](https://www.nature.com/articles/s41467-025-56095-3) | KPPR1 bacteremic dissemination Tn-seq. | 4.1b |
| [Jana et al. 2023 (AEM)](https://journals.asm.org/doi/10.1128/aem.00956-23) | Mobile-CRISPRi-seq for conditionally-essential Kp genes. | 4.1c |
| [Goodall et al. 2018 (mBio)](https://journals.asm.org/doi/10.1128/mbio.02096-17) | *E. coli* BW25113 TraDIS + Keio + PEC consensus essential set. | 4.2a |
| [ProteomeLM](https://github.com/Bitbol-Lab/ProteomeLM) | Whole-proteome transformer; `-Ess` head trained on OGEE v3. | 4.3a |
| [Geptop 2.0](http://cefg.uestc.cn/geptop) | Web server + standalone for phylogeny-aware essentiality. | 4.3b |
| [DeeplyEssential](https://github.com/ucrbioinfo/DeeplyEssential) | DNA + protein deep-NN essentiality predictor. | 4.3c |
| [Liao et al. 2011 (iYL1228)](https://pubmed.ncbi.nlm.nih.gov/21478289/) | Genome-scale metabolic reconstruction of *K. pneumoniae* MGH 78578; supports in silico knockouts. | 4.3d |

## Suggestions

_Audit findings from a 2026-05 literature review; not yet wired into the diagram or Tracks table. Single highest-leverage reframing: from binary essentiality to graded **vulnerability**. BacPROTACs don't knock out — they produce kinetically-controlled, often-incomplete depletion. A target that is "essential" in Tn-seq but tolerates 70% knockdown without growth defect is a poor BacPROTAC target; one that drops fitness sharply at 50% depletion is excellent._

### Add

- **Vulnerability score (graded) as new track 4.4** — refit the [Jana 2023](https://journals.asm.org/doi/10.1128/aem.00956-23) mobile-CRISPRi titration data already in §4.1c into a per-gene vulnerability index, modelled on the [Bosch & Rock 2021 *Cell* Mtb vulnerability index](https://pmc.ncbi.nlm.nih.gov/articles/PMC8382161/) (max fitness cost, sensitivity to partial knockdown, phenotypic lag). Reanalysis of in-hand data; biggest single shift in target ranking.
- **[Enterobacteriaceae TraDIS compendium](https://pubmed.ncbi.nlm.nih.gov/39207104/)** (Ghomi / Jung / Barquist 2024, [data](https://github.com/Gardner-BinfLab/Enterobacteriaceae-TraDIS)) — 13 high-density TIS libraries across *Escherichia / Salmonella / Klebsiella / Citrobacter / Enterobacter* with explicit modelling of essentiality turnover (~⅓ of essential genes switch across genera). Upgrades §4.2a from *E. coli*-only lift to graded consensus.
- **Broad-spectrum vulnerability transfer (new 4.2c)** — pool [Wang / Geisinger 2023 *eLife*](https://pubmed.ncbi.nlm.nih.gov/38126769/) (*A. baumannii* CRISPRi), [Poulsen 2019 *PNAS*](https://www.pnas.org/doi/10.1073/pnas.1900570116) (*P. aeruginosa* FiTnEss) and *Mtb* [Bosch 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8382161/) vulnerability scores via OrthoDB lift. Richer than Geptop 2.0's pure-phylogeny score.
- **Hypervirulent Kp Tn-seq in [*Galleria mellonella*](https://www.frontiersin.org/journals/cellular-and-infection-microbiology/articles/10.3389/fcimb.2025.1643224/full)** (Lin 2025; also [Insua 2021 MDR cKp](https://pubmed.ncbi.nlm.nih.gov/33512418/)). Breaks the KPPR1 monoculture in §4.1b — adds hvKp-specific in-vivo essentiality.
- **[OGEE v3](https://academic.oup.com/nar/article/49/D1/D998/5934414) (bacterial subset) + [DEG 15](http://origin.tubic.org/deg/) direct ortholog lookup** as new **4.3e**. Restrict to the bacterial subset (DEG-Bacteria, OGEE bacterial organisms). ProteomeLM-Ess uses OGEE; DeeplyEssential uses DEG — a direct lookup gives a non-ML interpretable baseline and audits when ProteomeLM is parroting its training set.
- **Synthetic-lethality / redundancy flag (new 4.5)** — informed by [Liu / van Opijnen 2024 *Nat Microbiol*](https://www.nature.com/articles/s41564-024-01759-x) CRISPRi–TnSeq genetic-interaction framework (*S. pneumoniae*). For BacPROTACs, partial degradation only becomes lethal if no buffering paralog exists. Approximate via [STRING](https://string-db.org/) / KEGG / paralog detection until a Kp dataset arrives.
- **Watch-item: [InducTn-seq](https://pubmed.ncbi.nlm.nih.gov/40148565/)** (Christen 2025) — inducible Tn mutagenesis bypassing the in-vivo bottleneck. No Kp dataset yet but the framework will produce one within project horizon.

### Upgrade

- **Integration framework**: replace implicit vote-counting with a weighted ensemble + calibration against the Jana vulnerability score as ground truth. Emit a graded 0–1 vulnerability rather than a binary call.
- **§4.1c (Jana CRISPRi)**: change the output from "list of conditionally-essential genes" to "per-gene vulnerability index" via depletion-curve refit.
- **§4.2a (*E. coli* lift)**: keep Keio ∩ PEC ∩ Goodall as the strict tier, overlay Enterobacteriaceae TraDIS for graded confidence.
- **§4.3a (ProteomeLM-Ess)**: log the OGEE training overlap so we can detect parroting vs genuine generalisation when comparing to the new §4.3e direct lookup.

### Skip

- [PATRIC / BV-BRC](https://www.bv-brc.org/) essentiality calls (funding instability), [KleTy](https://klety.dmicrobe.cn/) / [EnteroBase](https://enterobase.warwick.ac.uk/) (typing only), NPEpredictor / "EpilogeneEss" (no maintained tool), STRING-degree-centrality essentiality (weak vs ProteomeLM-Ess; correlated with conservation already captured), CEG / pan-bacterial compilations (subsumed by OGEE v3 + Enterobacteriaceae TraDIS), zebrafish-Kp (limited data, defer).
