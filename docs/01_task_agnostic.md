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

## Output schema

### 1.0 · Reference proteome

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `uniprot_accession` | string (PK) | one | UniProt AC, e.g. `A6T570`. |
| `locus_tag` | string | one | `KPHS_xxxxx` locus tag from the HS11286 assembly. |
| `gene_symbol` | string, nullable | one | Canonical gene symbol; blank for ~82% of HS11286. |
| `protein_name` | string | one | UniProt recommended protein name. |
| `sequence` | string | one | Amino-acid sequence. |
| `sequence_length` | int | one | Residue count. |

### 1.1a · PANTHER family / subfamily

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `panther_family_id` | string, nullable | one | `PTHRxxxxx`. |
| `panther_family_name` | string, nullable | one | Human-readable family name. |
| `panther_subfamily_id` | string, nullable | one | `PTHRxxxxx:SFx`. |
| `panther_subfamily_name` | string, nullable | one | Human-readable subfamily name. |

### 1.1b · InterPro domains

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `interpro_ids` | list[string] | many | Domain / family / site IPR accessions. |
| `interpro_names` | list[string] | many | Parallel names for the above. |
| `interpro_n_domains` | int | one | Count of distinct InterPro hits. |
| `pfam_ids` | list[string] | many | Pfam subset of InterPro hits. |

### 1.2a · PDB coverage

Implemented by `scripts/04c_pdb_coverage.py` (run per organism via `--organism kpneumoniae|ecoli`) → `output/results/<organism>/<prefix>_pdb_coverage.csv`. Source is **PDBe SIFTS**: per-protein covering structures come from the `best_structures` endpoint (batched POST; gives per-structure UniProt residue ranges + coverage + resolution + method), and apo/holo from `ligand_monomers/{pdb_id}` (404 = apo; otherwise bound monomers filtered against a water/ion/cryo-additive ignore-list). Raw responses cached under `data/processed/<organism>/pdb/{best_structures,ligands}/`. Coverage is **range-based** (UniProt-mapped residue segments). Most of the proteome has no experimental structure: **K. pneumoniae 30 / 5,728 (0.5%)**, **E. coli K-12 1,779 / 4,403 (40.4%)** — among those, median union coverage ≈ 1.00.

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `pdb_has_structure` | bool | one | Any experimental PDB structure maps to this UniProt. |
| `pdb_n_structures` | int | one | Number of distinct covering PDB entries. |
| `pdb_n_chains` | int | one | Number of covering (entry, chain) mappings. |
| `pdb_ids` | list[string] | many | `;`-joined covering PDB IDs. |
| `pdb_best_id`, `pdb_best_chain` | string | one | Best structure (max single-structure coverage, ties broken by resolution) and a representative chain. |
| `pdb_coverage_union` | float, [0,1] | one | Fraction of the UniProt sequence covered by the union of all covering structures. |
| `pdb_coverage_best_single` | float, [0,1] | one | Coverage of the single best structure (`≤ pdb_coverage_union`). |
| `pdb_best_resolution_A` | float, nullable | one | Best (lowest) resolution among covering structures, Å (null for NMR-only). |
| `pdb_best_method` | string | one | Experimental method of the best-resolution structure. |
| `pdb_has_holo` | bool | one | Any covering structure has a non-trivial bound ligand. |
| `pdb_holo_ids` | list[string] | many | `;`-joined covering PDB IDs that are holo. |
| `pdb_ligands` | list[string] | many | `;`-joined real (non-ignored) ligand chem-comp IDs across covering structures. |
| `pdb_apo_holo_per_structure` | string | many | `;`-joined `pdbid:LIG1|LIG2` or `pdbid:apo` per covering entry. |

### 1.2b · AlphaFold structural quality

Implemented by `scripts/04a_alphafold_structures.py` (run per organism via `--organism kpneumoniae|ecoli`) → `output/results/<organism>/<prefix>_alphafold_structure.csv`. Models + PAE are pulled fresh from the [AlphaFold DB API](https://alphafold.ebi.ac.uk/api/prediction/) and cached under `data/processed/<organism>/alphafold/{cif,pae,pred}/` for reuse by later structural tracks (§2 pockets, §3 degron exposure). **Monomer only** — AFDB holds single-chain models and does not predict oligomeric state, so monomer-vs-multimer is a separate track. Coverage: K. pneumoniae **5,727 / 5,728** (only `A0A0H3GQF9`, 3,163 aa, exceeds AFDB's ~2,700-residue limit); E. coli K-12 **4,371 / 4,403** (the 32 misses are demerged/recently-added accessions or oversized chains absent from the current AFDB release). Both proteomes show the same structure: median mean-pLDDT ≈ 91 and ≈ 45% multidomain.

Two confidence axes: **pLDDT** (local, read straight from the AFDB-precomputed fractions) and **PAE** (relative-position confidence, clustered into structural domains via the `pae_to_domains` graph approach — Louvain at resolution 0.3, calibrated so compact low-PAE chains resolve to one domain).

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `af_available` | bool | one | AFDB model exists for this accession. |
| `af_model_id`, `af_model_version` | string, int | one | e.g. `AF-<acc>-F1`, version 6. |
| `af_modeled_len` | int | one | Residues covered by the model. |
| `af_is_complex` | bool | one | Always False (AFDB monomer models) — documents the monomer-only caveat. |
| `af_mean_plddt` | float, [0,100] | one | Mean per-residue pLDDT. |
| `af_frac_high_plddt` | float, [0,1] | one | Fraction with pLDDT > 70 (= confident + very-high). |
| `af_frac_very_high_plddt` | float, [0,1] | one | Fraction with pLDDT ≥ 90. |
| `af_frac_confident_plddt` | float, [0,1] | one | Fraction with pLDDT 70–90. |
| `af_frac_low_plddt` | float, [0,1] | one | Fraction with pLDDT 50–70. |
| `af_frac_very_low_plddt` | float, [0,1] | one | Fraction with pLDDT < 50 — proxy for disorder. |
| `af_pae_max`, `af_pae_mean` | float | one | Max / mean Predicted Aligned Error (Å). |
| `af_n_domains` | int | one | PAE-clustered structural domains (≥10 residues each). |
| `af_largest_domain_frac` | float, [0,1] | one | Largest domain size / modeled length. |
| `af_interdomain_pae_mean` | float, nullable | one | Mean PAE between distinct domains (high ⇒ uncertain relative orientation); NaN if single-domain. |
| `af_is_multidomain` | bool | one | `af_n_domains ≥ 2` (46% of the proteome). |

### 1.3a · BV-BRC protein families

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `plfam_id` | string, nullable | one | Within-Kp PATtyFam local cluster ID. |
| `pgfam_id` | string, nullable | one | Global PATtyFam cluster ID. |

### 1.3b · Within-Kp pan-genome class

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `kp_pangenome_class` | categorical | one | `core` / `soft_core` / `shell` / `cloud`. |
| `kp_pangenome_frac` | float, [0,1] | one | Fraction of reference Kp strains carrying an ortholog. |
| `kp_n_strains_containing` | int | one | Numerator of the above for traceability. |

### 1.3c · Cross-species broad-spectrum

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `broad_spectrum_score` | float, [0,1] | one | Normalised phyletic breadth across ESKAPE-E pathogens. |
| `broad_spectrum_n_pathogens` | int | one | Pathogens (of ESKAPE-E) with an ortholog. |
| `broad_spectrum_pathogens` | list[string] | many | Which pathogens. |

### 1.3d · Selectivity vs human

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `human_ortholog_uniprot` | string, nullable | one | Closest human ortholog AC (blank if none). |
| `human_ortholog_pident` | float, [0,100], nullable | one | Sequence identity, %. |
| `human_ortholog_evalue` | float, nullable | one | BLAST E-value to closest human ortholog. |
| `is_selective_vs_human` | bool | one | `pident < 30%` or no ortholog (threshold TBD). |

### 1.4 · Bibliometric / popularity

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `popularity_tier` | categorical | one | `dark` / `studied` / `well_studied`. |
| `n_europepmc_hits` | int | one | Europe PMC mention count for the protein / gene. |
| `n_uniprot_pubs` | int | one | Publications cited by UniProt for this entry. |
| `uniprot_annotation_score` | int, 1–5 | one | UniProt's annotation-completeness score (stars). |

### 1.5 · ESM2 embeddings (sidecar — not joined into the main table)

| Column | Type | Cardinality | Description |
| --- | --- | --- | --- |
| `esm2_embedding` | float[1280] | one | Per-protein vector from ESM2-650M; stored separately (HDF5 / NPY keyed by `uniprot_accession`). |

> **Implementation note.** The embeddings actually computed for this proteome are **ESM-C 600M** (mean-pooled, 1152-d), stored at `data/processed/embeddings/kp_esmc600m_embeddings.npz` keyed by `uniprot_accession` (per the §1.5 → ESM-C suggestion below).

#### 1.5b · Absolute ESM Atlas coordinates

`scripts/03_esmatlas_coords.py` downloads the **absolute (x, y) coordinates each protein occupies in the [Biohub ESM Atlas](https://biohub.ai/esm/protein/atlas)** (ESMC + ESMFold2 global UMAP of the protein universe) and writes `output/results/kp_esmatlas_coords.csv` (`uniprot_accession`, `protein_hash`, `atlas_x`, `atlas_y`, `matched`, `source`, `atlas_api`) plus `output/plots/kp_esmatlas_projection.png`.

These are the atlas's **own** coordinates, not a local re-projection. The atlas keys every protein by `protein_hash = md5(sequence)` (verified), and its web app serves coordinates from a public, unauthenticated endpoint: `POST https://biohub.ai/esm/protein/api/v1alpha1/umap/coords/by-hash` with body `{"protein_hashes": [...≤10...]}`, returning `{"by_hash": {hash: [x, y]}}` in absolute world space (~0–1024; cf. the `umap/manifest` `transform`, `space_size=1024`). The script hashes each proteome sequence, batches the lookups (10/request, concurrent), and re-queries any misses — the atlas folds some proteins **on demand**, so a hash can be omitted on first request and appear once folded. Coverage: **5,727 / 5,728** proteins located; the lone miss (`A0A0H3GQF9`, 3,163 aa) is absent from the atlas (HTTP 404). `source = biohub_esm_atlas_absolute` in every matched row.

## Suggestions

- **[DeepLocPro 1.0](https://academic.oup.com/bioinformatics/article/40/12/btae677/7900293)** — bacteria-trained subcellular-localization (6 classes); anchor at §5, cross-link here.
- **[SignalP 6.0](https://www.nature.com/articles/s41587-021-01156-3) + LipoP** — five-class signal-peptide / lipoprotein detection (lipoprotein flag = "hard for BacPROTAC").
- **[Foldseek](https://www.nature.com/articles/s41587-023-01773-0) + [ProstT5](https://academic.oup.com/nargab/article/6/4/lqae150/7901286)** — extend §1.2 from coverage to structural-neighbour-with-ligand search; ProstT5 fills in for sequences without an AlphaFold model.
- **[PPanGGOLiN](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007732)** — implements the planned §1.3b within-Kp pan-genome partition (core / persistent / shell / cloud).
- **[eggNOG-mapper v2](https://academic.oup.com/mbe/article/38/12/5825/6379734)** — one pass for COG / KEGG / EC / GO and the human-ortholog flag (covers §1.3d).
- **[MobiDB-lite](https://pmc.ncbi.nlm.nih.gov/articles/PMC7779018/)** — IDR consensus alongside pLDDT bins in §1.2b. Bacterial proteomes have low IDR content, so few proteins trigger.
- **§1.5 ESM-2-650M → [SaProt-650M](https://github.com/westlake-repl/SaProt) or [ESM-C 600M](https://www.evolutionaryscale.ai/blog/esm-cambrian)** — SaProt's 3Di tokens are free given §1.2b structures.
- **§1.2a PDB coverage → coverage + Foldseek neighbours** — SIFTS misses fold-similar ligand-bound entries.
- **§1.3a BV-BRC PATtyFams: snapshot locally** — NIAID funding renewed Sept 2024 but treat as a third-party dependency.

