"""A browsable map of every deliverable table this project ships.

Ten axes each write one or two tables per species to `data/processed/<task>/`, 5,728 rows by
8 to 99 columns apiece. The schemas are documented, but across `CLAUDE.md` and ten `docs/*.md`
files -- the record, not something you can scan. This writes the missing overview: one
self-contained HTML page, one card per axis, built from the files on disk.

    output/plots/overview/tables.json   the payload
    output/plots/overview/tables.html   the page, JSON inlined, no network

**The column meanings are hand-written, in `COLUMNS` below.** There is no machine-readable
source for them, and deriving them from the data would be exactly the kind of plausible-wrong
this project avoids. The cost is that they drift, so **the script exits non-zero if a column on
disk has no entry, or an entry names a column that is gone** -- that reconciliation is what
keeps the page honest as the tables change.

**The sample rows are ONE FIXED PROTEIN SET, used in every table.** The first N rows are an
arbitrary slice and read as noise; a fixed set lets you follow the same proteins down the whole
page, which is also the canonical-row-order claim made visible. Resolution is by gene symbol,
so **an unresolved symbol is NAMED on the page, not dropped** -- on *S. aureus* the LPS/Lpt
entries are structurally absent, and that absence is worth seeing.

Env `gradi`. Reads through the `src/` loaders wherever one exists, so the page shows what an
axis's own consumers see rather than what `read_csv` happens to infer.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src import (  # noqa: E402
    degradability,
    essentiality,
    function,
    interest,
    ligandability,
    localization,
    orthology,
    pockets,
    projections,
    proteomes,
    studiedness,
)

SPECIES = ("kpneumoniae", "ecoli", "saureus")
SPECIES_LABEL = {
    "kpneumoniae": "K. pneumoniae",
    "ecoli": "E. coli",
    "saureus": "S. aureus",
}
SPECIES_SHORT = {"kpneumoniae": "Kp", "ecoli": "Ec", "saureus": "Sa"}
SPECIES_STRAIN = {
    "kpneumoniae": "HS11286 · UP000007841",
    "ecoli": "K-12 MG1655 · UP000000625",
    "saureus": "NCTC 8325 · UP000008816",
}

OUT_DIR = REPO_ROOT / "output" / "plots" / "overview"

# How many sample rows, and how we choose them. See the module docstring.
SAMPLE_SYMBOLS = ("gyrA", "gyrB", "clpP", "lpxC", "lptD", "secA", "ftsH", "rpoB", "rpsJ", "folA")
N_RANDOM = 2
SEED = 0


# --------------------------------------------------------------------------------------------
# What each column means. Hand-written, from CLAUDE.md and the `src/` loader docstrings; there
# is no machine-readable source. Reconciled against the files on every run -- see `_reconcile`.
# --------------------------------------------------------------------------------------------

AC = {"uniprot_ac": "UniProt accession. The canonical key in every table, and the row order every "
                    "other matrix in the project follows."}


COLUMNS: dict[str, dict[str, str]] = {
    "proteome": AC | {
        "proteome_evidence": "1-3, replacing is_reviewed, which was degenerate per species (Ec "
                             "and human 100% reviewed, Kp 7 of 5,728). 3 = the entry carries its "
                             "OWN identity (SwissProt-reviewed, or a gene symbol on the anchor "
                             "entry) AND a specific protein name · 2 = one of those · 1 = neither. "
                             "The two FILLED naming tiers do not count as own identity -- they are "
                             "inference. Ec and human have no level 1 because every entry is "
                             "curator-read. is_reviewed itself is in evidence/proteome_full_<sp>.tsv.",
        "gene_name": "gene symbol. Kp 63.4% / Sa 44.6% after three labelled filling tiers, from "
                     "18.4% / 28.2% raw. NEVER the join key -- use locus_tag.",
        "protein_name": "UniProt's recommended or submitted protein name.",
        "sequence": "the amino-acid sequence. The real join key to every external database.",
    },
    "projection": AC | {
        "tsne_x": "openTSNE coordinate over the 1,152-d ESM-C space (multiscale, cosine, PCA-50, "
                  "dof=0.8). PER-SPECIES, no shared frame.",
        "tsne_y": "the second coordinate. Never compute a cross-species distance on these -- and "
                  "(0, 0) is a real position, in the dense centre.",
    },
    "localization": AC | {
        "localization": "DeepLocPro's compartment, from sequence alone. The model always returns a "
                        "call and has no `unknown` class, so 100% coverage is a property of the "
                        "method, not evidence -- and `confidence` is no longer here to say which "
                        "calls are weak. 12-15% of them sit below 0.7; read "
                        "evidence/deeplocpro_<sp>.tsv before trusting one label.",
        "localization_evidence": "1-3. 3 = the two predictors CONCUR, DeepLocPro is confident "
                                 "(>=0.7), AND a curated GO cellular-component term agrees · 2 = "
                                 "one of those · 1 = neither, or a GO term CONTRADICTS the call "
                                 "(94 Kp / 195 Ec / 56 Sa). LEVEL 3 IS NOT 'EXPERIMENTALLY "
                                 "LOCALIZED' -- nothing here is an experiment; the GO term is the "
                                 "only signal that is not a sequence model. It agrees with "
                                 "DeepLocPro 90.8-94.8% where present (29.6-48.4% of proteins). "
                                 "The weak class sinks as it should: `extracellular` reaches 3 for "
                                 "0.4% on Kp against cytoplasm's 15.1%.",
        "cytoplasmic_fraction": "TMbed: fraction of residues on the cytoplasmic side. PREFER this "
                                "over the label where a choice is forced -- it corroborates "
                                "`extracellular`, the weakest class, from outside DeepLocPro.",
    },
    "degradability": AC | {
        "adep4_prob": "seed-averaged out-of-fold p(substrate), TabPFN-3.5 on ESM-C. One comparable "
                      "scale across all 13,020 proteins. RANK on this; 0.5 is the wrong threshold.",
        "onc212_prob": "the ONC212 probability. NOT independent of `adep4_prob` -- they correlate "
                       "at rho 0.89 while the labels agree at only 0.52.",
        "nn_similarity": "max ESM-C cosine to the S. aureus training set -- the price of the "
                         "extrapolation. A fifth of Kp sits in bands with no AUROC estimate at all.",
        "degradability_consensus": "0-1, the mean WITHIN-SPECIES percentile rank of the two "
                                   "activator probabilities. ADDS LITTLE: it tracks either one "
                                   "alone at rho 0.97, because the two are one opinion (rho 0.88). "
                                   "`nn_similarity` is deliberately not an input -- it measures "
                                   "reach, not degradability.",
        "degradability_evidence": "1-3. 3 = the protein was in an activator screen; 2 = predicted "
                                  "within a band where the model has a validated AUROC "
                                  "(nn_similarity >= 0.90); 1 = predicted beyond any validated "
                                  "band. KP AND EC CANNOT EXCEED 2 -- the screens are S. aureus "
                                  "only, so the cap is a consequence, not a rule. Here a 3 is ONE "
                                  "measurement, unlike the generic ladder in src/consensus.py.",
    },
    "essentiality": AC | {
        "geptop_ess": "Geptop 2.0 orthology+phylogeny score, 0-1, always present. CIRCULAR on Ec "
                      "and Sa (both are Geptop references, 53.2% / 58.3% self-derived).",
        "proteomelm_ess": "ProteomeLM-Ess, the authors' own head. An independent fifth opinion "
                          "(rho 0.44 with geptop) -- and its meaning differs per anchor: Ec "
                          "held-out, Sa in-training, Kp unseen species.",
        "essentiality_consensus": "0-1, the mean WITHIN-SPECIES percentile rank of the three "
                                  "predictors -- ranked first because they are not comparable as "
                                  "values. PREDICTORS ONLY: no measurement enters it. The three "
                                  "agree only loosely (rho 0.13-0.57), so this is a consensus of "
                                  "differing opinions. Never compare it across species.",
        "essentiality_evidence": "1-3. 3 = two or more experimental sources, unanimous, AND "
                                 "agreeing with the consensus · 2 = one source, or several that "
                                 "conflict · 1 = no measurement. NOT purely experimental -- a "
                                 "measurement the models contradict lands at 2, so do not read 2 "
                                 "as 'the experiment was weak'. Kp reaches 3 for 3,973 proteins "
                                 "and NOT ONE is measured on HS11286: all rest on >=2 screened "
                                 "K. pneumoniae strains.",
        "screens_ess_mean": "mean probability over the ten published-screen models. The headline "
                        "column's fallback since 2026-10-03 -- it reaches AUROC 0.89-0.96 on the "
                        "three measured Kp screens against Geptop's 0.59-0.81.",
    },
    "orthology": AC | {
        "has_human_ortholog": "selectivity liability. Under-detecting human homology would make a "
                              "target look MORE selective than it is, which is why the search runs "
                              "`--very-sensitive`.",
        "bacterial_panel_orthologs": "A FRACTION, 0-1, not a count -- the share of the bacterial "
                                    "panel sharing this protein's "
                                    "orthogroup. THE DENOMINATOR IS 28, not 26 -- the 26 tier-C "
                                    "comparators plus the three anchors, minus this protein's own "
                                    "species. Counted over SPECIES, never proteins, so a paralog "
                                    "pair does not inflate it. A 0 is MEASURED: OrthoFinder runs "
                                    "de novo on our own FASTAs, so it is a finding, not a miss.",
        "orthology_evidence": "1-3. 3 = placed in a grouping AND both OrthoFinder and RBH found a "
                              "bacterial ortholog AND the two do not conflict on the human call; "
                              "2 = placed but not corroborated; 1 = in NEITHER an OrthoFinder "
                              "orthogroup nor an OrthoDB group, so both columns beside it are "
                              "'could not look'. IT REQUIRES A POSITIVE FINDING, so S. aureus "
                              "reaches 3 for only 27.6% -- the lone Gram-positive among the "
                              "anchors, which is its biology, not our uncertainty. Nothing here "
                              "is an experiment: a 3 is not experimental corroboration.",
    },
    "ligands": AC | {
        "n_ligands_own": "POTENT (pChEMBL >= 6) distinct molecules on THIS protein, species-level.",
        "n_ligands_bacterial": "potent molecules over the bacterial pool. CONTAINS `n_ligands_own` -- "
                               "never sum the two. Distinct parent_molregno over the UNION of "
                               "homologous targets, never a per-target sum.",
        "n_ligands_human": "potent molecules on human targets. A LIABILITY, never summed in.",
        "n_assayed_own": "compounds anyone ASSAYED against this protein, whatever the outcome.",
        "n_assayed_bacterial": "the denominator. A 0 in `n_ligands_bacterial` against 158 assayed "
                               "is a measured discouragement; against 0 it is an open question. "
                               "NA, never 0, when the effort extract is missing.",
        "n_assayed_human": "assayed against human targets.",
        "best_pactivity_bacterial": "max pChEMBL over the bacterial pool. This is ChEMBL's "
                                   "`pchembl_value` renamed, not a new quantity.",
        "ligands_consensus": "0-1, how much POTENT precedent, ranked WITHIN the ~3% of proteins "
                             "that have any. EVERYTHING ELSE IS EXACTLY 0 -- a deliberate "
                             "departure from the other axes' plain percentile, because under "
                             "average-rank ties the 97% with nothing would read 0.49, mid-scale, "
                             "and 'nobody looked' would appear moderately ligandable. Excludes "
                             "`n_ligands_human` (a liability, points the other way) and "
                             "`n_assayed_*` (effort, not ligandability).",
        "ligands_evidence": "1-3, and it grades PROVENANCE, not outcome -- the consensus already "
                            "carries the outcome, so the two stay independent (rho 0.11-0.24 "
                            "within the evidence-bearing set). 3 = somebody assayed THIS protein; "
                            "2 = only a bacterial homolog was assayed, i.e. transfer, which this "
                            "axis showed CANNOT be calibrated; 1 = nothing in ChEMBL, so the 0 "
                            "beside it is an open question, never a measured negative. KP HAS 11 "
                            "AT LEVEL 3 and they are almost all beta-lactamases -- the only "
                            "K. pneumoniae proteins anyone has screened directly.",
    },
    "pockets": AC | {
        "p2rank_score": "best P2Rank score among ADMITTED pockets -- admitted means the lining "
                        "residues average pLDDT >= 70, the one place model confidence enters "
                        "(v1 applied it twice). On the AlphaFold v6 model. Prefer it to fpocket "
                        "-- but CONTROL FOR LENGTH before quoting any agreement: length alone "
                        "predicts a measured ligand at AUROC 0.65-0.67. Within length deciles, "
                        "against n_ligands_pdb > 0, P2Rank is 0.494 (Kp) / 0.561 (Ec) / 0.619 "
                        "(Sa) and fpocket 0.435 / 0.523 / 0.477 -- so ON THE ANCHOR the pocket "
                        "scores add nothing over protein size. A soft prior at best.",
        "fpocket_score": "best fpocket score, kept for comparison.",
        "n_ligands_pdb": "NOT the ligands axis's n_ligands -- these are SCAFFOLDS seen in a "
                         "structure, not MOLECULES with a measured potency. "
                         "MEASURED: non-redundant drug-like ligands seen bound to THIS "
                         "protein's own PDB structures (>=95% identity chains). Non-redundant = "
                         "distinct Bemis-Murcko generic scaffolds, the ChEMBL axis's definition.",
        "n_ligands_alphafill": "MODELLED: non-redundant drug-like ligands AlphaFill transplanted "
                               "onto the AlphaFold model, from ~30%-identity donors. NEVER add "
                               "this to n_ligands_pdb -- different evidence, 10x the reach.",
        "n_pdb_structures": "distinct PDB entries with a chain that IS this protein (>=95% "
                            "identity, matched by SEQUENCE -- the accession route reached 30 of "
                            "5,728 Kp proteins). Ligand or not: this is structural COVERAGE, the "
                            "question the ligand counts do not ask.",
        "af_plddt": "mean pLDDT of the AlphaFold model. NA means NO MODEL -- and an NA in the "
                    "pocket columns is 'could not look', not 'looked and found nothing'. A "
                    "protein WITH a model and no admitted pocket gets 0. Never fillna(0).",
    },
    "function": AC | {
        "cog_categories": "`;`-joined COG2024 category letters, vocabulary order. Empty = not "
                          "classified (Kp 20.9%). No evidence column beside it, because "
                          "`cogclassifier` vs `none` was 1:1 with empty-or-not on all three "
                          "species -- measured, not assumed.",
        "function_evidence": "1-3, and there is deliberately NO function_consensus -- 'how much "
                             "function' is not a quantity. 3 = BOTH schemes annotate it, the GO is "
                             "UniProt-CURATED and the COG is INFORMATIVE (a letter outside R/S) · "
                             "2 = annotated but not corroborated, NEVER 'badly annotated' (a "
                             "protein with curated GO and no COG hit caps here, and COG tops out "
                             "near 81.6% by NCBI's own curators) · 1 = neither scheme annotates. "
                             "NOT a fame measure: 210 E. coli proteins named 'Uncharacterized' sit "
                             "at 3, because UniProt leaves them unnamed while curating their class.",
        "goslim_terms": "`;`-joined GO ids from the 97-term goslim_prokaryote vocabulary. "
                        "Multi-label by design, built from the *_all columns.",
    },
    "studiedness": AC | {
        "n_papers_uniprot_own": "curated SwissProt references on THIS accession. A paper count, nothing "
                        "scaled. Near-dead on Kp -- 5,710 of 5,728 carry exactly one id, the genome "
                        "paper -- and that IS the measurement of darkness. Do not rank Kp or Sa on it.",
        "n_papers_uniprot_prokaryotic": "references on the best-cited prokaryotic SwissProt homolog. RANK ON "
                           "THIS. Beside `_own` it makes `dark in Klebsiella, famous in E. coli` "
                           "readable off one row.",
        "n_papers_pubtator_prokaryotic": "PubTator3 TEXT-MINED papers on the donor's gene SYMBOL. A "
                                    "third definition, NEVER summed with the other two. It "
                                    "transfers measurably better (0.4054 vs 0.3428 on a held-out "
                                    "E. coli control) but is NOT the ranking: that control is "
                                    "E. coli-only, and E. coli symbols are exactly the ones that "
                                    "entered human nomenclature. EMPTY means a donor exists with "
                                    "no symbol to look up -- not a measured zero.",
    },
}


# --------------------------------------------------------------------------------------------
# The cards. One per axis, in the order they are read, not the order they run -- execution
# order is deliberately unexpressed in this project.
# --------------------------------------------------------------------------------------------

TABLES = [
    dict(key="proteome", axis="proteomes", title="proteome_<sp>.tsv",
         path="data/processed/proteomes/proteome_<sp>.tsv",
         loader=proteomes.load, drop=("sequence",),
         question="Which proteins are we talking about?",
         read_first="This file defines THE ROW ORDER. Every other matrix in the project has the "
                    "same rows in the same order, which is what lets any two axes stack with no "
                    "join at all. Human exists here too (20,416 reviewed entries).",
         also="The provenance columns -- gene_name_source, gene_synonyms, refseq, geneid -- moved "
              "to evidence/proteome_full_<sp>.tsv via load_full(). geneid is NOT idle there: the "
              "literature axis keys NCBI counts on it. `sequence` stays HERE deliberately, because "
              "the house rule is map by sequence, not by accession. evidence/ also holds the "
              "naming audits, the locus tags and the xrefs; src/proteome_registry.tsv drives it."),
    dict(key="projection", axis="embeddings", title="projection_<sp>.tsv",
         path="data/processed/embeddings/projection_<sp>.tsv",
         loader=projections.load,
         question="Where does this protein sit in sequence space?",
         read_first="Coordinates are PER-SPECIES with no shared frame -- never compute a distance "
                    "across species on these columns; use the 1,152-d space. The control is "
                    "trustworthiness (k=10), floor 0.90, because a diverged t-SNE still returns "
                    "plausible numbers.",
         also="The real deliverables of this axis are three .npz matrices, not a table: "
              "embeddings_<sp>.npz (ESM-C, 1,152-d), prott5_<sp>.npz (1,024-d) and "
              "proteomelm_<sp>.npz. They have no accession column -- their rows are positional, "
              "so alignment is ONLY ever by order."),
    dict(key="function", axis="function", title="function_<sp>.tsv",
         path="data/processed/function/function_<sp>.tsv",
         loader=function.load,
         question="What does this protein do?",
         read_first="AN EMPTY LIST MEANS 'NOT ANNOTATED', NOT 'ABSENT'. On K. pneumoniae 1,506 "
                    "proteins (26.3%) carry no GO-slim term because nothing is known about them, "
                    "never because the function was ruled out. Coverage is NOT 100% and must not "
                    "be forced to be: NCBI's own curators reach 81.6% of E. coli K-12 with COG, so "
                    "there is no headroom, and raising the e-value buys noise.",
         also="Both columns are `;`-joined and MULTI-LABEL by design -- 34-52% of annotated "
              "proteins carry more than one slim term (max 11 on Kp). WHAT THIS FORM CANNOT SAY: "
              "a packed list cannot tell a term that is unannotated from one the organism "
              "structurally cannot reach (8 GO terms are eukaryote/plant, Sa has 19 because it is "
              "Gram-positive, COG `Y` is nuclear structure). Those stay as kept columns in "
              "evidence/{goslim,cog}_matrix_<sp>.tsv, which these columns re-expand to exactly. "
              "Nor does it say whether a GO term was UniProt-curated or inferred from an eggNOG "
              "orthogroup -- that is `goslim_source` in evidence/goslim_<sp>.tsv, and it "
              "separates 322 of 13,020 proteins."),
    dict(key="localization", axis="localization", title="localization_<sp>.tsv",
         path="data/processed/localization/localization_<sp>.tsv",
         loader=localization.load,
         question="Which compartment is this protein in, and how much of it faces the cytoplasm?",
         read_first="THE GRAM-POSITIVE TRAP: on S. aureus, `positive` mode does not merely mask the "
                    "two Gram-negative-only classes, it ADDS their probability mass into "
                    "Extracellular. So the absence of periplasm / outer_membrane on Sa is "
                    "STRUCTURAL, not missing data. E. coli K-12 is almost certainly in "
                    "DeepLocPro's training set, so Ec agreement is a sanity check, not validation.",
         also="Two predictors side by side, never reduced to one call: DeepLocPro answers WHICH "
              "compartment, TMbed HOW MUCH of the chain faces the cytoplasm, and where they "
              "disagree the disagreement is the information. `has_signal_peptide` ships in "
              "evidence/tmbed_<sp>.tsv, not here -- it says WHY a fraction is near zero (exported "
              "vs membrane-buried), which is what makes secreted proteins degradable and membrane "
              "proteins protected, so join it back before filtering on localization. "
              "No `evidence` column here, unlike "
              "every other axis: both predictors cover 100% by construction, so it was constant "
              "and said nothing. The per-predictor tables stay in evidence/, because the tracks "
              "run and resume independently."),
    dict(key="degradability", axis="degradability", title="degradability_<sp>.tsv",
         path="data/processed/degradability/degradability_<sp>.tsv",
         loader=degradability.load,
         question="Is this protein a substrate of activated, partnerless ClpP?",
         read_first="E. coli and K. pneumoniae rows are RANKING HYPOTHESES, NOT MEASUREMENTS -- the "
                    "model is trained on two S. aureus activator screens and applied outward. The "
                    "two probability columns are not independent evidence (rho 0.89 while the "
                    "labels agree at 0.52), and their top-100 lists share only 24/100 on Kp.",
         also="A mechanistic finding for filtering: hit rate is cytoplasm 0.180/0.275 but membrane "
              "0.031/0.105 and extracellular 0.049/0.346. 'Cytoplasmic only' is the WRONG filter; "
              "'not membrane' is closer. The MEASURED calls are not columns here -- they were "
              "empty for every Ec and Kp row -- they live in evidence/labels_saureus.tsv with "
              "their continuous log2FCs; src.degradability.measured() re-attaches them."),
    dict(key="essentiality", axis="essentiality", title="essentiality_<sp>.tsv",
         path="data/processed/essentiality/essentiality_<sp>.tsv",
         loader=essentiality.load,
         question="Can the organism live without it?",
         read_first="RANK WITHIN A SPECIES, NEVER ACROSS -- and this table ships THREE PREDICTORS "
                    "AND NO VERDICT, so which one you rank on is now your choice to make. A "
                    "geptop_ess of 0 has two meanings (confident non-essential vs no orthology "
                    "evidence); read geptop_evidence in geptop_<sp>.tsv before reading a zero. "
                    "The MEASURED calls are in deg_<sp>.tsv, where --rule any gives 695 E. coli "
                    "essentials and --rule all gives 205: a 3.4x spread from one choice.",
         also="Six per-source tables feed this one and stay beside it: geptop_, deg_, ogee_, "
              "proteomelm_ess_, screens_<sp>.tsv, plus dataset_registry.tsv. The MEASURED column "
              "lives in deg_<sp>.tsv, which carries both the `any` and `all` reads -- this "
              "summary could only have held whichever --rule picked."),
    dict(key="orthology", axis="orthology", title="orthology_<sp>.tsv",
         path="data/processed/orthology/orthology_<sp>.tsv",
         loader=orthology.load,
         question="Is this protein conserved across bacteria, and does a human have it?",
         read_first="A 0 HERE IS A MEASURED 0 -- OrthoFinder runs de novo on our own FASTAs and the "
                    "stage exits unless every protein is accounted for, so no bacterial ortholog "
                    "is a finding rather than a lookup miss. Median prop is 0.536 Kp / 0.571 Ec / "
                    "0.179 Sa -- that last is S. aureus's Gram-positive isolation against a mostly "
                    "Gram-negative panel, not a defect.",
         also="The 29-column dense table -- orthogroup, paralogs, per-species counts and "
              "identities -- is evidence/orthology_<sp>.tsv via load_dense(). orthologs.tsv and "
              "neighbors.tsv ship beside it and are SPARSE; orthodb_<sp>.tsv is a separate, stable "
              "grouping -- filter it on orthodb_confidence, not on identity."),
    dict(key="ligands", axis="ligands", title="ligands_<sp>.tsv",
         path="data/processed/ligands/ligands_<sp>.tsv",
         loader=ligandability.load,
         question="...and if not, did anyone look?",
         read_first="TWO QUESTIONS, NOT ONE. A 0 against 158 assayed compounds is a measured "
                    "discouragement; a 0 against 0 is an open question. Kp pyrH is the case to "
                    "remember -- 158 compounds against a 98.3%-identical target, not one potent. "
                    "Proteins with a potent ligand are ~2% of each proteome (Kp 113 / Ec 96 / "
                    "Sa 78) and that must not be forced upward -- it is a fact about how little "
                    "of the bacterial proteome anyone has screened.",
         also="`n_ligands_*` HERE IS NOT `n_ligands_*` IN POCKETS: these are distinct MOLECULES "
              "with a measured pChEMBL >= 6 over a taxonomic pool, and they NEST; the pockets ones "
              "are distinct SCAFFOLDS seen in a structure, with no potency, and they are DISJOINT. "
              "WHAT THIS TABLE CANNOT SAY: how many distinct SCAFFOLDS those compounds cover. On "
              "Kp the 5,286 potent compounds collapse to 1,593 Murcko scaffolds (3.3 per "
              "scaffold), which is the difference between 50 starting points and 50 analogues of "
              "one series -- that is `*_n_scaffolds` in evidence/chembl_<sp>.tsv, along with the "
              "identity bands and the match provenance. `pchembl_value` only exists for = "
              "relations on IC50/EC50/Ki/Kd/Potency, so MIC is absent BY CONSTRUCTION: this axis "
              "does not say 'has an antibiotic'."),
    dict(key="pockets", axis="pockets", title="pockets_<sp>.tsv",
         path="data/processed/pockets/pockets_<sp>.tsv",
         loader=pockets.load,
         question="Could a small molecule bind this fold at all?",
         read_first="Three kinds of column, never merged: p2rank/fpocket are PREDICTED on an "
                    "AlphaFold model, n_ligands_pdb and n_pdb_structures are MEASURED, and "
                    "n_ligands_alphafill is a third party's MODEL. Within length deciles the "
                    "pocket scores barely track where ligands are really found (P2Rank 0.49-0.62, "
                    "fpocket 0.44-0.52), so prefer the measured columns where they are non-zero. "
                    "Drug-likeness is built from published sources, not a denylist: QED >= 0.2 "
                    "and Ro3 were measured and REJECTED because both delete antibiotics."),
    dict(key="studiedness", axis="studiedness", title="studiedness_<sp>.tsv",
         path="data/processed/studiedness/studiedness_<sp>.tsv",
         loader=studiedness.load,
         question="How much is already known about it?",
         read_first="Wanted in BOTH directions -- an uncharacterised target is a risk and also the "
                    "novelty the collaboration is looking for. THE NUMBER IS A PAPER COUNT, nothing "
                    "scaled. A zero is an answer, not a gap; never impute it.",
         also="`evidence` is NOT here -- it ships in load_transfer(), so a 0 in n_papers_uniprot_prokaryotic is "
              "ambiguous in this table between `no_hit` (nothing among 575,748 curated entries "
              "resembles the protein -- the strongest novelty claim the axis makes) and "
              "`below_floor`. Read it there before calling a 0 novelty. The unpriced confound is "
              "essentiality: rho 0.38-0.50 with geptop_ess, and alone "
              "among the axes it survives stratification -- stacking the two double-counts."),
]

# The two matrices are EVIDENCE now, not deliverables -- `function_<sp>.tsv` carries the same
# content packed. Kept out of the cards; `src.matrices` still audits them.
MATRICES = []


# --------------------------------------------------------------------------------------------
# Measuring what is actually on disk
# --------------------------------------------------------------------------------------------

def _jsonable(v):
    """numpy and pandas scalars -> something `json.dumps` accepts, NaN -> None."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return None if np.isnan(v) else round(float(v), 4)
    if v is pd.NA or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return str(v)


def _dtype_label(s: pd.Series) -> str:
    if pd.api.types.is_bool_dtype(s):
        return "bool"
    if pd.api.types.is_integer_dtype(s):
        return "int"
    if pd.api.types.is_float_dtype(s):
        return "float"
    return "text"


def column_stats(s: pd.Series) -> dict:
    """Fill rate plus whatever summary the dtype supports. Fill is the number that matters:
    this project's completeness rule means every row exists, so an empty column is a real
    statement about what is unmeasured, never a short file."""
    n = len(s)
    # Several loaders read with `keep_default_na=False`, so an absent value arrives as "" rather
    # than NaN. Counting those as filled made `gene_name` read 100% on a page whose whole point
    # is coverage -- an empty string is an empty cell here, whatever dtype carries it.
    if not pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
        s = s.replace("", np.nan)
    filled = int(s.notna().sum())
    out = {
        "dtype": _dtype_label(s),
        "filled": filled,
        "fill_pct": round(100.0 * filled / n, 1) if n else 0.0,
        "n_unique": int(s.nunique(dropna=True)),
    }
    nonnull = s.dropna()
    if len(nonnull) == 0:
        return out
    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
        vals = pd.to_numeric(nonnull, errors="coerce").dropna()
        if len(vals):
            # A zero is not a missing value anywhere in this project, so report it separately
            # rather than letting it hide inside a fill rate of 100%.
            out |= {
                "min": _jsonable(vals.min()),
                "median": _jsonable(vals.median()),
                "max": _jsonable(vals.max()),
                "zero_pct": round(100.0 * float((vals == 0).mean()), 1),
            }
    elif out["n_unique"] == filled:
        # One value per row: an identifier. Listing its "top categories" would be six arbitrary
        # accessions at 0.0% each, which is noise dressed as a summary.
        out["is_key"] = True
    else:
        vc = nonnull.astype(str).value_counts()
        out["top"] = [{"value": str(k), "n": int(v), "pct": round(100.0 * v / n, 1)}
                      for k, v in vc.head(6).items()]
    return out


def _reconcile(key: str, cols: list[str]) -> list[str]:
    """The page is only honest if its descriptions still name the columns on disk. Returns the
    complaints; the caller aggregates them and exits non-zero. Deliberately both directions --
    a column that quietly disappeared is as wrong as one that quietly appeared."""
    known = COLUMNS[key]
    problems = []
    for c in cols:
        if c not in known:
            problems.append(f"{key}: column {c!r} is on disk with no description in COLUMNS")
    for c in known:
        if c not in cols:
            problems.append(f"{key}: COLUMNS describes {c!r}, which is no longer in the file")
    return problems


def sample_accessions(species: str, prot: pd.DataFrame) -> tuple[list[str], list[str]]:
    """One fixed protein set per species, used in EVERY table on the page -- so a reader follows
    the same proteins all the way down, which is the canonical-row-order claim made visible.

    Returns (accessions, unresolved symbols). The unresolved ones are RETURNED, not swallowed:
    on S. aureus the Lpt/LPS entries are structurally absent, and that absence is the point."""
    by_symbol = {}
    for ac, sym in zip(prot["uniprot_ac"], prot["gene_name"]):
        if isinstance(sym, str) and sym:
            by_symbol.setdefault(sym.lower(), ac)
    picked, missing = [], []
    for sym in SAMPLE_SYMBOLS:
        ac = by_symbol.get(sym.lower())
        if ac is None:
            missing.append(sym)
        elif ac not in picked:
            picked.append(ac)
    rng = np.random.default_rng(SEED)
    pool = [a for a in prot["uniprot_ac"] if a not in picked]
    picked += [str(a) for a in rng.choice(pool, size=min(N_RANDOM, len(pool)), replace=False)]
    return picked, missing


def build_table(spec: dict, species: str, sample: list[str], quiet: bool) -> tuple[dict, list[str]]:
    """One ordinary table, one species: schema, per-column stats, and the sample rows."""
    df = spec["loader"](species)
    path = REPO_ROOT / spec["path"].replace("<sp>", species)
    problems = _reconcile(spec["key"], list(df.columns))

    shown = [c for c in df.columns if c not in spec.get("drop", ())]
    cols = [{"name": c, "meaning": COLUMNS[spec["key"]].get(c, "")} | column_stats(df[c])
            for c in shown]

    idx = df.set_index("uniprot_ac")
    rows = []
    for ac in sample:
        if ac not in idx.index:
            continue
        r = idx.loc[ac]
        rows.append({"uniprot_ac": ac} |
                    {c: _jsonable(r[c]) for c in shown if c != "uniprot_ac"})

    if not quiet:
        print(f"    {path.name:34s} {len(df):>6,} x {len(df.columns):<3} "
              f"{path.stat().st_size / 1024:>7.0f} KB   {len(rows)} sample rows")
    return {
        "key": spec["key"], "axis": spec["axis"], "kind": "table",
        "title": spec["title"].replace("<sp>", SPECIES_SHORT[species]),
        "path": spec["path"].replace("<sp>", species),
        "question": spec["question"], "read_first": spec["read_first"],
        "also": spec.get("also", ""),
        "n_rows": len(df), "n_cols": len(df.columns),
        "size_kb": round(path.stat().st_size / 1024),
        "hidden_cols": list(spec.get("drop", ())),
        "columns": cols, "rows": rows,
    }, problems


def _vocabulary(kind: str) -> dict[str, str]:
    if kind == "goslim":
        t = function.load_goslim_terms()
        return {r.go_id: f"{r.name} ({r.aspect.split('_')[0]})" for r in t.itertuples()}
    cog = pd.read_csv(REPO_ROOT / "data/source/cdd/cog_func_category.tsv", sep="\t", header=None)
    return {str(r[0]): f"{r[3]}" for r in cog.itertuples(index=False)}


def build_matrix(spec: dict, species: str, sample: list[str], quiet: bool) -> dict:
    """A wide binary matrix. 97 columns of mostly zeros do not fit the layout the other cards
    use, so the terms become a counted, sorted list and the sample rows become chips -- only
    the terms those proteins actually carry. Always-zero columns are counted, not hidden:
    this project treats a structural zero as information."""
    df = spec["loader"](species)
    path = REPO_ROOT / spec["path"].replace("<sp>", species)
    vocab = _vocabulary(spec["vocab"])
    terms = [c for c in df.columns if c not in ("uniprot_ac", "evidence")]

    counts = df[terms].sum().astype(int)
    vterms = sorted(
        ({"term": t, "name": vocab.get(t, ""), "n": int(counts[t]),
          "pct": round(100.0 * counts[t] / len(df), 1)} for t in terms),
        key=lambda d: -d["n"],
    )
    idx = df.set_index("uniprot_ac")
    rows = []
    for ac in sample:
        if ac not in idx.index:
            continue
        r = idx.loc[ac]
        on = [t for t in terms if r[t] == 1]
        rows.append({"uniprot_ac": ac, "evidence": _jsonable(r.get("evidence")),
                     "terms": [{"term": t, "name": vocab.get(t, "")} for t in on]})

    all_zero = [v for v in vterms if v["n"] == 0]
    annotated = int((df[terms].sum(axis=1) > 0).sum())
    if not quiet:
        print(f"    {path.name:34s} {len(df):>6,} x {len(df.columns):<3} "
              f"{path.stat().st_size / 1024:>7.0f} KB   {annotated:,} annotated, "
              f"{len(all_zero)} always-zero columns")
    return {
        "key": spec["key"], "axis": spec["axis"], "kind": "matrix",
        "title": spec["title"].replace("<sp>", SPECIES_SHORT[species]),
        "path": spec["path"].replace("<sp>", species),
        "question": spec["question"], "read_first": spec["read_first"],
        "also": spec.get("also", ""),
        "n_rows": len(df), "n_cols": len(df.columns),
        "size_kb": round(path.stat().st_size / 1024),
        "n_terms": len(terms), "n_annotated": annotated,
        "pct_annotated": round(100.0 * annotated / len(df), 1),
        "n_always_zero": len(all_zero),
        "evidence": [{"value": str(k), "n": int(v)}
                     for k, v in df["evidence"].astype(str).value_counts().items()],
        "terms": vterms, "rows": rows,
    }


def unavailable_card(spec: dict, species: str, exc: Exception) -> dict:
    """A card for an axis whose loader will not run -- usually because the stage is mid-rewrite
    and its code expects columns the shipped table does not have yet. Named, never omitted."""
    return {
        "key": spec["key"], "axis": spec["axis"], "kind": "unavailable",
        "title": spec["title"].replace("<sp>", SPECIES_SHORT[species]),
        "path": spec["path"].replace("<sp>", species),
        "question": spec["question"],
        "error": f"{type(exc).__name__}: {exc}",
    }


def build(quiet: bool) -> dict:
    payload = {"species": [], "generated": pd.Timestamp.now().strftime("%Y-%m-%d"), "cards": {}}
    problems: list[str] = []
    unavailable: list[str] = []
    for sp in SPECIES:
        prot = proteomes.load(sp)
        sample, missing = sample_accessions(sp, prot)
        if not quiet:
            print(f"\n  {SPECIES_LABEL[sp]}  ({SPECIES_STRAIN[sp]})  {len(prot):,} proteins")
            print(f"    sample: {len(sample)} proteins"
                  + (f"   unresolved symbols: {', '.join(missing)}" if missing else ""))
        payload["species"].append({
            "key": sp, "label": SPECIES_LABEL[sp], "short": SPECIES_SHORT[sp],
            "strain": SPECIES_STRAIN[sp], "n_proteins": len(prot),
            "sample": sample, "unresolved": missing,
            "sample_names": {
                str(r.uniprot_ac): (r.gene_name if isinstance(r.gene_name, str) and r.gene_name
                                    else "")
                for r in prot[prot["uniprot_ac"].isin(sample)].itertuples()},
        })
        cards = []
        for spec in TABLES:
            # An axis mid-flight must not take the other nine down with it. A loader that RAISES
            # is a different failure from a schema that drifted: the first is loud and obvious,
            # the second is silent, which is why only the second exits non-zero. The card still
            # ships, carrying the error, so the page says what is missing instead of hiding it.
            try:
                card, probs = build_table(spec, sp, sample, quiet)
                problems += probs
            except Exception as exc:                                        # noqa: BLE001
                card = unavailable_card(spec, sp, exc)
                unavailable.append(f"{sp}/{spec['key']}: {type(exc).__name__}: {exc}")
                if not quiet:
                    print(f"    {spec['key']:<34} UNAVAILABLE -- {type(exc).__name__}: {exc}")
            cards.append(card)
        for spec in MATRICES:
            cards.append(build_matrix(spec, sp, sample, quiet))
        payload["cards"][sp] = cards

    if unavailable and not quiet:
        print("\n  UNAVAILABLE -- these axes could not be loaded and ship as a notice on the page:")
        for u in unavailable:
            print(f"    {u}")

    if problems:
        print("\nThe column descriptions no longer match the files:", file=sys.stderr)
        for p in sorted(set(problems)):
            print(f"  {p}", file=sys.stderr)
        print("\nFix COLUMNS in this script -- an overview that describes columns that are not "
              "there is worse than no overview.", file=sys.stderr)
        sys.exit(1)
    return payload


# --------------------------------------------------------------------------------------------
# The page. The markup lives in `tables_overview.html` beside this script -- it is a page, not a
# string, and editing it as one is the difference between tweaking a layout and editing Python.
# --------------------------------------------------------------------------------------------

TEMPLATE = Path(__file__).resolve().parent / "tables_overview.html"
PLACEHOLDER = "/*__PAYLOAD__*/null"


def render(payload: dict) -> str:
    """Inline the payload into the template. Self-contained: no network, no build step."""
    html = TEMPLATE.read_text()
    if PLACEHOLDER not in html:
        raise RuntimeError(f"{TEMPLATE} no longer contains {PLACEHOLDER!r}")
    # `</script>` inside a JSON string would close the tag early -- the one escape that matters.
    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    return html.replace(PLACEHOLDER, blob)


# --------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build a browsable map of every deliverable table the project ships.")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR,
                    help=f"where tables.json and tables.html go (default {OUT_DIR})")
    ap.add_argument("--json-only", action="store_true",
                    help="write the payload but not the page")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    if not args.quiet:
        print("=" * 92)
        print("tables_overview.py -- a browsable map of the project's deliverable tables")
        print(f"  in   data/processed/<task>/   {len(TABLES)} tables + {len(MATRICES)} matrices "
              f"x {len(SPECIES)} species")
        print(f"  out  {args.out_dir.relative_to(REPO_ROOT)}/tables.json, tables.html")
        print("=" * 92)

    payload = build(args.quiet)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "tables.json"
    json_path.write_text(json.dumps(payload, indent=1))
    if not args.quiet:
        print(f"\n  wrote {json_path.relative_to(REPO_ROOT)}  "
              f"({json_path.stat().st_size / 1024:.0f} KB)")

    if not args.json_only:
        html_path = args.out_dir / "tables.html"
        html_path.write_text(render(payload))
        if not args.quiet:
            print(f"  wrote {html_path.relative_to(REPO_ROOT)}  "
                  f"({html_path.stat().st_size / 1024:.0f} KB)")

    if not args.quiet:
        n = sum(len(v) for v in payload["cards"].values())
        print(f"\n  {n} cards over {len(payload['species'])} species. "
              f"Every column reconciled against the files on disk.")


if __name__ == "__main__":
    main()
