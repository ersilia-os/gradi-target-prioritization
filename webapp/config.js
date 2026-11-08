// config.js — single source of truth for the browser's weighted-score radar, table columns, views,
// gene-card panels, filters and Methods text. It describes the v2 tables exported by
// scripts/webapp/export.py (one row per protein; every axis keeps its own consensus + evidence).
// A column referencing a key absent from the loaded data is skipped automatically.
//
// RULES THIS FILE ENCODES (see DESIGN.md and CLAUDE.md):
//  * Each axis is shown as the axis ships it: `<axis>_consensus` (within-species 0-1 score,
//    never comparable across species) beside `<axis>_evidence` (1-3, provenance not outcome).
//  * The weighted-score radar is the owner's exception to the project's "no composite" rule: it combines
//    the consensus columns IN THE BROWSER ONLY (on by default, four axes equally, Novelty off) and is
//    never exported or stored in the data files. The result is labelled "Weighted score", not a project score.
//  * A missing value is "not annotated / could not look", never 0.

// ---- axis colours -----------------------------------------------------------
const AXIS_COLORS = {
  composite:     "var(--brand)",
  essentiality:  "var(--crimson)",
  ligandability: "var(--cobalt)",
  degradability: "var(--amber)",
  localization:  "var(--lime)",
  orthology:     "var(--turquoise)",
  studiedness:   "var(--orchid)",
  function:      "var(--periwinkle)",
};
// map a column key to its axis (for colouring)
function axisOf(key) {
  if (key === "__c") return "composite";
  if (/^(adep4|onc212)|degrad|nn_similarity/.test(key)) return "degradability";
  if (/essential|geptop|proteomelm|screens_ess|^deg_/.test(key)) return "essentiality";
  if (/ligand|n_drugs|n_assayed|pactivity|pocket|p2rank|fpocket|plddt|pdb|alphafill/.test(key)) return "ligandability";
  if (/localization|cytoplasmic_fraction/.test(key)) return "localization";
  if (/ortholog|orthology/.test(key)) return "orthology";
  if (/studied|n_papers/.test(key)) return "studiedness";
  if (/^cog|function|goslim|proteomes_evidence/.test(key)) return "function";
  return "composite";
}
function colColor(key) { return AXIS_COLORS[axisOf(key)] || "var(--brand)"; }

// No axis is a mock in v2. Kept as a no-op so the table code needs no special case.
const PROVISIONAL = {};
function isProvisional() { return false; }

// ---- table sections (grouped header band) ------------------------------------
const COLUMN_GROUPS = {
  Weighted: "var(--brand)", Essentiality: "var(--crimson)", "Ligand precedence": "var(--cobalt)",
  Degradability: "var(--amber)", Localization: "var(--lime)", Orthology: "var(--turquoise)",
  Studiedness: "var(--orchid)", Function: "var(--periwinkle)", Pockets: "var(--cobalt)",
};
const AXIS_GROUP_LABEL = {
  composite: "", essentiality: "Essentiality", ligandability: "Ligand precedence", degradability: "Degradability",
  localization: "Localization", orthology: "Orthology", studiedness: "Studiedness", function: "Function",
};
function columnGroup(col) {
  if (col && col.group) return col.group;
  return AXIS_GROUP_LABEL[axisOf(col.key)] || "";
}
function groupColor(name) { return COLUMN_GROUPS[name] || "var(--faint)"; }
const GROUP_DISPLAY = { Degradability: "Degrad.", "Ligand precedence": "Ligand prec.", Studiedness: "Studied", Localization: "Localiz." };
function groupDisplay(name) { return GROUP_DISPLAY[name] || name; }

// ---- function classes = first COG letter (vocabulary: data/source/cdd/cog_func_category.tsv) ----
// "none" = no COG hit (NOT annotated, not a class). Multi-letter COGs keep every letter in
// `cog_categories`; the filter and badge use the first, a display choice stated in the column help.
const FUNCTIONAL_CLASSES = [
  { id: "J", label: "Translation, ribosome",           abbr: "J", color: "var(--periwinkle)" },
  { id: "A", label: "RNA processing",                  abbr: "A", color: "var(--periwinkle)" },
  { id: "K", label: "Transcription",                   abbr: "K", color: "var(--orchid)" },
  { id: "L", label: "Replication and repair",          abbr: "L", color: "var(--fuchsia)" },
  { id: "B", label: "Chromatin structure",             abbr: "B", color: "var(--silver)" },
  { id: "D", label: "Cell cycle, division",            abbr: "D", color: "var(--tangerine)" },
  { id: "V", label: "Defense mechanisms",              abbr: "V", color: "var(--crimson)" },
  { id: "T", label: "Signal transduction",             abbr: "T", color: "var(--cobalt)" },
  { id: "M", label: "Cell wall, membrane, envelope",   abbr: "M", color: "var(--pink)" },
  { id: "N", label: "Cell motility",                   abbr: "N", color: "var(--turquoise)" },
  { id: "U", label: "Secretion, trafficking",          abbr: "U", color: "var(--cobalt)" },
  { id: "O", label: "Chaperones, protein turnover",    abbr: "O", color: "var(--amber)" },
  { id: "X", label: "Mobilome (prophages, transposons)", abbr: "X", color: "var(--silver)" },
  { id: "C", label: "Energy production",               abbr: "C", color: "var(--lime)" },
  { id: "G", label: "Carbohydrate metabolism",         abbr: "G", color: "var(--lime)" },
  { id: "E", label: "Amino acid metabolism",           abbr: "E", color: "var(--lime)" },
  { id: "F", label: "Nucleotide metabolism",           abbr: "F", color: "var(--lime)" },
  { id: "H", label: "Coenzyme metabolism",             abbr: "H", color: "var(--lime)" },
  { id: "I", label: "Lipid metabolism",                abbr: "I", color: "var(--lime)" },
  { id: "P", label: "Inorganic ion transport",         abbr: "P", color: "var(--lime)" },
  { id: "Q", label: "Secondary metabolites",           abbr: "Q", color: "var(--lime)" },
  { id: "R", label: "General function prediction only", abbr: "R", color: "var(--silver)" },
  { id: "S", label: "Function unknown",                abbr: "S", color: "var(--silver)" },
  { id: "none", label: "No COG hit (not annotated)",   abbr: "–", color: "var(--egray)" },
];
// plain-language names for the COG letters, grouped as the COG database groups them (the sidebar shows names, not letters)
const COG_SHORT = { J: "Translation", A: "RNA processing", K: "Transcription", L: "Replication", B: "Chromatin",
  D: "Cell division", V: "Defense", T: "Signaling", M: "Cell envelope", N: "Motility", U: "Secretion", O: "Chaperones, turnover", X: "Mobilome",
  C: "Energy", G: "Carbohydrates", E: "Amino acids", F: "Nucleotides", H: "Coenzymes", I: "Lipids", P: "Ion transport", Q: "Secondary",
  R: "General", S: "Unknown", none: "No COG hit" };
// one sidebar chip can stand for several codes: "Unknown" covers COG S (function unknown) and no COG hit at all
const COG_CHIP_CODES = { S: ["S", "none"] };
const COG_GROUPS = [
  ["Information storage", ["J", "A", "K", "L", "B"]],
  ["Cellular processes", ["D", "V", "T", "M", "N", "U", "O", "X"]],
  ["Metabolism", ["C", "G", "E", "F", "H", "I", "P", "Q"]],
  ["Poorly characterized", ["R", "S"]],
];
const FC_BY_ID = Object.fromEntries(FUNCTIONAL_CLASSES.map((c) => [c.id, c]));

// ---- weighted score (radar) ---------------------------------------------------
// Keys are the v2 consensus columns; `weight` is the default spoke length (0-100). Default: Degradability,
// Essentiality, Ligand precedence and Pockets at 100 each (a quarter apiece), Novelty off. `invert` uses 1 - value.

const COMPONENTS = [
  { key: "degradability_consensus", label: "Degradability", short: "Deg", weight: 100, on: true,
    help: "Rank of the predicted ClpP-activator susceptibility (ADEP4 and ONC212 models, trained on S. aureus screens, applied here). Within-species percentile; a ranking hypothesis for E. coli and K. pneumoniae." },
  { key: "essentiality_consensus", label: "Essentiality", short: "Ess", radar: "Essent.", weight: 100, on: true,
    help: "Mean within-species percentile of three predictors. On E. coli it is partly circular (built from the same screens as the measured calls)." },
  { key: "ligands_consensus", label: "Ligand precedence", short: "Lig", radar: "Ligand prec.", weight: 100, on: true,
    help: "Potent ChEMBL ligands for this protein or a bacterial homolog. Zero for the ~97% nobody has screened; that is 'not looked', not 'undruggable'." },
  { key: "pockets_consensus", label: "Pockets", short: "Pkt", weight: 100, on: true,
    help: "Predicted and observed pockets on their own 0-1 scales; a ligand seen in the protein's own PDB structure counts at least 0.8. Tracks protein length (rho about 0.65): use with care." },
  { key: "studiedness_consensus", label: "Novelty", short: "Nov", weight: 0, on: false, invert: true,
    help: "1 minus the studiedness rank: high = few papers on this protein or its closest curated homolog. Off by default." },
];

// ---- table columns -----------------------------------------------------------
// `fixed` columns belong to no group and are always shown, after all the groups (last column).
// type: score (0-1 heat) | pct (0-1 as %) | plddt | num | int | binary01 | bool | tier | text | class
const TABLE_COLUMNS = [
  { key: "length", label: "Length (aa)", type: "int", fixed: true, w: 116,
    desc: "Protein length in amino acids (UniProt reference proteome sequence). Degradability probabilities fall steeply with length (rho about -0.6 to -0.7), so read the degradability and pocket scores against it; the length bands in the filters keep proteins of similar size together." },
  // Degradability
  { key: "degradability_consensus", label: "Score", head: "Deg. score", type: "score", summary: true, heat: true, group: "Degradability",
    desc: "Degradability consensus (0-1): mean within-species percentile of the ADEP4 and ONC212 probabilities. Never compare across species. Tracks either probability at rho 0.97, so it adds little over them. Rows for E. coli and K. pneumoniae are predictions extrapolated from S. aureus screens." },
  { key: "adep4_prob", label: "ADEP4 prob.", type: "score", heat: true, group: "Degradability",
    desc: "Probability the protein is depleted by ADEP4-activated ClpP (TabPFN-3.5 on ESM-C and ProtT5 embeddings, trained on S. aureus; cluster-grouped CV AUROC 0.88). Strongly length-dependent: short proteins score higher." },
  { key: "onc212_prob", label: "ONC212 prob.", type: "score", heat: true, group: "Degradability",
    desc: "Same for the ONC212 screen (cluster-grouped CV AUROC 0.77). The two probabilities correlate at rho 0.89: one opinion, not two." },
  { key: "adep4_ortholog_hit", label: "ADEP4 orth.", type: "binary01", group: "Degradability",
    desc: "A measured S. aureus ortholog is depleted by ADEP4 (tick) or measured and not depleted (cross). Blank = no measured ortholog (about 80% of proteins). Not independent of the model: the ortholog's label was training data." },
  { key: "onc212_ortholog_hit", label: "ONC212 orth.", type: "binary01", group: "Degradability",
    desc: "Same for ONC212. Blank = no measured S. aureus ortholog." },
  { key: "nn_similarity", label: "NN cosine", type: "score", heat: false, group: "Degradability",
    desc: "Cosine of the nearest S. aureus protein in ESM-C space. Raw cosine, not relatedness: random pairs reach a median of 0.78." },

  // Essentiality
  { key: "essentiality_consensus", label: "Score", head: "Ess. score", type: "score", summary: true, heat: true, group: "Essentiality",
    desc: "Essentiality consensus (0-1): mean within-species percentile of Geptop, ProteomeLM-Ess and the screen-trained models. Partly circular on E. coli." },
  { key: "essentiality_evidence", label: "Evidence", type: "tier", group: "Essentiality",
    desc: "Measured essentiality. 3 = two or more independent screens agree with each other and with the consensus; 2 = one screen, or screens that conflict; 1 = no screen. Unanimous non-essential calls count, so most proteins reach 3. E. coli counts its three K-12 screens (one strain lineage); K. pneumoniae counts screens of other strains (DEG has none on HS11286)." },
  { key: "screens_ess_mean", label: "Screens mod.", type: "score", heat: false, group: "Essentiality",
    desc: "Mean probability of nine models trained on published screens (not nine independent votes: they read the same embedding). On E. coli this is not a validated prediction." },
  { key: "proteomelm_ess", label: "ProteomeLM", type: "score", heat: false, group: "Essentiality",
    desc: "p(essential) from the ProteomeLM essentiality head (E. coli K-12 is held out of its training set; K. pneumoniae is an unseen species). Comparable within a species only." },
  { key: "geptop_ess", label: "Geptop", type: "score", heat: false, group: "Essentiality",
    desc: "Geptop 2.0 score (orthology to known essential genes). A 0 is a confident non-essential call when orthologs exist and none is essential, otherwise 'no orthologs'." },

  // Ligand precedence
  { key: "ligands_consensus", label: "Score", head: "Prec. score", type: "score", summary: true, heat: true, group: "Ligand precedence",
    desc: "Ligand precedence consensus (0-1): percentile of potent ChEMBL ligand counts (own, bacterial homologs) and best potency, pinned to 0 for the ~97% with no evidence. Human ligands are excluded (a liability, not a benefit)." },
  { key: "ligands_evidence", label: "Evidence", type: "tier", group: "Ligand precedence",
    desc: "3 = somebody assayed this protein, or a drug is recorded against it; 2 = only a bacterial homolog was assayed; 1 = nothing in ChEMBL (an open question, not a negative)." },
  { key: "n_drugs", label: "Drugs", type: "int", group: "Ligand precedence",
    desc: "Drugs in ChEMBL's mechanism table whose target is this protein or one at least 95% identical. Target-level: a drug on a complex is credited to every member. 0 means none recorded, not undruggable. The 95% floor is strict between species: K. pneumoniae gyrA, parC, folA and ftsI (90-94% identical to E. coli) read 0." },
  { key: "n_ligands_own", label: "Lig. own", type: "int", group: "Ligand precedence",
    desc: "Distinct potent (pChEMBL >= 6) ChEMBL molecules for this protein (species-level match)." },
  { key: "n_ligands_bacterial", label: "Lig. bact.", type: "int", group: "Ligand precedence",
    desc: "Distinct potent molecules over all true-bacterial homologs at >= 40% identity (unique molecules, not a per-target sum). The identity floor controls coverage, not transfer reliability." },
  { key: "n_ligands_human", label: "Lig. human", type: "int", group: "Ligand precedence",
    desc: "Potent molecules against the closest human homolog: a selectivity liability, never added to the bacterial count." },
  { key: "best_pactivity_bacterial", label: "Best pActivity", type: "num", group: "Ligand precedence",
    desc: "Best potency among bacterial ligands (pChEMBL; 6 = 1 uM). Blank for about 97% of proteins." },
  { key: "n_assayed_own", label: "Assayed", type: "int", group: "Ligand precedence",
    desc: "Compounds assayed against this protein, active or not. 0 against 0 assayed is an open question; 0 potent against many assayed is a measured discouragement." },

  // Pockets / structure
  { key: "pockets_consensus", label: "Score", head: "Pocket score", type: "score", summary: true, heat: true, group: "Pockets",
    desc: "Pockets score (0-1): mean of the P2Rank probability, the fpocket druggability and the PDB and AlphaFill ligand counts (each as n/(n+1)); any protein with a drug-like ligand in its own PDB structures is floored at 0.8. Tracks protein length (rho about 0.65); within length deciles P2Rank is near chance." },
  { key: "pockets_evidence", label: "Evidence", type: "tier", group: "Pockets",
    desc: "Where the pocket information comes from. 3 = this protein has its own PDB structure; 2 = an AlphaFill transplant and no structure; 1 = predicted pockets only. Provenance, not pocket quality." },
  { key: "p2rank_score", label: "P2Rank", type: "num", group: "Pockets",
    desc: "Best P2Rank pocket probability-sum on the AlphaFold model (pLDDT >= 70 residues only)." },
  { key: "fpocket_score", label: "Fpocket", type: "num", group: "Pockets",
    desc: "Best fpocket druggability score on the AlphaFold model." },
  { key: "n_ligands_pdb", label: "PDB lig.", type: "int", group: "Pockets",
    desc: "Distinct drug-like ligand scaffolds seen in this protein's own PDB structures (measured). Fusion constructs can carry a partner's ligands." },
  { key: "n_ligands_alphafill", label: "AlphaFill lig.", type: "int", group: "Pockets",
    desc: "Distinct drug-like scaffolds transplanted from homologs (modelled, never summed with the PDB count)." },
  { key: "n_pdb_structures", label: "PDB entries", type: "int", group: "Pockets",
    desc: "PDB entries that are this protein (by sequence, >= 95% identity, >= 50% of the chain)." },
  { key: "af_plddt", label: "pLDDT", type: "plddt", group: "Pockets",
    desc: "Mean AlphaFold pLDDT (0-100). Low = disordered." },

  // Localization
  { key: "localization", label: "Class", type: "loc", w: 126, summary: true, group: "Localization",
    desc: "DeepLocPro 1.0 class: cytoplasm, cytoplasmic membrane, periplasm, outer membrane or extracellular. Always returns a call; 12-15% sit below 0.7 confidence. E. coli K-12 is probably in its training set." },
  { key: "cytoplasmic_fraction", label: "Cytoplasmic", type: "score", heat: true, group: "Localization",
    desc: "TMbed: fraction of the chain predicted cytoplasmic (1 = soluble cytoplasmic). Prefer it over the DeepLocPro label where they disagree." },
  { key: "localization_evidence", label: "Evidence", type: "tier", group: "Localization",
    desc: "3 = both predictors agree, DeepLocPro is confident (>= 0.7) and a curated GO cellular-component term agrees; 2 = one of those; 1 = neither or a GO term contradicts. Not experimental. Tracks GO annotation more than accuracy: 73% of E. coli inner-membrane proteins reach 3, 22% of cytoplasmic ones." },

  // Orthology
  { key: "has_human_ortholog", label: "Human ortholog", type: "binary01", w: 132, summary: true, group: "Orthology",
    desc: "A human ortholog by OrthoFinder or reciprocal best hit (either). A cross is a selectivity liability; a dash means none found; a question mark means the protein is in no orthogroup, so the call could not be made. Remote homologs count (E. coli yidC, gyrB sit near 25% identity); identity is in the project tables, not here. Sensitive to recall: E. coli secY has no human hit although Sec61a is its homolog." },
  { key: "bacterial_panel_orthologs", label: "Panel cons.", type: "ratio", group: "Orthology",
    desc: "Bacterial proteomes, out of 28, that have an ortholog (species counted, not proteins). Sorts as the fraction." },
  { key: "orthology_evidence", label: "Evidence", type: "tier", group: "Orthology",
    desc: "Quality of the orthology calls. 3 = placed in an orthogroup, OrthoFinder and reciprocal best hit both found a bacterial ortholog, and they do not conflict on the human call; 2 = placed, but one method is missing or they conflict; 1 = in no orthogroup, so both calls beside it could not be made." },

  // Studiedness
  { key: "studiedness_consensus", label: "Score", head: "Stud. score", type: "score", summary: true, heat: true, group: "Studiedness",
    desc: "Studiedness consensus (0-1): mean percentile of three paper counts. High = well studied. Proteins with no donor literature are pinned to 0 (S. aureus only)." },
  { key: "studiedness_evidence", label: "Evidence", type: "tier", group: "Studiedness",
    desc: "3 = a curated SwissProt entry at >= 95% identity (this protein's own literature); 2 = a homolog at >= 40%; 1 = nothing in SwissProt resembles it. Reflects the reference strain: 88% of E. coli K-12 reaches 3, 6% of K. pneumoniae HS11286, an uncurated proteome." },
  { key: "n_papers_uniprot_prokaryotic", label: "Papers", type: "int", group: "Studiedness",
    desc: "Curated UniProt references of the best-cited prokaryotic donor at >= 40% identity (the ranking count)." },
  { key: "n_papers_pubtator_prokaryotic", label: "PubTator", type: "int", group: "Studiedness",
    desc: "PubTator-annotated papers of the best PubTator donor, keyed on NCBI GeneID. Blank = no in-scope donor carried a GeneID (about 38% on K. pneumoniae)." },

  // Function
  { key: "cog", label: "COG", type: "class", summary: true, group: "Function",
    desc: "First COG category letter (COGclassifier 2.0 over CDD). 'none' = no hit (E. coli 84.5% classified, K. pneumoniae 79.1%). Multi-letter COGs (9-16%) are in the gene card." },
  { key: "function_evidence", label: "Evidence", type: "tier", group: "Function",
    desc: "3 = both COG and GO slim annotate it, the GO is UniProt-curated and the COG informative (not R or S); 2 = annotated but not corroborated; 1 = neither. 65-75% of proteins reach 3." },
  { key: "goslim_terms", label: "GO slim", type: "goslim", w: 360, group: "Function",
    desc: "GO-slim (prokaryote) terms by name. Empty = not annotated, never 'ruled out'." },
];

const LEADING_COL_DESC = {
  rank: "Row rank under the current filters and sort.",
  name: "Gene name (italic) and UniProt accession, the canonical identifier.",
  __c: "Weighted score: the weighted mean of the radar's axes (default: Degradability, Essentiality, Ligand precedence and Pockets equally; Novelty off), renormalised over the axes that have a value. NOT a project score and never exported. Switch it off in the sidebar to sort by a single axis.",
  __conf: "Evidence support for essentiality and ligand precedence: dots reflect the evidence levels (3 = full, 2 = partial, 1 = none). Provenance, not outcome.",
};

// ---- methods / provenance panel ------------------------------------------------
const METHODS = [
  { title: "Degradability", color: "var(--amber)", body:
    "TabPFN-3.5 models on ESM-C and ProtT5 embeddings, trained on two <i>S. aureus</i> ClpP-activator screens (ADEP4, ONC212) and applied to every protein. These are predictions, not measurements: no Gram-negative labels exist. Scores fall with protein length, and membrane proteins are protected from cytoplasmic ClpP. TabPFN is used under its non-commercial licence; this browser is a free, non-profit research resource.",
    refs: [
      ["TabPFN (Prior Labs)", "https://priorlabs.ai/"],
      ["ESM-C (EvolutionaryScale)", "https://www.evolutionaryscale.ai/"],
    ] },
  { title: "Essentiality", color: "var(--crimson)", body:
    "Consensus of Geptop 2.0, ProteomeLM-Ess and models trained on published screens. Measured calls (DEG) exist for E. coli mostly.",
    refs: [
      ["DEG, Database of Essential Genes", "http://origin.tubic.org/deg/public/index.php"],
      ["ProteomeLM", "https://github.com/Bitbol-Lab/ProteomeLM"],
    ] },
  { title: "Ligand precedence", color: "var(--cobalt)", body:
    "Counts potent (pChEMBL 6 or more) ChEMBL molecules for the protein or its bacterial homologs, plus recorded drugs. The majority have none, as expected.",
    refs: [
      ["ChEMBL", "https://www.ebi.ac.uk/chembl/"],
    ] },
  { title: "Pockets", color: "var(--cobalt)", body:
    "The pocket score averages P2Rank, Fpocket, PDB and AlphaFill ligands on their own scales, with a measured PDB ligand counting the most.",
    refs: [
      ["AlphaFold DB", "https://alphafold.ebi.ac.uk/"],
      ["P2Rank", "https://github.com/rdk/p2rank"],
      ["Fpocket", "https://github.com/Discngine/fpocket"],
    ] },
  { title: "Localization, function, orthology, studiedness", color: "var(--lime)", body:
    "Localization: DeepLocPro and TMbed, from sequence alone. Function: COG category and GO slim. Orthology: OrthoFinder on 30 proteomes (the four anchors and 26 bacterial comparators) plus reciprocal best hits. Studiedness: paper counts from UniProt and PubTator.",
    refs: [
      ["DeepLocPro", "https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/"],
      ["OrthoFinder", "https://github.com/davidemms/OrthoFinder"],
      ["COGclassifier", "https://github.com/moshi4/COGclassifier"],
    ] },
];

// ---- categorical filters (chip groups) -----------------------------------------
const CATEGORICAL_FILTERS = [
  { key: "localization", label: "Localization" },
  { key: "cog", label: "Function (COG)", meta: FC_BY_ID, order: FUNCTIONAL_CLASSES.map((c) => c.id) },
  { key: "essentiality_evidence", label: "Essentiality evidence" },
  { key: "ligands_evidence", label: "Ligand evidence" },
  { key: "pockets_evidence", label: "Pocket evidence" },
];
// compact abbreviations for categorical values (columns are fixed-width; full value on hover)
const EVID_ABBREV = { 3: "L3", 2: "L2", 1: "L1" };
const CAT_ABBREV = {
  localization: { cytoplasm: "Cyt", cytoplasmic_membrane: "IM", periplasm: "Peri", outer_membrane: "OM", extracellular: "Ext" },
  essentiality_evidence: EVID_ABBREV, ligands_evidence: EVID_ABBREV, pockets_evidence: EVID_ABBREV,
  localization_evidence: EVID_ABBREV, orthology_evidence: EVID_ABBREV, studiedness_evidence: EVID_ABBREV,
  function_evidence: EVID_ABBREV, proteomes_evidence: EVID_ABBREV,
};
// text shown on the gene card for an evidence level
const EVID_TEXT = { 3: "evidence 3", 2: "evidence 2", 1: "evidence 1" };

// boolean filters (tri-state: any / yes / no)
const BOOL_FILTERS = [
  { key: "has_human_ortholog", label: "Human ortholog" },
];

// ---- external links ---------------------------------------------------------------
function externalLinks(row) {
  const acc = row.uniprot_accession;
  const gene = (row.gene || "").trim();
  const links = [
    { label: "UniProt", icon: "", href: `https://www.uniprot.org/uniprotkb/${acc}/entry` },
    { label: "AlphaFold", icon: "", href: `https://alphafold.ebi.ac.uk/entry/${acc}` },
    { label: "InterPro", icon: "", href: `https://www.ebi.ac.uk/interpro/protein/UniProt/${acc}/` },
    { label: "STRING", icon: "", href: `https://string-db.org/cgi/network?identifiers=${encodeURIComponent(acc)}` },
    { label: "ChEMBL", icon: "", href: `https://www.ebi.ac.uk/chembl/g/#search_results/all/query=${encodeURIComponent(acc)}` },
  ];
  const org = (typeof state !== "undefined" && state.org) || "ec";
  const orgName = (ORGANISM_META[org] && ORGANISM_META[org].name) || "Escherichia coli";
  const term = encodeURIComponent(`${gene ? gene + " " : ""}${orgName}`.trim());
  links.push({ label: "PubMed", icon: "", href: `https://pubmed.ncbi.nlm.nih.gov/?term=${term}` });
  return links;
}

// ---- projection map: which fields you can colour points by -------------------------
const MAP_COLORS = [
  { key: "degradability_consensus", label: "Degradability" },
  { key: "essentiality_consensus", label: "Essentiality" },
  { key: "ligands_consensus", label: "Ligand precedence" },
  { key: "pockets_consensus", label: "Pockets" },
  { key: "studiedness_consensus", label: "Studiedness" },
  { key: "cytoplasmic_fraction", label: "Cytoplasmic fraction" },
  { key: "__c", label: "Weighted score" },
  { key: "cog", label: "Function (COG)" },
  { key: "localization", label: "Localization" },
  { key: "essentiality_evidence", label: "Essentiality evidence" },
  { key: "ligands_evidence", label: "Ligand evidence" },
];
const MAP_CAT_PALETTE = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
  "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#A0A0A0", "#6C5CE7"];
const EVID_COLORS = { 3: "var(--lime)", 2: "var(--amber)", 1: "var(--silver)" };
const MAP_CATEGORICAL = {
  cog: Object.fromEntries(FUNCTIONAL_CLASSES.map((c) => [c.id, c.color])),
  localization: { cytoplasm: "var(--lime)", cytoplasmic_membrane: "var(--amber)", periplasm: "var(--orange)",
    outer_membrane: "var(--crimson)", extracellular: "var(--fuchsia)" },
  essentiality_evidence: EVID_COLORS, ligands_evidence: EVID_COLORS,
};
function catLabel(field, v) {
  if (field === "cog") return (FC_BY_ID[v] && FC_BY_ID[v].label) || String(v);
  if (/_evidence$/.test(field)) return `evidence ${v}`;
  if (field === "localization" && typeof LOC_LABEL !== "undefined" && LOC_LABEL[v]) return LOC_LABEL[v];   // as the sidebar and table say it
  return String(v).replace(/_/g, " ");
}

// ---- per-page filter bars: intentionally empty; filters live in the sidebar ---------
const VIEW_FILTERS = {};

// ---- range filters (both ends) -------------------------------------------------
// Scores are within-species 0-1 values (not comparable between species). Length is log-scaled.
const RANGE_FILTERS = [
  { key: "degradability_consensus", label: "Degradability", min: 0, max: 1, step: 0.01, color: "var(--amber)" },
  { key: "essentiality_consensus",  label: "Essentiality",  min: 0, max: 1, step: 0.01, color: "var(--crimson)" },
  { key: "ligands_consensus",       label: "Ligand precedence", min: 0, max: 1, step: 0.01, color: "var(--cobalt)" },
  { key: "pockets_consensus",       label: "Pockets",       min: 0, max: 1, step: 0.01, color: "var(--cobalt)" },
  { key: "studiedness_consensus",   label: "Studiedness",   min: 0, max: 1, step: 0.01, color: "var(--orchid)" },
  { key: "length",                  label: "Length (aa)",   step: 1, scale: "log", color: "var(--periwinkle)" },
];

// ---- organisms ---------------------------------------------------------------------
const ORGANISM_META = {
  ec: { file: "data/ec.json", name: "E. coli",       italic: "E. coli",       strain: "K-12 MG1655", cssAccent: "var(--ec)" },
  kp: { file: "data/kp.json", name: "K. pneumoniae", italic: "K. pneumoniae", strain: "HS11286",     cssAccent: "var(--kp)" },
};

// ---- one single table; each group of columns expands and collapses ------------------
// Collapsed, a group shows its `summary` columns (its score and evidence, or its headline flag);
// expanded, it shows every column. Click a group's header band to toggle it.
const GROUP_ORDER = ["Degradability", "Essentiality", "Ligand precedence", "Pockets", "Localization", "Function", "Orthology", "Studiedness"];
// caveat shown above the table for every expanded group
const GROUP_NOTES = {
  Degradability: "Predicted from S. aureus ClpP-activator screens; a ranking hypothesis. Falls with length; membrane proteins are protected.",
  Essentiality: "Measured calls for E. coli K-12 only; the consensus is partly circular with them. None for K. pneumoniae.",
  "Ligand precedence": "Potent ChEMBL ligands. A zero means nobody looked, not undruggable.",
  Pockets: "Predicted and observed pockets on their own scales; a PDB ligand counts at least 0.8. Tracks protein length.",
  Localization: "Sequence-based predictions, not experiments.",
  Function: "First COG letter; a dash means no COG hit.",
  Orthology: "A human ortholog is a selectivity liability; conservation is a breadth benefit.",
  Studiedness: "Papers on the best curated donor; few papers can mean neglected or no donor.",
};

// ---- GO slim (prokaryote) term names, for the gene card ------------------------------------------------
// Fetched once from EBI QuickGO (2026-11-08) for every GO-slim id that occurs in webapp/data/*.json.
// [name, aspect]: P = biological process, F = molecular function, C = cellular component.
const GO_SLIM = {
  "GO:0000910": ["cytokinesis", "P"],
  "GO:0003677": ["DNA binding", "F"],
  "GO:0003723": ["RNA binding", "F"],
  "GO:0003774": ["cytoskeletal motor activity", "F"],
  "GO:0003924": ["GTPase activity", "F"],
  "GO:0005198": ["structural molecule activity", "F"],
  "GO:0005215": ["transporter activity", "F"],
  "GO:0005576": ["extracellular region", "C"],
  "GO:0005618": ["cell wall", "C"],
  "GO:0005694": ["chromosome", "C"],
  "GO:0005737": ["cytoplasm", "C"],
  "GO:0005840": ["ribosome", "C"],
  "GO:0005856": ["cytoskeleton", "C"],
  "GO:0005886": ["plasma membrane", "C"],
  "GO:0005975": ["carbohydrate metabolic process", "P"],
  "GO:0006091": ["generation of precursor metabolites and energy", "P"],
  "GO:0006260": ["DNA replication", "P"],
  "GO:0006281": ["DNA repair", "P"],
  "GO:0006310": ["DNA recombination", "P"],
  "GO:0006351": ["DNA-templated transcription", "P"],
  "GO:0006355": ["regulation of DNA-templated transcription", "P"],
  "GO:0006412": ["translation", "P"],
  "GO:0006457": ["protein folding", "P"],
  "GO:0006520": ["amino acid metabolic process", "P"],
  "GO:0006575": ["modified amino acid metabolic process", "P"],
  "GO:0006629": ["lipid metabolic process", "P"],
  "GO:0006790": ["sulfur compound metabolic process", "P"],
  "GO:0006810": ["transport", "P"],
  "GO:0006886": ["intracellular protein transport", "P"],
  "GO:0007010": ["cytoskeleton organization", "P"],
  "GO:0007059": ["chromosome segregation", "P"],
  "GO:0007155": ["cell adhesion", "P"],
  "GO:0007163": ["establishment or maintenance of cell polarity", "P"],
  "GO:0008092": ["cytoskeletal protein binding", "F"],
  "GO:0008289": ["lipid binding", "F"],
  "GO:0009101": ["glycoprotein biosynthetic process", "P"],
  "GO:0009279": ["cell outer membrane", "C"],
  "GO:0009975": ["cyclase activity", "F"],
  "GO:0012501": ["programmed cell death", "P"],
  "GO:0016070": ["RNA metabolic process", "P"],
  "GO:0016071": ["mRNA metabolic process", "P"],
  "GO:0016209": ["antioxidant activity", "F"],
  "GO:0016491": ["oxidoreductase activity", "F"],
  "GO:0016740": ["transferase activity", "F"],
  "GO:0016787": ["hydrolase activity", "F"],
  "GO:0016829": ["lyase activity", "F"],
  "GO:0016853": ["isomerase activity", "F"],
  "GO:0016874": ["ligase activity", "F"],
  "GO:0023052": ["signaling", "P"],
  "GO:0030030": ["cell projection organization", "P"],
  "GO:0030163": ["protein catabolic process", "P"],
  "GO:0033013": ["tetrapyrrole metabolic process", "P"],
  "GO:0036211": ["protein modification process", "P"],
  "GO:0036386": ["bacterial nucleoid DNA packaging", "P"],
  "GO:0040033": ["sRNA-mediated post-transcriptional gene silencing", "P"],
  "GO:0042178": ["xenobiotic catabolic process", "P"],
  "GO:0042254": ["ribosome biogenesis", "P"],
  "GO:0042597": ["periplasmic space", "C"],
  "GO:0042995": ["cell projection", "C"],
  "GO:0043226": ["organelle", "C"],
  "GO:0043590": ["bacterial nucleoid", "C"],
  "GO:0044183": ["protein folding chaperone", "F"],
  "GO:0044283": ["small molecule biosynthetic process", "P"],
  "GO:0044550": ["secondary metabolite biosynthetic process", "P"],
  "GO:0045182": ["translation regulator activity", "F"],
  "GO:0048870": ["cell motility", "P"],
  "GO:0050896": ["response to stimulus", "P"],
  "GO:0051604": ["protein maturation", "P"],
  "GO:0055086": ["nucleobase-containing small molecule metabolic process", "P"],
  "GO:0060089": ["molecular transducer activity", "F"],
  "GO:0060090": ["molecular adaptor activity", "F"],
  "GO:0061024": ["membrane organization", "P"],
  "GO:0065003": ["protein-containing complex assembly", "P"],
  "GO:0071554": ["cell wall organization or biogenesis", "P"],
  "GO:0071941": ["nitrogen cycle metabolic process", "P"],
  "GO:0072657": ["protein localization to membrane", "P"],
  "GO:0090729": ["toxin activity", "F"],
  "GO:0098542": ["defense response to other organism", "P"],
  "GO:0098754": ["detoxification", "P"],
  "GO:0098772": ["molecular function regulator activity", "F"],
  "GO:0140096": ["catalytic activity, acting on a protein", "F"],
  "GO:0140097": ["catalytic activity, acting on DNA", "F"],
  "GO:0140098": ["catalytic activity, acting on RNA", "F"],
  "GO:0140104": ["molecular carrier activity", "F"],
  "GO:0140110": ["transcription regulator activity", "F"],
  "GO:0140299": ["molecular sensor activity", "F"],
  "GO:0140313": ["molecular sequestering activity", "F"],
  "GO:0140657": ["ATP-dependent activity", "F"],
  "GO:1901135": ["carbohydrate derivative metabolic process", "P"]
};

// ---- protein-card badges ----------------------------------------------------------------------
// Every badge on the card comes from one rule below: a fixed list, each reading one column with a stated
// threshold, so the badges are the same for every protein and the Methods dialog can list them all.
// tone: good = helps a degrader target, bad = hurts, info = context. `test` returns the badge text or null.
// (nz, isYes, isNo, plural, fmtN, LOC_LABEL and PANEL_SIZE are defined in app.js and used only when a card opens.)
const FACT_TOP = 0.9, FACT_BOTTOM = 0.1;   // top / bottom decile of a within-species score
const FACT_RULES = [
  { id: "ess_measured", tone: "good", rule: "Called essential in at least one DEG screen on this strain",
    test: (r) => isNum(r.deg_n_datasets) && r.deg_n_datasets > 0 && nz(r.deg_n_essential) ? `Essential in ${r.deg_n_essential} of ${r.deg_n_datasets} DEG screens` : null },
  { id: "ess_not", tone: "bad", rule: "In DEG screens on this strain and called essential in none",
    test: (r) => isNum(r.deg_n_datasets) && r.deg_n_datasets > 0 && !nz(r.deg_n_essential) ? `Not essential in ${plural(r.deg_n_datasets, "DEG screen")}` : null },
  { id: "ess_pred", tone: "good", rule: "No DEG screen on this strain and essentiality score in the top decile",
    test: (r) => !(isNum(r.deg_n_datasets) && r.deg_n_datasets > 0) && isNum(r.essentiality_consensus) && r.essentiality_consensus >= FACT_TOP ? "Predicted essential (top 10%)" : null },
  { id: "deg_orth_yes", tone: "good", rule: "A measured S. aureus ortholog is depleted in the ADEP4 or ONC212 screen",
    test: (r) => [r.adep4_ortholog_hit, r.onc212_ortholog_hit].some(isYes) ? `S. aureus ortholog degraded (${["ADEP4", "ONC212"].filter((_, i) => isYes([r.adep4_ortholog_hit, r.onc212_ortholog_hit][i])).join(", ")})` : null },
  { id: "deg_orth_no", tone: "bad", rule: "Measured S. aureus orthologs exist and none is depleted",
    test: (r) => ![r.adep4_ortholog_hit, r.onc212_ortholog_hit].some(isYes) && [r.adep4_ortholog_hit, r.onc212_ortholog_hit].some(isNo) ? "S. aureus ortholog not degraded" : null },
  { id: "deg_top", tone: "good", rule: "Degradability score in the top decile",
    test: (r) => isNum(r.degradability_consensus) && r.degradability_consensus >= FACT_TOP ? "Top 10% degradability" : null },
  { id: "deg_bottom", tone: "bad", rule: "Degradability score in the bottom decile",
    test: (r) => isNum(r.degradability_consensus) && r.degradability_consensus <= FACT_BOTTOM ? "Bottom 10% degradability" : null },
  { id: "drugs", tone: "good", rule: "A drug in ChEMBL's mechanism table targets this protein or a >= 95% identical one",
    test: (r) => nz(r.n_drugs) ? plural(r.n_drugs, "drug") : null },
  { id: "lig_own", tone: "good", rule: "At least one potent (pChEMBL >= 6) ChEMBL ligand for this protein",
    test: (r) => nz(r.n_ligands_own) ? plural(r.n_ligands_own, "potent ligand") : null },
  { id: "lig_bact", tone: "good", rule: "No potent ligand for this protein, but some for bacterial homologs (>= 40% identity)",
    test: (r) => !nz(r.n_ligands_own) && nz(r.n_ligands_bacterial) ? `${plural(r.n_ligands_bacterial, "potent ligand")} on homologs` : null },
  { id: "assayed_none", tone: "bad", rule: "10 or more compounds assayed against this protein and none potent",
    test: (r) => !nz(r.n_ligands_own) && isNum(r.n_assayed_own) && r.n_assayed_own >= 10 ? `${fmtN(r.n_assayed_own)} assayed, none potent` : null },
  { id: "pdb_lig", tone: "good", rule: "A drug-like ligand in this protein's own PDB structures",
    test: (r) => nz(r.n_ligands_pdb) ? `${plural(r.n_ligands_pdb, "ligand")} in PDB structures` : null },
  { id: "pdb_struct", tone: "info", rule: "A PDB structure of this protein, with no drug-like ligand",
    test: (r) => !nz(r.n_ligands_pdb) && nz(r.n_pdb_structures) ? plural(r.n_pdb_structures, "PDB structure") : null },
  { id: "pdb_none", tone: "info", rule: "No PDB structure of this protein",
    test: (r) => isNum(r.n_pdb_structures) && r.n_pdb_structures === 0 ? "No PDB structure" : null },
  { id: "disordered", tone: "bad", rule: "AlphaFold model with mean pLDDT below 50",
    test: (r) => isNum(r.af_plddt) && r.af_plddt < 50 ? `Mostly disordered (pLDDT ${Math.round(r.af_plddt)})` : null },
  { id: "loc_cyto", tone: "good", rule: "DeepLocPro class is cytoplasm",
    test: (r) => r.localization === "cytoplasm" ? "Cytoplasmic" : null },
  { id: "loc_membrane", tone: "bad", rule: "DeepLocPro class is inner or outer membrane",
    test: (r) => r.localization === "cytoplasmic_membrane" || r.localization === "outer_membrane" ? `${LOC_LABEL[r.localization]} protein` : null },
  { id: "loc_other", tone: "info", rule: "DeepLocPro class is periplasm or extracellular",
    test: (r) => r.localization === "periplasm" || r.localization === "extracellular" ? LOC_LABEL[r.localization] : null },
  { id: "human_orth", tone: "bad", rule: "A human ortholog by OrthoFinder or reciprocal best hit",
    test: (r) => isYes(r.has_human_ortholog) ? "Human ortholog" : null },
  { id: "no_human_orth", tone: "good", rule: "Placed in an orthogroup and no human ortholog by either method",
    test: (r) => isNo(r.has_human_ortholog) && r.orthology_evidence >= 2 ? "No human ortholog" : null },
  { id: "human_lig", tone: "bad", rule: "Potent ligands on the closest human homolog",
    test: (r) => nz(r.n_ligands_human) ? `${plural(r.n_ligands_human, "potent ligand")} on human homolog` : null },
  { id: "conserved", tone: "info", rule: "Ortholog in at least 90% of the 28 bacterial proteomes",
    test: (r) => isNum(r.bacterial_panel_orthologs) && r.bacterial_panel_orthologs >= 0.9 ? `Conserved in ${Math.round(r.bacterial_panel_orthologs * PANEL_SIZE)}/${PANEL_SIZE} bacteria` : null },
  { id: "well_studied", tone: "info", rule: "30 or more papers on the best-cited homolog",
    test: (r) => isNum(r.n_papers_uniprot_prokaryotic) && r.n_papers_uniprot_prokaryotic >= 30 ? `Well studied (${fmtN(r.n_papers_uniprot_prokaryotic)} papers)` : null },
  { id: "unexplored", tone: "info", rule: "Nothing similar in SwissProt (studiedness evidence 1)",
    test: (r) => r.studiedness_evidence === 1 ? "Nothing similar in SwissProt" : null },
];
