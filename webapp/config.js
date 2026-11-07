// config.js — single source of truth for the browser's weight sliders, table columns, views,
// gene-card panels, filters and Methods text. It describes the v2 tables exported by
// scripts/webapp/export.py (one row per protein; every axis keeps its own consensus + evidence).
// A column referencing a key absent from the loaded data is skipped automatically.
//
// RULES THIS FILE ENCODES (see DESIGN.md and CLAUDE.md):
//  * Each axis is shown as the axis ships it: `<axis>_consensus` (within-species percentile rank,
//    never comparable across species) beside `<axis>_evidence` (1-3, provenance not outcome).
//  * The weight sliders are the owner's exception to the project's "no composite" rule: they combine
//    the consensus columns IN THE BROWSER ONLY, are OFF by default, and are never exported or stored
//    in the data files. The result is labelled "weighted (yours)", not a project score.
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
  D: "Cell division", V: "Defense", T: "Signaling", M: "Cell envelope", N: "Motility", U: "Secretion", O: "Chaperones", X: "Mobilome",
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
  { key: "essentiality_consensus", label: "Essentiality", short: "Ess", weight: 100, on: true,
    help: "Mean within-species percentile of three predictors. On E. coli it is partly circular (built from the same screens as the measured calls)." },
  { key: "ligands_consensus", label: "Ligand precedence", short: "Lig", weight: 100, on: true,
    help: "Potent ChEMBL ligands for this protein or a bacterial homolog. Zero for the ~97% nobody has screened; that is 'not looked', not 'undruggable'." },
  { key: "pockets_consensus", label: "Pockets", short: "Pkt", weight: 100, on: true,
    help: "Predicted and observed pockets on their own 0-1 scales; a ligand seen in the protein's own PDB structure counts at least 0.8. Tracks protein length (rho about 0.65): use with care." },
  { key: "studiedness_consensus", label: "Novelty", short: "Nov", weight: 0, on: false, invert: true,
    help: "1 minus the studiedness rank: high = few papers on this protein or its closest curated homolog. Off by default." },
];

// ---- table columns -----------------------------------------------------------
// `fixed` columns belong to no group and are always shown, next to Target.
// type: score (0-1 heat) | pct (0-1 as %) | plddt | num | int | binary01 | bool | tier | text | class
const TABLE_COLUMNS = [
  { key: "length", label: "Length (aa)", type: "int", fixed: true, w: 116,
    desc: "Protein length in amino acids (UniProt reference proteome sequence). Degradability probabilities fall steeply with length (rho about -0.6 to -0.7), so read the degradability and pocket ranks against it; the length bands in the filters keep proteins of similar size together." },
  // Degradability
  { key: "degradability_consensus", label: "Rank", type: "score", summary: true, heat: true, group: "Degradability",
    desc: "Degradability consensus (0-1): mean within-species percentile of the ADEP4 and ONC212 probabilities. Never compare across species. Tracks either probability at rho 0.97, so it adds little over them. Rows for E. coli and K. pneumoniae are predictions extrapolated from S. aureus screens." },
  { key: "adep4_prob", label: "ADEP4 prob.", type: "score", heat: true, group: "Degradability",
    desc: "Probability the protein is depleted by ADEP4-activated ClpP (TabPFN-3.5 on ESM-C and ProtT5 embeddings, trained on S. aureus; cluster-grouped CV AUROC 0.88). Strongly length-dependent: short proteins score higher." },
  { key: "onc212_prob", label: "ONC212 prob.", type: "score", heat: true, group: "Degradability",
    desc: "Same for the ONC212 screen (cluster-grouped CV AUROC 0.77). The two probabilities correlate at rho 0.89: one opinion, not two." },
  { key: "adep4_ortholog_hit", label: "ADEP4 ortholog", type: "binary01", group: "Degradability",
    desc: "A measured S. aureus ortholog is depleted by ADEP4 (tick) or measured and not depleted (cross). Blank = no measured ortholog (about 80% of proteins). Not independent of the model: the ortholog's label was training data." },
  { key: "onc212_ortholog_hit", label: "ONC212 ortholog", type: "binary01", group: "Degradability",
    desc: "Same for ONC212. Blank = no measured S. aureus ortholog." },
  { key: "nn_similarity", label: "NN cosine", type: "score", heat: false, group: "Degradability",
    desc: "Cosine of the nearest S. aureus protein in ESM-C space. Raw cosine, not relatedness: random pairs reach a median of 0.78." },

  // Essentiality
  { key: "essentiality_consensus", label: "Rank", type: "score", summary: true, heat: true, group: "Essentiality",
    desc: "Essentiality consensus (0-1): mean within-species percentile of Geptop, ProteomeLM-Ess and the screen-trained models. Partly circular on E. coli." },
  { key: "essentiality_evidence", label: "Evidence", type: "tier", group: "Essentiality",
    desc: "3 = two or more independent experimental sources, unanimous and agreeing with the consensus; 2 = one source, or several that conflict; 1 = no measurement. Provenance, not outcome. On K. pneumoniae a 3 rests on strains other than HS11286." },
  { key: "deg_essential_any", label: "Measured", type: "binary01", summary: true, group: "Essentiality",
    desc: "Called essential by at least one DEG screen on this exact strain (E. coli K-12 only). Blank on K. pneumoniae: DEG indexes no Klebsiella, so there is no measured call on HS11286." },
  { key: "deg_essential_all", label: "All screens", type: "binary01", group: "Essentiality",
    desc: "Called essential by every DEG screen on this strain. 'Any' and 'all' differ about 3.4-fold on E. coli. Blank on K. pneumoniae." },
  { key: "screens_ess_mean", label: "Screens model", type: "score", heat: false, group: "Essentiality",
    desc: "Mean probability of nine models trained on published screens (not nine independent votes: they read the same embedding). On E. coli this is not a validated prediction." },
  { key: "proteomelm_ess", label: "ProteomeLM", type: "score", heat: false, group: "Essentiality",
    desc: "p(essential) from the ProteomeLM essentiality head (E. coli K-12 is held out of its training set; K. pneumoniae is an unseen species). Comparable within a species only." },
  { key: "geptop_ess", label: "Geptop", type: "score", heat: false, group: "Essentiality",
    desc: "Geptop 2.0 score (orthology to known essential genes). A 0 is a confident non-essential call when orthologs exist and none is essential, otherwise 'no orthologs'." },

  // Ligand precedence
  { key: "ligands_consensus", label: "Rank", type: "score", summary: true, heat: true, group: "Ligand precedence",
    desc: "Ligand precedence consensus (0-1): percentile of potent ChEMBL ligand counts (own, bacterial homologs) and best potency, pinned to 0 for the ~97% with no evidence. Human ligands are excluded (a liability, not a benefit)." },
  { key: "ligands_evidence", label: "Evidence", type: "tier", group: "Ligand precedence",
    desc: "3 = somebody assayed this protein, or a drug is recorded against it; 2 = only a bacterial homolog was assayed; 1 = nothing in ChEMBL (an open question, not a negative)." },
  { key: "n_drugs", label: "# drugs", type: "int", group: "Ligand precedence",
    desc: "Drugs in ChEMBL's mechanism table whose target is this protein or one at least 95% identical. Target-level: a drug on a complex is credited to every member. 0 means none recorded, not undruggable." },
  { key: "n_ligands_own", label: "# lig. own", type: "int", group: "Ligand precedence",
    desc: "Distinct potent (pChEMBL >= 6) ChEMBL molecules for this protein (species-level match)." },
  { key: "n_ligands_bacterial", label: "# lig. bact.", type: "int", group: "Ligand precedence",
    desc: "Distinct potent molecules over all true-bacterial homologs at >= 40% identity (unique molecules, not a per-target sum). The identity floor controls coverage, not transfer reliability." },
  { key: "n_ligands_human", label: "# lig. human", type: "int", group: "Ligand precedence",
    desc: "Potent molecules against the closest human homolog: a selectivity liability, never added to the bacterial count." },
  { key: "best_pactivity_bacterial", label: "Best pActivity", type: "num", group: "Ligand precedence",
    desc: "Best potency among bacterial ligands (pChEMBL; 6 = 1 uM). Blank for about 97% of proteins." },
  { key: "n_assayed_own", label: "# assayed", type: "int", group: "Ligand precedence",
    desc: "Compounds assayed against this protein, active or not. 0 against 0 assayed is an open question; 0 potent against many assayed is a measured discouragement." },

  // Pockets / structure
  { key: "pockets_consensus", label: "Score", type: "score", summary: true, heat: true, group: "Pockets",
    desc: "Pockets score (0-1): mean of the P2Rank probability, the fpocket druggability and the PDB and AlphaFill ligand counts (each as n/(n+1)); any protein with a drug-like ligand in its own PDB structures is floored at 0.8. Tracks protein length (rho about 0.65); within length deciles P2Rank is near chance." },
  { key: "pockets_evidence", label: "Evidence", type: "tier", group: "Pockets",
    desc: "3 = this protein has its own PDB structure; 2 = an AlphaFill transplant and no structure; 1 = pocket scores only." },
  { key: "p2rank_score", label: "P2Rank", type: "num", group: "Pockets",
    desc: "Best P2Rank pocket probability-sum on the AlphaFold model (pLDDT >= 70 residues only)." },
  { key: "fpocket_score", label: "fpocket", type: "num", group: "Pockets",
    desc: "Best fpocket druggability score on the AlphaFold model." },
  { key: "n_ligands_pdb", label: "# PDB lig.", type: "int", group: "Pockets",
    desc: "Distinct drug-like ligand scaffolds seen in this protein's own PDB structures (measured). Fusion constructs can carry a partner's ligands." },
  { key: "n_ligands_alphafill", label: "# AlphaFill lig.", type: "int", group: "Pockets",
    desc: "Distinct drug-like scaffolds transplanted from homologs (modelled, never summed with the PDB count)." },
  { key: "n_pdb_structures", label: "# PDB", type: "int", group: "Pockets",
    desc: "PDB entries that are this protein (by sequence, >= 95% identity, >= 50% of the chain)." },
  { key: "af_plddt", label: "pLDDT", type: "plddt", group: "Pockets",
    desc: "Mean AlphaFold pLDDT (0-100). Low = disordered." },

  // Localization
  { key: "localization", label: "Class", type: "loc", w: 126, summary: true, group: "Localization",
    desc: "DeepLocPro 1.0 class: cytoplasm, cytoplasmic membrane, periplasm, outer membrane or extracellular. Always returns a call; 12-15% sit below 0.7 confidence. E. coli K-12 is probably in its training set." },
  { key: "cytoplasmic_fraction", label: "Cytoplasmic", type: "score", heat: true, group: "Localization",
    desc: "TMbed: fraction of the chain predicted cytoplasmic (1 = soluble cytoplasmic). Prefer it over the DeepLocPro label where they disagree." },
  { key: "localization_evidence", label: "Evidence", type: "tier", group: "Localization",
    desc: "3 = both predictors agree, DeepLocPro is confident (>= 0.7) and a curated GO cellular-component term agrees; 2 = one of those; 1 = neither or a GO term contradicts. Not experimental." },

  // Orthology
  { key: "has_human_ortholog", label: "Human ortholog", type: "binary01", w: 132, summary: true, group: "Orthology",
    desc: "A human ortholog by either OrthoFinder or reciprocal best hit (union). Tick = has one (a selectivity liability)." },
  { key: "bacterial_panel_orthologs", label: "Panel conserv.", type: "pct", group: "Orthology",
    desc: "Fraction of 28 bacterial proteomes with an ortholog (species counted, not proteins). A fraction, not a count." },
  { key: "orthology_evidence", label: "Evidence", type: "tier", group: "Orthology",
    desc: "3 = placed in a grouping, OrthoFinder and RBH both found a bacterial ortholog and do not conflict on the human call; 2 = placed but not corroborated; 1 = in no grouping (could not look)." },

  // Studiedness
  { key: "studiedness_consensus", label: "Rank", type: "score", summary: true, heat: true, group: "Studiedness",
    desc: "Studiedness consensus (0-1): mean percentile of three paper counts. High = well studied. Proteins with no donor literature are pinned to 0 (S. aureus only)." },
  { key: "studiedness_evidence", label: "Evidence", type: "tier", group: "Studiedness",
    desc: "3 = a curated SwissProt donor at >= 95% identity (this protein's own literature); 2 = a donor at >= 40%; 1 = nothing in SwissProt resembles it." },
  { key: "n_papers_uniprot_prokaryotic", label: "# papers", type: "int", group: "Studiedness",
    desc: "Curated UniProt references of the best-cited prokaryotic donor at >= 40% identity (the ranking count)." },
  { key: "n_papers_uniprot_own", label: "# papers (own)", type: "int", group: "Studiedness",
    desc: "Union of curated references over this protein's own entries (accession, identical sequence, same species at >= 95%)." },
  { key: "n_papers_pubtator_prokaryotic", label: "# PubTator", type: "int", group: "Studiedness",
    desc: "PubTator-annotated papers of the best PubTator donor, keyed on NCBI GeneID. Blank = no in-scope donor carried a GeneID (about 38% on K. pneumoniae)." },

  // Function
  { key: "cog", label: "COG", type: "class", summary: true, group: "Function",
    desc: "First COG category letter (COGclassifier 2.0 over CDD). 'none' = no hit (E. coli 84.5% classified, K. pneumoniae 79.1%). Multi-letter COGs (9-16%) are in the gene card." },
  { key: "function_evidence", label: "Evidence", type: "tier", group: "Function",
    desc: "3 = both COG and GO-slim annotate it, the GO is UniProt-curated and the COG informative (not R or S); 2 = annotated but not corroborated; 1 = neither." },
  { key: "goslim_terms", label: "GO slim", type: "text", group: "Function",
    desc: "GO-slim (prokaryote) terms, ';'-joined. Empty = not annotated, never 'ruled out'." },
  { key: "proteomes_evidence", label: "Naming evid.", type: "tier", group: "Function",
    desc: "3 = the entry carries its own identity (reviewed or a gene symbol) and a specific protein name; 2 = one of those; 1 = neither." },
];

const LEADING_COL_DESC = {
  rank: "Row rank under the current filters and sort.",
  name: "Gene name (italic) and UniProt accession, the canonical identifier.",
  __c: "Weighted score (yours): the weighted mean of the radar's axes (default: Degradability, Essentiality, Ligand precedence and Pockets equally; Novelty off), renormalised over the axes that have a value. NOT a project score and never exported. Switch it off in the sidebar to sort by a single axis.",
  __conf: "Evidence support for essentiality and ligand precedence: dots reflect the evidence levels (3 = full, 2 = partial, 1 = none). Provenance, not outcome.",
};

// ---- methods / provenance panel ------------------------------------------------
const METHODS = [
  { title: "Reading the table", color: "var(--brand)", body:
    "Each axis gives a <b>rank</b> (0-1, a percentile within that species) and an <b>evidence</b> level (1-3: how well the call is backed; provenance, not outcome). Ranks cannot be compared between species. A dash means not annotated or not measured, never zero. The optional weighted score is yours: it mixes ranks in your browser only.",
    refs: [] },
  { title: "Degradability", color: "var(--amber)", body:
    "TabPFN-3.5 models on ESM-C and ProtT5 embeddings, trained on two <i>S. aureus</i> ClpP-activator screens (ADEP4, ONC212) and applied to every protein here. <b>These are predictions, not measurements</b>: no Gram-negative labels exist. Scores fall with protein length, and membrane proteins are protected from cytoplasmic ClpP. TabPFN is used under its non-commercial licence; this browser is a free, non-profit research resource.",
    refs: [
      ["TabPFN (Prior Labs)", "https://priorlabs.ai/"],
      ["ESM-C (EvolutionaryScale)", "https://www.evolutionaryscale.ai/"],
    ] },
  { title: "Essentiality", color: "var(--crimson)", body:
    "Consensus of Geptop 2.0, ProteomeLM-Ess and models trained on published screens. Measured calls (DEG) exist for E. coli K-12 only, and the consensus is partly circular with them. K. pneumoniae HS11286 has no measured call.",
    refs: [
      ["DEG, Database of Essential Genes", "http://origin.tubic.org/deg/public/index.php"],
      ["ProteomeLM", "https://github.com/Bitbol-Lab/ProteomeLM"],
    ] },
  { title: "Ligand precedence and pockets", color: "var(--cobalt)", body:
    "Counts potent (pChEMBL 6 or more) ChEMBL molecules for the protein or bacterial homologs, plus recorded drugs. About 97% have none because nobody screened them. The pocket score averages P2Rank, fpocket, PDB and AlphaFill ligands on their own scales, with a measured PDB ligand counting at least 0.8; it tracks protein length.",
    refs: [
      ["ChEMBL", "https://www.ebi.ac.uk/chembl/"],
      ["AlphaFold DB", "https://alphafold.ebi.ac.uk/"],
      ["P2Rank", "https://github.com/rdk/p2rank"],
      ["fpocket", "https://github.com/Discngine/fpocket"],
    ] },
  { title: "Localization, function, orthology, studiedness", color: "var(--lime)", body:
    "Localization: DeepLocPro and TMbed, from sequence alone. Function: COG category and GO slim. Orthology: OrthoFinder on 30 proteomes plus reciprocal best hits. Studiedness: paper counts of the best curated donor protein.",
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
  { key: "degradability_consensus", label: "Degradability rank" },
  { key: "essentiality_consensus", label: "Essentiality rank" },
  { key: "ligands_consensus", label: "Ligand precedence rank" },
  { key: "pockets_consensus", label: "Pocket score" },
  { key: "studiedness_consensus", label: "Studiedness rank" },
  { key: "cytoplasmic_fraction", label: "Cytoplasmic fraction" },
  { key: "__c", label: "Weighted score (yours)" },
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
  return String(v).replace(/_/g, " ");
}

// ---- per-page filter bars: intentionally empty; filters live in the sidebar ---------
const VIEW_FILTERS = {};

// ---- range filters (both ends) -------------------------------------------------
// Ranks are within-species percentile-like values (0-1, not comparable between species). Length is log-scaled.
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
// Collapsed, a group shows its `summary` columns (its rank and evidence, or its headline flag);
// expanded, it shows every column. Click a group's header band to toggle it.
const GROUP_ORDER = ["Degradability", "Essentiality", "Ligand precedence", "Pockets", "Localization", "Function", "Orthology", "Studiedness"];
// caveat shown above the table for every expanded group
const GROUP_NOTES = {
  Degradability: "Predicted from S. aureus ClpP-activator screens; a ranking hypothesis. Falls with length; membrane proteins are protected.",
  Essentiality: "Measured calls for E. coli K-12 only; the consensus is partly circular with them. None for K. pneumoniae.",
  "Ligand precedence": "Potent ChEMBL ligands. A zero means nobody looked, not undruggable.",
  Pockets: "Largely a ranking by protein length.",
  Localization: "Sequence-based predictions, not experiments.",
  Function: "First COG letter; none = no COG hit.",
  Orthology: "A human ortholog is a selectivity liability; conservation is a breadth benefit.",
  Studiedness: "Papers on the best curated donor; few papers can mean neglected or no donor.",
};
