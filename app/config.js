// config.js — single source of truth for the browser's composite components,
// table columns, and detail-drawer field groups. Editing this file (not app.js)
// is how the prototype evolves: add a component or column here and the UI picks
// it up. A component/column referencing a key absent from the loaded data is
// automatically greyed out / skipped.

// ---- composite components (weight sliders) --------------------------------
// `key` must match a column emitted by scripts/08a_webapp_export.py.
// `available:false` => shown as a disabled "coming soon" slider regardless of data.
// Each axis has a distinct hue (Ersilia brand plum + NPG accents from stylia)
// so the table reads as a rich, colour-coded scorecard rather than one accent.
const AXIS_COLORS = {
  composite:     "var(--brand)",
  essentiality:  "var(--crimson)",
  ligandability: "var(--cobalt)",     // ligandability + structure
  orthology:     "var(--turquoise)",  // conservation + human-selectivity + selectivity class
  novelty:       "var(--orchid)",
};
// map a column/component key to its axis (for colouring)
function axisOf(key) {
  if (key === "__c") return "composite";
  if (/conservation|human_selective|selectivity|human_homolog|ortholog/.test(key)) return "orthology";
  if (key === "comp_breadth" || /entero_pct|bacteria_pct/.test(key)) return "essentiality";
  if (/essential|proteomelm|geptop|fba|evidence_experimental|evidence_transfer|evidence_predictor/.test(key)) return "essentiality";
  if (/ligand|chembl|bindingdb|pocket|fpocket|p2rank|alphafill|evidence_binding|evidence_structural|evidence_pocket|disorder|structure/.test(key)) return "ligandability";
  if (/novelty|popularity|homolog_|own_n_pubs|own_reviewed/.test(key)) return "novelty";
  if (/^af_|plddt|pdb_|domain/.test(key)) return "ligandability";
  return "composite";
}
function colColor(key) { return AXIS_COLORS[axisOf(key)] || "var(--brand)"; }

// ---- table sections (grouped header band) ---------------------------------
const COLUMN_GROUPS = {
  Composite:     "var(--brand)",
  Essentiality:  "var(--crimson)",
  Ligandability: "var(--cobalt)",
  Novelty:       "var(--orchid)",
  Orthology:     "var(--turquoise)",
  Annotation:    "var(--periwinkle)",
};
// axis → section-band label (fallback when a column has no explicit `group`)
const AXIS_GROUP_LABEL = {
  composite: "Composite", essentiality: "Essentiality", ligandability: "Ligandability",
  orthology: "Orthology", novelty: "Novelty",
};
function columnGroup(col) {
  if (col && col.group) return col.group;
  return AXIS_GROUP_LABEL[axisOf(col.key)] || "";
}
function groupColor(name) { return COLUMN_GROUPS[name] || "var(--faint)"; }

// ---- task-agnostic functional classes (heuristic, from the export) --------
// id must match scripts/08a_webapp_export.py FUNCTIONAL_CLASS_RULES.
const FUNCTIONAL_CLASSES = [
  { id: "transport",                label: "Transport",         abbr: "Tp",  color: "var(--cobalt)" },
  { id: "transcription_regulation", label: "Regulation",        abbr: "Rg",  color: "var(--orchid)" },
  { id: "signaling",                label: "Signaling",         abbr: "Sg",  color: "var(--periwinkle)" },
  { id: "oxidoreductase",           label: "Oxidoreductase",    abbr: "Ox",  color: "var(--crimson)" },
  { id: "transferase",              label: "Transferase",       abbr: "Tf",  color: "var(--tangerine)" },
  { id: "hydrolase_protease",       label: "Hydrolase/protease", abbr: "Hy", color: "var(--amber)" },
  { id: "lyase_isomerase_ligase",   label: "Lyase/isom./ligase", abbr: "Ly", color: "var(--lime)" },
  { id: "ribosomal_translation",    label: "Ribosomal/transl.", abbr: "Rb",  color: "var(--turquoise)" },
  { id: "dna_replication_repair",   label: "DNA repl./repair",  abbr: "DNA", color: "var(--fuchsia)" },
  { id: "cell_envelope",            label: "Cell envelope",     abbr: "Env", color: "var(--pink)" },
  { id: "uncharacterized",          label: "Uncharacterized",   abbr: "Un",  color: "var(--silver)" },
  { id: "other",                    label: "Other",             abbr: "Ot",  color: "var(--egray)" },
];
const FC_BY_ID = Object.fromEntries(FUNCTIONAL_CLASSES.map((c) => [c.id, c]));

// High-level scoring axes only. Secondary metrics (breadth, human-selectivity,
// novelty, conservation, structure) are table columns / tier filters, not weights.
const COMPONENTS = [
  { key: "comp_essentiality",    label: "Essentiality",  short: "Ess", weight: 50, on: true,
    help: "Required for survival (0–1). Higher = better target." },
  { key: "comp_ligandability",   label: "Ligandability", short: "Lig", weight: 50, on: true,
    help: "Small-molecule tractability (0–1). Higher = more druggable." },
  { key: "comp_degradability",   label: "Degradability", short: "Deg", weight: 0,  on: false, available: false,
    help: "Clp-protease susceptibility — axis not implemented in the pipeline yet." },
];

// ---- table columns --------------------------------------------------------
// type: score(0–1 heat) | pct(0–1 as %) | plddt(0–100) | num | int | bool | tier | text
// Fixed leading columns (rank, name, composite) are rendered by app.js; these are the rest.
const TABLE_COLUMNS = [
  { key: "comp_essentiality",    label: "Essentiality",   type: "score", heat: true,  default: true, group: "Essentiality",
    desc: "Essentiality component (0–1): how required the gene is for fitness/survival. From essentiality_score = 0.40·experimental + 0.20·ortholog-transfer + 0.40·predictors (ProteomeLM / Geptop / FBA). Higher = better target." },
  { key: "comp_breadth",         label: "Essential breadth", type: "score", heat: true, default: true, group: "Essentiality",
    desc: "Essential breadth (0–1): fraction of Enterobacteriaceae genomes in which the gene is essential (entero_pct_essential). Higher = essential across more pathogens." },
  { key: "comp_ligandability",   label: "Ligandability",  type: "score", heat: true,  default: true, group: "Ligandability",
    desc: "Ligandability component (0–1): small-molecule tractability. From ligandability_score = 0.45·binding + 0.30·structural + 0.25·pocket, penalised by disorder. Higher = more druggable." },
  { key: "structure_score",      label: "Structure",      type: "score", heat: true,  default: false, group: "Ligandability",
    desc: "Structure score (0–1) = pocket druggability (evidence_pocket): pLDDT-weighted fpocket/P2Rank consensus that a ligandable pocket exists on the AlphaFold model." },
  { key: "conservation_score",   label: "Conservation",   type: "score", heat: true,  default: false, group: "Orthology",
    desc: "Conservation (0–1): fraction of the ~24-species bacterial panel that carries an ortholog (presence/phyletic spread). Distinct from Essential breadth — a protein can be widely conserved yet not broadly essential." },
  { key: "comp_novelty",         label: "Novelty",        type: "score", heat: true,  default: false, group: "Novelty",
    desc: "Novelty / neglectedness (0–1) = 1 − bibliometric studiedness. High = under-studied protein — a more novel target." },
  { key: "comp_human_selective", label: "Human-selective", type: "score", heat: true,  default: true, group: "Orthology",
    desc: "Human-selectivity (0 or 1): 1 = no meaningful human ortholog (safer antibacterial target); 0 = has a human homolog. Derived from the selectivity class." },
  { key: "essentiality_tier",    label: "Ess. tier",      type: "tier",  default: true, group: "Essentiality",
    desc: "Discrete essentiality call: essential / likely_essential / non_essential (thresholded on the essentiality score + direct experimental hits)." },
  { key: "ligandability_tier",   label: "Lig. tier",      type: "tier",  default: true, group: "Ligandability",
    desc: "Discrete tractability call: tractable / partial / intractable (evidence-driven, not a pure score cutoff)." },
  { key: "selectivity",          label: "Selectivity",    type: "tier",  default: true, group: "Orthology",
    desc: "Spectrum × human-homology class: broad/narrow (essential across many pathogens?) crossed with selective/human_homolog (has a human ortholog?). broad_selective is the ideal." },
  { key: "has_hard_evidence",    label: "Hard evid.",     type: "bool",  default: false,
    desc: "Hard ligandability evidence: an experimentally measured ≤1 µM binder, or a genuine (own / ≥95%-identity) co-crystal structure — not just a predicted pocket." },
  { key: "experimentally_essential", label: "Exp. ess.",  type: "bool",  default: false,
    desc: "Called essential by at least one direct experimental screen (Tn-seq / CRISPRi), rather than only by prediction or cross-species transfer." },
  { key: "af_mean_plddt",        label: "pLDDT",          type: "plddt", default: false,
    desc: "Mean AlphaFold pLDDT (0–100): per-residue model confidence, used here as a proxy for how ordered/structured the protein is. Low pLDDT ≈ disorder." },
  { key: "pdb_has_structure",    label: "PDB",            type: "bool",  default: false,
    desc: "An experimental PDB structure exists for this protein (direct or a close sequence match)." },
  { key: "chembl_any_n_potent",  label: "ChEMBL potent",  type: "int",   default: false,
    desc: "Number of ChEMBL compounds tested against this target or an ortholog with potency ≤1 µM (pChEMBL ≥ 6)." },
  { key: "fpocket_max_drug_score", label: "Pocket",       type: "score", heat: false, default: false,
    desc: "Best fpocket druggability score across pockets detected on the AlphaFold model (0–1)." },
  { key: "disorder_frac",        label: "Disorder",       type: "pct",   default: false,
    desc: "Fraction of residues predicted disordered (very-low pLDDT). High disorder lowers the ligandability score." },
  { key: "family",               label: "Family",         type: "int",   default: false,
    desc: "ESM-C embedding cluster id: proteins whose language-model embeddings are similar share a family number (also the clustering behind the Map view)." },

  // --- task-agnostic annotation ---
  { key: "functional_class",     label: "Functional class", type: "class", default: false, group: "Annotation",
    desc: "Coarse functional class (heuristic, keyword-mapped from InterPro/PANTHER family names): transport, regulation, signaling, oxidoreductase, transferase, hydrolase/protease, lyase/isomerase/ligase, ribosomal/translation, DNA replication/repair, cell envelope, uncharacterized, other. An aid for filtering — use the family search for precision." },
  { key: "n_interpro_entries",   label: "InterPro entries", type: "int",  default: false, group: "Annotation",
    desc: "Number of InterPro signatures matched on this protein (families + superfamilies + domains)." },
  { key: "interpro_family_names", label: "InterPro family", type: "text", default: false,
    desc: "InterPro family name(s) — the fine-grained named family. Full text on hover / in the detail panel." },
  { key: "panther_family_names", label: "PANTHER family", type: "text",   default: false,
    desc: "PANTHER family name(s). Full text on hover / in the detail panel." },

  // --- essentiality-axis detail columns ---
  { key: "evidence_experimental", label: "Ess: experimental", type: "score", heat: true, default: false,
    desc: "Experimental essentiality evidence (0–1): strength of direct Tn-seq / CRISPRi essential calls (Kp screens + transferred E. coli screens)." },
  { key: "evidence_transfer",    label: "Ess: transfer",   type: "score", heat: true,  default: false,
    desc: "Cross-species transfer evidence (0–1): essentiality inferred from orthologs in related bacteria." },
  { key: "evidence_predictor",   label: "Ess: predictor",  type: "score", heat: true,  default: false,
    desc: "Predictor evidence (0–1): weighted mean of the available essentiality predictors (ProteomeLM 0.5 / Geptop 0.3 / FBA 0.2)." },
  { key: "proteomelm_ess_score", label: "ProteomeLM",      type: "score", heat: false, default: false,
    desc: "ProteomeLM essentiality probability (0–1): the primary embedding-based predictor (logistic head over ESM-C, trained on E. coli labels)." },
  { key: "geptop_score",         label: "Geptop",          type: "score", heat: false, default: false,
    desc: "Geptop 2.0 essentiality score (0–1): orthology to known essential genes across reference genomes (DIAMOND reimplementation)." },
  { key: "fba_essential",        label: "FBA lethal",      type: "bool",  default: false,
    desc: "Flux-balance analysis: single-gene deletion is lethal in the genome-scale metabolic model (iYL1228 Kp / iML1515 Ec)." },
  { key: "entero_pct_essential", label: "Entero % ess.",   type: "pct",   default: false,
    desc: "Fraction of Enterobacteriaceae genomes in which the gene is essential — the raw value behind the breadth component." },

  // --- ligandability-axis detail columns ---
  { key: "evidence_binding",     label: "Lig: binding",    type: "score", heat: true,  default: false,
    desc: "Binding evidence (0–1): measured bioactivity from ChEMBL / BindingDB, direct or via ortholog." },
  { key: "evidence_structural",  label: "Lig: structural", type: "score", heat: true,  default: false,
    desc: "Structural evidence (0–1): drug-like ligands seen in PDB co-crystals or transplanted by AlphaFill." },
  { key: "evidence_pocket",      label: "Lig: pocket",     type: "score", heat: true,  default: false,
    desc: "Pocket evidence (0–1): pLDDT-weighted fpocket / P2Rank druggable-pocket consensus." },
  { key: "chembl_any_n_compounds", label: "ChEMBL cmpds",  type: "int",   default: false,
    desc: "Total ChEMBL compounds tested against this target or an ortholog." },
  { key: "chembl_any_best_pchembl", label: "Best pChEMBL", type: "num",   default: false,
    desc: "Best pChEMBL potency recorded (higher = more potent; 6 ≈ 1 µM)." },
  { key: "bindingdb_any_n_potent", label: "BindingDB potent", type: "int", default: false,
    desc: "BindingDB compounds with affinity ≤1 µM (direct or ortholog)." },
  { key: "pocket_consensus_score", label: "Pocket consensus", type: "score", heat: false, default: false,
    desc: "Consensus druggable-pocket score (fpocket × P2Rank, pLDDT-weighted)." },
  { key: "alphafill_best_ligand", label: "AlphaFill lig.", type: "text",  default: false,
    desc: "Best drug-like ligand transplanted onto the AlphaFold model by AlphaFill, if any." },
  { key: "human_ligandable_family", label: "Human-drug fam.", type: "bool", default: false,
    desc: "Belongs to a protein family with established human-drug ligandability." },

  // --- structure detail columns ---
  { key: "af_n_domains",         label: "AF domains",      type: "int",   default: false,
    desc: "Number of structural domains segmented in the AlphaFold model." },
  { key: "af_is_multidomain",    label: "Multidomain",     type: "bool",  default: false,
    desc: "AlphaFold model contains more than one domain." },
  { key: "pdb_n_structures",     label: "# PDB",           type: "int",   default: false,
    desc: "Number of experimental PDB structures mapped to this protein." },
  { key: "pdb_best_resolution_A", label: "PDB res. (Å)",   type: "num",   default: false,
    desc: "Best experimental structure resolution in Ångström (lower = sharper)." },

  // --- cross-species detail columns ---
  { key: "ec_transfer_essential", label: "Ess. in E. coli", type: "bool", default: false,
    desc: "Essential in E. coli (via ortholog transfer) — supporting evidence for broad-spectrum essentiality." },
  { key: "n_ecoli_orthologs",    label: "# Ec orthologs",  type: "int",   default: false,
    desc: "Number of E. coli orthologs mapped to this protein." },

  // --- novelty / popularity (studiedness) ---
  { key: "popularity_tier",      label: "Studiedness",     type: "tier",  default: false, group: "Novelty",
    desc: "Bibliometric tier: dark (neglected) / studied / well_studied. Combines the protein's own publications with those of its orthologs." },
  { key: "own_n_pubs",           label: "# pubs",          type: "int",   default: false,
    desc: "Number of publications for this protein itself (EuropePMC + UniProt), before ortholog propagation." },
  { key: "best_homolog_gene",    label: "Best homolog",    type: "text",  default: false,
    desc: "Best-studied ortholog whose literature is propagated to this protein." },
  { key: "best_homolog_organism", label: "Homolog sp.",    type: "text",  default: false,
    desc: "Species of the best-studied ortholog." },
];

// ---- table VIEWS: column presets per axis ---------------------------------
// Each view fixes which columns the table shows (Target + Composite are always
// present). Selecting a view applies its column set; the Columns menu can then
// fine-tune. The special "map" view shows the projection scatter instead.
const TABLE_VIEWS = [
  { key: "overview", label: "Overview", cols: [
    "comp_essentiality", "comp_breadth", "essentiality_tier",
    "comp_ligandability", "structure_score", "ligandability_tier",
    "comp_novelty", "popularity_tier",
    "conservation_score", "comp_human_selective", "selectivity" ] },
  { key: "essentiality", label: "Essentiality", accent: true, cols: [
    "comp_essentiality", "essentiality_tier", "experimentally_essential",
    "evidence_experimental", "evidence_transfer", "evidence_predictor",
    "proteomelm_ess_score", "geptop_score", "fba_essential", "entero_pct_essential" ] },
  { key: "ligandability", label: "Ligandability", cols: [
    "comp_ligandability", "ligandability_tier", "has_hard_evidence",
    "evidence_binding", "evidence_structural", "evidence_pocket",
    "chembl_any_n_potent", "chembl_any_best_pchembl", "fpocket_max_drug_score",
    "pocket_consensus_score", "alphafill_best_ligand", "disorder_frac" ] },
  { key: "structure", label: "Structure", cols: [
    "af_mean_plddt", "af_n_domains", "af_is_multidomain",
    "pdb_has_structure", "pdb_n_structures", "pdb_best_resolution_A",
    "disorder_frac", "fpocket_max_drug_score" ] },
  { key: "crossspecies", label: "Cross-species", cols: [
    "comp_breadth", "comp_human_selective", "selectivity", "entero_pct_essential",
    "ec_transfer_essential", "n_ecoli_orthologs", "family" ] },
  { key: "novelty", label: "Novelty", cols: [
    "comp_novelty", "popularity_tier", "own_n_pubs", "best_homolog_gene",
    "best_homolog_organism", "comp_ligandability", "ligandability_tier" ] },
  { key: "annotation", label: "Annotation", cols: [
    "functional_class", "n_interpro_entries", "family",
    "pdb_has_structure", "af_is_multidomain", "af_mean_plddt", "disorder_frac" ] },
];

// tooltips for the fixed leading columns rendered by app.js
const LEADING_COL_DESC = {
  rank: "Row rank under the current composite weighting, filters and sort.",
  name: "Gene name (italic) and UniProt accession — the canonical identifier for the protein.",
  __c: "Weighted composite score (0–1): the per-protein mean of the enabled weight sliders at left, renormalised over the components that have a value for this protein. This is the ranking you tune.",
};

// ---- categorical filters (checkbox groups) --------------------------------
const CATEGORICAL_FILTERS = [
  { key: "functional_class",   label: "Functional class", meta: FC_BY_ID,
    order: FUNCTIONAL_CLASSES.map((c) => c.id) },
  { key: "essentiality_tier",  label: "Essentiality tier" },
  { key: "ligandability_tier", label: "Ligandability tier" },
  { key: "selectivity",        label: "Selectivity" },
  { key: "popularity_tier",    label: "Studiedness" },
];
// compact abbreviations for categorical values (columns are fixed-width; the
// full value is shown on hover via the cell title + the header tooltip).
const CAT_ABBREV = {
  essentiality_tier:  { essential: "Es", likely_essential: "Lk", non_essential: "No" },
  ligandability_tier: { tractable: "Tr", partial: "Pa", intractable: "In" },
  selectivity:        { broad_selective: "BS", narrow_selective: "NS",
                        broad_human_homolog: "BH", narrow_human_homolog: "NH" },
  popularity_tier:    { dark: "Dk", studied: "St", well_studied: "Ws" },
};

// boolean filters (tri-state: any / yes / no)
const BOOL_FILTERS = [
  { key: "has_hard_evidence",        label: "Hard ligand evidence" },
  { key: "experimentally_essential", label: "Experimentally essential" },
  { key: "pdb_has_structure",        label: "Has PDB structure" },
];

// ---- detail drawer field groups -------------------------------------------
const DETAIL_GROUPS = [
  { title: "Essentiality", fields: [
    ["essentiality_score", "Score", "score"], ["essentiality_tier", "Tier", "text"],
    ["evidence_experimental", "Experimental", "score"], ["evidence_transfer", "Transfer", "score"],
    ["evidence_predictor", "Predictor", "score"], ["experimentally_essential", "Exp. essential", "bool"],
    ["entero_pct_essential", "Enterobacteriaceae % ess.", "pct"], ["bacteria_pct_essential", "Bacteria % ess.", "pct"],
    ["proteomelm_ess_score", "ProteomeLM", "score"], ["geptop_score", "Geptop", "score"],
    ["fba_essential", "FBA essential", "bool"], ["evidence_sources", "Sources", "text"],
  ]},
  { title: "Ligandability", fields: [
    ["ligandability_score", "Score", "score"], ["ligandability_tier", "Tier", "text"],
    ["evidence_binding", "Binding", "score"], ["evidence_structural", "Structural", "score"],
    ["evidence_pocket", "Pocket", "score"], ["has_hard_evidence", "Hard evidence", "bool"],
    ["human_ligandable_family", "Human-ligandable family", "bool"], ["disorder_frac", "Disorder frac.", "pct"],
    ["chembl_any_n_compounds", "ChEMBL compounds", "int"], ["chembl_any_n_potent", "ChEMBL potent (≤1µM)", "int"],
    ["chembl_any_best_pchembl", "Best pChEMBL", "num"], ["bindingdb_any_n_potent", "BindingDB potent", "int"],
    ["pdb_lig_direct_pdb_ids", "PDB co-crystals", "text"], ["alphafill_best_ligand", "AlphaFill ligand", "text"],
    ["fpocket_max_drug_score", "fpocket drug score", "score"], ["p2rank_top_score", "P2Rank top", "num"],
    ["pocket_consensus_score", "Pocket consensus", "score"],
  ]},
  { title: "Structure", fields: [
    ["af_available", "AlphaFold model", "bool"], ["af_mean_plddt", "Mean pLDDT", "plddt"],
    ["af_n_domains", "Domains", "int"], ["af_is_multidomain", "Multidomain", "bool"],
    ["pdb_has_structure", "PDB structure", "bool"], ["pdb_n_structures", "# PDB", "int"],
    ["pdb_ids", "PDB IDs", "text"], ["pdb_best_resolution_A", "Best resolution (Å)", "num"],
    ["pdb_best_method", "Method", "text"], ["pdb_has_holo", "Has holo", "bool"],
  ]},
  { title: "Cross-species / selectivity", fields: [
    ["selectivity", "Selectivity class", "text"], ["comp_human_selective", "Human-selective", "score"],
    ["conservation_score", "Conservation (ortholog spread)", "score"],
    ["ec_transfer_essential", "E. coli ess. transfer", "bool"], ["n_ecoli_orthologs", "# E. coli orthologs", "int"],
    ["family", "ESM-C family cluster", "int"],
  ]},
  { title: "Novelty / studiedness", fields: [
    ["comp_novelty", "Novelty (1 − studied)", "score"], ["popularity_tier", "Studiedness tier", "text"],
    ["own_n_pubs", "Own publications", "int"], ["own_reviewed", "Swiss-Prot reviewed", "bool"],
    ["best_homolog_gene", "Best-studied homolog", "text"], ["best_homolog_organism", "Homolog species", "text"],
    ["best_homolog_n_pubs", "Homolog publications", "int"],
  ]},
  { title: "Annotation (task-agnostic)", fields: [
    ["functional_class", "Functional class", "class"], ["interpro_family_names", "InterPro family", "text"],
    ["interpro_superfamily_names", "InterPro superfamily", "text"], ["panther_family_names", "PANTHER family", "text"],
    ["panther_subfamily_names", "PANTHER subfamily", "text"], ["pfam_ids", "Pfam", "text"],
    ["n_interpro_entries", "InterPro entries", "int"],
  ]},
];

// ---- external links (built from the accession/urls) -----------------------
function externalLinks(row) {
  const acc = row.uniprot_accession;
  const links = [
    { label: "UniProt", href: `https://www.uniprot.org/uniprotkb/${acc}/entry` },
    { label: "AlphaFold", href: `https://alphafold.ebi.ac.uk/entry/${acc}` },
    { label: "InterPro", href: `https://www.ebi.ac.uk/interpro/protein/UniProt/${acc}/` },
  ];
  if (row.pdb_best_id || (row.pdb_ids && String(row.pdb_ids).trim())) {
    const pid = row.pdb_best_id || String(row.pdb_ids).split(/[;, ]/)[0];
    if (pid) links.push({ label: `PDB ${pid}`, href: `https://www.rcsb.org/structure/${pid}` });
  }
  if (row.pfam_ids && String(row.pfam_ids).trim()) {
    const pf = String(row.pfam_ids).split(/[;, ]/)[0];
    if (pf) links.push({ label: pf, href: `https://www.ebi.ac.uk/interpro/entry/pfam/${pf}/` });
  }
  return links;
}

// ---- projection map: which 0–1 fields you can colour points by ------------
// The map plots every protein at (tsne_x, tsne_y) from the ESM-C 2D projection.
const MAP_COLORS = [
  { key: "__c", label: "Composite" },
  { key: "comp_essentiality", label: "Essentiality" },
  { key: "comp_ligandability", label: "Ligandability" },
  { key: "comp_breadth", label: "Essential breadth" },
  { key: "comp_human_selective", label: "Human-selective" },
  { key: "comp_novelty", label: "Novelty / neglect" },
  { key: "conservation_score", label: "Conservation" },
  { key: "structure_score", label: "Structure" },
];

// ---- per-page filter bars (main area) -------------------------------------
// Each view shows only its relevant filters in a bar above the table.
const VIEW_FILTERS = {
  overview:     { presets: true },
  essentiality: { cats: ["essentiality_tier"], bools: ["experimentally_essential"] },
  ligandability:{ cats: ["ligandability_tier"], bools: ["has_hard_evidence"] },
  structure:    { bools: ["pdb_has_structure"] },
  crossspecies: { cats: ["selectivity"] },
  novelty:      { cats: ["popularity_tier"] },
  annotation:   { cats: ["functional_class"], family: true },
  map:          {},
};

// ---- tier-based selection -------------------------------------------------
// An alternative to the weighted composite: pick a band (Low / Med / High) per
// axis and filter targets to those bands. Binary axes use two labelled bands.
const TIER_AXES = [
  { key: "comp_essentiality",    label: "Essentiality" },
  { key: "comp_ligandability",   label: "Ligandability" },
  { key: "comp_breadth",         label: "Essential breadth" },
  { key: "comp_human_selective", label: "Human-selective", binary: true, hi: "Selective", lo: "Human homolog" },
  { key: "comp_novelty",         label: "Novelty / neglect" },
];
// bands are [lo, hi) on the 0–1 component value; high is inclusive of 1.
const TIER_BANDS = [
  { id: "low",  label: "Low",  lo: -0.0001, hi: 0.334, intensity: 22 },
  { id: "med",  label: "Med",  lo: 0.334,   hi: 0.667, intensity: 52 },
  { id: "high", label: "High", lo: 0.667,   hi: 1.0001, intensity: 90 },
];
const TIER_BANDS_BINARY = [
  { id: "low",  label: "lo", lo: -0.0001, hi: 0.5,    intensity: 22 },
  { id: "high", label: "hi", lo: 0.5,     hi: 1.0001, intensity: 90 },
];
function tierBands(axis) {
  if (!axis.binary) return TIER_BANDS;
  return TIER_BANDS_BINARY.map((b) => ({ ...b, label: b.id === "high" ? axis.hi : axis.lo }));
}

const ORGANISM_META = {
  kp: { file: "data/kp.json", name: "K. pneumoniae", italic: "K. pneumoniae", strain: "HS11286", cssAccent: "var(--kp)" },
  ec: { file: "data/ec.json", name: "E. coli",       italic: "E. coli",       strain: "K-12",     cssAccent: "var(--ec)" },
};
