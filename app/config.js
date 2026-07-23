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
  degradability: "var(--amber)",      // Clp-protease / degron susceptibility
  orthology:     "var(--turquoise)",  // conservation + human-selectivity + selectivity class
  novelty:       "var(--orchid)",
};
// map a column/component key to its axis (for colouring)
function axisOf(key) {
  if (key === "__c") return "composite";
  if (/degrad|degron|clp_trap|halflife|clp_class/.test(key)) return "degradability";
  if (/conservation|human_selective|selectivity|human_homolog|ortholog|closeness/.test(key)) return "orthology";
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
  Degradability: "var(--amber)",
  Novelty:       "var(--orchid)",
  Orthology:     "var(--turquoise)",
  Annotation:    "var(--periwinkle)",
};
// axis → section-band label (fallback when a column has no explicit `group`)
const AXIS_GROUP_LABEL = {
  composite: "", essentiality: "Essentiality", ligandability: "Ligandability",
  degradability: "Degradability", orthology: "Orthology", novelty: "Novelty",
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
  { key: "comp_degradability",   label: "Degradability", short: "Deg", weight: 0,  on: false,
    help: "Susceptibility to bacterial Clp-protease degradation (degron motifs + Clp-substrate evidence) — basis for BacPROTAC / targeted degradation. K. pneumoniae: real; E. coli: provisional mock. Off by default." },
  { key: "comp_novelty",         label: "Novelty", short: "Nov", weight: 0, on: false,
    help: "1 − studiedness. High = under-studied (novel) target. Off by default; turn on to reward neglected proteins." },
];

// ---- table columns --------------------------------------------------------
// type: score(0–1 heat) | pct(0–1 as %) | plddt(0–100) | num | int | bool | tier | text
// Fixed leading columns (rank, name, composite) are rendered by app.js; these are the rest.
const TABLE_COLUMNS = [
  { key: "comp_essentiality",    label: "Ess. score",     type: "score", heat: true,  default: true, group: "Essentiality",
    desc: "Essentiality component (0–1): how required the gene is for fitness/survival. From essentiality_score = 0.40·experimental + 0.20·ortholog-transfer + 0.40·predictors (ProteomeLM / Geptop / FBA). Higher = better target." },
  { key: "comp_breadth",         label: "Ess. breadth",   type: "score", heat: true, default: true, group: "Essentiality",
    desc: "Essential breadth (0–1): fraction of Enterobacteriaceae genomes in which the gene is essential (entero_pct_essential). Higher = essential across more pathogens." },
  { key: "comp_ligandability",   label: "Lig. score",     type: "score", heat: true,  default: false, group: "Ligandability",
    desc: "Ligandability component (0–1): small-molecule tractability. From ligandability_score = 0.45·binding + 0.30·structural + 0.25·pocket, penalised by disorder. Higher = more druggable." },
  { key: "evidence_binding",     label: "Liganded",       type: "score", heat: true,  default: false, group: "Ligandability",
    desc: "Liganded (0–1): strength of measured/known-ligand evidence (ChEMBL / BindingDB actives + PDB co-crystals), direct or via ortholog. High = it has been liganded." },
  { key: "evidence_pocket",      label: "Ligandable",     type: "score", heat: true,  default: false, group: "Ligandability",
    desc: "Ligandable (0–1): pocket druggability — pLDDT-weighted fpocket/P2Rank consensus that a druggable pocket exists on the AlphaFold model." },
  { key: "structure_score",      label: "Structure",      type: "score", heat: true,  default: false, group: "Ligandability",
    desc: "Alias of the pocket-druggability score (evidence_pocket)." },
  { key: "conservation_score",   label: "Conservation",   type: "score", heat: true,  default: false, group: "Orthology",
    desc: "Conservation (0–1): fraction of the ~24-species bacterial panel that carries an ortholog (presence/phyletic spread). Distinct from Essential breadth — a protein can be widely conserved yet not broadly essential." },
  { key: "human_closeness",      label: "Human closeness", type: "score", heat: true,  default: false, group: "Orthology",
    desc: "Human closeness (0–1): % sequence identity to the nearest human ortholog / 100. 0 = no human ortholog (maximally selective); high = close to a human protein (selectivity risk)." },
  { key: "comp_degradability",   label: "Deg. score",     type: "score", heat: true,  default: true, group: "Degradability",
    desc: "Degradability (0–1): susceptibility to bacterial Clp-protease degradation — degron motifs (C/N-terminal) modulated by structure, plus Clp-substrate/half-life evidence. High = better BacPROTAC / targeted-degradation candidate. K. pneumoniae: real; E. coli: provisional mock." },
  { key: "degradability_tier",   label: "Deg. tier",      type: "tier",  default: true, group: "Degradability",
    desc: "Discrete degradability call: high (labile) / medium (moderate) / low (stable), thresholded on the degradability score." },
  { key: "degron_score",         label: "Degron",         type: "score", heat: true,  default: false, group: "Degradability",
    desc: "Degron-motif strength (0–1): combined C-terminal (ssrA/MuA-like) and N-terminal (N-end rule / Flynn NM) degron features, downweighted when the terminus is buried (high pLDDT)." },
  { key: "degron_cterm",         label: "C-term degron",  type: "bool",  default: false, group: "Degradability",
    desc: "A C-terminal ssrA-like or MuA-like Clp-recognition degron motif is present." },
  { key: "degron_nterm",         label: "N-term degron",  type: "bool",  default: false, group: "Degradability",
    desc: "An N-terminal destabilising residue (N-end rule) or Flynn NM degron motif is present." },
  { key: "clp_trapped",          label: "Clp substrate",  type: "bool",  default: false, group: "Degradability",
    desc: "The E. coli ortholog is a documented Clp-protease substrate (trap experiments), transferred by orthology." },
  { key: "halflife_class",       label: "Half-life",      type: "text",  default: false, group: "Degradability",
    desc: "E. coli half-life class (fast / slow / stable) transferred via ortholog (Nagar 2021 pulsed-SILAC)." },
  { key: "degron_feature_count", label: "Degron feats",   type: "int",   default: false, group: "Degradability",
    desc: "Number of degron features detected (C-terminal + N-terminal motifs)." },
  { key: "comp_novelty",         label: "Novelty",        type: "score", heat: true,  default: false, group: "Novelty",
    desc: "Novelty / neglectedness (0–1) = 1 − bibliometric studiedness. High = under-studied protein — a more novel target." },
  { key: "comp_human_selective", label: "Selective",      type: "binary01", default: true, group: "Orthology",
    desc: "Selective (yes/no): yes = no meaningful human ortholog (safer antibacterial target); no = has a human homolog." },
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

  // --- ligandability-axis detail columns (evidence_binding/pocket defined above as Liganded/Ligandable) ---
  { key: "evidence_structural",  label: "Lig: structural", type: "score", heat: true,  default: false,
    desc: "Structural evidence (0–1): drug-like ligands seen in PDB co-crystals or transplanted by AlphaFill." },
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
    "evidence_binding", "evidence_pocket", "ligandability_tier",
    "comp_degradability", "degradability_tier",
    "comp_novelty", "popularity_tier",
    "conservation_score", "human_closeness", "comp_human_selective" ] },
  { key: "essentiality", label: "Essentiality", accent: true, cols: [
    "comp_essentiality", "essentiality_tier", "experimentally_essential",
    "evidence_experimental", "evidence_transfer", "evidence_predictor",
    "proteomelm_ess_score", "geptop_score", "fba_essential", "entero_pct_essential" ] },
  { key: "ligandability", label: "Ligandability", cols: [
    "comp_ligandability", "ligandability_tier", "has_hard_evidence",
    "evidence_binding", "evidence_structural", "evidence_pocket",
    "chembl_any_n_potent", "chembl_any_best_pchembl", "fpocket_max_drug_score",
    "pocket_consensus_score", "alphafill_best_ligand", "disorder_frac" ] },
  { key: "degradability", label: "Degradability", cols: [
    "comp_degradability", "degradability_tier", "degron_score", "degron_feature_count",
    "degron_cterm", "degron_nterm", "clp_trapped", "halflife_class" ] },
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
  { key: "degradability_tier", label: "Degradability tier" },
  { key: "selectivity",        label: "Selectivity" },
  { key: "popularity_tier",    label: "Studiedness" },
];
// compact abbreviations for categorical values (columns are fixed-width; the
// full value is shown on hover via the cell title + the header tooltip).
const CAT_ABBREV = {
  essentiality_tier:  { essential: "Es", likely_essential: "Lk", non_essential: "No" },
  ligandability_tier: { tractable: "Tr", partial: "Pa", intractable: "In" },
  degradability_tier: { high: "Hi", medium: "Me", low: "Lo" },
  selectivity:        { broad_selective: "BS", narrow_selective: "NS",
                        broad_human_homolog: "BH", narrow_human_homolog: "NH" },
  popularity_tier:    { dark: "Dk", studied: "St", well_studied: "Ws" },
};

// boolean filters (tri-state: any / yes / no)
const BOOL_FILTERS = [
  { key: "has_hard_evidence",        label: "Hard ligand evidence" },
  { key: "experimentally_essential", label: "Experimentally essential" },
  { key: "pdb_has_structure",        label: "Has PDB structure" },
  { key: "clp_trapped",              label: "Clp substrate (E. coli)" },
];

// ---- gene card (detail drawer) — per-axis panel model ---------------------
// A rich profile: each axis panel has an optional headline score, a tier badge,
// evidence bars (0–1), stat chips (counts / numbers) and boolean flags. Keys
// absent from the loaded data are skipped automatically by the renderer.
const CARD_AXES = [
  { key: "essentiality", title: "Essentiality", axis: "essentiality",
    headline: "comp_essentiality", tier: "essentiality_tier",
    blurb: "How required the gene is for survival and fitness.",
    bars: [
      ["evidence_experimental", "Experimental (Tn-seq / CRISPRi)"],
      ["evidence_transfer",     "Cross-species transfer"],
      ["evidence_predictor",    "Predictor consensus"],
    ],
    subbars: [
      ["proteomelm_ess_score", "ProteomeLM"],
      ["geptop_score",         "Geptop"],
    ],
    stats: [
      ["entero_pct_essential",   "Enterobacteriaceae", "pct"],
      ["bacteria_pct_essential", "Bacteria-wide",      "pct"],
    ],
    flags: [
      ["experimentally_essential", "Experimentally essential"],
      ["fba_essential",            "FBA-lethal knockout"],
      ["ec_transfer_essential",    "Essential in E. coli"],
    ],
    sources: "evidence_sources" },

  { key: "ligandability", title: "Ligandability", axis: "ligandability",
    headline: "ligandability_score", tier: "ligandability_tier",
    blurb: "Small-molecule tractability from ligand, structural and pocket evidence.",
    bars: [
      ["evidence_binding",    "Known ligands / binding"],
      ["evidence_structural", "Structural (co-crystal · AlphaFill)"],
      ["evidence_pocket",     "Predicted druggable pocket"],
    ],
    stats: [
      ["chembl_any_n_compounds",  "ChEMBL cmpds",    "int"],
      ["chembl_any_n_potent",     "Potent ≤1 µM",    "int"],
      ["chembl_any_best_pchembl", "Best pChEMBL",    "num"],
      ["bindingdb_any_n_potent",  "BindingDB ≤1 µM", "int"],
      ["fpocket_max_drug_score",  "fpocket",         "score"],
      ["pocket_consensus_score",  "Pocket consensus","score"],
    ],
    flags: [
      ["has_hard_evidence",       "Hard evidence (≤1 µM / co-crystal)"],
      ["human_ligandable_family", "Human-ligandable family"],
    ],
    chips: [
      ["pdb_lig_direct_pdb_ids", "Co-crystal PDBs", "pdb"],
      ["alphafill_best_ligand",  "AlphaFill ligand", "text"],
    ],
    penalty: ["disorder_frac", "Disorder (score penalty)"] },

  { key: "degradability", title: "Degradability", axis: "degradability",
    headline: "comp_degradability", tier: "degradability_tier", provisionalOrgs: ["ec"],
    blurb: "Susceptibility to bacterial Clp-protease degradation — the basis for BacPROTAC / targeted degradation.",
    bars: [ ["degron_score", "Degron-motif strength"] ],
    stats: [
      ["degron_feature_count", "Degron features", "int"],
      ["ecoli_halflife_min",   "E. coli half-life (min)", "num"],
    ],
    flags: [
      ["degron_cterm", "C-terminal degron"],
      ["degron_nterm", "N-terminal degron"],
      ["clp_trapped",  "E. coli Clp substrate"],
    ],
    text: [
      ["halflife_class",  "E. coli half-life class"],
      ["ecoli_clp_class", "Clp system"],
    ] },

  { key: "structure", title: "Structure", axis: "ligandability",
    plddt: "af_mean_plddt",
    blurb: "AlphaFold model confidence and experimental coverage.",
    flags: [
      ["af_available",      "AlphaFold model"],
      ["af_is_multidomain", "Multidomain fold"],
      ["pdb_has_structure", "Experimental PDB"],
      ["pdb_has_holo",      "Ligand-bound (holo)"],
    ],
    stats: [
      ["af_n_domains",         "Domains",        "int"],
      ["pdb_n_structures",     "# PDB",          "int"],
      ["pdb_best_resolution_A","Best res. (Å)",  "num"],
    ],
    chips: [ ["pdb_ids", "PDB structures", "pdb"] ],
    text:  [ ["pdb_best_method", "Best method"] ] },

  { key: "selectivity", title: "Selectivity & conservation", axis: "orthology",
    tierClass: "selectivity",
    blurb: "Human off-target risk versus phyletic spread across bacteria.",
    bars: [
      ["conservation_score", "Conservation (ortholog spread)"],
      ["human_closeness",    "Human closeness (identity) — lower is safer", true],
    ],
    stats: [
      ["n_orthologs",       "# orthologs",        "int"],
      ["n_ecoli_orthologs", "# E. coli orthologs","int"],
      ["family",            "ESM-C family",       "int"],
    ],
    flags: [
      ["comp_human_selective", "No human ortholog (selective)", "ge05"],
    ] },

  { key: "novelty", title: "Novelty & studiedness", axis: "novelty",
    headline: "comp_novelty", tier: "popularity_tier",
    blurb: "How under-studied the target is (bibliometric, ortholog-propagated).",
    stats: [
      ["own_n_pubs",          "Own publications", "int"],
      ["best_homolog_n_pubs", "Homolog pubs",     "int"],
    ],
    flags: [ ["own_reviewed", "Swiss-Prot reviewed"] ],
    homolog: ["best_homolog_gene", "best_homolog_organism"] },
];

// annotation panel (families / signatures) — rendered as chip clusters
const CARD_ANNOTATION = {
  title: "Annotation", axis: "annotation",
  classKey: "functional_class",
  chipGroups: [
    ["interpro_family_names",      "InterPro family"],
    ["interpro_superfamily_names", "InterPro superfamily"],
    ["interpro_domain_names",      "InterPro domain"],
    ["panther_family_names",       "PANTHER family"],
    ["panther_subfamily_names",    "PANTHER subfamily"],
    ["pfam_ids",                   "Pfam"],
  ],
};
const CARD_GROUP_COLORS = {
  essentiality: "var(--crimson)", ligandability: "var(--cobalt)", degradability: "var(--amber)",
  orthology: "var(--turquoise)", novelty: "var(--orchid)", annotation: "var(--periwinkle)",
};

// ---- external links (built from the accession/urls) -----------------------
function externalLinks(row) {
  const acc = row.uniprot_accession;
  const gene = (row.gene || "").trim();
  const links = [
    { label: "UniProt", icon: "🧬", href: `https://www.uniprot.org/uniprotkb/${acc}/entry` },
    { label: "AlphaFold", icon: "◈", href: `https://alphafold.ebi.ac.uk/entry/${acc}` },
    { label: "InterPro", icon: "▦", href: `https://www.ebi.ac.uk/interpro/protein/UniProt/${acc}/` },
    { label: "STRING", icon: "⚭", href: `https://string-db.org/cgi/network?identifiers=${encodeURIComponent(acc)}` },
  ];
  if (row.pdb_best_id || (row.pdb_ids && String(row.pdb_ids).trim())) {
    const pid = row.pdb_best_id || String(row.pdb_ids).split(/[;, ]/)[0];
    if (pid) links.push({ label: `PDB ${pid}`, icon: "◍", href: `https://www.rcsb.org/structure/${pid}` });
  }
  if (row.pfam_ids && String(row.pfam_ids).trim()) {
    const pf = String(row.pfam_ids).split(/[;, ]/)[0];
    if (pf) links.push({ label: pf, icon: "▦", href: `https://www.ebi.ac.uk/interpro/entry/pfam/${pf}/` });
  }
  // literature search (useful for the novelty angle) — gene + organism
  const org = (typeof state !== "undefined" && state.org) || "kp";
  const orgName = (ORGANISM_META[org] && ORGANISM_META[org].name) || "Klebsiella pneumoniae";
  const term = encodeURIComponent(`${gene ? gene + " " : ""}${orgName}`.trim());
  links.push({ label: "PubMed", icon: "🔎", href: `https://pubmed.ncbi.nlm.nih.gov/?term=${term}` });
  return links;
}

// ---- projection map: which 0–1 fields you can colour points by ------------
// The map plots every protein at (tsne_x, tsne_y) from the ESM-C 2D projection.
const MAP_COLORS = [
  { key: "__c", label: "Composite" },
  { key: "comp_essentiality", label: "Essentiality" },
  { key: "comp_ligandability", label: "Ligandability" },
  { key: "comp_degradability", label: "Degradability" },
  { key: "comp_breadth", label: "Essential breadth" },
  { key: "comp_human_selective", label: "Human-selective" },
  { key: "comp_novelty", label: "Novelty" },
  { key: "conservation_score", label: "Conservation" },
  { key: "human_closeness", label: "Human closeness" },
];

// ---- per-page filter bars (main area) -------------------------------------
// Each view shows only its relevant filters in a bar above the table.
const VIEW_FILTERS = {
  overview:     {},
  essentiality: { cats: ["essentiality_tier"], bools: ["experimentally_essential"] },
  ligandability:{ cats: ["ligandability_tier"], bools: ["has_hard_evidence"] },
  degradability:{ cats: ["degradability_tier"], bools: ["clp_trapped"] },
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
  { key: "comp_essentiality", label: "Essentiality", bands: [
      { id: "non",    label: "Non-ess.",  intensity: 22, test: (v, r) => r.essentiality_tier === "non_essential" },
      { id: "likely", label: "Likely",    intensity: 52, test: (v, r) => r.essentiality_tier === "likely_essential" },
      { id: "ess",    label: "Essential", intensity: 90, test: (v, r) => r.essentiality_tier === "essential" } ] },
  { key: "comp_breadth", label: "Essential breadth", bands: [
      { id: "low",  label: "Specific", lo: -0.0001, hi: 0.334,  intensity: 22 },
      { id: "med",  label: "Broad",    lo: 0.334,   hi: 0.667,  intensity: 52 },
      { id: "high", label: "Core",     lo: 0.667,   hi: 1.0001, intensity: 90 } ] },
  { key: "comp_ligandability", label: "Ligandability", bands: [
      { id: "none",   label: "Not ligandable",   intensity: 22,
        test: (v, r) => !r.has_hard_evidence && !(typeof r.evidence_pocket === "number" && r.evidence_pocket >= 0.5) },
      { id: "pocket", label: "Lig. pocket", intensity: 52,
        test: (v, r) => !r.has_hard_evidence && typeof r.evidence_pocket === "number" && r.evidence_pocket >= 0.5 },
      { id: "known",  label: "Known ligands",     intensity: 90,
        test: (v, r) => r.has_hard_evidence === true } ] },
  { key: "comp_degradability", label: "Degradability", bands: [
      { id: "stable",   label: "Stable",   intensity: 22, test: (v, r) => r.degradability_tier === "low" },
      { id: "moderate", label: "Moderate", intensity: 52, test: (v, r) => r.degradability_tier === "medium" },
      { id: "labile",   label: "Labile",   intensity: 90, test: (v, r) => r.degradability_tier === "high" } ] },
  { key: "comp_novelty", label: "Novelty", bands: [
      { id: "well",    label: "Well-studied", intensity: 22, test: (v, r) => r.popularity_tier === "well_studied" },
      { id: "studied", label: "Studied",      intensity: 52, test: (v, r) => r.popularity_tier === "studied" },
      { id: "dark",    label: "Dark",         intensity: 90, test: (v, r) => r.popularity_tier === "dark" } ] },
  { key: "conservation_score", label: "Conservation", bands: [
      { id: "low",  label: "Specific",     lo: -0.0001, hi: 0.334,  intensity: 22 },
      { id: "med",  label: "Intermediate", lo: 0.334,   hi: 0.667,  intensity: 52 },
      { id: "high", label: "Broad",        lo: 0.667,   hi: 1.0001, intensity: 90 } ] },
  { key: "comp_human_selective", label: "Selectivity", bands: [
      { id: "high", label: "Human selective", lo: 0.5,     hi: 1.0001, intensity: 90 },
      { id: "low",  label: "Human homologs",  lo: -0.0001, hi: 0.5,    intensity: 22 } ] },
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
function tierBands(axis) { return axis.bands || TIER_BANDS; }
// does a value/row fall in a band? numeric range by default, or a custom predicate.
function bandMatches(band, v, row) {
  if (band.test) return band.test(v, row);
  return typeof v === "number" && v >= band.lo && v < band.hi;
}

const ORGANISM_META = {
  kp: { file: "data/kp.json", name: "K. pneumoniae", italic: "K. pneumoniae", strain: "HS11286", cssAccent: "var(--kp)" },
  ec: { file: "data/ec.json", name: "E. coli",       italic: "E. coli",       strain: "K-12",     cssAccent: "var(--ec)" },
};
