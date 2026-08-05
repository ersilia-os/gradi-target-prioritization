"""Shared helpers for the subcellular-localization stage (docs §5.1, scripts/09*).

Deliberately thin, in the same spirit as `src/essentiality.py`: it **reuses**
`src/ligandability.py` for everything generic — organism config, proteome loading, the 03a
ortholog table, the DIAMOND-by-sequence engine and the output-path helpers — and adds only the
localization-specific pieces:

  * `LOC_CLASSES` — the one canonical class vocabulary every 09* track normalises onto, plus the
    per-source translation tables (`UNIPROT_*`, `DEEPLOCPRO_MAP`, `PSORTB_MAP`, `STEPDB_MAP`);
  * `classify_uniprot()` — keyword parse of the UniProt `cc_subcellular_location` free text, with
    a fall-back to the signal-peptide / transmembrane features when the text is silent;
  * `parse_tmbed()` — turns TMbed's per-residue `S/H/B/i/o` label string into the per-protein
    topology features the accessibility ladder needs (β-barrel call, cytoplasmic-domain fraction);
  * `lipoprotein_sorting()` — the Lol pathway **"+2 rule"**, which is what separates an
    outer-membrane lipoprotein from an inner-membrane one (PSORTb's best-known failure mode);
  * `clp_accessibility()` — the graded [0-1] ladder specified in
    `docs/05_expression_and_localization.md`, replacing the earlier coarse 1.0/0.5/0.0 stub.

Keying is always by UniProt accession (project convention). Run with the `gradi` conda env, except
the two predictor scripts (09c DeepLocPro, 09d TMbed) which need `gradi-loc` — see install.sh.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src import ligandability as L  # noqa: E402
from src import essentiality as E  # noqa: E402

# Re-export the generics so 09* scripts can reach everything through `from src import localization as LOC`.
ORGANISMS = L.ORGANISMS
REPO_ROOT = L.REPO_ROOT
acc_from_header = L.acc_from_header
load_accessions = L.load_accessions
load_genes = L.load_genes
load_orthologs = L.load_orthologs
ortholog_map = L.ortholog_map
proteome_fasta = L.proteome_fasta
proteome_tsv = L.proteome_tsv
results_dir = L.results_dir
processed_dir = L.processed_dir
run_diamond_blastp = L.run_diamond_blastp
load_diamond_hits = L.load_diamond_hits
transfer_ecoli_to_kp = E.transfer_ecoli_to_kp
ORG_DISPLAY = E.ORG_DISPLAY


# --------------------------------------------------------------------------- vocabulary
# The canonical class set. Every track normalises onto this; `membrane` is a *provisional*
# label for a UniProt "Membrane" call with no stated side, which 09g resolves using TMbed
# topology (β-barrel -> outer_membrane, otherwise -> inner_membrane).
LOC_CLASSES: tuple[str, ...] = (
    "cytoplasm",
    "inner_membrane",
    "periplasm",
    "outer_membrane",
    "extracellular",
    "cell_wall_surface",
    "membrane",
    "unknown",
)

# Short labels for plots and the webapp badge.
LOC_ABBREV: dict[str, str] = {
    "cytoplasm": "Cyt",
    "inner_membrane": "IM",
    "periplasm": "Peri",
    "outer_membrane": "OM",
    "extracellular": "Ext",
    "cell_wall_surface": "CW",
    "membrane": "Mem",
    "unknown": "?",
}

# DeepLocPro 1.0 emits these six strings (DeepLocPro/deeplocpro.py `labels`).
DEEPLOCPRO_MAP: dict[str, str] = {
    "cytoplasmic": "cytoplasm",
    "cytoplasmic membrane": "inner_membrane",
    "periplasmic": "periplasm",
    "outer membrane": "outer_membrane",
    "extracellular": "extracellular",
    "cell wall & surface": "cell_wall_surface",
}

# PSORTb 3.0 `Final_Localization` values (as served by PSORTdb).
PSORTB_MAP: dict[str, str] = {
    "cytoplasmic": "cytoplasm",
    "cytoplasmicmembrane": "inner_membrane",
    "cytoplasmic membrane": "inner_membrane",
    "periplasmic": "periplasm",
    "outermembrane": "outer_membrane",
    "outer membrane": "outer_membrane",
    "extracellular": "extracellular",
    "cellwall": "cell_wall_surface",
    "cell wall": "cell_wall_surface",
}

# STEPdb 2.0 subcellular classes (E. coli), matched as ordered substrings because the published
# CSV carries compound values ("Nucleoid, Integral Inner Membrane") and free text. First match
# wins, so the more specific patterns must come first.
#
# The important modelling choice: STEPdb resolves *which side* a peripheral membrane protein faces,
# and a peripheral inner-membrane protein facing the cytoplasm is as Clp-reachable as a soluble
# cytoplasmic one — so it folds to `cytoplasm`, not `inner_membrane`. Nucleoid and ribosomal
# proteins are likewise cytoplasmic.
STEPDB_PATTERNS: tuple[tuple[str, str], ...] = (
    ("peripheral inner membrane protein facing the periplasm", "periplasm"),
    ("peripheral inner membrane protein facing the cytoplasm", "cytoplasm"),
    ("peripheral outer membrane protein facing the extra-cellular", "extracellular"),
    ("outer membrane lipoprotein", "outer_membrane"),
    ("outer membrane b-barrel", "outer_membrane"),
    ("outer membrane", "outer_membrane"),
    ("inner membrane lipoprotein", "inner_membrane"),
    ("integral inner membrane", "inner_membrane"),
    ("inner membrane", "inner_membrane"),
    ("periplasmic", "periplasm"),
    ("extracellular", "extracellular"),
    ("cell surface", "cell_wall_surface"),
    ("nucleoid", "cytoplasm"),
    ("ribosomal", "cytoplasm"),
    ("cytoplasmic", "cytoplasm"),
)


def classify_stepdb(value: str) -> str:
    """Map a STEPdb subcellular-location string onto `LOC_CLASSES`."""
    s = (value or "").strip().lower()
    if not s:
        return "unknown"
    for pattern, cls in STEPDB_PATTERNS:
        if pattern in s:
            return cls
    return "unknown"

# ECO codes that mark a UniProt annotation as experimentally determined (rather than inferred by
# similarity / sequence model / automatic assertion). 0000269 = experimental evidence used in a
# manual assertion; 0007744 = combinatorial evidence used in a manual assertion.
UNIPROT_EXPERIMENTAL_ECO: frozenset[str] = frozenset({"ECO:0000269", "ECO:0007744"})
_ECO_RE = re.compile(r"ECO:\d{7}")


# --------------------------------------------------------------------------- paths
def localization_raw_dir(organism: str, *parts: str) -> Path:
    d = REPO_ROOT / "data" / "raw" / organism / "localization"
    for p in parts:
        d = d / p
    d.mkdir(parents=True, exist_ok=True)
    return d


def localization_processed_dir(organism: str, *parts: str) -> Path:
    return L.processed_dir(organism, "localization", *parts)


# --------------------------------------------------------------------------- sequences
def read_fasta(path: Path) -> dict[str, str]:
    """UniProt-style FASTA -> {accession: sequence}, keyed with `L.acc_from_header`."""
    seqs: dict[str, str] = {}
    acc, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(chunks)
                acc, chunks = L.acc_from_header(line[1:]), []
            elif acc is not None:
                chunks.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(chunks)
    return seqs


def write_fasta(seqs: dict[str, str], path: Path, width: int = 60) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for acc, seq in seqs.items():
            fh.write(f">{acc}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    return path


# --------------------------------------------------------------------------- UniProt parsing
def uniprot_evidence(subcell_raw: str) -> str:
    """`experimental` if any experimental ECO code is cited, `curated` if annotated at all."""
    if not subcell_raw or not subcell_raw.strip():
        return "none"
    codes = set(_ECO_RE.findall(subcell_raw))
    if codes & UNIPROT_EXPERIMENTAL_ECO:
        return "experimental"
    return "curated"


def classify_uniprot(subcell: str, has_signal: bool, n_tm: int) -> str:
    """Keyword-parse `cc_subcellular_location` onto `LOC_CLASSES`.

    Order matters: the more specific compartment wins, because UniProt CC text routinely chains
    several ("Cell inner membrane; Periplasmic side"). Falls back to the sequence features when
    the CC text is absent, and to `unknown` when there is nothing to go on at all.
    """
    s = (subcell or "").lower()

    def has(*keys: str) -> bool:
        return any(k in s for k in keys)

    if has("outer membrane", "outer-membrane"):
        return "outer_membrane"

    # Peripheral membrane proteins are placed by the side they sit on, not by the membrane.
    # This matters directly for Clp accessibility — a peripherally-attached protein on the
    # cytoplasmic face is as reachable as a soluble cytoplasmic one — and it keeps the UniProt
    # parse consistent with `classify_stepdb`, which draws the same distinction explicitly.
    # It also fixes chaperone-type entries such as DnaK, whose CC reads
    # "Cytoplasm. Cell inner membrane; Peripheral membrane protein."
    if has("peripheral membrane protein"):
        if has("cytoplasmic side"):
            return "cytoplasm"
        if has("periplasmic side"):
            return "periplasm"
        if has("cytoplasm", "cytosol"):
            return "cytoplasm"

    if has("periplasm"):
        return "periplasm"
    if has("fimbri", "pilus", "pili", "flagell", "capsule", "cell surface", "cell wall"):
        return "cell_wall_surface"
    if has("secreted", "extracellular"):
        return "extracellular"
    if has("inner membrane", "cytoplasmic membrane", "plasma membrane", "cell membrane"):
        return "inner_membrane"
    if has("cytoplasm", "cytosol"):
        return "cytoplasm"
    if has("membrane"):
        return "membrane"          # side unstated -> 09g resolves it from topology
    if n_tm and n_tm > 0:
        return "membrane"
    if has_signal:
        return "extracellular"     # exported, compartment unstated
    return "unknown"


# --------------------------------------------------------------------------- TMbed topology
def parse_tmbed(path: Path) -> dict[str, dict]:
    """Parse TMbed 3-line `--out-format=1` output into per-protein topology features.

    Label alphabet: `S` signal peptide, `H` transmembrane α-helix, `B` transmembrane β-strand,
    `i` non-transmembrane inside (cytoplasmic), `o` non-transmembrane outside (periplasmic).
    """
    out: dict[str, dict] = {}
    lines = [ln.rstrip("\n") for ln in open(path) if ln.strip()]
    for i in range(0, len(lines) - 2, 3):
        header, _seq, labels = lines[i], lines[i + 1], lines[i + 2]
        if not header.startswith(">"):
            continue
        out[L.acc_from_header(header[1:])] = tmbed_features(labels)
    return out


def _n_segments(labels: str, ch: str) -> int:
    return len([seg for seg in re.split(f"[^{ch}]+", labels) if seg])


def tmbed_features(labels: str) -> dict:
    """Per-protein features from one TMbed label string."""
    n = len(labels) or 1
    n_helix = _n_segments(labels, "H")
    n_strand = _n_segments(labels, "B")
    inside_runs = [seg for seg in re.split("[^i]+", labels) if seg]
    return {
        "n_tm_helix": n_helix,
        "n_tm_strand": n_strand,
        # A Gram-negative outer-membrane β-barrel needs several strands; one or two stray
        # predicted strands are noise, so require the canonical minimum of 8.
        "is_beta_barrel": bool(n_strand >= 8),
        "has_signal_peptide_tmbed": bool("S" in labels),
        "cyto_residue_fraction": round(labels.count("i") / n, 4),
        "cyto_longest_segment": max((len(r) for r in inside_runs), default=0),
        "tmbed_length": len(labels),
    }


# --------------------------------------------------------------------------- lipoprotein sorting
# Lol pathway "+2 rule": after signal-peptidase II cleavage the mature protein starts with the
# lipidated Cys. An Asp at position +2 relative to that Cys is the Lol avoidance signal and
# retains the lipoprotein in the inner membrane; anything else routes it to the outer membrane.
LOL_AVOIDANCE_RESIDUES: frozenset[str] = frozenset({"D"})


# Lipobox detection, used by 09e only when licensed SignalP 6.0 is unavailable.
# Motif is Prosite PS51257 ([LVI][ASTVI][GAS]C, the Cys becoming the mature N-terminus); on its own
# it is far too permissive, so we also require the two other defining features of a signal peptide:
# a positively charged n-region and a hydrophobic h-region immediately before the lipobox.
# Calibrated against the 99 UniProt lipid-anchor entries of E. coli K-12: with these three filters
# together, precision 0.74 / recall 0.82 (motif alone: precision 0.36).
_LIPOBOX_RE = re.compile(r"[LVI][ASTVI][GAS](C)")
LIPOBOX_SEARCH_WINDOW = 35     # signal peptides are short; the Cys is within the first ~35 residues
LIPOBOX_MIN_CYS_POS = 14       # ...and after a real n+h region, never right at the N-terminus
LIPOBOX_H_REGION = 12          # residues before the lipobox scored for hydrophobicity
LIPOBOX_MIN_HYDROPATHY = 1.0   # mean Kyte-Doolittle over the h-region
LIPOBOX_N_REGION = 8           # residues scanned for the positively charged n-region

# Kyte-Doolittle hydropathy.
KYTE_DOOLITTLE: dict[str, float] = dict(zip(
    "ARNDCQEGHILKMFPSTWYV",
    [1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5,
     3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2],
))


def mean_hydropathy(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(KYTE_DOOLITTLE.get(c, 0.0) for c in seq) / len(seq)


def find_lipobox(sequence: str) -> int | None:
    """1-based position of the lipidated Cys of a Sec/SPII signal peptide, or None."""
    m = _LIPOBOX_RE.search(sequence[:LIPOBOX_SEARCH_WINDOW])
    if not m:
        return None
    cys_idx = m.start(1)                       # 0-based index of the Cys
    if cys_idx + 1 < LIPOBOX_MIN_CYS_POS:
        return None
    h_region = sequence[max(0, cys_idx - LIPOBOX_H_REGION):cys_idx]
    if mean_hydropathy(h_region) < LIPOBOX_MIN_HYDROPATHY:
        return None
    if not any(c in "KR" for c in sequence[:LIPOBOX_N_REGION]):
        return None
    return cys_idx + 1


def lipoprotein_sorting(sequence: str, cleavage_pos: int) -> str | None:
    """Sort a Sec/SPII lipoprotein to `inner_membrane` or `outer_membrane`.

    `cleavage_pos` is the 1-based position of the lipidated cysteine (the first residue of the
    mature protein). Returns None when the sequence is too short to read position +2.
    """
    if not sequence or cleavage_pos < 1:
        return None
    plus_two_idx = cleavage_pos + 1  # 0-based index of the residue two after the Cys
    if plus_two_idx >= len(sequence):
        return None
    return "inner_membrane" if sequence[plus_two_idx] in LOL_AVOIDANCE_RESIDUES else "outer_membrane"


# --------------------------------------------------------------------------- accessibility ladder
# docs/05_expression_and_localization.md. `clp_accessibility` scores reachability by the
# *cytoplasmic* Clp machinery, which is the gating requirement for BacPROTAC-style degradation.
CYTO_DOMAIN_FRACTION_MIN = 0.30   # IM protein needs this much cytoplasm-facing sequence to score 0.6
CLP_IM_EXPOSED = 0.6
CLP_IM_BURIED = 0.2
CLP_IM_UNKNOWN_TOPOLOGY = 0.4     # documented midpoint when TMbed topology is missing
CLP_PERIPLASM = 0.2
MIN_CLP_ACCESSIBILITY = 0.5       # shortlist gate; mirrors src/degradability.py


def clp_accessibility(
    localization: str,
    cyto_residue_fraction: float | None = None,
    is_beta_barrel: bool = False,
) -> float | None:
    """Graded [0-1] Clp reachability, or None when the localization is unknown."""
    if is_beta_barrel:
        return 0.0
    if localization == "cytoplasm":
        return 1.0
    if localization in ("inner_membrane", "membrane"):
        if cyto_residue_fraction is None:
            return CLP_IM_UNKNOWN_TOPOLOGY
        return CLP_IM_EXPOSED if cyto_residue_fraction >= CYTO_DOMAIN_FRACTION_MIN else CLP_IM_BURIED
    if localization == "periplasm":
        return CLP_PERIPLASM
    if localization in ("outer_membrane", "extracellular", "cell_wall_surface"):
        return 0.0
    return None


# --------------------------------------------------------------------------- ortholog transfer
ECOLI_SPECIES_KEY = "Ecoli_K12_MG1655"


def transfer_categorical_ecoli_to_kp(ec_label_by_acc: dict[str, str]) -> dict[str, dict]:
    """Lift an E.-coli-keyed *categorical* localization onto HS11286 anchors via orthology.

    `E.transfer_ecoli_to_kp` reduces with max/mean, which is meaningless for a class label, so the
    transfer rule here is **donor consensus** instead.

    Note the 03a kp->ec orthology is OrthoFinder orthogroup membership only — `pident`, `coverage`
    and `bitscore` are all empty for these rows — so there is no identity to threshold on or to
    break ties with. Orthogroup membership is itself the evidence. Of the 3,179 Kp anchors with an
    E. coli ortholog only 81 have more than one, and for those we take the majority label and
    **drop outright ties** rather than pick arbitrarily: a wrong compartment silently becomes a
    wrong Clp-accessibility score, so abstaining is the cheaper error.

    Returns {kp_accession: {label, donors, n_donors, agreement}}.
    """
    orth = L.load_orthologs("kpneumoniae")
    sub = orth[orth["species"] == ECOLI_SPECIES_KEY].copy()
    sub["label"] = sub["target_uniprot"].map(ec_label_by_acc)
    sub = sub.dropna(subset=["label"])
    if sub.empty:
        return {}

    out: dict[str, dict] = {}
    for anchor, grp in sub.groupby("anchor_uniprot"):
        counts = grp["label"].value_counts()
        if len(counts) > 1 and counts.iloc[0] == counts.iloc[1]:
            continue  # donors disagree with no majority -> abstain
        label = counts.index[0]
        out[anchor] = {
            "label": label,
            "donors": ";".join(sorted(grp.loc[grp["label"] == label, "target_uniprot"])),
            "n_donors": int(len(grp)),
            "agreement": round(counts.iloc[0] / len(grp), 3),
        }
    return out


# --------------------------------------------------------------------------- merge bookkeeping
# 09g precedence, best first. `localization_source` records which one supplied the final call.
SOURCE_PRECEDENCE: tuple[str, ...] = (
    "uniprot_experimental",
    "stepdb",
    "ortholog_transfer",
    "uniprot_curated",
    "deeplocpro",
    "psortb",
)

# Evidence tier reported alongside the call.
SOURCE_EVIDENCE: dict[str, str] = {
    "uniprot_experimental": "experimental",
    "stepdb": "experimental",
    "ortholog_transfer": "experimental",
    "uniprot_curated": "curated",
    "deeplocpro": "predicted",
    "psortb": "predicted",
    "none": "none",
}
