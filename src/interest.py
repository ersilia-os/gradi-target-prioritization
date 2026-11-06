"""The GraDi consortium's proteins of interest, as the repository records them.

**There is no machine-readable list in this repo.** The consortium's targets are stated in prose, in
two places, and this module is the first time they are written down as data. Both sources:

1. `legacy/docs/03_degradability.md:420` — *"the ~30 consortium envelope targets (LPS biosynthesis,
   LptA-G, BamA/D, Lol/Lnt/LspA, Sec/YidC/LepB, FtsH)"*, in a passage warning that these are
   **almost all envelope**, so nearly all are reachable by activated ClpP *only* through the
   pre-export window, if at all.
2. `legacy/HISTORY.md:611` and `legacy/docs/degradability_report.md:181` — the unresolved conflict
   between the kick-off note ("targets in the periplasm") and the **v5 proposal, GyrA/GyrB,
   "cytosolic"**. Both were recorded as never settled.

So the panel below is **an expansion of prose into gene symbols**, not a file someone handed over.
Where the source names a range or a family it is expanded literally -- `LptA-G` becomes `lptA`..`lptG`,
`BamA/D` becomes `bamA` and `bamD` only (the source names those two, not the whole Bam complex) --
and where it names a pathway ("LPS biosynthesis") the canonical *lpx*/*kds* core is used, which is a
judgement and is flagged as one in `EXPANSION_NOTES`.

**Treat this as a draft to be corrected, not as the consortium's own list.** `PANEL` is a plain dict,
so replacing it with a curated file later changes nothing downstream.

Matching is by **gene symbol, case-insensitively**, because that is what the sources use. That is a
real limitation on this project's anchor: Kp `gene_name` covers 63.4% of the proteome (44.6% on Sa),
so a consortium target whose Kp ortholog has no symbol cannot be found this way and will show up in
`missing()`. Always read the coverage report before concluding a target is absent.
"""

from __future__ import annotations

import pandas as pd

# family label -> gene symbols. Order is envelope-outward-ish, then the cytosolic proposal last.
PANEL: dict[str, tuple[str, ...]] = {
    # "LPS biosynthesis" -- the source names the pathway, not the genes. This is the Raetz pathway
    # core (lipid A) plus the Kdo arm; see EXPANSION_NOTES.
    "LPS biosynthesis": (
        "lpxA", "lpxB", "lpxC", "lpxD", "lpxH", "lpxK", "lpxL", "lpxM",
        "kdsA", "kdsB", "kdsC", "kdsD", "waaA", "gmhA",
    ),
    # "LptA-G" -- expanded literally across the range.
    "LPS transport (Lpt)": ("lptA", "lptB", "lptC", "lptD", "lptE", "lptF", "lptG"),
    # "BamA/D" -- only the two the source names.
    "OM assembly (Bam)": ("bamA", "bamD"),
    # "Lol/Lnt/LspA" -- the lipoprotein sorting pathway.
    "Lipoprotein sorting": ("lolA", "lolB", "lolC", "lolD", "lolE", "lnt", "lspA"),
    # "Sec/YidC/LepB" -- the general secretion route plus insertase and signal peptidase I.
    "Secretion (Sec)": (
        "secA", "secB", "secD", "secE", "secF", "secG", "secY", "yajC", "yidC", "lepB",
    ),
    # Named on its own in the source, and separately flagged there as one of only three
    # experimentally essential proteases -- i.e. both a target and a candidate handle.
    "FtsH": ("ftsH",),
    # The v5 proposal, explicitly "cytosolic" and therefore the only part of this panel that
    # activated ClpP could reach directly.
    "v5 proposal (cytosolic)": ("gyrA", "gyrB"),
}

EXPANSION_NOTES: dict[str, str] = {
    "LPS biosynthesis": (
        "The source names a pathway, not genes. Expanded to the lipid A (Raetz) core plus the Kdo "
        "arm and waaA/gmhA. This is the least faithful entry in the panel -- a curated list would "
        "likely differ at the edges."
    ),
    "LPS transport (Lpt)": "'LptA-G' expanded literally as lptA..lptG.",
    "OM assembly (Bam)": (
        "'BamA/D' expanded to bamA and bamD only. bamB/C/E are deliberately NOT included: the "
        "source names two subunits, and adding the rest would be inventing scope."
    ),
    "Lipoprotein sorting": "'Lol/Lnt/LspA' expanded to lolA..lolE plus lnt and lspA.",
    "Secretion (Sec)": (
        "'Sec/YidC/LepB' expanded to the secA/B/D/E/F/G/Y core (plus yajC, which forms the "
        "SecDF-YajC complex) with yidC and lepB."
    ),
    "FtsH": "Named directly.",
    "v5 proposal (cytosolic)": "Named directly: GyrA/GyrB.",
}

# Which families the source itself describes as envelope, hence reachable by activated ClpP only
# through the pre-export window. Everything here is a caveat, not a prediction.
ENVELOPE_FAMILIES: tuple[str, ...] = (
    "LPS biosynthesis", "LPS transport (Lpt)", "OM assembly (Bam)",
    "Lipoprotein sorting", "Secretion (Sec)", "FtsH",
)

SOURCES = (
    "legacy/docs/03_degradability.md:420 (the ~30 consortium envelope targets); "
    "legacy/HISTORY.md:611 and legacy/docs/degradability_report.md:181 (the v5 GyrA/GyrB proposal)"
)


def genes() -> list[str]:
    """Every gene symbol in the panel, deduplicated, in panel order."""
    seen: dict[str, None] = {}
    for members in PANEL.values():
        for g in members:
            seen.setdefault(g, None)
    return list(seen)


def family_of() -> dict[str, str]:
    """`{gene_symbol_lower: family label}`."""
    return {g.lower(): fam for fam, members in PANEL.items() for g in members}


def annotate(df: pd.DataFrame, gene_col: str = "gene_name") -> pd.DataFrame:
    """Add `interest_family` (the panel family, or empty) and `is_interest` to a frame.

    Matches on lowercased gene symbol. A row whose `gene_col` is empty can never match, which is the
    coverage limitation described in the module docstring -- not evidence of absence.
    """
    fam = family_of()
    key = df[gene_col].fillna("").astype(str).str.strip().str.lower()
    out = df.copy()
    out["interest_family"] = key.map(fam).fillna("")
    out["is_interest"] = out["interest_family"] != ""
    return out


def coverage(df: pd.DataFrame, gene_col: str = "gene_name") -> pd.DataFrame:
    """Per family: how many panel genes were found in this frame, and which are missing.

    Print this before reading any panel figure. A missing gene means "no protein in this proteome
    carries that symbol", which on Kp and Sa is often a naming gap rather than a biological absence.
    """
    present = set(df[gene_col].fillna("").astype(str).str.strip().str.lower())
    rows = []
    for fam, members in PANEL.items():
        found = [g for g in members if g.lower() in present]
        missing = [g for g in members if g.lower() not in present]
        rows.append({
            "family": fam,
            "n_panel": len(members),
            "n_found": len(found),
            "found": ",".join(found),
            "missing": ",".join(missing),
            "is_envelope": fam in ENVELOPE_FAMILIES,
        })
    return pd.DataFrame(rows)


def missing(df: pd.DataFrame, gene_col: str = "gene_name") -> list[str]:
    """Panel genes with no matching protein in this frame."""
    present = set(df[gene_col].fillna("").astype(str).str.strip().str.lower())
    return [g for g in genes() if g.lower() not in present]
