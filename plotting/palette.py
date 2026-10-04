"""The deck's colours: stylia's **npg** palette, NOT stylia's ersilia palette.

**Standing instruction from the project owner (2026-10-04):** use stylia, but take colours from
`stylia.CategoricalPalette("npg")` -- ggsci's Nature Publishing Group palette -- rather than from
`stylia.NamedColors()`, which under `set_style("ersilia")` returns the ersilia plum/orange/mint set.

`set_style` and the palette are DIFFERENT KNOBS and only the second one is constrained here.
`set_style` controls typography, grid and spines; the colours come from this module. Every script
still calls `stylia.set_format("slide")` and `stylia.set_style(...)` for the former.

**This is the one place the deck's colours are defined, and that is a deliberate departure from
`scripts/plots/`,** where style constants are copy-pasted into each script by design. The reason
the exception is justified: a palette is no longer per-script taste once the owner has set it as a
rule, and a deck rendered half in one palette and half in another is a defect that no single script
can see. `plotting/filters.py` is the other shared module, for the same class of reason -- it
carries a claim, not a preference.

**The localization vocabulary is already npg and is NOT redefined here.** `src/localization.py`
pins `LOC_CLASS_COLOR` to #E64B35 / #00A087 / #3C5488 / #F39B7F / #8491B4 / #7E6148, which ARE
ggsci npg colours. Figures keep importing it from `src/` so the deck and the localization axis's own
plots cannot drift apart; redefining those six here would create exactly that drift.

Semantic names, not indices, because `NPG[5]` at a call site says nothing about what it is for and
silently changes meaning if the palette is ever reordered.

Run with the `gradi` env. Imported, never executed.
"""

from __future__ import annotations

import stylia
from matplotlib.colors import LinearSegmentedColormap

#: The raw palette, 10 colours, in stylia's npg order.
NPG = stylia.CategoricalPalette("npg").colors

# -- semantic roles ------------------------------------------------------------------------------
PRIMARY = NPG[5]      # steel blue  -- the main series, the thing the panel is about
SECONDARY = NPG[0]    # red         -- the contrasting series
TERTIARY = NPG[4]     # turquoise   -- a third series
QUATERNARY = NPG[2]   # amber       -- a fourth, rarely needed
ACCENT = NPG[8]       # magenta     -- a highlight ON TOP of a muted field, never a series colour
MUTED = NPG[9]        # grey        -- "everything else", denominators, backgrounds-with-data

#: Text, axis rules and reference lines. Not from the palette: it must stay legible in print.
INK = (0.17254901960784313, 0.24313725490196078, 0.3137254901960784)

#: A near-white field for points that are present but not the subject. NOT a data colour, which is
#: why it is not drawn from the categorical palette.
BACKDROP = "#E4E4E0"

#: A sequential ramp for continuous fields (heatmaps, percentiles), built FROM the palette rather
#: than taken from matplotlib's own set -- `BuPu`, `viridis` and friends are not npg, and a deck
#: that is npg everywhere except its heatmaps is the drift this module exists to prevent.
SEQUENTIAL = LinearSegmentedColormap.from_list("npg_sequential", ["#FFFFFF", PRIMARY])

#: Per-species, where a figure compares the three bacteria.
SPECIES_COLOR = {"kpneumoniae": PRIMARY, "ecoli": SECONDARY, "saureus": TERTIARY}

#: The three COG top-level groups plus the two "we do not know" cases. npg, and deliberately
#: desaturated for the last two so "poorly characterised" never reads as a finding.
COG_GROUP_COLOR = {
    "information": NPG[4],
    "cellular processes": NPG[5],
    "metabolism": NPG[0],
    "poorly characterised": NPG[9],
    "not classified": BACKDROP,
}

#: Studiedness evidence tiers, strongest first. A diverging ramp: the two tiers that both read 0 in
#: the shipped table are the warm pair, so the ambiguity is visible as colour.
TIER_COLOR = {
    "swissprot_direct": NPG[4],
    "swissprot_close": NPG[3],
    "swissprot_homolog": NPG[9],
    "below_floor": NPG[1],
    "no_hit": NPG[0],
}
