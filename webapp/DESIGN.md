---
version: alpha
name: Gram-negative target browser for protein degradation
description: A calm, table-first browser for per-protein evidence on E. coli and K. pneumoniae. Colour carries the axis, never a verdict.
colors:
  primary: "#6C5CE7"
  on-primary: "#FFFFFF"
  ink: "#2C3E50"
  muted: "#6B6675"
  faint: "#9A93A6"
  background: "#FAFAFC"
  surface: "#FFFFFF"
  surface-2: "#F4F4F8"
  border: "#E6E6EE"
  essentiality: "#E63946"
  degradability: "#FCBF49"
  ligandability: "#457B9D"
  localization: "#6BBF59"
  orthology: "#2EC4B6"
  studiedness: "#B05CC8"
  function: "#6C5CE7"
  ecoli: "#F4845F"
  kpneumoniae: "#457B9D"
  caution-ink: "#7A520A"
  caution-surface: "#FEF3DC"
typography:
  wordmark:
    fontFamily: Inter
    fontSize: 22px
    fontWeight: 700
    lineHeight: 1.1
    letterSpacing: -0.01em
  body:
    fontFamily: Inter
    fontSize: 13px
    fontWeight: 400
    lineHeight: 1.5
  label:
    fontFamily: ui-monospace
    fontSize: 10px
    fontWeight: 600
    letterSpacing: 0.06em
  number:
    fontFamily: ui-monospace
    fontSize: 12px
    fontWeight: 500
    lineHeight: 1.3
rounded:
  sm: 7px
  md: 10px
  lg: 14px
  full: 999px
spacing:
  xs: 4px
  sm: 8px
  md: 16px
  lg: 24px
components:
  button-primary:
    backgroundColor: "{colors.primary}"
    textColor: "{colors.on-primary}"
    typography: "{typography.body}"
    rounded: "{rounded.sm}"
    padding: 8px
  button-quiet:
    backgroundColor: "{colors.surface}"
    textColor: "{colors.ink}"
    typography: "{typography.body}"
    rounded: "{rounded.sm}"
    padding: 8px
  panel:
    backgroundColor: "{colors.surface}"
    textColor: "{colors.ink}"
    rounded: "{rounded.lg}"
    padding: 16px
  note-caution:
    backgroundColor: "{colors.caution-surface}"
    textColor: "{colors.caution-ink}"
    typography: "{typography.body}"
    rounded: "{rounded.md}"
    padding: 8px
  evidence-chip:
    backgroundColor: "{colors.surface-2}"
    textColor: "{colors.muted}"
    typography: "{typography.label}"
    rounded: "{rounded.full}"
    padding: 4px
---

# Gram-negative target browser for protein degradation

## Overview

A plain, table-first browser: one row per protein, one species at a time. It should read like a well-made data table, not a dashboard: no gradients, no emoji, no pills, no explanatory banners. Colour is used for one job, telling the axes apart. Caveats live in hover tooltips and the Methods dialog, never in the page body.

## Colors

The project's palette: periwinkle for the interface (active tab, selected filter, links, badges), and one hue per axis: amber degradability, crimson essentiality, cobalt ligandability and pockets, lime localization, periwinkle function, turquoise orthology, orchid studiedness. Axis hues appear in four places only: the tinted group header, the soft data bar behind each value, the dot beside a filter, and the section marker in the gene card. A thin strip of the axis colours runs along the top of the header. Text stays ink or grey.

## Typography

Inter throughout, sentence case, no letter-spaced capitals, no monospace except accessions. Numbers are tabular. Headings are 15px at most.

## Layout

A slim header (title, species switch, Methods, Ersilia, Code) over a left filter column and one table. The filter column is built for dragging, not clicking: **range sliders with both ends** for every score and for length (log scale), each with a small histogram of the species' values behind it; localization drawn as the cell itself (a quiet diagram in the Localization green, click a ring or its row, shade follows the other filters); ligandability as any / has / none plus a slider over the proteins that have evidence; a three-way human ortholog switch; Function as a grid of named COG tiles in their four groups, and Evidence level as 1-3 range sliders; nothing in the sidebar is collapsed or hidden. Search matches gene names, UniProt accessions and protein names only. The table has a sticky header and sticky Protein column. Each axis is a group of columns: collapsed it shows one column (its rank, or its class or flag), open it shows all of them, one group at a time. Length is its own always-visible column. Below 900px the filters become a drawer; below 700px the table keeps only a thin number column and the Protein column.

## Elevation & Depth

None. One-pixel borders only.

## Motion

Short, springy-eased (cubic-bezier .22,1,.36,1) and always opacity or transform, never layout. Segmented controls slide a thumb; table rows fade in and value bars grow when the list changes (not when a column group opens); the gene card slides in and its sections rise in turn; dialogs and tooltips fade; stars pop; loading shows a shimmering skeleton. Everything is switched off under prefers-reduced-motion.

## Shapes

6px radius on controls; no pills.

## Components

- **Value cell**: the number, right-aligned, over a soft bar in the axis colour that fills in proportion to the value. A missing value is a dash.
- **Evidence**: a plain 1, 2 or 3, green when 3, amber when 2, grey when 1. Provenance, not outcome.
- **Group header**: name with a chevron on a tint of the axis colour, underlined in it; hover shows the group's one-line caveat.
- **Gene card**: the AlphaFold model, then the same groups as the table, with its own open section (using the card never changes the table behind it). Up/down arrows (or keys) step through the filtered list without closing it; the table does not scroll or change page.
- **Weighted score**: optional, collapsed, in the sidebar; adds one column while on; never exported.

## Do's and Don'ts

- Do keep rank beside evidence wherever both are shown, and say "not measured" rather than showing zero.
- Do keep ranks within one species.
- Don't add banners, pills, emoji, shadows or gradients (the header strip is the one exception).
- Don't add a combined score other than the optional, labelled, browser-only weighted one.
- Don't show measured essentiality for K. pneumoniae: there is none.
- Don't use axis colours for decoration, or put long explanations in the interface.
- Do keep every control reachable with a thumb: 36px minimum on touch screens.
