# The v2 website (`webapp/`): "Gram-negative target browser for protein degradation"

Record of the 2026-11-06 build: a new browser for **E. coli** and **K. pneumoniae**, built from the v2 per-axis tables. The
v1 site was not touched while this was built (owner: "do not supersede the previous site"). On 2026-11-06 the owner reversed that ("I want the new site to be the default"): the new
browser is the root of Pages and the v1 browser is served at `/v1/` (see Revision 17).

## Decisions (owner)

1. **Reuse the existing app** (shell, table, gene card, 3D viewer, map) with a `webapp/DESIGN.md` as guardrail.
2. **Keep the weight sliders.** This is the owner's single exception to "no composite, ever". Mitigation: the sliders combine the
   `*_consensus` columns in the browser only, are OFF by default, the result is labelled "Weighted score (yours)", and it is never
   exported (CSV) or written to the data files. `Neglectedness` = 1 - `studiedness_consensus`.
3. **No shortlist focus.** The site is the final per-axis tables. Neither the E. coli ASKA 200 nor a Kp shortlist is shown; the
   quick selections are stated predicates instead.
4. **E. coli and K. pneumoniae only.**

## What was built

- `scripts/webapp/export.py` (env `gradi`): reads only the `src/` loaders (proteomes, function, localization, degradability,
  essentiality, ligands, pockets, orthology, studiedness, projection, DEG), keeps canonical order, writes `webapp/data/{ec,kp}.json`
  (columnar, 55 columns, E. coli 4,403 rows 1.51 MB, K. pneumoniae 5,728 rows 1.80 MB; v1 files were 4.0 and 4.4 MB). Adds `cog` (first COG
  letter, `none` = no hit) and renames `uniprot_ac`/`gene_name`/`protein_name` to the names the browser code reads. Verifies the written
  file against the frame (shape, column order, per-column non-null counts, 2,000 random cells). An independent re-read of the
  source TSVs agreed on every consensus/probability/count column checked, same NaNs, same row order.
- `webapp/config.js` rewritten for v2 (columns with definitions, eight expandable column groups, gene-card panels, filters, rank bands, quick selections,
  Methods text, group notes). `webapp/app.js` patched (about 40 edits, listed in git history once committed): the v1 composite
  components, tiers, transfer toggle, mock hatching, pocket-residue overlay and gene-symbol family search are gone.
- `webapp/DESIGN.md` (Google DESIGN.md, `version: alpha`) and `tools/check_design_md.py`. The offline check found one real defect on its
  first run: the caution note text (#9A6B12 on a pale amber) had contrast 4.25 against the AA minimum 4.5; fixed to #7A520A.
  `npx @google/design.md lint` was denied in this session (it needs a download), so the official linter has NOT been run.

## What the page says and refuses to say

- Every axis: rank (within-species percentile) beside evidence level (1-3, provenance not outcome). Ranks are never compared across species.
- Degradability on both species is extrapolated from *S. aureus* screens and falls with protein length; the view note says so.
- E. coli essentiality consensus is partly circular with the measured DEG calls. K. pneumoniae HS11286 has **no measured essentiality**:
  `deg_*` columns are exported empty and the table hides columns that are empty for the active species; the gene card shows a note
  instead of a row of zeros.
- A blank is "not annotated / not measured", never 0. The ligandability consensus is pinned to 0 for the ~97% with no evidence, so the
  gene card suppresses the percentile for tie blocks above 25%.

## Revision 2 (same day, owner feedback)

"Not clear that we need so many tabs": the nine view tabs are replaced by **one table whose column groups expand and collapse** (click the coloured
group header; "Expand all" toggles everything; state persists and is in the share link). The site is renamed **"Gram-negative target browser for
protein degradation"**. The weighted-score column now only appears while a slider is on. Pockets became their own group (it was part of
Ligandability) and each expanded group adds its caveat above the table. The browser checks below predate this revision except the expand/collapse
test: collapsed 19 columns, expanding Degradability adds 5, Expand all gives 48, K. pneumoniae hides the empty measured-essentiality column.

**Revision 3:** the gene card was rebuilt to match the single table: header, 3D model, then one section per column group (same order, colours and
open/closed state as the table; clicking a section also opens it in the table). The ring, contribution bars, per-axis bespoke panels and the
"Protein universe" mini-map were removed. Checked on E. coli `gyrB`: opening a card section adds it to the table; the model loads; no console errors.

**Revision 4 (sidebar):** rebuilt as search, a live "n of N proteins" count with Clear all, quick selections, then filter accordions with
active-count badges (Rank bands open; Localization, Evidence level by axis as L1-L3 chips, Human ortholog, Function as COG letters with names on hover);
the weighted score stays collapsed at the bottom. Previously the sidebar was a long unfolded list (24 COG name chips, unlabelled evidence chips).
Checked in Chrome: the Degrader preset gives 63 of 4,403, badges read 2/1/1, Clear all appears, no console errors.

**Revision 5:** the three quick selections (Degrader, Druggable, Neglected) were removed on the owner's instruction, together with their code. Filters are now
only the rank bands, localization, evidence level, human ortholog and COG letters, plus search. Revision 5 also made the table groups
one-open-at-a-time ("Expand all" removed), put Function before Orthology, and reworked the table (horizontal headers, sticky Target column and
header rows) and the responsive layout (off-canvas filters below 900 px).
Quick-selection counts quoted below (63 / 53 / 193) are historical.

**Responsiveness (revision 5), measured and checked.** Speed: on K. pneumoniae (5,728 rows) a filter or sort recompute takes about 17-25 ms and rendering a
100-row page about 13 ms (Chrome, this Mac), so interaction is not the bottleneck; the work went into layout. Table: explicit `<colgroup>` widths (the
colspan band row broke fixed layout), horizontal two-line headers, sticky Target column and two sticky header rows, widths chosen per breakpoint in
`colWidths()` and re-applied on resize. Breakpoints: above 1200 px full sidebar; 900-1200 px narrower sidebar; below 900 px the filters become an off-canvas
drawer opened by "Filters" with a "Show N proteins" button; below 700 px the toolbar buttons become icons and the table keeps only a slim # and Target column.
Checked at 1500, 820 and 500 px wide. A class-name clash (`.acc` for both the accordions and the accession text) put a border round every accession on phones; renamed to `.facc`.
Not checked: real touch devices, landscape phones, very wide screens, Safari/Firefox.

**Revision 6: UX audit and restyle (owner: "very extensive audit", "cannot look claudesque", "the header looks claudesque").** Findings and what changed:
header with a mono eyebrow, italic purple title and an amber preview chip -> a slim plain header bar; rainbow heat pills, coloured L1/L2/L3 badges, red crosses, emoji
links, soft shadows and 14 px radii -> ink on white, a 3 px bar in the axis colour under each value, plain 1/2/3 for evidence, 6 px radius, no shadows; a paragraph
banner above the table and notes in every filter and card section -> removed (one-line caveat in the group header tooltip, full text in Methods); the unexplained
three-dot evidence column -> removed; 15 columns collapsed -> 10 (one per group, plus Length and measured essentiality); repeated "Deg. rank/Ess. rank" labels ->
"Rank" under the group name; segmented "any" badges on every filter -> click again to clear; "Shortlist" (confusable with the removed shortlists) -> "Starred".
The stylesheet was rewritten from scratch (about 330 lines, was 600+ with v1 leftovers) and the Methods text cut by half. Protein length is now its own fixed
column and a filter (bands up to 100, 101-150, 151-300, over 300 aa). A scripted pass over both species (every column sorted with every group opened, every filter chip,
search, 30 random gene cards, every map colouring) found no errors. Open: header tooltips are still one to two sentences each; touch devices untested.

**Revision 7 (gene card, sidebar, responsiveness, colour).** Owner: continue the audit on card and sidebar, make the site extra responsive; then "too sad, I liked the
previous palette". Card: up/down buttons and arrow keys (also j/k) step through the filtered list without closing it, with a position counter ("3 of 4,403") that
crosses page boundaries; swipe right closes on touch; focus returns to where it was on close. Sidebar: the search box and count stay pinned while the filters scroll,
evidence filters are one inline row per axis, which sections are open is remembered, touch targets are 36-38 px below 900 px. Table: bottom pager, `touch-action`
and press states, sort cached so a filter change no longer re-sorts, rows built as one HTML string with delegated clicks, band counts memoised, dead code for
the old family search and presets removed. Timings on K. pneumoniae (5,728 rows): filter or sort click 15-19 ms (was 17-25), card open 0.4 ms, next protein 2.6 ms.
Colour: the earlier palette is back (periwinkle interface, one hue per axis) but used as tinted group headers, soft data bars behind values, coloured 3/2/1 evidence
and a thin axis-colour strip above the header, not as pills.

**Revision 8 (map defaults, nothing moves on click).** Map: "Show >=" now defaults per colouring from the data (ranks: their 80th percentile, so about the top fifth;
ligandability, zero for 97%: anything above zero; cytoplasmic fraction: the maximum), is remembered per colouring, and double-click resets it; the controls moved
onto the map so the toolbar no longer reflows; "Weighted score" is only offered as a colouring while it exists. Click-time movement, found with a scripted layout
diff (positions of the toolbar, header, sidebar and table before and after 25 actions): a 1 px press transform on every button (mine, removed); the table jumping
to the top when a column group was opened or a row clicked (now only on filter, sort and page changes); scrollIntoView on a plain row click (now only when stepping
with the card); the sidebar scroll snapping back when a filter panel was rebuilt (bands update in place, other rebuilds restore the scroll); scrollbars appearing and
shifting content (stable gutters); header rows changing height when a label wraps (fixed 50 px); sort arrows and counters changing width (reserved space); the pager
vanishing on the map (kept, invisible); the Starred button widening with its count. After the fixes every action leaves the toolbar, header and sidebar within 0.6 px and
the table scroll position unchanged. Secondary text colour darkened for contrast (#766f84).

**Revision 9 (critical aesthetic pass).** Findings and changes: the full-cell data bars made numbers straddle the fill edge and the table read as a wall of slabs ->
a slim 8 px bar in the axis colour followed by the number (colourful, no overlap); group headers were flush tinted blocks whose labels did not line up with the numbers
below -> tab-style headers (3 px coloured top rule, light tint, a gap between groups, ink label); "no" repeated thousands of times -> a quiet dash (kept as yes/no only
for human ortholog); a star on every row -> shown on hover or when starred (always on touch); rows 50 -> 46 px; the sort arrow made short headers wrap and shifted them
-> label and arrow are one unit; Length given a wider column; the open gene-card section repeated the rank and evidence already shown in its header line -> removed.
Kept on purpose: the periwinkle interface colour and the axis hues, per the owner's preference.

**Revision 10 (the card no longer moves the table).** Owner: the table moves when using the gene card. Cause: the card shared its open section with the table (opening a card section
opened that group in the table, re-rendering it), and the step buttons scrolled the table to the selected row and switched its page. The card now has its own open section, and
stepping leaves the table's scroll, page and columns alone. Checked by script with the table scrolled 500 px down and 200 px across: opening and switching card sections,
120 steps with the buttons, an arrow key and closing the card all left scroll position, page, open group and row positions unchanged.

**Revision 11 (motion).** Owner: transitions should feel dynamic and modern. Added, all opacity or transform so layout never moves: a thumb that slides behind the selected species and view
buttons; rows fading in with a short stagger and value bars growing when the list changes (filter, sort, page; not on opening a column group); an easing curve of
cubic-bezier(.22,1,.36,1) on the gene card, which slides in while its sections rise in turn (only when the card first opens, not while stepping); open accordions fade in; chevrons turn;
the star pops; dialogs scale in; the map and compare views fade in; a shimmering skeleton while a species first loads (a cached species no longer blanks the table);
hover tints fade. Disabled under prefers-reduced-motion. After the change, 8 click scenarios (band, sort, group, page, card open and close) still leave the toolbar, headers
and sidebar in place.

**Revision 12 (sidebar rebuilt around sliders).** Owner: the left bar needs drastic focus; ranges should be sliders with both ends, less clicking, search limited to names and accessions.
The three-band buttons per score and the length bands are replaced by dual-handle range sliders (drag either end, click the track, arrow keys with shift for big steps,
double-click the name to clear) with a square-root histogram of the species' values behind each so the shape of the data is visible (ligandability's spike at zero no longer hides
the rest). Length is on a log scale. Evidence level is also sliders (1 to 3), collapsed. Localization chips, a one-line human-ortholog switch and COG letters (collapsed)
remain categorical. Dragging recomputes the table every frame without replaying the fade or scrolling it; releasing commits. The sidebar went from about 25 clickable
bands and chips to 6 sliders and 1 chip row on screen. Search no longer reads COG or GO text. Verified by script: drag, track click, keyboard, reset, log length,
evidence sliders, URL round trip, clear all, and search scope (accession exact, gene prefix, name, no COG text).

**Revision 13 (sidebar audit, faster table).** Sidebar: each slider's histogram is now cross-filtered: the dark part of every bin is how many proteins remain once all the OTHER filters apply, the light part
the full distribution, so a slider shows what the rest of the filters have already removed (cost about 3 ms per update). The "Clear filters" link no longer changes the
row height by a tenth of a pixel. Table: measured on K. pneumoniae (5,728 rows, 100 per page), a filter or sort went from 24 ms to 10.6 ms and the page render from 13.4 ms to
3.6 ms, by (a) dropping the extra element inside every value cell (the bar is now a pseudo-element, about 17% fewer nodes), (b) rendering the first 28 rows at once and the rest on
the next frame, and (c) `content-visibility: auto` on rows so rows below the fold are not laid out until they scroll in (sticky columns checked, still correct). Dragging a
slider costs about 9 ms per frame. A scripted pass over every column of both species, a burst of ten re-renders and a slider drag found no errors and no layout movement.

**Revision 14 (localization circle, ligandability control).** Localization is now a drawing of the cell: concentric rings (centre cytoplasm, then inner membrane, periplasm, outer membrane,
extracellular) in the axis colours; click a ring or its legend row to select (several allowed); ring opacity and the counts follow the other filters. Ligandability: a plain 0-1 slider was
not the best option, because the score is pinned to 0 for the ~97% nobody screened, so the histogram was one spike and any drag mixed "no evidence" with "weak evidence". It is now a
three-way choice (any / has / none) plus a slider over only the proteins that have evidence (starting at 0.01); touching the slider implies "has". E. coli: 160 proteins have
ligand evidence, 4,243 none. Checked by script: ring and row clicks, cross-filtered counts (E. coli cytoplasm 2,601 -> 294 under an essentiality filter), any / has / none, slider drag,
URL round trip and clear. Other options considered for ligandability: a log slider on the number of potent ligands (needs a count column the consensus does not use), or a
histogram only; the three-way plus slider was the simplest that matches how the number is built.

**Revision 15 (circle aesthetics).** The first localization circle was a thick rainbow bullseye (five axis colours as pastel rings, plus a coloured dot per legend row) that fought with the rest of the
panel. Redrawn as a quiet diagram: one lavender hue, thin white membrane lines between the compartments, shade showing how many proteins remain, the selection in the interface colour, and a plain
text legend with counts that fits the column. Smaller (80 px).

**Revision 16 (sidebar: colour, COG names, nothing hidden).** Localization uses its own axis colour (green) instead of lavender, with no counts (compartments with nothing left just fade);
Function is a grid of named tiles in the four COG groups (Information storage, Cellular processes, Metabolism, Poorly characterized) with plain names instead of letters (tooltip: letter and full COG name),
no counts, tiles with nothing left dimmed; and nothing in the sidebar is collapsed any more: Function, Evidence level and Weighted score are always open. Verified: no `details` or hidden
elements remain in the sidebar; tile selection, badge, filter and reset work; no label is truncated.

**Revision 17 (release and deploy).** v1.0.0 was tagged on main (`dd6e76f`); `legacy/app/` was restored byte-for-byte from it; `webapp/data/*.json` became committed files (gitignore exception);
`pages.yml` now assembles `_site/` = `webapp/` at the root plus the v1 browser at `/v1/` and fails the deploy if a payload is empty; the header gained a "Previous version" link (removed 2026-11-06 at the owner's request; `/v1/` is still served, just not linked). A local dry run of the
assembly served both sites with no console errors (new: 100 rows, v1: 100 rows).

**Licence note.** TabPFN-3.5 weights are non-commercial (CLAUDE.md flags this must be settled before outputs ship). Owner's call, 2026-11-06: proceed, treat the use as non-commercial
(Ersilia is a non-profit; the browser is free), and say so in the Methods dialog and in the pull request.

## Verified in a browser (Chrome, local server, 2026-11-06)

All nine views render for both species with no console errors; the quick selections give 63 / 53 / 193 E. coli proteins, equal to an independent
pandas recomputation of the stated predicates; the weighted score equals a hand computation (0.6 degradability + 0.4 essentiality = 0.9864 for the top
row); the CSV export carries no weighted column; the URL hash round trip restores the state; each view opens sorted by its own headline column
(K. pneumoniae: KPC-2 tops ligandability, `rplE` essentiality, `crp` studiedness); the AlphaFold model of `clpP` loads live in the 3D viewer. Phone width
was only checked as far as the sidebar stacking.

## Gaps and not done

- **Not committed, not deployed.** `webapp/`, `scripts/webapp/`, `tools/check_design_md.py`, `docs/webapp.md` and the CLAUDE.md section are
  untracked or uncommitted. Publishing needs a decision: a new Pages path, a staging branch, or replacing `app/`.
- **3D viewer** fetches AlphaFold DB models live (`alphafold.ebi.ac.uk`); it was not exercised in the browser check (needs the network).
  The 32 E. coli ESMFold-modelled proteins and 1 K. pneumoniae protein have no AFDB model and show a message. v2 has pocket
  **scores** only, so the v1 pocket-residue overlay is gone.
- **K. pneumoniae annotation gaps** show as blanks (PubTator counts 62.5% filled, many accession-only names).
- **Not checked:** Safari/Firefox, keyboard navigation after the patches, the hash deep link round trip across species, print.
- The v1 columns with no v2 equivalent were dropped (family clusters, InterPro/PANTHER family names, popularity tiers, selectivity
  classes, Clp/degron features, FBA). The site no longer has a protein-family search.
- Top of the default degradability sort is dominated by short, uncharacterised proteins: a property of the model (length confound), stated
  in the view note, not a bug.

## Revision 18 (2026-11-06): state hardening, no "Previous version" link

A live link carried a stale hash (v1 weight keys, `map: "localization"`) and the site opened filtered and coloured by things the
sidebar does not show. `load()` and `applyHash()` now pass saved and shared state through `cleanWeights`, `cleanRanges`,
`cleanCats`, `cleanBools` and `cleanMapColor`, which keep only keys that exist in `config.js`. The "Previous version" header link
was removed; `/v1/` is still deployed by `pages.yml`.

## Revision 19 (2026-11-06): social preview

`index.html` now carries Open Graph and Twitter card tags (absolute URLs, `summary_large_image`) and a `description`. The card is `webapp/social-preview.png` (1200x630), rendered by `scripts/webapp/social_preview.py` (env `gradi`), which reads the protein counts from `webapp/data/*.json` so the card cannot disagree with the site. Text stays inside a 90 px safe margin and was checked at 320 px wide, the size chat apps show on a phone. Re-run the script if the payload size changes. Link scrapers cache the image, so a changed card may take time to refresh after a deploy.
