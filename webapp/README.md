# Gram-negative target browser for protein degradation

A dependency-free web app (vanilla HTML/CSS/JS, no backend) to browse the v2 per-axis tables for
**E. coli** K-12 and **K. pneumoniae** HS11286: degradability, essentiality, ligandability, pockets,
localization, orthology, studiedness and function. One row per protein; every axis shows its
within-species **rank** beside its **evidence level** (1-3).

**This does not replace the public site.** The deployed app (`app/` on `main`, `legacy/app/` on this branch)
and `.github/workflows/pages.yml` are untouched. Nothing here deploys until the owner decides how to publish it.

## Regenerate the data

```bash
~/miniconda3/envs/gradi/bin/python scripts/webapp/export.py        # writes webapp/data/{ec,kp}.json
```

It reads only the `src/` loaders, keeps the canonical row order, and verifies the written files
(shape, per-column non-null counts, 2,000 random cells). `webapp/data/*.json` is about 1.5-1.8 MB each.

## Serve

```bash
bash webapp/serve.sh        # PORT=8080 by default
```

## What is on the page

- **One table.** Each axis is a column group (degradability, essentiality, ligand precedence, pockets, localization, function, orthology, studiedness), with length as the last column. Collapsed, a group shows its score (and, for some, a headline flag); click its coloured header to expand all its columns. Both species show exactly the same columns; a column with no data for a species is faded and blank. Compare (starred proteins side by side; ✕ removes one, its name opens its card) and Map (t-SNE of ESM-C) are the only other modes. The table shows the recommended protein name; the full UniProt name is on hover and in the card.
- **Filters**: range sliders with both ends for each score and for length (log scale, histogram behind each), localization drawn as the cell, a "No human ortholog" on/off switch, and Function (COG) as compact chips, all on by default. Search covers gene names, UniProt accessions and protein names.
- **Weighted score (yours)**: a radar right after Localization. Drag each axis's dot along its spoke (centre = off);
  the default is Degradability, Essentiality, Ligand precedence and Pockets equally, with Novelty (1 − studiedness) off.
  The table opens sorted by it. It is computed in the browser only. It is the
  owner's single exception to the project's "no composite" rule, is labelled as the reader's own, and is never exported
  or written to the data files.
- **Gene card**: the AlphaFold model (live from the AlphaFold DB, coloured by pLDDT) and the same column groups as the table, each section opening on its own; the table behind never moves.

## Caveats the site states

Degradability is extrapolated from *S. aureus* screens (ranking hypotheses) and falls with protein length; E. coli
essentiality consensus is partly circular with its measured calls; K. pneumoniae HS11286 has no measured essentiality
(those columns are hidden there); pocket rank is largely protein length; a blank is "not annotated", never 0.

## Guardrail

`DESIGN.md` is the design contract (Google DESIGN.md format). Check it offline with
`~/miniconda3/envs/gradi/bin/python tools/check_design_md.py`; the official linter is `npx @google/design.md lint webapp/DESIGN.md`.

Method and decisions: `docs/webapp.md`.
