# Target-Prioritization Browser (prototype)

An OpenTargets-style, dependency-free web app to browse the GraDi per-protein
prioritization results for **K. pneumoniae** (HS11286) and **E. coli** (K-12),
with an **adjustable-weight composite score** and rich filtering.

Everything runs client-side (vanilla HTML/CSS/JS); there is no backend. The only
build step is a Python export that snapshots the pipeline CSVs into JSON.

## Live site

Deployed to GitHub Pages by `.github/workflows/pages.yml` on every push to `main`:
**https://ersilia-os.github.io/gradi-target-prioritization/**

The site data (`app/data/{kp,ec}.json`) is committed to the repo (it can't be regenerated in CI —
the pipeline inputs live in eosvc/S3, not git). To refresh the published data: re-run
`scripts/08a_webapp_export.py`, commit the updated `app/data/*.json`, and push to `main` → auto-redeploys.

One-time setup: repo **Settings → Pages → Source: GitHub Actions**.

## 1. Regenerate the data

Run after the pipeline outputs change:

```bash
conda run -n gradi python scripts/08a_webapp_export.py
```

This joins `output/results/{kpneumoniae,ecoli}/*_essentiality.csv`,
`*_ligandability.csv`, `*_esmc600m_projection.csv`, `*_alphafold_structure.csv`,
and `*_pdb_coverage.csv` on `uniprot_accession`, derives the 0–1 composite
components, and writes `app/data/kp.json` and `app/data/ec.json` (gitignored).

## 2. Serve it

```bash
bash app/serve.sh          # PORT=8080 by default
```

Open from your laptop (same Tailscale network as this Mac mini):

- `http://miquel-macmini:8080`  (MagicDNS), or
- `http://100.105.44.112:8080`  (Tailscale IP)

Locally on the Mac mini: `http://localhost:8080`.

## 3. Use it

- **Organism toggle** (top-left) switches between K. pneumoniae and E. coli.
- **Views** (top tabs) — pick the column set:
  - **Overview** — the high-level scorecard (both axes + tiers + selectivity).
  - **Essentiality / Ligandability / Degradability / Structure / Cross-species / Novelty /
    Annotation** — per-axis views that surface that axis's detailed evidence columns.
  - **◵ Map** — the ESM-C protein-universe projection (every protein at its
    `tsne_x/tsne_y`), coloured by any component via the **Colour by** dropdown. Current
    filters carve the universe: passing proteins are coloured, the rest fade to grey
    context. Hover for a tooltip; click a point to open its detail drawer.
- **Composite weights** — toggle each component on/off and set its 0–100 weight.
  The composite is a per-protein weighted mean over *enabled components that have a
  value*, renormalized per protein (missing components are dropped, not zero-filled).
  Essentiality + Ligandability are on by default; **Degradability** and **Novelty** are wired
  (off by default). *Expression/localization is not yet implemented.* An **Evidence** column
  shows how well each target is supported (measured vs inferred/predicted); a **Methods ⓘ**
  panel documents every score.
- **Orthology transfer** switch (top-left) recomputes essentiality without the E. coli→Kp
  cross-species transfer track.
- **Gene card** (click any row) — a per-protein profile: composite ring + rank + evidence,
  per-axis panels, an interactive **AlphaFold 3D structure** (coloured by pLDDT, top predicted
  pocket highlighted), a mini locator map, and annotation.
- **Provisional data** (E. coli degradability is a mock) is shown with hatched cells + a
  "provisional" flag, never as if it were real.
- **Presets** — "★ Prime targets" = essential ∧ tractable ∧ broad-selective;
  "◐ Neglected & druggable" = tractable ∧ dark (under-studied).
- **Filters** — text search (by **gene name or UniProt accession**), min-threshold
  sliders, categorical chips, tri-state booleans.
- **Table** — click any header to sort; **hover a header for a plain-language
  explanation** of that column; the **Columns** menu fine-tunes which columns show on top
  of the chosen view; per-axis cells are heat-shaded; **Export CSV** dumps the current
  filtered/sorted view.
- **Click a row** for the detail drawer (all evidence grouped by axis, composite
  contributions, and external links to UniProt / AlphaFold / InterPro / PDB).

State (view, weights, filters, columns, map colour, theme) persists in `localStorage`.

## Extending

- **Add a composite component or table column:** edit `config.js` (single source of
  truth) and add the underlying column to `scripts/08a_webapp_export.py`. Anything
  referencing a key absent from the data is auto-skipped / greyed out. Each column carries
  a `desc` string — that is the hover explanation.
- **Add a per-axis view:** add an entry to `TABLE_VIEWS` in `config.js` with a `cols` list;
  it appears as a new top tab automatically.
- **Add a map colour option:** add a `{key,label}` to `MAP_COLORS` in `config.js`.
- Files: `index.html` (shell) · `styles.css` (design system) · `config.js` (columns +
  components) · `app.js` (logic) · `serve.sh` (server) · `data/` (generated, gitignored).
