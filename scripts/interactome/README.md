# interactome — no code yet

**What does this protein interact with, and how central is it?** Feeds target prioritisation two
ways: a hub is more likely essential, and a degrader's effect propagates through partners.

Candidate sources, none yet evaluated:
  - **STRING** v12 — per-organism protein links with separate evidence channels. Keep the channels
    apart: `experiments` and `database` are evidence, `textmining` and `neighborhood` are inference,
    and merging them into one `combined_score` is what makes STRING look better than it is.
  - IntAct / BioGRID for curated binary interactions
  - operon / genomic-neighbourhood structure (already partly available from `orthology/`)

**Check accession coverage before building anything.** STRING keys on its own protein ids per
genome; whether HS11286 is present at all is an open question, and the project has hit the
"external database does not know the anchor" wall four times (CLAUDE.md § *Identifier convention*).
If it is absent, map **by sequence**, not by accession.

Deliverable: a **complete matrix** — degree / betweenness / channel-specific counts, one row per
protein, canonical row order. A protein with no recorded partner gets a measured 0, not an empty row.
