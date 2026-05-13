# Pipeline overview

Diagrams describing the GraDi target-prioritization pipeline. The pipeline is
documented as five sections, each in its own file under `docs/`:

1. [Task-agnostic per-protein annotation](./01_task_agnostic.md) — proteome →
   per-protein evidence consumed by all downstream axes.
2. [Ligandability assessment](./02_ligandability.md) — can a small-molecule
   recruiter engage the target?
3. [Degradability assessment](./03_degradability.md) — is the target a
   plausible Clp-protease substrate?
4. [Essentiality / vulnerability assessment](./04_essentiality.md) — is the
   target required for fitness or survival?
5. [Expression and localization](./05_expression_and_localization.md) — is
   the target present in meaningful amounts and physically reachable by the
   cytoplasmic Clp machinery?

A final-ranking section will be added once the corresponding workstream comes
online.

## Diagram style

Every diagram in this directory uses the canonical Mermaid style defined in
[`mermaid_style.md`](./mermaid_style.md), anchored on the Ersilia brand
palette. Seven node classes:

- `:::source` (purple) — anchor proteome / starting entity.
- `:::dataset` (yellow) — curated reference data file (e.g. Flynn 2003, Nagar 2021).
- `:::method` (blue) — script / compute step.
- `:::result` (mint, thicker border) — output, score, or sink.
- `:::tagnostic` (pink) — pre-computed input from the task-agnostic layer (Section 1).
- `:::stub` (orange, dashed border) — parser exists in the codebase but the data file is not yet staged in `data/raw/`.
- `:::planned` (gray, dashed border) — not yet implemented.
