# pockets — no v2 code yet

Structure-based druggability: AlphaFold models, fpocket / P2Rank pocket detection, and the
ray-traced cartoons that went with them.

**Nothing here has been ported to v2.** The v1 implementation is frozen under `legacy/scripts/`
(`06e_pockets.py`, `06m_pocket_*`, `06n_structure_snapshots.py`, `_06n_pymol_render.py`) and runs
from inside `legacy/` — see CLAUDE.md § *Legacy*.

Envs it will need, both already in `install.sh`: `gradi-pockets` (fpocket + `openjdk=17` for
P2Rank, osx-64 under Rosetta) and `gradi-pymol`.

When it is built, the deliverable is a **complete matrix** — one row per protein of the reference
proteome, in the canonical row order (CLAUDE.md § *Every axis ends in COMPLETE matrices*).
