# pockets — structural ligandability

Can a small molecule bind this protein's fold? Two answers that share no input, never merged:
**predicted** pockets on AlphaFold models (fpocket, P2Rank) and **measured** drug-like holo
structures in the PDB (BioLiP).

Run in order, all with the `gradi` env:

| script | does | writes |
|---|---|---|
| `structures.py` | AlphaFold DB v6 models, validated against the proteome sequence | `data/source/alphafold/<sp>/`, `evidence/alphafold_<sp>.tsv` |
| `esmfold.py` | ESMFold v1 for the 32 proteins AlphaFold DB lacks (re-run `structures.py` after) | `scratch/esmfold/<sp>/`, `evidence/esmfold_<sp>.tsv` |
| `predict.py` | fpocket 4.0 + P2Rank 2.5.1 (`gradi-pockets`, across a process boundary) | `evidence/pocket_list_<sp>.tsv` (long) |
| `holo.py` | drug-like ligand classes + DIAMOND vs BioLiP holo chains (`gradi-ortho`) | `evidence/{ligand_classes,holo_<sp>}.tsv` |
| `pdb_coverage.py` | experimental structures per protein, ligand or not (DIAMOND vs `pdb_seqres`) | `evidence/pdb_<sp>.tsv` |
| `alphafill_check.py` | measurement only: does AlphaFill add holo evidence? (no) | `evidence/alphafill_comparison.tsv` |
| `merge.py` | the deliverable + its checks | `pockets_<sp>.tsv` |

Deliverable, complete and canonical: `uniprot_ac · p2rank_score · fpocket_score ·
p2rank_n_pockets · holo_identity · pdb_n_structures · pdb_coverage · af_plddt`. Load through `src/pockets.py`.
Details, measurements and traps: `docs/pockets.md`.

The v1 cartoons (`legacy/scripts/06n_structure_snapshots.py`, `gradi-pymol`) are not ported.
