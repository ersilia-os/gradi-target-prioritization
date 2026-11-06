# pockets — structural ligandability

Can a small molecule bind this protein's fold? Three kinds of answer, never merged: **predicted**
pockets on AlphaFold models (fpocket, P2Rank), **measured** ligands in this protein's own PDB
structures, and ligands **modelled** onto it by AlphaFill.

Run in order, all with the `gradi` env:

| script | does | writes |
|---|---|---|
| `structures.py` | AlphaFold DB v6 models, validated against the proteome sequence | `data/source/alphafold/<sp>/`, `evidence/alphafold_<sp>.tsv` |
| `esmfold.py` | ESMFold v1 for the 32 proteins AlphaFold DB lacks (re-run `structures.py` after) | `scratch/esmfold/<sp>/`, `evidence/esmfold_<sp>.tsv` |
| `predict.py` | fpocket 4.0 + P2Rank 2.5.1 (`gradi-pockets`, across a process boundary) | `evidence/pocket_list_<sp>.tsv` (long) |
| `pdb_coverage.py` | which PDB chains ARE this protein (DIAMOND vs `pdb_seqres`, ≥ 95%) | `evidence/pdb_<sp>.tsv`, `evidence/pdb_chains_<sp>.tsv` |
| `alphafill.py` | fetch + flatten AlphaFill transplants; filters nothing | `data/source/alphafill/<sp>/`, `evidence/transplants_<sp>.tsv` |
| `holo.py` | the drug-likeness vocabulary, then both ligand counts | `evidence/{ligand_classes, ligands_pdb_<sp>, ligands_alphafill_<sp>, ligand_counts_<sp>}.tsv` |
| `merge.py` | the deliverable + its checks | `pockets_<sp>.tsv` |

Deliverable, complete and canonical: `uniprot_ac · p2rank_score · fpocket_score · n_ligands_pdb ·
n_ligands_alphafill · n_pdb_structures · af_plddt`. Load through `src/pockets.py`.
Details, measurements and traps: `docs/pockets.md`.

**We do not transfer ligands ourselves.** An earlier `holo.py` did — DIAMOND against BioLiP's holo
chains, a binding-site span test, and the best hit's % identity as a `holo_identity` column. That
was a hand-rolled AlphaFill and reached a fraction as far (184 Kp proteins against AlphaFill's
1,503). It and `alphafill_check.py`, whose "AlphaFill adds little" conclusion came from imposing an
identity floor AlphaFill does not use, were both deleted on 2026-10-03.

The v1 cartoons (`legacy/scripts/06n_structure_snapshots.py`, `gradi-pymol`) are not ported.
