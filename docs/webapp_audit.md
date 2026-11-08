# Browser audit, 2026-11-08

Every column of the table was checked against the shipped payloads (`webapp/data/{ec,kp}.json`), the orthology and ligand tables, and, for the biology, the literature. Numbers are from those files. Nothing here changes a pipeline stage; the fixes made in the browser are listed first, the suggestions for the pipeline after.

## Fixed in the browser (this round)

- **Human ortholog "no" was shown for proteins that were never looked at.** 110 E. coli and 291 K. pneumoniae proteins are in no orthogroup (orthology evidence 1), so `has_human_ortholog = False` there means "could not look". The cell now shows `?`, and the "No human ortholog" badge needs evidence of at least 2.
- **Papers "0" for proteins with no similar curated protein** (1,961 K. pneumoniae, studiedness evidence 1) now shows a dash: it is "no donor", not "a donor with no papers".
- **Ligand evidence 3 with no potent ligand** (104 E. coli proteins) now says "Assayed, none potent" in its hover.
- COG `O` is now "Chaperones, turnover" (it hid protein turnover, the category that matters for a degrader).
- Hover texts rewritten with the measured limitations below.

## Findings, by column

**Degradability.** Ranks on length: Spearman with length is -0.70 (E. coli) and -0.73 (K. pneumoniae), so the top of the default sort is short proteins. The two probabilities are one opinion (rho 0.89). `nn_similarity` is raw cosine (median 0.93-0.94), not relatedness. The ortholog flags are empty for 75-85% of proteins. ADEP4 and ONC212 as ClpP activators with similar proteome-wide degradation in *S. aureus* is supported by Jacques et al. 2020 (ONC212 and ADEP4 caused similar proteome-wide degradation profiles).

**Essentiality.** The evidence level does not discriminate: 84% of E. coli and 69% of K. pneumoniae proteins are at 3, because a unanimous non-essential call counts. The mean consensus by level is not monotone (E. coli 0.46 / 0.73 / 0.47 for levels 1 / 2 / 3). For K. pneumoniae a 3 rests on other strains (3,974 proteins, none measured on HS11286). E. coli counts three K-12 screens of one lineage, while DEG lists two datasets for the strain, so "n of m screens" in the card refers to DEG only.
*Suggestion:* split the dot into "measured essential" vs "measured non-essential" vs "predicted only", the three states a reader cares about.

**Ligand precedence.** Sound. The score is pinned to 0 for 96-97% of proteins by design. Evidence 3 mixes "has potent ligands" with "assayed, none potent" and "a drug is recorded". **`n_drugs` is understated for K. pneumoniae**: the 95% identity floor splits close orthologs. E. coli to K. pneumoniae identity is 92% (gyrA), 94% (parC, parE), 92% (folA), 92% (murA), 94% (ftsI/PBP3), and those read 0 drugs while E. coli reads 27, 19, 3, 1 and 73. At a 90% floor, 15 more K. pneumoniae proteins gain a drug (23 at 85%), among them gyrA, parC, parE, folA, thyA, fabI, murA, ftsI, mrdA and several ribosomal proteins.
*Suggestion:* transfer drugs between E. coli and K. pneumoniae through the OrthoFinder ortholog, or lower the floor to 90% for this one column and document it.

**Pockets.** Tracks length (rho about 0.65). Evidence is provenance (3 = own PDB structure): E. coli 1,893 at 3, K. pneumoniae 569. Not a statement about pocket quality; the hover now says so.

**Localization.** Classes match DeepLocPro 1.0 (the Gram-negative runs use cytoplasm, inner membrane, periplasm, outer membrane, extracellular; the tool also has a cell-wall class, unused here). The evidence level tracks GO annotation more than accuracy: 73% of E. coli inner-membrane proteins reach 3 against 22% of cytoplasmic ones, since a cytoplasmic protein often carries no curated cellular-component term. Extracellular is the weakest class (E. coli 128 of 188 at level 1).
*Suggestion:* report the DeepLocPro confidence beside the class (it exists in `evidence/deeplocpro_<sp>.tsv`) and let the evidence ladder stop at 2 for cytoplasmic proteins without a GO term instead of penalising them.

**Function.** 15% (E. coli) to 21% (K. pneumoniae) of proteins have no COG hit. Function evidence 3 for 75% / 65%, so it is a coverage statement more than a quality one. GO slim now shows names. "Naming evidence" was dropped from the table as requested; it remains in the proteome table.

**Orthology.** Evidence 3 is reached by 63% of E. coli proteins; it needs both methods to find a bacterial ortholog, so it describes conservation as much as reliability. **`has_human_ortholog` has a recall gap and a precision gap**: E. coli `secY` has no human hit at all in `neighbors.tsv`, although Sec61a is its homolog (the mammalian homolog of Sec61p is similar to bacterial SecY, Cell 1992); `ftsY` has a 31% identity, 54% coverage hit to the SRP receptor and is called "no". The other direction: `gyrB` is a human ortholog on a 25% identity hit over 36% of the human chain, `yidC` on 25% over 34%. The column mixes true orthologs (`atpD` 72%) with remote family members. For selectivity the useful quantity is human identity over the binding region, which this column does not carry.
*Suggestion:* ship the best human identity and coverage as columns, and call "ortholog" only above a stated floor, or rename the column "human homolog".

**Studiedness.** The reference strain drives the evidence: E. coli K-12 reaches 3 for 88% of proteins, K. pneumoniae HS11286 for 6%, because HS11286 is an uncurated proteome whose entries have few papers of their own. That makes K. pneumoniae look less studied than it is (`KPC-2` is not novel). "Papers" is the best-cited homolog's count at 40% identity or more, not this protein's; the own-papers column was removed from the table. The score and the evidence invert for E. coli (consensus 0.50 at level 3, 0.54 at level 2).

## Evidence levels, overall

Three of the seven ladders say little because most proteins sit at 3 (essentiality 84% / 69%, studiedness 88% in E. coli, function 75% / 65%); two measure annotation richness rather than reliability (localization, function); one is provenance only (pockets). The two that read best are ligands and degradability's absence of one. A single shared 1-2-3 vocabulary hides that each means something different, which the hover text now spells out.
