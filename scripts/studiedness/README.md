# studiedness — no code yet

**How much is already known about this protein?** A target nobody has ever characterised is a
different proposition from one with fifty papers, and the axis is needed in both directions: an
unstudied protein is a risk, but it can also be the novelty the collaboration is looking for.

Candidate sources, none yet evaluated:
  - PubMed / Europe PMC citation counts per gene (gene2pubmed, or UniProt's own reference list)
  - UniProt **protein existence** level (1 evidence at protein level ... 5 uncertain) and whether
    the entry is reviewed — already fetched by `proteomes/download.py`, so this part is nearly free
  - annotation score / number of curated GO terms with experimental evidence codes
  - structures in the PDB, ligands in ChEMBL (the `ligands/` axis already measures the latter)

**The dark-proteome warning applies hardest here.** Kp HS11286 is a TrEMBL proteome: a low count
may mean the protein is unstudied, or that the *strain's* accession is unstudied while the family is
well known. Any count must be taken over orthologs (`orthology/`), not the accession alone, or the
axis will simply re-measure how dark the anchor is.

Deliverable: a **complete matrix**, one row per protein, canonical row order.
