# Stage 00 — reference proteomes and the identified species tables

`scripts/proteomes/download.py` · registry `src/proteome_registry.tsv` ·
outputs `data/processed/proteomes/`

Four species, **one reference proteome each**. That proteome is the unit of analysis for every
downstream stage; nothing else is an analysis organism.

| species | proteome | strain | n | locus tags |
|---|---|---|---|---|
| *K. pneumoniae* | `UP000007841` | HS11286 | 5,728 | `KPHS_*` |
| *E. coli* | `UP000000625` | K-12 MG1655 | 4,403 | b-numbers |
| *S. aureus* | `UP000008816` | NCTC 8325 | 2,889 | `SAOUHSC_*` |
| *H. sapiens* | `UP000005640` | — | 20,416 | — |

Human **must** be fetched with `reviewed:true`; the unfiltered proteome is 147,506 TrEMBL-bloated
entries.

## Why these strains

**HS11286** is the only *K. pneumoniae* proteome UniProt flags "Reference and representative", so it
has the most complete annotation and cleanest cross-references. Rejected: **ATCC 43816 / KPPR1** (no
curated UniProt reference proteome — it is a mouse-essentiality workhorse, not an annotation
reference) and **MGH 78578** (the historical reference; holds most of the species' reviewed entries,
but reviewed coverage is ~0.5% species-wide and propagates by orthology anyway).

**NCTC 8325** over the alternatives, by the same rule. Both NCTC 8325 and USA300 (`UP000001939`,
2,607) are UniProt Reference proteomes; NCTC 8325 is the classical genetic reference, has the deeper
annotation (2,889 vs 2,607), and is what AureoWiki and most *S. aureus* resources key on. **COL**
(`UP000000530`, 2,675) is tempting because the ClpP-activator proteomics is keyed on `SACOL_*` tags —
but it is *not* a reference proteome, and UniProt serves **zero entries** for it. COL is therefore a
tier-D bridge strain, not the anchor.

*S. aureus* is new in v2. It is the organism the activated-ClpP proteomics (Conlon 2013 ADEP4,
Jacques 2020 ONC212) is native to — v1 could only reach that data by cross-phylum DIAMOND RBH at
median 42% identity — and the only one of the three bacteria with ClpC/McsB, the machine every
published BacPROTAC actually targets.

## The problem this stage solves

The anchors are badly under-named. Measured against UniProt release **2026_02** (10 June 2026):

| field | Kp HS11286 | Ec K-12 | Sa NCTC 8325 |
|---|---|---|---|
| **gene name (primary)** | **1,055 = 18.4%** | 4,402 = 100% | **815 = 28.2%** |
| gene synonyms | 36 = 0.6% | 2,003 = 45.5% | 85 = 2.9% |
| ordered locus name | 100% | 100% | 98.4% |
| RefSeq · GeneID | 100% · 100% | 96.1% · 95.1% | 95.7% · 94.7% |
| KEGG · STRING | 100% · 87.5% | 97.2% · 92.8% | 95.4% · 93.0% |
| eggNOG · BioCyc | **0% · 0%** | 93.0% · 97.2% | 79.2% · 1.3% |
| InterPro · Pfam · PANTHER | 87.1 · 85.3 · 71.2% | 96.0 · 95.0 · 78.2% | 82.7 · 81.3 · 65.2% |
| PDB · AlphaFoldDB | 0.5% · 99.7% | 40.2% · 99.3% | 8.7% · 99.8% |
| reviewed | 7 = 0.1% | 100% | 815 = 28.2% |

Coverage tracks curation almost exactly — S. aureus has 815 reviewed entries and 815 gene names. The
*mappings* are fine; it is specifically **names and synonyms** that are missing, and those are what
literature and database lookups key on. In v1 a gene-symbol join on the raw anchor lost ~27% of known
essentials (`legacy/HISTORY.md` trap 59).

### Routes tested

1. **NCBI RefSeq annotation — dead end, measured.** The HS11286 GFF sets `Name=KPHS_00010` for all
   5,867 genes: **zero real symbols**. S. aureus NCTC 8325 yields only 219. NCBI is worse than
   UniProt here, not better.
2. **The rest of the species — this is the answer.** The Kp species tree (`taxonomy_id:573 AND
   gene:*`) returns 26,956 entries, of which **15,395 carry a primary gene name** over 11,973
   distinct sequences. S. aureus returns 17,562, of which **12,318** are named over 5,599 distinct
   sequences. (Note `gene:*` matches *any* gene field, including ordered locus names — hence the gap
   between pool size and actual donors.)
3. **UniRef90 clustering — free, and worth as much as exact matching.** UniProt's ID-mapping service
   maps `UniProtKB_AC-ID → UniRef90`, giving ≥90%-identity grouping with no alignment step and no
   `gradi-ortho` dependency.

## Name filling — three tiers, all inside the species

```
anchor            UniProt's own name on the anchor entry
species_exact     identical sequence, same species
species_uniref90  same UniRef90 cluster, same species
none
```

**Measured yield:**

| | Kp HS11286 | Sa NCTC 8325 |
|---|---|---|
| `anchor` | 1,055 (18.4%) | 815 (28.2%) |
| `+ species_exact` | +1,259 (22.0%) | +326 (11.3%) |
| `+ species_uniref90` | +1,316 (23.0%) | +148 (5.1%) |
| **named total** | **3,630 (63.4%)** | **1,289 (44.6%)** |
| still unnamed | 2,098 (36.6%) | 1,600 (55.4%) |
| synonyms gained | 175 | 159 |

Kp improves **3.4×**. E. coli is unchanged at 100% — it is the control: enrichment must add *zero*
names there, and if it adds any the donor join is wrong.

**Cross-species naming is deliberately NOT done here.** Inheriting an *E. coli* ortholog's name would
reach ~3,179 more Kp anchors, but it is an inference rather than an annotation, it belongs downstream
of real orthology where it can be labelled as such, and it would drag a DIAMOND dependency into
stage 00. As built, this stage needs nothing beyond `requests` and `pandas`.

### There is always exactly one preferred name

`gene_name` is a single **preferred** name, never an unordered set. When several donors disagree, a
**total order** decides, so a tie can never leave the field empty:

1. a name supported by a **reviewed (Swiss-Prot) donor** beats an unreviewed one — curation wins;
2. a **real name beats a `y###` systematic placeholder** — enterobacterial convention is that `yhdE`,
   `yqiE`, `yjgQ` are assigned before a function is known and survive only as synonyms once the gene
   is characterised;
3. the **more frequently attested** name;
4. the **shorter** name, then alphabetical — a final deterministic tiebreak.

Measured: 319 of 2,575 Kp fills were contested and **all 319 resolved**; 44 of 474 for S. aureus.
**Zero proteins are left with candidates but no preferred name.** The deciding rules for Kp:
non-placeholder over y-name 55, more-attested 21, shortest-then-alphabetical 25, reviewed donor 6
(at the `species_exact` tier). Worked examples — all of which are the correct modern name:

```
rdgB over yggV     maf  over yhdE     nudF over yqiE
lptG over yjgQ     kefA over yjeP     pldB over ytpA     hutC over yvoA
```

`gene_name_candidates` retains every alternative **in preference order**, and `name_audit.tsv`
records which rule decided each contested case, so a wrong preferred name is traceable rather than
mysterious. Rule 4 is the only arbitrary one; it fires on genuine paralog pairs such as
`trkG`/`trkH`, where no convention picks a winner.

One mechanical normalization is applied before any of this: automated annotation numbers duplicated
loci as `phoE_2`, `manX_1`, and that suffix is not part of the gene name (underscore-digit endings are
not valid bacterial nomenclature). Stripping it before voting stops `phoE` and `phoE_2` from looking
like a genuine conflict.

## Outputs

**Four tables at the top level, one per species, and nothing else.** Everything that supports them
is split between `evidence/` and `scratch/`.

```
data/processed/proteomes/
  proteome_kpneumoniae.tsv    5,728 x 9
  proteome_ecoli.tsv          4,403 x 9
  proteome_saureus.tsv        2,889 x 9
  proteome_human.tsv         20,416 x 9
  evidence/
    locus_tags_<species>.tsv   locus_tag + locus_tag_all
    annotation_<species>.tsv   the wide xref layer, free in the same request
    name_audit.tsv             every name fill: donor, candidates, contested, rule
    locus_bridge_strains.tsv   tier-D strain tags, for joining published screens
    registry.tsv  manifest.tsv
  scratch/
    .uniref90_<species>.json   the clustering cache
```
```
data/source/uniprot/proteomes/<label>.{fasta,tsv}   as fetched, + <label>.SOURCE.md
```

**The split is by the directory contract's own test — would you cite or check it, or would you
delete it to reclaim space?** The audits, the registry, the manifest, the locus tags and the wide
xref layer are all things a reader may need to *check*: `name_audit.tsv` is what makes a wrong
preferred name traceable, `registry.tsv` carries the required `why` for every row, and
`annotation_<species>.tsv` is cited whenever a downstream stage asks whether an xref exists. So
they are `evidence/`. The UniRef90 clustering cache is the one thing here that is purely a saved
API round-trip — regenerable, worth nothing to a reader, deleted without loss — so it is
`scratch/`. (These were one `accessory/` directory before the September 2026 reorganisation; the
reasoning that put them together, that they are supporting detail for a finished table rather than
staging artifacts on the way to something else, is still right — it just does not distinguish the
two tiers, which is what the cite-or-delete test adds.)

### The 9 columns

```
uniprot_ac         primary key
is_reviewed        near-useless for Kp (0.1%) but real for Sa (28.2%)
gene_name          the PREFERRED name — exactly one, always
gene_name_source   anchor | species_exact | species_uniref90 | none
gene_synonyms      literature uses old names (e.g. yggV for rdgB)
protein_name       the only readable handle for the 37% Kp / 55% Sa still unnamed
sequence           kept in-table so a single file is self-contained
refseq  geneid     the bridges to NCBI and literature lookup
```

Identity and nothing else. Everything dropped from the first run's 17 columns was measured, not
guessed: `taxid` and `species` had exactly one distinct value per file (they existed only for a
stacked parquet); `sequence_md5` was verified identical to `md5(sequence)`, which stays; `length` is
`len(sequence)`; `gene_name_donor` + `gene_name_candidates` moved to `evidence/name_audit.tsv`.

**The locus tags moved to `evidence/locus_tags_<species>.tsv`** rather than being discarded, since
they are the join key for published bacterial data. Coverage there: Kp 5,728/5,728 · Ec 4,402/4,403 ·
Sa 2,844/2,889 · human 0/20,416. Join them back with `src.proteomes.with_locus_tags(species)`.

**Two files were deleted outright**, not moved: `proteins.parquet` (19 MB, verified a pure concat of
the four tables) and `id_bridge.tsv` (8 MB, a melt of columns already present). Between them they
carried **zero** new information. `src/proteomes.py` provides `load`, `load_all`, `load_locus_tags`,
`with_locus_tags`, `load_annotation` and `id_bridge` instead.

### `evidence/annotation_<species>.tsv`

`kegg` · `string` · `embl` · `eggnog` · `biocyc` · `interpro` · `pfam` · `panther` · `go_id` · `ec` ·
`protein_families` · `pdb` · `alphafolddb`.

These come down in the *same* stream request at zero marginal cost, and UniProt's InterPro xref is as
complete as v1's dedicated InterPro stage (87.1% vs 86.8% on Kp), so banking them now avoids
re-streaming later. **`eggnog` and `biocyc` are 0% on Kp** — nothing downstream should assume
otherwise.

## The registry

`src/proteome_registry.tsv`, 41 rows, four tiers. Every row carries a **required `why`** naming the
dataset or question that justifies it — that column is what keeps the registry from decaying into
clutter. All 30 UniProt `expected_n` values were validated against the live API.

| tier | rows | what | fetched by default |
|---|---|---|---|
| **A** anchors | 4 | the four species proteomes | yes |
| **B** name donors | 3 | same-species pools (Kp 26,956 · Ec 65,479 · Sa 17,562) | yes |
| **C** comparator panel | 26 | v1's curated orthology panel, **with the 17 species it left as bare taxids now pinned** to explicit proteome IDs so the panel cannot drift when UniProt reassigns a reference | `--tier C` |
| **D** bridge strains | 8 | strains published screens are keyed on; **NCBI only** | `--tier D` |

### Tier D and the reference-proteome-only finding

**UniProt serves UniProtKB entries only for proteomes it flags "Reference proteome". Every other
proteome record returns zero entries.** Verified: ECL8, KPNIH1, KPPR1, NJST258_1/2, Kp52.145,
BW25113, W3110, UTI89, COL, Mu50, Newman, N315, MW2 and TCH1516 all return **0**.

That zero-entry list is precisely the set of strains published screens key on, so they come from
**NCBI RefSeq** instead (`datasets v2alpha .../download?include_annotation_type=PROT_FASTA,GENOME_GFF`)
and exist **only as identifier bridges** — never as analysis organisms.

| strain | assembly | unlocks |
|---|---|---|
| Kp ECL8 | `GCF_000315385.1` | Eichelberger/Short 2024 TraDIS |
| Kp KPNIH1 | `GCF_000281535.2` | Ramage 2017, the 424-gene essential set |
| Kp KPPR1 / ATCC 43816 | `GCF_000742755.1` | Bachman 2015/2025, Mike & Bachman 2023 in-vivo (`VK055_*`) |
| Kp MGH 78578 | `GCF_000016305.1` | iYL1228 FBA model + Jana 2023 CRISPRi library (`KPN_*`) |
| Kp NJST258_1 | `GCF_000598005.1` | Cain 2017 |
| Ec BW25113 | `GCF_000750555.1` | Keio, RB-TnSeq, most CRISPRi — low priority, b-numbers map exactly onto MG1655 |
| **Sa COL** | **`GCA_000012045.1`** | Conlon 2013 + Jacques 2020 activator proteomics (`SACOL_*`) |
| Sa Mu50 | `GCF_000009665.1` | Jacques 2020's `Mu50_homolog` column |

Two traps, both verified:

- **ECL8's GFF carries `old_locus_tag=BN373_00001,KPNEcl8_00001`** — both legacy namespaces. This is
  the bridge v1 ran Prokka to reconstruct, so it should retire the `gradi-prokka` environment.
- **COL must come from GenBank, not RefSeq.** `GCF_000012045.1` re-tagged everything `SACOL_RS#####`
  and retains only **60** original `SACOL####` tags; `GCA_000012045.1` retains **2,711**.

The other 10 Enterobacteriaceae-TraDIS genomes (RH201207, EC958, NCTC13441, ICC168, Ty2, A130,
D23580, SL3261, SL1344, P125109) are deliberately **not** in the registry: the compendium is consumed
through the column names in its own `giant-tab_final.tsv`, so it never needed per-genome proteomes.

## Figures

`scripts/plots/proteomes.py` writes `output/plots/proteomes/gene_names.png`. Three panels, all
read from this stage's own outputs (`gene_name_source` and `evidence/name_audit.tsv`) — nothing is
recomputed.

**A — coverage by fill tier.** Stacked, so the first segment is what UniProt supplied unaided and the
distance to the grey is what the fill bought:

| species | anchor | + species_exact | + species_uniref90 | = named | unnamed |
|---|---|---|---|---|---|
| *K. pneumoniae* | **18.4%** | 22.0% | 23.0% | **63.4%** | 36.6% |
| *E. coli* | 100.0% | — | — | 100.0% | 0.0% |
| *S. aureus* | **28.2%** | 11.3% | 5.1% | **44.6%** | 55.4% |
| *H. sapiens* | 99.3% | — | — | 99.3% | 0.7% |

The two donor tiers more than **triple** Kp's naming and roughly halve Sa's gap; *E. coli* and human
need no help at all, which is exactly why neither would have revealed the problem. The grey is left
grey on purpose — cross-species naming is deferred to a stage that has real orthology and can label
it as inference.

**B — how the 3,049 fills were decided.** 2,686 had a single candidate. The remaining **363 were
contested**, and the total order settled them: non-placeholder over a `y###` name (135), more
attested (116), shortest-then-alphabetical (72), reviewed donor (40). Every one is a row in
`evidence/name_audit.tsv` with the alternatives and the rule that fired, so no fill is a black box.

**C — what you can actually join on**, and the panel that matters operationally:

| species | gene_name | locus_tag | refseq | geneid |
|---|---|---|---|---|
| *K. pneumoniae* | 63.4% | **100.0%** | 100.0% | 100.0% |
| *E. coli* | 100.0% | **100.0%** | 96.1% | 95.1% |
| *S. aureus* | 44.6% | **98.4%** | 95.7% | 94.7% |
| *H. sapiens* | 99.3% | **0.0%** | 93.5% | 93.8% |

This is the figure behind CLAUDE.md's rule to join on `locus_tag`, not on gene name: on the anchor
the name column would silently drop **36.6%** of the proteome, and on *S. aureus* **55.4%**.

Two corrections to the way this has been stated: `locus_tag` is at 100% on Kp and Ec but **98.4% on
*S. aureus***, not 100% — so a locus-tag join is near-lossless, not lossless. And human has **no
locus tags at all** (0.0%), which is expected but means the rule is a *bacterial* rule; use
`uniprot_ac` for human, as everywhere else.

## Running it

```bash
python scripts/proteomes/download.py                      # tiers A + B
python scripts/proteomes/download.py --dry-run            # show the plan, fetch nothing
python scripts/proteomes/download.py --tier C --tier D    # comparator panel + bridge strains
python scripts/proteomes/download.py --only saureus__nctc8325__UP000008816 --refresh
```

Verbose by default (`-q` to quiet): a banner, one line per registry row with running counts, and a
closing coverage table plus `gene_name_source` tally. Idempotent — anything already on disk whose row
count matches is skipped unless `--refresh`.

### Guards

- **The count is verified before the stream is trusted.** `X-Total-Results` is read from a `size=1`
  search with the *same* query, and a streamed row count that disagrees raises. This is the guard
  against v1's most common failure: a short or empty payload behind an HTTP 200.
- **UniProt gzips when `compressed=true`** and a naive save yields a file that reads as UTF-8 garbage.
  Decompression is explicit, and the header row is asserted to parse.
- **The ID-mapping protocol has two sharp edges**: `/idmapping/status/{job}` answers 200
  `{"jobStatus":"RUNNING"}` and then **303 with a `Location` header** — following that redirect
  blindly yields a 400 — and results **paginate at 500 rows** via `Link: rel="next"`. Reading only the
  first response silently keeps 500 clusters out of N.
- **A `gene_name_source` outside the four allowed values raises**, so a cross-species name cannot leak
  in while that tier is deferred.
- **Every protein with any candidate gets a preferred name.** A row with a non-empty
  `gene_name_candidates` and an empty `gene_name` is a bug, and is asserted against.

## Run log

### 2026-09-01 — first full run, tiers A + B

UniProt release **2026_02** (10 June 2026). Wall clock ~45 min, of which ~20 min was fetching and the
rest UniRef90 clustering. All seven registry rows fetched with **exact** count agreement.

| species | n | anchor | + species_exact | + species_uniref90 | **named** | synonyms gained |
|---|---|---|---|---|---|---|
| *K. pneumoniae* HS11286 | 5,728 | 1,055 (18.4%) | +1,259 (22.0%) | +1,316 (23.0%) | **3,630 (63.4%)** | 175 |
| *E. coli* K-12 | 4,403 | 4,402 (100.0%) | +0 | +0 | **4,402 (100.0%)** | 180 |
| *S. aureus* NCTC 8325 | 2,889 | 815 (28.2%) | +326 (11.3%) | +148 (5.1%) | **1,289 (44.6%)** | 159 |
| human | 20,416 | 20,281 (99.3%) | — | — | **20,281 (99.3%)** | — |

**The E. coli control passed:** enrichment added **exactly zero** names. Its one unnamed protein,
`P75688` / `b0309` ("Putative uncharacterized protein b0309", 70 aa, reviewed), has no gene name in
UniProt at all and no same-sequence donor that carries one — a real gap, not a join failure.

Combined: `proteins.parquet` 33,436 rows · `id_bridge.tsv` 227,447 rows over 4 namespaces ·
`name_audit.tsv` 3,049 fills of which 363 contested, **all resolved to a preferred name**.

Coverage of the identity table (Kp): `uniprot_ac`, `is_reviewed`, `locus_tag`, `protein_name`,
`length`, `sequence_md5`, `taxid`, `species`, `refseq`, `geneid` all **100%**; `gene_name` 63.4%,
`gene_synonyms` 3.6%. Human `locus_tag` is 0.3% — expected, human has no ordered locus names.

### What the run demonstrated

`legacy/HISTORY.md` trap 58 records that **`clpA` and `sspB` gene symbols are absent from HS11286**,
so v1 could only find them by a sequence-based census. Both are now resolved by `species_exact`:

```
clpA  A0A0H3GQX4  KPHS_17930  <- donor A6T6Y0  (MGH 78578)
sspB  A0A0H3GTU5  KPHS_47670  <- donor A6TEN6  (MGH 78578)
```

Every Gr-ADI target and every Clp-machinery component resolves by gene name in all three bacteria:
`clpP clpX clpA clpS lon hslU hslV ftsH dnaK acpP gyrA gyrB sspB smpB`. **`clpC` is present and named
in *S. aureus* only** (`Q2G0P5` / `SAOUHSC_00505`), confirming the asymmetry the degradability axis
turns on — the Enterobacteriaceae do not have it.

### Output simplified, same numbers (later on 2026-09-01)

The first run wrote 13 files and a 17-column table. Measured, three of those files carried **27 MB and
one column of new information** between them, so the shape was cut to **4 tables + a supporting
directory** (now split `evidence/` / `scratch/`) and the table to **9 columns**, not the 11 this
line claimed until 2026-10-03 — `head -1 proteome_ecoli.tsv` returns 9 and *Outputs* has always
said 9. Re-run after the change: every figure identical —
Kp 3,630 / Ec 4,402 / Sa 1,289 / human 20,281 named, 363 contested and all resolved.

A third bug surfaced while doing it: **human `locus_tag` was literally `";"`** and `locus_tag_all`
`"; ; ; ;"`, because human has no ordered locus names and the raw UniProt field is just separators.
Now a separator-only field yields an empty string rather than punctuation posing as an identifier.

### Two bugs found and fixed during the run

1. **`gene_name_donor` pointed at the protein itself.** The donor index recorded "the first accession
   sharing this sequence", and the anchor's own entry is nearly always in the pool — it matches
   `gene:*` on its locus tag while carrying no primary name. So every attribution named the anchor and
   the audit trail was useless. Now the index stores a supporting accession *per gene name*, preferring
   a reviewed one. Verified: **0 of 3,049 attributions are self-referential.**
2. **No UniRef90 cache.** Re-running the assembly cost ~10 min of pure re-mapping. Now cached per
   species in `.uniref90_<species>.json`, keyed on an md5 of the accession set so a changed proteome
   invalidates it. A re-run is now seconds.

Also added: the clustering pass is skipped entirely when the anchor has no unnamed proteins. Note
E. coli does **not** hit that branch — it has exactly one — so it clusters once, then reads its cache.

### Residual gaps worth knowing

- **37% of Kp and 55% of S. aureus still have no gene name.** For those, `protein_name` (100%) plus
  `locus_tag` (100% / 98.4%) are the only readable handles. Cross-species transfer from E. coli would
  reach further but is deliberately deferred to a stage that has real orthology.
- **A second S. aureus ClpP-family protein is unnamed**: `Q2G2S1` / `SAOUHSC_01536`, "ATP-dependent
  Clp protease proteolytic subunit". Exactly the kind of protein the degradability axis cares about,
  and it has no gene name — a concrete reminder that the residual gap is not harmless.
- **No entry in the NCTC 8325 proteome is annotated as a Lon protease**, so `lon` resolves in Kp and
  E. coli but not S. aureus. Treat that as a UniProt annotation gap for this strain rather than
  evidence of absence.
