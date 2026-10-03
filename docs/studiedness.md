# studiedness — how much is already known about this protein?

The axis is wanted in **both directions**: a target nobody has characterised is a risk, but it is
also the novelty the GraDi collaboration is looking for. `src.studiedness.novelty()` is the same
number read the other way.

**Deliverable:** `data/processed/studiedness/studiedness_<species>.tsv`, complete and canonical for
the three bacteria — `uniprot_ac · n_papers_own · n_papers_family · evidence`.

**The number is a paper count.** `n_papers_family` is the number of curated PubMed references on
the best-studied prokaryotic SwissProt homolog. Nothing is scaled, weighted or blended: a 5 is
five papers, a 0 is zero papers. See §6 for the 0–1 composite this replaced and why it went.

---

## 1. Why the anchors' own literature cannot carry the axis

Measured on UniProt 2026_03 (2026-09-21), fetching
`accession,reviewed,annotation_score,protein_existence,lit_pubmed_id`:

| | reviewed | PubMed ids per protein | distinct values | verdict |
|---|---|---|---|---|
| **Kp HS11286** | 0.1% | **5,710 of 5,728 have exactly 1** (the genome paper) | 7 | dead column |
| **Sa NCTC 8325** | 28.2% | **2,532 of 2,889 have none**; median 0 | 14 | nearly dead |
| **Ec K-12** | 100% | median 5, max 58 | 48 | rich |

99.7% of *K. pneumoniae* carries the identical count. This confirms v1's diagnosis
(`legacy/HISTORY.md:60` — "HS11286 is bibliometrically dark … pub counts flat at 1") and sharpens
it. So the usable number is **transferred from homologs**, which is also CLAUDE.md's standing rule:
*map by sequence, not by accession, when reaching an external database*.

`n_papers_own` still ships, because it *is* the measurement of darkness — and keeping it beside
`n_papers_family` is what makes "dark in *Klebsiella*, famous in *E. coli*" readable off one row.

**Correction to the old spec:** `scripts/studiedness/README.md` (now retired) claimed protein
existence was "already fetched by `proteomes/download.py`, so this part is nearly free". It was
not — `IDENTITY_FIELDS`/`ANNOTATION_FIELDS` at `scripts/proteomes/download.py:101-113` request
neither `protein_existence`, `annotation_score` nor `lit_pubmed_id`. `fetch.py` gets them.

---

## 2. The route, and the three that were measured against it

`transfer.py` DIAMONDs each proteome against **SwissProt (575,748 reviewed entries, ~14,000
species)** and transfers the literature of the best-characterised hit.

**Why not the four-species ortholog table already on disk?** Measured from `src.orthology.load`,
it reaches an E. coli or human ortholog for only **Kp 68.5% and Sa 46.6%** — *S. aureus* is
Gram-positive, so only 43.3% of it has an E. coli ortholog at all, and the free route structurally
abandons most of that proteome. The head-to-head ships in `evidence/route_comparison.tsv`; nothing
is merged.

**Why reviewed-only?** Curated literature is by definition where reviewed entries are; an
unreviewed entry has no curated reference list.

**Deferred, with the reason recorded so nobody re-derives it:**

- **OrthoDB group members** (92.7% Kp coverage) is the widest net but needs the ~46 GB re-stream.
  Revisit only if SwissProt coverage becomes the binding constraint.
- **NCBI Identical Protein Groups** would aggregate literature across strains sharing an exact
  sequence — but strain-identical proteins carry no *extra* literature, and the ≥95%
  `swissprot_direct` band already reaches the cases where a curated entry exists.
- **PubTator3** (`gene2pubtator3.gz`, 756 MB, text-mined rather than curated) is larger than
  gene2pubmed, but a format probe returns `P-glycoprotein`, `S100`, `mTOR`, `CD8` — visibly
  human-centric, and bacterial gene normalisation is weak.

---

## 2b. Donor scope: not eukaryotic, and all three scopes are measured

`--donor-scope prokaryotic` (the default) admits any donor whose UniProt taxonomic lineage does
**not** contain `Eukaryota (domain)` — Bacteria, Archaea **and phages**. `bacteria` is strict
Bacteria only; `any` lets a human, fungal or plant homolog donate. SwissProt splits
**337,294 Bacteria · 201,017 Eukaryota · 19,908 Archaea · 17,529 Viruses**.

**Why restrict at all.** Under `any`, five of K. pneumoniae's ten highest-scoring proteins took
human donors — HSPD1 (934 papers), HADHA, CTPS1, LONP1, AFG3L2 — so the claim became "this Kp
protein is well studied because its human mitochondrial homolog is". Reasonable for *is this
family understood?*, wrong for *is this a novel antibacterial target?*. It is also the mistake the
ligands axis measured and recorded: an unrestricted non-human bucket gave **424** apparent potent
Kp proteins against a true **175**. So the test is **lineage**, not organism name and not "is not
human" — strains are filed under their own taxids.

**Why NOT strict Bacteria, which is the obvious choice and is wrong.** Of the 93 proteins a strict
bacterial rule stranded with no donor at all, **70 (75%) lost theirs to a virus** — *Escherichia*
phage lambda and P1 — not to a eukaryote. A K. pneumoniae prophage protein whose best-characterised
relative is a lambda protein is emphatically not novel; lambda is among the most studied systems in
biology. Scoring it 0 fails in the one direction this axis must not fail. Admitting phages recovers
**62 of Kp's 77 stranded proteins and 12 of Sa's 16** while still excluding every eukaryotic donor.

**Measured, all three, every run** (`evidence/donor_scope_comparison.tsv`):

| species | scope | with a donor | lost vs `any` | different donor | rho vs `any` |
|---|---|---|---|---|---|
| Kp | **prokaryotic** | 3,767 | **15** | 215 | 0.989 |
| Kp | bacteria | 3,705 | 77 | 290 | 0.983 |
| Kp | any | 3,782 | 0 | — | 1.000 |
| Sa | **prokaryotic** | 1,613 | **4** | 129 | 0.991 |
| Sa | bacteria | 1,601 | 16 | 158 | 0.987 |
| Ec | prokaryotic / bacteria / any | 4,374 | 0 | 183 / 201 | 0.982 / 0.978 |

**The honest headline: the scope barely matters numerically.** The three scores correlate at
**rho 0.978–0.991** and only 215 of 5,728 Kp proteins change at all. The human donors were
visually striking and numerically almost irrelevant, because the log scale saturates above
P_REF = 204 papers — Kp `groEL` scores 1.000 from human HSPD1 (934 papers) *and* 1.000 from
E. coli `groEL` (246). The restriction is kept because it makes the number mean the right thing,
not because it moves it.

**The held-out control could not arbitrate, and `prokaryotic` wins on the tie-break:**

| scope | scored | spearman | pearson |
|---|---|---|---|
| **prokaryotic** | 2,945 (66.9%) | 0.328 | **0.350** |
| bacteria | 2,920 (66.3%) | 0.322 | 0.344 |
| any | 2,954 (67.1%) | **0.339** | 0.348 |

Spearman spans 0.017 across all three — a tie — and the two metrics *disagree* about the winner,
which is itself the signal that the difference is noise. `prokaryotic` has the best pearson, is
within noise on spearman, and strands 15 Kp proteins instead of 77.

**The evidence tiers are scope-consistent.** Under a restricted scope `below_floor` means "has an
in-scope hit below the floor" and `no_hit` means "no in-scope hit at all". Reusing the
unrestricted hit set would file a protein whose only curated relative is eukaryotic as though it
had a distant prokaryotic one — a much weaker novelty claim than the tier implies.

**Nothing is discarded.** `n_papers_family_prokaryotic`, `_bacteria` and `_any` all ship in
`evidence/transfer_<species>.tsv`, with `donor_scope` naming the shipped one and `donor_ac_any` /
`donor_organism_any` recording the donor the unrestricted scope would have chosen. Switching the
deliverable to another scope is a column swap, not a re-run.

---

## 3. NCBI gene2pubmed: it does not rescue the anchors, it sharpens the donors

Joined on the `geneid` column **already in `proteome_<species>.tsv`** (Kp 100% / Ec 95.1% /
Sa 94.7%) — no id-mapping round-trip. gene2pubmed is 83,245,982 rows over 43,745,777 distinct
GeneIDs; `gene2pubmed.py` reduces it once to the 273,002 GeneIDs this project can reach (anchors
plus every SwissProt donor), of which 247,961 (90.8%) carry at least one paper.

| | ≥1 paper | distinct values | max | UniProt `lit_pubmed_id` |
|---|---|---|---|---|
| **Ec** | 4,175 (94.8%) | **193** | **513** | 48 distinct, max 58 — **~4× coarser** |
| Kp | 5,342 (93.3%) | 12 | 11 | 7 distinct — both dead |
| Sa | 427 (14.8%) | 15 | 17 | UniProt reaches more: 357 |

Two conclusions:

1. **NCBI is dark on Kp and Sa too**, through a different door. It is not a fix for the own-literature
   problem and the axis does not pretend otherwise.
2. **It does sharpen the donor ranking** — 189 distinct values against UniProt's 47, and on
   *S. aureus* the union would lift coverage from 12.4% to 18.7% with 184 proteins getting a
   first paper.

**But it is NOT the shipped number** — that changed on 2026-09-22. `n_papers_family` counts
UniProt curated references only, because "papers a curator read and used" is **one consistent
definition applied to every row**, whereas `max(curated, gene-linked)` silently switched
definition per protein: Kp `rpoB` took NCBI's 350 while ~11% of donors took SwissProt's number.
gene2pubmed remains a **measured alternative** in `evidence/` — the `interpro2go` precedent — and
both components ship in `evidence/transfer_<species>.tsv`. See §6.

---

## 4. Unknome: usable, with one trap that would have zeroed a whole proteome

Unknome v3 (Rocha, Jayaram, Stevens et al., *PLoS Biol* 21(8):e3002222, 2023) scores each
orthologous cluster by how much is annotated about its best-annotated member.

**USE THE CLUSTER TABLE, NEVER THE PER-PROTEIN TABLE.** The per-protein download
(`/download/prot_tsv`, 276 MB) reads **`knownness = 0.000` for all 1,882 *S. aureus* entries** —
including `clpP`, `rpoB`, `ftsZ`, `gyrB`, `dnaA`, `rplB` — while those same proteins' clusters carry
real values. Sa `clpP` sits in `UKP00027`, the *same cluster* as Ec `clpP`, which reads **10.6**.
Measured 2026-09-21. Joining the protein table would have silently zeroed a proteome while looking
like a completed run. The cluster table (8.4 MB, 15,589 rows) is the one to use.

**The join is through PANTHER, for all three species.** Unknome covers 143 species; E. coli (83333)
and S. aureus (93061) are both there and are exactly our anchor taxids, but **there is no
Klebsiella**. The `panther` xref is already on disk from stage 00, so nothing is downloaded but the
cluster table.

Coverage is capped by that xref: **Kp 71.2% · Ec 78.2% · Sa 65.1%**, and essentially every protein
that has a PANTHER family joins (4,079 of 4,081 on Kp). That cap is why Unknome is *evidence* and
not the axis.

**Two caveats before quoting a knownness.** It counts **GO terms, not papers** — which is what makes
it a genuinely independent second opinion, and also means the two are not interchangeable. And the
clusters span 143 mostly-eukaryote species, so a bacterial protein's value can be set by a plant
homolog: FtsZ's `best_known_protein` is `FTSZ1_ARATH`. Read `unknome_best_known_protein`.

A protein with no PANTHER family gets an **empty** `unknome_knownness`, never a zero — a zero would
claim the family is unstudied, a much stronger statement than "no PANTHER family".

---

## 5. REJECTED — training ESM-C on Unknome knownness

Considered 2026-09-21 and rejected. **Do not re-add it.** (The `interpro2go` precedent: a rejected
alternative is recorded, not silently dropped.)

1. **The effective training set is 15,589, not ~1.9M.** Knownness is a *cluster* property — every
   member of `UKP00027` carries ClpP's 10.6 — so the protein rows are a few thousand distinct
   labels, massively replicated. CV would have to be cluster-grouped or the score is pure leakage,
   the trap `essentiality/paralog_clusters.py` exists to prevent (ungrouped 0.8743 → grouped 0.8705).
2. **It is a lossy approximation of a lookup we can do exactly.** ESM-C encodes family, so
   sequence→knownness is "find similar proteins, average their knownness". DIAMOND does that
   exactly and *names* the donor, its organism and its identity. The opposite of the essentiality
   situation, where Kp has no measurements at all and prediction is the only option.
3. **Mean-reversion points the wrong way.** A regressor gives a genuinely uncharacterised protein
   that merely resembles a known family a middling score — systematically hiding exactly the novel
   proteins the axis exists to surface.
4. **Circularity.** Degradability, essentiality and projections already run on ESM-C; a fourth
   ESM-C-derived column would correlate with them for reasons unrelated to studiedness — the trap
   CLAUDE.md records at stage 04 ("the embedding already encodes them").

---

## 6. The number, and the composite score it replaced

`n_papers_family` = **curated PubMed references on the best-studied prokaryotic SwissProt
homolog**. `n_papers_own` = the same count on this accession. Two integers, one definition, no
scaling. `src.studiedness.scaled()` derives a 0–1 version on the fly for anyone combining this
axis with the others — deliberately **not stored**, so there is one source of truth on disk.

### REJECTED — the 0–1 composite that shipped first

Until 2026-09-22 the axis shipped
`0.6 × min(1, log1p(n)/log1p(204)) + 0.4 × (annotation_score − 1)/4`. The project owner could not
interpret it, and three measurements showed that was the design's fault, not a documentation gap.
**Do not reintroduce it.**

1. **The same value meant different things.** At 0.5: **13 papers** if the donor's annotation
   score was 3, **83** if it was 1, **1** if it was 5. A protein with *zero* papers scored 0.4
   when its donor was annotation-5.
2. **The two halves double-counted.** UniProt's annotation score is largely a function of how much
   is known, so it correlated with the paper count at **r = 0.64–0.70**.
3. **The weights did the opposite of what the code claimed.** Literature was weighted 0.6
   "because it is the quantity the axis is named for", but the annotation component's spread was
   nearly double (sd 0.33–0.36 against 0.17–0.20), so **annotation swung the score more**. The
   0.6/0.4 split was invented and never measured.

Dropping the annotation term moved the ranking by spearman 0.86–0.91 — a real change, not a
relabelling.

### There is no ceiling — an earlier claim in this project was wrong

An earlier version of this document said SwissProt reference counts "saturate at 58, a
curation-practice ceiling". **That was false.** Verified directly against UniProt: the TSV export
is not truncated (human TP53 returns **225** PubMed ids), counts reach **225 overall and 119
among prokaryotic entries**, and the 58 is simply *E. coli* GroEL (`P0A6F5`) happening to be the
most-curated donor selected in all three species.

### The dynamic range is small, and the ties are in the right place

Bacterial proteins carry few curated references — median donor 4–6, 37–48 distinct values. What
matters for a prioritized shortlist is *where* the ties fall:

| | distinct | median (scored) | top of the ranking | largest non-zero tie |
|---|---|---|---|---|
| Kp | 46 | 6 | 58p ×1 · 56p ×2 · 50p ×1 | 466 proteins at 4 papers (12.4% of scored) |
| Ec | 48 | 6 | 58p ×1 · 56p ×2 · 50p ×1 | 625 at 4 papers (14.3%) |
| Sa | 37 | 4 | 58p ×1 · 46p ×2 · 38p ×1 | 305 at 2 papers (18.9%) |

The median column is over *scored* proteins only; §9's deliverable table gives the whole-proteome
median, which on Kp and Sa is dragged down by the unscored tail (Kp 4, Sa 1).

**The top is near-unique and the ties sit in the poorly-studied bulk** — exactly the right shape.
A `MAX_TIE_FRACTION` guard (25% of the *scored* proteins, zero excluded) exits the stage if that
ever stops being true.

### Why SwissProt-curated only, and not the union with NCBI

`gene2pubmed` is larger for **87–94% of donors, median 2.8×** (see §3), and it was in the shipped
number until 2026-09-22. It is now a measured alternative rather than the count, because
"references a UniProt curator read and used" is **one consistent definition** applied to every
row, whereas `max(curated, gene-linked)` silently switched definitions per protein — for Kp
`rpoB` the shipped number came from NCBI (350) and for ~11% of donors from SwissProt. Both
components still ship in `evidence/transfer_<species>.tsv`.

## 7. Traps, each hit at least once during the build

**1. `--max-target-seqs` silently biased the whole axis, and it was the worst bug here.** DIAMOND
returns the top *k* hits **by bitscore**, i.e. the CLOSEST relatives — but this stage wants the
**best-cited** homolog, which for a conserved protein is a distant model-organism entry. At k=50
the candidate list for *S. aureus* GroEL was 50 Staphylococci and Bacilli, the chosen donor was
*B. subtilis* GroEL with **9 papers**, and *E. coli* GroEL (246 papers) never appeared at all.
The bias is systematic and points the wrong way: it understated studiedness **precisely for the
most conserved, most studied families**. Every spot-check gene showed `n_candidates == 50`, which
is what exposed it.

Fixed in two stages, both measured: k=500 globally, then a **deep second pass at k=5,000 over
only the 1,403 queries (12.6%) still capped** — raising k globally to 5,000 would be 65M rows.
After the deep pass **0 of 1,403 remain capped**, so every query now has its complete homolog
list. Effect on the held-out control, at every floor — **measured while the axis still shipped the
0–1 composite, so the absolute values are superseded** (the paper-count control reads 0.2909 /
0.3280 / 0.2507 at the same three floors; see §8). What this table establishes is the *ordering*,
which the method change does not touch: k=500 beats k=50 and the deep pass beats k=500, at every
floor.

| floor | k=50 | k=500 | k=500 + deep |
|---|---|---|---|
| 25% | 0.3265 | 0.4248 | 0.4537 |
| **40%** | 0.4666 | 0.5338 | **0.5436** |
| 95% | 0.5901 | 0.5922 | 0.5922 |

Sa `groEL` went from the 85.8th percentile to the 98.2nd. `n_candidates` ships per protein and
the run prints the capped fraction — if it is ever large again, raise `--deep-targets`.

**2. Anchoring the log scale at P95 made half of E. coli an unrankable tie.** See §6.

**3. Preferring the closest band over the most-cited donor is worse, not better.** It looks like
the safer rule and was tried: held-out control **0.4426 against 0.4666** at the 40% floor
(composite-era numbers, like trap 1's — the comparison is the finding, not the absolute level). The
quantity being carried is literature volume, and the best-cited homolog predicts "this family is
studied" better than the nearest one does.

**4. Unknome's per-protein table is zero for all of *S. aureus*.** See §4.

**5. "Restrict donors to Bacteria" is the obvious rule and it is wrong** — it strands prophage
proteins whose best-characterised relative is a lambda or P1 protein. Only looking at *which*
donors were lost exposed it: 70 of 93 went to viruses, not eukaryotes. See §2b.

**6. A biological control calibrated on an assumption instead of a measurement.** The spot check
first demanded the top decile and failed on `rplB` — but ribosomal protein L2 genuinely has ~40
papers against RNA polymerase's 350. It is well studied *as part of the ribosome*, not as itself,
and 77th–90th percentile is the axis being right. The bar is now the 70th percentile, set to
catch a **broken join** (which scatters the eight genes to random percentiles; all eight clearing
the 70th by chance is p ≈ 7e-5) rather than to enforce an intuition. This is the same lesson
CLAUDE.md already records for the essentiality controls, where a 0.50 separation bar would have
failed the gold standard.

---

## 8. Controls

**The E. coli held-out transfer control** — the one that makes the axis believable. *E. coli* is
the only anchor whose own literature is real, so it is the only place transfer can be tested at
all. Every **Escherichia** taxid (90 of them, 24,423 SwissProt entries) is struck out of the donor
pool and `n_papers_family` recomputed, then correlated against E. coli's own measured
`n_papers_own`.

**Measured under the shipped `prokaryotic` scope: spearman 0.328, pearson 0.350 (on log1p, since
both sides are skewed integers), over 2,945 of 4,403 proteins (66.9%) that still find a
non-Escherichia donor.** Floor **0.25**, frozen in `transfer.py`; the run exits non-zero below it.
All three donor scopes are scored (see §2b) and land in `evidence/control_ecoli_summary.tsv`
(prokaryotic 0.328 · bacteria 0.322 · any 0.339 — still a tie).

**It reads lower than it used to, and that is a circularity being removed, not a regression.**
While the axis shipped the 0–1 composite this control read 0.5411. But the blend put a
0.4-weighted annotation term on *both* sides, and own-vs-donor annotation score correlates at
**spearman 0.950** — so a large part of that number was annotation agreeing with itself rather
than literature transferring. Counting papers only gives **0.328**, which is the honest figure.
The floor was re-derived accordingly; 0.25 leaves headroom while staying far above chance, since
the guard exists to catch a broken join (which would read ~0), not to certify a particular rho.

**DO NOT MAXIMISE IT** — but note that the shipped floor is nevertheless its optimum. The control
correlates two deliberately different quantities: a protein with 3 papers of its own whose homolog
has 300 *should* score low on `own` and high on `family`, and that gap is the entire point of the
axis. So the control has a ceiling well below 1 and is a **floor-check, not an objective** — it
exists to catch a broken join (which reads ~0), not to certify a particular rho.

**The identity floor is 40%, and the control's own optimum agrees with it**
(`evidence/floor_sensitivity.tsv`, recomputed every run):

| identity floor | Kp with a donor | Ec | Sa | control spearman | control n |
|---|---|---|---|---|---|
| 20% | 78.0% | 99.3% | 70.3% | 0.2804 | 3,571 |
| 25% | 77.4% | 99.3% | 68.8% | 0.2909 | 3,532 |
| 30% | 73.8% | 99.3% | 65.4% | 0.3149 | 3,369 |
| **40%** | **65.8%** | **99.3%** | **55.8%** | **0.3280** | **2,945** |
| 60% | 56.5% | 99.3% | 40.5% | 0.2767 | 2,411 |
| 95% | 19.6% | 99.3% | 35.9% | 0.2507 | 1,787 |

**rho PEAKS at the shipped 40% floor and falls away on both sides, reaching its LOWEST value
(0.2507) at 95%.** So 40% is not a compromise between coverage and the control — it is the
control's best point, and it is also CLAUDE.md's standing transfer threshold, which is the reason
it was tried first. The coverage trade-off either side is real and monotonic: 25% buys Kp 77.4%
at rho 0.2909; 40% gives Kp 65.8% at 0.3280; 60% drops to Kp 56.5% and 0.2767, paying 9 points of
Kp coverage for a *worse* control.

**An earlier version of this document (and of CLAUDE.md) claimed the opposite** — that rho rose
monotonically with the floor, reaching 0.5922 at 95%, so that maximising it would drive the floor
to 95% where `family` collapses onto `own`. That was **composite-era residue**: those numbers came
from the 0–1 blend, whose annotation term correlated with itself across the join (own-vs-donor
annotation score, spearman 0.950) and rose with identity for that circular reason. Counting papers
only, the circularity is gone and the shape inverts. The *caution* survives the numbers and is why
the paragraph is kept — the control is still not a thing to maximise — but the empirical claim
that it conflicts with the shipped floor is withdrawn.

**Decoy calibration** — composition-preserving shuffles of our own 13,020 sequences, searched
against the same database. **0 of 13,020 match at every floor from 20% to 95%.** So the floor is
not defending against spurious homology; DIAMOND's e-value already does.

**Spot checks** — eight named workhorses, which is the only control that catches a join to the
wrong accessions (that would still produce a well-formed table). Measured on the shipped run:
All eight land above the 70th-percentile bar in all three species. Kp `groEL` reads
**1 paper of its own, 58 on its family** (donor *E. coli* GroEL, `P0A6F5`) — the axis's whole
point in one row.

**Completeness and canonical order** — enforced through `M.reindex()`, which refuses to invent a
missing protein. `python -m src.matrices` reports all three studiedness tables canonical.

---

## 9. Results

### The deliverable

**The five tiers and their thresholds.** `evidence` is nested on identity: **`swissprot_direct`
(≥95%) ⊃ `swissprot_close` (≥60%) ⊃ `swissprot_homolog` (≥40%)**, then the two zeros —
**`below_floor`** (an in-scope hit exists but under the 40% floor) and **`no_hit`** (nothing
in-scope in SwissProt at all). Both score 0 and the axis refuses to invent a number for either,
but they are different claims: `no_hit` is the strongest novelty statement this axis makes, while
`below_floor` merely has a too-distant curated relative. Without the split, a third of
*K. pneumoniae* is one undifferentiated tie.

**The hard ceiling is DIAMOND's own hit rate, and it is independent of the floor**: a SwissProt
hit *of any kind* exists for only **81.6% of Kp and 73.0% of Sa** (99.4% of Ec) — i.e.
`direct + close + homolog + below_floor`. No choice of floor can take the axis above that; moving
the floor only shuffles proteins between the scored tiers and `below_floor`.

| species | n | with a donor | direct | close | homolog | below_floor | no_hit | median own | median family | max family | distinct family |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Kp** | 5,728 | 3,767 (65.8%) | 364 | 2,491 | 912 | 903 | 1,058 | **1** | **4** | 58 | 46 |
| **Ec** | 4,403 | 4,374 (99.3%) | 3,855 | 154 | 365 | 4 | 25 | 5 | 6 | 58 | 48 |
| **Sa** | 2,889 | 1,613 (55.8%) | 332 | 225 | 1,056 | 493 | 783 | **0** | **1** | 58 | 37 |

**`family − own` is the axis working, and the contrast is the clearest thing here: +3 papers
(median) on K. pneumoniae against +0 on E. coli**, whose own literature already *is* its
family's. On Kp the median protein has **1** paper of its own and **4** on its family; the top of
the ranking reaches 58.

Largest non-zero tie: Kp 466 proteins at 4 papers (12.4% of scored) · Ec 625 at 4 (14.3%) ·
Sa 305 at 2 (18.9%) — all under the 25% guard, and all in the poorly-studied bulk rather than at
the top.

### Loading it

**Load through `src/studiedness.py`**, never by path — the evidence tables are the ones that move.

| | |
|---|---|
| `load(species)` · `load_all()` | the deliverable, one species or all three stacked |
| `novelty(species)` | the same number read the other way round |
| `scaled(n_papers)` | the 0–1 version, derived on the fly and **deliberately not stored** |
| `definition()` | what the shipped number counts, as data |
| `load_own` · `load_transfer` | the two components in `evidence/` |
| `load_unknome` · `load_gene2pubmed` | the two measured alternatives |
| `load_route_comparison` · `load_floor_sensitivity` · `load_donor_scope_comparison` | the three sweeps |
| `control()` | the E. coli held-out summary |
| `manifest()` | provenance and row counts |

### Route comparison (`evidence/route_comparison.tsv`) — nothing is merged

| | SwissProt DIAMOND | free 4-species orthologs |
|---|---|---|
| Kp | 65.8% | **68.5%** |
| Ec | **99.3%** | 22.3% |
| Sa | **55.8%** | 46.6% |

**Read this honestly: the free route reaches slightly more of *K. pneumoniae*.** The two are not
like-for-like — the ortholog route applies no identity or coverage floor at all, while SwissProt
requires ≥40% identity and ≥50% coverage both ways. SwissProt wins decisively on the other two
(and on E. coli by 4.5×), and it is the only route that yields a *literature count* rather than
just the existence of an ortholog. The free route is measured because the *Identifier mapping*
rule says to measure every route; it is not used.

### Unknome agreement (`evidence/unknome_agreement.tsv`)

| species | n both | spearman vs `n_papers_family` | median family where Unknome present / absent |
|---|---|---|---|
| Kp | 4,079 | +0.305 | 5 papers / 0 (1,649 absent) |
| Ec | 3,444 | +0.398 | 6 / 3 (959 absent) |
| Sa | 1,882 | +0.288 | 3 / 0 (1,007 absent) |

**Verdict: correlated but not redundant — keep it as evidence, do not promote it.** At rho
0.29–0.40 it agrees on direction while measuring a genuinely different quantity (GO annotation
depth, not papers), so it is a usable second opinion for a protein where the literature route is
weak. It cannot become the axis: coverage caps at the PANTHER xref, and the proteins it cannot
reach are exactly the ones with `n_papers_family = 0` — it is silent precisely where an
independent opinion would be most valuable.

### The consortium's own panel (`output/plots/studiedness/interest_panel.png`)

**The GraDi envelope targets are not novel.** `src/interest.py`'s panel sits at the **80th
percentile (median) on Kp, 74th on Ec, 69th on Sa** — well studied, as would be expected of
LPS/Lpt/Bam/Sec biology. Read `interest.coverage()` first: the panel matches by gene symbol, and
only 43 of its genes are named in Kp and 11 in Sa.

---

## 10. Run log

| date | what |
|---|---|
| 2026-09-21 | Axis built. `fetch.py` (274 s, all payloads verified against `x-total-results`), `gene2pubmed.py` (83.2M rows streamed in 71 s), `unknome.py` (seconds), `transfer.py` (~12 min incl. two DIAMOND passes), `merge.py` (seconds). UniProt release **2026_03 (02-September-2026)**; SwissProt fasta dated 2026-09-03; gene2pubmed dated 2026-09-21. |
| 2026-09-22 | **The 0-1 composite was replaced by a plain paper count** after the project owner could not interpret it; see §6. The control floor was re-derived (0.45 -> 0.25) because dropping the blend also removed a circularity in the control. Figures became survival curves. |
| 2026-09-21 | Donor scope added. SwissProt metadata refetched with `lineage` (432 s, 30.1 MB). Three scopes computed every run; `prokaryotic` shipped. The strict-Bacteria rule was implemented first and rejected on the phage measurement in §2b. |

**Order matters:** `fetch.py` → `gene2pubmed.py` → `unknome.py` → `transfer.py` → `merge.py`.
`transfer.py` reads the gene2pubmed counts cache for donor scoring, and will run without it while
printing that donor literature is UniProt-only and ~4× coarser. **~20 min cold**, and every step
caches, so a re-run is a re-parse.

**Environment: the `gradi` env throughout** — no process boundary on this axis. DIAMOND is
borrowed from **`gradi-ortho`** via `GRADI_DIAMOND_BIN`, the standing pattern (there is no
osx-arm64 DIAMOND build, and installing it into `gradi` would flip the env to osx-64 and take
ESM-C down with it).

**CLI.** All five scripts take `--species` · `--refresh` · `--dry-run` · `-q`. `transfer.py`, which
does the work, adds:

| flag | |
|---|---|
| `--limit N` | smoke test — **writes only to `scratch/`**, so it can never clobber a full run |
| `--donor-scope {prokaryotic,bacteria,any}` | §2b; all three are computed every run regardless |
| `--max-targets` · `--deep-targets` | the two DIAMOND k values behind trap 1 (500, then 5,000) |
| `--threads` · `--sensitivity` | DIAMOND tuning |
| `--no-control` | skips the E. coli held-out control and its non-zero exit |

**Figures:** `scripts/plots/studiedness.py` (stylia) → `output/plots/studiedness/studiedness.png`,
`control.png`, `interest_panel.png`. **Run plot scripts one at a time** — stylia clears the
matplotlib font cache at import, so two concurrent plot scripts delete each other's and die with a
`FileNotFoundError` naming `~/.matplotlib`.

**Deletable after a run:** `data/source/uniprot/literature/uniprot_sprot.fasta.gz` (89.5 MB),
`data/source/ncbi/gene/gene2pubmed.gz` (273.8 MB) and everything in
`data/processed/studiedness/scratch/` (805 MB, mostly the DIAMOND hit tables). All are public and
re-derivable from the `SOURCE.md` files — **never upload them to eosvc.**
