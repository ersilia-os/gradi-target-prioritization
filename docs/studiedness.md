# studiedness — how much is already known about this protein?

The axis is wanted in **both directions**: a target nobody has characterised is a risk, but it is
also the novelty the GraDi collaboration is looking for. `src.studiedness.novelty()` is the same
number read the other way.

**Deliverable:** `data/processed/studiedness/studiedness_<species>.tsv`, complete and canonical for
the three bacteria — `uniprot_ac · studiedness_own · studiedness_family · evidence`.

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

`studiedness_own` still ships, because it *is* the measurement of darkness — and keeping it beside
`studiedness_family` is what makes "dark in *Klebsiella*, famous in *E. coli*" readable off one row.

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
| Kp | **prokaryotic** | 3,767 | **15** | 254 | 0.990 |
| Kp | bacteria | 3,705 | 77 | 329 | 0.985 |
| Kp | any | 3,782 | 0 | — | 1.000 |
| Sa | **prokaryotic** | 1,613 | **4** | 152 | 0.990 |
| Sa | bacteria | 1,601 | 16 | 179 | 0.986 |
| Ec | prokaryotic / bacteria / any | 4,374 | 0 | 224 / 240 | 0.986 / 0.983 |

**The honest headline: the scope barely matters numerically.** The three scores correlate at
**rho 0.983–0.990** and only 254 of 5,728 Kp proteins change score at all. The human donors were
visually striking and numerically almost irrelevant, because the log scale saturates above
P_REF = 204 papers — Kp `groEL` scores 1.000 from human HSPD1 (934 papers) *and* 1.000 from
E. coli `groEL` (246). The restriction is kept because it makes the number mean the right thing,
not because it moves it.

**The held-out control could not arbitrate, and `prokaryotic` wins on the tie-break:**

| scope | scored | spearman | pearson |
|---|---|---|---|
| **prokaryotic** | 2,945 (66.9%) | 0.5411 | **0.5351** |
| bacteria | 2,920 (66.3%) | 0.5375 | 0.5319 |
| any | 2,954 (67.1%) | **0.5436** | 0.4998 |

Spearman spans 0.006 across all three — a tie — and the two metrics *disagree* about the winner,
which is itself the signal that the difference is noise. `prokaryotic` has the best pearson, is
within noise on spearman, and strands 15 Kp proteins instead of 77.

**The evidence tiers are scope-consistent.** Under a restricted scope `below_floor` means "has an
in-scope hit below the floor" and `no_hit` means "no in-scope hit at all". Reusing the
unrestricted hit set would file a protein whose only curated relative is eukaryotic as though it
had a distant prokaryotic one — a much weaker novelty claim than the tier implies.

**Nothing is discarded.** `studiedness_family_prokaryotic`, `_bacteria` and `_any` all ship in
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
2. **It materially sharpens the DONOR ranking**, which is what `studiedness_family` consumes — and
   the union earns its keep on *S. aureus* as well, taking it from 12.4% to **18.7%** with 184
   proteins getting their first paper.

So `n_pubs` is the **union** of the two sources (the `max`, not the sum — they overlap heavily and
adding them would double-count). Both components ship, so the union can be undone.

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

## 6. The scale

Both columns are
`0.6 * min(1, log1p(n_pubs) / log1p(P_REF)) + 0.4 * (annotation_score - 1) / 4`, with **one fixed
global `P_REF`** over all 575,748 reviewed entries (`evidence/scale.tsv`).

**Not a per-proteome percentile** — that is the documented `geptop_score` trap, where a
proteome-relative score means 0.6 in Kp is not 0.6 in Sa. Here the three species *and* the two
columns are directly comparable.

**`P_REF` is the 99th percentile (204 papers), and P95 was measured and rejected.** The SwissProt
distribution is extremely skewed:

| P50 | P75 | P90 | P95 | P99 | P99.9 | max |
|---|---|---|---|---|---|---|
| 1 | 3 | 12 | 38 | **204** | 1,195 | 20,402 |

Anchoring at P95 clipped the top and left E. coli with a median `studiedness_own` of **exactly
1.000** — unrankable ties over half a proteome, the same failure CLAUDE.md records for v1's
under-regularised logistic in stage 04. At P99 only ~1% of SwissProt saturates, and a
**`saturated > 20%` guard** exits `transfer.py` if that ever changes.

---

---

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
list. Effect on the held-out control, at every floor:

| floor | k=50 | k=500 | k=500 + deep |
|---|---|---|---|
| 25% | 0.3265 | 0.4248 | 0.4537 |
| **40%** | 0.4666 | 0.5338 | **0.5436** |
| 95% | 0.5901 | 0.5922 | 0.5922 |

Sa `groEL` went from the 85.8th percentile to the 98.2nd. `n_candidates` ships per protein and
the run prints the capped fraction — if it is ever large again, raise `--deep-targets`.

**2. Anchoring the log scale at P95 made half of E. coli an unrankable tie.** See §6.

**3. Preferring the closest band over the most-cited donor is worse, not better.** It looks like
the safer rule and was tried: held-out control **0.4426 against 0.4666** at the 40% floor. The
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
pool and `studiedness_family` recomputed, then correlated against E. coli's own measured
`studiedness_own`.

**Measured under the shipped `prokaryotic` scope: spearman 0.5411, pearson 0.5351, over 2,945 of
4,403 proteins (66.9%) that still find a non-Escherichia donor.** Floor 0.45, frozen in
`transfer.py`; the run exits non-zero below it. All three donor scopes are scored (see §2b) and
land in `evidence/control_ecoli_summary.tsv`.

**DO NOT MAXIMISE IT.** The control correlates two deliberately different quantities — a protein
with 3 papers of its own whose human homolog has 300 *should* score low on `own` and high on
`family`, and that gap is the entire point of the axis. Worse, rho rises monotonically with the
floor (0.5922 at 95%) for a nearly circular reason: the exclusion removes *Escherichia* only, so
at a high floor the surviving donors are largely *Salmonella* and *Shigella* near-duplicates whose
publication counts track E. coli's because they are effectively the same proteins. Tuning the
floor to maximise rho would drive it to 95%, where `family` collapses onto `own` and the axis
stops doing anything.

**Decoy calibration** — composition-preserving shuffles of our own 13,020 sequences, searched
against the same database. **0 of 13,020 match at every floor from 20% to 95%.** So the floor is
not defending against spurious homology; DIAMOND's e-value already does.

**Spot checks** — eight named workhorses, which is the only control that catches a join to the
wrong accessions (that would still produce a well-formed table). Measured on the shipped run:
`rpoB`, `gyrB`, `ftsZ`, `clpP`, `dnaA`, `secA`, `groEL` all land at **99.0–99.5th percentile in
all three species**; `rplB` at 85.5 (Ec) / 89.1 (Kp) / 94.1 (Sa) for the reason in §7.

**Completeness and canonical order** — enforced through `M.reindex()`, which refuses to invent a
missing protein. `python -m src.matrices` reports all three studiedness tables canonical.

---

## 9. Results

### The deliverable

| species | n | with a donor | direct | close | homolog | below_floor | no_hit | median own | median family | family − own |
|---|---|---|---|---|---|---|---|---|---|---|
| **Kp** | 5,728 | 3,767 (65.8%) | 364 | 2,489 | 914 | 903 | 1,058 | 0.078 | 0.398 | **+0.292** |
| **Ec** | 4,403 | 4,374 (99.3%) | 3,908 | 122 | 344 | 4 | 25 | 0.589 | 0.613 | +0.000 |
| **Sa** | 2,889 | 1,613 (55.8%) | 248 | 268 | 1,097 | 493 | 783 | 0.000 | 0.224 | +0.046 |

Distinct `studiedness_family` values: Kp 325 · Ec 327 · Sa 213 — no saturation at either end.

**`family − own` is the transfer's contribution and it behaves exactly as the axis predicts:**
**+0.292 on K. pneumoniae**, where the anchor is dark and the family is not, and **+0.000 on
E. coli**, where the protein's own literature already *is* its family's. That contrast is the
clearest single number in the stage.

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

| species | n both | spearman vs `studiedness_family` | median family where Unknome present / absent |
|---|---|---|---|
| Kp | 4,079 | +0.371 | 0.598 / 0.000 (1,649 absent) |
| Ec | 3,444 | +0.501 | 0.670 / 0.302 (959 absent) |
| Sa | 1,882 | +0.328 | 0.381 / 0.000 (1,007 absent) |

**Verdict: correlated but not redundant — keep it as evidence, do not promote it.** At rho
0.35–0.52 it agrees on direction while measuring a genuinely different quantity (GO annotation
depth, not papers), so it is a usable second opinion for a protein where the literature route is
weak. It cannot become the axis: coverage caps at the PANTHER xref, and the proteins it cannot
reach are exactly the ones with `studiedness_family = 0` — it is silent precisely where an
independent opinion would be most valuable.

### The consortium's own panel (`output/plots/studiedness/interest_panel.png`)

**The GraDi envelope targets are not novel.** `src/interest.py`'s panel sits at the **86th
percentile (median) on Kp, 82nd on Ec, 74th on Sa** — well studied, as would be expected of
LPS/Lpt/Bam/Sec biology. Read `interest.coverage()` first: the panel matches by gene symbol, and
only 43 of its genes are named in Kp and 11 in Sa.

---

## 10. Run log

| date | what |
|---|---|
| 2026-09-21 | Axis built. `fetch.py` (274 s, all payloads verified against `x-total-results`), `gene2pubmed.py` (83.2M rows streamed in 71 s), `unknome.py` (seconds), `transfer.py` (~12 min incl. two DIAMOND passes), `merge.py` (seconds). UniProt release **2026_03 (02-September-2026)**; SwissProt fasta dated 2026-09-03; gene2pubmed dated 2026-09-21. |
| 2026-09-21 | Donor scope added. SwissProt metadata refetched with `lineage` (432 s, 30.1 MB). Three scopes computed every run; `prokaryotic` shipped. The strict-Bacteria rule was implemented first and rejected on the phage measurement in §2b. |

**Order matters:** `fetch.py` → `gene2pubmed.py` → `unknome.py` → `transfer.py` → `merge.py`.
`transfer.py` reads the gene2pubmed counts cache for donor scoring, and will run without it while
printing that donor literature is UniProt-only and ~4× coarser.

**Deletable after a run:** `data/source/uniprot/literature/uniprot_sprot.fasta.gz` (89.5 MB),
`data/source/ncbi/gene/gene2pubmed.gz` (273.8 MB) and everything in
`data/processed/studiedness/scratch/` (805 MB, mostly the DIAMOND hit tables). All are public and
re-derivable from the `SOURCE.md` files — **never upload them to eosvc.**
