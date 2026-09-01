# Localization — access routes and identifier bridges (docs §5.1)

Everything `scripts/09*` fetches, how to fetch it again, and the traps found on the way. Companion
to `docs/localization_log.md` (decisions and results). Same role as
`docs/degradability_downloads.md` plays for §3.

All caches live under `data/` and are eosvc-tracked, not Git. Every fetcher is idempotent: if the
cache file exists and is non-empty it is reused, so re-running a script is cheap.

---

## A. UniProt subcellular location — 09a

```
https://rest.uniprot.org/uniprotkb/stream
  ?compressed=false&format=tsv
  &query=proteome:UP000007841          # HS11286  (UP000000625 = E. coli K-12)
  &fields=accession,cc_subcellular_location,ft_signal,ft_transmem,ft_lipid
```

No cache file — it is a fast stream, re-fetched each run and written straight to
`output/results/<org>/<prefix>_loc_uniprot.csv`.

Parsing notes:

- **Evidence split.** The CC text carries inline ECO codes. `ECO:0000269` (experimental evidence,
  manual assertion) and `ECO:0007744` (combinatorial, manual assertion) mark a call as
  `experimental`; anything else annotated is `curated`. This is the split that makes 09f worth
  doing — see the counts in the log.
- **`ft_lipid`** looks like `LIPID 21; /note="S-diacylglycerol cysteine"; /evidence=…`. The number
  is the 1-based position of the lipidated Cys **in the precursor** (signal peptide included), which
  is exactly what `lipoprotein_sorting()` expects. `"diacylglycerol cysteine"` in the note is the
  Sec/SPII lipoprotein signature.
- **Compartment order in the CC text does not imply priority.** See the peripheral-membrane trap in
  the log — `classify_uniprot()` handles it explicitly.

## B. PSORTdb — precomputed PSORTb 3.0 — 09b

Cache: `data/raw/<org>/localization/psortdb/assembly_<id>.tab`

**Resolving the assembly id.** Not discoverable from the browse URLs (`?query=` is ignored — the
pager returns all 638 pages regardless). POST the full genome name to the search endpoint and read
the `assembly_id` out of the resulting HTML:

```bash
curl -sL -X POST "https://db.psort.org/search/results" \
  --data-urlencode "organism=Klebsiella pneumoniae subsp. pneumoniae HS11286" \
  --data "dataset=c" | grep -oE 'assembly_id=[0-9]+' | sort -u
```

Already resolved, and hard-coded in `09b_psortdb.py`:

| Organism | PSORTdb `assembly_id` | RefSeq assembly |
|---|---|---|
| K. pneumoniae HS11286 | `446671` | GCF_000240185.1 |
| E. coli K-12 MG1655 | `449203` | GCF_000005845.2 |

**The download.** A per-genome tab file, ~3 MB:

```bash
curl -sL -X POST "http://db.psort.org/search/results/download?assembly=446671&id=" \
  -H "Content-Type: application/x-www-form-urlencoded" --data "id=" -o assembly_446671.tab
```

Traps, all hit during the build:

- **The POST needs a body.** With no body the server returns `411 Invalid Request`. `--data "id="`
  is enough.
- **The search endpoint is slow** — a genome-name query took over 2 minutes at one point. Cache the
  resolved id (they are in the table above) and never look it up in a loop.
- **Do not use the bulk download.** `Computed-Gram_negative-PSORTdb-3.00.tab.tar.gz` is 1.2 GB
  gzipped and **18.7 GB uncompressed as a single TSV** with every Gram-negative genome concatenated
  and no genome column in the header — you would have to stream-scan the whole thing filtering on a
  RefSeq accession set. The per-genome route above returns the same rows in seconds.
- **Columns.** `SeqID` is `ref|YP_005221101.1`. The call is `Final_Localization`
  (`Cytoplasmic` / `CytoplasmicMembrane` / `Periplasmic` / `OuterMembrane` / `Extracellular` /
  `Cellwall` / `Unknown`), scored 0–10 in `Final_Score`. Roughly 29% of Kp rows come back `Unknown`,
  which is the behaviour DeepLocPro was chosen to avoid.

## C. UniProt ID mapping — RefSeq → UniProt — 09b

Cache: `data/raw/<org>/localization/psortdb/assembly_<id>_refseq2uniprot.tsv`

Three-step async job; `from=RefSeq_Protein`, `to=UniProtKB`, chunked at 5,000 ids:

```
POST https://rest.uniprot.org/idmapping/run                      -> {"jobId": …}
GET  https://rest.uniprot.org/idmapping/status/<jobId>           -> poll while jobStatus RUNNING|NEW
GET  https://rest.uniprot.org/idmapping/uniprotkb/results/stream/<jobId>?format=tsv
```

Result columns are `From` / `Entry`. Mapping quality is excellent — **every** mapped Kp accession
landed inside the reference proteome (5,728/5,728, 100%); E. coli reaches 93.4%.

`09b` asserts on this: below `MIN_COVERAGE = 0.50` it exits with the remedy in the message rather
than quietly writing a thin table.

> **A fallback that was deliberately not written.** The plan called for a DIAMOND best-hit fallback
> (the house rule for the dark Kp proteome). At 100%/93.4% accession coverage it would be dead code
> that has never executed, so it is a loud assert plus a documented remedy instead. If RefSeq ever
> drifts, DIAMOND the anchor proteome against the genome's RefSeq FASTA with
> `L.run_diamond_blastp` and transfer by best hit.

## D. STEPdb 2.0 — E. coli experimental localization — 09f

Cache: `data/raw/ecoli/localization/stepdb/stepdb_k12_basic_proteome.csv`

```
http://stepdb.eu/static/downs/E%20coli%20K-12%20strain%20MG1655%20basic%20proteome.csv
```

**Pick the right file.** `stepdb.eu/info/downloads` offers a dozen CSV/XLS exports and only this one
carries per-protein subcellular classes keyed by UniProt accession. In particular
`E coli K-12 MG1655 strain - Protein location and pseudogenes.csv` sounds right and is **not** —
it is a pseudogene table, `@`-delimited, with padded junk columns.

Columns used: `Accession (UniProt)` (already UniProt — no mapping needed),
`STEPdb Sub-cellular Location (Full Name)`, `Annotation References`.

Traps:

- **Semicolon-delimited, with semicolons inside free-text cells.** A minority of rows shift columns
  as a result. `classify_stepdb()` matches ordered substrings and returns `unknown` for anything it
  does not recognise, so a shifted row is dropped rather than mis-assigned.
- **Compound classes** — `"Nucleoid, Integral Inner Membrane"` — so first-match-wins ordering
  matters; the more specific patterns come first in `STEPDB_PATTERNS`.
- **There is no evidence flag.** A non-empty `Annotation References` is used as the proxy for bench
  support (2,114/3,897 rows).

## E. Model weights — 09c, 09d

Downloaded automatically on first run, into the `gradi-loc` env, not into `data/`:

| Model | Size | Pulled by | Landing |
|---|---|---|---|
| ESM-2 650M (`esm2_t33_650M_UR50D`) | ~2.5 GB | `fair-esm` | `~/.cache/torch/hub/checkpoints/` |
| ProtT5-XL-U50 encoder | ~2.25 GB | `transformers` | `<site-packages>/tmbed/models/t5/` |

DeepLocPro's 20 classifier checkpoints (1.3 MB each) ship inside the pip package. Per-protein and
per-shard prediction caches live under `data/processed/<org>/localization/{deeplocpro,tmbed}/` and
are fully regenerable — they are the bulk of the ~53 MB this axis writes.

## F. SignalP 6.0 — 09e — manual, licensed

Not fetchable programmatically. Accept the academic terms at
<https://services.healthtech.dtu.dk/services/SignalP-6.0/> (current package version 6.0i, `fast` and
`slow_sequential` variants; `fast` is a distilled model and is what 09e expects), install the
package, then point the script at the binary:

```bash
export SIGNALP6_BIN=/path/to/signalp6
python scripts/09e_signalp.py --organism kpneumoniae   # then re-run 09g, 09h, 08a
```

Without it, 09e falls back to the calibrated lipobox and records `signalp_source=lipobox` so the two
are always distinguishable downstream. The direct `sw_request` CGI URL pattern seen in search results
returns `Error 24` — go through the service page.

## G. Identifier bridges in this axis

| From | To | Route |
|---|---|---|
| RefSeq `YP_*` / `NP_*` (PSORTdb) | UniProt | UniProt ID-mapping job (§C) |
| UniProt (STEPdb) | UniProt | none needed — STEPdb is already UniProt-keyed |
| E. coli UniProt | Kp UniProt | 03a OrthoFinder orthogroups, donor consensus (see log) |
| FASTA header | UniProt | `L.acc_from_header` — handles `sp\|P12345\|…`, `tr\|…\|…` and bare accessions |

## H. Stale artifacts to clean

The pre-rework track left two things behind that nothing reads any more:

- `data/raw/<org>/localization/<prefix>_localization.tsv` — the old merged TSV, superseded by
  `output/results/<org>/<prefix>_localization.csv`.
- `data/raw/<org>/localization/psortb/` — empty; the Docker PSORTb run never produced output.

Safe to delete when convenient; both are eosvc-tracked, so removing them is a data-tree change.
