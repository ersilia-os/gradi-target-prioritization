"""Essentiality negatives -- complete proteomes for the DEG corpus (stage 07, part 2).

DEG ships POSITIVES ONLY: 26,619 essential genes, no non-essential list (38 of 66 datasets report a
non-essential *count* and ship no genes). A classifier needs both classes, and ProteomeLM needs the
WHOLE proteome regardless -- it contextualises every protein against every other, so unlabelled
proteins still have to be in the forward pass even though they contribute nothing to the loss.

So: fetch each retained dataset's complete proteome from the RefSeq/EMBL replicon accessions DEG
records, then label it. Essential = matched a DEG positive. Non-essential = everything else, which is
only legitimate because part 1 already excluded the non-genome-wide screens.

The join is BY SEQUENCE, the house rule, and here it is exact-match first: DEG's sequences and the
replicon translations should both descend from the same RefSeq annotation. The measured exact-match
rate per dataset is the control -- a low rate means DEG's vintage annotation has drifted from the
current one, and that is reported per dataset rather than papered over.

Two edge cases, both flagged and neither silently dropped:
  * DEG1053 B. cenocepacia K56-2 records `LAUA00000000`, a WGS *master* accession, and efetch
    returns zero CDS for it.
  * DEG1037 S. pyogenes MGAS5448 records no replicon at all (the field is "-").

The replicon field itself is not clean: separators seen are ",", ";", ", " and bare spaces, so
accessions are extracted by regex rather than split on any one delimiter.

Output
  data/source/ncbi/deg_proteomes/<deg_dataset_id>.faa    as fetched, one file per dataset
  data/processed/essentiality/evidence/
      proteome_join.tsv       one row per dataset -- CDS count, matched, join rate, edge-case flag
      labeled_proteins.tsv    one row per protein per dataset: ids, sequence, essential 0/1

Run
  python scripts/essentiality/deg_proteomes.py
  python scripts/essentiality/deg_proteomes.py --limit 3      # smoke test, writes scratch/smoke_*
  python scripts/essentiality/deg_proteomes.py --refresh      # re-download every proteome
  python scripts/essentiality/deg_proteomes.py --dry-run
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

RAW_DIR = REPO_ROOT / "data" / "source" / "ncbi" / "deg_proteomes"
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "essentiality" / "evidence"
SCRATCH_DIR = REPO_ROOT / "data" / "processed" / "essentiality" / "scratch"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# NCBI allows 3 requests/second without an API key; stay under it.
REQUEST_PAUSE_S = 0.5
# An accession: NC_000964(.3), NZ_CP008957, AM747720, CP029332, LAUA00000000, AE014133.
# Validate whole tokens rather than scanning: a scanning regex anchored on \b silently drops
# NZ_-prefixed RefSeq accessions, because "_" is a word character so the boundary never fires.
# That cost E. coli O157:H7 -- the largest single dataset, 1,071 positives -- on the first run.
ACCESSION_RE = re.compile(r"^(?:[A-Z]{2}_)?[A-Z]{0,4}\d{5,}(?:\.\d+)?$")
ACCESSION_SPLIT_RE = re.compile(r"[,;\s]+")
# Below this, DEG's annotation has drifted from the current replicon and the labels are unreliable.
MIN_JOIN_RATE = 0.80
# Exact sequence match recovers only 59-95% of positives: DEG's vintage annotation and the current
# replicon differ by point substitutions (same length, different residues -- measured, not assumed;
# it is not a start-codon or prefix/suffix artefact, and locus_tag is 0% on the older datasets).
# So fall back to DIAMOND at the house threshold, the same >=95% identity stage 04 and 06 use.
DIAMOND_MIN_PIDENT = 95.0
DIAMOND_MIN_COVERAGE = 80.0
DEFAULT_DIAMOND_DIR = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def ensure_diamond() -> str:
    """Locate the DIAMOND binary, or exit with the fix. Same borrow pattern as stages 04/05/06."""
    override = os.environ.get("GRADI_DIAMOND_BIN")
    if shutil.which("diamond") is None:
        for d in ([Path(override)] if override else [DEFAULT_DIAMOND_DIR]):
            if (d / "diamond").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("diamond")
    if exe is None:
        sys.exit(
            "diamond not found.\n"
            f"  looked on PATH and in {override or DEFAULT_DIAMOND_DIR}\n"
            "  it lives in the `gradi-ortho` env (osx-64, no arm64 build); point GRADI_DIAMOND_BIN\n"
            "  at a directory containing it. Exact sequence matching recovers as little as 59% of\n"
            "  DEG's positives against the current replicon annotation, so without DIAMOND this\n"
            "  stage would label hundreds of genuinely essential genes as non-essential."
        )
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return f"{exe}  ({(out.stdout or out.stderr).strip().splitlines()[0]})"


def diamond_recover(unmatched: list[str], prot: pd.DataFrame) -> dict[str, tuple[str, float]]:
    """Map DEG positives that failed exact match onto CDS of the same proteome.

    Returns {deg_sequence: (cds_sequence, pident)} for best hits clearing the identity and
    coverage floors. Within one organism this is near-lossless; it is still a threshold, so the
    identity distribution is reported.
    """
    if not unmatched or prot.empty:
        return {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        q, db = tmp / "q.faa", tmp / "t.faa"
        q.write_text("".join(f">q{i}\n{s}\n" for i, s in enumerate(unmatched)))
        db.write_text("".join(f">t{i}\n{s}\n" for i, s in enumerate(prot["sequence"])))
        subprocess.run(["diamond", "makedb", "--in", str(db), "-d", str(tmp / "t"), "--quiet"],
                       check=True, capture_output=True)
        res = tmp / "hits.tsv"
        subprocess.run(
            ["diamond", "blastp", "-q", str(q), "-d", str(tmp / "t"), "-o", str(res),
             "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore",
             "--max-target-seqs", "1", "--very-sensitive", "--quiet"],
            check=True, capture_output=True)
        out: dict[str, tuple[str, float]] = {}
        targets = list(prot["sequence"])
        best: dict[str, tuple[float, str, float]] = {}
        for line in res.read_text().splitlines():
            qid, sid, pid, qcov, bits = line.split("\t")
            pid, qcov, bits = float(pid), float(qcov), float(bits)
            if pid < DIAMOND_MIN_PIDENT or qcov < DIAMOND_MIN_COVERAGE:
                continue
            if qid not in best or bits > best[qid][0]:
                best[qid] = (bits, sid, pid)
        for qid, (_, sid, pid) in best.items():
            out[unmatched[int(qid[1:])]] = (targets[int(sid[1:])], pid)
    return out


def accessions(field: str) -> list[str]:
    """DEG's replicon field uses ',', ';', ', ' and bare spaces. Split, then validate each token."""
    if not field or str(field) == "-":
        return []
    return [t for t in ACCESSION_SPLIT_RE.split(str(field).upper()) if ACCESSION_RE.match(t)]


def efetch_cds_aa(accs: list[str], dest: Path, refresh: bool) -> str:
    if dest.exists() and not refresh and dest.stat().st_size > 0:
        return dest.read_text()
    url = f"{EFETCH}?db=nuccore&id={','.join(accs)}&rettype=fasta_cds_aa&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "gradi/2.0 (essentiality stage 07)"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                text = resp.read().decode("utf8", errors="replace")
            break
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == 2:
                say(f"      efetch failed after 3 attempts: {exc}")
                return ""
            time.sleep(2 ** attempt)
    time.sleep(REQUEST_PAUSE_S)
    if text.count(">"):
        dest.write_text(text)
    return text


def parse_cds_fasta(text: str) -> pd.DataFrame:
    """efetch fasta_cds_aa headers look like
    >lcl|NC_000964.3_prot_NP_387882.1_1 [gene=dnaA] [locus_tag=BSU_00010] [protein=...]"""
    rows, header, buf = [], None, []

    def flush() -> None:
        if header is None:
            return
        seq = "".join(buf).strip().rstrip("*")
        tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", header))
        rows.append({
            "fasta_id": header.split()[0].lstrip(">"),
            "locus_tag": tags.get("locus_tag", ""),
            "protein_id": tags.get("protein_id", ""),
            "gene": tags.get("gene", ""),
            "protein_name": tags.get("protein", ""),
            "pseudo": "pseudo" in tags,
            "sequence": seq,
        })

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            header, buf = line, []
        elif header is not None:
            buf.append(line.strip())
    flush()
    return pd.DataFrame(rows)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, help="only the first N datasets (smoke test)")
    ap.add_argument("--refresh", action="store_true", help="re-download every proteome")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    pre = "smoke_" if args.limit else ""

    ds = pd.read_csv(EVIDENCE_DIR / "deg_datasets.tsv", sep="\t")
    gn = pd.read_csv(EVIDENCE_DIR / "deg_genes.tsv", sep="\t", low_memory=False)
    keep = ds[ds.retained].copy()
    if args.limit:
        keep = keep.head(args.limit)

    rule("=")
    say("STAGE 07 part 2 - complete proteomes, so the corpus has a negative class")
    rule("=")
    say(f"  datasets      {len(keep)} retained of {len(ds)}   ({keep.species.nunique()} species)")
    say(f"  positives     {int(gn.trainable.sum()):,} trainable DEG essential genes")
    say("  negatives     complete proteome minus positives (genome-wide screens only)")
    say("  source        NCBI efetch db=nuccore rettype=fasta_cds_aa, by DEG's replicon accessions")
    say("  join          exact sequence match, then DIAMOND >=95% id for the drift")
    say(f"  floor         a dataset below {MIN_JOIN_RATE:.0%} join rate is flagged, not trusted")
    say(f"  raw           {RAW_DIR.relative_to(REPO_ROOT)}/<deg_dataset_id>.faa")
    say(f"  out           {EVIDENCE_DIR.relative_to(REPO_ROOT)}/{pre}proteome_join.tsv, "
        f"{pre}labeled_proteins.tsv")
    rule("=")
    if args.dry_run:
        say("dry run -- nothing fetched, nothing written.")
        return

    say(f"  diamond       {ensure_diamond()}")
    rule("=")
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    pos_by_ds = {k: set(v) for k, v in
                 gn[gn.trainable].groupby("deg_dataset_id")["sequence"].apply(list).items()}

    rule()
    say("FETCH + JOIN")
    rule()
    say(f"  {'dataset':<9} {'CDS':>6} {'DEG+':>6} {'exact':>6} {'+dmnd':>6} {'rate':>6}  organism")
    joins, labeled = [], []
    for i, (_, r) in enumerate(keep.iterrows(), 1):
        accs = accessions(r.replicon)
        pos = pos_by_ds.get(r.deg_dataset_id, set())
        note = ""
        if not accs:
            note = "no_replicon_accession"
            text = ""
        else:
            text = efetch_cds_aa(accs, RAW_DIR / f"{r.deg_dataset_id}.faa", args.refresh)
            if not text.count(">"):
                note = "efetch_returned_no_cds"

        prot = parse_cds_fasta(text) if text else pd.DataFrame(
            columns=["fasta_id", "locus_tag", "protein_id", "gene", "protein_name",
                     "pseudo", "sequence"])
        n_exact = n_diamond = 0
        pidents: list[float] = []
        if len(prot):
            pset = set(prot["sequence"])
            exact = pos & pset
            n_exact = len(exact)
            # DEG's vintage annotation drifts from the current replicon; recover the rest by homology.
            recovered: dict[str, tuple[str, float]] = {}
            if pos - pset:
                recovered = diamond_recover(sorted(pos - pset), prot)
                pidents = [pid for _, pid in recovered.values()]
            hit_seqs = exact | {cds for cds, _ in recovered.values()}
            n_diamond = len({cds for cds, _ in recovered.values()} - exact)
            prot["deg_dataset_id"] = r.deg_dataset_id
            prot["species"] = r.species
            prot["essential"] = prot["sequence"].isin(hit_seqs).astype(int)
            prot["match_method"] = [
                "exact" if q in exact else ("diamond" if q in hit_seqs else "")
                for q in prot["sequence"]]
            labeled.append(prot)
        matched = n_exact + n_diamond
        rate = matched / len(pos) if pos else float("nan")
        joins.append({
            "deg_dataset_id": r.deg_dataset_id, "organism": r.organism, "species": r.species,
            "replicon_field": r.replicon, "accessions": ";".join(accs),
            "n_cds": len(prot), "n_deg_positives": len(pos),
            "n_exact": n_exact, "n_diamond": n_diamond, "n_matched": matched,
            "diamond_median_pident": round(float(pd.Series(pidents).median()), 2) if pidents else None,
            "join_rate": round(rate, 4) if pos else None,
            "n_essential_labeled": int(prot["essential"].sum()) if len(prot) else 0,
            "note": note,
        })
        flag = ""
        if note:
            flag = f"  <- {note}"
        elif pos and rate < MIN_JOIN_RATE:
            flag = f"  <- LOW JOIN ({rate:.1%})"
        rate_s = f"{rate:>5.1%}" if pos else f"{'n/a':>6}"
        say(f"  [{i:>2}/{len(keep)}] {r.deg_dataset_id:<9} {len(prot):>6,} {len(pos):>6,} "
            f"{n_exact:>6,} {n_diamond:>6,} {rate_s}  {r.organism[:36]}{flag}")

    jdf = pd.DataFrame(joins)
    lab = (pd.concat(labeled, ignore_index=True) if labeled
           else pd.DataFrame(columns=["deg_dataset_id", "sequence", "essential"]))

    rule()
    say("SUMMARY")
    rule()
    ok = jdf[(jdf.note == "") & (jdf.join_rate.notna())]
    say(f"  datasets fetched            {int((jdf.n_cds > 0).sum())}/{len(jdf)}")
    say(f"  edge cases                  {int((jdf.note != '').sum())}"
        f"   {sorted(set(jdf.loc[jdf.note != '', 'note']))}")
    if len(ok):
        say(f"  join rate  median {ok.join_rate.median():.1%}   min {ok.join_rate.min():.1%}   "
            f"max {ok.join_rate.max():.1%}")
        say(f"  recovered by exact match {int(ok.n_exact.sum()):,}  "
            f"(+{int(ok.n_diamond.sum()):,} by DIAMOND, "
            f"median identity {ok.diamond_median_pident.median():.1f}%)")
        low = ok[ok.join_rate < MIN_JOIN_RATE]
        say(f"  below the {MIN_JOIN_RATE:.0%} floor          {len(low)}")
        for _, r in low.sort_values("join_rate").iterrows():
            say(f"      {r.deg_dataset_id}  {r.join_rate:.1%}  "
                f"{r.n_matched:,}/{r.n_deg_positives:,}  {r.organism[:44]}")
    say("")
    say(f"  proteins labeled            {len(lab):,}")
    if len(lab):
        say(f"  essential                   {int(lab.essential.sum()):,} "
            f"({100 * lab.essential.mean():.2f}%)")
        say(f"  non-essential               {int((lab.essential == 0).sum()):,}")
        say(f"  distinct species            {lab.species.nunique()}")
        say(f"  pseudogenes (context only)  {int(lab.pseudo.sum()):,}")

    rule()
    say("OUTPUTS")
    rule()
    for label, df in ((f"{pre}proteome_join", jdf), (f"{pre}labeled_proteins", lab)):
        p = (SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{label}.tsv"
        df.to_csv(p, sep="\t", index=False)
        say(f"  {p.relative_to(REPO_ROOT)}  ({len(df):,} rows x {df.shape[1]} cols, "
            f"{p.stat().st_size:,} bytes)")

    # Every DEG positive that matched must be labeled essential exactly once per dataset.
    if len(lab):
        per = lab.groupby("deg_dataset_id")["essential"].sum()
        exp = jdf.set_index("deg_dataset_id")["n_matched"]
        bad = [d for d in per.index if per[d] < exp[d]]
        if bad:
            sys.exit(f"FATAL {len(bad)} datasets label fewer essentials than matched: {bad[:5]}")
        say(f"  integrity: every matched positive is labeled, across {len(per)} datasets")

    rule("=")
    say("stage 07 part 2 complete.")
    rule("=")


if __name__ == "__main__":
    main()
