"""Essentiality labels -- the training corpus for ProteomeLM-Ess (stage 07, part 1).

ProteomeLM-Ess, as published, is trained on OGEE v3. **OGEE v3 is unreachable**: the host
`v3.ogee.info` aborts the TLS handshake under every configuration tested -- LibreSSL, OpenSSL 3.5
and real Chrome (`ERR_SSL_PROTOCOL_ERROR`), with SNI, without SNI, pinned to TLS 1.2, ALPN disabled,
and by direct IP -- so it is server-side and no client-side workaround exists. The OGEE v3 paper
names no FTP, mirror or deposit, and there is no GitHub copy (only human subsets). v1 recorded the
same failure in May (legacy/HISTORY.md trap 29).

So the corpus is **DEG**, which is OGEE's own prokaryotic upstream, and which is alive. This script
builds the label table from it and nothing else: fetch, parse, flag, report. Training happens in
`07_essentiality_proteomelm.py`.

What DEG gives, and why it is usable:
  * 66 bacterial datasets over 42 species, with RefSeq replicon accessions per dataset
  * 26,619 essential genes, each with **its protein sequence** (`DEG10.aa`, 1:1 by DEG gene id)
  * the sequences are the point -- we join to our proteomes BY SEQUENCE, never on DEG's identifier
    columns, which are patchy: gene symbol 63.1%, GI 54.5%, COG 48.4%, UniProt AC 45.7%

What DEG does NOT give: the non-essential list. 38 of 66 datasets report a non-essential *count* but
ship no genes. Negatives therefore come from each organism's complete proteome minus the positives,
which is only valid for a genome-wide screen -- hence the two exclusions below.

Three exclusions, FLAGGED not dropped (every row is emitted; `retained` says what training may use):
  * not genome-wide -- antisense RNA, MATT, insertion-duplication, transposon-hybridisation. These
    did not test every gene, so "absent from the list" does not mean non-essential.
  * condition-specific -- tobramycin, murine pneumonia, bile, cholesterol, kanamycin. Conditional
    fitness, not core essentiality; these are why P. aeruginosa PAO1 swings 117 -> 336 -> 551
    essential genes across its three datasets. Detected from the condition string AND, crucially,
    from DEG's own per-gene notes: a dataset whose genes are tagged "Condition-Specific Essential
    Gene" or "Recovered from blood/cerebrospinal fluid/meninges" is a conditional screen whatever
    its condition string says. That data-driven rule is what catches DEG1058 (S. suis), whose
    condition reads "Columbia blood base agar" but whose paper is titled "Identification of
    CONDITIONALLY essential genes for Streptococcus suis infection in pigs".
  * placeholder sequences -- **326 of DEG's 26,619 "sequences" are the literal English string
    "Not available now."** The join to DEG10.aa succeeds for them, so a join-coverage check passes
    at 100% while the payload is junk; only validating the amino-acid alphabet catches it. E. coli
    O157:H7 (DEG1056) is 12.0% placeholder and E. coli is a held-out anchor, so this one matters.
    A further 52 carry an INTERNAL "*" stop -- frameshifted or mis-annotated translations,
    mostly N. gonorrhoeae MS11 -- and are rejected too; a TRAILING "*"/"$" is merely a stop
    marker and is stripped. Short proteins are NOT filtered: rpmJ (50S ribosomal L36) is a
    genuine 37-aa essential gene.

Consolidation to one label set per species is DEFERRED on purpose. 66 datasets over 42 species means
S. aureus appears 7x, P. aeruginosa 4x, Salmonella 4x, E. coli 4x, M. tuberculosis 3x. The paper
consolidated OGEE's 127 studies into 87 taxids, so we will too, but the union-vs-intersection rule
gets chosen against measured counts rather than guessed. This script keeps every dataset separate.

Output
  data/source/deg/                 the three bulk files as fetched, + SOURCE.md
  data/processed/essentiality/evidence/
      deg_datasets.tsv    66 rows -- one per DEG dataset, every source column + the flags
      deg_genes.tsv       26,619 rows -- one per essential gene, every source column + sequence
      deg_manifest.tsv    counts, bytes, checksums, fetch time

Run
  python scripts/essentiality/labels.py
  python scripts/essentiality/labels.py --refresh      # re-download the bulk files
  python scripts/essentiality/labels.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import sys
import urllib.request
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

RAW_DIR = REPO_ROOT / "data" / "source" / "deg"
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
DEG_BASE = "http://tubic.org/deg/public/download"
DEG_FILES = {
    "datasets": "deg_bacteria.csv.zip",      # 66 dataset records
    "genes": "deg_annotation_p.csv.zip",     # 26,619 essential genes
    "sequences": "DEG10.aa.gz",              # 26,619 protein sequences
}

# deg_bacteria.csv is ';'-quoted with no header. Column order verified by inspection.
DATASET_COLS = [
    "organism", "reference", "pmid", "replicon", "n_essential", "n_nonessential",
    "reference_2", "pmid_2", "reference_3", "pmid_3",
    "method", "condition", "deg_dataset_id", "deg_date",
]

# deg_annotation_p.csv, likewise headerless.
GENE_COLS = [
    "deg_dataset_id", "deg_gene_id", "gene_symbol", "gi", "cog",
    "function_class", "description", "organism", "replicon", "condition",
    "locus_tag_raw", "go_terms", "uniprot_ac", "deg_gene_note",
]

# Sequence sanity. DEG serves this exact string where it has no sequence.
PLACEHOLDER_SEQUENCE = "Not available now."
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

# DEG's own per-gene notes that mark a screen as conditional rather than core-essentiality.
CONDITIONAL_NOTE_PREFIXES = ("Condition-Specific", "Recovered from")

# Methods that did NOT test every gene, so "absent from the essential list" != non-essential.
NOT_GENOME_WIDE = {
    "Antisense RNA",
    "MATT",
    "Insertion-duplication",
    "Insertion-duplication and allelic replacement",
    "Transposon mutagenesis followed by hybridization",
}

# A condition string is core-essentiality only if it is plain rich medium. Anything naming a drug,
# a host, or a required nutrient is a conditional-fitness screen.
CONDITION_DISQUALIFIERS = (
    "tobramycin", "kanamycin", "bile", "cholesterol", "model", "murine",
    "required for", "optimal growth",
)

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def sha16(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()[:16]


def fetch(name: str, refresh: bool) -> bytes:
    """Download one DEG bulk file, verifying the payload against Content-Length."""
    fname = DEG_FILES[name]
    dest = RAW_DIR / fname
    if dest.exists() and not refresh:
        data = dest.read_bytes()
        say(f"  {fname:<26} cached   {len(data):>10,} bytes  sha={sha16(data)}")
        return data

    url = f"{DEG_BASE}/{fname}"
    req = urllib.request.Request(url, headers={"User-Agent": "gradi/2.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        declared = resp.headers.get("Content-Length")
        data = resp.read()
    # An HTTP 200 is not evidence of data -- assert the payload.
    if declared is not None and int(declared) != len(data):
        sys.exit(f"FATAL {fname}: Content-Length {int(declared):,} != {len(data):,} received")
    if len(data) < 1024:
        sys.exit(f"FATAL {fname}: {len(data)} bytes is too small to be real")
    dest.write_bytes(data)
    say(f"  {fname:<26} fetched  {len(data):>10,} bytes  sha={sha16(data)}"
        f"  (Content-Length {'ok' if declared else 'absent'})")
    return data


def read_csv_from_zip(data: bytes, expect_suffix: str) -> list[list[str]]:
    """DEG ships ';'-delimited, quoted, headerless CSV inside a zip (plus __MACOSX cruft)."""
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        names = [n for n in z.namelist()
                 if n.endswith(expect_suffix) and not n.startswith("__MACOSX")]
        if len(names) != 1:
            sys.exit(f"FATAL expected exactly one {expect_suffix} in the zip, found {names}")
        text = z.read(names[0]).decode("utf8", errors="replace")
    return list(csv.reader(io.StringIO(text), delimiter=";"))


def read_fasta(data: bytes) -> dict[str, str]:
    text = gzip.decompress(data).decode("utf8", errors="replace")
    seqs: dict[str, str] = {}
    name, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def frame(rows: list[list[str]], cols: list[str], label: str) -> pd.DataFrame:
    widths = Counter(len(r) for r in rows)
    if list(widths) != [len(cols)]:
        sys.exit(f"FATAL {label}: expected {len(cols)} columns throughout, saw {dict(widths)}")
    df = pd.DataFrame(rows, columns=cols)
    return df.apply(lambda s: s.str.strip() if s.dtype == object else s)


def is_condition_ok(condition: str) -> bool:
    c = condition.lower()
    if not any(m in c for m in ("rich medium", "todd-hewitt", "columbia blood")):
        return False
    return not any(d in c for d in CONDITION_DISQUALIFIERS)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="re-download the DEG bulk files")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STAGE 07 part 1 - essentiality labels, from DEG")
    rule("=")
    say(f"  corpus        DEG bacteria ({DEG_BASE})")
    say("  why not OGEE  v3.ogee.info aborts the TLS handshake for every client incl. Chrome")
    say("  positives     essential genes + their protein sequences (join downstream BY SEQUENCE)")
    say("  negatives     NOT shipped by DEG -- complete proteome minus positives, genome-wide only")
    say("  exclusions    flagged, never dropped: non-genome-wide methods, condition-specific screens")
    say("  consolidation deferred -- every dataset kept separate, one row per (dataset, gene)")
    say(f"  raw           {RAW_DIR.relative_to(REPO_ROOT)}")
    say(f"  out           {EVIDENCE_DIR.relative_to(REPO_ROOT)}/deg_datasets.tsv, deg_genes.tsv")
    rule("=")
    if args.dry_run:
        say("dry run -- nothing fetched, nothing written.")
        return

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    rule()
    say("FETCH")
    rule()
    blobs = {k: fetch(k, args.refresh) for k in DEG_FILES}

    rule()
    say("PARSE")
    rule()
    ds = frame(read_csv_from_zip(blobs["datasets"], "deg_bacteria.csv"), DATASET_COLS, "datasets")
    gn = frame(read_csv_from_zip(blobs["genes"], "deg_annotation_p.csv"), GENE_COLS, "genes")
    seqs = read_fasta(blobs["sequences"])
    say(f"  datasets   {len(ds):>7,} rows x {ds.shape[1]} cols")
    say(f"  genes      {len(gn):>7,} rows x {gn.shape[1]} cols")
    say(f"  sequences  {len(seqs):>7,} records")

    # The sequence join is 1:1 by DEG gene id, and everything downstream depends on it.
    gn["sequence"] = gn["deg_gene_id"].map(seqs)
    missing = int(gn["sequence"].isna().sum())
    if missing:
        sys.exit(f"FATAL {missing:,} of {len(gn):,} genes have no sequence in DEG10.aa")
    say(f"  sequence join: {len(gn):,}/{len(gn):,} (100.0%) by DEG gene id")

    # A successful join is NOT a valid payload. Validate the alphabet.
    gn["sequence"] = gn["sequence"].str.strip().str.rstrip("*$")
    gn["sequence_valid"] = gn["sequence"].map(
        lambda s: bool(s) and s != PLACEHOLDER_SEQUENCE and not (set(s.upper()) - AA_ALPHABET))
    gn["seq_len"] = gn["sequence"].str.len().where(gn["sequence_valid"])
    n_bad = int((~gn["sequence_valid"]).sum())
    n_ph = int((gn["sequence"] == PLACEHOLDER_SEQUENCE).sum())
    say(f"  sequence VALID: {int(gn.sequence_valid.sum()):,}/{len(gn):,} "
        f"({100 * gn.sequence_valid.mean():.1f}%)  median length {gn.seq_len.median():.0f} aa")
    say(f"  rejected {n_bad:,}: {n_ph:,} are the placeholder {PLACEHOLDER_SEQUENCE!r}; "
        f"{n_bad - n_ph:,} carry an INTERNAL '*' stop, i.e. a frameshifted or mis-annotated "
        "translation (a trailing stop is stripped, not rejected)")
    if n_bad:
        worst = gn[~gn.sequence_valid].groupby("deg_dataset_id").size().sort_values(ascending=False)
        tot = gn.groupby("deg_dataset_id").size()
        say("    worst-affected datasets:")
        for dsid, n in worst.head(5).items():
            say(f"      {dsid}  {n:>4}/{tot[dsid]:<5} ({100 * n / tot[dsid]:.1f}%)")

    # locus_tag_raw packs two identifiers; split them out rather than leaving a compound string.
    lt = gn["locus_tag_raw"].str.extract(r"locus_tag:([^;]+)")[0]
    gn["locus_tag"] = lt
    gn["gi_2"] = gn["locus_tag_raw"].str.extract(r"gi:(\d+)")[0]
    say(f"  locus_tag parsed from locus_tag_raw: {int(lt.notna().sum()):,} "
        f"({100 * lt.notna().mean():.1f}%)")

    # DEG's per-gene notes betray conditional screens the condition string does not.
    note = gn["deg_gene_note"].fillna("")
    gn["note_is_conditional"] = note.str.startswith(CONDITIONAL_NOTE_PREFIXES)
    cond_ds = set(gn.loc[gn.note_is_conditional, "deg_dataset_id"].unique())
    say(f"  datasets whose GENE NOTES mark them conditional: "
        f"{sorted(cond_ds) if cond_ds else 'none'}")

    rule()
    say("FLAG")
    rule()
    ds["n_essential_int"] = pd.to_numeric(ds["n_essential"], errors="coerce").astype("Int64")
    ds["n_nonessential_int"] = pd.to_numeric(ds["n_nonessential"], errors="coerce").astype("Int64")
    ds["has_nonessential_list"] = ds["n_nonessential_int"].fillna(0) > 0
    ds["is_genome_wide"] = ~ds["method"].isin(NOT_GENOME_WIDE)
    ds["is_core_condition"] = (ds["condition"].map(is_condition_ok)
                               & ~ds["deg_dataset_id"].isin(cond_ds))
    ds["retained"] = ds["is_genome_wide"] & ds["is_core_condition"]
    ds["exclusion_reason"] = [
        "" if r else ("not_genome_wide" if not g else "condition_specific")
        for r, g in zip(ds["retained"], ds["is_genome_wide"])
    ]
    ds["species"] = ds["organism"].str.split().str[:2].str.join(" ")

    say(f"  genome-wide method     {int(ds.is_genome_wide.sum()):>3}/{len(ds)}")
    say(f"  core (rich-medium)     {int(ds.is_core_condition.sum()):>3}/{len(ds)}")
    say(f"  RETAINED               {int(ds.retained.sum()):>3}/{len(ds)}"
        f"   over {ds.loc[ds.retained, 'species'].nunique()} species")
    say(f"  ships a non-ess list   {int(ds.has_nonessential_list.sum()):>3}/{len(ds)}"
        "   (the rest need a complete proteome for negatives)")
    say("")
    say("  excluded datasets, with the reason:")
    for _, r in ds[~ds.retained].iterrows():
        say(f"    {r.deg_dataset_id}  {r.exclusion_reason:<18} {r.method[:34]:<34} "
            f"cond={r.condition[:26]!r:<28} {r.organism[:34]}")

    gn = gn.merge(
        ds[["deg_dataset_id", "species", "method", "retained", "exclusion_reason",
            "has_nonessential_list", "replicon", "pmid"]].rename(
                columns={"replicon": "dataset_replicon", "pmid": "dataset_pmid"}),
        on="deg_dataset_id", how="left", validate="many_to_one")
    if gn["retained"].isna().any():
        sys.exit("FATAL some genes reference a DEG dataset id absent from the index")

    # DEG's index states n_essential; the annotation file ships rows. They do not always agree.
    actual = gn.groupby("deg_dataset_id").size().rename("n_rows")
    chk = ds.set_index("deg_dataset_id")[["organism", "n_essential_int"]].join(actual)
    chk["delta"] = chk["n_rows"] - chk["n_essential_int"]
    off = chk[chk["delta"] != 0]
    say("")
    if len(off):
        say(f"  index vs shipped rows: {len(off)} of {len(chk)} datasets DISAGREE "
            "(DEG's own inconsistency; the shipped rows are what we use)")
        for dsid, r in off.iterrows():
            say(f"    {dsid}  index={r.n_essential_int}  shipped={r.n_rows}  "
                f"delta={r.delta:+d}  {str(r.organism)[:40]}")
    else:
        say(f"  index vs shipped rows: all {len(chk)} datasets agree")

    # The same protein listed twice inside one screen would double-weight it in training.
    dup = gn[gn.sequence_valid].groupby(["deg_dataset_id", "sequence"]).size()
    say(f"  duplicate (dataset, sequence) pairs: {int((dup > 1).sum())}"
        " -- deduplicate before training, keep both rows here for provenance")

    rule()
    say("RETAINED CORPUS, BY SPECIES")
    rule()
    gn["trainable"] = gn["retained"].astype(bool) & gn["sequence_valid"]
    say("  trainable = retained dataset AND a valid protein sequence")
    say(f"  {'species':<40} {'sets':>4} {'essential':>10} {'unique seq':>11}")
    keep = gn[gn.trainable]
    for sp, g in sorted(keep.groupby("species"), key=lambda kv: -len(kv[1])):
        say(f"  {sp[:40]:<40} {g.deg_dataset_id.nunique():>4} {len(g):>10,} "
            f"{g.sequence.nunique():>11,}")
    rule()
    say(f"  {'TOTAL':<40} {keep.deg_dataset_id.nunique():>4} {len(keep):>10,} "
        f"{keep.sequence.nunique():>11,}")
    say(f"  {'excluded / invalid':<40} {(~ds.retained).sum():>4} {len(gn) - len(keep):>10,}")

    rule()
    say("OUTPUTS")
    rule()
    paths = {}
    for label, df in (("deg_datasets", ds), ("deg_genes", gn)):
        p = EVIDENCE_DIR / f"{label}.tsv"
        df.to_csv(p, sep="\t", index=False)
        paths[label] = p
        say(f"  {p.relative_to(REPO_ROOT)}  ({len(df):,} rows x {df.shape[1]} cols, "
            f"{p.stat().st_size:,} bytes)")

    man = pd.DataFrame([{
        "file": f, "bytes": len(b), "sha256_16": sha16(b),
        "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    } for f, b in ((DEG_FILES[k], blobs[k]) for k in DEG_FILES)])
    man.to_csv(EVIDENCE_DIR / "deg_manifest.tsv", sep="\t", index=False)
    say(f"  {(EVIDENCE_DIR / 'deg_manifest.tsv').relative_to(REPO_ROOT)}")

    # Integrity re-read: the file on disk is what the rest of the stage will consume.
    back = pd.read_csv(paths["deg_genes"], sep="\t", low_memory=False)
    if len(back) != len(gn) or back["sequence"].isna().any():
        sys.exit("FATAL deg_genes.tsv did not round-trip")
    if int(back["trainable"].sum()) != int(gn["trainable"].sum()):
        sys.exit("FATAL trainable flag did not round-trip")
    bad_back = back.loc[back["trainable"], "sequence"].map(
        lambda x: bool(set(str(x).upper()) - AA_ALPHABET)).sum()
    if bad_back:
        sys.exit(f"FATAL {bad_back} trainable rows carry a non-amino-acid sequence")
    say(f"  re-read check: {len(back):,} rows, {int(back.trainable.sum()):,} trainable, "
        "every trainable sequence alphabet-clean")

    rule("=")
    say("stage 07 part 1 complete.")
    rule("=")


if __name__ == "__main__":
    main()
