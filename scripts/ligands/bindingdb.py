"""Does BindingDB add anything ChEMBL does not? -- a measurement, not a deliverable.

Stage 06, part 1b. BindingDB is the obvious second bioactivity source, and v1 built it as a
co-equal track (`legacy/scripts/06b_bindingdb_bioactivity.py`) reporting 93 potent Kp proteins
against ChEMBL's 175 -- **without ever measuring the overlap**, so that 93 may have been almost
entirely redundant. This script measures it before anything is promoted.

The approach is CLAUDE.md's eggNOG lesson: test the yield before paying for the database. There it
was "eggNOG places 73% of unannotated proteins in a group, but only 2-18% of those groups carry any
GO" -- a ceiling estimated from xref coverage that was far too optimistic. The same discipline
applies here, and the measurement is unusually cheap because BindingDB carries **the target chain
sequence inline** (column 39) and **the ligand's InChIKey** (column 4). So both halves of the join
that matters are already in the file: no accession mapping, and structural compound identity.

What is measured
----------------
Per species, at the same pChEMBL/pAff cutoff and the same identity bands as the ChEMBL track:

  * proteins with potent bacterial evidence from BindingDB alone
  * proteins that gain their FIRST potent ligand -- the only number that can justify the database
  * scaffolds genuinely new to the union, by InChIKey then by Murcko generic scaffold

Two deliberate restrictions
---------------------------
**Single-chain targets only.** A multichain BindingDB entry is a complex, and attributing its
ligand to each chain is the same over-claim the ChEMBL track refuses for PROTEIN COMPLEX. The
skipped fraction is reported rather than hidden.

**The bacterial genus map is a hard dependency, and fails loudly.** BindingDB records an organism
string and no taxonomy id, so the superkingdom has to come from somewhere. It is derived from the
ChEMBL target table this stage already wrote (`organism_class.l1`, 1,493 bacterial taxa). v1 did
the same thing but *warned and returned an empty set* when the file was missing, which silently
emptied its bacterial bucket. Here a missing or empty map exits non-zero, and the number of
BindingDB organisms that could not be classified is printed every run.

Run with the `gradi` env, after scripts/ligands/chembl.py:
    python scripts/ligands/bindingdb.py
    python scripts/ligands/bindingdb.py --limit-rows 200000   # smoke test
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "ligands"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
BDB_ROOT = REPO_ROOT / "data" / "raw" / "other" / "bindingdb"

DEFAULT_DIAMOND_BIN = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

# Column indices, 0-based, in BindingDB_All.tsv. The file is RAGGED -- measured 50 / 62 / 74 / 86
# columns in the first 200k rows, 50 base columns plus 12 per additional target chain -- so it must
# be streamed with the csv module and indexed positionally. pd.read_csv with a fixed header breaks.
C_SMILES, C_INCHIKEY, C_ORGANISM = 1, 3, 7
C_KI, C_IC50, C_KD, C_EC50 = 8, 9, 10, 11
C_LIGAND_CHEMBL = 32
C_N_CHAINS, C_SEQUENCE = 37, 38

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def ensure_diamond() -> str:
    override = os.environ.get("GRADI_DIAMOND_BIN")
    if shutil.which("diamond") is None:
        for d in [Path(override)] if override else [DEFAULT_DIAMOND_BIN]:
            if (d / "diamond").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("diamond")
    if exe is None:
        sys.exit("FAILED: no `diamond` on PATH (see scripts/ligands/chembl.py).")
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip()
    return f"{exe}  ({ver})"


def find_tsv() -> Path:
    hits = sorted(BDB_ROOT.glob("BindingDB_All.tsv"))
    if not hits:
        sys.exit(
            f"FAILED: no BindingDB_All.tsv under {BDB_ROOT.relative_to(REPO_ROOT)}.\n"
            "  unzip BindingDB_All_202504_tsv.zip in that directory (6.45 GB extracted)."
        )
    return hits[0]


def bacterial_genera() -> set[str]:
    """Genus names classified as Bacteria, from the ChEMBL target table this stage wrote.

    A hard dependency by design: v1 warned and returned an empty set here, which silently emptied
    its bacterial bucket rather than failing.
    """
    oc_path = EVIDENCE_DIR / "organism_class.tsv"
    if not oc_path.exists():
        sys.exit(
            f"FAILED: {oc_path.relative_to(REPO_ROOT)} is missing, and it is where the bacterial\n"
            "  genus map comes from. Run scripts/ligands/chembl.py first. (v1's bug was to\n"
            "  warn and carry on with an empty map, which silently emptied the bacterial bucket.)"
        )
    oc = pd.read_csv(oc_path, sep="\t", dtype=str, keep_default_na=False)
    # l3 is the genus-bearing rank in organism_class; l2/l1 back it up for the coarse entries.
    genera: set[str] = set()
    for col in ("l3", "l2"):
        genera |= {
            str(v).split()[0].lower()
            for v in oc.loc[oc.l1 == "Bacteria", col].dropna().unique()
            if str(v).strip()
        }
    try:
        tgt = L.load_targets()
        genera |= {
            str(o).split()[0].lower()
            for o in tgt.loc[tgt.superkingdom == "Bacteria", "organism"].dropna().unique()
            if str(o).strip()
        }
    except FileNotFoundError:
        pass
    if not genera:
        sys.exit("FAILED: the bacterial genus map is empty -- refusing to run with no restriction.")
    return genera


def parse_nm(raw: str) -> float | None:
    """A BindingDB affinity in nM. `>` values are censored non-binders and are discarded."""
    v = (raw or "").strip()
    if not v or v.startswith(">"):
        return None
    v = v.lstrip("<=~ ").strip()
    try:
        nm = float(v)
    except ValueError:
        return None
    return nm if nm > 0 else None


def stream(path: Path, limit_rows: int | None) -> tuple[pd.DataFrame, dict[str, str], dict]:
    """One pass: best pAff per (target sequence, ligand InChIKey), plus the sequences."""
    csv.field_size_limit(50_000_000)
    best: dict[tuple[str, str], float] = {}
    chembl_id: dict[str, str] = {}
    seqs: dict[str, str] = {}
    organism: dict[str, str] = {}
    stats = {"rows": 0, "multichain": 0, "no_affinity": 0, "kept": 0, "short_seq": 0}

    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        next(reader, None)
        for row in reader:
            stats["rows"] += 1
            if limit_rows and stats["rows"] > limit_rows:
                break
            if len(row) <= C_SEQUENCE:
                continue
            try:
                n_chains = int(row[C_N_CHAINS] or 1)
            except ValueError:
                n_chains = 1
            if n_chains != 1:
                stats["multichain"] += 1
                continue
            seq = (row[C_SEQUENCE] or "").strip().upper().replace("*", "")
            if len(seq) < 30 or not seq.isalpha():
                stats["short_seq"] += 1
                continue
            nm = None
            for col in (C_KI, C_KD, C_IC50, C_EC50):
                nm = parse_nm(row[col] if col < len(row) else "")
                if nm is not None:
                    break
            if nm is None:
                stats["no_affinity"] += 1
                continue
            # sha1, NOT the builtin hash(): Python randomises string hashing per process, so a
            # builtin-hash id would differ between the run that writes the cache and the run that
            # reads it, and the FASTA would stop matching the ligand table.
            key_seq = hashlib.sha1(seq.encode()).hexdigest()[:16]
            seqs.setdefault(key_seq, seq)
            organism.setdefault(key_seq, (row[C_ORGANISM] or "").strip())
            ik = (row[C_INCHIKEY] or "").strip()
            if not ik:
                continue
            cid = (row[C_LIGAND_CHEMBL] or "").strip() if C_LIGAND_CHEMBL < len(row) else ""
            if cid:
                chembl_id.setdefault(ik, cid)
            paff = 9.0 - math.log10(nm)
            k = (key_seq, ik)
            if paff > best.get(k, -99):
                best[k] = paff
            stats["kept"] += 1

    df = pd.DataFrame(
        [(s, i, p) for (s, i), p in best.items()], columns=["seq_id", "inchikey", "paff"]
    )
    df["organism"] = df.seq_id.map(organism)
    df["ligand_chembl_id"] = df.inchikey.map(chembl_id).fillna("")
    return df, seqs, stats


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--species", nargs="+", default=list(L.SPECIES), choices=list(L.SPECIES))
    ap.add_argument("--pchembl", type=float, default=L.PCHEMBL_HEADLINE)
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit-rows", type=int, help="only the first N rows of the TSV (smoke test)")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc)

    rule("=")
    say("STAGE 06 - part 1b: does BindingDB add anything ChEMBL does not?")
    rule("=")
    say("  purpose  : a MEASUREMENT. Nothing is promoted into the deliverable by this script.")
    say(f"  potency  : pAff >= {args.pchembl}, the same line as the ChEMBL track")
    say(f"  bands    : same as ChEMBL -- close >={L.CLOSE_PIDENT}, remote >={L.REMOTE_PIDENT}")
    say("  targets  : single-chain only; multichain is a complex, and per-chain attribution")
    say("             would be the over-claim the ChEMBL track refuses")
    if args.dry_run:
        say("\n  --dry-run: nothing read, nothing written.")
        return

    say()
    rule()
    say("TOOLS AND INPUTS")
    rule()
    say(f"  diamond         : {ensure_diamond()}")
    tsv = find_tsv()
    say(f"  bindingdb       : {tsv.relative_to(REPO_ROOT)}  ({tsv.stat().st_size / 1e9:.2f} GB)")
    genera = bacterial_genera()
    say(f"  genus map       : {len(genera):,} bacterial genera, from scratch/chembl_targets.tsv")

    say()
    rule()
    say("STREAM - one pass over a ragged 6.5 GB TSV")
    rule()
    cache = SCRATCH_DIR / ("smoke_bindingdb_ligands.tsv" if args.limit_rows else "bindingdb_ligands.tsv")
    fasta = SCRATCH_DIR / ("smoke_bindingdb_targets.faa" if args.limit_rows else "bindingdb_targets.faa")
    if not args.refresh and cache.exists() and fasta.exists():
        bdb = pd.read_csv(cache, sep="\t", dtype=str, keep_default_na=False)
        bdb["paff"] = pd.to_numeric(bdb.paff, errors="coerce")
        say(f"  cached          : {cache.name}")
    else:
        bdb, seqs, stats = stream(tsv, args.limit_rows)
        say(f"  rows read       : {stats['rows']:>10,}")
        say(f"  multichain      : {stats['multichain']:>10,} skipped "
            f"({100 * stats['multichain'] / max(stats['rows'], 1):.1f}%)")
        say(f"  no affinity     : {stats['no_affinity']:>10,} skipped "
            f"(no Ki/Kd/IC50/EC50, or a censored '>' value)")
        say(f"  measurements    : {stats['kept']:>10,} -> {len(bdb):,} (sequence, ligand) pairs")
        with open(fasta, "w") as fh:
            for sid, seq in seqs.items():
                if sid in set(bdb.seq_id):
                    fh.write(f">{sid}\n{seq}\n")
        bdb.to_csv(cache, sep="\t", index=False)
    say(f"  targets         : {bdb.seq_id.nunique():>10,} distinct sequences")
    say(f"  ligands         : {bdb.inchikey.nunique():>10,} distinct InChIKeys")

    known = (bdb.organism.fillna("").str.split().str[0].str.lower()).isin(genera)
    unclassified = bdb.loc[~known, "organism"].nunique()
    say(f"  bacterial       : {bdb.loc[known, 'seq_id'].nunique():>10,} sequences "
        f"({unclassified:,} organism strings could not be classified -- reported, not hidden)")

    say()
    rule()
    say("MAP - DIAMOND against the BindingDB target sequences")
    rule()
    tmp = Path(tempfile.mkdtemp(prefix="gradi06b_"))
    dbp = tmp / "bdb"
    subprocess.run(["diamond", "makedb", "--in", str(fasta), "-d", str(dbp), "--quiet"], check=True)
    hits = {}
    for sp in args.species:
        q = tmp / f"{sp}.faa"
        pro = P.load(sp)
        with open(q, "w") as fh:
            for ac, seq in zip(pro.uniprot_ac, pro.sequence):
                fh.write(f">{ac}\n{seq}\n")
        raw = tmp / f"{sp}.m8"
        subprocess.run(
            ["diamond", "blastp", "-q", str(q), "-d", str(dbp), "-o", str(raw),
             "--very-sensitive", "--id", "25", "--evalue", "1e-5", "--max-target-seqs", "100",
             "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore",
             "--quiet", "--threads", str(args.threads)],
            check=True,
        )
        h = pd.read_csv(raw, sep="\t",
                        names=["uniprot_ac", "seq_id", "pident", "qcov", "scov", "bitscore"])
        h = h.sort_values("bitscore", ascending=False).drop_duplicates(["uniprot_ac", "seq_id"])
        h["seq_id"] = h.seq_id.astype(str)
        hits[sp] = h
        say(f"  {sp:14s}: {len(h):>8,} hits  {h.uniprot_ac.nunique():>6,} proteins")

    say()
    rule()
    say("MARGINAL GAIN - the only number that can justify the database")
    rule()
    potent_bdb = bdb[(bdb.paff >= args.pchembl) & known]
    scaf = L.load_scaffolds()
    chembl_keys = set(scaf.standard_inchi_key.dropna())
    rows = []
    for sp in args.species:
        chembl = L.load_chembl(sp).set_index("uniprot_ac")
        h = hits[sp]
        h = h[(h.qcov >= L.MIN_QCOV) & (h.scov >= L.MIN_SCOV) & (h.pident >= L.REMOTE_PIDENT)]
        pairs = h[["uniprot_ac", "seq_id"]].merge(potent_bdb, on="seq_id", how="inner")
        per = pairs.groupby("uniprot_ac").inchikey.nunique()
        with_bdb = set(per[per > 0].index)
        with_chembl = set(chembl.index[chembl.remote_n_compounds > 0])
        gained = with_bdb - with_chembl
        novel_ik = set(pairs.inchikey) - chembl_keys
        rows.append({
            "species": sp,
            "chembl_proteins": len(with_chembl),
            "bindingdb_proteins": len(with_bdb),
            "overlap": len(with_bdb & with_chembl),
            "first_ligand_gained": len(gained),
            "bindingdb_ligands": pairs.inchikey.nunique(),
            "ligands_not_in_chembl": len(novel_ik),
            "pchembl": args.pchembl,
            "run_at": started.isoformat(timespec="seconds"),
        })
        say(f"  {sp:14s}: ChEMBL {len(with_chembl):>4,} proteins   BindingDB {len(with_bdb):>4,}"
            f"   overlap {len(with_bdb & with_chembl):>4,}"
            f"   NEW {len(gained):>4,}")
        say(f"  {'':14s}  ligands {pairs.inchikey.nunique():>6,}, of which "
            f"{len(novel_ik):>6,} carry an InChIKey ChEMBL does not have "
            f"({100 * len(novel_ik) / max(pairs.inchikey.nunique(), 1):.1f}%)")
    gain = pd.DataFrame(rows)
    pre = "smoke_" if args.limit_rows else ""
    gain.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}bindingdb_gain.tsv", sep="\t", index=False)

    say()
    rule()
    say("VERDICT")
    rule()
    total_new = int(gain.first_ligand_gained.sum())
    total_chembl = int(gain.chembl_proteins.sum())
    pct = 100 * total_new / max(total_chembl, 1)
    say(f"  {total_new} proteins across all species gain their first potent ligand from BindingDB,")
    say(f"  against {total_chembl} that ChEMBL already covers -- {pct:.1f}%.")
    say()
    if pct < 10:
        say("  READ THIS AS: BindingDB is largely redundant here. ChEMBL 37 already ingests")
        say("  BindingDB patent bioactivity (13,835 assays, 2.68M activities), which is most of")
        say("  the overlap. Not promoted into the deliverable.")
    else:
        say("  READ THIS AS: a material gain. Promoting it into the deliverable is defensible --")
        say("  see docs/ligands.md before doing so.")

    shutil.rmtree(tmp, ignore_errors=True)
    say()
    rule("=")
    say(f"stage 06 (bindingdb) complete in "
        f"{(datetime.now(timezone.utc) - started).seconds / 60:.1f} min.")
    rule("=")


if __name__ == "__main__":
    main()
