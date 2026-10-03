"""How many experimental structures exist for this protein, and how much of it do they cover?

    data/processed/pockets/evidence/pdb_<species>.tsv        one row per protein, COMPLETE
    data/processed/pockets/evidence/pdb_chains_<species>.tsv LONG: one row per matched PDB chain
    data/processed/pockets/scratch/pdb_seqres/               unique-sequence FASTA, DIAMOND db, hits

    uniprot_ac · n_pdb_structures · pdb_n_chains · pdb_coverage · pdb_best_identity · pdb_ids

Ligand or not, it does not matter -- this is structural COVERAGE, the question `holo.py` does not
ask. Requested by the project owner on 2026-10-03.

* `n_pdb_structures` -- distinct PDB entries with at least one chain that IS this protein.
* `pdb_coverage` -- fraction (0-1) of this protein's residues covered by the union of those chains'
  alignments.

"IS THIS PROTEIN" = >= 95% identity, by sequence
--------------------------------------------------
`src.ligandability.DIRECT_PIDENT` -- the house "essentially this protein, possibly another strain"
rule. Matched by DIAMOND against every protein chain in the wwPDB's `pdb_seqres.txt`, NOT by
accession: v1's SIFTS accession route reached **30 of 5,728 Kp proteins (0.5%)** because HS11286 is
a dark TrEMBL proteome whose accessions the PDB never uses (`legacy/HISTORY.md:54`).

A chain counts only if **>= 50% of the PDB chain is aligned** (`src.ligandability.MIN_SCOV`): a
10-residue peptide co-crystallised with something else, or a fusion construct whose bulk is another
protein, is not a structure of this one. There is NO query-coverage floor -- a structure of one
domain is a structure, and `pdb_coverage` says how much it covers.

**Coverage is of the deposited SEQUENCE (SEQRES), not of resolved residues.** Disordered loops and
tails missing from the electron density still count as covered. Observed-residue coverage would
need the coordinates (SIFTS residue mappings) for every chain; SEQRES is the cheap, standard proxy
and overstates coverage where constructs carry unresolved termini.

`pdb_chains_<species>.tsv` (`uniprot_ac · pdb · chain · identity · chain_coverage`) is what lets
`holo.py` count
the ligands bound to this protein's OWN structures with no alignment of its own: these chains ARE
this protein, so their BioLiP ligand rows attach directly by `(pdb, chain)`.

**`chain_coverage` is how much of the PDB CHAIN our protein accounts for, and it is the fusion
warning.** At 1.0 the chain is this protein. Well below it, the chain is a construct in which this
protein is only the larger part -- an MBP fusion, say -- and a ligand bound to the *other* part is
still a BioLiP row on that chain. `MIN_SCOV = 50` lets those through by design (it exists to admit
domain constructs), so read this column before trusting a ligand count on a known crystallisation
chaperone. See `docs/pockets.md`.

Run with the `gradi` env; DIAMOND from `gradi-ortho` (`GRADI_DIAMOND_BIN` overrides). ~5 min.
  python scripts/pockets/pdb_coverage.py
  python scripts/pockets/pdb_coverage.py --species saureus -q
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402
from src.precedents import diamond_bin  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
SEQRES = REPO_ROOT / "data" / "source" / "wwpdb" / "pdb_seqres.txt.gz"
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
WORK = TASK_DIR / "scratch" / "pdb_seqres"
MIN_PIDENT = L.DIRECT_PIDENT
MIN_SCOV = L.MIN_SCOV
HIT_COLS = ["qseqid", "sseqid", "pident", "qstart", "qend", "sstart", "send", "qlen", "slen"]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def build_subjects() -> dict[str, list[str]]:
    """Unique protein sequences of pdb_seqres -> the chains carrying each. Writes the FASTA."""
    by_seq: dict[str, list[str]] = {}
    head, buf = None, []

    def flush():
        if head and "mol:protein" in head:
            by_seq.setdefault("".join(buf), []).append(head[1:].split()[0])

    with gzip.open(SEQRES, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                head, buf = line.strip(), []
            else:
                buf.append(line.strip())
    flush()
    ids = {}
    with open(WORK / "seqres_unique.faa", "w") as out:
        for i, (seq, chains) in enumerate(by_seq.items()):
            out.write(f">u{i}\n{seq}\n")
            ids[f"u{i}"] = chains
    return ids


def run_diamond(species: str, threads: int) -> pd.DataFrame:
    dmnd = WORK / "seqres_unique.dmnd"
    if not dmnd.exists():
        subprocess.run([diamond_bin(), "makedb", "--in", str(WORK / "seqres_unique.faa"),
                        "-d", str(dmnd), "--quiet"], check=True)
    q = WORK / f"query_{species}.faa"
    prot = P.load(species)
    with open(q, "w") as fh:
        for a, s in zip(prot["uniprot_ac"], prot["sequence"]):
            fh.write(f">{a}\n{s}\n")
    out = WORK / f"hits_{species}.tsv"
    if not out.exists() or out.stat().st_mtime < dmnd.stat().st_mtime:
        # High identity only, so --fast suffices; -k large because one protein can match
        # thousands of distinct construct sequences (lysozyme, DHFR, beta-lactamases).
        subprocess.run([diamond_bin(), "blastp", "-q", str(q), "-d", str(dmnd), "-o", str(out),
                        "--fast", "--id", str(MIN_PIDENT), "-k", "20000", "-e", "1e-5",
                        "--threads", str(threads), "--quiet", "--outfmt", "6", *HIT_COLS],
                       check=True)
    return pd.read_csv(out, sep="\t", header=None, names=HIT_COLS)


def score(species: str, hits: pd.DataFrame,
          chains_of: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(per-protein summary, long table of every matched PDB chain)."""
    prot = P.load(species)[["uniprot_ac", "sequence"]]
    hits = hits.assign(
        scov=100 * (hits["send"] - hits["sstart"] + 1) / hits["slen"])
    h = hits[(hits["pident"] >= MIN_PIDENT) & (hits["scov"] >= MIN_SCOV)]
    g = dict(tuple(h.groupby("qseqid")))
    rows, chain_rows = [], []
    for acc, seq in zip(prot["uniprot_ac"], prot["sequence"]):
        d = g.get(acc)
        if d is None:
            rows.append({"uniprot_ac": acc, "n_pdb_structures": 0, "pdb_n_chains": 0,
                         "pdb_coverage": 0.0, "pdb_best_identity": np.nan, "pdb_ids": ""})
            continue
        covered = np.zeros(len(seq), dtype=bool)
        for s, e in zip(d["qstart"], d["qend"]):
            covered[s - 1:e] = True
        chains = [c for u in d["sseqid"] for c in chains_of[u]]
        best = d.sort_values("pident", ascending=False).drop_duplicates("sseqid")
        for u, pid, sc in zip(best["sseqid"], best["pident"], best["scov"]):
            for c in chains_of[u]:
                pdb_id, _, ch = c.partition("_")
                chain_rows.append({"uniprot_ac": acc, "pdb": pdb_id, "chain": ch,
                                   "identity": float(pid), "chain_coverage": round(sc / 100, 4)})
        entries = sorted({c.split("_")[0] for c in chains})
        rows.append({"uniprot_ac": acc, "n_pdb_structures": len(entries),
                     "pdb_n_chains": len(chains), "pdb_coverage": round(float(covered.mean()), 4),
                     "pdb_best_identity": float(d["pident"].max()),
                     "pdb_ids": ";".join(entries)})
    long = pd.DataFrame(chain_rows,
                        columns=["uniprot_ac", "pdb", "chain", "identity", "chain_coverage"])
    return M.reindex(pd.DataFrame(rows), species), long.drop_duplicates()


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say("pockets/pdb_coverage.py -- experimental structures per protein, by sequence")
    say(f"  in : {SEQRES.relative_to(REPO_ROOT)} (every protein chain in the PDB)")
    say("  out: data/processed/pockets/evidence/pdb_<species>.tsv")
    say(f"  a chain IS this protein at >= {MIN_PIDENT:.0f}% identity with >= {MIN_SCOV:.0f}% of the "
        "PDB chain aligned; coverage = union of aligned residues (SEQRES)")
    say("=" * 92)
    WORK.mkdir(parents=True, exist_ok=True)
    chains_of = build_subjects()
    n_chains = sum(len(v) for v in chains_of.values())
    say(f"  {n_chains:,} protein chains, {len(chains_of):,} distinct sequences")

    summary = []
    for sp in args.species:
        hits = run_diamond(sp, args.threads)
        out, chains = score(sp, hits, chains_of)
        path = EVIDENCE_DIR / f"pdb_{sp}.tsv"
        out.to_csv(path, sep="\t", index=False)
        chains.to_csv(EVIDENCE_DIR / f"pdb_chains_{sp}.tsv", sep="\t", index=False)
        has = out["n_pdb_structures"] > 0
        summary.append({
            "species": sp, "n": len(out), "with_structure": int(has.sum()),
            "pct": 100 * has.mean(),
            "median_n_when>0": float(out.loc[has, "n_pdb_structures"].median()),
            "median_cov_when>0": float(out.loc[has, "pdb_coverage"].median()),
            "cov>=0.9": int((out["pdb_coverage"] >= 0.9).sum()),
        })
        say(f"  [{sp}] {len(hits):,} raw hits -> {path.relative_to(REPO_ROOT)} "
            f"({len(chains):,} matched chains)")

    say("\n" + "-" * 92)
    say("SUMMARY -- proteins with an experimental structure of their own (>= 95% identity)")
    say(pd.DataFrame(summary).to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()
