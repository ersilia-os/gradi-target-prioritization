"""Ligand precedent for ANY protein sequence: three counts, from cached ChEMBL extracts.

    (a) n_ligands_exact      ligands on an EXACT match -- UniProt accession, or identical sequence
    (b) n_ligands_bacteria   UNIQUE ligands reachable across BACTERIAL targets by identity
    (c) n_ligands_human      ligands on HUMAN orthologs

This is a QUERY TOOL over arbitrary input, which is what separates it from
`scripts/ligands/chembl.py`. That script answers "what does our proteome have" and needs the 30.5 GB
ChEMBL dump; this answers "what about this sequence" for any sequence, from three cached files
totalling 82 MB, and returns in about a second.

    chembl_targets.faa    8,469 sequences, headers = component_id   the DIAMOND subject database
    chembl_targets.tsv    9,347 targets, with organism + superkingdom
    chembl_ligands.tsv    2,591,526 rows: tid -> parent_molregno, pchembl, track

"UNIQUE" IS THE WHOLE POINT OF (b). `parent_molregno` is ChEMBL's `molecule_hierarchy` parent, so
salts are already collapsed; counting DISTINCT `parent_molregno` over the union of every homologous
target is what makes (b) a count of MOLECULES rather than of target-compound pairs. One compound
tested against three homologs counts ONCE. Summing per-target counts inflates it, and the wider the
identity band the worse that gets -- which is exactly the regime this tool is built for.

IDENTITY ALONE CANNOT TRANSFER A LIGAND COUNT, and the case that proves it is real. K. pneumoniae
`A0A0H3GWM6` is 99.2% identical to E. coli `P0ADG7`, which carries 32 compounds at pChEMBL 8.82 --
and it correctly gets nothing. The hit exists (component 1947 = CHEMBL3630, 99.2% id) but runs
`qcov 100.0 / scov 26.6`: the Kp entry is a 130-aa FRAGMENT against a 488-aa IMP dehydrogenase.
Both coverage floors are enforced here, imported from `src.ligandability`, so a fragment can never
inherit a whole enzyme's pharmacology.

(c) IS A LIABILITY, NOT A PRECEDENT. A ligand-bearing human ortholog says the fold is druggable AND
that hitting it may be dangerous. It is returned in its own columns and must never be summed into
(b). `clpP` is the case to remember: 106 human compounds against 61 bacterial.

WHAT THESE COUNTS DO NOT SAY. `pchembl` is 100% populated in the extract, so every count is of
POTENCY-MEASURABLE ligands -- `=` relations on IC50/EC50/Ki/Kd/Potency in nM. MIC and %-inhibition
are absent by construction, so **this tool does not say "has an antibiotic"**, and a ribosomal
protein looking empty here is a fact about the assay type, not about the biology.

    from src import precedents as PR
    PR.count({"my_protein": "MKTAYIAKQR..."})
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

from src import ligandability as L

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRATCH = REPO_ROOT / "data" / "processed" / "ligands" / "scratch"
TARGETS_FAA = SCRATCH / "chembl_targets.faa"
TARGETS_TSV = SCRATCH / "chembl_targets.tsv"
LIGANDS_TSV = SCRATCH / "chembl_ligands.tsv"
DEFAULT_DIAMOND_BIN = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

HUMAN = "Homo sapiens"
BACTERIA = "Bacteria"
DEFAULT_MIN_IDENTITY = L.REMOTE_PIDENT      # 40.0 -- the house floor for transfer

HIT_COLS = ["query", "component_id", "pident", "ppos", "length",
            "qlen", "slen", "qcov", "scov", "evalue", "bitscore"]

OUT_COLS = [
    "id", "n_ligands_exact", "n_ligands_bacteria", "n_ligands_human",
    "exact_route", "exact_target",
    "n_targets_bacteria", "best_pident_bacteria", "best_pchembl_bacteria",
    "n_targets_human", "best_pident_human", "best_pchembl_human",
    "n_ligands_bacteria_complex",
]

_CACHE: dict = {}


def diamond_bin() -> str:
    """DIAMOND from `gradi-ortho` (`GRADI_DIAMOND_BIN` overrides) -- it has no osx-arm64 build."""
    if shutil.which("diamond"):
        return "diamond"
    d = Path(os.environ.get("GRADI_DIAMOND_BIN", DEFAULT_DIAMOND_BIN))
    if (d / "diamond").exists():
        return str(d / "diamond")
    raise FileNotFoundError(
        f"no diamond on PATH or at {d}. It lives in `gradi-ortho` (osx-64, Rosetta); do NOT "
        "install it into `gradi`. Set GRADI_DIAMOND_BIN to override.")


def _tables() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """(targets, ligands, component_id -> sequence). Loaded once per process."""
    if "t" not in _CACHE:
        for p in (TARGETS_FAA, TARGETS_TSV, LIGANDS_TSV):
            if not p.exists():
                raise FileNotFoundError(
                    f"{p} missing -- it is written by scripts/ligands/chembl.py. The cached "
                    "extracts are all this tool needs; the 30.5 GB dump is not required.")
        t = pd.read_csv(TARGETS_TSV, sep="\t")
        t["component_id"] = t["component_id"].astype(str)
        lg = pd.read_csv(LIGANDS_TSV, sep="\t")
        seqs, cur, buf = {}, None, []
        for line in TARGETS_FAA.read_text().splitlines():
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if cur:
            seqs[cur] = "".join(buf)
        _CACHE["t"], _CACHE["l"], _CACHE["s"] = t, lg, seqs
    return _CACHE["t"], _CACHE["l"], _CACHE["s"]


def _diamond(seqs: dict[str, str], threads: int) -> pd.DataFrame:
    """All hits of every input against the ChEMBL target database, best HSP per pair.

    `--id 25` is a permissive PREFILTER, not the decision: the identity floor and both coverage
    floors are applied by the caller, so one search serves any `min_identity`.
    """
    with tempfile.TemporaryDirectory() as td:
        q, db, out = Path(td) / "q.faa", Path(td) / "db", Path(td) / "hits.tsv"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        subprocess.run([diamond_bin(), "makedb", "--in", str(TARGETS_FAA), "-d", str(db),
                        "--quiet"], check=True)
        subprocess.run([diamond_bin(), "blastp", "-q", str(q), "-d", str(db), "-o", str(out),
                        "--very-sensitive", "--id", "25", "--evalue", "1e-5",
                        "--max-target-seqs", "500", "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "ppos", "length", "qlen", "slen",
                        "qcovhsp", "scovhsp", "evalue", "bitscore",
                        "--quiet", "--threads", str(threads)], check=True)
        if not out.exists() or out.stat().st_size == 0:
            return pd.DataFrame(columns=HIT_COLS)
        d = pd.read_csv(out, sep="\t", names=HIT_COLS)
    d["component_id"] = d["component_id"].astype(str)
    # several HSPs per pair -> keep the best-scoring, as chembl.py does
    return d.sort_values("bitscore", ascending=False).drop_duplicates(["query", "component_id"])


def count(sequences: dict[str, str], accessions: dict[str, str] | None = None,
          min_identity: float = DEFAULT_MIN_IDENTITY, min_pchembl: float | None = None,
          threads: int = 4) -> pd.DataFrame:
    """Ligand precedent for each input sequence. One row per input, in input order.

    `sequences`  id -> protein sequence
    `accessions` id -> UniProt accession, optional. Only route (a) uses it; arbitrary input
                 usually has none, which is why the exact-sequence route exists.
    """
    # DIAMOND truncates `qseqid` at the first whitespace, so an id with a space would silently
    # return 0 rather than raising -- measured: `{"my prot A": folA}` gave 0, `{"lc": folA}` gave
    # the real count. Refuse rather than mislead.
    bad = [k for k in sequences if any(ch.isspace() for ch in str(k))]
    if bad:
        raise ValueError(
            f"sequence ids must not contain whitespace (DIAMOND truncates them): {bad[:5]}")
    t, lg, tseq = _tables()
    accessions = accessions or {}
    if min_pchembl is not None:
        lg = lg[lg["pchembl"] >= min_pchembl]

    # ONE COMPONENT MAPS TO MANY TIDS -- 9,347 rows over 8,469 distinct component_ids, 490 of them
    # with more than one tid and up to 15. `dict(zip(...))` keeps only the LAST and silently drops
    # the rest: measured, that understated 21 proteins by 3,324 ligands and made E. coli gyrA read
    # ZERO single-protein ligands, because component 166 carries tid 53 (SINGLE PROTEIN, 117
    # ligands) and tid 104721 (PROTEIN COMPLEX) and the complex row came last. That is the v1
    # "GyrA/GyrB look unliganded" error re-entering through a different door -- the very thing the
    # track split exists to prevent. `scripts/ligands/chembl.py:481` does it correctly with a
    # one-to-many merge; this now matches.
    comp2tids: dict[str, list] = {}
    for cid, tid in zip(t["component_id"], t["tid"]):
        comp2tids.setdefault(cid, []).append(tid)
    by_tid = lg.groupby("tid")
    single = {k: set(g.loc[g["track"] == "single", "parent_molregno"]) for k, g in by_tid}
    cplx = {k: set(g.loc[g["track"] == "complex", "parent_molregno"]) for k, g in by_tid}
    # best pChEMBL PER TRACK. Reading it across both tracks put a non-null
    # `best_pchembl_bacteria` next to `n_ligands_bacteria == 0` on 14 Kp rows -- a row that says
    # "no ligands, best potency 9.68" is incoherent on its face.
    best_single = lg[lg["track"] == "single"].groupby("tid")["pchembl"].max().to_dict()
    kingdom = dict(zip(t["component_id"], t["superkingdom"].astype(str)))
    organism = dict(zip(t["component_id"], t["organism"].astype(str)))
    # Identical sequences are shared by several components (57 of them, 2-4 each). Keeping one
    # would decide the winner by FASTA file order and silently understate 10 exact rows -- Kp
    # `bla` 24 against a true 241. Take them all; they are the same protein.
    seq2comps: dict[str, list] = {}
    for c, q in tseq.items():
        seq2comps.setdefault(q, []).append(c)
    acc2comps: dict[str, list] = {}
    for a, c in zip(t["accession"].astype(str), t["component_id"]):
        if a != "nan":
            acc2comps.setdefault(a, []).append(c)

    hits = _diamond(sequences, threads) if sequences else pd.DataFrame(columns=HIT_COLS)
    # ASSERT, AND FAIL LOUDLY. A whole proteome returning no hit at all means the search broke, not
    # that nothing binds -- and without this a --species run would write a table of zeros and exit
    # 0. Measured baseline: 6-8% of an anchor proteome hits something, so 0 of >=100 is impossible.
    if len(sequences) >= 100 and hits.empty:
        raise RuntimeError(
            f"DIAMOND returned no hit for any of {len(sequences):,} sequences. Expected ~6-8% to "
            "hit something; this is a broken search, not an empty result.")
    # BOTH coverage floors, imported not restated -- a 130-aa fragment must not inherit a 488-aa
    # enzyme's ligands. See the module docstring for the case that proves it.
    hits = hits[(hits["qcov"] >= L.MIN_QCOV) & (hits["scov"] >= L.MIN_SCOV)
                & (hits["pident"] >= min_identity)]
    by_query = dict(list(hits.groupby("query"))) if len(hits) else {}

    def ligands_of(comps) -> tuple[set, set]:
        """UNION of compounds over components -> DISTINCT molecules, counted once.

        Walks EVERY tid a component maps to, not just one. See the comment on `comp2tids`.
        """
        sg, cx = set(), set()
        for cid in comps:
            for tid in comp2tids.get(cid, ()):
                sg |= single.get(tid, set())
                cx |= cplx.get(tid, set())
        return sg, cx

    rows = []
    for qid, seq in sequences.items():
        # ---- (a) exact: accession first, then identical sequence
        route, comps = "none", []
        ac = accessions.get(qid)
        if ac and ac in acc2comps:
            route, comps = "accession", acc2comps[ac]
        elif seq in seq2comps:
            route, comps = "sequence", seq2comps[seq]
        ex_single, ex_cplx = ligands_of(comps)

        h = by_query.get(qid, pd.DataFrame(columns=HIT_COLS))
        bact = h[h["component_id"].map(kingdom).eq(BACTERIA)]
        hum = h[h["component_id"].map(organism).eq(HUMAN)]
        b_single, b_cplx = ligands_of(bact["component_id"])
        h_single, _ = ligands_of(hum["component_id"])

        def bp(sub):
            """Best pChEMBL over the SINGLE track only, matching what the count reports."""
            v = [best_single.get(tid) for c in sub["component_id"]
                 for tid in comp2tids.get(c, ())]
            v = [x for x in v if x is not None and pd.notna(x)]
            return max(v) if v else pd.NA

        rows.append({
            "id": qid,
            # SINGLE track, like (b) and (c). Mixing tracks here put 713 next to a
            # bacterial 0 on E. coli gyrA -- three columns that cannot be compared.
            "n_ligands_exact": len(ex_single),
            "n_ligands_bacteria": len(b_single),
            "n_ligands_human": len(h_single),
            "exact_route": route,
            "exact_target": ";".join(str(c) for c in comps) if comps else pd.NA,
            "n_targets_bacteria": len(bact),
            "best_pident_bacteria": round(bact["pident"].max(), 1) if len(bact) else pd.NA,
            "best_pchembl_bacteria": bp(bact),
            "n_targets_human": len(hum),
            "best_pident_human": round(hum["pident"].max(), 1) if len(hum) else pd.NA,
            "best_pchembl_human": bp(hum),
            "n_ligands_bacteria_complex": len(b_cplx),
        })
    out = pd.DataFrame(rows, columns=OUT_COLS)
    # TYPE THE COLUMNS EXPLICITLY. Building the frame from dicts containing `pd.NA` makes every
    # column that has a missing value OBJECT dtype, so `best_pchembl_bacteria` came out holding
    # the STRINGS '4.36' and '' -- `>= 6` then raises TypeError and a sort orders '9.02' above
    # '10.1'. Nullable Int64/Float64 keep "no measurement" distinguishable from zero while still
    # comparing and sorting numerically.
    for c in ("n_ligands_exact", "n_ligands_bacteria", "n_ligands_human",
              "n_targets_bacteria", "n_targets_human", "n_ligands_bacteria_complex"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    for c in ("best_pident_bacteria", "best_pchembl_bacteria",
              "best_pident_human", "best_pchembl_human"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Float64")
    out["exact_target"] = out["exact_target"].astype("string")
    return out


def for_accession(accession: str) -> pd.DataFrame:
    """Convenience: look the sequence up in our own proteomes, then `count` it."""
    from src import proteomes as P
    for sp in ("kpneumoniae", "ecoli", "saureus", "human"):
        d = P.load(sp)
        r = d[d["uniprot_ac"] == accession]
        if len(r):
            return count({accession: r["sequence"].iloc[0]}, {accession: accession})
    raise KeyError(f"{accession} is in none of the four reference proteomes; pass --sequence")
