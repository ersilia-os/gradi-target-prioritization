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

# "Exact" means THIS PROTEIN, and a protein does not stop being itself between strains. The first
# version required an identical sequence (or an accession match), so one substitution in another
# K. pneumoniae isolate demoted the same enzyme to a "homolog" and its ligands left the exact
# count. `chembl.py` has had a `species` bucket for this reason since v2; `precedents.py` never
# got one. So the exact count is the UNION of three routes -- accession, identical sequence, and
# same species at >= EXACT_PIDENT -- all three of which name the same protein.
#
# 95 is `L.DIRECT_PIDENT`, documented there as "essentially this protein, possibly another strain",
# which is exactly the claim being made. It is NOT set by the transfer calibration:
# `scripts/ligands/transfer_calibration.py` measured P(potent | neighbour potent) as FLAT at
# 0.84-0.98 from 25% to 100% identity, so no identity threshold is an accuracy threshold here.
# This one is a statement about protein identity, not about how far evidence travels.
EXACT_PIDENT = L.DIRECT_PIDENT              # 95.0

HIT_COLS = ["query", "component_id", "pident", "ppos", "length",
            "qlen", "slen", "qcov", "scov", "evalue", "bitscore"]

OUT_COLS = [
    "id", "n_ligands_exact", "n_ligands_bacteria", "n_ligands_human",
    "exact_route", "exact_target",
    "n_targets_bacteria", "best_pident_bacteria", "best_pchembl_bacteria",
    "n_targets_human", "best_pident_human", "best_pchembl_human",
    "n_ligands_bacteria_complex",
    "n_compounds_assayed_bacteria", "n_compounds_potent_bacteria", "hit_rate_bacteria",
    "precedent_evidence",
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

        # SUPPLEMENT: bacterial targets the pChEMBL filter hid from the FASTA entirely.
        # `chembl.py` builds chembl_targets.faa from the tids that survived that filter, so a
        # target whose compounds were all assayed and none measurable has no sequence -- DIAMOND
        # can never hit it, and no protein could ever be told "somebody screened your homolog and
        # nothing came out". That made the `screened_clean` tier structurally unreachable, i.e. a
        # category that would have shipped permanently empty while looking meaningful. 326 targets.
        # They carry NO rows in chembl_ligands.tsv by construction, so adding them cannot change
        # any ligand count -- only `n_targets_bacteria`, `best_pident_bacteria`, and the tier.
        for extra, faa in ((SCRATCH / "chembl_effort_targets.tsv",
                            SCRATCH / "chembl_effort_targets.faa"),):
            if extra.exists() and faa.exists():
                e = pd.read_csv(extra, sep="\t")
                e["component_id"] = e["component_id"].astype(str)
                t = pd.concat([t, e[[c for c in e.columns if c in t.columns]]], ignore_index=True)
                cur, buf = None, []
                for line in faa.read_text().splitlines():
                    if line.startswith(">"):
                        if cur:
                            seqs.setdefault(cur, "".join(buf))
                        cur, buf = line[1:].split()[0], []
                    else:
                        buf.append(line.strip())
                if cur:
                    seqs.setdefault(cur, "".join(buf))

        _CACHE["t"], _CACHE["l"], _CACHE["s"] = t, lg, seqs
    return _CACHE["t"], _CACHE["l"], _CACHE["s"]


def binomial(organism: str) -> str:
    """The two-word binomial -- `chembl.py`'s own same-species rule, reused not restated.

    ChEMBL files strains under their own names: `Escherichia coli`, `Escherichia coli K-12` and
    `Escherichia coli (strain K12)` are one species. Genus alone would sweep in
    `Klebsiella aerogenes`; the full string would split a species into strains, which is the whole
    thing the species route exists to stop doing.
    """
    if not isinstance(organism, str) or not organism.strip():
        return ""
    return " ".join(organism.split()[:2]).lower()


def _effort() -> tuple[dict, bool]:
    """tid -> set of ASSAYED parent compounds, from `scripts/ligands/effort.py`.

    Sets, not counts, because the denominator must be UNIONED over the homology pool exactly as the
    numerator is: a compound assayed against three homologous targets is one compound. Summing
    per-target counts would inflate the denominator precisely where the pool is widest, which is
    where the hit rate matters most.

    Optional by design -- it needs the 30.5 GB dump restored, while everything else here runs off
    82 MB of cached extracts. Absent, the effort columns are NA and the tier degrades to what the
    old two-way split could say. NA is honest; a zero would claim nobody ever screened the protein.
    """
    if "e" not in _CACHE:
        path = SCRATCH / "chembl_assayed.tsv"
        if not path.exists():
            _CACHE["e"] = ({}, False)
        else:
            a = pd.read_csv(path, sep="\t")
            _CACHE["e"] = (a.groupby("tid")["parent_molregno"].apply(set).to_dict(), True)
    return _CACHE["e"]


def _diamond(seqs: dict[str, str], threads: int) -> pd.DataFrame:
    """All hits of every input against the ChEMBL target database, best HSP per pair.

    `--id 25` is a permissive PREFILTER, not the decision: the identity floor and both coverage
    floors are applied by the caller, so one search serves any `min_identity`.

    **The subject database is built from `_tables()`'s sequences, not from `TARGETS_FAA` directly.**
    It used to read the file, which silently excluded the 326 effort-only targets that `_tables()`
    concatenates -- so they were reachable by exact match and invisible to DIAMOND, and the
    `screened_clean` tier stayed empty while looking implemented. One source of subject sequences,
    or the two drift apart with nothing raising.
    """
    _, _, subject = _tables()
    with tempfile.TemporaryDirectory() as td:
        q, db, out = Path(td) / "q.faa", Path(td) / "db", Path(td) / "hits.tsv"
        subj = Path(td) / "subject.faa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        subj.write_text("".join(f">{k}\n{v}\n" for k, v in subject.items()))
        subprocess.run([diamond_bin(), "makedb", "--in", str(subj), "-d", str(db),
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
          organisms: dict[str, str] | str | None = None,
          min_identity: float = DEFAULT_MIN_IDENTITY, min_pchembl: float | None = None,
          exact_pident: float = EXACT_PIDENT, threads: int = 4) -> pd.DataFrame:
    """Ligand precedent for each input sequence. One row per input, in input order.

    `sequences`  id -> protein sequence
    `accessions` id -> UniProt accession, optional. Only the accession route uses it; arbitrary
                 input usually has none, which is why the other two routes exist.
    `organisms`  id -> organism name, or ONE name for every input. Enables the species route.
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
    # One organism for every input, or one per input. A bare string is the common case: a whole
    # proteome is one species.
    org_of = ({k: organisms for k in sequences} if isinstance(organisms, str)
              else dict(organisms or {}))
    lg_all = lg                              # unfiltered -- the potency tier needs the real values
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

    assayed_by_tid, have_effort = _effort()
    potent_by_tid = (lg_all[(lg_all["track"] == "single")
                            & (lg_all["pchembl"] >= L.PCHEMBL_HEADLINE)]
                     .groupby("tid")["parent_molregno"].apply(set).to_dict())

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
        h = by_query.get(qid, pd.DataFrame(columns=HIT_COLS))
        bact = h[h["component_id"].map(kingdom).eq(BACTERIA)]
        hum = h[h["component_id"].map(organism).eq(HUMAN)]

        # ---- exact: the UNION of three routes, each of which names THIS protein.
        # Not an if/elif chain: an accession match and a same-species 99% match can be different
        # ChEMBL components carrying different compounds, and both are this protein's evidence.
        routes, comps = [], set()
        ac = accessions.get(qid)
        if ac and ac in acc2comps:
            routes.append("accession")
            comps.update(acc2comps[ac])
        if seq in seq2comps:
            routes.append("sequence")
            comps.update(seq2comps[seq])
        want = binomial(org_of.get(qid, ""))
        if want and len(bact):
            same = bact[(bact["pident"] >= exact_pident)
                        & bact["component_id"].map(lambda c: binomial(organism.get(c, "")) == want)]
            if len(same):
                routes.append(f"species_{exact_pident:g}")
                comps.update(same["component_id"])
        route = ";".join(routes) if routes else "none"
        comps = sorted(comps)
        ex_single, ex_cplx = ligands_of(comps)

        b_single, b_cplx = ligands_of(bact["component_id"])
        h_single, _ = ligands_of(hum["component_id"])

        # ---- effort: distinct compounds ASSAYED across the same bacterial pool, and the tier.
        b_assayed, b_potent = set(), set()
        for cid in bact["component_id"]:
            for tid in comp2tids.get(cid, ()):
                b_assayed |= assayed_by_tid.get(tid, set())
                b_potent |= potent_by_tid.get(tid, set())
        n_assayed = len(b_assayed) if have_effort else pd.NA
        n_potent = len(b_potent)
        hit_rate = (len(b_potent & b_assayed) / len(b_assayed)) if (have_effort and b_assayed) else pd.NA

        # A zero is an ANSWER, not a gap -- the studiedness axis's `no_hit` vs `below_floor` rule
        # applied here. `screened_clean` (somebody tried, nothing measurable came out) and
        # `never_screened` (nobody opened it) are opposite evidence and were the same zero across
        # ~96% of every proteome until now.
        if len(b_single):
            tier = "liganded"
        elif not len(bact):
            tier = "no_homolog"
        elif not have_effort:
            tier = "unknown_effort"
        elif b_assayed:
            tier = "screened_clean"
        else:
            tier = "never_screened"

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
            "n_compounds_assayed_bacteria": n_assayed,
            "n_compounds_potent_bacteria": n_potent,
            "hit_rate_bacteria": hit_rate,
            "precedent_evidence": tier,
        })
    out = pd.DataFrame(rows, columns=OUT_COLS)
    # TYPE THE COLUMNS EXPLICITLY. Building the frame from dicts containing `pd.NA` makes every
    # column that has a missing value OBJECT dtype, so `best_pchembl_bacteria` came out holding
    # the STRINGS '4.36' and '' -- `>= 6` then raises TypeError and a sort orders '9.02' above
    # '10.1'. Nullable Int64/Float64 keep "no measurement" distinguishable from zero while still
    # comparing and sorting numerically.
    for c in ("n_ligands_exact", "n_ligands_bacteria", "n_ligands_human",
              "n_targets_bacteria", "n_targets_human", "n_ligands_bacteria_complex",
              "n_compounds_assayed_bacteria", "n_compounds_potent_bacteria"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Int64")
    for c in ("best_pident_bacteria", "best_pchembl_bacteria",
              "best_pident_human", "best_pchembl_human", "hit_rate_bacteria"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("Float64")
    for c in ("exact_target", "exact_route", "precedent_evidence"):
        out[c] = out[c].astype("string")
    return out


def for_accession(accession: str) -> pd.DataFrame:
    """Convenience: look the sequence up in our own proteomes, then `count` it."""
    from src import proteomes as P
    for sp in ("kpneumoniae", "ecoli", "saureus", "human"):
        d = P.load(sp)
        r = d[d["uniprot_ac"] == accession]
        if len(r):
            return count({accession: r["sequence"].iloc[0]}, {accession: accession},
                         organisms=L.SPECIES_ORGANISM.get(sp))
    raise KeyError(f"{accession} is in none of the four reference proteomes; pass --sequence")
