"""Ligandability from measured bioactivity -- can a small molecule bind this protein?

Stage 06, part 1. Reads the ChEMBL 37 SQLite dump, keeps every activity with a pChEMBL value on a
single-protein or protein-complex target, maps our three bacterial proteomes onto those targets by
DIAMOND sequence search, and counts the NON-REDUNDANT ligands reachable at three levels of
evolutionary distance. Writes one dense row per protein.

Why this axis exists
--------------------
Stages 00-05 answer identity, function, localization, degradability and orthology. None of them
answers whether a protein can be drugged -- which, for a BacPROTAC programme, decides whether a
degrader is even conceivable. This is the strongest available evidence: somebody measured a real
compound against a real protein.

Why homology transfer is the whole game
---------------------------------------
The anchor proteome is dark TrEMBL and its accessions are essentially absent from ChEMBL. Measured
here: ChEMBL 37 holds **21** K. pneumoniae, **225** E. coli and **80** S. aureus single-protein
targets in total, against proteomes of 5,728 / 4,403 / 2,889. Exact lookup is not a strategy. So
the axis is built on CLAUDE.md's house rule -- map by sequence, never by accession -- and every
count is "ligands measured against something resembling this protein", with the resemblance
reported next to it.

Same-species matching is by organism NAME, not tax_id
------------------------------------------------------
ChEMBL files strains under their own taxids. Measured: `tax_id = 562` finds 65 E. coli
single-protein targets, while `organism LIKE 'Escherichia coli%'` finds 225 -- K-12 lives at 83333.
A taxid match would have undercounted the best-represented of our three species threefold. The
prefix is the two-word binomial so that `Klebsiella aerogenes` is not swept in by a genus match.

One confidence gate cannot serve both target types
---------------------------------------------------
ChEMBL's `confidence_score` is target-type specific: 9/8 are direct/homologous SINGLE PROTEIN,
7/6 are direct/homologous PROTEIN COMPLEX subunits, 5/4 are family-level. Two measured
consequences, both load-bearing:

  * `>= 8` removes **0 of 3,271,336** single-protein rows -- it is entirely subsumed by requiring
    a pChEMBL value at all. It is kept as an explicit guard, but it must not be described as
    quality work it did not do.
  * `>= 8` returns **exactly zero** protein complexes. Applying one gate to both tracks would have
    shipped empty `complex_*` columns that read as a real biological zero.

That matters more than it sounds: **DNA gyrase is a PROTEIN COMPLEX in ChEMBL** (E. coli 713
compounds, S. aureus 491), as is topoisomerase IV (200) and Mtb ClpP1P2. GyrA/GyrB are in the
consortium's own interest panel, and a SINGLE-PROTEIN-only rule would have made them look
unliganded.

Non-redundant means scaffolds, not molregnos
---------------------------------------------
v1 counted `COUNT(DISTINCT molregno)`, so a compound and its hydrochloride salt counted twice. Here
compounds are collapsed to their `molecule_hierarchy` parent, and the headline diversity number is
distinct **Bemis-Murcko generic scaffolds**: 104 potent compounds may be one optimised analog series
or 40 independent chemical ideas, and only the second is strong evidence that a protein is
ligandable. An acyclic compound has no Murcko scaffold and joins one shared empty-string bucket --
never dropped, never counted as zero.

Counts are unions over a pool, not one best hit
------------------------------------------------
v1 picked one target per bucket by `max(n_potent, best_pchembl)` and reported its counts. That
under-counts a protein resembling several liganded targets and biases toward whichever target was
screened hardest. Here each bucket's count is the union of distinct parent compounds over every
target in the pool. The closest target is still reported, as provenance.

Output
------
    data/processed/ligands/evidence/chembl_<species>.tsv
        uniprot_ac, direct_hit, best_target, best_pident, best_organism, n_compounds_tested,
        {species,close,remote}_{n_compounds,n_scaffolds,best_pchembl},
        human_{n_compounds,n_scaffolds,best_pchembl,best_pident},
        allorg_{n_compounds,n_scaffolds}, complex_n_compounds, complex_best_target

    scratch/chembl_targets.tsv       component_id x tid bridge
    scratch/chembl_ligands.tsv       tid x parent_molregno x best pchembl
    evidence/chembl_hits.tsv          the raw DIAMOND join, unfiltered by band
    scratch/chembl_targets.faa       the DIAMOND target database
    evidence/scaffolds.tsv            parent_molregno -> smiles -> generic scaffold
    evidence/cutoff_sensitivity.tsv   every bucket count at pChEMBL >= 5 / 6 / 7
    evidence/control.tsv              v1's numbers against v2's
    evidence/manifest.tsv

Run with the `gradi` env (DIAMOND is borrowed from `gradi-ortho`; `GRADI_DIAMOND_BIN` overrides):
    python scripts/ligands/chembl.py                  # full run
    python scripts/ligands/chembl.py --limit 200      # smoke test -> scratch/smoke_*
    python scripts/ligands/chembl.py --refresh        # discard the SQL + DIAMOND caches
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "ligands"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
CHEMBL_ROOT = REPO_ROOT / "data" / "raw" / "other" / "chembl"

DEFAULT_DIAMOND_BIN = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------- external tools


def ensure_diamond() -> str:
    """Put DIAMOND on PATH, from `gradi-ortho` if it is not there already."""
    override = os.environ.get("GRADI_DIAMOND_BIN")
    if shutil.which("diamond") is None:
        for d in [Path(override)] if override else [DEFAULT_DIAMOND_BIN]:
            if (d / "diamond").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("diamond")
    if exe is None:
        sys.exit(
            "FAILED: no `diamond` on PATH.\n"
            "  it lives in the `gradi-ortho` env (osx-64, no arm64 build); point GRADI_DIAMOND_BIN\n"
            "  at that env's bin/ directory, or `conda create -n gradi-ortho -c bioconda diamond`."
        )
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip()
    return f"{exe}  ({ver})"


def find_db(required: bool) -> Path | None:
    """The ChEMBL SQLite dump, or None when every cache it feeds is already present.

    The dump is deleted after each run (30.5 GB, public, re-derivable), and the caches exist so a
    re-run does not need it back. So its absence is only fatal when something actually has to be
    re-extracted -- otherwise this returns None and the run proceeds from cache.
    """
    hits = sorted(CHEMBL_ROOT.glob("**/chembl_*.db"))
    if not hits:
        if not required:
            return None
        sys.exit(
            f"FAILED: no ChEMBL SQLite database under {CHEMBL_ROOT.relative_to(REPO_ROOT)}.\n"
            "  it is deleted after each run (30.5 GB, public, re-derivable). Restore it with:\n"
            "    tar -xzf data/raw/other/chembl/chembl_37_sqlite.tar.gz -C data/raw/other/chembl\n"
            "  see data/raw/other/chembl/SOURCE.md."
        )
    db = hits[0]
    assert_version(db)
    return db


def assert_version(db: Path) -> None:
    """The dump must be the version this stage CLAIMS, and nothing here proved that before.

    `CHEMBL_VERSION` is a hard-coded literal and `find_db` globs `chembl_*.db`, so extracting a
    different release beside this one would have been consumed silently and labelled 37 in every
    manifest. ChEMBL's SQLite ships a `version` table; that is the authority. The filename is the
    fallback, because it is what the glob actually matched.

    **`version` is NOT one row.** It holds 11 -- every upstream resource the release was built
    from (Bioassay Ontology, EFO, MeSH, UBERON, RDKit, InChI, COCONUT, GO, Swiss-Prot 2025_03,
    ChEMBL_Structure_Pipeline) beside the release itself, and `ChEMBL_37` is not first. A
    `fetchone()` returns `Bioassay Ontology 2.0`. Nor does `LIKE 'ChEMBL_%'` disambiguate, because
    `ChEMBL_Structure_Pipeline 1.2.0` matches it too; only `ChEMBL_<digits>` exactly is the release.
    """
    want = f"ChEMBL_{L.CHEMBL_VERSION}"
    seen, route = None, "filename"
    try:
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        try:
            names = [str(r[0]) for r in con.execute("SELECT name FROM version")]
        finally:
            con.close()
        releases = [n for n in names if re.fullmatch(r"ChEMBL_\d+", n)]
        if len(releases) == 1:
            seen, route = releases[0], "version table"
        else:
            sys.exit(
                f"FAILED: {db.name} has {len(releases)} release rows in its `version` table "
                f"({releases or 'none'}); expected exactly one. Rows seen: {names}"
            )
    except sqlite3.Error:
        m = re.search(r"chembl_(\d+)", db.name)
        seen = f"ChEMBL_{m.group(1)}" if m else db.name
    if seen.lower() != want.lower():
        sys.exit(
            f"FAILED: {db.name} is {seen} ({route}), but this stage is pinned to {want}.\n"
            f"  either restore the {want} dump or change CHEMBL_VERSION in src/ligandability.py "
            "and re-extract -- the cached extracts are version-specific."
        )
    say(f"  db version      : {seen}  (from the {route})")


# ---------------- extract


def _activity_where(kind: str) -> str:
    """The SQL predicate for one target track. See the module docstring on confidence_score."""
    if kind == "single":
        types = ",".join(f"'{t}'" for t in L.SINGLE_TYPES)
        conf = L.MIN_CONFIDENCE_SINGLE
    else:
        types = ",".join(f"'{t}'" for t in L.COMPLEX_TYPES)
        conf = L.MIN_CONFIDENCE_COMPLEX
    assays = ",".join(f"'{t}'" for t in L.ASSAY_TYPES)
    return f"""
        td.target_type IN ({types})
        AND a.confidence_score >= {conf}
        AND a.assay_type IN ({assays})
        AND act.pchembl_value IS NOT NULL
        AND act.standard_relation = '='
        AND act.data_validity_comment IS NULL
        AND (act.potential_duplicate IS NULL OR act.potential_duplicate = 0)
    """


def extract_chembl(db: Path | None, refresh: bool) -> tuple[pd.DataFrame, pd.DataFrame, Path]:
    """One pass over the dump: the target bridge, the ligand table and the DIAMOND FASTA."""
    t_path, l_path, f_path = (
        SCRATCH_DIR / "chembl_targets.tsv",
        SCRATCH_DIR / "chembl_ligands.tsv",
        SCRATCH_DIR / "chembl_targets.faa",
    )
    if not refresh and t_path.exists() and l_path.exists() and f_path.exists():
        say(f"  cached          : {t_path.name}, {l_path.name}, {f_path.name}")
        return L._read(t_path), L._read(l_path), f_path

    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    ligands = []
    for kind in ("single", "complex"):
        sql = f"""
            SELECT a.tid                                        AS tid,
                   COALESCE(mh.parent_molregno, act.molregno)   AS parent_molregno,
                   MAX(act.pchembl_value)                       AS pchembl,
                   SUM(CASE WHEN a.assay_type='B' THEN 1 ELSE 0 END) AS n_binding,
                   COUNT(*)                                     AS n_activities
              FROM activities act
              JOIN assays a            ON a.assay_id  = act.assay_id
              JOIN target_dictionary td ON td.tid      = a.tid
              LEFT JOIN molecule_hierarchy mh ON mh.molregno = act.molregno
             WHERE {_activity_where(kind)}
             GROUP BY 1, 2
        """
        df = pd.read_sql_query(sql, con)
        df["track"] = kind
        say(
            f"  {kind:8s}        : {len(df):>9,} (tid, compound) pairs  "
            f"{df.tid.nunique():>6,} targets  {df.parent_molregno.nunique():>9,} parent compounds"
        )
        ligands.append(df)
    lig = pd.concat(ligands, ignore_index=True)

    tids = ",".join(str(t) for t in sorted(lig.tid.unique()))
    tgt = pd.read_sql_query(
        f"""
        SELECT DISTINCT tc.component_id                AS component_id,
               td.tid                                  AS tid,
               td.chembl_id                            AS target_chembl_id,
               td.target_type                          AS target_type,
               td.tax_id                               AS tax_id,
               td.organism                             AS organism,
               td.pref_name                            AS pref_name,
               cs.accession                            AS accession,
               COALESCE(oc.l1, 'Unclassified')         AS superkingdom
          FROM target_dictionary td
          JOIN target_components  tc ON tc.tid = td.tid
          JOIN component_sequences cs ON cs.component_id = tc.component_id
          LEFT JOIN organism_class oc ON oc.tax_id = td.tax_id
         WHERE td.tid IN ({tids}) AND cs.sequence IS NOT NULL
        """,
        con,
    )
    say(
        f"  targets         : {len(tgt):>9,} (component, target) rows  "
        f"{tgt.component_id.nunique():>6,} distinct sequences"
    )

    # The FULL organism classification, not just the taxa our targets happen to use. It is the
    # only superkingdom source available to the BindingDB track (which has organism strings and no
    # taxonomy ids), and deriving that map from the hit targets alone covers 78 genera against
    # 1,493 taxa -- turning BindingDB's bacterial count into a lower bound of unknown tightness.
    oc = pd.read_sql_query("SELECT tax_id, l1, l2, l3 FROM organism_class", con)
    oc.to_csv(EVIDENCE_DIR / "organism_class.tsv", sep="\t", index=False)
    say(f"  organism_class  : {len(oc):>9,} taxa  "
        f"({int((oc.l1 == 'Bacteria').sum()):,} Bacteria) -> organism_class.tsv")

    seqs = pd.read_sql_query(
        f"SELECT component_id, sequence FROM component_sequences "
        f"WHERE component_id IN ({','.join(str(c) for c in sorted(tgt.component_id.unique()))})",
        con,
    )
    con.close()

    with open(f_path, "w") as fh:
        for cid, seq in zip(seqs.component_id, seqs.sequence):
            fh.write(f">{cid}\n{seq}\n")
    say(f"  fasta           : {len(seqs):>9,} records -> {f_path.name}")

    tgt.to_csv(t_path, sep="\t", index=False)
    lig.to_csv(l_path, sep="\t", index=False)
    return L._read(t_path), L._read(l_path), f_path


# ---------------- map by sequence


def write_query_fasta(species: str, limit: int | None, tmp: Path) -> tuple[Path, int]:
    """The proteome FASTA, built from proteome_<species>.tsv so it cannot drift from stage 00."""
    df = P.load(species)
    if limit:
        df = df.head(limit)
    path = tmp / f"{species}.faa"
    with open(path, "w") as fh:
        for ac, seq in zip(df.uniprot_ac, df.sequence):
            fh.write(f">{ac}\n{seq}\n")
    return path, len(df)


HIT_COLS = [
    "uniprot_ac",
    "component_id",
    "pident",
    "ppos",
    "alnlen",
    "qlen",
    "slen",
    "qcov",
    "scov",
    "evalue",
    "bitscore",
]


def diamond_map(species: str, faa: Path, target_faa: Path, threads: int, refresh: bool,
                limit: int | None, tmp: Path) -> pd.DataFrame:
    """DIAMOND blastp, one proteome against the ChEMBL target database.

    No identity or coverage filter is applied here beyond a permissive prefilter: the bands are
    applied in pandas so that `chembl_hits.tsv` stays the honest raw join and re-banding never
    needs DIAMOND again.
    """
    cache = SCRATCH_DIR / "hits" / f"{'smoke_' if limit else ''}hits_{species}.tsv.gz"
    cache.parent.mkdir(parents=True, exist_ok=True)
    if not refresh and cache.exists():
        say(f"  {species:14s}: cached {cache.name}")
        return pd.read_csv(cache, sep="\t")

    dbp = tmp / "chembl"
    if not dbp.with_suffix(".dmnd").exists():
        subprocess.run(
            ["diamond", "makedb", "--in", str(target_faa), "-d", str(dbp), "--quiet"], check=True
        )
    raw = tmp / f"{species}.m8"
    subprocess.run(
        [
            "diamond", "blastp", "-q", str(faa), "-d", str(dbp), "-o", str(raw),
            "--very-sensitive", "--id", "25", "--evalue", "1e-5",
            "--max-target-seqs", "100", "--outfmt", "6",
            "qseqid", "sseqid", "pident", "ppos", "length", "qlen", "slen",
            "qcovhsp", "scovhsp", "evalue", "bitscore",
            "--quiet", "--threads", str(threads),
        ],
        check=True,
    )
    df = pd.read_csv(raw, sep="\t", names=HIT_COLS)
    # multiple HSPs per pair -> keep the best-scoring one
    df = df.sort_values("bitscore", ascending=False).drop_duplicates(["uniprot_ac", "component_id"])
    df.to_csv(cache, sep="\t", index=False)
    say(
        f"  {species:14s}: {len(df):>8,} hits  "
        f"{df.uniprot_ac.nunique():>6,} proteins with any hit at >=25% id"
    )
    return df


# ---------------- scaffolds


def _scaffold_chunk(smiles: list[str]) -> list[str]:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
    out = []
    for smi in smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                out.append("")
                continue
            core = MurckoScaffold.GetScaffoldForMol(mol)
            if core.GetNumAtoms() == 0:
                out.append("")  # acyclic: one shared bucket, not a missing value
                continue
            out.append(Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(core)))
        except Exception:
            out.append("")
    return out


def compute_scaffolds(db: Path | None, molregnos: np.ndarray, threads: int,
                      refresh: bool) -> pd.DataFrame:
    """Bemis-Murcko generic scaffolds for every parent compound in play.

    Computed for the whole ligand table rather than only what the current bands reach, so that
    re-banding never needs the 30 GB dump back.
    """
    path = EVIDENCE_DIR / "scaffolds.tsv"
    if not refresh and path.exists():
        cached = L._read(path)
        # NOT a subset test against `molregnos`: ~1,500 parent compounds carry no structure in
        # ChEMBL at all, so they can never appear here and a strict subset check would never pass,
        # silently forcing a re-derivation (and demanding the 30 GB dump) on every run. The cache
        # is derived from chembl_ligands.tsv, which is itself cached alongside it; what is reported
        # instead is how much of the request it covers.
        have = set(cached.parent_molregno.dropna().astype(int))
        missing = len({int(m) for m in molregnos} - have)
        say(f"  cached          : {path.name} ({len(cached):,} compounds; {missing:,} of "
            f"{len(molregnos):,} requested have no structure in ChEMBL)")
        return cached

    if db is None:
        sys.exit(
            "FAILED: scaffolds must be computed but the ChEMBL database is not on disk.\n"
            "  restore it -- see data/raw/other/chembl/SOURCE.md -- or drop --refresh."
        )
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    ids = ",".join(str(int(m)) for m in sorted(molregnos))
    # standard_inchi_key and the molecule's own ChEMBL id are carried so that a later source
    # (BindingDB) can be compared against this one WITHOUT the 30 GB dump being restored. The
    # InChIKey is the better join of the two: it is structural, so it matches compounds that
    # carry no cross-reference.
    smi = pd.read_sql_query(
        f"SELECT cs.molregno AS parent_molregno, cs.canonical_smiles, cs.standard_inchi_key,"
        f"       md.chembl_id AS molecule_chembl_id "
        f"  FROM compound_structures cs "
        f"  LEFT JOIN molecule_dictionary md ON md.molregno = cs.molregno "
        f" WHERE cs.molregno IN ({ids}) AND cs.canonical_smiles IS NOT NULL",
        con,
    )
    con.close()
    say(f"  structures      : {len(smi):>9,} of {len(molregnos):,} parent compounds have SMILES")

    values = smi.canonical_smiles.tolist()
    chunks = [values[i : i + 5000] for i in range(0, len(values), 5000)]
    with mp.Pool(threads) as pool:
        results = pool.map(_scaffold_chunk, chunks)
    smi["scaffold"] = [s for chunk in results for s in chunk]

    n_acyclic = int((smi.scaffold == "").sum())
    say(
        f"  scaffolds       : {smi.scaffold.nunique():>9,} distinct  "
        f"({n_acyclic:,} acyclic -> one shared bucket)"
    )
    smi.to_csv(path, sep="\t", index=False)
    return smi


# ---------------- aggregate


def _pool_masks(tgt: pd.DataFrame, species: str) -> dict[str, pd.Series]:
    """Which ChEMBL targets belong to which bucket's pool, before any identity band."""
    single = tgt.target_type.isin(L.SINGLE_TYPES)
    bact = tgt.superkingdom == "Bacteria"
    binomial = L.SPECIES_ORGANISM[species]
    same = tgt.organism.fillna("").str.startswith(binomial)
    # `close` is a SUPERSET of `species` by construction. The bands alone do not nest -- a
    # same-species paralog at 45% identity is in `species` but below the 60% `close` cut -- and a
    # non-monotonic ladder is a footgun for anything downstream that subtracts one bucket from
    # another. The consequence, stated plainly: a same-species target between 40 and 60% identity
    # is counted as `close`. Being the same organism is treated as its own form of closeness.
    return {
        "direct": single & bact,
        "species": single & same,
        "close": single & (bact | same),
        "remote": single & bact,
        "human": single & (tgt.tax_id.astype(str) == "9606"),
        "allorg": single,
        "complex": tgt.target_type.isin(L.COMPLEX_TYPES),
    }


BAND = {
    "direct": L.DIRECT_PIDENT,
    "species": L.REMOTE_PIDENT,
    "close": L.CLOSE_PIDENT,
    "remote": L.REMOTE_PIDENT,
    "human": L.REMOTE_PIDENT,
    "allorg": L.REMOTE_PIDENT,
    "complex": L.REMOTE_PIDENT,
}


def aggregate(species: str, hits: pd.DataFrame, tgt: pd.DataFrame, lig: pd.DataFrame,
              scaf: dict[int, str], accessions: list[str], cutoff: float) -> pd.DataFrame:
    """One row per protein: union counts per bucket at one pChEMBL cutoff."""
    hits = hits[(hits.qcov >= L.MIN_QCOV) & (hits.scov >= L.MIN_SCOV)].copy()

    # component -> target rows, so one hit can reach several targets (and a complex several ways)
    bridge = tgt[["component_id", "tid", "target_chembl_id", "target_type", "organism",
                  "superkingdom", "tax_id"]]
    masks = _pool_masks(bridge, species)

    potent = lig[lig.pchembl >= cutoff]

    out = pd.DataFrame({"uniprot_ac": accessions})
    joined = hits.merge(bridge, on="component_id", how="inner")

    # provenance: the closest bacterial single-protein target
    prov = joined[
        joined.target_type.isin(L.SINGLE_TYPES) & (joined.superkingdom == "Bacteria")
    ]
    prov = prov.sort_values("pident", ascending=False).drop_duplicates("uniprot_ac")
    out = out.merge(
        prov[["uniprot_ac", "target_chembl_id", "pident", "organism"]].rename(
            columns={"target_chembl_id": "best_target", "pident": "best_pident",
                     "organism": "best_organism"}
        ),
        on="uniprot_ac", how="left",
    )
    direct = set(
        joined.loc[
            (joined.pident >= L.DIRECT_PIDENT) & joined.target_type.isin(L.SINGLE_TYPES),
            "uniprot_ac",
        ]
    )
    out["direct_hit"] = out.uniprot_ac.isin(direct)

    same_tids = set(bridge.loc[masks["species"], "tid"])
    for bucket, mask in masks.items():
        keep_tids = set(bridge.loc[mask, "tid"])
        band = joined.pident >= BAND[bucket]
        if bucket == "close":
            # nested by construction: same-species evidence enters `close` at the remote floor
            band = band | (joined.tid.isin(same_tids) & (joined.pident >= L.REMOTE_PIDENT))
        sub = joined[joined.tid.isin(keep_tids) & band]
        pairs = sub[["uniprot_ac", "tid"]].drop_duplicates().merge(
            potent[["tid", "parent_molregno", "pchembl"]], on="tid", how="inner"
        )
        if len(pairs):
            pairs["scaffold"] = pairs.parent_molregno.map(scaf)
            grp = pairs.groupby("uniprot_ac")
            agg = pd.DataFrame({
                f"{bucket}_n_compounds": grp.parent_molregno.nunique(),
                f"{bucket}_n_scaffolds": grp.scaffold.nunique(dropna=True),
                f"{bucket}_best_pchembl": grp.pchembl.max(),
            }).reset_index()
        else:
            agg = pd.DataFrame(columns=["uniprot_ac", f"{bucket}_n_compounds",
                                        f"{bucket}_n_scaffolds", f"{bucket}_best_pchembl"])
        out = out.merge(agg, on="uniprot_ac", how="left")

    # the identity to the closest human target -- provenance for the liability block, and the
    # quantity stage 05 measures independently against the real human proteome
    htids = set(bridge.loc[masks["human"], "tid"])
    hsub = joined[joined.tid.isin(htids) & (joined.pident >= L.REMOTE_PIDENT)]
    if len(hsub):
        hbest = hsub.sort_values("pident", ascending=False).drop_duplicates("uniprot_ac")
        out = out.merge(
            hbest[["uniprot_ac", "pident"]].rename(columns={"pident": "human_best_pident"}),
            on="uniprot_ac", how="left",
        )
    else:
        out["human_best_pident"] = np.nan

    # compounds tested at any potency over the bacterial pool -- screened-and-clean differs from
    # never-screened, and only this column can tell them apart
    keep = set(bridge.loc[masks["remote"], "tid"])
    sub = joined[joined.tid.isin(keep) & (joined.pident >= L.REMOTE_PIDENT)]
    tested = (
        sub[["uniprot_ac", "tid"]].drop_duplicates()
        .merge(lig[["tid", "parent_molregno"]], on="tid", how="inner")
        .groupby("uniprot_ac").parent_molregno.nunique().rename("n_compounds_tested").reset_index()
    )
    out = out.merge(tested, on="uniprot_ac", how="left")

    # the closest liganded complex, as provenance for the complex columns
    ctids = set(bridge.loc[masks["complex"], "tid"])
    csub = joined[joined.tid.isin(ctids) & (joined.pident >= L.REMOTE_PIDENT)]
    if len(csub):
        cbest = csub.sort_values("pident", ascending=False).drop_duplicates("uniprot_ac")
        out = out.merge(
            cbest[["uniprot_ac", "target_chembl_id"]].rename(
                columns={"target_chembl_id": "complex_best_target"}),
            on="uniprot_ac", how="left",
        )
    else:
        out["complex_best_target"] = ""

    out = out.drop(columns=["complex_n_scaffolds", "complex_best_pchembl"], errors="ignore")
    # every bucket sits inside `remote`, so a tighter one may never exceed it

    for col in L.DELIVERABLE_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    for col in out.columns:
        if col.endswith(("_n_compounds", "_n_scaffolds", "_tested")):
            out[col] = out[col].fillna(0).astype(int)
    for col in ("best_target", "best_organism", "complex_best_target"):
        out[col] = out[col].fillna("")
    return out[list(L.DELIVERABLE_COLUMNS)]


# ---------------- main


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--species", nargs="+", default=list(L.SPECIES), choices=list(L.SPECIES))
    ap.add_argument("--pchembl", type=float, default=L.PCHEMBL_HEADLINE,
                    help="headline potency cutoff (default 6.0 = 1 uM). NOT a knob to tune.")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, help="only the first N proteins per species (smoke test)")
    ap.add_argument("--refresh", action="store_true",
                    help="discard the SQL extract, the DIAMOND hits and the scaffold cache")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    pre = "smoke_" if args.limit else ""
    started = datetime.now(timezone.utc)

    rule("=")
    say("STAGE 06 - ligandability, part 1: measured bioactivity from ChEMBL")
    rule("=")
    say(f"  source   : ChEMBL {L.CHEMBL_VERSION} SQLite")
    say(f"  assays   : type {'/'.join(L.ASSAY_TYPES)}, pchembl_value present, relation '=', valid, not duplicate")
    say(f"  targets  : {', '.join(L.SINGLE_TYPES)} (confidence >={L.MIN_CONFIDENCE_SINGLE})"
        f" + {', '.join(L.COMPLEX_TYPES)} (confidence >={L.MIN_CONFIDENCE_COMPLEX})")
    say(f"  excluded : {', '.join(L.EXCLUDED_TYPES)} -- confidence 4-5, target ambiguous within the group")
    say(f"  potency  : headline pChEMBL >= {args.pchembl}, sweep {L.PCHEMBL_CUTOFFS}")
    say(f"  bands    : direct >={L.DIRECT_PIDENT} - close >={L.CLOSE_PIDENT} - remote >={L.REMOTE_PIDENT}"
        f"  (qcov >={L.MIN_QCOV}, scov >={L.MIN_SCOV})")
    say(f"  ligands  : parent molregno, then Bemis-Murcko generic scaffold")
    say(f"  species  : {', '.join(args.species)}   (human is the reference pool, not a query)")
    say(f"  output   : {EVIDENCE_DIR.relative_to(REPO_ROOT)}/chembl_<species>.tsv")

    if args.dry_run:
        say("\n  --dry-run: nothing extracted, nothing written.")
        for sp in args.species:
            say(f"    would write {(EVIDENCE_DIR / f'chembl_{sp}.tsv').relative_to(REPO_ROOT)}")
        return

    say()
    rule()
    say("TOOLS")
    rule()
    say(f"  diamond         : {ensure_diamond()}")
    cached = all((SCRATCH_DIR / f).exists() for f in (
        "chembl_targets.tsv", "chembl_ligands.tsv", "chembl_targets.faa",
        "scaffolds.tsv", "organism_class.tsv"))
    db = find_db(required=args.refresh or not cached)
    if db is None:
        say("  chembl          : not on disk -- running from cache (see "
            "data/raw/other/chembl/SOURCE.md)")
    else:
        say(f"  chembl          : {db.relative_to(REPO_ROOT)}  ({db.stat().st_size / 1e9:.1f} GB)")

    say()
    rule()
    say("EXTRACT - one pass over the dump")
    rule()
    tgt, lig, target_faa = extract_chembl(db, args.refresh)
    n_b = int(lig.n_binding.fillna(0).sum())
    n_a = int(lig.n_activities.fillna(0).sum())
    say(f"  assay split     : {n_b:,} binding of {n_a:,} activities "
        f"({100 * n_b / max(n_a, 1):.1f}% B, rest F -- F is where bacterial enzymology lives)")
    for name, mask in [("true Bacteria", tgt.superkingdom == "Bacteria"),
                       ("human", tgt.tax_id.astype(str) == "9606")]:
        say(f"  {name:15s} : {tgt.loc[mask, 'tid'].nunique():>6,} targets")

    say()
    rule()
    say("SCAFFOLDS - collapsing salts, then counting chemical ideas")
    rule()
    scaf_df = compute_scaffolds(db, lig.parent_molregno.dropna().unique(), args.threads,
                                args.refresh)
    scaf = dict(zip(scaf_df.parent_molregno.astype(int), scaf_df.scaffold))

    say()
    rule()
    say("MAP - DIAMOND, our proteomes against the ChEMBL target sequences")
    rule()
    tmp = Path(tempfile.mkdtemp(prefix="gradi06_"))
    hits_by_sp, expected = {}, {}
    for sp in args.species:
        faa, n = write_query_fasta(sp, args.limit, tmp)
        expected[sp] = n
        hits_by_sp[sp] = diamond_map(sp, faa, target_faa, args.threads, args.refresh,
                                     args.limit, tmp)

    all_hits = pd.concat(
        [h.assign(species=sp) for sp, h in hits_by_sp.items()], ignore_index=True
    )
    all_hits = all_hits[["species"] + HIT_COLS]
    all_hits.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}chembl_hits.tsv", sep="\t", index=False)

    say()
    rule()
    say("COVERAGE GUARDS - what qcov/scov cost, since v1 used qcov alone")
    rule()
    for sp in args.species:
        h = hits_by_sp[sp]
        band = h[h.pident >= L.REMOTE_PIDENT]
        q = band[band.qcov >= L.MIN_QCOV]
        both = q[q.scov >= L.MIN_SCOV]
        say(f"  {sp:14s}: {len(band):>6,} hits >={L.REMOTE_PIDENT:.0f}% id  ->  "
            f"{len(q):>6,} after qcov>={L.MIN_QCOV:.0f}  ->  {len(both):>6,} after scov>="
            f"{L.MIN_SCOV:.0f}   (scov drops {len(q) - len(both):,}, "
            f"{100 * (len(q) - len(both)) / max(len(q), 1):.1f}%)")

    say()
    rule()
    say("AGGREGATE")
    rule()
    rows, frames = [], {}
    for sp in args.species:
        df = P.load(sp)
        if args.limit:
            df = df.head(args.limit)
        accs = df.uniprot_ac.tolist()
        out = aggregate(sp, hits_by_sp[sp], tgt, lig, scaf, accs, args.pchembl)
        frames[sp] = out
        path = (SCRATCH_DIR / f"smoke_chembl_{sp}.tsv") if args.limit \
            else (EVIDENCE_DIR / f"chembl_{sp}.tsv")
        out.to_csv(path, sep="\t", index=False)
        n_any = int((out.remote_n_compounds > 0).sum())
        n_direct = int(out.direct_hit.sum())
        say(f"  {sp:14s}: {len(out):>6,} proteins  "
            f"{n_direct:>4,} direct-hit  "
            f"direct {int((out.direct_n_compounds > 0).sum()):>4,}  "
            f"species {int((out.species_n_compounds > 0).sum()):>4,}  "
            f"close {int((out.close_n_compounds > 0).sum()):>4,}  "
            f"remote {n_any:>4,}  "
            f"human {int((out.human_n_compounds > 0).sum()):>4,}  "
            f"complex {int((out.complex_n_compounds > 0).sum()):>4,}")
        rows.append({
            "species": sp, "n": len(out), "n_expected": expected[sp],
            "n_direct_hit": n_direct,
            "n_direct": int((out.direct_n_compounds > 0).sum()), "n_species": int((out.species_n_compounds > 0).sum()),
            "n_close": int((out.close_n_compounds > 0).sum()), "n_remote": n_any,
            "n_human": int((out.human_n_compounds > 0).sum()),
            "n_complex": int((out.complex_n_compounds > 0).sum()),
            "pchembl": args.pchembl, "chembl_version": L.CHEMBL_VERSION,
            "path": str(path.relative_to(REPO_ROOT)),
            "run_at": started.isoformat(timespec="seconds"),
        })
    pd.DataFrame(rows).to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}manifest.tsv", sep="\t", index=False)

    say()
    rule()
    say("CUTOFF SENSITIVITY - reported, never optimised")
    rule()
    sens = []
    for cut in L.PCHEMBL_CUTOFFS:
        for sp in args.species:
            df = P.load(sp)
            if args.limit:
                df = df.head(args.limit)
            o = (frames[sp] if cut == args.pchembl else
                 aggregate(sp, hits_by_sp[sp], tgt, lig, scaf, df.uniprot_ac.tolist(), cut))
            sens.append({"pchembl": cut, "species": sp,
                         "n_species": int((o.species_n_compounds > 0).sum()),
                         "n_close": int((o.close_n_compounds > 0).sum()),
                         "n_remote": int((o.remote_n_compounds > 0).sum()),
                         "n_human": int((o.human_n_compounds > 0).sum())})
    sens_df = pd.DataFrame(sens)
    sens_df.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}cutoff_sensitivity.tsv", sep="\t", index=False)
    for cut in L.PCHEMBL_CUTOFFS:
        s = sens_df[sens_df.pchembl == cut]
        say(f"  pChEMBL >={cut:.0f}    : " + "  ".join(
            f"{r.species[:2]} species {r.n_species:>4,} remote {r.n_remote:>4,}"
            for r in s.itertuples()))

    say()
    rule()
    say("CONTROL - v1's pre-registered numbers")
    rule()
    ctl = []
    if args.limit:
        say(f"  skipped: --limit {args.limit} is a subset of the proteome, so v1's whole-proteome"
            f" counts are not comparable.")
    for sp, v1 in (("kpneumoniae", L.V1_KP_POTENT), ("ecoli", L.V1_EC_POTENT)):
        if sp not in frames or args.limit:
            continue
        v2 = int((frames[sp].remote_n_compounds > 0).sum())
        ok = abs(v2 - v1) <= L.TOLERANCE_FRACTION * v1
        ctl.append({"quantity": f"{sp} proteins with a potent bacterial ligand",
                    "v1": v1, "v2": v2, "within_tolerance": ok})
        say(f"  {sp:14s}: v1 {v1:>5,}   v2 {v2:>5,}   "
            f"{'ok' if ok else 'OUT OF TOLERANCE'}  (+-{L.TOLERANCE_FRACTION:.0%})")
    ctl.append({"quantity": "v1 uncorrected (no identity floor, non-human bucket)",
                "v1": L.V1_KP_POTENT_UNCORRECTED, "v2": "", "within_tolerance": ""})
    pd.DataFrame(ctl).to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}control.tsv", sep="\t", index=False)

    say()
    rule()
    say("SPOT CHECKS")
    rule()
    fails: list[str] = []
    ac, want_target, want_pident = L.V1_SPOT_KP_BLA
    if "kpneumoniae" in frames:
        h = hits_by_sp["kpneumoniae"]
        row = h[h.uniprot_ac == ac].merge(tgt, on="component_id").sort_values(
            "pident", ascending=False)
        if len(row) and row.iloc[0].accession == want_target and row.iloc[0].pident >= want_pident:
            say(f"  {ac} -> {want_target} at {row.iloc[0].pident:.1f}% id  "
                f"({row.iloc[0].organism})  -- v1's validation case holds")
        else:
            got = f"{row.iloc[0].accession} at {row.iloc[0].pident:.1f}%" if len(row) else "nothing"
            fails.append(f"{ac} should map to {want_target} at {want_pident}%, got {got}")

    for sp in args.species:
        pro = P.load(sp)
        if args.limit:
            pro = pro.head(args.limit)
        named = frames[sp].merge(pro[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left")
        for gene, expectation in (("clpP", "human"), ("ftsZ", "no human")):
            r = named[named.gene_name == gene]
            if not len(r):
                continue
            hn = int(r.iloc[0].human_n_compounds)
            bn = int(r.iloc[0].remote_n_compounds)
            say(f"  {sp:14s} {gene:5s}: human {hn:>5,} compounds   bacterial {bn:>5,}")
            if gene == "clpP" and expectation == "human" and hn == 0 and sp == "kpneumoniae":
                fails.append(
                    "Kp clpP has zero human evidence -- stage 05 measures a 56.3% human CLPP "
                    "ortholog and ONC212 is a characterised human ClpP ligand, so the human "
                    "bucket is broken"
                )

    say()
    rule()
    say("CROSS-CHECK - human identity, against stage 05's own DIAMOND run")
    rule()
    try:
        from src import orthology as O

        nb = O.load_neighbors()
        nb = nb[nb.target_species == "human"].sort_values("pident", ascending=False)
        nb = nb.drop_duplicates("query_ac")[["query_ac", "pident"]].rename(
            columns={"query_ac": "uniprot_ac", "pident": "stage05_human_pident"})
        for sp in args.species:
            m = frames[sp].merge(nb, on="uniprot_ac", how="inner")
            m = m[m.human_best_pident.notna() & m.stage05_human_pident.notna()]
            if len(m) < 5:
                say(f"  {sp:14s}: only {len(m)} proteins comparable -- not reported")
                continue
            r = float(np.corrcoef(m.human_best_pident, m.stage05_human_pident)[0, 1])
            md = float((m.human_best_pident - m.stage05_human_pident).abs().median())
            say(f"  {sp:14s}: n={len(m):>4,}  r={r:.3f}  median |diff| {md:.1f} pp"
                f"   (independent routes to the same quantity)")
    except FileNotFoundError:
        say("  stage 05 not built -- skipped (it is a soft dependency)")

    say()
    rule()
    say("ASSERTIONS")
    rule()
    for sp in args.species:
        o = frames[sp]
        if len(o) != expected[sp]:
            fails.append(f"{sp} has {len(o)} rows, expected {expected[sp]}")
        if not (o.remote_n_compounds >= o.close_n_compounds).all():
            fails.append(f"{sp}: remote < close for some protein -- the buckets are not nested")
        if not (o.close_n_compounds >= o.species_n_compounds).all():
            fails.append(f"{sp}: close < species for some protein -- the buckets are not nested")
        if not (o.remote_n_compounds >= o.direct_n_compounds).all():
            fails.append(f"{sp}: remote < direct for some protein -- the buckets are not nested")
        for b in L.BUCKETS + ("human", "allorg"):
            if not (o[f"{b}_n_scaffolds"] <= o[f"{b}_n_compounds"]).all():
                fails.append(f"{sp}: {b}_n_scaffolds exceeds {b}_n_compounds")
        say(f"  {sp:14s}: rows ok, buckets nested, scaffolds <= compounds")
    ratio_rows = []
    for sp in args.species:
        o = frames[sp]
        m = o.remote_n_compounds > 0
        if m.any():
            r = (o.loc[m, "remote_n_compounds"] / o.loc[m, "remote_n_scaffolds"].clip(lower=1)).median()
            ratio_rows.append(f"{sp} {r:.1f}")
    say(f"  compounds per scaffold (median, remote): {', '.join(ratio_rows)}"
        f"   -- near 1.0 everywhere would mean the scaffold step is not working")

    say()
    rule()
    say("OUTPUTS")
    rule()
    seen = set()
    for r in rows:
        seen.add(REPO_ROOT / r["path"])
        say(f"  {r['path']:<54} {(REPO_ROOT / r['path']).stat().st_size / 1e3:>8.1f} kB")
    for f in sorted([*EVIDENCE_DIR.glob("*"), *SCRATCH_DIR.glob("*")]):
        keep = f.name.startswith(pre) or f.name == "scaffolds.tsv" if pre else \
            not f.name.startswith("smoke_")
        if f.is_file() and f not in seen and keep:
            say(f"  {str(f.relative_to(REPO_ROOT)):<54} {f.stat().st_size / 1e3:>8.1f} kB")

    shutil.rmtree(tmp, ignore_errors=True)
    if fails:
        sys.exit("FAILED:\n  - " + "\n  - ".join(fails))
    say()
    rule("=")
    say(f"stage 06 (chembl) complete in {(datetime.now(timezone.utc) - started).seconds / 60:.1f} min.")
    rule("=")


if __name__ == "__main__":
    main()
