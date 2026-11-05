"""Which proteins are the recorded TARGET of an approved drug? -- clinical validation, not potency.

Ligandability, part 1c. `chembl.py` and `ligands.py` count compounds with a measured pChEMBL, so
they cannot see a target whose drugs were characterised any other way: MIC at the whole-organism
level, covalent acylation, a ribosome stored as a protein-nucleic-acid complex. The result is that
ribosomal proteins, PBPs and RNA polymerase read as empty while being the best-validated
antibacterial targets there are. ChEMBL's `drug_mechanism` table records, for every approved drug,
the target it acts on, whatever the assay -- so it says what the potency counts cannot.

Reads    ChEMBL REST API (network)                  mechanisms, targets and components, molecules
         data/processed/proteomes/                  the three bacterial proteomes
Writes   data/processed/ligands/evidence/mechanisms_<species>.tsv     one row per protein
         data/processed/ligands/evidence/mechanisms_audit.tsv          what could and could not be mapped
         data/processed/ligands/scratch/mechanisms/                   cached API pages

**`n_approved_drugs_own` is the deliverable column** (`ligands.py` reads this file and exits if it
is missing). It counts DISTINCT approved drugs (ChEMBL `max_phase` 4, collapsed to the parent
molecule) whose recorded mechanism targets a protein that THIS protein is >= 95% identical to, with
both coverages >= 50% -- the house "this is the same protein" rule (`L.DIRECT_PIDENT`, `L.MIN_QCOV`,
`L.MIN_SCOV`, imported). Unlike the `exact` ligand count it does not require the same species: a drug
recorded against a generic "Bacterial 70S ribosome" built from E. coli subunits is a drug for a
K. pneumoniae ribosomal protein at 99% identity, and that is the point of the column.

**What this CANNOT see, and the audit says by how much.** Only mechanism targets that carry component
sequences can be mapped. A target stored as a PROTEIN FAMILY or SELECTIVITY GROUP (many beta-lactam
mechanisms are "penicillin-binding proteins") has none, so its drugs are unmappable -- they are
counted in `mechanisms_audit.tsv` by ATC class and never silently dropped. A drug that acts on a
non-protein (vancomycin on D-Ala-D-Ala, polymyxin on LPS) correctly maps to nothing. **A 0 here is
"no approved drug is recorded against this protein", NOT "not druggable".**

Run with the `gradi` env (DIAMOND is borrowed from `gradi-ortho`; `GRADI_DIAMOND_BIN` overrides):
    python scripts/ligands/mechanisms.py [--refresh] [-q]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import matrices as M  # noqa: E402
from src import precedents as PR  # noqa: E402
from src import proteomes as P  # noqa: E402

API = "https://www.ebi.ac.uk/chembl/api/data/"
OUT_DIR = REPO_ROOT / "data" / "processed" / "ligands"
EVIDENCE_DIR = OUT_DIR / "evidence"
CACHE_DIR = OUT_DIR / "scratch" / "mechanisms"
SPECIES = ("kpneumoniae", "ecoli", "saureus")
APPROVED = 4                      # ChEMBL max_phase
BATCH = 80                        # ids per `__in` request, well inside URL limits
# systemic antibacterials and antimycobacterials, plus the topical classes (retapamulin, finafloxacin,
# ozenoxacin carry D06/S02 codes, not J01, and read as "non-antibacterial" without them)
ANTIBACTERIAL_ATC = ("J01", "J04", "A07AA", "D06A", "D06BA", "G01AA", "R02AB", "S01AA", "S02AA")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def get_json(endpoint: str, params: dict) -> dict:
    url = API + endpoint + ".json?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception:
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("unreachable")


def fetch_all(endpoint: str, key: str, params: dict, cache: str, refresh: bool) -> list[dict]:
    """Every row of a paged endpoint, cached as one JSON file."""
    path = CACHE_DIR / f"{cache}.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    rows, offset = [], 0
    while True:
        d = get_json(endpoint, {**params, "limit": 1000, "offset": offset})
        rows += d[key]
        total = d["page_meta"]["total_count"]
        offset += 1000
        if offset >= total:
            break
    if len(rows) != total:
        sys.exit(f"FATAL {endpoint}: fetched {len(rows)} rows but the API reports {total}")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows))
    return rows


def fetch_in(endpoint: str, key: str, field: str, ids: list[str], cache: str,
             refresh: bool) -> list[dict]:
    """Rows for a list of ids through `<field>__in`, batched and cached."""
    path = CACHE_DIR / f"{cache}.json"
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    rows = []
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        d = get_json(endpoint, {f"{field}__in": ",".join(chunk), "limit": 1000})
        rows += d[key]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows))
    return rows


def main() -> int:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refresh", action="store_true", help="discard the cached API pages")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("mechanisms.py -- approved-drug mechanism targets, mapped onto the three bacterial proteomes")
    rule("=")
    say(f"  rule   max_phase >= {APPROVED}, parent molecule, identity >= {L.DIRECT_PIDENT:.0f}%, "
        f"coverage q>={L.MIN_QCOV:.0f} s>={L.MIN_SCOV:.0f} (imported, not restated)")
    rule()

    # ---- 1. every approved-drug mechanism
    mech = pd.DataFrame(fetch_all("mechanism", "mechanisms", {}, "mechanisms", a.refresh))
    say(f"  {len(mech):,} mechanism rows in ChEMBL")
    mech["max_phase"] = pd.to_numeric(mech["max_phase"], errors="coerce")
    appr = mech[mech["max_phase"] >= APPROVED].copy()
    appr["drug"] = appr["parent_molecule_chembl_id"].fillna(appr["molecule_chembl_id"])
    say(f"  {len(appr):,} are approved ({appr.drug.nunique():,} distinct drugs); "
        f"{int(appr.target_chembl_id.isna().sum()):,} rows carry no target at all")
    appr = appr[appr["target_chembl_id"].notna()]

    # ---- 2. the molecules, for names and the ATC class
    mol = pd.DataFrame(fetch_in("molecule", "molecules", "molecule_chembl_id",
                                sorted(appr.drug.unique()), "molecules", a.refresh))
    name = dict(zip(mol["molecule_chembl_id"], mol["pref_name"].fillna(mol["molecule_chembl_id"])))
    atc = {r.molecule_chembl_id: (r.atc_classifications or []) for r in mol.itertuples()}
    antibacterial = {k for k, v in atc.items() if any(c.startswith(ANTIBACTERIAL_ATC) for c in v)}
    say(f"  {len(antibacterial):,} of those drugs carry a systemic-antibacterial ATC code "
        f"({'/'.join(ANTIBACTERIAL_ATC)})")

    # ---- 3. the targets and their component sequences
    tgt = fetch_in("target", "targets", "target_chembl_id", sorted(appr.target_chembl_id.unique()),
                   "targets", a.refresh)
    comp_of, ttype = {}, {}
    for t in tgt:
        ttype[t["target_chembl_id"]] = t["target_type"]
        comp_of[t["target_chembl_id"]] = [c["accession"] for c in t["target_components"]
                                          if c.get("accession")]
    accs = sorted({x for v in comp_of.values() for x in v})
    comps = fetch_in("target_component", "target_components", "accession", accs, "components",
                     a.refresh)
    seq = {c["accession"]: c["sequence"] for c in comps if c.get("sequence")}
    say(f"  {len(tgt):,} targets, {len(accs):,} component accessions, {len(seq):,} with a sequence")

    # which approved mechanisms can be mapped at all
    appr["target_type"] = appr.target_chembl_id.map(ttype)
    appr["n_seq_components"] = appr.target_chembl_id.map(
        lambda t: sum(1 for x in comp_of.get(t, []) if x in seq))
    appr["mappable"] = appr["n_seq_components"] > 0
    appr["antibacterial"] = appr.drug.isin(antibacterial)
    audit = (appr.groupby(["target_type", "mappable", "antibacterial"])
             .agg(mechanisms=("drug", "size"), drugs=("drug", "nunique")).reset_index())
    rule()
    say("  APPROVED MECHANISMS by target type -- the mappable share is what this column can see")
    say(audit.to_string(index=False))
    ab = appr[appr.antibacterial]
    say(f"  antibacterial drugs: {ab.drug.nunique():,} total, {ab[ab.mappable].drug.nunique():,} with "
        f"a mappable target, {ab[~ab.mappable].drug.nunique():,} with none (families, groups, "
        f"non-protein targets)")
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    audit.to_csv(EVIDENCE_DIR / "mechanisms_audit.tsv", sep="\t", index=False)

    # ---- 4. DIAMOND the three proteomes against the component sequences
    rows = appr[appr.mappable]
    drug_of_comp: dict[str, set] = {}
    for r in rows.itertuples():
        for x in comp_of[r.target_chembl_id]:
            if x in seq:
                drug_of_comp.setdefault(x, set()).add(r.drug)
    subject = {x: seq[x] for x in drug_of_comp}
    rule()
    say(f"  {len(subject):,} component sequences are the subject database")
    ev_rows = {}
    with tempfile.TemporaryDirectory() as td:
        subj, db = Path(td) / "subject.faa", Path(td) / "db"
        subj.write_text("".join(f">{k}\n{v}\n" for k, v in subject.items()))
        subprocess.run([PR.diamond_bin(), "makedb", "--in", str(subj), "-d", str(db), "--quiet"],
                       check=True)
        for sp in SPECIES:
            d = P.load(sp)
            q, out = Path(td) / f"{sp}.faa", Path(td) / f"{sp}.tsv"
            q.write_text("".join(f">{k}\n{v}\n" for k, v in zip(d.uniprot_ac, d.sequence.astype(str))))
            subprocess.run([PR.diamond_bin(), "blastp", "-q", str(q), "-d", str(db), "-o", str(out),
                            "--very-sensitive", "--id", str(L.DIRECT_PIDENT - 5), "--evalue", "1e-5",
                            "--max-target-seqs", "100", "--outfmt", "6", "qseqid", "sseqid",
                            "pident", "qcovhsp", "scovhsp", "bitscore", "--quiet",
                            "--threads", str(a.threads)], check=True)
            h = (pd.read_csv(out, sep="\t", names=["q", "s", "pident", "qcov", "scov", "bits"])
                 if out.exists() and out.stat().st_size else
                 pd.DataFrame(columns=["q", "s", "pident", "qcov", "scov", "bits"]))
            h = h[(h.pident >= L.DIRECT_PIDENT) & (h.qcov >= L.MIN_QCOV) & (h.scov >= L.MIN_SCOV)]
            per = {}
            for r in h.itertuples():
                per.setdefault(r.q, []).append((r.s, r.pident))
            recs = []
            for ac in d.uniprot_ac:
                hits = per.get(ac, [])
                drugs = sorted({x for s, _ in hits for x in drug_of_comp[s]})
                recs.append({"uniprot_ac": ac,
                             "n_approved_drugs_own": len(drugs),
                             "n_antibacterial_drugs_own": sum(1 for x in drugs if x in antibacterial),
                             "approved_drugs": ";".join(name.get(x, x) for x in drugs),
                             "best_pident": max((p for _, p in hits), default=float("nan"))})
            t = M.reindex(pd.DataFrame(recs), sp)
            t.to_csv(EVIDENCE_DIR / f"mechanisms_{sp}.tsv", sep="\t", index=False)
            ev_rows[sp] = t
            say(f"  {sp:<12} {len(t):>5,} proteins   with an approved drug "
                f"{int((t.n_approved_drugs_own > 0).sum()):>4,} "
                f"({100 * (t.n_approved_drugs_own > 0).mean():.1f}%)   with an antibacterial "
                f"{int((t.n_antibacterial_drugs_own > 0).sum()):>4,}")

    # ---- 5. control: the join must reach classes that are certainly there
    rule()
    ec = ev_rows["ecoli"].merge(P.load("ecoli")[["uniprot_ac", "gene_name"]], on="uniprot_ac")
    probe = {"gyrA": "quinolones", "gyrB": "aminocoumarins", "folA": "trimethoprim",
             "rpsL": "aminoglycosides (30S)", "rplC": "50S", "rpoB": "rifamycins"}
    say("  CONTROL (E. coli): proteins that should carry an approved drug")
    for g, why in probe.items():
        r = ec[ec.gene_name == g]
        n = int(r.n_approved_drugs_own.iloc[0]) if len(r) else -1
        say(f"    {g:<6} {n:>3} approved drugs   ({why})")
    if int((ec.n_approved_drugs_own > 0).sum()) == 0:
        sys.exit("FAILED control: no E. coli protein carries an approved drug -- the join is broken")
    say(f"\n  wrote evidence/mechanisms_<species>.tsv x3 and evidence/mechanisms_audit.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
