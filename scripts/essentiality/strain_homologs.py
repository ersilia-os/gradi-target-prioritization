"""Which anchor proteins have a MEASURED counterpart in a screened strain of their own species?

    data/processed/essentiality/evidence/strain_homologs_<species>.tsv   one row per protein

K. pneumoniae is the reason this exists. It is absent from DEG and from OGEE, so the anchor has no
measured essentiality of its own -- exact-sequence matching against ALL of DEG reaches 29 of 5,728
Kp proteins (0.5%). But three published TraDIS screens WERE run on other K. pneumoniae strains,
and a protein does not stop being itself between strains. This records, per anchor protein,
whether such a counterpart was measured and what it said.

THIS IS A FLAG, NOT A LABEL TRANSFER
--------------------------------------
`CLAUDE.md` records that transferring the screens' measured labels onto the anchor was considered
and REJECTED on the owner's instruction: exact-sequence transfer recovers only 13.7-71.6% of rows,
so it would silently mislabel thousands of proteins. Nothing here transfers a label. The measured
call never enters `essentiality_consensus`, which stays pure prediction. What it feeds is
`essentiality_evidence` -- the COUNT of independent measurements and whether they agree -- so the
rejected thing stays rejected. Read `n_strains_measured` as "how many screened relatives covered
this protein", never as "this protein's essentiality".

MATCHING: THE HOUSE `exact` RULE
----------------------------------
DIAMOND, >= `DIRECT_PIDENT` (95%) identity and >= `MIN_SCOV` (50%) of the SUBJECT aligned -- the
same rule `src/precedents.py` uses for "this protein across strains" and `studiedness/transfer.py`
for `_own`. Best hit per (anchor protein, strain); ties broken by bitscore then identity.

There is deliberately NO query-coverage floor, for the reason `pockets/holo.py` documents: a
measured gene that covers one domain of a longer anchor protein is still a measurement of that
domain. `MIN_SCOV` is what stops a fragment counting.

Run in `gradi`; DIAMOND is borrowed from `gradi-ortho` (`GRADI_DIAMOND_BIN` overrides).
    python scripts/essentiality/strain_homologs.py
    python scripts/essentiality/strain_homologs.py --species kpneumoniae -q
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import strains as S  # noqa: E402
from src.precedents import diamond_bin  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
TRAINING_DIR = OUT_DIR / "training_sets"

# screen training set -> the strain proteome its `key` column is keyed on. Only same-species
# screens belong here: a measurement in another species is not "this protein measured elsewhere",
# it is a different protein, and `geptop`/`screens_ess_mean` already carry cross-species signal.
SAME_SPECIES_SCREENS: dict[str, dict[str, str]] = {
    "kpneumoniae": {
        "essential_kpneumoniae_ecl8_tradis": "kpneumoniae__ecl8__GCA_000315385.1",
        "essential_kpneumoniae_rh201207_tradis": "rh201207_bruchmann",
        "essential_kpneumoniae_atcc43816_tradis": "kpneumoniae__kppr1__GCF_000742755.1",
    },
    # E. coli and S. aureus have DEG on the exact anchor strain (96.6% / 92.7%), so a proxy adds
    # nothing there. Left empty deliberately rather than omitted, so the asymmetry is visible.
    "ecoli": {},
    "saureus": {},
}

COLUMNS = ["uniprot_ac", "n_strains_measured", "n_strains_essential", "strains_measured"]
VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def best_hits(query: pd.DataFrame, subject: pd.DataFrame, threads: int) -> pd.DataFrame:
    """Best DIAMOND hit per query above the house `exact` floors. (uniprot_ac, key, pident)."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        qf, sf = td / "q.faa", td / "s.faa"
        for path, df in ((qf, query), (sf, subject)):
            path.write_text("".join(f">{a}\n{s}\n" for a, s in
                                    zip(df["uniprot_ac"], df["sequence"])))
        subprocess.run([diamond_bin(), "makedb", "--in", str(sf), "-d", str(td / "db"),
                        "--quiet"], check=True)
        out = td / "hits.tsv"
        subprocess.run(
            [diamond_bin(), "blastp", "-q", str(qf), "-d", str(td / "db"), "-o", str(out),
             "--outfmt", "6", "qseqid", "sseqid", "pident", "scovhsp", "bitscore",
             "--very-sensitive", "--max-target-seqs", "5", "--threads", str(threads),
             "--quiet"], check=True)
        if not out.stat().st_size:
            return pd.DataFrame(columns=["uniprot_ac", "key", "pident"])
        h = pd.read_csv(out, sep="\t", header=None,
                        names=["uniprot_ac", "key", "pident", "scovhsp", "bitscore"])
    h = h[(h.pident >= L.DIRECT_PIDENT) & (h.scovhsp >= L.MIN_SCOV)]
    h = h.sort_values(["uniprot_ac", "bitscore", "pident"], ascending=[True, False, False])
    return h.drop_duplicates("uniprot_ac")[["uniprot_ac", "key", "pident"]]


def build(species: str, threads: int) -> pd.DataFrame:
    prot = P.load(species)[["uniprot_ac", "sequence"]]
    screens = SAME_SPECIES_SCREENS[species]
    per_strain: dict[str, pd.Series] = {}

    for screen, label in screens.items():
        tpath = TRAINING_DIR / f"{screen}.tsv"
        if not tpath.exists():
            say(f"    {screen}: training set absent -- skipped (run essentiality/screens.py)")
            continue
        lab = pd.read_csv(tpath, sep="\t", dtype={"key": str})
        calls = lab.set_index("key")["label"].astype(float)
        sub = S.load_frame(label)
        hits = best_hits(prot, sub, threads)
        # a hit only counts if the matched strain gene was actually MEASURED by that screen
        hits["call"] = hits["key"].map(calls)
        hits = hits[hits["call"].notna()]
        per_strain[screen] = hits.set_index("uniprot_ac")["call"]
        say(f"    {screen:40s} {len(sub):>5} proteins · {len(calls):>5} measured · "
            f"{len(hits):>5} anchor proteins matched "
            f"({100 * len(hits) / len(prot):.1f}%), {int((hits['call'] == 1).sum()):>4} essential")

    out = pd.DataFrame({"uniprot_ac": prot["uniprot_ac"]})
    if per_strain:
        wide = pd.DataFrame(per_strain).reindex(out["uniprot_ac"])
        out["n_strains_measured"] = wide.notna().sum(axis=1).to_numpy()
        out["n_strains_essential"] = (wide == 1).sum(axis=1).to_numpy()
        out["strains_measured"] = [
            ";".join(sorted(c for c in wide.columns if pd.notna(r[c]))) for _, r in wide.iterrows()]
    else:
        out["n_strains_measured"] = 0
        out["n_strains_essential"] = 0
        out["strains_measured"] = ""
    return M.reindex(out[COLUMNS], species)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SAME_SPECIES_SCREENS),
                    choices=list(SAME_SPECIES_SCREENS))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("essentiality/strain_homologs.py - measured counterparts in screened strains")
    rule("=")
    say(f"  rule     : >= {L.DIRECT_PIDENT:g}% identity, >= {L.MIN_SCOV:g}% of the SUBJECT aligned")
    say("  output   : a COUNT of measurements, never a transferred label")
    rule()

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    for sp in args.species:
        say(f"  {sp}")
        df = build(sp, args.threads)
        path = EVIDENCE_DIR / f"strain_homologs_{sp}.tsv"
        df.to_csv(path, sep="\t", index=False)
        n = int((df.n_strains_measured > 0).sum())
        say(f"    -> {path.relative_to(REPO_ROOT)}  {len(df):,} rows, "
            f"{n:,} with a measured counterpart ({100 * n / len(df):.1f}%)"
            + (f", {int((df.n_strains_measured >= 2).sum()):,} in 2+ strains" if n else ""))
    rule("=")


if __name__ == "__main__":
    main()
