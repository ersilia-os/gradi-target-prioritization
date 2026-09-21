"""Published essentiality screens -> one tidy training set per endpoint.

Each screen becomes its OWN endpoint. They are never merged, for the same reason stage 04 never
merges its two ClpP activators: these measure different things on different organisms with
different assays, and the base rate alone spans 0.061-0.111 here (15x across the wider DEG corpus).
A merged label would be a number with no referent.

WHY THIS SCRIPT EXISTS: v1 had all of this data and joined it BY GENE SYMBOL, losing most of it --
ECL8 954/5,165, Ramage 212/424, KPPR1 32 of 3,791. Its own retrospective calls that "the single
most expensive recurring cost in v1", and v2's stage 00 exists to fix it. So every join here is by
a STABLE IDENTIFIER (b-number) or by sequence, never by symbol -- with one measured exception,
noted below.

    THE ONE PLACE A SYMBOL JOIN IS SAFE. Gene-name coverage is 100% on E. coli and 18.4% -> 63.4%
    on K. pneumoniae. The ~27% loss trap 59 documents is a KLEBSIELLA problem. Goodall 2018 ships
    gene symbols and nothing else, and on E. coli that costs essentially nothing -- which this
    script proves every run by reporting the join rate rather than assuming it.

Six endpoints:

    keio_ess       E. coli K-12, PEC/Keio arrayed single-gene knockout      b-number
    goodall_ess    E. coli BW25113, TraDIS                                  symbol -> b-number
    bw25113_ess    E. coli BW25113, TraDIS reprocessed through DESeq        BW25113_ -> b-number
    conservation   how conserved is essentiality across Gammaproteobacteria b-number
    bn373_ess      K. pneumoniae ECL8, TraDIS/DESeq                         BN373_ old_locus_tag
    kpnih1_ess     K. pneumoniae KPNIH1, Tn-seq                             KPNIH1_ old_locus_tag

The two Kp endpoints need their strain proteome from tier D (`download.py --tier D`) and an
ESM-C run (`esmc.py --strain`). `--list` shows what is built and what was measured and dropped.

Run with the `gradi` env:
    python scripts/essentiality/screens.py
    python scripts/essentiality/screens.py --screen keio_ess --dry-run
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

RAW = REPO_ROOT / "data" / "raw"
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"

COMPENDIUM = RAW / "other" / "essentiality" / "enterobacteriaceae_tradis" / "giant-tab_final.tsv"
KEIO_B = "Locus: Escherichia coli BW25113 (Keio)"      # b-numbers
BW_LOCUS = "Locus: Escherichia coli BW25113"            # BW25113_#### tags

# The clade whose conservation label we ship. Conservation is strongly BIMODAL -- a gene is
# essential nearly everywhere or nearly nowhere -- so the middle class is thin at narrow clades and
# fattens as the clade widens. Measured over the 4,301 b-numbered rows at 80/20 cuts:
#     Enterobacteriaceae   core 264  mid 100  non 3937     <- 100 is too thin to learn
#     Gammaproteobacteria  core 205  mid 530  non 3566     <- chosen
#     Bacteria             core 133  mid 557  non 3611
# The raw 0-100 value ships beside the 3-class call so the binning can change without a re-run.
CLADE_COLUMN = "Gammaproteobacteria %essential"
CORE_CUT, NON_CUT = 80.0, 20.0

# Join-rate floors. Below these something is structurally wrong (a renamed column, a changed
# identifier namespace) and the run must stop rather than ship a quietly decimated training set.
# Set from v1's measured symbol-join rates: whatever we do must BEAT what v1 already achieved.
FLOORS = {"keio_ess": 0.95, "goodall_ess": 0.90, "bw25113_ess": 0.90,
          "conservation": 0.90, "bn373_ess": 0.90, "kpnih1_ess": 0.99}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# --------------------------------------------------------------------------- the E. coli key maps

def ecoli_bnumber_map() -> dict[str, str]:
    """b-number -> uniprot_ac. E. coli locus tags ARE b-numbers (4,402 of 4,403)."""
    lt = P.load_locus_tags("ecoli")
    m = lt[lt["locus_tag"].astype(str).str.match(r"^b\d{4}$", na=False)]
    return dict(zip(m["locus_tag"], m["uniprot_ac"]))


def ecoli_symbol_map() -> dict[str, str]:
    """normalised gene symbol -> uniprot_ac, from EVERY name E. coli carries, not just the primary.

    `gene_name` holds one preferred name; `gene_name_all` holds the alternatives. A screen may use
    either, so both are indexed -- the preferred name wins a collision.
    """
    df = P.load("ecoli")
    out: dict[str, str] = {}
    for col in ("gene_name_all", "gene_name"):          # primary LAST so it overwrites
        if col not in df.columns:
            continue
        for ac, names in zip(df["uniprot_ac"], df[col].fillna("")):
            for n in re.split(r"[;,\s]+", str(names)):
                n = n.strip().lower()
                if n:
                    out[n] = ac
    return out


def compendium(usecols: list[str] | None = None) -> pd.DataFrame:
    if not COMPENDIUM.exists():
        sys.exit(f"FATAL missing {COMPENDIUM}. It is the Enterobacteriaceae-TraDIS compendium "
                 "(PMID 39207104), staged by v1 under data/raw/other/essentiality/.")
    return pd.read_csv(COMPENDIUM, sep="\t", low_memory=False, usecols=usecols)


# --------------------------------------------------------------------------------- the screens

def screen_keio() -> tuple[pd.DataFrame, str]:
    """PEC's curation of the Keio arrayed knockout collection. The only non-transposon assay here.

    Column 9 is `Class(1:essential 2:noessential 3:unknown)`; column 3 packs synonyms including the
    b-number (`b0001,ECK0001,JW4367`). Unknown (n=5) is DROPPED, never made negative.
    """
    p = RAW / "ecoli" / "essentiality" / "pec" / "PECData.dat"
    d = pd.read_csv(p, sep="\t", header=None, low_memory=False, skiprows=1,
                    names=[f"c{i}" for i in range(13)])
    d = d[pd.to_numeric(d["c1"], errors="coerce") == 1]          # feature type 1 = protein-coding
    cls = pd.to_numeric(d["c9"], errors="coerce")
    d = d.assign(label=np.where(cls == 1, 1, np.where(cls == 2, 0, np.nan)))
    d["bnum"] = d["c3"].astype(str).str.extract(r"\b(b\d{4})\b")[0]
    bmap = ecoli_bnumber_map()
    d["uniprot_ac"] = d["bnum"].map(bmap)
    return d[["uniprot_ac", "label", "bnum", "c2"]].rename(columns={"c2": "source_id"}), "b-number"


def screen_goodall() -> tuple[pd.DataFrame, str]:
    """Goodall 2018 TraDIS on BW25113. Ships gene symbols and nothing else.

    `Unclear` (155) is DROPPED. It is the authors saying the insertion profile does not separate --
    treating it as non-essential would manufacture negatives out of ignorance.

    TWO ROUTES, stable one first. The compendium's `Gene_Keio` column is itself a symbol -> b-number
    index built by its authors, so going symbol -> b-number -> accession lands on a STABLE
    identifier; our own proteome's symbol index is the fallback for what it misses. Measured on the
    4,158 labelled rows: compendium alone 96.4%, our symbols alone 93.2%, **union 98.6%** -- so the
    fallback is worth keeping, and the symbol-only join this screen invites is worth avoiding.
    """
    p = RAW / "ecoli" / "essentiality" / "goodall2018_Ec_BW25113" / "mbo001183726st4.xlsx"
    x = pd.ExcelFile(p)
    d = x.parse(x.sheet_names[0], header=1)
    ess = d["Essential"] == True            # noqa: E712 -- real booleans in the sheet
    non = d["Non-essential"] == True        # noqa: E712
    d = d.assign(label=np.where(ess, 1, np.where(non, 0, np.nan)))

    sym = d["Gene"].astype(str).str.strip().str.lower()
    c = compendium(["Gene_Keio", KEIO_B])
    c = c[c[KEIO_B].astype(str).str.match(r"^b\d{4}$", na=False)]
    sym2b = {str(s).strip().lower(): b for s, b in zip(c["Gene_Keio"], c[KEIO_B])
             if str(s) not in ("nan", "-", "")}
    bmap = ecoli_bnumber_map()
    via_b = sym.map(sym2b).map(bmap)
    d["uniprot_ac"] = via_b.fillna(sym.map(ecoli_symbol_map()))
    d["join_route"] = np.where(via_b.notna(), "symbol->b-number", "symbol->proteome")
    return (d[["uniprot_ac", "label", "Gene", "Insertion Index Score", "join_route"]]
            .rename(columns={"Gene": "source_id", "Insertion Index Score": "insertion_index"}),
            "symbol->b-number (+sym)")


def screen_bw25113() -> tuple[pd.DataFrame, str]:
    """The compendium's own DESeq reprocessing of BW25113 -- same pipeline as its other 11 genomes.

    Keyed on `BW25113_####`, bridged to b-numbers through the compendium's two parallel Locus
    columns. `Essentiality` is `Reduced` / `Unchanged`; there is no third class.
    """
    p = COMPENDIUM.parent / "BW25113.out.DESeq.tsv"
    d = pd.read_csv(p, sep="\t")
    bridge = compendium([KEIO_B, BW_LOCUS])
    bridge = bridge[(bridge[KEIO_B].astype(str).str.match(r"^b\d{4}$", na=False))
                    & (bridge[BW_LOCUS].astype(str) != "-")]
    tag2b = dict(zip(bridge[BW_LOCUS].astype(str), bridge[KEIO_B].astype(str)))
    d["bnum"] = d["id"].astype(str).map(tag2b)
    d["uniprot_ac"] = d["bnum"].map(ecoli_bnumber_map())
    d["label"] = np.where(d["Essentiality"].astype(str) == "Reduced", 1,
                          np.where(d["Essentiality"].astype(str) == "Unchanged", 0, np.nan))
    return (d[["uniprot_ac", "label", "id", "log2FoldChange", "padj"]]
            .rename(columns={"id": "source_id"}), "BW25113_ -> b-number")


def screen_conservation() -> tuple[pd.DataFrame, str]:
    """How conserved is this gene's essentiality across the clade -- core / mid / non.

    NOT a per-organism essentiality call: it is the fraction of the clade's genomes in which the
    gene is essential, which is a property of the gene family. That makes it the one endpoint here
    a per-protein embedding is straightforwardly the right feature for -- conservation is a family
    property, and the embedding encodes the family.

    The label is the compendium's own computation over 12 TraDIS genomes; we only bin it. Ships the
    raw percentage too, so the cut can be revisited without re-running anything.
    """
    d = compendium([KEIO_B, CLADE_COLUMN])
    d = d[d[KEIO_B].astype(str).str.match(r"^b\d{4}$", na=False)].copy()
    v = pd.to_numeric(d[CLADE_COLUMN], errors="coerce")
    d["pct_essential"] = v
    d["label"] = np.where(v >= CORE_CUT, 2, np.where(v < NON_CUT, 0, 1))   # 0 non, 1 mid, 2 core
    d.loc[v.isna(), "label"] = np.nan
    d["uniprot_ac"] = d[KEIO_B].astype(str).map(ecoli_bnumber_map())
    return (d[["uniprot_ac", "label", KEIO_B, "pct_essential"]]
            .rename(columns={KEIO_B: "source_id"}), "b-number")


# ------------------------------------------------------------------- K. pneumoniae screens
#
# These key on a STRAIN's locus tags, so they join to that strain's own proteome and are embedded
# there -- no label ever crosses a species boundary. The NCBI annotation puts `old_locus_tag` on
# the *gene* feature and `protein_id` on the *CDS*, linked by ID/Parent, so both must be walked.
#
# WHY `ecl8_ess` IS ABSENT, measured rather than assumed. Eichelberger 2024's own `ecl8_#####` tags
# annotate an assembly that is NOT the deposited one, and all three routes fail:
#     locus tags   475 of 5,165 numeric suffixes shared, and gene names disagree where they overlap
#     coordinates  0 of 5,074 CDS match the deposited assembly's start/end
#     gene symbols the deposited ECL8 annotation carries symbols for just 144 CDS -> 3.7%
# `bn373_ess` is the SAME organism and the SAME TraDIS assay, uniformly reprocessed by the
# compendium and keyed on deposited tags, so it carries the measurement instead. What changes is
# the analysis (DESeq vs the authors' bimodal fit) and the positive count (562 vs 373).

def _ncbi_tag_to_protein(label: str) -> pd.DataFrame:
    """(old_locus_tag, current locus_tag, protein_id) for one NCBI assembly."""
    gff = REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "ncbi" / f"{label}.gff"
    if not gff.exists():
        sys.exit(f"FATAL {gff} missing -- run scripts/proteomes/download.py --tier D --only {label}")
    from urllib.parse import unquote
    genes, cds = {}, []
    for line in gff.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] in ("gene", "pseudogene") and a.get("ID"):
            genes[a["ID"]] = (unquote(a.get("locus_tag", "")), unquote(a.get("old_locus_tag", "")))
        elif f[2] == "CDS" and a.get("protein_id"):
            cds.append({"parent": a.get("Parent", ""), "protein_id": a["protein_id"]})
    d = pd.DataFrame(cds).drop_duplicates("parent")
    d["locus_tag"] = d["parent"].map(lambda p: genes.get(p, ("", ""))[0])
    d["old_locus_tag"] = d["parent"].map(lambda p: genes.get(p, ("", ""))[1])
    return d


def _tag_map(label: str) -> dict[str, str]:
    """Every tag form -> protein_id, so a screen can key on either the old or current namespace.

    USE THE GENBANK (GCA) ASSEMBLY, NOT REFSEQ (GCF). The screens key on the submitter's original
    locus tags, and PGAP's RefSeq re-annotation drops some of them outright. Measured on these two:

        BN373 / ECL8    GenBank 4,930/5,048 (97.7%)   RefSeq 4,767 (94.4%)   RefSeq adds 0
        KPNIH1          GenBank   424/424  (100.0%)   RefSeq   418 (98.6%)   v1 symbol join: 212

    263 of RefSeq's 281 misses were tags ABSENT from its annotation entirely, not proteinless
    genes. GenBank is a strict superset here. CLAUDE.md records the same trap for S. aureus COL
    (GCF retains 60 SACOL tags, GCA retains 2,711) -- it generalises to any strain whose published
    screen predates the current RefSeq annotation.
    """
    d = _ncbi_tag_to_protein(label)
    m: dict[str, str] = {}
    for col in ("locus_tag", "old_locus_tag"):
        for tags, pid in zip(d[col], d["protein_id"]):
            for t in str(tags).split(","):
                t = t.strip()
                if t:
                    m[t] = pid
    return m


ECL8 = "kpneumoniae__ecl8__GCA_000315385.1"
KPNIH1 = "kpneumoniae__kpnih1__GCA_000281535.2"


def screen_bn373() -> tuple[pd.DataFrame, str]:
    """The compendium's DESeq call on K. pneumoniae ECL8 -- THE Klebsiella endpoint.

    `Essentiality` is `Reduced` / `Unchanged`, no third class. Keyed on `BN373_#####`, which is the
    deposited assembly's `old_locus_tag`.
    """
    d = pd.read_csv(COMPENDIUM.parent / "BN373.out.DESeq.tsv", sep="\t")
    d["uniprot_ac"] = d["id"].astype(str).map(_tag_map(ECL8))
    d["label"] = np.where(d["Essentiality"].astype(str) == "Reduced", 1,
                          np.where(d["Essentiality"].astype(str) == "Unchanged", 0, np.nan))
    return (d[["uniprot_ac", "label", "id", "log2FoldChange", "padj"]]
            .rename(columns={"id": "source_id"}), "BN373_ old_locus_tag")


def screen_kpnih1() -> tuple[pd.DataFrame, str]:
    """Ramage 2017's KPNIH1 essential set, as re-tabulated in Jana 2023 s0001.

    A POSITIVES-ONLY LIST: 419 `+` against 5 `-`. It becomes a training set only by the
    construction DEG itself uses -- take the whole proteome, mark the listed genes essential, treat
    the rest as non-essential. That is standard and it is also AN ASSUMPTION WE ADD: these
    negatives are "absent from a published list", which is weaker evidence than a measured
    non-essential call. `label_evidence` records which is which so no consumer can forget.
    """
    x = pd.ExcelFile(RAW / "kpneumoniae" / "essentiality" / "jana2023_crispri" /
                     "aem.00956-23-s0001.xlsx")
    k = x.parse("in vitro screening-KPNIH1", header=1)
    tagcol = [c for c in k.columns if "ocus" in str(c)][0]
    esscol = [c for c in k.columns if "Tnseq" in str(c)][0]
    tmap = _tag_map(KPNIH1)
    listed = {tmap[t] for t in k[tagcol].astype(str).str.strip() if t in tmap
              and str(k.loc[k[tagcol].astype(str).str.strip() == t, esscol].iloc[0]).strip() != "–"}
    n_listed = len(k)
    if len(listed) < 0.95 * n_listed:
        sys.exit(f"FAILED kpnih1_ess: only {len(listed)}/{n_listed} listed genes mapped to a "
                 "protein. v1 managed 212/424 by gene symbol; a tag join must do far better.")
    say(f"    {len(listed)}/{n_listed} listed essential genes mapped "
        f"({100 * len(listed) / n_listed:.1f}%)  -- v1's symbol join reached 212/424")
    allp = _ncbi_tag_to_protein(KPNIH1)["protein_id"].dropna().unique()
    d = pd.DataFrame({"uniprot_ac": allp})
    d["label"] = d["uniprot_ac"].isin(listed).astype(int)
    d["source_id"] = d["uniprot_ac"]
    d["label_evidence"] = np.where(d["label"] == 1, "listed_essential", "absent_from_list")
    return d[["uniprot_ac", "label", "source_id", "label_evidence"]], "KPNIH1_ old_locus_tag"


SCREENS = {
    "keio_ess":     (screen_keio,         "ecoli", "E. coli K-12 (PEC/Keio)", "arrayed knockout"),
    "goodall_ess":  (screen_goodall,      "ecoli", "E. coli BW25113",         "TraDIS"),
    "bw25113_ess":  (screen_bw25113,      "ecoli", "E. coli BW25113",         "TraDIS/DESeq"),
    "conservation": (screen_conservation, "ecoli", "Gammaproteobacteria",     "clade conservation"),
    "bn373_ess":    (screen_bn373,        ECL8,    "Kp ECL8",                 "TraDIS/DESeq"),
    "kpnih1_ess":   (screen_kpnih1,       KPNIH1,  "Kp KPNIH1",               "Tn-seq (list)"),
}

# Measured and dropped -- see the block above. Kept here so nobody re-derives the dead ends.
UNJOINABLE = {
    "ecl8_ess": "Eichelberger 2024 annotates a different ECL8 assembly: tags 475/5165, "
                "coordinates 0/5074, symbols 3.7%. Superseded by bn373_ess (same organism, "
                "same assay, 94.4% join).",
}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--screen", nargs="+", choices=list(SCREENS), default=list(SCREENS))
    ap.add_argument("--list", action="store_true", help="show every screen and exit")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("scripts/essentiality/screens.py -- published screens -> one training set per endpoint")
    rule("=")
    say(f"  in       {COMPENDIUM.parent.relative_to(REPO_ROOT)}/ + data/raw/ecoli/essentiality/")
    say(f"  out      {EVIDENCE_DIR.relative_to(REPO_ROOT)}/screen_<endpoint>.tsv")
    say(f"  clade    {CLADE_COLUMN}  (core >= {CORE_CUT:.0f}, non < {NON_CUT:.0f})")
    say("  join     stable identifiers only -- v1 used gene symbols and lost most of the data")

    if a.list:
        rule()
        say("BUILT")
        for k, (_, sp, org, assay) in SCREENS.items():
            say(f"  {k:14s} {org:24s} {assay}")
        say("\nMEASURED AND DROPPED -- kept here so nobody re-derives the dead ends")
        for k, why in UNJOINABLE.items():
            say(f"  {k:14s} {why}")
        return

    rows, tables = [], {}
    rule()
    say("SCREENS")
    rule()
    for name in a.screen:
        fn, species, org, assay = SCREENS[name]
        df, key = fn()
        n_raw = len(df)
        labelled = df[df["label"].notna()].copy()
        joined = labelled[labelled["uniprot_ac"].notna()].copy()
        # one row per protein: a screen may list paralogs separately; keep the strongest call
        joined = (joined.sort_values("label", ascending=False)
                        .drop_duplicates("uniprot_ac", keep="first"))
        rate = len(joined) / len(labelled) if len(labelled) else 0.0
        vc = joined["label"].value_counts().to_dict()
        say(f"  {name:14s} {org:24s} {assay}")
        say(f"    {'rows in file':<22} {n_raw:6d}   labelled {len(labelled):6d}  "
            f"(dropped {n_raw - len(labelled)} unclear/unknown)")
        say(f"    {'joined by ' + key:<22} {len(joined):6d}   {100 * rate:5.1f}% of labelled")
        say(f"    {'label counts':<22} " + "  ".join(f"{int(k)}={v}" for k, v in sorted(vc.items())))
        if name != "conservation":
            pos = int(joined["label"].sum())
            say(f"    {'base rate':<22} {pos / len(joined):.4f}  ({pos} positive)")
        floor = FLOORS.get(name, 0.0)
        if rate < floor:
            sys.exit(f"FAILED {name}: join rate {rate:.3f} is below the {floor:.2f} floor. "
                     "Something structural changed -- do NOT ship a decimated training set.")
        rows.append({"screen": name, "organism": org, "assay": assay, "join_key": key,
                     "n_file": n_raw, "n_labelled": len(labelled), "n_joined": len(joined),
                     "join_rate": round(rate, 4),
                     "n_positive": int(joined["label"].sum()) if name != "conservation" else pd.NA,
                     "base_rate": round(float(joined["label"].mean()), 4) if name != "conservation" else pd.NA})
        tables[name] = joined
        say("")

    rule()
    say("SHAPES  (X will be n_joined x 1152 for ESM-C/ProteomeLM, x 1024 for ProtT5)")
    rule()
    for r in rows:
        say(f"  {r['screen']:14s} {r['n_joined']:6d} x 1152   base "
            f"{r['base_rate'] if r['base_rate'] is not pd.NA else 'n/a (3-class)'}")

    if a.dry_run:
        rule("=")
        say("dry run: nothing written.")
        return

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    rule()
    say("OUTPUTS")
    rule()
    for name, df in tables.items():
        out = EVIDENCE_DIR / f"screen_{name}.tsv"
        df.to_csv(out, sep="\t", index=False)
        say(f"  {out.relative_to(REPO_ROOT)}  ({len(df)} rows x {df.shape[1]} cols)")
    audit = pd.DataFrame(rows)
    audit["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    audit.to_csv(EVIDENCE_DIR / "screen_join_audit.tsv", sep="\t", index=False)
    say(f"  {(EVIDENCE_DIR / 'screen_join_audit.tsv').relative_to(REPO_ROOT)}")
    rule("=")
    say("screens complete.")
    rule("=")


if __name__ == "__main__":
    main()
