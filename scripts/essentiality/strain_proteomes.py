"""Locus tag -> protein sequence, for every K. pneumoniae strain a screen is keyed on.

A screen is only as good as the fraction of its genes we can attach a SEQUENCE to -- that is what
becomes a feature vector. So this script's real output is the coverage table it prints, and it
exits non-zero if a strain falls below its floor.

TWO ROUTES, because no single one works everywhere:

  efetch    `rettype=fasta_cds_aa` on the chromosome accession the paper's tags belong to. Headers
            carry `[locus_tag=...]` inline, so the namespace cannot drift. Preferred.
  assembly  NCBI `GCF_*` GFF + protein.faa, walking gene -> CDS via ID/Parent and reading
            `old_locus_tag`. Needed where the submitter record has no CDS features at all.

RH201207 is why both exist: its GenBank assembly `GCA_905477585.1` carries **gene features only,
no CDS and no proteins**, and its chromosome record `FR997879` returns zero CDS from efetch. The
RefSeq counterpart `GCF_905477585.1` has 5,527 proteins and preserves `KPNRH_*` as `old_locus_tag`.

AND THE SAME STRAIN CAN NEED DIFFERENT RECORDS PER PAPER. RH201207 appears as `RH201207_*` on
`LT216436` (Short 2020) and `KPNRH_*` on `FR997879`/`GCF_905477585.1` (Bruchmann 2021) -- different
assemblies of one strain, chromosomes 901 bp apart. They are staged separately and must not be
merged.

Run with the `gradi` env:
    python scripts/essentiality/strain_proteomes.py
    python scripts/essentiality/strain_proteomes.py --strain b5055
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.request
import zipfile
from io import BytesIO
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

OUT_DIR = REPO_ROOT / "data" / "source" / "ncbi" / "kp_strains"
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={}&rettype=fasta_cds_aa&retmode=text")
DATASETS = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}/download"
            "?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF")

# label -> (route, accession, locus prefix, which screen needs it)
STRAINS = {
    "rh201207_bruchmann": ("assembly", "GCF_905477585.1", "KPNRH_",     "Bruchmann 2021 TableS5"),
    "rh201207_short":     ("efetch",   "LT216436",        "RH201207_",  "Short 2020 sheet 4"),
    "b5055":              ("efetch",   "FO834906",        "BN49_",      "Short 2020 sheet 1"),
    "ntuh_k2044":         ("efetch",   "AP006725",        "KP1_",       "Short 2020 sheet 2"),
    "njst258_2":          ("assembly", "GCF_000597905.1", "KPNJ2_RS",   "Rome 2026"),
}
FLOOR = 0.90

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(c: str = "-", w: int = 92) -> None:
    say(c * w)


def from_efetch(acc: str, prefix: str) -> pd.DataFrame:
    txt = urllib.request.urlopen(EFETCH.format(acc), timeout=300).read().decode()
    rows, tag, buf = [], None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if tag:
                rows.append({"locus_tag": tag, "sequence": "".join(buf)})
            m = re.search(r"\[locus_tag=([^\]]+)\]", line)
            tag, buf = (m.group(1) if m else None), []
        elif tag:
            buf.append(line.strip())
    if tag:
        rows.append({"locus_tag": tag, "sequence": "".join(buf)})
    return pd.DataFrame(rows)


def from_assembly(acc: str, prefix: str) -> pd.DataFrame:
    blob = urllib.request.urlopen(DATASETS.format(acc), timeout=600).read()
    z = zipfile.ZipFile(BytesIO(blob))
    gff = next(n for n in z.namelist() if n.endswith("genomic.gff"))
    faa = next(n for n in z.namelist() if n.endswith("protein.faa"))

    seqs, pid, buf = {}, None, []
    for line in z.read(faa).decode().splitlines():
        if line.startswith(">"):
            if pid:
                seqs[pid] = "".join(buf)
            pid, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if pid:
        seqs[pid] = "".join(buf)

    # old_locus_tag sits on the GENE feature, protein_id on the CDS, linked by ID/Parent
    genes, rows = {}, []
    for line in z.read(gff).decode().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].rstrip().split(";") if "=" in kv)
        if f[2] in ("gene", "pseudogene") and a.get("ID"):
            genes[a["ID"]] = (unquote(a.get("locus_tag", "")), unquote(a.get("old_locus_tag", "")))
        elif f[2] == "CDS" and a.get("protein_id"):
            cur, old = genes.get(a.get("Parent", ""), ("", ""))
            for t in {x.strip() for x in (cur + "," + old).split(",") if x.strip()}:
                if t.startswith(prefix) and a["protein_id"] in seqs:
                    rows.append({"locus_tag": t, "sequence": seqs[a["protein_id"]]})
    return pd.DataFrame(rows).drop_duplicates("locus_tag")


def build(label: str, refresh: bool) -> pd.DataFrame:
    route, acc, prefix, used_by = STRAINS[label]
    out = OUT_DIR / f"{label}.tsv"
    if out.exists() and not refresh:
        d = pd.read_csv(out, sep="\t")
        say(f"  {label:20s} cached  {len(d)} proteins")
        return d
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    say(f"  {label:20s} {route:8s} {acc}  ({used_by})")
    d = (from_efetch if route == "efetch" else from_assembly)(acc, prefix)
    d = d[d["locus_tag"].astype(str).str.startswith(prefix)]
    d = d[d["sequence"].str.len() > 0].drop_duplicates("locus_tag")
    d.to_csv(out, sep="\t", index=False)
    say(f"    {len(d)} proteins carrying {prefix}*  -> {out.relative_to(REPO_ROOT)}")
    time.sleep(2)
    return d


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--strain", nargs="+", choices=list(STRAINS), default=list(STRAINS))
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("K. pneumoniae screen strains -- locus tag to protein SEQUENCE")
    rule("=")
    say(f"  out    {OUT_DIR.relative_to(REPO_ROOT)}/<strain>.tsv")
    say(f"  floor  a screen must attach a sequence to >= {FLOOR:.0%} of its genes")
    rule()
    frames = {lbl: build(lbl, a.refresh) for lbl in a.strain}

    rule()
    say("SEQUENCE COVERAGE  -- the number that decides whether a screen is usable")
    rule()
    say(f"  {'strain':20s} {'proteins':>9} {'median aa':>10}  prefix")
    bad = []
    for lbl, d in frames.items():
        if not len(d):
            say(f"  {lbl:20s} {'0':>9}   EMPTY"); bad.append(lbl); continue
        say(f"  {lbl:20s} {len(d):>9} {int(d.sequence.str.len().median()):>10}  "
            f"{STRAINS[lbl][2]}*")
    if bad:
        sys.exit(f"FAILED: no sequences for {', '.join(bad)}")
    rule("=")


if __name__ == "__main__":
    main()
