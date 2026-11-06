"""Fetch published essentiality screens, and record exactly what failed and how to get it by hand.

A screen we cannot download is not a dead end -- it is a TASK FOR A HUMAN, and this script's real
output is that task list. Every dataset ends in `evidence/screen_fetch_status.tsv` with a status, the
routes tried, and -- when automation fails -- a `manual` field saying precisely what to click and
where to put the file. v1 lost real datasets to fetch failures that were never written down; the
423-gene Ramage set survived only because a later paper happened to re-tabulate it.

FETCH LADDER, cheapest first:
  1. direct URL (publisher CDN, or a lab/Zenodo/GitHub mirror)
  2. Europe PMC supplementary ZIP -- `.../{PMCID}/supplementaryFiles`, which bypasses the
     Cloudflare 403s that ASM, PNAS, Nature and Cell return to non-browser clients
  3. NCBI OA tarball for anything in PMC Open Access
  4. give up, record `manual` instructions, keep going

Nothing here overwrites an existing file; a dataset already on disk is reported `cached` and skipped.

Run with the `gradi` env:
    python scripts/essentiality/fetch_screens.py --list
    python scripts/essentiality/fetch_screens.py --dry-run
    python scripts/essentiality/fetch_screens.py --only cain2017_njst258
"""

from __future__ import annotations

import argparse
import io
import sys
import time
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

RAW = REPO_ROOT / "data" / "raw"
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "essentiality" / "evidence"
UA = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                    "(KHTML, like Gecko) Chrome/125.0 Safari/537.36"}
TIMEOUT = 120
MIN_BYTES = 2000        # anything smaller is an error page, not a supplement

# An HTTP 200 is NOT evidence of data -- the house rule, and it bit this script on its first run:
# a Nature article URL returned 407 KB of HTML and passed a size check as if it were a supplement.
# A payload must look like a data file by content-type AND not start like markup.
DATA_TYPES = ("application/vnd.ms-excel", "application/vnd.openxmlformats", "text/csv",
              "text/tab-separated-values", "application/zip", "application/gzip",
              "application/octet-stream", "text/plain", "application/excel",
              "application/x-excel", "application/msexcel", "application/x-gzip",
              "application/binary", "application/force-download")


def looks_like_data(content: bytes, content_type: str) -> tuple[bool, str]:
    head = content[:512].lstrip().lower()
    if head.startswith((b"<!doctype", b"<html", b"<?xml version=\"1.0\" encoding=\"utf-8\"?><!doctype")):
        return False, "HTML page, not a data file"
    if b"<html" in head:
        return False, "HTML page, not a data file"
    ct = (content_type or "").split(";")[0].strip().lower()
    if ct and not any(ct.startswith(d) for d in DATA_TYPES):
        return False, f"content-type {ct}"
    return True, ct or "unknown"


@dataclass
class Screen:
    key: str
    organism: str
    strain: str
    assay: str
    readout: str                 # "binary" | "continuous" | "both"
    citation: str
    why: str                     # why this screen is worth having
    dest: str                    # directory under data/raw/
    urls: list[str] = field(default_factory=list)
    pmcid: str = ""
    doi: str = ""
    manual: str = ""             # what a human must do if every route fails


# Seeded with the two we already know are missing. The research pass appends to this list; every
# entry carries `manual` BEFORE it is tried, so a failure always leaves usable instructions.
SCREENS: list[Screen] = [
    # RETRACTED ENTRY -- kept as a warning, not as a task.
    #
    # `cain2017_njst258` was in this list on the strength of a note that turned out to be wrong in
    # BOTH respects: the DOI I had (10.1038/srep44322) belongs to an unrelated neurovascular paper,
    # and NO NJST258 TraDIS screen was ever published. NJST258_1/2 exist only as assemblies. Amy
    # Cain's Klebsiella library is RH201207 (Sci Rep 2017, srep42483), whose sole supplement is a
    # 6-page PDF carrying no per-gene table at all -- so it is not a data source either.
    #
    # Do not re-add either without checking the DOI resolves to the paper you think it does.

    Screen(
        key="nichols2011_chemgen",
        organism="Escherichia coli", strain="K-12 BW25113", assay="chemical genomics",
        readout="continuous",
        citation="Nichols RJ et al. 2011, Cell 144:143-156 (phenotypic landscape of a bacterial cell)",
        why="~4,000 gene x 324 condition fitness matrix -- by far the widest conditional-fitness "
            "resource for E. coli, and the only one that would let an endpoint ask 'essential under "
            "WHICH stress'. v1 failed to fetch it (PMC serves the supplement behind an image "
            "CAPTCHA) and judged it redundant with RB-TnSeq; that judgement was never tested.",
        dest="ecoli/essentiality/nichols2011_chemgen",
        doi="10.1016/j.cell.2010.11.052", pmcid="PMC3060659",
        # The Elsevier CDN is UNPROTECTED; cell.com and the PMC landing page are the gated parts.
        # The PII is S0092867410013747 -- not derivable from the DOI, which is why this sat
        # "unobtainable" since v1. Verified: 13,938,077 bytes, sheet TableS2-FinalData,
        # 3,981 rows x 325 columns of continuous S-scores keyed ECK####-GENENAME.
        urls=["https://ars.els-cdn.com/content/image/1-s2.0-S0092867410013747-mmc2.xls"],
        manual="Should not be needed any more -- the ars.els-cdn.com URL above works. If it ever "
               "stops: do NOT use the PMC route, which serves a reCAPTCHA page as HTTP 200 (21 KB "
               "of HTML that passes a naive size check). Open the article in a browser instead.",
    ),
    Screen(
        key="fitnessbrowser_2024_db",
        organism="48 bacteria + archaea (incl. E. coli BW25113; NO K. pneumoniae)",
        strain="various", assay="RB-TnSeq", readout="continuous",
        citation="Price MN, Deutschbauer AM et al. -- February 2024 release of the Fitness Browser "
                 "(figshare 10.6084/m9.figshare.25236931), 7,552 genome-wide fitness experiments",
        why="The single widest CONTINUOUS fitness resource in existence, and CC BY 4.0. DEG gives "
            "us 49 binary endpoints; this gives graded fitness over 7,552 experiments AND ships "
            "`aaseqs` -- protein sequences for every organism -- so it is directly embeddable with "
            "no identifier join at all. Contains `Keio` = E. coli BW25113 (168 experiments, 3,789 "
            "genes), the strain most of our E. coli screens use. v1 held only a 3,789-gene summary "
            "and a 280-condition slice of it; the full matrix was never retained.\n\n"
            "TWO THINGS IT IS NOT. (1) There is NO K. pneumoniae here. The one Klebsiella, `Koxy`, "
            "is `Klebsiella | michiganensis | M5al` per the database's own Organism table (taxid "
            "290337) and its species assignment is unsettled across sources -- a comparator, never "
            "a Kp label. (2) RB-TnSeq CANNOT see truly essential genes: a gene with no surviving "
            "insertions has no fitness value, so it is ABSENT rather than extreme. This is a "
            "FITNESS resource, and an essentiality endpoint built from it would be measuring the "
            "wrong thing.",
        dest="other/essentiality/fitness_browser_2024",
        doi="10.6084/m9.figshare.25236931",
        urls=["https://ndownloader.figshare.com/files/44580544",   # aaseqs.gz, 43 MB
              "https://ndownloader.figshare.com/files/44580529",   # db.StrainFitness.Keio.gz
              "https://ndownloader.figshare.com/files/44580535"],  # db.StrainFitness.Koxy.gz
        manual="figshare article 25236931. If ndownloader links fail, open "
               "https://figshare.com/articles/dataset/25236931 in a browser and take `aaseqs.gz`, "
               "`feba.db.gz` (2.3 GB, the whole sqlite database) and the per-organism "
               "`db.StrainFitness.*.gz` files you want. Put them in "
               "data/raw/other/essentiality/fitness_browser_2024/.",
    ),
    Screen(
        key="choe2025_ecoli_tnseq",
        organism="Escherichia coli", strain="K-12 BW25113", assay="Tn-seq", readout="continuous",
        citation="Choe D et al. 2025, iScience -- E. coli essentiality across 13 conditions",
        why="4,198 genes x 13 conditions with a continuous `ER` (essentiality ratio), and it joins "
            "to our E. coli proteome at 99.1% on b-numbers. The widest CONDITION-RESOLVED E. coli "
            "essentiality table we have found that is not chemical genomics -- it answers "
            "'essential in which medium', which every binary screen we hold cannot.",
        dest="ecoli/essentiality/choe2025_ecoli",
        pmcid="PMC12063145", doi="10.1016/j.isci.2025.112438",
        manual="Europe PMC ZIP for PMC12063145 is the working route (mmc3.xlsx). NOTE the file has "
               "a TWO-ROW header: row 1 is `Name, b number, ..., ER`, row 2 carries the 13 "
               "condition names -- read with header=[0,1] or the columns come out wrong.",
    ),
    Screen(
        key="choe2023_ecoli_tnseq",
        organism="Escherichia coli", strain="K-12 BW25113", assay="Tn-seq", readout="continuous",
        citation="Choe D et al. 2023, mSystems -- E. coli essentiality in LB and M9",
        why="4,498 genes with `IPKM` (insertions per kb per million) in TWO media, LB and M9 "
            "glucose. Two media is the minimum needed to separate 'essential' from 'essential in "
            "rich medium', which is the confound sitting under every single-condition screen we "
            "hold.",
        dest="ecoli/essentiality/choe2023_ecoli",
        pmcid="PMC9948719", doi="10.1128/msystems.01011-22",
        manual="Europe PMC ZIP for PMC9948719, Table S1. Carries BOTH media as separate column "
               "blocks (`Insertion, IPKM, ec Insertion, ecIPKM, Essentiality` per medium).",
    ),
    Screen(
        key="bruchmann2021_kp",
        organism="Klebsiella pneumoniae", strain="RH201207 + ATCC 43816", assay="TraDIS",
        readout="both",
        citation="Bruchmann S et al. 2021, Nucleic Acids Res -- Kp TraDIS, two strains",
        why="TWO Klebsiella genetic backgrounds we hold nothing for, each with a three-state "
            "essentiality call AND logFC/q, plus a built-in cross-strain locus map and COG letters. "
            "TableS5 = RH201207 (5,390 rows, `KPNRH_` on FR997879), TableS6 = ATCC 43816 (5,217, "
            "`VK055_RS` on NZ_CP009208.1). Plain CSV inside a nested zip.",
        dest="kpneumoniae/essentiality/bruchmann2021",
        pmcid="PMC7981267",
        manual="Europe PMC zip for PMC7981267, then unzip the INNER "
               "`ftab009_supplemental_files.zip` to reach TableS5.csv / TableS6.csv.",
    ),
    Screen(
        key="short2020_kp4",
        organism="Klebsiella pneumoniae", strain="B5055, NTUH-K2044, ATCC 43816, RH201207",
        assay="TraDIS", readout="both",
        citation="Short FL et al. 2020, Infect Immun -- Kp TraDIS across four strains",
        why="FOUR complete Klebsiella proteomes with three paired contrasts each -- the widest "
            "single Kp resource found. Sheets are per strain: B5055 (`BN49_` on FO834906), "
            "NTUH-K2044 (`KP1_` on AP006725), ATCC 43816 (`VK055_` on CP009208 GenBank, NOT the "
            "_RS form), RH201207 (`RH201207_` on LT216436).",
        dest="kpneumoniae/essentiality/short2020",
        pmcid="PMC7375759",
        manual="Europe PMC zip for PMC7375759 -> IAI.00043-20-s0002.xlsx.",
    ),
    Screen(
        key="gray2024_bridge",
        organism="Klebsiella pneumoniae", strain="ECL8 / KPNIH1 / RH201207 / ATCC 43816",
        assay="TraDIS", readout="binary",
        citation="Gray J, Eichelberger K, Short FL et al. 2024, eLife 88971",
        why="The FOUR-STRAIN BRIDGE (385 / 642 / 608 / 476 essential) that lets the Kp screens be "
            "stacked, plus the raw ECL8 three-state table. Read the caveats: the paper's headline "
            "'427 essential' for ECL8 is a POST-CURATION list, while the raw call is 373/4551/241 "
            "-- different labels. And its RH201207 column is `KPNRH_` on FR997879, which is a "
            "DIFFERENT ASSEMBLY from Jana's `RH201207_` on LT216436 (chromosomes differ by 901 bp). "
            "Do not merge them.",
        dest="kpneumoniae/essentiality/gray2024_elife",
        urls=["https://cdn.elifesciences.org/articles/88971/elife-88971-fig3-data1-v1.xlsx"],
        manual="Direct eLife CDN; unprotected.",
    ),
    Screen(
        key="rome2026_st258",
        organism="Klebsiella pneumoniae", strain="NJST258_2 (ST258)", assay="TnSeq/TRANSIT",
        readout="both",
        citation="Rome K et al. 2026, Antimicrob Agents Chemother -- ST258 TnSeq",
        why="A SECOND ST258 background beside KPNIH1, with TRANSIT ES/GD/NE/GA calls and log2FC "
            "(5,125 chromosome genes). Tags `KPNJ2_RS` on GCF_000597905.1. CAVEAT: grown in "
            "iron-depleted medium, so its essential set is condition-shifted relative to rich-medium "
            "screens -- a difference to model, not to average away.",
        dest="kpneumoniae/essentiality/rome2026",
        pmcid="PMC13436387",
        manual="Europe PMC zip for PMC13436387 -> aac.00197-26-s0004.xlsx.",
    ),
]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def _write(dest: Path, name: str, content: bytes) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    p = dest / name
    p.write_bytes(content)
    return p


def try_direct(s: Screen, dest: Path) -> tuple[bool, str, list[str]]:
    got = []
    for url in s.urls:
        try:
            r = requests.get(url, headers=UA, timeout=TIMEOUT, allow_redirects=True)
        except Exception as exc:                                     # noqa: BLE001
            say(f"      direct {url[:70]} -> {type(exc).__name__}")
            continue
        if r.status_code != 200:
            say(f"      direct {url[:70]} -> HTTP {r.status_code}")
            continue
        if len(r.content) < MIN_BYTES:
            say(f"      direct {url[:70]} -> {len(r.content)} bytes (an error page, not data)")
            continue
        ok, why = looks_like_data(r.content, r.headers.get("Content-Type", ""))
        if not ok:
            say(f"      direct {url[:70]} -> REJECTED: {why}")
            continue
        name = url.split("/")[-1].split("?")[0].split("#")[0] or "download.bin"
        got.append(_write(dest, name, r.content).name)
    return bool(got), "direct", got


def try_curl(s: Screen, dest: Path) -> tuple[bool, str, list[str]]:
    """curl -L, for hosts `requests` cannot follow cleanly.

    figshare is the case that forced this: `ndownloader` answers 202 while it prepares the file,
    then 302s to an S3 URL signed with `X-Amz-Expires=10`. Ten seconds is long enough for curl's
    single-process redirect follow and short enough that a two-step client loses the race.
    """
    import subprocess
    got = []
    for url in s.urls:
        name = url.rstrip("/").split("/")[-1].split("?")[0] or "download.bin"
        dest.mkdir(parents=True, exist_ok=True)
        tmp = dest / f".{name}.part"
        # NO --retry here: it makes curl re-issue the request, and figshare's redirect target is
        # an S3 URL signed with X-Amz-Expires=10, so a retry races its own signature and lands on
        # 202/403. A plain single-pass -L follows it inside the window. Retry the WHOLE call
        # instead, below.
        # TWO User-Agents, short one FIRST, and the order is not arbitrary -- measured on figshare:
        #     full Chrome UA  -> HTTP 202, 0 bytes, every time
        #     "Mozilla/5.0"   -> HTTP 200, the file
        #     no UA at all    -> HTTP 200, the file
        # i.e. the browser-like header is what BREAKS it, the opposite of the usual publisher case
        # where a bare client gets a 403. Publishers and repositories want opposite things, so try
        # both rather than assuming either.
        code = "?"
        for attempt, agent in enumerate(["Mozilla/5.0", UA["User-Agent"], "Mozilla/5.0"]):
            r = subprocess.run(["curl", "-sL", "-A", agent, "-o", str(tmp),
                                "-w", "%{http_code}", url], capture_output=True, text=True)
            code = (r.stdout or "?").strip()
            if code == "200" and tmp.exists() and tmp.stat().st_size >= MIN_BYTES:
                break
            time.sleep(2 * (attempt + 1))
        if code != "200" or not tmp.exists() or tmp.stat().st_size < MIN_BYTES:
            say(f"      curl {url[:60]} -> HTTP {code}, "
                f"{tmp.stat().st_size if tmp.exists() else 0} bytes")
            tmp.unlink(missing_ok=True)
            continue
        ok, why = looks_like_data(tmp.read_bytes()[:512], "")
        if not ok:
            say(f"      curl {url[:60]} -> REJECTED: {why}")
            tmp.unlink(missing_ok=True)
            continue
        # figshare serves the real filename only in the redirect; recover it when we can
        final = dest / name
        tmp.rename(final)
        say(f"      curl {final.name}  {final.stat().st_size / 1e6:.1f} MB")
        got.append(final.name)
    return bool(got), "curl", got


def try_europepmc(s: Screen, dest: Path) -> tuple[bool, str, list[str]]:
    """Europe PMC ships the publisher's own supplementary files as one ZIP.

    This is the route that bypasses the Cloudflare 403s ASM, PNAS, Nature and Cell return to
    non-browser clients -- documented in CLAUDE.md under the essentiality corpus.
    """
    if not s.pmcid:
        return False, "europepmc", []
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{s.pmcid}/supplementaryFiles"
    try:
        r = requests.get(url, headers=UA, timeout=TIMEOUT)
    except Exception as exc:                                         # noqa: BLE001
        say(f"      europepmc -> {type(exc).__name__}")
        return False, "europepmc", []
    if r.status_code != 200 or len(r.content) < MIN_BYTES:
        say(f"      europepmc -> HTTP {r.status_code}, {len(r.content)} bytes")
        return False, "europepmc", []
    try:
        z = zipfile.ZipFile(io.BytesIO(r.content))
    except zipfile.BadZipFile:
        say("      europepmc -> not a zip")
        return False, "europepmc", []
    got = []
    for n in z.namelist():
        data = z.read(n)
        if len(data) >= MIN_BYTES:
            got.append(_write(dest, Path(n).name, data).name)
    return bool(got), "europepmc", got


def fetch(s: Screen, dry: bool) -> dict:
    dest = RAW / s.dest
    # A PLACEHOLDER left by a previous failed fetch is NOT data -- v1 left several, and counting
    # them as "cached" is how a missing dataset stays missing silently.
    DATA_EXT = {".xlsx", ".xls", ".csv", ".tsv", ".txt", ".zip", ".gz", ".dat", ".faa", ".fa"}
    existing = ([p.name for p in dest.glob("*")
                 if p.is_file() and p.suffix.lower() in DATA_EXT
                 and not p.name.upper().startswith("PLACEHOLDER")]
                if dest.exists() else [])
    say(f"  {s.key:24s} {s.organism} {s.strain}")
    say(f"    {s.assay}, readout={s.readout}")
    if existing:
        say(f"    cached: {len(existing)} file(s) already present -- skipping")
        return dict(screen=s.key, status="cached", route="-", n_files=len(existing),
                    files=";".join(existing[:6]), manual="", **_meta(s))
    if dry:
        say(f"    would try: {len(s.urls)} direct URL(s)"
            + (f" then Europe PMC {s.pmcid}" if s.pmcid else ""))
        return dict(screen=s.key, status="dry-run", route="-", n_files=0, files="",
                    manual=s.manual, **_meta(s))

    for fn in (try_direct, try_curl, try_europepmc):
        ok, route, got = fn(s, dest)
        if ok:
            say(f"    OK via {route}: {len(got)} file(s) -> {dest.relative_to(REPO_ROOT)}")
            return dict(screen=s.key, status="ok", route=route, n_files=len(got),
                        files=";".join(got[:6]), manual="", **_meta(s))

    say("    FAILED every automated route -- manual instructions recorded")
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "PLACEHOLDER.md").write_text(
        f"# {s.key} -- NOT DOWNLOADED\n\n{s.citation}\n\nDOI: {s.doi}\nPMCID: {s.pmcid}\n\n"
        f"## Why we want it\n\n{s.why}\n\n## How to get it by hand\n\n{s.manual}\n")
    return dict(screen=s.key, status="MANUAL", route="-", n_files=0, files="",
                manual=s.manual, **_meta(s))


def _meta(s: Screen) -> dict:
    return dict(organism=s.organism, strain=s.strain, assay=s.assay, readout=s.readout,
                citation=s.citation, doi=s.doi, pmcid=s.pmcid, why=s.why,
                dest=f"data/raw/{s.dest}")


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", nargs="+", help="fetch just these keys")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("scripts/essentiality/fetch_screens.py -- published screens, and what needs a human")
    rule("=")
    say(f"  out      {EVIDENCE_DIR.relative_to(REPO_ROOT)}/screen_fetch_status.tsv")
    say("  ladder   direct URL -> Europe PMC supplementaryFiles -> record MANUAL instructions")

    todo = [s for s in SCREENS if not a.only or s.key in a.only]
    if a.list:
        rule()
        for s in todo:
            say(f"  {s.key:24s} {s.organism:24s} {s.assay:18s} {s.readout}")
            say(f"      {s.citation}")
        return

    rule()
    say("FETCH")
    rule()
    rows = [fetch(s, a.dry_run) for s in todo]
    say("")

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df["checked_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    out = EVIDENCE_DIR / "screen_fetch_status.tsv"
    if out.exists() and not a.dry_run:
        prev = pd.read_csv(out, sep="\t")
        df = pd.concat([prev[~prev["screen"].isin(df["screen"])], df])
    if not a.dry_run:
        df.to_csv(out, sep="\t", index=False)

    rule()
    say("STATUS")
    rule()
    for _, r in df.iterrows():
        say(f"  {r['status']:8s} {r['screen']:24s} {r['organism']}")
    n_manual = int((df["status"] == "MANUAL").sum())
    if n_manual:
        rule()
        say(f"{n_manual} SCREEN(S) NEED A HUMAN -- instructions in "
            f"{out.relative_to(REPO_ROOT)} and in each PLACEHOLDER.md:")
        for _, r in df[df["status"] == "MANUAL"].iterrows():
            say(f"\n  {r['screen']}  ({r['citation']})")
            say(f"    {r['manual']}")
    rule("=")


if __name__ == "__main__":
    main()
