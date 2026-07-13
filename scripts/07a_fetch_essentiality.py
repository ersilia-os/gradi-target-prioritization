"""Fetch the K. pneumoniae / E. coli essentiality source tables (docs §4.1/§4.2).

Robust, resumable, manifest-driven downloader. Most of the primary essentiality data lives in
per-paper supplementary spreadsheets that are gated to non-browser clients (ASM returns 403; PMC
fronts binaries with a reCAPTCHA). We get around that with a candidate-URL ladder, tried in order:

  (1) publisher CDN direct  (eLife / Springer / PLoS — no gate)
  (2) Europe PMC `supplementaryFiles` zip  (bypasses the ASM 403 & PMC reCAPTCHA for OA articles)
  (3) NCBI OA-package `.tar.gz`  (resolved via the oa.fcgi service)
  (4) give up gracefully  -> write a PLACEHOLDER + a row in fetch_status.tsv, so a genuinely gated
      dataset never blocks the rest of the pipeline (downstream treats it as "present, 0 rows").

Downloaded bundles (zip/tar) are opened in-memory and only the spreadsheet members are extracted
into the dataset's raw dir. Everything is idempotent: a dataset whose output dir already holds a
spreadsheet is skipped. ECL8 (Eichelberger/Short 2024) is already on disk under data/raw/legacy/;
it is re-staged (copied) into the organism-first tree rather than re-fetched.

Output:
  data/raw/<organism>/essentiality/<key>/*.xlsx|xls|csv|tsv
  data/raw/<organism>/essentiality/fetch_status.tsv   (dataset, status, route, http, note)

Run with the `gradi` conda env interpreter. Network-only; independent of the compute steps.
"""

from __future__ import annotations

import argparse
import io
import shutil
import sys
import tarfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124 Safari/537.36")
HEADERS = {"User-Agent": UA, "Connection": "close"}
SPREADSHEET_EXTS = (".xlsx", ".xls", ".csv", ".tsv")
EUROPEPMC_SUPP = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
NCBI_OA = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"

# The Gardner-BinfLab Enterobacteriaceae-TraDIS compendium (Goodall/Gardner, PMID 39207104) — the
# spec's recommended §4.2a/§4.2c upgrade. GitHub raw is reliably fetchable headless. `giant-tab_final`
# is a fully ortholog-cluster-aligned essentiality matrix keyed by E. coli Keio b-numbers, carrying
# the curated EcoGene essential call (299 genes) AND graded "Enterobacteriaceae/Bacteria %essential"
# columns spanning E. coli, Klebsiella (incl. ECL8), Salmonella, Citrobacter, Enterobacter.
TRADIS_RAW = "https://raw.githubusercontent.com/Gardner-BinfLab/Enterobacteriaceae-TraDIS/master"
TRADIS_FILES = [
    "data/giant-tab_final.tsv",
    "data/pastml_reconstruction/table_only_ess_genes.tsv",
    "results/Tn-seq/BW25113.out.DESeq.tsv",
    "results/Tn-seq/BN373.out.DESeq.tsv",   # BN373 = K. pneumoniae ECL8 reference annotation
]


@dataclass
class Dataset:
    key: str                       # output subdir name
    organism: str                  # which organism tree it belongs under
    paper: str                     # short citation
    flavor: str                    # essentiality flavor(s) it feeds
    strain: str
    pmcid: str = ""                # drives the Europe PMC + NCBI-OA routes
    direct_urls: list[str] = field(default_factory=list)  # publisher-CDN candidates (route 1)
    open_access: bool = True
    note: str = ""


# --- the manifest (verified routes; see docs/04_essentiality.md + the acquisition research) ---
MANIFEST: list[Dataset] = [
    Dataset("eichelberger2024_ECL8", "kpneumoniae",
            "Eichelberger/Short 2024 (eLife)", "in_vitro_essential;urine;serum", "ECL8",
            pmcid="PMC11349299",
            direct_urls=[f"https://cdn.elifesciences.org/articles/88971/elife-88971-fig{n}-data1-v1.xlsx"
                         for n in (1, 4, 6)],
            note="already downloaded under legacy; re-staged from there"),
    Dataset("goodall2018_Ec_BW25113", "ecoli",
            "Goodall 2018 (mBio)", "in_vitro_essential_Ec", "E. coli BW25113",
            pmcid="PMC5821084", note="NCBI OA tarball route (Europe PMC 404s)"),
    Dataset("bachman2015_KPPR1", "kpneumoniae",
            "Bachman 2015 (mBio)", "in_vivo_lung", "KPPR1",
            pmcid="PMC4462621"),
    Dataset("cain2017_NJST258", "kpneumoniae",
            "Cain 2017 (Sci Rep)", "in_vitro_essential", "NJST258",
            pmcid="PMC5309761"),
    Dataset("mikebachman2023_KPPR1", "kpneumoniae",
            "Mike & Bachman 2023 (PLoS Pathog)", "in_vivo_blood;spleen;liver;lung", "KPPR1",
            pmcid="PMC10381055"),
    Dataset("bachman2025_KPPR1", "kpneumoniae",
            "Bachman 2025 (Nat Commun)", "in_vivo_dissemination", "KPPR1",
            pmcid="PMC11742683",  # open access; Europe PMC supp zip works
            direct_urls=["https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-56095-3"
                         "/MediaObjects/41467_2025_56095_MOESM3_ESM.xlsx"]),
    Dataset("ramage2017_KPNIH1", "kpneumoniae",
            "Ramage 2017 (J Bacteriol)", "in_vitro_essential", "KPNIH1/MKP103",
            pmcid="PMC5637181", open_access=False,
            note="not OA; but the 424-gene KPNIH1 essential set is re-tabulated in Jana 2023 s0001"),
    Dataset("jana2023_crispri", "kpneumoniae",
            "Jana/Zhu 2023 (AEM)", "crispri_library;in_vivo;kpnih1_essential", "Kp Mobile-CRISPRi",
            pmcid="PMC10617577", open_access=False,
            note="ASM 403s non-browser clients; s0001-s0003.xlsx fetched via an authenticated Chrome "
                 "session (evaluate_script fetch with cookies). 870-gene CRISPRi library + in-vivo "
                 "KPPR1 (ratio/p) + KPNIH1 essential (Ramage) + Bachman in-vivo, all in s0001."),
]


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=2, max=30))
def _get(url: str, timeout=(10, 120)) -> requests.Response:
    r = requests.get(url, headers=HEADERS, stream=True, timeout=timeout)
    r.raise_for_status()
    return r


def _resolve_oa_tarball(pmcid: str) -> str | None:
    """Ask the NCBI OA service for the .tar.gz href; rewrite ftp:// -> https://."""
    try:
        r = requests.get(NCBI_OA.format(pmcid=pmcid), headers=HEADERS, timeout=40)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for link in root.iter("link"):
            href = link.get("href", "")
            if href.endswith(".tar.gz"):
                return href.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
    except Exception as exc:  # noqa: BLE001
        print(f"    [oa] resolve failed: {exc}", flush=True)
    return None


def _looks_like(data: bytes, kind: str) -> bool:
    if kind == "zip":
        return data[:2] == b"PK"
    if kind == "gz":
        return data[:2] == b"\x1f\x8b"
    if kind == "xlsx":
        return data[:2] == b"PK"
    if kind == "xls":
        return data[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
    return False


def _is_spreadsheet_name(name: str) -> bool:
    return name.lower().endswith(SPREADSHEET_EXTS)


def _extract_spreadsheets(data: bytes, dest: Path) -> list[str]:
    """Extract spreadsheet members from a zip or tar.gz byte blob; return written filenames."""
    written: list[str] = []
    if _looks_like(data, "zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            for m in zf.namelist():
                if _is_spreadsheet_name(m) and not m.endswith("/"):
                    out = dest / Path(m).name
                    out.write_bytes(zf.read(m))
                    written.append(out.name)
    elif _looks_like(data, "gz"):
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            for m in tf.getmembers():
                if m.isfile() and _is_spreadsheet_name(m.name):
                    f = tf.extractfile(m)
                    if f is not None:
                        out = dest / Path(m.name).name
                        out.write_bytes(f.read())
                        written.append(out.name)
    return written


def _save_single(data: bytes, url: str, dest: Path) -> list[str]:
    """Save a single (non-archive) spreadsheet download."""
    name = url.split("?")[0].rstrip("/").split("/")[-1]
    if not _is_spreadsheet_name(name):
        name = name + ".xlsx" if _looks_like(data, "xlsx") else name
    out = dest / name
    out.write_bytes(data)
    return [out.name]


def _has_data(dest: Path) -> bool:
    return dest.exists() and any(_is_spreadsheet_name(p.name) for p in dest.iterdir() if p.is_file())


def _restage_ecl8(dest: Path) -> list[str]:
    """Copy the already-downloaded ECL8 xlsx from data/raw/legacy/ into the organism tree."""
    src = E.legacy_essentiality_dir("literature", "eichelberger2024_ECL8")
    written = []
    if src.exists():
        for p in sorted(src.glob("*.xlsx")):
            out = dest / p.name
            shutil.copy2(p, out)
            written.append(out.name)
    return written


def _fetch_one(ds: Dataset) -> dict:
    dest = E.essentiality_raw_dir(ds.organism, ds.key)
    if _has_data(dest):
        return {"dataset": ds.key, "status": "ok", "route": "cached",
                "http": "", "note": "already present"}

    # ECL8 is re-staged from legacy, not re-fetched.
    if ds.key == "eichelberger2024_ECL8":
        w = _restage_ecl8(dest)
        if w:
            return {"dataset": ds.key, "status": "ok", "route": "legacy",
                    "http": "", "note": f"restaged {len(w)} xlsx"}

    # candidate URL ladder
    candidates: list[tuple[str, str]] = [("cdn", u) for u in ds.direct_urls]
    if ds.pmcid and ds.open_access:
        candidates.append(("europepmc", EUROPEPMC_SUPP.format(pmcid=ds.pmcid)))
        oa = _resolve_oa_tarball(ds.pmcid)
        if oa:
            candidates.append(("ncbi_oa", oa))

    last_http = ""
    for route, url in candidates:
        try:
            print(f"    [{route}] {url}", flush=True)
            r = _get(url)
            last_http = str(r.status_code)
            data = r.content
            if _looks_like(data, "zip") or _looks_like(data, "gz"):
                w = _extract_spreadsheets(data, dest)
            elif _looks_like(data, "xlsx") or _looks_like(data, "xls"):
                w = _save_single(data, url, dest)
            else:
                print(f"      rejected: not a spreadsheet/archive ({data[:16]!r})", flush=True)
                continue
            if w:
                return {"dataset": ds.key, "status": "ok", "route": route,
                        "http": last_http, "note": f"got {len(w)} file(s): {','.join(w[:4])}"}
            print("      archive held no spreadsheet members", flush=True)
        except Exception as exc:  # noqa: BLE001
            last_http = str(getattr(getattr(exc, "response", None), "status_code", "") or "")
            print(f"      failed: {exc}", flush=True)

    # placeholder fallback — never blocks the pipeline
    (dest / "PLACEHOLDER.txt").write_text(
        f"{ds.paper} ({ds.strain}) — {ds.flavor}\n"
        f"open_access={ds.open_access}; pmcid={ds.pmcid}\n"
        f"Automated fetch failed/gated. {ds.note}\n"
        "Drop the supplementary spreadsheet here manually and re-run the 07b/07c parser.\n"
    )
    return {"dataset": ds.key, "status": "placeholder", "route": "none",
            "http": last_http, "note": ds.note or "gated / not fetchable"}


def _fetch_tradis_compendium() -> dict:
    """Download the key Enterobacteriaceae-TraDIS files into data/raw/other/essentiality/."""
    dest = E.REPO_ROOT / "data" / "raw" / "other" / "essentiality" / "enterobacteriaceae_tradis"
    dest.mkdir(parents=True, exist_ok=True)
    got = []
    for rel in TRADIS_FILES:
        out = dest / Path(rel).name
        if out.exists() and out.stat().st_size > 0:
            got.append(out.name)
            continue
        try:
            r = _get(f"{TRADIS_RAW}/{rel}", timeout=(10, 180))
            out.write_bytes(r.content)
            got.append(out.name)
            print(f"    [github] {rel} -> {len(r.content)} bytes", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"    [github] {rel} FAILED: {exc}", flush=True)
    return {"dataset": "enterobacteriaceae_tradis", "organism": "other",
            "status": "ok" if got else "placeholder", "route": "github",
            "http": "", "note": f"{len(got)}/{len(TRADIS_FILES)} files"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", help="fetch only this dataset key (default: all)")
    args = ap.parse_args()

    rows = []
    if not args.only or args.only == "enterobacteriaceae_tradis":
        print("[other] enterobacteriaceae_tradis — Goodall/Gardner TraDIS compendium", flush=True)
        rows.append(_fetch_tradis_compendium())
        print(f"  -> {rows[-1]['status']} ({rows[-1]['route']}) {rows[-1]['note']}", flush=True)

    todo = [d for d in MANIFEST if not args.only or d.key == args.only]
    for ds in todo:
        print(f"[{ds.organism}] {ds.key} — {ds.paper}", flush=True)
        rows.append(_fetch_one(ds))
        print(f"  -> {rows[-1]['status']} ({rows[-1]['route']}) {rows[-1]['note']}", flush=True)

    # one status table per organism touched (org carried on each row)
    status = pd.DataFrame(rows)
    org_of = {d.key: d.organism for d in MANIFEST}
    org_of["enterobacteriaceae_tradis"] = "other"
    status["organism"] = status["dataset"].map(org_of)
    for org, sub in status.groupby("organism"):
        base = (E.REPO_ROOT / "data" / "raw" / org / "essentiality") if org == "other" \
            else E.essentiality_raw_dir(org)
        base.mkdir(parents=True, exist_ok=True)
        out = base / "fetch_status.tsv"
        sub.drop(columns=["organism"]).to_csv(out, sep="\t", index=False)
        print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)}", flush=True)

    n_ok = (status["status"] == "ok").sum()
    print(f"\nDone: {n_ok}/{len(status)} datasets fetched; "
          f"{(status['status']=='placeholder').sum()} placeholders.", flush=True)


if __name__ == "__main__":
    main()
