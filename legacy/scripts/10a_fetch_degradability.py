"""Fetch the primary degradability source tables (docs §3, tracks 3.2b/3.2c/3.3b).

The legacy prototype never obtained any of these: `data/raw/legacy/clp_substrates/SOURCE.md` admits
its Flynn table is a 45-row hand-curated stand-in ("not a verbatim copy of Flynn 2003 Tables 1-2")
and its Nagar table is 35 rows against ~1,149 in the real Data Set S2. Everything downstream was
bottlenecked there: only 21 K. pneumoniae proteins carried any Clp-trap flag and 15 a half-life class.

Route reality, probed directly (see docs/degradability_references.md):

  * **Nagar 2021** (mSystems, PMC7857536) is CC-BY and fully automated — the Europe PMC
    `supplementaryFiles` endpoint returns a zip carrying `sd001` (1,602-protein raw pulsed-SILAC time
    course), `sd002` (1,149 proteins with `label` + `half_life`) and `sd003` (1,149 `Stability`
    classes), all keyed by UniProt accession. This one fetch is the axis's main unlock.
  * **Every Clp-trap paper is gated.** Europe PMC answers "is not open access" and NCBI-OA answers
    `idIsNotOpenAccess` for Lunge 2020 (PMC7363115), Bhat 2013 (PMC3681837) and Graham 2013
    (PMC3807464); Flynn 2003 (Cell Press), Neher 2006 (Cell Press), Feng 2013 (ACS) and Ziemski 2021
    (Wiley) have no PMC record at all. They need an authenticated browser session — the same
    chrome-devtools `evaluate_script` same-origin fetch already used for the gated essentiality
    tables (Jana 2023 / Goodall / Hawkins / Rousset 2021).
  * **Flynn 2003 Tables 1-2 are MAIN-TEXT tables, not a paywalled supplement** — `SOURCE.md` is wrong
    about this. The fetch target is the article HTML, parsed into a TSV; Table 2 carries the five
    authoritative CM/NM recognition-signal consensuses that docs §3.1a/§3.1b claim to encode.

So this script does what it can headlessly and then reports precisely what is missing. `--stage`
ingests files a browser session has already dropped into the landing directory (or any directory you
point at), so the Chrome pass and the automated pass converge on the same layout and the parsers in
10c/10d never need to know which route a file arrived by.

Landing zones (eosvc-tracked, gitignored):
  data/raw/ecoli/degradability/<key>/     — E. coli-native evidence (Flynn, Neher, Nagar, RegulonDB)
  data/raw/other/degradability/<key>/     — cross-species trap papers (Bhat, Lunge, Graham, Feng, …)
  …/fetch_status.tsv                      — one row per dataset: status, route, http, note

Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import io
import shutil
import subprocess
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
from src import degradability as D  # noqa: E402

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124 Safari/537.36")
HEADERS = {"User-Agent": UA, "Connection": "close"}
DATA_EXTS = (".xlsx", ".xls", ".csv", ".tsv", ".txt", ".html", ".json", ".pdf", ".docx")
# `.txt` is a legitimate data extension here (RegulonDB ships .txt), so the placeholder marker has to
# be excluded by name — otherwise `_has_data()` reads its own "not obtained" note as a success and the
# status table reports coverage the axis does not have.
PLACEHOLDER_NAME = "PLACEHOLDER.txt"
EUROPEPMC_SUPP = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
EUROPEPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
NCBI_OA = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"


@dataclass
class Dataset:
    key: str                        # landing subdirectory name
    organism: str                   # which tree it belongs under ("ecoli" / "other")
    paper: str                      # short citation
    track: str                      # which docs §3 track it feeds
    pmcid: str = ""                 # drives the Europe PMC + NCBI-OA routes
    direct_urls: list[str] = field(default_factory=list)   # publisher-CDN candidates (route 1)
    open_access: bool = True        # False -> skip the automated PMC routes entirely, they 403/deny
    chrome_url: str = ""            # what to open in an authenticated browser session
    filename: str = ""              # name to save a single-file download under (figshare et al. serve
                                    # opaque numeric paths, so the URL basename carries no provenance)
    note: str = ""


# --- the manifest ------------------------------------------------------------------------------
MANIFEST: list[Dataset] = [
    Dataset("nagar2021_halflife", "ecoli",
            "Nagar 2021 (mSystems)", "3.2c half-life",
            pmcid="PMC7857536", open_access=True,
            note="CC-BY; Europe PMC supplementaryFiles works headlessly. sd002 = label + half_life "
                 "(1,149 UniProt-keyed); sd003 = Stability class; sd001 = raw SILAC time course."),
    Dataset("flynn2003_clpxp_trap", "ecoli",
            "Flynn 2003 (Mol Cell)", "3.2b trap + 3.1a/3.1b motif consensuses",
            open_access=False,
            chrome_url="https://www.cell.com/molecular-cell/fulltext/S1097-2765(03)00060-1",
            note="No PMC record. Tables 1-2 are MAIN-TEXT tables (SOURCE.md wrongly calls them a "
                 "paywalled supplement). Table 1 = ClpXP trap census; Table 2 = the CM1/CM2/NM1-3 "
                 "consensuses. Save the article HTML here and re-run with --stage."),
    Dataset("neher2006_clpxp_sos", "ecoli",
            "Neher 2006 (Mol Cell)", "3.2b trap (SOS/DNA-damage condition)",
            open_access=False,
            chrome_url="https://www.cell.com/molecular-cell/fulltext/S1097-2765(06)00168-7",
            note="Not in the docs §3 paper list; a second E. coli ClpXP substrate-trap set from the "
                 "same lab under DNA damage, which roughly doubles the E. coli trap positives."),
    Dataset("calloni2012_dnak", "ecoli",
            "Calloni 2012 (Cell Rep)", "3.4 chaperone dependence + assembly state",
            open_access=False,
            chrome_url="https://www.cell.com/cell-reports/fulltext/S2211-1247(11)00017-9",
            note="SUPPLIED MANUALLY — ScienceDirect `attachment/mmc{1,2,3}.xls` all answer 403, so there is "
                 "no automated route. 8-sheet legacy BIFF workbook (needs `xlrd>=2.0.1`, NOT openpyxl); "
                 "**the header row differs per sheet** (Table S1 at index 7, Table S8 at 8) so detect it per "
                 "sheet. Table S1 = the 674-interactor DnaK set (688 rows) and also carries `GroEL class`, "
                 "`Oligomeric state`, `Sol (%)` and `Localization`. Keys are `EG` (EcoGene) + `Gene Name`, no "
                 "UniProt — join via src/essentiality.gene_aliases_to_uniprot. See its SOURCE.md. Value caveat: "
                 "eSOL already gives quantitative chaperone dependence for 788 proteins, so treat this as "
                 "corroborating/annotation, and as an independent route to `assembly_state` for E. coli."),
    # --- E. coli feature layer, added 2026-08-06 -----------------------------------------------
    # All routes probed on that date. Two lessons, in tension, and both cost time:
    #   1. "The journal is paywalled" does NOT imply "the data is unreachable". Cappelletti 2021 is
    #      Cell Press yet sits in the PMC OA subset, and Springer serves supplementary workbooks free
    #      from `media.springernature.com` even for paywalled articles (Conlon 2013, Gupta, Mateus).
    #      So always probe before declaring something gated.
    #   2. But VERIFY THE PMCID BY ITS TITLE before trusting a 200. `to2021_refoldability` was first
    #      entered with PMC8382223, which is an unrelated Anal Chem paper; the fetch "succeeded" and
    #      landed the wrong supplement. A 200 from the wrong record looks exactly like a 200 from the
    #      right one. `_fetch_fulltext()` caches the article XML for OA datasets, which makes this
    #      checkable after the fact — use it.
    Dataset("gupta2024_turnover", "ecoli",
            "Gupta 2024 (Nat Commun)", "3.2c turnover + protease attribution",
            pmcid="PMC11246515", open_access=True,
            direct_urls=["https://media.springernature.com/original/springer-static/esm/"
                         "art%3A10.1038%2Fs41467-024-49920-8/MediaObjects/"
                         "41467_2024_49920_MOESM4_ESM.xlsx"],
            filename="gupta2024_TableS1_turnover.xlsx",
            note="THE PRIMARY LABEL SET. Table S1 = half-lives for ~3,200 proteins x 13-14 conditions "
                 "INCLUDING the protease-KO panel (clpP / lon / hslV / triple / smpB / ftsH), keyed by "
                 "UniProt accession (`sp|P00350|6PGD_ECOLI`). A per-protein change in half-life on "
                 "protease deletion is a LABEL, not a feature. The Europe PMC zip carries the other "
                 "Supp Data files (protease assignment, ~600 in-vivo N-termini) if the CDN route is used."),
    Dataset("cappelletti2021_lipms", "ecoli",
            "Cappelletti 2021 (Cell)", "3.4 proteolytic accessibility (LiP-MS)",
            pmcid="PMC7836100", open_access=True,
            note="The most mechanistically apt single column for this axis — limited proteolysis measures "
                 "exactly the accessible/flexible-region property AAA+ engagement requires. ~1,900 proteins, "
                 "8 carbon sources, peptide resolution. Strain BW25113, the SAME as Mateus 2018, so the two "
                 "join cleanly. Cell Press but IN the PMC OA subset — no browser needed."),
    Dataset("mateus2018_tpp", "ecoli",
            "Mateus 2018 (Mol Syst Biol)", "3.4 thermal stability",
            pmcid="PMC6056769", open_access=True,
            note="CC-BY. Dataset EV1 = 1,738 fitted apparent Tm (of 1,831 identified), strain BW25113. "
                 "MUST-HANDLE CAVEAT: Tm follows a cell-surface -> cytoplasm high-to-low gradient, so it "
                 "partly encodes compartment and has to be RESIDUALISED ON LOCALIZATION before use, or it "
                 "will smuggle the 09g localization call back into the score."),
    Dataset("mateus2020_tpp121", "ecoli",
            "Mateus 2020 (Nature)", "3.4 conformational plasticity",
            pmcid="PMC7612278", open_access=True,
            direct_urls=["https://media.springernature.com/original/springer-static/esm/"
                         "art%3A10.1038%2Fs41586-020-3002-5/MediaObjects/"
                         "41586_2020_3002_MOESM5_ESM.xlsx"],
            filename="mateus2020_SuppData_tpp121strains.xlsx",
            note="TPP across 121 strains (mostly Keio deletions); Supp Data 5 = abundance + thermal-stability "
                 "score for 1,764 proteins. The feature worth deriving is STABILITY-SCORE VARIANCE across the "
                 "121 perturbations — a conformational-plasticity measure arguably more informative than any "
                 "single Tm. MOESM6 is 60 MB and not needed for that."),
    Dataset("to2021_refoldability", "ecoli",
            "To 2021 (JACS)", "3.4 refoldability / kinetic stability",
            pmcid="PMC8650709", open_access=False,
            direct_urls=["https://www.biorxiv.org/content/biorxiv/early/2020/10/22/"
                         "2020.08.28.273110/DC1/embed/media-1.pdf?download=true"],
            filename="to2021_biorxiv_SI.pdf",
            chrome_url="https://pubs.acs.org/doi/10.1021/jacs.1c03270",
            note="396 of 1,198 proteins (33%) NON-REFOLDABLE after 2 h — a kinetic-stability/metastability "
                 "proxy complementary to Tm; protein- AND domain-level tables. "
                 "!! CORRECTION 2026-08-06: an earlier version of this entry used pmcid PMC8382223 and claimed "
                 "the paper was in the PMC OA subset. **PMC8382223 is a DIFFERENT PAPER** ('Evaluation of the "
                 "Higher Order Structure of Biotherapeutics Embedded in Hydrogels', Anal Chem) — the 200 OK "
                 "came from an unrelated article and it fetched `ac1c01850_si_001.pdf`, which was deleted. The "
                 "real record is PMC8650709 and it is **NOT open access**. Lesson: verify a PMCID by its title "
                 "before trusting a 200. "
                 "The MACHINE-READABLE TABLES ARE ACS-GATED: bioRxiv 2020.08.28.273110 (v1 and v2) ships only "
                 "`media-1.pdf`, and NSF-PAR `purl/10311595` is the accepted-manuscript PDF. So the direct_url "
                 "here gets the open SI **PDF only** — the per-protein table still needs a browser session or "
                 "text extraction from the PDF."),
    Dataset("to2022_chaperone_refold", "ecoli",
            "To 2022 (PNAS)", "3.4 chaperone-assisted refolding",
            pmcid="PMC9704704", open_access=True,
            chrome_url="https://www.pnas.org/doi/10.1073/pnas.2210536119",
            note="The follow-up separating INTRINSIC from CHAPERONE-RESCUED refoldability in a cytosol-like "
                 "milieu. Pairs with to2021_refoldability. "
                 "ROUTE NOTE 2026-08-06: the Europe PMC supplementaryFiles zip for PMC9704704 contains "
                 "**figure images only** (12 jpg/gif, no data) — verified by listing the archive. So the "
                 "endpoint returns 200 with a real zip and still yields nothing usable. Needs the PNAS "
                 "DCSupplemental route via a browser session. LOW PRIORITY: eSOL already supplies quantitative "
                 "chaperone dependence for 788 proteins."),
    Dataset("mackrell2026_boncat", "ecoli",
            "MacKrell 2026 (PNAS)", "3.2c turnover (independent validation)",
            pmcid="PMC12974527", open_access=True,
            note="BONCAT + TMT + FAIMS, time-resolved. THE ONLY SOURCE COVERING STATIONARY PHASE (1,339 "
                 "proteins) as well as exponential (1,810); 56 and 88 pronouncedly unstable. Best independent "
                 "validation set for whatever we fit on Gupta/Nagar. Supp S4 has ML predictions for "
                 "unmeasured proteins — do NOT mix those into a label set."),
    Dataset("schmidt2016_abundance", "ecoli",
            "Schmidt 2016 (Nat Biotechnol)", "3.5 abundance / resynthesis burden",
            pmcid="PMC4888949", open_access=True,
            note="2,359 proteins in MOLECULES PER CELL across 22 conditions (~55% of ORFs, >95% of proteome "
                 "mass). Supplies the `resynthesis_burden` modifier and the abundance term Won 2024 measured "
                 "at r = -0.69. Free author manuscript carries the tables."),
    Dataset("gyorkei2022_solubility", "ecoli",
            "Gyorkei 2022 (Sci Rep)", "3.4 in-vivo solubility limits",
            pmcid="PMC9023497", open_access=True,
            note="2,577 cytosolic proteins (ASKA GFP library, K-12 AG1): per-protein solubility threshold "
                 "plus a 3-class label {soluble, rapidly aggregating, slowly aggregating}. In-VIVO, which "
                 "complements eSOL's cell-free PURE-system measurement."),
    Dataset("niwa2022_lon_clpxp", "ecoli",
            "Niwa 2022 (Molecules)", "3.3 protease attribution (Lon vs ClpXP vs HslUV)",
            pmcid="PMC9228906", open_access=True,
            note="The best public side-by-side of which protease degrades which E. coli protein; Lon is "
                 "dominant for the ~80 obligate GroE substrates. Regime A evidence, so provenance/cross-check "
                 "rather than a scored column under the ClpP decision."),
    Dataset("esol_solubility", "ecoli",
            "Niwa 2009 + Niwa 2012 (PNAS)", "3.4 solubility + chaperone dependence",
            open_access=True,
            direct_urls=["https://dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip"],
            note="TWO PAPERS IN ONE FILE, verified 2026-08-06. `esol.csv`, 4,132 rows keyed by JW_ID / "
                 "B number / gene name: `Solubility (%)` for 3,173 proteins (Niwa 2009, chaperone-free PURE "
                 "system, strongly bimodal) AND `Minus/TF/GroE/KJE Sol (%)` + yields for EXACTLY 788 proteins "
                 "(Niwa 2012 chaperone effects — matching that paper's own '788 proteins x 4 conditions'). "
                 "So Niwa 2012 needs NO separate supplementary fetch; its PNAS DCSupplemental route 403s "
                 "anyway. NBDC LSDB mirror — the original tanpaku.org / tp-esol.genes.nig.ac.jp site is dead."),
    # MobiDB is fetched by `_fetch_mobidb()` rather than from this MANIFEST: it is a database bulk
    # export, not a paper supplement, and it is per-organism rather than per-publication. That
    # route also validates the payload harder than the generic ladder can (header must start
    # `acc` AND >=1000 rows, because a wrong proteome id still returns a header-only HTTP 200).
    Dataset("bhat2013_caulobacter", "other",
            "Bhat 2013 (Mol Microbiol)", "3.3b trap (C. crescentus ClpP)",
            pmcid="PMC3681837", open_access=False,
            chrome_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3681837/",
            note="In PMC but NOT in the OA subset; the PMC article page still serves /bin/ "
                 "supplementary to a browser session."),
    Dataset("lunge2020_mtb_clpc1", "other",
            "Lunge 2020 (JBC)", "3.3b ClpC1-regulated (M. tuberculosis)",
            pmcid="PMC7363115", open_access=False,
            chrome_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC7363115/",
            note="Data Set S2 = the 219 ClpC1-regulated proteins. Knockdown-derived, not a trap, so "
                 "10d weights it below the true traps — and ClpC1 is absent from both focal organisms."),
    Dataset("graham2013_saureus_clpc", "other",
            "Graham 2013 (J Bacteriol)", "3.3b trap (S. aureus ClpC)",
            pmcid="PMC3807464", open_access=False,
            chrome_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3807464/",
            note="In PMC but not the OA subset. Spec suggestion: extend the §3.3b pool with this."),
    Dataset("feng2013_saureus_clpp", "other",
            "Feng 2013 (J Proteome Res)", "3.3b trap (S. aureus ClpXP + ClpCP)",
            open_access=False,
            chrome_url="https://pubs.acs.org/doi/10.1021/pr300394r",
            note="ACS, no PMC record at all — the hardest of the set."),
    Dataset("ziemski2021_mtb_clpcp", "other",
            "Ziemski 2021 (FEBS J)", "3.3b interaction screen (Mtb ClpCP)",
            open_access=False,
            chrome_url="https://febs.onlinelibrary.wiley.com/doi/10.1111/febs.15335",
            note="Wiley, no PMC record. Interaction screen (~67 candidates), weakest inference in "
                 "the pool; 10d weights it 0.40."),

    # --- ACTIVATED-ClpP (partnerless) evidence -------------------------------------------------
    # These two are the ONLY datasets that measure the criterion the Gr-ADI proposal actually gates
    # on: "susceptibility to degradation by the activated ClpP of the target organism only (i.e. in
    # the absence of an unfoldase partner)". Every other entry above is unfoldase-DEPENDENT evidence
    # (a ClpXP/ClpAP/ClpCP trap), which is a different mechanism — see src/degradability.PARTNERLESS.
    # The v5 Final proposal picked its first two targets (DnaK, AcpP) off precisely these tables.
    Dataset("jacques2020_onc212_saureus", "other",
            "Jacques 2020 (Genetics)", "3.3c activated-ClpP proteomics (ONC212, S. aureus)",
            pmcid="PMC7153937", open_access=True,
            direct_urls=["https://ndownloader.figshare.com/files/21767271"],
            filename="jacques2020_TableS3_saureus_onc212.xlsx",
            chrome_url="https://doi.org/10.25386/genetics.11873841",
            note="Table S3 via the GSA figshare deposit (10.25386/genetics.11873841) — fully open, so "
                 "this one needs no browser. 1,620 S. aureus proteins x 30 uM ONC212. Two readouts: "
                 "`24H_log2_fold-change` (abundance loss) and "
                 "`10-40_minutes_non-tryptic_peptides_log2_fold-change` — the early non-tryptic "
                 "peptides are direct CLEAVAGE products and are the better evidence of degradation, "
                 "since 24 h abundance also moves with growth arrest and regulon change. Keyed by "
                 "RefSeq (`Gene`, ODV*.1) with UniProt in `Mu50_homolog` and SACOL locus tags in "
                 "`COL_homolog`. NOTE the Europe PMC supplementaryFiles endpoint for this PMCID "
                 "returns only figure images — the data lives on figshare only. Table S2 (human "
                 "NALM-6) is deliberately not fetched; it is the human-mitochondrial-ClpP off-target "
                 "picture, not target evidence."),
    Dataset("conlon2013_adep4_saureus", "other",
            "Conlon 2013 (Nature)", "3.3c activated-ClpP proteomics (ADEP4, S. aureus)",
            pmcid="PMC4031760", open_access=False,
            direct_urls=["https://media.springernature.com/original/springer-static/esm/"
                         "art%3A10.1038%2Fnature12790/MediaObjects/"
                         "41586_2013_BFnature12790_MOESM97_ESM.xlsx"],
            filename="conlon2013_ADEP4_MRSA_TableS1-S2.xlsx",
            chrome_url="https://www.nature.com/articles/nature12790",
            note="The ADEP4 reference proteome; Jacques 2020 benchmarks against it and reports the "
                 "same signature (200 most-depleted enriched for ribosome-related functions, "
                 "i.e. nascent chains). ROUTE NOTE: the article is paywalled and Europe PMC answers "
                 "'is not open access', but the supplementary workbook itself is FREE on "
                 "`media.springernature.com/original/springer-static/esm/...` — note the host. The "
                 "`static-content.springer.com` mirror answers 403 to bare curl and PMC now fronts "
                 "PMC4031760 with reCAPTCHA, so both of those routes look like a paywall and are not. "
                 "Two sheets, both keyed by SACOL locus tag (strain COL), `Average` = log2 FC "
                 "(ADEP4/control) with an adjusted p-value: 'Table S1' (1,715 proteins, fully tryptic "
                 "peptides) = ABUNDANCE loss; 'Table S2' (2,382 rows, partially tryptic peptides, one "
                 "terminus made by an endogenous protease) = direct CLEAVAGE evidence. SACOL is also "
                 "Jacques 2020's `COL_homolog`, so the two activator papers share a join key."),
]

# RegulonDB sigma-factor regulons — the spec's sigma32/sigmaS suggestion. Carried as ANNOTATION
# ONLY, never scored: the correlation with Clp substrate status is confounded with function (the
# proteostasis network's clients are heat-shock/stationary-phase regulators), and scoring it would
# re-import the exact bias that produced the legacy top-10 (rpoS, rpoH, dps, lexA, …).
REGULONDB_FILES = {
    "network_tf_gene.txt": "https://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt",
    "network_sigma_gene.txt": "https://regulondb.ccg.unam.mx/menu/download/datasets/files/network_sigma_gene.txt",
}


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=2, max=20))
def _get(url: str, timeout=(10, 120)) -> requests.Response:
    r = requests.get(url, headers=HEADERS, stream=True, timeout=timeout)
    r.raise_for_status()
    # `raise_for_status()` passes 2xx, but 202 Accepted means "not ready yet, ask again" — figshare's
    # ndownloader prepares large deposits asynchronously and answers 202 with a ZERO-LENGTH body until
    # the file is warm. Banking that as success writes an empty file and reports coverage we do not
    # have, which is the same class of silent-gap bug `_is_html_shell()` guards against. Raise so
    # tenacity backs off and retries; an empty 200 is equally unusable.
    if r.status_code == 202 or not r.content:
        raise requests.HTTPError(f"{r.status_code} with {len(r.content)}-byte body (not ready)",
                                 response=r)
    return r


def _dest_dir(ds: Dataset) -> Path:
    return D.xbac_raw_dir(ds.key) if ds.organism == "other" else D.degradability_raw_dir(ds.organism, ds.key)


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
    if kind in ("zip", "xlsx"):
        return data[:2] == b"PK"
    if kind == "gz":
        return data[:2] == b"\x1f\x8b"
    if kind == "xls":
        return data[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
    return False


def _is_ooxml(data: bytes) -> bool:
    """True when a PK blob is itself an Office document (xlsx/docx), not an archive OF documents.

    Both share the `PK` magic, so a magic-byte test alone cannot tell "a zip containing Table_S2.xlsx"
    from "Table_S2.xlsx". Getting this wrong is silent: the archive branch opens the workbook, finds no
    member whose name ends in a data extension (the real members are `xl/worksheets/sheet1.xml`), and
    reports the dataset as gated even though the bytes were served with HTTP 200. Office Open XML
    always carries a top-level `[Content_Types].xml`, so that is the discriminator.
    """
    if not _looks_like(data, "zip"):
        return False
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            names = zf.namelist()
        return "[Content_Types].xml" in names or any(n.startswith("xl/") for n in names)
    except zipfile.BadZipFile:
        return False


def _is_data_name(name: str) -> bool:
    return name != PLACEHOLDER_NAME and name.lower().endswith(DATA_EXTS)


def _extract_data_members(data: bytes, dest: Path, _depth: int = 0) -> list[str]:
    """Extract data members from a zip or tar.gz byte blob; return written filenames.

    Recurses ONE level into nested archives. MDPI (and some other publishers) ship the whole
    supplement as a single inner `…-s001.zip` inside the Europe PMC bundle, and since `.zip` is not a
    data extension the outer pass would report "archive held no data members" for a bundle that
    actually contains everything — a silent false negative. Depth is capped at 1 so a zip-bomb or a
    self-referential archive cannot spin.
    """
    written: list[str] = []
    if _looks_like(data, "zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            for m in zf.namelist():
                if m.endswith("/"):
                    continue
                if _is_data_name(m):
                    out = dest / Path(m).name
                    out.write_bytes(zf.read(m))
                    written.append(out.name)
                elif _depth == 0 and m.lower().endswith((".zip", ".tar.gz", ".tgz")):
                    inner = zf.read(m)
                    got = _extract_data_members(inner, dest, _depth + 1)
                    if got:
                        print(f"      [nested] {Path(m).name} -> {len(got)} file(s)", flush=True)
                    written.extend(got)
    elif _looks_like(data, "gz"):
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            for m in tf.getmembers():
                if m.isfile() and _is_data_name(m.name):
                    f = tf.extractfile(m)
                    if f is not None:
                        out = dest / Path(m.name).name
                        out.write_bytes(f.read())
                        written.append(out.name)
    return written


def _save_single(data: bytes, url: str, dest: Path, filename: str = "") -> list[str]:
    name = filename or url.split("?")[0].rstrip("/").split("/")[-1] or "download"
    if not _is_data_name(name):
        name = f"{name}.xlsx" if _looks_like(data, "xlsx") else f"{name}.txt"
    out = dest / name
    out.write_bytes(data)
    return [out.name]


# Payloads that are neither a spreadsheet nor an archive, but are still legitimate data here: MobiDB
# serves TSV, and some supplements are only ever published as a PDF (To 2021's bioRxiv SI, Graham 2013).
RAW_EXTS = (".tsv", ".csv", ".txt", ".pdf", ".json")


def _wants_raw(ds: Dataset) -> bool:
    """True when this dataset declares a non-spreadsheet payload via its `filename`.

    Opt-in by design: without it, any HTML error page that slipped past `_is_html_shell` could be
    banked as data for every dataset. Requiring the manifest to name the expected extension keeps the
    permissive path narrow.
    """
    return ds.filename.lower().endswith(RAW_EXTS)


def _has_data(dest: Path) -> bool:
    return dest.exists() and any(
        _is_data_name(p.name) and p.stat().st_size > 0 for p in dest.iterdir() if p.is_file()
    )


def _write_placeholder(ds: Dataset, dest: Path, last_http: str) -> None:
    """A placeholder that says exactly how to obtain the file — never a silent gap."""
    (dest / "PLACEHOLDER.txt").write_text(
        f"{ds.paper} — docs §{ds.track}\n"
        f"open_access={ds.open_access}; pmcid={ds.pmcid or 'none'}; last_http={last_http or 'n/a'}\n"
        f"{ds.note}\n\n"
        "NOT OBTAINED by the automated ladder. To supply it:\n"
        f"  1. open {ds.chrome_url or '<publisher page>'} in an authenticated browser session\n"
        "  2. download the supplementary/article table(s) into this directory, or anywhere, then\n"
        f"     python scripts/10a_fetch_degradability.py --stage <dir> --only {ds.key}\n"
        "This dataset is listed under 'not obtained' in docs/degradability_references.md until then.\n"
    )


def _fetch_one(ds: Dataset) -> dict:
    dest = _dest_dir(ds)
    if _has_data(dest):
        names = sorted(p.name for p in dest.iterdir() if _is_data_name(p.name))
        return {"dataset": ds.key, "status": "ok", "route": "cached", "http": "",
                "note": f"already present: {','.join(names[:4])}"}

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
            # Office documents must be tested BEFORE the archive branch — see `_is_ooxml()`.
            if _is_ooxml(data) or _looks_like(data, "xls"):
                w = _save_single(data, url, dest, ds.filename)
            elif _looks_like(data, "zip") or _looks_like(data, "gz"):
                w = _extract_data_members(data, dest)
            elif _wants_raw(ds) and not _is_html_shell(data):
                # Neither a spreadsheet nor an archive, but still real data: plain TSV/CSV, or a
                # PDF-only supplement. Accepted only when the dataset OPTS IN by naming the extension
                # in `filename`, so a stray HTML error page can never be banked as data for the
                # others — and `_is_html_shell` still guards the opt-in case.
                w = _save_single(data, url, dest, ds.filename)
            else:
                print(f"      rejected: not a spreadsheet/archive ({data[:24]!r})", flush=True)
                continue
            if w:
                return {"dataset": ds.key, "status": "ok", "route": route, "http": last_http,
                        "note": f"got {len(w)} file(s): {','.join(sorted(w)[:4])}"}
            print("      archive held no data members", flush=True)
        except Exception as exc:  # noqa: BLE001
            last_http = str(getattr(getattr(exc, "response", None), "status_code", "") or "")
            print(f"      failed: {exc}", flush=True)

    if not candidates:
        print("    [skip] no automated route exists (gated publisher, no PMC OA record)", flush=True)
    _write_placeholder(ds, dest, last_http)
    return {"dataset": ds.key, "status": "placeholder", "route": "none", "http": last_http,
            "note": "gated — needs the authenticated-browser route"}


def _fetch_fulltext(ds: Dataset) -> None:
    """Cache the Europe PMC full-text XML for an OA paper (methods/growth conditions provenance).

    Needed for the Nagar growth correction: the pulsed-SILAC half-lives conflate proteolysis with
    growth dilution, so `src/degradability.NAGAR_DOUBLING_MIN` has to be traceable to the paper.
    """
    if not (ds.pmcid and ds.open_access):
        return
    out = _dest_dir(ds) / f"{ds.pmcid}_fulltext.xml"
    if out.exists() and out.stat().st_size > 0:
        return
    try:
        r = _get(EUROPEPMC_FULLTEXT.format(pmcid=ds.pmcid))
        if r.content[:5] == b"<?xml" and b"errorBean" not in r.content[:400]:
            out.write_bytes(r.content)
            print(f"    [europepmc] full text -> {out.name} ({len(r.content)} bytes)", flush=True)
    except Exception as exc:  # noqa: BLE001
        print(f"    [europepmc] full text failed: {exc}", flush=True)


def _is_html_shell(data: bytes) -> bool:
    """True for an HTML error/SPA-catch-all page masquerading as a data download.

    RegulonDB's rebuilt site serves its single-page-app `index.html` with HTTP 200 for *any* unknown
    path, so a status-code check alone reports success for a 1.6 kB `<!doctype html>` shell. PMC's
    download interstitial behaves the same way. Anything HTML and small is not a dataset.
    """
    head = data[:512].lstrip().lower()
    return head.startswith((b"<!doctype html", b"<html")) and len(data) < 200_000


def _curl(url: str, out: Path, timeout: int = 180) -> bool:
    """Fetch via the system curl. Needed only for RegulonDB, whose TLS stack negotiates in a way
    OpenSSL 3 (via requests) rejects with an SSLError while macOS curl accepts it."""
    try:
        r = subprocess.run(["curl", "-sSL", "-m", str(timeout), "-A", UA, "-o", str(out), url],
                           capture_output=True, text=True, timeout=timeout + 20)
        if r.returncode != 0 or not out.exists() or out.stat().st_size == 0:
            return False
        if _is_html_shell(out.read_bytes()):
            out.unlink()
            print(f"      rejected: HTML shell, not data ({url})", flush=True)
            return False
        return True
    except Exception:  # noqa: BLE001
        return False


def _fetch_regulondb() -> dict:
    """RegulonDB sigma/TF-gene networks (annotation only — never scored; see the module docstring)."""
    dest = D.degradability_raw_dir("ecoli", "regulondb")
    got, failed = [], []
    for name, url in REGULONDB_FILES.items():
        out = dest / name
        if out.exists() and out.stat().st_size > 0:
            got.append(name)
            continue
        try:
            r = _get(url, timeout=(10, 180))
            if _is_html_shell(r.content):
                raise RuntimeError("HTML shell, not data")
            out.write_bytes(r.content)
            got.append(name)
            print(f"    [regulondb] {name} -> {len(r.content)} bytes", flush=True)
        except Exception as exc:  # noqa: BLE001
            if _curl(url, out):
                got.append(name)
                print(f"    [regulondb:curl] {name} -> {out.stat().st_size} bytes", flush=True)
                continue
            failed.append(name)
            print(f"    [regulondb] {name} FAILED: {type(exc).__name__}", flush=True)
    if not got:
        (dest / "PLACEHOLDER.txt").write_text(
            "RegulonDB sigma32/sigmaS regulon membership (docs §3 suggestion).\n"
            "Automated download failed; carried as annotation only, so this does not gate the axis.\n"
        )
    return {"dataset": "regulondb", "status": "ok" if got else "placeholder", "route": "regulondb",
            "http": "", "note": f"{len(got)}/{len(REGULONDB_FILES)} files"
                                + (f"; failed: {','.join(failed)}" if failed else "")}


def _fetch_mobidb() -> dict:
    """MobiDB bulk, per proteome — the disorder + domain-architecture layer (docs §3.1, track 10d).

    One GET per organism supplies, at 100% of both proteomes: 8 independent disorder predictors
    (dis465, disHL, glo, th_50, espX, iups, espN, iupl) plus a `priority` consensus, an
    AlphaFold-derived disorder track, Pfam/Gene3D/merged domain boundaries (the input to
    `two_domain_architecture`), LIP (linear interacting peptide) regions, low-complexity, and TM /
    signal-peptide tracks. Long format: acc / feature / start..end / content_fraction /
    content_count / length.

    Two reasons this is worth its own route rather than a MANIFEST entry: it is a database bulk
    export rather than a paper supplement, and it is per-organism rather than per-publication.

    **The endpoint answers 405 to HEAD — GET only.** It is also slow to first byte for a whole
    proteome, hence the long read timeout.
    """
    got, failed = [], []
    for org, (proteome_id, prefix) in D.ORGANISMS.items():
        if org not in ("kpneumoniae", "ecoli"):
            continue
        # ORGANISMS stores the FASTA stem ("UP000007841_HS11286"); MobiDB wants the bare UPID.
        upid = proteome_id.split("_")[0]
        dest = D.degradability_raw_dir(org, "mobidb")
        out = dest / f"{prefix}_mobidb_{upid}.tsv"
        if out.exists() and out.stat().st_size > 0:
            got.append(out.name)
            continue
        url = f"https://mobidb.org/api/download?proteome={upid}&format=tsv"
        try:
            r = _get(url, timeout=(10, 300))
            if _is_html_shell(r.content):
                raise RuntimeError("HTML shell, not data")
            # A valid export is a TSV whose header starts with the accession column AND carries rows.
            # A wrong proteome id still returns a well-formed 61-byte header-only file with HTTP 200,
            # so the header check alone would bank an empty dataset as a success.
            if not r.content[:3].lower().startswith(b"acc"):
                raise RuntimeError(f"unexpected header {r.content[:40]!r}")
            n_rows = r.content.count(b"\n") - 1
            if n_rows < 1000:
                raise RuntimeError(f"only {n_rows} data rows — wrong proteome id or empty export")
            out.write_bytes(r.content)
            got.append(out.name)
            print(f"    [mobidb] {org} {upid} -> {len(r.content):,} bytes, {n_rows:,} rows",
                  flush=True)
        except Exception as exc:  # noqa: BLE001
            failed.append(org)
            print(f"    [mobidb] {org} FAILED: {type(exc).__name__}: {exc}", flush=True)
            (dest / "PLACEHOLDER.txt").write_text(
                f"MobiDB bulk for {org} ({upid}) — disorder predictors + domain boundaries.\n"
                f"GET {url}  (HEAD returns 405; use GET)\n"
                "Automated download failed. 10d falls back to pLDDT-derived disorder only, and\n"
                "`two_domain_architecture` cannot be computed without it.\n"
            )
    return {"dataset": "mobidb", "status": "ok" if not failed else "placeholder", "route": "mobidb",
            "http": "", "note": f"{len(got)}/2 proteomes"
                                + (f"; failed: {','.join(failed)}" if failed else "")}


def _stage(stage_dir: Path, only: str | None) -> list[dict]:
    """Ingest browser-downloaded files from `stage_dir` into the right landing directories.

    Matching is by dataset key appearing in the filename, so a browser download named e.g.
    `lunge2020_mtb_clpc1_DataSetS2.xlsx` lands automatically. Files that match nothing are reported
    rather than guessed at.
    """
    rows: list[dict] = []
    if not stage_dir.exists():
        print(f"[stage] {stage_dir} does not exist", flush=True)
        return rows
    files = [p for p in sorted(stage_dir.rglob("*")) if p.is_file() and _is_data_name(p.name)]
    targets = [d for d in MANIFEST if not only or d.key == only]
    claimed: set[Path] = set()
    for ds in targets:
        dest = _dest_dir(ds)
        moved = []
        for p in files:
            if p in claimed:
                continue
            hay = f"{p.parent.name}/{p.name}".lower()
            if ds.key.lower() in hay or (only and only == ds.key and stage_dir in p.parents):
                out = dest / p.name
                shutil.copy2(p, out)
                claimed.add(p)
                moved.append(out.name)
        if moved:
            ph = dest / "PLACEHOLDER.txt"
            if ph.exists():
                ph.unlink()
            print(f"[stage] {ds.key} <- {len(moved)} file(s): {','.join(moved[:4])}", flush=True)
            rows.append({"dataset": ds.key, "status": "ok", "route": "chrome_staged", "http": "",
                         "note": f"staged {len(moved)} file(s): {','.join(sorted(moved)[:4])}"})
    unclaimed = [p.name for p in files if p not in claimed]
    if unclaimed:
        print(f"[stage] unmatched (rename to contain a dataset key): {', '.join(unclaimed[:8])}",
              flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="operate on this dataset key only (default: all)")
    ap.add_argument("--stage", type=Path,
                    help="ingest browser-downloaded files from this directory instead of fetching")
    args = ap.parse_args()

    rows: list[dict] = []
    if args.stage:
        rows = _stage(args.stage, args.only)
    else:
        if not args.only or args.only == "mobidb":
            print("[both] mobidb — disorder predictors + domain boundaries (100% of both proteomes)",
                  flush=True)
            rows.append(_fetch_mobidb())
            print(f"  -> {rows[-1]['status']} ({rows[-1]['route']}) {rows[-1]['note']}", flush=True)
        if not args.only or args.only == "regulondb":
            print("[ecoli] regulondb — sigma32/sigmaS regulon membership (annotation only)", flush=True)
            rows.append(_fetch_regulondb())
            print(f"  -> {rows[-1]['status']} ({rows[-1]['route']}) {rows[-1]['note']}", flush=True)
        for ds in [d for d in MANIFEST if not args.only or d.key == args.only]:
            print(f"[{ds.organism}] {ds.key} — {ds.paper} (docs §{ds.track})", flush=True)
            rows.append(_fetch_one(ds))
            print(f"  -> {rows[-1]['status']} ({rows[-1]['route']}) {rows[-1]['note']}", flush=True)
            if rows[-1]["status"] == "ok":
                _fetch_fulltext(ds)

    if not rows:
        print("nothing to do", flush=True)
        return

    org_of = {d.key: d.organism for d in MANIFEST}
    org_of["regulondb"] = "ecoli"
    org_of["mobidb"] = "ecoli"   # written per-organism; the status row is filed once, under ecoli
    status = pd.DataFrame(rows)
    status["organism"] = status["dataset"].map(org_of).fillna("other")
    for org, sub in status.groupby("organism"):
        base = D.xbac_raw_dir() if org == "other" else D.degradability_raw_dir(org)
        out = base / "fetch_status.tsv"
        # merge with any previous run so a --stage pass does not erase the automated results
        if out.exists():
            prev = pd.read_csv(out, sep="\t")
            sub = (pd.concat([prev, sub.drop(columns=["organism"])], ignore_index=True)
                     .drop_duplicates("dataset", keep="last"))
        else:
            sub = sub.drop(columns=["organism"])
        sub.to_csv(out, sep="\t", index=False)
        print(f"[{org}] wrote {out.relative_to(D.REPO_ROOT)}", flush=True)

    n_ok = int((status["status"] == "ok").sum())
    gated = sorted(status.loc[status["status"] == "placeholder", "dataset"])
    print(f"\nDone: {n_ok}/{len(status)} obtained.", flush=True)
    if gated:
        print(f"NOT OBTAINED ({len(gated)}) — need the authenticated-browser route: "
              f"{', '.join(gated)}", flush=True)


if __name__ == "__main__":
    main()
