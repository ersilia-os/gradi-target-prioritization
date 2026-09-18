"""Check that every path CLAUDE.md claims actually exists.

CLAUDE.md is the map. A wrong path in it is worse than no map, because it is believed — and the
two ways it goes wrong are both silent: a bulk rename that catches a case it should not (smoke
files are `scratch/`, never `evidence/`), and a path rewrite that doubles a segment
(`scripts/localization/localization/workers/`). Both read perfectly.

Run it after editing CLAUDE.md, and after any restructure:

    python tools/check_claude_md.py            # non-zero exit if anything is stale

What it checks, and what it deliberately does not:
  - `scripts/<task>/<name>.py`, `src/<name>.py`, `docs/<name>.md`, `tools/<name>.py` must exist.
  - `data/...` and `output/...` paths must exist, EXCEPT when they contain a `<placeholder>`,
    `{brace}` or `*`, which are patterns -- those are globbed and need at least one match.
  - Deleted-after-use dumps (ChEMBL, eggNOG, OrthoDB) are exempt: they are documented as removed,
    so their absence is correct, not stale.
  - Prose like `legacy/...` is skipped: that tree is frozen and indexed elsewhere.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DOC = REPO / "CLAUDE.md"

# documented as deleted after the run, or fetched on demand -- absence is correct
EXEMPT = (
    "data/source/eggnog", "data/raw/other/chembl", "data/raw/other/bindingdb",
    "data/source/orthodb/odb", "data/processed/legacy",
)

# Paths CLAUDE.md quotes precisely BECAUSE they are wrong -- cautionary examples in the
# "Maintaining this file" section, and the pre-restructure layout it describes as history.
# They must stay broken; that is the point of quoting them.
NEGATIVE_EXAMPLES = {
    "scripts/localization/localization/workers/",   # a doubled segment from a bad path rewrite
    "docs/02_function.md",                          # the old numbered doc name, named as forbidden
    "scripts/NN_*.py",                              # the layout this repo moved away from
}

PATH_RE = re.compile(r"`([a-z][A-Za-z0-9_./<>{}*-]*"
                     r"(?:\.py|\.md|\.tsv|\.npz|\.sh|\.txt|\.obo|\.h5|\.faa|\.pkl|/))`")


def looks_like_path(s: str) -> bool:
    return "/" in s or s.endswith((".py", ".md", ".tsv", ".npz"))


def main() -> int:
    text = DOC.read_text()
    seen: set[str] = set()
    stale: list[str] = []
    checked = 0

    for m in PATH_RE.finditer(text):
        raw = m.group(1)
        if raw in seen or not looks_like_path(raw):
            continue
        seen.add(raw)
        if raw in NEGATIVE_EXAMPLES:
            continue
        if raw.startswith("legacy/") or any(raw.startswith(e) for e in EXEMPT):
            continue
        # only check paths rooted at a directory we own
        if not raw.split("/")[0] in {"scripts", "src", "docs", "tools", "data", "output"}:
            continue
        checked += 1
        if any(c in raw for c in "<>{}*"):
            pattern = re.sub(r"<[^>]*>|\{[^}]*\}", "*", raw).rstrip("/")
            if not list(REPO.glob(pattern)):
                stale.append(f"{raw}   (no match for glob {pattern})")
        elif not (REPO / raw).exists():
            stale.append(raw)

    line_of = {}
    for i, line in enumerate(text.splitlines(), 1):
        for s in stale:
            key = s.split("   ")[0]
            if f"`{key}`" in line and key not in line_of:
                line_of[key] = i

    print(f"checked {checked} distinct paths in CLAUDE.md")
    if stale:
        print(f"\n{len(stale)} STALE:")
        for s in stale:
            key = s.split("   ")[0]
            print(f"  CLAUDE.md:{line_of.get(key, '?')}  {s}")
        return 1
    print("all good -- every path CLAUDE.md names exists")
    return 0


if __name__ == "__main__":
    sys.exit(main())
