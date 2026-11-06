"""Offline check of webapp/DESIGN.md against the parts of the Google DESIGN.md spec that matter here.

Checks, with no network and no npm: front matter parses; `name` and `colors.primary` exist; every
`{token.path}` reference resolves; component backgroundColor/textColor pairs meet WCAG AA 4.5:1;
`##` sections are in the spec's order with no duplicates. It is NOT the official linter
(`npx @google/design.md lint`, which needs a download); run that too where npm is allowed.

Run with the `gradi` env:  python tools/check_design_md.py [path]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
ORDER = ["Overview", "Colors", "Typography", "Layout", "Elevation & Depth", "Shapes", "Components", "Do's and Don'ts"]
ALIASES = {"Brand & Style": "Overview", "Layout & Spacing": "Layout", "Elevation": "Elevation & Depth"}


def lum(hex_):
    h = hex_.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4  # noqa: E731
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(a, b):
    la, lb = sorted((lum(a), lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "webapp" / "DESIGN.md"
    text = path.read_text(encoding="utf-8")
    m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    assert m, "no front matter"
    fm = yaml.safe_load(m.group(1))
    body = text[m.end():]
    errors = []
    if not fm.get("name"):
        errors.append("missing name")
    if "primary" not in fm.get("colors", {}):
        errors.append("missing colors.primary")

    def resolve(ref):
        node = fm
        for part in ref.split("."):
            if not isinstance(node, dict) or part not in node:
                return None
            node = node[part]
        return node

    refs = re.findall(r"\{([a-z0-9_.-]+)\}", m.group(1))
    for r in refs:
        if resolve(r) is None:
            errors.append(f"broken reference {{{r}}}")
    for name, comp in fm.get("components", {}).items():
        bg, fg = comp.get("backgroundColor"), comp.get("textColor")
        if bg and fg:
            bgv = resolve(bg.strip("{}")) if bg.startswith("{") else bg
            fgv = resolve(fg.strip("{}")) if fg.startswith("{") else fg
            c = contrast(bgv, fgv)
            print(f"  contrast {name:16s} {fgv} on {bgv}: {c:.2f}" + ("" if c >= 4.5 else "  < 4.5 AA"))
            if c < 4.5:
                errors.append(f"contrast {name} {c:.2f}")
    heads = [ALIASES.get(h.strip(), h.strip()) for h in re.findall(r"^## (.+)$", body, re.M)]
    if len(set(heads)) != len(heads):
        errors.append("duplicate ## heading")
    idx = [ORDER.index(h) for h in heads if h in ORDER]
    if idx != sorted(idx):
        errors.append(f"sections out of order: {heads}")
    unknown = [h for h in heads if h not in ORDER]
    print(f"  tokens: {len(fm.get('colors', {}))} colors, {len(fm.get('typography', {}))} typography, "
          f"{len(fm.get('components', {}))} components; {len(refs)} references; sections {heads}")
    if unknown:
        print(f"  note: sections not in the spec: {unknown}")
    for e in errors:
        print("ERROR", e)
    print("OK" if not errors else f"{len(errors)} problem(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
