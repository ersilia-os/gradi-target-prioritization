#!/usr/bin/env python3
"""Render every ```mermaid``` code block in docs/*.md to a PNG under assets/.

A file with one block produces ``assets/<stem>.png``. A file with N blocks
produces ``assets/<stem>_1.png`` … ``<stem>_N.png``.

Uses dockerised mermaid-cli pinned to 10.9.1 — the :latest tag currently ships
without Chrome and the renderer fails. Run from anywhere; output paths are
resolved relative to the repo root.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
TMP = ROOT / "tmp" / "mmd"
OUT = ROOT / "assets"
IMAGE = "minlag/mermaid-cli:10.9.1"
WIDTH = "1800"

BLOCK_RE = re.compile(r"^```mermaid\n(.*?)^```$", re.MULTILINE | re.DOTALL)

# Lines that frame a diagram but don't add any renderable nodes / edges.
_FRAMING_PREFIXES = (
    "classDef ",
    "flowchart",
    "graph ",
    "subgraph",
    "direction ",
    "end",
)


def has_renderable_content(block: str) -> bool:
    """Return True if the block has at least one node or edge declaration.

    Skeleton / legend blocks that contain only the init magic, the flowchart
    header, classDef declarations and comments render to a near-empty PNG;
    they should be skipped.
    """
    for raw in block.splitlines():
        line = raw.strip()
        if not line or line.startswith("%%"):
            continue
        if line.startswith(_FRAMING_PREFIXES):
            continue
        return True
    return False


def relpath(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


def render(mmd: Path, png: Path, cfg: Path) -> None:
    cmd = [
        "docker", "run", "--rm",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "-e", "HOME=/home/mermaidcli",
        "-v", f"{ROOT}:/data",
        IMAGE,
        "-p", f"/data/{relpath(cfg)}",
        "-i", f"/data/{relpath(mmd)}",
        "-o", f"/data/{relpath(png)}",
        "-w", WIDTH,
        "-b", "white",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(f"mermaid-cli failed on {mmd}")


def main() -> int:
    if shutil.which("docker") is None:
        print("docker not on PATH", file=sys.stderr)
        return 1
    try:
        subprocess.run(
            ["docker", "info"], check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        print("docker daemon is not running", file=sys.stderr)
        return 1

    TMP.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    cfg = TMP / "puppeteer-config.json"
    if not cfg.exists():
        cfg.write_text('{"args":["--no-sandbox","--disable-setuid-sandbox"]}\n')

    n = 0
    for md in sorted(DOCS.glob("*.md")):
        blocks = BLOCK_RE.findall(md.read_text())
        if not blocks:
            continue
        for i, block in enumerate(blocks, 1):
            suffix = "" if len(blocks) == 1 else f"_{i}"
            mmd = TMP / f"{md.stem}{suffix}.mmd"
            png = OUT / f"{md.stem}{suffix}.png"
            if not has_renderable_content(block):
                print(f"  {md.name} block {i}: skeleton only, skipping")
                continue
            mmd.write_text(block)
            render(mmd, png, cfg)
            print(f"  {md.name} → {png.name} ({png.stat().st_size // 1024} KB)")
            n += 1

    print(f"\nRendered {n} diagram(s) to {OUT}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
