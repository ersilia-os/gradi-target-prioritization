"""Render the link-preview card (webapp/social-preview.png, 1200x630).

Counts are read from the committed payloads, never typed, so the card cannot
disagree with the site. Text is sized to stay legible when a chat app shrinks the
card to ~300 px wide, and everything sits inside a 90 px safe margin because phone
apps crop the edges of a 1.91:1 image. Env: gradi (needs Pillow and Inter).
"""
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = REPO_ROOT / "webapp" / "social-preview.png"
FONTS = Path.home() / "Library" / "Fonts"
W, H, M = 1200, 630, 90
SCALE = 2  # draw at 2x, downsample for clean edges

INK, MUTED, BG, BORDER = "#2c3e50", "#6b6675", "#fafafc", "#e6e6ee"
STRIP = ["#fcbf49", "#e63946", "#457b9d", "#6bbf59", "#6c5ce7", "#2ec4b6", "#b05cc8"]
ORG = {"ec": ("E. coli", "#f4845f"), "kp": ("K. pneumoniae", "#457b9d")}


def font(weight, size):
    return ImageFont.truetype(str(FONTS / f"Inter_28pt-{weight}.ttf"), size * SCALE)


def main():
    counts = {k: json.load(open(REPO_ROOT / "webapp" / "data" / f"{k}.json"))["n"] for k in ORG}
    im = Image.new("RGB", (W * SCALE, H * SCALE), BG)
    d = ImageDraw.Draw(im)
    s = SCALE

    seg = W * s / len(STRIP)
    for i, c in enumerate(STRIP):
        d.rectangle([round(i * seg), 0, round((i + 1) * seg), 10 * s], fill=c)

    d.text((M * s, 120 * s), "Gram-negative target browser", font=font("SemiBold", 72), fill=INK)
    d.text((M * s, 205 * s), "for protein degradation", font=font("SemiBold", 72), fill=INK)
    d.text((M * s, 320 * s), "Proteome-wide, evidence-graded ranking of candidate targets",
           font=font("Regular", 32), fill=MUTED)

    x = M
    for k, (name, col) in ORG.items():
        label = f"{name}  {counts[k]:,} proteins"
        f = font("Medium", 34)
        w = d.textlength(label, font=f) / s
        d.rounded_rectangle([x * s, 405 * s, (x + w + 64) * s, 473 * s], radius=10 * s,
                            fill="#fff", outline=BORDER, width=2 * s)
        d.ellipse([(x + 24) * s, 429 * s, (x + 44) * s, 449 * s], fill=col)
        d.text(((x + 56) * s, 439 * s), label, font=f, fill=INK, anchor="lm")
        x += w + 64 + 20

    d.line([M * s, 520 * s, (W - M) * s, 520 * s], fill=BORDER, width=2 * s)
    d.text((M * s, 568 * s), "GraDi  ·  Ersilia Open Source Initiative", font=font("Medium", 28),
           fill=MUTED, anchor="lm")

    im = im.resize((W, H), Image.LANCZOS)
    im.save(OUT, optimize=True)
    print(f"wrote {OUT.relative_to(REPO_ROOT)} {W}x{H} {OUT.stat().st_size/1024:.0f} KB; n = {counts}")


if __name__ == "__main__":
    main()
