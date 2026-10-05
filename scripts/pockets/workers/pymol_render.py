"""PyMOL worker: one ray-traced cartoon per protein, with its top pocket highlighted.

    pymol -cq scripts/pockets/workers/pymol_render.py -- <job.json> <outdir>

**Runs under `gradi-pymol`, never imported.** `from pymol import cmd` only exists in that env, so
this file is a worker by the layout rule: an entry point under a different interpreter. The caller
is `plotting/cavities.py`, which reaches across the process boundary with `conda run`.

job.json = {"jobs": [{"acc", "label", "cif", "pocket": [resi, ...]}, ...]}

Adapted from the v1 renderer `legacy/scripts/_06n_pymol_render.py`, which is frozen. The colouring
is the official AlphaFold per-residue pLDDT palette, so a reader who has seen an AlphaFold page
reads the confidence here the same way, with the pocket over the top in npg teal.

**The model is coloured by pLDDT, which is NOT pocket confidence.** An orange loop means the fold
there is uncertain; it says nothing about the pocket. That distinction belongs in the caption.
"""
import json
import os
import sys

from pymol import cmd

args = sys.argv[1:]
if len(args) < 2:
    raise SystemExit("usage: pymol -cq pymol_render.py -- <job.json> <outdir>")
job_path, outdir = args[0], args[1]
with open(job_path) as fh:
    data = json.load(fh)
os.makedirs(outdir, exist_ok=True)

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("surface_quality", 1)
# Without this, anything the near clipping plane cuts renders as a BLACK hole, which reads as a
# void in the protein. Measured while choosing this idiom: a tight zoom on the pocket produced
# large black regions that looked like rendering damage.
cmd.set("ray_interior_color", "grey70")
# DEPTH CUE OFF. With a COMMON camera distance (see the two-pass block below) the smaller
# proteins sit deeper in the fog range and render visibly washed out -- gyrB came out pale beside
# irp2 purely because it is smaller, which is exactly the misreading the common scale exists to
# prevent. Fog is a depth cue for one molecule, not a fair comparison across several.
cmd.set("depth_cue", 0)
cmd.set("ray_trace_fog", 0)
cmd.set("fog", 0)

# ALPHAFOLD pLDDT CARTOON + a TRANSPARENT pocket surface (owner, 2026-10-05). The cartoon uses
# the official AlphaFold per-residue palette, so a reader who has seen an AlphaFold page reads the
# model confidence here the same way.
#
# **pLDDT IS NOT POCKET CONFIDENCE.** An orange loop means the fold there is uncertain; it says
# nothing about whether the pocket is real. Keep that out of the caption's way.
#
# The pocket is MAGENTA rather than the v1 teal: against a blue/cyan pLDDT cartoon, teal is nearly
# the same hue and the pocket disappears into it. Magenta is the one colour in the npg set that
# separates cleanly from the whole pLDDT ramp.
cmd.set_color("afvl", [1.00, 0.490, 0.271])    # <50  very low
cmd.set_color("afl", [1.00, 0.859, 0.075])     # 50-70 low
cmd.set_color("afc", [0.396, 0.796, 0.953])    # 70-90 confident
cmd.set_color("afvh", [0.00, 0.325, 0.839])    # >90  very high
cmd.set_color("pocket", [0.800, 0.000, 0.470])  # npg magenta, PAL.ACCENT

# TWO PASSES, FOR A COMMON SCALE. Rendering each protein with its own `zoom` fills every frame
# equally, so an 805-aa protein comes out LARGER than a 1,048-aa one and the picture contradicts
# the lengths printed in the figure's titles. Pass 1 orients each molecule and records the camera
# distance its own zoom chose; pass 2 re-renders every molecule at the LARGEST of those distances,
# so physical size is comparable across panels and the biggest protein looks biggest.
def setup(j):
    """Load, colour, select the pocket and orient. Returns the pocket selection name or None."""
    cmd.delete("all")
    cmd.load(j["cif"], "m")
    cmd.hide("everything", "m")
    cmd.show("cartoon", "m")
    # Painted low-to-high so each band overwrites the one below it; the order matters. AlphaFold
    # models carry pLDDT in the B-factor column, so `b>` is a pLDDT test HERE -- that is true of
    # these models, not of a crystal structure.
    cmd.color("afvl", "m")
    cmd.color("afl", "m and b>50")
    cmd.color("afc", "m and b>70")
    cmd.color("afvh", "m and b>90")
    resi = j.get("pocket") or []
    if resi:
        sel = "+".join(str(r) for r in resi)
        cmd.select("pock", f"m and resi {sel}")
        # TRANSPARENT SURFACE. An opaque one interpenetrates the cartoon -- ribbon slabs cut
        # across the pocket and read as broken geometry. At 0.4 the ribbon passes through cleanly
        # and the cavity still reads as a volume rather than as painted residues.
        cmd.show("surface", "pock")
        cmd.color("pocket", "pock")
        cmd.set("transparency", 0.4, "pock")
        # Orient on the POCKET so it faces the camera; zoom on the WHOLE molecule so nothing
        # clips. Zooming on the pocket puts the camera inside the protein and the render goes
        # black where the near plane cuts.
        cmd.orient("pock")
    else:
        cmd.orient("m")
    cmd.zoom("m", 4)
    cmd.clip("slab", 500)
    return bool(resi)


# Pass 1: what distance does each molecule need?
distances = []
for j in data["jobs"]:
    setup(j)
    distances.append(cmd.get_view()[11])
far = min(distances)          # most negative = camera furthest back = the largest molecule

written = []
for j in data["jobs"]:
    setup(j)
    view = list(cmd.get_view())
    view[11] = far            # one camera distance for every panel
    cmd.set_view(view)
    out = os.path.join(outdir, f"{j['acc']}.png")
    cmd.png(out, width=1400, height=1400, dpi=300, ray=1)
    written.append(out)

cmd.delete("all")
print(f"rendered {len(written)} structures", flush=True)
