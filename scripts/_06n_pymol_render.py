"""PyMOL render helper for 06n (runs under the `gradi-pymol` env, invoked by the orchestrator):

    pymol -cq scripts/_06n_pymol_render.py -- <job.json> <outdir>

<job.json> = {"jobs": [{"acc","label","cif","consensus","plddt","pocket":[resi,...]}, ...]}.
Renders one ray-traced cartoon PNG per target: AlphaFold model coloured by per-residue pLDDT
(official palette) with the top P2Rank pocket as green sticks + a translucent green surface.
"""
from pymol import cmd
import json, os, sys

args = sys.argv[1:]
job, outdir = args[0], args[1]
data = json.load(open(job))
os.makedirs(outdir, exist_ok=True)

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("surface_quality", 1)
# official AlphaFold pLDDT palette
cmd.set_color("afvl", [1.00, 0.490, 0.271])   # <50 orange
cmd.set_color("afl",  [1.00, 0.859, 0.075])   # 50-70 yellow
cmd.set_color("afc",  [0.396, 0.796, 0.953])  # 70-90 cyan
cmd.set_color("afvh", [0.00, 0.325, 0.839])   # 90+ blue
cmd.set_color("pkgreen", [0.0, 0.627, 0.529])  # #00A087 pocket

for i, j in enumerate(data["jobs"]):
    cmd.delete("all")
    cmd.load(j["cif"], "m")
    cmd.hide("everything", "m")
    cmd.show("cartoon", "m")
    cmd.color("afvl", "m")                # very-low (orange); higher pLDDT bands overwrite
    cmd.color("afl", "m and b>50")
    cmd.color("afc", "m and b>70")
    cmd.color("afvh", "m and b>90")
    if j.get("pocket"):
        resi = "+".join(str(r) for r in j["pocket"])
        cmd.select("pock", f"m and resi {resi}")
        cmd.show("sticks", "pock and not (name C+N+O)")
        cmd.color("pkgreen", "pock")
        cmd.set("stick_radius", 0.22)
        cmd.show("surface", "pock")
        cmd.set("transparency", 0.5)
    cmd.orient("m")
    cmd.zoom("m", buffer=3)
    cmd.ray(1000, 800)
    cmd.png(os.path.join(outdir, f"{i:02d}_{j['acc']}.png"), dpi=150)
    print("rendered", j["label"])
