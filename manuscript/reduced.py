import sys, time, json
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import omega_length, omega_panel
from pinnbuck.solver import configure
import numpy as np
cases = [("Isotropic",285.5,0.75,"SS"),("Cross-ply",150.0,1.5,"SS"),("Unsym",150.0,1.5,"SS"),("Sym2",150.0,1.5,"SS"),("Unsym",150.0,1.5,"CC")]
out = {}
for lay, bs, beta, bc in cases:
    with configure(OMEGA_SIDE_BC="periodic") as m:
        if lay == "Isotropic":
            pan = m.make_omega_panel(bs, "iso", [], []); L = beta * bs
        else:
            s = validated_omega_sample(lay, bs, beta, bc, "periodic"); pan = omega_panel(s); L = omega_length(s)
        sec = m.omega_gbt_section(pan)
        full = m.omega_gbt_solve(sec, L, bc, return_mode=True)
        hm = m.hierarchical_modes(sec)
        res = {"full": (full["N"], sec.nm)}
        for nc in (20, 40, 60, 100, 150):
            for mode in ("a", "b3"):
                kw = dict(coupled_shear=True, kfac=(1.0,) if mode == "a" else (0.5, 1.0, 2.0))
                r = m.truncated_section(sec, nc, 0, 0, 0, hm=hm, **kw)
                N = m.omega_gbt_solve(r, L, bc, pre=full["pre"])["N"]
                res[f"{mode}_{nc}"] = (N, r.nm, 100 * (N - full["N"]) / full["N"])
        out[f"{lay}|{bs}|{beta}|{bc}"] = res
        print(lay, bs, beta, bc, {k: (round(v[0], 3), v[1], round(v[2], 2) if len(v) > 2 else None) for k, v in res.items()}, flush=True)
json.dump(out, open("results/reduced.json", "w"), indent=1)
