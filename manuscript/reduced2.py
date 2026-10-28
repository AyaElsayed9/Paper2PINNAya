import sys, json
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import omega_length, omega_panel
from pinnbuck.solver import configure
cases = [("Isotropic",150.0,1.5,"SS"),("Cross-ply",150.0,1.5,"SS"),("Sym2",150.0,1.5,"SS"),("Unsym",150.0,1.5,"SS"),("Unsym",150.0,1.5,"CC"),("Sym1",285.5,1.5,"CC")]
out = {}
for lay, bs, beta, bc in cases:
    with configure(OMEGA_SIDE_BC="periodic") as m:
        if lay == "Isotropic":
            pan = m.make_omega_panel(bs, "iso", [], []); L = beta * bs
        else:
            s = validated_omega_sample(lay, bs, beta, bc, "periodic"); pan = omega_panel(s); L = omega_length(s)
        sec = m.omega_gbt_section(pan); full = m.omega_gbt_solve(sec, L, bc); hm = m.hierarchical_modes(sec)
        res = dict(full=[full["N"], sec.nm], nip=sec.nip)
        for nc in (20, 40, 60, 100, 150):
            r = m.truncated_section(sec, nc, None, None, None, hm=hm)
            N = m.omega_gbt_solve(r, L, bc, pre=full["pre"])["N"]; res[f"c{nc}"] = [N, r.nm]
        r = m.truncated_section(sec, None, None, 0, 0, hm=hm)
        N = m.omega_gbt_solve(r, L, bc, pre=full["pre"])["N"]; res["noM"] = [N, r.nm]
        out[f"{lay}|{bs}|{beta}|{bc}"] = res
        print(lay, bs, beta, bc, {k: v for k, v in res.items()}, flush=True)
json.dump(out, open("results/reduced2.json", "w"), indent=1)
