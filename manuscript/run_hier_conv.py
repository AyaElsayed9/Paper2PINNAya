"""Dense convergence curve of N_cr with the number n_C of hierarchical conventional modes
(all shear and Mindlin modes kept, pre-buckling state of the complete basis) -- Figure 10."""
import sys, json, os
from multiprocessing import Pool
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import omega_length, omega_panel
from pinnbuck.solver import configure
cases = [("Isotropic",150.0,1.5,"SS"),("Cross-ply",150.0,1.5,"SS"),("Sym2",150.0,1.5,"SS"),("Unsym",150.0,1.5,"SS"),("Unsym",150.0,1.5,"CC"),("Sym1",285.5,1.5,"CC")]
NC = [4, 8, 12, 16, 20, 25, 30, 40, 50, 60, 80, 100, 125, 150, 200]
def one(c):
    lay, bs, beta, bc = c
    with configure(OMEGA_SIDE_BC="periodic") as m:
        if lay == "Isotropic":
            pan = m.make_omega_panel(bs, "iso", [], []); L = beta * bs
        else:
            s = validated_omega_sample(lay, bs, beta, bc, "periodic"); pan = omega_panel(s); L = omega_length(s)
        sec = m.omega_gbt_section(pan); full = m.omega_gbt_solve(sec, L, bc); hm = m.hierarchical_modes(sec)
        res = dict(full=full["N"], nip=sec.nip, nm=sec.nm, nc={})
        for nc in NC:
            if nc >= sec.nip: continue
            r = m.truncated_section(sec, nc, None, None, None, hm=hm)
            res["nc"][nc] = m.omega_gbt_solve(r, L, bc, pre=full["pre"])["N"]
        print(lay, bc, res, flush=True)
        return f"{lay}|{bs}|{beta}|{bc}", res
if __name__ == "__main__":
    out = {}
    for c in cases:
        k, r = one(c); out[k] = r
        json.dump(out, open("results/hier_conv.json", "w"), indent=1)
