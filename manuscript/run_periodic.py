"""Periodic-only GBT runs for the manuscript: N_cr over the aspect ratio for all
lay-ups, pitches and end conditions, with the energy-based modal participation
(template Eqs. 50-51: rho_g = U_g / sum U, U_g = energy of family g alone).
Appends one JSON line per case to results/periodic_curves.jsonl (resumable)."""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample          # noqa: E402
from pinnbuck.oracle import omega_length, omega_panel             # noqa: E402
from pinnbuck.solver import configure                             # noqa: E402

OUT = Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
FILE = OUT / "periodic_curves.jsonl"
LAYUPS = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
BETAS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]


def energy_participation(m, sec, lb, C):
    """Energy-based participation of the families.
    rho_g  = a_g^T (K a) / a^T K a   (additive energy partition, sums to 1)
    self_g = a_g^T K_gg a_g / sum_h a_h^T K_hh a_h  (self energies only)."""
    Kcs = m._gbt_Kcs(sec)
    X = {(d, e): (lb.B[d].T @ sp.diags(lb.wx) @ lb.B[e]).toarray() for d in range(3) for e in range(3)}
    grp = np.array(["S" if g == "A" else g for g in sec.group])
    KA = np.zeros_like(C)
    for de, A in Kcs.items():
        KA += A @ C @ X[de].T
    tot = float(np.sum(C * KA))
    rho, selfe = {}, {}
    for g in ("G", "D", "L", "T", "S", "M"):
        k = np.nonzero(grp == g)[0]
        if not len(k):
            continue
        rho[g] = float(np.sum(C[k] * KA[k])) / tot
        Cg = C[k]
        selfe[g] = sum(float(np.sum(A[np.ix_(k, k)] * (Cg @ X[de] @ Cg.T))) for de, A in Kcs.items())
    st = sum(selfe.values())
    return rho, {g: v / st for g, v in selfe.items()}


def one(args):
    lay, bs, beta, bc = args
    t0 = time.perf_counter()
    with configure(OMEGA_SIDE_BC="periodic") as m:
        if lay == "Isotropic":
            pan = m.make_omega_panel(bs, "iso", [], [])
            L = beta * bs
        else:
            s = validated_omega_sample(lay, bs, beta, bc, "periodic")
            pan = omega_panel(s); L = omega_length(s)
        sec = m.omega_gbt_section(pan)
        r = m.omega_gbt_solve(sec, L, bc, nev=2, return_mode=True)
        rho, rho_self = energy_participation(m, sec, r["lb"], r["coef"])
    return dict(layup=lay, bs=bs, beta=beta, L=L, bc=bc, N=r["N"], eigs=r["eigs"], ndof=r["ndof"],
                nmodes=r["nmodes"], rho=rho, rho_self=rho_self, p_disp=r["participation"], time=time.perf_counter() - t0)


if __name__ == "__main__":
    done = set()
    if FILE.exists():
        for ln in FILE.read_text().splitlines():
            d = json.loads(ln); done.add((d["layup"], d["bs"], d["beta"], d["bc"]))
    todo = [(l, bs, b, bc) for bc in ("SS", "CC") for bs in (150.0, 285.5) for l in LAYUPS for b in BETAS
            if (l, bs, b, bc) not in done]
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        todo = todo[:1]
    print(len(todo), "cases", flush=True)
    with Pool(int(os.environ.get("NPROC", "2"))) as pool, FILE.open("a") as f:
        for d in pool.imap_unordered(one, todo):
            f.write(json.dumps(d) + "\n"); f.flush()
            print(d["layup"], d["bs"], d["beta"], d["bc"], round(d["N"], 4), f'{d["time"]:.1f}s', flush=True)
