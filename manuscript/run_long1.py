"""One long periodic cell (b_s = 150 mm, SS) in its own process (memory): N_cr and energy
participation; pss=0 drops the transverse pre-buckling force S_ss from G (mechanism check).
Appends to results/long_curves.jsonl."""
import json, sys, time, resource
from pathlib import Path
import numpy as np
from run_periodic import energy_participation
import fp_common as P
lay, beta, pss = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
hx = float(sys.argv[4]) if len(sys.argv) > 4 else None
t0 = time.perf_counter()
m, sec, L = P.case(lay, 150.0, beta, "SS")
r = m.omega_gbt_solve(sec, L, "SS", hx=hx, nev=2, return_mode=True, pss_scale=pss)
rho, rho_self = energy_participation(m, sec, r["lb"], r["coef"])
d = dict(layup=lay, bs=150.0, beta=beta, L=L, bc="SS", pss=pss, hx=hx, N=r["N"], eigs=r["eigs"], ndof=r["ndof"], nmodes=r["nmodes"],
         rho=rho, p_disp=r["participation"], time=time.perf_counter() - t0,
         maxrss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
with (Path(__file__).parent / "results" / "long_curves.jsonl").open("a") as f:
    f.write(json.dumps(d) + "\n")
print(lay, beta, "pss", pss, round(d["N"], 3), {g: round(100 * v, 1) for g, v in rho.items()}, f'{d["time"]:.0f}s', f'{d["maxrss_MB"]:.0f}MB', flush=True)
