"""Time the four long cells of Table 7 on this computer (same run as run_long1.py:
cross-section analysis + 2 eigenvalues with mode + energy participation).

Usage (from the manuscript folder):   python time_table7.py            # 2 threads, 3 repeats
                                       python time_table7.py 2 1        # threads, repeats
Results are printed and saved to results/pc_timing.json. N_cr must equal the paper values;
only the time depends on the computer."""
import os, sys
THREADS = sys.argv[1] if len(sys.argv) > 1 else "2"            # 2 cores, as in the Abaqus runs
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = THREADS                                     # must be set before numpy is imported
import json, time, platform
from pathlib import Path
import numpy as np
from run_periodic import energy_participation
import fp_common as P

REPEATS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
CELLS = [("Isotropic", 8.0), ("Sym1", 8.0), ("Unsym", 12.0), ("Cross-ply", 14.0)]
try:
    import pypardiso  # noqa: F401
    solver = "pypardiso"
except ImportError:
    solver = "SuperLU (pip install pypardiso for faster solves)"
print(f"threads={THREADS}  repeats={REPEATS}  sparse solver: {solver}  cpu: {platform.processor() or platform.machine()}", flush=True)

out = []
for lay, beta in CELLS:
    ts = []
    for k in range(REPEATS):
        t0 = time.perf_counter()
        m, sec, L = P.case(lay, 150.0, beta, "SS")
        r = m.omega_gbt_solve(sec, L, "SS", nev=2, return_mode=True)
        energy_participation(m, sec, r["lb"], r["coef"])
        ts.append(time.perf_counter() - t0)
        print(f"  {lay:9s} beta={beta:4.1f}  run {k + 1}: N_cr = {r['N']:.4f} N/mm  ndof = {r['ndof']}  time = {ts[-1]:.1f} s", flush=True)
    out.append(dict(layup=lay, beta=beta, N=r["N"], ndof=r["ndof"], times=ts, median=float(np.median(ts))))
    print(f"{lay:9s} beta={beta:4.1f}  median time = {np.median(ts):.1f} s", flush=True)

res = dict(threads=int(THREADS), solver=solver, cpu=platform.processor(), python=platform.python_version(),
           numpy=np.__version__, cells=out)
(Path(__file__).parent / "results" / "pc_timing.json").write_text(json.dumps(res, indent=1))
print("saved results/pc_timing.json")
