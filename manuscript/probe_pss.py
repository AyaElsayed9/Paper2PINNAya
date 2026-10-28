import sys, time
sys.argv = sys.argv[:1] + sys.argv[1:]
import fp_common as P
lay, beta = sys.argv[1], float(sys.argv[2])
m, sec, L = P.case(lay, 150.0, beta, "SS")
for ps in (1.0, 0.0):
    t = time.time(); r = m.omega_gbt_solve(sec, L, "SS", pss_scale=ps)
    print(lay, beta, "pss_scale", ps, "N", round(r["N"], 3), "eigs", [round(e, 3) for e in r["eigs"]], f"{time.time()-t:.0f}s", flush=True)
