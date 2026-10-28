import sys, time, json
import fp_common as P
lay, beta, hx = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
m, sec, L = P.case(lay, 150.0, beta, "SS")
t = time.time(); r = m.omega_gbt_solve(sec, L, "SS", hx=hx)
print(json.dumps(dict(layup=lay, beta=beta, hx=hx, N=r["N"], eigs=r["eigs"], t=time.time() - t)), flush=True)
