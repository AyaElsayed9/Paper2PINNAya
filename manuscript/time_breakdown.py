import sys, time, json
from run_periodic import energy_participation
import fp_common as P
lay, beta = sys.argv[1], float(sys.argv[2])
t0 = time.perf_counter(); m, sec, L = P.case(lay, 150.0, beta, "SS"); t1 = time.perf_counter()
r = m.omega_gbt_solve(sec, L, "SS", nev=1); t2 = time.perf_counter()
r2 = m.omega_gbt_solve(sec, L, "SS", nev=2, return_mode=True); t3 = time.perf_counter()
rho, _ = energy_participation(m, sec, r2["lb"], r2["coef"]); t4 = time.perf_counter()
print(json.dumps(dict(lay=lay, beta=beta, ndof=r["ndof"], N=round(r["N"], 3), section=round(t1-t0, 1), solve_nev1=round(t2-t1, 1),
                      solve_nev2_mode=round(t3-t2, 1), participation=round(t4-t3, 1))), flush=True)
