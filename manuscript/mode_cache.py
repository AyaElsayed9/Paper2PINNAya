"""GBT critical mode of one SS/CC cell (b_s = 150 mm) -> figp/longmode_<lay>_<beta>_<pss>.pkl (same format as figures_paper)."""
import sys, pickle
import numpy as np
import fp_common as P
lay, beta, pss = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
hx = float(sys.argv[4]) if len(sys.argv) > 4 else None
f = P.CACHE / f"longmode_{lay}_{beta}_{pss}.pkl"
if not f.exists():
    m, sec, L = P.case(lay, 150.0, beta, "SS")
    r = m.omega_gbt_solve(sec, L, "SS", hx=hx, nev=1, return_mode=True, pss_scale=pss)
    xs = np.linspace(0, L, 161)
    U, u = P.nodal_disp(sec, r["lb"], np.asarray(r["coef"]), xs)
    d = dict(L=L, xs=xs, U=U, u=u, N=r["N"], X=sec.X.copy(),
             S=[(s_["a"], s_["b"], sec.wpart[s_["w"]], s_["n"]) for s_ in sec.S], sides=sec.sides)
    f.write_bytes(pickle.dumps(d))
print(lay, beta, pss, pickle.loads(f.read_bytes())["N"])
