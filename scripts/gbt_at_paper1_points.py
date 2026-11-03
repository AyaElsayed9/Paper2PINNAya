"""Full GBT (solver V32, symmetry side edges, SS ends, cross-ply lay-up of Paper 1) at the 22 Abaqus
points of Paper 1 (Elsayed et al., Compos. Struct. 399 (2026) 120959, Table 2). The Abaqus and PINN
values are copied from that table (not computed here)."""
import json, sys, time
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import solve

P1 = [  # region, bs, beta, N_FEA, N_PINN(Paper 1), m_FEA
 ("Regular",150,1.69,30.837,30.865,3),("Regular",150,2.8,30.84,30.982,5),("Regular",205,2.53,13.813,13.558,4),
 ("Regular",205,2.94,13.766,13.597,5),("Regular",170,2.25,22.144,21.881,4),("Regular",285.5,2.53,6.281,6.206,4),
 ("Regular",136,2.11,39.813,41.222,4),("Regular",150,2.11,31.018,31.073,4),("Transition",205,0.86,16.459,16.018,2),
 ("Transition",205,1.41,14.346,14.021,2),("Transition",150,0.86,34.389,34.479,2),("Transition",150,1.27,31.605,31.562,2),
 ("Transition",170,1.41,23.529,22.937,2),("Transition",285.5,1.41,6.425,6.364,2),("Transition",150,1.83,31.178,31.453,3),
 ("Transition",205,2.11,14.145,13.87,4),("Boundary",124.25,0.72,55.099,55.398,1),("Boundary",124.25,1.55,49.43,54.768,3),
 ("Boundary",124.25,2.53,49.642,54.763,5),("Boundary",124.25,3.5,49.714,55.014,6),("Boundary",285.5,3.5,6.34,6.256,6),
 ("Boundary",205,3.5,13.775,13.689,6)]
out = sys.argv[1]
for reg, bs, beta, nf, npinn, mf in P1:
    t0 = time.time()
    r = solve(validated_omega_sample("Cross-ply", float(bs), float(beta), "SS", "symmetry"), "gbt", "standard")
    rec = dict(region=reg, bs=bs, beta=beta, N_FEA=nf, N_PINN_P1=npinn, m_FEA=mf, ok=r["ok"],
               N_GBT=r.get("value"), m_GBT=r.get("m"), err=r.get("error"), t=round(time.time() - t0, 1))
    open(out, "a").write(json.dumps(rec) + "\n"); print(rec, flush=True)
