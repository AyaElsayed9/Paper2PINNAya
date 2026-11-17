"""Fixed cost of the check split finer (17 Nov 2026): section, pre-buckling (x4 coarse), other operators,
Ritz space, energy metric (block Cholesky), one LOBPCG step. Warm-up design not counted."""
import sys, time, pickle
from pathlib import Path
import numpy as np, torch
ROOT = Path("/home/claude/paper2pinnaya")
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio
from pinnbuck.designspace import sample_design
from pinnbuck.families import family_context
from pinnbuck.operators import build_operator
from pinnbuck.ritzspace import RitzSpace, EnergyMetric
from pinnbuck.dictionary import dictionary_kw
from pinnbuck.spectral import lobpcg
torch.set_num_threads(2)
tag = sys.argv[1]; side, lay = tag.split("_")[:2]
S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")
sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
epinn.INPUT_COLS = [0, 1]
rows = []
for i, s in enumerate(S[:int(sys.argv[2]) + 1]):
    T = {}
    t = time.perf_counter()
    with family_context(s) as fc:
        T["section"] = time.perf_counter() - t; t = time.perf_counter()
        fc.pre(coarsen=4); T["pre_x4"] = time.perf_counter() - t; t = time.perf_counter()
        r = fc.solve(nev=3, return_mode=False); T["full_gbt"] = time.perf_counter() - t
    t = time.perf_counter()
    op = build_operator(s, sb, lam_full=float("nan"), with_shared_reference=False, pre_coarsen=4)
    op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32); op.info["dictionary"] = "d2"
    T["build_operator"] = time.perf_counter() - t; t = time.perf_counter()
    rs = RitzSpace(op, raw=True, **dictionary_kw(op)); T["ritz_space"] = time.perf_counter() - t; t = time.perf_counter()
    em = EnergyMetric(rs, groups="labels", ridge=1e-8); T["metric"] = time.perf_counter() - t
    stamps = [time.perf_counter()]
    lobpcg(rs, R=8, iters=6, X0=None, precond=em.solve, tol=0, seed=1, callback=lambda it, mu: stamps.append(time.perf_counter()))
    T["step"] = (stamps[-1] - stamps[1]) / (len(stamps) - 2)
    T["other_ops"] = T["build_operator"] - T["section"] - T["pre_x4"]
    if i == 0:
        continue                                            # warm-up
    rows.append(T)
    print(f"{s.sid} bs {s.params['bs']:.0f} beta {s.params['beta']:.2f} nz {len(op.z)} sizes {len(em.sizes)} blocks max {max(em.sizes)}: " +
          "  ".join(f"{k} {v:.2f}" for k, v in T.items()), flush=True)
    del op, rs, em
print("mean: " + "  ".join(f"{k} {np.mean([r[k] for r in rows]):.2f}" for k in rows[0]))
