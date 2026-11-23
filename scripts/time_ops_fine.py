"""Online operator time of the certified check (fine pre-buckling state, pre_coarsen = 1), 2 test designs per expert,
and the full GBT time of the same designs on this machine (23 Nov 2026, for the cost study)."""
import pickle, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import sample_design
from pinnbuck.families import family_context
from pinnbuck.operators import build_operator
to, tf = [], []
for side in ("symmetry", "periodic"):
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
    for lay in ("crossply", "quasiiso", "zerorich", "pm45"):
        for s in sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")[:2]:
            t = time.perf_counter(); build_operator(s, sb, lam_full=float("nan"), with_shared_reference=False, pre_coarsen=1)
            to.append(time.perf_counter() - t)
            t = time.perf_counter()
            with family_context(s) as fc:
                fc.solve(nev=3, return_mode=False)
            tf.append(time.perf_counter() - t)
            print(f"{side} {lay} {s.sid}: operators (fine pre) {to[-1]:.2f} s, full GBT {tf[-1]:.2f} s", flush=True)
print(f"mean: operators {np.mean(to):.2f} s, full GBT {np.mean(tf):.2f} s, ratio {np.mean(to) / np.mean(tf):.2f}")
