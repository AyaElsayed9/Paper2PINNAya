"""Long periodic cells (local -> global transition): N_cr and energy participation for
beta = 3 ... 20, b_s = 150 mm, all lay-ups, SS ends.  Sequential (memory), resumable."""
import json, sys, time, resource
from pathlib import Path
from run_periodic import one
F = Path(__file__).parent / "results" / "long_curves.jsonl"
LAY = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
BETAS = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0]
BC = sys.argv[1] if len(sys.argv) > 1 else "SS"
done = set()
if F.exists():
    for l in F.read_text().splitlines():
        d = json.loads(l); done.add((d["layup"], d["bs"], d["beta"], d["bc"]))
for beta in BETAS:
    for lay in LAY:
        k = (lay, 150.0, beta, BC)
        if k in done: continue
        d = one(k)
        d["maxrss_MB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        with F.open("a") as f: f.write(json.dumps(d) + "\n")
        print(lay, beta, BC, round(d["N"], 3), {g: round(100 * v, 1) for g, v in d["rho"].items()}, f'{d["time"]:.0f}s',
              f'{d["maxrss_MB"]:.0f}MB', flush=True)
