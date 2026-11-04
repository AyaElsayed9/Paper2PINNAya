"""GBT (frozen solver V32, standard settings) for every Abaqus case of cases.csv, computed and committed
BEFORE the Abaqus results exist, so that the comparison is blind.
    python gbt_blind.py   -> gbt_predictions.csv"""
import csv, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import validated_omega_sample   # noqa: E402
from pinnbuck.oracle import solve                         # noqa: E402

here = Path(__file__).parent
out = here / "gbt_predictions.csv"
done = {(r["set"], r["bs"], r["beta"], r["layup"]) for r in csv.DictReader(open(out))} if out.exists() else set()
new = not out.exists()
with open(out, "a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(["set", "sides", "bc", "bs", "beta", "layup", "N_GBT", "m_GBT", "t_s"])
    for r in csv.DictReader(open(here / "cases.csv")):
        if (r["set"], r["bs"], r["beta"], r["layup"]) in done:
            continue
        side = "periodic" if r["sides"] == "PBC" else "symmetry"
        t0 = time.time()
        s = solve(validated_omega_sample(r["layup"], float(r["bs"]), float(r["beta"]), r["bc"], side), "gbt", "standard")
        w.writerow([r["set"], r["sides"], r["bc"], r["bs"], r["beta"], r["layup"],
                    f"{s['value']:.5f}" if s.get("ok") else "FAILED", s.get("m"), round(time.time() - t0, 1)])
        f.flush(); print(r, s.get("value"), s.get("m"), flush=True)
