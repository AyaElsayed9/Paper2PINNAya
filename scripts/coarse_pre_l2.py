"""How much does the coarse pre-buckling state change full-GBT N_cr on the 45 level-2 test designs? (10 Nov 2026)
For coarsen k = 2, 4: time of the coarse pre-buckling state and dN = N(full GBT, coarse pre) / N(full GBT) - 1.
results/fast_check/level2_coarse_pre.csv"""
import csv, glob, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import opio
from pinnbuck.designspace import Sample
from pinnbuck.families import family_context
full = {r["sid"]: float(r["N_full"]) for r in csv.DictReader(open(ROOT / "results/fast_check/level2_industrial.csv"))}
out = ROOT / "results/fast_check/level2_coarse_pre.csv"
done = {r["sid"] for r in csv.DictReader(open(out))} if out.exists() else set()
with open(out, "a", newline="") as fh:
    w = csv.writer(fh)
    if not done:
        w.writerow(["sid", "bs", "beta", "skin", "N_full", "t_pre_fine", "t_pre2", "dN2_pct", "t_pre4", "dN4_pct"])
    for f in sorted(glob.glob(str(ROOT / "data/operators/periodic_industrial_SS_v2/test/*.pkl"))):
        sd = opio.load(f).info["sample"]
        if sd["sid"] in done:
            continue
        s = Sample(family=sd["family"], params=sd["params"], layups=sd["layups"], fixed=sd["fixed"], space=sd["space"])
        row = [s.sid, s.params["bs"], s.params["beta"], "/".join(f"{a:g}" for a in sd["layups"]["skin"]), full[s.sid]]
        with family_context(s) as fc:
            t = time.perf_counter(); fc.pre(); row.append(time.perf_counter() - t)
            for k in (2, 4):
                t = time.perf_counter(); pre = fc.pre(coarsen=k); row.append(time.perf_counter() - t)
                N = float(fc.m.omega_gbt_solve(fc.sec, fc.L, fc.bc, hx=fc.hx, nev=3, pre=pre)["N"])
                row.append(100 * (N / full[s.sid] - 1))
        w.writerow(row); fh.flush(); print(row, flush=True)
