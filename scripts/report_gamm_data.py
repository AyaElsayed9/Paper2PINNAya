"""Numbers and chart data for the GAMM 2027 results report (25 Nov 2026); all from logged result files."""
import csv, glob, json
from pathlib import Path
import numpy as np
out = {}
tags = ["symmetry_crossply", "periodic_crossply", "symmetry_quasiiso", "periodic_quasiiso",
        "symmetry_zerorich", "periodic_zerorich", "symmetry_pm45", "periodic_pm45"]
ex = []
lab = {}
pool5 = []
for t in tags:
    R = list(csv.DictReader(open(f"results/error_indicator/{t}_SS_v2.csv")))
    T = {r["sid"]: float(r["N_full"]) for r in R if r["N_full"]}
    row = {"expert": t}
    for m in sorted({r["model"] for r in R}):
        rr = [r for r in R if r["model"] == m and r["sid"] in T]
        e = lambda k: np.array([100 * abs(float(r[k]) / T[r["sid"]] - 1) for r in rr])
        h = e("N_head")
        d = {"head_mean": h.mean(), "head_max": h.max()}
        if not m.startswith("gp"):
            for s in (0, 5, 12):
                d[f"s{s}_mean"] = e(f"N{s}").mean(); d[f"s{s}_max"] = e(f"N{s}").max()
            if m in ("physics_L0", "hybrid_L120"):
                pool5 += [(float(r["eta5"]), float(r["N5"]) / T[r["sid"]] - 1, m) for r in rr]
        row[m] = d
        if t == "symmetry_crossply":
            lab[m] = d
    ex.append(row)
out["experts"] = ex
out["labels_crossply_symmetry"] = lab
# convergence from the cost study traces (value after stop carried forward)
conv = {}
for st in ("random", "nn10", "nn120", "net_L0", "net_L120"):
    E = []
    for t in tags:
        for r in csv.DictReader(open(f"results/cost_study/{t}_SS_v2.csv")):
            if r["start"] != st: continue
            N = json.loads(r["N"]); Nf = float(r["N_full"])
            N = (N + [N[-1]] * 13)[:13]
            E.append([100 * (v / Nf - 1) for v in N])
    E = np.array(E); conv[st] = {"median": np.median(E, 0).round(3).tolist(), "p90": np.quantile(E, 0.9, 0).round(3).tolist(), "n": len(E)}
out["convergence"] = conv
out["estimate_step5"] = [[round(100 * a, 4), round(100 * b, 4), m] for a, b, m in pool5]
out["cost_tol5"] = list(csv.DictReader(open("results/cost_study/summary_tol5.csv")))
json.dump(out, open("results/report_gamm/report_data.json", "w"), indent=1, default=float)
for r in ex:
    a, b = r["physics_L0"], r["hybrid_L120"]
    print(f"{r['expert']:18s} L0 head {a['head_mean']:5.2f}/{a['head_max']:5.2f} check12 {a['s12_mean']:.2f}/{a['s12_max']:.2f} | "
          f"L120 head {b['head_mean']:.2f}/{b['head_max']:.2f} check12 {b['s12_mean']:.2f}/{b['s12_max']:.2f} | GP120 {r['gp_L120']['head_mean']:.2f}/{r['gp_L120']['head_max']:.2f}")
print({m: (round(d['head_mean'], 2), round(d['head_max'], 2)) for m, d in lab.items()})
for st, c in conv.items(): print(st, c["median"][:13])
