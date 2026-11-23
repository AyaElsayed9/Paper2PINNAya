"""Summary of the 5 % cost study (scripts/cost_study.py) with the GP and full-GBT baselines.

Accuracy target: every delivered N_cr within TOL of full GBT.
  full GBT          one full solve per design (time: measured full-GBT times of the same designs, fast_check csv)
  GP, NL labels     mean of the GP (inputs b_s, beta; kernel of expert_tables.py) on the first NL labelled designs;
                    'GP + trigger': designs whose +-2 std band is wider than TOL get a full GBT solve
  check from a start (random, nn10, nn120, net_L0, net_L120): stop when the label-free estimate eta <= TOL; if not
                    reached in 40 steps -> full GBT for that design
Costs: full GBT solves (labels offline + solves online), offline time (labels at the measured full-GBT time,
network training measured from the logs, operators of the training designs), online time per design (operators of
the design + steps).  Break-even = number of designs from which a method is cheaper in total time than full GBT.
"""
import csv
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "scripts")
TOL = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
T_OPS = float(sys.argv[2]) if len(sys.argv) > 2 else float("nan")     # online operator time per design (fine pre), s
T_OPS_TRAIN = T_OPS                                                     # operators of one training design, s


def train_hours(tag, m):
    v = [int(x) for x in re.findall(r"^outer\s+\d+:.*\((\d+) s\)", open(f"runs/experts/{tag}_SS_v2_{m}.log").read(), re.M)]
    return sum(a for a, b in zip(v, v[1:] + [0]) if b < a) / 3600


def gp(tag, nl, X):
    from error_indicator import gp_mean_std
    return gp_mean_std(f"{tag}_SS_v2", nl, X)


rows_out = []
for f in sorted(glob.glob("results/cost_study/*_SS_v2.csv")):
    tag = Path(f).stem.replace("_SS_v2", "")
    R = list(csv.DictReader(open(f)))
    fc = {r["sid"]: r for r in csv.DictReader(open(f"results/fast_check/{tag}_SS_v2.csv"))}
    sids = sorted({r["sid"] for r in R})
    if any(sum(1 for r in R if r["sid"] == s) < 5 for s in sids) or len(sids) < 40:
        print(f"{tag}: incomplete ({len(sids)} designs), skipped"); continue
    Nf = {r["sid"]: float(r["N_full"]) for r in R}
    t_full = np.mean([float(fc[s]["t_full"]) for s in sids])
    t_prep = np.mean([float(fc[s]["t_prep"]) for s in sids])
    n = len(sids)
    res = [dict(expert=tag, method="full GBT", calls=n, offline_s=0.0, online_s=t_full, err_mean=0.0, err_max=0.0, guaranteed="yes")]
    X = np.array([[float(next(r for r in R if r["sid"] == s)["bs"]), float(next(r for r in R if r["sid"] == s)["beta"])] for s in sids])
    N = np.array([Nf[s] for s in sids])
    for nl in (10, 20, 40, 80, 120):
        m, sd = gp(tag, nl, X)
        e = 100 * np.abs(m / N - 1)
        res.append(dict(expert=tag, method=f"GP {nl} labels", calls=nl, offline_s=nl * t_full, online_s=0.0,
                        err_mean=e.mean(), err_max=e.max(), guaranteed="no"))
        trig = 2 * sd > np.log(1 + TOL)
        e2 = np.where(trig, 0.0, e)
        res.append(dict(expert=tag, method=f"GP {nl} + trigger", calls=nl + int(trig.sum()), offline_s=nl * t_full,
                        online_s=trig.mean() * t_full, err_mean=e2.mean(), err_max=e2.max(), guaranteed="no"))
    off = {"random": 0.0, "nn10": 10 * t_full, "nn120": 120 * t_full,
           "net_L0": 152 * T_OPS_TRAIN + 3600 * train_hours(tag, "physics_L0"),
           "net_L120": 152 * T_OPS_TRAIN + 120 * t_full + 3600 * train_hours(tag, "hybrid_L120")}
    lab = {"random": 0, "nn10": 10, "nn120": 120, "net_L0": 0, "net_L120": 120}
    for st in ("random", "nn10", "nn120", "net_L0", "net_L120"):
        errs, online, calls, steps = [], [], 0, []
        for s in sids:
            r = next(x for x in R if x["sid"] == s and x["start"] == st)
            Ns, etas, ts = json.loads(r["N"]), json.loads(r["eta"]), json.loads(r["t"])
            k = next((i for i, v in enumerate(etas) if v <= TOL), None)
            if k is None:
                calls += 1; errs.append(0.0); online.append(T_OPS + t_prep + float(r["t_start"]) + ts[-1] + t_full)
                steps.append(len(ts) - 1)
            else:
                errs.append(100 * (Ns[k] / Nf[s] - 1)); online.append(T_OPS + t_prep + float(r["t_start"]) + ts[k]); steps.append(k)
        errs = np.array(errs)
        res.append(dict(expert=tag, method=f"check from {st}", calls=lab[st] + calls, offline_s=off[st],
                        online_s=np.mean(online), err_mean=np.abs(errs).mean(), err_max=np.abs(errs).max(),
                        guaranteed="yes" if errs.min() >= -1e-9 else "no", steps_mean=np.mean(steps), steps_max=max(steps)))
    for o in res:
        gain = t_full - o["online_s"]
        o["break_even_designs"] = (o["offline_s"] / gain) if gain > 0 else float("inf")
    rows_out += res

if rows_out:
    keys = ["expert", "method", "calls", "offline_s", "online_s", "err_mean", "err_max", "guaranteed", "steps_mean",
            "steps_max", "break_even_designs"]
    with open(f"results/cost_study/summary_tol{int(100 * TOL)}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for o in rows_out:
            w.writerow({k: (f"{o[k]:.4g}" if isinstance(o.get(k), float) else o.get(k, "")) for k in keys})
    print(f"target {100 * TOL:.0f} %, per expert: method | full-GBT solves for 40 designs | offline h | online s/design | "
          f"error mean / max % | above full GBT | steps | break-even designs")
    for o in rows_out:
        print(f"{o['expert']:18s} {o['method']:18s} {o['calls']:4d} {o['offline_s'] / 3600:6.2f} h {o['online_s']:6.2f} s "
              f"{o['err_mean']:5.2f}/{o['err_max']:5.2f} % {o['guaranteed']:3s} "
              f"{(str(round(o['steps_mean'], 1)) + '/' + str(o['steps_max'])) if 'steps_mean' in o else '':8s} "
              f"{o['break_even_designs']:.0f}")
