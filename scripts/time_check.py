"""Clean timing of every step, on this machine, with nothing else running (9 Nov 2026).

Per test design (crossply / symmetry expert, first --n test designs, same samples as the operator files):
  full_gbt   : full GBT as used for the labels / reference (cross-section analysis + pre-buckling + eigen solve, nev 3)
  ops        : what the physics check needs for a NEW design: cross-section analysis + pre-buckling state + the
               operators in the shared basis (no eigen solve)
  prep       : energy space of the design (epinn.prepare) + energy metric
  step       : one check step (block LOBPCG, R = 8), mean of 12 steps started from the network
  net        : one network forward pass (modes + head)
  gp         : one GP prediction (GP fitted once on the 120 labels; fit time reported separately)
Machine: reported from /proc (cores, memory).  results/timing/check_timing.csv and .txt
"""
import argparse
import csv
import os
import pickle
import platform
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio                         # noqa: E402
from pinnbuck.designspace import sample_design           # noqa: E402
from pinnbuck.families import family_context             # noqa: E402
from pinnbuck.operators import build_operator            # noqa: E402
from pinnbuck.ritzspace import EnergyMetric              # noqa: E402
from pinnbuck.spectral import lobpcg                     # noqa: E402


def machine():
    cpu = [l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")]
    mem = [l for l in open("/proc/meminfo") if l.startswith("MemTotal")][0].split()[1]
    return f"{len(cpu)} CPU cores ({cpu[0] if cpu else '?'}), {int(mem) / 2**20:.1f} GB RAM, {platform.system()}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="symmetry_crossply_SS_v2")
    ap.add_argument("--model", default="runs/experts/symmetry_crossply_SS_v2_hybrid_L120.pt")
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    side, lay = a.tag.split("_")[:2]
    S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")[:a.n]
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
    ck = torch.load(ROOT / a.model, weights_only=False)
    epinn.INPUT_COLS = None if ck["args"].get("inputs", "all") == "all" else [0, 1]

    # GP on the 120 labels (as in scripts/expert_tables.py)
    sys.path.insert(0, str(ROOT / "scripts"))
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    sc = lambda bs, be: [2 * (bs - 124.25) / (285.5 - 124.25) - 1, 2 * (be - 0.75) / (3.5 - 0.75) - 1]
    Xtr = np.array([sc(float(r["bs"]), float(r["beta"])) for r in tab]); ytr = np.log([float(r["N_full"]) for r in tab])
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
    t = time.perf_counter()
    gp = GaussianProcessRegressor(C(1.0) * RBF(np.ones(2), (1e-2, 1e3)) + WhiteKernel(1e-6, (1e-10, 1e-2)),
                                  normalize_y=True, n_restarts_optimizer=1, random_state=0).fit(Xtr, ytr)
    t_gpfit = time.perf_counter() - t

    rows, net = [], None
    for s in S:
        r = dict(sid=s.sid)
        t = time.perf_counter()                                            # full GBT
        with family_context(s) as fc:
            rf = fc.solve(nev=3, return_mode=False)
        r["full_gbt"] = time.perf_counter() - t; N_full = float(rf["N"])
        t = time.perf_counter()                                            # operators for the check (no eigen solve)
        op = build_operator(s, sb, lam_full=float("nan"), with_shared_reference=False)
        op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32)
        op.info["dictionary"] = "d2"                                      # as build_operators.py for space v2
        r["ops"] = time.perf_counter() - t
        op.lam_full = N_full
        t = time.perf_counter()
        d = epinn.prepare(op); em = EnergyMetric(d.rs, groups="labels", ridge=1e-8)
        r["prep"] = time.perf_counter() - t
        if net is None:
            net = epinn.net_from_ckpt(ck, d)
        t = time.perf_counter()
        with torch.no_grad():
            Y, logN = net(torch.as_tensor(d.x[None]))
        r["net"] = time.perf_counter() - t
        t = time.perf_counter()
        gp.predict(np.array([sc(s.params["bs"], s.params["beta"])]))
        r["gp"] = time.perf_counter() - t
        X0 = d.rs.unflat(d.fm.color(Y[0, :, d.sel].double()))
        stamps = [time.perf_counter()]
        lam, _, h = lobpcg(d.rs, R=8, iters=12, X0=X0, precond=em.solve, tol=0, seed=1,
                           callback=lambda it, mu: stamps.append(time.perf_counter()))
        dts = np.diff(stamps)
        r["step"] = float(np.mean(dts[1:])) if len(dts) > 1 else float(dts[0])
        r["first_step"] = float(dts[0])
        r["err5"] = 100 * (h[4] / N_full - 1); r["err12"] = 100 * (h[11] / N_full - 1)
        r["N_full"] = N_full
        rows.append(r)
        print(f"{s.sid}: full GBT {r['full_gbt']:.1f} s | check: operators {r['ops']:.1f} s + energy space {r['prep']:.2f} s"
              f" + step {r['step']:.3f} s (first {r['first_step']:.3f}) | net {1e3 * r['net']:.1f} ms | GP {1e3 * r['gp']:.2f} ms",
              flush=True)
        del d, em, op
    out = ROOT / "results/timing"; out.mkdir(parents=True, exist_ok=True)
    with open(out / "check_timing.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    m = {k: float(np.mean([r[k] for r in rows])) for k in rows[0] if k != "sid"}
    lines = [f"Machine: {machine()}; torch threads {a.threads}; nothing else running.",
             f"{len(rows)} test designs of {a.tag}; mean times:",
             f"  full GBT (N_cr, reference)           {m['full_gbt']:.2f} s",
             f"  check: operators of the new design   {m['ops']:.2f} s",
             f"  check: energy space + metric         {m['prep']:.2f} s",
             f"  check: one step                      {m['step']:.3f} s (first step {m['first_step']:.3f} s)",
             f"  check, 5 steps in total              {m['ops'] + m['prep'] + m['first_step'] + 4 * m['step']:.2f} s",
             f"  check, 12 steps in total             {m['ops'] + m['prep'] + m['first_step'] + 11 * m['step']:.2f} s",
             f"  network forward pass                 {1e3 * m['net']:.1f} ms",
             f"  GP prediction                        {1e3 * m['gp']:.2f} ms (fit on 120 labels once: {t_gpfit:.2f} s)"]
    open(out / "check_timing.txt", "w").write("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
