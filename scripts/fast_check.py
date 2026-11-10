"""Fast physics check of a NEW design vs full GBT, timed per design on this machine (10 Nov 2026, plan B).

Online, per test design (nothing of the design is precomputed):
  full_gbt : the reference, full GBT exactly as for the labels (cross-section + pre-buckling + eigen solve, nev 3)
  ops      : cross-section analysis + pre-buckling state on a 4x coarser longitudinal mesh (interpolated to the fine
             Gauss points, `families._Omega.pre`) + operators in the shared basis.  No factorisation of the fine
             full-GBT system.
  prep     : energy space + energy metric (built ONCE; it is both the whitening and the preconditioner)
  steps    : block LOBPCG (R = 8) from a start, times and errors after 3 / 5 / 8 / 12 steps
Starts (offline per expert: the 120 GBT labels in the fixed layout, the trained network):
  nn120    : modes of the nearest labelled design (scaled b_s, beta)
  net      : the expert network (hybrid, 120 labels)
Every check value is an upper bound of the critical load of the design's Ritz space (with the coarse pre-buckling
state); the error columns are vs full GBT with the solver's own pre-buckling state, so they include both
approximations.  results/fast_check/<tag>.csv (resumable) and the summary at the end.

    python scripts/fast_check.py symmetry_crossply_SS_v2 --n 40
"""
import argparse
import csv
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio                         # noqa: E402
from pinnbuck.designspace import sample_design           # noqa: E402
from pinnbuck.families import family_context             # noqa: E402
from pinnbuck.operators import build_operator            # noqa: E402
from pinnbuck.spectral import lobpcg                     # noqa: E402
from check_start_baselines import scaled                 # noqa: E402

STEPS = (3, 5, 8, 12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--coarsen", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    side, lay = a.tag.split("_")[:2]
    S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")[:a.n]
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))

    # ---- offline, once per expert
    t = time.perf_counter()
    ck = torch.load(ROOT / f"runs/experts/{a.tag}_hybrid_L120.pt", weights_only=False)
    epinn.INPUT_COLS = None if ck["args"].get("inputs", "all") == "all" else [0, 1]
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab])
    bank = []
    for r in tab:
        op = opio.load(ROOT / f"data/operators/{a.tag}/train/{r['sid']}.pkl")
        sel = epinn.design_sel(op, epinn.layout_labels(op)); del op
        bank.append((sel, torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{r['sid']}.pt", weights_only=False)["T"]))
    t_off = time.perf_counter() - t
    print(f"{a.tag}: offline (labels in the fixed layout + network) {t_off:.0f} s", flush=True)

    out = ROOT / "results/fast_check"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / f"{a.tag}.csv"
    done = {r["sid"] for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    cols = ["sid", "bs", "beta", "N_full", "t_full", "t_ops", "t_prep"] + \
        [f"{s}_{k}" for s in ("nn120", "net") for k in ["t_start"] + [f"t{j}" for j in STEPS] + [f"e{j}" for j in STEPS]]
    net = None
    with open(fcsv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if not done:
            w.writeheader()
        for s in S:
            if s.sid in done:
                continue
            r = dict(sid=s.sid, bs=s.params["bs"], beta=s.params["beta"])
            t = time.perf_counter()                                        # reference
            with family_context(s) as fc:
                N_full = float(fc.solve(nev=3, return_mode=False)["N"])
            r["t_full"] = time.perf_counter() - t; r["N_full"] = N_full
            t = time.perf_counter()                                        # fast check, online
            op = build_operator(s, sb, lam_full=float("nan"), with_shared_reference=False, pre_coarsen=a.coarsen)
            op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32); op.info["dictionary"] = "d2"
            r["t_ops"] = time.perf_counter() - t
            op.lam_full = N_full                                           # only for the error columns
            t = time.perf_counter()
            d = epinn.prepare(op)
            r["t_prep"] = time.perf_counter() - t
            if net is None:
                net = epinn.net_from_ckpt(ck, d)
            for name in ("nn120", "net"):
                t = time.perf_counter()
                if name == "nn120":
                    j = int(np.argmin(((Xtr - scaled(s.params["bs"], s.params["beta"])) ** 2).sum(1)))
                    sel, T = bank[j]
                    full = torch.zeros(T.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=T.dtype)
                    full[:, sel] = T; Y = full[:, d.sel].double()
                else:
                    with torch.no_grad():
                        Y = net(torch.as_tensor(d.x[None]))[0][0, :, d.sel].double()
                X0 = d.rs.unflat(d.fm.color(Y))
                r[f"{name}_t_start"] = time.perf_counter() - t
                stamps = [time.perf_counter()]
                _, _, h = lobpcg(d.rs, R=8, iters=max(STEPS), X0=X0, precond=d.em.solve, tol=0, seed=1,
                                 callback=lambda it, mu: stamps.append(time.perf_counter()))
                for j in STEPS:
                    r[f"{name}_t{j}"] = stamps[j] - stamps[0]
                    r[f"{name}_e{j}"] = 100 * (h[j - 1] / N_full - 1)
            w.writerow(r); fh.flush()
            base = r["t_ops"] + r["t_prep"]
            print(f"{s.sid} bs {s.params['bs']:.0f} beta {s.params['beta']:.2f}: full GBT {r['t_full']:.1f} s | check "
                  f"{r['t_ops']:.1f} + {r['t_prep']:.1f} s + steps | nn120: 5 steps {base + r['nn120_t5']:.1f} s "
                  f"{r['nn120_e5']:+.2f} %, 12 steps {base + r['nn120_t12']:.1f} s {r['nn120_e12']:+.2f} % | "
                  f"net: 5 steps {r['net_e5']:+.2f} %, 12 steps {r['net_e12']:+.2f} %", flush=True)
            del d, op

    rows = list(csv.DictReader(open(fcsv)))
    f = lambda k: np.array([float(x[k]) for x in rows])
    base = f("t_ops") + f("t_prep")
    print(f"\n{a.tag}: {len(rows)} test designs, mean (max) on this machine")
    print(f"  full GBT                       {f('t_full').mean():6.2f} s ({f('t_full').max():.1f})")
    print(f"  check: operators (coarse pre)  {f('t_ops').mean():6.2f} s   energy space + metric {f('t_prep').mean():.2f} s")
    for name in ("nn120", "net"):
        for j in STEPS:
            tt = base + f(f"{name}_t_start") + f(f"{name}_t{j}"); e = np.abs(f(f"{name}_e{j}"))
            print(f"  {name:6s} {j:2d} steps: time {tt.mean():5.2f} s ({tt.max():.1f}) = {100 * tt.mean() / f('t_full').mean():3.0f} % "
                  f"of full GBT | |error| {e.mean():.2f} / {e.max():.2f} %, below full GBT {int((f(f'{name}_e{j}') < 0).sum())}")


if __name__ == "__main__":
    main()
