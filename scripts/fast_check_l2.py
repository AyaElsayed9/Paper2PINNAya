"""Fast physics check of a NEW design vs full GBT for the LEVEL-2 network (one PINN, 38 inputs, all industrial
lay-ups, periodic, SS), timed per design on this machine (10 Nov 2026, plan B).  Same pipeline as
scripts/fast_check.py; the test designs are the 45 level-2 test designs (data/operators/periodic_industrial_SS_v2/test,
only their design definition `info["sample"]` is used: everything else is recomputed online).

Starts (offline: the 518 labelled training designs of the network, their inputs, layouts and labels):
  nn518 : modes of the nearest labelled design in the network's standardised 38 inputs
  knn2  : modes of the 2 nearest labelled designs (16 vectors)
  net   : the level-2 network (hybrid, 518 labels)
Check values are upper bounds of the Ritz-space critical load with the coarse pre-buckling state; errors are vs full
GBT (solver's own pre-buckling state).  results/fast_check/level2_industrial.csv (resumable) + summary.

    python scripts/fast_check_l2.py
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
from pinnbuck import epinn, opio                         # noqa: E402
from pinnbuck.designspace import Sample                  # noqa: E402
from pinnbuck.families import family_context             # noqa: E402
from pinnbuck.operators import build_operator            # noqa: E402
from pinnbuck.spectral import lobpcg                     # noqa: E402

STEPS = (3, 5, 8, 12)
STARTS = ("nn518", "knn2", "net")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="runs/level2/industrial_hybrid_L518.pt")
    ap.add_argument("--files", default="runs/level2/files_industrial.txt")
    ap.add_argument("--n-train", type=int, default=518)
    ap.add_argument("--labels", default="data/exact_targets/periodic_mixed_SS_v2/train")
    ap.add_argument("--test", default="data/operators/periodic_industrial_SS_v2/test")
    ap.add_argument("--n", type=int, default=45)
    ap.add_argument("--coarsen", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    tests = sorted((ROOT / a.test).glob("*.pkl"))[:a.n]
    S = []
    for f in tests:
        sd = opio.load(f).info["sample"]
        s = Sample(family=sd["family"], params=sd["params"], layups=sd["layups"], fixed=sd["fixed"], space=sd["space"])
        assert s.sid == sd["sid"]
        S.append(s)
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / "basis_periodic.pkl", "rb"))

    # ---- offline, once: network + inputs / layout / labels of the 518 labelled designs
    t = time.perf_counter()
    ck = torch.load(ROOT / a.model, weights_only=False)
    epinn.INPUT_COLS = None if ck["args"].get("inputs", "all") == "all" else [0, 1]
    files = [ROOT / l.strip() for l in open(ROOT / a.files) if l.strip()][:a.n_train]
    X, BANK = [], []
    for f in files:
        op = opio.load(f); X.append(epinn.design_inputs(op)); sel = epinn.design_sel(op); sid = op.sid; del op
        BANK.append((sel, torch.load(ROOT / a.labels / f"{sid}.pt", weights_only=False)["T"]))
    X = np.stack(X).astype(float); mu, sg = X.mean(0), X.std(0); sg[sg < 1e-6] = 1.0
    Z = (X - mu) / sg
    print(f"offline (518 labels in the fixed layout + network) {time.perf_counter() - t:.0f} s", flush=True)

    def label_start(d, j):
        sel, T = BANK[j]
        full = torch.zeros(T.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=T.dtype)
        full[:, sel] = T
        return full[:, d.sel]

    out = ROOT / "results/fast_check"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / "level2_industrial.csv"
    done = {r["sid"] for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    cols = ["sid", "bs", "beta", "N_full", "t_full", "t_ops", "t_prep"] + \
        [f"{s}_{k}" for s in STARTS for k in ["t_start"] + [f"t{j}" for j in STEPS] + [f"e{j}" for j in STEPS]]
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
            for name in STARTS:
                t = time.perf_counter()
                if name == "net":
                    with torch.no_grad():
                        Y = net(torch.as_tensor(d.x[None]))[0][0, :, d.sel].double()
                else:
                    z = (np.asarray(d.x, float) - mu) / sg
                    order = np.argsort(((Z - z) ** 2).sum(1))
                    Y = torch.cat([label_start(d, int(j)) for j in order[:1 if name == "nn518" else 2]]).double()
                X0 = d.rs.unflat(d.fm.color(Y))
                r[f"{name}_t_start"] = time.perf_counter() - t
                stamps = [time.perf_counter()]
                _, _, h = lobpcg(d.rs, R=8, iters=max(STEPS), X0=X0, precond=d.em.solve, tol=0, seed=1,
                                 callback=lambda it, mu_: stamps.append(time.perf_counter()))
                for j in STEPS:
                    r[f"{name}_t{j}"] = stamps[j] - stamps[0]
                    r[f"{name}_e{j}"] = 100 * (h[j - 1] / N_full - 1)
            w.writerow(r); fh.flush()
            base = r["t_ops"] + r["t_prep"]
            print(f"{s.sid} bs {s.params['bs']:.0f} beta {s.params['beta']:.2f}: full GBT {r['t_full']:.1f} s | check "
                  f"{r['t_ops']:.1f} + {r['t_prep']:.1f} s + steps | " + " | ".join(
                      f"{n}: 5 st {base + r[f'{n}_t5']:.1f} s {r[f'{n}_e5']:+.2f} %, 12 st {r[f'{n}_e12']:+.2f} %"
                      for n in STARTS), flush=True)
            del d, op

    rows = list(csv.DictReader(open(fcsv)))
    f = lambda k: np.array([float(x[k]) for x in rows])
    base = f("t_ops") + f("t_prep")
    print(f"\nlevel 2: {len(rows)} test designs, mean (max) on this machine")
    print(f"  full GBT                       {f('t_full').mean():6.2f} s ({f('t_full').max():.1f})")
    print(f"  check: operators (coarse pre)  {f('t_ops').mean():6.2f} s   energy space + metric {f('t_prep').mean():.2f} s")
    for name in STARTS:
        for j in STEPS:
            tt = base + f(f"{name}_t_start") + f(f"{name}_t{j}"); e = np.abs(f(f"{name}_e{j}"))
            print(f"  {name:6s} {j:2d} steps: time {tt.mean():5.2f} s ({tt.max():.1f}) = {100 * tt.mean() / f('t_full').mean():3.0f} % "
                  f"of full GBT | |error| {e.mean():.2f} / {e.max():.2f} %, below full GBT {int((f(f'{name}_e{j}') < 0).sum())}")


if __name__ == "__main__":
    main()
