"""Cost study for a 5 % accuracy target (Aya, 23 Nov 2026; step 2 of the Paper-2 plan, "GBT-call study").

For each expert, the same 40 unseen test designs (stored operators, fine pre-buckling state = the certified check):
the physics check (block LOBPCG, R = 8, energy metric, as in all checks) is run from five starts and stopped by the
label-free estimate eta (relative residual, scripts/error_indicator.py: eta >= true error in 100 % of 3800 cases):
  random    no labels, no network
  nn10      modes of the nearest of the first 10 labelled training designs (10 full GBT solves offline)
  nn120     the same with 120 labelled designs
  net_L0    physics-only network (0 labels)
  net_L120  hybrid network (120 labels)
The trace is kept until eta <= 2 % or 40 steps, so tolerances of 5 % (the target), 3 % and 2 % can be read from it.
Per start and design: N and eta after every step, the step times.  Truth (full GBT) only to score.
results/cost_study/<tag>.csv (resumable); summary in scripts/cost_study_summary.py.

    python scripts/cost_study.py symmetry_crossply_SS_v2 --budget 540
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio                                   # noqa: E402
from pinnbuck.designspace import sample_design                     # noqa: E402
from pinnbuck.spectral import cat, dot, lincomb, rayleigh_ritz     # noqa: E402
from check_start_baselines import scaled                           # noqa: E402

TOL_STOP, MAX_STEPS = 0.02, 40


def trace(rs, X0, T, R=8):
    """LOBPCG as spectral.lobpcg; after every step (0 = start): N of the best vector, its relative residual, time."""
    t0 = time.perf_counter()
    X = X0
    KX, GX = rs.apply_K(X), rs.apply_G(X)
    mu, C = rayleigh_ritz(X, KX, GX, R)
    X, KX, GX = lincomb(X, C), lincomb(KX, C), lincomb(GX, C)
    P = KP = GP = None
    Ns, etas, ts = [], [], []
    for it in range(MAX_STEPS + 1):
        res = [g - mu[:, None, None] * k for g, k in zip(GX, KX)]
        W = T(res)
        eps = float(torch.sqrt(torch.clamp(dot([r[:1] for r in res], [w[:1] for w in W])[0, 0], min=0.0)))
        Ns.append(1.0 / float(mu[0])); etas.append(eps / float(mu[0])); ts.append(time.perf_counter() - t0)
        if etas[-1] <= TOL_STOP or it == MAX_STEPS:
            break
        KW, GW = rs.apply_K(W), rs.apply_G(W)
        if P is None:
            S, KS, GS = cat(X, W), cat(KX, KW), cat(GX, GW)
        else:
            S, KS, GS = cat(X, W, P), cat(KX, KW, KP), cat(GX, GW, GP)
        mu, C = rayleigh_ritz(S, KS, GS, R)
        Xn, KXn, GXn = lincomb(S, C), lincomb(KS, C), lincomb(GS, C)
        Cp = C.clone(); Cp[:R] = 0.0
        P, KP, GP = lincomb(S, Cp), lincomb(KS, Cp), lincomb(GS, Cp)
        X, KX, GX = Xn, KXn, GXn
    return Ns, etas, ts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--budget", type=float, default=540.0)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    t_start = time.time()
    side, lay = a.tag.split("_")[:2]
    S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")
    cks = {m: torch.load(ROOT / f"runs/experts/{a.tag}_{m}.pt", weights_only=False) for m in ("physics_L0", "hybrid_L120")}
    epinn.INPUT_COLS = None if cks["physics_L0"]["args"].get("inputs", "all") == "all" else [0, 1]
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab])
    cache = {}

    def nn_start(d, NL):
        j = int(np.argmin(((Xtr[:NL] - scaled(d.params["bs"], d.params["beta"])) ** 2).sum(1)))
        sid = tab[j]["sid"]
        if sid not in cache:
            op = opio.load(ROOT / f"data/operators/{a.tag}/train/{sid}.pkl")
            sel = epinn.design_sel(op, epinn.layout_labels(op)); del op
            cache[sid] = (sel, torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{sid}.pt", weights_only=False)["T"])
        sel, T = cache[sid]
        full = torch.zeros(T.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=T.dtype)
        full[:, sel] = T
        return full[:, d.sel].double()

    out = ROOT / "results/cost_study"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / f"{a.tag}.csv"
    done = {(r["sid"], r["start"]) for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    cols = ["sid", "bs", "beta", "N_full", "start", "t_start", "N", "eta", "t"]
    nets = {}
    starts = ["random", "nn10", "nn120", "net_L0", "net_L120"]
    with open(fcsv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if not done:
            w.writeheader()
        for s in S:
            todo = [x for x in starts if (s.sid, x) not in done]
            if not todo:
                continue
            if time.time() - t_start > a.budget:
                print("budget reached (resumable)", flush=True); break
            d = epinn.prepare(opio.load(ROOT / f"data/operators/{a.tag}/test/{s.sid}.pkl"))
            for name in todo:
                t = time.perf_counter()
                if name == "random":
                    X0 = d.em.solve(d.rs.random(8, torch.Generator().manual_seed(1)))
                else:
                    if name.startswith("nn"):
                        Y = nn_start(d, int(name[2:]))
                    else:
                        m = {"net_L0": "physics_L0", "net_L120": "hybrid_L120"}[name]
                        if m not in nets:
                            nets[m] = epinn.net_from_ckpt(cks[m], d)
                        with torch.no_grad():
                            Y = nets[m](torch.as_tensor(d.x[None]))[0][0][:, d.sel].double()
                    X0 = d.rs.unflat(d.fm.color(Y))
                t_st = time.perf_counter() - t
                Ns, etas, ts = trace(d.rs, X0, d.em.solve)
                w.writerow(dict(sid=s.sid, bs=s.params["bs"], beta=s.params["beta"], N_full=f"{d.lam_full:.6f}", start=name,
                                t_start=f"{t_st:.3f}", N=json.dumps([round(v, 6) for v in Ns]),
                                eta=json.dumps([round(v, 6) for v in etas]), t=json.dumps([round(v, 3) for v in ts])))
                fh.flush()
                k5 = next((i for i, v in enumerate(etas) if v <= 0.05), None)
                msg = (f"{k5} steps, error {100 * (Ns[k5] / d.lam_full - 1):+.2f} %" if k5 is not None else "not reached")
                print(f"{s.sid} {name:8s} eta <= 5 %: {msg}", flush=True)
            del d


if __name__ == "__main__":
    main()
