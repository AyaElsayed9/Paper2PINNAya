"""What does the network add to the physics check?  (9 Nov 2026, Aya: "do what we need before writing")

The check = block LOBPCG (energy metric, R vectors) in the design's Ritz space; every value is an upper bound of
the Ritz-space critical load.  Here the SAME check is started from different starting modes on the same test
designs, and the error vs full GBT is recorded after each step:

  random      : random start (no network, no labels) -- what the check costs without any learning
  nn<NL>      : the GBT modes of the nearest labelled training design (in scaled b_s, beta) among the first NL
                training designs, copied into the test design's coordinates (no network; uses NL labels)
  net_<name>  : the trained network's predicted modes (as in eval_epinn.py)

Nothing is tuned; the test designs, steps and seeds are those of eval_epinn.py.  Results per design are appended
to results/check_baselines/<tag>.csv (resumable), the summary is printed at the end.

    python scripts/check_start_baselines.py symmetry_crossply_SS_v2 --models physics_L0,hybrid_L120
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio                 # noqa: E402
from pinnbuck.ritzspace import EnergyMetric      # noqa: E402
from pinnbuck.spectral import lobpcg             # noqa: E402


def scaled(bs, beta):
    return np.array([2 * (bs - 124.25) / (285.5 - 124.25) - 1, 2 * (beta - 0.75) / (3.5 - 0.75) - 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--models", default="physics_L0,hybrid_L120")
    ap.add_argument("--nn", default="10,120")
    ap.add_argument("--steps", default="3,5,8,12,20,30")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    steps = [int(s) for s in a.steps.split(",")]; S = max(steps)
    out = ROOT / "results/check_baselines"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / f"{a.tag}.csv"
    done = set()
    if fcsv.exists():
        done = {(r["sid"], r["start"]) for r in csv.DictReader(open(fcsv))}

    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]  # sample order
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab])
    nn_cache = {}

    def nn_start(d, NL, rank=0):
        x = scaled(d.params["bs"], d.params["beta"])
        j = int(np.argsort(((Xtr[:NL] - x) ** 2).sum(1))[rank])
        sid = tab[j]["sid"]
        if sid not in nn_cache:
            op = opio.load(ROOT / f"data/operators/{a.tag}/train/{sid}.pkl")
            sel = epinn.design_sel(op, epinn.layout_labels(op)); del op
            T = torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{sid}.pt", weights_only=False)["T"]
            nn_cache[sid] = (sel, T)
        sel, T = nn_cache[sid]
        full = torch.zeros(T.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=T.dtype)
        full[:, sel] = T
        Y = full[:, d.sel]
        return Y, sid

    nets = {}
    ar0 = torch.load(ROOT / f"runs/experts/{a.tag}_{a.models.split(',')[0]}.pt", weights_only=False)["args"]
    epinn.INPUT_COLS = None if ar0.get("inputs", "all") == "all" else [0, 1]   # before prepare() (sets d.x)
    files = sorted((ROOT / f"data/operators/{a.tag}/test").glob("*.pkl"))
    with open(fcsv, "a", newline="") as fh:
        w = csv.writer(fh)
        if not done:
            w.writerow(["sid", "bs", "beta", "N_full", "start", "start_from", "time_s"] + [f"s{s}" for s in range(1, S + 1)])
        for f in files:
            d = epinn.prepare(opio.load(f))
            em = EnergyMetric(d.rs, groups="labels", ridge=1e-8)
            starts = [("random", None, "-")]
            for NL in [int(v) for v in a.nn.split(",")]:
                Y, sid = nn_start(d, NL); starts.append((f"nn{NL}", Y, sid))
            for m in a.models.split(","):
                if m not in nets:
                    ck = torch.load(ROOT / f"runs/experts/{a.tag}_{m}.pt", weights_only=False); ar = ck["args"]
                    epinn.INPUT_COLS = None if ar.get("inputs", "all") == "all" else [0, 1]
                    nets[m] = epinn.net_from_ckpt(ck, d)
                with torch.no_grad():
                    Y, _ = nets[m](torch.as_tensor(d.x[None]))
                starts.append((f"net_{m}", Y[0, :, d.sel], m))
            # combined starts (9 Nov evening): the first Rayleigh-Ritz step keeps the best R of all vectors
            byname = {n: Y for n, Y, _ in starts}
            starts.append(("knn4_10", torch.cat([nn_start(d, 10, k)[0] for k in range(4)]), "4 nearest of 10"))
            starts.append(("knn2_120", torch.cat([nn_start(d, 120, k)[0] for k in range(2)]), "2 nearest of 120"))
            for m in a.models.split(","):
                starts.append((f"nn10+net_{m}", torch.cat([byname["nn10"], byname[f"net_{m}"].double()]), m))
            for name, Y, src in starts:
                if (d.sid, name) in done:
                    continue
                X0 = None if Y is None else d.rs.unflat(d.fm.color(Y.double()))
                t0 = time.perf_counter()
                _, _, h = lobpcg(d.rs, R=8, iters=S, X0=X0, precond=em.solve, tol=0, seed=1)
                dt = time.perf_counter() - t0
                e = [100 * (v / d.lam_full - 1) for v in h]
                w.writerow([d.sid, d.params["bs"], d.params["beta"], d.lam_full, name, src, f"{dt:.1f}"] + [f"{v:.4f}" for v in e])
                fh.flush()
                print(f"{d.sid} {name:18s} " + "  ".join(f"s{s} {e[s - 1]:+.2f}" for s in steps) + f"  ({dt:.0f} s)", flush=True)
            del d, em

    rows = list(csv.DictReader(open(fcsv)))
    print(f"\n{a.tag}: |error| vs full GBT after s steps, mean / max over the test designs")
    for name in sorted({r["start"] for r in rows}):
        rr = [r for r in rows if r["start"] == name]
        print(f"  {name:18s} ({len(rr)} designs) " + "  ".join(
            f"s{s} {np.mean([abs(float(r[f's{s}'])) for r in rr]):.2f}/{np.max([abs(float(r[f's{s}'])) for r in rr]):.2f}"
            for s in steps))


if __name__ == "__main__":
    main()
