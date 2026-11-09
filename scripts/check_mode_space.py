"""How good can predicted mode shapes be?  (9 Nov 2026, after the start baselines)

For each test design, the Ritz value (no LOBPCG steps) of:
  nn<NL>   : the 8 GBT modes of the nearest labelled training design (as in check_start_baselines.py)
  knn<K>   : the modes of the K nearest labelled designs together (8K vectors)
  pod<r>   : the r leading POD vectors of all 120 x 8 training label modes (fixed layout, whitened coordinates)
  net_<m>  : the trained network's 8 modes
The pod<r> values say how many shared vectors a network would need to predict (r coefficients per mode instead
of ~11865), i.e. whether a reduced-basis output layer is a realistic way to better mode shapes.
results/check_baselines/<tag>_modespace.csv
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio                 # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from check_start_baselines import scaled          # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--models", default="physics_L0,hybrid_L120")
    ap.add_argument("--pod", default="8,16,32,64,128,256")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab])
    ar0 = torch.load(ROOT / f"runs/experts/{a.tag}_{a.models.split(',')[0]}.pt", weights_only=False)["args"]
    epinn.INPUT_COLS = None if ar0.get("inputs", "all") == "all" else [0, 1]

    # all training label modes in the fixed layout
    rows, n_lay = [], 0
    for r in tab:
        op = opio.load(ROOT / f"data/operators/{a.tag}/train/{r['sid']}.pkl")
        sel = epinn.design_sel(op, epinn.layout_labels(op)); del op
        T = torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{r['sid']}.pt", weights_only=False)["T"]
        rows.append((sel, T)); n_lay = max(n_lay, int(sel.max()) + 1)
    L = torch.zeros(120, 8, n_lay, dtype=torch.float64)
    for i, (sel, T) in enumerate(rows):
        L[i][:, sel] = T
    U, S, _ = torch.linalg.svd(L.reshape(-1, n_lay).T, full_matrices=False)       # (n_lay, 960)
    print(f"POD of {L.shape[0] * 8} label modes: energy kept by r vectors " + ", ".join(
        f"{r}: {float((S[:r] ** 2).sum() / (S ** 2).sum()):.4f}" for r in [int(v) for v in a.pod.split(",")]), flush=True)

    nets = {}
    out = ROOT / f"results/check_baselines/{a.tag}_modespace.csv"
    res = []
    for f in sorted((ROOT / f"data/operators/{a.tag}/test").glob("*.pkl")):
        d = epinn.prepare(opio.load(f))
        x = scaled(d.params["bs"], d.params["beta"])
        order = np.argsort(((Xtr - x) ** 2).sum(1))
        r = dict(sid=d.sid)

        def err(Y):
            mu, _ = epinn.ritz(d, Y)
            return 100 * (1.0 / float(mu[0]) / d.lam_full - 1)
        for NL in (10, 120):
            j = [i for i in order if i < NL][0]
            r[f"nn{NL}"] = err(L[j][:, d.sel])
        for K in (2, 4):
            r[f"knn{K}"] = err(L[order[:K]].reshape(-1, n_lay)[:, d.sel])
        for p in [int(v) for v in a.pod.split(",")]:
            r[f"pod{p}"] = err(U[:, :p].T[:, d.sel].contiguous())
        for m in a.models.split(","):
            if m not in nets:
                ck = torch.load(ROOT / f"runs/experts/{a.tag}_{m}.pt", weights_only=False); ar = ck["args"]
                net = epinn.build_net(ar.get("arch", "linear"), d.x[None].repeat(2, 0), d.rs, R=ar["R"], k=ar["k"],
                                      width=ar["width"], depth=ar["depth"])
                net.load_state_dict(ck["state"]); net.eval(); nets[m] = net
            with torch.no_grad():
                Y, _ = nets[m](torch.as_tensor(d.x[None]))
            r[f"net_{m}"] = err(Y[0, :, d.sel])
        res.append(r)
        print(d.sid + " " + "  ".join(f"{k} {v:+.2f}" for k, v in r.items() if k != "sid"), flush=True)
        del d
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    print(f"\n{a.tag}: Ritz value of the start (no steps), |error| vs full GBT, mean / max")
    for k in res[0]:
        if k != "sid":
            v = np.abs([r[k] for r in res]); print(f"  {k:18s} {v.mean():7.2f} / {v.max():7.2f} %")


if __name__ == "__main__":
    main()
