"""Expert data tables + plain-data baseline (9 Nov 2026).

(1) Writes every expert's GBT numbers to small CSV files that are kept in git (the big label / operator files are
    not): results/expert_tables/<tag>_train.csv (152 designs in sample order: sid, b_s, beta, full-GBT N, the
    label's 3 Ritz loads, label time) and <tag>_test.csv (40 designs: sid, b_s, beta, full-GBT N).
(2) Baseline for the label study: a Gaussian process (same kernel as scripts/a_surrogate.py, ARD RBF + white
    noise, log N) on the SAME first NL labelled training designs as the PINN (NL = 10, 20, 40, 120), inputs b_s
    and beta, errors on the same 40 test designs vs full GBT.  Data only, no physics; nothing is tuned.
    results/expert_tables/gp_baseline.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, "src")
from pinnbuck import opio
from pinnbuck.designspace import sample_design

OUT = Path("results/expert_tables"); OUT.mkdir(parents=True, exist_ok=True)
LAYS = ["crossply", "quasiiso", "zerorich", "pm45"]; SIDES = ["symmetry", "periodic"]


def gp_fit_predict(Xt, yt, Xs, seed=0):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
    k = C(1.0) * RBF(length_scale=np.ones(Xt.shape[1]), length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(1e-6, (1e-10, 1e-2))
    gp = GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=1, random_state=seed).fit(Xt, yt)
    return gp.predict(Xs)


def scaled(S):
    return np.array([[2 * (s.params["bs"] - 124.25) / (285.5 - 124.25) - 1,
                      2 * (s.params["beta"] - 0.75) / (3.5 - 0.75) - 1] for s in S])


rows = []
for side in SIDES:
    for lay in LAYS:
        tag = f"{side}_{lay}_SS_v2"
        cat = {"side": side, "bc": "SS"}
        Str = sample_design("omega", 152, seed=11, categorical=cat, layups=lay, space="v2")
        Ste = sample_design("omega", 40, seed=99, categorical=cat, layups=lay, space="v2")
        tr = []
        for s in Str:
            r = torch.load(f"data/exact_targets/{tag}/train/{s.sid}.pt", weights_only=False)
            tr.append((s, float(r["N_full"]), [float(x) for x in r["N_ritz"][:3]], float(r.get("time", np.nan))))
        te = []
        for s in Ste:
            op = opio.load(f"data/operators/{tag}/test/{s.sid}.pkl")
            te.append((s, float(op.lam_full))); del op
        with open(OUT / f"{tag}_train.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["sid", "bs", "beta", "N_full", "N_ritz1", "N_ritz2", "N_ritz3", "label_time_s", "role"])
            for i, (s, nf, nr, t) in enumerate(tr):
                w.writerow([s.sid, s.params["bs"], s.params["beta"], f"{nf:.6f}", *[f"{x:.6f}" for x in nr], f"{t:.1f}",
                            "train" if i < 120 else "validation"])
        with open(OUT / f"{tag}_test.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["sid", "bs", "beta", "N_full"])
            for s, nf in te:
                w.writerow([s.sid, s.params["bs"], s.params["beta"], f"{nf:.6f}"])
        Xtr = scaled([t[0] for t in tr[:120]]); ytr = np.log([t[1] for t in tr[:120]])
        Xte = scaled([t[0] for t in te]); Nte = np.array([t[1] for t in te])
        for nl in (10, 20, 40, 120):
            e = np.abs(np.exp(gp_fit_predict(Xtr[:nl], ytr[:nl], Xte)) / Nte - 1) * 100
            rows.append(dict(expert=tag, n_labels=nl, mean=round(e.mean(), 3), max=round(e.max(), 3)))
            print(f"{tag:28s} GP NL {nl:3d}: test mean {e.mean():.2f} %  max {e.max():.2f} %", flush=True)
with open(OUT / "gp_baseline.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
