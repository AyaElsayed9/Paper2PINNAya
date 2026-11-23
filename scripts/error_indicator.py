"""Label-free physics error indicator from the network's modes vs the GP's own uncertainty (23 Nov 2026, step 1).

Story question (Aya, 23 Nov): the GP gives only a number (black box); the network gives mode shapes, so the physics
can be checked.  This script measures how well a label-free physics indicator tells how wrong an answer is.

Per expert (4 lay-ups x 2 side conditions), per test design (the same 40 unseen designs as all expert tables), per
network (physics-only L0 and hybrid L120; crossply / symmetry also L10, L20, L40), from the stored test operators:
  head            the network's instant N_cr (no operators)
  check after s steps (s = 0 = Ritz of the predicted modes, 1, 3, 5, 12): block LOBPCG (R = 8, energy metric as
                  preconditioner, as in all earlier checks).  N_s is an upper bound of the design's Ritz-space N_cr.
  indicators (no labels): eta_s = ||G x - mu K x||_{M^-1} / mu of the best vector (x^T K x = 1, M ~ K the energy
                  metric), i.e. the relative residual; est_s = eta_s^2 / relgap, the quadratic eigenvalue-error
                  estimate (relgap = (mu1 - mu2) / mu1 of the current Ritz values, floored at 1e-3);
                  dlt_s = N_{s-1} / N_s - 1, the change of the last step.
  truth           full GBT N (stored with every test design), only to score the indicators.
GP (same first 120 / NL labelled training designs, inputs b_s and beta, log N, kernel of expert_tables.py):
  mean and predictive std -> its own error indicator (relative std).
Output: results/error_indicator/<tag>.csv (resumable, one row per design and model) and the summary script
scripts/error_indicator_summary.py.

    python scripts/error_indicator.py symmetry_crossply_SS_v2 --budget 520
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
from pinnbuck import epinn, opio                                   # noqa: E402
from pinnbuck.designspace import sample_design                     # noqa: E402
from pinnbuck.spectral import cat, dot, lincomb, rayleigh_ritz     # noqa: E402

STEPS = (0, 1, 3, 5, 12)


def lobpcg_trace(rs, X0, T, R=8, iters=12):
    """The block LOBPCG of spectral.lobpcg (same steps), but it also returns, after every step, the critical load,
    the relative residual and the relative gap of the best Ritz vector."""
    X = X0
    KX, GX = rs.apply_K(X), rs.apply_G(X)
    mu, C = rayleigh_ritz(X, KX, GX, R)
    X, KX, GX = lincomb(X, C), lincomb(KX, C), lincomb(GX, C)
    P = KP = GP = None
    out = []
    for it in range(iters + 1):
        res = [g - mu[:, None, None] * k for g, k in zip(GX, KX)]
        W = T(res)
        eps = float(torch.sqrt(torch.clamp(dot([r[:1] for r in res], [w[:1] for w in W])[0, 0], min=0.0)))
        m0, m1 = float(mu[0]), float(mu[1])
        out.append((1.0 / m0, eps / m0, max((m0 - m1) / m0, 1e-3)))
        if it == iters:
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
    return out


def gp_mean_std(tag, nl, Xs):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{tag}_train.csv")))[:nl]
    Xt = np.array([[float(r["bs"]), float(r["beta"])] for r in tab]); yt = np.log([float(r["N_full"]) for r in tab])
    sc = lambda X: np.c_[2 * (X[:, 0] - 124.25) / (285.5 - 124.25) - 1, 2 * (X[:, 1] - 0.75) / (3.5 - 0.75) - 1]
    k = C(1.0) * RBF(length_scale=np.ones(2), length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(1e-6, (1e-10, 1e-2))
    gp = GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=1, random_state=0).fit(sc(Xt), yt)
    m, s = gp.predict(sc(Xs), return_std=True)
    return np.exp(m), s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--budget", type=float, default=520.0, help="stop after this many seconds (resumable)")
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    t_start = time.time()
    side, lay = a.tag.split("_")[:2]
    S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")[:a.n]
    models = ["physics_L0", "hybrid_L120"] + (["hybrid_L10", "hybrid_L20", "hybrid_L40"] if a.tag == "symmetry_crossply_SS_v2" else [])
    cks = {m: torch.load(ROOT / f"runs/experts/{a.tag}_{m}.pt", weights_only=False) for m in models}
    inp = {m: cks[m]["args"].get("inputs", "all") for m in models}
    assert len(set(map(str, inp.values()))) == 1, inp
    epinn.INPUT_COLS = None if inp[models[0]] == "all" else [0, 1]

    Xs = np.array([[s.params["bs"], s.params["beta"]] for s in S])
    gp = {nl: gp_mean_std(a.tag, nl, Xs) for nl in ([10, 20, 40, 120] if a.tag == "symmetry_crossply_SS_v2" else [120])}

    out = ROOT / "results/error_indicator"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / f"{a.tag}.csv"
    done = {(r["sid"], r["model"]) for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    cols = ["sid", "bs", "beta", "N_full", "model", "N_head"] + \
        [f"{k}{s}" for s in STEPS for k in ("N", "eta", "gap")] + ["t_prep", "t_steps"]
    nets = {}
    with open(fcsv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if not done:
            w.writeheader()
            for nl, (gm, gs) in gp.items():                     # GP rows: mean, predictive std in eta0
                for i, s in enumerate(S):
                    w.writerow(dict(sid=s.sid, bs=s.params["bs"], beta=s.params["beta"], N_full="", model=f"gp_L{nl}",
                                    N_head=f"{gm[i]:.6f}", eta0=f"{gs[i]:.6f}"))
        for s in S:
            todo = [m for m in models if (s.sid, m) not in done]
            if not todo:
                continue
            if time.time() - t_start > a.budget:
                print("budget reached, stopping (resumable)", flush=True); break
            t = time.perf_counter()
            d = epinn.prepare(opio.load(ROOT / f"data/operators/{a.tag}/test/{s.sid}.pkl"))
            t_prep = time.perf_counter() - t
            for m in todo:
                if m not in nets:
                    nets[m] = epinn.net_from_ckpt(cks[m], d)
                t = time.perf_counter()
                with torch.no_grad():
                    Y, logN = nets[m](torch.as_tensor(d.x[None]))
                    Y = Y[0][:, d.sel].double()
                tr = lobpcg_trace(d.rs, d.rs.unflat(d.fm.color(Y)), d.em.solve, R=cks[m]["args"]["R"], iters=max(STEPS))
                r = dict(sid=s.sid, bs=s.params["bs"], beta=s.params["beta"], N_full=f"{d.lam_full:.6f}", model=m,
                         N_head=f"{float(torch.exp(logN[0, 0])):.6f}", t_prep=f"{t_prep:.2f}",
                         t_steps=f"{time.perf_counter() - t:.2f}")
                for st in STEPS:
                    r[f"N{st}"], r[f"eta{st}"], r[f"gap{st}"] = (f"{v:.6g}" for v in tr[st])
                w.writerow(r); fh.flush()
                e = lambda v: 100 * (float(v) / d.lam_full - 1)
                print(f"{s.sid} {m:12s} head {e(r['N_head']):+6.2f} %  Ritz {e(r['N0']):+6.2f} % eta {float(r['eta0']):.3f} | "
                      f"5 steps {e(r['N5']):+.2f} % eta {float(r['eta5']):.4f} | 12 steps {e(r['N12']):+.3f} % eta "
                      f"{float(r['eta12']):.5f}", flush=True)
            del d


if __name__ == "__main__":
    main()
