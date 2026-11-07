"""Option A (7 Nov 2026): plain data surrogate for N_cr from GBT labels only (no physics in the loss).

Question: how good is the fast answer when a small model learns ONLY N_cr from the GBT labels, with the
same 38 inputs as the PINN (scaled geometry, flags, lamination parameters, closed-form plate/column loads,
plate wave content)?  If it is clearly better than the PINN head (4.6 % mean on validation, round 20 of
B2 L1000), the fast answer should be split from the mode-shape network (Paper 1 structure: fast tool +
physics check).  If it also stalls around 4 %, the limit is the design space, not our network.

Data: v2 periodic SS mixed lay-ups; train = first N of the 1000 train designs (sample order, seed 11),
validation = last 64 (model choice / early stopping), test = the 100 v2 test designs (seed 99), errors vs
full GBT (N_full).  Training targets = the GBT labels used by the PINN (N_ritz, 400 LOBPCG iterations).

    OMP_NUM_THREADS=1 python3 scripts/a_surrogate.py [--ns 50,100,200,400,1000] [--out results/a_surrogate]
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import sample_design  # noqa: E402
from pinnbuck.net import physics_features  # noqa: E402
from pinnbuck.epinn import wave_features  # noqa: E402

EX = ROOT / "data" / "exact_targets" / "periodic_mixed_SS_v2"


def load(seed, n, sub):
    S = sample_design("omega", n, seed=seed, categorical={"side": "periodic", "bc": "SS"}, layups="mixed", space="v2")
    X, yl, yf, ind, sid = [], [], [], [], []
    for s in S:
        d = torch.load(EX / sub / f"{s.sid}.pt", weights_only=False)
        X.append(np.r_[s.features(), physics_features(s), wave_features(s)])
        yl.append(np.log(float(d["N_ritz"][0]))); yf.append(d["N_full"]); sid.append(s.sid)
        ind.append(_industrial()[s.sid])
    return np.array(X, float), np.array(yl), np.array(yf), np.array(ind), sid


_IND = None
def _industrial():
    """sid -> industrial lay-up flag, from the label summaries (results/b1_targets)."""
    global _IND
    if _IND is None:
        import csv
        _IND = {}
        for sub in ("train", "test"):
            with open(ROOT / "results" / "b1_targets" / f"periodic_mixed_SS_v2_{sub}.csv") as fh:
                _IND.update({r["sid"]: r["industrial"] == "1" for r in csv.DictReader(fh)})
    return _IND


def err(logpred, Nfull):
    return 100 * np.abs(np.exp(logpred) / Nfull - 1)


def stats(e):
    return dict(mean=float(e.mean()), median=float(np.median(e)), p95=float(np.percentile(e, 95)), max=float(e.max()))


# ------------------------------------------------------------------ models
def fit_mlp(Xt, yt, Xv, yv, seed, width=128, depth=3, wd=1e-4, epochs=4000, lr=3e-3):
    torch.manual_seed(seed)
    layers, d = [], Xt.shape[1]
    for _ in range(depth):
        layers += [torch.nn.Linear(d, width), torch.nn.SiLU()]; d = width
    net = torch.nn.Sequential(*layers, torch.nn.Linear(d, 1)).double()
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs, eta_min=lr * 1e-2)
    Xt_, yt_, Xv_, yv_ = map(torch.as_tensor, (Xt, yt, Xv, yv))
    best, state = 1e9, None
    for ep in range(epochs):
        opt.zero_grad(); loss = ((net(Xt_)[:, 0] - yt_) ** 2).mean(); loss.backward(); opt.step(); sch.step()
        if ep % 20 == 0 or ep == epochs - 1:
            with torch.no_grad():
                v = (torch.exp((net(Xv_)[:, 0] - yv_).abs()) - 1).mean().item()
            if v < best:
                best, state = v, {k: x.clone() for k, x in net.state_dict().items()}
    net.load_state_dict(state)
    return lambda X: net(torch.as_tensor(X)).detach().numpy()[:, 0]


def fit_gp(Xt, yt, Xv, yv, seed):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
    k = C(1.0) * RBF(length_scale=np.ones(Xt.shape[1]), length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(1e-6, (1e-10, 1e-2))
    gp = GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=1, random_state=seed).fit(Xt, yt)
    return gp.predict


def fit_gbr(Xt, yt, Xv, yv, seed):
    from sklearn.ensemble import HistGradientBoostingRegressor
    m = HistGradientBoostingRegressor(max_iter=2000, learning_rate=0.03, max_leaf_nodes=15, l2_regularization=1e-3,
                                      random_state=seed, early_stopping=False).fit(Xt, yt)
    return m.predict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="50,100,200,400,1000")
    ap.add_argument("--models", default="mlp5,gp,gbr")
    ap.add_argument("--out", default=str(ROOT / "results" / "a_surrogate"))
    a = ap.parse_args()
    torch.set_num_threads(1)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    Xa, yla, yfa, inda, sida = load(11, 1064, "train")
    Xte, ylte, yfte, indte, sidte = load(99, 100, "test")
    print(f"features {Xa.shape[1]}, train+val {len(Xa)} (industrial {inda.sum()}), test {len(Xte)} "
          f"(industrial {indte.sum()}), loaded in {time.perf_counter() - t0:.0f} s", flush=True)
    Xv, ylv, yfv = Xa[1000:], yla[1000:], yfa[1000:]
    mu, sd = Xa[:1000].mean(0), Xa[:1000].std(0) + 1e-12
    z = lambda X: (X - mu) / sd
    rows = []
    for n in [int(x) for x in a.ns.split(",")]:
        Xt, ylt = Xa[:n], yla[:n]
        for name in a.models.split(","):
            t1 = time.perf_counter()
            if name.startswith("mlp"):
                E = int(name[3:] or 1)
                fs = [fit_mlp(z(Xt), ylt, z(Xv), ylv, seed=s) for s in range(E)]
                f = lambda X, fs=fs: np.mean([g(z(X)) for g in fs], 0)
                f1 = lambda X, fs=fs: fs[0](z(X))
            elif name == "gp":
                g = fit_gp(z(Xt), ylt, z(Xv), ylv, 0); f = lambda X, g=g: g(z(X)); f1 = None
            elif name == "gbr":
                g = fit_gbr(Xt, ylt, Xv, ylv, 0); f = g; f1 = None
            tt = time.perf_counter() - t1
            t2 = time.perf_counter(); pte = f(Xte); tp = (time.perf_counter() - t2) / len(Xte)
            ev, et = err(f(Xv), yfv), err(pte, yfte)
            r = dict(model=name, n=n, fit_s=round(tt, 1), pred_ms=round(1e3 * tp, 3), val=stats(ev), test=stats(et),
                     test_industrial=stats(et[indte]) if indte.any() else None,
                     test_general=stats(et[~indte]) if (~indte).any() else None,
                     test_bias=float(np.mean(100 * (np.exp(pte) / yfte - 1))))
            if f1 is not None:
                r["test_single"] = stats(err(f1(Xte), yfte))
            rows.append(r)
            print(f"N {n:5d} {name:5s} fit {tt:7.1f} s | val mean {r['val']['mean']:.2f} max {r['val']['max']:.2f} | "
                  f"test mean {r['test']['mean']:.2f} p95 {r['test']['p95']:.2f} max {r['test']['max']:.2f} bias "
                  f"{r['test_bias']:+.2f} | ind {r['test_industrial']['mean'] if r['test_industrial'] else float('nan'):.2f} "
                  f"gen {r['test_general']['mean'] if r['test_general'] else float('nan'):.2f}"
                  + (f" | single net test mean {r['test_single']['mean']:.2f} max {r['test_single']['max']:.2f}" if f1 else ""),
                  flush=True)
            if n == max(int(x) for x in a.ns.split(",")):
                np.savetxt(out / f"test_pred_{name}_N{n}.csv",
                           np.c_[np.exp(pte), yfte, err(pte, yfte), indte],
                           delimiter=",", header="N_pred,N_full,abs_err_pct,industrial", comments="",
                           fmt=["%.6f", "%.6f", "%.4f", "%d"])
    (out / "results.json").write_text(json.dumps(rows, indent=1))
    print("A_SURROGATE_DONE", flush=True)


if __name__ == "__main__":
    main()
