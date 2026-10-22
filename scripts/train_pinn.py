"""Train the parametric GBT deep-Ritz network (fast PINN) for omega panels.

    python scripts/train_pinn.py --side symmetry --steps 4000 --batch 8 --threads 8

Data: data/operators/<side>/train/*.pkl (physics only, no labels) and
      data/operators/<side>/test/*.pkl  (with the full GBT N_cr for the error).
Output: runs/<name>/model.pt, results.json, and a printed summary:
  lam_R    Rayleigh value of the predicted mode on the panel's GBT operators (upper bound)
  lam_head direct prediction of the network (no GBT at all: the fast path)
both compared with the full GBT (the validated model).
"""
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.pinn import GBTDeepRitz, GBTWaveRitz, TorchPanel, evaluate, mode_scale, train   # noqa: E402


def files(d):
    return sorted(Path(d).glob("*.pkl"))


def load(f):
    return pickle.load(open(f, "rb"))


def summary(res, key):
    e = np.array([100 * (r[key] / r["lam_full"] - 1) for r in res if np.isfinite(r["lam_full"])])
    return f"mean |err| {np.abs(e).mean():5.2f} %   max |err| {np.abs(e).max():5.2f} %   bias {e.mean():+5.2f} %" if e.size else "-"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--J", type=int, default=24, help="Fourier order along z")
    ap.add_argument("--R", type=int, default=6,
                    help="wave model: number of wavenumbers (2R candidates); fourier model: candidates")
    ap.add_argument("--width", type=int, default=256)
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--steps", type=int, default=4000)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--model", default="wave", choices=["wave", "fourier"])
    ap.add_argument("--teacher-steps", type=int, default=1000,
                    help="first steps on the full-GBT teacher modes only (0 = physics from the start)")
    ap.add_argument("--w-teach-late", type=float, default=0.1, help="teacher weight in the physics stage")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--name", default=None)
    ap.add_argument("--float64", action="store_true", help="train in double precision (default float32)")
    a = ap.parse_args()
    torch.set_num_threads(a.threads); torch.manual_seed(a.seed)
    base = ROOT / "data" / "operators" / a.side
    f_tr, f_te = files(base / "train"), files(base / "test")
    ms = mode_scale([load(f) for f in f_tr[:64]])           # Jacobi scaling of the shared modes
    dt = torch.float64 if a.float64 else torch.float32
    J = a.J if a.model == "fourier" else 0
    P_tr = [TorchPanel(load(f), J, mode_scale=ms, dtype=dt) for f in f_tr]   # operators freed on the fly
    P_te = [TorchPanel(load(f), J, mode_scale=ms, dtype=dt) for f in f_te]
    nm, n_in = P_tr[0].nm, P_tr[0].x.numel()
    print(f"{len(P_tr)} training panels, {len(P_te)} test panels, {nm} shared modes", flush=True)
    model = (GBTWaveRitz(n_in, nm, R=a.R, width=a.width, depth=a.depth) if a.model == "wave"
             else GBTDeepRitz(n_in, nm, a.J, R=a.R, width=a.width, depth=a.depth)).to(dt)
    t0 = time.perf_counter()

    def cb(it, tot, rec):
        err = [100 * (r[1] / P_tr[r[0]].lam_full - 1) for r in rec if np.isfinite(P_tr[r[0]].lam_full)]
        c2 = [float(r[3][0]) for r in rec if r[3] is not None]
        print(f"step {it:5d}  loss {tot:+.4f}  batch: Ritz error vs full GBT {np.round(err, 1)} %"
              + (f"  teacher cos2 {np.mean(c2):.4f}" if c2 else "") + f"  {time.perf_counter() - t0:.0f} s", flush=True)
        if it % (5 * 200) == 0 or it == a.steps:
            print("   test:", "Rayleigh", summary(evaluate(model, P_te), "lam_R"), flush=True)

    train(model, P_tr, steps=a.steps, batch=a.batch, lr=a.lr, log_every=200, seed=a.seed, callback=cb,
          teacher_steps=a.teacher_steps, w_teach_late=a.w_teach_late)
    t_train = time.perf_counter() - t0
    del P_tr, P_te
    model = model.double()                                  # final evaluation in double precision
    ev64 = lambda fl: [evaluate(model, [TorchPanel(load(f), J, mode_scale=ms)])[0] for f in fl]
    res_tr, res_te = ev64(f_tr), ev64(f_te)
    P_te = [TorchPanel(load(f), J, mode_scale=ms) for f in f_te]
    # speed of the fast path: network only
    x = torch.stack([p.x for p in P_te])
    with torch.no_grad():
        t = time.perf_counter()
        for _ in range(100):
            model.hL(model.body(x))
        t_pred = (time.perf_counter() - t) / 100 / len(P_te)
    name = a.name or f"{a.side}_J{a.J}_R{a.R}_s{a.steps}"
    out = ROOT / "runs" / name; out.mkdir(parents=True, exist_ok=True)
    torch.save(dict(state=model.state_dict(), args=vars(a), mode_scale=ms, nm=nm, n_in=n_in),
               out / "model.pt")
    json.dump(dict(args=vars(a), train=res_tr, test=res_te, t_train=t_train, t_pred_per_design=t_pred),
              open(out / "results.json", "w"), indent=1)
    print(f"\ntraining time {t_train / 60:.1f} min; prediction {1e3 * t_pred:.3f} ms per design (network only)")
    print("test panels vs full GBT:")
    print("  Rayleigh value of predicted mode:", summary(res_te, "lam_R"))
    print("  direct N_cr prediction (fast)   :", summary(res_te, "lam_head"))
    print(f"results in {out}")


if __name__ == "__main__":
    main()
