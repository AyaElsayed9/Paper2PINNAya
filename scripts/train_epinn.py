"""Train the physics-only energy-whitened Ritz PINN (`pinnbuck/epinn.py`) on omega panels.

    python scripts/train_epinn.py --train data/operators/symmetry_industrial/train \
        --test data/operators/symmetry_industrial/test --bc CC --method ritzng --outer 30

No labels are used in training: the loss is the GBT energy (Rayleigh-Ritz) of the
predicted modes.  The full-GBT N_cr stored with the operators only reports the error.
One network per end condition (the coordinates of the modes differ between SS and CC).
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
from pinnbuck import epinn   # noqa: E402


def files_of(folder, bc, nmax):
    out = []
    for f in sorted(Path(folder).glob("*.pkl")):
        op = pickle.load(open(f, "rb"))
        if op.params["bc"] == bc:
            out.append(f)
        if len(out) >= nmax:
            break
    return out


def report(tag, res):
    e = np.array([r["N_ritz"] / r["N_full"] - 1 for r in res]); h = np.array([r["N_head"] / r["N_full"] - 1 for r in res])
    return (f"{tag}: Ritz N (upper bound) vs full GBT: mean {100 * np.mean(np.abs(e)):.2f} %, max {100 * np.max(np.abs(e)):.2f} %"
            f"  |  head N: mean {100 * np.mean(np.abs(h)):.2f} %, max {100 * np.max(np.abs(h)):.2f} %")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--bc", default="CC", choices=["SS", "CC"])
    ap.add_argument("--n-train", type=int, default=64)
    ap.add_argument("--method", default="ritzng", choices=["ritzng", "adam", "ritzng+adam"])
    ap.add_argument("--arch", default="decoder", choices=["decoder", "linear"])
    ap.add_argument("--cache", type=int, default=40, help="energy spaces kept in memory")
    ap.add_argument("--n-eval-train", type=int, default=16)
    ap.add_argument("--outer", type=int, default=30)
    ap.add_argument("--inner", type=int, default=20)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--R", type=int, default=8)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--width", type=int, default=256)
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    torch.set_num_threads(a.threads); torch.manual_seed(a.seed)
    t0 = time.perf_counter()
    tr = epinn.DesignPool(files_of(a.train, a.bc, a.n_train), cache=a.cache)
    te = epinn.DesignPool(files_of(a.test, a.bc, 10000), cache=0)
    n = tr.get(0).rs.n
    print(f"{a.bc}: {len(tr)} training designs, {len(te)} test designs, n = {n} per mode; "
          f"inputs ready in {time.perf_counter() - t0:.0f} s", flush=True)
    net = epinn.build_net(a.arch, tr.X.numpy(), tr.get(0).rs, R=a.R, k=a.k, width=a.width, depth=a.depth)
    ev_idx = list(range(min(a.n_eval_train, len(tr))))
    print(f"network: {sum(p.numel() for p in net.parameters()) / 1e6:.1f} M parameters", flush=True)
    t1 = time.perf_counter()

    def ev(it):
        print("   " + report(f"[{it}] train", epinn.evaluate(net, tr, idx=ev_idx)), flush=True)
        print("   " + report(f"[{it}] test ", epinn.evaluate(net, te)) + f"   ({time.perf_counter() - t1:.0f} s)", flush=True)

    lg = lambda s: print(s + f"   ({time.perf_counter() - t1:.0f} s)", flush=True)
    if a.method == "ritzng+adam":
        opt = epinn.train_ritz_ng(net, tr, outer=a.outer, inner=a.inner, batch=a.batch, lr=a.lr, k=a.k, seed=a.seed,
                                  log=lg, eval_fn=ev, eval_every=5)
        epinn.train_adam(net, tr, steps=a.steps, batch=a.batch, lr=a.lr, k=a.k, seed=a.seed, log=lg,
                         eval_fn=ev, eval_every=200, opt=opt)
    elif a.method == "ritzng":
        epinn.train_ritz_ng(net, tr, outer=a.outer, inner=a.inner, batch=a.batch, lr=a.lr, k=a.k, seed=a.seed,
                            log=lambda s: print(s + f"   ({time.perf_counter() - t1:.0f} s)", flush=True),
                            eval_fn=ev, eval_every=5)
    else:
        epinn.train_adam(net, tr, steps=a.steps, batch=a.batch, lr=a.lr, k=a.k, seed=a.seed,
                         log=lambda s: print(s + f"   ({time.perf_counter() - t1:.0f} s)", flush=True),
                         eval_fn=ev, eval_every=200)
    res = dict(train=epinn.evaluate(net, tr, idx=ev_idx), test=epinn.evaluate(net, te), args=vars(a))
    print(report("final train", res["train"])); print(report("final test ", res["test"]))
    for r in res["test"]:
        print(f"   test {r['sid']} bs {r['params']['bs']:.0f} beta {r['params']['beta']:.2f}: full {r['N_full']:.4f}  "
              f"Ritz {r['N_ritz']:.4f} ({100 * (r['N_ritz'] / r['N_full'] - 1):+.2f} %)  head {r['N_head']:.4f} "
              f"({100 * (r['N_head'] / r['N_full'] - 1):+.2f} %)")
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        torch.save(dict(state=net.state_dict(), args=vars(a)), a.out)
        json.dump(res, open(str(a.out).replace(".pt", ".json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
