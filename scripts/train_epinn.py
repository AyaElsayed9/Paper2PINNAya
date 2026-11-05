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
from pinnbuck import epinn
from pinnbuck import opio                        # noqa: E402   # noqa: E402


def files_of(folder, bc, nmax):
    out = []
    for f in sorted(Path(folder).glob("*.pkl")):
        op = opio.load(f)
        if op.params["bc"] == bc:
            out.append(f)
        if len(out) >= nmax:
            break
    return out


def report(tag, res):
    e = np.array([r["N_ritz"] / r["N_full"] - 1 for r in res]); h = np.array([r["N_head"] / r["N_full"] - 1 for r in res])
    return (f"{tag}: Ritz N (upper bound) vs full GBT: mean {100 * np.mean(np.abs(e)):.2f} %, max {100 * np.max(np.abs(e)):.2f} %"
            f"  |  head N: mean {100 * np.mean(np.abs(h)):.2f} %, max {100 * np.max(np.abs(h)):.2f} %, "
            f"bias {100 * np.mean(h):+.2f} %, above {int(np.sum(h > 0))}/{len(h)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--bc", default="CC", choices=["SS", "CC"])
    ap.add_argument("--n-train", type=int, default=64)
    ap.add_argument("--method", default="ritzng", choices=["ritzng", "adam", "ritzng+adam", "ngfull", "ngfull+adam"])
    ap.add_argument("--epochs", type=int, default=10, help="fit passes per outer iteration (ngfull)")
    ap.add_argument("--arch", default="decoder", choices=["decoder", "linear", "hybrid"])
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
    ap.add_argument("--init", default=None, help="start from these weights (transfer, as in Paper 1)")
    ap.add_argument("--oracle", default=None, help="folder of training designs with an oracle N_cr (full GBT)")
    ap.add_argument("--gamma", type=float, default=0.0, help="weight of the Paper 1 data term SmoothL1(N/N_oracle, 1)")
    ap.add_argument("--n-oracle-batch", type=int, default=4)
    ap.add_argument("--mode", default="physics", choices=["physics", "data", "hybrid"],
                    help="three-way test: physics-only, data-only (exact modes and loads), hybrid (ngfull only)")
    ap.add_argument("--exact", default=None, help="folder of exact targets (scripts/exact_targets.py)")
    ap.add_argument("--val", type=int, default=0,
                    help="hold out the last VAL training designs (with exact targets) as validation; "
                         "the saved weights are chosen on them, never on the test designs")
    ap.add_argument("--lr-final", type=float, default=None, help="cosine decay of the learning rate to this value")
    ap.add_argument("--ema", type=float, default=None, help="moving average of the weights (decay per step)")
    ap.add_argument("--no-test-log", action="store_true", help="do not look at the test designs during training")
    a = ap.parse_args()
    torch.set_num_threads(a.threads); torch.manual_seed(a.seed)
    t0 = time.perf_counter()
    all_files = files_of(a.train, a.bc, a.n_train)
    if a.val:
        assert a.exact, "validation needs exact targets"
        all_files, val_files = all_files[:-a.val], all_files[-a.val:]
    tr = epinn.DesignPool(all_files, cache=a.cache)
    te = epinn.DesignPool(files_of(a.test, a.bc, 10000), cache=0)
    va = None
    if a.val:
        va = epinn.DesignPool(val_files, cache=0)
        for i, f in enumerate(va.files):
            r = torch.load(Path(a.exact) / f"{Path(f).stem}.pt", weights_only=False)
            va.lam_full[i] = float(r["N_ritz"][0])
        print(f"{len(va)} validation designs held out", flush=True)
    n = tr.get(0).rs.n
    print(f"{a.bc}: {len(tr)} training designs, {len(te)} test designs, n = {n} per mode; "
          f"inputs ready in {time.perf_counter() - t0:.0f} s", flush=True)
    net = epinn.build_net(a.arch, tr.X.numpy(), tr.get(0).rs, R=a.R, k=a.k, width=a.width, depth=a.depth)
    ev_idx = list(range(min(a.n_eval_train, len(tr))))
    if a.init:
        net.load_state_dict(torch.load(a.init, weights_only=False)["state"])
        print(f"weights from {a.init}", flush=True)
    oracle = {}
    if a.oracle:
        names = {f.name for f in Path(a.oracle).glob("*.pkl")}
        for i, f in enumerate(tr.files):
            if Path(f).name in names and np.isfinite(tr.lam_full[i]):
                oracle[i] = tr.lam_full[i]
        ev_idx = sorted(oracle)[:a.n_eval_train] or ev_idx
        print(f"{len(oracle)} oracle designs (full GBT), gamma = {a.gamma}", flush=True)
    exact = None
    if a.exact:
        exact = []
        for i, f in enumerate(tr.files):
            r = torch.load(Path(a.exact) / f"{Path(f).stem}.pt", weights_only=False)
            exact.append((r["T"].double(), torch.log(torch.as_tensor(r["N_ritz"][:a.k])).float()))
            if not np.isfinite(tr.lam_full[i]):
                tr.lam_full[i] = float(r["N_ritz"][0])        # exact load in the Ritz space (reporting)
        print(f"exact targets for {len(exact)} training designs, mode = {a.mode}", flush=True)
    print(f"network: {sum(p.numel() for p in net.parameters()) / 1e6:.1f} M parameters", flush=True)
    t1 = time.perf_counter()

    best = dict(score=np.inf)
    best_modes = dict(score=np.inf)      # weights chosen for the mode shapes (validation Ritz bound)

    def swap(sd):
        old = {k: v.detach().clone() for k, v in net.state_dict().items()}
        net.load_state_dict(sd, strict=False)
        return old

    def score(res):
        return float(np.mean([abs(r["N_head"] / r["N_full"] - 1) for r in res]))

    def ev(it):
        print("   " + report(f"[{it}] train", epinn.evaluate(net, tr, idx=ev_idx)), flush=True)
        if va is not None:
            cands = [("net", None)] + ([("ema", net.ema_state)] if getattr(net, "ema_state", None) else [])
            for name, sd in cands:
                old = swap(sd) if sd is not None else None
                rv = epinn.evaluate(net, va)
                print("   " + report(f"[{it}] val {name}", rv), flush=True)
                if score(rv) < best["score"]:
                    best.update(score=score(rv), it=it, name=name,
                                state={k: v.detach().clone() for k, v in net.state_dict().items()})
                sm = float(np.mean([r["N_ritz"] / r["N_full"] - 1 for r in rv]))
                if sm < best_modes["score"]:
                    best_modes.update(score=sm, it=it, name=name,
                                      state={k: v.detach().clone() for k, v in net.state_dict().items()})
                if old is not None:
                    net.load_state_dict(old)
        if not a.no_test_log:
            print("   " + report(f"[{it}] test ", epinn.evaluate(net, te)), flush=True)
        print(f"   ({time.perf_counter() - t1:.0f} s)", flush=True)

    lg = lambda s: print(s + f"   ({time.perf_counter() - t1:.0f} s)", flush=True)
    if a.method.startswith("ngfull"):
        opt = epinn.train_ritz_ng_full(net, tr, outer=a.outer, epochs=a.epochs, batch=a.batch, lr=a.lr, k=a.k,
                                       seed=a.seed, log=lg, eval_fn=ev, eval_every=2, oracle=oracle, gamma=a.gamma,
                                       mode=a.mode, exact=exact, log_every=1 if a.mode != "data" else 5,
                                       lr_final=a.lr_final, ema=a.ema)
        if a.method.endswith("+adam"):
            epinn.train_adam(net, tr, steps=a.steps, batch=a.batch, lr=a.lr, k=a.k, seed=a.seed, log=lg,
                             eval_fn=ev, eval_every=200, opt=opt)
    elif a.method == "ritzng+adam":
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
                         oracle=oracle, gamma=a.gamma, n_oracle_batch=a.n_oracle_batch,
                         log=lambda s: print(s + f"   ({time.perf_counter() - t1:.0f} s)", flush=True),
                         eval_fn=ev, eval_every=100)
    state_last = {k: v.detach().clone() for k, v in net.state_dict().items()}
    if va is not None:
        ev("end")
        print(f"chosen on validation: {best['name']} weights of evaluation {best['it']} "
              f"(validation head mean {100 * best['score']:.2f} %); for the modes: {best_modes['name']} weights of "
              f"evaluation {best_modes['it']} (validation Ritz mean {100 * best_modes['score']:.2f} %)", flush=True)
        net.load_state_dict(best["state"])
    res = dict(train=epinn.evaluate(net, tr, idx=ev_idx), test=epinn.evaluate(net, te), args=vars(a),
               chosen=dict(it=str(best.get("it")), name=best.get("name"), val_head=best.get("score"),
                           modes_it=str(best_modes.get("it")), modes_name=best_modes.get("name"),
                           val_ritz=best_modes.get("score")))
    print(report("final train", res["train"])); print(report("final test ", res["test"]))
    for r in res["test"]:
        print(f"   test {r['sid']} bs {r['params']['bs']:.0f} beta {r['params']['beta']:.2f}: full {r['N_full']:.4f}  "
              f"Ritz {r['N_ritz']:.4f} ({100 * (r['N_ritz'] / r['N_full'] - 1):+.2f} %)  head {r['N_head']:.4f} "
              f"({100 * (r['N_head'] / r['N_full'] - 1):+.2f} %)")
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        torch.save(dict(state=net.state_dict(), state_last=state_last, state_modes=best_modes.get("state"),
                        state_ema=getattr(net, "ema_state", None), args=vars(a), chosen=res["chosen"]), a.out)
        json.dump(res, open(str(a.out).replace(".pt", ".json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
