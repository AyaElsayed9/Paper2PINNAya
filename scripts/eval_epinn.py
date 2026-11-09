"""Evaluate a trained energy-whitened PINN on test designs, with the physics certificate.

    python scripts/eval_epinn.py runs/t3_cc259.pt --test data/operators/symmetry_industrial/test

Per design: N_cr of the head (instant, no operators), the Ritz upper bound of the
predicted modes, and the certificate: block LOBPCG with the energy metric started from the
predicted modes (each step ~0.3 s once the operators exist).  Every Ritz value is an upper
bound of the GBT critical load in the Ritz space; the relative change of the last step and
the energy residual tell how far it is from convergence.
"""
import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn
from pinnbuck import opio                        # noqa: E402                       # noqa: E402
from pinnbuck.ritzspace import EnergyMetric      # noqa: E402
from pinnbuck.spectral import lobpcg             # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--test", default=str(ROOT / "data/operators/symmetry_industrial/test"))
    ap.add_argument("--steps", default="3,5,8")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--modes", default="state", help="weights for the mode shapes / certificate start: state (the "
                    "head's, default), state_modes (chosen on the validation Ritz bound), state_last, state_ema")
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    ck = torch.load(a.model, weights_only=False); ar = ck["args"]
    epinn.INPUT_COLS = None if ar.get("inputs", "all") == "all" else [0, 1]
    steps = [int(s) for s in a.steps.split(",")]
    files = [f for f in sorted(Path(a.test).glob("*.pkl")) if opio.load(f).params["bc"] == ar["bc"]]
    rows = []
    net = None
    for f in files:
        d = epinn.prepare(opio.load(f))
        if net is None:
            net = epinn.net_from_ckpt(ck, d)
            net_m = net
            if a.modes != "state":                     # second copy of the network for the mode shapes
                import copy
                net_m = copy.deepcopy(net); net_m.load_state_dict(ck[a.modes], strict=False); net_m.eval()
                print(f"head: weights 'state'; mode shapes: weights '{a.modes}'", flush=True)
        with torch.no_grad():
            _, logN = net(torch.as_tensor(d.x[None]))
            Y, _ = net_m(torch.as_tensor(d.x[None]))
            Y = Y[:, :, d.sel]                                     # the design's part of the fixed layout
            mu, _ = epinn.ritz(d, Y[0])
        em = EnergyMetric(d.rs, groups="labels", ridge=1e-8)
        lam, _, h = lobpcg(d.rs, R=ar["R"], iters=max(steps), X0=d.rs.unflat(d.fm.color(Y[0].double())),
                           precond=em.solve, tol=0, seed=1)
        e = lambda v: 100 * (v / d.lam_full - 1)
        r = dict(head=e(float(torch.exp(logN[0, 0]))), ritz=e(1.0 / float(mu[0])), **{f"s{s}": e(h[s - 1]) for s in steps})
        rows.append(r)
        print(f"{d.sid} bs {d.params['bs']:.0f} beta {d.params['beta']:.2f}: full GBT {d.lam_full:.3f} | head {r['head']:+.2f} % | "
              f"network Ritz {r['ritz']:+.2f} % | certificate " + "  ".join(f"{s} steps {r[f's{s}']:+.2f} %" for s in steps),
              flush=True)
    print(f"\n{len(rows)} test designs ({ar['bc']}), |error| vs full GBT, mean / max:")
    for key in ["head", "ritz"] + [f"s{s}" for s in steps]:
        v = np.abs([r[key] for r in rows])
        print(f"  {key:5s} {np.mean(v):6.2f} % / {np.max(v):6.2f} %")


if __name__ == "__main__":
    main()
