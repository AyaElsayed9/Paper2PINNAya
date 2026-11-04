"""Step B0: are the shared bases and the Ritz space still accurate in the extended design space
(b_s from 124.25 mm as Paper 1, beta up to 3.5)?

    python scripts/check_extended_space.py --side periodic --workers 2

For every check design three GBT values of the lowest buckling load:
    full    full GBT (the validated model, all cross-section modes)
    shared  GBT in the shared POD basis (241 modes)
    ritz    converged eigenvalue in the network's Ritz space (shared basis x longitudinal dictionary),
            which is the label / the best value the network can reach
and the errors shared/full - 1 and ritz/full - 1.  Designs: corners and edges of the new region,
a few random (Sobol) points inside it, and control points inside the old region (b_s >= 150, beta <= 3)
where the errors are known (0.01-0.02 % basis, 0.05-0.3 % Ritz).  Lay-ups: the named validated ones
and random industrial ones.  Results: results/b0_check/<side>.jsonl (resumable).
"""
import argparse
import json
import multiprocessing as mp
import os
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def designs(side):
    import numpy as np
    from scipy.stats import qmc
    from pinnbuck.designspace import FAMILIES, Sample, industrial_layups, validated_omega_sample
    out = []
    # corners / edges of the new region and controls, named lay-ups
    pts = [(124.25, 0.75, "Cross-ply", "SS"), (124.25, 1.5, "Cross-ply", "SS"), (124.25, 3.5, "Cross-ply", "SS"),
           (124.25, 1.5, "Unsym", "SS"), (124.25, 3.5, "Sym1", "SS"), (124.25, 0.75, "Sym2", "CC"),
           (124.25, 3.5, "Unsym", "CC"), (285.5, 3.5, "Cross-ply", "SS"), (285.5, 3.5, "Sym1", "CC"),
           (150.0, 1.5, "Sym1", "SS"), (150.0, 1.5, "Unsym", "CC"), (285.5, 1.5, "Unsym", "SS")]
    for bs, beta, lay, bc in pts:
        out.append(("corner" if (bs < 150 or beta > 3) else "control", validated_omega_sample(lay, bs, beta, bc, side)))
    # random points inside the new region, industrial lay-ups
    rng = np.random.default_rng(7)
    ind8, ind7 = industrial_layups(8), industrial_layups(7)
    sob = qmc.Sobol(d=2, scramble=True, seed=31).random(8)
    fx = dict(FAMILIES["omega"]["fixed"])
    for i, u in enumerate(sob):
        if i < 6:
            bs, beta = 124.25 + u[0] * (150 - 124.25), 0.75 + u[1] * (3.5 - 0.75)
        else:
            bs, beta = 150 + u[0] * (285.5 - 150), 3.0 + u[1] * 0.5
        lay = {"skin": list(ind8[rng.integers(len(ind8))]), "stringer": list(ind7[rng.integers(len(ind7))])}
        out.append(("random", Sample("omega", dict(bs=round(float(bs), 2), beta=round(float(beta), 3),
                                                   bc=("SS", "CC")[i % 2], side=side), lay, fx)))
    return out


def _one(args):
    kind, sample, side = args
    import numpy as np
    import torch
    torch.set_num_threads(1)
    from pinnbuck.operators import build_operator
    sys.path.insert(0, str(ROOT / "scripts"))
    from exact_targets import exact_modes
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
    t = time.perf_counter()
    op = build_operator(sample, sb, lam_full=None, with_shared_reference=True)
    t_op = time.perf_counter() - t
    t = time.perf_counter()
    ex = exact_modes(op, R=4)
    t_ritz = time.perf_counter() - t
    p = sample.params
    return dict(kind=kind, sid=sample.sid, bs=p["bs"], beta=p["beta"], bc=p["bc"], side=side, layups=sample.layups,
                full=op.lam_full, shared=op.lam_shared, ritz=float(ex["N_ritz"][0]), ritz_iters=ex["iters"],
                err_shared=100 * (op.lam_shared / op.lam_full - 1), err_ritz=100 * (ex["N_ritz"][0] / op.lam_full - 1),
                nm=op.nm, L=op.L, t_op=t_op, t_ritz=t_ritz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="periodic", choices=["symmetry", "periodic"])
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[k] = "1"
    out = ROOT / "results" / "b0_check" / f"{a.side}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    have = {json.loads(l)["sid"] for l in open(out)} if out.exists() else set()
    todo = [(k, s, a.side) for k, s in designs(a.side) if s.sid not in have]
    print(f"{a.side}: {len(have)} done, {len(todo)} to go", flush=True)
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for r in pool.imap_unordered(_one, todo):
            with open(out, "a") as f:
                f.write(json.dumps(r) + "\n")
            print(f"{r['kind']:7s} bs {r['bs']:7.2f} beta {r['beta']:5.2f} {r['bc']} full {r['full']:9.4f} "
                  f"shared {r['err_shared']:+7.3f} %  ritz {r['err_ritz']:+7.3f} % ({r['ritz_iters']} it)  "
                  f"{r['t_op']:.0f}+{r['t_ritz']:.0f} s", flush=True)


if __name__ == "__main__":
    main()
