"""GBT operators of omega panels in the shared POD basis (training / test data of the PINN).

    python scripts/build_operators.py --side symmetry --set train --n 256 --seed 11 --workers 6
    python scripts/build_operators.py --side symmetry --set test  --n 64  --seed 99 --workers 6 --full

Needs data/shared_basis/basis_<side>.pkl (scripts/make_shared_basis.py).
One file per panel: data/operators/<side>/<set>/<sid>.pkl (resumable).
--full also runs the full GBT (the validated model): its N_cr gives the error of the PINN and
its three lowest modes (projected on the shared basis) are the teacher of the first
training stage.
"""
import argparse
import multiprocessing as mp
import os
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _one(args):
    sample, side, full, out = args
    import numpy as np
    from pinnbuck.operators import build_operator
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
    t = time.perf_counter()
    op = build_operator(sample, sb, lam_full=None if full else float("nan"), with_shared_reference=False)
    op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32)
    op.info["time"] = time.perf_counter() - t
    pickle.dump(op, open(out, "wb"))
    return sample.sid, op.info["time"], op.lam_full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--set", default="train")
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(k, "1")
    from pinnbuck.designspace import sample_design
    d = ROOT / "data" / "operators" / a.side / a.set; d.mkdir(parents=True, exist_ok=True)
    todo = [(s, a.side, a.full, d / f"{s.sid}.pkl")
            for s in sample_design("omega", a.n, seed=a.seed, categorical={"side": a.side})
            if not (d / f"{s.sid}.pkl").exists()]
    print(f"{a.side}/{a.set}: {len(todo)} panels to build with {a.workers} workers", flush=True)
    t0 = time.perf_counter()
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for k, (sid, t, lam) in enumerate(pool.imap_unordered(_one, todo), 1):
            print(f"[{k}/{len(todo)}] {sid} {t:.1f} s  N_full {lam:.4f}  (elapsed {time.perf_counter() - t0:.0f} s)",
                  flush=True)


if __name__ == "__main__":
    main()
