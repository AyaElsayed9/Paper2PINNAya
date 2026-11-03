"""Exact GBT modes of the training designs in the network's coordinates (three-way test, 3 Nov 2026).

For every design: block LOBPCG with the energy metric in the design's Ritz space, run to convergence
(R = 8 vectors), stored as unit whitened vectors (the coordinates the network predicts), together with
the converged Ritz loads and the full-GBT loads stored in the operator file.

    python scripts/exact_targets.py data/operators/periodic_industrial_SS/train data/exact_targets/periodic_industrial_SS \
        --part 0 --parts 2
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn                       # noqa: E402
from pinnbuck.ritzspace import EnergyMetric      # noqa: E402
from pinnbuck.spectral import lobpcg             # noqa: E402


def exact_modes(op, R=8, iters=80, tol=1e-10):
    d = epinn.prepare(op)
    em = EnergyMetric(d.rs, groups="labels", ridge=1e-8)
    N, X, hist = lobpcg(d.rs, R=R, iters=iters, precond=em.solve, tol=tol, seed=1)
    Y = d.fm.whiten(d.rs.flat(X))
    Y = Y / Y.norm(dim=1, keepdim=True)
    return dict(sid=d.sid, T=Y.double(), N_ritz=np.asarray(N, float), N_full=float(op.lam_full),
                eigs_full=np.asarray(op.teacher_eigs, float), iters=len(hist))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--part", type=int, default=0); ap.add_argument("--parts", type=int, default=1)
    ap.add_argument("--R", type=int, default=8); ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    out = Path(a.dst); out.mkdir(parents=True, exist_ok=True)
    files = sorted(Path(a.src).glob("*.pkl"))[a.part::a.parts]
    for f in files:
        if (out / f"{f.stem}.pt").exists():
            continue
        t0 = time.perf_counter()
        r = exact_modes(pickle.load(open(f, "rb")), R=a.R)
        torch.save(r, out / f"{f.stem}.pt")
        print(f"{r['sid']}: Ritz N1 {r['N_ritz'][0]:.4f}, full GBT {r['N_full']:.4f} "
              f"({100 * (r['N_ritz'][0] / r['N_full'] - 1):+.3f} %), {r['iters']} iterations, {time.perf_counter() - t0:.1f} s",
              flush=True)


if __name__ == "__main__":
    main()
