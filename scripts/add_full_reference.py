"""Add the full-GBT N_cr (the validated model) to operator files that lack it.

    python scripts/add_full_reference.py --side symmetry --set train --seed 11 --n 256 --workers 6

Training panels do not need it (the physics loss has no labels); it is used to report
the training error and for the data term of the active learning.
"""
import argparse
import math
import multiprocessing as mp
import os
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _one(args):
    f, sample = args
    from pinnbuck.basis import NINT_REF, hx_of
    from pinnbuck.oracle import omega_length, omega_panel
    from pinnbuck.solver import configure
    op = pickle.load(open(f, "rb"))
    if not math.isfinite(op.lam_full):
        p = sample.params
        with configure(OMEGA_SIDE_BC=p["side"], OMEGA_GBT_NINT=list(NINT_REF)) as m:
            sec = m.omega_gbt_section(omega_panel(sample)); L = omega_length(sample)
            op.lam_full = m.omega_gbt_solve(sec, L, p["bc"], hx=hx_of(sec, L), nev=2)["N"]
        pickle.dump(op, open(f, "wb"))
    return f.name, op.lam_full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--set", default="train")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(k, "1")
    from pinnbuck.designspace import sample_design
    smp = {s.sid: s for s in sample_design("omega", a.n, seed=a.seed, categorical={"side": a.side})}
    d = ROOT / "data" / "operators" / a.side / a.set
    todo = [(f, smp[f.stem]) for f in sorted(d.glob("*.pkl")) if f.stem in smp]
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for k, (name, lam) in enumerate(pool.imap_unordered(_one, todo), 1):
            print(f"[{k}/{len(todo)}] {name} N_full {lam:.4f}", flush=True)


if __name__ == "__main__":
    main()
