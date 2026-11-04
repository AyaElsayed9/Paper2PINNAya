"""Step B0 follow-up: why is the Ritz value 5.3 % too high for the long clamped panel
(b_s 285.5, beta 3.5, CC, L ~ 1000 mm)?

Two candidate causes, tested separately on a few long panels:
  (1) LOBPCG not converged (all B0 runs stopped at the 80-iteration cap)  -> rerun with 400 iterations;
  (2) the longitudinal dictionary is too short (longest half-wave 350 mm)  -> add the next two
      half-waves of the same geometric series (593 and 1006 mm), fixed physical lengths, so the
      dictionary labels stay the same for every design.

    python scripts/b0_long_cc_probe.py --side periodic
Results: results/b0_check/long_cc_probe_<side>.jsonl
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

from pinnbuck.designspace import validated_omega_sample        # noqa: E402
from pinnbuck.dictionary import HALF_WAVES                     # noqa: E402
from pinnbuck.operators import build_operator                  # noqa: E402
from pinnbuck.ritzspace import RitzSpace                       # noqa: E402
from pinnbuck.spectral import lobpcg                           # noqa: E402
from pinnbuck.ritzspace import EnergyMetric                     # noqa: E402

HW_LONG = tuple(HALF_WAVES) + (593.0, 1006.0)


def ritz(op, half_waves, iters):
    rs = RitzSpace(op, raw=True, half_waves=half_waves)
    em = EnergyMetric(rs, groups="labels", ridge=1e-8)
    t = time.perf_counter()
    N, X, hist = lobpcg(rs, R=4, iters=iters, precond=em.solve, tol=1e-10, seed=1)
    return float(N[0]), len(hist), rs.n, time.perf_counter() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="periodic")
    a = ap.parse_args()
    torch.set_num_threads(2)
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{a.side}.pkl", "rb"))
    cases = [("Sym1", 285.5, 3.5, "CC"), ("Sym1", 285.5, 3.0, "CC"), ("Unsym", 285.5, 3.0, "CC"),
             ("Sym1", 285.5, 2.0, "CC"), ("Cross-ply", 285.5, 3.5, "SS"), ("Unsym", 124.25, 3.5, "CC")]
    out = ROOT / "results" / "b0_check" / f"long_cc_probe_{a.side}.jsonl"
    done = {(j["lay"], j["bs"], j["beta"], j["bc"]) for j in map(json.loads, open(out))} if out.exists() else set()
    for lay, bs, beta, bc in cases:
        if (lay, bs, beta, bc) in done:
            continue
        smp = validated_omega_sample(lay, bs, beta, bc, a.side)
        op = build_operator(smp, sb, lam_full=None, with_shared_reference=False)
        r = dict(lay=lay, bs=bs, beta=beta, bc=bc, side=a.side, L=op.L, full=op.lam_full)
        for tag, hw, it in (("base80", HALF_WAVES, 80), ("base400", HALF_WAVES, 400), ("long400", HW_LONG, 400)):
            N, n_it, n, dt = ritz(op, hw, it)
            r[tag] = dict(N=N, err=100 * (N / op.lam_full - 1), iters=n_it, n=n, t=dt)
            print(f"{lay:9s} bs {bs:6.2f} beta {beta:4.2f} {bc} L {op.L:6.0f}  {tag:8s} N {N:8.4f} "
                  f"({r[tag]['err']:+7.3f} %)  {n_it:3d} it  n {n}  {dt:.0f} s", flush=True)
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
        del op


if __name__ == "__main__":
    main()
