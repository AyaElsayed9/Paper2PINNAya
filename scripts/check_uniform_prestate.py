"""16 Nov 2026: how much of the pre-buckling state does the check really need?

The fixed cost of the check is mostly the pre-buckling state of the new design (full-GBT static solve, 2.4-4.3 s of
~6.4 s).  Test: replace the computed membrane forces P(z, s) by their mid-length values at every z (a state that
needs no solve along the length) and converge the Ritz-space eigenproblem again (LOBPCG from the exact modes).
Error vs the converged load with the true state.  Labelled training designs of one expert.
    python scripts/check_uniform_prestate.py symmetry_crossply_SS_v2 --n 20
-> results/check_baselines/<tag>_uniform_prestate.csv"""
import argparse, csv, sys
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio
from pinnbuck.spectral import lobpcg


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--iters", type=int, default=60); ap.add_argument("--inputs", default="geom")
    a = ap.parse_args(); torch.set_num_threads(2)
    epinn.INPUT_COLS = None if a.inputs == "all" else [0, 1]
    files = [l.strip() for l in open(ROOT / f"runs/experts/{a.tag}_files.txt") if l.strip()][:a.n]
    rows = []
    for f in files:
        op = opio.load(ROOT / f)
        ex = torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{Path(f).stem}.pt", weights_only=False)
        P = np.array(op.P, dtype=float); nz = P.shape[1]; mid = nz // 2
        dev = [float(np.abs(P[i] - P[i][mid][None]).max() / max(np.abs(P[0]).max(), 1e-300)) for i in range(3)]
        op.P = np.repeat(P[:, mid:mid + 1, :], nz, axis=1).astype(op.P.dtype)
        d = epinn.prepare(op)
        N, _, h = lobpcg(d.rs, R=8, iters=a.iters, X0=d.rs.unflat(d.fm.color(ex["T"].double())), precond=d.em.solve, tol=1e-9, seed=1)
        e = 100 * (float(N[0]) / float(ex["N_ritz"][0]) - 1)
        rows.append(dict(sid=op.sid, bs=round(op.params["bs"], 1), beta=round(op.params["beta"], 3), N_true=round(float(ex["N_ritz"][0]), 4),
                         N_uniform=round(float(N[0]), 4), err_pct=round(e, 3), iters=len(h),
                         max_dev_Pxx=round(dev[0], 3), max_dev_Pss=round(dev[1], 3), max_dev_Pxs=round(dev[2], 3)))
        print(rows[-1], flush=True); del op, d
    with open(ROOT / f"results/check_baselines/{a.tag}_uniform_prestate.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    e = np.array([r["err_pct"] for r in rows])
    print(f"\n{a.tag}: N_cr with the mid-length pre-buckling state at every z vs the true state, {len(rows)} designs: "
          f"mean |e| {np.abs(e).mean():.2f} %, max {np.abs(e).max():.2f} %, range {e.min():+.2f} .. {e.max():+.2f} %")


if __name__ == "__main__":
    main()
