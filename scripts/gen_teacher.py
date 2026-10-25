"""GBT teacher data for the buckling network (Paper 1 style oracle, here the validated GBT).

    python scripts/gen_teacher.py --side symmetry --set train --n 512 --seed 11 --workers 6
    python scripts/gen_teacher.py --side symmetry --set test  --n 64  --seed 99 --workers 6 --keep-ops
    python scripts/gen_teacher.py --family c --set train --n 256 --seed 11 --workers 6   (members)

Per design (Sobol sample of the omega design space): full GBT N_cr of the three lowest
modes, the modes in the network coordinates (`teacher.py`: projected and refined), and
the design vector.  One .npz per design in data/teacher/<side>/<set>/ (resumable).
--keep-ops also stores the panel operators (needed for the physics check; ~12 MB each).
Omega needs data/shared_basis/basis_<side>.pkl; members use their complete nodal basis.
About 20-60 s per design on one core.
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

J = 24


def _one(args):
    sample, side, out, keep_ops = args          # side: omega side condition or member family
    import numpy as np
    import torch
    torch.set_num_threads(1)
    from pinnbuck.operators import build_operator
    from pinnbuck.spectral import FourierPanel
    from pinnbuck.teacher import project, rayleigh, refine
    t0 = time.perf_counter()
    if sample.family == "omega":
        sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
    else:
        from pinnbuck.families import member_basis
        sb = member_basis(sample.family)
    op = build_operator(sample, sb, lam_full=None, with_shared_reference=False, n_teacher=3)
    t_gbt = time.perf_counter() - t0
    fp = FourierPanel(op, J=J)
    A0 = project(fp, op.teacher.astype(float))
    lam_proj = rayleigh(fp, A0).numpy()
    lam, M, hist = refine(fp, A0)
    np.savez_compressed(out, sid=sample.sid, params=json.dumps(sample.params), layups=json.dumps(sample.layups),
                        x=op.features.astype(np.float32), L=op.L, bc=op.bc, side=side,
                        lam_full=np.asarray(op.teacher_eigs, float), lam_proj=lam_proj, lam=lam.numpy(),
                        M=M.numpy().astype(np.float64), J=J, nm=op.nm, hist=np.asarray(hist))
    if keep_ops:
        op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32)
        pickle.dump(op, open(str(out).replace(".npz", ".op.pkl"), "wb"))
    return sample.sid, float(op.teacher_eigs[0]), float(lam[0]), float(lam_proj[0]), t_gbt, time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default="omega", choices=["omega", "c", "z", "l", "t", "i", "box"])
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--set", default="train")
    ap.add_argument("--n", type=int, default=512)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--keep-ops", action="store_true")
    a = ap.parse_args()
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(k, "1")
    from pinnbuck.designspace import sample_design
    tag = a.side if a.family == "omega" else a.family
    d = ROOT / "data" / "teacher" / tag / a.set; d.mkdir(parents=True, exist_ok=True)
    smp = (sample_design("omega", a.n, seed=a.seed, categorical={"side": a.side}) if a.family == "omega"
           else sample_design(a.family, a.n, seed=a.seed))
    todo = [(s, tag, d / f"{s.sid}.npz", a.keep_ops) for s in smp if not (d / f"{s.sid}.npz").exists()]
    print(f"{tag}/{a.set}: {len(todo)} designs to compute with {a.workers} workers", flush=True)
    t0 = time.perf_counter()
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for k, (sid, lf, lr, lp, tg, tt) in enumerate(pool.imap_unordered(_one, todo), 1):
            print(f"[{k}/{len(todo)}] {sid}  N_full {lf:.4f}  network-space {lr:.4f} ({100 * (lr / lf - 1):+.2f} %, "
                  f"projection only {100 * (lp / lf - 1):+.2f} %)  {tt:.0f} s  (elapsed {time.perf_counter() - t0:.0f} s)",
                  flush=True)


if __name__ == "__main__":
    main()
