"""Full-GBT buckling modes of training panels for the shared POD basis.

    python scripts/build_pod_snapshots.py --side symmetry --n 32 --seed 2026 --workers 2

Writes data/shared_basis/snap_<side>.npz (resumable: panels already in the file are skipped).
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


def _one(sample):
    from pinnbuck.basis import panel_snapshots
    t = time.perf_counter()
    s = panel_snapshots(sample, nev=3)
    s["params"] = sample.params; s["time"] = time.perf_counter() - t
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    from pinnbuck.designspace import sample_design
    out = ROOT / "data" / "shared_basis" / f"snap_{a.side}.pkl"
    out.parent.mkdir(parents=True, exist_ok=True)
    snaps = pickle.load(open(out, "rb")) if out.exists() else []
    have = {s["sid"] for s in snaps}
    todo = [s for s in sample_design("omega", a.n, seed=a.seed, categorical={"side": a.side}) if s.sid not in have]
    print(f"{len(have)} panels done, {len(todo)} to go", flush=True)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[k] = str(a.threads)
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for s in pool.imap_unordered(_one, todo):
            snaps.append(s)
            pickle.dump(snaps, open(out, "wb"))
            print(f"[{len(snaps)}] {s['sid']} {s['params']} N_cr {s['eigs'][0]:.3f}  {s['time']:.1f} s", flush=True)


if __name__ == "__main__":
    main()
