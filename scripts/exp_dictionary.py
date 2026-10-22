"""Experiment: the longitudinal dictionary (bulk waves + end layers, physical units)
with exact cross-section freedom, against the full GBT.

    python scripts/exp_dictionary.py --side symmetry --set test

Needs operator files with the full-GBT reference (scripts/build_operators.py --full).
Each panel takes 2-4 minutes (dense Ritz problem of about 6000 unknowns).
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pinnbuck.operators as O                     # noqa: E402
from pinnbuck.dictionary import dictionary_classes  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--set", default="test")
    ap.add_argument("--max", type=int, default=8)
    a = ap.parse_args()
    kg = O.kg_svd
    O.kg_svd = lambda op, tol=1e-10: kg(op, tol=1e-6)   # pre-stress SVD at 1e-6 (dense Ritz cost)
    for f in sorted((ROOT / "data" / "operators" / a.side / a.set).glob("*.pkl"))[:a.max]:
        op = pickle.load(open(f, "rb"))
        t = time.perf_counter()
        B = dictionary_classes(op)
        lam = O.ritz_lambda(op, B)
        n = sum(c.Pi.shape[1] * b.shape[2] for c, b in zip(op.classes, B))
        print(f"{op.params}  full {op.lam_full:.4f}  dictionary {lam:.4f}  ({100 * (lam / op.lam_full - 1):+.2f} %)"
              f"  unknowns {n}  {time.perf_counter() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()
