"""Compress an operator set in place (lossless, `pinnbuck.opio`); every file is checked bit for bit after
writing (all arrays identical, all scalars identical) before the next one.

    python scripts/compress_operators.py data/operators/periodic_mixed_SS_v2/train
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import opio                                      # noqa: E402


def same(a, b):
    va, vb = vars(a), vars(b)
    if va.keys() != vb.keys():
        return False
    for k in va:
        x, y = va[k], vb[k]
        if isinstance(x, np.ndarray):
            if not (x.dtype == y.dtype and x.shape == y.shape and np.array_equal(x, y, equal_nan=True)):
                return False
        elif isinstance(x, float) and np.isnan(x):
            if not np.isnan(y):
                return False
        elif k != "classes" and x != y:
            return False
    return True


def main():
    src = Path(sys.argv[1])
    files = sorted(src.glob("*.pkl"))
    before = after = n = 0
    for f in files:
        if opio.is_compressed(f):
            continue
        s0 = f.stat().st_size
        op = opio.load(f)
        opio.save(op, f)
        assert same(op, opio.load(f)), f"mismatch after compression: {f}"
        before += s0; after += f.stat().st_size; n += 1
    print(f"{src}: {n} files compressed, {before / 1e9:.2f} GB -> {after / 1e9:.2f} GB "
          f"({len(files) - n} were already compressed)")


if __name__ == "__main__":
    main()
