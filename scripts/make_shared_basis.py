"""Shared POD basis from the snapshot modes (scripts/build_pod_snapshots.py).

    python scripts/make_shared_basis.py --side symmetry --n-conv 120 --n-shear 60 --n-mindlin 30

Writes data/shared_basis/basis_<side>.pkl, used by build_operators.py and the PINN.
"""
import argparse
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.basis import pod_basis, pod_energy   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--n-conv", type=int, default=120)
    ap.add_argument("--n-shear", type=int, default=60)
    ap.add_argument("--n-mindlin", type=int, default=30)
    a = ap.parse_args()
    snaps = pickle.load(open(ROOT / "data" / "shared_basis" / f"snap_{a.side}.pkl", "rb"))
    sb = pod_basis(snaps, a.side, a.n_conv, a.n_shear, a.n_mindlin)
    out = ROOT / "data" / "shared_basis" / f"basis_{a.side}.pkl"
    pickle.dump(sb, open(out, "wb"))
    nm = sb.nodal_C.shape[1] + 1 + sb.Ts.shape[1] + sb.nodal_mx.shape[1] + sb.nodal_ms.shape[1]
    print(f"{out.name}: {len(snaps)} training panels, {nm} shared modes "
          f"({sb.nodal_C.shape[1]} conventional, 1 axial, {sb.Ts.shape[1]} shear, "
          f"{sb.nodal_mx.shape[1]} + {sb.nodal_ms.shape[1]} Mindlin)")
    for key, name in (("F", "conventional"), ("H", "shear hats"), ("X", "Mindlin psi"), ("S", "Mindlin phi")):
        e = pod_energy(snaps, key); n = {"F": a.n_conv, "H": a.n_shear}.get(key, a.n_mindlin)
        print(f"  {name:13s} snapshot energy not captured: {1 - e[min(n, e.size) - 1]:.1e}")


if __name__ == "__main__":
    main()
