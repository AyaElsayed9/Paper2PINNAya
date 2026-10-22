"""Quick check of the shared GBT basis (a few minutes on a laptop).

    python scripts/check_shared_basis.py
    python scripts/check_shared_basis.py --side periodic --layup Sym1

Uses the committed basis data/shared_basis/basis_<side>.pkl (241 shared modes) on the
validated Paper 1 panels (b_s = 150 and 285.5 mm, beta = 1.5) for SS and CC ends:
full GBT (the validated model, all modes) against GBT restricted to the shared basis.
Expected: differences of a few hundredths of a percent.
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.basis import NINT_REF, hx_of, shared_section          # noqa: E402
from pinnbuck.designspace import validated_omega_sample              # noqa: E402
from pinnbuck.oracle import omega_length, omega_panel                # noqa: E402
from pinnbuck.solver import configure                                # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--layup", default="Unsym", choices=["Cross-ply", "Sym1", "Sym2", "Unsym"])
    a = ap.parse_args()
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{a.side}.pkl", "rb"))
    print(f"side {a.side}, lay-up {a.layup}")
    print(f"{'bs':>6s} {'bc':>3s} {'full GBT':>10s} {'shared':>10s} {'diff':>8s} {'modes':>11s} {'time':>12s}")
    for bc in ("SS", "CC"):
        for bs in (150.0, 285.5):
            s = validated_omega_sample(a.layup, bs, 1.5, bc, a.side)
            with configure(OMEGA_SIDE_BC=a.side, OMEGA_GBT_NINT=list(NINT_REF)) as m:
                sec = m.omega_gbt_section(omega_panel(s)); L = omega_length(s); hx = hx_of(sec, L)
                t = time.perf_counter()
                full = m.omega_gbt_solve(sec, L, bc, hx=hx, nev=1)
                t_full = time.perf_counter() - t
                w = shared_section(m, sec, sb)
                t = time.perf_counter()
                red = m.omega_gbt_solve(w, L, bc, hx=hx, nev=1, pre=full["pre"])
                t_red = time.perf_counter() - t
            print(f"{bs:6.1f} {bc:>3s} {full['N']:10.4f} {red['N']:10.4f} {100 * (red['N'] / full['N'] - 1):+7.3f}% "
                  f"{sec.nm:4d} -> {w.nm:3d} {t_full:5.1f} / {t_red:4.1f} s", flush=True)


if __name__ == "__main__":
    main()
