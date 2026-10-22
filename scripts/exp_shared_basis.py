"""Experiment: accuracy of a shared GBT basis on unseen omega panels.

    python scripts/exp_shared_basis.py --side symmetry --n-test 16 --seed 77

Needs the snapshots of scripts/build_pod_snapshots.py.  For every TEST panel (a Sobol
set with another seed than the training panels, plus the validated panels of Paper 1):
  full     : GBT with all modes (the validated model)
  own80    : the panel's own hierarchical basis, 80 conventional modes (panel-specific)
  pod<n>s<k>m<j>: POD basis of the training buckling modes: n conventional + k shear-hat
             + j Mindlin-psi + j Mindlin-phi modes (conventional block energy-weighted)
All reduced solves use the full-basis pre-buckling state of the panel; errors are
relative to `full` (a reduced basis can only over-estimate).
Writes data/shared_basis/exp_<side>.jsonl and prints a summary table.
"""
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.basis import (NINT_REF, build_shared_basis, hx_of, pod_basis,   # noqa: E402
                            pod_energy, shared_section)
from pinnbuck.designspace import VALIDATED_LAYUPS, sample_design, validated_omega_sample  # noqa: E402
from pinnbuck.oracle import omega_length, omega_panel                          # noqa: E402
from pinnbuck.solver import configure                                          # noqa: E402

CFGS = ((40, 40, 20), (60, 60, 30), (80, 80, 40), (120, 60, 30))   # (n conv, n shear, n Mindlin per block)
TAG = lambda c: f"pod{c[0]}s{c[1]}m{c[2]}"


def test_panels(side, n, seed):
    out = sample_design("omega", n, seed=seed, categorical={"side": side})
    for lay in VALIDATED_LAYUPS["omega"]:                                       # Paper 1 validation panels
        for bs in (150.0, 285.5):
            out.append(validated_omega_sample(lay, bs, 1.5, "SS", side))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--n-test", type=int, default=12)
    ap.add_argument("--seed", type=int, default=77)
    ap.add_argument("--cfg", nargs="*", default=None, help="basis sizes n_conv,n_shear,n_mindlin (e.g. 120,60,30)")
    ap.add_argument("--no-own", action="store_true", help="skip the panel's own 80-mode basis")
    a = ap.parse_args()
    snaps = pickle.load(open(ROOT / "data" / "shared_basis" / f"snap_{a.side}.pkl", "rb"))
    print(f"{len(snaps)} training panels; share of the snapshot energy NOT captured by n modes:")
    for key, name in (("F", "conventional"), ("H", "shear hats"), ("X", "Mindlin psi"), ("S", "Mindlin phi")):
        e = pod_energy(snaps, key)
        print(f"  {name:13s} " + "  ".join(f"n={n}: {1 - e[min(n, e.size) - 1]:.1e}" for n in (10, 20, 40, 60, 80)))
    cfgs = CFGS if a.cfg is None else [tuple(int(v) for v in c.split(",")) for c in a.cfg]
    bases = {c: pod_basis(snaps, a.side, c[0], c[1], c[2]) for c in cfgs}
    out = ROOT / "data" / "shared_basis" / f"exp_{a.side}.jsonl"
    if out.exists():                                   # reuse the full-GBT values of earlier runs
        pass
    done = {json.loads(l)["sid"] for l in open(out)} if out.exists() else set()
    with open(out, "a") as f:
        for s in test_panels(a.side, a.n_test, a.seed):
            if s.sid in done:
                continue
            p = s.params; rec = dict(sid=s.sid, params=p, layups=s.layups)
            with configure(OMEGA_SIDE_BC=a.side, OMEGA_GBT_NINT=list(NINT_REF)) as m:
                sec = m.omega_gbt_section(omega_panel(s)); L = omega_length(s); hx = hx_of(sec, L)
                t = time.perf_counter()
                rec["full"] = m.omega_gbt_solve(sec, L, p["bc"], hx=hx, nev=2)["N"]
                rec["t_full"] = time.perf_counter() - t
                pre = m.omega_gbt_solve(sec, L, p["bc"], hx=hx, pre_only=True)["pre"]
                if not a.no_own:
                    own = m.truncated_section(sec, None, 0, 0, 0, hm=m.omega_gbt_modes(sec),
                                              coupled_shear=True, cidx=list(range(80)))
                    rec["own80"] = m.omega_gbt_solve(own, L, p["bc"], hx=hx, nev=2, pre=pre)["N"]
                for c in cfgs:
                    w = shared_section(m, sec, bases[c])
                    rec[TAG(c)] = m.omega_gbt_solve(w, L, p["bc"], hx=hx, nev=2, pre=pre)["N"]
                    rec["nm_" + TAG(c)] = int(w.nm)
            f.write(json.dumps(rec) + "\n"); f.flush()
            keys = [k for k in ["own80"] + [TAG(c) for c in cfgs] if k in rec]
            print(s.sid, {k: p[k] for k in ("bs", "beta", "bc")},
                  " ".join(f"{k} {100 * (rec[k] / rec['full'] - 1):+.2f}%" for k in keys), flush=True)
    recs = [json.loads(l) for l in open(out)]
    print(f"\nsummary over {len(recs)} test panels (error vs full GBT, %):")
    for k in ["own80"] + [TAG(c) for c in cfgs]:
        e = np.array([100 * (r[k] / r["full"] - 1) for r in recs if k in r])
        if not e.size:
            continue
        print(f"  {k:6s} mean {e.mean():6.2f}  max {e.max():6.2f}  min {e.min():6.2f}")


if __name__ == "__main__":
    main()
