"""Step B1 check: is the full-GBT label converged along the length?

The end-layer probe (`b1_end_layer_probe.py`) gave Ritz values up to 0.15 % BELOW the full GBT value when the
end splines were finer than the solver's longitudinal span hx.  The Ritz space is not a subspace of the
solver's spline space (other functions along z), so this can mean either
  (a) the full GBT (span hx, as the solver chooses it) is not converged along z at the 0.1 % level, or
  (b) an artefact of the Ritz integration at the solver's Gauss points.
Here the full GBT itself (frozen solver, only its hx argument changed) is rerun with hx/2 and hx/4.
If N drops by about the same amount, (a) is true.

    python scripts/b1_hx_convergence.py data/operators/periodic_mixed_SS_v2/test --sids 040508c6cbc1 ...
Results: results/b1_end_layer_probe/hx_convergence.jsonl (resumable)
"""
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pinnbuck.designspace import Sample            # noqa: E402
from pinnbuck import opio                            # noqa: E402
from pinnbuck.families import family_context        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--sids", nargs="+", required=True)
    ap.add_argument("--factors", nargs="+", type=float, default=[1.0, 0.5, 0.25])
    a = ap.parse_args()
    out = ROOT / "results" / "b1_end_layer_probe" / "hx_convergence.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {(j["sid"], j["factor"]) for j in map(json.loads, open(out))} if out.exists() else set()
    for sid in a.sids:
        op = opio.load(Path(a.src) / f"{sid}.pkl")
        d = dict(op.info["sample"]); d.pop("sid")
        smp = Sample(**d)
        assert smp.sid == sid, (smp.sid, sid)
        stored = op.lam_full
        del op
        for fac in a.factors:
            if (sid, fac) in done:
                continue
            with family_context(smp) as fc:
                hx = fc.hx * fac
                t = time.perf_counter()
                r = fc.m.omega_gbt_solve(fc.sec, fc.L, fc.bc, hx=hx, nev=1)
                dt = time.perf_counter() - t
            N = float(r["N"])
            rec = dict(sid=sid, factor=fac, hx=hx, L=fc.L, bs=smp.params["bs"], N=N, stored=stored,
                       diff_vs_stored=100 * (N / stored - 1), t=dt)
            print(f"{sid} hx x{fac:<5} ({hx:6.2f} mm)  N {N:9.5f}  vs label {rec['diff_vs_stored']:+.4f} %  {dt:.0f} s",
                  flush=True)
            with open(out, "a") as f:
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
