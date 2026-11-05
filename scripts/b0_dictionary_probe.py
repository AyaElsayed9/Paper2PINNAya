"""Step B0 follow-up 2: dictionary variants on the worst check designs of both side conditions.

The first probe (`b0_long_cc_probe.py`) showed that the 5.3 % Ritz error of the long clamped panel comes from
the dictionary (longest half-wave 350 mm), not from LOBPCG.  Here, on every B0 check design with a Ritz error
above 0.3 % plus the long clamped corner and two controls, four Ritz spaces are compared at 400 iterations:
    base       half-waves 25 .. 350 mm (6, geometric)                    -- current
    long       + 593, 1006 mm (next two of the same series)
    longshort  + 14.75 mm (the one before 25 mm) + 593, 1006 mm
and base with the old 80-iteration cap.

    python scripts/b0_dictionary_probe.py --side symmetry
Results: results/b0_check/dictionary_probe_<side>.jsonl
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from b0_long_cc_probe import HW_LONG, ritz                     # noqa: E402
from pinnbuck.designspace import FAMILIES, Sample              # noqa: E402
from pinnbuck.dictionary import HALF_WAVES                     # noqa: E402
from pinnbuck.operators import build_operator                  # noqa: E402

HW_LONGSHORT = (14.75,) + tuple(HW_LONG)
VARIANTS = (("base80", HALF_WAVES, 80), ("base400", HALF_WAVES, 400), ("long400", HW_LONG, 400),
            ("longshort400", HW_LONGSHORT, 400))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry")
    ap.add_argument("--min-err", type=float, default=0.3)
    a = ap.parse_args()
    torch.set_num_threads(2)
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{a.side}.pkl", "rb"))
    chk = [json.loads(l) for l in open(ROOT / "results" / "b0_check" / f"{a.side}.jsonl")]
    pick = [r for r in chk if r["err_ritz"] > a.min_err or r["kind"] == "control"]
    out = ROOT / "results" / "b0_check" / f"dictionary_probe_{a.side}.jsonl"
    done = {json.loads(l)["sid"] for l in open(out)} if out.exists() else set()
    for c in pick:
        if c["sid"] in done:
            continue
        smp = Sample("omega", dict(bs=c["bs"], beta=c["beta"], bc=c["bc"], side=a.side), c["layups"],
                     dict(FAMILIES["omega"]["fixed"]))
        assert smp.sid == c["sid"], (smp.sid, c["sid"])
        op = build_operator(smp, sb, lam_full=None, with_shared_reference=False)
        r = dict(sid=c["sid"], kind=c["kind"], bs=c["bs"], beta=c["beta"], bc=c["bc"], side=a.side, L=op.L,
                 full=op.lam_full, check_err_ritz=c["err_ritz"])
        for tag, hw, it in VARIANTS:
            N, n_it, n, dt = ritz(op, hw, it)
            r[tag] = dict(N=N, err=100 * (N / op.lam_full - 1), iters=n_it, n=n, t=dt)
            print(f"{c['kind']:7s} bs {c['bs']:7.2f} beta {c['beta']:5.2f} {c['bc']} L {op.L:5.0f}  {tag:12s} "
                  f"N {N:9.4f} ({r[tag]['err']:+7.3f} %)  {n_it:3d} it  n {n}  {dt:.0f} s", flush=True)
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
        del op


if __name__ == "__main__":
    main()
