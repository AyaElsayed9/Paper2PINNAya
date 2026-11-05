"""Step B1 check: do the end layers of the dictionary need to be longer for general lay-ups?

On operators that are already built (full GBT inside), the converged Ritz value (400 iterations) with the d2
waves and three end-layer choices:
    e8     8 splines over b_s           (current)
    e12   12 splines over 1.5 b_s       (same spline spacing, longer zone)
    e16   16 splines over 2.5 b_s       (longer zone again)
(the zone is capped at L/2).

    python scripts/b1_end_layer_probe.py data/operators/periodic_mixed_SS_v2/test --n 24 [--set spacing]
Results: results/b1_end_layer_probe/<set>.jsonl (resumable)
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

from b0_long_cc_probe import ritz                              # noqa: E402
from pinnbuck.dictionary import DICTIONARIES                   # noqa: E402
from pinnbuck import opio                                      # noqa: E402
from pinnbuck.designspace import industrial_layups             # noqa: E402

VARIANTS = (("e8", {}), ("e12", dict(n_end=12, end_length=1.5)), ("e16", dict(n_end=16, end_length=2.5)))
# second probe (--set spacing): the same zones but the spline spacing never below c * hx (solver span), so the
# splines are integrated correctly at the solver's Gauss points (`dictionary.end_spacing`)
VARIANTS_SPACING = (("e8s15", dict(min_spacing=1.5)),
                    ("e16s15", dict(n_end=16, end_length=2.5, min_spacing=1.5)),
                    ("e16s2", dict(n_end=16, end_length=2.5, min_spacing=2.0)))
SETS = {"zone": (VARIANTS, ""), "spacing": (VARIANTS_SPACING, "_spacing")}


def is_industrial(lay):
    ind8 = {tuple(map(float, x)) for x in industrial_layups(8)}
    ind7 = {tuple(map(float, x)) for x in industrial_layups(7)}
    return tuple(map(float, lay["skin"])) in ind8 and tuple(map(float, lay["stringer"])) in ind7


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--set", default="zone", choices=list(SETS))
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    src = Path(a.src)
    variants, suffix = SETS[a.set]
    out = ROOT / "results" / "b1_end_layer_probe" / f"{src.parent.name}_{src.name}{suffix}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {json.loads(l)["sid"] for l in open(out)} if out.exists() else set()
    files = sorted(src.glob("*.pkl"))[:a.n]
    for f in files:
        if f.stem in done:
            continue
        op = opio.load(f)
        smp = op.info["sample"]
        r = dict(sid=f.stem, bs=op.params["bs"], beta=op.params["beta"], bc=op.params["bc"], L=op.L,
                 full=op.lam_full, industrial=is_industrial(smp["layups"]), layups=smp["layups"])
        for tag, kw in variants:
            N, n_it, n, dt = ritz(op, DICTIONARIES["d2"], 400, **kw)
            r[tag] = dict(N=N, err=100 * (N / op.lam_full - 1), iters=n_it, n=n, t=dt)
        print(f"{f.stem} bs {r['bs']:6.1f} beta {r['beta']:4.2f} L {op.L:5.0f} {'ind' if r['industrial'] else 'gen'}  "
              + "  ".join(f"{t} {r[t]['err']:+.3f} % (n {r[t]['n']})" for t, _ in variants), flush=True)
        with open(out, "a") as fo:
            fo.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()
