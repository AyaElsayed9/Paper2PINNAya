"""Step B1 check: why do a few designs have a large Ritz error (test set: 5cc7be690865 +1.53 %, 6cc2e9904ad4 +0.47 %)?

Two possible causes:
  (a) the longitudinal dictionary (11 half-waves, ratio 1.45 between neighbours) misses the mode's wavelength;
  (b) the shared cross-section basis (POD, 241 modes) misses the mode's cross-section shape.
Here the Ritz value is recomputed with denser half-wave lists (same end layers e16s15).  If it falls towards the
full GBT value, (a) is the cause; if it stays, (b) is.

    python scripts/b1_outlier_probe.py data/operators/periodic_mixed_SS_v2/test --sids 5cc7be690865 6cc2e9904ad4
Results: results/b1_targets/outlier_probe.jsonl (resumable)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from b0_long_cc_probe import ritz                              # noqa: E402
from pinnbuck import opio                                      # noqa: E402
from pinnbuck.dictionary import DICTIONARIES, END_LAYERS       # noqa: E402

VARIANTS = {"d2": DICTIONARIES["d2"],
            "hw21": tuple(np.geomspace(25, 1006, 21)),
            "hw31": tuple(np.geomspace(25, 1006, 31)),
            "hw21_from12": tuple(np.geomspace(12, 1006, 23))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--sids", nargs="+", required=True)
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    a = ap.parse_args()
    torch.set_num_threads(1)
    out = ROOT / "results" / "b1_targets" / "outlier_probe.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {(j["sid"], j["variant"]) for j in map(json.loads, open(out))} if out.exists() else set()
    for sid in a.sids:
        op = opio.load(Path(a.src) / f"{sid}.pkl")
        for v in a.variants:
            if (sid, v) in done:
                continue
            N, it, n, dt = ritz(op, VARIANTS[v], 400, **END_LAYERS["d2"])
            r = dict(sid=sid, variant=v, half_waves=list(map(float, VARIANTS[v])), N=N, full=op.lam_full,
                     err=100 * (N / op.lam_full - 1), iters=it, n=n, t=dt, bs=op.params["bs"], L=op.L)
            print(f"{sid} {v:12s} N {N:9.4f}  full {op.lam_full:9.4f}  {r['err']:+.3f} %  n {n}  {dt:.0f} s",
                  flush=True)
            with open(out, "a") as f:
                f.write(json.dumps(r) + "\n")
        del op


if __name__ == "__main__":
    main()
