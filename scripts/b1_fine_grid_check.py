"""Step B1 check: are the end-layer gains real or an integration artefact?

On the solver's own grid (span hx, 4 Gauss points per span) the end-layer probe gave Ritz values up to 0.15 %
below the full GBT label for short panels, where the end splines are finer than hx.  Here the operator of a
design is rebuilt on a grid `factor` times finer (only hx changes: finer pre-buckling field, more Gauss points;
the frozen solver is not touched) and the same Ritz spaces (d2 waves, end layers e8 / e12 / e16, 400
iterations) are compared with the full GBT on that grid.  On a fine enough grid the Ritz value must lie above
the full GBT value (Rayleigh-Ritz), so a negative difference on the solver grid that turns positive here is
an artefact of the coarse integration.

    python scripts/b1_fine_grid_check.py data/operators/periodic_mixed_SS_v2/test --sids 040508c6cbc1 --factor 0.5
Results: results/b1_end_layer_probe/fine_grid.jsonl (resumable).  Needs about 3 GB per design for long panels.
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

import pinnbuck.families as fam                                # noqa: E402
from b0_long_cc_probe import ritz                              # noqa: E402
from b1_end_layer_probe import VARIANTS                        # noqa: E402
from pinnbuck.designspace import Sample                        # noqa: E402
from pinnbuck.dictionary import DICTIONARIES                   # noqa: E402
from pinnbuck.operators import build_operator                  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--sids", nargs="+", required=True)
    ap.add_argument("--factor", type=float, default=0.5)
    a = ap.parse_args()
    torch.set_num_threads(1)
    out = ROOT / "results" / "b1_end_layer_probe" / "fine_grid.jsonl"
    done = {(j["sid"], j["factor"]) for j in map(json.loads, open(out))} if out.exists() else set()
    orig = fam.hx_of
    fam.hx_of = lambda sec, L: orig(sec, L) * a.factor
    for sid in a.sids:
        if (sid, a.factor) in done:
            continue
        op0 = pickle.load(open(Path(a.src) / f"{sid}.pkl", "rb"))
        d = dict(op0.info["sample"]); d.pop("sid"); smp = Sample(**d)
        assert smp.sid == sid
        side = smp.params["side"]
        label = op0.lam_full
        del op0
        sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
        op = build_operator(smp, sb, lam_full=None, with_shared_reference=True)
        r = dict(sid=sid, factor=a.factor, bs=smp.params["bs"], L=op.L, label=label, full=op.lam_full,
                 shared=op.lam_shared)
        for tag, kw in VARIANTS:
            r[tag] = ritz(op, DICTIONARIES["d2"], 400, **kw)[0]
        r["vs_own_full"] = {k: 100 * (r[k] / r["full"] - 1) for k in ("shared",) + tuple(t for t, _ in VARIANTS)}
        r["vs_label"] = {k: 100 * (r[k] / label - 1) for k in ("full", "shared") + tuple(t for t, _ in VARIANTS)}
        print(sid, "vs own full %", {k: round(v, 4) for k, v in r["vs_own_full"].items()},
              "| full vs label %", round(r["vs_label"]["full"], 4), flush=True)
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
        del op, sb


if __name__ == "__main__":
    main()
