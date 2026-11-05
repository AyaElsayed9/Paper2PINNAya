"""Step B1 check, cause (b) of the outliers: does the shared cross-section basis (POD, 241 modes) hold the mode?

Rebuilds the operator of a design with `with_shared_reference=True`: the frozen solver is run once more with
the cross-section restricted to the shared basis but with its own full spline space along the length.  This
value (lam_shared) is the best any longitudinal dictionary can reach; if it is close to the full GBT label, the
shared basis is fine and the gap is longitudinal (dictionary or iterations).

    python scripts/b1_shared_basis_check.py data/operators/periodic_mixed_SS_v2/test --sids 5cc7be690865
Results: results/b1_targets/shared_basis_check.jsonl (resumable).  Needs about 3 GB per design.
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import opio                                      # noqa: E402
from pinnbuck.designspace import Sample                        # noqa: E402
from pinnbuck.operators import build_operator                  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("--sids", nargs="+", required=True)
    a = ap.parse_args()
    out = ROOT / "results" / "b1_targets" / "shared_basis_check.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    done = {json.loads(l)["sid"] for l in open(out)} if out.exists() else set()
    for sid in a.sids:
        if sid in done:
            continue
        op0 = opio.load(Path(a.src) / f"{sid}.pkl")
        d = dict(op0.info["sample"]); d.pop("sid"); smp = Sample(**d)
        assert smp.sid == sid
        label = op0.lam_full
        del op0
        sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{smp.params['side']}.pkl", "rb"))
        op = build_operator(smp, sb, lam_full=label, with_shared_reference=True)
        r = dict(sid=sid, label=label, shared=op.lam_shared, shared_err=100 * (op.lam_shared / label - 1),
                 bs=smp.params["bs"], L=op.L)
        print(f"{sid}  label {label:.4f}  shared basis {op.lam_shared:.4f}  ({r['shared_err']:+.3f} %)", flush=True)
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
        del op, sb


if __name__ == "__main__":
    main()
