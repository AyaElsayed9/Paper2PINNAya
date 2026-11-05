"""Step B1 summary of the exact targets of a v2 set (no operator file is loaded, so it can run beside the jobs).

For every design with a target file: Ritz error of the network's shape space (N_ritz[0] vs the full-GBT label),
lay-up type (industrial / general), b_s, beta, L, the full-GBT time (from the build log) and the LOBPCG time.

    python scripts/b1_targets_summary.py --set test  --n 100  --seed 99
    python scripts/b1_targets_summary.py --set train --n 1064 --seed 11
Writes results/b1_targets/periodic_mixed_SS_v2_<set>.csv and prints the statistics.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck.designspace import sample_design                 # noqa: E402
from b1_end_layer_probe import is_industrial                   # noqa: E402


def build_times(log):
    t = {}
    if log.exists():
        for m in re.finditer(r"\] (\w{12}) ([\d.]+) s  N_full", log.read_text()):
            t[m.group(1)] = float(m.group(2))
    return t


def stats(x):
    x = np.asarray(x, float)
    if not len(x):
        return "n 0"
    return (f"n {len(x):4d}  mean {x.mean():.4f}  median {np.median(x):.4f}  p95 {np.percentile(x, 95):.4f}  "
            f"max {x.max():.4f}  min {x.min():.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="test")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=99)
    ap.add_argument("--side", default="periodic")
    a = ap.parse_args()
    tag = f"{a.side}_mixed_SS_v2"
    tdir = ROOT / "data" / "exact_targets" / tag / a.set
    bt = build_times(ROOT / "logs" / f"b1_build_{a.side}_SS_{a.set}.log")
    rows = []
    for s in sample_design("omega", a.n, seed=a.seed, categorical={"side": a.side, "bc": "SS"},
                           layups="mixed", space="v2"):
        f = tdir / f"{s.sid}.pt"
        if not f.exists():
            continue
        d = torch.load(f, weights_only=False)
        p = s.params
        rows.append(dict(sid=s.sid, industrial=int(is_industrial(s.layups)), bs=p["bs"],
                         beta=p.get("beta", float("nan")), N_full=d["N_full"], N_ritz=float(d["N_ritz"][0]),
                         ritz_err_pct=100 * (float(d["N_ritz"][0]) / d["N_full"] - 1), n_ritz=d["n"],
                         t_full_gbt=bt.get(s.sid, float("nan")), t_lobpcg=d["time"]))
    out = ROOT / "results" / "b1_targets" / f"{tag}_{a.set}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    e = np.array([r["ritz_err_pct"] for r in rows]); ind = np.array([r["industrial"] for r in rows], bool)
    print(f"{tag}/{a.set}: {len(rows)} of {a.n} designs have targets  -> {out.relative_to(ROOT)}")
    print("Ritz error % (all)       ", stats(e))
    print("Ritz error % (industrial)", stats(e[ind]))
    print("Ritz error % (general)   ", stats(e[~ind]))
    print("negative Ritz errors:", int((e < 0).sum()), " above 0.2 %:", int((e > 0.2).sum()),
          " above 0.5 %:", int((e > 0.5).sum()))
    bs = np.array([r["bs"] for r in rows])
    for lo, hi in ((0, 200), (200, 300), (300, 1e9)):
        m = (bs >= lo) & (bs < hi)
        print(f"  b_s {lo:>4}-{hi if hi < 1e9 else 'max':>4} mm:", stats(e[m]))
    tf = np.array([r["t_full_gbt"] for r in rows]); tl = np.array([r["t_lobpcg"] for r in rows])
    print("full GBT time s (build, 1 thread)", stats(tf[np.isfinite(tf)]))
    print("LOBPCG target time s             ", stats(tl))
    worst = sorted(rows, key=lambda r: -r["ritz_err_pct"])[:5]
    print("worst:", [(r["sid"], round(r["ritz_err_pct"], 3), round(r["bs"], 1), r["industrial"]) for r in worst])


if __name__ == "__main__":
    main()
