"""Start baselines for the level-2 network (one PINN, 38 inputs, all industrial lay-ups; 9 Nov 2026 evening).

Same test as scripts/check_start_baselines.py, but the nearest labelled design is searched in the network's own
standardised 38 inputs (mean / std of the training designs), among the first NL labelled training designs.
The question: with many inputs, is copying the nearest labelled design still better than the network?
results/check_baselines/level2_<model>.csv
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio                 # noqa: E402
from pinnbuck.ritzspace import EnergyMetric      # noqa: E402
from pinnbuck.spectral import lobpcg             # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="runs/level2/industrial_hybrid_L518.pt")
    ap.add_argument("--files", default="runs/level2/files_industrial.txt")
    ap.add_argument("--n-train", type=int, default=518)
    ap.add_argument("--labels", default="data/exact_targets/periodic_mixed_SS_v2/train")
    ap.add_argument("--test", default="data/operators/periodic_industrial_SS_v2/test")
    ap.add_argument("--nn", default="50,518")
    ap.add_argument("--steps", default="3,5,8,12,20,30")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    steps = [int(s) for s in a.steps.split(",")]; S = max(steps)
    ck = torch.load(ROOT / a.model, weights_only=False); ar = ck["args"]
    epinn.INPUT_COLS = None if ar.get("inputs", "all") == "all" else [0, 1]
    files = [ROOT / l.strip() for l in open(ROOT / a.files) if l.strip()][:a.n_train]
    t0 = time.perf_counter()
    X, SEL, SID = [], [], []
    for f in files:                                      # inputs and layout of every labelled training design
        op = opio.load(f); X.append(epinn.design_inputs(op)); SEL.append(epinn.design_sel(op)); SID.append(op.sid); del op
    X = np.stack(X).astype(float); m, s = X.mean(0), X.std(0); s[s < 1e-6] = 1.0
    Z = (X - m) / s
    print(f"{len(files)} training designs read in {time.perf_counter() - t0:.0f} s", flush=True)

    def nn_start(d, NL, rank=0):
        z = (np.asarray(d.x, float) - m) / s
        j = int(np.argsort(((Z[:NL] - z) ** 2).sum(1))[rank])
        T = torch.load(ROOT / a.labels / f"{SID[j]}.pt", weights_only=False)["T"]
        full = torch.zeros(T.shape[0], int(max(SEL[j].max(), d.sel.max())) + 1, dtype=T.dtype)
        full[:, SEL[j]] = T
        return full[:, d.sel], SID[j], float(np.sqrt(((Z[j] - z) ** 2).sum()))

    out = ROOT / "results/check_baselines"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / f"level2_{Path(a.model).stem}.csv"
    done = {(r["sid"], r["start"]) for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    net = None
    with open(fcsv, "a", newline="") as fh:
        w = csv.writer(fh)
        if not done:
            w.writerow(["sid", "N_full", "start", "start_from", "distance", "time_s"] + [f"s{k}" for k in range(1, S + 1)])
        for f in sorted((ROOT / a.test).glob("*.pkl")):
            d = epinn.prepare(opio.load(f))
            if net is None:
                net = epinn.net_from_ckpt(ck, d)
            em = EnergyMetric(d.rs, groups="labels", ridge=1e-8)
            starts = [("random", None, "-", np.nan)]
            for NL in [int(v) for v in a.nn.split(",")]:
                Y, sid, dist = nn_start(d, NL); starts.append((f"nn{NL}", Y, sid, dist))
            starts.append(("knn2_518", torch.cat([nn_start(d, a.n_train, k)[0] for k in range(2)]), "2 nearest", np.nan))
            with torch.no_grad():
                Y, _ = net(torch.as_tensor(d.x[None]))
            starts.append(("net", Y[0, :, d.sel], Path(a.model).stem, np.nan))
            for name, Y, src, dist in starts:
                if (d.sid, name) in done:
                    continue
                X0 = None if Y is None else d.rs.unflat(d.fm.color(Y.double()))
                t1 = time.perf_counter()
                _, _, h = lobpcg(d.rs, R=8, iters=S, X0=X0, precond=em.solve, tol=0, seed=1)
                dt = time.perf_counter() - t1
                e = [100 * (v / d.lam_full - 1) for v in h]
                w.writerow([d.sid, d.lam_full, name, src, f"{dist:.3f}", f"{dt:.1f}"] + [f"{v:.4f}" for v in e]); fh.flush()
                print(f"{d.sid} {name:10s} " + "  ".join(f"s{k} {e[k - 1]:+.2f}" for k in steps) + f"  ({dt:.0f} s)", flush=True)
            del d, em
    rows = list(csv.DictReader(open(fcsv)))
    print("\nlevel 2: |error| vs full GBT after s steps, mean / max over the test designs")
    for name in sorted({r["start"] for r in rows}):
        rr = [r for r in rows if r["start"] == name]
        print(f"  {name:10s} ({len(rr)} designs) " + "  ".join(
            f"s{k} {np.mean([abs(float(r[f's{k}'])) for r in rr]):.2f}/{np.max([abs(float(r[f's{k}'])) for r in rr]):.2f}"
            for k in steps))


if __name__ == "__main__":
    main()
