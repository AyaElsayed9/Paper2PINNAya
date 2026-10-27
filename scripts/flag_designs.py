"""Supervisor of the active learning (Paper 1, section 2.5, adapted): flag the training
designs where the PINN is least trustworthy, so that the oracle (full GBT) is called there.

    python scripts/flag_designs.py runs/p1_ss_dec.pt --train data/operators/periodic_industrial_SS/train --frac 0.1

Physics indicator (no labels needed): the relative gap between the Rayleigh-Ritz load of the
predicted modes (an upper bound) and the head,  s = (N_Ritz - N_head) / N_head.  A large gap
means the predicted modes are far from an eigenmode (the energy is not at a minimum) or the
head disagrees with its own physics, the analogue of the suspicion check and the ensemble
CoV of Paper 1.  Writes <model>.flags.json with the scores and the selected designs.
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--train", required=True)
    ap.add_argument("--frac", type=float, default=0.1)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    ck = torch.load(a.model, weights_only=False); ar = ck["args"]
    files = [f for f in sorted(Path(a.train).glob("*.pkl")) if pickle.load(open(f, "rb")).params["bc"] == ar["bc"]]
    net, rows = None, []
    for f in files:
        d = epinn.prepare(pickle.load(open(f, "rb")))
        if net is None:
            net = epinn.build_net(ar.get("arch", "linear"), d.x[None].repeat(2, 0), d.rs, R=ar["R"], k=ar["k"],
                                  width=ar["width"], depth=ar["depth"])
            net.load_state_dict(ck["state"]); net.eval()
        with torch.no_grad():
            Y, logN = net(torch.as_tensor(d.x[None]))
            mu, _ = epinn.ritz(d, Y[0])
        nr, nh = 1.0 / float(mu[0]), float(torch.exp(logN[0, 0]))
        rows.append(dict(file=str(f), sid=d.sid, N_ritz=nr, N_head=nh, score=(nr - nh) / nh))
    rows.sort(key=lambda r: -r["score"])
    nsel = max(1, int(round(a.frac * len(rows))))
    sel = [r["file"] for r in rows[:nsel]]
    out = Path(str(a.model).replace(".pt", ".flags.json"))
    json.dump(dict(rows=rows, selected=sel), open(out, "w"), indent=1)
    s = np.array([r["score"] for r in rows])
    print(f"{len(rows)} designs; gap Ritz vs head: median {100 * np.median(s):.1f} %, 90 % {100 * np.quantile(s, 0.9):.1f} %, "
          f"max {100 * s.max():.1f} %; {nsel} flagged (gap >= {100 * rows[nsel - 1]['score']:.1f} %) -> {out}")


if __name__ == "__main__":
    main()
