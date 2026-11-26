"""New PINN (cross-ply / symmetry expert) and its GBT check at the 22 Abaqus points of Paper 1.

Aya (26 Nov 2026): compare like with like in the GAMM abstract: Paper 1 PINN vs Abaqus (2.53 % mean) against the
new PINN vs Abaqus, instant and after the GBT check.  Abaqus and Paper-1 PINN values are copied from Paper 1, Table 2
(Elsayed et al., Compos. Struct. 399 (2026) 120959), exactly as in scripts/gbt_at_paper1_points.py; nothing tuned.

Per point: GBT operators of the design (cross-ply of Paper 1, symmetry side, SS, design space v2, dictionary d2, full
GBT on the same operator), then for the experts symmetry_crossply_SS_v2 physics_L0 (0 labels) and hybrid_L120
(120 labels): instant N (head) and the certified check (block LOBPCG, energy metric, R = 8, as in all checks),
read at eta <= 5 % (the tool's stopping rule; traced up to 40 steps as in the cost study) and after 12 steps.  Note: point (124.25, 0.72) lies just below the
expert range beta >= 0.75 (extrapolation).

    python scripts/pinn_at_paper1_points.py            -> results_p1/pinn_p1.jsonl (resumable)
"""
import json
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn                                            # noqa: E402
from pinnbuck.designspace import FAMILIES, LAYUP_CATALOG, Sample      # noqa: E402
from pinnbuck.operators import build_operator                         # noqa: E402
import cost_study                                                     # noqa: E402

import ast                                                            # noqa: E402
_src = open(ROOT / "scripts/gbt_at_paper1_points.py").read()             # same Paper 1 Table 2 list, not re-typed
P1 = ast.literal_eval(_src[_src.index("P1 = [") + 5:_src.index("]\nout")+1])
TAG = "symmetry_crossply_SS_v2"
OUT = ROOT / "results_p1" / "pinn_p1.jsonl"


def main():
    torch.set_num_threads(int(os.environ.get("THREADS", "2")))
    cost_study.TOL_STOP, cost_study.MAX_STEPS = 0.0, 40            # trace to 40 steps (as the cost study); eta per step
    done = {(r["bs"], r["beta"]) for r in map(json.loads, open(OUT))} if OUT.exists() else set()
    sb = pickle.load(open(ROOT / "data/shared_basis/basis_symmetry.pkl", "rb"))
    cks = {m: torch.load(ROOT / f"runs/experts/{TAG}_{m}.pt", weights_only=False) for m in ("physics_L0", "hybrid_L120")}
    epinn.INPUT_COLS = [0, 1]
    sk, st = LAYUP_CATALOG["crossply"]
    nets = {}
    for reg, bs, beta, n_fea, n_p1, m_fea in P1:
        if (bs, beta) in done:
            continue
        s = Sample("omega", dict(bs=float(bs), beta=float(beta), bc="SS", side="symmetry"),
                   {"skin": [float(a) for a in sk], "stringer": [float(a) for a in st]},
                   dict(FAMILIES["omega"]["fixed"]), "v2")
        t = time.perf_counter()
        op = build_operator(s, sb, lam_full=None)
        op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32)
        op.info["dictionary"] = "d2"
        t_op = time.perf_counter() - t
        d = epinn.prepare(op)
        rec = dict(region=reg, bs=bs, beta=beta, N_FEA=n_fea, N_PINN_P1=n_p1, N_GBT=float(d.lam_full), t_op=round(t_op, 1))
        for m, ck in cks.items():
            if m not in nets:
                nets[m] = epinn.net_from_ckpt(ck, d)
            with torch.no_grad():
                Y, logN = nets[m](torch.as_tensor(d.x[None]))
                Y = Y[0][:, d.sel].double()
            Ns, etas, ts = cost_study.trace(d.rs, d.rs.unflat(d.fm.color(Y)), d.em.solve)
            k5 = next((i for i, v in enumerate(etas) if v <= 0.05), None)
            rec[m] = dict(head=float(torch.exp(logN[0, 0])), N_eta5=Ns[k5] if k5 is not None else None, steps_eta5=k5,
                          N12=Ns[12], eta12=etas[12], N=[round(v, 6) for v in Ns], eta=[round(v, 6) for v in etas])
        OUT.parent.mkdir(exist_ok=True)
        open(OUT, "a").write(json.dumps(rec) + "\n")
        e = lambda v: float("nan") if v is None else 100 * (v / n_fea - 1)
        print(f"{reg:10s} bs {bs:7.2f} beta {beta:4.2f}  Abaqus {n_fea:7.3f} | GBT {e(rec['N_GBT']):+5.2f} % | "
              + " | ".join(f"{m}: head {e(r['head']):+6.2f} %, check(eta<=5%) {e(r['N_eta5']):+5.2f} % "
                           f"in {r['steps_eta5']} st, 12 st {e(r['N12']):+5.2f} %" for m, r in rec.items() if isinstance(r, dict)),
              flush=True)
        del d, op


if __name__ == "__main__":
    main()
