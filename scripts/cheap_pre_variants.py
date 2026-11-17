"""17 Nov 2026: cheaper pre-buckling states for the check, accuracy vs time.
For each design: the check's operators with different pre-buckling states, the Ritz-space critical load
converged (LOBPCG to tol 1e-10), compared with the same Ritz space using the solver's own fine state (x1).
  full_x1   : solver state, fine mesh (reference)          full_x4 : current check (coarse x4, cubic interp.)
  red_x1/x4 : the state computed in the SHARED reduced section w (241 modes instead of 570)
  uni_x4    : mid-length values of full_x4 at every z (a cross-section-like uniform state)
    python pre_variants.py <tag> <n>   -> /home/claude/paper2pinnaya/results/cheap_check/pre_variants_<tag>.csv"""
import csv, sys, time, pickle
from pathlib import Path
import numpy as np, torch
ROOT = Path("/home/claude/paper2pinnaya")
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn
from pinnbuck.designspace import sample_design
from pinnbuck.families import family_context
from pinnbuck.operators import PanelOperator, shared_section, section_operators
from pinnbuck.ritzspace import RitzSpace, EnergyMetric
from pinnbuck.dictionary import dictionary_kw
from pinnbuck.spectral import lobpcg
from scipy.interpolate import CubicSpline
torch.set_num_threads(2)
tag = sys.argv[1]; side, lay = tag.split("_")[:2]; n = int(sys.argv[2])
S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")[:n]
sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
out = Path(f"/home/claude/paper2pinnaya/results/cheap_check/pre_variants_{tag}.csv")
done = {r["sid"] for r in csv.DictReader(open(out))} if out.exists() else set()


def pre_on(m, sec_like, L, bc, hx, zf, coarsen):
    hc = coarsen * hx
    P = m.omega_gbt_solve(sec_like, L, bc, hx=hc, pre_only=True)["pre"]
    if coarsen == 1:
        return P
    zc = np.asarray(m._LongSpline(L, m._sfs_breaks(L, hc, hc, L)).xg)
    return tuple(CubicSpline(zc, np.asarray(p), axis=0)(zf) for p in P)


def converged(op):
    rs = RitzSpace(op, raw=True, **dictionary_kw(op)); em = EnergyMetric(rs, groups="labels", ridge=1e-8)
    N, _, h = lobpcg(rs, R=8, iters=120, X0=None, precond=em.solve, tol=1e-9, seed=1)
    return float(N[0]), len(h)


with open(out, "a", newline="") as fh:
    wr = None
    for s in S:
        if s.sid in done:
            continue
        r = dict(sid=s.sid, bs=round(s.params["bs"], 1), beta=round(s.params["beta"], 3))
        with family_context(s) as fc:
            m, sec, L, lb = fc.m, fc.sec, fc.L, fc.lb; zf = np.asarray(lb.xg)
            w = shared_section(m, sec, sb); cls = fc.classes(w); full = "T" in sec.families
            pres = {}
            for name, sl, c in (("full_x1", sec, 1), ("full_x4", sec, 4), ("full_x8", sec, 8), ("full_x16", sec, 16)):
                t = time.perf_counter(); pres[name] = pre_on(m, sl, L, s.params["bc"], fc.hx, zf, c); r[f"t_{name}"] = round(time.perf_counter() - t, 3)
            mid = len(zf) // 2
            pres["uni_x4"] = tuple(np.repeat(np.asarray(p)[mid:mid + 1], len(zf), axis=0) for p in pres["full_x4"])
            ops = {}
            for name, pre in pres.items():
                K, G, P = section_operators(w, pre, full)
                ops[name] = PanelOperator(sid=s.sid, params=dict(s.params), features=s.features(), L=float(L), bc=s.params["bc"],
                                          side=s.params.get("side"), nm=int(w.nm), K=K, G=G.astype(np.float32), P=P.astype(np.float32),
                                          z=zf, wz=np.asarray(lb.wx), classes=cls, lam_full=float("nan"), lam_shared=float("nan"),
                                          info=dict(hx=fc.hx, nspan=lb.nspan, sample=s.to_dict(), dictionary="d2"))
        for name, op in ops.items():
            N, it = converged(op); r[f"N_{name}"] = N; r[f"it_{name}"] = it
        for name in ops:
            r[f"e_{name}"] = round(100 * (r[f"N_{name}"] / r["N_full_x1"] - 1), 4)
        if wr is None:
            wr = csv.DictWriter(fh, fieldnames=list(r))
            if not done:
                wr.writeheader()
        wr.writerow(r); fh.flush()
        print(s.sid, {k: v for k, v in r.items() if k.startswith(("t_", "e_", "it_"))}, flush=True)
rows = list(csv.DictReader(open(out)))
print(f"\n{tag}: {len(rows)} designs")
for name in ("full_x4", "full_x8", "full_x16", "uni_x4"):
    e = np.array([float(x[f"e_{name}"]) for x in rows])
    tt = np.array([float(x.get(f"t_{name}", "nan")) for x in rows]) if name != "uni_x4" else np.array([float(x["t_full_x4"]) for x in rows])
    print(f"  {name:8s} error vs fine state: mean |e| {np.abs(e).mean():.3f} %, max {np.abs(e).max():.3f} % (range {e.min():+.3f} .. {e.max():+.3f})   pre time {tt.mean():.2f} s")
print(f"  full_x1  pre time {np.mean([float(x['t_full_x1']) for x in rows]):.2f} s")
