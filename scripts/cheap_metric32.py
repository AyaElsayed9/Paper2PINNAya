"""17 Nov 2026: energy metric (preconditioner + whitening) in float32 vs float64: time and check accuracy.
The Rayleigh-Ritz values stay in float64, so the upper-bound property is unchanged; only speed and the
quality of the preconditioner / start change.  nn120 and net starts, crossply / symmetry test designs."""
import csv, sys, time, pickle, math
from pathlib import Path
import numpy as np, torch
ROOT = Path("/home/claude/paper2pinnaya")
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio
from pinnbuck.designspace import sample_design
from pinnbuck.operators import build_operator
from pinnbuck.ritzspace import RitzSpace, EnergyMetric, FlatMetric
from pinnbuck.dictionary import dictionary_kw
from pinnbuck.spectral import lobpcg
from check_start_baselines import scaled
torch.set_num_threads(2)
tag = sys.argv[1]; n = int(sys.argv[2]); side, lay = tag.split("_")[:2]
S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")[:n]
sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
ck = torch.load(ROOT / f"runs/experts/{tag}_hybrid_L120.pt", weights_only=False)
epinn.INPUT_COLS = None if ck["args"].get("inputs", "all") == "all" else [0, 1]
tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{tag}_train.csv")))[:120]
Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab])
full_rows = {r["sid"]: r for r in csv.DictReader(open(ROOT / f"results/fast_check/{tag}.csv"))}

_orig_chol = torch.linalg.cholesky
class Metric32(EnergyMetric):
    def __init__(self, *a, **k):
        torch.linalg.cholesky = lambda M: _orig_chol(M.float())
        try:
            super().__init__(*a, **k)
        finally:
            torch.linalg.cholesky = _orig_chol
    def solve(self, Rv):
        out = [torch.zeros_like(r) for r in Rv]
        for b, off, Lc in self.bins:
            self._scatter(out, torch.cholesky_solve(self._gather(Rv, b, off).T.float(), Lc).T.double(), b, off)
        return out

net = None; res = []
for s in S:
    if s.sid not in full_rows:
        continue
    N_full = float(full_rows[s.sid]["N_full"])
    op = build_operator(s, sb, lam_full=float("nan"), with_shared_reference=False, pre_coarsen=4)
    op.G = op.G.astype(np.float32); op.P = op.P.astype(np.float32); op.info["dictionary"] = "d2"; op.lam_full = N_full
    d = epinn.prepare(op)                                           # float64 metric (as in fast_check)
    if net is None:
        net = epinn.net_from_ckpt(ck, d)
    j = int(np.argmin(((Xtr - scaled(s.params["bs"], s.params["beta"])) ** 2).sum(1)))
    sid = tab[j]["sid"]; o2 = opio.load(ROOT / f"data/operators/{tag}/train/{sid}.pkl")
    sel = epinn.design_sel(o2, epinn.layout_labels(o2)); del o2
    T = torch.load(ROOT / f"data/exact_targets/{tag}/train/{sid}.pt", weights_only=False)["T"]
    full = torch.zeros(T.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=T.dtype); full[:, sel] = T; Ynn = full[:, d.sel].double()
    with torch.no_grad():
        Ynet = net(torch.as_tensor(d.x[None]))[0][0, :, d.sel].double()
    r = dict(sid=s.sid)
    for prec in ("f64", "f32"):
        t = time.perf_counter()
        em = EnergyMetric(d.rs, groups="labels", ridge=1e-8) if prec == "f64" else Metric32(d.rs, groups="labels", ridge=1e-8)
        fm = FlatMetric(d.rs, em); r[f"t_metric_{prec}"] = time.perf_counter() - t
        for name, Y in (("nn120", Ynn), ("net", Ynet)):
            X0 = d.rs.unflat(fm.color(Y))
            stamps = [time.perf_counter()]
            _, _, h = lobpcg(d.rs, R=8, iters=12, X0=X0, precond=em.solve, tol=0, seed=1, callback=lambda it, mu: stamps.append(time.perf_counter()))
            r[f"t_step_{prec}_{name}"] = (stamps[12] - stamps[0]) / 12
            for k in (5, 12):
                r[f"e{k}_{prec}_{name}"] = 100 * (h[k - 1] / N_full - 1)
    res.append(r)
    print(s.sid, {k: round(v, 3) for k, v in r.items() if k != "sid"}, flush=True)
    del d, op
print(f"\n{tag}, {len(res)} designs (mean, |error| mean / max):")
for prec in ("f64", "f32"):
    tm = np.mean([r[f't_metric_{prec}'] for r in res]); ts = np.mean([r[f't_step_{prec}_nn120'] for r in res])
    line = f"  {prec}: metric {tm:.2f} s, step {ts:.3f} s"
    for name in ("nn120", "net"):
        for k in (5, 12):
            e = np.abs([r[f'e{k}_{prec}_{name}'] for r in res]); line += f" | {name} s{k} {e.mean():.2f}/{e.max():.2f}"
    print(line)
print("  below full GBT (any):", sum(r[k] < 0 for r in res for k in r if k.startswith("e")))
