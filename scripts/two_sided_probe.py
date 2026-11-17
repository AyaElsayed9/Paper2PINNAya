"""17 Nov 2026 (plan: a lower bound beside the check's upper bound). Kato-Temple estimate after s check steps.
Pencil G x = mu K x, N = 1/mu, critical = largest mu.  For the check's mode-1 Ritz pair (theta, x), x^T K x = 1:
residual r = G x - theta K x, eps^2 = r^T K^-1 r (K^-1 by preconditioned CG, the check's metric as preconditioner).
Temple: if nu >= mu_2 and nu < theta, then mu_1 <= theta + eps^2 / (theta - nu)  ->  N_1 >= 1 / (theta + eps^2/(theta - nu)).
nu = the check's own Ritz value of mode 2 (an ESTIMATE of mu_2 from below, so the lower value is an estimate, not a
proof).  Compared with the converged Ritz-space load (the quantity the check bounds) and with full GBT.
    python scripts/two_sided_probe.py <tag> --n 10 --steps 5,12   -> results/two_sided/<tag>.csv"""
import argparse, csv, sys, time
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio
from pinnbuck.spectral import lobpcg


def pcg(rs, em, b, tol=1e-10, maxit=400):
    """K y = b in list format, preconditioned with em.solve; returns y and iterations"""
    flat, unflat = rs.flat, rs.unflat
    bf = flat(b); x = torch.zeros_like(bf); r = bf.clone(); z = flat(em.solve(unflat(r))); p = z.clone(); rz = (r * z).sum()
    nb = float(bf.norm())
    for it in range(maxit):
        Ap = flat(rs.apply_K(unflat(p))); alpha = rz / (p * Ap).sum()
        x = x + alpha * p; r = r - alpha * Ap
        if float(r.norm()) < tol * nb:
            return unflat(x), it + 1
        z = flat(em.solve(unflat(r))); rz_new = (r * z).sum(); p = z + (rz_new / rz) * p; rz = rz_new
    return unflat(x), maxit


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--steps", default="5,12"); a = ap.parse_args(); torch.set_num_threads(2)
    ck = torch.load(ROOT / f"runs/experts/{a.tag}_hybrid_L120.pt", weights_only=False)
    epinn.INPUT_COLS = None if ck["args"].get("inputs", "all") == "all" else [0, 1]
    out = ROOT / "results/two_sided"; out.mkdir(parents=True, exist_ok=True); fcsv = out / f"{a.tag}.csv"
    rows = []; net = None
    for f in sorted((ROOT / f"data/operators/{a.tag}/test").glob("*.pkl"))[:a.n]:
        op = opio.load(f); d = epinn.prepare(op)
        if net is None:
            net = epinn.net_from_ckpt(ck, d)
        with torch.no_grad():
            Y = net(torch.as_tensor(d.x[None]))[0][0, :, d.sel].double()
        Nconv, _, hc = lobpcg(d.rs, R=8, iters=150, X0=d.rs.unflat(d.fm.color(Y)), precond=d.em.solve, tol=1e-12, seed=1)
        Nconv = float(Nconv[0])
        for s in [int(v) for v in a.steps.split(",")]:
            N, X, h = lobpcg(d.rs, R=8, iters=s, X0=d.rs.unflat(d.fm.color(Y)), precond=d.em.solve, tol=0, seed=1)
            x = [t[:1] for t in X]
            Kx, Gx = d.rs.apply_K(x), d.rs.apply_G(x)
            kk = float(d.rs.flat(x) @ d.rs.flat(Kx).T); gg = float(d.rs.flat(x) @ d.rs.flat(Gx).T)
            sc = 1 / np.sqrt(kk); x = [t * sc for t in x]; Kx = [t * sc for t in Kx]; Gx = [t * sc for t in Gx]
            theta = gg / kk
            r = [g - theta * k for g, k in zip(Gx, Kx)]
            t0 = time.perf_counter(); y, it = pcg(d.rs, d.em, r); t_cg = time.perf_counter() - t0
            eps2 = float(d.rs.flat(r) @ d.rs.flat(y).T)
            nu = 1.0 / float(N[1])
            low = 1.0 / (theta + eps2 / (theta - nu)) if theta > nu else float("nan")
            row = dict(sid=op.sid, steps=s, N_up=1 / theta, N_low_est=low, N_ritz_conv=Nconv, N_full=float(op.lam_full),
                       up_vs_conv=100 * (1 / theta / Nconv - 1), low_vs_conv=100 * (low / Nconv - 1),
                       width_pct=100 * (1 / theta - low) / Nconv, gap12_pct=100 * (float(N[1]) * theta - 1),
                       bracket_ok=bool(low <= Nconv * (1 + 1e-12) <= 1 / theta * (1 + 1e-12)), cg_it=it, t_cg=round(t_cg, 2))
            rows.append(row)
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        del d, op
    with open(fcsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for s in sorted({r["steps"] for r in rows}):
        rr = [r for r in rows if r["steps"] == s]
        print(f"{a.tag}, {s} steps: upper bound {np.mean([r['up_vs_conv'] for r in rr]):+.3f} %, lower estimate "
              f"{np.mean([r['low_vs_conv'] for r in rr]):+.3f} % (min {min(r['low_vs_conv'] for r in rr):+.3f}), width mean "
              f"{np.mean([r['width_pct'] for r in rr]):.3f} % max {max(r['width_pct'] for r in rr):.3f} %, bracket holds "
              f"{sum(r['bracket_ok'] for r in rr)}/{len(rr)}, CG {np.mean([r['cg_it'] for r in rr]):.0f} it, {np.mean([r['t_cg'] for r in rr]):.2f} s")


if __name__ == "__main__":
    main()
