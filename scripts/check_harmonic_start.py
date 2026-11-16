"""Step 1 (16 Nov 2026): starts for the physics check built in LENGTH-HARMONIC coordinates.

Finding of scripts/check_mode_structure.py: an exact GBT mode is ~98 % ONE harmonic sin(j pi z / L) times a
cross-section vector, while on the dictionary's fixed wavelengths it is spread over 6-10 groups.  So a mode is
carried from a labelled design to a new one much better as "cross-section vectors per harmonic" than as dictionary
coefficients.  Starts compared here (same check, same 40 test designs, same seed as check_start_baselines.py):

  nn<NL>            : old copy of the dictionary coefficients of the nearest labelled design (reference, recomputed)
  hs<NL> (stretch)  : the neighbour's modes stretched to the new length (same harmonic numbers j)
  hw<NL> (wave)     : the neighbour's modes with the same half-wavelengths (j' = round(j L / L_nb))
  hsw<NL>           : both together (16 vectors; the first Rayleigh-Ritz step keeps the best 8)

No network yet: this measures what a network that predicts cross-section vectors per harmonic could reach.
    python scripts/check_harmonic_start.py symmetry_crossply_SS_v2
-> results/check_baselines/<tag>_harmonic.csv (resumable)"""
import argparse, csv, sys, time
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio
from pinnbuck.ritzspace import EnergyMetric
from pinnbuck.spectral import lobpcg
from check_start_baselines import scaled


def harmonics(z, L, J):
    """(nz, 2J + 1): sin(j pi z / L), j = 1..J, then cos(j pi z / L), j = 0..J."""
    j = np.arange(1, J + 1)[None, :]; zz = np.asarray(z)[:, None]
    return np.concatenate([np.sin(j * np.pi * zz / L), np.ones_like(zz), np.cos(j * np.pi * zz / L)], 1)


def decompose(d, T, J=None):
    """Harmonic cross-section vectors of whitened modes T (R, n) of design d: (R, 2J + 1, nm)."""
    op = d.rs.op; z, wz, L = np.asarray(op.z), np.asarray(op.wz), op.L
    J = J or int(min(80, max(8, np.ceil(L / 15.0))))
    Q = d.rs.field(d.rs.unflat(d.fm.color(T.double())))[:, 0].numpy()            # (R, nz, nm)
    B = harmonics(z, L, J) * np.sqrt(wz)[:, None]
    C = np.stack([np.linalg.lstsq(B, q * np.sqrt(wz)[:, None], rcond=None)[0] for q in Q])
    return C, J, L


def compose(d, C, J, jmap):
    """Whitened start vectors (R, n) of design d from harmonic vectors C (R, 2J + 1, nm) with harmonic numbers
    mapped by jmap (array of the new j for j = 1..J); projected on the design's dictionary by least squares."""
    op = d.rs.op; z, wz, L = np.asarray(op.z), np.asarray(op.wz), op.L
    jn = np.asarray(jmap); zz = z[:, None]
    H = np.concatenate([np.sin(jn[None, :] * np.pi * zz / L), np.ones_like(zz), np.cos(jn[None, :] * np.pi * zz / L)], 1)
    Q = np.einsum("zh,rhm->rzm", H, C)                                             # (R, nz, nm)
    sw = np.sqrt(wz)[:, None]; A = []
    for Bc, Pi in zip(d.rs.B, d.rs.Pi):
        B0 = Bc[0].numpy() * sw; P = Pi.numpy()                                    # (nz, nb), (nm, k)
        rhs = np.einsum("rzm,mk->rzk", Q, P) * sw[None]                            # (R, nz, k)
        X = np.linalg.lstsq(B0, rhs.transpose(1, 0, 2).reshape(len(z), -1), rcond=1e-10)[0]     # (nb, R k)
        A.append(torch.as_tensor(X.reshape(B0.shape[1], Q.shape[0], P.shape[1]).transpose(1, 2, 0).copy()))
    Y = d.fm.whiten(d.rs.flat(A))
    return Y / Y.norm(dim=1, keepdim=True).clamp(min=1e-300)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--nn", default="10,120")
    ap.add_argument("--steps", default="0,1,2,3,5,8,12"); ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=2); ap.add_argument("--inputs", default="geom")
    ap.add_argument("--variants", default="nn,hs", help="nn, hs, hw, hsw (hw: 3-design smoke test 16 Nov was far worse, "
                                                         "rounding the harmonic numbers destroys the modes)")
    a = ap.parse_args(); torch.set_num_threads(a.threads)
    epinn.INPUT_COLS = None if a.inputs == "all" else [0, 1]
    steps = [int(s) for s in a.steps.split(",")]; S = max(steps)
    fcsv = ROOT / f"results/check_baselines/{a.tag}_harmonic.csv"
    done = {(r["sid"], r["start"]) for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab]); cache = {}

    def neighbour(d, NL):
        x = scaled(d.params["bs"], d.params["beta"]); sid = tab[int(np.argmin(((Xtr[:NL] - x) ** 2).sum(1)))]["sid"]
        if sid not in cache:
            dn = epinn.prepare(opio.load(ROOT / f"data/operators/{a.tag}/train/{sid}.pkl"))
            T = torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{sid}.pt", weights_only=False)["T"]
            C, J, Ln = decompose(dn, T); cache[sid] = (dn.sel, T, C, J, Ln); del dn
        return sid, cache[sid]

    files = sorted((ROOT / f"data/operators/{a.tag}/test").glob("*.pkl"))[:a.n]
    with open(fcsv, "a", newline="") as fh:
        w = csv.writer(fh)
        if not done:
            w.writerow(["sid", "bs", "beta", "N_full", "start", "start_from", "n_vec", "t_start_s", "t_check_s"] + [f"s{s}" for s in range(0, S + 1)])
        for f in files:
            d = epinn.prepare(opio.load(f)); em = d.em
            for NL in [int(v) for v in a.nn.split(",")]:
                sid, (sel, T, C, J, Ln) = neighbour(d, NL)
                full = torch.zeros(T.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=T.dtype); full[:, sel] = T
                t0 = time.perf_counter()
                j0 = np.arange(1, J + 1)
                Ys = compose(d, C, J, j0); ts = time.perf_counter() - t0
                Yw = compose(d, C, J, np.maximum(1, np.rint(j0 * d.rs.op.L / Ln)))
                cand = [(f"nn{NL}", full[:, d.sel], 0.0), (f"hs{NL}", Ys, ts), (f"hw{NL}", Yw, ts),
                        (f"hsw{NL}", torch.cat([Ys, Yw]), 2 * ts)]
                for name, Y, tt in [c for c in cand if c[0].rstrip("0123456789") in a.variants.split(",")]:
                    if (d.sid, name) in done:
                        continue
                    mu0, _ = epinn.ritz(d, Y)
                    X0 = d.rs.unflat(d.fm.color(Y.double())); t1 = time.perf_counter()
                    _, _, h = lobpcg(d.rs, R=8, iters=S, X0=X0, precond=em.solve, tol=0, seed=1)
                    e = [100 * (v / d.lam_full - 1) for v in [1.0 / float(mu0[0])] + list(h)]
                    w.writerow([d.sid, d.params["bs"], d.params["beta"], d.lam_full, name, sid, Y.shape[0], f"{tt:.2f}",
                                f"{time.perf_counter() - t1:.1f}"] + [f"{v:.4f}" for v in e]); fh.flush()
                    print(f"{d.sid} {name:8s} " + "  ".join(f"s{s} {e[s]:+.2f}" for s in steps), flush=True)
            del d, em
    rows = list(csv.DictReader(open(fcsv)))
    print(f"\n{a.tag}: |error| vs full GBT after s steps (s0 = the start itself), mean / max over the test designs")
    for name in sorted({r["start"] for r in rows}):
        rr = [r for r in rows if r["start"] == name]
        print(f"  {name:8s} ({len(rr)} designs) " + "  ".join(
            f"s{s} {np.mean([abs(float(r[f's{s}'])) for r in rr]):.2f}/{np.max([abs(float(r[f's{s}'])) for r in rr]):.2f}" for s in steps))


if __name__ == "__main__":
    main()
