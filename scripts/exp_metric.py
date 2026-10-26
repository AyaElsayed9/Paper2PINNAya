"""Experiment: quality of the energy metric (wavenumber-bin blocks of K) on the dictionary
Ritz space, and the eigen-iteration it allows.

    python scripts/exp_metric.py --folder data/operators/symmetry_industrial/test --max 2 [--dense]

Per design: size of the Ritz space; check of the Kronecker K and pointwise G against the
direct z-integrals; LOBPCG (block of 6) with the metric as preconditioner: iterations and
N_cr against the full GBT; --dense also gives the condition number of M^-1 K.
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.ritzspace import EnergyMetric, RitzSpace   # noqa: E402
from pinnbuck.spectral import dot, lobpcg                 # noqa: E402


def direct_energies(rs, A):
    """U and U_G of the vectors A from the z-integrals of the operator (no Kronecker)."""
    op = rs.op; nm = op.nm
    Q = rs.field(A)                                                    # (R, 3, nz, nm)
    Qc = Q.permute(0, 2, 1, 3).reshape(Q.shape[0], Q.shape[2], 3 * nm)
    K = torch.as_tensor(op.K)
    U = 0.5 * torch.einsum("z,rzi,ij,rzj->r", rs.wz, Qc, K, Qc)
    G = torch.as_tensor(op.G.astype(float))                            # (np, nc, 2, 3nm)
    g = torch.einsum("pcai,rzi->rzpca", G, Qc)
    P = torch.as_tensor(op.P.astype(float))                            # (3, nz, np)
    dens = (P[0][None, :, :, None] * g[..., 0] ** 2 + P[1][None, :, :, None] * g[..., 1] ** 2
            + 2 * P[2][None, :, :, None] * g[..., 0] * g[..., 1]).sum((-1, -2))
    UG = 0.5 * torch.einsum("z,rz->r", rs.wz, dens)
    return U, UG


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", default=str(ROOT / "data/operators/symmetry_industrial/test"))
    ap.add_argument("--max", type=int, default=2)
    ap.add_argument("--dense", action="store_true")
    ap.add_argument("--ratio", type=float, default=1.35)
    ap.add_argument("--raw", action="store_true", help="fixed dictionary coordinates, label groups")
    ap.add_argument("--ridge", type=float, default=1e-12)
    a = ap.parse_args()
    torch.set_num_threads(2)
    for f in sorted(Path(a.folder).glob("*.pkl"))[:a.max]:
        op = pickle.load(open(f, "rb"))
        t0 = time.perf_counter()
        rs = RitzSpace(op, raw=a.raw)
        t1 = time.perf_counter()
        gen = torch.Generator().manual_seed(0)
        A = rs.random(3, gen)
        U, UG = direct_energies(rs, A)
        Uk = 0.5 * torch.diagonal(dot(A, rs.apply_K(A))); Ug = 0.5 * torch.diagonal(dot(A, rs.apply_G(A)))
        em = EnergyMetric(rs, ratio=a.ratio, groups="labels" if a.raw else None, ridge=a.ridge)
        t2 = time.perf_counter()
        print(f"{f.stem} {op.params['bc']} bs {op.params['bs']:.0f} beta {op.params['beta']:.2f} L {op.L:.0f}: "
              f"n = {rs.n}, classes {rs.k}, functions {rs.nb}; space {t1 - t0:.1f} s, metric {t2 - t1:.1f} s "
              f"(blocks {len(em.sizes)}, largest {max(em.sizes)})")
        print(f"   check K: {float((Uk / U - 1).abs().max()):.1e}   G: {float((Ug / UG - 1).abs().max()):.1e}")
        for name, pre in (("metric", em.solve),):
            hist = []
            t = time.perf_counter()
            lam, X, h = lobpcg(rs, R=6, iters=200, precond=pre, tol=1e-9, seed=1)
            el = time.perf_counter() - t
            h = np.asarray(h); lam1 = float(lam[0])
            it01 = int(np.argmax(h / lam1 - 1 < 1e-3)) + 1
            print(f"   LOBPCG [{name}]: {len(h)} iterations ({el:.0f} s), within 0.1 % after {it01}; "
                  f"N = {lam1:.4f} vs full GBT {op.lam_full:.4f} ({100 * (lam1 / op.lam_full - 1):+.3f} %), "
                  f"next {', '.join(f'{float(x):.3f}' for x in lam[1:4])}")
        if a.dense:
            t = time.perf_counter()
            K = rs.dense_K()
            Lfull = torch.zeros_like(K)
            o = 0
            # assemble C^-1 action column-wise: M^-1/2 K M^-1/2 spectrum = eig of C^-T K C^-1
            I = torch.eye(rs.n, dtype=K.dtype)
            Ci = rs.flat(em.color(rs.unflat(I)))                       # rows: C^-1 e_i  -> (n, n) = (C^-1)^T
            S = Ci @ K @ Ci.T
            ev = torch.linalg.eigvalsh(0.5 * (S + S.T))
            evK = torch.linalg.eigvalsh(K)
            q = np.quantile(ev.numpy(), [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0])
            print(f"   cond(K) = {float(evK[-1] / max(float(evK[0]), 1e-300)):.2e},  cond(M^-1 K) = {float(ev[-1] / ev[0]):.2e}  "
                  f"(eig quantiles 0.1/1/10/50/90/99/100 %: {' '.join(f'{x:.3g}' for x in q)}; {time.perf_counter() - t:.0f} s)")


if __name__ == "__main__":
    main()
