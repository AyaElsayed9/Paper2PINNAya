"""GBT teacher modes in the network coordinates (shared modes x Fourier in z).

For every design the full GBT (the validated model) gives the three lowest buckling
modes.  They are
  1. projected onto the network representation: shared GBT modes (`basis.py`) times the
     longitudinal functions b_c(zeta) {cos j pi zeta, sin j pi zeta} of every end-condition
     class (`spectral.FourierPanel`), by least squares on values, slopes and curvatures;
  2. refined by a few block-LOBPCG iterations on the exact GBT energies in that same
     representation, so the targets are exact eigenvectors of the network's own
     discretisation (projection alone leaves up to ~1 % in N_cr on end-localised modes).

A mode is stored as M (nm x nf): q(z) = sum_c Pi_c Pi_c^T M F_c(z)^T, the same layout for
every end condition (the classes only change the functions F_c).
"""
from __future__ import annotations

import numpy as np
import torch

from .spectral import FourierPanel, WavePrecond, cat, dot, lincomb, lobpcg, rayleigh_ritz

L_SCALE = 20.0          # [mm] weight of slopes and curvatures in the projection


def project(fp: FourierPanel, T: np.ndarray):
    """Teacher fields T (k, 3, nz, nm) -> coefficient list per class (k, k_c, nf)."""
    w = np.sqrt(fp.op.wz)[:, None]
    out = []
    for F, Pi in zip(fp.F, fp.Pi):
        Fn, Pn = F.numpy(), Pi.numpy()
        A = np.concatenate([Fn[d] * w * L_SCALE ** d for d in range(3)], 0)
        X = []
        for k in range(T.shape[0]):
            Y = np.concatenate([(T[k, d] @ Pn) * w * L_SCALE ** d for d in range(3)], 0)
            X.append(np.linalg.lstsq(A, Y, rcond=None)[0].T)          # (k_c, nf)
        out.append(torch.as_tensor(np.stack(X), dtype=fp.dtype))
    return out


def to_modes(fp: FourierPanel, A):
    """Class coefficients -> M (k, nm, nf)."""
    return sum(torch.einsum("mc,kcf->kmf", Pi, Ac) for Pi, Ac in zip(fp.Pi, A))


def from_modes(fp: FourierPanel, M):
    """M (k, nm, nf) -> class coefficients."""
    return [torch.einsum("mc,kmf->kcf", Pi, M) for Pi in fp.Pi]


def rayleigh(fp: FourierPanel, A):
    """Buckling factors of the individual modes (no Ritz combination)."""
    KA, GA = fp.apply_K(A), fp.apply_G(A)
    return (dot(A, KA).diagonal() / dot(A, GA).diagonal())


def refine(fp: FourierPanel, A0, iters=12, extra=3, seed=0):
    """Block LOBPCG started from the projected modes (+ `extra` random vectors).
    Returns the buckling factors (k,) and the refined modes (k, nm, nf), K-normalised."""
    k = A0[0].shape[0]
    T = WavePrecond(fp)
    gen = torch.Generator().manual_seed(seed)
    X0 = cat(A0, T(fp.random(extra, gen)))
    lam, X, hist = lobpcg(fp, R=k + extra, iters=iters, X0=X0, precond=T)
    Xk = [x[:k] for x in X]
    return lam[:k], to_modes(fp, Xk), hist
