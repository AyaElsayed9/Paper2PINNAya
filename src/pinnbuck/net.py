"""Fast buckling network of a family (omega first): design -> N_cr of the three lowest
modes and the modes themselves, taught by the validated GBT (Paper 1 workflow with the
GBT as oracle), with a physics check of every prediction.

Coordinates of a mode (the same for every design of the family):
  M (nm x nf): shared GBT modes x longitudinal functions (`teacher.py`).
  Output space y = whitened M:  y = W^-1 M D, with
    W = Kbar(k_ref)^-1/2  (mean cross-section stiffness of a sin(k_ref z) shape) and
    D = diag(1 + (k_j / k_ref)^2)  (bending energy of the longitudinal function j),
  so that Euclidean distances in y are close to energy-norm distances of the modes.
  A PCA of the training modes in y (`ModeCodec`) gives n_pc coordinates.

Losses:
  subspace  1 - |Q_pred^T Q_teacher|_F^2 / k   (principal angles of the k-mode spaces:
            independent of signs, scaling and of the order of near-equal modes);
  N_cr      squared error of log N_i (weights 1, 1/2, 1/2).

Physics check (`physics_check`): on the GBT operators of the design, the Rayleigh-Ritz
value of the predicted modes (an upper bound of the GBT N_cr in this representation)
and the energy residual eta = |G x - mu K x|_T / mu (T: wavenumber-block inverse of K)
of the Ritz mode; together with the ensemble spread they decide when the oracle is
called (replacing the fixed -15 % / +3 % suspicion window of Paper 1).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch import nn

K_REF = math.pi / 30.0


def load_teacher(folder, family="omega"):
    """Teacher records; x = design vector + closed-form physics features (omega)."""
    from .designspace import FAMILIES, Sample
    recs = []
    for f in sorted(Path(folder).glob("*.npz")):
        d = np.load(f, allow_pickle=False)
        params, lay = json.loads(str(d["params"])), json.loads(str(d["layups"]))
        x = d["x"].astype(np.float32)
        if family == "omega":
            x = np.r_[x, physics_features(Sample("omega", params, lay, dict(FAMILIES["omega"]["fixed"])))]
        recs.append(dict(sid=str(d["sid"]), params=params, layups=lay, x=x.astype(np.float32),
                         L=float(d["L"]), bc=str(d["bc"]), lam_full=d["lam_full"], lam=d["lam"],
                         M=d["M"].astype(np.float64), J=int(d["J"]), file=str(f)))
    return recs


def cross_section_whitening(ops, k=K_REF, ridge=1e-8):
    """Symmetric Kbar^-1/2 and Kbar^1/2 of the mean cross-section stiffness at wavenumber k."""
    nm = ops[0].nm
    blk = lambda o, d, e: o.K[d * nm:(d + 1) * nm, e * nm:(e + 1) * nm]
    Kb = sum(blk(o, 0, 0) + k ** 2 * (blk(o, 1, 1) - blk(o, 0, 2) - blk(o, 2, 0)) + k ** 4 * blk(o, 2, 2)
             for o in ops) / len(ops)
    Kb = 0.5 * (Kb + Kb.T); Kb += ridge * np.trace(Kb) / nm * np.eye(nm)
    ev, V = np.linalg.eigh(Kb); ev = np.maximum(ev, 1e-12 * ev.max())
    return (V / np.sqrt(ev)) @ V.T, (V * np.sqrt(ev)) @ V.T


def fourier_wavenumbers(J, L):
    return np.r_[np.arange(J + 1), np.arange(1, J + 1)] * math.pi / L


@dataclass
class ModeCodec:
    """Whitening + PCA of the mode tensors."""
    Wm: np.ndarray            # Kbar^-1/2 ... applied as y = Wh^{1/2}... (see encode)
    Wp: np.ndarray            # Kbar^+1/2
    basis: np.ndarray         # (n_y, n_pc) PCA basis in the whitened space
    sv: np.ndarray
    J: int

    def scale(self, L):
        return 1.0 + (fourier_wavenumbers(self.J, L) / K_REF) ** 2

    def encode(self, M, L):
        """M (k, nm, nf) -> PCA coefficients (k, n_pc).  y = Kbar^1/2 M D."""
        y = np.einsum("mn,knf->kmf", self.Wp, M) * self.scale(L)[None, None, :]
        return y.reshape(y.shape[0], -1) @ self.basis

    def decode(self, c, L, nm):
        """PCA coefficients (k, n_pc) -> M (k, nm, nf)."""
        y = (c @ self.basis.T).reshape(c.shape[0], nm, -1) / self.scale(L)[None, None, :]
        return np.einsum("mn,knf->kmf", self.Wm, y)

    @staticmethod
    def fit(recs, Wm, Wp, n_pc=256):
        Y = []
        for r in recs:
            sc = 1.0 + (fourier_wavenumbers(r["J"], r["L"]) / K_REF) ** 2
            y = np.einsum("mn,knf->kmf", Wp, r["M"].astype(float)) * sc[None, None, :]
            y /= np.linalg.norm(y.reshape(y.shape[0], -1), axis=1)[:, None, None]
            Y.append(y.reshape(y.shape[0], -1))
        Y = np.concatenate(Y, 0)
        _, sv, Vt = np.linalg.svd(Y, full_matrices=False)
        n = min(n_pc, Vt.shape[0])
        return ModeCodec(Wm, Wp, Vt[:n].T.copy(), sv[:n], recs[0]["J"])


def physics_features(sample):
    """Classical closed-form buckling estimates used as extra network inputs (log N/mm):
    long simply supported orthotropic plates for the skin bay, the stringer web and the
    crown, N = pi^2 / b^2 (2 sqrt(D11 D22) + 2 (D12 + 2 D66)), and the Euler column load
    per unit width of the whole panel (pi^2 EI / (w L^2), x 4 for clamped ends)."""
    from .oracle import omega_length, omega_panel
    pan = omega_panel(sample); nd = pan["nodes"]; L = omega_length(sample)
    plate = lambda D, b: math.pi ** 2 / b ** 2 * (2 * math.sqrt(D[0, 0] * D[1, 1]) + 2 * (D[0, 1] + 2 * D[2, 2]))
    best = {}
    ys, EA, walls = [], [], []
    for w in pan["walls"]:
        P0, P1 = np.asarray(nd[w["i"]], float), np.asarray(nd[w["j"]], float)
        b = float(np.hypot(*(P1 - P0))); A, B, D, _ = w["ABD"]
        part = w.get("part", "")
        if part in ("skin", "web", "crown"):
            best[part] = max(best.get(part, 0.0), b) if part != "skin" else max(best.get(part, 0.0), b)
            best[part + "_D"] = D
        walls.append((P0, P1, b, A, D))
    feats = [plate(best[p + "_D"], best[p]) for p in ("skin", "web", "crown")]
    ybar = sum(b * A[0, 0] * 0.5 * (P0[1] + P1[1]) for P0, P1, b, A, D in walls) / sum(b * A[0, 0] for P0, P1, b, A, D in walls)
    EI = 0.0
    for P0, P1, b, A, D in walls:
        y0, y1 = P0[1] - ybar, P1[1] - ybar
        EI += A[0, 0] * b * (y0 * y0 + y0 * y1 + y1 * y1) / 3.0 + D[0, 0] * b * abs(P1[0] - P0[0]) / max(b, 1e-12)
    width = float(np.ptp([p[0] for p in nd]))
    fac = 4.0 if sample.params["bc"] == "CC" else 1.0
    feats.append(fac * math.pi ** 2 * EI / (width * L ** 2))
    return np.log(np.asarray(feats, np.float32))


class BucklingNet(nn.Module):
    def __init__(self, n_in, n_pc, k=3, width=512, depth=4):
        super().__init__()
        layers, d = [], n_in
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.SiLU()]; d = width
        self.body = nn.Sequential(*layers)
        self.modes = nn.Linear(width, k * n_pc)
        self.lam = nn.Linear(width, k)
        self.k, self.n_pc = k, n_pc

    def forward(self, x):
        h = self.body(x)
        return self.modes(h).reshape(*x.shape[:-1], self.k, self.n_pc), self.lam(h)   # coefficients, log N


def subspace_loss(Cp, Ct):
    """Cp, Ct (B, k, n): 1 - |Q_p^T Q_t|_F^2 / k per sample."""
    Qp, _ = torch.linalg.qr(Cp.transpose(1, 2)); Qt, _ = torch.linalg.qr(Ct.transpose(1, 2))
    return 1.0 - (Qp.transpose(1, 2) @ Qt).pow(2).sum((1, 2)) / Cp.shape[1]


def physics_check(fp, M_pred, iters=(0, 2, 5)):
    """Physics check of a prediction on the GBT operators of the design.

    Starting from the predicted modes, block-LOBPCG iterations on the exact GBT energies
    (network representation, wavenumber-block preconditioner) give after `it` iterations
    the Ritz value N_it, an upper bound of the GBT critical load in this representation,
    and the relative energy residual eta_it = |G x - mu K x|_T / mu of the Ritz mode
    (T: wavenumber-block inverse of K).  Raw network modes carry small errors in stiff
    directions that inflate N_0; a few iterations remove them, so N_it and eta_it measure
    how far the prediction is from a true buckling mode.  Returns {it: (N_it, eta_it)}."""
    from .spectral import WavePrecond, cat, lincomb, lobpcg
    from .teacher import from_modes
    A = from_modes(fp, torch.as_tensor(M_pred, dtype=fp.dtype))
    T = WavePrecond(fp)
    out = {}
    for it in iters:
        if it == 0:
            from .spectral import rayleigh_ritz
            KA, GA = fp.apply_K(A), fp.apply_G(A)
            mu, C = rayleigh_ritz(A, KA, GA, 1)
            X, KX, GX = lincomb(A, C), lincomb(KA, C), lincomb(GA, C)
        else:
            lam, X, _ = lobpcg(fp, R=A[0].shape[0], iters=it, X0=A, precond=T, tol=0.0)
            X = [x[:1] for x in X]
            KX, GX = fp.apply_K(X), fp.apply_G(X)
            mu = (sum((x * g).sum() for x, g in zip(X, GX)) / sum((x * k).sum() for x, k in zip(X, KX))).reshape(1)
            nrm = torch.sqrt(sum((x * k).sum() for x, k in zip(X, KX)))
            X = [x / nrm for x in X]; KX = [k / nrm for k in KX]; GX = [g / nrm for g in GX]
        res = [g - mu[:, None, None] * k for g, k in zip(GX, KX)]
        eta = math.sqrt(abs(float(sum((r * t).sum() for r, t in zip(res, T(res)))))) / abs(float(mu[0]))
        out[it] = (1.0 / float(mu[0]), eta)
    return out
