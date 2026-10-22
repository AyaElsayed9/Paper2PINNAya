"""Spectral (Fourier-in-z) GBT operators of a panel in the shared basis, with a
wavenumber-block preconditioner and a block LOBPCG eigen-iteration.

Unknowns: for every end-condition class c (`operators.py`) a coefficient tensor
A_c (k_c x nf) on the longitudinal functions F_c = b_c(zeta) * {cos j pi zeta (j = 0..J),
sin j pi zeta (j = 1..J)}; the mode is q(z) = sum_c Pi_c A_c F_c(z)^T.

    U   = 1/2 a^T K a,      U_G = 1/2 a^T G a      (K SPD, G indefinite)

`apply_K`, `apply_G` return the gradients K a, G a (dual tensors of the same shape),
computed at the solver's Gauss points from the exact GBT operators.

Preconditioner (`WavePrecond`): the exact diagonal blocks of K for every wavenumber j
(functions cos j pi zeta and sin j pi zeta of all classes together, 2 nm unknowns).  For a
prismatic panel the energy couples different wavenumbers only through the ends and
the end factors, so these blocks carry most of K; their inverse is a spectral GBT
preconditioner (one cross-section problem per wavenumber).
"""
from __future__ import annotations

import math

import numpy as np
import torch

from .operators import PanelOperator, bubble, fourier_basis

_DIR_D = {0: (1, 2), 1: (0, 1)}       # nonzero amplitude-derivative orders of the z- and s-gradients


class FourierPanel:
    def __init__(self, op: PanelOperator, J: int = 24, dtype=torch.float64, basis="bubble", sets=None):
        """basis = "bubble": every class uses cos/sin j pi zeta times its end factor.
        basis = "harmonic": every class uses one trigonometric family that meets its end
        conditions (`sets`: class name -> family, see `harmonic_sets`)."""
        t = lambda a: torch.as_tensor(np.ascontiguousarray(a), dtype=dtype)
        self.op = op; self.J = J; self.dtype = dtype
        self.L, self.nm, self.sid, self.params = op.L, op.nm, op.sid, op.params
        self.lam_full = op.lam_full
        zeta = op.z / op.L
        self.wz = t(op.wz); self.nz = op.z.size
        self.F, self.Pi, self.jc = [], [], []
        if basis == "bubble":
            Fz = fourier_basis(J, zeta, op.L)                          # (3, nz, nf)
            for c in op.classes:
                bb = bubble(zeta, c.left, c.right, op.L)
                b0, b1, b2 = bb[0][:, None], bb[1][:, None], bb[2][:, None]
                F = np.stack([b0 * Fz[0], b1 * Fz[0] + b0 * Fz[1], b2 * Fz[0] + 2 * b1 * Fz[1] + b0 * Fz[2]])
                self.F.append(t(F)); self.Pi.append(t(c.Pi)); self.jc.append(np.r_[np.arange(J + 1), np.arange(1, J + 1)])
        else:
            sets = sets or harmonic_sets(op.bc)
            for c in op.classes:
                F, jc = trig_family(sets[c.name], J, zeta, op.L)
                self.F.append(t(F)); self.Pi.append(t(c.Pi)); self.jc.append(jc)
        self.k = [c.Pi.shape[1] for c in op.classes]
        nm = op.nm
        self.K = t(op.K.reshape(3, nm, 3, nm).transpose(0, 2, 1, 3))   # (3, 3, nm, nm)
        npnt, nc = op.G.shape[:2]
        G5 = op.G.astype(float).reshape(npnt, nc, 2, 3, nm)
        self.Ga = {(a, d): t(G5[:, :, a, d, :].reshape(npnt * nc, nm)) for a in (0, 1) for d in _DIR_D[a]}
        P = op.P.astype(float)                                         # (3, nz, np), includes w_p
        self.wP = [t(np.repeat(P[i], nc, axis=1)) for i in range(3)]   # (nz, np*nc)
        self.nfs = [F.shape[2] for F in self.F]

    # ---- vectors are lists over classes of tensors (R, k_c, nf) -----------------------
    def zeros(self, R):
        return [torch.zeros(R, k, nf, dtype=self.dtype) for k, nf in zip(self.k, self.nfs)]

    def random(self, R, gen=None):
        return [torch.randn(R, k, nf, dtype=self.dtype, generator=gen) for k, nf in zip(self.k, self.nfs)]

    def field(self, A):
        """Q (R, 3, nz, nm): amplitude derivatives at the Gauss points."""
        Q = 0.0
        for Ac, F, Pi in zip(A, self.F, self.Pi):
            Q = Q + torch.einsum("rkf,dzf,mk->rdzm", Ac, F, Pi)
        return Q

    def dual(self, D):
        """Coefficient-space gradient of sum_z w_z D(z) . Q(z) for D (R, 3, nz, nm)."""
        Dw = D * self.wz[None, None, :, None]
        return [torch.einsum("rdzm,dzf,mk->rkf", Dw, F, Pi) for F, Pi in zip(self.F, self.Pi)]

    def apply_K(self, A):
        Q = self.field(A)
        KQ = torch.einsum("demn,rezn->rdzm", self.K, Q)
        return self.dual(KQ)

    def apply_G(self, A):
        Q = self.field(A)
        g = {a: sum(Q[:, d] @ self.Ga[(a, d)].T for d in _DIR_D[a]) for a in (0, 1)}   # (R, nz, X)
        h0 = self.wP[0] * g[0] + self.wP[2] * g[1]
        h1 = self.wP[1] * g[1] + self.wP[2] * g[0]
        GQ = torch.zeros_like(Q)
        for (a, d), Gm in self.Ga.items():
            GQ[:, d] = GQ[:, d] + (h0 if a == 0 else h1) @ Gm
        return self.dual(GQ)


def harmonic_sets(bc):
    """Trigonometric family of every end-condition class.  SS: phi-type (zero ends) sin,
    psi-type (free) cos, axial (zero at z = 0) quarter-wave sin; CC: clamped
    conventional modes sin(j pi zeta) sin(pi zeta), axial and shear quarter-wave sin,
    Mindlin rotations (zero ends) sin."""
    if bc == "SS":
        return dict(conv_fixed="sin", conv_anchor="qsin", conv_free="sincos", axial="qsin",
                    shear="cos", mx="cos", ms="sincos")
    return dict(conv="clamp", axial="qsin", shear="qsin", mx="sin", ms="sin")


def trig_family(kind, J, zeta, L):
    """Functions (3, nz, n) [value, d/dz, d2/dz2] and their wavenumber index (n,)."""
    out, jc = [], []
    def add(f, f1, f2, j):
        out.append(np.stack([f, f1 / L, f2 / L ** 2])); jc.append(j)
    for j in range(0, J + 1):
        k = j * np.pi
        if kind in ("cos", "sincos"):
            add(np.cos(k * zeta), -k * np.sin(k * zeta), -k * k * np.cos(k * zeta), j)
        if kind in ("sin", "sincos") and j >= 1:
            add(np.sin(k * zeta), k * np.cos(k * zeta), -k * k * np.sin(k * zeta), j)
        if kind == "qsin" and j >= 1:
            q = (j - 0.5) * np.pi
            add(np.sin(q * zeta), q * np.cos(q * zeta), -q * q * np.sin(q * zeta), j)
        if kind == "clamp" and j >= 1:                        # sin(j pi z) sin(pi z): zero value and slope
            a, b = np.sin(k * zeta), np.sin(np.pi * zeta)
            a1, b1 = k * np.cos(k * zeta), np.pi * np.cos(np.pi * zeta)
            a2, b2 = -k * k * a, -np.pi ** 2 * b
            add(a * b, a1 * b + a * b1, a2 * b + 2 * a1 * b1 + a * b2, j)
    F = np.stack(out, axis=2)                                  # (3, nz, n)
    return F, np.array(jc)


def dot(A, B):
    """(R_A, R_B) matrix of inner products."""
    return sum(torch.einsum("rkf,skf->rs", a, b) for a, b in zip(A, B))


def lincomb(A, C):
    """Columns combined: A (R vectors) times C (R, R') -> R' vectors."""
    return [torch.einsum("rkf,rs->skf", a, C) for a in A]


def cat(*Vs):
    return [torch.cat(parts, 0) for parts in zip(*Vs)]


class WavePrecond:
    """Inverse of the wavenumber-diagonal blocks of K (all classes, cos and sin of one j)."""

    def __init__(self, fp: FourierPanel, shift=1e-8):
        self.fp = fp; self.shift = shift
        nC = len(fp.k); J = fp.J
        # Gram of all longitudinal functions of all classes and derivative orders
        Fs = [F.permute(0, 2, 1) for F in fp.F]                          # (3, nf_c, nz)
        Zc = {(c, c2): torch.einsum("dfz,z,egz->dfeg", Fs[c], fp.wz, Fs[c2]) for c in range(nC) for c2 in range(nC)}
        KP = [[torch.einsum("mk,demn,nl->dekl", fp.Pi[c], fp.K, fp.Pi[c2]) for c2 in range(nC)] for c in range(nC)]
        # geometric stiffness with the z-averaged pre-stress (wavenumber-diagonal part of G)
        L = float(fp.wz.sum())
        Pm = [(fp.wz[:, None] * w).sum(0) / L for w in fp.wP]                # (X,)
        Gm = torch.zeros(3, 3, fp.nm, fp.nm, dtype=fp.dtype)
        for (a, d), Ga in fp.Ga.items():
            for (b, e), Gb in fp.Ga.items():
                wt = Pm[0] if (a, b) == (0, 0) else Pm[1] if (a, b) == (1, 1) else Pm[2]
                Gm[d, e] += Ga.T @ (wt[:, None] * Gb)
        GP = [[torch.einsum("mk,demn,nl->dekl", fp.Pi[c], Gm, fp.Pi[c2]) for c2 in range(nC)] for c in range(nC)]
        self.parts = []
        for j in range(J + 1):
            idx = [(c, f) for c in range(nC) for f in np.nonzero(fp.jc[c] == j)[0]]
            if not idx:
                continue
            n = [fp.k[c] for c, f in idx]; off = np.r_[0, np.cumsum(n)]
            Mk = torch.zeros(off[-1], off[-1], dtype=fp.dtype); Mg = torch.zeros_like(Mk)
            for i, (c, f) in enumerate(idx):
                for i2, (c2, f2) in enumerate(idx):
                    zz = Zc[(c, c2)][:, f, :, f2]
                    Mk[off[i]:off[i + 1], off[i2]:off[i2 + 1]] = torch.einsum("dekl,de->kl", KP[c][c2], zz)
                    Mg[off[i]:off[i + 1], off[i2]:off[i2 + 1]] = torch.einsum("dekl,de->kl", GP[c][c2], zz)
            self.parts.append((idx, off, 0.5 * (Mk + Mk.T), 0.5 * (Mg + Mg.T)))
        self.set_shift(0.0)

    def set_shift(self, sigma):
        """Factor K_j - sigma G_j (SPD for sigma below the critical load of the block)."""
        self.sigma = sigma; self.blocks = []
        for idx, off, Mk, Mg in self.parts:
            M = Mk - sigma * Mg
            M = M + self.shift * torch.diagonal(Mk).mean() * torch.eye(M.shape[0], dtype=M.dtype)
            Lc, info = torch.linalg.cholesky_ex(M)
            if int(info) != 0:                                          # not SPD: fall back to K_j
                Lc = torch.linalg.cholesky(Mk + self.shift * torch.diagonal(Mk).mean()
                                           * torch.eye(M.shape[0], dtype=M.dtype))
            self.blocks.append((idx, off, Lc))

    def __call__(self, Rv):
        out = [torch.zeros_like(r) for r in Rv]
        R = Rv[0].shape[0]
        for idx, off, Lc in self.blocks:
            b = torch.cat([Rv[c][:, :, f] for c, f in idx], 1)          # (R, n)
            x = torch.cholesky_solve(b.T, Lc).T
            for i, (c, f) in enumerate(idx):
                out[c][:, :, f] = x[:, off[i]:off[i + 1]]
        return out


def rayleigh_ritz(S, KS, GS, nkeep, drop=1e-12):
    """Ritz pairs (largest mu) of G x = mu K x on span(S); returns mu (nkeep,), coefficients."""
    Ks = dot(S, KS); Gs = dot(S, GS)
    Ks = 0.5 * (Ks + Ks.T); Gs = 0.5 * (Gs + Gs.T)
    ev, V = torch.linalg.eigh(Ks)
    keep = ev > drop * ev.max()
    W = V[:, keep] / torch.sqrt(ev[keep])
    mu, Y = torch.linalg.eigh(W.T @ Gs @ W)
    order = torch.argsort(mu, descending=True)[:nkeep]
    return mu[order], W @ Y[:, order]


def lobpcg(fp: FourierPanel, R=6, iters=40, X0=None, precond=None, tol=1e-7, callback=None, seed=0,
           shift_every=0, shift_frac=0.9):
    """Block LOBPCG for the largest mu of G x = mu K x (N_cr = 1 / mu_max)."""
    gen = torch.Generator().manual_seed(seed)
    T = precond or WavePrecond(fp)
    X = X0 if X0 is not None else T(fp.random(R, gen))
    KX, GX = fp.apply_K(X), fp.apply_G(X)
    mu, C = rayleigh_ritz(X, KX, GX, R)
    X, KX, GX = lincomb(X, C), lincomb(KX, C), lincomb(GX, C)
    P = KP = GP = None
    hist = []
    for it in range(iters):
        if shift_every and it % shift_every == 0 and it > 0:           # shift-invert flavour
            T.set_shift(shift_frac / float(mu[0]))
        res = [g - mu[:, None, None] * k for g, k in zip(GX, KX)]      # G x - mu K x
        W = T(res); KW, GW = fp.apply_K(W), fp.apply_G(W)
        if P is None:
            S, KS, GS = cat(X, W), cat(KX, KW), cat(GX, GW)
        else:
            S, KS, GS = cat(X, W, P), cat(KX, KW, KP), cat(GX, GW, GP)
        mu_new, C = rayleigh_ritz(S, KS, GS, R)
        Xn, KXn, GXn = lincomb(S, C), lincomb(KS, C), lincomb(GS, C)
        # implicit P: the part of the new X outside the old X
        Cp = C.clone(); Cp[:R] = 0.0
        P, KP, GP = lincomb(S, Cp), lincomb(KS, Cp), lincomb(GS, Cp)
        X, KX, GX = Xn, KXn, GXn
        rel = float((mu_new[0] - mu[0]).abs() / mu_new[0].abs())
        mu = mu_new
        hist.append(1.0 / float(mu[0]))
        if callback:
            callback(it, mu)
        if rel < tol:
            break
    return 1.0 / mu, X, hist
