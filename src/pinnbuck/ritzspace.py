"""Energy space of one design for the physics-only PINN: dictionary Ritz space,
Kronecker stiffness, pointwise geometric stiffness and the energy-whitening factor.

Unknowns.  For every end-condition class c (`operators.py`) a coefficient matrix
A_c (k_c x nb_c) on longitudinal functions B_c(z) built from the physical dictionary
(`dictionary.py`: bulk waves 25-350 mm + end layers), so that the mode is

    q(z) = sum_c Pi_c A_c B_c(z)^T,        U = 1/2 a^T K a,   U_G = 1/2 a^T G a.

Kronecker stiffness.  The panel is prismatic, so K is a sum of Kronecker products,
K[(c,f),(c',g)] = sum_{d,e} (Pi_c^T K_de Pi_c')[k,l] Z^{cc'}_de[f,g],
Z^{cc'}_de = int B_c^(d) B_c'^(e) dz; K a costs ~10^7 flops (`apply_K`).  G is applied
pointwise along z with the full-GBT pre-buckling state (`apply_G`, as in spectral.py).

Longitudinal functions.  Per class the dictionary is orthonormalised and then rotated
to the eigenfunctions of the class "beam pencil" Z_22 t = w Z_00 t, so every function
has a wavenumber kappa = w^(1/4) and Z_00, Z_22 are diagonal in each class.

Energy whitening (`EnergyMetric`).  K is approximated by its blocks over wavenumber bins
(all classes and cross-section modes of the functions whose kappa lies in one bin;
couplings between bins are neglected), K ~ M = C^T C with C block upper-triangular.
In the coordinates y = C a the strain energy is ~ |y|^2 / 2 for every design, so
plain gradient steps on y are energy natural-gradient steps (Mueller & Zeinhofer 2023,
metric = energy Hessian K), and M^-1 is the preconditioner of the eigen-iterations.
"""
from __future__ import annotations

import math

import numpy as np
import torch

from .dictionary import END_LENGTH, HALF_WAVES, N_END, dictionary
from .operators import PanelOperator, orthonormalise

_DIR_D = {0: (1, 2), 1: (0, 1)}       # amplitude-derivative orders entering the z- and s-gradients


def raw_functions(op: PanelOperator, half_waves=HALF_WAVES, n_end=N_END, end_length=END_LENGTH):
    """Per class: the admissible dictionary functions themselves (fixed labels for every
    design of an end condition), their labels and wavenumbers [1/mm]."""
    D, labels = dictionary(op, half_waves=half_waves, n_end=n_end, end_length=end_length)
    le = min(end_length * op.params["bs"], 0.5 * op.L)
    kw = lambda lab: math.pi / float(lab[4:-1]) if lab.startswith("wave") else math.pi / (2 * le / n_end)
    out, labs, kap = [], [], []
    for B in D:
        keep = np.abs(B).sum(axis=(0, 1)) > 0
        out.append(B[:, :, keep]); labs.append([l for l, k_ in zip(labels, keep) if k_])
        kap.append(np.array([kw(l) for l in labs[-1]]))
    return out, labs, kap


def class_functions(op: PanelOperator, half_waves=HALF_WAVES, n_end=N_END, end_length=END_LENGTH, tol=1e-10):
    """Per class: longitudinal functions (3, nz, nb_c), orthonormal in Z_00 and rotated to
    the eigenfunctions of Z_22 t = w Z_00 t, and their wavenumbers kappa_c (nb_c,)."""
    D, _ = dictionary(op, half_waves=half_waves, n_end=n_end, end_length=end_length)
    out, kap = [], []
    for B in D:
        B = B[:, :, np.abs(B).sum(axis=(0, 1)) > 0]
        B = orthonormalise(B, op.wz, op.L, tol=tol)                   # non-redundant (H2)
        Z0 = (B[0] * op.wz[:, None]).T @ B[0]; Z2 = (B[2] * op.wz[:, None]).T @ B[2]
        ev0, V0 = np.linalg.eigh(0.5 * (Z0 + Z0.T))
        T0 = V0 / np.sqrt(np.maximum(ev0, 1e-300))                    # Z0-orthonormal
        w, V = np.linalg.eigh(T0.T @ (0.5 * (Z2 + Z2.T)) @ T0)
        T = T0 @ V
        out.append(np.einsum("dzf,fg->dzg", B, T))
        kap.append(np.maximum(w, 0.0) ** 0.25)
    return out, kap


class RitzSpace:
    """Operators of one design on its dictionary Ritz space (vectors: lists over classes of
    tensors (R, k_c, nb_c))."""

    def __init__(self, op: PanelOperator, dtype=torch.float64, raw=False, **kw):
        """raw=True: the admissible dictionary functions themselves (same coordinates for
        every design with the same end condition; `labels`), else the design's
        orthonormal beam-pencil functions."""
        t = lambda a: torch.as_tensor(np.ascontiguousarray(a), dtype=dtype)
        self.op, self.dtype, self.raw = op, dtype, raw
        self.sid, self.params, self.L, self.nm = op.sid, op.params, op.L, op.nm
        self.lam_full = op.lam_full
        if raw:
            Bn, self.labels, kap = raw_functions(op, **kw)
        else:
            Bn, kap = class_functions(op, **kw); self.labels = None
        self.B = [t(b) for b in Bn]; self.kappa = kap
        self.Pi = [t(c.Pi) for c in op.classes]
        self.k = [c.Pi.shape[1] for c in op.classes]; self.nb = [b.shape[2] for b in Bn]
        self.n = sum(k * nb for k, nb in zip(self.k, self.nb))
        self.wz = t(op.wz)
        nm = op.nm; K4 = op.K.reshape(3, nm, 3, nm)
        nC = len(self.k)
        # Kronecker terms of K per class pair (exactly zero cross-section blocks skipped)
        self.KZ = {}
        for c in range(nC):
            for c2 in range(nC):
                terms = []
                for d in range(3):
                    for e in range(3):
                        Kc = op.classes[c].Pi.T @ K4[d, :, e, :] @ op.classes[c2].Pi
                        if not np.any(Kc):
                            continue
                        Z = (Bn[c][d] * op.wz[:, None]).T @ Bn[c2][e]
                        terms.append((t(Kc), t(Z)))
                self.KZ[(c, c2)] = terms
        npnt, nc = op.G.shape[:2]
        G5 = op.G.astype(float).reshape(npnt, nc, 2, 3, nm)
        # pointwise gradient matrices, compressed to their nonzero rows and columns
        self.Ga = {}
        for a in (0, 1):
            for d in _DIR_D[a]:
                Gm = G5[:, :, a, d, :].reshape(npnt * nc, nm)
                rows = np.nonzero(np.abs(Gm).sum(1) > 0)[0]; cols = np.nonzero(np.abs(Gm).sum(0) > 0)[0]
                self.Ga[(a, d)] = (torch.as_tensor(rows), torch.as_tensor(cols), t(Gm[np.ix_(rows, cols)]))
        self.nX = npnt * nc
        P = op.P.astype(float)
        self.wP = [t(np.repeat(P[i], nc, axis=1)) for i in range(3)]          # (nz, np*nc)
        self.features = op.features

    # ---- vectors -------------------------------------------------------------------
    def zeros(self, R):
        return [torch.zeros(R, k, nb, dtype=self.dtype) for k, nb in zip(self.k, self.nb)]

    def random(self, R, gen=None):
        return [torch.randn(R, k, nb, dtype=self.dtype, generator=gen) for k, nb in zip(self.k, self.nb)]

    def flat(self, A):
        return torch.cat([a.reshape(a.shape[0], -1) for a in A], 1)

    def unflat(self, X):
        out, o = [], 0
        for k, nb in zip(self.k, self.nb):
            out.append(X[:, o:o + k * nb].reshape(-1, k, nb)); o += k * nb
        return out

    def field(self, A):
        """Q (R, 3, nz, nm)."""
        return sum(torch.einsum("rkf,dzf,mk->rdzm", a, B, Pi) for a, B, Pi in zip(A, self.B, self.Pi))

    def dual(self, D):
        Dw = D * self.wz[None, None, :, None]
        return [torch.einsum("rdzm,dzf,mk->rkf", Dw, B, Pi) for B, Pi in zip(self.B, self.Pi)]

    # ---- operators -----------------------------------------------------------------
    def apply_K(self, A):
        out = []
        for c in range(len(self.k)):
            acc = 0.0
            for c2, a2 in enumerate(A):
                for Kc, Z in self.KZ[(c, c2)]:
                    acc = acc + torch.einsum("kl,rlg,fg->rkf", Kc, a2, Z)
            out.append(acc)
        return out

    def apply_G(self, A):
        Q = self.field(A)
        R, _, nz, nm = Q.shape
        g = {}
        for (a, d), (rows, cols, Gm) in self.Ga.items():
            v = Q[:, d][..., cols] @ Gm.T                                     # (R, nz, n_rows)
            if rows.numel() == self.nX:
                g[a] = v if a not in g else g[a] + v
            else:
                full = torch.zeros(R, nz, self.nX, dtype=Q.dtype).index_copy(2, rows, v)
                g[a] = full if a not in g else g[a] + full
        h = {0: self.wP[0] * g[0] + self.wP[2] * g[1], 1: self.wP[1] * g[1] + self.wP[2] * g[0]}
        parts = {d: [] for d in range(3)}
        for (a, d), (rows, cols, Gm) in self.Ga.items():
            ha = h[a] if rows.numel() == self.nX else h[a].index_select(2, rows)
            parts[d].append(torch.zeros(R, nz, nm, dtype=Q.dtype).index_add(2, cols, ha @ Gm))
        GQ = torch.stack([sum(parts[d]) if parts[d] else torch.zeros(R, nz, nm, dtype=Q.dtype) for d in range(3)], 1)
        return self.dual(GQ)

    def dense_K(self):
        """Dense K (n x n), for checks."""
        n = self.n; I = torch.eye(n, dtype=self.dtype)
        cols = [self.flat(self.apply_K(self.unflat(I[i:i + 512]))) for i in range(0, n, 512)]
        K = torch.cat(cols, 0)
        return 0.5 * (K + K.T)


class EnergyMetric:
    """Block approximation M = C^T C of K over wavenumber bins (see module doc).
    `whiten(A)` = C a,  `color(Y)` = C^-1 y,  `solve(R)` = M^-1 r."""

    def __init__(self, rs: RitzSpace, ratio=1.35, kmin=None, max_block=2500, groups=None, ridge=1e-12):
        """groups: None = bins of log kappa; "labels" (raw spaces) = one block per dictionary
        wavelength (sin + cos) and per end layer (fixed for every design)."""
        self.rs = rs
        nC = len(rs.k)
        if groups == "labels":
            key = lambda lab: lab[:-1] if lab.startswith("wave") else lab[:4]
            names = sorted({key(l) for labs in rs.labels for l in labs})
            bins = [[(c, f) for c in range(nC) for f, l in enumerate(rs.labels[c]) if key(l) == nm_]
                    for nm_ in names]
        else:
            kmin = kmin or math.pi / 2000.0
            items = [(c, f, float(rs.kappa[c][f])) for c in range(nC) for f in range(rs.nb[c])]
            lk = np.array([math.log(max(k_, kmin)) for _, _, k_ in items])
            order = np.argsort(lk, kind="stable")
            # greedy bins: width log(ratio) in log kappa, capped in size
            bins, cur, start, size = [], [], None, 0
            for i in order:
                c, f, _ = items[i]
                if cur and (lk[i] - start > math.log(ratio) or size + rs.k[c] > max_block):
                    bins.append(cur); cur, size = [], 0
                if not cur:
                    start = lk[i]
                cur.append((c, f)); size += rs.k[c]
            if cur:
                bins.append(cur)
        self.bins = []
        for b in bins:
            b = sorted(b)                                         # grouped by class, then function
            n = [rs.k[c] for c, _ in b]; off = np.r_[0, np.cumsum(n)]
            M = torch.zeros(off[-1], off[-1], dtype=rs.dtype)
            cls = sorted({c for c, _ in b})
            fs = {c: [f for c2, f in b if c2 == c] for c in cls}
            start = {c: int(off[[i for i, (c2, _) in enumerate(b) if c2 == c][0]]) for c in cls}
            for c in cls:
                for c2 in cls:
                    terms = rs.KZ[(c, c2)]
                    if not terms:
                        continue
                    Kt = torch.stack([Kc for Kc, _ in terms])                      # (T, k, k2)
                    fi = torch.as_tensor(fs[c]); gi = torch.as_tensor(fs[c2])
                    Zt = torch.stack([Z[fi][:, gi] for _, Z in terms])             # (T, nf, ng)
                    blk = torch.einsum("tkl,tfg->fkgl", Kt, Zt).reshape(len(fs[c]) * rs.k[c], len(fs[c2]) * rs.k[c2])
                    r0, c0 = start[c], start[c2]
                    M[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]] = blk
            M = 0.5 * (M + M.T)
            M = M + ridge * torch.diagonal(M).mean() * torch.eye(M.shape[0], dtype=M.dtype)
            Lc = torch.linalg.cholesky(M)                 # M = L L^T, C = L^T
            self.bins.append((b, off, Lc))
        self.sizes = [Lc.shape[0] for _, _, Lc in self.bins]

    def _gather(self, A, b, off):
        return torch.cat([A[c][:, :, f] for c, f in b], 1)                     # (R, n_bin)

    def _scatter(self, out, X, b, off):
        for i, (c, f) in enumerate(b):
            out[c][:, :, f] = X[:, off[i]:off[i + 1]]

    def whiten(self, A):
        out = [torch.zeros_like(a) for a in A]
        for b, off, Lc in self.bins:
            self._scatter(out, self._gather(A, b, off) @ Lc, b, off)          # y^T = a^T L
        return out

    def color(self, Y):
        out = [torch.zeros_like(y) for y in Y]
        for b, off, Lc in self.bins:
            x = torch.linalg.solve_triangular(Lc.T, self._gather(Y, b, off).T, upper=True).T
            self._scatter(out, x, b, off)
        return out

    def solve(self, Rv):
        out = [torch.zeros_like(r) for r in Rv]
        for b, off, Lc in self.bins:
            self._scatter(out, torch.cholesky_solve(self._gather(Rv, b, off).T, Lc).T, b, off)
        return out


class _RitzForms(torch.autograd.Function):
    """(K_R, G_R) = (X^T K X, X^T G X) for the flat candidate matrix X (R, n).  The backward
    pass reuses K X and G X from the forward pass (d tr(S^T X^T M X) = M X (S + S^T)),
    so one energy evaluation costs one K and one G product."""

    @staticmethod
    def forward(ctx, X, rs):
        A = rs.unflat(X)
        KX = rs.flat(rs.apply_K(A)); GX = rs.flat(rs.apply_G(A))
        ctx.save_for_backward(KX, GX)
        return X @ KX.T, X @ GX.T

    @staticmethod
    def backward(ctx, gK, gG):
        KX, GX = ctx.saved_tensors
        return (gK + gK.T) @ KX + (gG + gG.T) @ GX, None


def ritz_forms(rs: RitzSpace, X):
    """Stiffness and geometric matrices (R x R) of the candidate rows of X (R, n)."""
    Kr, Gr = _RitzForms.apply(X, rs)
    return 0.5 * (Kr + Kr.T), 0.5 * (Gr + Gr.T)


def ritz_values(Kr, Gr, drop=1e-10):
    """Ritz values mu (descending) of Gr x = mu Kr x on the range of Kr; N = 1 / mu."""
    dsc = 1.0 / torch.sqrt(torch.clamp(torch.diagonal(Kr), min=1e-300))
    Ks = Kr * dsc[:, None] * dsc[None, :]; Gs = Gr * dsc[:, None] * dsc[None, :]
    ev, V = torch.linalg.eigh(Ks)
    keep = ev > drop * ev.max()
    W = V[:, keep] / torch.sqrt(ev[keep])
    mu, Y = torch.linalg.eigh(W.T @ Gs @ W)
    return mu.flip(0), (W @ Y.flip(1)) * dsc[:, None]


class FlatMetric:
    """EnergyMetric on flat vectors (R, n), differentiable: whitened coordinates y = C a and
    back, a = C^-1 y, with C block-triangular over the bins (one permutation + per-bin
    triangular products)."""

    def __init__(self, rs: RitzSpace, em: EnergyMetric):
        offs = np.r_[0, np.cumsum([k * nb for k, nb in zip(rs.k, rs.nb)])]
        perm, self.blocks = [], []
        for b, off, Lc in em.bins:
            start = len(perm)
            for c, f in b:                    # flat index of (c, k, f) = offs[c] + k * nb_c + f
                perm += list(offs[c] + np.arange(rs.k[c]) * rs.nb[c] + f)
            self.blocks.append((start, len(perm), Lc))
        self.perm = torch.as_tensor(np.asarray(perm))
        self.inv = torch.empty_like(self.perm); self.inv[self.perm] = torch.arange(self.perm.numel())

    def whiten(self, X):
        Xp = X[:, self.perm]
        return torch.cat([Xp[:, s:e] @ Lc for s, e, Lc in self.blocks], 1)[:, self.inv]

    def color(self, Y):
        Yp = Y[:, self.perm]
        parts = [torch.linalg.solve_triangular(Lc.T.to(Y.dtype), Yp[:, s:e].T, upper=True).T for s, e, Lc in self.blocks]
        return torch.cat(parts, 1)[:, self.inv]
