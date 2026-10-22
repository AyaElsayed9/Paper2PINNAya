"""Parametric GBT deep-Ritz network (the fast PINN of Paper 2).

One network per member family and side condition maps the design vector p (geometry,
end condition, lamination parameters) to R candidate buckling modes of the panel,

    phi_r(z) = sum_c Pi_c Pi_c^T U_r(p) b_c(zeta) psi_r(z; p),

q = phi: amplitudes of the shared GBT modes (`basis.py`); Pi_c, b_c: end-condition
classes and their end factors (`operators.py`), so every end condition of the solver
holds exactly; U_r: cross-section shape, psi_r: longitudinal shape of candidate r.
The buckling factor is the Rayleigh-Ritz value over the candidates (exact eigenvalue
of the R x R problem), an upper bound of the GBT critical load of that panel.  The
network also returns N_cr directly (head) for millisecond predictions.

Longitudinal shapes:
  * GBTWaveRitz (default): psi = e(zeta) sin(kappa z) and e(zeta) cos(kappa z) with a
    physical wavenumber kappa_r(p) [1/mm] and a smooth envelope e (cos j pi zeta, j <= 4).
    The critical half-wavelength depends smoothly on the cross-section, not on L, so the
    network learns a smooth map; the number of half-waves follows from L.
  * GBTDeepRitz: psi_r = sum_j C_rj f_j(zeta) (Fourier, cos and sin j pi zeta).

Losses (`panel_loss`): physics = maximise mu = U_G / U of the Ritz mode (-mu/|mu|, the
gradient of -log mu; U_G is indefinite, so U / U_G is not used); optional teacher =
1 - K-cosine^2 between the candidate span and the lowest full-GBT modes; head =
(log N_head - log N_Ritz)^2.
"""
from __future__ import annotations

import math

import numpy as np
import torch
from torch import nn

from .operators import PanelOperator, bubble, fourier_basis


class TorchPanel:
    """Operators of one panel as torch tensors."""

    def __init__(self, op: PanelOperator, J: int = 24, k0: float = math.pi / 30.0, mode_scale=None,
                 dtype=torch.float64, device="cpu"):
        t = lambda a: torch.as_tensor(np.ascontiguousarray(a), dtype=dtype, device=device)
        self.sid = op.sid; self.L = op.L; self.nm = op.nm; self.k0 = k0; self.dtype = dtype
        self.x = t(op.features)
        self.K = t(op.K); self.wz = t(op.wz)
        npnt, nc, _, n3 = op.G.shape
        G5 = op.G.reshape(npnt, nc, 2, 3, op.nm)
        # Ga[a][d]: (np*nc, nm) gradient direction a (z, s), amplitude derivative d (contiguous)
        self.Ga = [[t(G5[:, :, a, d, :].reshape(npnt * nc, op.nm)) for d in range(3)] for a in range(2)]
        self.P = t(op.P)                                              # (3, nz, np)
        self.Kcol = [self.K[:, e * op.nm:(e + 1) * op.nm].contiguous() for e in range(3)]
        wP = self.wz[None, :, None] * self.P
        self.wP = [wP[t_][:, :, None].expand(-1, -1, nc).reshape(-1).contiguous() for t_ in range(3)]
        self.nc = nc; self.npnt = npnt
        zeta = op.z / op.L
        self.z = t(op.z); self.zeta = t(zeta)
        self.bub = [t(bubble(zeta, c.left, c.right, op.L)) for c in op.classes]   # (3, nz) per class
        if J:                                                         # Fourier set (GBTDeepRitz)
            Fz = fourier_basis(J, zeta, op.L)
            kj = np.r_[np.arange(J + 1), np.arange(1, J + 1)] * np.pi / op.L
            self.F0 = t(Fz / (1.0 + (kj / k0) ** 2)[None, None, :])   # equal bending energy per function
        ms = np.ones(op.nm) if mode_scale is None else np.asarray(mode_scale)
        self.Pi = [t(ms[:, None] * c.Pi) for c in op.classes]       # q = diag(ms) Pi Pi^T U ...
        self.PiT = [t(c.Pi.T) for c in op.classes]
        self.lam_full = op.lam_full; self.lam_shared = op.lam_shared
        self.KT = self.tt = None
        if getattr(op, "teacher", None) is not None:              # K-weighted teacher modes
            Tq = t(op.teacher)                                         # (k, 3, nz, nm)
            k_, _, nz_, nm_ = Tq.shape
            Qc = Tq.permute(0, 2, 1, 3).reshape(k_, nz_, 3 * nm_)
            KTc = Qc @ self.K
            self.tt = torch.einsum("z,kzx,kzx->k", self.wz, KTc, Qc)
            self.KT = KTc.reshape(k_, nz_, 3, nm_).permute(0, 2, 1, 3).contiguous()   # (k, 3, nz, nm)
            self.teacher_eigs = np.asarray(op.teacher_eigs)
        self.params = op.params


def times_bubble_t(psi, bb):
    """(b psi)^(d), d = 0, 1, 2 for psi (3, nz, R), b (3, nz)."""
    b0, b1, b2 = bb[0][:, None], bb[1][:, None], bb[2][:, None]
    return torch.stack([b0 * psi[0], b1 * psi[0] + b0 * psi[1], b2 * psi[0] + 2 * b1 * psi[1] + b0 * psi[2]])


def mode(tp: TorchPanel, U, psi0):
    """Per class: longitudinal functions with the end factor (3, nz, R) and the projected
    cross-section vectors (nm, R)."""
    psi, Uc = [], []
    for bb, Pi, PiT in zip(tp.bub, tp.Pi, tp.PiT):
        psi.append(times_bubble_t(psi0, bb))
        Uc.append(Pi @ (PiT @ U))
    return psi, Uc


def energies(tp: TorchPanel, U, psi0):
    """U and U_G of the single mode sum_r phi_r (candidates added, not Ritz-combined)."""
    psi, Uc = mode(tp, U, psi0)
    nz = tp.wz.numel(); nm = tp.nm
    Q = sum(torch.einsum("dzr,mr->dzm", p_, u_) for p_, u_ in zip(psi, Uc))   # (3, nz, nm)
    Qc = Q.permute(1, 0, 2).reshape(nz, 3 * nm)
    Ue = 0.5 * torch.sum(tp.wz * torch.sum((Qc @ tp.K) * Qc, dim=1))
    Ps = torch.cat([p_[d] for p_ in psi for d in range(3)], dim=1)   # (nz, n_terms)
    g0, g1 = [(Ps @ torch.cat([tp.Ga[a][d] @ u_ for u_ in Uc for d in range(3)], dim=1).T
               ).reshape(nz, tp.npnt, tp.nc) for a in range(2)]
    dens = (tp.P[0] * (g0 * g0).sum(-1) + tp.P[1] * (g1 * g1).sum(-1) + 2 * tp.P[2] * (g0 * g1).sum(-1))
    UG = 0.5 * torch.sum(tp.wz * dens.sum(-1))
    return Ue, UG


def ritz_forms(tp: TorchPanel, U, psi0, return_mode=False):
    """Stiffness and geometric matrices (R x R) of the R candidate modes."""
    psi, Uc = mode(tp, U, psi0)
    nz = tp.wz.numel(); nm = tp.nm; R = U.shape[-1]
    nC = len(Uc)
    # KR[r,s] = sum_{c,d,c',e} (Uc_c^T K_de Uc_c')[r,s] * int psi_c[d](z,r) psi_c'[e](z,s) dz
    Ua = torch.cat(Uc, dim=1)                                         # (nm, nC R)
    Y = torch.stack([tp.Kcol[e] @ Ua for e in range(3)], 0)          # (3e, 3nm, nC R)
    A = torch.stack([torch.stack([Ua.T @ Y[e, d * nm:(d + 1) * nm] for e in range(3)], 0)
                     for d in range(3)], 0).reshape(3, 3, nC, R, nC, R)
    Ps = torch.stack([torch.stack(list(p_), 0) for p_ in psi], 0)      # (nC, 3, nz, R)
    Pf = Ps.permute(0, 1, 3, 2).reshape(nC * 3 * R, nz)              # rows (c, d, r)
    Z = ((Pf * tp.wz) @ Pf.T).reshape(nC, 3, R, nC, 3, R)
    KR = torch.einsum("decrfs,cdrfes->rs", A, Z)
    # g_r(z, p, c, a) = sum_{class, d} psi_c[d](z, r) (G_ad Uc)(p, c, r): one batched GEMM over r
    Ps = torch.stack([p_[d] for p_ in psi for d in range(3)], dim=0).permute(2, 1, 0)   # (R, nz, T)
    g0, g1 = [torch.bmm(Ps, torch.stack([tp.Ga[a][d] @ u_ for u_ in Uc for d in range(3)], 0).permute(2, 0, 1)
                        ).reshape(R, -1) for a in range(2)]            # (R, nz*np*nc)
    wP = tp.wP
    G01 = (g0 * wP[2]) @ g1.T
    GR = (g0 * wP[0]) @ g0.T + (g1 * wP[1]) @ g1.T + G01 + G01.T
    KR, GR = 0.5 * (KR + KR.T), 0.5 * (GR + GR.T)
    return (KR, GR, psi, Uc) if return_mode else (KR, GR)


def teacher_cos2(tp: TorchPanel, KR, psi, Uc, jitter=1e-8):
    """Squared K-cosine between each teacher mode and its best approximation in the span
    of the candidate modes: b^T KR^-1 b / <t, t>_K  (1 = teacher mode reproduced)."""
    b = 0.0
    for p_, u_ in zip(psi, Uc):
        for d in range(3):
            b = b + torch.einsum("z,zr,kzr->rk", tp.wz, p_[d], tp.KT[:, d] @ u_)   # (R, k)
    dsc = 1.0 / torch.sqrt(torch.clamp(torch.diagonal(KR), min=1e-30))
    Ks = KR * dsc[:, None] * dsc[None, :]
    Lc = torch.linalg.cholesky(Ks + jitter * torch.eye(Ks.shape[0], dtype=Ks.dtype))
    y = torch.linalg.solve_triangular(Lc, b * dsc[:, None], upper=False)
    return (y * y).sum(0) / tp.tt


def ritz_mu(KR, GR, jitter=1e-8):
    """Largest mu of GR x = mu KR x (Jacobi-scaled Cholesky); returns mu and the combination x."""
    dsc = 1.0 / torch.sqrt(torch.clamp(torch.diagonal(KR), min=1e-30))
    Ks = KR * dsc[:, None] * dsc[None, :]; Gs = GR * dsc[:, None] * dsc[None, :]
    Lc = torch.linalg.cholesky(Ks + jitter * torch.eye(Ks.shape[0], dtype=Ks.dtype))
    X = torch.linalg.solve_triangular(Lc, Gs, upper=False)          # L^-1 Gs
    A = torch.linalg.solve_triangular(Lc, X.T, upper=False)         # L^-1 Gs L^-T
    ev, V = torch.linalg.eigh(0.5 * (A + A.T))
    y = torch.linalg.solve_triangular(Lc.T, V[:, -1:], upper=True)[:, 0]
    return ev[-1], y * dsc


def _mlp(n_in, width, depth):
    layers, d = [], n_in
    for _ in range(depth):
        layers += [nn.Linear(d, width), nn.SiLU()]; d = width
    return nn.Sequential(*layers)


class GBTDeepRitz(nn.Module):
    """Fourier longitudinal shapes (zeta = z / L)."""

    def __init__(self, n_in: int, nm: int, J: int, R: int = 8, width: int = 256, depth: int = 4):
        super().__init__()
        self.nm, self.R, self.nf = nm, R, 2 * J + 1
        self.body = _mlp(n_in, width, depth)
        self.hU = nn.Linear(width, nm * R)
        self.hC = nn.Linear(width, R * self.nf)
        self.hL = nn.Linear(width, 1)

    def candidates(self, tp):
        h = self.body(tp.x)
        U = self.hU(h).reshape(self.nm, self.R)
        C = self.hC(h).reshape(self.R, self.nf)
        return U, tp.F0 @ C.T, self.hL(h).squeeze(-1)                  # psi0 (3, nz, R), log N_cr


class GBTWaveRitz(nn.Module):
    """Wave longitudinal shapes: e_r(zeta) [sin, cos](kappa_r z), R wavenumbers -> 2R candidates.
    The wavenumbers start log-spaced between half-wavelengths of 12 and 400 mm."""

    def __init__(self, n_in: int, nm: int, R: int = 6, width: int = 256, depth: int = 4, n_env: int = 5,
                 half_wave=(12.0, 400.0)):
        super().__init__()
        self.nm, self.R, self.ne = nm, R, n_env
        self.body = _mlp(n_in, width, depth)
        self.hU = nn.Linear(width, nm * 2 * R)
        self.hK = nn.Linear(width, R)
        self.hE = nn.Linear(width, R * n_env)
        self.hL = nn.Linear(width, 1)
        with torch.no_grad():
            self.hK.weight.mul_(0.01)
            self.hK.bias.copy_(torch.log(math.pi / torch.logspace(math.log10(half_wave[0]),
                                                                  math.log10(half_wave[1]), R)))
            self.hE.weight.mul_(0.01); self.hE.bias.zero_()
            self.hE.bias.view(R, n_env)[:, 0] = 1.0

    def candidates(self, tp):
        h = self.body(tp.x)
        U = self.hU(h).reshape(self.nm, 2 * self.R)
        kap = torch.exp(self.hK(h))                                    # (R,) [1/mm]
        a = self.hE(h).reshape(self.R, self.ne)
        j = torch.arange(self.ne, dtype=tp.z.dtype) * math.pi
        cz = torch.cos(tp.zeta[:, None] * j); sz = torch.sin(tp.zeta[:, None] * j)
        e0 = cz @ a.T; e1 = (-sz * j) @ a.T / tp.L; e2 = (-cz * j * j) @ a.T / tp.L ** 2   # (nz, R)
        ph = tp.z[:, None] * kap[None, :]
        s, c = torch.sin(ph), torch.cos(ph)
        pre = 1.0 / (1.0 + (kap / tp.k0) ** 2)                         # equal bending energy scaling
        ps = torch.stack([e0 * s, e1 * s + e0 * kap * c, e2 * s + 2 * e1 * kap * c - e0 * kap ** 2 * s]) * pre
        pc = torch.stack([e0 * c, e1 * c - e0 * kap * s, e2 * c - 2 * e1 * kap * s - e0 * kap ** 2 * c]) * pre
        return U, torch.cat([ps, pc], dim=2), self.hL(h).squeeze(-1)


def panel_loss(model, tp: TorchPanel, w_phys=1.0, w_teach=0.0, eps=1e-4):
    """Loss of one panel: physics (Rayleigh-Ritz over the candidates, -mu/|mu|, i.e. the
    gradient of -log mu), teacher (1 - K-cosine^2 to the lowest full-GBT modes; first mode
    weight 1, others 1/2) and the N_cr head (log N_cr = log of the network's own Ritz value)."""
    U, psi0, loglam = model.candidates(tp)
    KR, GR, psi, Uc = ritz_forms(tp, U, psi0, return_mode=True)
    mu, _ = ritz_mu(KR, GR)
    loss = w_phys * (-mu / torch.clamp(mu.detach().abs(), min=eps))
    cos2 = None
    if w_teach and tp.KT is not None:
        cos2 = teacher_cos2(tp, KR, psi, Uc)
        wk = torch.ones_like(cos2) * 0.5; wk[0] = 1.0
        loss = loss + w_teach * torch.sum(wk * (1.0 - cos2)) / torch.sum(wk)
    lam_R = 1.0 / float(mu.detach())
    if lam_R > 0:
        loss = loss + (loglam - math.log(lam_R)) ** 2
    return loss, lam_R, float(torch.exp(loglam.detach())), (None if cos2 is None else cos2.detach())


def mode_scale(ops, k0=math.pi / 30.0):
    """Jacobi scaling of the shared modes: 1 / sqrt(mean over the panels of the diagonal
    of K(k0) = K00 + k0^2 (K11 - K02 - K20) + k0^4 K22) (stiffness of a sin(k0 z) shape)."""
    acc = 0.0
    for op in ops:
        nm = op.nm; B = lambda d, e: np.diag(op.K[d * nm:(d + 1) * nm, e * nm:(e + 1) * nm])
        acc = acc + B(0, 0) + k0 ** 2 * (B(1, 1) - B(0, 2) - B(2, 0)) + k0 ** 4 * B(2, 2)
    acc = acc / len(ops)
    return 1.0 / np.sqrt(np.maximum(acc, 1e-12 * acc.max()))


def train(model, panels, steps=2000, batch=8, lr=1e-3, log_every=100, seed=0, callback=None,
          teacher_steps=0, w_teach_late=0.1):
    """Adam with cosine decay.  The first `teacher_steps` steps train the candidate
    subspace on the teacher modes (w_phys = 0); afterwards the physics loss drives the
    training and the teacher term keeps the weight w_teach_late (0 = physics only)."""
    gen = np.random.default_rng(seed)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=lr * 0.02)
    hist = []
    for it in range(1, steps + 1):
        early = it <= teacher_steps
        w_phys, w_teach = (0.0, 1.0) if early else (1.0, w_teach_late)
        idx = gen.choice(len(panels), size=min(batch, len(panels)), replace=False)
        opt.zero_grad()
        tot = 0.0; rec = []
        for i in idx:
            loss, lamR, lamH, cos2 = panel_loss(model, panels[i], w_phys=w_phys, w_teach=w_teach)
            (loss / len(idx)).backward(); tot += float(loss.detach()) / len(idx)
            rec.append((i, lamR, lamH, cos2))
        opt.step(); sched.step()
        hist.append(tot)
        if callback and (it % log_every == 0 or it == steps):
            callback(it, tot, rec)
    return hist


@torch.no_grad()
def evaluate(model, panels):
    out = []
    for tp in panels:
        U, psi0, loglam = model.candidates(tp)
        mu, _ = ritz_mu(*ritz_forms(tp, U, psi0))
        out.append(dict(sid=tp.sid, params=tp.params, lam_R=float(1.0 / mu), lam_head=float(torch.exp(loglam)),
                        lam_full=tp.lam_full, lam_shared=tp.lam_shared))
    return out
