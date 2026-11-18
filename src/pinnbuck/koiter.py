"""Koiter initial post-buckling of the omega panel in the validated GBT model (18 Nov 2026).

The frozen solver (solver/aya_gbt_v32.py, omega_gbt_solve) is NOT changed.  `assemble` repeats its assembly
(stiffness K, end conditions C, end loads f, geometric stiffness Kg) by calling the solver's own helpers, and the
buckling load is checked against omega_gbt_solve.

Energy of the GBT model with the quadratic (von Karman type) membrane strains that the solver already uses for Kg:
    eps = eps_L(u) + eps_Q(u),  eps_Q,xx = 1/2 sum_c g_c,x^2,  eps_Q,ss = 1/2 sum_c g_c,s^2,  eps_Q,xs = sum_c g_c,x g_c,s
    (g_c,a = sec.Gd: gradient of displacement component c in direction a; the same comps / terms as _gbt_kg_tensors)
    Pi(u) = 1/2 u'K u + C3(u) + C4(u) - lam f'u,   C3 = int N(u) . eps_Q(u),   C4 = 1/2 int eps_Q' A eps_Q
    with N(u) = Cp[:, :3, :] eps_L(u) (membrane resultants of the linear strain, tension +) and A = Cp[:, :3, :3].
With linear pre-buckling, the second variation at lam is K - lam Kg (Kg = the solver's geometric stiffness, built
from the compression-positive first-order resultants), so the eigenproblem is exactly the solver's.

Koiter expansion, u = xi u1 + xi^2 u2, u1 the critical mode, u2 G-orthogonal to u1 (G11 = u1'Kg u1):
    lam / lam_c = 1 + a xi + b xi^2,   a = 3 C3(u1) / (lam_c G11)
    (K - lam_c Kg) u2 = -r,  r = grad C3(u1)  (bordered solve, u2' Kg u1 = 0)
    E4 = C4(u1) + 1/2 r'u2,   b = 4 E4 / (lam_c G11)        (a = 0 case)
    end-shortening  D = lam e0 + xi^2 e2,  e0 = f'K^-1 f / P,  e2 = f'u2 / P   (P = total end load per unit lam)
    stable (b > 0): post-buckling stiffness / pre-buckling stiffness = e0 / (e0 + e2 / (b lam_c))
Mode-shaped initial imperfection xi_bar (same normalisation as u1):  xi (1 - lam/lam_c) + b xi^3 = (lam/lam_c) xi_bar.
u1 is normalised to a largest in-plane nodal displacement of 1 mm, so xi and xi_bar are in mm.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch


def assemble(fc):
    """K, Kg, end reduction Cm, load f of the omega GBT model of a family context `fc` (families._Omega),
    following omega_gbt_solve line by line (SS / CC ends, ss_xfree=True, pss_scale=1)."""
    m29 = fc.sec.m29; mm = fc.sec.mm; sec = fc.sec; bc = fc.bc.upper(); L = float(fc.L)
    lb = fc.lb
    nb, nm, nip = lb.nb, sec.nm, sec.nip; ndof = nm * nb
    K = None
    Kcs = m29._gbt_Kcs(sec)
    for (d, e), A in Kcs.items():
        X = (lb.B[d].T @ sp.diags(lb.wx) @ lb.B[e]).tocsr()
        T = sp.kron(sp.csr_matrix(A), X, format="csr")
        K = T if K is None else K + T
    U = sec.U
    rows_end = {}; fixed = set()
    def fix(m, i):
        fixed.add(m * nb + i)
    ip = np.arange(nip)
    if bc == "SS":
        Ry = U[1::2, :]
        Ry = Ry[np.abs(Ry).sum(axis=1) > 0]
        r0 = Ry
        if sec.side == "periodic":
            r0 = np.vstack([Ry, U[2 * sec.sides[1], :][None, :]])
        rows_end[0] = r0; rows_end[nb - 1] = Ry
        fix(sec.o_ax, 0)
    else:
        for m in range(nip):
            for i in (0, 1, nb - 2, nb - 1):
                fix(m, i)
        fix(sec.o_ax, 0)
        for j in range(sec.ns):
            fix(sec.o_s + j, 0)
        for m in range(sec.o_mx, sec.nm):
            fix(m, 0); fix(m, nb - 1)
    kd = sum(np.abs(np.diag(A)) for A in Kcs.values())
    gz = np.abs(sec.Gd).sum(axis=(0, 1, 2, 4))
    for m in np.nonzero((kd <= 1e-12 * kd.max()) & (gz == 0))[0]:
        for i in range(nb):
            fix(int(m), i)
    rr, cc, vv = [], [], []; nq = 0; done = set()
    for i, R in rows_end.items():
        act = [m for m in ip if (m * nb + i) not in fixed]
        Rr = R[:, act]
        Z = sla.null_space(Rr) if Rr.size else np.eye(len(act))
        for jj in range(Z.shape[1]):
            for kk, m in enumerate(act):
                if abs(Z[kk, jj]) > 0:
                    rr.append(m * nb + i); cc.append(nq + jj); vv.append(Z[kk, jj])
        nq += Z.shape[1]
        done |= {m * nb + i for m in ip}
    for dof in range(ndof):
        if dof in fixed or dof in done:
            continue
        rr.append(dof); cc.append(nq); vv.append(1.0); nq += 1
    Cm = sp.csr_matrix((vv, (rr, cc)), shape=(ndof, nq))
    f = np.zeros(ndof)
    (B0, BL), (D0, DL) = lb.end[0], lb.end[1]
    w_q = sec.wp * sec.qp
    iuk = w_q @ sec.ukp; iwk = (w_q * sec.ep) @ sec.VW[:, 1, :]; iuj = w_q @ sec.ujp
    for m in range(nip):
        f[m * nb:(m + 1) * nb] += (iuk[m] - iwk[m]) * (D0 - DL)
    f[sec.o_ax * nb:(sec.o_ax + 1) * nb] += float(w_q.sum()) * (B0 - BL)
    for j in range(sec.ns):
        m = sec.o_s + j; f[m * nb:(m + 1) * nb] += iuj[j] * (B0 - BL)
    Ptot = float(w_q.sum())
    Kr = (Cm.T @ K @ Cm).tocsc(); Kr = 0.5 * (Kr + Kr.T)
    solve, release = mm._osfe_factor(Kr)
    full = "T" in sec.families
    a0 = Cm @ solve(Cm.T @ f); C = a0.reshape(nm, nb)
    amp = np.stack([(lb.B[d] @ C.T) for d in range(3)], axis=2)
    strain = np.einsum("pamd,xmd->xpa", sec.E, amp, optimize=True)
    Nres = -np.einsum("pca,xpa->xpc", sec.Cp[:, :3, :], strain, optimize=True)
    pre = (Nres[..., 0], Nres[..., 1], Nres[..., 2])
    Kg = m29._omega_kg_assemble(sec, lb, pre, full, ndof)
    Gr = (Cm.T @ Kg @ Cm).tocsc(); Gr = 0.5 * (Gr + Gr.T)
    return dict(sec=sec, lb=lb, K=K, Kg=Kg, Cm=Cm, f=f, Ptot=Ptot, Kr=Kr, Gr=Gr, solve=solve, release=release,
                a0=a0, full=full, nm=nm, nb=nb, ndof=ndof)


def buckle(S, nev=3):
    Kinv = spla.LinearOperator(S["Kr"].shape, matvec=S["solve"], dtype=float)
    mu, V = spla.eigsh(S["Gr"], k=nev, M=S["Kr"], Minv=Kinv, which="LA", ncv=max(2 * nev + 16, 20),
                       tol=1e-10, maxiter=20000)
    ok = mu > 0; lam = 1.0 / mu[ok]; V = V[:, ok]; o = np.argsort(lam)
    return lam[o], V[:, o]


class Fields:
    """Strain fields of a coefficient vector a (full, mode-major m*nb+i), in torch for gradients."""

    def __init__(self, S):
        sec, lb = S["sec"], S["lb"]
        self.nm, self.nb = S["nm"], S["nb"]
        dt = torch.float64
        self.B = [torch.as_tensor(lb.B[d].toarray(), dtype=dt) for d in range(3)]          # (nx, nb)
        self.wx = torch.as_tensor(np.asarray(lb.wx), dtype=dt)
        self.wp = torch.as_tensor(sec.wp, dtype=dt)
        comps = [0, 1, 2] if S["full"] else [1, 2]
        self.full = S["full"]
        # keep only modes with strain / gradient support (others are all zero)
        self.E = torch.as_tensor(sec.E, dtype=dt)                                              # (np, 8, nm, 3)
        self.G = torch.as_tensor(sec.Gd[:, comps], dtype=dt)                                   # (np, nc, 2, nm, 3)
        self.Cp = torch.as_tensor(sec.Cp, dtype=dt)                                            # (np, 8, 8)
        self.U = torch.as_tensor(sec.U, dtype=dt); self.nip = sec.nip

    def amp(self, a):
        C = a.reshape(self.nm, self.nb)
        return torch.stack([self.B[d] @ C.T for d in range(3)], dim=2)                       # (nx, nm, 3)

    def lin(self, A):
        return torch.einsum("pamd,xmd->xpa", self.E, A)                                        # (nx, np, 8)

    def grad(self, A):
        return torch.einsum("pcamd,xmd->xpca", self.G, A)                                      # (nx, np, nc, 2)

    @staticmethod
    def quad(g, h=None):
        """eps_Q (xx, ss, xs) of gradient g, or its symmetric bilinear form with h."""
        if h is None:
            h = g
        xx = 0.5 * (g[..., 0] * h[..., 0]).sum(-1)
        ss = 0.5 * (g[..., 1] * h[..., 1]).sum(-1)
        xs = 0.5 * ((g[..., 0] * h[..., 1]).sum(-1) + (g[..., 1] * h[..., 0]).sum(-1))
        if g is not h:
            pass
        return torch.stack([xx, ss, xs], dim=-1)                                               # (nx, np, 3)

    def W(self):
        return self.wx[:, None] * self.wp[None, :]

    def C3(self, a):
        A = self.amp(a); N = torch.einsum("pca,xpa->xpc", self.Cp[:, :3, :], self.lin(A))
        return (self.W() * (N * self.quad(self.grad(A))).sum(-1)).sum()

    def C4(self, a):
        q = self.quad(self.grad(self.amp(a)))
        return 0.5 * (self.W() * torch.einsum("xpi,pij,xpj->xp", q, self.Cp[:, :3, :3], q)).sum()

    def maxdisp(self, a):
        A = self.amp(a)[:, :self.nip, 0]                                                       # (nx, nip)
        return float((A @ self.U.T).abs().max())


def koiter(S, lam, V, mode=0):
    """Koiter coefficients of buckling mode `mode` (0 = critical)."""
    Cm, Kg, f, P = S["Cm"], S["Kg"], S["f"], S["Ptot"]
    F = Fields(S)
    v = V[:, mode]; a1 = Cm @ v
    s = 1.0 / F.maxdisp(torch.as_tensor(a1)); a1 = a1 * s; v = v * s                        # max nodal disp 1 mm
    lc = float(lam[mode])
    G11 = float(a1 @ (Kg @ a1))
    t = torch.as_tensor(a1, dtype=torch.float64).requires_grad_(True)
    c3 = F.C3(t); (r,) = torch.autograd.grad(c3, t)
    C3v = float(c3); r = r.numpy(); C4v = float(F.C4(torch.as_tensor(a1)))
    a_coef = 3.0 * C3v / (lc * G11)
    rr = Cm.T @ r
    Kt = (S["Kr"] - lc * S["Gr"]).tocsc()
    gv = S["Gr"] @ v
    Bd = sp.bmat([[Kt, sp.csc_matrix(gv[:, None])], [sp.csc_matrix(gv[None, :]), None]], format="csc")
    try:
        import pypardiso
        sol = pypardiso.spsolve(Bd, np.concatenate([-rr, [0.0]]))
    except Exception:
        sol = spla.spsolve(Bd, np.concatenate([-rr, [0.0]]))
    u2 = sol[:-1]; res = np.linalg.norm(Bd @ sol - np.concatenate([-rr, [0.0]])) / max(np.linalg.norm(rr), 1e-300)
    E4 = C4v + 0.5 * float(rr @ u2)
    b_coef = 4.0 * E4 / (lc * G11)
    e0 = float(f @ S["a0"]) / P
    e2 = float(f @ (Cm @ u2)) / P
    ratio = e0 / (e0 + e2 / (b_coef * lc)) if b_coef > 0 else float("nan")
    return dict(lam_c=lc, G11=G11, C3=C3v, C4=C4v, E4=E4, a=a_coef, b=b_coef, e0=e0, e2=e2,
                stiff_ratio=ratio, solve_res=res, u1=a1, u2=Cm @ u2, lmult=float(sol[-1]))


def imperfect_path(k, xi_bar, lam_max_ratio=2.0, n=400, e0=None, e2=None):
    """Load-shortening path of the one-mode Koiter model with a mode-shaped imperfection xi_bar [mm]:
    xi (1 - r) + a xi^2 + b xi^3 = r xi_bar,  r = lam / lam_c.  Returns lam [N/mm] and shortening
    D = lam e0 + (xi^2 + 2 xi xi_bar) e2 (load-weighted mean, or of one node when e0, e2 are given)."""
    lc, a, b = k["lam_c"], k["a"], k["b"]
    e0 = k["e0"] if e0 is None else e0; e2 = k["e2"] if e2 is None else e2
    out_l, out_d = [], []
    for xi in np.linspace(0.0, 6.0, n * 4):
        r = (xi + a * xi ** 2 + b * xi ** 3) / (xi + xi_bar) if xi + xi_bar > 0 else 0.0
        if r > lam_max_ratio:
            break
        out_l.append(r * lc); out_d.append(r * lc * e0 + (xi ** 2 + 2 * xi * xi_bar) * e2)
    return np.array(out_l), np.array(out_d)
