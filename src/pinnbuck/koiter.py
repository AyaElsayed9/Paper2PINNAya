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
    (K - lam_c Kg) u2 = -r,  r = grad C3(u1)  (constrained PCG, u2' Kg u1 = 0)
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
    Kcs = m29._gbt_Kcs(sec); Kblk = []
    for (d, e), A in Kcs.items():
        X = (lb.B[d].T @ sp.diags(lb.wx) @ lb.B[e]).tocsr()
        Kblk.append((np.asarray(A), X))
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
    del K, Kg                                    # 19 Nov: full matrices dropped (memory); K a via the blocks

    def Kmul(a):
        C = a.reshape(nm, nb)
        return sum(A @ C @ X.T for A, X in Kblk).ravel()
    return dict(sec=sec, lb=lb, Kmul=Kmul, Cm=Cm, f=f, Ptot=Ptot, Kr=Kr, Gr=Gr, solve=solve, release=release,
                a0=a0, full=full, nm=nm, nb=nb, ndof=ndof)


def constrained_solve(S, lam_c, Vr, rhs, tol=1e-10, maxit=2000):
    """w with (Kr - lam_c Gr) w = rhs - Gr Vr mu,  Vr' Gr w = 0  (Vr: buckling modes, reduced coords, columns).
    PCG with the Kr factor as preconditioner: on the Gr-complement of the modes K - lam_c G is positive definite,
    and both the operator and the preconditioner keep the iterates there (19 Nov, replaces the bordered solve)."""
    Kr, Gr, Ksolve = S["Kr"], S["Gr"], S["solve"]
    GV = Gr @ Vr; VGV = Vr.T @ GV
    def projR(y):                                 # remove the GV part: V' y = 0
        return y - GV @ np.linalg.solve(VGV, Vr.T @ y)
    def projC(x):                                 # Gr-orthogonal complement of the modes
        return x - Vr @ np.linalg.solve(VGV, GV.T @ x)
    b = projR(rhs); x = np.zeros_like(b); r = b.copy(); z = projC(Ksolve(r)); p_ = z.copy(); rz = r @ z
    nb_ = np.linalg.norm(b)
    for it in range(maxit):
        Ap = Kr @ p_ - lam_c * (Gr @ p_)
        al = rz / (p_ @ Ap); x = x + al * p_; r = projR(r - al * Ap)
        if np.linalg.norm(r) < tol * nb_:
            return projC(x), it + 1
        z = projC(Ksolve(r)); rzn = r @ z; p_ = z + (rzn / rz) * p_; rz = rzn
    return projC(x), maxit


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
    Cm, f, P = S["Cm"], S["f"], S["Ptot"]
    F = Fields(S)
    v = V[:, mode]; a1 = Cm @ v
    s = 1.0 / F.maxdisp(torch.as_tensor(a1)); a1 = a1 * s; v = v * s                        # max nodal disp 1 mm
    lc = float(lam[mode])
    G11 = float(v @ (S["Gr"] @ v))
    t = torch.as_tensor(a1, dtype=torch.float64).requires_grad_(True)
    c3 = F.C3(t); (r,) = torch.autograd.grad(c3, t)
    C3v = float(c3); r = r.numpy(); C4v = float(F.C4(torch.as_tensor(a1)))
    a_coef = 3.0 * C3v / (lc * G11)
    rr = Cm.T @ r
    u2, res = constrained_solve(S, lc, v[:, None], -rr)
    E4 = C4v + 0.5 * float(rr @ u2)
    b_coef = 4.0 * E4 / (lc * G11)
    e0 = float(f @ S["a0"]) / P
    e2 = float(f @ (Cm @ u2)) / P
    ratio = e0 / (e0 + e2 / (b_coef * lc)) if b_coef > 0 else float("nan")
    return dict(lam_c=lc, G11=G11, C3=C3v, C4=C4v, E4=E4, a=a_coef, b=b_coef, e0=e0, e2=e2,
                stiff_ratio=ratio, cg_iters=res, u1=a1, u2=Cm @ u2)


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


# ---------------------------------------------------------------------------------------------------------------
# Multi-mode (interaction) Koiter, 19 Nov 2026.  Modes u_i (i < M) with loads lam_i close to lam_1:
#   Pi(xi) = sum_i 1/2 G_ii (lam_i - lam) xi_i^2 - lam sum_i G_ii xibar_i xi_i + sum_ijk c_ijk xi_i xi_j xi_k
#            + sum_ijkl B_ijkl xi_i xi_j xi_k xi_l
#   c_ijk = c(u_i, u_j, u_k) (C3(u) = c(u, u, u));  second-order fields w_ij: (K - lam_1 Kg) w_ij = -r_ij,
#   r_ij = 3 c(u_i, u_j, .),  w_ij Kg-orthogonal to every u_k;  B_ijkl = Q_ijkl + 1/2 r_ij' w_kl,
#   Q_ijkl = 1/2 int q(u_i,u_j)' A q(u_k,u_l)  (q = symmetric bilinear form of eps_Q).
#   Shortening D = lam e0 + sum_ij [(xi_i + xibar_i)(xi_j + xibar_j) - xibar_i xibar_j] e_ij,  e_ij = f'w_ij / P.
# ---------------------------------------------------------------------------------------------------------------
def koiter_multi(S, lam, V, M):
    Cm, f, P = S["Cm"], S["f"], S["Ptot"]
    F = Fields(S); W = F.W(); dt = torch.float64
    U = []; Vr = []
    for i in range(M):
        a = Cm @ V[:, i]; s = 1.0 / F.maxdisp(torch.as_tensor(a)); U.append(a * s); Vr.append(V[:, i] * s)
    Gd = np.array([float(Vr[i] @ (S["Gr"] @ Vr[i])) for i in range(M)])
    T = [torch.as_tensor(u, dtype=dt) for u in U]
    A = [F.amp(t) for t in T]; g = [F.grad(x) for x in A]
    Nl = [torch.einsum("pca,xpa->xpc", F.Cp[:, :3, :], F.lin(x)) for x in A]
    q = {(i, j): F.quad(g[i], g[j]) for i in range(M) for j in range(M)}
    c = np.zeros((M, M, M)); Q = np.zeros((M, M, M, M))
    for i in range(M):
        for j in range(M):
            for k in range(M):
                c[i, j, k] = float((W * ((Nl[i] * q[j, k]).sum(-1) + (Nl[j] * q[i, k]).sum(-1)
                                          + (Nl[k] * q[i, j]).sum(-1))).sum()) / 3.0
                for l in range(M):
                    Q[i, j, k, l] = 0.5 * float((W * torch.einsum("xpi,pij,xpj->xp", q[i, j], F.Cp[:, :3, :3], q[k, l])).sum())
    Vm = np.stack(Vr, 1); its = []
    r = {}; w = {}
    for i in range(M):
        for j in range(i, M):
            d = torch.zeros(len(U[0]), dtype=dt, requires_grad=True)
            Ad = F.amp(d); gd = F.grad(Ad); Nd = torch.einsum("pca,xpa->xpc", F.Cp[:, :3, :], F.lin(Ad))
            val = (W * ((Nd * q[i, j]).sum(-1) + (Nl[i] * F.quad(g[j], gd)).sum(-1)
                        + (Nl[j] * F.quad(g[i], gd)).sum(-1))).sum()
            # 3 c(u_i,u_j,d) = int [N(d) q_ij + N_i q(u_j,d) + N_j q(u_i,d)],  q = quad = symmetric bilinear eps_Q
            (rv,) = torch.autograd.grad(val, d)
            rr = Cm.T @ rv.numpy()
            sol, it = constrained_solve(S, float(lam[0]), Vm, -rr); its.append(it)
            del d, Ad, gd, Nd, val, rv
            r[i, j] = r[j, i] = rr; w[i, j] = w[j, i] = sol
    B = np.zeros((M, M, M, M))
    for i in range(M):
        for j in range(M):
            for k in range(M):
                for l in range(M):
                    B[i, j, k, l] = Q[i, j, k, l] + 0.5 * float(r[i, j] @ w[k, l])
    import itertools
    B = sum(np.transpose(B, p) for p in itertools.permutations(range(4))) / 24.0     # fully symmetric quartic form
    e0 = float(f @ S["a0"]) / P
    E = np.array([[float(f @ (Cm @ w[i, j])) / P for j in range(M)] for i in range(M)])
    return dict(M=M, lam=np.asarray(lam[:M], float), G=Gd, c=c, Q=Q, B=B, e0=e0, E=E, cg_iters=its,
                U=U, W={k: Cm @ v for k, v in w.items()})


def multi_path(km, xibar, lam_max_ratio=1.9, ds=0.01, nmax=20000, e0=None, E=None):
    """Equilibrium path of the reduced M-mode model with imperfection xibar (vector, mm) by pseudo-arc-length
    continuation in (xi, lam).  Returns lam, xi (n, M), shortening D."""
    lamv, Gd, c, B = km["lam"], km["G"], km["c"], km["B"]
    e0 = km["e0"] if e0 is None else e0; E = km["E"] if E is None else E
    M = len(lamv); xb = np.asarray(xibar, float); sc = lamv[0]
    def R(x, l):          # dPi/dxi, scaled by lam_1 G_ii
        return (Gd * (lamv - l) * x - l * Gd * xb + 3 * np.einsum("ijk,j,k->i", c, x, x)
                + 4 * np.einsum("ijkl,j,k,l->i", B, x, x, x)) / (sc * Gd)
    def J(x, l):
        Jx = (np.diag(Gd * (lamv - l)) + 6 * np.einsum("ijk,k->ij", c, x) + 12 * np.einsum("ijkl,k,l->ij", B, x, x)) / (sc * Gd)[:, None]
        Jl = (-Gd * x - Gd * xb) / (sc * Gd)
        return Jx, Jl
    x = np.zeros(M); l = 0.0
    # predictor direction: start along increasing lam
    t = np.zeros(M + 1); t[-1] = 1.0
    out = [(l, x.copy())]
    for _ in range(nmax):
        Jx, Jl = J(x, l); Ja = np.hstack([Jx, Jl[:, None]])
        # tangent: null vector of Ja, oriented like the previous tangent
        _, _, Vt = np.linalg.svd(np.vstack([Ja, np.zeros((1, M + 1))]))
        tn = Vt[-1]; tn = tn if tn @ t >= 0 else -tn; t = tn
        y = np.concatenate([x, [l / sc]]) + ds * t
        for _ in range(30):
            xx, ll = y[:-1], y[-1] * sc
            Rv = R(xx, ll); Jx, Jl = J(xx, ll)
            Ja = np.vstack([np.hstack([Jx, (Jl * sc)[:, None]]), t[None, :]])
            rhs = np.concatenate([Rv, [t @ (y - np.concatenate([x, [l / sc]])) - ds]])
            dy = np.linalg.solve(Ja, -rhs); y = y + dy
            if np.abs(dy).max() < 1e-12:
                break
        x, l = y[:-1], y[-1] * sc
        out.append((l, x.copy()))
        if l > lam_max_ratio * sc or np.abs(x).max() > 8:
            break
    L = np.array([o[0] for o in out]); X = np.array([o[1] for o in out])
    Xt = X + xb[None, :]
    D = L * e0 + np.einsum("ni,ij,nj->n", Xt, E, Xt) - float(xb @ E @ xb)
    return L, X, D
