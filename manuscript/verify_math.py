"""Independent checks of the manuscript equations against the solver (V32).

1. Explicit modal matrices (B, C, D1, D2, F1, H1 and the shear blocks) written
   out in the paper == generic A^(de) = int E^(d)T C E^(e) ds of the code.
2. Classical limit: conventional modes only, uniform compression, classical SS
   ends -> B-spline Ritz eigenvalue == min_m of the sine eigenproblem (signature curve).
3. Periodicity of every mode family at the side nodes.
"""
import sys
import numpy as np
import scipy.linalg as sla

sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample          # noqa: E402
from pinnbuck.oracle import omega_panel                           # noqa: E402
from pinnbuck.solver import configure                             # noqa: E402


def get(lay, side="periodic", bs=150.0):
    s = validated_omega_sample(lay, bs, 1.5, "SS", side)
    with configure(OMEGA_SIDE_BC=side) as m:
        sec = m.omega_gbt_section(omega_panel(s))
    return m, sec


def explicit_blocks(sec):
    """Paper Eqs. for the conventional block, from the sampled profiles."""
    E, C, w = sec.E, sec.Cp, sec.wp
    I = slice(0, sec.nip)
    vs = E[:, 1, I, 0]; wss = -E[:, 4, I, 0]          # v_,s ; w_,ss
    gam = E[:, 2, I, 1]; ws = -0.5 * E[:, 5, I, 1]    # u_,s + v ; w_,s
    u = E[:, 0, I, 2]; wv = -E[:, 3, I, 2]            # u ; w
    A = C[:, :3, :3]; B = C[:, :3, 3:6]; D = C[:, 3:6, 3:6]
    def I2(c, f, g):
        return np.einsum("p,p,pk,pl->kl", w, c, f, g)
    # strain order in E: 0 e_xx, 1 e_ss, 2 g_xs | B/D index 0 xx, 1 ss, 2 xs
    Bb = I2(A[:, 1, 1], vs, vs) - I2(B[:, 1, 1], vs, wss) - I2(B[:, 1, 1], wss, vs) + I2(D[:, 1, 1], wss, wss)
    Cb = I2(A[:, 0, 0], u, u) - I2(B[:, 0, 0], u, wv) - I2(B[:, 0, 0], wv, u) + I2(D[:, 0, 0], wv, wv)
    D1 = (I2(A[:, 2, 2], gam, gam) - 2 * I2(B[:, 2, 2], gam, ws) - 2 * I2(B[:, 2, 2], ws, gam)
          + 4 * I2(D[:, 2, 2], ws, ws))
    D2 = I2(A[:, 1, 0], vs, u) - I2(B[:, 1, 0], vs, wv) - I2(B[:, 1, 0], wss, u) + I2(D[:, 1, 0], wss, wv)
    F1 = (I2(A[:, 1, 2], vs, gam) - 2 * I2(B[:, 1, 2], vs, ws) - I2(B[:, 1, 2], wss, gam)
          + 2 * I2(D[:, 1, 2], wss, ws))
    H1 = (I2(A[:, 2, 0], gam, u) - I2(B[:, 2, 0], gam, wv) - 2 * I2(B[:, 2, 0], ws, u)
          + 2 * I2(D[:, 2, 0], ws, wv))
    return {(0, 0): Bb, (2, 2): Cb, (1, 1): D1, (0, 2): D2, (0, 1): F1, (1, 2): H1}


def check_blocks(m, sec, tag):
    Kcs = m._gbt_Kcs(sec)
    ex = explicit_blocks(sec)
    I = slice(0, sec.nip)
    worst = 0.0
    for de, X in ex.items():
        R = Kcs[de][I, I]
        err = np.abs(R - X).max() / max(np.abs(R).max(), 1e-30)
        worst = max(worst, err)
        print(f"  {tag} block {de}: rel. max diff {err:.2e}")
    # shear blocks: Dcheck1 = int A66 u_j,s u_j,s ; Hcheck1 = int A16 u_j,s u_j ; Ccheck = int A11 u_j u_j
    J = slice(sec.o_s, sec.o_s + sec.ns)
    ujs = sec.E[:, 2, J, 0]; uj = sec.E[:, 0, J, 1]; w = sec.wp; A = sec.Cp[:, :3, :3]
    for de, X in {(0, 0): np.einsum("p,p,pk,pl->kl", w, A[:, 2, 2], ujs, ujs),
                  (0, 1): np.einsum("p,p,pk,pl->kl", w, A[:, 2, 0], ujs, uj),
                  (1, 1): np.einsum("p,p,pk,pl->kl", w, A[:, 0, 0], uj, uj)}.items():
        R = Kcs[de][J, J]
        err = np.abs(R - X).max() / max(np.abs(R).max(), 1e-30)
        worst = max(worst, err)
        print(f"  {tag} shear block {de}: rel. max diff {err:.2e}")
    return worst


def classical_limit(m, sec, L, M=12):
    """Ritz (B-splines, conventional modes, phi(0)=phi(L)=0) vs sine eigenproblem."""
    nip = sec.nip
    Kcs = m._gbt_Kcs(sec)
    keep = np.arange(nip)
    # uniform compression P_xx = 1 (compression positive) at every cross-section point
    npnt = sec.E.shape[0]
    # sine: K(k) = B + k^2 (D1 - D2 - D2^T) + k^4 C ; G(k) = k^2 X_VW + k^4 X_U
    Bb = Kcs[(0, 0)][np.ix_(keep, keep)]; D1 = Kcs[(1, 1)][np.ix_(keep, keep)]
    D2 = Kcs[(0, 2)][np.ix_(keep, keep)]; Cb = Kcs[(2, 2)][np.ix_(keep, keep)]
    Gd = sec.Gd; w = sec.wp
    XVW = (np.einsum("p,pk,pl->kl", w, Gd[:, 1, 0, keep, 1], Gd[:, 1, 0, keep, 1])
           + np.einsum("p,pk,pl->kl", w, Gd[:, 2, 0, keep, 1], Gd[:, 2, 0, keep, 1]))
    XU = np.einsum("p,pk,pl->kl", w, Gd[:, 0, 0, keep, 2], Gd[:, 0, 0, keep, 2])
    lam_s = []
    for mm in range(1, M + 1):
        k = mm * np.pi / L
        Kk = Bb + k * k * (D1 - D2 - D2.T) + k ** 4 * Cb
        Gk = k * k * XVW + k ** 4 * XU
        Kk = 0.5 * (Kk + Kk.T); Gk = 0.5 * (Gk + Gk.T)
        mu = sla.eigh(Gk, Kk + 1e-14 * np.abs(Kk).max() * np.eye(nip), eigvals_only=True)
        lam_s.append(1.0 / mu.max())
    lam_sine = min(lam_s); m_cr = int(np.argmin(lam_s)) + 1
    # Ritz with B-splines, same matrices
    hx = min(15.0, max(2.0, np.ptp(sec.X[:, 0]) / 40.0), L / 10.0) / 2.0
    lb = m._LongSpline(L, m._sfs_breaks(L, hx, hx, L)); nb = lb.nb
    import scipy.sparse as sp, scipy.sparse.linalg as spla
    Wx = lb.wx
    X = {(d, e): (lb.B[d].T @ sp.diags(Wx) @ lb.B[e]).tocsr() for d in range(3) for e in range(3)}
    K = None
    for (d, e), A in Kcs.items():
        T = sp.kron(sp.csr_matrix(A[np.ix_(keep, keep)]), X[(d, e)], format="csr")
        K = T if K is None else K + T
    G = sp.kron(sp.csr_matrix(XVW), X[(1, 1)]) + sp.kron(sp.csr_matrix(XU), X[(2, 2)])
    free = np.array([i for i in range(nip * nb) if (i % nb) not in (0, nb - 1)])
    Kr = K[free][:, free].tocsc(); Gr = G[free][:, free].tocsc()
    Kr = 0.5 * (Kr + Kr.T); Gr = 0.5 * (Gr + Gr.T)
    # rotation-only / empty directions: tiny regularisation of K
    pass
    lu = spla.splu(Kr.tocsc())
    Minv = spla.LinearOperator(Kr.shape, matvec=lu.solve)
    OP = spla.LinearOperator(Kr.shape, matvec=lambda x: lu.solve(Gr @ x))
    mu = spla.eigs(OP, k=6, which="LR", tol=1e-12, maxiter=50000)[0]
    lam_ritz = 1.0 / mu.real.max()
    return lam_sine, m_cr, lam_ritz


def periodicity(sec):
    sL, sR = sec.sides; N = sec.N
    dU = np.abs(sec.U[2 * sL:2 * sL + 2] - sec.U[2 * sR:2 * sR + 2]).max()
    dT = np.abs(sec.Th[sL] - sec.Th[sR]).max()
    dW = np.abs(sec.Uw[sL] - sec.Uw[sR]).max()
    dS = np.abs(sec.Us[sL] - sec.Us[sR]).max()
    dO = np.abs(sec.Om[2 * sL:2 * sL + 2] - sec.Om[2 * sR:2 * sR + 2]).max()
    dM = np.abs(sec.Ms[sL] - sec.Ms[sR]).max()
    return dict(U=dU, theta=dT, warping=dW, shear=dS, Omega=dO, b=dM)


if __name__ == "__main__":
    for lay in ("Cross-ply", "Sym2", "Unsym"):
        m, sec = get(lay)
        print(lay, "nm", sec.nm, "nip", sec.nip, "periodicity", {k: f"{v:.1e}" for k, v in periodicity(sec).items()})
        print("  worst explicit-vs-code", f"{check_blocks(m, sec, lay):.2e}")
    for lay in ("Cross-ply",):
        try:
            m, sec = get(lay)
        except KeyError:
            continue
        for L in (225.0, 450.0):
            ls, mcr, lr = classical_limit(m, sec, L)
            print(f"classical limit {lay} L={L}: sine {ls:.4f} (m={mcr}), Ritz {lr:.4f}, diff {100*(lr-ls)/ls:+.3f} %")
