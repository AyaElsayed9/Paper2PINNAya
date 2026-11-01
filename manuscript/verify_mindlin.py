"""Check of the general tensor formulas added to the paper:
 (a) Xi^(de) = int E^(d)T C E^(e) ds with the Mindlin columns written as in the paper
     (E0: zeta_s -> k_ss, xi_s -> k_xs, xi -> g_xz, zeta -> g_sz; E1: xi -> k_xx, zeta -> k_xs),
     built independently from the nodal vectors Omega and the nodal zeta functions;
 (b) X^(de) = int [Sxx Gx^T Gx + Sss Gs^T Gs + Sxs (Gx^T Gs + Gs^T Gx)] ds with the shear-mode
     gradients written as in the paper (u = u_j varphi_j), for random membrane-force fields."""
import sys, numpy as np
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import omega_panel
from pinnbuck.solver import configure

for lay in ("Cross-ply", "Sym2", "Unsym"):
    s = validated_omega_sample(lay, 150.0, 1.5, "SS", "periodic")
    with configure(OMEGA_SIDE_BC="periodic") as m:
        sec = m.omega_gbt_section(omega_panel(s))
        gp, gw = np.polynomial.legendre.leggauss(4); xi = 0.5 * (gp + 1)
        nm, N = sec.nm, sec.N
        Em = np.zeros_like(sec.E)                    # paper columns for Mindlin modes only
        Gs_ = np.zeros_like(sec.Gd)                  # paper gradients for shear modes only
        p = 0
        for sg in sec.S:
            a, b, L = sg["a"], sg["b"], sg["L"]
            for q in range(4):
                h = np.array([1 - xi[q], xi[q]])
                # xi_r = t . Omega (linear), zeta_r = nodal hat (linear)
                ta = sg["t"] @ sec.Om[2*a:2*a+2]; tb = sg["t"] @ sec.Om[2*b:2*b+2]
                X = h[0]*ta + h[1]*tb; Xs = (tb - ta)/L
                MX = slice(sec.o_mx, sec.o_mx + sec.nmx)
                Em[p, 3, MX, 1] = X; Em[p, 5, MX, 0] = Xs; Em[p, 6, MX, 0] = X
                Z = h[0]*sec.Ms[a] + h[1]*sec.Ms[b]; Zs = (sec.Ms[b] - sec.Ms[a])/L
                MS = slice(sec.o_ms, sec.o_ms + sec.nms)
                Em[p, 4, MS, 0] = Zs; Em[p, 5, MS, 1] = Z; Em[p, 7, MS, 0] = Z
                J = slice(sec.o_s, sec.o_s + sec.ns)
                uj = h[0]*sec.Us[a] + h[1]*sec.Us[b]; ujs = (sec.Us[b] - sec.Us[a])/L
                Gs_[p, 0, 0, J, 1] = uj; Gs_[p, 0, 1, J, 0] = ujs
                p += 1
        Mi = np.r_[np.arange(sec.o_mx, sec.o_mx+sec.nmx), np.arange(sec.o_ms, sec.o_ms+sec.nms)]
        Jx = np.arange(sec.o_s, sec.o_s + sec.ns)
        Epap = sec.E.copy(); Epap[:, :, Mi, :] = Em[:, :, Mi, :]
        # (a) stiffness tensors, all blocks that involve a Mindlin mode
        Kc = m._gbt_Kcs(sec); worst = 0
        for d in range(3):
            for e in range(3):
                A = np.einsum("p,pai,pab,pbj->ij", sec.wp, Epap[:, :, :, d], sec.Cp, Epap[:, :, :, e], optimize=True)
                R = Kc.get((d, e), np.zeros((nm, nm)))
                sel = np.zeros((nm, nm), bool); sel[Mi, :] = True; sel[:, Mi] = True
                if np.abs(A[sel]).max() + np.abs(R[sel]).max() == 0: continue
                worst = max(worst, np.abs(A[sel] - R[sel]).max() / max(np.abs(R).max(), 1e-30))
        print(f"{lay}: Mindlin stiffness tensors, paper vs code: rel. max diff {worst:.1e}")
        # (b) geometric tensors with random S fields, shear-mode blocks
        rng = np.random.default_rng(1); npn = sec.E.shape[0]
        Pxx, Pss, Pxs = (rng.normal(size=(1, npn)) for _ in range(3))
        T = m._gbt_kg_tensors(sec, Pxx, Pss, Pxs, full=True)
        Gp = sec.Gd.copy(); Gp[:, :, :, Jx, :] = Gs_[:, :, :, Jx, :]
        worst = 0
        for d in range(3):
            for e in range(3):
                Gx_d, Gs_d = Gp[:, :, 0, :, d], Gp[:, :, 1, :, d]; Gx_e, Gs_e = Gp[:, :, 0, :, e], Gp[:, :, 1, :, e]
                w = sec.wp
                X = (np.einsum("p,pci,pcj->ij", w*Pxx[0], Gx_d, Gx_e, optimize=True) + np.einsum("p,pci,pcj->ij", w*Pss[0], Gs_d, Gs_e, optimize=True)
                     + np.einsum("p,pci,pcj->ij", w*Pxs[0], Gx_d, Gs_e, optimize=True) + np.einsum("p,pci,pcj->ij", w*Pxs[0], Gs_d, Gx_e, optimize=True))
                R = T[(d, e)][0] if (d, e) in T else np.zeros((nm, nm))
                sel = np.zeros((nm, nm), bool); sel[Jx, :] = True; sel[:, Jx] = True
                if np.abs(R).max() == 0 and np.abs(X).max() == 0: continue
                worst = max(worst, np.abs(X - R).max() / max(np.abs(R).max(), 1e-30))
                # Mindlin columns of the code gradients are zero
                assert np.abs(sec.Gd[:, :, :, Mi, :]).max() == 0
        print(f"{lay}: geometric tensors (all blocks), paper vs code: rel. max diff {worst:.1e}; Mindlin gradients zero")
