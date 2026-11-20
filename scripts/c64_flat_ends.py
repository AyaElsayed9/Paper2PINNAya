"""20 Nov 2026, comment C64, exact version of the uniform-end-strain check: the loaded ends are kept PLANE (the
warping at each end is restricted to a + b Y + c Z, i.e. a rigid end plate that may rotate), the total end force is
applied to the axial mode, and the end forces distribute themselves according to the wall stiffnesses.  This is what
an Abaqus model with the end edges coupled to a reference point (load on the RP) does.  Solver unchanged: the
assembly of koiter.assemble + a penalty on the end constraints; the pre-buckling state and Kg are recomputed.
Compared quantity: critical mean axial force per unit width P_cr / W, against the paper's equal line loads.
Variants: 'plane' (end plate may rotate, default) and 'flat' (end plate may only translate); for each, the effect of
the load distribution alone (buckling with the paper's end conditions) and with the plane end also in the buckling mode.
    python scripts/c64_flat_ends.py     -> results/c64_flat_ends.csv"""
import csv, sys
from pathlib import Path
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.families import family_context
from pinnbuck import koiter as KO


def end_warping_rows(S, which):
    """sparse (N x ndof): nodal warping u(s) at the end x = 0 (which 0) or x = L (which 1)"""
    sec, lb, nb, nm = S["sec"], S["lb"], S["nb"], S["nm"]
    Bv = lb.end[0][which]; Dv = lb.end[1][which]; N = sec.Uw.shape[0]
    M = np.zeros((N, nm, nb))
    M[:, :sec.nip, :] += sec.Uw[:, :, None] * Dv[None, None, :]
    M[:, sec.o_ax, :] += Bv[None, :]
    M[:, sec.o_s:sec.o_s + sec.ns, :] += sec.Us[:, :, None] * Bv[None, None, :]
    return sp.csr_matrix(M.reshape(N, nm * nb))


def solve_variant(fc, S, mode, kappa_rel=1e4):
    sec, Cm = S["sec"], S["Cm"]; X = sec.X; N = X.shape[0]
    basis = [np.ones(N)] + ([X[:, 0], X[:, 1]] if mode == "plane" else [])
    Q, _ = np.linalg.qr(np.stack(basis, 1)); Pp = np.eye(N) - Q @ Q.T            # removes a + bY + cZ
    A = sp.vstack([sp.csr_matrix(Pp) @ end_warping_rows(S, 0), sp.csr_matrix(Pp) @ end_warping_rows(S, 1)]).tocsr()
    Ar = (A @ Cm).tocsr()
    kap = kappa_rel * abs(S["Kr"].diagonal()).max()
    Kp = (S["Kr"] + kap * (Ar.T @ Ar)).tocsc()
    lb = S["lb"]; nb = S["nb"]; f = np.zeros(S["ndof"])
    (B0, BL) = lb.end[0]; Ptot = S["Ptot"]
    f[sec.o_ax * nb:(sec.o_ax + 1) * nb] = Ptot * (B0 - BL)                    # same total force as the paper's load
    fsolve, frel = sec.mm._osfe_factor(Kp); a0 = Cm @ fsolve(Cm.T @ f)
    C = a0.reshape(S["nm"], nb)
    amp = np.stack([(lb.B[d] @ C.T) for d in range(3)], axis=2)
    strain = np.einsum("pamd,xmd->xpa", sec.E, amp, optimize=True)
    Nres = -np.einsum("pca,xpa->xpc", sec.Cp[:, :3, :], strain, optimize=True)
    Kg = sec.m29._omega_kg_assemble(sec, lb, (Nres[..., 0], Nres[..., 1], Nres[..., 2]), S["full"], S["ndof"])
    Gr = (Cm.T @ Kg @ Cm).tocsc(); Gr = 0.5 * (Gr + Gr.T)
    Kinv = spla.LinearOperator(Kp.shape, matvec=fsolve, dtype=float)
    mu = spla.eigsh(Gr, k=2, M=Kp, Minv=Kinv, which="LA", tol=1e-10, return_eigenvectors=False)
    lam_rig = 1.0 / mu.max()                       # plane ends in pre-buckling AND in the buckling mode
    Kinv0 = spla.LinearOperator(Kp.shape, matvec=S["solve"], dtype=float)
    mu0 = spla.eigsh(Gr, k=2, M=S["Kr"], Minv=Kinv0, which="LA", tol=1e-10, return_eigenvectors=False)
    lam = 1.0 / mu0.max()                          # only the load distribution of a plane end (buckling ends as paper)
    frel()
    # end-force shares per wall part (pre-buckling resultant n_x at the first Gauss line)
    part = np.array([sec.wpart[sg["w"]] for sg in sec.S for _ in range(4)])
    nx0 = Nres[Nres.shape[0] // 2, :, 0]; share = {p: float((sec.wp * nx0)[part == p].sum() / (sec.wp * nx0).sum()) for p in ("skin", "foot", "web", "crown")}
    return lam, share, lam_rig


def main():
    rows = []
    for lay in ("Cross-ply", "Sym1", "Sym2", "Unsym"):
        for bs in (150.0, 285.5):
            for beta in (1.5,):
                s = validated_omega_sample(lay, bs, beta, "SS", "periodic")
                with family_context(s) as fc:
                    S = KO.assemble(fc); lam0, _ = KO.buckle(S, 2); W = 2 * bs; P = S["Ptot"]
                    res = {}
                    for mode in ("plane", "flat"):
                        lam, share, viol = solve_variant(fc, S, mode); res[mode] = (lam, share, viol)
                    lam_k10 = solve_variant(fc, S, "plane", kappa_rel=1e5)[0] if (lay, bs) == ("Cross-ply", 150.0) else float("nan")
                    sh0 = None
                    S["release"]()
                row = dict(layup=lay, bs=bs, beta=beta, Pcr_W_equal_load=float(lam0[0]) * P / W, P_total_per_W=P / W)
                for mode in ("plane", "flat"):
                    lam, share, lam_rig = res[mode]
                    row[f"diff_load_{mode}_pct"] = 100 * (lam / lam0[0] - 1)            # distribution effect only
                    row[f"diff_rigid_{mode}_pct"] = 100 * (lam_rig / lam0[0] - 1)       # + plane ends in the buckling mode
                    row[f"skin_share_{mode}"] = share["skin"]; row[f"web_share_{mode}"] = share["web"]
                    row[f"Pcr_W_rigid_{mode}"] = lam_rig * P / W                        # = Abaqus RP eigenvalue (flat)
                row["kappa_check_pct"] = 100 * (lam_k10 / res["plane"][0] - 1)
                rows.append(row); print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    with open(ROOT / "results/c64_flat_ends.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
