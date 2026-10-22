"""GBT operators of one panel in the shared basis (input of the parametric PINN).

The GBT energies of a panel are integrals along z of quadratic forms in the amplitude
derivatives Q(z) = [q(z); q'(z); q''(z)] (3 nm values, derivative-major):

    U   = 1/2 int Q^T K Q dz                       K: (3nm, 3nm), constant along z (prismatic)
    U_G = 1/2 int sum_p w_p [P_zz g_z.g_z + P_ss g_s.g_s + 2 P_zs g_z.g_s] dz,
          g_a(z, p) = G_a(p) Q(z)                  G: (np, nc, 2, 3nm) displacement gradients

with the full-GBT pre-buckling resultants P_t(z, p) at the solver's Gauss points.
Both are exactly the solver's integrands (`verify_against_solver`: 1e-13).

For a Ritz solve the z-dependence of P is split by an SVD, P_t = sum_s alpha_ts(z)
beta_ts(p), truncated at 1e-10 (exact to round-off; `kg_svd`).

End conditions of the solver (omega_gbt_solve) as longitudinal CLASSES of modes,
each a set of directions Pi (nm x k) with a condition at z = 0 and at z = L:
  SS: conventional directions with zero U_Y at the nodes (rows R) zero at both ends;
      periodic: the U_X anchor of the right side node zero at z = 0; axial mode zero
      at z = 0; shear and Mindlin modes free.
  CC: conventional modes clamped (value and slope) at both ends; axial and shear modes
      zero at z = 0; Mindlin rotations zero at both ends.
In the shared basis these directions are the same for every panel (R acts on the
shared nodal shapes)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla

from .basis import NINT_REF, SharedBasis, hx_of, shared_section
from .designspace import Sample
from .oracle import omega_length, omega_panel
from .solver import configure

SVD_TOL = 1e-10


@dataclass
class LongClass:
    name: str
    Pi: np.ndarray           # (nm, k) directions in the shared mode coordinates
    left: str                # "free" | "zero" | "clamp"
    right: str


@dataclass
class PanelOperator:
    sid: str
    params: dict
    features: np.ndarray
    L: float
    bc: str
    side: str
    nm: int
    K: np.ndarray            # (3nm, 3nm)
    G: np.ndarray            # (np, nc, 2, 3nm)
    P: np.ndarray            # (3, nz, np) pre-buckling resultants times w_p (compression > 0)
    z: np.ndarray            # (nz,) solver Gauss points
    wz: np.ndarray           # (nz,)
    classes: list
    lam_full: float          # full GBT (all modes), the validated model
    lam_shared: float        # solver with the shared basis (spline in z): best the PINN can reach
    info: dict = field(default_factory=dict)
    teacher: np.ndarray = None       # (k, 3, nz, nm) lowest full-GBT modes in the shared coordinates
    teacher_eigs: np.ndarray = None  # their buckling factors

    def nbytes(self):
        return self.K.nbytes + self.G.nbytes + self.P.nbytes


def _stack(T):
    """(..., nm, 3) [mode, derivative] -> (..., 3nm) derivative-major."""
    return np.moveaxis(T, -1, -2).reshape(*T.shape[:-2], -1)


def section_operators(w, pre, full):
    """K, pointwise G and the w_p-weighted pre-buckling resultants of a (reduced) section."""
    nm = w.nm
    E = _stack(w.E)                                              # (np, a, 3nm)
    Y = np.matmul(w.Cp, E) * w.wp[:, None, None]
    K = E.reshape(-1, 3 * nm).T @ Y.reshape(-1, 3 * nm); K = 0.5 * (K + K.T)
    comps = (0, 1, 2) if full else (1, 2)
    G = _stack(w.Gd[:, comps])                                   # (np, nc, 2, 3nm)
    P = np.stack([np.asarray(pre[0]), np.asarray(pre[1]) if full else 0 * pre[0],
                  np.asarray(pre[2]) if full else 0 * pre[0]]) * w.wp[None, None, :]
    return K, G, P


def kg_svd(op, tol=1e-10):
    """Kg_k (3nm, 3nm) and alpha_k(z) with U_G = 1/2 sum_k int alpha_k Q^T Kg_k Q dz."""
    terms = [(0, 0, 0, 1.0), (1, 1, 1, 1.0), (2, 0, 1, 2.0)]   # (P_t, dir a, dir b, factor)
    Kg, alpha = [], []
    for t, a, b, fac in terms:
        P = op.P[t]
        if not np.any(P):
            continue
        Uz, sv, Vp = np.linalg.svd(P, full_matrices=False)
        r = int(np.sum(sv > tol * sv[0]))
        Ga = op.G[:, :, a, :]; Gb = op.G[:, :, b, :]
        for s_ in range(r):
            M = np.einsum("p,pcm,pcn->mn", Vp[s_] * sv[s_], Ga, Gb, optimize=True) * fac
            Kg.append(0.5 * (M + M.T)); alpha.append(Uz[:, s_])
    return np.array(Kg), np.array(alpha)


def long_classes(m, sec, w, bc):
    """End-condition classes in the shared coordinates of `w` (see module doc)."""
    nm = w.nm; nip = w.nip; T = w.T_basis
    Tc = T[:sec.nip, :nip]                                       # shared conv -> full conv
    I = np.eye(nm)
    blk = dict(conv=np.arange(0, nip), axial=np.array([w.o_ax]), shear=np.arange(w.o_s, w.o_s + w.ns),
               mx=np.arange(w.o_mx, w.o_mx + w.nmx), ms=np.arange(w.o_ms, w.o_ms + w.nms))
    out = []
    if bc == "SS":
        Ry = sec.U[1::2, :]; Ry = Ry[np.abs(Ry).sum(axis=1) > 0]
        R = Ry @ Tc
        Rr = sla.orth(R.T, rcond=1e-10)                          # constrained directions (nip, k)
        Z = sla.null_space(R, rcond=1e-10)                       # free directions
        out.append(LongClass("conv_fixed", I[:, blk["conv"]] @ Rr, "zero", "zero"))
        if sec.side == "periodic":
            a = (sec.U[2 * sec.sides[1], :] @ Tc) @ Z            # anchor within the free directions
            if np.linalg.norm(a) > 1e-12:
                an = a / np.linalg.norm(a)
                out.append(LongClass("conv_anchor", I[:, blk["conv"]] @ (Z @ an[:, None]), "zero", "free"))
                Z = Z @ sla.null_space(an[None, :])
        out.append(LongClass("conv_free", I[:, blk["conv"]] @ Z, "free", "free"))
        out.append(LongClass("axial", I[:, blk["axial"]], "zero", "free"))
        for k in ("shear", "mx", "ms"):
            out.append(LongClass(k, I[:, blk[k]], "free", "free"))
    elif bc == "CC":
        out.append(LongClass("conv", I[:, blk["conv"]], "clamp", "clamp"))
        out.append(LongClass("axial", I[:, blk["axial"]], "zero", "free"))
        out.append(LongClass("shear", I[:, blk["shear"]], "zero", "free"))
        for k in ("mx", "ms"):
            out.append(LongClass(k, I[:, blk[k]], "zero", "zero"))
    else:
        raise ValueError(bc)
    return [c for c in out if c.Pi.shape[1]]


def build_operator(sample: Sample, sb: SharedBasis, lam_full: float | None = None,
                   with_shared_reference: bool = True, n_teacher: int = 3) -> PanelOperator:
    """lam_full = None runs the full GBT (validated model) and stores its n_teacher lowest
    modes as `teacher`; pass a number (or nan) to skip it."""
    p = sample.params
    with configure(OMEGA_SIDE_BC=p["side"], OMEGA_GBT_NINT=list(NINT_REF)) as m:
        sec = m.omega_gbt_section(omega_panel(sample)); L = omega_length(sample); hx = hx_of(sec, L)
        lb = m._LongSpline(L, m._sfs_breaks(L, hx, hx, L))
        rfull = None
        if lam_full is None:
            rfull = m.omega_gbt_solve(sec, L, p["bc"], hx=hx, nev=n_teacher, return_mode=True)
            lam_full = rfull["N"]
        pre = m.omega_gbt_solve(sec, L, p["bc"], hx=hx, pre_only=True)["pre"]
        w = shared_section(m, sec, sb)
        teacher = teig = None
        if rfull is not None:                  # lowest modes, least-squares in the shared coordinates
            Tp = np.linalg.pinv(w.T_basis)                               # (nm, nm_full)
            teacher = np.stack([np.stack([Tp @ (lb.B[d] @ c.T).T for d in range(3)]).transpose(0, 2, 1)
                                for c in rfull["coefs"]]).astype(np.float32)  # (k, 3, nz, nm)
            teig = np.asarray(rfull["eigs"])
        full = "T" in sec.families
        K, G, P = section_operators(w, pre, full)
        classes = long_classes(m, sec, w, p["bc"])
        lam_sh = (m.omega_gbt_solve(w, L, p["bc"], hx=hx, nev=2, pre=pre)["N"]
                  if with_shared_reference else float("nan"))
    return PanelOperator(sid=sample.sid, params=dict(p), features=sample.features(), L=float(L), bc=p["bc"],
                         side=p["side"], nm=int(w.nm), K=K, G=G, P=P, z=np.asarray(lb.xg),
                         wz=np.asarray(lb.wx), classes=classes, lam_full=float(lam_full),
                         lam_shared=float(lam_sh), info=dict(hx=hx, nspan=lb.nspan, sample=sample.to_dict()),
                         teacher=teacher, teacher_eigs=teig)


# ---------------------------------------------------------------------------------
# longitudinal B-spline basis in zeta = z / L (fixed number of spans for every panel)
def zeta_basis(nspan: int, zeta: np.ndarray, L: float):
    """Uniform cubic B-splines on [0, 1]: values and first two z-derivatives at zeta.
    Returns (3, nz, nb)."""
    from scipy.interpolate import BSpline
    t = np.r_[[0.0] * 3, np.linspace(0.0, 1.0, nspan + 1), [1.0] * 3]
    nb = nspan + 3
    out = np.zeros((3, zeta.size, nb))
    for i in range(nb):
        f = BSpline(t, np.eye(nb)[i], 3, extrapolate=False)
        out[0, :, i] = f(zeta); out[1, :, i] = f.derivative(1)(zeta) / L; out[2, :, i] = f.derivative(2)(zeta) / L ** 2
    return np.nan_to_num(out)


def keep_columns(nb: int, left: str, right: str):
    drop = {"free": 0, "zero": 1, "clamp": 2}
    return np.arange(drop[left], nb - drop[right])


def bubble(zeta, left, right, L):
    """End factor b(zeta) and its first two z-derivatives: zeta^a (1 - zeta)^b with
    a, b = 0 (free), 1 (zero value), 2 (zero value and slope).  (3, nz)"""
    a = {"free": 0, "zero": 1, "clamp": 2}[left]; b = {"free": 0, "zero": 1, "clamp": 2}[right]
    P = np.polynomial.Polynomial([0, 1]) ** a * np.polynomial.Polynomial([1, -1]) ** b
    return np.stack([P(zeta), P.deriv(1)(zeta) / L, P.deriv(2)(zeta) / L ** 2])


def fourier_basis(J: int, zeta, L):
    """cos(j pi zeta), j = 0..J and sin(j pi zeta), j = 1..J with z-derivatives: (3, nz, 2J+1).
    Period 2 on [0, 1], so any smooth function on [0, 1] converges spectrally."""
    j = np.arange(J + 1); k = j * np.pi
    c, s = np.cos(np.outer(zeta, k)), np.sin(np.outer(zeta, k[1:]))
    kz = k / L; ks = kz[1:]
    return np.stack([np.hstack([c, s]), np.hstack([-kz * np.sin(np.outer(zeta, k)), ks * np.cos(np.outer(zeta, k[1:]))]),
                     np.hstack([-kz ** 2 * c, -ks ** 2 * s])])


def times_bubble(F, bb):
    """Product rule: (b f)^(d) for F (3, nz, nf), bb (3, nz)."""
    b0, b1, b2 = bb[0][:, None], bb[1][:, None], bb[2][:, None]
    return np.stack([b0 * F[0], b1 * F[0] + b0 * F[1], b2 * F[0] + 2 * b1 * F[1] + b0 * F[2]])


def orthonormalise(F, wz, L, tol=1e-12):
    """Non-redundant H2-orthonormal combination of the functions F (3, nz, nf)."""
    Gm = sum(L ** (2 * d) * (F[d] * wz[:, None]).T @ F[d] for d in range(3))
    ev, V = np.linalg.eigh(Gm); keep = ev > tol * ev.max()
    T = V[:, keep] / np.sqrt(ev[keep])
    return np.einsum("dzf,fg->dzg", F, T)


def class_bases(op, kind="fourier", n=24):
    """Longitudinal basis of every end-condition class at op.z: (3, nz, nf_c).
    The Fourier set (cos and sin) is redundant on [0, 1]; for a Ritz solve it is made
    non-redundant (H2-orthonormal) per class."""
    zeta = op.z / op.L; out = []
    for c in op.classes:
        if kind == "fourier":
            out.append(times_bubble(fourier_basis(n, zeta, op.L), bubble(zeta, c.left, c.right, op.L)))
        elif kind == "fourier_orth":
            out.append(orthonormalise(times_bubble(fourier_basis(n, zeta, op.L),
                                                   bubble(zeta, c.left, c.right, op.L)), op.wz, op.L))
        else:                                                    # (3, nz, nb) spline given as n
            out.append(n[:, :, keep_columns(n.shape[2], c.left, c.right)])
    return out


# ---------------------------------------------------------------------------------
# Ritz eigen-solve on a longitudinal basis: the best a network on that basis can reach
def ritz_matrices(op: PanelOperator, Bz):
    """Ritz K, Kg for the class-wise tensor basis Pi_c (x) basis_c.
    Bz: (3, nz, nb) spline functions with z-derivatives at op.z (end classes by dropping
    coefficients) or a list of per-class bases from `class_bases`."""
    nm = op.nm
    if isinstance(Bz, np.ndarray):
        Bz = class_bases(op, "spline", Bz)
    blocks = [(c.Pi, B) for c, B in zip(op.classes, Bz)]
    # Phi(z): (nz, 3nm, ndof) would be large; assemble block by block instead
    def form(Kmat, wfun):
        n = [Pi.shape[1] * Bc.shape[2] for Pi, Bc in blocks]; off = np.r_[0, np.cumsum(n)]
        A = np.zeros((off[-1], off[-1]))
        for i, (Pi, Bi) in enumerate(blocks):
            for j, (Pj, Bj) in enumerate(blocks):
                if j < i:
                    continue
                S = 0.0
                for d in range(3):
                    for e in range(3):
                        Kde = Kmat[d * nm:(d + 1) * nm, e * nm:(e + 1) * nm]
                        if not np.any(Kde):
                            continue
                        Kc = Pi.T @ Kde @ Pj                                    # (ki, kj)
                        Zc = (Bi[d] * wfun[:, None]).T @ Bj[e]                   # (nbi, nbj)
                        S = S + np.kron(Kc, Zc)
                A[off[i]:off[i + 1], off[j]:off[j + 1]] = S
                if j > i:
                    A[off[j]:off[j + 1], off[i]:off[i + 1]] = S.T
        return A
    Kr = form(op.K, op.wz)
    Kg, alpha = kg_svd(op)
    Gr = sum(form(Kg[k], op.wz * alpha[k]) for k in range(len(Kg)))
    return 0.5 * (Kr + Kr.T), 0.5 * (Gr + Gr.T)


def ritz_lambda(op: PanelOperator, Bz, drop=1e-15):
    """Lowest positive buckling factor of the Ritz problem (restricted to the range of K;
    use a non-redundant basis: the soft buckling directions have small K eigenvalues)."""
    Kr, Gr = ritz_matrices(op, Bz)
    ev, V = np.linalg.eigh(Kr)
    keep = ev > drop * ev.max()
    W = V[:, keep] / np.sqrt(ev[keep])                         # K-orthonormal range basis
    mu = np.linalg.eigvalsh(W.T @ Gr @ W)
    return float(1.0 / mu[-1])
