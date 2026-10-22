"""A GBT mode basis shared by all panels of a family (omega).

A parametric network can only learn amplitudes if output k means the same physical
mode on every panel.  The solver's hierarchical modes do not guarantee this: they
are eigenvectors of each panel's own stiffness, so their order and sign change
with the lay-up and the geometry.

Every panel is built on the SAME node layout (`NINT_REF`: the same number of
intermediate nodes in every wall), so node i is the same point of the section for
every b_s, and the conventional GBT basis of any panel spans the same space of
nodal fields [U (2N); theta (N)].  A reduced mode is therefore defined once by its
nodal shape and is expressed exactly in each panel's own complete GBT basis; the
panel's own laminate then gives the warping, the energies and the pre-buckling
state.  Only the reduced subspace is shared, never the stiffness.

Two ways to choose the shared nodal shapes:
  * `build_shared_basis` - the first n hierarchical modes of ONE reference panel
    (quasi-isotropic, b_s = 217.75).  Kept for the record: it fails (+80-150 % on
    random lay-ups, `scripts/exp_shared_basis.py`) because the hierarchy depends on
    the laminate.
  * `pod_basis` - proper orthogonal decomposition of the buckling modes of a set of
    training panels (full GBT eigenmodes, no Abaqus): the n nodal shapes that best
    represent the critical modes over the design space.  A Galerkin projection on a
    subspace, so the eigenvalue stays an upper bound of the full-GBT one.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .designspace import FAMILIES, Sample
from .oracle import omega_length, omega_panel
from .solver import configure

QUASI_ISO_SKIN = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]
QUASI_ISO_STRINGER = [0.0, 45.0, -45.0, 90.0, -45.0, 45.0, 0.0]
NINT_REF = [6, 1, 9, 1, 12, 1, 9, 1, 6, 5, 4, 5, 5, 4, 5]   # default layout at b_s = 150 mm
ELL_THETA = 5.0          # length [mm] that scales nodal rotations in the POD (about the node spacing)
K_METRIC = np.pi / np.array([25.0, 60.0, 150.0])   # wavenumbers [1/mm] of the energy metric (local .. global)


def reference_sample(side: str, bs: float = 217.75) -> Sample:
    fx = dict(FAMILIES["omega"]["fixed"])
    return Sample("omega", dict(bs=bs, beta=1.5, bc="SS", side=side),
                  {"skin": QUASI_ISO_SKIN, "stringer": QUASI_ISO_STRINGER}, fx)


@dataclass
class SharedBasis:
    side: str
    n_conv: int
    nodal_C: np.ndarray      # (3N, n_conv): [U; theta] of the shared conventional modes
    Ts: np.ndarray           # (ns, k_s): shared shear fields in the (geometric) hat coordinates
    Us_ref: np.ndarray       # hat basis of the layout (layout check)
    lam_ref: np.ndarray      # reference eigenvalues / POD singular values
    nint: list
    kind: str = "reference"
    nodal_mx: np.ndarray = None   # (2N, k_mx): shared Mindlin rotation-vector fields (nodal)
    nodal_ms: np.ndarray = None   # (N, k_ms)


def hx_of(sec, L):
    return min(15.0, max(2.0, np.ptp(sec.X[:, 0]) / 40.0), L / 10.0)        # as the solver


def build_shared_basis(side: str, n_conv: int, bs_ref: float = 217.75, nint=NINT_REF) -> SharedBasis:
    ref = reference_sample(side, bs_ref)
    with configure(OMEGA_SIDE_BC=side, OMEGA_GBT_NINT=list(nint)) as m:
        sec = m.omega_gbt_section(omega_panel(ref))
        hm = m.omega_gbt_modes(sec)
        n = min(int(n_conv), hm["C"][1].shape[1])
        work = m.truncated_section(sec, None, 0, 0, 0, hm=hm, coupled_shear=True, cidx=list(range(n)))
        Tc = hm["C"][1][:, :n]
        nodal = np.vstack([sec.U @ Tc, sec.Th @ Tc])
        Ts = work.T_basis[sec.o_s:sec.o_s + sec.ns, work.o_s:work.o_s + work.ns]
        return SharedBasis(side, n, nodal, Ts, sec.Us.copy(), hm["C"][0][:n], list(nint))


def nodal_metric(m, sec):
    """Energy metric of the conventional nodal fields [U; ELL_THETA*theta] of a panel:
    sum over K_METRIC of the cross-section stiffness K(k) = K00 + k^2 (K11 - K02 - K20)
    + k^4 K22 of a sin(k z) shape, each normalised by its trace, mapped to nodal fields."""
    A = m._gbt_Kcs(sec); ip = np.arange(sec.nip)
    Z = lambda d, e: A[(d, e)][np.ix_(ip, ip)] if (d, e) in A else 0.0
    Kc = 0.0
    for k in K_METRIC:
        Kk = Z(0, 0) + k ** 2 * (Z(1, 1) - Z(0, 2) - Z(2, 0)) + k ** 4 * Z(2, 2)
        Kk = 0.5 * (Kk + Kk.T); Kc = Kc + Kk / np.trace(Kk)
    S3p = np.linalg.pinv(np.vstack([sec.U, ELL_THETA * sec.Th]))      # (nip, 3N)
    return S3p.T @ Kc @ S3p


def panel_snapshots(sample: Sample, nev: int = 3, nint=NINT_REF) -> dict:
    """Full-GBT buckling modes of one panel as nodal fields (one column per spline
    coefficient): conventional [U; ELL_THETA*theta] and shear-hat amplitudes.  Each
    mode is scaled to unit Frobenius norm of its conventional field.  The Mindlin
    blocks are stored as NODAL rotation fields (Om c, Ms c): the number of Mindlin
    coordinates varies by panel (empty round-off vectors normal to a web), the nodal
    fields do not."""
    p = sample.params
    with configure(OMEGA_SIDE_BC=p["side"], OMEGA_GBT_NINT=list(nint)) as m:
        sec = m.omega_gbt_section(omega_panel(sample)); L = omega_length(sample)
        r = m.omega_gbt_solve(sec, L, p["bc"], hx=hx_of(sec, L), nev=nev, return_mode=True)
        S3 = np.vstack([sec.U, ELL_THETA * sec.Th])
        F, H, X, S = [], [], [], []
        for c in r["coefs"]:
            f = S3 @ c[:sec.nip]
            s = np.linalg.norm(f) or 1.0
            F.append(f / s); H.append(c[sec.o_s:sec.o_s + sec.ns] / s)
            X.append(sec.Om @ c[sec.o_mx:sec.o_mx + sec.nmx] / s)          # nodal rotation vectors (2N)
            S.append(sec.Ms @ c[sec.o_ms:sec.o_ms + sec.nms] / s)          # nodal rotations (N)
        return dict(sid=sample.sid, eigs=r["eigs"], F=np.hstack(F), H=np.hstack(H), X=np.hstack(X),
                    S=np.hstack(S), Us=sec.Us.copy(), Om=sec.Om.copy(), Ms=sec.Ms.copy(), M=nodal_metric(m, sec))


def _pod(snaps, key, k):
    A = np.hstack([s[key] for s in snaps])
    if A.shape[0] == 0 or k == 0:
        return np.zeros((A.shape[0], 0)), np.zeros(0)
    U, sv, _ = np.linalg.svd(A, full_matrices=False)
    k = min(int(k), int(np.sum(sv > 1e-12 * sv[0])))
    return U[:, :k].copy(), sv[:k]


def pod_basis(snaps: list[dict], side: str, n_conv: int, n_shear: int | None = None,
              n_mindlin: int = 0, nint=NINT_REF) -> SharedBasis:
    """POD of the snapshot fields, block by block (conventional, shear hats, Mindlin
    psi and phi rotations); n_shear = None uses the same count as n_conv; n_mindlin
    modes are taken in each of the two Mindlin blocks."""
    if "M" in snaps[0]:                                                 # energy-weighted POD
        Mb = sum(sn["M"] for sn in snaps) / len(snaps); Mb = 0.5 * (Mb + Mb.T)
        ev, V = np.linalg.eigh(Mb); ev = np.maximum(ev, 0.0) + 1e-6 * ev.max() / Mb.shape[0]
        Wh, Wi = (V * np.sqrt(ev)) @ V.T, (V / np.sqrt(ev)) @ V.T
        A = Wh @ np.hstack([sn["F"] for sn in snaps])
        U_, sc, _ = np.linalg.svd(A, full_matrices=False)
        k = min(int(n_conv), int(np.sum(sc > 1e-12 * sc[0])))
        Uc, sc = Wi @ U_[:, :k], sc[:k]
    else:
        Uc, sc = _pod(snaps, "F", n_conv)
    Uh, _ = _pod(snaps, "H", n_conv if n_shear is None else n_shear)
    Ux, _ = _pod(snaps, "X", n_mindlin) if "X" in snaps[0] else (None, None)
    Um, _ = _pod(snaps, "S", n_mindlin) if "S" in snaps[0] else (None, None)
    nodal = Uc.copy(); N = nodal.shape[0] // 3
    nodal[2 * N:] /= ELL_THETA                                          # back to rotations
    return SharedBasis(side, nodal.shape[1], nodal, Uh, snaps[0]["Us"], sc, list(nint), kind="pod",
                       nodal_mx=Ux, nodal_ms=Um)


def pod_energy(snaps: list[dict], key: str = "F"):
    """Cumulative captured share of the snapshot 'energy' (squared singular values;
    conventional block in the energy metric when available)."""
    A = np.hstack([s[key] for s in snaps])
    if key == "F" and "M" in snaps[0]:
        Mb = sum(sn["M"] for sn in snaps) / len(snaps); ev, V = np.linalg.eigh(0.5 * (Mb + Mb.T))
        A = (V * np.sqrt(np.maximum(ev, 0) + 1e-6 * ev.max() / Mb.shape[0])) @ V.T @ A
    sv = np.linalg.svd(A, compute_uv=False)
    return np.cumsum(sv ** 2) / np.sum(sv ** 2)


def shared_section(m, sec, sb: SharedBasis, own_shear: bool = False):
    """Reduced copy of `sec` (a panel built with the same node layout) whose
    conventional modes have the nodal shapes of `sb`.  Shear modes: the shared hat
    fields of `sb` (own_shear=False) or the panel's own condensed warping of the
    shared conventional modes (own_shear=True; panel-specific, for comparison)."""
    if sec.Us.shape != sb.Us_ref.shape or not np.allclose(sec.Us, sb.Us_ref):
        raise ValueError("node layout / side condition differ from the shared basis")
    S3 = np.vstack([sec.U, sec.Th])                                  # (3N, nip) complete basis
    Tc, *_ = np.linalg.lstsq(S3, sb.nodal_C, rcond=None)
    resid = np.linalg.norm(S3 @ Tc - sb.nodal_C) / np.linalg.norm(sb.nodal_C)
    if resid > 1e-8:
        raise ValueError(f"shared nodal shapes not representable (residual {resid:.2e})")
    if own_shear:
        hm = dict(m.omega_gbt_modes(sec)); hm["C"] = (np.zeros(Tc.shape[1]), Tc)
        return m.truncated_section(sec, None, 0, 0, 0, hm=hm, coupled_shear=True, cidx=list(range(Tc.shape[1])))
    def onto(B, F, name):
        if F is None or F.shape[1] == 0:
            return np.zeros((B.shape[1], 0))
        T, *_ = np.linalg.lstsq(B, F, rcond=None)
        res = np.linalg.norm(B @ T - F) / np.linalg.norm(F)
        if res > 1e-8:
            raise ValueError(f"shared {name} fields not representable (residual {res:.2e})")
        return T
    Tmx = onto(sec.Om, sb.nodal_mx, "Mindlin psi"); Tms = onto(sec.Ms, sb.nodal_ms, "Mindlin phi")
    hm = {"C": (np.zeros(Tc.shape[1]), Tc), "S": (np.zeros(sb.Ts.shape[1]), sb.Ts),
          "MX": (np.zeros(Tmx.shape[1]), Tmx), "MS": (np.zeros(Tms.shape[1]), Tms)}
    return m.truncated_section(sec, None, None, None, None, hm=hm)
