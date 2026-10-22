"""GBT physics engine of the PINN (omega panel).

For one panel `build_omega(sample, n_conv)` runs the validated GBT code once and
returns, as torch tensors, everything the physics loss needs:

  E   (np, 8, nm, 3)     generalised strains of mode k for amplitude derivative d
  C   (np, 8, 8)         wall stiffness [A B 0; B D 0; 0 0 As] at the section points
  wp  (np,)              section quadrature weights
  G   (np, 3, 2, nm, 3)  displacement gradients g_{c,a} of mode k, derivative d
  P   (3, nz, np)        pre-buckling resultants P_zz, P_ss, P_zs (compression > 0)
  z, wz (nz,)            quadrature along the panel (the solver's own Gauss points)

with the energies (pointwise form, identical to the solver's K and K_G):
  U   = 1/2 sum_z wz sum_p wp  eps^T C eps,               eps = sum_d E^(d) q^(d)
  U_G = 1/2 sum_z wz sum_p wp  sum_c (Pzz g_cz^2 + Pss g_cs^2 + 2 Pzs g_cz g_cs)

The pre-buckling state always comes from the FULL basis (one sparse linear
solve); the buckling basis is either the full basis (n_conv=None, empty modes
removed) or the first n_conv hierarchical conventional modes with the shear
fields they activate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

from .designspace import Sample
from .oracle import omega_length, omega_panel
from .solver import configure, solver

DTYPE = torch.float64


@dataclass
class EngineData:
    E: torch.Tensor
    C: torch.Tensor
    wp: torch.Tensor
    G: torch.Tensor
    P: torch.Tensor
    z: torch.Tensor
    wz: torch.Tensor
    L: float
    bc: str
    side: str
    groups: dict            # index arrays: conventional, axial, shear, mindlin_psi, mindlin_phi
    Ry: torch.Tensor        # (r, n_conv) U_Y of the conventional modes at all nodes (SS end rows)
    Rx: torch.Tensor        # (1, n_conv) U_X of the conventional modes at s_R (periodic anchor) or empty
    lam_ref: float          # eigenvalue of the same (reduced or full) GBT model, from the solver
    lam_full: float | None  # eigenvalue of the full GBT model (if computed)
    nm: int
    info: dict


def _t(x):
    return torch.as_tensor(np.ascontiguousarray(x), dtype=DTYPE)


def energies(d: EngineData, q: torch.Tensor):
    """U and U_G for amplitude derivatives q: (nz, nm, 3) [value, 1st, 2nd derivative]."""
    eps = torch.einsum("pamd,zmd->zpa", d.E, q)
    U = 0.5 * torch.einsum("z,p,zpa,pab,zpb->", d.wz, d.wp, eps, d.C, eps)
    g = torch.einsum("pcamd,zmd->zpca", d.G, q)
    Pzz, Pss, Pzs = d.P[0], d.P[1], d.P[2]
    dens = Pzz[:, :, None] * g[..., 0] ** 2 + Pss[:, :, None] * g[..., 1] ** 2 + 2 * Pzs[:, :, None] * g[..., 0] * g[..., 1]
    UG = 0.5 * torch.einsum("z,p,zpc->", d.wz, d.wp, dens)
    return U, UG


def build_omega(sample: Sample, n_conv: int | None = 80, with_full_reference: bool = False) -> EngineData:
    p = sample.params; bc = p["bc"]; side = p["side"]
    with configure(OMEGA_SIDE_BC=side) as m:
        pan = omega_panel(sample); L = omega_length(sample)
        sec = m.omega_gbt_section(pan)
        hx = min(15.0, max(2.0, np.ptp(sec.X[:, 0]) / 40.0), L / 10.0)        # as the solver
        lb = m._LongSpline(L, m._sfs_breaks(L, hx, hx, L))
        pre = m.omega_gbt_solve(sec, L, bc, hx=hx, pre_only=True)["pre"]      # full-basis state
        if n_conv is None:
            work = sec
            kd = sum(np.abs(np.diag(A)) for A in m._gbt_Kcs(sec).values())
            gz = np.abs(sec.Gd).sum(axis=(0, 1, 2, 4))
            keep = np.nonzero(~((kd <= 1e-12 * kd.max()) & (gz == 0)))[0]  # drop empty modes
            r = m.omega_gbt_solve(sec, L, bc, hx=hx, nev=2)
            lam_ref = r["N"]
        else:
            hm = m.omega_gbt_modes(sec)
            ncv = hm["C"][1].shape[1]
            idx = list(range(min(int(n_conv), ncv)))
            work = m.truncated_section(sec, None, 0, 0, 0, hm=hm, coupled_shear=True, cidx=idx)
            keep = np.arange(work.nm)
            r = m.omega_gbt_solve(work, L, bc, hx=hx, nev=2, pre=pre)
            lam_ref = r["N"]
        lam_full = None
        if with_full_reference:
            lam_full = m.omega_gbt_buckling(pan, L, bc, nmodes=2)[0][0][0]
        nip = work.nip; o_ax = work.o_ax; o_s = work.o_s; ns = work.ns
        o_mx = getattr(work, "o_mx", o_s + ns); o_ms = getattr(work, "o_ms", o_mx); nm_all = work.nm
        U = work.U                                                       # (2N, nip)
        Ry = U[1::2, :]; Ry = Ry[np.abs(Ry).sum(axis=1) > 0]
        Rx = U[2 * work.sides[1], :][None, :] if side == "periodic" else np.zeros((0, nip))
        E = work.E[:, :, keep, :]; G = work.Gd[:, :, :, keep, :]
        pos = {k: i for i, k in enumerate(keep)}
        def grp(a, b):
            return np.array([pos[k] for k in range(a, b) if k in pos], int)
        groups = dict(conventional=grp(0, nip), axial=grp(o_ax, o_ax + 1), shear=grp(o_s, o_s + ns),
                      mindlin_psi=grp(o_mx, o_ms), mindlin_phi=grp(o_ms, nm_all))
        P = np.stack([pre[0], pre[1], pre[2]])                           # (3, nz, np)
        info = dict(hx=hx, nb=lb.nb, nspan=lb.nspan, n_conv=n_conv, nm_total=int(work.nm),
                    nm_used=len(keep), npnt=int(work.E.shape[0]), nz=int(lb.xg.size))
        return EngineData(E=_t(E), C=_t(work.Cp), wp=_t(work.wp), G=_t(G), P=_t(P), z=_t(lb.xg), wz=_t(lb.wx),
                          L=float(L), bc=bc, side=side, groups={k: torch.as_tensor(v) for k, v in groups.items()},
                          Ry=_t(Ry), Rx=_t(Rx), lam_ref=float(lam_ref), lam_full=lam_full, nm=len(keep), info=info)


def solver_mode_check(sample: Sample, n_conv: int | None = 80) -> dict:
    """Consistency test of the engine: the solver's own buckling mode, inserted into
    the engine energies, must give U/U_G = lambda of the solver."""
    d = build_omega(sample, n_conv)
    p = sample.params
    with configure(OMEGA_SIDE_BC=p["side"]) as m:
        pan = omega_panel(sample); L = omega_length(sample)
        sec = m.omega_gbt_section(pan)
        hx = d.info["hx"]
        lb = m._LongSpline(L, m._sfs_breaks(L, hx, hx, L))
        pre = m.omega_gbt_solve(sec, L, p["bc"], hx=hx, pre_only=True)["pre"]
        if n_conv is None:
            r = m.omega_gbt_solve(sec, L, p["bc"], hx=hx, nev=2, return_mode=True)
            kd = sum(np.abs(np.diag(A)) for A in m._gbt_Kcs(sec).values())
            gz = np.abs(sec.Gd).sum(axis=(0, 1, 2, 4))
            keep = np.nonzero(~((kd <= 1e-12 * kd.max()) & (gz == 0)))[0]
        else:
            hm = m.omega_gbt_modes(sec)
            work = m.truncated_section(sec, None, 0, 0, 0, hm=hm, coupled_shear=True,
                                       cidx=list(range(min(n_conv, hm["C"][1].shape[1]))))
            r = m.omega_gbt_solve(work, L, p["bc"], hx=hx, nev=2, return_mode=True, pre=pre)
            keep = np.arange(work.nm)
        Cm = r["coef"][keep]                                             # (nm, nb)
        q = np.stack([lb.B[k] @ Cm.T for k in range(3)], axis=2)        # (nz, nm, 3)
    U, UG = energies(d, _t(q))
    return dict(lam_solver=float(r["N"]), lam_engine=float(U / UG), rel=float(U / UG / r["N"] - 1))
