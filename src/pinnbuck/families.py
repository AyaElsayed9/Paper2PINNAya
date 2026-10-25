"""Family adapters: one interface to the GBT of the omega panel and of the members
(C, Z, L, T, I, box), used by `operators.build_operator` and the teacher data.

    with family_context(sample) as fc:
        fc.sec, fc.L, fc.lb              cross-section analysis, length, longitudinal splines
        fc.solve(nev, return_mode)       full GBT (validated model): {"N", "eigs", "coefs"}
        fc.pre()                         first-order pre-buckling resultants at the Gauss points
        fc.solve_reduced(w, pre, nev)    GBT restricted to a reduced section w
        fc.classes(w)                    end-condition classes in the coordinates of w

Shared bases:
  omega    POD of training buckling modes (`basis.pod_basis`), 241 modes;
  members  the cross-section of a member family is fixed in the design space (H, B, t), so
           the COMPLETE nodal GBT basis is the same for every design: nodal fields
           [U; theta] of all conventional modes (3N), all N-1 warping hats and all Mindlin
           rotations.  No truncation: the shared basis is exact (`member_basis`).
"""
from __future__ import annotations

from contextlib import contextmanager

import numpy as np

from .basis import NINT_REF, SharedBasis, hx_of
from .designspace import Sample
from .oracle import member_cfg, omega_length, omega_panel
from .solver import configure

MEMBERS = ("c", "z", "l", "t", "i", "box")


class _Omega:
    def __init__(self, m, sample):
        p = sample.params; self.m, self.bc = m, p["bc"]
        self.sec = m.omega_gbt_section(omega_panel(sample)); self.L = omega_length(sample)
        self.hx = hx_of(self.sec, self.L)
        self.lb = m._LongSpline(self.L, m._sfs_breaks(self.L, self.hx, self.hx, self.L))

    def solve(self, nev=3, return_mode=True):
        r = self.m.omega_gbt_solve(self.sec, self.L, self.bc, hx=self.hx, nev=nev, return_mode=return_mode)
        return dict(N=r["N"], eigs=r["eigs"], coefs=r.get("coefs"))

    def pre(self):
        return self.m.omega_gbt_solve(self.sec, self.L, self.bc, hx=self.hx, pre_only=True)["pre"]

    def solve_reduced(self, w, pre, nev=2):
        return self.m.omega_gbt_solve(w, self.L, self.bc, hx=self.hx, nev=nev, pre=pre)["N"]

    def classes(self, w):
        from .operators import long_classes
        return long_classes(self.m, self.sec, w, self.bc)


class _Member:
    def __init__(self, m, sample):
        p = sample.params; self.m, self.bc = m, p["bc"]
        _, fam, _, seg, hx = m._gbt_opts({})
        self.sec = m.GBTSection(member_cfg(sample), seg, fam, m.GBT_NMIN)
        if not hasattr(self.sec, "Ms"):
            self.sec.Ms = np.eye(self.sec.nms)                 # Mindlin phi: nodal rotations
        self.sec.m29 = m; self.sec.mm = m
        self.L = float(p["L"]); self.hx = min(hx, self.L / 20.0)
        zone = max(np.ptp(self.sec.X[:, 0]), np.ptp(self.sec.X[:, 1]))
        self.lb = m._LongSpline(self.L, m._sfs_breaks(self.L, self.hx, 0.5 * self.hx, zone))

    def solve(self, nev=3, return_mode=True):
        r = self.m.gbt_fsm_solve(self.sec, self.L, self.bc, method="spline", hx=self.hx, nev=nev,
                                 return_mode=return_mode)
        return dict(N=r["P"], eigs=r["eigs"], coefs=r.get("coefs"))

    def pre(self):
        return self.m.gbt_fsm_solve(self.sec, self.L, self.bc, method="spline", hx=self.hx, pre_only=True)["pre"]

    def solve_reduced(self, w, pre, nev=2):
        return self.m.gbt_fsm_solve(w, self.L, self.bc, method="spline", hx=self.hx, nev=nev, pre=pre)["P"]

    def classes(self, w):
        """Member ends (gbt_fsm_solve): SS - conventional modes held at both ends, axial
        anchor at z = 0, shear and Mindlin free; CC - conventional clamped, axial anchor,
        shear and Mindlin held at both ends."""
        from .operators import LongClass
        I = np.eye(w.nm)
        blk = dict(conv=np.arange(0, w.nip), axial=np.array([w.o_ax]), shear=np.arange(w.o_s, w.o_s + w.ns),
                   mx=np.arange(w.o_mx, w.o_mx + w.nmx), ms=np.arange(w.o_ms, w.o_ms + w.nms))
        if self.bc == "SS":
            spec = dict(conv=("zero", "zero"), axial=("zero", "free"), shear=("free", "free"),
                        mx=("free", "free"), ms=("free", "free"))
        else:
            spec = dict(conv=("clamp", "clamp"), axial=("zero", "free"), shear=("zero", "zero"),
                        mx=("zero", "zero"), ms=("zero", "zero"))
        return [LongClass(k, I[:, blk[k]], *spec[k]) for k in blk if len(blk[k])]


@contextmanager
def family_context(sample: Sample):
    if sample.family == "omega":
        with configure(OMEGA_SIDE_BC=sample.params["side"], OMEGA_GBT_NINT=list(NINT_REF)) as m:
            yield _Omega(m, sample)
    else:
        with configure() as m:
            yield _Member(m, sample)


def member_basis(family: str, reference: Sample | None = None) -> SharedBasis:
    """Complete nodal GBT basis of a member family (exact for every design of the family)."""
    from .designspace import sample_design
    ref = reference or sample_design(family, 1, seed=0)[0]
    with family_context(ref) as fc:
        sec = fc.sec
        nodal = np.vstack([sec.U, sec.Th])                     # (3N, nip), full rank
        return SharedBasis(side=family, n_conv=nodal.shape[1], nodal_C=nodal, Ts=np.eye(sec.ns),
                           Us_ref=sec.Us.copy(), lam_ref=np.zeros(0), nint=[], kind="complete",
                           nodal_mx=sec.Om.copy(), nodal_ms=sec.Ms.copy())
