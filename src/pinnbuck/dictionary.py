"""Longitudinal dictionary of the GBT buckling modes: bulk waves + end layers.

The buckling modes of the omega panels are of two kinds (`docs/ROADMAP.md`, 2c):
waves in the bulk of the panel and modes localised at the loaded ends, where the load
introduction concentrates the pre-buckling stress.  Both are described in PHYSICAL
coordinates, independent of the panel length L:

  * bulk waves   sin(kappa_w z), cos(kappa_w z)   for fixed half-wavelengths (25 .. 350 mm),
                 times the end factor b_c of every end-condition class (`operators.bubble`);
  * end layers   cubic B-splines in the distance from each end over a length l_e = b_s,
                 the first one or two dropped for classes with zero value / clamped ends.

A mode is q(z) = sum_c Pi_c Pi_c^T sum_f u_f(p) phi_f^c(z): one cross-section vector
u_f per dictionary function f.  For the Rayleigh-Ritz problem the vectors u_f only have
to span the right space, the amplitudes and phases follow exactly from the Ritz step.
With exact cross-section vectors (all nm per function) this dictionary reproduces the
full GBT within 0.05-0.27 % on the hardest (end-localised) panels.
"""
from __future__ import annotations

import math

import numpy as np
import torch
from scipy.interpolate import BSpline

from .operators import PanelOperator, bubble, orthonormalise

HALF_WAVES = tuple(np.geomspace(25.0, 350.0, 6))   # [mm]
N_END = 8                                          # end-layer splines per end
END_LENGTH = 1.0                                   # l_e / b_s
_DROP = {"free": 0, "zero": 1, "clamp": 2}


def end_splines(z, L, le, ne, side):
    """(3, nz, ne): cubic B-splines in the distance s from the end (s = z or L - z);
    function 0 has value 1 at the end, function 1 the slope; z-derivatives."""
    s = z if side == "L" else L - z
    h = le / ne
    t = np.r_[0.0, 0.0, 0.0, np.arange(0, ne + 4) * h]
    out = []
    for i in range(ne):
        c = np.zeros(len(t) - 4); c[i] = 1.0
        b = BSpline(t, c, 3, extrapolate=False)
        v, d1, d2 = (np.nan_to_num(b(s)), np.nan_to_num(b.derivative(1)(s)), np.nan_to_num(b.derivative(2)(s)))
        if side == "R":
            d1 = -d1
        out.append(np.stack([v, d1, d2]))
    return np.stack(out, axis=2)


def waves(z, half_waves):
    """(3, nz, 2 nw): sin and cos of the physical wavenumbers pi / half_wave."""
    out = []
    for hw in half_waves:
        k = math.pi / hw
        out.append(np.stack([np.sin(k * z), k * np.cos(k * z), -k * k * np.sin(k * z)]))
        out.append(np.stack([np.cos(k * z), -k * np.sin(k * z), -k * k * np.cos(k * z)]))
    return np.stack(out, axis=2)


def dictionary(op: PanelOperator, half_waves=HALF_WAVES, n_end=N_END, end_length=END_LENGTH):
    """Per class c: functions (3, nz, nD) of the common dictionary (zero where a function
    is not admissible for the class), and labels of the nD functions."""
    z, L, bs = op.z, op.L, op.params["bs"]
    le = min(end_length * bs, 0.5 * L)
    W = waves(z, half_waves)
    labels = [f"wave{hw:.0f}{sc}" for hw in half_waves for sc in ("s", "c")]
    labels += [f"endL{i}" for i in range(n_end)] + [f"endR{i}" for i in range(n_end)]
    EL, ER = end_splines(z, L, le, n_end, "L"), end_splines(z, L, le, n_end, "R")
    out = []
    for c in op.classes:
        bb = bubble(z / L, c.left, c.right, L)
        b0, b1, b2 = bb[0][:, None], bb[1][:, None], bb[2][:, None]
        Wb = np.stack([b0 * W[0], b1 * W[0] + b0 * W[1], b2 * W[0] + 2 * b1 * W[1] + b0 * W[2]])
        El, Er = EL.copy(), ER.copy()
        El[:, :, :_DROP[c.left]] = 0.0; Er[:, :, :_DROP[c.right]] = 0.0
        out.append(np.concatenate([Wb, El, Er], axis=2))
    return out, labels


def dictionary_classes(op: PanelOperator, **kw):
    """Per-class bases for `operators.ritz_lambda` (admissible functions only, orthonormalised)."""
    D, _ = dictionary(op, **kw)
    res = []
    for B in D:
        keep = np.abs(B).sum(axis=(0, 1)) > 0
        res.append(orthonormalise(B[:, :, keep], op.wz, op.L, tol=1e-10))
    return res
