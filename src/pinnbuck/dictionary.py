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

HALF_WAVES = tuple(np.geomspace(25.0, 350.0, 6))   # [mm]  dictionary "d1" (design space v1, steps A and before)
# "d2" (design space v2, step B0): 11 half-waves 25 .. 1006 mm (ratio 1.45 instead of 1.70).  The longest
# clamped panels of v2 (L ~ 1000 mm) need the long waves (Ritz error 5.3 % -> 0.01 %) and the short panels
# (b_s ~ 125 mm) the closer spacing (0.3 -> 0.1 %); results/b0_check/final_probe_*.jsonl, RESEARCH_LOG 5 Nov.
DICTIONARIES = {"d1": HALF_WAVES, "d2": tuple(float(x) for x in np.geomspace(25.0, 1006.0, 11))}


# End layers of each dictionary version.  d2 (step B1): 16 splines over 2.5 b_s, spacing never below 1.5 hx.
# On a correctly integrated grid the longer zone is a real gain (e8 -> e16 0.10 % mean, general lay-ups
# 0.13 %); the spacing limit keeps the splines integrable at the solver's Gauss points (no Ritz value below the
# label).  24 test designs: mean/max 0.076/0.160 -> 0.053/0.116 % (industrial), 0.123/0.256 -> 0.078/0.141 %
# (general); results/b1_end_layer_probe/*, RESEARCH_LOG 5 Nov (B1).
# Revised the same day (B1, outlier probe): 24 splines instead of 16 (e24s15).  The worst test designs (tight
# pitch, general lay-ups, short-wave modes localised at both ends) need finer end splines: 5cc7be690865
# +1.53 -> +0.62 %, 6cc2e9904ad4 +0.47 -> +0.18 %; normal designs unchanged (there the 1.5 hx limit is already
# reached).  e32 = e24 everywhere (limit reached).  results/b1_targets/outlier_probe.jsonl.  The first test
# targets (e16s15) are kept in data/exact_targets/periodic_mixed_SS_v2_e16s15 for comparison.
END_LAYERS = {"d1": dict(n_end=8, end_length=1.0, min_spacing=0.0),
              "d2": dict(n_end=24, end_length=2.5, min_spacing=1.5)}


def version_of(op) -> str:
    return getattr(op, "info", {}).get("dictionary", "d1")


def half_waves_of(op) -> tuple:
    """Half-waves of the dictionary an operator was built for (op.info["dictionary"], default "d1")."""
    return DICTIONARIES[version_of(op)]


def dictionary_kw(op) -> dict:
    """All dictionary settings of an operator's version (half-waves and end layers), for RitzSpace."""
    v = version_of(op)
    return dict(half_waves=DICTIONARIES[v], **END_LAYERS[v])


N_END = 8                                          # end-layer splines per end
END_LENGTH = 1.0                                   # l_e / b_s
MIN_SPACING = 0.0                                  # smallest spline spacing / solver span hx (0: no limit)
_DROP = {"free": 0, "zero": 1, "clamp": 2}


def end_spacing(op, n_end, end_length, min_spacing=MIN_SPACING):
    """Spline spacing h of the end layers and the number of splines that fit in half the panel.
    Zone l_e = min(end_length b_s, L/2), h = l_e / n_end, but never below min_spacing * hx (hx: the
    solver's longitudinal span, where the pre-buckling field is known).  Splines finer than about the
    solver grid are integrated wrongly at its Gauss points (Ritz values below the label; RESEARCH_LOG
    5 Nov, B1); with the limit, splines that would reach past L/2 are left out (zero, dropped later)."""
    le = min(end_length * op.params["bs"], 0.5 * op.L)
    h = max(le / n_end, min_spacing * float(op.info.get("hx", 0.0)) if min_spacing else 0.0)
    n_fit = min(n_end, int(math.floor(0.5 * op.L / h * (1 + 1e-9))))
    return h, n_fit


def end_splines(z, L, le, ne, side, h=None):
    """(3, nz, ne): cubic B-splines in the distance s from the end (s = z or L - z);
    function 0 has value 1 at the end, function 1 the slope; z-derivatives.
    Spacing h (default le / ne)."""
    s = z if side == "L" else L - z
    h = le / ne if h is None else h
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


def dictionary(op: PanelOperator, half_waves=HALF_WAVES, n_end=N_END, end_length=END_LENGTH,
               min_spacing=MIN_SPACING):
    """Per class c: functions (3, nz, nD) of the common dictionary (zero where a function
    is not admissible for the class or does not fit in half the panel), and labels of the nD functions."""
    z, L = op.z, op.L
    h, n_fit = end_spacing(op, n_end, end_length, min_spacing)
    W = waves(z, half_waves)
    labels = [f"wave{hw:.0f}{sc}" for hw in half_waves for sc in ("s", "c")]
    labels += [f"endL{i}" for i in range(n_end)] + [f"endR{i}" for i in range(n_end)]
    EL, ER = end_splines(z, L, None, n_end, "L", h=h), end_splines(z, L, None, n_end, "R", h=h)
    EL[:, :, n_fit:] = 0.0; ER[:, :, n_fit:] = 0.0
    out = []
    for c in op.classes:
        bb = bubble(z / L, c.left, c.right, L)
        b0, b1, b2 = bb[0][:, None], bb[1][:, None], bb[2][:, None]
        Wb = np.stack([b0 * W[0], b1 * W[0] + b0 * W[1], b2 * W[0] + 2 * b1 * W[1] + b0 * W[2]])
        El, Er = EL.copy(), ER.copy()
        El[:, :, :_DROP[c.left]] = 0.0; Er[:, :, :_DROP[c.right]] = 0.0
        out.append(np.concatenate([Wb, El, Er], axis=2))
    return out, labels


def dictionary_classes(op: PanelOperator, **kw):  # noqa: D401
    """Per-class bases for `operators.ritz_lambda` (admissible functions only, orthonormalised)."""
    D, _ = dictionary(op, **kw)
    res = []
    for B in D:
        keep = np.abs(B).sum(axis=(0, 1)) > 0
        res.append(orthonormalise(B[:, :, keep], op.wz, op.L, tol=1e-10))
    return res
