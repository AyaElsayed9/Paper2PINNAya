"""Illustration for Aya (22 Nov 2026, C61): cubic Hermite elements (C1) vs cubic B-splines with simple knots (C2).
Same test deflection on 5 elements: both approximate w well, but the curvature w'' of the Hermite interpolant
jumps at the element ends, while the cubic spline's curvature is continuous."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline, CubicSpline

w = lambda x: np.sin(np.pi * x) + 0.3 * np.sin(3 * np.pi * x)
dw = lambda x: np.pi * np.cos(np.pi * x) + 0.9 * np.pi * np.cos(3 * np.pi * x)
d2w = lambda x: -np.pi**2 * np.sin(np.pi * x) - 2.7 * np.pi**2 * np.sin(3 * np.pi * x)
xn = np.linspace(0, 1, 6)                                   # 5 elements / spans
H = CubicHermiteSpline(xn, w(xn), dw(xn))                   # C1: value + slope shared at nodes
S = CubicSpline(xn, w(xn), bc_type=((1, dw(0)), (1, dw(1))))  # cubic spline = cubic B-splines, simple knots: C2
x = np.linspace(0, 1, 2001)
BLUE, ORANGE, INK, MUTED = "#2a6fdb", "#e07b39", "#222222", "#9a9a9a"
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6), dpi=200)
for a in ax:
    for xv in xn[1:-1]:
        a.axvline(xv, color=MUTED, lw=0.6, ls=(0, (3, 3)), zorder=0)
    a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=8)
    a.set_xlabel("x / L", fontsize=9)
ax[0].plot(x, w(x), color=MUTED, lw=3, alpha=0.5, label="exact")
ax[0].plot(x, H(x), color=BLUE, lw=1.6, label="cubic Hermite (C1)")
ax[0].plot(x, S(x), color=ORANGE, lw=1.6, ls=(0, (5, 2)), label="cubic B-spline (C2)")
ax[0].set_title("(a) deflection w: both close to the exact curve", fontsize=9.5, loc="left")
ax[0].legend(fontsize=7.5, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.18))
for k in range(5):                                          # plot Hermite curvature element by element (jumps visible)
    xs = np.linspace(xn[k], xn[k + 1], 300)
    ax[1].plot(xs, H(xs, 2), color=BLUE, lw=1.8, label="cubic Hermite (C1)" if k == 0 else None)
ax[1].plot(x, S(x, 2), color=ORANGE, lw=1.6, ls=(0, (5, 2)), label="cubic B-spline (C2)")
for xv in xn[1:-1]:
    lo, hi = H(xv - 1e-9, 2), H(xv + 1e-9, 2)
    ax[1].plot([xv, xv], [lo, hi], color=BLUE, lw=0.8, ls=":")
    ax[1].annotate("jump", (xv, 0.5 * (lo + hi)), xytext=(6, 0), textcoords="offset points", fontsize=7.5, color=BLUE, va="center")
ax[1].set_title("(b) curvature w'': Hermite jumps at element ends", fontsize=9.5, loc="left")
ax[1].legend(fontsize=7.5, frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.18))
fig.text(0.5, -0.10, "dashed grey lines: element ends (Hermite) = knots (B-spline); same nodal data for both", ha="center", fontsize=8, color=INK)
fig.tight_layout(); fig.savefig("results/hermite_vs_bspline.png", bbox_inches="tight")
for xv in xn[1:-1]:
    print(f"x={xv:.3f}: Hermite w'' left {H(xv-1e-9,2):+.2f} right {H(xv+1e-9,2):+.2f} | spline left {S(xv-1e-9,2):+.2f} right {S(xv+1e-9,2):+.2f}")
