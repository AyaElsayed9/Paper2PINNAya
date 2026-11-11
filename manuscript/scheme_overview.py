"""Scheme 1: flowchart of the periodic GBT-Ritz analysis (11 Nov 2026, fourth version).
Drawn in the style of the flowchart of the authors' previous paper (Fig. 3 of [22]): plain rectangles with thin black
borders, light blue process boxes, light green decision diamonds, an elliptic end node, elbow connectors, equations
inside the boxes, Yes/No labels and small drawings.  The drawings are computed (Unsym cell, b_s = 150 mm, beta = 1.5,
SS, the cell of Figures 6 and 7); the mode count is read from the solver's section.
"""
import pickle
import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import fp_common as P
import figures_paper as FP
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Ellipse, FancyArrowPatch, Polygon, Rectangle

LAY_, BS, BETA, BC = "Unsym", 150.0, 1.5, "SS"
h = FP._hier_mode(LAY_, BS, BETA, BC)
md = pickle.loads((P.CACHE / f"mode_{LAY_}_{BS}_{BETA}_{BC}.pkl").read_bytes())
m_, sec, _L = P.case(LAY_, BS, BETA, BC)
NTOT = sec.nm
plt.rcParams.update({"text.usetex": False, "font.family": "Carlito", "mathtext.fontset": "stix", "font.size": 7})

BLUE = "#dceaf7"; GREEN = "#d9f2d0"; PINK = "#f8cef0"; K = "black"; GREY = "#b4b4b4"; RED = "#c0392b"
X = h["X"]; S = h["S"]; sL, sR = h["sides"]
rho = 100 * h["rho"]; dom = list(np.argsort(-np.abs(rho[:h["nip"]]))[:3])

W, H = 160.0, 113.0
fig = plt.figure(figsize=(W * P.MM, H * P.MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")


def box(x, y, w, h_, fc="white", lw=0.8):
    ax.add_patch(Rectangle((x, y), w, h_, fc=fc, ec=K, lw=lw, zorder=2))


def txt(x, y, s, fs=7, ha="center", va="center", **kw):
    ax.text(x, y, s, fontsize=fs, ha=ha, va=va, zorder=5, **kw)


def diamond(cx, cy, w, h_, lines, fs=6.6):
    ax.add_patch(Polygon([(cx - w / 2, cy), (cx, cy + h_ / 2), (cx + w / 2, cy), (cx, cy - h_ / 2)], closed=True, fc=GREEN,
                         ec=K, lw=0.8, zorder=2))
    txt(cx, cy, lines, fs=fs, linespacing=1.15)


def line(pts, head=True):
    pts = np.asarray(pts, float)
    if len(pts) > 2:
        ax.plot(pts[:-1, 0], pts[:-1, 1], color=K, lw=0.7, zorder=3, solid_joinstyle="miter")
    if head:
        ax.add_patch(FancyArrowPatch(tuple(pts[-2]), tuple(pts[-1]), arrowstyle="-|>", mutation_scale=6, lw=0.7, color=K,
                                     shrinkA=0, shrinkB=0, zorder=3))
    else:
        ax.plot(pts[-2:, 0], pts[-2:, 1], color=K, lw=0.7, zorder=3)


def inset(x, y, w, h_):
    a = fig.add_axes([x / W, y / H, w / W, h_ / H]); a.patch.set_alpha(0); a.axis("off"); return a


def section(a, U=None, amp=22.0, col=K, lw=0.8, ms=2.6, nodes=False):
    segs = [[X[p], X[q]] for p, q, part in S if part != "link"]
    a.add_collection(LineCollection(segs, colors=GREY if U is not None else col, linewidths=0.5 if U is not None else lw))
    Xd = X if U is None else X + amp * U / np.hypot(*U.T).max()
    if U is not None:
        a.add_collection(LineCollection([[Xd[p], Xd[q]] for p, q, part in S if part != "link"], colors=col, linewidths=lw))
    if nodes:
        a.plot(X[:, 0], X[:, 1], "o", ms=0.9, mfc="white", mec=K, mew=0.3, zorder=4)
    for n in (sL, sR):
        a.plot(*Xd[n], "s", ms=ms, mfc="white", mec=RED, mew=0.8, zorder=6)
    a.set_xlim(-320, 20); a.set_ylim(-34, 52); a.set_aspect("equal")
    return Xd


# ---------------------------------------------------------------------------------------------- row 1
box(2, 92, 21, 20)
txt(12.5, 107.5, "Panel data", fs=7.4)
txt(12.5, 99.5, "$b_s$, $\\beta$, $L=\\beta b_s$\nlay-ups (ABD)\nend conditions\nload $n_0$", fs=6.6, linespacing=1.3)
line([(23, 102), (28, 102)])

box(28, 92, 38, 20)
txt(47, 108.6, "Periodic cell (two pitches)", fs=7.2)
a = inset(29, 93.0, 36, 13.0); Xd = section(a, nodes=True)
a.text(X[sL, 0], -12, "$s_L$", ha="center", va="top", fontsize=6.5, color=RED)
a.text(X[sR, 0], -12, "$s_R$", ha="center", va="top", fontsize=6.5, color=RED)
line([(66, 102), (71, 102)])

box(71, 92, 39, 20, fc=BLUE)
txt(90.5, 108.6, "Discretization (Section 2.1.3)", fs=7.2)
txt(90.5, 101.6, f"{sec.N} nodes, {len(sec.S)} segments", fs=6.6)
txt(90.5, 96.6, "$w(\\eta)=H_1 w_q+b_pH_2\\theta_q+H_3 w_{q'}+b_pH_4\\theta_{q'}$", fs=6.4)
line([(110, 102), (115, 102)])

box(115, 92, 43, 20, fc=BLUE)
txt(136.5, 108.6, "GBT mode families", fs=7.2)
txt(136.5, 104.2, "(side edges still free)", fs=6.4, style="italic")
txt(136.5, 97.6, "G global, D distortional, L local,\nT transverse extension,\nS shear, M Mindlin", fs=6.4, linespacing=1.25)

# ---------------------------------------------------------------------------------------------- periodic modes
ax.add_patch(Rectangle((28, 39.5), 130, 49, fc=PINK, ec=K, lw=0.8, zorder=1, alpha=0.45))
txt(30, 85.6, "Periodic modes: every mode takes equal values at $s_L$ and $s_R$ (Section 2.2)", fs=7.2, ha="left")
# feed from the families box (outside the frame, left)
line([(136.5, 92), (136.5, 90.0), (25, 90.0), (25, 44.5)], head=False)
for yy in (71.0, 53.5, 44.5):
    line([(25, yy), (34, yy)])

# lane 1: in-plane modes, family by family
box(34, 66, 24, 10, fc=BLUE)
txt(46, 71, "In-plane modes\nG, D, L, T\n(each family)", fs=6.4, linespacing=1.15)
line([(58, 71), (62, 71)])
diamond(73, 71, 22, 12.0, "$\\mathbf{R}\\,\\mathbf{v}=\\mathbf{0}$ ?\n(Eq. (12))", fs=6.2)
box(91, 74, 17, 6.0); txt(99.5, 77, "keep $\\mathbf{v}$", fs=6.4)
line([(73, 77), (73, 79.6), (88, 79.6), (88, 77), (91, 77)]); txt(75, 80.8, "Yes", fs=6.4, ha="left")
box(90, 60.2, 39, 10.6)
txt(109.5, 65.5, "combine within the family:\n$\\mathbf{V}'_g\\mathbf{Z}_g$ with $\\mathbf{R}\\mathbf{V}'_g\\mathbf{Z}_g=\\mathbf{0}$\nmissing periodic fields $\\rightarrow$ T or D", fs=6.1, linespacing=1.3)
line([(73, 65), (73, 63.4), (87, 63.4), (87, 65.5), (90, 65.5)]); txt(75, 62.2, "No", fs=6.4, ha="left")
line([(108, 77), (131, 77), (131, 71), (133, 71)])
line([(129, 65.5), (131, 65.5), (131, 71.5)], head=False)

# lane 2: warping
box(34, 49.0, 24, 9, fc=BLUE)
txt(46, 53.5, "Warping $u_k$\nof G and D", fs=6.4, linespacing=1.2)
line([(58, 53.5), (62, 53.5)])
box(62, 49.0, 66, 9)
txt(95, 55.7, "least shear strain energy, cell closed through $s_L$, $s_R$:", fs=6.4)
txt(95, 51.4, "$\\min\\ \\Pi_w(\\mathbf{u})$  with  $u_{s_L}=u_{s_R}$,  $\\mathbf{c}_{ax}^T\\mathbf{u}=0$   (Eqs. (9), (13))", fs=6.4)
line([(128, 53.5), (133, 53.5)])

# lane 3: shear and Mindlin
box(34, 41.0, 24, 7, fc=BLUE)
txt(46, 44.5, "Shear S, Mindlin M", fs=6.4)
line([(58, 44.5), (62, 44.5)])
box(62, 41.0, 66, 7)
txt(95, 44.5, "merge the side-node functions:  $u=u_{s_L}+u_{s_R}$,  $\\Omega=\\Omega_{s_L}+\\Omega_{s_R}$  (Eq. (14))", fs=6.4)
line([(128, 44.5), (133, 44.5)])

# periodic basis with three computed modes
box(133, 41.0, 23, 42.5, fc="white")
txt(144.5, 79.6, f"periodic basis\n{NTOT} modes", fs=6.6, linespacing=1.15)
for j, k in enumerate(dom):
    a = inset(133.6, 62.0 - 9.6 * j, 21.8, 9.6)
    section(a, h["U"][:, k].reshape(-1, 2), amp=26.0, col=[FP.COL["Isotropic"], FP.COL["Cross-ply"], FP.COL["Sym1"]][j], lw=0.7, ms=1.8)

# ---------------------------------------------------------------------------------------------- row 3, member analysis
line([(144.5, 41.0), (144.5, 35)])
box(126, 22, 32, 13, fc=BLUE)
txt(142, 31.4, "Modal stiffness tensors", fs=6.8)
txt(142, 25.8, "$\\Xi_{de}=\\int_b\\mathbf{E}_d^T\\mathbf{C}_{lam}\\mathbf{E}_e\\,ds$  (Eq. (26))", fs=6.3)
line([(126, 28.5), (121, 28.5)])

box(89, 22, 32, 13, fc=BLUE)
txt(105, 31.4, "Pre-buckling state", fs=6.8)
txt(105, 25.6, "$\\delta U+\\delta\\Pi_{ext}=0\\ \\rightarrow\\ S_{xx},S_{ss},S_{xs}$", fs=6.3)
line([(89, 28.5), (84, 28.5)])

box(52, 22, 32, 13, fc=BLUE)
txt(68, 31.4, "Geometric tensors", fs=6.8)
txt(68, 25.8, "$\\mathbf{X}_{de}$ from $S_{xx},S_{ss},S_{xs}$  (Eq. (29))", fs=6.3)
line([(52, 28.5), (47, 28.5)])

box(13, 22, 34, 13, fc=BLUE)
txt(30, 31.4, "Ritz solution (cubic B-splines)", fs=6.8)
txt(30, 25.8, "$(\\mathbf{K}+\\lambda\\,\\mathbf{G})\\,\\mathbf{a}=\\mathbf{0}$,  smallest $\\lambda>0$", fs=6.3)

line([(30, 22), (30, 18.0)])
diamond(30, 12, 26, 12.0, "other length or\nend condition?", fs=6.4)
line([(43, 12), (105, 12), (105, 22)]); txt(46, 13.6, "Yes (cross-section reused)", fs=6.4, ha="left")
line([(17, 12), (8, 12), (8, 5.8)], head=True); txt(12.5, 13.4, "No", fs=6.4, ha="center")
ax.add_patch(Ellipse((17, 3.3), 30, 6.2, fc="white", ec=K, lw=0.8, zorder=2))
txt(17, 3.3, "$N_{cr}$, buckling mode, $\\rho_g$", fs=6.6)
a = inset(117, 3.0, 40, 17.5)
U = md["U"].copy()
if U.flat[np.argmax(np.abs(U))] < 0:
    U = -U
FP._mode3d_ortho(a, dict(md, U=U), amp=14.0, azim=-35.0, elev=28.0, cmap=plt.get_cmap("RdBu_r"), norm=TwoSlopeNorm(0, -1, 1))
txt(137, 1.0, f"example: critical mode of the {LAY_} cell, $N_{{cr}}$ = {md['N']:.2f} N/mm", fs=6.0, ha="center")
P.save(fig, "scheme01_overview")
print("N", md["N"], "modes", NTOT, [int(k) + 1 for k in dom])
