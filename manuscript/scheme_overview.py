"""Scheme 1: the periodic GBT-Ritz analysis in pictures (11 Nov 2026; Aya: not paragraphs in boxes).
Every picture is drawn from computed data of ONE cell of the paper (Unsym, b_s = 150 mm, beta = 1.5, SS: the cell of
Figures 6 and 7), with the drawing functions and colours of the paper figures:
  1  the discretised cross-section (solver section: walls, nodes, side nodes)
  2  a mode of the cell with free side edges that is not periodic: the rigid rotation about X (built from the
     nodal coordinates exactly as in the solver, U = omega x (X - centroid); Section 2.2.2 names it)
  3  three periodic modes: the dominant hierarchical conventional modes of the critical mode (Figure 7)
  4  the computed pre-buckling forces across the skin at mid-length (Figure 3 data)
  5  the amplitude functions of the same three modes (Figure 7e data)
  6  the critical mode and N_cr (Figure 6 data)
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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK, MUTED, RED, GREY = FP.INK, FP.MUTED, "#c0392b", "#c9d0d8"
LAY_, BS, BETA, BC = "Unsym", 150.0, 1.5, "SS"

h = FP._hier_mode(LAY_, BS, BETA, BC)                       # hierarchical modes, amplitudes, participations
pre = FP._pre(LAY_, BS, BETA, BC)                           # pre-buckling forces
md = pickle.loads((P.CACHE / f"mode_{LAY_}_{BS}_{BETA}_{BC}.pkl").read_bytes())      # critical mode (Figure 6)
P.apply_style()
X = h["X"]; S = h["S"]; sL, sR = h["sides"]
rho = 100 * h["rho"]; dom = list(np.argsort(-np.abs(rho[:h["nip"]]))[:3])
CM = [FP.COL["Isotropic"], FP.COL["Cross-ply"], FP.COL["Sym1"]]

W, H = 160.0, 128.0
fig = plt.figure(figsize=(W * P.MM, H * P.MM))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis("off")


def axes(x, y, w, h_):
    return fig.add_axes([x / W, y / H, w / W, h_ / H])


def band(y, h_, label):
    bg.add_patch(FancyBboxPatch((1.5, y), W - 3, h_, boxstyle="round,pad=0,rounding_size=2", fc="#f6f7f9", ec="none"))
    bg.text(4, y + h_ - 2.2, label, ha="left", va="top", fontsize=7.6, color=MUTED)


def arrow(x0, y0, x1, y1, col=INK):
    bg.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9, lw=0.9, color=col, shrinkA=0, shrinkB=0))


def step(x, y, n, title, new=False):
    bg.text(x, y, rf"\textbf{{{n}}}", ha="center", va="center", fontsize=8, color="white",
            bbox=dict(boxstyle="circle,pad=0.22", fc=RED if new else INK, ec="none"))
    bg.text(x + 3.2, y, rf"\textbf{{{title}}}", ha="left", va="center", fontsize=7.8, color=RED if new else INK)


def section(ax, U=None, amp=20.0, color=INK, lw=1.0, squares=True, nodes=False):
    segs0 = [[X[a], X[b]] for a, b, p_ in S if p_ != "link"]
    ax.add_collection(LineCollection(segs0, colors=GREY if U is not None else color, linewidths=0.6 if U is not None else lw))
    Xd = X if U is None else X + amp * U / np.hypot(*U.T).max()
    if U is not None:
        ax.add_collection(LineCollection([[Xd[a], Xd[b]] for a, b, p_ in S if p_ != "link"], colors=color, linewidths=lw,
                                         capstyle="round"))
    if nodes:
        ax.plot(X[:, 0], X[:, 1], "o", ms=1.3, mfc="white", mec=INK, mew=0.35, zorder=4)
    if squares:
        for n in (sL, sR):
            ax.plot(*Xd[n], "s", ms=4.2, mfc="white", mec=RED, mew=0.9, zorder=6)
    ax.set_xlim(-318, 18); ax.set_aspect("equal"); ax.axis("off")
    return Xd


# ------------------------------------------------------------------ cross-section analysis
band(70.5, 56.0, "Cross-section analysis: once per cross-section (Sections 2.1, 2.2)")
step(6.5, 117.0, 1, "Cell")
a1 = axes(3, 96, 46, 18); a1.set_ylim(-22, 52); section(a1, nodes=True)
a1.text(X[sL, 0], -9, r"$s_L$", ha="center", va="top", fontsize=7.5, color=RED); a1.text(X[sR, 0], -9, r"$s_R$", ha="center", va="top", fontsize=7.5, color=RED)
bg.text(26, 93.5, rf"{len(X)} nodes; displacements, rotations" + "\n" + r"and warping at the nodes (Eq. (8))", ha="center", va="top", fontsize=6.6,
        color=INK, linespacing=1.25)

step(57.5, 117.0, 2, "Modes, side edges free")
a2 = axes(54, 92, 46, 22); a2.set_ylim(-34, 62)
c = X.mean(axis=0); Urot = np.c_[-(X[:, 1] - c[1]), X[:, 0] - c[0]]            # rigid rotation, as in the solver
Xd = section(a2, Urot, amp=24.0)
a2.annotate("", xy=Xd[sL], xytext=X[sL], arrowprops=dict(arrowstyle="-|>", lw=0.7, color=RED, mutation_scale=6))
a2.annotate("", xy=Xd[sR], xytext=X[sR], arrowprops=dict(arrowstyle="-|>", lw=0.7, color=RED, mutation_scale=6))
bg.text(77, 91.0, r"families G, D, L, T, S, M (Section 2.1.3)" + "\n" + r"e.g. rigid rotation: $s_L \neq s_R$, \textbf{not periodic}",
        ha="center", va="top", fontsize=6.6, color=INK, linespacing=1.25)

step(110.5, 117.0, 3, "Periodic modes (new)", new=True)
for j, k in enumerate(dom):
    a3 = axes(108, 101.0 - 12.5 * j, 46, 12.5); a3.set_ylim(-30, 46)
    section(a3, h["U"][:, k].reshape(-1, 2), amp=20.0, color=CM[j], lw=1.0)
    a3.text(-150, 44, rf"mode {k + 1}", ha="center", va="top", fontsize=6.3, color=CM[j])
bg.text(131, 75.4, r"every mode: $s_L = s_R$ (Eqs. (12)--(14))", ha="center", va="top", fontsize=6.6, color=RED)
arrow(49.5, 104, 54.5, 104); arrow(101, 104, 106, 104, col=RED)

# ------------------------------------------------------------------ Eq. (1): modes x amplitudes
bg.text(80, 65.6, r"displacement field $=\ \sum_k\ $ \textit{periodic mode} $k\,(s)\ \times\ $ \textit{amplitude} $\phi_k(x)$"
        r"\quad (Eq. (1)) $\ \Rightarrow\ $ periodic for any amplitudes (Proposition 1)", ha="center", va="center", fontsize=7.4, color=INK)
arrow(131, 70.0, 131, 68.0, col=RED); arrow(80, 63.2, 80, 60.6)

# ------------------------------------------------------------------ member analysis
band(13.5, 47.0, "Member analysis: for every length and end condition (Section 2.3)")
step(6.5, 51.5, 4, "Pre-buckling state")
a4 = axes(10.5, 23.5, 36, 23)
i = int(np.argmin(np.abs(pre["xg"] - pre["L"] / 2))); sk = np.abs(pre["pts"][:, 1]) < 1e-9; o = np.argsort(pre["pts"][sk, 0])
Ys = pre["pts"][sk, 0][o]
a4.plot(Ys, pre["Pxx"][i, sk][o], color=INK, lw=0.9); a4.plot(Ys, pre["Pss"][i, sk][o], color=RED, lw=0.9)
a4.set_xlim(-300, 0); a4.set_ylim(0, None); a4.set_xlabel(r"$Y$ [mm]", fontsize=6.5, labelpad=1); a4.tick_params(labelsize=5.8, length=2)
a4.set_xticks([-300, -150, 0]); FP.grid(a4)
yx = float(np.median(pre["Pxx"][i, sk])); ys = float(np.median(pre["Pss"][i, sk]))
a4.text(-150, yx + 0.22, r"$-S_{xx}/n_0$", ha="center", va="bottom", fontsize=6.3, color=INK)
a4.text(-150, ys + 0.10, r"$-S_{ss}/n_0$", ha="center", va="bottom", fontsize=6.3, color=RED)
bg.text(27, 17.2, r"computed, not assumed (Eq. (31))", ha="center", va="top", fontsize=6.6, color=INK)

step(57.5, 51.5, 5, "Ritz solution")
a5 = axes(61.5, 23.5, 36, 23)
A = h["amp"]; Am = np.abs(A[:, dom]).max(); sg = np.sign(A[np.argmax(np.abs(A[:, dom[0]])), dom[0]])
for j, k in enumerate(dom):
    a5.plot(h["xs"] / h["L"], sg * A[:, k] / Am, color=CM[j], lw=1.0)
a5.axhline(0, color=INK, lw=0.4); a5.set_xlim(0, 1); a5.set_xlabel(r"$x/L$", fontsize=6.5, labelpad=1)
a5.set_ylabel(r"$\phi_k(x)$", fontsize=6.5, labelpad=1); a5.tick_params(labelsize=5.8, length=2); a5.set_yticks([-1, 0, 1]); FP.grid(a5)
bg.text(78, 17.2, r"B-splines; $(\bm{K} + \lambda\,\bm{G})\,\bm{a} = \bm{0}$ (Eq. (38))", ha="center", va="top", fontsize=6.6, color=INK)

step(108.5, 51.5, 6, "Critical load and mode")
a6 = axes(104, 18.5, 52, 30)
U = md["U"].copy()
if U.flat[np.argmax(np.abs(U))] < 0:
    U = -U
FP._mode3d_ortho(a6, dict(md, U=U), amp=14.0, azim=-35.0, elev=28.0, cmap=plt.get_cmap("RdBu_r"), norm=TwoSlopeNorm(0, -1, 1))
bg.text(131, 17.2, rf"$N_{{cr}} = {md['N']:.2f}$ N/mm; mode families (Eq. (41))", ha="center", va="top", fontsize=6.6, color=INK)
arrow(48.5, 35, 53.5, 35); arrow(99.5, 35, 104.0, 35)

bg.text(80, 9.0, rf"Example: {LAY_} cell, $b_s = {BS:g}$ mm, $\beta = {BETA:g}$, {BC} (Figures 6 and 7). Modes {', '.join(str(k + 1) for k in dom)} "
        rf"are its dominant hierarchical modes ({rho[dom].sum():.1f}\,\% of the strain energy).", ha="center", va="center", fontsize=6.6, color=MUTED)
P.save(fig, "scheme01_overview")
print("dominant modes", [int(k) + 1 for k in dom], "rho", [round(float(rho[k]), 1) for k in dom], "N", md["N"])
