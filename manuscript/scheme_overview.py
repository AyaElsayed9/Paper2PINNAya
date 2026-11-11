"""Scheme 1: flowchart of the periodic GBT-Ritz analysis (11 Nov 2026, third version; Aya: a scientific flowchart,
attractive, complete, and showing where the complexity is).

Flowchart elements: terminals (input / output), processes, one decision, three parallel branches, a loop.
Colour: red = new in this study (periodic modes); blue tag = differs from classical GBT (the four differences named
in Section 2); family chips in the colours of Figure 7.  The small pictures are computed results of the Unsym cell
(b_s = 150 mm, beta = 1.5, SS: the cell of Figures 6 and 7), drawn with the functions of the paper figures; the mode
counts are read from the solver's section (b_s = 150 mm; Table 2).  Only the hat functions in branch 3c are a
schematic (of Eq. (14)).
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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

INK, MUTED, GREY = FP.INK, "#5d6875", "#c9d0d8"
RED, RED_BG = "#c0392b", "#fbeeec"
BLUE, BLUE_BG = "#2a78d6", "#eef4fc"
BOX = "#ffffff"; LANE = "#f3f5f8"
LAY_, BS, BETA, BC = "Unsym", 150.0, 1.5, "SS"

h = FP._hier_mode(LAY_, BS, BETA, BC)
pre = FP._pre(LAY_, BS, BETA, BC)
md = pickle.loads((P.CACHE / f"mode_{LAY_}_{BS}_{BETA}_{BC}.pkl").read_bytes())
m_, sec, _L = P.case(LAY_, BS, BETA, BC)
grp = np.array(sec.group)
NF = {g: int((grp == g).sum()) for g in "GDLT"}; NF["S"] = int((grp == "S").sum()) + 1; NF["M"] = int((grp == "M").sum())
NTOT = sec.nm; NC = sec.nip
assert sum(NF.values()) == NTOT and NC == 3 * sec.N - 3
P.apply_style()
X = h["X"]; S = h["S"]; sL, sR = h["sides"]
rho = 100 * h["rho"]; dom = list(np.argsort(-np.abs(rho[:h["nip"]]))[:3])
CM = [FP.COL["Isotropic"], FP.COL["Cross-ply"], FP.COL["Sym1"]]

W, H = 160.0, 196.0
fig = plt.figure(figsize=(W * P.MM, H * P.MM))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis("off")


def axes(x, y, w, h_):
    a = fig.add_axes([x / W, y / H, w / W, h_ / H]); a.patch.set_alpha(0); return a


def rbox(x, y, w, h_, fc=BOX, ec=INK, lw=0.6, r=1.4, ls="-", z=2):
    bg.add_patch(FancyBboxPatch((x, y), w, h_, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, ls=ls, zorder=z))


def t(x, y, s, fs=6.6, c=INK, ha="left", va="center", **kw):
    bg.text(x, y, s, fontsize=fs, color=c, ha=ha, va=va, zorder=5, **kw)


def num(x, y, n, c=INK):
    bg.text(x, y, rf"\textbf{{{n}}}", ha="center", va="center", fontsize=7.6, color="white", zorder=6,
            bbox=dict(boxstyle="circle,pad=0.2", fc=c, ec="none"))


def arr(pts, c=INK, lw=0.8):
    """poly-line connector with an arrow head at the end"""
    pts = np.asarray(pts, float)
    if len(pts) > 2:
        bg.plot(pts[:-1, 0], pts[:-1, 1], color=c, lw=lw, zorder=3, solid_capstyle="round", solid_joinstyle="round")
    bg.add_patch(FancyArrowPatch(tuple(pts[-2]), tuple(pts[-1]), arrowstyle="-|>", mutation_scale=7.5, lw=lw, color=c, shrinkA=0,
                                 shrinkB=0, zorder=3))


def chip(x, y, lab, col, w=5.2):
    rbox(x, y - 1.9, w, 3.8, fc=col, ec="none", r=0.9, z=4)
    t(x + w / 2, y - 0.05, rf"\textbf{{{lab}}}", fs=6.2, c="white", ha="center")


def tag(x, y, s="differs from classical GBT"):
    bg.text(x, y, s, fontsize=5.6, color=BLUE, ha="right", va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", fc=BLUE_BG, ec=BLUE, lw=0.4))


def section(ax, U=None, amp=20.0, color=INK, lw=0.9, nodes=False, ms=3.4):
    segs0 = [[X[a], X[b]] for a, b, p_ in S if p_ != "link"]
    ax.add_collection(LineCollection(segs0, colors=GREY if U is not None else color, linewidths=0.5 if U is not None else lw))
    Xd = X if U is None else X + amp * U / np.hypot(*U.T).max()
    if U is not None:
        ax.add_collection(LineCollection([[Xd[a], Xd[b]] for a, b, p_ in S if p_ != "link"], colors=color, linewidths=lw, capstyle="round"))
    if nodes:
        ax.plot(X[:, 0], X[:, 1], "o", ms=1.0, mfc="white", mec=INK, mew=0.3, zorder=4)
    for n in (sL, sR):
        ax.plot(*Xd[n], "s", ms=ms, mfc="white", mec=RED, mew=0.8, zorder=6)
    ax.set_xlim(-320, 20); ax.set_aspect("equal"); ax.axis("off")
    return Xd


# ====================================================================== input
rbox(34, 187.5, 92, 7.5, fc=INK, ec="none", r=3.7)
t(80, 191.25, r"\textbf{INPUT}\quad cross-section $\cdot$ lay-ups $\cdot$ length $L$ $\cdot$ end conditions $\cdot$ load $n_0$", fs=7, c="white", ha="center")
arr([(80, 187.5), (80, 183.6)])

# ====================================================================== A  cross-section analysis
rbox(1.5, 62.5, 157, 120.5, fc=LANE, ec="none", r=2.5, z=0)
t(4.5, 179.6, r"\textbf{A\quad CROSS-SECTION ANALYSIS}", fs=7.4, c=MUTED)
t(155.5, 179.6, r"once per cross-section, reused for every length (Sections 2.1, 2.2)", fs=6.4, c=MUTED, ha="right")

# 1 discretize
rbox(5, 160, 72, 16.5); num(9, 172.6, 1); t(12.5, 172.6, r"\textbf{Discretize the cell}", fs=7.2)
t(8, 167.6, rf"walls $\rightarrow$ {len(sec.S)} segments, {sec.N} nodes"); t(8, 163.6, r"nodal $\bm{U}$, $\bm{\Theta}$, warping (Eq. (8))")
a1 = axes(45.5, 160.5, 31, 12.5); a1.set_ylim(-26, 52); section(a1, nodes=True, lw=0.7, ms=2.8)
a1.text(X[sL, 0], -10, r"$s_L$", ha="center", va="top", fontsize=6, color=RED); a1.text(X[sR, 0], -10, r"$s_R$", ha="center", va="top", fontsize=6, color=RED)
arr([(77, 168.25), (83, 168.25)])

# 2 mode families
rbox(83, 160, 72, 16.5); num(87, 172.6, 2); t(90.5, 172.6, r"\textbf{Mode families, side edges free}", fs=7.2)
tag(153.5, 172.6)
x0 = 86.0
for g in "GDLT":
    chip(x0, 166.6, g, P.FAMC[g]); x0 += 6.0
t(x0 + 0.3, 166.6, r"conventional"); x0 = 133.0
for g in "SM":
    chip(x0, 166.6, g, P.FAMC[g]); x0 += 6.0
t(86, 162.4, r"global, distortional, local, transv.\ extension $\cdot$ shear, Mindlin (Section 2.1.3)", fs=5.9, c=MUTED)
arr([(119, 160), (119, 156.6)], c=RED)

# 3 periodic modes (new)
rbox(5, 84, 150, 72.5, fc=RED_BG, ec=RED, lw=1.0, r=2.0, z=1)
num(9, 152.4, 3, c=RED); t(12.5, 152.4, r"\textbf{Make every single mode periodic}", fs=7.4, c=RED)
t(153, 152.4, r"\textbf{NEW}\quad equal values at $s_L$ and $s_R$, mode by mode (Proposition 1)", fs=6.4, c=RED, ha="right")

# --- branch a: conventional modes
bx, bw = 7.5, 62.0
rbox(bx, 86.5, bw, 61.5, ec=RED, lw=0.5)
x0 = bx + 2.5
for g in "GDLT":
    chip(x0, 144.2, g, P.FAMC[g], w=4.6); x0 += 5.2
t(x0 + 0.5, 144.2, r"\textbf{in-plane modes}, family by family", fs=6.6)
t(bx + 2.5, 139.8, r"periodic condition: $\bm{R}\,[\bm{U};\bm{\Theta}] = \bm{0}$ \ (3 equations, Eq. (12))", fs=6.2)
cx_, cy_ = bx + 20, 130.0
bg.add_patch(Polygon([(cx_ - 15, cy_), (cx_, cy_ + 5.6), (cx_ + 15, cy_), (cx_, cy_ - 5.6)], closed=True, fc="white", ec=INK, lw=0.6, zorder=3))
t(cx_, cy_ + 0.9, r"mode already", fs=6.0, ha="center"); t(cx_, cy_ - 1.8, r"periodic?", fs=6.0, ha="center")
arr([(cx_, 137.6), (cx_, cy_ + 5.6)])
# yes
rbox(bx + 40.5, 126.6, 19, 6.8, fc=BOX); t(bx + 50, 131.3, r"\textbf{keep it}", fs=6.2, ha="center"); t(bx + 50, 128.5, r"(step 1)", fs=5.6, c=MUTED, ha="center")
arr([(cx_ + 15, cy_), (bx + 40.5, cy_)]); t(cx_ + 16.2, cy_ + 1.5, r"yes", fs=5.8, c=MUTED)
# no
rbox(bx + 3, 108.2, 34, 11.6, fc=BOX)
t(bx + 20, 117.2, r"\textbf{combine within its family}", fs=6.2, ha="center")
t(bx + 20, 113.9, r"so that the side-edge", fs=5.9, ha="center"); t(bx + 20, 111.2, r"differences cancel (step 2)", fs=5.9, ha="center")
arr([(cx_, cy_ - 5.6), (cx_, 119.8)]); t(cx_ + 1.2, 122.4, r"no", fs=5.8, c=MUTED)
# complete
rbox(bx + 3, 93.8, 56.5, 9.6, fc=BOX)
t(bx + 31.2, 100.6, r"\textbf{add the periodic deformations still missing}", fs=6.2, ha="center")
t(bx + 31.2, 96.9, r"stretching a segment $\rightarrow$ T, otherwise $\rightarrow$ D (step 3)", fs=5.9, ha="center")
arr([(cx_, 108.2), (cx_, 103.4)]); arr([(bx + 50, 126.6), (bx + 50, 103.4)])
t(bx + 31, 90.0, rf"$\Rightarrow$ {NC} modes $= 3N-3$: nothing lost, nothing spurious", fs=6.1, c=RED, ha="center")

# --- branch b: warping
bx2, bw2 = 72.0, 40.5
rbox(bx2, 86.5, bw2, 61.5, ec=RED, lw=0.5)
chip(bx2 + 2.5, 144.2, "G", P.FAMC["G"], w=4.6); chip(bx2 + 7.7, 144.2, "D", P.FAMC["D"], w=4.6)
t(bx2 + 13.4, 144.2, r"\textbf{warping} $u_k$", fs=6.6)
tag(bx2 + bw2 - 1.5, 139.6, "closed contour")
rbox(bx2 + 2.5, 125.5, bw2 - 5, 10.5, fc=BOX)
t(bx2 + bw2 / 2, 132.6, r"$u_k(s_L) = u_k(s_R)$:", fs=6.1, ha="center"); t(bx2 + bw2 / 2, 129.6, r"the cell is a contour", fs=5.9, ha="center")
t(bx2 + bw2 / 2, 127.0, r"closed through its side edges", fs=5.9, ha="center")
rbox(bx2 + 2.5, 108.2, bw2 - 5, 12.3, fc=BOX)
t(bx2 + bw2 / 2, 117.6, r"\textbf{least shear strain energy}", fs=6.2, ha="center")
t(bx2 + bw2 / 2, 114.2, r"Eq. (9) $\rightarrow$ linear system", fs=5.9, ha="center"); t(bx2 + bw2 / 2, 111.2, r"Eq. (13)", fs=5.9, ha="center")
arr([(bx2 + bw2 / 2, 125.5), (bx2 + bw2 / 2, 120.5)])
t(bx2 + bw2 / 2, 101.6, r"open part: $\gamma_k = 0$ (Vlasov)", fs=5.9, ha="center")
t(bx2 + bw2 / 2, 98.4, r"closed part: shear flow kept", fs=5.9, ha="center")
t(bx2 + bw2 / 2, 90.0, r"$\Rightarrow$ unique periodic warping", fs=6.1, c=RED, ha="center")
arr([(bx2 + bw2 / 2, 108.2), (bx2 + bw2 / 2, 104.2)])

# --- branch c: shear and Mindlin
bx3, bw3 = 115.0, 37.5
rbox(bx3, 86.5, bw3, 61.5, ec=RED, lw=0.5)
chip(bx3 + 2.5, 144.2, "S", P.FAMC["S"], w=4.6); chip(bx3 + 7.7, 144.2, "M", P.FAMC["M"], w=4.6)
t(bx3 + 13.4, 144.2, r"\textbf{nodal functions}", fs=6.6)
t(bx3 + bw3 / 2, 139.2, r"one function per node", fs=5.9, ha="center")
ah = axes(bx3 + 3, 120.5, bw3 - 6, 15); xs_ = np.linspace(0, 1, 9)
for k_, c_ in ((0, RED), (8, RED)):                           # hat functions of the two side nodes (schematic of Eq. (14))
    y_ = np.maximum(0, 1 - np.abs(np.arange(9) - k_)); ah.plot(xs_, y_, color=c_, lw=0.9)
for k_ in (3, 4, 5):
    y_ = np.maximum(0, 1 - np.abs(np.arange(9) - k_)); ah.plot(xs_, y_, color=GREY, lw=0.6)
ah.plot([0, 1], [0, 0], color=INK, lw=0.5); ah.plot([0, 1], [1, 1], "s", ms=3.0, mfc="white", mec=RED, mew=0.8)
ah.text(0, 1.22, r"$s_L$", ha="center", va="bottom", fontsize=5.8, color=RED); ah.text(1, 1.22, r"$s_R$", ha="center", va="bottom", fontsize=5.8, color=RED)
ah.set_xlim(-0.08, 1.08); ah.set_ylim(-0.1, 1.6); ah.axis("off")
rbox(bx3 + 2.5, 104.5, bw3 - 5, 12.3, fc=BOX)
t(bx3 + bw3 / 2, 113.9, r"\textbf{merge the two}", fs=6.2, ha="center"); t(bx3 + bw3 / 2, 110.9, r"\textbf{side-node functions}", fs=6.2, ha="center")
t(bx3 + bw3 / 2, 107.3, r"into one, Eq. (14)", fs=5.9, ha="center")
arr([(bx3 + bw3 / 2, 120.5), (bx3 + bw3 / 2, 116.8)])
t(bx3 + bw3 / 2, 97.4, rf"$\Rightarrow$ {NF['S']} + {NF['M']} modes", fs=6.1, c=RED, ha="center")
t(bx3 + bw3 / 2, 90.0, r"(axial mode in S)", fs=5.7, c=MUTED, ha="center")

# result of step 3
arr([(80, 84), (80, 80.6)], c=RED)
rbox(5, 65, 150, 15.6, ec=INK)
t(8, 76.4, rf"\textbf{{Periodic basis: {NTOT} modes}}", fs=7.0)
t(8, 72.4, r"every mode, every truncated set and every", fs=6.1); t(8, 69.5, r"Ritz trial function is periodic", fs=6.1)
num(64.5, 76.4, 4); t(68, 76.4, r"\textbf{Stiffness tensors}", fs=7.0)
t(62, 72.4, r"cross-section integrals of the", fs=6.1); t(62, 69.5, r"modes, Eqs. (17)--(26)", fs=6.1)
bg.plot([59, 59], [66.5, 79.1], color=GREY, lw=0.5, zorder=4)
for j, k in enumerate(dom):
    am = axes(96.5 + 19.3 * j, 65.6, 19, 11.0); am.set_ylim(-34, 52)
    section(am, h["U"][:, k].reshape(-1, 2), amp=22.0, color=CM[j], lw=0.8, ms=2.2)
    t(106 + 19.3 * j, 78.2, rf"mode {k + 1}", fs=5.6, c=CM[j], ha="center")

# ====================================================================== B  member analysis
arr([(80, 62.5), (80, 58.6)])
rbox(1.5, 12.5, 157, 46.0, fc=LANE, ec="none", r=2.5, z=0)
t(4.5, 55.2, r"\textbf{B\quad MEMBER ANALYSIS}", fs=7.4, c=MUTED)
t(155.5, 55.2, r"repeated for every length $L$ and end condition (Section 2.3)", fs=6.4, c=MUTED, ha="right")
bw_ = 47.5
for k_, x_ in enumerate((5, 56.25, 107.5)):
    rbox(x_, 15.5, bw_, 36.0)
    if k_:
        arr([(x_ - 3.75, 33.5), (x_, 33.5)])
# 5 pre-buckling
num(9, 47.6, 5); t(12.5, 47.6, r"\textbf{Pre-buckling state}", fs=7.0); tag(50.5, 43.6, "computed, not assumed")
a4 = axes(13.5, 24.0, 35, 15.5)
i = int(np.argmin(np.abs(pre["xg"] - pre["L"] / 2))); sk = np.abs(pre["pts"][:, 1]) < 1e-9; o = np.argsort(pre["pts"][sk, 0]); Ys = pre["pts"][sk, 0][o]
a4.plot(Ys, pre["Pxx"][i, sk][o], color=INK, lw=0.8); a4.plot(Ys, pre["Pss"][i, sk][o], color=RED, lw=0.8)
a4.set_xlim(-300, 0); a4.set_ylim(0, 2.4); a4.set_xticks([-300, 0]); a4.set_yticks([0, 1, 2]); a4.tick_params(labelsize=5.2, length=1.6, pad=1.5)
a4.set_xlabel(r"$Y$ [mm]", fontsize=5.6, labelpad=-4)
a4.text(-150, 1.12, r"$-S_{xx}/n_0$", ha="center", va="bottom", fontsize=5.6, color=INK); a4.text(-150, 0.14, r"$-S_{ss}/n_0$", ha="center", va="bottom", fontsize=5.6, color=RED)
t(28.75, 18.4, r"linear solution under $n_0$, Eqs. (30), (31)", fs=5.9, ha="center")
# 6 geometric
num(60.25, 47.6, 6); t(63.75, 47.6, r"\textbf{Geometric tensors}", fs=7.0)
t(80, 40.6, r"second-order work of", fs=6.1, ha="center"); t(80, 37.4, r"$S_{xx}$, $S_{ss}$, $S_{xs}$", fs=6.4, ha="center")
t(80, 34.0, r"on the buckling mode", fs=6.1, ha="center")
t(80, 27.4, r"$S_{ss}$ comes from the constant", fs=5.9, c=MUTED, ha="center"); t(80, 24.6, r"cell width (Section 3.2)", fs=5.9, c=MUTED, ha="center")
t(80, 18.4, r"Eqs. (27)--(29)", fs=5.9, ha="center")
# 7 Ritz
num(111.5, 47.6, 7); t(115, 47.6, r"\textbf{Ritz solution}", fs=7.0)
a5 = axes(116, 24.0, 35, 15.5); A = h["amp"]; Am = np.abs(A[:, dom]).max(); sg = np.sign(A[np.argmax(np.abs(A[:, dom[0]])), dom[0]])
for j, k in enumerate(dom):
    a5.plot(h["xs"] / h["L"], sg * A[:, k] / Am, color=CM[j], lw=0.8)
a5.axhline(0, color=INK, lw=0.35); a5.set_xlim(0, 1); a5.set_ylim(-1.15, 1.15); a5.set_xticks([0, 1]); a5.set_yticks([-1, 0, 1])
a5.tick_params(labelsize=5.2, length=1.6, pad=1.5); a5.set_xlabel(r"$x/L$", fontsize=5.6, labelpad=-4)
t(131.25, 43.2, r"amplitudes $\phi_k(x)$ of modes " + ", ".join(str(k + 1) for k in dom), fs=5.9, ha="center")
t(131.25, 18.4, r"B-splines $\cdot$ $(\bm{K} + \lambda\,\bm{G})\,\bm{a} = \bm{0}$, Eqs. (37), (38)", fs=5.9, ha="center")

# ====================================================================== output
arr([(131.25, 15.5), (131.25, 9.6)])
rbox(64, 1.5, 94.5, 8.0, fc=INK, ec="none", r=4.0)
t(111.25, 5.5, rf"\textbf{{OUTPUT}}\quad $N_{{cr}}$ ({md['N']:.2f} N/mm here) $\cdot$ mode $\cdot$ family participations (Eq. (41))", fs=6.8, c="white", ha="center")
# legend
rbox(2.5, 2.2, 3.4, 3.4, fc=RED_BG, ec=RED, lw=0.8, r=0.6); t(7.2, 3.9, r"new in this study", fs=5.9)
bg.text(31.5, 3.9, r"tag", fontsize=5.4, color=BLUE, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", fc=BLUE_BG, ec=BLUE, lw=0.4))
t(35.0, 3.9, r"differs from classical GBT", fs=5.9)
t(2.5, 8.6, rf"pictures: computed, {LAY_} cell, $b_s = {BS:g}$ mm, $\beta = {BETA:g}$, {BC}", fs=5.7, c=MUTED)
P.save(fig, "scheme01_overview")
print("modes", NF, NTOT, "NC", NC, "dominant", [int(k) + 1 for k in dom], "N", md["N"])
