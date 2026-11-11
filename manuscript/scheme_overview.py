"""Scheme 1: overview of the periodic GBT-Ritz analysis (11 Nov 2026; Navid/Aya: the theory must be easy to follow).
Only the steps and equation numbers of the manuscript; drawn in the style of the paper figures."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import fp_common as P
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

INK, MUTED, RED, FILL, FILL_NEW = "#1f2a36", "#5b6775", "#c0392b", "#f3f5f8", "#fdf1ef"
P.apply_style()
W, H = 160.0, 123.0                                    # mm
fig = plt.figure(figsize=(W * P.MM, H * P.MM)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")


def box(x, y, w, h, title, lines, new=False, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=FILL_NEW if new else FILL,
                                ec=RED if new else INK, lw=0.9 if new else 0.5))
    ax.text(x + 2.2, y + h - 2.4, (r"\textbf{%s}" % title) if bold else title, ha="left", va="top", fontsize=7.6, color=INK)
    for k, s in enumerate(lines):
        ax.text(x + 2.2, y + h - 6.6 - 3.7 * k, s, ha="left", va="top", fontsize=6.9, color=INK)


def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK, shrinkA=0, shrinkB=0,
                                                                   mutation_scale=7))


def group(x, y, w, h, label):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5", fc="none", ec=MUTED, lw=0.5, ls=(0, (3, 2))))
    ax.text(x + 2, y + h + 1.0, label, ha="left", va="bottom", fontsize=7.4, color=MUTED)


# input
box(22, 113.5, 116, 8.5, "Input", [])
ax.text(35, 117.7, r"cross-section, lay-ups, length $L$, end conditions (SS or CC), end load $n_0$", ha="left", va="center", fontsize=6.9, color=INK)
arrow(80, 113.5, 80, 109.4)

# cross-section analysis
group(3, 41.5, 154, 64.0, r"Cross-section analysis (Sections 2.1 and 2.2): once per cross-section, reused for every length")
box(6, 89.5, 72, 13.5, "1  Discretization",
    [r"walls $\rightarrow$ segments and nodes; nodal displacements $\bm{U}$,", r"rotations $\bm{\Theta}$ and warping (Section 2.1.3, Eq. (8))"])
box(82, 89.5, 72, 13.5, "2  Mode families of the cell with free side edges",
    [r"G, D, L, T (conventional), S (shear), M (Mindlin)", r"(Section 2.1.3)"])
arrow(78, 96.25, 82, 96.25)
arrow(118, 89.5, 118, 85.7)
box(6, 58.0, 148, 27.5, r"3  Periodic modes (new): every single mode takes equal values at the side edges $s_L$ and $s_R$", [], new=True)
bw = 46.0
for k, (ttl, l1, l2, l3) in enumerate((
        ("conventional modes", r"each family restricted to its", r"periodic subspace, $\bm{R}\,[\bm{U};\bm{\Theta}] = \bm{0}$", r"(Section 2.2.2, Eq. (12))"),
        ("warping of G and D modes", r"least shear strain energy on the", r"cell closed through its side edges", r"(Section 2.2.3, Eqs. (9), (13))"),
        ("shear and Mindlin modes", r"functions of the two side nodes", r"merged into one", r"(Section 2.2.4, Eq. (14))"))):
    x = 9 + k * (bw + 2.0)
    ax.add_patch(FancyBboxPatch((x, 60.0), bw, 19.5, boxstyle="round,pad=0,rounding_size=1.0", fc="white", ec=RED, lw=0.5))
    ax.text(x + 2, 77.5, r"\textit{%s}" % ttl, ha="left", va="top", fontsize=6.9, color=INK)
    for j_, s_ in enumerate((l1, l2, l3)):
        ax.text(x + 2, 73.5 - 3.7 * j_, s_, ha="left", va="top", fontsize=6.7, color=INK)
arrow(80, 58.0, 80, 54.2)
box(6, 43.2, 148, 11.0, "4  Modal stiffness tensors",
    [r"cross-section integrals of the periodic modes with the laminate stiffnesses (Section 2.3.1, Eqs. (17)--(26))"])
arrow(80, 41.5, 80, 36.0)

# member analysis
group(3, 13.0, 154, 19.0, r"Member analysis (Section 2.3): for every length and end condition")
bw = 47.5
for k, (ttl, lines) in enumerate((
        ("5  Pre-buckling state", [r"linear solution under the end load:", r"$S_{xx}$, $S_{ss}$, $S_{xs}$ (Eqs. (30), (31))"]),
        ("6  Geometric tensors", [r"second-order work of the", r"pre-buckling forces (Eqs. (27)--(29))"]),
        ("7  Ritz solution", [r"cubic B-splines, end conditions,", r"$(\bm{K} + \lambda\,\bm{G})\,\bm{a} = \bm{0}$ (Eqs. (37), (38))"]))):
    x = 6 + k * (bw + 2.75)
    box(x, 15.5, bw, 14.0, ttl, lines)
    if k:
        arrow(x - 2.75, 22.5, x, 22.5)
arrow(80, 13.0, 80, 9.2)
box(26, 0.6, 108, 8.5, "Output", [])
ax.text(41, 4.8, r"critical load $N_{cr}$, buckling mode, participation of the mode families (Eq. (41))", ha="left", va="center",
        fontsize=6.9, color=INK)
P.save(fig, "scheme01_overview")
