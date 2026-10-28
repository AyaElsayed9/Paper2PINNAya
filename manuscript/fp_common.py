"""Shared helpers and the publication style of the figures (figures_paper.py)."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample          # noqa: E402
from pinnbuck.oracle import omega_length, omega_panel             # noqa: E402
from pinnbuck.solver import configure                             # noqa: E402

H = Path(__file__).parent
OUT = H / "fig"; OUT.mkdir(exist_ok=True)
CACHE = H / "figp"; CACHE.mkdir(exist_ok=True)

MM = 1 / 25.4
FULL = 160 * MM                 # figure width [in]; inserted at 15.5 cm
STYLE = {
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{mathpazo}\usepackage{amsmath}\usepackage{bm}",
    "font.family": "serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "lines.linewidth": 1.0, "lines.markersize": 4, "legend.frameon": False, "legend.handlelength": 1.8,
    "axes.titlepad": 4, "savefig.dpi": 600, "pdf.fonttype": 42, "axes.unicode_minus": True,
    "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
}


def apply_style():
    """The solver module sets its own rcParams on import: restore the paper style."""
    plt.rcdefaults(); plt.rcParams.update(STYLE)


apply_style()
INK = "#1b2733"; MUTED = "#8a94a0"; GRID = "#e3e7ec"
LAY = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
COL = dict(zip(LAY, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]))
MK = dict(zip(LAY, ["o", "s", "^", "D", "v"]))
FAMC = {"G": "#2a78d6", "D": "#eb6834", "L": "#1baf7a", "T": "#eda100", "S": "#e87ba4", "M": "#008300"}
FAMN = {"G": "global (G)", "D": "distortional (D)", "L": "local (L)", "T": "transverse extension (T)",
        "S": "shear (S)", "M": "Mindlin (M)"}
PART = {"skin": "#1b2733", "foot": "#2a78d6", "web": "#eb6834", "crown": "#1baf7a", "link": "#9aa3ad"}


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight", pad_inches=0.02, facecolor="white")
    plt.close(fig)
    from PIL import Image
    im = Image.open(OUT / f"{name}.png"); im.convert("RGB").save(OUT / f"{name}.png", optimize=True, dpi=(600, 600))


def grid(ax):
    ax.grid(True, color=GRID, lw=0.4, zorder=0); ax.set_axisbelow(True)


def case(lay, bs, beta, bc):
    """(module, section, L) of a periodic cell (same set-up as run_periodic.py)."""
    cm = configure(OMEGA_SIDE_BC="periodic"); m = cm.__enter__()
    if lay == "Isotropic":
        pan = m.make_omega_panel(bs, "iso", [], []); L = beta * bs
    else:
        s = validated_omega_sample(lay, bs, beta, bc, "periodic"); pan = omega_panel(s); L = omega_length(s)
    sec = m.omega_gbt_section(pan)
    apply_style()
    return m, sec, L


def points(sec):
    """Coordinates (np, 2) of the cross-section Gauss points (4 per segment, order of sec.S)."""
    gp = np.polynomial.legendre.leggauss(4)[0]; xi = 0.5 * (gp + 1.0)
    return np.array([sec.X[s["a"]] + q * (sec.X[s["b"]] - sec.X[s["a"]]) for s in sec.S for q in xi])


def seg_lines(sec, X=None):
    X = sec.X if X is None else X
    return np.array([[X[s["a"]], X[s["b"]]] for s in sec.S])


def amplitudes(lb, C, xs, d=0):
    from scipy.interpolate import BSpline
    spl = BSpline(lb.t, np.eye(lb.nb), 3)
    if d:
        spl = spl.derivative(d)
    return spl(np.clip(xs, 0, lb.breaks[-1])) @ C.T      # (nx, nm)


def nodal_disp(sec, lb, C, xs):
    """Nodal in-plane displacements (nx, N, 2) and warping (nx, N) of a mode."""
    A0 = amplitudes(lb, C, xs); A1 = amplitudes(lb, C, xs, 1); nip = sec.nip
    U = (A0[:, :nip] @ sec.U.T).reshape(len(xs), sec.N, 2)
    u = A1[:, :nip] @ sec.Uw.T + A0[:, [sec.o_ax]]
    if sec.ns:
        u = u + A0[:, sec.o_s:sec.o_s + sec.ns] @ sec.Us.T
    return U, u


def panel_label(ax, s, x=-0.01, y=1.0, **kw):
    ax.text(x, y, s, transform=ax.transAxes, ha="left", va="bottom", fontsize=8, **kw)
