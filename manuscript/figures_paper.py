"""Publication figures of the manuscript (vector PDF + 600 dpi PNG in fig/).
All data come from the solver (V32) and the result files; nothing is tuned.
    python3 figures_paper.py [fig01 fig02 ...]"""
import json
import pickle
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
import matplotlib.ticker as mt

import fp_common as P
from fp_common import FULL, MM, INK, MUTED, LAY, COL, MK, FAMC, FAMN, PART, save, grid, panel_label

KREF = lambda bs: np.pi / bs      # reference wavenumber of the hierarchical ordering (= solver default)


def cached(name, fn):
    f = P.CACHE / f"{name}.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    d = fn(); P.apply_style(); f.write_bytes(pickle.dumps(d)); return d


def draw_sec(ax, sec, X=None, color=None, lw=0.9, z=2, alpha=1.0, ls="-"):
    segs = P.seg_lines(sec, X)
    cols = [color or PART[sec.wpart[s["w"]]] for s in sec.S]
    lss = [(0, (2, 1.2)) if sec.wpart[s["w"]] == "link" and color is None else ls for s in sec.S]
    lc = LineCollection(segs, colors=cols, linewidths=lw, zorder=z, alpha=alpha, linestyles=lss, capstyle="round")
    ax.add_collection(lc)
    return lc


# =====================================================================================
# Figure 1 -- panel, periodic cell, load, GBT discretisation
# =====================================================================================
def fig01():
    m, sec, L = P.case("Cross-ply", 150.0, 1.5, "SS")
    X = sec.X; W = np.ptp(X[:, 0]); x0, x1 = X[:, 0].min(), X[:, 0].max()
    fig = plt.figure(figsize=(FULL, 118 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.55, 0.62, 1.05], hspace=0.12)

    # (a) oblique view of the cell with the load and the periodic side edges -----------------
    a = fig.add_subplot(gs[0]); c, sn = 0.62, 0.36          # oblique projection of the x axis
    pr = lambda XX, YY, xx: (XX + c * xx, YY + sn * xx)
    Ldraw = L
    # walls as filled quadrilaterals, back to front (painter): larger Y and larger x first
    quads = []
    for s in sec.S:
        A_, B_ = X[s["a"]], X[s["b"]]; part = sec.wpart[s["w"]]
        if part == "link":
            continue
        pts = [pr(A_[0], A_[1], 0), pr(B_[0], B_[1], 0), pr(B_[0], B_[1], Ldraw), pr(A_[0], A_[1], Ldraw)]
        depth = 0.5 * (A_[1] + B_[1]) + (0.0 if part in ("skin", "foot") else 1.0)
        quads.append((depth, part, pts))
    quads.sort(key=lambda q: q[0])
    fillc = {"skin": "#dfe7f1", "foot": "#c9d8ec", "web": "#f6dccd", "crown": "#d3eee2"}
    for _d, part, pts in quads:
        a.add_patch(plt.Polygon(pts, closed=True, fc=fillc[part], ec="none", lw=0, zorder=1 + 0.01 * _d))
    for xx, lw in ((0.0, 0.8), (Ldraw, 0.8)):
        for s in sec.S:
            if sec.wpart[s["w"]] == "link":
                continue
            A_, B_ = X[s["a"]], X[s["b"]]
            a.plot(*zip(pr(A_[0], A_[1], xx), pr(B_[0], B_[1], xx)), color=INK, lw=lw, zorder=5, solid_capstyle="round")
    for n in sec.junctions + sec.tips:                    # longitudinal fold lines
        if abs(X[n, 1] - 0.938) < 0.01:
            continue
        a.plot(*zip(pr(X[n, 0], X[n, 1], 0), pr(X[n, 0], X[n, 1], Ldraw)), color=INK, lw=0.35, zorder=4)
    sL, sR = sec.sides
    for n in (sL, sR):
        a.plot(*zip(pr(X[n, 0], 0, 0), pr(X[n, 0], 0, Ldraw)), color="#c0392b", lw=1.3, zorder=6)
    xm, ym = pr(X[sL, 0], 0, 0.3 * Ldraw)
    a.text(xm - 11, ym, r"side edge $s_L$", ha="right", va="center", fontsize=7.5, color="#c0392b")
    a.text(X[sR, 0] + 6, -13, r"side edge $s_R$", ha="left", va="center", fontsize=7.5, color="#c0392b")
    # end load: arrows along the loaded edges (compression, pointing into the cell)
    for xx, sgn in ((0.0, 1.0), (Ldraw, -1.0)):
        for Xa in np.linspace(x0 + 12, x1 - 12, 11):
            p0 = pr(Xa, 0, xx - sgn * 34); p1 = pr(Xa, 0, xx - sgn * 3)
            a.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="-|>", lw=0.5, color=INK, mutation_scale=5), zorder=7)
    a.text(*pr(x0 - 4, 0, -38), r"$n_0$", ha="right", va="center", fontsize=8)
    a.text(*pr(x1 + 4, 0, Ldraw + 38), r"$n_0$", ha="left", va="center", fontsize=8)
    # axes triad
    from matplotlib.patches import FancyArrowPatch
    o = np.array([x0 - 30, 52.0])
    for d_, lab in ((np.array([28.0, 0]), r"$X$"), (np.array([0, 26.0]), r"$Y$"), (np.array([c, sn]) * 46, r"$x$")):
        a.add_patch(FancyArrowPatch(tuple(o), tuple(o + d_), arrowstyle="-|>", mutation_scale=6, lw=0.6, color=INK, zorder=9))
        a.text(*(o + d_ * 1.12 + np.array([1.5, 0])), lab, fontsize=8, ha="left", va="center")
    # dimension of the cell length
    p0, p1 = pr(x1 + 22, 0, 0), pr(x1 + 22, 0, Ldraw)
    a.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="<|-|>", lw=0.5, mutation_scale=5, color=MUTED))
    a.text(0.5 * (p0[0] + p1[0]) + 10, 0.5 * (p0[1] + p1[1]) - 10, r"$L = \beta\, b_s$", fontsize=8, ha="left", va="top")
    p0, p1 = pr(x0, 0, -55), pr(x1, 0, -55)
    a.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="<|-|>", lw=0.5, mutation_scale=5, color=MUTED))
    a.text(0.5 * (p0[0] + p1[0]), p0[1] - 4, r"$W = 2b_s$", fontsize=8, ha="center", va="top")
    a.set_aspect("equal"); a.axis("off")
    a.set_xlim(x0 - 75, x1 + c * Ldraw + 95); a.set_ylim(-62, sn * Ldraw + 40)
    panel_label(a, r"(a) periodic cell of the stiffened panel: compressive end load $n_0$ and periodic side edges", x=0.0, y=0.98)

    # (b) infinite panel and the cell ------------------------------------------------------
    b = fig.add_subplot(gs[1])
    for k in (-1, 0, 1):
        Xs = X.copy(); Xs[:, 0] += k * W
        draw_sec(b, sec, Xs, color=INK if k == 0 else "#b6bec8", lw=0.9 if k == 0 else 0.7)
    b.axvspan(x0, x1, color="#eef3f8", zorder=0, lw=0)
    for xx in (x0, x1):
        b.axvline(xx, color="#c0392b", lw=0.6, ls=(0, (3, 2)), zorder=1)
    for xs in (x0 - W - 18, x1 + W + 18):
        b.text(xs, 12, r"$\cdots$", ha="center", va="center", fontsize=11)
    st = [X[sec.walls[i][0], 0] for i, p in enumerate(sec.wpart) if p == "crown"]
    b.annotate("", xy=(-75 - 150, 38), xytext=(-75, 38), arrowprops=dict(arrowstyle="<|-|>", lw=0.5, mutation_scale=5, color=MUTED))
    b.text(-150, 40, r"$b_s$", ha="center", va="bottom", fontsize=8)
    b.set_xlim(x0 - W - 40, x1 + W + 40); b.set_ylim(-6, 62); b.set_aspect("equal"); b.axis("off")
    panel_label(b, r"(b) infinitely wide panel with identical stringers at the pitch $b_s$; the cell contains two stringers", x=0.0, y=0.98)

    # (c) GBT discretisation ---------------------------------------------------------------
    d = fig.add_subplot(gs[2])
    draw_sec(d, sec, lw=1.3)
    nat = set(sec.junctions + sec.tips)
    inter = [n for n in range(sec.N) if n not in nat]
    d.plot(X[inter, 0], X[inter, 1], "o", ms=1.7, mfc="w", mec=INK, mew=0.45, zorder=4)
    natl = sorted(nat)
    d.plot(X[natl, 0], X[natl, 1], "o", ms=3.0, mfc=INK, mec=INK, zorder=5)
    for n, lab in ((sL, r"$s_L$"), (sR, r"$s_R$")):
        d.plot(*X[n], "s", ms=5.5, mfc="w", mec="#c0392b", mew=1.0, zorder=6)
        d.text(X[n, 0], -5.5, lab, ha="center", va="top", fontsize=8, color="#c0392b")
    # local axes on the right crown
    o = np.array([-150.0, 0.0]); t = np.array([1.0, 0.0]); nrm = np.array([0.0, 1.0])
    from matplotlib.patches import FancyArrowPatch
    o = o + np.array([0.0, 0.0])
    for d_ in (24 * t, 13 * nrm):
        d.add_patch(FancyArrowPatch(tuple(o), tuple(o + d_), arrowstyle="-|>", mutation_scale=6, lw=0.7, color=INK, zorder=9))
    d.text(*(o + 24 * t + np.array([0, 2.2])), r"$s,\ v$", fontsize=7.5, ha="center", va="bottom")
    d.text(*(o + 13 * nrm + np.array([0, 0.8])), r"$z,\ w$", fontsize=7.5, ha="center", va="bottom")
    d.plot(*o, "o", ms=3.5, mfc="w", mec=INK, mew=0.6, zorder=7); d.plot(*o, ".", ms=1.5, color=INK, zorder=8)
    d.text(o[0] - 3, o[1] + 1.5, r"$x,\ u$", fontsize=7.5, ha="right", va="bottom")
    hs = [Line2D([], [], color=PART[k], lw=1.6, ls=(0, (2, 1.2)) if k == "link" else "-") for k in PART] + \
         [Line2D([], [], ls="", marker="o", ms=3.0, mfc=INK, mec=INK),
          Line2D([], [], ls="", marker="o", ms=1.9, mfc="w", mec=INK, mew=0.45)]
    d.legend(hs, ["skin", "foot (skin + flange)", "web", "crown", "link wall (offset tie)", "natural node", "intermediate node"],
             ncol=7, loc="lower center", bbox_to_anchor=(0.5, -0.36), fontsize=6.8, handlelength=1.5, columnspacing=1.1, frameon=False)
    d.set_xlim(x0 - 8, x1 + 8); d.set_ylim(-9, 42); d.set_aspect("equal"); d.axis("off")
    panel_label(d, rf"(c) GBT discretisation of the cell ($b_s = 150$ mm): {sec.N} nodes, {len(sec.walls)} walls, {len(sec.S)} segments", x=0.0, y=1.0)
    save(fig, "fig01_geometry")


# =====================================================================================
# Figure 2 -- hierarchical periodic conventional modes (in-plane shape + warping)
# =====================================================================================
def fig02():
    m, sec, L = P.case("Cross-ply", 150.0, 1.5, "SS")
    lam, T = m.hierarchical_modes(sec, kref=KREF(150.0))["C"]
    o = np.argsort(lam)
    g0 = np.array(sec.group[:sec.nip])
    fig, axs = plt.subplots(4, 2, figsize=(FULL, 92 * MM), gridspec_kw=dict(hspace=0.42, wspace=0.06))
    norm = TwoSlopeNorm(0, -1, 1); cmap = plt.get_cmap("RdBu_r")
    sL, sR = sec.sides
    out = []
    for k, (ax, i) in enumerate(zip(axs.ravel(), o[:8])):
        U = (sec.U @ T[:, i]).reshape(-1, 2); uw = sec.Uw @ T[:, i]
        sc = 20.0 / max(np.hypot(*U.T).max(), 1e-12)
        Xd = sec.X + sc * U
        draw_sec(ax, sec, color="#c9d0d8", lw=0.6, z=1)
        segs = P.seg_lines(sec, Xd)
        um = np.abs(uw).max()
        if um > 1e-9 * max(np.hypot(*U.T).max(), 1e-30):
            val = np.array([0.5 * (uw[s["a"]] + uw[s["b"]]) for s in sec.S]) / um
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=1.25, zorder=3, capstyle="round")
            lc.set_array(val); ax.add_collection(lc)
        else:
            ax.add_collection(LineCollection(segs, colors=INK, linewidths=1.1, zorder=3, capstyle="round"))
        for n in (sL, sR):
            ax.plot(*Xd[n], "s", ms=3.2, mfc="w", mec="#c0392b", mew=0.8, zorder=5)
        lam_s = r"0" if lam[i] < 1e-12 else rf"{lam[i]:.2e}".replace("e-0", r"\times10^{-").replace("e+0", r"\times10^{") + "}"
        ax.set_title(rf"mode {k + 1}: $\Lambda = {lam_s}$", fontsize=7.5, pad=1)
        ax.set_xlim(-326, 26); ax.set_ylim(-26, 56); ax.set_aspect("equal"); ax.axis("off")
        out.append(dict(mode=k + 1, lam=float(lam[i]),
                        dU=float(np.abs(U[sL] - U[sR]).max()), du=float(abs(uw[sL] - uw[sR]))))
    cax = fig.add_axes([0.35, 0.035, 0.3, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label(r"warping displacement $u_k(s)/|u_k|_{\max}$", fontsize=7.5); cb.ax.tick_params(labelsize=6.5, length=2)
    cb.outline.set_linewidth(0.4)
    save(fig, "fig02_modes")
    return out


# =====================================================================================
# Figure 3 -- pre-buckling membrane forces
# =====================================================================================
def _pre(lay, bs, beta, bc):
    def run():
        m, sec, L = P.case(lay, bs, beta, bc)
        r = m.omega_gbt_solve(sec, L, bc, return_mode=True)
        pts = P.points(sec)
        part = np.array([sec.wpart[s["w"]] for s in sec.S for _q in range(4)])
        return dict(L=L, xg=np.asarray(r["lb"].xg), pts=pts, part=part, Pxx=r["pre"][0], Pss=r["pre"][1], Pxs=r["pre"][2],
                    N=r["N"])
    return cached(f"pre_{lay}_{bs}_{beta}_{bc}", run)


def fig03():
    fig = plt.figure(figsize=(FULL, 112 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.22)
    out = {}
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
    for lay in LAY:
        d = _pre(lay, 150.0, 1.5, "SS")
        i = int(np.argmin(np.abs(d["xg"] - d["L"] / 2)))
        sk = np.abs(d["pts"][:, 1]) < 1e-9
        o = np.argsort(d["pts"][sk, 0]); Xs = d["pts"][sk, 0][o]
        me = (LAY.index(lay) * 7, 36)
        kw = dict(color=COL[lay], lw=0.9, marker=MK[lay], ms=2.6, mfc="w", mew=0.6, markevery=me)
        ax_a.plot(Xs, d["Pxx"][i, sk][o], label=lay, **kw)
        ax_b.plot(Xs, d["Pss"][i, sk][o], label=lay, **kw)
        out[lay] = dict(Sxx_skin_mid=float(np.median(d["Pxx"][i, sk & (d["part"] == "skin")])),
                        Sss_skin_mid=float(np.median(d["Pss"][i, sk & (d["part"] == "skin")])))
    for ax in (ax_a, ax_b):
        d = _pre("Isotropic", 150.0, 1.5, "SS")
        ft = d["pts"][d["part"] == "foot"]
        for X0 in np.unique(np.round(ft[:, 0] / 25) * 25):
            pass
        segs = sorted(set((round(a, 3)) for a in ft[:, 0]))
        # foot bands: contiguous runs of foot points along X
        xs = np.sort(ft[np.abs(ft[:, 1]) < 1e-9, 0]); br = np.nonzero(np.diff(xs) > 5)[0]
        for st, en in zip(np.r_[0, br + 1], np.r_[br, len(xs) - 1]):
            ax.axvspan(xs[st] - 1, xs[en] + 1, color="#eef1f5", lw=0, zorder=0)
        ax.set_xlim(-300, 0); ax.set_xlabel(r"$X$ [mm]"); ax.xaxis.set_major_locator(mt.MultipleLocator(50))
        grid(ax)
    ax_a.set_ylabel(r"$-S_{xx}/n_0$"); ax_b.set_ylabel(r"$-S_{ss}/n_0$")
    ax_a.set_ylim(0, 2.3); ax_a.text(-75, 2.12, "foot", ha="center", fontsize=6.5, color=MUTED)
    ax_a.text(-225, 2.12, "foot", ha="center", fontsize=6.5, color=MUTED)
    panel_label(ax_a, r"(a) skin and feet at $x = L/2$, SS")
    panel_label(ax_b, r"(b) skin and feet at $x = L/2$, SS")
    # (c), (d): S_xx and S_ss along x in the middle of the skin bay between the stringers
    ax_c = fig.add_subplot(gs[1, 0]); ax_d = fig.add_subplot(gs[1, 1])
    for lay in LAY:
        for bc, ls in (("SS", "-"), ("CC", (0, (3.5, 1.6)))):
            d = _pre(lay, 150.0, 1.5, bc)
            sk = (np.abs(d["pts"][:, 1]) < 1e-9) & (d["part"] == "skin")
            j = np.nonzero(sk)[0][np.argmin(np.abs(d["pts"][sk, 0] + 150.0))]
            xr = d["xg"] / d["L"]
            kw = dict(color=COL[lay], lw=0.9, ls=ls, marker=MK[lay], ms=2.6, mfc="w" if bc == "SS" else COL[lay], mew=0.6,
                      markevery=(LAY.index(lay) * 3 + (0 if bc == "SS" else 9), 18))
            ax_c.plot(xr, d["Pxx"][:, j], **kw)
            ax_d.plot(xr, d["Pss"][:, j], **kw)
            out[f"{lay}_{bc}_bay_end"] = (float(d["Pxx"][0, j]), float(d["Pss"][0, j]), float(d["Pss"][len(xr) // 2, j]))
    for ax in (ax_c, ax_d):
        ax.set_xlim(0, 1); ax.set_xlabel(r"$x/L$"); grid(ax)
    ax_c.set_ylabel(r"$-S_{xx}/n_0$"); ax_d.set_ylabel(r"$-S_{ss}/n_0$")
    panel_label(ax_c, r"(c) middle of the skin bay, $X = -b_s$")
    panel_label(ax_d, r"(d) middle of the skin bay, $X = -b_s$")
    hs = [Line2D([], [], color=COL[l], lw=1.0, marker=MK[l], ms=2.8, mfc="w", mew=0.6) for l in LAY] + [
        Line2D([], [], color=INK, lw=0.9, marker="o", ms=2.6, mfc="w", mew=0.6),
        Line2D([], [], color=INK, lw=0.9, ls=(0, (3.5, 1.6)), marker="o", ms=2.6, mew=0.6)]
    fig.legend(hs, LAY + ["SS (open markers)", "CC (dashed, filled)"], ncol=7, loc="upper center",
               bbox_to_anchor=(0.5, 0.975), fontsize=6.8, handlelength=1.8, columnspacing=1.2)
    save(fig, "fig03_prebuckling")
    return out


# =====================================================================================
# Figures 4-6 -- critical loads and accuracy
# =====================================================================================
def _load():
    """GBT curves, Abaqus SS (results table) and, when available, Abaqus CC
    (src/abaqus_cc_periodic.csv = omega_validation_eigs_CC_PBC.csv of the CC script)."""
    import csv
    cur = [json.loads(l) for l in (P.H / "results" / "periodic_curves.jsonl").read_text().splitlines()]
    gbt = {(c["layup"], c["bs"], c["beta"], c["bc"]): c["N"] for c in cur}
    ref = []
    for r in csv.DictReader(open(P.H / "src" / "omega_results_V32.csv")):
        if r["bc"] == "SS":
            ref.append(dict(layup=r["layup"], bc="SS", bs=float(r["bs"]), beta=float(r["beta"]), abq=float(r["abq_per"])))
    fcc = P.H / "src" / "abaqus_cc_periodic.csv"
    if fcc.exists():
        for r in csv.DictReader(open(fcc)):
            if r["Ncr"] not in ("", "NA"):
                ref.append(dict(layup=r["layup"], bc="CC", bs=float(r["bs"]), beta=float(r["beta"]), abq=float(r["Ncr"])))
    for r in ref:
        r["gbt"] = gbt[(r["layup"], r["bs"], r["beta"], r["bc"])]
    return cur, ref


def _curves(bc, name):
    """N_cr versus beta: rows SS (and CC), columns b_s = 150 and 285.5 mm; Abaqus markers where available."""
    cur, ref = _load()
    bcs = ["SS", "CC"] if bc is None else [bc]
    fig, axs = plt.subplots(len(bcs), 2, figsize=(FULL, (58 * len(bcs) + 8) * MM), squeeze=False,
                            gridspec_kw=dict(wspace=0.18, hspace=0.38))
    k = 0
    for row, b_ in enumerate(bcs):
        for col, bs in enumerate((150.0, 285.5)):
            ax = axs[row][col]
            for lay in LAY:
                d = sorted([c for c in cur if c["layup"] == lay and c["bs"] == bs and c["bc"] == b_], key=lambda c: c["beta"])
                ax.plot([c["beta"] for c in d], [c["N"] for c in d], "-", color=COL[lay], lw=1.0, marker=".", ms=3.2, zorder=3)
                rr = [r for r in ref if r["layup"] == lay and r["bs"] == bs and r["bc"] == b_]
                ax.plot([r["beta"] for r in rr], [r["abq"] for r in rr], MK[lay], ls="none", mfc="none", mec=COL[lay], ms=6.0,
                        mew=0.9, zorder=4)
            ax.set_yscale("log"); ax.set_xlim(0.4, 3.1)
            ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}"))
            ax.yaxis.set_minor_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}" if f"{v:g}"[0] in "2357" and v < 1000 else ""))
            ax.tick_params(axis="y", which="minor", labelsize=6.3)
            grid(ax); ax.grid(True, which="minor", axis="y", color=P.GRID, lw=0.3)
            if row == len(bcs) - 1:
                ax.set_xlabel(r"aspect ratio $\beta = L/b_s$")
            if col == 0:
                ax.set_ylabel(r"$N_{cr}$ [N/mm]")
            ends = "simply supported" if b_ == "SS" else "clamped"
            panel_label(ax, rf"({'abcd'[k]}) {ends} loaded edges, $b_s = {bs:g}$ mm")
            k += 1
    hs = [Line2D([], [], color=COL[l], lw=1.0, marker=MK[l], ms=5, mfc="none", mew=0.8) for l in LAY] + \
         [Line2D([], [], color=INK, lw=1.0, marker=".", ms=3.2),
          Line2D([], [], ls="none", marker="o", mfc="none", mec=INK, ms=5.5, mew=0.8)]
    fig.subplots_adjust(top=0.90 if len(bcs) == 2 else 0.82)
    fig.legend(hs, LAY + ["GBT", "Abaqus"], ncol=7, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=6.8,
               columnspacing=1.2, handlelength=1.8)
    save(fig, name)


def fig04():
    _curves(None, "fig04_ncr")


def fig05():
    """Deviation of GBT from Abaqus, case by case (SS; CC once the Abaqus CC results exist)."""
    cur, ref = _load()
    bcs = [b for b in ("SS", "CC") if any(r["bc"] == b for r in ref)]
    fig, axs = plt.subplots(1, len(bcs), figsize=(FULL if len(bcs) == 2 else 0.6 * FULL, 66 * MM),
                            gridspec_kw=dict(wspace=0.12), sharey=True, squeeze=False)
    out = {}
    for k, (ax, bc) in enumerate(zip(axs[0], bcs)):
        devs = []
        for g, lay in enumerate(LAY):
            rr = sorted([r for r in ref if r["layup"] == lay and r["bc"] == bc], key=lambda r: (r["bs"], r["beta"]))
            for q, r in enumerate(rr):
                x = g + (q - 2.5) * 0.12 + (0.06 if r["bs"] > 200 else -0.06)
                e = 100 * (r["gbt"] / r["abq"] - 1)
                fill = COL[lay] if r["bs"] < 200 else "w"
                ax.plot(x, e, MK[lay], mfc=fill, mec=COL[lay], ms=4.6, mew=0.8, zorder=4); devs.append(e)
            if g:
                ax.axvline(g - 0.5, color=P.GRID, lw=0.6, zorder=0)
        ax.axhspan(-1, 1, color="#f1f4f7", zorder=0, lw=0); ax.axhline(0, color=INK, lw=0.5, zorder=1)
        ax.set_xticks(range(5)); ax.set_xticklabels(LAY); ax.tick_params(axis="x", which="both", length=0)
        lim = max(1.75, 0.25 * np.ceil(4 * (max(np.abs(devs)) + 0.2)))
        ax.set_xlim(-0.55, 4.55); ax.set_ylim(-lim, lim); ax.yaxis.set_major_locator(mt.MultipleLocator(0.5))
        ax.tick_params(axis="x", top=False)
        mean_abs = float(np.mean(np.abs(devs))); mx = float(np.max(np.abs(devs)))
        out[bc] = dict(mean_abs=mean_abs, max_abs=mx, n=len(devs))
        panel_label(ax, (rf"({'ab'[k]}) " if len(bcs) == 2 else "") +
                    ("simply supported loaded edges" if bc == "SS" else "clamped loaded edges"))
        ax.text(0.99, 0.97, rf"mean $|\Delta|$ = {mean_abs:.2f}\,\%, max = {mx:.2f}\,\%", transform=ax.transAxes,
                ha="right", va="top", fontsize=6.8, bbox=dict(fc="white", ec="none", pad=1.0))
    axs[0][0].set_ylabel(r"deviation $\Delta$ of GBT from Abaqus [\%]")
    hs = [Line2D([], [], ls="none", marker="o", mfc=INK, mec=INK, ms=4.6), Line2D([], [], ls="none", marker="o", mfc="w", mec=INK, ms=4.6, mew=0.8),
          plt.Rectangle((0, 0), 1, 1, fc="#f1f4f7")]
    fig.legend(hs, [r"$b_s = 150$ mm", r"$b_s = 285.5$ mm", r"$\pm 1\,\%$"], ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 1.07), fontsize=6.8, columnspacing=1.6)
    save(fig, "fig05_accuracy")
    return out




# =====================================================================================
# Figure 7 -- critical buckling modes (3D view and skin deflection)
# =====================================================================================
def _mode(lay, bs, beta, bc, nx=91):
    def run():
        m, sec, L = P.case(lay, bs, beta, bc)
        r = m.omega_gbt_solve(sec, L, bc, nev=1, return_mode=True)
        xs = np.linspace(0, L, nx)
        U, u = P.nodal_disp(sec, r["lb"], np.asarray(r["coef"]), xs)
        return dict(L=L, xs=xs, U=U, u=u, N=r["N"], X=sec.X.copy(), S=[(s_["a"], s_["b"], sec.wpart[s_["w"]], s_["n"]) for s_ in sec.S],
                    sides=sec.sides)
    return cached(f"mode_{lay}_{bs}_{beta}_{bc}", run)


def _mode3d(ax, d, amp=22.0, c=0.62, sn=0.36, cmap=None, norm=None):
    """Oblique parallel projection (X + c x, Y + sn x) of the deformed cell, painter's algorithm
    along the projection direction, colour = normal displacement, Lambert shading of the walls."""
    X, U, xs = d["X"], d["U"], d["xs"]
    Um = np.abs(U).max(); Ui = U * amp / Um
    light = np.array([0.35, 1.0, -0.45]); light /= np.linalg.norm(light)
    polys, cols, depth = [], [], []
    for (a, b, part, n) in d["S"]:
        if part == "link":
            continue
        shade = 0.62 + 0.38 * abs(np.r_[n, 0.0] @ light)
        for i in range(len(xs) - 1):
            P4 = []
            for (node, xi) in ((a, i), (b, i), (b, i + 1), (a, i + 1)):
                Xd = X[node] + Ui[xi, node]
                P4.append((Xd[0] + c * xs[xi], Xd[1] + sn * xs[xi]))
            wv = 0.25 * sum(n @ U[xi_, nd] for nd, xi_ in ((a, i), (b, i), (b, i + 1), (a, i + 1))) / Um
            rgb = np.array(cmap(norm(wv))[:3]) * shade
            xm = 0.5 * (xs[i] + xs[i + 1]); Xm = 0.5 * (X[a] + X[b])
            polys.append(P4); cols.append(np.clip(rgb, 0, 1)); depth.append(c * Xm[0] + sn * Xm[1] - xm)
    o = np.argsort(depth)
    pc = PolyCollection([polys[k] for k in o], facecolors=[cols[k] for k in o], edgecolors=[cols[k] for k in o],
                        linewidths=0.2, rasterized=True, zorder=2)
    ax.add_collection(pc)
    # fold lines of the stringers (crown and web edges) and the end sections
    nodes = sorted({nd for (a, b, part, n) in d["S"] if part in ("crown",) for nd in (a, b)
                    if nd in {q for (a2, b2, p2, _n) in d["S"] if p2 == "web" for q in (a2, b2)}})
    for nd in nodes:
        P_ = X[nd] + Ui[:, nd]
        ax.plot(P_[:, 0] + c * xs, P_[:, 1] + sn * xs, color=INK, lw=0.35, zorder=3)
    for xi in (0, len(xs) - 1):
        for (a, b, part, n) in d["S"]:
            if part == "link":
                continue
            p0 = X[a] + Ui[xi, a]; p1 = X[b] + Ui[xi, b]
            ax.plot([p0[0] + c * xs[xi], p1[0] + c * xs[xi]], [p0[1] + sn * xs[xi], p1[1] + sn * xs[xi]],
                    color=INK, lw=0.5, zorder=3 if xi else 4)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(X[:, 0].min() - 10, X[:, 0].max() + c * xs[-1] + 10)
    ax.set_ylim(-amp - 5, X[:, 1].max() + sn * xs[-1] + amp + 5)


def _mode3d_ortho(ax, d, amp=22.0, azim=-35.0, elev=28.0, cmap=None, norm=None, xstride=1, edge_lw=0.55):
    """Orthographic 3D view (azimuth about the vertical Y axis, elevation) of the deformed cell, painter's
    algorithm along the viewing direction, colour = displacement normal to the walls, Lambert shading."""
    X, U, xs = d["X"], d["U"], d["xs"]
    Um = np.abs(U).max(); Ui = U * amp / Um
    a, e = np.radians(azim), np.radians(elev)
    def proj(Xc, Yc, zc):                     # world: X across, Y up, z = x along the stringers
        h = np.cos(a) * Xc + np.sin(a) * zc
        dpt = -np.sin(a) * Xc + np.cos(a) * zc
        return h, np.cos(e) * Yc + np.sin(e) * dpt, -np.sin(e) * Yc + np.cos(e) * dpt
    light = np.array([0.4, 1.0, -0.6]); light /= np.linalg.norm(light)
    idx = list(range(0, len(xs), xstride)) + ([len(xs) - 1] if (len(xs) - 1) % xstride else [])
    polys, cols, depth = [], [], []
    for (pa, pb, part, n) in d["S"]:
        if part == "link":
            continue
        for i0, i1 in zip(idx[:-1], idx[1:]):
            P3 = []
            for (node, xi) in ((pa, i0), (pb, i0), (pb, i1), (pa, i1)):
                Xd = X[node] + Ui[xi, node]
                P3.append((Xd[0], Xd[1], xs[xi]))
            P3 = np.array(P3)
            nrm = np.cross(P3[2] - P3[0], P3[3] - P3[1]); nrm /= max(np.linalg.norm(nrm), 1e-30)
            shade = 0.55 + 0.45 * abs(nrm @ light)
            wv = 0.25 * sum(n @ U[xi_, nd] for nd, xi_ in ((pa, i0), (pb, i0), (pb, i1), (pa, i1))) / Um
            rgb = np.clip(np.array(cmap(norm(wv))[:3]) * shade, 0, 1)
            h, v, dp = proj(P3[:, 0], P3[:, 1], P3[:, 2])
            polys.append(np.c_[h, v]); cols.append(rgb); depth.append(dp.mean())
    o = np.argsort(depth)[::-1]                # far (large depth) first
    ax.add_collection(PolyCollection([polys[k] for k in o], facecolors=[cols[k] for k in o],
                                     edgecolors=[cols[k] for k in o], linewidths=0.15, rasterized=True, zorder=2))
    # outline: end cross-sections and the fold lines of the stringers
    for xi in (0, len(xs) - 1):
        for (pa, pb, part, n) in d["S"]:
            if part == "link":
                continue
            P_ = np.array([X[pa] + Ui[xi, pa], X[pb] + Ui[xi, pb]])
            h, v, _ = proj(P_[:, 0], P_[:, 1], np.full(2, xs[xi]))
            ax.plot(h, v, color=INK, lw=edge_lw, zorder=3, solid_capstyle="round")
    web = {q for (pa, pb, part, n) in d["S"] if part == "web" for q in (pa, pb)}
    fold = sorted({q for (pa, pb, part, n) in d["S"] if part == "crown" for q in (pa, pb)} & web)
    for q in fold:
        P_ = X[q] + Ui[:, q]
        h, v, _ = proj(P_[:, 0], P_[:, 1], xs)
        ax.plot(h, v, color=INK, lw=0.3, alpha=0.8, zorder=3)
    for q in d["sides"]:
        P_ = X[q] + Ui[:, q]
        h, v, _ = proj(P_[:, 0], P_[:, 1], xs)
        ax.plot(h, v, color="#c0392b", lw=0.8, zorder=3)
    ax.set_aspect("equal"); ax.axis("off"); ax.autoscale_view()


def _skinmap(b, Xg, zg, wn, cmap, norm, tri=False):
    """Skin deflection map with contour lines (GBT on its grid, Abaqus on its node triangulation)."""
    lv = [-0.75, -0.5, -0.25, 0.25, 0.5, 0.75]
    if tri:
        import matplotlib.tri as mtri
        T = mtri.Triangulation(Xg, zg)
        pc = b.tripcolor(T, wn, cmap=cmap, norm=norm, shading="gouraud", rasterized=True)
        b.tricontour(T, wn, levels=lv, colors=INK, linewidths=0.3, alpha=0.6)
        b.tricontour(T, wn, levels=[0.0], colors=INK, linewidths=0.45, linestyles="dashed")
    else:
        pc = b.pcolormesh(Xg, zg, wn, cmap=cmap, norm=norm, shading="gouraud", rasterized=True)
        b.contour(Xg, zg, wn, levels=lv, colors=INK, linewidths=0.3, alpha=0.6)
        b.contour(Xg, zg, wn, levels=[0.0], colors=INK, linewidths=0.45, linestyles="dashed")
    for Xs in (-75.0, -225.0):
        b.axvline(Xs, color=MUTED, lw=0.5, ls=(0, (1, 1.5)))
    b.set_aspect("equal"); b.set_xlim(-300, 0)
    b.xaxis.set_major_locator(mt.MultipleLocator(50)); b.yaxis.set_major_locator(mt.MultipleLocator(50))
    return pc


def fig06():
    import abq_modes as AQ
    fig = plt.figure(figsize=(FULL, 168 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.78, 1.0, 1.0], hspace=0.16, wspace=0.10)
    norm = TwoSlopeNorm(0, -1, 1); cmap = plt.get_cmap("RdBu_r")
    out = {}
    for k, bc in enumerate(("SS", "CC")):
        d = _mode("Unsym", 150.0, 1.5, bc)
        sk = [n for n in range(len(d["X"])) if abs(d["X"][n, 1]) < 1e-9]; sk = np.array(sorted(sk, key=lambda n: d["X"][n, 0]))
        w = d["U"][:, sk, 1]
        if w.flat[np.argmax(np.abs(w))] < 0:              # sign: largest skin deflection positive
            d["U"] = -d["U"]; w = -w
        a = fig.add_subplot(gs[0, k]); _mode3d_ortho(a, d, amp=18.0, azim=-38.0, elev=30.0, cmap=cmap, norm=norm)
        panel_label(a, rf"({'ab'[k]}) GBT, {bc}: $N_{{cr}} = {d['N']:.2f}$ N/mm", x=0.0, y=0.98)
        b = fig.add_subplot(gs[1, k])
        wn = w / np.abs(w).max()
        pc = _skinmap(b, d["X"][sk, 0], d["xs"], wn, cmap, norm)
        b.set_ylim(0, d["L"]); b.tick_params(labelbottom=False)
        if k == 0:
            b.set_ylabel(r"$x$ [mm]")
        panel_label(b, rf"({'cd'[k]}) GBT, {bc}: skin deflection")
        # Abaqus mode of the same cell
        ab = AQ.read_abq("MODE1_val_bs150_beta1p5_Unsym_PBC.csv" if bc == "SS" else "MODE1_val_bs150_beta1p5_Unsym_CC_PBC.csv")
        mac, sg = AQ.mac(d, ab)
        mk = ab["inst"] == "INST_SKIN"
        wa = sg * ab["U"][mk, 1]; wa = wa / np.abs(wa).max()
        c = fig.add_subplot(gs[2, k])
        _skinmap(c, ab["xyz"][mk, 0], ab["xyz"][mk, 2], wa, cmap, norm, tri=True)
        c.set_ylim(0, d["L"]); c.set_xlabel(r"$X$ [mm]")
        if k == 0:
            c.set_ylabel(r"$x$ [mm]")
        panel_label(c, rf"({'ef'[k]}) Abaqus, {bc}: $N_{{cr}} = {ab['eig']:.2f}$ N/mm; MAC {mac:.4f}")
        out[bc] = dict(N=d["N"], abq=ab["eig"], mac=mac, edge=float(np.abs(wn[:, 0] - wn[:, -1]).max()))
    cax = fig.add_axes([0.935, 0.12, 0.012, 0.5])
    cb = fig.colorbar(pc, cax=cax); cb.ax.tick_params(labelsize=6.5, length=2); cb.outline.set_linewidth(0.4)
    cb.set_label(r"normal displacement $w/|w|_{\max}$", fontsize=7.5)
    save(fig, "fig06_buckling_modes")
    return out


# =====================================================================================
# Figure 8 -- mode decomposition: families over beta, hierarchical modes of one critical mode
# =====================================================================================
def _hier_mode(lay, bs, beta, bc):
    def run():
        import scipy.sparse as sp
        m, sec, L = P.case(lay, bs, beta, bc)
        hm = m.hierarchical_modes(sec)
        r = m.truncated_section(sec, None, None, None, None, hm=hm)
        full = m.omega_gbt_solve(sec, L, bc)
        res = m.omega_gbt_solve(r, L, bc, pre=full["pre"], return_mode=True)
        C = np.asarray(res["coef"]); lb = res["lb"]
        Kcs = m._gbt_Kcs(r)
        Xb = {(d_, e): (lb.B[d_].T @ sp.diags(lb.wx) @ lb.B[e]).toarray() for d_ in range(3) for e in range(3)}
        KA = np.zeros_like(C)
        for de, A in Kcs.items():
            KA += A @ C @ Xb[de].T
        rho = np.sum(C * KA, axis=1) / np.sum(C * KA)
        xs = np.linspace(0, L, 201)
        amp = P.amplitudes(lb, C, xs)[:, :r.nip]                    # (nx, nip)
        scale = np.array([np.hypot(*r.U[:, k].reshape(-1, 2).T).max() for k in range(r.nip)])
        return dict(N_full=full["N"], N_hier=res["N"], rho=rho, nip=r.nip, o_ax=r.o_ax, o_mx=r.o_mx, lam=hm["C"][0],
                    xs=xs, L=L, amp=amp * scale[None, :], U=r.U.copy(), X=sec.X.copy(),
                    S=[(s_["a"], s_["b"], sec.wpart[s_["w"]]) for s_ in sec.S], sides=sec.sides)
    return cached(f"hier_{lay}_{bs}_{beta}_{bc}", run)


def fig07():
    cur, _ref = _load()
    fig = plt.figure(figsize=(FULL, 150 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.95, 0.42], hspace=0.5)
    top = outer[0].subgridspec(2, 3, height_ratios=[0.42, 1.0], hspace=0.06, wspace=0.12)
    cases = [("Cross-ply", 285.5, "SS"), ("Unsym", 150.0, "SS"), ("Unsym", 150.0, "CC")]
    fams = ["G", "D", "T", "S", "M"]
    out = {}
    for k, (lay, bs, bc) in enumerate(cases):
        d = sorted([c for c in cur if c["layup"] == lay and c["bs"] == bs and c["bc"] == bc], key=lambda c: c["beta"])
        be = [c["beta"] for c in d]
        aL = fig.add_subplot(top[0, k]); aO = fig.add_subplot(top[1, k], sharex=aL)
        aL.plot(be, [100 * c["rho"]["L"] for c in d], "-o", color=FAMC["L"], ms=2.8, mfc="w", mew=0.7, lw=0.9)
        aL.set_ylim(95.5, 101); aL.yaxis.set_major_locator(mt.MultipleLocator(2))
        for g, mk in zip(fams, ["o", "s", "^", "D", "v"]):
            aO.plot(be, [100 * c["rho"].get(g, 0.0) for c in d], "-", marker=mk, color=FAMC[g], ms=2.8, mfc="w", mew=0.7, lw=0.9)
        aO.set_ylim(-2.2, 4.2); aO.axhline(0, color=INK, lw=0.4)
        for a_ in (aL, aO):
            grid(a_)
        aL.spines["bottom"].set_visible(False); aO.spines["top"].set_visible(False)
        aL.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False); aO.tick_params(axis="x", top=False)
        kw = dict(transform=aL.transAxes, color=INK, clip_on=False, lw=0.5)
        aL.plot([-0.015, 0.015], [-0.04, 0.04], **kw); aL.plot([0.985, 1.015], [-0.04, 0.04], **kw)
        kw["transform"] = aO.transAxes
        aO.plot([-0.015, 0.015], [0.985, 1.015], **kw); aO.plot([0.985, 1.015], [0.985, 1.015], **kw)
        aO.set_xlabel(r"$\beta$"); aO.set_xlim(0.4, 3.1)
        if k == 0:
            aL.set_ylabel(r"$\rho_L$ [\%]"); aO.set_ylabel(r"$\rho_g$ [\%]")
        else:
            aL.tick_params(labelleft=False); aO.tick_params(labelleft=False)
        panel_label(aL, rf"({'abc'[k]}) {lay}, $b_s = {bs:g}$ mm, {bc}")
    hs = [Line2D([], [], color=FAMC[g], lw=0.9, marker=mk, ms=3, mfc="w", mew=0.7) for g, mk in zip(["L"] + fams, ["o"] * 2 + ["s", "^", "D", "v"])]
    fig.legend(hs, [FAMN[g] for g in ["L"] + fams], ncol=6, loc="upper center", bbox_to_anchor=(0.5, 0.965), fontsize=6.8,
               columnspacing=1.1, handlelength=1.8)

    # (d) participation of the hierarchical conventional modes, (e) amplitude functions
    h = _hier_mode("Unsym", 150.0, 1.5, "SS")
    out["N_full"], out["N_hier"] = h["N_full"], h["N_hier"]
    rho = 100 * h["rho"]; nshow = 60
    mid = outer[1].subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.22)
    ad = fig.add_subplot(mid[0])
    dom = list(np.argsort(-np.abs(rho[:h["nip"]]))[:4])
    cb_ = [COL["Isotropic"], COL["Cross-ply"], COL["Sym1"], COL["Sym2"]]
    colors = [cb_[dom.index(k)] if k in dom else "#b9c1ca" for k in range(nshow)]
    ad.bar(np.arange(1, nshow + 1), rho[:nshow], width=0.72, color=colors, lw=0, zorder=3)
    for j, k in enumerate(dom):
        ad.text(k + 1, rho[k] + 1.2, rf"{k + 1}", ha="center", va="bottom", fontsize=6.8)
    ad.set_xlim(0, nshow + 1); ad.set_ylim(-2, 56); ad.axhline(0, color=INK, lw=0.4)
    ad.set_xlabel(r"hierarchical conventional mode $k$"); ad.set_ylabel(r"$\rho_k$ [\%]")
    ad.xaxis.set_major_locator(mt.MultipleLocator(10)); ad.xaxis.set_minor_locator(mt.MultipleLocator(1))
    ad.tick_params(axis="x", which="minor", length=1.2)
    grid(ad)
    s20 = rho[:20].sum(); sS = rho[h["o_ax"]:h["o_mx"]].sum(); sM = rho[h["o_mx"]:].sum()
    ad.text(0.98, 0.95, rf"modes 1--20: {s20:.1f}\,\%" + "\n" + rf"all conventional: {rho[:h['nip']].sum():.2f}\,\%" + "\n"
            + rf"shear (incl. axial): {sS:.3f}\,\%" + "\n" + r"Mindlin: 0 (identically)", transform=ad.transAxes, ha="right", va="top", fontsize=6.6,
            linespacing=1.3)
    panel_label(ad, r"(d) participation $\rho_k$ of the hierarchical modes")
    ae = fig.add_subplot(mid[1])
    A = h["amp"]; Am = np.abs(A[:, dom]).max()
    sg = np.sign(A[np.argmax(np.abs(A[:, dom[0]])), dom[0]])        # one common sign: relative phases kept
    for j, k in enumerate(dom):
        ae.plot(h["xs"] / h["L"], sg * A[:, k] / Am, color=cb_[j], lw=1.0, label=rf"mode {k + 1}")
    ae.axhline(0, color=INK, lw=0.4); ae.set_xlim(0, 1); ae.set_ylim(-1.1, 1.5)
    ae.set_xlabel(r"$x/L$"); ae.set_ylabel(r"$\phi_k(x)$ (scaled)")
    ae.set_yticks([-1, -0.5, 0, 0.5, 1])
    grid(ae); lg = ae.legend(ncol=4, loc="upper center", fontsize=6.4, columnspacing=0.8, handlelength=1.2, frameon=True,
                             borderaxespad=0.6, borderpad=0.25)
    lg.get_frame().set(facecolor="white", edgecolor="none", alpha=1.0)
    panel_label(ae, r"(e) amplitude functions of the dominant modes")
    # (f) cross-section shapes of the dominant modes
    bot = outer[2].subgridspec(1, 4, wspace=0.05)
    for j, k in enumerate(dom):
        ax = fig.add_subplot(bot[j]); Uk = h["U"][:, k].reshape(-1, 2); Uk = Uk / np.hypot(*Uk.T).max()
        Xd = h["X"] + 22 * Uk
        segs0 = [[h["X"][a], h["X"][b]] for a, b, p_ in h["S"] if p_ != "link"]
        segs1 = [[Xd[a], Xd[b]] for a, b, p_ in h["S"] if p_ != "link"]
        ax.add_collection(LineCollection(segs0, colors="#c9d0d8", linewidths=0.5))
        ax.add_collection(LineCollection(segs1, colors=cb_[j], linewidths=0.9))
        ax.set_xlim(-305, 5); ax.set_ylim(-26, 56); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(rf"mode {k + 1}: $\rho_k = {rho[k]:.1f}\,\%$", fontsize=7, pad=1)
        out[f"mode{k + 1}"] = float(rho[k])
    fig.text(0.125, outer[2].get_position(fig).y1 + 0.012, r"(f) in-plane shapes of the dominant hierarchical modes", fontsize=8)
    save(fig, "fig07_participation")
    out["sum20"] = float(s20); out["S"] = float(sS); out["M"] = float(sM)
    return out


# =====================================================================================
# Figure 9 -- convergence with the number of hierarchical conventional modes
# =====================================================================================
def fig10():
    hc = json.loads((P.H / "results" / "hier_conv.json").read_text())
    fig, ax = plt.subplots(figsize=(122 * MM, 70 * MM))
    out = {}
    for key, r in hc.items():
        lay, bs, beta, bc = key.split("|")
        nc = sorted(int(k) for k in r["nc"]); dev = [100 * (r["nc"][str(k)] / r["full"] - 1) for k in nc]
        ls = "-" if bc == "SS" else (0, (3.5, 1.6))
        lab = rf"{lay}, {bc}" + (rf", $b_s = {float(bs):g}$ mm" if float(bs) != 150.0 else "")
        ax.plot(nc, dev, ls=ls, color=COL[lay], marker=MK[lay], ms=3.4, mfc="w" if bc == "SS" else COL[lay], mew=0.7, lw=0.9,
                label=lab)
        out[key] = dict(min_dev=min(dev), n_lt1=next((k for k, d_ in zip(nc, dev) if d_ < 1.0), None),
                        n_lt01=next((k for k, d_ in zip(nc, dev) if d_ < 0.1), None))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"number of hierarchical conventional modes $n_C$"); ax.set_ylabel(r"deviation of $N_{cr}$ from the complete basis [\%]")
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.set_xlim(3.5, 230)
    for yv in (1.0, 0.1):
        ax.axhline(yv, color=MUTED, lw=0.5, ls=(0, (1, 1.5)), zorder=1)
    grid(ax)
    ax.legend(loc="upper right", fontsize=6.8, ncol=1, handlelength=2.4)
    save(fig, "fig10_hierarchical")
    return out


# =====================================================================================
# Figures 9, 10 -- long cells: mode transition and the effect of the fixed cell width
# =====================================================================================
def _long():
    """All SS curves of b_s = 150 mm: short cells (periodic_curves) + long cells (long_curves);
    pss = 1: complete pre-buckling state, pss = 0: S_ss removed from G."""
    cur = [json.loads(l) for l in (P.H / "results" / "periodic_curves.jsonl").read_text().splitlines()]
    lng = [json.loads(l) for l in (P.H / "results" / "long_curves.jsonl").read_text().splitlines()]
    out = {}
    for c in cur:
        if c["bs"] == 150.0 and c["bc"] == "SS":
            out[(c["layup"], c["beta"], 1.0)] = c
    for c in lng:
        out[(c["layup"], c["beta"], float(c.get("pss", 1.0)))] = c     # later lines win
    return out


def _abq_long():
    import csv
    f = P.H / "src" / "abaqus_long.csv"          # layup,beta,sides,Ncr  (sides = PBC | PBC_FREE)
    return [dict(layup=r["layup"], beta=float(r["beta"]), sides=r["sides"], N=float(r["Ncr"]))
            for r in csv.DictReader(open(f))] if f.exists() else []


def fig08():
    D = _long(); A = _abq_long()
    fig = plt.figure(figsize=(FULL, 128 * MM))
    gs = fig.add_gridspec(2, 5, height_ratios=[1.35, 0.8], hspace=0.45, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    out = {}
    for lay in LAY:
        for pss, ls in ((1.0, "-"), (0.0, (0, (3.5, 1.6)))):
            d = sorted([v for (l, b, p_), v in D.items() if l == lay and p_ == pss], key=lambda c: c["beta"])
            if len(d) < 2:
                continue
            ax.plot([c["beta"] for c in d], [c["N"] for c in d], ls=ls, color=COL[lay], lw=1.0, marker=MK[lay], ms=3.2,
                    mfc="w" if pss == 1.0 else COL[lay], mew=0.7, zorder=3)
            out[f"{lay}_{int(pss)}"] = [(c["beta"], round(c["N"], 3), {g: round(100 * c["rho"].get(g, 0), 1) for g in "GDL"}) for c in d]
        for r in A:
            if r["layup"] == lay:
                ax.plot(r["beta"], r["N"], "*" if r["sides"] == "PBC" else "X", color=COL[lay], ms=8, mec=INK, mew=0.5, zorder=5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.45, 22)
    for a_ in (ax.xaxis, ax.yaxis):
        a_.set_major_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.xaxis.set_minor_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}" if f"{v:g}"[0] in "25" else ""))
    ax.yaxis.set_minor_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}" if f"{v:g}"[0] in "25" else ""))
    ax.tick_params(which="minor", labelsize=6.3)
    ax.set_xlabel(r"aspect ratio $\beta = L/b_s$"); ax.set_ylabel(r"critical load $N_{cr}$ [N/mm]")
    grid(ax); ax.grid(True, which="minor", color=P.GRID, lw=0.3)
    panel_label(ax, r"(a) simply supported loaded edges, $b_s = 150$ mm")
    hs = [Line2D([], [], color=COL[l], lw=1.0, marker=MK[l], ms=3.5, mfc="w", mew=0.7) for l in LAY] + \
         [Line2D([], [], color=INK, lw=1.0), Line2D([], [], color=INK, lw=1.0, ls=(0, (3.5, 1.6)))]
    lab = LAY + [r"GBT, complete pre-buckling state", r"GBT, $S_{ss}$ omitted"]
    if A:
        hs += [Line2D([], [], ls="none", marker="*", color=MUTED, mec=INK, mew=0.5, ms=8),
               Line2D([], [], ls="none", marker="X", color=MUTED, mec=INK, mew=0.5, ms=7)]
        lab += ["Abaqus, constant width", "Abaqus, free lateral expansion"]
    ax.legend(hs, lab, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.07), fontsize=6.6, columnspacing=1.0, handlelength=2.2)
    for k, lay in enumerate(LAY):
        b = fig.add_subplot(gs[1, k])
        d = sorted([v for (l, bb, p_), v in D.items() if l == lay and p_ == 1.0], key=lambda c: c["beta"])
        be = [c["beta"] for c in d]
        for g, mk in (("L", "o"), ("D", "s"), ("G", "^")):
            b.plot(be, [100 * c["rho"].get(g, 0) for c in d], "-", marker=mk, color=FAMC[g], ms=2.4, mfc="w", mew=0.6, lw=0.8)
        b.set_xscale("log"); b.set_xlim(0.45, 22); b.set_ylim(-5, 105)
        b.xaxis.set_major_formatter(mt.FuncFormatter(lambda v, _p: f"{v:g}"))
        b.set_xlabel(r"$\beta$"); grid(b)
        if k:
            b.tick_params(labelleft=False)
        else:
            b.set_ylabel(r"$\rho_g$ [\%]")
        panel_label(b, rf"({'bcdef'[k]}) {lay}")
    hs = [Line2D([], [], color=FAMC[g], lw=0.8, marker=mk, ms=2.8, mfc="w", mew=0.6) for g, mk in (("L", "o"), ("D", "s"), ("G", "^"))]
    fig.legend(hs, [FAMN["L"], FAMN["D"], FAMN["G"]], ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.025), fontsize=6.8)
    save(fig, "fig08_long_cells")
    return out


def fig09():
    """Critical modes of long cells: S_ss-driven long-wave mode, the same cell without S_ss, global mode."""
    import abq_modes as AQ
    cases = [("Isotropic", 8.0, 1.0, r"(a) Isotropic, $\beta = 8$, complete pre-buckling state", "MODE1_val_bs150_beta8_Isotropic_PBC.csv"),
             ("Isotropic", 8.0, 0.0, r"(b) Isotropic, $\beta = 8$, $S_{ss}$ omitted", None),
             ("Unsym", 12.0, 1.0, r"(c) Unsym, $\beta = 12$, complete pre-buckling state", "MODE1_val_bs150_beta12_Unsym_PBC.csv")]
    fig = plt.figure(figsize=(FULL, 122 * MM))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 0.62], hspace=0.12, wspace=0.06)
    norm = TwoSlopeNorm(0, -1, 1); cmap = plt.get_cmap("RdBu_r")
    out = {}
    for k, (lay, beta, pss, ttl, afile) in enumerate(cases):
        def run(lay=lay, beta=beta, pss=pss):
            m, sec, L = P.case(lay, 150.0, beta, "SS")
            r = m.omega_gbt_solve(sec, L, "SS", hx=10.0 if beta >= 14 else None, nev=1, return_mode=True, pss_scale=pss)
            xs = np.linspace(0, L, 161)
            U, u = P.nodal_disp(sec, r["lb"], np.asarray(r["coef"]), xs)
            return dict(L=L, xs=xs, U=U, u=u, N=r["N"], X=sec.X.copy(),
                        S=[(s_["a"], s_["b"], sec.wpart[s_["w"]], s_["n"]) for s_ in sec.S], sides=sec.sides)
        d = cached(f"longmode_{lay}_{beta}_{pss}", run)
        sk = np.array(sorted([n for n in range(len(d["X"])) if abs(d["X"][n, 1]) < 1e-9], key=lambda n: d["X"][n, 0]))
        if d["U"].flat[np.argmax(np.abs(d["U"]))] < 0:
            d["U"] = -d["U"]
        a = fig.add_subplot(gs[k, 0]); _mode3d_ortho(a, d, amp=40.0, azim=-62.0, elev=24.0, cmap=cmap, norm=norm, xstride=2)
        panel_label(a, ttl + rf": $N_{{cr}} = {d['N']:.2f}$ N/mm", x=0.0, y=0.92)
        # deformed cross-section at the section of largest displacement
        b = fig.add_subplot(gs[k, 1])
        i = int(np.argmax(np.abs(d["U"]).max(axis=(1, 2))))
        Ui = d["U"][i] / np.abs(d["U"][i]).max()
        segs0 = [[d["X"][a_], d["X"][b_]] for a_, b_, p_, _n in d["S"] if p_ != "link"]
        segs1 = [[d["X"][a_] + 20 * Ui[a_], d["X"][b_] + 20 * Ui[b_]] for a_, b_, p_, _n in d["S"] if p_ != "link"]
        b.add_collection(LineCollection(segs0, colors="#c9d0d8", linewidths=0.5))
        b.add_collection(LineCollection(segs1, colors=INK, linewidths=0.9))
        if afile:                                     # Abaqus mode, same section (nearest node plane)
            ab = AQ.read_abq(afile); mac, sg = AQ.mac(d, ab)
            zc = ab["xyz"][np.argmin(np.abs(ab["xyz"][:, 2] - d["xs"][i])), 2]
            sel = np.abs(ab["xyz"][:, 2] - zc) < 1e-3
            Ua = sg * ab["U"][sel][:, :2]; Ua = Ua / np.abs(Ua).max()
            b.plot(ab["xyz"][sel, 0] + 20 * Ua[:, 0], ab["xyz"][sel, 1] + 20 * Ua[:, 1], "o", ms=1.1, mfc="#c0392b",
                   mec="none", zorder=4)
            b.text(0.5, -0.04, rf"line: GBT; dots: Abaqus ({ab['eig']:.2f} N/mm, MAC {mac:.4f})", transform=b.transAxes,
                   ha="center", va="top", fontsize=6.5, color="#c0392b")
            out[f"{lay}_{beta}_abq"] = dict(N=ab["eig"], mac=mac)
        elif pss == 0.0:
            fr = [r for r in AQ.csv.DictReader(open(P.H / "src" / "abaqus_long.csv")) if r["layup"] == lay and r["sides"] == "PBC_FREE"
                  and float(r["beta"]) == beta]
            if fr:
                b.text(0.5, -0.04, rf"Abaqus with free lateral expansion: {float(fr[0]['Ncr']):.2f} N/mm (local mode)",
                       transform=b.transAxes, ha="center", va="top", fontsize=6.5, color="#c0392b")
        b.set_xlim(-310, 10); b.set_ylim(-28, 58); b.set_aspect("equal"); b.axis("off")
        b.set_title(rf"cross-section at $x/L = {d['xs'][i] / d['L']:.2f}$", fontsize=7, pad=1)
        out[f"{lay}_{beta}_{pss}"] = d["N"]
    cax = fig.add_axes([0.16, 0.075, 0.30, 0.010])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label(r"normal displacement $w/|w|_{\max}$", fontsize=7.5); cb.ax.tick_params(labelsize=6.5, length=2)
    cb.outline.set_linewidth(0.4)
    save(fig, "fig09_long_modes")
    return out


if __name__ == "__main__":
    want = sys.argv[1:] or sorted(k for k in dir() if k.startswith("fig") and k[3:].isdigit())
    fd = P.H / "results" / "figdata.json"
    allr = json.loads(fd.read_text()) if fd.exists() else {}
    for k in want:
        r = globals()[k]()
        if r is not None:
            allr[k] = r
            print(k, json.dumps(r, default=float)[:1500])
    fd.write_text(json.dumps(allr, indent=1, default=float))
