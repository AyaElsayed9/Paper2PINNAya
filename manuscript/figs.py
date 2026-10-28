"""Figures of the manuscript (data: solver V32 and the V32 results table)."""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample          # noqa: E402
from pinnbuck.oracle import omega_length, omega_panel             # noqa: E402
from pinnbuck.solver import configure                             # noqa: E402

OUT = Path(__file__).parent / "fig"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 8.5, "axes.linewidth": 0.6,
                     "mathtext.fontset": "dejavuserif"})
PART_COL = {"skin": "#1F3A5F", "foot": "#2E7D32", "web": "#B26A00", "crown": "#B3261E", "link": "#7A7A7A"}
INK = "#1F3A5F"


def section(lay, bs, beta, bc, side):
    s = validated_omega_sample(lay, bs, beta, bc, side)
    with configure(OMEGA_SIDE_BC=side) as m:
        sec = m.omega_gbt_section(omega_panel(s))
    return m, sec, s


def draw_section(ax, sec, U=None, scale=0.0, lw=1.2, nodes=True, col=None, alpha=1.0):
    X = sec.X.copy()
    if U is not None:
        X = X + scale * U.reshape(-1, 2)
    for sg in sec.S:
        a, b, w = sg["a"], sg["b"], sg["w"]
        c = col or PART_COL[sec.wpart[w]]
        ax.plot([X[a, 0], X[b, 0]], [X[a, 1], X[b, 1]], color=c, lw=lw, alpha=alpha, solid_capstyle="round")
    if nodes:
        ax.plot(X[:, 0], X[:, 1], "o", ms=1.6, color="k", alpha=alpha)
    return X


# ---------------------------------------------------------------- Fig. 1: geometry and periodic cell
def fig_geometry():
    m, sec, s = section("Cross-ply", 150.0, 1.5, "SS", "periodic")
    W = np.ptp(sec.X[:, 0])
    fig, ax = plt.subplots(2, 1, figsize=(7.0, 3.6), gridspec_kw=dict(height_ratios=[1, 1.35]))
    a = ax[0]
    for k in (-1, 0, 1):
        for sg in sec.S:
            A_, B_ = sec.X[sg["a"]], sec.X[sg["b"]]
            a.plot([A_[0] + k * W, B_[0] + k * W], [A_[1], B_[1]], color=INK if k == 0 else "#9AA7B8", lw=1.1 if k == 0 else 0.8)
    x0, x1 = sec.X[:, 0].min(), sec.X[:, 0].max()
    a.axvspan(x0, x1, color="#EEF3F8", zorder=-1)
    for xx in (x0, x1):
        a.axvline(xx, color="k", lw=0.5, ls="--")
    a.annotate("", xy=(x1, 44), xytext=(x0, 44), arrowprops=dict(arrowstyle="<->", lw=0.6))
    a.text(0.5 * (x0 + x1), 47, r"periodic cell, width $W = 2b_s$", ha="center", va="bottom", fontsize=8)
    a.text(x0 - W * 1.08, 2, r"$\cdots$", ha="right", va="center", fontsize=12); a.text(x1 + W * 1.08, 2, r"$\cdots$", ha="left", va="center", fontsize=12)
    a.set_xlim(x0 - W * 1.2, x1 + W * 1.2); a.set_ylim(-6, 58); a.set_aspect("equal"); a.axis("off")
    a.set_title("(a) infinitely wide omega-stringer-stiffened panel and its periodic cell", fontsize=8.5, loc="left")
    b = ax[1]
    draw_section(b, sec, lw=1.4)
    sL, sR = sec.sides
    for i, lab in ((sL, r"$s_L$"), (sR, r"$s_R$")):
        b.plot(*sec.X[i], "s", ms=5, mfc="w", mec="k", mew=0.9)
        b.text(sec.X[i, 0], -7, lab, ha="center", va="top", fontsize=9)
    b.annotate("", xy=(sec.X[sR, 0], -15), xytext=(sec.X[sL, 0], -15),
               arrowprops=dict(arrowstyle="<->", lw=0.6, ls="--"))
    b.text(0.5 * (sec.X[sL, 0] + sec.X[sR, 0]), -14, r"$\mathbf{u}(s_L,x)=\mathbf{u}(s_R,x)$ for every mode",
           ha="center", va="bottom", fontsize=8, bbox=dict(fc="w", ec="none", pad=0.5))
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=c, lw=2) for c in PART_COL.values()]
    b.legend(hs, ["skin", "foot (skin + flange)", "web", "crown", "link wall (offset tie)"], ncol=5, fontsize=7,
             loc="upper center", bbox_to_anchor=(0.5, 1.13), frameon=False)
    b.set_aspect("equal"); b.set_ylim(-22, 45); b.axis("off")
    b.text(0.0, 0.0, "", transform=b.transAxes)
    b.set_title(r"(b) GBT discretisation of the cell ($b_s = 150$ mm): natural and intermediate nodes", fontsize=8.5,
                loc="left", y=1.12)
    fig.tight_layout(); fig.savefig(OUT / "fig1_geometry.png", dpi=300); plt.close(fig)
    return sec


# ---------------------------------------------------------------- Fig. 2: hierarchical periodic modes
def fig_modes():
    m, sec, s = section("Cross-ply", 150.0, 1.5, "SS", "periodic")
    lam, T = m.hierarchical_modes(sec, kref=2 * np.pi / 150.0)["C"]
    o = np.argsort(lam)
    fig, axs = plt.subplots(4, 2, figsize=(7.0, 4.3))
    for ax, i in zip(axs.ravel(), o[:8]):
        U = sec.U @ T[:, i]
        sc = 18.0 / max(np.abs(U).max(), 1e-12)
        draw_section(ax, sec, col="#C9D2DD", nodes=False, lw=0.8)
        Xd = draw_section(ax, sec, U=U, scale=sc, col=INK, nodes=False, lw=1.1)
        sL, sR = sec.sides
        ax.plot(*Xd[sL], "s", ms=3.5, mfc="w", mec="#B3261E", mew=0.9); ax.plot(*Xd[sR], "s", ms=3.5, mfc="w", mec="#B3261E", mew=0.9)
        k = list(o).index(i) + 1
        ax.set_title(rf"mode {k}:  $\Lambda$ = {lam[i]:.3g}", fontsize=8)
        ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "fig2_modes.png", dpi=300); plt.close(fig)
    # side-edge check for every conventional mode of the basis
    sL, sR = sec.sides
    dU = np.abs(sec.U[2 * sL:2 * sL + 2] - sec.U[2 * sR:2 * sR + 2]).max(); dT = np.abs(sec.Th[sL] - sec.Th[sR]).max()
    dW = np.abs(sec.Uw[sL] - sec.Uw[sR]).max()
    return dict(dU=float(dU), dTh=float(dT), dW=float(dW), counts={g: sec.group.count(g) for g in set(sec.group)})


# ---------------------------------------------------------------- Fig. 3: validation
def fig_validation(rows):
    LAY = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
    MK = dict(zip(LAY, ["o", "s", "^", "D", "v"]))
    CO = dict(zip(LAY, ["#1F3A5F", "#2E7D32", "#B26A00", "#B3261E", "#6A3D9A"]))
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.55))
    a = ax[0]
    for r in rows:
        e = 100 * (r["gbt_sym"] / r["abq_sym"] - 1)
        a.plot(r["abq_sym"], e, MK[r["layup"]], mfc="none" if r["bc"] == "SS" else CO[r["layup"]], mec=CO[r["layup"]], ms=4, mew=0.8)
    a.set_xscale("log"); a.axhline(0, color="k", lw=0.5); a.set_ylim(-2.2, 1.2)
    a.set_xlabel(r"Abaqus $N_{cr}$ [N/mm]"); a.set_ylabel("GBT vs Abaqus [%]")
    a.set_title("(a) symmetry sides, 60 cases", fontsize=8.5)
    b = ax[1]
    for r in rows:
        if r["abq_per"] == "":
            continue
        e = 100 * (r["gbt_per"] / float(r["abq_per"]) - 1)
        b.plot(float(r["abq_per"]), e, MK[r["layup"]], mfc="none", mec=CO[r["layup"]], ms=4, mew=0.8)
    b.set_xscale("log"); b.axhline(0, color="k", lw=0.5); b.set_ylim(-2.2, 1.2)
    b.set_xlabel(r"Abaqus $N_{cr}$ [N/mm]"); b.set_title("(b) periodic sides, 30 SS cases", fontsize=8.5)
    c = ax[2]
    for r in rows:
        if r["abq_per"] == "":
            continue
        sa = 100 * (float(r["abq_per"]) / r["abq_sym"] - 1); sg = 100 * (r["gbt_per"] / r["gbt_sym"] - 1)
        c.plot(sa, sg, MK[r["layup"]], mfc="none", mec=CO[r["layup"]], ms=4, mew=0.8)
    c.plot([-0.3, 2.3], [-0.3, 2.3], "k-", lw=0.5)
    c.set_xlabel("Abaqus shift [%]"); c.set_ylabel("GBT shift [%]")
    c.set_title("(c) shift periodic − symmetry", fontsize=8.5)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls="", marker=MK[l], mec=CO[l], mfc="none", ms=4) for l in LAY]
    fig.legend(hs, LAY, ncol=5, loc="lower center", fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.27, top=0.9, wspace=0.5); fig.savefig(OUT / "fig3_validation.png", dpi=300); plt.close(fig)


def fig_cc_shift(rows):
    LAY = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.4))
    for j, lay in enumerate(LAY):
        for k, bc in enumerate(("SS", "CC")):
            g = [100 * (r["gbt_per"] / r["gbt_sym"] - 1) for r in rows if r["layup"] == lay and r["bc"] == bc]
            st = [100 * (r["strip_per"] / r["strip_sym"] - 1) for r in rows if r["layup"] == lay and r["bc"] == bc]
            x = j + (k - 0.5) * 0.36
            col = "#9AA7B8" if bc == "SS" else "#1F3A5F"
            ax.plot([x, x], [min(g), max(g)], color=col, lw=5, solid_capstyle="butt")
            ax.plot([x + 0.08] * len(st), st, "x", color="#B3261E", ms=3, mew=0.7)
    ax.set_xticks(range(5)); ax.set_xticklabels(LAY, fontsize=7.5)
    ax.set_ylabel("periodic − symmetry [%]"); ax.axhline(0, color="k", lw=0.5)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([], [], color="#9AA7B8", lw=5), Line2D([], [], color="#1F3A5F", lw=5),
               Line2D([], [], ls="", marker="x", color="#B3261E")], ["GBT, SS (range)", "GBT, CC (range)", "spline FSM"],
              fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig4_shift_ss_cc.png", dpi=300); plt.close(fig)


# ---------------------------------------------------------------- Fig. 5: buckling modes, symmetry vs periodic
def skin_field(lay, bs, beta, bc, side):
    s = validated_omega_sample(lay, bs, beta, bc, side)
    with configure(OMEGA_SIDE_BC=side) as m:
        sec = m.omega_gbt_section(omega_panel(s)); L = omega_length(s)
        r = m.omega_gbt_solve(sec, L, bc, nev=1, return_mode=True)
    C = np.asarray(r["coefs"][0]); lb = r["lb"]
    V = (lb.B[0] @ C[:sec.nip].T)                     # (nz, nip) conventional amplitudes at the Gauss points
    skin = np.nonzero(np.abs(sec.X[:, 1]) < 1e-9)[0]; skin = skin[np.argsort(sec.X[skin, 0])]
    UY = V @ sec.U[2 * skin + 1].T                     # (nz, n_skin)
    UY = UY / np.abs(UY).max()
    return sec.X[skin, 0], np.asarray(lb.xg), UY, r["N"], r["participation"], sec


def fig_buckling():
    out = {}
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.3), sharey=True)
    for a, side in zip(ax, ("symmetry", "periodic")):
        X, z, w, N, part, sec = skin_field("Unsym", 150.0, 1.5, "SS", side)
        if np.mean(w) < 0:
            w = -w
        pc = a.pcolormesh(X, z, w, cmap="RdBu_r", vmin=-1, vmax=1, shading="gouraud", rasterized=True)
        a.set_title(f"({'a' if side == 'symmetry' else 'b'}) {side}: $N_{{cr}}$ = {N:.2f} N/mm", fontsize=8.5)
        a.set_xlabel(r"$X$ [mm]"); a.set_aspect("equal")
        out[side] = dict(N=N, part=part, edge_diff=float(np.abs(w[:, 0] - w[:, -1]).max()))
    ax[0].set_ylabel(r"$x$ [mm]")
    cb = fig.colorbar(pc, ax=ax, shrink=0.85, pad=0.02); cb.set_label(r"skin deflection $w/|w|_{max}$")
    fig.savefig(OUT / "fig5_buckling_modes.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    return out


plt.rcParams.update({"xtick.direction": "out", "ytick.direction": "out", "axes.labelsize": 8.5, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5})

if __name__ == "__main__":
    rows = list(csv.DictReader(open(Path(__file__).parent / "src" / "omega_results_V32.csv")))
    for r in rows:
        for k in ("bs", "beta", "L", "abq_sym", "gbt_sym", "strip_sym", "err_sym", "gbt_per", "strip_per"):
            r[k] = float(r[k])
    fig_geometry()
    print("modes:", fig_modes())
    fig_validation(rows)
    fig_cc_shift(rows)
    print("buckling:", fig_buckling())
