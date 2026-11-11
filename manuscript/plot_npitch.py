"""Modes of the n-pitch cells (11 Nov 2026): which transverse pattern (Bloch wavenumber) does every mode have?
Measured, not assumed: for a real mode u of a cell of n pitches, c = <u(Y + b_s), u(Y)> / <u, u> = cos(k_B b_s)
(skin deflection U_Z on a uniform grid over the cell, shifted by one pitch, all sections).  c = -1: alternating (pi); -0.5: 2pi/3;
0: pi/2; +1: the same in every bay (0).   -> results/npitch_modes.json, fig/check_npitch_<lay-up>.png/.pdf (style and drawing
functions of the paper figures, cf. Figure 9; every number in the figure is read from the computed modes)"""
import json, sys
from pathlib import Path
import numpy as np
H = Path(__file__).parent; sys.path.insert(0, str(H))
import fp_common as P
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from run_npitch import npitch_panel, configure, validated_omega_sample, omega_panel

BS, BETA, HX = 150.0, 8.0, 15.0


def modes(lay, n, nev=4):
    with configure(OMEGA_SIDE_BC="periodic") as m:
        pan2 = m.make_omega_panel(BS, "iso", [], []) if lay == "Isotropic" else omega_panel(
            validated_omega_sample(lay, BS, BETA, "SS", "periodic"))
        sec = m.omega_gbt_section(npitch_panel(m, pan2, n)); L = BETA * BS
        r = m.omega_gbt_solve(sec, L, "SS", hx=HX, nev=nev, return_mode=True)
        xs = np.linspace(0, L, 161); X = sec.X.copy(); W = n * BS
        # skin deflection U_Z(Y) on a uniform grid over the cell (linear interpolation between the skin nodes)
        sk = np.array(sorted([k for k in range(len(X)) if abs(X[k, 1]) < 1e-9], key=lambda k: X[k, 0]))
        ng = 60 * n; Yg = -W + W * np.arange(ng) / ng
        out = []
        for j in range(nev):
            U, _u = P.nodal_disp(sec, r["lb"], np.asarray(r["coefs"][j]), xs)       # (nx, N, 2)
            w = np.array([np.interp(Yg, X[sk, 0], U[i, sk, 1]) for i in range(len(xs))])     # (nx, ng)
            c = float(np.sum(np.roll(w, -60, axis=1) * w) / np.sum(w * w))
            out.append(dict(N=float(r["eigs"][j]), c=c, U=U))
        S = [(s_["a"], s_["b"], sec.wpart[s_["w"]], s_["n"]) for s_ in sec.S]
    return dict(X=X, S=S, xs=xs, modes=out, W=W, L=L, sides=sec.sides, lay=lay, n=n)


def figure(d):
    """Paper style (figures_paper.fig09): deformed cell and cross-section at the largest displacement, one row per
    computed mode.  Titles: eigenvalue, measured c and k_B b_s = arccos(c) of that mode."""
    import figures_paper as FP
    P.apply_style()
    nm = len(d["modes"])
    fig = plt.figure(figsize=(P.FULL, (16 + 36 * nm) * P.MM))
    gs = fig.add_gridspec(nm, 2, width_ratios=[1.0, 0.95], hspace=0.16, wspace=0.05)
    norm = TwoSlopeNorm(0, -1, 1); cmap = plt.get_cmap("RdBu_r"); X = d["X"]
    for k, m_ in enumerate(d["modes"]):
        U = m_["U"].copy()
        if U.flat[np.argmax(np.abs(U))] < 0:
            U = -U
        dd = dict(d, U=U)
        a = fig.add_subplot(gs[k, 0]); FP._mode3d_ortho(a, dd, amp=40.0, azim=-62.0, elev=24.0, cmap=cmap, norm=norm, xstride=2)
        kb = np.arccos(np.clip(m_["c"], -1, 1)) / np.pi
        P.panel_label(a, rf"({chr(97 + k)}) mode {k + 1}: $N = {m_['N']:.2f}$ N/mm, $c = {m_['c']:+.3f}$, "
                         rf"$k_B b_s = {kb:.2f}\,\pi$", x=0.0, y=0.95)
        b = fig.add_subplot(gs[k, 1])
        i = int(np.argmax(np.abs(U).max(axis=(1, 2)))); Ui = U[i] / np.abs(U[i]).max()
        b.add_collection(LineCollection([[X[p], X[q]] for p, q, part, _n in d["S"] if part != "link"], colors="#c9d0d8", linewidths=0.5))
        b.add_collection(LineCollection([[X[p] + 20 * Ui[p], X[q] + 20 * Ui[q]] for p, q, part, _n in d["S"] if part != "link"],
                                        colors=FP.INK, linewidths=0.9))
        b.set_xlim(-d["W"] - 10, 10); b.set_ylim(-28, 58); b.set_aspect("equal"); b.axis("off")
        b.set_title(rf"cross-section at $X/L = {d['xs'][i] / d['L']:.2f}$", fontsize=7, pad=1)
    cax = fig.add_axes([0.16, 0.045, 0.30, 0.008])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label(r"normal displacement $w/|w|_{\max}$", fontsize=7.5); cb.ax.tick_params(labelsize=6.5, length=2)
    cb.outline.set_linewidth(0.4)
    P.save(fig, f"check_npitch_{d['lay']}_n{d['n']}")


if __name__ == "__main__":
    import pickle, subprocess
    if len(sys.argv) == 4:                       # one case per process (memory)
        d = modes(sys.argv[1], int(sys.argv[2])); pickle.dump(d, open(sys.argv[3], "wb")); sys.exit(0)
    res = {}
    for lay in ("Isotropic", "Sym1"):
        for n in (2, 3, 4):
            f = P.CACHE / f"npitch_{lay}_{n}.pkl"
            if not f.exists():
                subprocess.run([sys.executable, __file__, lay, str(n), str(f)], check=True)
            d = pickle.load(open(f, "rb"))
            res[f"{lay}_n{n}"] = [dict(N=m_["N"], c=m_["c"], kB_bs_over_pi=float(np.arccos(np.clip(m_["c"], -1, 1)) / np.pi))
                                  for m_ in d["modes"]]
            print(lay, "pitches", n, [(round(m_["N"], 3), round(m_["c"], 3)) for m_ in d["modes"]], flush=True)
            if n == 4:
                figure(d)
    json.dump(res, open(H / "results" / "npitch_modes.json", "w"), indent=1)
