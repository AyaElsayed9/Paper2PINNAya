"""Modes of the n-pitch cells (11 Nov 2026): which transverse pattern (Bloch wavenumber) does every mode have?
Measured, not assumed: for a real mode u of a cell of n pitches, c = <u(Y + b_s), u(Y)> / <u, u> = cos(k_B b_s)
(skin deflection U_Z on a uniform grid over the cell, shifted by one pitch, all sections).  c = -1: alternating (pi); -0.5: 2pi/3;
0: pi/2; +1: the same in every bay (0).   -> results/npitch_modes.json, fig/check_npitch_modes.png"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
H = Path(__file__).parent; sys.path.insert(0, str(H))
import fp_common as P
from run_npitch import npitch_panel, configure, validated_omega_sample, omega_panel

BS, BETA, HX = 150.0, 8.0, 15.0


def modes(lay, n, nev=4):
    with configure(OMEGA_SIDE_BC="periodic") as m:
        pan2 = m.make_omega_panel(BS, "iso", [], []) if lay == "Isotropic" else omega_panel(
            validated_omega_sample(lay, BS, BETA, "SS", "periodic"))
        sec = m.omega_gbt_section(npitch_panel(m, pan2, n)); L = BETA * BS
        r = m.omega_gbt_solve(sec, L, "SS", hx=HX, nev=nev, return_mode=True)
        xs = np.linspace(0, L, 121); X = sec.X.copy(); W = n * BS
        # skin deflection U_Z(Y) on a uniform grid over the cell (linear interpolation between the skin nodes)
        sk = np.array(sorted([k for k in range(len(X)) if abs(X[k, 1]) < 1e-9], key=lambda k: X[k, 0]))
        ng = 60 * n; Yg = -W + W * np.arange(ng) / ng
        out = []
        for j in range(nev):
            U, _u = P.nodal_disp(sec, r["lb"], np.asarray(r["coefs"][j]), xs)       # (nx, N, 2)
            w = np.array([np.interp(Yg, X[sk, 0], U[i, sk, 1]) for i in range(len(xs))])     # (nx, ng)
            c = float(np.sum(np.roll(w, -60, axis=1) * w) / np.sum(w * w))
            out.append(dict(N=float(r["eigs"][j]), c=c, U=U))
        S = [(s_["a"], s_["b"]) for s_ in sec.S if sec.wpart[s_["w"]] != "link"]
    return dict(X=X, S=S, xs=xs, modes=out, W=W)


def name(c):
    for v, s in ((-1, r"alternating ($k_B b_s=\pi$)"), (-0.5, r"every 3 bays ($2\pi/3$)"), (0, r"every 4 bays ($\pi/2$)"),
                 (1, r"same in every bay ($0$)")):
        if abs(c - v) < 0.05:
            return s
    return "mixed"


if __name__ == "__main__":
    import pickle, subprocess, tempfile
    if len(sys.argv) == 4:                       # one case per process (memory)
        d = modes(sys.argv[1], int(sys.argv[2])); pickle.dump(d, open(sys.argv[3], "wb")); sys.exit(0)
    res = {}; keep = None; tmp = tempfile.mkdtemp()
    for lay in ("Isotropic", "Sym1"):
        for n in (2, 3, 4):
            f = f"{tmp}/{lay}_{n}.pkl"
            subprocess.run([sys.executable, __file__, lay, str(n), f], check=True)
            d = pickle.load(open(f, "rb"))
            res[f"{lay}_n{n}"] = [dict(N=m_["N"], c=m_["c"], pattern=name(m_["c"])) for m_ in d["modes"]]
            print(lay, "pitches", n, [(round(m_["N"], 3), round(m_["c"], 3)) for m_ in d["modes"]], flush=True)
            if lay == "Isotropic" and n == 4:
                keep = d
    json.dump(res, open(H / "results" / "npitch_modes.json", "w"), indent=1)
    d = keep; X = d["X"]
    fig, axs = plt.subplots(4, 1, figsize=(7.2, 7.4))
    for ax, m_, k in zip(axs, d["modes"], range(4)):
        U = m_["U"]; i = int(np.argmax(np.abs(U).max(axis=(1, 2)))); Ui = U[i] / np.abs(U[i]).max()
        ax.add_collection(LineCollection([[X[a], X[b]] for a, b in d["S"]], colors="#c9d0d8", linewidths=0.6))
        ax.add_collection(LineCollection([[X[a] + 22 * Ui[a], X[b] + 22 * Ui[b]] for a, b in d["S"]], colors="#1f2a36", linewidths=1.1))
        for y in np.arange(-d["W"], 1, BS):
            ax.axvline(y, color="#c0392b", lw=0.5, ls=":")
        ax.set_xlim(-d["W"] - 10, 10); ax.set_ylim(-32, 62); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"mode {k + 1}: N = {m_['N']:.2f} N/mm,  shift-by-one-pitch value c = {m_['c']:+.3f}:  {name(m_['c'])}",
                     fontsize=8.5, loc="left")
    fig.suptitle("Isotropic, $b_s$ = 150 mm, $\\beta$ = 8, SS: 4-pitch cell, cross-section at the largest displacement\n"
                 "(grey: undeformed, dark: deformed, dotted: one pitch)", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(H / "fig" / "check_npitch_modes.png", dpi=200)
