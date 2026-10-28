"""Periodic-only result figures: N_cr over beta (SS with Abaqus, CC with spline FSM),
energy-based mode decomposition over beta, and the buckling modes (SS, CC)."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import figs as F   # shared helpers (section, skin_field, styles)

H = Path(__file__).parent
OUT = H / "fig"
LAY = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
COL = {"Isotropic": "#1F3A5F", "Cross-ply": "#2E7D32", "Sym1": "#B26A00", "Sym2": "#B3261E", "Unsym": "#6A3D9A"}
MK = {"Isotropic": "o", "Cross-ply": "s", "Sym1": "^", "Sym2": "D", "Unsym": "v"}
FAM = [("L", "local (L)", "#1F3A5F"), ("D", "distortional (D)", "#B3261E"), ("G", "global (G)", "#2E7D32"),
       ("T", "transverse ext. (T)", "#B26A00"), ("S", "shear (S)", "#6A3D9A"), ("M", "Mindlin (M)", "#7A7A7A")]


def load():
    cur = [json.loads(l) for l in (H / "results" / "periodic_curves.jsonl").read_text().splitlines()]
    ref = list(csv.DictReader(open(H / "src" / "omega_results_V32.csv")))
    for r in ref:
        for k in ("bs", "beta", "gbt_per", "strip_per"):
            r[k] = float(r[k])
        r["abq_per"] = float(r["abq_per"]) if r["abq_per"] not in ("", None) else None
    return cur, ref


def fig_curves(cur, ref, bc, fname):
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
    for a, bs in zip(ax, (150.0, 285.5)):
        for lay in LAY:
            d = sorted([c for c in cur if c["layup"] == lay and c["bs"] == bs and c["bc"] == bc], key=lambda c: c["beta"])
            if not d:
                continue
            a.plot([c["beta"] for c in d], [c["N"] for c in d], "-", color=COL[lay], lw=1.1, label=f"{lay} (GBT)")
            rr = [r for r in ref if r["layup"] == lay and r["bs"] == bs and r["bc"] == bc]
            if bc == "SS":
                a.plot([r["beta"] for r in rr], [r["abq_per"] for r in rr], MK[lay], mfc="none", mec=COL[lay], ms=5.5, mew=0.9)
            else:
                a.plot([r["beta"] for r in rr], [r["strip_per"] for r in rr], MK[lay], mfc="none", mec=COL[lay], ms=5.5, mew=0.9)
        a.set_yscale("log"); a.set_xlabel(r"aspect ratio $\beta = L/b_s$")
        a.set_title(f"({'a' if bs == 150 else 'b'}) $b_s$ = {bs:g} mm", fontsize=8.5)
        a.grid(alpha=0.25, which="both")
    ax[0].set_ylabel(r"$N_{cr}$ [N/mm]")
    ref_lbl = "Abaqus (markers)" if bc == "SS" else "spline FSM (markers)"
    h, l = ax[0].get_legend_handles_labels()
    h.append(plt.Line2D([], [], ls="none", marker="o", mfc="none", mec="k")); l.append(ref_lbl)
    fig.legend(h, l, loc="lower center", ncol=6, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT / fname, dpi=300, bbox_inches="tight"); plt.close(fig)


def fig_decomp(cur, fname):
    cases = [("Cross-ply", "SS"), ("Unsym", "SS"), ("Sym2", "CC"), ("Unsym", "CC")]
    fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.4), sharey=True)
    for k, (a, (lay, bc)) in enumerate(zip(ax, cases)):
        d = sorted([c for c in cur if c["layup"] == lay and c["bs"] == 150.0 and c["bc"] == bc], key=lambda c: c["beta"])
        b = [c["beta"] for c in d]
        for g, lbl, col in FAM:
            a.plot(b, [100 * c["rho"].get(g, 0.0) for c in d], "-o", ms=2.5, lw=1.0, color=col, label=lbl)
        a.set_title(f"({'abcd'[k]}) {lay}, {bc}", fontsize=8.5); a.set_xlabel(r"$\beta$")
        a.grid(alpha=0.25); a.set_ylim(-5, 105)
    ax[0].set_ylabel(r"participation $\rho_g$ [%]")
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT / fname, dpi=300, bbox_inches="tight"); plt.close(fig)


def fig_modes_ss_cc(fname):
    out = {}
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.3), sharey=True)
    for a, bc in zip(ax, ("SS", "CC")):
        X, z, w, N, part, sec = F.skin_field("Unsym", 150.0, 1.5, bc, "periodic")
        if np.mean(w) < 0:
            w = -w
        pc = a.pcolormesh(X, z, w, cmap="RdBu_r", vmin=-1, vmax=1, shading="gouraud", rasterized=True)
        a.set_title(f"({'a' if bc == 'SS' else 'b'}) {bc}: $N_{{cr}}$ = {N:.2f} N/mm", fontsize=8.5)
        a.set_xlabel(r"$X$ [mm]"); a.set_aspect("equal")
        out[bc] = dict(N=N, edge_diff=float(np.abs(w[:, 0] - w[:, -1]).max()))
    ax[0].set_ylabel(r"$x$ [mm]")
    cb = fig.colorbar(pc, ax=ax, shrink=0.85, pad=0.02); cb.set_label(r"skin deflection $w/|w|_{max}$")
    fig.savefig(OUT / fname, dpi=300, bbox_inches="tight"); plt.close(fig)
    return out


if __name__ == "__main__":
    cur, ref = load()
    fig_curves(cur, ref, "SS", "fig3_ncr_ss.png")
    fig_curves(cur, ref, "CC", "fig4_ncr_cc.png")
    print(fig_modes_ss_cc("fig5_buckling_modes.png"))
