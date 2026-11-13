"""Figure 10 (11 Nov 2026): buckling modes of the four-pitch cells, b_s = 150 mm, beta = 8, SS.
Data: GBT modes of plot_npitch.py (figp/npitch_<lay>_4.pkl, spans 15 mm); Bloch wavenumber of every mode from the
measured c = cos(k_B b_s) (skin deflection shifted by one pitch); Abaqus eigenvalues of the four-pitch shell models
(src/abaqus_npitch/omega_validation_eigs_PBC_n4.csv, Aya's run 11 Nov 2026).
Also writes results/npitch_table.json (the values of Table 6)."""
import csv
import json
import pickle
import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import fp_common as P
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

H = Path(__file__).parent
P.apply_style()
abq = {r["layup"]: [float(r[f"m{k}_eig"]) for k in range(1, 7)]
       for r in csv.DictReader(open(H / "src" / "abaqus_npitch" / "omega_validation_eigs_PBC_n4.csv"))}
cell = {(lay, n): pickle.load(open(P.CACHE / f"npitch_{lay}_{n}.pkl", "rb")) for lay in ("Isotropic", "Sym1") for n in (2, 3, 4)}
AQ_CHECK = json.load(open(H / "results" / "npitch_abq_modes.json"))
KB = [(-1.0, r"\pi"), (-0.5, r"2\pi/3"), (0.0, r"\pi/2"), (1.0, r"0")]


def kb_of(c):
    return min(KB, key=lambda t: abs(t[0] - c))[1]


# ---- table data: for every lay-up and wavenumber, the GBT load (lowest mode with that wavenumber) in each cell
table = {}
for lay in ("Isotropic", "Sym1"):
    rows = {}
    for n in (2, 3, 4):
        for m in cell[(lay, n)]["modes"]:
            k = kb_of(m["c"]); assert abs(np.cos({r"\pi": np.pi, r"2\pi/3": 2 * np.pi / 3, r"\pi/2": np.pi / 2, r"0": 0}[k]) - m["c"]) < 1e-3
            rows.setdefault(k, {}).setdefault(n, []).append(m["N"])
    # Abaqus four-pitch eigenvalues 1-4: wavenumbers pi, pi/2 (double), 0, measured from the Abaqus modes (13 Nov 2026)
    assert [kb_of(q_["c_abaqus"]) for q_ in AQ_CHECK[lay]] == [r"\pi", r"\pi/2", r"\pi/2", r"0"]
    g4 = [kb_of(m["c"]) for m in cell[(lay, 4)]["modes"]]
    assert g4 == [r"\pi", r"\pi/2", r"\pi/2", r"0"], g4
    a4 = abq[lay][:4]; assert abs(a4[1] - a4[2]) < 1e-9 * a4[1]
    table[lay] = {k: dict(cells=sorted(v), N_gbt=min(min(x) for x in v.values()), N_abq=None) for k, v in rows.items()}
    table[lay][r"\pi"]["N_abq"] = a4[0]; table[lay][r"\pi/2"]["N_abq"] = a4[1]; table[lay][r"0"]["N_abq"] = a4[3]
json.dump(table, open(H / "results" / "npitch_table.json", "w"), indent=1)

# ---- figure: rows = Bloch wavenumber (pi, pi/2, 0), columns = Isotropic | Sym1 (deformed cell + cross-section),
#      drawn like Figure 9 (figures_paper._mode3d_ortho)
import figures_paper as FP
from matplotlib.colors import TwoSlopeNorm
P.apply_style()
FW, FH = 160.0, 126.0
fig = plt.figure(figsize=(FW * P.MM, FH * P.MM))
ax_ = lambda x, y, w, h_: fig.add_axes([x / FW, y / FH, w / FW, h_ / FH])
norm = TwoSlopeNorm(0, -1, 1); cmap = plt.get_cmap("RdBu_r")
for c, lay in enumerate(("Isotropic", "Sym1")):
    d = cell[(lay, 4)]; X = d["X"]
    for r, k in enumerate((0, 1, 3)):                  # pi, pi/2 (first of the double pair), 0
        x0, y0 = 1.0 + 80.0 * c, FH - 38.0 * (r + 1)
        m = d["modes"][k]; U = m["U"].copy()
        if U.flat[np.argmax(np.abs(U))] < 0:
            U = -U
        kb = kb_of(m["c"]); t = table[lay][kb]
        fig.text(x0 / FW, (y0 + 35.5) / FH, rf"({'abcdef'[2 * r + c]}) {lay}, $k_B b_s = {kb}$", fontsize=8, va="top")
        a3 = ax_(x0, y0 + 1.0, 44.0, 29.0)
        FP._mode3d_ortho(a3, dict(d, U=U), amp=40.0, azim=-62.0, elev=24.0, cmap=cmap, norm=norm, xstride=2)
        b = ax_(x0 + 45.0, y0 + 6.0, 33.0, 12.0)
        i = int(np.argmax(np.abs(U).max(axis=(1, 2)))); Ui = U[i] / np.abs(U[i]).max()
        segs0 = [[X[a], X[b_]] for a, b_, part, _n in d["S"] if part != "link"]
        segs1 = [[X[a] + 40 * Ui[a], X[b_] + 40 * Ui[b_]] for a, b_, part, _n in d["S"] if part != "link"]
        b.add_collection(LineCollection(segs0, colors="#c9d0d8", linewidths=0.4))
        b.add_collection(LineCollection(segs1, colors=P.INK, linewidths=0.7))
        for y in np.arange(-d["W"], 1, 150.0):
            b.axvline(y, color="#c0392b", lw=0.3, ls=(0, (2, 2)), zorder=0)
        b.set_xlim(-d["W"] - 12, 12); b.set_ylim(-50, 80); b.set_aspect("equal"); b.axis("off")
        dbl = " (double)" if k == 1 else ""
        a_dbl = dbl
        mac_txt = ""
        # Abaqus modes of the four-pitch shell model (Aya: mode 1 on 11 Nov, modes 2-4 on 13 Nov 2026); their wavenumbers
        # are measured from the Abaqus skin deflection itself (check_npitch_modes.py -> results/npitch_abq_modes.json)
        import abq_modes as AQ
        AQ.FE = H / "src" / "abaqus_npitch" / "fe_modes"
        ka = {0: 1, 1: 2, 3: 4}[k]
        ab = AQ.read_abq(f"MODE{ka}_val_bs150_beta8_{lay}_PBC_n4.csv"); assert abs(ab["eig"] - t["N_abq"]) < 1e-6
        qa = AQ_CHECK[lay][ka - 1]; assert kb_of(qa["c_abaqus"]) == kb, (lay, ka, qa["c_abaqus"], kb)
        if k == 1:                                     # double pair: share of the Abaqus mode in the plane of the GBT pair
            mac_txt = "\n" + rf"subspace MAC {min(AQ_CHECK[lay][1]['share_in_gbt_pi2_plane'], AQ_CHECK[lay][2]['share_in_gbt_pi2_plane']):.4f}"
            table[lay][kb]["subspace_mac"] = [AQ_CHECK[lay][1]["share_in_gbt_pi2_plane"], AQ_CHECK[lay][2]["share_in_gbt_pi2_plane"]]
        else:
            mc, sg = AQ.mac(dict(d, U=U), ab)
            zc = ab["xyz"][np.argmin(np.abs(ab["xyz"][:, 2] - d["xs"][i])), 2]
            sel = np.abs(ab["xyz"][:, 2] - zc) < 1e-3
            Ua = sg * ab["U"][sel][:, :2]; Ua = Ua / np.abs(Ua).max()
            o_ = np.argsort(ab["xyz"][sel, 0])[::3]                     # every third node: the GBT line stays visible
            b.plot(ab["xyz"][sel, 0][o_] + 40 * Ua[o_, 0], ab["xyz"][sel, 1][o_] + 40 * Ua[o_, 1], "o", ms=0.9, mfc="#c0392b",
                   mec="none", zorder=4)
            mac_txt = "\n" + rf"(dots), MAC {mc:.4f}"
            table[lay][kb]["mac_skin"] = mc
        a_lab = "Abaqus"
        fig.text((x0 + 61.5) / FW, (y0 + 26.0) / FH, rf"GBT: {m['N']:.2f} N/mm{dbl}" + "\n" + rf"{a_lab}: {t['N_abq']:.2f} N/mm{a_dbl}{mac_txt}",
                 ha="center", va="top", fontsize=6.8, linespacing=1.35)
        fig.text((x0 + 61.5) / FW, (y0 + 4.5) / FH, rf"cross-section at $X/L = {d['xs'][i] / d['L']:.2f}$", ha="center", va="top",
                 fontsize=6.2, color=P.MUTED)
cax = fig.add_axes([0.30, 4.5 / FH, 0.40, 1.2 / FH])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label(r"normal displacement $w/|w|_{\max}$", fontsize=7.5); cb.ax.tick_params(labelsize=6.5, length=2)
cb.set_ticks([-1, -0.5, 0, 0.5, 1]); cb.outline.set_linewidth(0.4)
P.save(fig, "fig10_npitch")
json.dump(table, open(H / "results" / "npitch_table.json", "w"), indent=1)
print(json.dumps(table, indent=1))
