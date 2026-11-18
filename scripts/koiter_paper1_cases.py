"""18 Nov 2026: Koiter initial post-buckling (src/pinnbuck/koiter.py) for the four Riks cases of Paper 1
(Compos Struct 399 (2026) 120959, Table 5): cross-ply, symmetry sides, SS, imperfection = mode 1 scaled to 0.01 mm.
Writes results/koiter/paper1_cases.json and the GBT paths at three edge nodes next to the digitized Abaqus curves
(docs/paper1_riks/case<k>_abaqus_digitized.csv) -> results/koiter/paper1_cases_paths.csv / .png"""
import json, sys
from pathlib import Path
import numpy as np, torch
from scipy.interpolate import BSpline
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.families import family_context
from pinnbuck import koiter as KO

CASES = {1: (205.0, 0.86), 2: (205.0, 1.67), 3: (150.0, 0.50), 4: (205.0, 1.77)}
XI_BAR = 0.01


def node_shortening(S, a, L):
    sec, lb = S["sec"], S["lb"]
    C = a.reshape(S["nm"], S["nb"]); spl = BSpline(lb.t, np.eye(lb.nb), 3); xs = np.array([0.0, L])
    A0 = spl(xs) @ C.T; A1 = spl.derivative()(xs) @ C.T
    u = sec.Uw @ A1[:, :sec.nip].T + A0[:, sec.o_ax][None, :] + sec.Us @ A0[:, sec.o_s:sec.o_s + sec.ns].T
    return u[:, 0] - u[:, 1]


def main():
    torch.set_num_threads(2); out = ROOT / "results/koiter"; out.mkdir(parents=True, exist_ok=True)
    res = {}; curves = {}
    for c, (bs, beta) in CASES.items():
        s = validated_omega_sample("Cross-ply", bs, beta, "SS", "symmetry")
        with family_context(s) as fc:
            S = KO.assemble(fc); lam, V = KO.buckle(S, 4); k = KO.koiter(S, lam, V)
            sh0 = node_shortening(S, S["a0"], fc.L); sh2 = node_shortening(S, k["u2"], fc.L); S["release"]()
            X = S["sec"].X
        rn = sh0 / (sh0 + sh2 / (k["b"] * k["lam_c"])) if k["b"] > 0 else np.full(len(sh0), np.nan)
        i_mid = int(np.argmin(np.hypot(X[:, 0] - 0.0, X[:, 1])))            # skin node on the symmetry edge x = 0
        i_crown = int(np.argmax(X[:, 1]))
        r = dict(bs=bs, beta=beta, L=bs * beta, lam=[float(v) for v in lam], gap12_pct=100 * (lam[1] / lam[0] - 1),
                 a=k["a"], b=k["b"], e0=k["e0"], e2=k["e2"], ratio_mean=k["stiff_ratio"],
                 ratio_node_min=float(np.nanmin(rn)), ratio_node_max=float(np.nanmax(rn)),
                 ratio_edge_x0=float(rn[i_mid]), ratio_crown=float(rn[i_crown]),
                 pre_slope=1.0 / k["e0"])
        res[c] = r
        print(c, {q: (round(v, 4) if isinstance(v, float) else v) for q, v in r.items()}, flush=True)
        dig = np.loadtxt(ROOT / f"docs/paper1_riks/case{c}_abaqus_digitized.csv", delimiter=",")
        lc = k["lam_c"]
        def slope(lo, hi):                               # Abaqus slope between two load levels (fractions of lam_c)
            m = (dig[:, 1] >= lo * lc) & (dig[:, 1] <= hi * lc)
            return float(np.polyfit(dig[m, 0], dig[m, 1], 1)[0]) if m.sum() > 5 else float("nan")
        r.update(abq_pre_slope=slope(0.2, 0.7), abq_post_slope_1=slope(1.08, 1.3), abq_post_slope_2=slope(1.4, 1.8))
        r["abq_ratio_1"] = r["abq_post_slope_1"] / r["abq_pre_slope"]; r["abq_ratio_2"] = r["abq_post_slope_2"] / r["abq_pre_slope"]
        curves[c] = (dig, {}, k)
        for name, (e0, e2) in dict(load_mean=(k["e0"], k["e2"]), skin_edge_x0=(sh0[i_mid], sh2[i_mid]),
                                   crown=(sh0[i_crown], sh2[i_crown])).items():
            if k["b"] > 0:
                l, d = KO.imperfect_path(k, XI_BAR, lam_max_ratio=1.9, e0=e0, e2=e2); curves[c][1][name] = (d, l)
        print("   Abaqus slopes pre %.1f post %.1f / %.1f -> ratio %.3f / %.3f" % (r["abq_pre_slope"], r["abq_post_slope_1"],
              r["abq_post_slope_2"], r["abq_ratio_1"], r["abq_ratio_2"]), flush=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))
    sty = dict(load_mean=("-", "GBT Koiter: load-weighted mean shortening"), skin_edge_x0=("--", "GBT Koiter: skin node on the side edge"),
               crown=(":", "GBT Koiter: stringer crown node"))
    for ax, (c, (dig, cv, k)) in zip(axs.ravel(), curves.items()):
        ax.plot(dig[:, 0], dig[:, 1], color="#d62728", lw=2.2, label="Abaqus Riks (Paper 1, digitized)")
        for name, (d, l) in cv.items():
            ax.plot(d, l, color="#1f4e79", ls=sty[name][0], lw=1.6, label=sty[name][1])
        ax.axhline(k["lam_c"], color="gray", lw=0.8, ls="-.", label="GBT critical load")
        xm = dig[:, 0].max() * 1.05; ax.set_xlim(0, xm); ax.set_ylim(0, max(dig[:, 1].max(), k["lam_c"]) * 1.12)
        bs, beta = CASES[c]
        ax.set_title(f"Case {c}: b_s {bs:g} mm, beta {beta:g}, gap modes 1-2 {res[c]['gap12_pct']:.1f} %, b = {k['b']:.3f}", fontsize=9)
        ax.set_xlabel("end-shortening [mm]"); ax.set_ylabel("N [N/mm]"); ax.grid(alpha=0.3)
    axs[0, 1].legend(fontsize=7, loc="lower right"); fig.tight_layout(); fig.savefig(out / "paper1_cases_paths.png", dpi=140)
    with open(out / "paper1_cases_paths.csv", "w") as fh:
        fh.write("case,curve,D_mm,N_Nmm\n")
        for c, (dig, cv, k) in curves.items():
            for name, (d, l) in cv.items():
                for x, y in zip(d, l):
                    fh.write(f"{c},{name},{x:.6f},{y:.5f}\n")
    json.dump(res, open(out / "paper1_cases.json", "w"), indent=1)


if __name__ == "__main__":
    main()
