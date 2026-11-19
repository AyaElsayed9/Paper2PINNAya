"""19 Nov 2026: multi-mode (interaction) Koiter vs one-mode Koiter vs Paper 1's Riks, shortening at the Abaqus
monitor (skin node midway between the stringers).  Modes within GAP % of the critical load are coupled.
-> results/koiter/paper1_multi.json, paper1_multi_paths.png/.csv"""
import json, sys
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.families import family_context
from pinnbuck import koiter as KO
from koiter_paper1_cases import CASES, XI_BAR, node_shortening

GAP = 5.0


def midbay(X):
    top = X[X[:, 1] > 0.9 * X[:, 1].max(), 0]; cut = np.median(top)
    xm = 0.5 * (top[top < cut].mean() + top[top >= cut].mean())
    sk = np.nonzero(X[:, 1] < 1e-9)[0]
    return int(sk[np.argmin(np.abs(X[sk, 0] - xm))])


def main():
    torch.set_num_threads(2); out = ROOT / "results/koiter"; res = {}; curves = {}
    for c, (bs, beta) in CASES.items():
        s = validated_omega_sample("Cross-ply", bs, beta, "SS", "symmetry")
        with family_context(s) as fc:
            S = KO.assemble(fc); lam, V = KO.buckle(S, 6)
            M = int(np.sum(lam <= lam[0] * (1 + GAP / 100)))
            k1 = KO.koiter_multi(S, lam, V, 1); km = KO.koiter_multi(S, lam, V, M)
            i = midbay(S["sec"].X)
            e0n = node_shortening(S, S["a0"], fc.L)[i]
            def En(k):
                return np.array([[node_shortening(S, k["W"][a, b], fc.L)[i] for b in range(k["M"])] for a in range(k["M"])])
            E1, Em = En(k1), En(km)
            S["release"]()
        L1, X1, D1 = KO.multi_path(k1, [XI_BAR], lam_max_ratio=1.9, e0=e0n, E=E1)
        Lm, Xm, Dm = KO.multi_path(km, [XI_BAR] + [0.0] * (M - 1), lam_max_ratio=1.9, e0=e0n, E=Em)
        dig = np.loadtxt(ROOT / f"docs/paper1_riks/case{c}_abaqus_digitized.csv", delimiter=",")
        cmp = {}
        for name, (D, L) in dict(one=(D1, L1), multi=(Dm, Lm)).items():
            o = np.argsort(D); rows = []
            for fr in (0.5, 0.75, 1.0):
                d = fr * dig[:, 0].max(); na = np.interp(d, dig[:, 0], dig[:, 1]); ng = np.interp(d, D[o], L[o])
                rows.append(100 * (ng / na - 1))
            cmp[name] = rows
        res[c] = dict(bs=bs, beta=beta, lam=lam[:M].tolist(), M=M, b_one=4 * k1["B"][0, 0, 0, 0] / (lam[0] * k1["G"][0]),
                      peak_one=float(L1.max()), peak_multi=float(Lm.max()), abaqus_peak=float(dig[:, 1].max()),
                      dN_pct_one=cmp["one"], dN_pct_multi=cmp["multi"], xi_end_multi=Xm[-1].tolist())
        curves[c] = (dig, (D1, L1), (Dm, Lm), lam[0], M)
        print(c, {q: (np.round(v, 3).tolist() if isinstance(v, (list, np.ndarray)) else round(v, 4) if isinstance(v, float) else v)
                  for q, v in res[c].items()}, flush=True)
    json.dump(res, open(out / "paper1_multi.json", "w"), indent=1)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.5))
    with open(out / "paper1_multi_paths.csv", "w") as fh:
        fh.write("case,model,D_mm,N_Nmm\n")
        for ax, (c, (dig, one, mul, lc, M)) in zip(axs.ravel(), curves.items()):
            ax.plot(dig[:, 0], dig[:, 1], color="#d62728", lw=2.2, label="Abaqus Riks (Paper 1, digitized)")
            ax.plot(one[0], one[1], color="#1f4e79", ls="--", lw=1.5, label="GBT Koiter, one mode")
            ax.plot(mul[0], mul[1], color="#2ca02c", lw=1.8, label="GBT Koiter, coupled modes within 5 %")
            ax.axhline(lc, color="gray", lw=0.8, ls="-.", label="GBT critical load")
            ax.set_xlim(0, dig[:, 0].max() * 1.05); ax.set_ylim(0, max(dig[:, 1].max(), lc) * 1.15); ax.grid(alpha=0.3)
            bs, beta = CASES[c]
            ax.set_title(f"Case {c}: b_s {bs:g} mm, beta {beta:g} ({M} coupled modes)", fontsize=9)
            ax.set_xlabel("end-shortening at mid-bay [mm]"); ax.set_ylabel("N [N/mm]")
            for name, (D, L) in dict(one=one, multi=mul).items():
                for x, y in zip(D, L):
                    fh.write(f"{c},{name},{x:.6f},{y:.5f}\n")
    axs[0, 0].legend(fontsize=7.5, loc="lower right"); fig.tight_layout(); fig.savefig(out / "paper1_multi_paths.png", dpi=140)


if __name__ == "__main__":
    main()
