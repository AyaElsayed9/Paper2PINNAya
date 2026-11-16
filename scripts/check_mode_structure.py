"""Step 1 (16 Nov 2026): what do the exact GBT modes look like, and in which coordinates are they simple?

The network predicts ~12000 coefficients per mode on FIXED physical wavelengths (25 .. 1006 mm) + end splines.
Test of another description: harmonics of the panel length, q(z) = sum_j sin(j pi z / L) v_j + cos(j pi z / L) w_j
(v_j, w_j: cross-section vectors in the shared basis).  For every labelled design and its first modes: share of the
mode (L2 over z and cross-section amplitudes) carried by the best 1, 2, 3, 5 harmonics, and the leading harmonic j*.
For comparison: number of dictionary coefficient groups (wavelength or end layer) needed for the same share.
    python scripts/check_mode_structure.py symmetry_crossply_SS_v2 --n 40
-> results/mode_structure/<tag>.csv"""
import argparse, csv, sys
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import epinn, opio


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--modes", type=int, default=3); ap.add_argument("--inputs", default="geom")
    a = ap.parse_args(); torch.set_num_threads(2)
    epinn.INPUT_COLS = None if a.inputs == "all" else [0, 1]
    files = [l.strip() for l in open(ROOT / f"runs/experts/{a.tag}_files.txt") if l.strip()][:a.n]
    rows = []
    for f in files:
        op = opio.load(ROOT / f); d = epinn.prepare(op)
        ex = torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{Path(f).stem}.pt", weights_only=False)
        T = ex["T"][:a.modes]
        A = d.rs.unflat(d.fm.color(T)); Q = d.rs.field(A)[:, 0].numpy()          # (R, nz, nm) amplitudes along z
        z, wz, L = np.asarray(op.z), np.asarray(op.wz), op.L
        J = int(min(60, max(8, np.ceil(L / 20.0))))
        S = np.stack([np.sin(j * np.pi * z / L) for j in range(1, J + 1)], 1)     # (nz, J)
        Cc = np.stack([np.cos(j * np.pi * z / L) for j in range(0, J + 1)], 1)    # (nz, J + 1)
        Bh = np.concatenate([S, Cc], 1); sw = np.sqrt(wz)[:, None]
        # dictionary groups (whitened coordinates): energy share per wavelength / end layer
        key = lambda lab: lab[:-1] if lab.startswith("wave") else lab[:4]
        gk = []
        for c, labs in enumerate(d.rs.labels):
            for kk in range(d.rs.k[c]):
                gk += [key(l) for l in labs]
        gk = np.array(gk)
        for r in range(T.shape[0]):
            q = Q[r]; tot = float(((q * sw) ** 2).sum())
            coef, *_ = np.linalg.lstsq(Bh * sw, q * sw, rcond=None)                # (2J + 1, nm)
            fit = Bh @ coef; res_all = float((((q - fit) * sw) ** 2).sum()) / tot
            # share per harmonic j (sin j + cos j together), from separate projections on each pair (not orthogonal
            # under the quadrature only to round-off for integer harmonics on [0, L])
            sh = np.zeros(J + 1)
            for j in range(0, J + 1):
                cols = ([j - 1] if j >= 1 else []) + [J + j]
                Bj = Bh[:, cols] * sw; cj, *_ = np.linalg.lstsq(Bj, q * sw, rcond=None)
                sh[j] = float(((Bj @ cj) ** 2).sum()) / tot
            order = np.argsort(-sh)
            def best(m):                                                           # joint fit with the m best harmonics
                cols = []
                for j in order[:m]:
                    cols += ([j - 1] if j >= 1 else []) + [J + j]
                Bm = Bh[:, cols] * sw; cm, *_ = np.linalg.lstsq(Bm, q * sw, rcond=None)
                return float(((Bm @ cm) ** 2).sum()) / tot
            y = T[r].numpy() ** 2; y = y / y.sum()
            gs = np.sort([y[gk == g].sum() for g in np.unique(gk)])[::-1]; cg = np.cumsum(gs)
            rows.append(dict(sid=op.sid, bs=round(op.params["bs"], 1), beta=round(op.params["beta"], 3), L=round(L, 1), mode=r + 1,
                             N=round(float(ex["N_ritz"][r]), 4), jstar=int(order[0]), half_wave=round(L / max(int(order[0]), 1), 1),
                             h1=round(best(1), 4), h2=round(best(2), 4), h3=round(best(3), 4), h5=round(best(5), 4),
                             h_all=round(1 - res_all, 4), dict_groups_90=int(np.searchsorted(cg, 0.90) + 1),
                             dict_groups_99=int(np.searchsorted(cg, 0.99) + 1), dict_top1=round(float(gs[0]), 3)))
        del op, d
        print(rows[-a.modes], flush=True)
    out = ROOT / f"results/mode_structure/{a.tag}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    R = np.array([[r[k] for k in ("h1", "h2", "h3", "h5", "h_all", "dict_groups_90", "dict_groups_99", "dict_top1")] for r in rows if r["mode"] == 1])
    print("mode 1, mean over designs: best 1/2/3/5 harmonics %.3f %.3f %.3f %.3f, all harmonics %.4f | dictionary groups for 90 / 99 %%: %.1f / %.1f, top group %.2f"
          % tuple(R.mean(0)))
    print("mode 1, worst design: best 1/2/3/5 harmonics %.3f %.3f %.3f %.3f" % tuple(R[:, :4].min(0)))


if __name__ == "__main__":
    main()
