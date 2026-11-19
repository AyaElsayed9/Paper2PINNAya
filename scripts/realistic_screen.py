"""19 Nov 2026: stability screen of a realistic design space, anchored on the NASA / Politecnico hat-stiffened panel
(Bisagni, Vescovini & Davila: L 720 mm, pitch 152 mm, h 30 mm, web 25 deg, 15 mm flanges, skin [45/90/-45/0]s 1 mm,
stringer [-45/0/45/0/45/0/-45] 0.875 mm) = the geometry of Paper 1 and of the GBT paper.
Grid: the four expert lay-ups x b_s {124.25, 152, 205} x h {20, 25, 30, 35} x L {500, 750} mm (frame pitch),
symmetry sides, SS.  Per design: lowest 8 GBT loads with family type, global load (G + axial + shear-warping
modes only), coupled-mode Koiter (modes within 10 % of N1, at most 5), paths with mode-shaped imperfections
+-0.1 / +-0.5 / +-1.0 mm to 1.5 N1: limit point or not; one-mode post-buckling stiffness ratio (load-weighted).
    python scripts/realistic_screen.py   -> results/koiter/realistic_screen.csv"""
import csv, sys, time
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck.designspace import Sample, LAYUP_CATALOG, FAMILIES
from pinnbuck.families import family_context
from pinnbuck import koiter as KO
from long_panel_scan import parts
from interaction_scan import global_load


def main():
    torch.set_num_threads(2); out = ROOT / "results/koiter/realistic_screen.csv"
    done = set(); rows = []
    if out.exists():
        import pandas as pd
        old = pd.read_csv(out); rows = old.to_dict("records"); done = {(r["lay"], r["bs"], r["h"], r["L"]) for r in rows}
    for lay in ("crossply", "quasiiso", "zerorich", "pm45"):
        sk, st = LAYUP_CATALOG[lay]
        for h in (20.0, 25.0, 30.0, 35.0):
            for bs in (124.25, 152.0, 205.0):
                for L in (500.0, 750.0):
                    if (lay, bs, h, L) in done:
                        continue
                    fixed = dict(FAMILIES["omega"]["fixed"]); fixed["h"] = h
                    s = Sample("omega", dict(bs=bs, beta=L / bs, bc="SS", side="symmetry"),
                               {"skin": [float(x) for x in sk], "stringer": [float(x) for x in st]}, fixed)
                    t0 = time.perf_counter()
                    with family_context(s) as fc:
                        S = KO.assemble(fc); lam, V = KO.buckle(S, 8)
                        types = []
                        for i in range(len(lam)):
                            p = parts(S, S["Cm"] @ V[:, i]); d = max(p, key=p.get); types.append(f"{d}{p[d]:.2f}")
                        Ng = global_load(S)
                        M = int(min(5, np.sum(lam <= lam[0] * 1.10)))
                        k1 = KO.koiter(S, lam, V); km = KO.koiter_multi(S, lam, V, M); S["release"]()
                    row = dict(lay=lay, bs=bs, h=h, L=L, beta=L / bs, N1=float(lam[0]), type1=types[0],
                               gap12=100 * (lam[1] / lam[0] - 1), N_global=float(Ng), r_glob=float(Ng / lam[0]), M=M,
                               b_one=k1["b"], a_one=k1["a"], stiff_ratio=k1["stiff_ratio"],
                               modes=" ".join(f"{l:.2f}:{t}" for l, t in zip(lam, types)))
                    for xb in (0.1, 0.5, 1.0):
                        pk = []
                        for sg in (1, -1):
                            L_, X_, D_ = KO.multi_path(km, [sg * xb] + [0.0] * (M - 1), lam_max_ratio=1.5, ds=0.02)
                            i = int(np.argmax(L_)); lim = bool(i < len(L_) - 1 and L_[-1] < L_[i] - 1e-6)
                            pk.append(L_[i] / lam[0] if lim else np.nan)
                        row[f"peak_{xb}"] = float(np.nanmin(pk)) if np.isfinite(pk).any() else np.nan
                    row["t_s"] = round(time.perf_counter() - t0, 1); rows.append(row)
                    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items() if k != "modes"}, flush=True)
                    with open(out, "w", newline="") as fh:
                        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
