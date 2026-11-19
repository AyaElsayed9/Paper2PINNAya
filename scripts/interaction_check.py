"""19 Nov 2026: full GBT model + coupled-mode Koiter at the local-global crossing found by interaction_scan.py.
Designs: h 15 mm, b_s 124.25, crossply and quasiiso, betas around the crossing (first pass 0.9-1.1 beta*).  Per design: lowest modes with
family type, coupled modes (types mixed within 10 %), coupled-Koiter path with mode-shaped imperfections
(0.1, 0.5 mm, both signs) -> maximum load / N_cr if there is a limit point.  -> results/koiter/interaction_check.csv"""
import csv, sys, time
from pathlib import Path
import numpy as np, pandas as pd, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck.designspace import Sample, LAYUP_CATALOG, FAMILIES
from pinnbuck.families import family_context
from pinnbuck import koiter as KO
from long_panel_scan import parts


def main():
    torch.set_num_threads(2); scan = pd.read_csv(ROOT / "results/koiter/interaction_scan.csv")
    out = ROOT / "results/koiter/interaction_check.csv"; rows = []
    # first pass (0.9-1.1 beta* of the scan) showed a lower mixed global mode for quasiiso -> betas set by hand
    plan = dict(quasiiso=[3.6, 4.0, 4.4], crossply=[5.3, 5.9, 6.5])
    for lay, betas in plan.items():
        r0 = scan[(scan.lay == lay) & (scan.h == 15.0) & (scan.bs == 124.25)].iloc[0]
        sk, st = LAYUP_CATALOG[lay]; fixed = dict(FAMILIES["omega"]["fixed"]); fixed["h"] = 15.0
        for beta in betas:
            f = beta / float(r0.beta_star)
            s = Sample("omega", dict(bs=124.25, beta=beta, bc="SS", side="symmetry"),
                       {"skin": [float(x) for x in sk], "stringer": [float(x) for x in st]}, fixed)
            t0 = time.perf_counter()
            with family_context(s) as fc:
                S = KO.assemble(fc); lam, V = KO.buckle(S, 8)
                types = []
                for i in range(len(lam)):
                    p = parts(S, S["Cm"] @ V[:, i]); d = max(p, key=p.get); types.append(f"{d}{p[d]:.2f}")
                M = int(np.sum(lam <= lam[0] * 1.10)); M = min(M, 4)
                km = KO.koiter_multi(S, lam, V, M); S["release"]()
            row = dict(lay=lay, h=15.0, bs=124.25, beta=beta, f_of_beta_star=f, L=124.25 * beta, M=M,
                       modes=" ".join(f"{l:.2f}:{t}" for l, t in zip(lam, types)))
            for xb in (0.1, 0.5):
                for sg in (1, -1):
                    L_, X_, D_ = KO.multi_path(km, [sg * xb] + [0.0] * (M - 1), lam_max_ratio=1.5, ds=0.02)
                    i = int(np.argmax(L_)); lim = bool(i < len(L_) - 1 and L_[-1] < L_[i] - 1e-6)
                    row[f"peak_{xb}_{'p' if sg > 0 else 'm'}"] = float(L_[i] / lam[0]) if lim else np.nan
            row["t_s"] = round(time.perf_counter() - t0, 1); rows.append(row)
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
            with open(out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
