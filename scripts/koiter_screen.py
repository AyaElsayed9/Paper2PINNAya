"""19 Nov 2026: stability screen of the Paper 2 design space with coupled-mode Koiter (modes within 5 %).
For each test design: one-mode b, number of coupled modes, and whether the coupled path with a mode-shaped
imperfection (0.1 / 0.5 mm) reaches a maximum load below 1.5 N_cr (= imperfection-sensitive).
    python scripts/koiter_screen.py --per 5      -> results/koiter/screen.csv"""
import argparse, csv, sys, time
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck import opio, koiter as KO
from pinnbuck.families import family_context
from pinnbuck.designspace import Sample

TAGS = [f"{s}_{l}_SS_v2" for l in ("crossply", "quasiiso", "zerorich", "pm45") for s in ("symmetry", "periodic")]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--per", type=int, default=5); ap.add_argument("--gap", type=float, default=5.0)
    a = ap.parse_args(); torch.set_num_threads(2)
    out = ROOT / "results/koiter/screen.csv"; rows = []
    for tag in TAGS:
        files = sorted((ROOT / f"data/operators/{tag}/test").glob("*.pkl"))
        for f in files[:: max(1, len(files) // a.per)][:a.per]:
            op = opio.load(f); smp = Sample(**{k: v for k, v in op.info["sample"].items() if k != "sid"})
            t0 = time.perf_counter()
            with family_context(smp) as fc:
                S = KO.assemble(fc); lam, V = KO.buckle(S, 6)
                M = int(np.sum(lam <= lam[0] * (1 + a.gap / 100)))
                k1 = KO.koiter_multi(S, lam, V, 1); km = KO.koiter_multi(S, lam, V, M); S["release"]()
            b1 = 4 * k1["B"][0, 0, 0, 0] / (lam[0] * k1["G"][0]); a1 = 3 * k1["c"][0, 0, 0] / (lam[0] * k1["G"][0])
            row = dict(tag=tag, sid=op.sid, bs=float(smp.params["bs"]), beta=float(smp.params["beta"]),
                       L=float(fc.L), N1=float(lam[0]), gap12=100 * (lam[1] / lam[0] - 1), M=M, a_one=a1, b_one=b1)
            for xb in (0.1, 0.5):
                for sgn in (1, -1):
                    L_, X_, D_ = KO.multi_path(km, [sgn * xb] + [0.0] * (M - 1), lam_max_ratio=1.5, ds=0.02)
                    i = int(np.argmax(L_)); peak = L_[i] / lam[0]
                    limit = bool(i < len(L_) - 1 and L_[-1] < L_[i] - 1e-6)
                    row[f"peak_{xb}_{'p' if sgn > 0 else 'm'}"] = peak if limit else np.nan
            row["t_s"] = round(time.perf_counter() - t0, 1)
            rows.append(row); print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
            with open(out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
