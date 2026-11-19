"""19 Nov 2026: where does the critical mode stop being local?  GBT buckling of longer omega panels (beta 3.5-10)
for the four expert lay-ups, symmetry sides, SS.  Per design: the lowest 6 loads and the family energy
participation (G global, D distortional, L local, other) of each mode -> results/koiter/long_scan.csv"""
import csv, sys, time
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import Sample, LAYUP_CATALOG, FAMILIES
from pinnbuck.families import family_context
from pinnbuck import koiter as KO


def parts(S, a):
    sec = S["sec"]; nb = S["nb"]; Ka = S["K"] @ a; tot = float(a @ Ka)
    grp = np.array(sec.group); out = {}
    for g in ("G", "D", "L"):
        k = np.nonzero(grp == g)[0]
        idx = (k[:, None] * nb + np.arange(nb)[None, :]).ravel()
        out[g] = float(a[idx] @ Ka[idx]) / tot
    out["other"] = 1 - sum(out.values())
    return out


def main():
    torch.set_num_threads(2); out = ROOT / "results/koiter/long_scan.csv"; rows = []
    for lay in ("crossply", "quasiiso", "zerorich", "pm45"):
        sk, st = LAYUP_CATALOG[lay]
        for bs in (124.25, 205.0, 285.5):
            for beta in (3.5, 5.0, 7.0, 10.0):
                s = Sample("omega", dict(bs=bs, beta=beta, bc="SS", side="symmetry"),
                           {"skin": [float(x) for x in sk], "stringer": [float(x) for x in st]}, dict(FAMILIES["omega"]["fixed"]))
                t0 = time.perf_counter()
                with family_context(s) as fc:
                    S = KO.assemble(fc); lam, V = KO.buckle(S, 6)
                    P = [parts(S, S["Cm"] @ V[:, i]) for i in range(len(lam))]
                    S["release"]()
                row = dict(lay=lay, bs=bs, beta=beta, L=bs * beta)
                for i, (l, p) in enumerate(zip(lam, P)):
                    row[f"N{i+1}"] = float(l); dom = max(p, key=p.get)
                    row[f"type{i+1}"] = f"{dom}{p[dom]:.2f}"
                row["t_s"] = round(time.perf_counter() - t0, 1)
                rows.append(row); print(row, flush=True)
                with open(out, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
