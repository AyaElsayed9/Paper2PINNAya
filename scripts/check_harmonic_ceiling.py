"""Step 1 (16 Nov 2026): ceiling of a length-harmonic description of the modes.

For labelled designs: the design's OWN exact modes (8), cut down to their m strongest harmonics of the panel length
(sin + cos of the same j count as one), projected back on the dictionary, then Ritz value (s0) and the check.
If one or a few harmonics already give a Ritz value within ~1 %, a network that predicts cross-section vectors per
harmonic has a high ceiling; if not, the small non-harmonic part (end zones) matters and it has not.
    python scripts/check_harmonic_ceiling.py symmetry_crossply_SS_v2 --n 20
-> results/check_baselines/<tag>_harmonic_ceiling.csv"""
import argparse, csv, sys
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio
from pinnbuck.spectral import lobpcg
from check_harmonic_start import decompose, compose


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--keep", default="1,2,3,5,all"); ap.add_argument("--inputs", default="geom")
    a = ap.parse_args(); torch.set_num_threads(2)
    epinn.INPUT_COLS = None if a.inputs == "all" else [0, 1]
    files = [l.strip() for l in open(ROOT / f"runs/experts/{a.tag}_files.txt") if l.strip()][:a.n]
    rows = []; steps = [0, 1, 2, 3, 5]
    for f in files:
        d = epinn.prepare(opio.load(ROOT / f))
        ex = torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{Path(f).stem}.pt", weights_only=False)
        T = ex["T"]; Nref = float(ex["N_ritz"][0])                 # reference = converged load in the same Ritz space
        C, J, L = decompose(d, T)
        e2 = (C ** 2).sum(2)                                       # (R, 2J + 1)
        share = np.concatenate([e2[:, J:J + 1], e2[:, :J] + e2[:, J + 1:]], 1)      # j = 0..J
        for keep in a.keep.split(","):
            Ck = C.copy()
            if keep != "all":
                for r in range(C.shape[0]):
                    top = np.argsort(-share[r])[:int(keep)]; mask = np.zeros(2 * J + 1, bool)
                    for j in top:
                        mask[J + j] = True
                        if j >= 1:
                            mask[j - 1] = True
                    Ck[r][~mask] = 0.0
            Y = compose(d, Ck, J, np.arange(1, J + 1))
            mu0, _ = epinn.ritz(d, Y)
            _, _, h = lobpcg(d.rs, R=8, iters=5, X0=d.rs.unflat(d.fm.color(Y.double())), precond=d.em.solve, tol=0, seed=1)
            e = [100 * (v / Nref - 1) for v in [1.0 / float(mu0[0])] + list(h)]
            rows.append(dict(sid=d.sid, bs=round(d.params["bs"], 1), beta=round(d.params["beta"], 3), keep=keep, **{f"s{s}": round(e[s], 4) for s in steps}))
            print(rows[-1], flush=True)
        del d
    with open(ROOT / f"results/check_baselines/{a.tag}_harmonic_ceiling.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n{a.tag}: error vs the converged Ritz-space load, mean / max over {len(files)} labelled designs")
    for keep in a.keep.split(","):
        rr = [r for r in rows if r["keep"] == keep]
        print(f"  harmonics kept per mode: {keep:>3s}  " + "  ".join(f"s{s} {np.mean([abs(r[f's{s}']) for r in rr]):.2f}/{np.max([abs(r[f's{s}']) for r in rr]):.2f}" for s in steps))


if __name__ == "__main__":
    main()
