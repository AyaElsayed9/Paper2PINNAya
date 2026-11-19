"""19 Nov 2026: where can local and global buckling interact?  For each lay-up, stringer pitch b_s and stringer
height h: the local plateau load (lowest GBT mode at beta = 3, full model) and the global (column) load from the GBT
model restricted to the global (G) deformation modes plus the axial / shear warping modes, at L = 3 b_s and 6 b_s,
fitted by the Euler-with-shear form 1/N = L^2/a + 1/N_s; the crossing length L* solves N_glob(L*) = N_loc, and
beta* = L*/b_s.  A check run at L* follows for a few designs.
-> results/koiter/interaction_scan.csv"""
import csv, sys, time
from pathlib import Path
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import Sample, LAYUP_CATALOG, FAMILIES
from pinnbuck.families import family_context
from pinnbuck import koiter as KO


def global_load(S):
    sec, nb, Cm = S["sec"], S["nb"], S["Cm"]
    grp = np.array(sec.group); keep_modes = np.nonzero(np.isin(grp, ["G", "A", "S"]))[0]
    full = np.zeros(S["ndof"], bool)
    for m in keep_modes:
        full[m * nb:(m + 1) * nb] = True
    cols = np.nonzero(np.asarray(abs(Cm[~full, :]).sum(axis=0)).ravel() == 0)[0]   # reduced dofs living on kept modes only
    K = S["Kr"][cols][:, cols].tocsc(); G = S["Gr"][cols][:, cols].tocsc()
    lu = spla.splu(K); Kinv = spla.LinearOperator(K.shape, matvec=lu.solve, dtype=float)
    mu = spla.eigsh(G, k=2, M=K, Minv=Kinv, which="LA", tol=1e-9, return_eigenvectors=False)
    return 1.0 / mu.max()


def main():
    torch.set_num_threads(2); out = ROOT / "results/koiter/interaction_scan.csv"; rows = []
    for lay in ("crossply", "quasiiso", "zerorich", "pm45"):
        sk, st = LAYUP_CATALOG[lay]
        for h in (30.0, 20.0, 15.0):
            for bs in (124.25, 205.0, 285.5):
                fixed = dict(FAMILIES["omega"]["fixed"]); fixed["h"] = h
                s = Sample("omega", dict(bs=bs, beta=3.0, bc="SS", side="symmetry"),
                           {"skin": [float(x) for x in sk], "stringer": [float(x) for x in st]}, fixed)
                t0 = time.perf_counter()
                with family_context(s) as fc:
                    S = KO.assemble(fc); lam, V = KO.buckle(S, 2); Ng3 = global_load(S); S["release"]()
                s6 = Sample("omega", dict(bs=bs, beta=6.0, bc="SS", side="symmetry"),
                            {"skin": [float(x) for x in sk], "stringer": [float(x) for x in st]}, fixed)
                with family_context(s6) as fc:
                    S = KO.assemble(fc); Ng6 = global_load(S); S["release"]()
                L3, L6 = 3 * bs, 6 * bs
                # 1/N = L^2/a + 1/Ns  from two points
                inv_a = (1 / Ng6 - 1 / Ng3) / (L6 ** 2 - L3 ** 2); inv_Ns = 1 / Ng3 - inv_a * L3 ** 2
                Nl = float(lam[0]); rhs = 1 / Nl - inv_Ns
                Ls = float(np.sqrt(rhs / inv_a)) if rhs > 0 and inv_a > 0 else float("nan")
                row = dict(lay=lay, h=h, bs=bs, N_local=Nl, N_global_3bs=float(Ng3), N_global_6bs=float(Ng6),
                           N_shear=float(1 / inv_Ns) if inv_Ns > 0 else float("inf"), L_star=Ls, beta_star=Ls / bs,
                           t_s=round(time.perf_counter() - t0, 1))
                rows.append(row); print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
                with open(out, "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
