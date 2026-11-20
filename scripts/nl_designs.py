"""20 Nov 2026: geometrically non-linear GBT paths (src/pinnbuck/nlgbt.py, basis enrichment to an energy-norm error
<= 2 %) for the designs flagged by the realistic screen.  Imperfection: mode 1 shape, 0.5 mm.
    python scripts/nl_designs.py <name>    -> results/koiter/nl_<name>.npz / .log"""
import sys, time
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck.designspace import Sample, LAYUP_CATALOG, FAMILIES
from pinnbuck.families import family_context
from pinnbuck import koiter as KO, nlgbt as NL
from koiter_paper1_cases import node_shortening
from koiter_multi_paper1 import midbay

DESIGNS = {"qi_152_20_750": ("quasiiso", 152.0, 20.0, 750.0),      # mixed global mode 1.4 % below the local cluster
           "pm45_152_20_750": ("pm45", 152.0, 20.0, 750.0),         # mixed global mode far below local
           "qi_124_20_500": ("quasiiso", 124.25, 20.0, 500.0)}      # the one screen limit point (0.1 mm only)


def main():
    name = sys.argv[1]; lay, bs, h, L = DESIGNS[name]; xb = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    torch.set_num_threads(2); sk, st = LAYUP_CATALOG[lay]; fixed = dict(FAMILIES["omega"]["fixed"]); fixed["h"] = h
    s = Sample("omega", dict(bs=bs, beta=L / bs, bc="SS", side="symmetry"),
               {"skin": [float(x) for x in sk], "stringer": [float(x) for x in st]}, fixed)
    t0 = time.time()
    with family_context(s) as fc:
        S = KO.assemble(fc); lam, V = KO.buckle(S, 8)
        Phi, km = NL.basis(S, lam, V, M_modes=8, M_pairs=5)
        F = KO.Fields(S); u1 = S["Cm"] @ V[:, 0]; u1 = u1 / F.maxdisp(torch.as_tensor(u1))
        out, R = NL.solve_path(S, Phi, xb * u1, 1.6 * lam[0], tol_err=0.02, rounds=4, ds_max=2.0)
        Phi = out["basis"]; i = midbay(S["sec"].X)
        dmid = np.array([node_shortening(S, Phi[:, k], fc.L)[i] for k in range(Phi.shape[1])])
        dmean = (Phi.T @ S["f"]) / S["Ptot"]
        wmax = np.array([F.maxdisp(torch.as_tensor(Phi @ y)) for y in out["y"][::max(1, len(out["y"]) // 60)]])
        S["release"]()
    Lp = out["lam"]; k = int(np.argmax(Lp)); limit = bool(k < len(Lp) - 1 and Lp[-1] < Lp[k] - 1e-6)
    np.savez(ROOT / f"results/koiter/nl_{name}.npz", lam=Lp, D_mid=out["y"] @ dmid, D_mean=out["y"] @ dmean,
             wmax=wmax, N_lin=lam, err=np.array(out["res"]), basis=Phi.shape[1], rounds=out["round"], xibar=xb)
    print(f"{name}: N1 {lam[0]:.3f}, path to {Lp.max():.3f} ({Lp.max()/lam[0]:.3f} N1), limit point {limit}, "
          f"basis {Phi.shape[1]}, max error {max(e for _, e in out['res']):.4f}, {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
