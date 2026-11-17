import sys, time
from pathlib import Path
import numpy as np, scipy.sparse as sp
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import sample_design
from pinnbuck.families import family_context
tag = sys.argv[1]; side, lay = tag.split("_")[:2]
S = sample_design("omega", 40, seed=99, categorical={"side": side, "bc": "SS"}, layups=lay, space="v2")
for s in S[:4]:
    with family_context(s) as fc:
        m, sec, L = fc.m, fc.sec, fc.L; m29 = sec.m29
        hc = 4 * fc.hx
        t = time.perf_counter(); Kcs = m29._gbt_Kcs(sec); t1 = time.perf_counter() - t
        lb = m29._LongSpline(L, m29._sfs_breaks(L, hc, hc, L))
        t = time.perf_counter(); K = None
        for (d, e), A in Kcs.items():
            X = (lb.B[d].T @ sp.diags(lb.wx) @ lb.B[e]).tocsr(); T = sp.kron(sp.csr_matrix(A), X, format="csr"); K = T if K is None else K + T
        t2 = time.perf_counter() - t
        t = time.perf_counter(); fc.pre(coarsen=4); t3 = time.perf_counter() - t
        t = time.perf_counter(); Kcs2 = m29._gbt_Kcs(sec); t4 = time.perf_counter() - t
        nz_frac = {k: float((np.abs(A) > 0).mean()) for k, A in Kcs.items()}
        print(f"{s.sid} L {L:.0f} nb {lb.nb} nm {sec.nm}: Kcs {t1:.2f} s (again {t4:.2f}), kron assembly {t2:.2f} s, nnz K {K.nnz/1e6:.1f} M, whole pre x4 {t3:.2f} s; keys {sorted(Kcs)}; density {[round(v,2) for v in nz_frac.values()]}", flush=True)
