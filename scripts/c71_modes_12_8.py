"""C71 check: what really separates hierarchical modes 12 and 8 (Unsym, b_s 150, beta 1.5, SS).

Prints: in-plane shape correlation, crown share of the displacement, and the split of the
transverse (B-matrix) energy into membrane stretching, plate bending and coupling per wall type,
checked against the solver's own B quadratic form; then how the two modes combine in the
critical mode. Run from manuscript/:  python3 ../scripts/c71_modes_12_8.py
"""
import sys
import numpy as np

sys.argv = ["x"]
sys.path.insert(0, ".")
import figures_paper as F  # noqa: E402

P = F.P
h = F._hier_mode("Unsym", 150.0, 1.5, "SS")
m, sec, L = P.case("Unsym", 150.0, 1.5, "SS")
hm = m.hierarchical_modes(sec)
lam, T = hm["C"]
A = hm["_A"]; ix = np.arange(sec.nip)
B = A[(0, 0)][np.ix_(ix, ix)]
X = sec.X
nz = lambda t: t / np.hypot(*(sec.U @ t).reshape(-1, 2).T).max()      # max in-plane displacement = 1


def split(t):
    U = (sec.U @ t).reshape(-1, 2); th = sec.Th @ t
    mem, ben, cpl, crown = {}, {}, {}, 0.0
    for s_ in sec.S:
        a, b = s_["a"], s_["b"]; p = sec.wpart[s_["w"]]; abd = sec.walls[s_["w"]][2]
        d = X[b] - X[a]; l = np.hypot(*d); e = d / l; n = np.array([-e[1], e[0]])
        es = (U[b] - U[a]) @ e / l
        wa, wb = U[a] @ n, U[b] @ n
        k0 = 6 * (wb - wa) / l ** 2 - (4 * th[a] + 2 * th[b]) / l
        k1 = -6 * (wb - wa) / l ** 2 + (2 * th[a] + 4 * th[b]) / l
        mem[p] = mem.get(p, 0) + abd[1, 1] * es ** 2 * l
        ben[p] = ben.get(p, 0) + abd[4, 4] * l * (k0 ** 2 + k0 * k1 + k1 ** 2) / 3
        cpl[p] = cpl.get(p, 0) - 2 * abd[1, 4] * es * l * (k0 + k1) / 2
        if p == "crown":
            crown = max(crown, np.hypot(*U[a]), np.hypot(*U[b]))
    return mem, ben, cpl, crown


t12, t8 = nz(T[:, 11]), nz(T[:, 7])
u12, u8 = sec.U @ t12, sec.U @ t8
print("in-plane shape correlation 12 vs 8: %.4f" % (u12 @ u8 / np.linalg.norm(u12) / np.linalg.norm(u8)))
for nm, t in (("12", t12), ("8", t8)):
    mem, ben, cpl, cr = split(t)
    M, Bn, C = sum(mem.values()), sum(ben.values()), sum(cpl.values())
    print(f"mode {nm}: lam={lam[int(nm) - 1]:.3f}  solver B={t @ B @ t:.3f} = membrane {M:.3f} (webs {mem['web']:.3f})"
          f" + bending {Bn:.3f} + coupling {C:.3f};  crown max displacement {cr:.2f}")
Aa = h["amp"]; i = np.argmax(np.abs(Aa[:, 11])); c12, c8 = Aa[i, 11], Aa[i, 7]
for nm, t in (("12 part", c12 * t12), ("8 part", c8 * t8), ("12 + 8", c12 * t12 + c8 * t8)):
    mem, *_ = split(t)
    print(f"critical mode, {nm:7s}: web membrane energy {mem['web']:.2e}, max displacement "
          f"{np.hypot(*(sec.U @ t).reshape(-1, 2).T).max():.4f}")
