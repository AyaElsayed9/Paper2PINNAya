"""Check the explicit conventional-shear interaction matrices of the paper against the code."""
import sys; sys.path.insert(0, '.')
import numpy as np
import verify_math as V
for lay in ("Unsym", "Sym2"):
    m, sec = V.get(lay)
    K = m._gbt_Kcs(sec); E, C, w = sec.E, sec.Cp, sec.wp
    I = slice(0, sec.nip); J = slice(sec.o_s, sec.o_s + sec.ns)
    vs = E[:, 1, I, 0]; wss = -E[:, 4, I, 0]; gam = E[:, 2, I, 1]; ws = -0.5 * E[:, 5, I, 1]
    u = E[:, 0, I, 2]; wv = -E[:, 3, I, 2]
    ujs = E[:, 2, J, 0]; uj = E[:, 0, J, 1]
    A = C[:, :3, :3]; B = C[:, :3, 3:6]
    f = lambda c, a, b: np.einsum("p,p,pk,pl->kl", w, c, a, b)
    ex = {(0, 0): f(A[:, 1, 2], vs, ujs) - f(B[:, 1, 2], wss, ujs),          # F~
          (0, 1): f(A[:, 1, 0], vs, uj) - f(B[:, 1, 0], wss, uj),            # D~2
          (1, 0): f(A[:, 2, 2], gam, ujs) - 2 * f(B[:, 2, 2], ws, ujs),      # D~1
          (1, 1): f(A[:, 2, 0], gam, uj) - 2 * f(B[:, 2, 0], ws, uj),        # H~1
          (2, 0): f(A[:, 0, 2], u, ujs) - f(B[:, 0, 2], wv, ujs),            # H~3
          (2, 1): f(A[:, 0, 0], u, uj) - f(B[:, 0, 0], wv, uj)}              # C~
    for de, X in ex.items():
        R = K[de][I, J]
        print(lay, de, f"{np.abs(R - X).max() / max(np.abs(R).max(), 1e-30):.1e}")
    for de in [(0, 2), (1, 2), (2, 2)]:
        if de in K:
            print(lay, de, "conv-shear block should be zero:", np.abs(K[de][I, J]).max())
