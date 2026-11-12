import sys, json; sys.path.insert(0, '/home/claude/paper2pinnaya/manuscript')
import numpy as np, fp_common as P
from run_periodic import energy_participation
out = []
for beta in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]:
    m, sec, L = P.case('Unsym', 150.0, beta, 'SS')
    r = m.omega_gbt_solve(sec, L, 'SS', nev=1, return_mode=True)
    lb, C = r['lb'], np.asarray(r['coef'])
    rho, _ = energy_participation(m, sec, lb, C)
    xs = np.linspace(0, L, 801)
    U, u = P.nodal_disp(sec, lb, C, xs)[:2]
    skin = [i for i in range(sec.N) if abs(sec.X[i, 1]) < 1e-9]
    ib = min(skin, key=lambda i: abs(sec.X[i, 0] + 75.0))      # middle of the skin bay
    w = U[:, ib, 1]; w = w / np.abs(w).max()
    # number of half-waves: sign changes of w along x (ignore tiny values at the ends)
    s = np.sign(w[np.abs(w) > 0.02]); hw = int(np.sum(s[1:] != s[:-1])) + 1
    # symmetry about mid-length: w(L-x) = +w(x) (symmetric) or -w(x) (antisymmetric)
    sym = np.abs(w - w[::-1]).max(); anti = np.abs(w + w[::-1]).max()
    # global (rigid) translations: amplitude functions of the 2 G modes
    A = P.amplitudes(lb, C, xs)[:, :sec.nip]
    g = [k for k, f in enumerate(sec.group) if f == 'G']
    Ag = A[:, g]; Aall = np.abs(A).max()
    mean_g = np.abs(np.trapezoid(Ag, xs, axis=0) / L).max() / Aall          # average (net) translation along the length
    max_g = np.abs(Ag).max() / Aall
    out.append(dict(beta=beta, N=round(r['N'], 3), rhoG=round(100 * rho['G'], 2), halfwaves=hw,
                    shape='symmetric' if sym < anti else 'antisymmetric', sym_err=round(min(sym, anti), 3),
                    G_max=round(max_g, 4), G_mean=round(mean_g, 4)))
    print(out[-1], flush=True)
json.dump(out, open('/home/claude/paper2pinnaya/manuscript/results/unsym_G_jump.json', 'w'), indent=1)
