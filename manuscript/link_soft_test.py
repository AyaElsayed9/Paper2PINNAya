"""Sensitivity of N_cr to the soft factor of the link walls (A11, D11, D66 = factor x skin).
Builds a patched copy of the solver (the solver file itself is not modified) and writes results/link_soft.json."""
import sys, importlib.util, json, os
src = open('/home/claude/paper2pinnaya/solver/aya_gbt_v32.py').read()
for line in ["CL[0, 0] = 1e-3 * A0[0, 0]", "CL[3, 3] = 1e-3 * D0[0, 0]", "CL[5, 5] = 1e-3 * D0[2, 2]"]:
    assert src.count(line) == 1; src = src.replace(line, line.replace("1e-3", "LINK_SOFT"))
src = src.replace("OMEGA_GBT_LINK_K = 1.0e3", "OMEGA_GBT_LINK_K = 1.0e3\nLINK_SOFT = 1e-3", 1)
open('/tmp/aya_soft_copy.py', 'w').write(src)
sys.argv = ['aya_soft_copy.py']
spec = importlib.util.spec_from_file_location('aya_soft_copy', '/tmp/aya_soft_copy.py'); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
m.OMEGA_SIDE_BC = 'periodic'
out = {}
for name, mt, sk, st in [('Isotropic', 'iso', [], []), ('Unsym', 'cfrp', [0, 30, 60, 90] * 2, [0, 30, 60, 90, 0, 30, 60])]:
    out[name] = {}
    for soft in [1e-2, 1e-3, 1e-4, 0.0]:
        m.LINK_SOFT = soft
        sec = m.omega_gbt_section(m.make_omega_panel(150.0, mt, sk, st))
        out[name][f"{soft:g}"] = round(m.omega_gbt_solve(sec, 225.0, 'SS', nev=1)['N'], 4)
print(out); json.dump(out, open('results/link_soft.json', 'w'), indent=1)
