"""16 Nov 2026: how much does the constant cell width (transverse pre-buckling force S_ss) change N_cr in the
Paper-2 design space?  Full GBT, periodic, SS, with S_ss (pss_scale = 1, as all Paper-2 labels) and without
(pss_scale = 0, the paper's approximation of a laterally free panel: within 2.4 % of Abaqus for a local mode,
Section 3.5).  Four expert lay-ups x b_s 150 / 285.5 x beta 1.5.  -> results/fixed_width_effect.csv"""
import csv, sys, time
from pathlib import Path
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / "manuscript")); sys.path.insert(0, str(R / "src"))
import fp_common as P
from pinnbuck.designspace import Sample, FAMILIES, VALIDATED_LAYUPS
from pinnbuck.oracle import omega_panel, omega_length
LAY = {"crossply": VALIDATED_LAYUPS["omega"]["Cross-ply"],
       "quasiiso": ([45, -45, 0, 90, 90, 0, -45, 45], [45, -45, 0, 90, 0, -45, 45]),
       "zerorich": ([45, 0, -45, 0, 0, -45, 0, 45], [45, 0, -45, 0, -45, 0, 45]),
       "pm45": VALIDATED_LAYUPS["omega"]["Sym1"]}
rows = []
cm = P.configure(OMEGA_SIDE_BC="periodic"); m = cm.__enter__()
for lay, (sk, st) in LAY.items():
    for bs in (150.0, 285.5):
        s = Sample("omega", dict(bs=bs, beta=1.5, bc="SS", side="periodic"),
                   {"skin": [float(a) for a in sk], "stringer": [float(a) for a in st]}, dict(FAMILIES["omega"]["fixed"]))
        sec = m.omega_gbt_section(omega_panel(s)); L = omega_length(s)
        N = [m.omega_gbt_solve(sec, L, "SS", pss_scale=ps)["N"] for ps in (1.0, 0.0)]
        rows.append(dict(layup=lay, bs=bs, beta=1.5, N_fixed_width=round(N[0], 4), N_no_Sss=round(N[1], 4),
                         change_pct=round(100 * (N[1] / N[0] - 1), 2)))
        print(rows[-1], flush=True)
with open(R / "results" / "fixed_width_effect.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
