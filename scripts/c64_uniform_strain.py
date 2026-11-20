"""20 Nov 2026, Prof. Mittelstedt's comment C64 on the GBT paper: the end load "should be dictated by the
longitudinal strain".  The paper applies the same line load (1 N/mm) to skin, webs and crowns (2 N/mm on a foot =
skin + flange), as Paper 1 and the earlier Abaqus models.  Here the alternative "uniform end strain" distribution:
each wall's line load proportional to its axial membrane stiffness A11 (foot: its laminate), same lever arms.
Solver unchanged: only the section's load intensities sec.qp are set before the solve.  Compared quantity: the
critical mean axial force per unit panel width, P_cr / W.   -> results/c64_uniform_strain.csv"""
import csv, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.families import family_context


def main():
    rows = []
    for lay in ("Cross-ply", "Sym1", "Sym2", "Unsym"):
        for side in ("periodic", "symmetry"):
            for bs in (150.0, 285.5):
                for beta in (1.5, 3.0):
                    s = validated_omega_sample(lay, bs, beta, "SS", side)
                    with family_context(s) as fc:
                        sec = fc.sec; W = 2 * bs
                        q0 = sec.qp.copy()
                        r0 = fc.solve(nev=1); P0 = float(r0["N"]) * float((sec.wp * q0).sum()) / W
                        A11 = sec.Cp[:, 0, 0]
                        part = np.array([sec.wpart[sg["w"]] for sg in sec.S for _ in range(4)])
                        skin = part == "skin"
                        qs = A11 / np.median(A11[skin])                                 # skin walls keep 1 N/mm
                        sec.qp = qs
                        r1 = fc.solve(nev=1); P1 = float(r1["N"]) * float((sec.wp * qs).sum()) / W
                        sec.qp = q0
                    row = dict(layup=lay, side=side, bs=bs, beta=beta, Ncr_uniform_load=float(r0["N"]),
                               Pcr_W_uniform_load=P0, Pcr_W_uniform_strain=P1, diff_pct=100 * (P1 / P0 - 1),
                               q_web=float(np.median(qs[part == "web"])), q_foot=float(np.median(qs[part == "foot"])))
                    rows.append(row); print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    with open(ROOT / "results/c64_uniform_strain.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
