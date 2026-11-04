"""Abaqus cases for the extended design space of Paper 2 (4 Nov 2026).

Positions from a scrambled Sobol sequence with a fixed seed (not chosen by looking at any GBT or
network result); lay-ups assigned in a fixed rotation.  Run once: python make_cases.py > cases.csv
  A  periodic, SS : 6 cells with b_s in [124.25, 150] mm, beta in [0.75, 3.5]
                    + 2 cells with b_s in [150, 285.5] mm, beta in [3.0, 3.5]  (beta above the old 3.0)
  B  symmetry, SS : 6 cells with b_s in [124.25, 150] mm, beta in [0.75, 3.5], lay-ups other than
                    cross-ply (Paper 1 has only cross-ply there)
  C  periodic, CC : 4 cells (optional), b_s in [124.25, 150] mm
The first cell of A and of B is a control (b_s 150, beta 1.5, Sym1): already run in Abaqus for the
periodic-GBT paper, it shows that the run on Aya's computer reproduces the earlier model.
"""
from scipy.stats import qmc


def pts(n, lo, hi, seed):
    u = qmc.Sobol(d=2, scramble=True, seed=seed).random_base2(3)[:n]
    return qmc.scale(u, lo, hi)


rows = [("A", "PBC", "SS", 150.0, 1.5, "Sym1"),      # control: one of the 30 validated periodic cells
        ("B", "SYM", "SS", 150.0, 1.5, "Sym1")]      # control: one of the 60 validated symmetry cells
for (bs, beta), lay in zip(pts(6, [124.25, 0.75], [150.0, 3.5], 2026), ["Cross-ply", "Sym1", "Sym2", "Unsym"] * 2):
    rows.append(("A", "PBC", "SS", bs, beta, lay))
for (bs, beta), lay in zip(pts(2, [150.0, 3.0], [285.5, 3.5], 2027), ["Unsym", "Sym2"]):
    rows.append(("A", "PBC", "SS", bs, beta, lay))
for (bs, beta), lay in zip(pts(6, [124.25, 0.75], [150.0, 3.5], 2028), ["Sym1", "Sym2", "Unsym"] * 2):
    rows.append(("B", "SYM", "SS", bs, beta, lay))
for (bs, beta), lay in zip(pts(4, [124.25, 0.75], [150.0, 3.5], 2029), ["Cross-ply", "Sym1", "Sym2", "Unsym"]):
    rows.append(("C", "PBC", "CC", bs, beta, lay))
print("set,sides,bc,bs,beta,layup")
for s, sd, bc, bs, beta, lay in rows:
    print(f"{s},{sd},{bc},{bs:.2f},{beta:.2f},{lay}")
