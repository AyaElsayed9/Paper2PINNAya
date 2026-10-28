"""Run Aya's own delivered solver file (outputs/AYA-GBT-V32.py) on the cases shown in the
figures and compare with the values used in the manuscript (repo copy of V32)."""
import json, sys, importlib.util
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
import pinnbuck.solver as PS
PS.SOLVER_PATH = __import__("pathlib").Path("/mnt/user-data/outputs/AYA-GBT-V32.py")
from pinnbuck.solver import configure
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import omega_length, omega_panel
cur = {(c["layup"], c["bs"], c["beta"], c["bc"]): c["N"] for c in map(json.loads, open("results/periodic_curves.jsonl"))}
cases = [("Isotropic", 150.0, 1.5, "SS"), ("Unsym", 150.0, 1.5, "SS"), ("Unsym", 150.0, 1.5, "CC"),
         ("Sym1", 285.5, 3.0, "CC"), ("Cross-ply", 285.5, 0.75, "SS")]
for lay, bs, beta, bc in cases:
    with configure(OMEGA_SIDE_BC="periodic") as m:
        assert m.__file__.endswith("AYA-GBT-V32.py")
        if lay == "Isotropic":
            pan = m.make_omega_panel(bs, "iso", [], []); L = beta * bs
        else:
            s = validated_omega_sample(lay, bs, beta, bc, "periodic"); pan = omega_panel(s); L = omega_length(s)
        N = m.omega_gbt_solve(m.omega_gbt_section(pan), L, bc)["N"]
    ref = cur[(lay, bs, beta, bc)]
    print(f"{lay:10s} bs={bs:6.1f} beta={beta:4.2f} {bc}: user file {N:.6f}  manuscript {ref:.6f}  diff {abs(N-ref)/ref:.1e}")
