"""Sequential timing of GBT and of the spline FSM (periodic), b_s = 150 mm, beta = 1.5, all lay-ups, SS and CC."""
import json, sys, time
from pathlib import Path
sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.oracle import omega_length, omega_panel
from pinnbuck.solver import configure
out = dict(times=[], gbt=[], cases=[], strip_N=[], gbt_N=[], strip_ndof=[])
for lay in ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]:
    for bc in ("SS", "CC"):
        with configure(OMEGA_SIDE_BC="periodic") as m:
            if lay == "Isotropic":
                mk = lambda: m.make_omega_panel(150.0, "iso", [], []); L = 225.0
            else:
                s = validated_omega_sample(lay, 150.0, 1.5, bc, "periodic"); mk = lambda: omega_panel(s); L = omega_length(s)
            t = time.perf_counter(); modes, info = m.omega_spline_eigs(mk(), L, bc, nmodes=2); ts = time.perf_counter() - t
            t = time.perf_counter(); sec = m.omega_gbt_section(mk()); r = m.omega_gbt_solve(sec, L, bc); tg = time.perf_counter() - t
        out["times"].append(ts); out["gbt"].append(tg); out["cases"].append([lay, bc])
        out["strip_N"].append(float(modes[0][0])); out["gbt_N"].append(r["N"]); out["strip_ndof"].append(info.get("ndof"))
        print(lay, bc, f"strip {modes[0][0]:.3f} {ts:.1f}s  gbt {r['N']:.3f} {tg:.1f}s ndof {info.get('ndof')}", flush=True)
Path(__file__).parent.joinpath("results", "strip_timing.json").write_text(json.dumps(out, indent=1))
