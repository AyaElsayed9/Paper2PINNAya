"""Discretisation and link-wall checks with periodic side edges."""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/claude/paper2pinnaya/src")
from pinnbuck.designspace import validated_omega_sample          # noqa: E402
from pinnbuck.oracle import omega_length, omega_panel             # noqa: E402
from pinnbuck.solver import configure                             # noqa: E402

OUT = Path(__file__).parent / "results" / "periodic_conv.json"


def panel(m, lay, bs, beta, bc):
    if lay == "Isotropic":
        return m.make_omega_panel(bs, "iso", [], []), beta * bs
    s = validated_omega_sample(lay, bs, beta, bc, "periodic")
    return omega_panel(s), omega_length(s)


def solve(lay, bs, beta, bc, seg=25.0, hfac=1.0, link_k=1.0e3, stiff_d66=False):
    with configure(OMEGA_SIDE_BC="periodic") as m:
        old = (m.OMEGA_GBT_SEG, m.OMEGA_GBT_LINK_K)
        m.OMEGA_GBT_SEG, m.OMEGA_GBT_LINK_K = seg, link_k
        orig = m.omega_gbt_geometry
        if stiff_d66:
            def geo(mm, pan, link_k=1.0e3):
                X, walls, sides = orig(mm, pan, link_k)
                sk = next(w for w in pan["walls"] if w["part"] == "skin")
                D0 = sk["ABD"][2]
                for w in walls:
                    if w[6] == "link":
                        w[2][5, 5] = link_k * D0[2, 2]
                return X, walls, sides
            m.omega_gbt_geometry = geo
        try:
            pan, L = panel(m, lay, bs, beta, bc)
            sec = m.omega_gbt_section(pan)
            hx = min(15.0, max(2.0, 2 * bs / 40.0), L / 10.0) * hfac
            return m.omega_gbt_solve(sec, L, bc, hx=hx)["N"]
        finally:
            m.omega_gbt_geometry = orig
            m.OMEGA_GBT_SEG, m.OMEGA_GBT_LINK_K = old


if __name__ == "__main__":
    res = {}
    cases = [("Isotropic", 150.0, 1.5, "SS"), ("Sym2", 150.0, 3.0, "SS"), ("Unsym", 150.0, 1.5, "CC"),
             ("Sym1", 285.5, 3.0, "CC"), ("Sym1", 285.5, 1.5, "CC")]
    for c in cases:
        t = time.perf_counter()
        r = dict(default=solve(*c), spans_half=solve(*c, hfac=0.5), seg40=solve(*c, seg=40.0))
        if c == ("Sym1", 285.5, 1.5, "CC"):
            r["spans_quarter"] = solve(*c, hfac=0.25)
        res["|".join(map(str, c))] = r
        print(c, r, f"{time.perf_counter() - t:.0f} s", flush=True)
    c = ("Isotropic", 150.0, 1.5, "SS")
    res["linkk"] = {str(k): solve(*c, link_k=k) for k in (1e2, 1e3, 1e4)}
    res["link_stiffD66"] = solve(*c, stiff_d66=True)
    print(res["linkk"], res["link_stiffD66"], flush=True)
    OUT.write_text(json.dumps(res, indent=1))
