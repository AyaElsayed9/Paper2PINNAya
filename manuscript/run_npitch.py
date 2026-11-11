"""Cell of n stringer pitches (Navid, 11 Nov 2026): a real periodic cell of n pitches contains the Bloch
wavenumbers k_B b_s = 2 pi m / n.  n = 2: 0 and pi (the paper).  n = 4 adds pi/2.  Does the critical load of the
long-wave mode (Section 3.5) drop?  The solver is unchanged: only the panel (cross-section) is built here with n
pitches, exactly like solver.make_omega_panel (which is the case n = 2; checked below).

    python run_npitch.py Isotropic 8 4 15  (lay-up, beta, pitches, hx)        -> results/npitch.jsonl
"""
import json, sys, time, resource
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from run_periodic import configure, validated_omega_sample, omega_panel, energy_participation   # noqa: E402

F = Path(__file__).parent / "results" / "npitch.jsonl"


def npitch_panel(m, pan2, n):
    """the panel of make_omega_panel (n = 2) rebuilt with n pitches; same walls, ABD, loads"""
    g = pan2["info"]; bs = g["bs"]; b1 = g["b1"]; ht = g["ht"]; clear = g["clear"]
    w2 = pan2["walls"]; P2 = pan2["nodes"]
    ABD = {p: next(w["ABD"] for w in w2 if w["part"] == p) for p in ("skin", "foot", "web", "crown")}
    crown = next(w for w in w2 if w["part"] == "crown"); h = crown["yload"] - clear
    foot = next(w for w in w2 if w["part"] == "foot"); b3 = None
    c0 = -bs / 2.0
    xf = sorted(abs(P2[w[k], 0] - c0) for w in [(foot["i"], foot["j"])] for k in (0, 1))
    # first foot found is at the far left stringer; use distances to the nearest stringer centre
    cs = [-(k + 0.5) * bs for k in range(2)]
    d = sorted(min(abs(P2[k, 0] - c) for c in cs) for k in (foot["i"], foot["j"]))
    b3, W = d[0], n * bs
    centers = [-(k + 0.5) * bs for k in range(n)]
    xs = [-W, 0.0]
    for c in centers:
        xs += [c - b3 - b1, c - b3, c + b3, c + b3 + b1]
    xs = sorted(set(round(x, 9) for x in xs)); nodes = [(x, 0.0) for x in xs]
    idx = {round(x, 9): k for k, x in enumerate(xs)}; walls = []
    add = lambda i, j, A, part, edges, yload: walls.append(dict(i=i, j=j, ABD=A, part=part, edges=edges, yload=yload))
    feet = set()
    for c in centers:
        feet.add((round(c - b3 - b1, 9), round(c - b3, 9))); feet.add((round(c + b3, 9), round(c + b3 + b1, 9)))
    for k in range(len(xs) - 1):
        if (round(xs[k], 9), round(xs[k + 1], 9)) in feet:
            add(k, k + 1, ABD["foot"], "foot", 2, 0.5 * clear)
        else:
            add(k, k + 1, ABD["skin"], "skin", 1, 0.0)
    cells = []
    for c in centers:
        rl, rr = idx[round(c - b3, 9)], idx[round(c + b3, 9)]
        nodes.append((c - ht, h + clear)); tl = len(nodes) - 1
        nodes.append((c + ht, h + clear)); tr = len(nodes) - 1
        add(rl, tl, ABD["web"], "web", 1, 0.5 * (clear + h + clear)); add(tl, tr, ABD["web"], "crown", 1, h + clear)
        add(tr, rr, ABD["web"], "web", 1, 0.5 * (clear + h + clear)); cells.append([rl, tl, tr, rr])
    info = dict(g); info.update(W=W, cells=cells)
    return dict(nodes=np.asarray(nodes, float), walls=walls, sides=[idx[round(-W, 9)], idx[0.0]], info=info)


def one(lay, beta, n, hx, bs=150.0, bc="SS"):
    t0 = time.perf_counter()
    with configure(OMEGA_SIDE_BC="periodic") as m:
        if lay == "Isotropic":
            pan2 = m.make_omega_panel(bs, "iso", [], []); L = beta * bs
        else:
            s = validated_omega_sample(lay, bs, beta, bc, "periodic"); pan2 = omega_panel(s); L = beta * bs
        pan = npitch_panel(m, pan2, n)
        if n == 2:                                    # the rebuilt panel must be the solver's panel
            assert np.allclose(pan["nodes"], pan2["nodes"]) and len(pan["walls"]) == len(pan2["walls"])
            for a, b in zip(pan["walls"], pan2["walls"]):
                assert (a["i"], a["j"], a["part"], a["edges"]) == (b["i"], b["j"], b["part"], b["edges"])
                assert a["ABD"] is b["ABD"] or all(np.allclose(u, v) for u, v in zip(a["ABD"], b["ABD"]))
                assert abs(a["yload"] - b["yload"]) < 1e-12
        sec = m.omega_gbt_section(pan)
        r = m.omega_gbt_solve(sec, L, bc, hx=hx, nev=4, return_mode=True)   # same hx for every n (default depends on W)
        rho, _ = energy_participation(m, sec, r["lb"], r["coef"])
    return dict(layup=lay, bs=bs, beta=beta, bc=bc, npitch=n, hx=hx, N=r["N"], eigs=[float(e) for e in r["eigs"]], ndof=r["ndof"],
                nmodes=r["nmodes"], rho=rho, time=time.perf_counter() - t0,
                maxrss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)


if __name__ == "__main__":
    lay, beta, n, hx = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    d = one(lay, beta, n, hx)
    with F.open("a") as f:
        f.write(json.dumps(d) + "\n")
    print(lay, beta, "pitches", n, "hx", hx, "N", round(d["N"], 4), "eigs", [round(e, 3) for e in d["eigs"]], "ndof", d["ndof"],
          {g: round(100 * v, 1) for g, v in d["rho"].items()}, f'{d["time"]:.0f}s {d["maxrss_MB"]:.0f}MB', flush=True)
