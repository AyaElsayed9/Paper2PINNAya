"""17 Nov 2026 (Aya: "in my first paper I got the right wave number"): number of half-waves m of the critical mode.
m of a mode = the harmonic sin/cos(j pi z / L) carrying the largest share (L2 over z, all shared amplitudes) of its
amplitude field (exact GBT modes are ~98 % one harmonic, check_mode_structure.py).  For each test design:
  gbt        : full-GBT mode 1 (stored teacher); also m of modes 2 and 3
  net_L0 / net_L120 : lowest Ritz vector of the network's 8 candidates (instant, no check)
  nn120      : lowest Ritz vector of the copied modes of the nearest of 120 labelled designs
  check12    : mode 1 of the physics check after 12 steps from the net_L120 start (what the tool reports)
    python scripts/half_waves.py <tag>   -> results/half_waves/<tag>.csv"""
import argparse, csv, sys
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio
from pinnbuck.ritzspace import ritz_forms, ritz_values
from pinnbuck.spectral import lobpcg
from check_start_baselines import scaled


def m_of(Q, z, wz, L):
    """dominant harmonic of the amplitude field Q (nz, nm); returns (m, share)"""
    J = int(min(80, max(8, np.ceil(L / 15.0)))); sw = np.sqrt(wz)[:, None]; q = Q * sw; tot = float((q ** 2).sum())
    best = (0, 0.0)
    for j in range(1, J + 1):
        B = np.stack([np.sin(j * np.pi * z / L), np.cos(j * np.pi * z / L)], 1) * sw
        c, *_ = np.linalg.lstsq(B, q, rcond=None); sh = float(((B @ c) ** 2).sum()) / tot
        if sh > best[1]:
            best = (j, sh)
    return best


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=40)
    a = ap.parse_args(); torch.set_num_threads(2)
    cks = {m: torch.load(ROOT / f"runs/experts/{a.tag}_{m}.pt", weights_only=False) for m in ("physics_L0", "hybrid_L120")}
    epinn.INPUT_COLS = None if cks["hybrid_L120"]["args"].get("inputs", "all") == "all" else [0, 1]
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab]); bank = {}; nets = {}
    out = ROOT / "results/half_waves"; out.mkdir(parents=True, exist_ok=True); fcsv = out / f"{a.tag}.csv"
    done = {r["sid"] for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    cols = ["sid", "bs", "beta", "L", "m_gbt1", "m_gbt2", "m_gbt3", "gap12_pct", "share_gbt1"] + \
           [f"m_{p}" for p in ("net_L0", "net_L120", "nn120", "check12")] + [f"share_{p}" for p in ("net_L0", "net_L120", "nn120", "check12")]
    with open(fcsv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if not done:
            w.writeheader()
        for f in sorted((ROOT / f"data/operators/{a.tag}/test").glob("*.pkl"))[:a.n]:
            op = opio.load(f)
            if op.sid in done:
                continue
            d = epinn.prepare(op); z, wz, L = np.asarray(op.z), np.asarray(op.wz), op.L
            T = np.asarray(op.teacher, dtype=float)
            r = dict(sid=op.sid, bs=round(d.params["bs"], 1), beta=round(d.params["beta"], 3), L=round(L, 1),
                     gap12_pct=round(100 * (op.teacher_eigs[1] / op.teacher_eigs[0] - 1), 3))
            for i in range(3):
                mm, sh = m_of(T[i][0], z, wz, L); r[f"m_gbt{i + 1}"] = mm
                if i == 0:
                    r["share_gbt1"] = round(sh, 3)
            Ys = {}
            for name in ("physics_L0", "hybrid_L120"):
                if name not in nets:
                    nets[name] = epinn.net_from_ckpt(cks[name], d)
                with torch.no_grad():
                    Ys["net_L0" if name == "physics_L0" else "net_L120"] = nets[name](torch.as_tensor(d.x[None]))[0][0, :, d.sel].double()
            j = int(np.argmin(((Xtr - scaled(d.params["bs"], d.params["beta"])) ** 2).sum(1))); sid = tab[j]["sid"]
            if sid not in bank:
                o2 = opio.load(ROOT / f"data/operators/{a.tag}/train/{sid}.pkl"); sel = epinn.design_sel(o2, epinn.layout_labels(o2)); del o2
                bank[sid] = (sel, torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{sid}.pt", weights_only=False)["T"])
            sel, Tn = bank[sid]
            full = torch.zeros(Tn.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=Tn.dtype); full[:, sel] = Tn
            Ys["nn120"] = full[:, d.sel].double()
            for key, Y in Ys.items():
                X = d.fm.color(Y); Kr, Gr = ritz_forms(d.rs, X); mu, Cc = ritz_values(Kr, Gr)
                Q = d.rs.field(d.rs.unflat((Cc[:, 0] @ X)[None]))[0, 0].detach().numpy()
                r[f"m_{key}"], sh = m_of(Q, z, wz, L); r[f"share_{key}"] = round(sh, 3)
            N, Xc, h = lobpcg(d.rs, R=8, iters=12, X0=d.rs.unflat(d.fm.color(Ys["net_L120"])), precond=d.em.solve, tol=0, seed=1)
            Q = d.rs.field(Xc)[0, 0].detach().numpy()
            r["m_check12"], sh = m_of(Q, z, wz, L); r["share_check12"] = round(sh, 3)
            w.writerow(r); fh.flush()
            print({k: r[k] for k in ("sid", "L", "m_gbt1", "m_gbt2", "gap12_pct", "m_net_L0", "m_net_L120", "m_nn120", "m_check12")}, flush=True)
            del d, op
    rows = list(csv.DictReader(open(fcsv)))
    print(f"\n{a.tag}: m (half-waves) right vs full-GBT mode 1, of {len(rows)} designs")
    for p in ("net_L0", "net_L120", "nn120", "check12"):
        ok = sum(r[f"m_{p}"] == r["m_gbt1"] for r in rows); ok3 = sum(r[f"m_{p}"] in (r["m_gbt1"], r["m_gbt2"], r["m_gbt3"]) for r in rows)
        off = [abs(int(r[f"m_{p}"]) - int(r["m_gbt1"])) for r in rows]
        print(f"  {p:9s} exact {ok}/{len(rows)}, one of modes 1-3 {ok3}/{len(rows)}, off by 1: {sum(o == 1 for o in off)}, by >1: {sum(o > 1 for o in off)}")


if __name__ == "__main__":
    main()
