"""17 Nov 2026 (Aya: "are the modes ok? is it the modes that make the PINN better than a GP?"):
quality of the predicted critical mode as INFORMATION, not as a start of the check.

For each test design of an expert, the mode-1 shape of each predictor is its lowest Ritz vector (Rayleigh-Ritz on its
R = 8 candidates in the design's Ritz space).  Compared with the full-GBT mode 1 (the stored teacher, projected on the
shared coordinates):
  MAC1    : MAC of the in-plane nodal displacements (U_X, U_Y of all nodes, all Gauss points along z)
  MACbest : best MAC with any of the three lowest full-GBT modes (close modes can swap)
  span3   : share of the predicted mode lying in the span of the three lowest full-GBT modes
  family  : energy participation of the GBT families (G, D, L, T, S, M; additive, as Eq. 41 of the GBT paper),
            dominant family of predicted vs full GBT, and half the L1 distance of the two participations
Predictors: net_L0 (physics only), net_L120 (hybrid, 120 labels), nn10 / nn120 (modes of the nearest labelled design).
    python scripts/mode_quality.py symmetry_crossply_SS_v2 [--n 40]   -> results/mode_quality/<tag>.csv"""
import argparse, csv, sys
from pathlib import Path
import numpy as np, torch
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from pinnbuck import epinn, opio
from pinnbuck.designspace import Sample
from pinnbuck.families import family_context
from pinnbuck.basis import shared_section
from pinnbuck.ritzspace import ritz_forms, ritz_values
from check_start_baselines import scaled

FAM = ("G", "D", "L", "T", "S", "M")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=40)
    a = ap.parse_args(); torch.set_num_threads(2)
    side = a.tag.split("_")[0]
    import pickle
    sb = pickle.load(open(ROOT / "data" / "shared_basis" / f"basis_{side}.pkl", "rb"))
    nets = {}
    for m in ("physics_L0", "hybrid_L120"):
        nets[m] = torch.load(ROOT / f"runs/experts/{a.tag}_{m}.pt", weights_only=False)
    epinn.INPUT_COLS = None if nets["hybrid_L120"]["args"].get("inputs", "all") == "all" else [0, 1]
    tab = list(csv.DictReader(open(ROOT / f"results/expert_tables/{a.tag}_train.csv")))[:120]
    Xtr = np.array([scaled(float(r["bs"]), float(r["beta"])) for r in tab]); bank = {}
    out = ROOT / "results/mode_quality"; out.mkdir(parents=True, exist_ok=True)
    fcsv = out / f"{a.tag}.csv"
    done = {(r["sid"], r["pred"]) for r in csv.DictReader(open(fcsv))} if fcsv.exists() else set()
    files = sorted((ROOT / f"data/operators/{a.tag}/test").glob("*.pkl"))[:a.n]
    cols = ["sid", "bs", "beta", "pred", "MAC1", "MACbest", "span3", "dom_pred", "dom_gbt", "dist"] + \
           [f"rho_{g}_pred" for g in FAM] + [f"rho_{g}_gbt" for g in FAM]
    built = {}
    with open(fcsv, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if not done:
            w.writeheader()
        for f in files:
            op = opio.load(f); d = epinn.prepare(op)
            if all((op.sid, p) in done for p in ("net_L0", "net_L120", "nn10", "nn120")):
                continue
            s = Sample(**{k: v for k, v in op.info["sample"].items() if k != "sid"}); assert s.sid == op.sid, (s.sid, op.sid)
            with family_context(s) as fc:
                m, sec = fc.m, fc.sec
                wsec = shared_section(m, sec, sb)
                Tb = np.asarray(wsec.T_basis)                       # (nm_full, nm_shared)
                Kcs = m._gbt_Kcs(sec)
                grp = np.array(["S" if g == "A" else g for g in sec.group])
            wz = np.asarray(op.wz); sw = np.sqrt(wz)[:, None]
            Cn = sb.nodal_C; nN = Cn.shape[0] // 3
            def nodal(Q):                                           # Q (3, nz, nm) shared -> (nz, 2N) in-plane nodal
                return (Q[0][:, :sb.n_conv] @ Cn[:2 * nN].T) * sw
            def rho(Q):
                qf = np.einsum("fm,dzm->dzf", Tb, Q)                # full amplitudes and derivatives
                KQ = {dd: np.zeros_like(qf[0]) for dd in range(3)}
                for (dd, e), A in Kcs.items():
                    KQ[dd] += qf[e] @ A.T                           # (nz, f)
                tot = sum(float(np.sum(wz[:, None] * qf[dd] * KQ[dd])) for dd in range(3))
                r = {}
                for g in FAM:
                    k = np.nonzero(grp == g)[0]
                    r[g] = 0.0 if not len(k) else sum(float(np.sum(wz[:, None] * qf[dd][:, k] * KQ[dd][:, k])) for dd in range(3)) / tot
                return r
            T = np.asarray(op.teacher, dtype=float)                  # (3 modes, 3, nz, nm)
            E = [nodal(T[i]) .ravel() for i in range(T.shape[0])]
            Eb = np.stack(E, 1)
            r_g = rho(T[0]); dom_g = max(r_g, key=r_g.get)
            preds = {}
            for name in ("physics_L0", "hybrid_L120"):
                key = "net_L0" if name == "physics_L0" else "net_L120"
                if key not in built:
                    built[key] = epinn.net_from_ckpt(nets[name], d)
                with torch.no_grad():
                    preds[key] = built[key](torch.as_tensor(d.x[None]))[0][0, :, d.sel].double()
            for NL in (10, 120):
                j = int(np.argmin(((Xtr[:NL] - scaled(d.params["bs"], d.params["beta"])) ** 2).sum(1))); sid = tab[j]["sid"]
                if sid not in bank:
                    o2 = opio.load(ROOT / f"data/operators/{a.tag}/train/{sid}.pkl"); sel = epinn.design_sel(o2, epinn.layout_labels(o2)); del o2
                    bank[sid] = (sel, torch.load(ROOT / f"data/exact_targets/{a.tag}/train/{sid}.pt", weights_only=False)["T"])
                sel, Tn = bank[sid]
                full = torch.zeros(Tn.shape[0], int(max(sel.max(), d.sel.max())) + 1, dtype=Tn.dtype); full[:, sel] = Tn
                preds[f"nn{NL}"] = full[:, d.sel].double()
            for key, Y in preds.items():
                if (op.sid, key) in done:
                    continue
                X = d.fm.color(Y)
                Kr, Gr = ritz_forms(d.rs, X); mu, Cc = ritz_values(Kr, Gr)
                x1 = (Cc[:, 0] @ X)[None]
                Q = d.rs.field(d.rs.unflat(x1))[0].detach().numpy()          # (3, nz, nm)
                v = nodal(Q).ravel()
                mac = [float((v @ e) ** 2 / ((v @ v) * (e @ e))) for e in E]
                coef, *_ = np.linalg.lstsq(Eb, v, rcond=None); span = float(1 - np.sum((v - Eb @ coef) ** 2) / np.sum(v ** 2))
                r_p = rho(Q); dom_p = max(r_p, key=r_p.get)
                row = dict(sid=op.sid, bs=round(d.params["bs"], 1), beta=round(d.params["beta"], 3), pred=key, MAC1=round(mac[0], 5),
                           MACbest=round(max(mac), 5), span3=round(span, 5), dom_pred=dom_p, dom_gbt=dom_g,
                           dist=round(0.5 * sum(abs(r_p[g] - r_g[g]) for g in FAM), 5),
                           **{f"rho_{g}_pred": round(r_p[g], 5) for g in FAM}, **{f"rho_{g}_gbt": round(r_g[g], 5) for g in FAM})
                w.writerow(row); fh.flush()
                print(op.sid, key, {k: row[k] for k in ("MAC1", "MACbest", "span3", "dom_pred", "dom_gbt", "dist")}, flush=True)
            del op, d
    rows = list(csv.DictReader(open(fcsv)))
    print(f"\n{a.tag}: {len({r['sid'] for r in rows})} designs | mean (min) | dominant family right | participation distance mean (max)")
    for key in ("net_L0", "net_L120", "nn10", "nn120"):
        rr = [r for r in rows if r["pred"] == key]
        g = lambda c: np.array([float(r[c]) for r in rr])
        print(f"  {key:9s} MAC1 {g('MAC1').mean():.3f} ({g('MAC1').min():.3f})  MACbest {g('MACbest').mean():.3f} ({g('MACbest').min():.3f})  "
              f"span3 {g('span3').mean():.3f} ({g('span3').min():.3f})  dominant {sum(r['dom_pred'] == r['dom_gbt'] for r in rr)}/{len(rr)}  "
              f"dist {g('dist').mean():.3f} ({g('dist').max():.3f})")


if __name__ == "__main__":
    main()
