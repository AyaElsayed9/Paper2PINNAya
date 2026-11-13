"""13 Nov 2026: Abaqus modes 1-4 of the four-pitch cells (Aya's export, PANEL_MODES=1,2,3,4) vs GBT.
For every Abaqus mode: (1) its own Bloch wavenumber, measured from the Abaqus skin deflection alone:
c = sum w(Y) w(Y - b_s) / sum w^2 over the whole periodic width (c = cos(k_B b_s): -1 pi, 0 pi/2, +1 zero);
(2) MAC with each of the four GBT modes; (3) for the double pair, the share of the Abaqus mode that lies in the
plane of the two GBT pi/2 modes (least squares). Writes results/npitch_abq_modes.json."""
import json, pickle, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import LinearNDInterpolator
sys.path.insert(0, str(Path(__file__).parent))
import fp_common as P
import abq_modes as AQ
H = Path(__file__).parent
AQ.FE = H / "src" / "abaqus_npitch" / "fe_modes"
BS = 150.0
out = {}
for lay in ("Isotropic", "Sym1"):
    d = pickle.load(open(P.CACHE / f"npitch_{lay}_4.pkl", "rb")); W = d["W"]
    G = [dict(d, U=m["U"]) for m in d["modes"]]
    res = []
    for k in range(1, 5):
        f = AQ.FE / f"MODE{k}_val_bs150_beta8_{lay}_PBC_n4.csv"
        if not f.exists():
            f = AQ.FE / f"MODE{k}_bs150_beta8_{lay}_PBC_n4.csv"          # name written by the n-pitch run script
        if not f.exists():
            print("missing", f.name); continue
        ab = AQ.read_abq(f.name); m = ab["inst"] == "INST_SKIN"
        x, z, w = ab["xyz"][m, 0], ab["xyz"][m, 2], ab["U"][m, 1]
        xs_ = x - BS; xs_ = np.where(xs_ < x.min() - 1e-6, xs_ + W, xs_)       # periodic shift by one pitch
        ws = LinearNDInterpolator(np.c_[x, z], w)(np.c_[xs_, z]); ok = np.isfinite(ws)
        c = float(np.sum(w[ok] * ws[ok]) / np.sum(w[ok] ** 2))
        macs = [AQ.mac(g, ab)[0] for g in G]
        pts = np.c_[z, x]; wg = np.c_[[AQ.gbt_skin_interp(G[j])(pts) for j in (1, 2)]].T
        coef = np.linalg.lstsq(wg, w, rcond=None)[0]; share = float(1 - np.sum((w - wg @ coef) ** 2) / np.sum(w ** 2))
        res.append(dict(mode=k, eig=ab["eig"], c_abaqus=round(c, 4), mac_gbt=[round(v, 6) for v in macs], share_in_gbt_pi2_plane=round(share, 6)))
        print(lay, res[-1])
    out[lay] = res
json.dump(out, open(H / "results" / "npitch_abq_modes.json", "w"), indent=1)
