"""Abaqus mode-1 exports (Aya's fe_modes) vs GBT critical modes: skin deflection MAC."""
import csv, pickle, json, sys
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
H = Path(__file__).parent
FE = H / "src" / "fe_modes"

def read_abq(name):
    rows = list(csv.reader(open(FE / name)))
    ev = float(rows[0][-1].split("=")[1])
    inst = np.array([r[0] for r in rows[1:]]); a = np.array([[float(v) for v in r[2:8]] for r in rows[1:]])
    return dict(eig=ev, inst=inst, xyz=a[:, :3], U=a[:, 3:6])

def gbt_skin_interp(d):
    X = d["X"]; sk = np.array(sorted([n for n in range(len(X)) if abs(X[n, 1]) < 1e-9], key=lambda n: X[n, 0]))
    w = d["U"][:, sk, 1]                                    # (nx, nsk)
    return RegularGridInterpolator((d["xs"], X[sk, 0]), w, bounds_error=False, fill_value=None)

def mac(d, ab):
    m = ab["inst"] == "INST_SKIN"
    pts = np.c_[ab["xyz"][m, 2], ab["xyz"][m, 0]]           # (z, X)
    wg = gbt_skin_interp(d)(pts); wa = ab["U"][m, 1]
    return float((wg @ wa) ** 2 / ((wg @ wg) * (wa @ wa))), float(np.sign(wg @ wa))

if __name__ == "__main__":
    out = {}
    cases = [("mode_Unsym_150.0_1.5_SS", "MODE1_val_bs150_beta1p5_Unsym_PBC.csv"),
             ("mode_Unsym_150.0_1.5_CC", "MODE1_val_bs150_beta1p5_Unsym_CC_PBC.csv"),
             ("longmode_Isotropic_8.0_1.0", "MODE1_val_bs150_beta8_Isotropic_PBC.csv"),
             ("longmode_Sym1_8.0_1.0", "MODE1_val_bs150_beta8_Sym1_PBC.csv"),
             ("longmode_Unsym_12.0_1.0", "MODE1_val_bs150_beta12_Unsym_PBC.csv"),
             ("longmode_Cross-ply_14.0_1.0", "MODE1_val_bs150_beta14_Crossmply_PBC.csv")]
    for gk, af in cases:
        d = pickle.loads((H / "figp" / f"{gk}.pkl").read_bytes()); ab = read_abq(af)
        mc, sg = mac(d, ab)
        out[gk] = dict(abq_file=af, abq_eig=ab["eig"], gbt_N=d["N"], mac_skin=mc)
        print(f"{gk:32s} GBT {d['N']:9.4f}  Abaqus {ab['eig']:9.4f}  skin MAC {mc:.4f}")
    (H / "results" / "abq_mode_mac.json").write_text(json.dumps(out, indent=1))
