"""Physics-reinforced training (Paper 1 hybrid in parametric form): GBT teacher on labelled
designs + Rayleigh-Ritz physics on unlabelled designs (only their GBT energies, ~5 s each).

    python scripts/train_prnn.py --side symmetry --epochs 100 --w-phys 0.1 --unlab 200

Labelled:   data/teacher/<side>/train/*.npz  (scripts/gen_teacher.py)
Unlabelled: data/operators/<side>/unlab/*.pkl (scripts/build_operators.py, no --full)
Test:       data/teacher/<side>/test/*.npz + *.op.pkl

Losses per step:
  teacher  subspace loss (principal angles) + log N_cr (labelled batch)
  physics  -log mu_Ritz of the predicted modes on the design's GBT energy (unlabelled batch);
           mu = 1/N, so this lowers the Rayleigh-Ritz buckling load of the predicted mode
           space towards the true critical load (an upper bound).
--w-phys 0 gives the purely taught network for comparison.
"""
import argparse
import math
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.designspace import Sample                                    # noqa: E402
from pinnbuck.net import (BucklingNet, ModeCodec, cross_section_whitening,  # noqa: E402
                          load_teacher, physics_check, physics_features, subspace_loss)
from pinnbuck.spectral import FourierPanel, rayleigh_ritz                  # noqa: E402
from pinnbuck.teacher import from_modes                                     # noqa: E402


class TorchCodec:
    def __init__(self, codec, nm):
        self.B = torch.as_tensor(codec.basis); self.Wm = torch.as_tensor(codec.Wm); self.c = codec; self.nm = nm

    def decode(self, C, L):
        y = (C.double() @ self.B.T).reshape(C.shape[0], self.nm, -1) / torch.as_tensor(self.c.scale(L))[None, None, :]
        return torch.einsum("mn,knf->kmf", self.Wm, y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--ubatch", type=int, default=4)
    ap.add_argument("--unlab", type=int, default=200)
    ap.add_argument("--w-phys", type=float, default=0.1)
    ap.add_argument("--n-pc", type=int, default=128)
    ap.add_argument("--k", type=int, default=3, help="predicted modes (subspace size)")
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--check", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    torch.set_num_threads(a.threads); torch.manual_seed(a.seed)
    base = ROOT / "data" / "teacher" / a.side
    tr, te = load_teacher(base / "train"), load_teacher(base / "test")
    ops_w = [pickle.load(open(f, "rb")) for f in sorted((base / "test").glob("*.op.pkl"))[:16]]
    Wm, Wp = cross_section_whitening(ops_w)
    codec = ModeCodec.fit(tr, Wm, Wp, a.n_pc); nm = tr[0]["M"].shape[1]; J = tr[0]["J"]
    tc = TorchCodec(codec, nm)
    # unlabelled designs: features from the stored sample, GBT energies streamed from disk
    ufiles = sorted((ROOT / "data" / "operators" / a.side / "unlab").glob("*.pkl"))[:a.unlab]
    ux = []
    for f in ufiles:
        op = pickle.load(open(f, "rb")); sd = op.info["sample"]
        s = Sample(sd["family"], sd["params"], sd["layups"], sd["fixed"])
        ux.append(np.r_[op.features, physics_features(s)].astype(np.float32))
    ux = np.stack(ux) if ux else np.zeros((0, tr[0]["x"].size), np.float32)
    Xtr = np.stack([r["x"] for r in tr]); mu_x, sd_x = Xtr.mean(0), Xtr.std(0) + 1e-6
    nx = lambda X: torch.tensor((X - mu_x) / sd_x, dtype=torch.float32)
    Xtr_t, Xte_t, Xu_t = nx(Xtr), nx(np.stack([r["x"] for r in te])), nx(ux)
    k = a.k
    Ctr = torch.tensor(np.stack([codec.encode(r["M"][:k], r["L"]) for r in tr]), dtype=torch.float32)
    Ltr = torch.tensor(np.log(np.stack([r["lam_full"] for r in tr])), dtype=torch.float32)
    print(f"{len(tr)} labelled, {len(ufiles)} unlabelled, {len(te)} test designs; w_phys = {a.w_phys}", flush=True)
    net = BucklingNet(Xtr_t.shape[1], a.n_pc, k=k, width=256, depth=3)
    opt = torch.optim.Adam(net.parameters(), lr=a.lr, weight_decay=1e-5)
    nsteps = a.epochs * max(1, len(tr) // a.batch)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=nsteps, eta_min=a.lr * 0.02)
    w = torch.tensor([1.0, 0.5, 0.5])[:k]
    gen = np.random.default_rng(a.seed); t0 = time.perf_counter(); fp_cache = {}
    def fp_of(i):
        if i not in fp_cache:
            if len(fp_cache) > 64:
                fp_cache.pop(next(iter(fp_cache)))
            fp_cache[i] = FourierPanel(pickle.load(open(ufiles[i], "rb")), J=J)
        return fp_cache[i]
    step = 0
    for ep in range(a.epochs):
        perm = gen.permutation(len(tr)); lp_hist = []
        for i in range(0, len(tr) - a.batch + 1, a.batch):
            idx = perm[i:i + a.batch]
            C, lg = net(Xtr_t[idx])
            loss = subspace_loss(C, Ctr[idx]).mean() + ((lg - Ltr[idx]) ** 2 * w).sum(1).mean()
            if a.w_phys and len(ufiles):
                lp = 0.0
                for j in gen.choice(len(ufiles), size=a.ubatch, replace=False):
                    fp = fp_of(int(j))
                    Cu, _ = net(Xu_t[j:j + 1])
                    M = tc.decode(Cu[0], fp.L)
                    A = from_modes(fp, M)
                    mu, _ = rayleigh_ritz(A, fp.apply_K(A), fp.apply_G(A), 1)
                    lp = lp - torch.log(torch.clamp(mu[0], min=1e-12)) / a.ubatch
                lp_hist.append(float(lp))
                loss = loss + a.w_phys * lp
            opt.zero_grad(); loss.backward(); opt.step(); sched.step(); step += 1
        if (ep + 1) % max(1, a.epochs // 10) == 0:
            with torch.no_grad():
                _, lt = net(Xte_t)
            err = 100 * (np.exp(lt[:, 0].numpy()) / np.array([r["lam_full"][0] for r in te]) - 1)
            ph = f"  physics: mean log N_Ritz (unlabelled) {np.mean(lp_hist):.2f}" if lp_hist else ""
            print(f"epoch {ep + 1:4d}  test N_cr head mean |err| {np.abs(err).mean():5.2f} %  max {np.abs(err).max():5.1f} %"
                  f"{ph}  ({time.perf_counter() - t0:.0f} s)", flush=True)
    # physics check of the predicted modes on test designs
    with torch.no_grad():
        Cte, lt = net(Xte_t)
    rows = []
    for i, r in enumerate(te[:a.check]):
        opf = Path(r["file"].replace(".npz", ".op.pkl"))
        if not opf.exists():
            continue
        fp = FourierPanel(pickle.load(open(opf, "rb")), J=J)
        chk = physics_check(fp, tc.decode(Cte[i], r["L"]).numpy(), iters=(0, 2, 5))
        rows.append([100 * (v[0] / r["lam_full"][0] - 1) for v in chk.values()])
    rows = np.array(rows)
    print("physics check on test designs (Ritz N_cr error %, median over designs) after 0 / 2 / 5 LOBPCG iterations:",
          np.round(np.median(rows, 0), 2))
    err = 100 * (np.exp(lt[:, 0].numpy()) / np.array([r["lam_full"][0] for r in te]) - 1)
    print(f"final: test N_cr head mean |err| {np.abs(err).mean():.2f} %  max {np.abs(err).max():.1f} %")


if __name__ == "__main__":
    main()
