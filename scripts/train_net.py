"""Train the fast buckling network of the omega family on GBT teacher data.

    python scripts/train_net.py --side symmetry --epochs 400 --ensemble 5 --threads 8

Data: data/teacher/<side>/train/*.npz and test/*.npz (+ test/*.op.pkl for the physics
check), from scripts/gen_teacher.py.  Whitening uses the operators of the test folder
only through the mean cross-section stiffness (no labels).
Output: runs/net_<side>/ (models, codec, results.json) and a summary:
  N_cr (head)      the fast prediction, no GBT involved
  N_cr (physics)   Ritz value of the predicted modes on the design's GBT operators
  residual eta     energy residual of that mode (the physics check)
  ensemble CoV     spread of the ensemble heads
"""
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pinnbuck.net import (BucklingNet, ModeCodec, cross_section_whitening, load_teacher,   # noqa: E402
                          physics_check, subspace_loss)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", default="symmetry", choices=["symmetry", "periodic"])
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--n-pc", type=int, default=256)
    ap.add_argument("--width", type=int, default=512)
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--ensemble", type=int, default=5)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--check", type=int, default=16, help="test designs with the physics check")
    a = ap.parse_args()
    torch.set_num_threads(a.threads)
    base = ROOT / "data" / "teacher" / a.side
    tr, te = load_teacher(base / "train"), load_teacher(base / "test")
    ops_w = [pickle.load(open(f, "rb")) for f in sorted((base / "test").glob("*.op.pkl"))[:16]]
    Wm, Wp = cross_section_whitening(ops_w)
    codec = ModeCodec.fit(tr, Wm, Wp, a.n_pc)
    nm = tr[0]["M"].shape[1]
    enc = lambda r: codec.encode(r["M"].astype(float), r["L"])
    Xtr = np.stack([r["x"] for r in tr]); Xte = np.stack([r["x"] for r in te])
    mu_x, sd_x = Xtr.mean(0), Xtr.std(0) + 1e-6                   # input standardisation
    Xtr = torch.tensor((Xtr - mu_x) / sd_x, dtype=torch.float32); Xte = torch.tensor((Xte - mu_x) / sd_x, dtype=torch.float32)
    Ctr = torch.tensor(np.stack([enc(r) for r in tr]), dtype=torch.float32)
    Ltr = torch.tensor(np.log(np.stack([r["lam_full"] for r in tr])), dtype=torch.float32)
    # PCA representation error of the teacher modes (energy-like norm)
    rep = np.mean([1 - np.sum(enc(r) ** 2) / np.sum((np.einsum("mn,knf->kmf", Wp, r["M"].astype(float))
                   * codec.scale(r["L"])[None, None, :]) ** 2) for r in te])
    print(f"{len(tr)} training / {len(te)} test designs, {nm} shared modes, {a.n_pc} PCA coordinates "
          f"(test modes outside the PCA space: {100 * rep:.2f} % of the energy-like norm)", flush=True)
    w = torch.tensor([1.0, 0.5, 0.5])
    models = []
    t0 = time.perf_counter()
    for e in range(a.ensemble):
        torch.manual_seed(e)
        net = BucklingNet(Xtr.shape[1], a.n_pc, width=a.width, depth=a.depth)
        opt = torch.optim.Adam(net.parameters(), lr=a.lr)
        steps = a.epochs * max(1, len(tr) // a.batch)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=a.lr * 0.01)
        gen = torch.Generator().manual_seed(100 + e)
        for ep in range(a.epochs):
            perm = torch.randperm(len(tr), generator=gen)
            for i in range(0, len(tr) - a.batch + 1, a.batch):
                idx = perm[i:i + a.batch]
                C, lg = net(Xtr[idx])
                loss = subspace_loss(C, Ctr[idx]).mean() + ((lg - Ltr[idx]) ** 2 * w).sum(1).mean()
                opt.zero_grad(); loss.backward(); opt.step(); sched.step()
            if (ep + 1) % max(1, a.epochs // 5) == 0:
                with torch.no_grad():
                    _, lt = net(Xte)
                    err = 100 * (torch.exp(lt[:, 0]).numpy() / np.array([r["lam_full"][0] for r in te]) - 1)
                print(f"  model {e}  epoch {ep + 1}  loss {float(loss):.4f}  test N_cr head: mean |err| "
                      f"{np.abs(err).mean():.2f} %  max {np.abs(err).max():.2f} %  ({time.perf_counter() - t0:.0f} s)",
                      flush=True)
        models.append(net.eval())
    # ---- evaluation ---------------------------------------------------------------
    with torch.no_grad():
        outs = [m(Xte) for m in models]
    lam_heads = np.stack([torch.exp(o[1][:, 0]).numpy() for o in outs])     # (E, n_te)
    lam_head = lam_heads.mean(0); cov = lam_heads.std(0) / lam_head
    ref = np.array([r["lam_full"][0] for r in te])
    # ---- reliability: split-conformal calibration of the ensemble (replaces -15 % / +3 %)
    conf = None
    cal_dir = base / "calib"
    if cal_dir.exists() and any(cal_dir.glob("*.npz")) and a.ensemble > 1:
        ca = load_teacher(cal_dir)
        Xca = torch.tensor((np.stack([r["x"] for r in ca]) - mu_x) / sd_x, dtype=torch.float32)
        with torch.no_grad():
            lc = np.stack([m(Xca)[1][:, 0].numpy() for m in models])        # log N
        sc = np.abs(lc.mean(0) - np.log([r["lam_full"][0] for r in ca])) / (lc.std(0) + 1e-3)
        lt = np.log(lam_heads)
        conf = {}
        for alpha in (0.10, 0.05):
            n = len(sc); q = np.quantile(sc, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n))
            half = q * (lt.std(0) + 1e-3)                                     # half-width in log N
            inside = np.abs(lt.mean(0) - np.log(ref)) <= half
            conf[alpha] = dict(q=float(q), coverage=float(inside.mean()), half_width_pct=float(100 * np.mean(np.expm1(half))))

    res = []
    from pinnbuck.spectral import FourierPanel
    for i, r in enumerate(te):
        row = dict(sid=r["sid"], params=r["params"], lam_full=float(ref[i]), lam_head=float(lam_head[i]),
                   cov=float(cov[i]))
        opf = Path(r["file"].replace(".npz", ".op.pkl"))
        if i < a.check and opf.exists():
            op = pickle.load(open(opf, "rb"))
            fp = FourierPanel(op, J=r["J"])
            Mp = codec.decode(outs[0][0][i].numpy().astype(np.float64), r["L"], nm)
            chk = physics_check(fp, Mp)
            row.update(check={str(k): v for k, v in chk.items()})
            row.update(lam_phys=chk[max(chk)][0], eta=chk[max(chk)][1])
        res.append(row)
    t_pred = time.perf_counter()
    with torch.no_grad():
        for _ in range(200):
            models[0](Xte)
    t_pred = (time.perf_counter() - t_pred) / 200 / len(te)
    out = ROOT / "runs" / f"net_{a.side}"; out.mkdir(parents=True, exist_ok=True)
    for e, m in enumerate(models):
        torch.save(m.state_dict(), out / f"model{e}.pt")
    pickle.dump(dict(codec=codec, mu_x=mu_x, sd_x=sd_x), open(out / "codec.pkl", "wb"))
    json.dump(dict(args=vars(a), results=res, conformal=conf), open(out / "results.json", "w"), indent=1)
    eh = np.array([100 * (r["lam_head"] / r["lam_full"] - 1) for r in res])
    if conf:
        for alpha, c in conf.items():
            print(f"conformal {100 * (1 - alpha):.0f} % interval (calibrated on {len(ca)} designs): coverage on test "
                  f"{100 * c['coverage']:.0f} %, mean half-width {c['half_width_pct']:.2f} %")
    print(f"\nprediction time {1e3 * t_pred:.3f} ms per design (network only)")
    print(f"N_cr (head, ensemble mean) vs full GBT: mean |err| {np.abs(eh).mean():.2f} %  max {np.abs(eh).max():.2f} %")
    chk = [r for r in res if "lam_phys" in r]
    if chk:
        ep_ = np.array([100 * (r["lam_phys"] / r["lam_full"] - 1) for r in chk])
        print(f"N_cr (physics check, 5 iterations from the predicted modes) on {len(chk)} designs: "
              f"mean {ep_.mean():+.2f} %  max {ep_.max():+.2f} %  min {ep_.min():+.2f} %")
        print("  per design: head err %, ensemble CoV %, then physics check after 0/2/5 LOBPCG iterations: err % (eta)")
        for r in chk:
            c = "  ".join(f"{100 * (v[0] / r['lam_full'] - 1):+8.2f} ({v[1]:.1e})" for v in r["check"].values())
            print(f"   {r['params']['bc']} bs {r['params']['bs']:6.1f} beta {r['params']['beta']:.2f}:"
                  f"  {100 * (r['lam_head'] / r['lam_full'] - 1):+6.2f}  {100 * r['cov']:5.2f}   {c}")
    print(f"results in {out}")


if __name__ == "__main__":
    main()
