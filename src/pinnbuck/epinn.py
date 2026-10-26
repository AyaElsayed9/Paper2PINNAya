"""Physics-only parametric PINN of the omega panels: energy-whitened Ritz network trained by
an energy natural-gradient step with Ritz line search ("projected LOBPCG").

Output space.  For every end condition the network returns R candidate modes y_r in the
energy-whitened coordinates of the dictionary Ritz space (`ritzspace.py`): the physical
amplitudes are a = C^-1 y with C^T C ~ K (the design's own stiffness, wavenumber blocks),
so |y|^2 ~ 2 U for every design.  The coordinates are fixed (same dictionary labels and
cross-section modes for every design), only C changes with the design.

Loss.  Physics only: the Ritz values mu_1 >= mu_2 >= ... of G x = mu K x over the R
candidates (N_i = 1/mu_i, upper bounds of the GBT critical loads in the Ritz space), with
the Ky Fan sum of the k lowest loads (`physics_loss`), plus the N_cr head, which learns the
network's own Ritz values (no external labels anywhere).

Optimiser (`ritz_ng_step`).  In whitened coordinates the energy natural gradient
(Mueller & Zeinhofer 2023; metric = energy Hessian) is the plain residual
r = G~ y - mu y.  For each design of the batch the best update of the candidates in
span{y, r, p} (p: previous change) follows exactly from a Rayleigh-Ritz step (the line
search of ENGD, in the form of LOBPCG); the network is then moved towards these updated
candidates by a few inner steps on the subspace distance (Gauss-Newton projection of the
function-space step onto the network).  No converged eigen-solution is ever computed:
every target is one physics step away from the network's current prediction.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import torch
from torch import nn

from .ritzspace import EnergyMetric, FlatMetric, RitzSpace, ritz_forms, ritz_values


# ----------------------------------------------------------------------------- designs
def wave_features(sample):
    """Closed-form wave content of the plates (skin bay, web, crown as long orthotropic
    plates): log critical half-wavelength lambda = b (D11 / D22)^(1/4) and log of the
    number of half-waves L / lambda that fit in the panel."""
    from .oracle import omega_length, omega_panel
    pan = omega_panel(sample); nd = pan["nodes"]; L = omega_length(sample)
    best = {}
    for w in pan["walls"]:
        part = w.get("part", "")
        if part in ("skin", "web", "crown"):
            P0, P1 = np.asarray(nd[w["i"]], float), np.asarray(nd[w["j"]], float)
            b = float(np.hypot(*(P1 - P0)))
            if b >= best.get(part, (0.0, None))[0]:
                best[part] = (b, w["ABD"][2])
    lam = [b * (D[0, 0] / D[1, 1]) ** 0.25 for b, D in (best[p] for p in ("skin", "web", "crown"))]
    return np.log(np.r_[lam, [L / l for l in lam]]).astype(np.float32)


def design_inputs(op):
    """Design vector (scaled geometry, flags, lamination parameters) + log closed-form
    plate / column loads (`net.physics_features`) + plate wave content (`wave_features`)."""
    from .designspace import FAMILIES, Sample
    from .net import physics_features
    s = op.info["sample"]
    smp = Sample(s["family"], s["params"], s["layups"], s.get("fixed") or dict(FAMILIES["omega"]["fixed"]))
    return np.r_[op.features, physics_features(smp), wave_features(smp)].astype(np.float32)


@dataclass
class Design:
    sid: str
    params: dict
    bc: str
    x: np.ndarray
    rs: RitzSpace
    fm: FlatMetric
    lam_full: float
    P: torch.Tensor = None              # previous change of the candidates (LOBPCG momentum)
    hist: list = field(default_factory=list)


def prepare(op, ridge=1e-8) -> Design:
    rs = RitzSpace(op, raw=True)
    fm = FlatMetric(rs, EnergyMetric(rs, groups="labels", ridge=ridge))
    return Design(op.sid, op.params, op.params["bc"], design_inputs(op), rs, fm, float(op.lam_full))


def _color_T(fm: FlatMetric, V):
    """C^-T v (dual vectors to whitened coordinates)."""
    Vp = V[:, fm.perm]
    parts = [torch.linalg.solve_triangular(Lc.to(V.dtype), Vp[:, s:e].T, upper=False).T for s, e, Lc in fm.blocks]
    return torch.cat(parts, 1)[:, fm.inv]


FlatMetric.color_T = _color_T


# ----------------------------------------------------------------------------- network
class Scaler(nn.Module):
    def __init__(self, X):
        super().__init__()
        X = torch.as_tensor(np.asarray(X), dtype=torch.float32)
        sd = X.std(0); sd[sd < 1e-6] = 1.0
        self.register_buffer("m", X.mean(0)); self.register_buffer("s", sd)

    def forward(self, x):
        return (x - self.m) / self.s


class EnergyRitzPINN(nn.Module):
    """design vector -> R whitened candidate modes (R, n) and log N_1..k."""

    def __init__(self, X_train, n_out, R=8, k=3, width=256, depth=4):
        super().__init__()
        self.R, self.k, self.n = R, k, n_out
        self.scale = Scaler(X_train)
        layers, d = [], X_train.shape[1]
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.SiLU()]; d = width
        self.body = nn.Sequential(*layers)
        self.hY = nn.Linear(width, R * n_out)
        self.hN = nn.Linear(width, k)
        with torch.no_grad():
            self.hN.weight.mul_(0.01); self.hN.bias.fill_(math.log(10.0))

    def forward(self, x):
        h = self.body(self.scale(x))
        return self.hY(h).reshape(-1, self.R, self.n), self.hN(h)


def function_codes(rs: RitzSpace):
    """Physical code of every longitudinal function of the fixed layout (class c, label f):
    [log half-wave / 100 mm, sin, cos, left end, right end, end-spline index / 8, class one-hot]."""
    nC = len(rs.k); codes, slots = [], []
    for c, labs in enumerate(rs.labels):
        for f, lab in enumerate(labs):
            if lab.startswith("wave"):
                v = [math.log(float(lab[4:-1]) / 100.0), float(lab[-1] == "s"), float(lab[-1] == "c"), 0.0, 0.0, 0.0]
            else:
                v = [0.0, 0.0, 0.0, float(lab[3] == "L"), float(lab[3] == "R"), int(lab[4:]) / 8.0]
            codes.append(v + list(np.eye(nC)[c])); slots.append((c, f))
    return np.asarray(codes, np.float32), slots


class FunctionDecoderPINN(nn.Module):
    """Nonlinear decoder (parameter encoder + coordinate encoder, as in P2INN): for design
    features h and the code e_f of every longitudinal function f, an MLP gives the
    cross-section vectors of the R candidates on that function.  Unlike a linear output
    layer (outputs confined to a fixed span; Kolmogorov n-width barrier), the outputs can
    follow modes whose wavelength and localisation change with the design."""

    def __init__(self, X_train, rs: RitzSpace, R=8, k=3, width=256, depth=4, dec_width=256, dec_depth=2):
        super().__init__()
        self.R, self.k, self.n = R, k, rs.n
        self.scale = Scaler(X_train)
        layers, d = [], X_train.shape[1]
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.SiLU()]; d = width
        self.body = nn.Sequential(*layers)
        codes, slots = function_codes(rs)
        self.register_buffer("codes", torch.as_tensor(codes))
        self.kmax = max(rs.k)
        offs = np.r_[0, np.cumsum([kc * nb for kc, nb in zip(rs.k, rs.nb)])]
        slot_of = {cf: i for i, cf in enumerate(slots)}
        gidx = np.zeros(rs.n, np.int64)
        for c in range(len(rs.k)):
            for kk in range(rs.k[c]):
                for f in range(rs.nb[c]):
                    gidx[offs[c] + kk * rs.nb[c] + f] = slot_of[(c, f)] * self.kmax + kk
        self.register_buffer("gidx", torch.as_tensor(gidx))
        self.inH = nn.Linear(width, dec_width); self.inE = nn.Linear(codes.shape[1], dec_width, bias=False)
        dec = []
        for _ in range(dec_depth - 1):
            dec += [nn.SiLU(), nn.Linear(dec_width, dec_width)]
        self.dec = nn.Sequential(*dec, nn.SiLU(), nn.Linear(dec_width, R * self.kmax))
        self.hN = nn.Linear(width, k)
        with torch.no_grad():
            self.hN.weight.mul_(0.01); self.hN.bias.fill_(math.log(10.0))

    def forward(self, x):
        h = self.body(self.scale(x))                                       # (B, w)
        z = self.inH(h)[:, None, :] + self.inE(self.codes)[None, :, :]     # (B, n_f, dw)
        out = self.dec(z).reshape(h.shape[0], -1, self.R, self.kmax)       # (B, n_f, R, kmax)
        out = out.permute(0, 2, 1, 3).reshape(h.shape[0], self.R, -1)      # (B, R, n_f * kmax)
        return out[:, :, self.gidx], self.hN(h)


class HybridPINN(FunctionDecoderPINN):
    """Decoder output plus a linear head on the same design features: the linear part
    learns the common span fast, the decoder adds what lies outside a fixed span."""

    def __init__(self, X_train, rs, R=8, k=3, width=256, depth=4, lin_rank=64):
        super().__init__(X_train, rs, R=R, k=k, width=width, depth=depth)
        self.lin_in = nn.Linear(width, lin_rank)
        self.lin_out = nn.Linear(lin_rank, R * rs.n)
        with torch.no_grad():
            self.lin_out.weight.mul_(0.1); self.lin_out.bias.zero_()

    def forward(self, x):
        h = self.body(self.scale(x))
        z = self.inH(h)[:, None, :] + self.inE(self.codes)[None, :, :]
        out = self.dec(z).reshape(h.shape[0], -1, self.R, self.kmax).permute(0, 2, 1, 3).reshape(h.shape[0], self.R, -1)
        lin = self.lin_out(self.lin_in(h)).reshape(h.shape[0], self.R, self.n)
        return out[:, :, self.gidx] + lin, self.hN(h)


def build_net(arch, X_train, rs, R=8, k=3, width=256, depth=4):
    if arch == "hybrid":
        return HybridPINN(X_train, rs, R=R, k=k, width=width, depth=depth)
    if arch == "decoder":
        return FunctionDecoderPINN(X_train, rs, R=R, k=k, width=width, depth=depth)
    return EnergyRitzPINN(X_train, rs.n, R=R, k=k, width=width, depth=depth)


# ----------------------------------------------------------------------------- physics
def ritz(d: Design, Y):
    """Ritz values mu (descending) and combinations of the whitened candidates Y (R, n)."""
    X = d.fm.color(Y.double())
    return ritz_values(*ritz_forms(d.rs, X))


def physics_loss(d: Design, Y, logN, k=3, w_head=1.0):
    """Ky Fan sum of the k largest Ritz values (each normalised: gradient of -sum log mu_i)
    and the N_cr head trained on the network's own Ritz values."""
    mu, _ = ritz(d, Y)
    kk = min(k, int((mu > 0).sum()))
    loss = -(mu[:kk] / mu[:kk].detach().abs()).sum()
    tgt = torch.log(1.0 / mu[:kk].detach()).float()
    loss = loss + w_head * ((logN[:kk] - tgt) ** 2).sum()
    return loss, mu.detach()


@torch.no_grad()
def ritz_ng_target(d: Design, Y, R):
    """One energy natural-gradient step with Ritz line search for design d: the best R
    vectors of span{Ritz(Y), residuals, previous change} (whitened coordinates)."""
    Y = Y.double()
    X = d.fm.color(Y)
    A = d.rs.unflat(X)
    KX = d.rs.flat(d.rs.apply_K(A)); GX = d.rs.flat(d.rs.apply_G(A))
    Kr = X @ KX.T; Gr = X @ GX.T
    mu, C = ritz_values(0.5 * (Kr + Kr.T), 0.5 * (Gr + Gr.T))
    m = min(R, C.shape[1])
    Yr, KXr, GXr = C[:, :m].T @ Y, C[:, :m].T @ KX, C[:, :m].T @ GX
    W = d.fm.color_T(GXr - mu[:m, None] * KXr)                            # whitened residuals
    W = W / W.norm(dim=1, keepdim=True).clamp(min=1e-300)
    blocks = [Yr, W] + ([d.P] if d.P is not None else [])
    S = torch.cat(blocks, 0)
    XS = d.fm.color(S)
    extra = XS[m:]
    AE = d.rs.unflat(extra)
    KS = torch.cat([KXr, d.rs.flat(d.rs.apply_K(AE))], 0); GS = torch.cat([GXr, d.rs.flat(d.rs.apply_G(AE))], 0)
    Ks = XS @ KS.T; Gs = XS @ GS.T
    mu2, C2 = ritz_values(0.5 * (Ks + Ks.T), 0.5 * (Gs + Gs.T))
    T = C2[:, :R].T @ S
    return T, mu2[:R], mu[:m]


def subspace_fit_loss(Y, T, w):
    """sum_i w_i (1 - |P_Y t_i|^2), t_i unit target vectors, P_Y the orthogonal projector
    onto span(Y) (whitened coordinates ~ energy metric)."""
    Q, _ = torch.linalg.qr(Y.double().T)                                  # (n, R)
    Tn = T / T.norm(dim=1, keepdim=True)
    c = Tn @ Q                                                            # (R_T, R)
    return (w * (1.0 - (c * c).sum(1))).sum()


# ----------------------------------------------------------------------------- training
class DesignPool:
    """Training designs streamed from operator files: inputs x, references (report only)
    and the LOBPCG momentum are kept for all; the energy spaces (~60 MB each) are built on
    demand (0.4 s) and kept in an LRU cache of `cache` designs."""

    def __init__(self, files, cache=48, designs=None):
        from collections import OrderedDict
        import pickle
        self.files = list(files); self.cache = cache; self._c = OrderedDict()
        self.x, self.lam_full, self.sid, self.params = [], [], [], []
        if designs is not None:                                  # in-memory designs
            for i, d in enumerate(designs):
                self._c[i] = d
            self.cache = max(cache, len(designs))
            src = designs
        else:
            src = [pickle.load(open(f, "rb")) for f in self.files]
        for o in src:
            if isinstance(o, Design):
                self.x.append(o.x); self.lam_full.append(o.lam_full); self.sid.append(o.sid); self.params.append(o.params)
            else:
                self.x.append(design_inputs(o)); self.lam_full.append(float(o.lam_full))
                self.sid.append(o.sid); self.params.append(o.params)
        self.X = torch.as_tensor(np.stack(self.x))
        self.P = [None] * len(self.x); self.hist = [[] for _ in self.x]

    def __len__(self):
        return len(self.x)

    def get(self, i) -> Design:
        import pickle
        if i in self._c:
            self._c.move_to_end(i); return self._c[i]
        d = prepare(pickle.load(open(self.files[i], "rb")))
        self._c[i] = d
        while len(self._c) > self.cache:
            self._c.popitem(last=False)
        return d

    def designs(self, idx):
        return [self.get(i) for i in idx]


def _as_pool(designs):
    return designs if isinstance(designs, DesignPool) else DesignPool([], designs=designs)


def _train_errors(pool):
    e = [pool.hist[i][-1] / pool.lam_full[i] - 1 for i in range(len(pool)) if pool.hist[i] and math.isfinite(pool.lam_full[i])]
    return np.abs(np.asarray(e)) if e else np.array([np.nan])


def train_ritz_ng(net, designs, outer=40, inner=20, batch=8, lr=1e-3, k=3, seed=0, log=print,
                  w_tail=0.25, eval_fn=None, eval_every=5, opt=None):
    """Outer loop: physics targets (one natural-gradient + Ritz step per design); inner loop:
    Adam on the subspace distance to the targets and on the N_cr head."""
    pool = _as_pool(designs)
    gen = np.random.default_rng(seed)
    opt = opt or torch.optim.Adam(net.parameters(), lr=lr)
    R = net.R
    w = torch.full((R,), w_tail, dtype=torch.float64); w[:k] = 1.0
    for it in range(1, outer + 1):
        idx = gen.permutation(len(pool))
        for b0 in range(0, len(idx), batch):
            bi = idx[b0:b0 + batch]
            with torch.no_grad():
                Y0, _ = net(pool.X[bi])
            T, logNt, cur = [], [], []
            for j, i in enumerate(bi):
                d = pool.get(i); d.P = pool.P[i]
                t, mu_t, mu_c = ritz_ng_target(d, Y0[j], R)
                T.append(t); logNt.append(torch.log(1.0 / mu_t[:k].clamp(min=1e-12)).float()); cur.append(mu_c)
            for _ in range(inner):
                Y, logN = net(pool.X[bi])
                loss = sum(subspace_fit_loss(Y[j], T[j], w) for j in range(len(bi))) / len(bi)
                loss = loss + sum(((logN[j] - logNt[j]) ** 2).sum() for j in range(len(bi))) / len(bi)
                opt.zero_grad(); loss.backward(); opt.step()
            with torch.no_grad():
                Y1, _ = net(pool.X[bi])
            for j, i in enumerate(bi):
                pool.P[i] = (Y1[j].double() - Y0[j].double())
                pool.hist[i].append(1.0 / float(cur[j][0]))
        e = _train_errors(pool)
        log(f"outer {it:3d}: training designs, N of the network before the step vs full GBT: "
            f"mean {100 * np.mean(e):.2f} %, max {100 * np.max(e):.2f} %")
        if eval_fn and (it % eval_every == 0 or it == outer):
            eval_fn(it)
    return opt


def train_ritz_ng_full(net, designs, outer=20, epochs=10, batch=16, lr=1e-3, k=3, seed=0, log=print,
                       w_tail=0.25, eval_fn=None, eval_every=2, opt=None):
    """Energy natural gradient with full projection: every outer iteration computes the
    physics step (Ritz line search, `ritz_ng_target`) for ALL designs from the current
    network, then fits the network to all these targets for `epochs` passes (the
    least-squares projection of the function-space step onto the network)."""
    pool = _as_pool(designs)
    gen = np.random.default_rng(seed)
    opt = opt or torch.optim.Adam(net.parameters(), lr=lr)
    R = net.R
    w = torch.full((R,), w_tail, dtype=torch.float64); w[:k] = 1.0
    nd = len(pool)
    for it in range(1, outer + 1):
        with torch.no_grad():
            Y0 = torch.cat([net(pool.X[i:i + 64])[0] for i in range(0, nd, 64)], 0)
        T, logNt = [None] * nd, [None] * nd
        for i in range(nd):
            d = pool.get(i); d.P = pool.P[i]
            t, mu_t, mu_c = ritz_ng_target(d, Y0[i], R)
            T[i] = t / t.norm(dim=1, keepdim=True); logNt[i] = torch.log(1.0 / mu_t[:k].clamp(min=1e-12)).float()
            pool.hist[i].append(1.0 / float(mu_c[0]))
        fit0 = fit1 = 0.0
        for ep in range(epochs):
            perm = gen.permutation(nd); tot = 0.0
            for b0 in range(0, nd, batch):
                bi = perm[b0:b0 + batch]
                Y, logN = net(pool.X[bi])
                lf = sum(subspace_fit_loss(Y[j], T[i], w) for j, i in enumerate(bi)) / len(bi)
                loss = lf + sum(((logN[j] - logNt[i]) ** 2).sum() for j, i in enumerate(bi)) / len(bi)
                opt.zero_grad(); loss.backward(); opt.step()
                tot += float(lf) * len(bi)
            if ep == 0:
                fit0 = tot / nd
            fit1 = tot / nd
        with torch.no_grad():
            for i in range(0, nd, 64):
                Y1 = net(pool.X[i:i + 64])[0]
                for j in range(Y1.shape[0]):
                    pool.P[i + j] = Y1[j].double() - Y0[i + j].double()
        e = _train_errors(pool)
        log(f"outer {it:3d}: training designs, N of the network before the step vs full GBT: mean {100 * np.mean(e):.2f} %, "
            f"max {100 * np.max(e):.2f} %; fit to the physics step: {fit0:.4f} -> {fit1:.4f}")
        if eval_fn and (it % eval_every == 0 or it == outer):
            eval_fn(it)
    return opt


def train_adam(net, designs, steps=1000, batch=8, lr=1e-3, k=3, seed=0, log=print, log_every=50,
               eval_fn=None, eval_every=200, opt=None, cosine=True):
    """Adam on the physics loss itself (Ky Fan Ritz sum + head)."""
    pool = _as_pool(designs)
    gen = np.random.default_rng(seed)
    opt = opt or torch.optim.Adam(net.parameters(), lr=lr)
    for g in opt.param_groups:
        g["lr"] = lr
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=0.05 * lr) if cosine else None
    for it in range(1, steps + 1):
        bi = gen.choice(len(pool), size=min(batch, len(pool)), replace=False)
        Y, logN = net(pool.X[bi])
        opt.zero_grad()
        for j, i in enumerate(bi):
            loss, mu = physics_loss(pool.get(i), Y[j], logN[j], k=k)
            (loss / len(bi)).backward(retain_graph=j < len(bi) - 1)
            pool.hist[i].append(1.0 / float(mu[0]))
        opt.step()
        if sched:
            sched.step()
        if it % log_every == 0:
            e = _train_errors(pool)
            log(f"step {it:5d}: training designs (last visit) vs full GBT: mean {100 * np.mean(e):.2f} %, "
                f"max {100 * np.max(e):.2f} %")
        if eval_fn and it % eval_every == 0:
            eval_fn(it)
    return opt


@torch.no_grad()
def evaluate(net, designs, k=3, idx=None):
    """Per design: Ritz N_cr of the predicted candidates (certified upper bound in the Ritz
    space), N_cr of the head (instant), full GBT."""
    pool = _as_pool(designs)
    idx = range(len(pool)) if idx is None else idx
    out = []
    for i in idx:
        d = pool.get(i)
        Y, logN = net(pool.X[i:i + 1])
        mu, _ = ritz(d, Y[0])
        out.append(dict(sid=d.sid, bc=d.bc, params=d.params, N_ritz=float(1.0 / mu[0]) if mu[0] > 0 else float("inf"),
                        N_head=float(torch.exp(logN[0, 0])), N_full=d.lam_full))
    return out
