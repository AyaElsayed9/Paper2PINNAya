"""Geometrically non-linear GBT path of the omega panel (20 Nov 2026), reduced-basis form.

Same energy as src/pinnbuck/koiter.py (the solver's quadratic membrane strains), solved along the whole path:
    Pi(a) = 1/2 a'K a + int N(a).e(a) + 1/2 int e(a)' A e(a) - lam f'a,
    e(a) = eps_Q(a) + 2 B(a, a_bar)          (a_bar = initial geometric imperfection, stress free)
    N(a) = Cp[:, :3, :] eps_L(a)
The displacement is sought in a basis Phi = [pre-buckling state, lowest buckling modes, Koiter second-order fields
w_ij] (Koiter-Newton type reduced basis, cf. Liang, Abdalla & Gurdal 2013); the reduced equations are solved
exactly by Newton with pseudo-arc-length continuation.  The residual of the FULL equations, ||C'(dPi/da)||, is
reported along the path, so the basis error is visible; the basis can be enriched with K^-1 (full residual).
Fields of the basis vectors are linear in the coefficients and are precomputed once.
"""
from __future__ import annotations

import numpy as np
import torch

from . import koiter as KO

DT = torch.float64


class Reduced:
    def __init__(self, S, Phi, abar=None):
        """S from koiter.assemble; Phi (ndof, n) full-coordinate basis (admissible: columns = Cm @ q)."""
        self.S = S; self.F = F = KO.Fields(S)
        self.Phi = Phi; n = Phi.shape[1]
        W = F.W(); self.W = W
        Ns, gs = [], []
        with torch.no_grad():
            for k in range(n):
                A = F.amp(torch.as_tensor(Phi[:, k], dtype=DT))
                Ns.append(torch.einsum("pca,xpa->xpc", F.Cp[:, :3, :], F.lin(A)))
                gs.append(F.grad(A))
        self.N = torch.stack(Ns, 0)                  # (n, nx, np, 3)
        self.g = torch.stack(gs, 0)                  # (n, nx, np, nc, 2)
        self.Am = F.Cp[:, :3, :3]
        Kphi = np.stack([S["Kmul"](Phi[:, k]) for k in range(n)], 1)
        self.K = torch.as_tensor(Phi.T @ Kphi, dtype=DT); self.K = 0.5 * (self.K + self.K.T)
        self.f = torch.as_tensor(Phi.T @ S["f"], dtype=DT)
        if abar is None:
            self.gbar = torch.zeros_like(self.g[0])
        else:
            with torch.no_grad():
                self.gbar = F.grad(F.amp(torch.as_tensor(abar, dtype=DT)))

    def energy(self, y):
        """internal energy (without the load term) of a = Phi y"""
        N = torch.tensordot(y, self.N, dims=1); g = torch.tensordot(y, self.g, dims=1)
        e = KO.Fields.quad(g) + 2.0 * KO.Fields.quad(g, self.gbar)
        lin = 0.5 * y @ (self.K @ y)
        c3 = (self.W * (N * e).sum(-1)).sum()
        c4 = 0.5 * (self.W * torch.einsum("xpi,pij,xpj->xp", e, self.Am, e)).sum()
        return lin + c3 + c4

    def grad(self, y):
        yt = torch.as_tensor(y, dtype=DT).requires_grad_(True)
        (g,) = torch.autograd.grad(self.energy(yt), yt)
        return g.numpy()

    def grad_hess(self, y):
        yt = torch.as_tensor(y, dtype=DT)
        g = torch.autograd.functional.jacobian(self.energy, yt, create_graph=False)
        H = torch.autograd.functional.hessian(self.energy, yt)
        return g.numpy(), H.numpy()

    def full_residual(self, y, lam):
        """|| C'(dPi/da - lam f) || / || lam C'f || of the full model at a = Phi y"""
        S, F = self.S, self.F; Cm = S["Cm"]
        a = torch.as_tensor(self.Phi @ y, dtype=DT).requires_grad_(True)
        A = F.amp(a); N = torch.einsum("pca,xpa->xpc", F.Cp[:, :3, :], F.lin(A)); g = F.grad(A)
        e = KO.Fields.quad(g) + 2.0 * KO.Fields.quad(g, self.gbar)
        nl = (self.W * (N * e).sum(-1)).sum() + 0.5 * (self.W * torch.einsum("xpi,pij,xpj->xp", e, self.Am, e)).sum()
        (gn,) = torch.autograd.grad(nl, a)
        an = a.detach().numpy()
        R = Cm.T @ (S["Kmul"](an) + gn.numpy() - lam * S["f"])
        # error estimate in the energy norm: delta = K^-1 R relative to the displacement itself
        dl = S["solve"](R); q = Cm.T @ an
        en_a = np.sqrt(max(an @ S["Kmul"](an), 1e-300)); en_d = np.sqrt(max(dl @ (S["Kr"] @ dl), 0.0))
        return float(en_d / en_a), Cm @ dl


def basis(S, lam, V, M_modes=8, M_pairs=5):
    """[pre-buckling a0, lowest M_modes modes (max nodal displacement 1 mm), second-order fields w_ij of the
    lowest M_pairs modes] -> orthonormalised in the K metric (full coordinates)"""
    km = KO.koiter_multi(S, lam, V, M_pairs)
    cols = [S["a0"] / np.abs(S["a0"]).max()]
    F = KO.Fields(S)
    for i in range(min(M_modes, V.shape[1])):
        a = S["Cm"] @ V[:, i]; cols.append(a / F.maxdisp(torch.as_tensor(a)))
    for (i, j), w in km["W"].items():
        if i <= j:
            cols.append(w / max(np.abs(w).max(), 1e-300))
    P = np.stack(cols, 1)
    # K-orthonormalise (Gram-Schmidt twice), drop dependent columns
    out = []
    for k in range(P.shape[1]):
        v = P[:, k].copy()
        for _ in range(2):
            for u in out:
                v -= (u @ S["Kmul"](v)) * u
        nv = np.sqrt(max(v @ S["Kmul"](v), 0.0))
        if nv > 1e-8 * np.sqrt(max(P[:, k] @ S["Kmul"](P[:, k]), 1e-300)):
            out.append(v / nv)
    return np.stack(out, 1), km


def enrich(Phi, S, extra):
    """append correction vectors (full coords) to a K-orthonormal basis"""
    out = [Phi[:, k] for k in range(Phi.shape[1])]
    for v in extra:
        v = v.copy()
        for _ in range(2):
            for u in out:
                v -= (u @ S["Kmul"](v)) * u
        nv = np.sqrt(max(v @ S["Kmul"](v), 0.0))
        if nv > 1e-10:
            out.append(v / nv)
    return np.stack(out, 1)


def solve_path(S, Phi, abar, lam_max, tol_err=0.01, rounds=4, ds=0.05, ds_max=0.5, check_every=8, log=print):
    """path with basis enrichment: rerun until the energy-norm error estimate along the path is below tol_err"""
    for rd in range(rounds + 1):
        R = Reduced(S, Phi, abar=abar)
        out = path(R, lam_max, ds=ds, check_every=check_every, ds_max=ds_max)
        errs = out["res"]; emax = max(e for _, e in errs) if errs else 0.0
        log(f"round {rd}: basis {Phi.shape[1]}, steps {len(out['lam'])}, lam_end {out['lam'][-1]:.3f}, max error {emax:.4f}")
        out["basis"] = Phi; out["round"] = rd
        if emax <= tol_err or rd == rounds:
            return out, R
        bad = [i for i, e in errs if e > tol_err]
        pick = bad[:: max(1, len(bad) // 6)][:6]
        extra = [R.full_residual(out["y"][i], out["lam"][i])[1] for i in pick]
        Phi = enrich(Phi, S, extra)


def path(R: Reduced, lam_max, ds=0.05, nmax=400, tol=1e-9, check_every=10, verbose=False, ds_max=2.0):
    """pseudo-arc-length continuation of  dE/dy - lam f = 0  from (0, 0).  Returns dict of arrays."""
    n = R.Phi.shape[1]; f = R.f.numpy()
    sc_l = 1.0 / max(lam_max, 1e-12)              # lam scaled so that lam_max ~ 1
    y = np.zeros(n); lam = 0.0
    # initial tangent: linear response direction
    y1 = np.linalg.solve(R.K.numpy(), f); t = np.concatenate([y1, [1.0 / sc_l]]); t /= np.linalg.norm(np.concatenate([y1, [1.0]]))
    t = np.concatenate([y1, [1.0]]); t /= np.linalg.norm(t)
    Ys, Ls, Res = [y.copy()], [0.0], []
    for step in range(nmax):
        yp = y + ds * t[:-1]; lp = lam + ds * t[-1] / sc_l
        ok = False; H = None
        for it in range(30):
            # modified Newton (20 Nov): the Hessian is rebuilt only at the first iteration and when convergence slows
            if H is None or it in (6, 12, 18, 24):
                g, H = R.grad_hess(yp)
            else:
                g = R.grad(yp)
            r = g - lp * f
            c = t[:-1] @ (yp - y) + t[-1] * (lp - lam) * sc_l - ds
            J = np.zeros((n + 1, n + 1)); J[:n, :n] = H; J[:n, n] = -f / sc_l; J[n, :n] = t[:-1]; J[n, n] = t[-1]
            d = np.linalg.solve(J, -np.concatenate([r, [c]]))
            yp = yp + d[:-1]; lp = lp + d[-1] / sc_l
            if np.linalg.norm(d) < tol * max(1.0, np.linalg.norm(yp)):
                ok = True; nit = it + 1; break
        if not ok:
            ds *= 0.5
            if ds < 1e-5:
                break
            continue
        tn = np.concatenate([yp - y, [(lp - lam) * sc_l]]); tn /= np.linalg.norm(tn)
        t = tn; y = yp; lam = lp
        Ys.append(y.copy()); Ls.append(lam)
        if nit <= 6:
            ds = min(ds * 1.4, ds_max)            # adaptive arc length
        elif nit > 14:
            ds *= 0.7
        if check_every and len(Ls) % check_every == 0:
            Res.append((len(Ls) - 1, R.full_residual(y, lam)[0]))
        if verbose:
            print(step, round(lam, 4), flush=True)
        if lam > lam_max or lam < -0.05 * lam_max:
            break
    return dict(y=np.array(Ys), lam=np.array(Ls), res=Res)
