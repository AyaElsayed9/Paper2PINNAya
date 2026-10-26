"""Energy space of the physics-only PINN (ritzspace.py) against direct z-integrals of the
operators (pytest -m solver)."""
import numpy as np
import pytest
import torch

from pinnbuck.basis import build_shared_basis
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.operators import build_operator
from pinnbuck.ritzspace import EnergyMetric, FlatMetric, RitzSpace, ritz_forms
from pinnbuck.spectral import dot, lobpcg

pytestmark = pytest.mark.solver


@pytest.fixture(scope="module", params=["SS", "CC"])
def rs(request):
    sb = build_shared_basis("symmetry", 20)
    op = build_operator(validated_omega_sample("Cross-ply", 200.0, 1.5, request.param, "symmetry"), sb,
                        lam_full=float("nan"), with_shared_reference=True)
    return RitzSpace(op, raw=True)


def _direct(rs, A):
    op = rs.op; nm = op.nm
    Q = rs.field(A); Qc = Q.permute(0, 2, 1, 3).reshape(Q.shape[0], Q.shape[2], 3 * nm)
    U = 0.5 * torch.einsum("z,rzi,ij,rzj->r", rs.wz, Qc, torch.as_tensor(op.K), Qc)
    g = torch.einsum("pcai,rzi->rzpca", torch.as_tensor(op.G.astype(float)), Qc)
    P = torch.as_tensor(op.P.astype(float))
    dens = (P[0][None, :, :, None] * g[..., 0] ** 2 + P[1][None, :, :, None] * g[..., 1] ** 2
            + 2 * P[2][None, :, :, None] * g[..., 0] * g[..., 1]).sum((-1, -2))
    return U, 0.5 * torch.einsum("z,rz->r", rs.wz, dens)


def test_kronecker_K_and_compressed_G_equal_direct_integrals(rs):
    A = rs.random(3, torch.Generator().manual_seed(0))
    U, UG = _direct(rs, A)
    assert torch.allclose(0.5 * torch.diagonal(dot(A, rs.apply_K(A))), U, rtol=1e-10)
    assert torch.allclose(0.5 * torch.diagonal(dot(A, rs.apply_G(A))), UG, rtol=1e-9)


def test_ritz_forms_gradient(rs):
    X = torch.randn(3, rs.n, dtype=torch.float64, requires_grad=True)
    W1, W2 = torch.randn(3, 3, dtype=torch.float64), torch.randn(3, 3, dtype=torch.float64)
    f = lambda X: (W1 * ritz_forms(rs, X)[0]).sum() + (W2 * ritz_forms(rs, X)[1]).sum()
    g, = torch.autograd.grad(f(X), X)
    v = torch.randn_like(X); e = 1e-6
    with torch.no_grad():
        fd = (f(X + e * v) - f(X - e * v)) / (2 * e)
    assert abs(float(fd) / float((g * v).sum()) - 1) < 1e-6


def test_metric_whitening_roundtrip_and_eigen_iteration(rs):
    em = EnergyMetric(rs, groups="labels", ridge=1e-8); fm = FlatMetric(rs, em)
    Y = torch.randn(2, rs.n, dtype=torch.float64)
    assert torch.allclose(fm.whiten(fm.color(Y)), Y, rtol=1e-8, atol=1e-8)
    lam, _, hist = lobpcg(rs, R=6, iters=40, precond=em.solve, tol=1e-10, seed=1)
    # the dictionary Ritz space reproduces the solver in the shared basis
    assert abs(float(lam[0]) / rs.op.lam_shared - 1) < 3e-3
