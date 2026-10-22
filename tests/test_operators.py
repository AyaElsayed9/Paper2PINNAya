"""Shared-basis operators and PINN energies against the solver (pytest -m solver)."""
import numpy as np
import pytest
import torch

from pinnbuck.basis import NINT_REF, build_shared_basis, hx_of
from pinnbuck.designspace import validated_omega_sample
from pinnbuck.operators import build_operator, class_bases, ritz_lambda, ritz_matrices
from pinnbuck.oracle import omega_length, omega_panel
from pinnbuck.pinn import TorchPanel, energies, ritz_forms, ritz_mu
from pinnbuck.solver import configure

pytestmark = pytest.mark.solver


@pytest.fixture(scope="module", params=["SS", "CC"])
def op(request):
    sb = build_shared_basis("symmetry", 20)
    smp = validated_omega_sample("Unsym", 150.0, 0.75, request.param, "symmetry")
    return smp, build_operator(smp, sb, lam_full=float("nan"))


def test_operators_reproduce_solver(op):
    """Ritz on the solver's own spline with the compressed operators = solver eigenvalue."""
    smp, o = op
    with configure(OMEGA_SIDE_BC="symmetry", OMEGA_GBT_NINT=list(NINT_REF)) as m:
        sec = m.omega_gbt_section(omega_panel(smp)); L = omega_length(smp); hx = hx_of(sec, L)
        lb = m._LongSpline(L, m._sfs_breaks(L, hx, hx, L))
        Bz = np.stack([lb.B[d].toarray() for d in range(3)])
    assert abs(ritz_lambda(o, Bz) / o.lam_shared - 1) < 1e-9


def test_pinn_energies_match_ritz_forms(op):
    """Energies of a separated mode in pinn.py = quadratic forms of the Ritz matrices."""
    _, o = op
    J, R = 6, 3
    tp = TorchPanel(o, J, k0=1e12)                     # no Fourier scaling, unit mode scale
    g = np.random.default_rng(0)
    U = g.standard_normal((o.nm, R)); C = g.standard_normal((R, 2 * J + 1))
    Ue, UG = energies(tp, torch.as_tensor(U), tp.F0 @ torch.as_tensor(C).T)
    Kr, Gr = ritz_matrices(o, class_bases(o, "fourier", J))
    a = np.concatenate([((c.Pi.T @ U) @ C).ravel() for c in o.classes])
    assert abs(float(Ue) / (0.5 * a @ Kr @ a) - 1) < 1e-10
    assert abs(float(UG) / (0.5 * a @ Gr @ a) - 1) < 1e-8


def test_ritz_forms_consistent(op):
    """Candidate-mode matrices: diagonal = 2 x energies of each candidate; Ritz value
    not above the best single candidate."""
    _, o = op
    J, R = 6, 3
    tp = TorchPanel(o, J, k0=1e12)
    g = np.random.default_rng(1)
    U = torch.as_tensor(g.standard_normal((o.nm, R))); C = torch.as_tensor(g.standard_normal((R, 2 * J + 1)))
    KR, GR = ritz_forms(tp, U, tp.F0 @ C.T)
    for r in range(R):
        Ue, UG = energies(tp, U[:, r:r + 1], tp.F0 @ C[r:r + 1].T)
        assert abs(float(KR[r, r]) / (2 * float(Ue)) - 1) < 1e-10
        assert abs(float(GR[r, r]) / (2 * float(UG)) - 1) < 1e-8
    mu, x = ritz_mu(KR, GR)
    assert float(mu) >= max(float(GR[r, r] / KR[r, r]) for r in range(R)) - 1e-12
    assert abs(float(x @ GR @ x) / float(x @ KR @ x) / float(mu) - 1) < 1e-8
