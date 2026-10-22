"""The wrappers must reproduce the validated V32 results (slow: real solves).
Run with:  pytest -m solver"""
import pytest

from pinnbuck.designspace import validated_member_sample, validated_omega_sample
from pinnbuck.engine import solver_mode_check
from pinnbuck.oracle import solve

pytestmark = pytest.mark.solver


@pytest.mark.parametrize("kind,mesh,ref", [("gbt", "standard", 32.2603), ("strip", "standard", 32.1153)])
def test_omega_unsym_periodic(kind, mesh, ref):
    s = validated_omega_sample("Unsym", 150.0, 1.5, "SS", "periodic")
    r = solve(s, kind, mesh)
    assert r["ok"], r.get("error")
    assert abs(r["value"] / ref - 1) < 2e-5


@pytest.mark.parametrize("sec,lay,L,bc,kind,ref", [
    ("c", "Sym1", 500.0, "SS", "gbt", 2674.849),
    ("c", "Sym1", 500.0, "SS", "strip", 2650.08),
    ("l", "Unsym", 200.0, "CC", "gbt", 458.751),
])
def test_members(sec, lay, L, bc, kind, ref):
    r = solve(validated_member_sample(sec, lay, L, bc), kind, "standard")
    assert r["ok"], r.get("error")
    assert abs(r["value"] / ref - 1) < 5e-5


@pytest.mark.parametrize("n_conv", [80, None])
def test_engine_energy_equals_solver(n_conv):
    s = validated_omega_sample("Unsym", 150.0, 1.5, "SS", "periodic")
    r = solver_mode_check(s, n_conv)
    assert abs(r["rel"]) < 1e-9
