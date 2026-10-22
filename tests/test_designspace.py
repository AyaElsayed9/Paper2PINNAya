import numpy as np
import pytest

from pinnbuck.designspace import FAMILIES, lamination_parameters, sample_design, sample_layup


def test_lamination_parameters_unidirectional():
    xi = lamination_parameters([0.0] * 8, 0.125)
    assert np.allclose(xi[[0, 2]], 1.0) and np.allclose(xi[[1, 3]], 0.0)   # A
    assert np.allclose(xi[4:8], 0.0)                                         # B (symmetric)
    assert np.allclose(xi[[8, 10]], 1.0)                                     # D


def test_lamination_parameters_symmetric_has_no_coupling():
    xi = lamination_parameters([45, -45, -45, 45, 45, -45, -45, 45], 0.125)
    assert np.allclose(xi[4:8], 0.0)
    assert abs(xi[0]) < 1e-12 and np.isclose(xi[2], -1.0)                    # [+-45]: cos2t = 0, cos4t = -1


def test_lamination_parameters_bounds():
    rng = np.random.default_rng(1)
    for _ in range(200):
        xi = lamination_parameters(sample_layup(8, rng, symmetric=False), 0.125)
        assert np.all(np.abs(xi) <= 1 + 1e-12)


def test_sampling_is_deterministic_and_in_range():
    a = sample_design("omega", 32, seed=7); b = sample_design("omega", 32, seed=7)
    assert [s.sid for s in a] == [s.sid for s in b]
    lo, hi = FAMILIES["omega"]["continuous"]["bs"]
    assert all(lo <= s.params["bs"] <= hi for s in a)
    assert len({s.sid for s in a}) == 32
    assert all(len(s.features()) == len(a[0].features()) for s in a)


def test_symmetric_layups_are_symmetric():
    rng = np.random.default_rng(3)
    for n in (4, 7, 8):
        lay = sample_layup(n, rng, symmetric=True)
        assert lay == lay[::-1]
