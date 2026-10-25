"""Fast checks of the network utilities (no solver)."""
import numpy as np
import torch

from pinnbuck.net import ModeCodec, subspace_loss


def test_codec_roundtrip():
    g = np.random.default_rng(0)
    nm, J = 7, 3; nf = 2 * J + 1
    A = g.standard_normal((nm, nm)); Kb = A @ A.T + nm * np.eye(nm)
    ev, V = np.linalg.eigh(Kb)
    Wm, Wp = (V / np.sqrt(ev)) @ V.T, (V * np.sqrt(ev)) @ V.T
    recs = [dict(M=g.standard_normal((3, nm, nf)), L=300.0 + 50 * i, J=J) for i in range(30)]
    codec = ModeCodec.fit(recs, Wm, Wp, n_pc=nm * nf)            # complete PCA: exact
    M = recs[4]["M"]
    assert np.allclose(codec.decode(codec.encode(M, recs[4]["L"]), recs[4]["L"], nm), M, atol=1e-10)


def test_subspace_loss_invariances():
    g = torch.Generator().manual_seed(0)
    C = torch.randn(2, 3, 20, generator=g, dtype=torch.float64)
    R = torch.randn(3, 3, generator=g, dtype=torch.float64)       # any invertible mixing
    assert torch.allclose(subspace_loss(C, R @ C), torch.zeros(2, dtype=torch.float64), atol=1e-12)
    assert float(subspace_loss(C, torch.randn(2, 3, 20, generator=g, dtype=torch.float64)).min()) > 0.1
