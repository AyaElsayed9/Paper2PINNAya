"""Access to the frozen, validated solver (solver/aya_gbt_v32.py).

The solver is imported once as a module and never modified here.  Its
module-level switches (side condition, solver choice, mesh) are set through
`configure()`, which restores a known state before every call so that results
never depend on the order of earlier calls.
"""
from __future__ import annotations

import hashlib
import importlib.util
import sys
from contextlib import contextmanager
from pathlib import Path

SOLVER_PATH = Path(__file__).resolve().parents[2] / "solver" / "aya_gbt_v32.py"
_MODULE = None


def solver():
    """The solver module (imported on first use)."""
    global _MODULE
    if _MODULE is None:
        argv = sys.argv
        sys.argv = [str(SOLVER_PATH)]          # the solver reads flags from argv at import
        try:
            spec = importlib.util.spec_from_file_location("aya_gbt_v32", SOLVER_PATH)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
        finally:
            sys.argv = argv
        _MODULE = mod
    return _MODULE


def solver_hash() -> str:
    """SHA-256 of the solver file: stored with every result for provenance."""
    return hashlib.sha256(SOLVER_PATH.read_bytes()).hexdigest()[:16]


_DEFAULTS = dict(OMEGA_SIDE_BC="symmetry", OMEGA_PERIODIC_MACRO=False, OMEGA_MESH="standard",
                 OMEGA_MESH_FINE=False, OMEGA_GBT_MODES=None, OMEGA_GBT_NINT=None,
                 OMEGA_CC_END="edge-load", GOVERNING_SOLVER="gbt", GBT_FAMILIES="CSTMR",
                 MESH_STRIP=None, MESH_SPAN=None)


@contextmanager
def configure(**switches):
    """Set solver switches for the duration of a block, then restore defaults.

    Example: ``with configure(OMEGA_SIDE_BC="periodic", OMEGA_MESH="fast"): ...``
    """
    m = solver()
    unknown = set(switches) - set(_DEFAULTS)
    if unknown:
        raise KeyError(f"unknown solver switches: {sorted(unknown)}")
    state = {**_DEFAULTS, **switches}
    old = {k: getattr(m, k) for k in _DEFAULTS}
    for k, v in state.items():
        setattr(m, k, v)
    try:
        yield m
    finally:
        for k, v in old.items():
            setattr(m, k, v)
