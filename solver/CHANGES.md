# Changes to the frozen solver (aya_gbt_v32.py)

The validated V32 solver is edited only to fix bugs that do not change any validated
result. Every fix is listed here. After each fix the regression tests pass
(`pytest -m solver`, which reproduces the validated N_cr values).

| date | function | bug | effect on validated results |
|---|---|---|---|
| 2026-10 | `_sfs_breaks` | Rounding could put the last span breakpoint past L. Spline knots were then out of order, and some lengths crashed with "Knots must be in a non-decreasing order". Fix: `b[0], b[-1] = 0, L`. | none (identical values) |
| 2026-10 | `truncated_section` | The reduced copy of a section inherited two caches of its parent: the Kg sparsity pattern `_kg_pattern` and the modes `_hm`. A reduced solve after a full solve on the same section then failed with an index error. Fix: drop both caches in the copy. | none (validated runs never truncate a solved section) |
| 2026-10 | `gbt_fsm_solve` (members) | Not a bug: new optional outputs for the pinnbuck operators. `pre_only=True` returns the first-order pre-buckling state, `pre=` imposes one, and `return_mode` also returns all computed modes (`coefs`). | none (defaults unchanged) |
