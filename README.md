# Paper2PINNAya — physics-informed learning of buckling of thin-walled composite sections

Paper 2 of the PINN buckling project (target journal: *Composite Structures*).
It extends Paper 1 (deep energy method with active learning for the omega-stiffened panel)
to the cross-sections of the validated GBT code and to imperfection sensitivity.

| role | model | where |
|---|---|---|
| physics engine (loss of the PINN) | shear-deformable GBT in a shared POD mode basis, side conditions built into the modes | `src/pinnbuck/basis.py`, `operators.py` |
| fast PINN | one network per family: design → N_cr and buckling mode (GBT deep-Ritz) | `src/pinnbuck/pinn.py` |
| oracle (training data, active learning) | shell-type spline finite strip | `src/pinnbuck/oracle.py` |
| blind validation | Abaqus (never used for training) | user-supplied results |

The physics comes from the frozen, validated solver `solver/aya_gbt_v32.py`
(AYA-GBT-V32; omega: 0.58 % mean error vs 60 Abaqus cases with symmetry edges, 0.57 % vs 30 periodic cases;
members: 144 Abaqus cases). It is imported as a library; it is edited only for bug fixes that leave every validated value unchanged
(listed in `solver/CHANGES.md`). Its SHA-256 is stored with every result.

## Layout
```
solver/aya_gbt_v32.py     frozen validated solver (do not edit)
src/pinnbuck/
  solver.py               imports the solver, sets its switches safely
  designspace.py          families, validated ranges, lay-ups, lamination parameters, Sobol sampling
  oracle.py               strip / GBT reference solves, resumable parallel batches
  engine.py               GBT operators + pre-buckling state (per-panel hierarchical basis; checks)
  basis.py                shared GBT mode basis of a family (POD of training buckling modes)
  operators.py            GBT operators of a panel in the shared basis, end-condition classes, Ritz checks
  pinn.py                 parametric GBT deep-Ritz network, physics loss, training
configs/                  one YAML file per experiment
scripts/                  command-line entry points (long runs)
tests/                    unit tests and regression against the validated values
notebooks/                exploration and figures only (import from src/, no core logic)
docs/                     design space, roadmap
data/                     generated datasets (JSON lines)
```

## Install (Windows laptop or Linux)
```
python -m venv .venv
.venv\Scripts\activate            (Linux: source .venv/bin/activate)
pip install -e .[fast,test]
```
`pypardiso` is optional (same results, 2–3x faster solves).

## Check the installation
```
pytest -m "not solver"            fast tests (seconds)
pytest -m solver                  regression against the validated values (a few minutes)
```

## Generate data
```
python scripts/generate_oracle_data.py configs/omega_pilot.yaml --workers 4
```
Runs are resumable: finished samples are skipped when the command is repeated.
