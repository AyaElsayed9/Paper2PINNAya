# Abaqus runs for Paper 2: extended design space (4 Nov 2026)

**Why.** Paper 2 covers the same design space as Paper 1: b_s from 124.25 mm and β up to 3.5, with
both side conditions. GBT has been compared with Abaqus there only for symmetry with the cross-ply
lay-up (Paper 1). These cells fill the gaps:

| set | side edges | ends | cells | why |
|---|---|---|---|---|
| A | periodic | SS | 1 control + 8 | b_s < 150 mm and β > 3.0 had no Abaqus check for periodic |
| B | symmetry (as Paper 1) | SS | 1 control + 6 | tight pitch with other lay-ups than cross-ply |
| C (optional) | periodic | CC | 4 | for the clamped network later |

The positions come from a Sobol sequence with a fixed seed (`make_cases.py`), not from any GBT or network
result. The lay-ups are the named ones of the Abaqus script (Cross-ply, Sym1, Sym2, Unsym). The first cell
of A and B (b_s 150, β 1.5, Sym1) was already run for the periodic-GBT paper. It checks that the run
reproduces the earlier model.

**Blind comparison.** The GBT values (`gbt_predictions.csv`) were computed and committed to git **before**
the Abaqus runs. The Abaqus numbers are used only for the comparison, never to change anything.

## How to run (Windows, PowerShell)
1. Put `run_extension.ps1` in the same folder as `Omegastiffener_PBC_v2.py` (the same script as for the
   periodic paper; a copy is in `manuscript/abaqus/`).
2. Open PowerShell in that folder and run
   `powershell -ExecutionPolicy Bypass -File .\run_extension.ps1`
3. The results go to `Downloads\abaqus_paper2_extension\A_periodic_SS`, `...\B_symmetry_SS` and
   `...\C_periodic_CC`. Send back the CSV files of every folder:
   `omega_validation_eigs*.csv` (N_cr and the half-wave numbers) and `omega_validation_times*.csv`.
   The `fe_modes` folders are not needed.

Set C can be skipped: delete its four lines in the `.ps1`, or stop Abaqus after set B.
