# Roadmap of Paper 2

Working title: *GBT-informed deep energy learning of buckling and imperfection sensitivity of
thin-walled composite sections* (Composite Structures).

## Contributions
1. A deep-energy PINN whose physics is shear-deformable GBT: the network predicts a small set of GBT amplitude
   functions instead of displacement fields, side conditions (symmetry, periodic) hold exactly by construction,
   and one framework covers the omega panel and the member sections.
2. Lay-up generality through lamination parameters, with a spline-strip oracle replacing Abaqus inside the
   active-learning loop (Abaqus kept for blind validation).
3. Imperfection sensitivity of the critical modes by Koiter's initial post-buckling analysis in the same GBT
   basis, giving knock-down factors for mode-shaped imperfections with uncertainty.

## Phases
| phase | content | status |
|---|---|---|
| 0 | frozen solver, package, design space, engine and oracle wrappers, tests | done (this commit) |
| 1 | oracle datasets: omega (pilot → ~500), then members (~150 per section) | pilot running |
| 2 | PINN per geometry (Paper 1 workflow): GBT energy, reduced basis, warm start from the neighbour, ensemble + strip oracle | next |
| 3 | parametric surrogate on the PINN dataset (inputs: geometry + lamination parameters) and its uncertainty | |
| 4 | Koiter imperfection sensitivity in GBT (see below) | |
| 5 | blind Abaqus validation (new random points), figures, manuscript | |

## Phase 2 — PINN per geometry
* Network: (ζ = z/L) → amplitudes of the n reduced GBT modes (n_conv = 80 → about 160 modes; 1–3.4 % from the full
  model; full-model eigenvalue stored as reference).
* End conditions by output transforms; normalisation as in Paper 1 (RMS amplitude = 1).
* Loss: maximise the Rayleigh quotient μ = U_G / U (U_G is indefinite, so U / U_G is not used) + normalisation;
  data term (smooth L1 on λ / λ_oracle) only when the oracle is called.
* Active learning: suspicion check against the neighbour, deep ensemble, CoV > 5 % → strip oracle.
* Check: the engine reproduces the solver eigenvalue to 1e-12 when given the solver's mode (`tests/`).

## Phase 4 — Koiter imperfection sensitivity in GBT
For a critical mode u₁ at λ_c (normalised), the post-buckling path and imperfection sensitivity follow from
the cubic and quartic terms of the potential energy (membrane strains with the quadratic terms
½(v,z² + w,z²) etc. in the GBT kinematics):
* second-order field u₂: one linear solve with the stiffness K and a load built from the cubic term,
* coefficients a (cubic) and b (quartic) of λ/λ_c = 1 + a ξ + b ξ² + …,
* imperfection sensitivity (mode-shaped imperfection ξ̄):
  a ≠ 0 (asymmetric): (1 − λ_s/λ_c)² = 4 |a ξ̄| λ_s/λ_c for the imperfection sign that lowers the load;
  a = 0, b < 0 (unstable symmetric): (1 − λ_s/λ_c)^{3/2} = (3√3/2) √(−b) |ξ̄| λ_s/λ_c (Koiter/Hutchinson).
* Close eigenvalues (found for several lay-ups) need the multi-mode (interactive) form.
Validation: Abaqus Riks with mode-shaped imperfections on a subset (the workflow of Paper 1).
Uncertainty: imperfection amplitude ξ̄ and material/ply-angle scatter propagated by Monte Carlo through the
surrogate → distribution of the limit load.

## Data the user supplies
* Blind Abaqus results at the random points listed by `scripts/` (phase 5), Riks runs for phase 4.
