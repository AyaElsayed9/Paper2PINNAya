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
| 1 | oracle datasets: omega (pilot → ~500), then members (~150 per section) | pilot done (16 omega samples) |
| 2 | shared GBT basis (done), operators (done), fast parametric PINN (open, see below) | in progress |
| 3 | parametric surrogate on the PINN dataset (inputs: geometry + lamination parameters) and its uncertainty | |
| 4 | Koiter imperfection sensitivity in GBT (see below) | |
| 5 | blind Abaqus validation (new random points), figures, manuscript | |

## Phase 2: shared GBT basis and the fast PINN (in progress)

### 2a. Shared mode basis: done
A network for the whole family needs outputs that mean the same mode on every panel.
* Every panel uses the same node layout (`basis.NINT_REF`), so the nodal fields are comparable.
* Reference-panel modes failed. With the first 80 hierarchical modes of one quasi-isotropic panel, the error on random lay-ups was +20 to +130 %, because the mode order depends on the laminate.
* **POD of buckling modes** (`basis.pod_basis`):
  * Snapshots: the 3 lowest full-GBT modes of 64 training panels, Sobol seed 2026, no Abaqus.
  * Block structure: conventional (energy-weighted), shear hats and Mindlin, with the Mindlin block as nodal fields.
  * Chosen size: 120 conventional + 1 axial + 60 shear + 30 + 30 Mindlin, so 241 modes in total.
* **Error vs full GBT** on unseen panels (Sobol seed 77 plus the Paper 1 validated lay-ups, `scripts/exp_shared_basis.py`):

  | side | test panels | mean | max |
  |---|---|---|---|
  | symmetry | 19 | 0.02 % | 0.07 % |
  | periodic | 20 | 0.01 % | 0.03 % |

  The panel's own 80-mode basis reaches 1.6 % mean, 2.9 % max.
* Smaller bases fail on some panels: 80 conventional modes reach 2 % on one panel, and the Mindlin modes are needed.

### 2b. Operators: done
* `operators.build_operator` gives K, the pointwise G and the pre-buckling state of a panel in the shared basis, plus the end-condition classes.
* Check: the Ritz solve on the solver's spline reproduces the solver eigenvalue to 1e-13 (`tests/test_operators.py`).
* A Fourier basis along z (J = 16 to 24, cos and sin) adds at most 0.1 %.

### 2c. Parametric network: open problem
Built in `pinn.py`:
* R candidate modes, each a cross-section vector × a longitudinal function (Fourier, or a wave sin/cos(κz) × envelope), with exact end conditions.
* Rayleigh–Ritz over the candidates (upper bound) plus a direct N_cr head, so prediction is about 0.2 ms per design.

Findings:
* **One panel:** physics-only training converges to 0.2 % of the GBT value in about 1 min.
* **Many panels (8 to 32), physics only:** convergence is very slow. After 2000 Adam steps the training panels are still 15–140 % above the GBT value.
* **Why:** the lowest modes are close together (e.g. 16.02 / 16.08 / 16.74) and spread over many wavelengths. The Rayleigh landscape is therefore flat, and gradient training of high-dimensional mode outputs is slow.
* **Teacher stage** (K-cosine to the full-GBT modes): also slow across panels.
* **Kantorovich form** (one wavenumber + envelope, cross-section solved exactly): 1–13 % error, so the separable form is too restrictive.
* **Pre-stress rank:** an SVD rank of 5 changes N_cr by up to ±0.5 %.

Decision needed (options in the message of this commit): hybrid physics-informed surrogate (GBT labels + energy consistency), physics fine-tuning per design, or further research on physics-only amortised training.

### 2d. Research log: label-free training of one network for the family (Oct 2026)
* **Physics of the panels.**
  * The load introduction (stringer foot eccentricity) raises the pre-buckling stress near the loaded ends to up to 5× the mid-length value.
  * Many critical modes are localised at the ends; on one panel, 78 % of the mode lies in the two end fifths.
  * The lowest modes are close together (e.g. 16.02 / 16.08 / 16.74).
* **Representations along z** (exact cross-section, Ritz, vs full GBT):

  | representation | error |
  |---|---|
  | single wave (Kantorovich) | 1–13 % |
  | wavenumber blocks, signature-curve shapes | 4–55 % |
  | Fourier J = 16–24 (L-dependent) | ≤ 0.1 % |
  | **dictionary** (`dictionary.py`): 6 bulk wavelengths 25–350 mm (sin, cos) + 8 cubic end-layer splines per end over l_e = b_s, physical units, independent of L | **0.05–0.27 %** on the hardest panels (`scripts/exp_dictionary.py`) |

* **Training.** Every gradient-based method converges slowly on this problem:
  * Rayleigh–Ritz loss on separable candidates, single panel: 4–12 % after 600–1000 steps.
  * Neural LOBPCG targets (`spectral.py`), wavenumber-block preconditioner: about 150 iterations to 0.2 %.
  * Energy-whitened outputs.
  * L-BFGS on free variables: 7.8 % after 300 iterations.
  * Only an exact stiffness solve converges fast (5 iterations), and it costs about one GBT solve per panel.
* **Conclusion:** a label-free amortised network costs at least about one exact solve per training panel. The GBT oracle gives these modes at the same cost, so it is used as the physics teacher (as Abaqus was in Paper 1, but automatic).

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
