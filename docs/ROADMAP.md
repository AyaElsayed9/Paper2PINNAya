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

### 2e. Physics-only PINN in energy-whitened coordinates (Oct 2026, current direction)
Decision: the network stays a PINN (physics loss only, no labels); GBT gives the operators and
the test references. Omega with industrial lay-ups first (`designspace.industrial_layups`:
symmetric, balanced, 0/±45/90, ≥ 10 % per direction: 24 skin × 12 stringer sequences), then
further experts (lay-up groups, periodic sides, members).

Diagnosis of 2c/2d: the obstacle was the conditioning, not the physics. The GBT stiffness in the
mode coordinates has a condition number of 1e9–1e19, so every gradient method crawls.

New pieces (`ritzspace.py`, `epinn.py`):
* **Ritz space**: shared modes × dictionary functions with FIXED labels (same coordinates for
  every design of one end condition), K as a sum of Kronecker products, compressed pointwise G;
  both equal the direct z-integrals to 1e-14 (`tests/test_ritzspace.py`).
* **Energy metric**: exact blocks of K per dictionary wavelength and per end layer, M = CᵀC.
  In y = C a the energy is ≈ |y|²: median eigenvalue of M⁻¹K = 1, 99 % below ≈ 5 (cond(K) 1e9–1e19).
  The block-preconditioned eigen-iteration reaches 0.1 % in 9–13 iterations (previously ~150).
* **Network**: design vector + closed-form plate loads → R = 8 candidate modes in y, exact
  Rayleigh–Ritz layer (N = upper bound), N head trained on the network's own Ritz values.
* **Optimiser** (energy natural gradient with Ritz line search, "projected LOBPCG"): per design
  one physics step span{y, residual, previous change} → Rayleigh–Ritz → the network follows it
  with 20 inner Adam steps on the subspace distance. ENGD (Müller & Zeinhofer 2023) with the energy
  metric = plain residual in y; LOBPCG supplies the exact step length.

Results (CC, symmetry, industrial lay-ups; errors vs full GBT, no labels in training):
| test | result |
|---|---|
| free variables, one design, Adam 1000 steps | raw 4.6 %, whitened 0.40 % |
| network, 7 designs (Test 1), 40 outer steps, 4 min, 1 core | mean 0.23 %, max 1.10 % (Ritz and head) |
| same network, Adam on the physics loss, 1000 steps, 21 min | mean 0.04 %, max 0.13 % |
| generalisation: 19 training designs → 7 unseen | head 10 %, Ritz 42 % (train 1 %) |
| 259 unlabelled training designs → 22 unseen (linear output, 10 passes) | test: head 9.4 %, Ritz 10.6 %; train 8.5 % |
| decoder output (P2INN-style), 10 passes | test: head 9.5 %, Ritz 20.6 %; train 19.7 % (still falling) |
| **decoder, full-projection natural gradient (`ngfull`), 20 outer, 40 min, 2 cores** | **22 unseen designs: head 5.0 % mean / 10.9 % max; Ritz 12.3 % / 25.9 %; train 9.3 %; still falling** |
| certificate from these predictions, 3 / 5 / 8 / 12 steps | 2.8 / 1.9 / 1.1 / 0.54 % mean; 6.7 / 5.1 / 3.1 / 1.7 % max |

Periodic side edges, SS ends (the case validated in Verification-Aya-2-10_SS-PBC), industrial lay-ups
(`build_operators.py --side periodic --layups industrial --bc SS`): shared basis vs full GBT on 32
test designs 0.007 % mean, 0.036 % max; same fixed layout for every design (n = 6564, 7 classes
incl. the periodic anchor); energy checks 1e-14; dictionary space +0.03 to +0.08 %; eigen-iteration
with the metric within 0.1 % in 7–11 steps.

Periodic SS result (decoder, `ngfull`, 30 outer steps, 320 unlabelled designs, 95 min on 2 cores),
32 unseen designs, |error| vs full GBT mean / max: head 3.65 / 9.77 %; Ritz of the predicted modes
10.7 / 34.2 %; certificate 3 / 5 / 8 / 12 steps: 1.78 / 1.01 / 0.46 / 0.22 % mean, 8.3 / 5.0 / 2.4 / 1.24 % max.

Paper 1 data term (Elsayed et al., Compos. Struct. 399 (2026) 120959, Eq. 18:
gamma * SmoothL1(N / N_FEA, 1), only where the supervisor calls the oracle), periodic SS.
Supervisor (`scripts/flag_designs.py`): gap between the Ritz bound and the head; the 32 (10 %)
most suspicious training designs get the full GBT. All runs start from the 30-outer periodic model.
| continuation | oracle designs: Ritz / head | 32 unseen designs: Ritz / head |
|---|---|---|
| start | — | 10.7 / 3.65 % |
| Adam 300 steps, data term on N_Ritz (gamma 1000) | 7.8 / 8.9 % | 11.0 / 8.3 % |
| Adam 300 steps, no data term | 10.6 / 10.6 % | 9.5 / 7.1 % |
| ngfull 10 outer, data term on the head (gamma 100) | 12.3 / 3.3 % | 10.3 / 3.66 % |
| ngfull 10 outer, no data term | 11.9 / 5.3 % | 10.2 / 3.65 % |
Findings: the data term corrects the designs it is applied to (as in Paper 1) but does not
transfer to unseen designs; Adam on the physics loss degrades the head (it learns the Ritz upper
bound of imperfect modes); 10 more ngfull steps do not improve unseen designs: plateau at about
3.6 % (head) and 10 % (Ritz of the predicted modes), also on training designs, so the limit is the
representation of the modes by the network (capacity), not data or optimiser steps.

Open: (i) the projection step (network fit to the physics step) limits the speed: the fit loss
falls only 0.98 → 0.30 over 20 outer steps; (ii) longer training (laptop); (iii) more candidates
(Ritz mixture of experts: the Rayleigh–Ritz layer itself selects the best candidates, so extra
candidates act as experts without a gating network); (iv) experts per mode family or geometry
range; (v) SS and periodic sides with the same scripts (`--bc SS`).

Diagnosis of the generalisation plateau (analysis only, converged modes not used in training):
the whitened modes of unseen designs lie 24 % (energy) outside the span of the modes of 45
other designs (Kolmogorov n-width barrier), so a linear output layer, whose outputs live in a
fixed span, cannot generalise. Hence `FunctionDecoderPINN`: design encoder + code of every
longitudinal function (wavelength, sin/cos, end layer, class) → MLP → cross-section vectors.

Physics certificate (`scripts/eval_epinn.py`, replaces the −15 %/+3 % check of Paper 1):
operators of the design (~10 s) + block LOBPCG with the energy metric started from the
predicted modes (~0.3 s per step); every value is an upper bound in the Ritz space.

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
