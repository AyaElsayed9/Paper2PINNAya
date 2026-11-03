# Research log of Paper 2 — every try, result and decision (chronological)

Sources: commit history, `docs/ROADMAP.md`, the reports *Paper2_status_report_simple.pdf* (26 Oct 2026)
and *Paper2_PINN_details.pdf* (27 Oct 2026), run files in `runs/`. Errors are always against full GBT
(solver V32) on designs the network has not seen, unless stated otherwise.

## Fixed requirements from Aya (do not change without asking)
* The network itself must be fast: design in, N_cr (and mode) out in milliseconds, for any design.
  No per-design training and no separate "extension" network as in Paper 1.
* It must stay a **PINN**: the loss is physics (GBT energy). Labels only as a supporting term.
* **No Abaqus in Paper 2** (decision 2 Nov 2026): the validated GBT code is the physics and the
  reference. (Older notes that mention Abaqus blind validation are superseded.)
* Motivation against Paper 1: remove the reduced-bending-stiffness (RBS) approximation of the plate
  energy and the expensive Abaqus data.
* Solver V32 stays frozen (only bug fixes that leave all validated values unchanged, `solver/CHANGES.md`).

## 22 Oct 2026 — foundations
1. **Phase 0** (`c62d27c`): frozen solver, package `pinnbuck`, design space (validated ranges, lamination
   parameters, Sobol sampling), oracle wrappers (strip/GBT), engine (operators, reproduces solver eigenvalue
   to 1e-12), 5 unit + 7 regression tests.
2. **Bug fix** (`83894f6`): spline breakpoints could pass L for odd lengths → ends set exactly to 0, L.
   Pilot: 16 random omega designs, GBT vs strip 0.45 % mean, 1.6 % max.
3. **Shared mode basis** (`fa5433b`).
   * Try: modes of one reference panel (quasi-isotropic) → **failed**, +20 to +130 % on random lay-ups
     (mode order changes with the laminate).
   * Solution: energy-weighted **POD** of the 3 lowest full-GBT modes of 64 panels on a fixed node layout:
     **241 modes** (120 conventional, 1 axial, 60 shear, 30 + 30 Mindlin). Error 0.02 % (symmetry), 0.01 %
     (periodic). A panel's own 80 modes: 1.6 % mean. Members: complete nodal basis, exact.
   * Operators K, G, pre-buckling state in this basis, checked to 1e-13.
4. **First physics-only network** (Rayleigh–Ritz loss on candidate modes + N head): one panel 0.2 % in 1 min;
   8–32 panels: **still 15–140 % after 2000 Adam steps → too slow**.
5. **Research on the length direction** (`761dee9`). Findings: load enters through the eccentric foot, so
   pre-buckling forces near the ends reach 5× the mid-length value; many modes are end-localised (78 % of one
   mode in the two end fifths); lowest loads very close (16.02 / 16.08 / 16.74).
   | description along z | error |
   |---|---|
   | one wave + envelope (Kantorovich) | 1–13 % |
   | one shape per wavenumber (signature curve) | 4–55 % |
   | Fourier 16–24 terms | ≤ 0.1 % |
   | **dictionary**: 6 wavelengths 25–350 mm (sin, cos) + 8 end splines per end | **0.05–0.27 %** |
   Training tries that all stayed slow: neural LOBPCG targets, wavenumber preconditioner (~150 iterations),
   teacher modes as guide, L-BFGS on free variables (7.8 % after 300 iterations). Only an exact stiffness
   solve converged fast (5 iterations) = the cost of a GBT solve.
   **Conclusion then:** label-free training of one network for the family costs about one GBT solve per design.

## 25 Oct 2026 — Attempt 2: GBT as teacher (supervised), later abandoned
* `dcda8f2`: teacher modes (full GBT projected on shared modes × Fourier, polished by LOBPCG, 0.01–0.2 %),
  stored in float64 (float32 destroys the energy of stiff components); network → log N1..3 + 3-mode subspace
  (principal-angle loss); closed-form plate/column loads as inputs; ensemble of 5 + split-conformal intervals
  (to replace the −15 %/+3 % window); member engine (C, Z, L, T, I, box) through `families.py`.
* Results (symmetry, 160 training, 32 test designs): hand formula 20 % mean / 91 % max; **network 10–16 %**;
  network + energy on 200 unlabelled designs 16.1 % (no help); Gaussian process 9.3 % / 30 %; intervals honest
  but wide (±46–81 %, coverage 84–88 %). Error falls slowly with data (40 → 160: 15 → 10 %).
* Plan at that time: overnight laptop run with ~4000 GBT teacher designs. (Status unknown — not in the repo.)
* `58b7426`: physics-reinforced experiment (teacher + Rayleigh–Ritz on unlabelled designs): no gain.

## 26 Oct 2026 — Decision: back to a pure PINN, with the conditioning fixed (current direction)
* **Diagnosis:** the obstacle was conditioning, not physics: cond(K) = 1e9–1e19 in mode coordinates.
* `7d4e89d` **Energy whitening** (`ritzspace.py`, `epinn.py`): Ritz space = shared modes × dictionary with
  fixed labels (same coordinates for every design of one end condition), Kronecker K, compressed G (both
  1e-14 exact); metric M = CᵀC from exact blocks of K → in y = C a the energy ≈ |y|². Eigen-iteration 0.1 %
  in 7–13 steps (was 150). Free variables, one design: raw 4.6 % → **whitened 0.40 %**.
* Design space for this stage: **industrial lay-ups** (symmetric, balanced, 0/±45/90, ≥ 10 % each direction:
  24 skin × 12 stringer sequences); CC with symmetry sides first.
* Optimisation is solved: 7 designs, Adam on physics loss → 0.04 % mean / 0.13 % max.
* **Generalisation history (CC, symmetry):**
  | network / training | unseen designs: head / Ritz bound |
  |---|---|
  | 19 training → 7 unseen | 10 % / 42 % (train 1 %: overfitting) |
  | 259 training, linear output layer | 9.4 % / 10.6 % (train 8.5 %) |
  | 259, function decoder (P2INN-style), Adam | 9.5 % / 20.6 % |
  | 259, decoder + `ngfull` (20 outer, 40 min) | **5.0 % / 12.3 %** (max 10.9 / 25.9 %) |
  | certificate from these predictions, 3/5/8/12 steps | 2.8 / 1.9 / 1.1 / **0.54 %** mean |
* `6a5c96b` **Why the decoder:** analysis showed the modes of unseen designs lie 15–24 % (energy) outside the
  span of the modes of 45 other designs (Kolmogorov n-width barrier) → a linear output layer cannot generalise.
  `FunctionDecoderPINN` = design encoder + code of each longitudinal function → MLP → cross-section vectors.
* `1a0d3d7` **`ngfull`** optimiser: per outer step a physics step for every design (residual → Rayleigh–Ritz on
  span{modes, residuals, previous change} = one LOBPCG step), then 10 epochs of Adam fitting the network to it.
* `27f7f06`, `6d048e5` **Periodic sides, SS ends** (the case of the periodic-GBT paper): shared basis 0.007 %
  mean; dictionary +0.03–0.08 %; n = 6564 unknowns, 7 classes. 320 training → 32 unseen, 30 outer, 95 min:
  **head 3.65 % / 9.77 % max; Ritz 10.7 % / 34.2 %; certificate 12 steps 0.22 % / 1.24 %.**
  Fit loss ≈ 0.2 → about 5 % of each physics target is missing from the predicted subspace.

## 27 Oct 2026 — Paper 1 data term + supervisor
* `b796c65`: data term γ·SmoothL1(N/N_oracle, 1) (Paper 1, Eq. 18) on the 10 % most suspicious designs
  (gap between Ritz bound and head), GBT as oracle, warm start from the periodic model.
* Result: corrects the designs it is applied to, **does not transfer** to unseen designs (head 3.65–8.3 %).
  Adam on the physics loss degrades the head (it learns the Ritz bound of imperfect modes).
* **Plateau:** ≈ 3.6 % head, ≈ 10 % Ritz, also on training designs → the limit is how well the network
  represents the modes (capacity / projection), not data or optimiser steps. Target for the paper: 1–2 %.

## 2–3 Nov 2026
* Periodic-GBT paper finished (manuscript in `manuscript/`). Decision: Paper 2 uses the GBT code, no Abaqus.
* Open next steps (from 27 Oct, not started):
  1. projection quality: more epochs per physics step, larger decoder, hybrid output;
  2. more candidates (R = 16): the Rayleigh–Ritz layer picks the best, so extra candidates act as experts;
  3. experts per lay-up group / geometry range;
  4. longer laptop run: 512 designs, 60 outer steps (curves still falling);
  5. periodic CC, members; then Koiter (phase 4).
* Honest cost note: the certificate needs the design's operators (5–10 s) ≈ one GBT solve (~10 s); the
  instant head (<1 ms) is the fast answer, the certificate is a check.

## 3 Nov 2026 — Aya's plan clarified, Paper 1 read, first test of the plan
* Aya's plan (her words summarised in `docs/PAPER1_SUMMARY.md`): Paper 1's framework with (1) GBT as the
  physics engine instead of the 1D RBS plate energy, (2) the GBT code instead of Abaqus as solver/oracle,
  (3) a reasonable check instead of the −15 %/+3 % neighbour window, (4) the lay-up as an input,
  (5) a fast tool.
* Paper 1 read from the uploaded PDF; facts in `docs/PAPER1_SUMMARY.md` (one cross-ply lay-up, per-panel
  PINN, RBS for the foot, max 10.8 % vs Abaqus at b_s = 124.25 mm, MLP surrogate on b_s and β).
* **Test of point (1): full GBT at the 22 Abaqus points of Paper 1, Table 2** (cross-ply, symmetry side,
  SS; `scripts/gbt_at_paper1_points.py`, `results_p1/gbt_p1.jsonl`; Abaqus values copied from the table):
  | points | Paper 1 PINN vs Abaqus: mean / max | GBT vs Abaqus: mean / max |
  |---|---|---|
  | all 22 | 2.53 % / 10.80 % | **0.54 % / 0.97 %** |
  | tight pitch b_s = 124.25 mm (4) | 8.08 % / 10.80 % | **0.89 % / 0.97 %** |
  | other 18 | 1.29 % / 3.54 % | **0.46 % / 0.95 %** |
  Half-wave number: GBT matches Abaqus in 22/22 (Paper 1: 20/22). GBT is 0.1–1.0 % below Abaqus everywhere.
  **Conclusion: the GBT engine removes the RBS error of Paper 1.** Note: b_s = 124.25 and 136 mm lie
  outside the 150–285.5 mm range of the earlier GBT validation; these Paper 1 points extend it.
* Mapping of the plan onto the existing code: (1) and (2) are what `epinn.py` already does (GBT
  Rayleigh–Ritz loss; GBT as oracle in the data term of 27 Oct); (3) physics certificate / supervisor;
  (4) lamination parameters as inputs (industrial lay-ups); (5) one network for all designs, < 1 ms.
  The open problem remains the accuracy of the fast network on unseen lay-ups (3.65 % mean, target 1–2 %).
