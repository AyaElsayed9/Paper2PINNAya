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
* **No order dependence** (3 Nov 2026): Paper 1 worked best when panels were asked in sweep order (transfer
  learning from the neighbour, −15 %/+3 % window vs the previous panel). In Paper 2 every prediction and
  every check must depend only on the design itself. (The current network has no batch normalisation or
  neighbour input, so this holds.)
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

## 3 Nov 2026 (afternoon) — error analysis of the best network and revised plan
* Analysis of `runs/p1_ss_dec.json` (periodic SS, 32 unseen designs) against the stored full-GBT values:
  * the head is **biased high**: above GBT in 27 of 32 designs, signed mean **+3.13 %** (unconservative).
    Cause: the head is trained on the network's own Ritz values, which are upper bounds of imperfect modes.
  * after removing the mean bias the scatter is still 2.79 % mean, 6.64 % max.
  * |error| is not correlated with b_s (0.0), the gap between the two lowest loads (−0.03) or N_cr (0.04);
    weakly with β (0.2). Several of the largest errors share the skin [−45/45/90/0]s.
* Honest correction: every operator file already stores the exact GBT N_cr and the 3 lowest modes
  (`lam_full`, `teacher`, `teacher_eigs`), because building the operators costs about one GBT solve.
  So "label-free" training saves no computing time here; the exact GBT targets are free.
  The supervised trial of 25 Oct (10–16 %) used the old network (linear output, no whitening, 160 designs)
  and was never repeated with whitening + decoder + ngfull.
* Revised plan (replaces the order proposed in the morning):
  1. hybrid loss as in Paper 1 but with exact GBT targets on all training designs (N_cr and mode subspace)
     plus the Rayleigh–Ritz physics loss; removes the head bias and measures the capacity directly;
  2. only if training designs still plateau: R = 16 candidates as branch experts (N_cr = min over mode
     branches, each branch smooth in the design) and a larger decoder;
  3. more designs where the error is large (active sampling with GBT), ensemble + certificate as the check;
  4. L-BFGS polish and curriculum last.

## 3 Nov 2026 (evening) — three-way test set up: physics-only vs data-only vs hybrid
* **Correction of the afternoon entry:** the exact GBT results are stored only for the test and oracle
  designs. The 320 training operators were built without the full GBT (`lam_full = nan`, about 2.5 s each).
  A label costs about 13 s per design (LOBPCG to convergence in the Ritz space, within 0.03–0.11 % of full
  GBT on 3 test designs), i.e. about 5× the physics data. Label-free training therefore does save time
  (≈ 80 % per design), although labels for a few hundred designs remain cheap (≈ 35 min on 2 cores).
* `scripts/exact_targets.py`: exact modes (R = 8, unit whitened vectors) and loads of the 320 training
  designs → `data/exact_targets/periodic_industrial_SS/`.
* `epinn.train_ritz_ng_full(mode=...)`: "physics" (unchanged), "data" (fit to exact modes, head on exact
  loads, no physics step), "hybrid" (physics-step targets + exact modes, head on exact loads);
  `train_epinn.py --mode --exact`. Reports now include the head bias and the count above GBT.
* Settings identical to the physics-only run `runs/p1_ss_dec` (decoder, R = 8, 30 outer × 10 epochs,
  batch 16, lr 1e-3, seed 0, 320 designs), which is reused as the physics-only arm.
* 19:20: the data-only run was killed by the out-of-memory killer after 45 min (data + hybrid together
  > 8 GB on this machine); nothing to do with the method. `logs/run_threeway2.sh` restarts it alone after
  the hybrid run and computes the certificates meanwhile. Hybrid at outer 16 of 30 (unseen designs):
  head 1.47 % mean / 5.14 % max, bias +0.26 % (15/32 above GBT); Ritz bound 13.4 % (still falling).

## 3 Nov 2026 (night) — three-way test result (SS periodic, 320 training → 32 unseen designs, one seed each)
| | physics-only (`p1_ss_dec`) | data-only (`w3_data`) | hybrid (`w3_hybrid`) |
|---|---|---|---|
| head (instant) mean / max | 3.65 / 9.77 % | 2.78 / 6.25 % | 3.94 / 11.02 % |
| head bias, designs above GBT | +3.13 %, 27/32 | −1.71 %, 8/32 | −1.08 %, 10/32 |
| Ritz bound of predicted modes | 10.68 / 34.17 % | 18.92 / 59.93 % | 14.09 / 37.27 % |
| certificate 3 / 5 / 8 / 12 steps (mean) | 1.78 / 1.01 / 0.46 / 0.22 % | 0.80 / 0.37 / 0.21 / 0.14 % | 0.79 / 0.37 / 0.19 / 0.13 % |
| certificate 12 steps, max | 1.24 % | 0.86 % | 0.70 % |
| Ritz, training designs (16) | not logged | 11.80 % | 9.06 % |
| cost | 95 min (2 threads) | 30 min (1 thread) + labels ≈ 35 min (2 cores) | 144 min (1 thread, first 45 min shared) + labels |
* Labels remove the high bias of the head (the cause found in the afternoon), in both data and hybrid.
* The head error jumps between checkpoints: hybrid 1.47 % at outer 16, 3.94 % at outer 30; data 2.0–3.05 %
  over outer 8–30. With one seed the head difference data vs hybrid (2.78 vs 3.94 %) is within this noise,
  so no winner on the head.
* Mode shapes (Ritz bound on unseen designs): physics-only best (10.7 %), hybrid 14.1 %, data-only worst
  (18.9 %, max 60 %, train/test gap 11.8 → 18.9 %). The physics loss makes the predicted shapes generalise;
  data alone overfits the shapes.
* Certificate: labels help it converge faster (3 steps 0.8 % instead of 1.8 %; 12 steps 0.13–0.14 %, max
  0.70 % for hybrid). The predicted block spans the right directions better although its Ritz value is worse.
* Ensemble (free, from the cert logs): mean of the three heads 1.61 % mean / 5.52 % max, bias +0.11 %;
  data + hybrid only: 2.59 / 7.65 %. Head errors of data and hybrid are uncorrelated (r = −0.04), so the
  error is mostly network variance, not a systematic physics error. Honest caveat: part of the 3-net gain
  is the +3 % and −1 to −2 % biases cancelling. The spread of the three heads does NOT predict the error
  (r = 0.02), so as in Paper 1 the ensemble spread is not a reliable flag; the certificate stays the check.
* Decision proposed to Aya: hybrid as the main training (best certificate, bias gone, shapes still guided by
  physics), then a deep ensemble of 5 seeds (Paper 1 idea) to remove the variance of the head; certificate as
  the check. Not yet started.

## 4 Nov 2026 — decision: quality first (Aya: "I don't care about the time, I care about the quality")
* Method for Paper 2: **hybrid** (GBT energy loss + exact GBT targets). Reasons: it stays a PINN (Aya's
  rule), labels come from GBT (no Abaqus, no RBS error) and are cheap, best certificate. Data-only is not a
  PINN and has the worst shapes: kept only as a baseline. Physics-only is the label-free baseline.
* Plan (each step only kept if the validation numbers improve):
  A. clean training: cosine learning-rate decay, weight average (EMA), 32 validation designs held out of
     the 320, weights chosen on validation, test designs not looked at during training; rerun all three arms.
  B. more training designs with exact GBT targets (320 → about 1000), sampled more where the validation
     error is large; a new final test set (about 100 designs) opened only at the end.
  C. more capacity for the shapes if training error stays high: R = 16 candidates, wider decoder, experts.
  D. ensemble of 5 seeds; test also stacking the candidates of all members in the Rayleigh–Ritz layer.
  E. paper evaluation: new test set, SS and CC, lay-ups, Paper 1's 22 Abaqus points (unchanged numbers),
     timing on Aya's laptop, order independence.
* Step A code: `train_ritz_ng_full(lr_final, ema)`, `train_epinn.py --val --lr-final --ema --no-test-log`
  (saves chosen, last and EMA weights). Smoke test OK. `logs/run_v2.sh` started 12:02: hybrid → data →
  physics, 2 threads, one at a time, with certificates. Note: the physics-only arm also uses the 32
  labelled validation designs, but only to choose its weights.

## 4 Nov 2026 (noon) — Aya's requirements for the paper
1. **Design space like Paper 1**: b_s 124.25–285.5 mm (the tight pitch of Paper 1 included), β up to 3.5,
   and **both side conditions**: symmetry (as Paper 1) and periodic. Paper 1's lay-up (cross-ply) is not in
   the "industrial" lay-up set (no ±45 plies), so the new sets mix general and industrial lay-ups.
   Validation: GBT vs Abaqus at b_s = 124.25 mm exists for symmetry/SS (Paper 1 points, 0.89 % mean);
   periodic below 150 mm has no Abaqus check yet (told Aya; optional few Abaqus runs by her).
2. **Keep every result** (model comparisons, convergence and tuning studies) for the paper:
   `runs/` (weights, logs, json, certificates) and `logs/*.log` are now in git; `scripts/collect_results.py`
   writes `results/experiments.csv` and `results/EXPERIMENTS.md` (one row per run); the exact labels are
   archived in `results/labels/*.jsonl`. The operator files (11 GB) stay out of git: too large, and they
   can be rebuilt from the seeds with `scripts/build_operators.py`.
3. **Count the GBT calls** as Paper 1 counted the Abaqus calls (Paper 1: 144 panels, 84 % fast PINN,
   10.4 % re-solved, 5.6 % Abaqus). Here a "GBT call" = a full GBT buckling solve (label, ≈ 13 s);
   the physics data of a design (GBT matrices, ≈ 2.5 s) is counted separately as GBT time.
* Step B becomes:
  B0 extend the design space; check both shared bases at b_s = 124.25 mm and β = 3.5, rebuild if the
     error is above about 0.1 %;
  B1 operators for SS, symmetry and periodic: about 1000 training + 64 validation + 100 test designs each
     (disk: about 11 MB per design; the old sets are deleted after step A, K maybe stored in float32);
  B2 label study: hybrid with the physics on all 1000 designs and GBT labels on 0 / 50 / 100 / 200 / 400 /
     1000 random designs, and active selection (labels added in rounds where an error indicator is largest;
     the head–3-step-certificate gap is the first candidate, because the ensemble spread did not work);
     plot validation error vs number of GBT calls and GBT time;
  B3 in use: how many new designs need the certificate / a GBT call (the Paper 1 84 / 10.4 / 5.6 % table).
