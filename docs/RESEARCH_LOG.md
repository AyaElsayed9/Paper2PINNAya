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

## 4 Nov 2026 (12:30) — Abaqus cells for the extended design space (to be run by Aya)
* `abaqus_tasks/2026-11-04_extension/`: 20 cells. Positions from a fixed-seed Sobol sequence, named lay-ups
  of the Abaqus script:
  * A: periodic/SS, 1 control + 8 (b_s 124.25–150, β to 3.5; two with β 3.0–3.5 at larger b_s);
  * B: symmetry/SS, 1 control + 6 (tight pitch, Sym1/Sym2/Unsym);
  * C: periodic/CC, 4 (optional).
  The control cell is b_s 150, β 1.5, Sym1, already in the validated sets.
* Files: `run_extension.ps1` (PowerShell, uses `Omegastiffener_PBC_v2.py` with PANEL_SIDES = PBC / SYM).
* **Blind:** GBT V32 values of all 20 cells (`gbt_predictions.csv`) computed and committed before any Abaqus
  result exists. Control: periodic 17.315, symmetry 17.312.

## 4 Nov 2026 (14:00) — step A, hybrid arm (`runs/v2_hybrid`)
Settings: 288 training + 32 validation designs (the last 32 training files), 30 outer steps, cosine LR
1e-3 → 1e-5, EMA 0.995, weights chosen on validation (head mean), test never looked at during training.
| run | head mean / max, bias | Ritz bound | certificate 3 / 12 steps (mean, 12-step max) |
|---|---|---|---|
| w3_hybrid (old, 320 designs, constant LR, last weights) | 3.94 / 11.02 %, −1.08 % | 14.1 % | 0.79 / 0.13 % (0.70) |
| v2_hybrid (chosen: EMA weights, outer 8) | **1.42 / 4.35 %, +0.55 %** | 20.1 % | 1.81 / 0.20 % (0.98) |
Findings (honest):
* The head is much better: 1.42 % mean, 4.35 % max, almost no bias (19 of 32 above GBT).
* The **head overfits after outer ≈ 8–10**: training head 1.3 → 0.77 %, validation head 1.0 → 1.9 %.
  The mode shapes do not overfit: validation Ritz keeps falling (20.8 % at outer 8 → 14.5 % at 30).
* So choosing the weights by the head alone takes early weights with poorer modes, and the certificate is
  slower (3 steps: 1.81 % vs 0.79 %). At 12 steps both are 0.1–0.2 %.
* The EMA weights are steadier and better than the raw weights on validation (outer 8: 1.02 vs 1.54 %).
* The validation set is small (32 designs) and the minimum over 15 evaluations is a little optimistic;
  the test value (1.42 %) confirms it is not far off.
Ideas from this (not done yet): choose the head weights and the mode weights separately (two weight sets of
the same network, or freeze the head early); regularise the head; more designs (step B1: about 1000),
which is the most direct fix for the overfitting.
* 4 Nov 14:00: the machine restarted and stopped the data-only arm at outer 4 (log kept as
  `runs/v2_data_killed_by_restart.log`). 16:03: data and physics arms restarted with the same settings
  (`logs/run_v2b.sh`).

## 4 Nov 2026 (16:35) — step A, data-only arm (`runs/v2_data`), same settings as v2_hybrid
| run | head mean / max, bias | Ritz bound | certificate 3 / 12 steps (mean, 12-step max) |
|---|---|---|---|
| w3_data (old) | 2.78 / 6.25 %, −1.71 % | 18.9 % | 0.80 / 0.14 % (0.86) |
| v2_data (chosen: EMA, outer 10) | **1.21 / 4.79 %, +0.05 %** | 22.8 % | 1.10 / 0.14 % (0.82) |
* Same pattern as the hybrid arm: the head is best at outer ≈ 10 (validation 0.96 %) and then overfits
  (training head 0.6 %, validation 1.7 % at outer 30). Validation Ritz still falls after that (22.5 → 18.6 %),
  but stays above the hybrid's (14.5 % at outer 30): the physics step gives better mode shapes.
* Without the physics step one outer step is much faster (whole run about 30 min instead of about 110 min).
* The head difference hybrid 1.42 % vs data 1.21 % is within the noise of one seed and 32 test designs.

## 4 Nov 2026 (18:20) — step A, physics-only arm (`runs/v2_physics`) and summary of step A
Same settings (288 training + 32 validation designs, cosine LR, EMA, choice on validation, one seed, SS
periodic industrial set, 32 unseen test designs).
| arm | chosen | head mean / max, bias (above) | Ritz bound | cert. 3 steps | cert. 12 steps mean / max | time |
|---|---|---|---|---|---|---|
| physics-only | net @ 30 | 4.48 / 11.89 %, +4.47 % (31/32) | **11.95 %** | 1.92 % | 0.22 / 1.07 % | 104 min |
| data-only | EMA @ 10 | **1.21 / 4.79 %, +0.05 %** (13/32) | 22.8 % | 1.10 % | 0.14 / 0.82 % | 26 min |
| hybrid | EMA @ 8 | 1.42 / 4.35 %, +0.55 % (19/32) | 20.1 % | 1.81 % | 0.20 / 0.98 % | 110 min |
Old runs (constant LR, last weights): physics 3.65 % / +3.13 % / Ritz 10.7 %; data 2.78 %; hybrid 3.94 %.
Findings (honest):
* The physics-only head is biased high on almost every design (+4.5 %). Its target is the Ritz value after
  the physics step, which is an upper bound, so the head learns an upper bound. It was still falling slowly
  at outer 30 (validation 5.2 → 4.8 % from outer 20 to 30), and it does not overfit. It is worse than the
  old physics run (3.65 %), which had 320 designs instead of 288 (and the same 30 outer steps).
* Labels fix the head: 1.2–1.4 % mean, max under 5 %, no bias. But with labels the head overfits after
  outer 8–10.
* The physics gives the best mode shapes (Ritz 12 % against 20–23 %). The weights chosen for the head in
  the labelled arms are early weights, so their modes are poorer.
* After 12 certificate steps all three are at 0.14–0.22 % mean, below 1.1 % max: the certificate makes the
  final answer good in every case; the head and the modes decide how many steps are needed.
Decision: keep hybrid as the method, and fix the conflict between "early weights for the head" and "late
weights for the modes": (1) a separate choice for the head and for the modes (two weight sets of the same
network, the cheap fix); (2) more designs (step B1, about 1000), which should reduce the overfitting;
(3) seeds (step D) before any claim that hybrid and data-only differ (1.42 vs 1.21 % is within the noise of
one seed and 32 test designs). Next: step B0 (extend the design space as Paper 1, check the shared bases).

## 4 Nov 2026 (evening) — step B0 started
* Code: `designspace.SPACES` (v1 = old ranges; v2 = b_s 124.25–285.5 mm, β 0.75–3.5), `Sample.space` (input
  scaling; not part of the design id), `sample_design(..., layups="mixed", space="v2")` (each design industrial
  with probability 1/2, else general ply-by-ply). Checked: the old sets are reproduced exactly (same ids).
  `build_operators.py --space v2 --layups mixed`.
* Code: `train_epinn.py` also keeps the weights with the best validation Ritz bound (`state_modes`);
  `eval_epinn.py --modes state_modes|state_last|state_ema` takes the mode shapes (certificate start) from
  other weights than the head.
* `scripts/check_extended_space.py` (running): full GBT vs shared basis vs converged Ritz space at 20 designs per
  side condition (corners of the new region, Sobol points with industrial lay-ups, controls in the old region).
* Incident (19:00): the check with 2 parallel workers ran out of memory (each worker about 3 GB; the
  memory limit killed one, the process pool then waited forever with 3 periodic designs left). 17 results
  were already saved; the first log is kept as `logs/b0_check_periodic_run1_oom.log`. Restarted with one
  worker (resumes from the saved results). Lesson: one design per process for the full GBT + Ritz check.
* First finding: shared basis error at most 0.144 % (Unsym, b_s 124.25, β 1.5), Ritz error 0.003–0.27 %,
  except **b_s 285.5, β 3.5, CC (L ≈ 1000 mm): Ritz +5.29 %** with the shared basis fine (+0.002 %). Two
  possible causes: (1) LOBPCG not converged (every run stopped at the 80-iteration cap); (2) the dictionary
  is too short (longest half-wave 350 mm). `scripts/b0_long_cc_probe.py` tests both on six panels: 400
  iterations, and the dictionary extended by the next two half-waves of the same series (593, 1006 mm).
* Head/modes split on the step-A checkpoints (same 32 test designs; head from the chosen weights, mode
  shapes from the weights at the end of training, so no choice on the test set; `runs/v2_*_cert_modes_*.log`):

  | arm | mode weights | Ritz mean/max | cert 3 mean/max | cert 5 | cert 12 mean/max |
  |---|---|---|---|---|---|
  | hybrid | chosen (ema@8) | 20.10/49.50 % | 1.81/6.84 % | 0.92 % | 0.20/0.98 % |
  | hybrid | final EMA | 14.67/36.74 % | 0.94/4.19 % | 0.44 % | 0.14/0.78 % |
  | hybrid | final weights | 14.69/36.93 % | 0.94/4.18 % | 0.44 % | 0.14/0.78 % |
  | data | chosen (ema@10) | 22.80/61.80 % | 1.10/4.84 % | 0.49 % | 0.14/0.82 % |
  | data | final EMA | 19.77/68.12 % | 0.88/4.10 % | 0.40 % | 0.14/0.85 % |

  The split halves the 3-step certificate error of hybrid (1.81 → 0.94 %) and makes hybrid's modes the best of
  the labelled arms. Decision: use separate weights for head and modes from now on (`state_modes`, already in
  `train_epinn.py`).

## 5 Nov 2026 (morning) — step B0 results
(The machine restarted during the night; both the check and the first probe had finished before.)
* B0 check, 20 designs per side (`results/b0_check/{periodic,symmetry}.jsonl`):

  | side | shared basis mean / max | Ritz space (80 it) mean / max | worst Ritz designs |
  |---|---|---|---|
  | periodic | 0.023 / 0.144 % | 0.37 / 5.29 % | 285.5/3.5 CC 5.29 %; 125.4/3.37 SS 0.42 %; 285.5/3.5 SS 0.27 % |
  | symmetry | 0.028 / 0.116 % | 0.45 / 5.33 % | 285.5/3.5 CC 5.33 %; 124.25/3.5 CC 0.87 %; 125.4/3.37 SS 0.65 % |

  The shared bases (built from b_s 150–285.5) stay accurate in the new region (max 0.14 %, at b_s 124.25):
  **no rebuild needed.** The problem is the Ritz space.
* Probe 1 (`b0_long_cc_probe.py`, periodic, `results/b0_check/long_cc_probe_periodic.jsonl`):

  | panel | L | 80 it | 400 it | + half-waves 593, 1006 mm |
  |---|---|---|---|---|
  | Sym1 285.5/3.5 CC | 999 | +5.292 % | +5.292 % (stopped at 93 it, converged) | **+0.012 %** |
  | Sym1 285.5/3.0 CC (old corner) | 856 | +0.122 % | +0.119 % | +0.012 % |
  | Unsym 285.5/3.0 CC | 856 | +0.026 % | +0.022 % | +0.008 % |
  | Sym1 285.5/2.0 CC | 571 | +0.027 % | +0.023 % | +0.017 % |
  | Cross-ply 285.5/3.5 SS | 999 | +0.268 % | +0.266 % | +0.068 % |
  | Unsym 124.25/3.5 CC | 435 | +0.050 % | +0.029 % | +0.030 % |

  Cause found: **the dictionary was too short** (longest half-wave 350 mm), not LOBPCG. With the next two
  half-waves of the same series the long clamped panel drops from 5.29 % to 0.012 %, and the old-space corner
  (285.5/3.0 CC, 0.12 %) also improves. Cost: Ritz size 6087 → 7051 (+16 %). The 80-iteration cap costs at
  most 0.02 % here.
* Probe 2 running (`b0_dictionary_probe.py`): all check designs with Ritz error > 0.3 % and the controls, both
  sides, with base / long / long + short (14.75 mm) dictionaries, to see whether the short panels at b_s ≈ 125
  (0.4–0.9 %) need a shorter wave or only more iterations.
* Probe 2 was stopped by the machine after its first design per side (no error in the log; jobs started
  detached from a session are not kept alive). Restarted as a tracked job; it resumes from the saved results.
* First probe-2 result, symmetry 124.25/3.5 CC: 80 it +0.874 %, 400 it +0.702 %, long +0.653 %, long + short
  +0.618 %. So for this short clamped panel neither more iterations nor the extra waves remove the error;
  something else limits it. Probe 3 (`b0_dictionary_probe.py --set fine`, queued) tests a denser wave series
  (15 half-waves 14.75–1006 mm, ratio 1.35 instead of 1.70) and longer, finer end layers (12 splines over
  1.5 b_s) on the designs with Ritz error > 0.4 %.
* Code: versioned dictionary. `dictionary.DICTIONARIES` = d1 (old, 25–350 mm) and d2 (+593, 1006 mm, may still
  change after probe 3); the operator remembers its dictionary (`op.info["dictionary"]`, set by
  `build_operators.py --dictionary`, default d2 for space v2); `epinn.prepare` uses it, so all old operator sets
  and results stay exactly as they were (d1).
* Probe 2 done (`results/b0_check/dictionary_probe_{symmetry,periodic}.jsonl`; Ritz error vs full GBT):

  | side | panel | 80 it | 400 it | long (d2) | long + 14.75 mm |
  |---|---|---|---|---|---|
  | sym | 124.25/3.5 CC | 0.874 | 0.702 | 0.653 | 0.618 |
  | sym | 285.5/3.5 SS | 0.300 | 0.267 | 0.076 | 0.076 |
  | sym | 285.5/3.5 CC | 5.331 | 5.330 | **0.014** | 0.014 |
  | sym | 125.4/3.37 SS | 0.651 | 0.280 | 0.282 | 0.277 |
  | sym | controls 150/1.5 SS, CC; 285.5/1.5 SS | 0.125, 0.092, 0.104 | 0.123, 0.040, 0.088 | 0.123, 0.042, 0.087 | 0.112, 0.025, 0.085 |
  | per | 285.5/3.5 CC | 5.292 | 5.292 | **0.012** | 0.012 |
  | per | 125.4/3.37 SS | 0.423 | 0.305 | 0.305 | 0.301 |
  | per | controls (same) | 0.117, 0.012, 0.090 | 0.116, 0.009, 0.075 | 0.116, 0.009, 0.074 | 0.104, 0.008, 0.073 |

  Findings: (1) the long waves fix the long panels and change nothing elsewhere (good: no harm); (2) the extra
  short wave (14.75 mm) gains at most 0.04 % for +7 % size: not worth it; (3) for the short panels (b_s ≈ 125)
  the 80-iteration cap was costing up to 0.37 % (0.65 → 0.28 %): the labels (exact targets) of the v2 sets
  will use 400 iterations; (4) one panel is still at 0.65 % (symmetry 124.25/3.5 CC): probe 3 (running) looks
  at wave spacing and end layers for it.
* Probe 3 done (`results/b0_check/fine_probe_*.jsonl`, 400 iterations; size n of the Ritz space in brackets):

  | side | panel | long d2 (7051) | dense 15 waves (10425) | long + end 12 (8979) | dense + end 12 (12353) |
  |---|---|---|---|---|---|
  | sym | 124.25/3.5 CC | 0.653 | 0.359 | 0.488 | 0.348 |
  | sym | 285.5/3.5 CC | 0.014 | 0.013 | 0.013 | 0.011 |
  | sym | 125.4/3.37 SS | 0.282 | 0.075 | 0.086 | 0.069 |
  | per | 285.5/3.5 CC | 0.012 | 0.011 | 0.011 | 0.009 |
  | per | 125.4/3.37 SS | 0.305 | 0.080 | 0.092 | 0.074 |

  (n for SS designs is about 480 larger.) The wave spacing is what limits the short panels: the denser
  series (ratio 1.35 instead of 1.70) brings them from 0.3 % to 0.08 % and the worst corner from 0.65 % to
  0.36 %; finer end layers add little on top. Price: Ritz space +48 %.
* Final B0 check started (`--set final`): all 20 check designs per side at 400 iterations with three candidate
  dictionaries: long (6+2 waves, ratio 1.70), mid (11 waves 25–1006 mm, ratio 1.45), dense (15 waves
  14.75–1006 mm, ratio 1.35). The choice is made on these 40 check designs only (no test designs involved).
* Final B0 check done (`results/b0_check/final_probe_*.jsonl`, all 20 check designs per side, 400 iterations;
  Ritz error vs full GBT, mean / max, and Ritz size n):

  | side | long (6+2 waves) | mid (11 waves, ratio 1.45) | dense (15 waves, ratio 1.35) |
  |---|---|---|---|
  | symmetry | 0.104 / 0.653 % (n 7337) | 0.077 / 0.429 % (n 8783) | 0.064 / 0.359 % (n 10711) |
  | periodic | 0.068 / 0.305 % | 0.053 / 0.220 % | 0.045 / 0.214 % |

  Smallest values −0.003 to −0.011 %: the Ritz value can be a little below the full GBT because the two use
  different longitudinal discretisations (not an error of the bound in the shared basis).
  **Decision: d2 = mid, 11 half-waves 25–1006 mm.** It takes most of the gain (worst 0.65 → 0.43 %, short
  SS panels 0.3 → 0.1 %); dense adds only 0.01–0.07 % for another +22 % size; the rest of the error at the
  worst points (0.2 % at 124.25/1.5 SS, 0.12 % at the 150/1.5 control) does not move with the waves at all.
  `scripts/exact_targets.py --iters` added (400 for the v2 sets).
* B0 closed: shared bases kept; dictionary d2; label iterations 400.

## 5 Nov 2026 (midday) — step B1 started, and an open question on general lay-ups
* Disk: deleted the old train operator sets (`periodic_industrial_SS/train`, `symmetry_industrial/train`, 9.3 GB;
  rebuildable from seeds 11, n 320 / 512) and the pip cache; kept both old test sets (step-A comparisons).
* Smoke test (2 v2 designs, `data/operators/periodic_mixed_SS_v2/smoke`): build and exact targets work, the
  operators carry `dictionary = d2`. One design with a GENERAL lay-up (skin [30/60/45/60]s, stringer
  [−60/−30/60/60/60/−30/−60], b_s 167.6, β 2.85, L 477) has Ritz error +0.376 % (400 it) / +0.363 % (1000 it)
  with the shared basis at only +0.041 %. The dense waves do not help (+0.42 % at 400 it); longer end layers
  do: 12 splines over 1.5 b_s +0.323 %, 16 splines over 2.5 b_s (capped at L/2) +0.193 %. So for general
  lay-ups part of the mode sits in a long transition zone near the ends that the 8 end splines over b_s
  cannot follow. Open: how often this happens. The operator files do not depend on the dictionary (it is only
  applied when a design is prepared), so the B1 build goes on, and the distribution of the Ritz error over
  all ~1000 B1 designs (exact targets) will decide whether the end layers need to change before training.
* Build of v2 sets started (`logs/run_b1_build.sh`): periodic SS, mixed lay-ups, full GBT for every design;
  test 100 (seed 99), train 1064 = 1000 + 64 validation (seed 11); no overlap of ids; 170 train designs
  have b_s < 150 and 192 have β > 3 (the new region).
* Incidents: (1) jobs started detached from the session are stopped a few minutes after a reply (the first
  build stopped after 24 designs, as probe 2 did); long jobs now run as tracked jobs that resume where they
  stopped. (2) My mistake: a check for partly written files deleted the 24 good files, because the check
  could not import the package (it reported every file as unreadable). 24 designs (about 4 min) rebuilt. Fix:
  `build_operators.py` now writes each file atomically (`.tmp` then rename), so a stopped job can never leave
  a partial file.

## 5 Nov 2026 (afternoon) — B1 end-layer probe: the gains of longer end layers are partly an artefact
* End-layer probe (`scripts/b1_end_layer_probe.py`, first 24 test designs of `periodic_mixed_SS_v2/test`,
  d2 waves, 400 iterations; `results/b1_end_layer_probe/periodic_mixed_SS_v2_test.jsonl`). Ritz error vs the
  full GBT label, mean / max / min:

  | lay-ups | e8 (8 over b_s, current) | e12 (12 over 1.5 b_s) | e16 (16 over 2.5 b_s) |
  |---|---|---|---|
  | industrial (13) | +0.076 / +0.160 / +0.023 % | +0.011 / +0.087 / −0.074 % (5 below 0) | −0.035 / +0.054 / −0.150 % (9 below 0) |
  | general (11) | +0.123 / +0.256 / +0.041 % | +0.087 / +0.186 / +0.037 % | +0.059 / +0.141 / +0.016 % |

  Ritz size n 8974 / 10902 / 12830; mean LOBPCG time 53 / 66 / 112 s.
  So the general lay-ups are a little worse with e8 (max 0.26 % here; the smoke design was 0.36 %), not a
  problem of the size seen in B0 for the long panels.
* The values BELOW the full GBT label (down to −0.15 %) were checked, because a Ritz value should not be lower
  than a converged reference. They appear only for short panels where the end splines (spacing L/(2 n_end)
  after the L/2 cap) are finer than the solver span hx (≈ b_s/20, max 15 mm). Two checks on design
  040508c6cbc1 (b_s 282.9, β 1.19, L 336, hx 14.15 mm, e16 −0.150 %):
  1. Full GBT with only the solver's hx argument changed (`scripts/b1_hx_convergence.py`,
     `results/b1_end_layer_probe/hx_convergence.jsonl`): hx 6.30873 (= stored label, reproduced to 1e-8 %),
     hx/2 6.30646 (−0.036 %), hx/4 6.30520 (−0.056 %). The label is converged along z to about 0.06–0.08 %
     for this short panel (slow convergence, most likely the pre-buckling end zones). hx/8 did not fit in
     memory next to the running build (killed, nothing else affected).
  2. The operator rebuilt on a 2× finer grid (hx/2: finer pre-buckling field and twice the Gauss points;
     frozen solver untouched) and the three Ritz spaces again: e8 / e12 / e16 = +0.072 / +0.048 / +0.026 %
     ABOVE the full GBT on that grid (shared basis alone +0.005 %); vs the original label −0.0095 % for e16
     instead of −0.150 %. So on the solver grid the e16 value was 0.14 % too low: an integration artefact
     of splines finer than the grid, not a better mode. (Script for this check: `scripts/b1_fine_grid_check.py`;
     this first run was done with a scratch copy of it, numbers above.)
* Meaning: (a) the negative e12/e16 values on industrial lay-ups are artefacts and cannot be counted as gains;
  (b) for general lay-ups part of the e16 gain may be real, but the same artefact can be in it; (c) the
  remaining e8 gap (about 0.08–0.12 % mean) is of the same size as the label's own longitudinal
  discretisation error (≈ 0.06–0.08 % on the checked design). Honest limit for the paper: the full GBT label
  at the solver's default hx carries about 0.05–0.1 % discretisation error on short panels.
* Decision (provisional): keep d2 with the e8 end layers (no change). Final check after the build, when memory
  is free: `b1_fine_grid_check.py --factor 0.5` on about 6 designs (3 general with the largest e8 errors,
  3 industrial). If on the fine grid e16 (or e12) is still clearly better for general lay-ups (> 0.05 % on
  average) a version d2e (longer end zone with spacing never below 2 hx) is considered before the exact
  targets. The fine-grid run for the general design 3265ef5ed33c (L 327, 42 spans) was killed for memory
  while the build ran; it waits for the build to finish.
* 14:30 the machine restarted (third time today); build stopped at train 338/1064, no partial files (atomic
  write works). Resumed as a tracked job; `run_b1_build.sh` now appends to the test log instead of
  overwriting it. Disk: 16 GB free, the train set needs about 9 GB more (≈ 12.8 MB per design), so before
  the symmetry v2 sets the operator storage has to be made smaller (or sets removed after use).

## 5 Nov 2026 (late afternoon) — B1 build done; operator files compressed
* B1 build finished: `periodic_mixed_SS_v2/test` 100 and `/train` 1064 designs, full GBT for every design,
  no failed design (0 errors / NaN in both logs). Wall time of the train set about 3.6 h on 2 workers over
  three runs (two machine restarts), 13.3 MB per file → disk 6.6 GB free.
* Operator files are now stored compressed (`src/pinnbuck/opio.py`: zlib on the pickle, lossless, 13.3 → 5.8 MB,
  0.05 s to unpack; `load` also reads the old plain files). `build_operators.py` writes compressed files;
  `epinn`, `exact_targets`, `train_epinn`, `eval_epinn` and the B1 probes read through `opio.load`.
  `scripts/compress_operators.py` converts a set in place and checks every file bit for bit after writing
  (tested on 55 copies first; `epinn.prepare` works on the compressed files). Train set running now; test set
  after the fine-grid check (which reads it).
* Fine-grid end-layer check started (`scripts/b1_fine_grid_check.py --factor 0.5`, 6 test designs: general
  3265ef5ed33c, 1d51e6f13ec3, 23665413e421; industrial 040508c6cbc1, 29027bd61ab1, 11ef6acdbbb8).
* Train operators compressed: 1064 files, 14.25 → 6.17 GB, all checked bit for bit; disk 15 GB free.
* Fine-grid check done (`results/b1_end_layer_probe/fine_grid.jsonl`; operator rebuilt with hx/2, Ritz vs the
  full GBT on the same fine grid, so this is a true Rayleigh–Ritz gap; all values positive):

  | design | lay-up | e8 coarse vs label | e8 fine | e12 fine | e16 fine | shared basis fine | label vs fine full |
  |---|---|---|---|---|---|---|---|
  | 040508c6cbc1 (b_s 283, L 336) | ind | +0.023 | +0.072 | +0.048 | +0.026 | +0.005 | +0.036 |
  | 29027bd61ab1 (b_s 170, L 250) | ind | +0.085 | +0.139 | +0.085 | +0.052 | +0.012 | +0.049 |
  | 11ef6acdbbb8 (b_s 124, L 204) | ind | +0.160 | +0.194 | +0.144 | +0.115 | +0.083 | +0.033 |
  | 3265ef5ed33c (b_s 161, L 327) | gen | +0.256 | +0.348 | +0.239 | +0.164 | +0.018 | +0.091 |
  | 1d51e6f13ec3 (b_s 203, L 354) | gen | +0.228 | +0.313 | +0.209 | +0.148 | +0.046 | +0.069 |
  | 23665413e421 (b_s 127, L 362) | gen | +0.192 | +0.238 | +0.229 | +0.184 | +0.073 | +0.043 |

  (all in %). Findings: (1) on a correctly integrated grid the order is always e16 < e12 < e8, so the longer
  end zone is a real gain, mean e8 → e16 0.10 % (general lay-ups 0.13 %); (2) on the solver grid ALL the
  variants read too low by 0.03–0.1 % (e8 too), and the fine splines of e12/e16 on short panels by up to
  0.19 %; (3) the label itself is 0.03–0.09 % above the fine-grid full GBT. Rule set before the check
  (> 0.05 % gain for general lay-ups on the fine grid) is met.
* So a spacing limit was added to the dictionary (`dictionary.end_spacing`, option `min_spacing` in units of
  the solver span hx; default 0 = old behaviour, reproduced to 4e-10 on a test design): the zone is
  min(end_length b_s, L/2), the spacing max(zone / n_end, min_spacing hx), and splines that do not fit in
  L/2 are left out (fixed labels, missing functions are dropped as for inadmissible ones). Probe 2
  (`--set spacing`, same 24 designs): e8s15 (8 over b_s, ≥ 1.5 hx), e16s15 (16 over 2.5 b_s, ≥ 1.5 hx),
  e16s2 (≥ 2 hx). Decision rule: the variant with the lowest mean error that has no value below the label.
  Test operators compressed before the probe.
* Probe 2 done (`results/b1_end_layer_probe/periodic_mixed_SS_v2_test_spacing.jsonl`, same 24 designs, solver
  grid, Ritz error vs label, mean / max / min, mean Ritz size n, mean LOBPCG time at 2 threads):

  | lay-ups | e8 (old) | e8s15 | e16s15 | e16s2 |
  |---|---|---|---|---|
  | industrial (13) | +0.076 / 0.160 | +0.077 / 0.160 / +0.023, n 8900, 33 s | **+0.053 / 0.116 / +0.012**, n 10198, 40 s | +0.088 / 0.158 / +0.027, n 9048, 35 s |
  | general (11) | +0.123 / 0.256 | +0.123 / 0.256 / +0.041, n 8974, 40 s | **+0.078 / 0.141 / +0.026**, n 11954, 55 s | +0.106 / 0.214 / +0.031, n 10902, 49 s |

  No variant has a value below the label any more. e16s2 is too coarse for the short panels (fewer splines
  fit in L/2). e8s15 = e8 except one design (its old splines were finer than 1.5 hx).
* **Decision: d2 end layers = e16s15** (16 splines over 2.5 b_s, spacing ≥ 1.5 hx), by the rule set before
  probe 2 (lowest mean, nothing below the label). Cost: Ritz size +14 % (industrial) to +33 % (general), LOBPCG
  time +20–40 %. Stored as the version: `dictionary.END_LAYERS["d2"]`, `dictionary_kw(op)`; `epinn.prepare`
  uses all settings of the operator's version (d1 sets unchanged).
* `exact_targets.py`: targets now record the dictionary version, its settings, n, labels and the time;
  written atomically. Smoke on 040508c6cbc1: N1 6.31013 (probe with R 4: 6.31017), +0.022 % vs label, 42 s.
* Exact targets started for test 100 and train 1064 (`logs/run_b1_targets.sh`, 2 parts in parallel, R 8,
  400 iterations).
* 5 Nov: the target job stopped at the 2 h limit (test 100 done, train 30); restarted (resumable). Speed about
  55 s per design with 2 parts, so train needs about 16 h more.
* Test set (100) summary (`scripts/b1_targets_summary.py`, `results/b1_targets/periodic_mixed_SS_v2_test.csv`),
  Ritz error of the network's shape space vs the full-GBT label: all mean 0.080 %, median 0.048 %, p95 0.151 %;
  industrial (45) mean 0.050 / max 0.167 %; general (55) mean 0.104 / max **1.535 %**. No negative values.
  3 designs above 0.2 %, all general lay-ups with b_s < 150 mm and a high N (30–40 N/mm, local modes):
  5cc7be690865 +1.535 % (b_s 136.7, β 3.30, skin [60,45,90,45,60,30,60,60], unbalanced), 6cc2e9904ad4
  +0.465 %, 4a8f2d5f7b80 +0.266 %. Full-GBT time per design (build, 1 thread) mean 18 s (6–37 s); LOBPCG
  target mean 111 s (39–183 s, 1 thread each).
* Outlier probe started (`scripts/b1_outlier_probe.py`): denser half-wave lists (21, 31, 23 from 12 mm) on
  the two worst designs, to tell a dictionary gap from a cross-section-basis gap.
* Outlier probe, rounds 1–2 (`results/b1_targets/outlier_probe.jsonl`, `shared_basis_check.jsonl`), Ritz error vs label:

  | design | d2 (400 it) | d2 1500 it | 21 half-waves | 31 half-waves | 31, 1500 it | shared basis, full splines in z |
  |---|---|---|---|---|---|---|
  | 5cc7be690865 | +1.532 | +1.503 | +1.388 | +1.049 | +0.891 | +0.124 |
  | 6cc2e9904ad4 | +0.467 | +0.447 | +0.431 | +0.371 | +0.313 | +0.040 |
  | 4a8f2d5f7b80 | +0.266 | | | | | +0.181 |

  So (c) iterations: not the cause (400 → 1500 changes 0.02–0.03 %); (b) shared cross-section basis: small for
  5cc7 / 6cc2 (0.12 / 0.04 %) but most of 4a8f's gap (0.18 of 0.27 %); (a) the longitudinal functions are
  the main cause for 5cc7 / 6cc2, and denser fixed half-waves close it only slowly. Note the waves are
  multiplied by the global SS end factor ζ(1−ζ), so a fixed half-wave that does not fit L is spoilt at the
  far end. Round 3: half-waves L/m (the panel's own SS harmonics; design dependent, diagnostic only).
* Round 3 (half-waves L/m, the panel's own SS harmonics): 5cc7 +1.532 → +0.894 % (with d2 added +0.837 %),
  6cc2 +0.467 → +0.274 %, 4a8f no change (+0.265 %, its gap is the cross-section basis). Better, but still far
  above the shared-basis floor (0.12 / 0.04 %). Suspect: the global SS end factor ζ(1−ζ) multiplies every wave,
  so even the exact harmonic sin(mπz/L) enters as ζ(1−ζ) sin(mπz/L) and the uniform amplitude of an
  orthotropic-like mode cannot be built.
* Round 4: local end factor (1 in the bulk, C1 ramp over l_env from each end: zero class 1−(1−s/l)², clamp
  3(s/l)²−2(s/l)³), l_env = b_s (env1) or 2.5 b_s (env25), with the d2 half-waves and e16s15 end layers; the 3
  outliers plus 3 normal designs (040508c6cbc1, 3265ef5ed33c, 11ef6acdbbb8) to check nothing gets worse.
* Round 4 result: the local end factor does NOT help (5cc7 +1.487 % with l = b_s, env25 identical to d2 because
  l is capped at L/2 and then spans the same space; normal designs unchanged within 0.01 %). Hypothesis of the
  global parabola rejected; the end splines already supply that freedom.
* Mode analysis (teacher mode 1 of the operator, projection on sin(mπz/L)): 5cc7 and 6cc2 have modes
  localised at BOTH ends (amplitude about 0.02–0.1 of the maximum in the middle third, non-zero at the ends
  for the free classes) made of short waves (5cc7: m = 11, 9, 13, half-waves 41, 50, 35 mm; 6cc2: m = 8, 10,
  6, half-waves 63, 50, 84 mm). The normal design 040508c6cbc1 has a bulk mode (m = 3, 1). So the hard case is
  end-localised local buckling with short waves on tight-pitch general lay-ups.
* Round 5: finer / longer end layers, spacing still ≥ 1.5 hx: e24s15, e32s15 (2.5 b_s), e32s15_l4 (4 b_s).
* Round 5, first design: e24s15 brings 5cc7 from +1.532 to +0.619 % (e32s15 and e32s15_l4 are identical,
  because the spacing limit 1.5 hx is reached). Train targets paused at 77 of 1064 (21:30) until the end-layer
  setting is decided, so no labels are made with a setting that may change.
* Round 5 result (Ritz error vs label, %):

  | design | e16s15 (d2 before) | e24s15 | e32s15 | e32s15, zone 4 b_s |
  |---|---|---|---|---|
  | 5cc7be690865 (outlier) | +1.532 | **+0.619** | +0.619 | +0.619 |
  | 6cc2e9904ad4 (outlier) | +0.467 | **+0.176** | +0.176 | +0.176 |
  | 4a8f2d5f7b80 (outlier, basis-limited) | +0.266 | +0.266 | +0.266 | +0.266 |
  | 040508c6cbc1 | +0.023 | +0.023 | +0.023 | +0.023 |
  | 3265ef5ed33c | +0.113 | +0.113 | +0.113 | +0.113 |
  | 11ef6acdbbb8 | +0.116 | +0.116 | +0.116 | +0.116 |

  Finer end splines help the end-localised short-wave modes; on the normal (short) designs the 1.5 hx limit is
  already reached, so nothing changes there (same Ritz size). e32 = e24 everywhere (limit reached).
  Remaining gap of 5cc7 (0.62 % vs 0.12 % shared-basis floor) is an honest limit for now; it may close
  further with half-waves below 25 mm or harmonics, both tried only with e16 (round 3).
* **Decision: d2 end layers = e24s15** (24 splines over 2.5 b_s, spacing ≥ 1.5 hx). Cost only on panels with
  L/2 > about 2.5 b_s / 24 · 16 … i.e. the longer ones (5cc7: n 12830 → 15722, LOBPCG 100 → 190 s at 1 thread).
  The first targets (test 100 and train 77 with e16s15) are kept in
  `data/exact_targets/periodic_mixed_SS_v2_e16s15/` (logs `logs/b1_targets_*_e16s15.log`, summary
  `results/b1_targets/e16s15/`) for the comparison over all 100 test designs; test and train rerun with e24s15.
* 6 Nov: test targets with e24s15 done (100; one machine restart and one 2 h stop, resumed). Ritz error vs label:
  all mean 0.064 % (e16s15: 0.080), median 0.041, p95 0.148, max 0.619 % (was 1.535); industrial mean 0.047 /
  max 0.167 %; general mean 0.078 / max 0.619 %; none negative; 2 above 0.2 % (5cc7 0.62, 4a8f 0.27
  basis-limited). Per design e24 − e16: 31 better by > 0.005 points, 65 identical Ritz space, largest
  increase +0.005 points (LOBPCG noise at 400 iterations), largest gain −0.92. LOBPCG time ×1.03 on average
  (max ×1.60), 3.27 h in total for 100 designs at 1 thread each. Train targets running (e24s15).

## 6 Nov 2026 — B2 pilot started early (Aya: "catch anything wrong early")
* Label job switched to 1 worker (`logs/run_b1_targets.sh`, PARTS=1); the other CPU trains a pilot network on
  the designs that already have labels, so problems show up now and not after all 1064 labels.
* **Problem found by the first smoke test** (and fixed before the pilot): with the spacing limit of the d2
  end layers, splines that do not fit in half the panel are left out, so the number of functions n differs
  between designs (v2 test: 8010–15722), but the network had ONE output layout (that of the first design).
  Fix (`epinn.Layout`, `layout_labels`, `master_layout`, `design_sel`): the network now works in a fixed
  layout with every function any design can have (all 24 end splines per end; n = 16686 per mode for v2) and
  every design uses its own part (`Design.sel`, `DesignPool.sel`); all training loops, `evaluate` and
  `eval_epinn.py` select it. `prepare` checks that the predicted labels equal the dictionary's. For d1 sets
  `sel` is the identity: the step-A hybrid model reproduces its stored test values to 1e-15.
* Training can now continue after the 2 h limit or a restart (`train_ritz_ng_full(ckpt=...)`,
  `train_epinn.py --ckpt`): weights, optimiser, EMA, random generator, physics-step directions, histories and
  the validation choice are saved after every outer iteration. Smoke test (12 v2 designs, 4 outer
  iterations, stopped after 2 and continued): final weights bit-identical to an uninterrupted run.
* `train_epinn.py --files`: fixed list of operator files. `DesignPool` reads the files one at a time (it held
  all of them in memory at once, about 4 GB for 255 v2 designs).
* Pilot (`logs/run_b2_pilot.sh`, `runs/b2/pilot_hybrid.*`): hybrid, the 255 train designs labelled at 09:40
  (`runs/b2/pilot_files.txt`; last 32 = validation, 223 training), step-A settings (30 outer, 10 epochs,
  batch 16, lr 1e-3 → 1e-5 cosine, EMA 0.995, R 8, k 3), 1 thread; test = the 100 v2 test designs, looked at
  only at the end.
* Pilot first start (09:25) killed by the memory limit after reading the inputs: 4.8 GB, mostly the cache of
  16 energy spaces (v2 designs are larger). The outer loop visits all designs in order, so an LRU cache
  smaller than the set gives no reuse anyway: restarted with `--cache 2` (no speed lost).
* Note on timings (6 Nov, Aya asked): all computations and times in this log are from the 2-CPU, 8 GB
  cloud machine of the working sessions, not from a laptop. For the paper, the timing step (E) has to be
  rerun on Aya's own laptop, or the text must name the 2-core 8 GB machine.
* **Pilot result (finished 6 Nov, 30 outer iterations, about 4.5 h with 3 restarts after the 2 h limit, all
  resumed from the checkpoint):** chosen head = net of outer 24 (val head 5.58 %), modes = EMA of outer 18
  (val Ritz 31.1 %). Test (100 v2 designs, unseen): head mean 7.75 %, max 29.8 %, bias +0.03 %, 44/100
  above the GBT label; Ritz bound mean 34.3 %, max 124 %. Training designs at outer 24: head 1.22 %, Ritz
  22.9 %. **Honest reading: far worse than step A (hybrid head 1.42 %, Ritz 20.1 %), and a large gap
  between training (1.2 %) and unseen designs (5.6–7.8 %) = overfitting.** Causes to separate: (i) v2 is a
  much wider design space than step A (b_s 124–280 mm, β to 3.5, general lay-ups, n up to 15722 per mode)
  with fewer designs (223 vs 288); (ii) same network size (0.5 M parameters) and settings as step A. Next:
  certificate on the test set (running), then the label-count study (B2) will show how much more labels fix.
* **Pilot certificate** (`scripts/eval_epinn.py ... --steps 3,5,8,12`, `results/b2/pilot_hybrid_cert.log`),
  100 v2 test designs, |error| vs full GBT mean / max: head 7.75 / 29.8 %, network Ritz 34.3 / 124 %,
  3 steps 2.44 / 20.2 %, 5 steps 1.16 / 9.73 %, 8 steps 0.46 / 3.24 %, 12 steps 0.21 / 1.30 %. After 12 steps
  all 100 values are above the label (upper bounds, min +0.03 %), 10 of 100 above 0.5 %. Reading: even with
  the weak pilot network the physics steps reach the step-A level (0.20 % at 12 steps); the fast head is the
  part that needs more labels / capacity.
* Pilot done; label job back to 2 workers (PARTS=2), 420 of 1064 train labels at the switch.

## 7 Nov 2026 — B1 train labels done (1064 / 1064)
* Exact targets for all 1064 v2 train designs (1000 train + 64 val; e24s15, LOBPCG R 8, 400 iterations), no
  failures; 2 machine restarts and several 2 h limits, all resumed (atomic files, nothing partial).
  Summary `results/b1_targets/periodic_mixed_SS_v2_train.csv` (`scripts/b1_targets_summary.py --set train
  --n 1064 --seed 11`).
* Ritz error of the network's shape space vs the full-GBT label: all mean 0.061 %, median 0.040, p95 0.18,
  max 0.49 % (71f4b0f32eba, general, b_s 142); industrial (553) 0.044 / max 0.22 %; general (511) 0.079 /
  max 0.49 %; b_s < 200 mm 0.081 %, b_s 200–300 mm 0.043 %; 41 above 0.2 %, none above 0.5 %. Same picture as
  the test set (mean 0.064 %): the worst cases are again general lay-ups with small b_s.
* **Honest note: 7 of 1064 values are slightly BELOW the label** (−0.0001 to −0.020 %; 5 industrial, 2
  general, b_s 185–281 mm). The test set had none. The shared-basis × dictionary space is not a subspace of
  the full GBT discretisation (different longitudinal functions on the same Gauss grid), so it is not a
  strict upper bound of the label; the largest gap (−0.02 %) is well inside the label's own discretisation
  error found in the fine-grid check (≈ 0.05–0.1 % on short panels). Not changed, kept as they are; to be
  stated in the paper as a limit ("the Ritz value bounds the continuum problem from above only up to the
  label's discretisation error").
* Times (2-CPU 8 GB cloud machine, 1 thread): full GBT build mean 24.7 s (5–85 s), LOBPCG target mean 111 s
  (30–269 s).

## 7 Nov 2026 — B2 label study started
* Code: hybrid training can now leave designs without a GBT label (`train_epinn.py --n-labels N`: only the
  first N training designs in file order get their exact target; `epinn.train_ritz_ng_full`: exact[i] = None
  -> physics-step target and head = its Ritz loads, as in physics mode). With all designs labelled the loop is
  unchanged. Smoke test (8 train, 3 labelled, 4 val): runs, "exact targets for 3 of 8".
* File list `runs/b2/files_train1000_val64.txt`: all 1064 v2 train designs in SAMPLE order (seed 11, so the
  order is random); last 64 = validation (always labelled, used only to choose the weights); "N labels" =
  the first N of the 1000 train designs (the 64 validation designs come on top). Test = the 100 v2 test designs, not
  looked at during training.
* Launcher `logs/run_b2_labels.sh` (NL = number of labelled train designs; NL = 0 -> physics mode), step-A
  settings, 1 thread, resumable. Expected about 18–20 h per run on this machine (time per outer iteration
  grows with the number of designs: pilot 223 designs 530 s).
* Order: first the two ends (NL = 1000 all labels; NL = 0 physics only), then 50 / 100 / 200 / 400, then the
  active selection; the second run is started once the memory of the first is known (8 GB for both).
* Mistake, caught after 30 min: the first start used NL = 936 (I had counted the 64 validation designs inside
  the 1000; the log says "936 of 1000 training designs"). Stopped before its first outer iteration, log kept
  in `runs/b2/aborted/hybrid_L936.log`, restarted as NL = 1000. Memory of that start before the first outer
  iteration: 3.4 GB (inputs ready in 504 s).
* NL = 1000 run: outer 1 took 2193 s at 1 thread (training N vs full 5272 % before the first step, as in the
  pilot); memory after outer 1: 4.5 GB, only about 3 GB left. A second training run (about 4 GB) or an
  operator build (about 3 GB per worker) does NOT fit next to it, so the B2 runs go ONE after the other, each
  with both CPUs: stopped after outer 1 (checkpoint) and continued with `THREADS=2` (launcher option). The
  thread count only changes the order of floating-point sums, so outer 2 on is not bit-identical to a
  1-thread run; noted, no effect expected on the results.
* NL = 1000 progress (2 h limit every 4–5 outer iterations, always resumed from the checkpoint; about 1500 s
  per outer iteration at 2 threads). Validation (64 designs, EMA weights), head mean / max and Ritz bound mean:
  outer 4: 4.34 % / – , 24.7 %; outer 8: 4.39 / 20.9 %, 17.8 %; outer 12: 4.58 / 24.5 %, 17.0 %; outer 16:
  4.62 / 24.5 %, 18.0 %; outer 20: 4.68 / 23.6 %, 19.9 %. Training designs (16 tracked) at outer 20: head 1.28 %
  mean / 2.93 % max. **More labels help (pilot 7.75 → 4.7 %), but the head stalls from outer 8 on and the
  train/val gap stays (1.3 vs 4.7 %).**

## 7 Nov 2026 — option A: plain data surrogate for N_cr (Aya: "go ahead try option a")
* Question: is the head's stall a limit of our network, or of the design space? Test: small models that learn
  ONLY log N_cr from the GBT labels, same 38 inputs as the PINN, no physics. `scripts/a_surrogate.py`:
  train = first N of the 1000 v2 train designs (sample order, as in B2), validation = the same 64, test = the
  100 v2 test designs, errors vs full GBT. Models: MLP (3 × 128, SiLU, AdamW, cosine, early stop on val),
  5-member average (mlp5); Gaussian process (ARD RBF + white noise, sklearn); gradient boosting (gbr).
  Results `results/a_surrogate/` (run.log, results.json, test predictions at N = 1000).
* Test mean / max (%), industrial / general mean:

  | N labels | mlp5 | GP | gbr |
  |---|---|---|---|
  | 50 | 18.5 / 112 (9.3 / 26.0) | 13.6 / 76 (7.8 / 18.3) | 26.4 / 173 |
  | 100 | 10.4 / 51 (5.1 / 14.7) | 10.6 / 50 (2.4 / 17.3) | 14.9 / 62 |
  | 200 | 7.9 / 41 (2.4 / 12.4) | 6.6 / 57 (1.5 / 10.8) | 10.4 / 35 |
  | 400 | 7.8 / 37 (1.7 / 12.8) | 8.0 / 36 (0.9 / 13.7) | 8.8 / 38 |
  | 1000 | 5.0 / 29.7 (1.0 / 8.3) | 5.2 / 27.9 (0.8 / 8.7) | 6.7 / 27.0 |

  At N = 1000 by lay-up type (test 45 industrial, 55 general): GP industrial mean 0.77 %, p95 2.67 %, max 4.21 %,
  none above 5 %; general mean 8.74 %, p95 22.9 %, max 27.9 %, 35 of 55 above 5 %. Training on the 518
  industrial designs only (`results/a_surrogate/industrial_only.log`): GP industrial test 0.77 / 4.25 %, mlp5
  1.06 / 7.86 % (same as with mixed training).
* **Reading (honest):** a plain data model with all 1000 labels is no better than the PINN head (5.0–5.2 % vs
  4.7 % on validation, max about 28 % for both). So the stall is NOT a weakness of our network: the general
  ply-by-ply lay-ups are the hard part. They are random stacks from the angle set, often with few or no 0°
  plies, unbalanced and unsymmetric (worst test designs e.g. skin [60/30/60/60/60/60/30/60]). Industrial
  lay-ups (design rules) are learnt well: 0.8 % mean, 4.2 % max, like Paper 1's MLP (all within 5 %).
* Options for Aya (decision needed): (1) design space for the fast answer = industrial lay-ups (realistic),
  general lay-ups kept as a stress test of the physics certificate (0.2 % after 12 steps for any lay-up);
  (2) fast answer from a separate data head (GP / MLP ensemble, Paper 1 structure), PINN for the mode shapes
  and the certificate; (3) keep general lay-ups and use the certificate (few physics steps) for them.

## 7 Nov 2026 (afternoon) — decision: lay-up experts (Aya)
* Aya: "build different networks for every lay-up and then do integration between them … one network cannot
  get all of them at once and the inputs are too much … make it more realistic but keep the same ideas (PINN)".
  Her choices (questions answered 12:5x): **4 standard lay-ups** and **selector + physics mix**.
* Lay-up catalogue (`designspace.LAYUP_CATALOG`, skin 8 / stringer 7 plies):
  crossply = Paper 1 [0/90/0/90]s / [0/90/0/90/0/90/0] (Abaqus-validated); quasiiso [45/-45/0/90]s /
  [45/-45/0/90/0/-45/45]; zerorich [45/0/-45/0]s / [45/0/-45/0/-45/0/45]; pm45 = the Abaqus-validated "Sym1"
  [45/-45/-45/45]s / [45/-45/45/-45/45/-45/45] (used instead of a new ±45 lay-up because it has Abaqus checks
  in the periodic-GBT paper; note its stringer has 4 × 45 and 3 × −45 plies, so it is not exactly balanced).
* One expert per lay-up AND side condition (symmetry and periodic have different shared GBT bases, so they
  cannot share one output layer) → 8 small PINNs. Inputs: scaled b_s and β only (`train_epinn.py --inputs geom`;
  `epinn.INPUT_COLS`, also read by `eval_epinn.py`). Same GBT physics (ngfull), hybrid with GBT labels, same
  certificate. Sets per expert: train file 152 (seed 11, last 32 = validation), test 40 (seed 99); the same
  (b_s, β) points for every lay-up. Pipeline `logs/run_expert.sh` (build → labels → train → certificate,
  every stage resumable).
* Integration (to build after the first experts): selector = the expert of the lay-up and side condition;
  physics mix = Rayleigh–Ritz over the candidate modes of ALL experts with the GBT energy of the actual design,
  then the certificate; tested on a lay-up that no expert has seen.
* The one-network B2 run (NL = 1000) paused during outer 23 (checkpoint = outer 21; the kill came during the
  outer-22 validation step, before its checkpoint, so outer 22 (about 37 min) is lost). It stays as the
  "one network for every lay-up" baseline and can be resumed. The rest of the B2 label study on the general
  lay-up space is dropped; the label (GBT-call) study is done per expert instead (cheap: 2 inputs).
  Small fix: `run_b2_labels.sh` wrote B2_TRAIN_DONE also after a kill; removed from the log, now only on success.
* Started 13:12: expert 1 = crossply, symmetry, SS (Paper 1's case, so its 22 Abaqus points are the final check).
* Expert 1 progress: operators built 13:12–13:35 (40 test + 152 train, full GBT 8–24 s per design). GBT labels
  (LOBPCG R 8, 400 iterations, 2 workers) 13:35–15:55, one restart after the 2 h limit; all 152 labels agree with
  full GBT within 0.04 % (most within 0.01 %).
* Training first attempt (15:55) was killed by the memory limit (6 GB, `--cache 24` energy spaces kept in
  memory). `run_expert.sh` now uses `--cache ${CACHE:-8}`; restarted with cache 8 (only speed changes, not
  the result: the cache only keeps energy spaces in memory instead of reading them again).
* **Expert 1 result (crossply, symmetry, SS; 120 train designs all labelled, 32 validation, 40 unseen test;
  inputs b_s, β; 30 outer rounds, 6446 s on the 2-CPU machine).** Test, |error| vs full GBT, mean / max:
  head N 1.14 % / 4.92 % (all 40 within 5 %, bias +0.01 %); network Ritz (upper bound from the predicted
  modes) 24.9 % / 49.3 %; physics certificate 3 steps 3.79 % / 9.40 %, 5 steps 1.78 % / 5.50 %, 8 steps
  0.67 % / 2.90 %, 12 steps 0.26 % / 1.09 % (all certificate values above full GBT, as an upper bound must be).
  Reading: the fast answer is now at Paper 1's level (one lay-up, 2 inputs) — compare the one-network general
  space (4.6 % mean / 24 % max on validation). Honest limit: the predicted mode shapes are rougher than in the
  general-space run (Ritz 25 % vs 17–20 %), so the certificate needs about 8–12 steps for < 1–3 %; the worst
  designs are at the tight pitch (b_s 124–141 mm, Ritz up to +49 %). Validation had one design at 16 % head
  error during training (round 28) while all 40 test designs stay within 5 % — to look at in the integration
  stage. Hardest certificate cases: b_s 137 β 3.30 (+1.09 % after 12 steps), b_s 141 β 2.31 (+0.83 %).
* Started expert 2 (crossply, periodic, SS) 17:5x.
* Expert 2 (crossply, periodic): operators built 17:5x–18:17; labels 95/152 when the 2 h limit stopped the
  job (resumed). Physics check: for the cross-ply lay-up, full GBT gives the SAME N_cr for the periodic and the
  symmetry side condition on all 192 (b_s, β) points (largest difference 0.005 %), i.e. the critical mode is
  symmetric about the panel edges for this orthotropic (D16 = D26 = 0) lay-up. Kept the periodic expert anyway:
  its mode basis is the periodic one (needed for the periodic certificate and the selector), and two independent
  networks that must give the same load are a useful consistency check. The side conditions are expected to
  differ for lay-ups with bending–twisting coupling (quasiiso, zerorich, pm45).
* **Expert 2 result (crossply, periodic, SS)** — training resumed once after the 2 h limit (outer 15), done 22:57.
  Test (40), |error| vs full GBT, mean / max: head 1.18 % / 4.75 % (all within 5 %); network Ritz 23.4 % / 53.3 %;
  certificate 3 steps 2.67 % / 8.66 %, 5 steps 1.11 % / 4.69 %, 8 steps 0.29 % / 2.21 %, 12 steps 0.11 % / 0.79 %.
  Consistency check with expert 1 (same 40 geometries, same full-GBT load): the two independent networks'
  fast answers differ by 0.77 % mean / 2.44 % max; their errors are strongly correlated (r = 0.84), so averaging
  the two hardly helps (1.08 % / 4.84 %) — the remaining head error comes from the shared data/labels and the
  2-input map, not from random training noise.
* Disk per expert ≈ 1.0 GB operators + 0.1 GB labels + 0.1 GB checkpoint; 11 GB free → enough for the 6 left.
* Started expert 3 (quasiiso, symmetry, SS) 22:5x.
* **Expert 3 result (quasiiso, symmetry, SS)** — labels within about 0.06 % of full GBT; two restarts after the
  2 h limit (labels, training outer 20); done 08 Nov 03:45. Test (40), |error| vs full GBT, mean / max:
  head 0.55 % / 1.74 %; network Ritz 14.0 % / 45.3 %; certificate 3 steps 2.00 % / 5.42 %, 5 steps 1.18 % / 2.69 %,
  8 steps 0.60 % / 1.55 %, 12 steps 0.28 % / 1.02 %. Best expert so far for the fast answer and the mode shapes
  (Ritz 14 % vs 23–25 % for cross-ply); the 12-step certificate is similar to expert 1.
* Started expert 4 (quasiiso, periodic, SS) 03:4x.
* Expert 4 (quasiiso, periodic): operators built, 100/152 labels at the 2 h limit (resumed 06:0x). Side
  conditions for quasiiso (full GBT, same 192 points): periodic vs symmetry N_cr differ by 0.20 % mean, 1.02 % max
  (−0.23 … +1.02 %; 124 of 192 points differ by more than 0.1 %) — small but no longer zero, unlike cross-ply.

## 8 Nov 2026 — decision: no unseen lay-up
* Aya (8 Nov, 08:48): the paper covers only the lay-ups the networks were trained on (the 4 catalogue lay-ups ×
  2 sides); no test on an unseen lay-up. Integration = selector (pick the expert of the given lay-up and side)
  + physics check (certificate). Keep all data (labels, operators, checkpoints, logs, certificates) for the
  paper tables at the end.
* Honest cost note repeated for the paper: the head answer needs no GBT (<1 ms); the certificate builds the
  design's GBT matrices (no full buckling solve), which costs about as much as one GBT solve, so it is a
  check for selected designs, not the fast answer. Certificate timing still to be measured properly.

## 8 Nov 2026 — expert 4 (quasiiso / periodic) done
* Training: 30 outer rounds (two 2 h-limit restarts, resumed from checkpoint), finished 08:55; certificate 08:58.
  Test (40 designs), |error| vs full GBT, mean / max: head 0.56 / 1.49 %, network Ritz 16.0 / 47.2 %,
  certificate 3 steps 1.37 / 5.82 %, 5 steps 0.63 / 2.14 %, 8 steps 0.26 / 1.16 %, 12 steps 0.13 / 0.56 %.
  Final train head 0.80 / 4.24 %; validation (EMA, round 28) 0.84 / 7.19 %.
* Same 40 test geometries as expert 3: full GBT periodic vs symmetry differ 0.20 % mean, 0.79 % max (as on the
  training points); head errors of experts 3 and 4 correlate r = 0.59.
* Started expert 5: zerorich / symmetry (08 Nov ~09:00).

## 8 Nov 2026 — decision: level 2 = one PINN for all industrial lay-ups (Aya: "yes u should")
* Why: Option A (7 Nov) showed that many inputs are not the problem; the random ply-by-ply lay-ups are. On the
  industrial lay-ups (symmetric, balanced, 0/±45/90, >= 10 % per direction; 24 skin × 12 stringer sequences)
  a plain GP with 38 inputs already reached 0.77 % mean / 4.2 % max on the 45 industrial test designs.
* Level 1 (the 8 experts, 2 inputs each) stays the core; level 2 shows that the same GBT-PINN + certificate
  scales to about 38 inputs on realistic lay-ups. Runs after the experts (one training run fits in memory).
* Set-up (`logs/run_level2.sh`): no new GBT runs; B2 operators and labels (periodic SS v2). Train = the 518
  industrial designs of the 1000 B2 training designs, validation = the 35 industrial of the 64 B2 validation
  designs (file list `runs/level2/files_industrial.txt`), test = the 45 industrial of the 100 v2 test designs
  (symlinks in `data/operators/periodic_industrial_SS_v2/test`). Same settings as the experts but `--inputs all`.
  First run NL = 518 (all labelled); then a label study (e.g. NL = 0 / 120) for the GBT-call table.
* Limit (honest): level 2 is periodic side only (operators/labels exist only for periodic); a symmetry version
  would need about 550 new GBT builds + labels (≈ 1 day here).
* Comparison for the paper: level-2 PINN head vs the Option A GP/MLP on the same 45 industrial test designs
  (GP 0.77 / 4.21 %), plus the certificate.
* Idea noted for later (not tested): choose the GBT labels with the physics check (label the designs where
  the head and the certificate disagree most) to cut GBT calls.

## 8 Nov 2026 — expert 5 (zerorich / symmetry) done
* Build 09:0x–09:30, labels 09:30–11:45 (one 2 h-limit restart), training 11:45–13:37 (one restart at outer 21),
  certificate 13:40. Test (40 designs), |error| vs full GBT, mean / max: head 0.75 / 3.10 %, network Ritz
  19.8 / 47.2 %, certificate 3 steps 3.28 / 8.26 %, 5 steps 1.88 / 4.72 %, 8 steps 0.83 / 2.97 %,
  12 steps 0.33 / 1.09 %. Final train head 0.88 / 3.46 %; validation (EMA, end) 1.01 / 6.35 %.
* The certificate converges slower than for quasiiso (3 steps 3.3 % vs 1.4–2.0 %), as for crossply; after
  12 steps all five experts are at 0.1–0.3 % mean, about 1 % max.
* Disk: 7.0 GB free. Started expert 6: zerorich / periodic (08 Nov ~13:45).

## 8 Nov 2026 — expert 6 (zerorich / periodic) done
* Labels done 16:30 (one 2 h-limit restart), training 16:30–18:26 (one restart at outer 19), certificate 18:29.
  Test (40 designs), |error| vs full GBT, mean / max: head 0.67 / 2.14 %, network Ritz 18.0 / 47.1 %,
  certificate 3 steps 1.45 / 4.79 %, 5 steps 0.60 / 2.46 %, 8 steps 0.21 / 0.97 %, 12 steps 0.11 / 0.38 %.
  Final train head 0.88 / 3.54 %; validation (EMA, end) 0.78 / 4.71 %.
* Same 40 test geometries as expert 5: full GBT periodic vs symmetry differ 0.28 % mean, 1.82 % max (25 of 40
  above 0.1 %), so the side condition matters for zerorich as for quasiiso (crossply: identical).
* Disk: 5.8 GB free (enough for experts 7 and 8, about 1.2 GB each). Started expert 7: pm45 / symmetry (18:3x).

## 8 Nov 2026 — expert 7 (pm45 / symmetry) done
* Build 18:3x–18:56, labels 18:56–21:24 (one 2 h-limit restart), training 21:24–23:16 (one restart at outer 18),
  certificate 23:19. Test (40 designs), |error| vs full GBT, mean / max: head 0.24 / 0.87 %, network Ritz
  15.9 / 49.3 %, certificate 3 steps 1.21 / 5.23 %, 5 steps 0.48 / 2.10 %, 8 steps 0.21 / 0.70 %,
  12 steps 0.14 / 0.44 %. Final train head 0.22 / 0.42 %; validation (EMA, end) 0.26 / 1.78 %.
* Best expert so far (reason not checked; possibly a smoother N_cr over (b_s, β) for this lay-up).
* Disk: 4.6 GB free (expert 8 needs about 1.2 GB). Started expert 8: pm45 / periodic (23:2x).

## 8 Nov 2026 (23:37) — GBT-call study planned (Aya asked: "like the one in the last paper")
* Paper 1 reference: 144 panels: 84 % fast PINN, 10.4 % re-solved PINN, 5.6 % Abaqus; 3.31 h vs all-Abaqus 1.64 h.
* Paper 2 plan, after all levels: (1) training cost = GBT labels per network (label study NL = 0 / 10 / 20 /
  40 / 120 per expert; 0 / 120 / 518 for level 2); (2) use cost = a dense sweep per network (e.g. 1000 points of
  (b_s, β), as in Paper 1) through a tiered check: network answer → cheap trigger (seed-ensemble spread, as
  Paper 1's CoV) → physics check (3–12 steps) for flagged designs → full GBT only if the check does not settle.
  Report % per tier, total time vs all-GBT, and accuracy vs full GBT on every point.
* Honest point to keep: the physics check costs about one GBT solve, so it must not run on every design,
  otherwise the time is no better than all-GBT (Paper 1's hybrid was also slower than all-Abaqus). The trigger
  decides the saving; its miss rate vs full GBT is measured on all sweep points.

## 9 Nov 2026 — expert 8 (pm45 / periodic) done: all 8 experts finished
* Labels 23:44–02:08 (one 2 h-limit restart), training 02:08–04:04 (restart at outer 19), certificate 04:08.
  Test (40 designs), |error| vs full GBT, mean / max: head 0.29 / 0.97 %, network Ritz 21.2 / 76.7 %,
  certificate 3 steps 1.07 / 4.79 %, 5 steps 0.40 / 2.17 %, 8 steps 0.19 / 0.67 %, 12 steps 0.14 / 0.44 %.
  Final train head 0.26 / 0.85 %; validation (EMA, end) 0.25 / 0.81 %.
* Same 40 test geometries as expert 7: full GBT periodic vs symmetry differ only 0.05 % mean, 0.12 % max
  (2 of 40 above 0.1 %) for pm45 (crossply: identical; quasiiso 0.20 / 0.79 %; zerorich 0.28 / 1.82 %).
* All 8 experts, test head mean / max: crossply 1.14 / 4.92, 1.18 / 4.75; quasiiso 0.55 / 1.74, 0.56 / 1.49;
  zerorich 0.75 / 3.10, 0.67 / 2.14; pm45 0.24 / 0.87, 0.29 / 0.97 % (symmetry, periodic). Certificate
  12 steps: 0.11–0.33 % mean, 0.38–1.09 % max.

## 9 Nov 2026 — level 2 started; out of memory with cache 8
* Level 2 (`NL=518 bash logs/run_level2.sh`) started 04:08 via `logs/run_night.sh`; outer rounds take about
  775 s (553 designs, 2 threads). Killed during outer 3/4 (no Python error, "Killed" = memory limit, as with
  expert 1 at cache 24): the 553-design pool with `--cache 8` is too large for the about 6 GB limit (B2 used
  cache 2 for 1064 designs). Resumed from the checkpoint with `CACHE=2` (only speed changes, not results).
  Expected: about 27 more rounds × 13 min ≈ 6 h plus 2 h-limit restarts.

## 9 Nov 2026 09:00 — level 2 at outer 19: validation got worse after round 8 (overfitting)
* Restarts at the 2 h limit: during outer 12 (06:55) and during outer 20 (09:00; checkpoint after outer 19).
* Validation (35 industrial designs, EMA head) mean / max by evaluation round: 2: 5.41 / 20.15 %,
  4: 1.89 / 8.14, 6: 1.61 / 6.82, **8: 1.60 / 4.91 (best)**, 10: 1.80 / 4.70, 12: 2.00 / 6.08,
  14: 2.24 / 5.53, 16: 2.23 / 5.69, 18: 2.45 / 7.37 %. Train head (16 designs) keeps improving
  (round 18: 1.10 / 4.61 %), so the network is overfitting the 518 training designs.
* `train_epinn.py` keeps the best-on-validation weights (now round 8 EMA), so the final model is not hurt by
  the later rounds; it can only change if validation gets below 1.60 % in rounds 20–30.
* Honest note: 1.60 % (validation) is about twice the GP of Option A (0.77 / 4.21 % on the 45 industrial test
  designs; different set, so the fair comparison is the test run after training). Not restarted at 09:00:
  waiting for Aya to choose (continue rounds 20–30, about 2.5 h, or finish now with the round-8 weights).

## 9 Nov 2026 09:45 — level 2 finished with the round-8 weights; label study started; GP baseline for the experts
* Aya: finish level 2 with the best-on-validation weights (round 8) and start the next step. `run_level2.sh` got
  `OUTER` (default 30); `OUTER=19` continues the checkpoint of outer 19 with no more rounds, runs the final
  evaluation, keeps the best-on-validation weights and saves. Then its certificate on the 45 industrial test designs.
* Chain `logs/run_label_study.sh` (resumable): level 2 finish → expert crossply / symmetry with NL = 0, 10, 20,
  40 labels of 120 → NL = 0 (physics only) for the other 7 experts. `run_expert.sh`: NL = 0 → `--mode physics`,
  tag `<side>_<lay>_SS_v2_physics_L0`.
* Aya asked why labels matter if the 1064-label run "did not work": that run (B2, NL = 1000) was on the MIXED
  space (half random ply-by-ply lay-ups); there even a GP with 1000 labels had 5.2 % mean / 28 % max, so it says
  nothing about labels on realistic lay-ups. More labels did help there (physics-only pilot 7.75 % → 4.7 %).
* **New baseline (`scripts/expert_tables.py`, `results/expert_tables/gp_baseline.csv`):** Gaussian process (same
  kernel as option A, log N, inputs b_s and β) on the SAME first NL training designs as the PINN, errors on the
  same 40 test designs vs full GBT, mean / max (%):

  | expert | GP 10 | GP 20 | GP 40 | GP 120 | PINN head 120 |
  |---|---|---|---|---|---|
  | crossply / sym | 3.32 / 18.30 | 2.29 / 14.90 | 1.68 / 7.68 | 0.48 / 3.65 | 1.14 / 4.92 |
  | quasiiso / sym | 0.98 / 5.05 | 0.58 / 3.26 | 0.60 / 3.10 | 0.42 / 2.25 | 0.55 / 1.74 |
  | zerorich / sym | 1.59 / 6.13 | 1.22 / 5.30 | 0.72 / 4.04 | 0.66 / 2.93 | 0.75 / 3.10 |
  | pm45 / sym | 1.61 / 14.29 | 0.49 / 5.99 | 0.07 / 0.50 | 0.01 / 0.17 | 0.24 / 0.87 |
  | crossply / per | 3.32 / 18.30 | 2.29 / 14.90 | 1.31 / 5.18 | 0.48 / 3.65 | 1.18 / 4.75 |
  | quasiiso / per | 1.06 / 4.89 | 0.62 / 3.48 | 0.52 / 2.87 | 0.46 / 2.47 | 0.56 / 1.49 |
  | zerorich / per | 1.53 / 5.63 | 1.95 / 11.43 | 0.67 / 3.24 | 0.69 / 2.77 | 0.67 / 2.14 |
  | pm45 / per | 1.62 / 14.30 | 0.48 / 5.95 | 0.07 / 0.51 | 0.01 / 0.17 | 0.24 / 0.97 |

* **Reading (honest):** with all 120 labels and only 2 inputs, a plain GP gives the critical load as well as or
  better than the PINN head (better mean on 6 of 8; PINN better max on quasiiso and zerorich / periodic; GP far
  better on pm45, 0.01 vs 0.24 %). So the paper cannot claim that the PINN head beats a data model with many
  labels. What only the PINN gives: the buckling mode, the physics check (the GP answer has no check), and
  training without labels or with few labels. With 10–20 labels the GP max error is 5–18 %, so the label study
  (NL = 0 / 10 / 20 / 40) decides whether the physics helps there. A few-label win (if it comes) and the
  physics check are the honest case for the PINN.
* The same script writes every expert's GBT numbers (152 training + 40 test designs: b_s, β, full-GBT N, label
  Ritz loads, label time) to `results/expert_tables/*.csv`, kept in git, so the numbers survive if the big label
  and operator files are deleted.
* Literature plan written: `docs/LITERATURE_PLAN.md` (6 topics, one comparison table).
* Disk: 2.9 GB free (data 18 GB). No new operators are needed for the label study.

## 9 Nov 2026 09:52 — level 2 result (industrial lay-ups, 38 inputs, NL 518, round-8 weights)
* Chosen on validation: EMA weights of evaluation 8 (val head 1.60 %); modes: EMA of evaluation 14.
* Test (45 industrial designs, unseen), |error| vs full GBT, mean / max: **head 1.35 / 8.40 %** (bias −0.23 %,
  25/45 above), network Ritz 16.86 / 34.68 %, certificate 3 steps 1.50 / 3.98 %, 5 steps 0.70 / 2.79 %,
  8 steps 0.31 / 2.13 %, **12 steps 0.14 / 0.74 %** (`results/level2/industrial_hybrid_L518_cert.log`).
* Honest comparison on the SAME 45 designs: GP of option A (1000 mixed labels, 518 industrial among them)
  0.77 / 4.21 %; GP on the 518 industrial only 0.77 / 4.25 %. The PINN head is worse (1.35 / 8.40 %); after
  the physics check (12 steps) it is better than the GP (0.14 / 0.74 %), but the check costs about one GBT
  matrix build per design.
* Old data sets (before the v2 space) exported to small CSV files in `results/old_sets/` (labels: sid, full-GBT
  N, Ritz loads, label time; operators: b_s, β, lay-ups, full-GBT N, GBT time), so their numbers stay in git if
  the big files are deleted. `periodic_industrial_SS/oracle` holds only broken links (train set deleted 2 Nov).
* Label study chain running: crossply / symmetry physics-only (NL 0) started 09:52; round 2 head 55 % (the
  physics-only start is slow, as in the pilot).

## 9 Nov 2026 10:10 — deleted old files (Aya: "ok delete them"); certificate cost corrected
* Deleted (numbers kept in git, see 09:52): the 8 expert restart files `runs/experts/*_hybrid_L120_ckpt.pt`
  (finished, final models in git), `runs/b2/pilot_hybrid_ckpt.pt`, and the pre-v2 sets
  `data/operators/{periodic_industrial_SS,symmetry_industrial}`, `data/exact_targets/{periodic_industrial_SS,
  periodic_mixed_SS_v2_e16s15,smoke_v2}`. Free disk 2.9 → 4.9 GB. Kept: `runs/b2/hybrid_L1000_ckpt.pt` (D, the
  only way to continue the one-network baseline; Aya's choice pending), all v2 expert / mixed operators and labels.
* Correction (honest): I wrote "the certificate costs about one GBT solve". From the logs it is less, but not
  measured cleanly: the certificate reads the design's GBT matrices in the network's Ritz space (building them
  without the full GBT eigen solve took about 2.5 s per design, 5 Nov) and then needs about 0.3 s per step
  (`eval_epinn.py` docstring), so about 2.5 + 12 × 0.3 ≈ 6 s with 12 steps, vs full GBT 8–24 s per expert design
  (mean 18 s on the mixed set). Rough estimate ≈ 3× cheaper than full GBT; to be timed cleanly (nothing else
  running) with a split into matrix build and steps.
* Plan to cut the check cost: (1) check only flagged designs (tiered, as Paper 1: network → ensemble trigger →
  check → full GBT); (2) fewer steps where 5 are enough; (3) profile the 2.5 s build and see whether the
  matrices can be prepared once per expert (offline) so that each design only needs a cheap online step
  (as in reduced-basis methods; literature topic 2).
* Aya (10:19): add point 3 (prepare the check once per expert, cheap online step per design) to the plan; done in
  `docs/STATUS.md` and `docs/LITERATURE_PLAN.md` (cheaper physics check: clean timing, profile the build,
  offline / online split, fewer steps). Starts after the label study (the CPUs are busy with training now).

### 9 Nov 2026, 11:40: label study, interim (crossply / symmetry, NL 0 = physics only)

* The chain was stopped by the 2 h limit at round 27 of 30 (checkpoint after round 27). Restarted
  `logs/run_label_study.sh` (Aya: "keep running"); level 2 is skipped, the expert resumes at round 28.
* Physics-only progress (validation, 32 designs, vs full GBT): round 6 head 35.1 / 73.9 %, round 24 13.6 / 35.5 %,
  round 26 13.4 / 34.8 % (mean / max). Every head value is **above** full GBT (bias +13 %, 32/32 above). The Ritz
  upper bound is 26 / 56 %.
* Honest reading so far: without labels the network learns slowly and stays far from full GBT; GP with 10 labels
  is 3.3 / 18.3 % and the hybrid PINN with 120 labels 1.14 / 4.92 %. The final test numbers (and the certificate,
  which corrects the network with physics steps) follow when the run ends. Possible reason for the bias: the
  physics loss only pushes the Ritz quotient down step by step (an upper bound), and 30 rounds are not enough
  to converge; to check after the run (not a change to the run).

### 9 Nov 2026, 11:55: label study, NL 0 finished (crossply / symmetry, physics only)

40 test designs, |error| vs full GBT, mean / max (`results/experts/symmetry_crossply_SS_v2_physics_L0_cert.log`):

| | head | network Ritz | check 3 steps | 5 steps | 8 steps | 12 steps |
|---|---|---|---|---|---|---|
| NL 0 (physics only) | 11.93 / 24.09 | 23.73 / 59.32 | 4.09 / 11.97 | 1.92 / 6.45 | 0.73 / 2.87 | 0.28 / 1.39 |
| NL 120 (hybrid) | 1.14 / 4.92 | 24.90 / 49.33 | 3.79 / 9.40 | 1.78 / 5.50 | 0.67 / 2.90 | 0.26 / 1.09 |
| GP, 10 labels | 3.32 / 18.30 | | | | | |
| GP, 40 labels | 1.68 / 7.68 | | | | | |

* Honest reading:
  1. The network's own number (head) is poor without labels: 11.9 % mean, always above full GBT. A GP with
     10 labels is better. So for the fast number alone, labels are needed.
  2. The physics check after the network is almost as good without labels as with 120 labels (5 steps 1.92 vs
     1.78 %, 12 steps 0.28 vs 0.26 %). The mode shapes, which the check uses, come from the physics; the labels
     mainly correct the head. NL 0 + 5-step check (1.92 / 6.45 %) is about as good as GP with 40 labels, with no
     labelled design, but it costs the check time per design (about 6 s, not yet timed cleanly).
* Deleted the NL 0 restart file (`_ckpt.pt`, 103 MB); the trained network and log are kept. NL 10 is running.

### 9 Nov 2026, 12:50: Navid: "something is wrong with the paper formulation, do not use it for the ML"

* Aya reports that Navid says the periodic-GBT formulation should not be used for the ML, although its results
  are correct. His exact reason is not known yet. Candidates (my reading of `manuscript/part2.md`, not
  confirmed):
  - the mesh changes with b_s (segments at most b_s/25: 92 nodes at b_s = 150 mm, 78 at 285.5 mm), so the
    number and meaning of the modes jump across the design space; a network needs a smooth map from design to
    answer;
  - the GBT modes themselves are not smooth in the design (unit-maximum scaling, family-by-family projection,
    hierarchical re-ordering at one reference wavenumber, mode crossings);
  - the link walls use factors 10^3 / 10^-3 (penalty-like), and thin walls with Mindlin and transverse-extension
    modes make K very ill-conditioned, which hurts gradient training;
  - the computed pre-buckling state makes G depend on a linear solve for every design;
  - the cell of width 2 b_s holds only Bloch waves 0 and pi (bay-repeating and alternating modes).
* Our ML already avoids the first two (shared POD basis), and the labels are full-GBT values, so the trained
  numbers are not affected by mode scaling. Decision: ask Navid for his exact point; nothing is stopped.
* 12:47 Aya: Navid means the periodic-BC formulation of GBT itself. Most likely point: the condition ties all six
  DOFs of the side edges (u(-W) = u(0)), so the cell width is fixed and the Poisson expansion of the skin is
  restrained. The paper says this (Section 3.2: S_ss ≈ (A12/A11) S_xx, 0.03-0.04 n0 cross-ply, 0.32 n0
  isotropic, 0.86-0.89 n0 angle-ply; Sections 3.5 and 4: "free lateral expansion" is future work). The solver
  has free expansion (constant jumps of U1, U3) only in the strip solver (`OMEGA_PERIODIC_MACRO`, raises an error
  in GBT). The "symmetry" side condition (U1 = 0 on both sides) also fixes the width. Abaqus used the same
  constraint, which is why the results agree. All Paper 2 labels use the fixed width. Proposed check: strip
  solver, fixed vs free width, a few designs per lay-up at beta <= 3.5. Waiting for Navid's exact point.

### 9 Nov 2026, 13:58: label study, NL 10 finished (crossply / symmetry)

40 test designs, mean / max vs full GBT: head 10.45 / 26.14 % (39/40 above), network Ritz 22.34 / 51.95 %,
check 3 steps 4.02 / 11.73 %, 5 steps 1.91 / 6.23 %, 8 steps 0.71 / 2.88 %, 12 steps 0.28 / 1.42 %.
(NL 0: head 11.93 / 24.09 %; NL 120: 1.14 / 4.92 %; GP with the same 10 labels: 3.32 / 18.30 %.)

* Honest reading: 10 labels barely improve the head (11.9 -> 10.5 %); the GP with the same 10 labels is 3x better.
  In the hybrid loss, the head of the 110 unlabelled designs is fitted to the physics step, whose value is still
  far above full GBT (upper bound, slow convergence), and this pulls the head up; 10 labels cannot correct it.
  The check is unchanged (mode shapes come from the physics).
* Not changed during the study (no tuning). Idea for later, to test on validation only: fit the head of
  unlabelled designs to the refined check value, or give the labelled designs more weight.
* Deleted the NL 10 restart file. NL 20 is running.

### 9 Nov 2026, 16:00: label study, NL 20 finished (crossply / symmetry)

40 test designs, mean / max vs full GBT: head 9.34 / 26.31 % (40/40 above), network Ritz 22.38 / 51.28 %,
check 3 steps 3.86 / 11.25 %, 5 steps 1.90 / 7.29 %, 8 steps 0.69 / 2.65 %, 12 steps 0.27 / 1.24 %.
GP with the same 20 labels: 2.29 / 14.90 %. Head so far: NL 0 11.93, NL 10 10.45, NL 20 9.34, NL 120 1.14 % mean.
The head improves only slowly with labels and stays biased high; the check is the same for every NL.
Deleted the NL 20 restart file. NL 40 is running.

### 9 Nov 2026, 18:05: label study, NL 40 finished; crossply / symmetry part complete

40 test designs, mean / max vs full GBT (head = network number; check after 5 / 12 physics steps):

| NL | head | check 5 | check 12 | GP head, same labels |
|---|---|---|---|---|
| 0 | 11.93 / 24.09 | 1.92 / 6.45 | 0.28 / 1.39 | - |
| 10 | 10.45 / 26.14 | 1.91 / 6.23 | 0.28 / 1.42 | 3.32 / 18.30 |
| 20 | 9.34 / 26.31 | 1.90 / 7.29 | 0.27 / 1.24 | 2.29 / 14.90 |
| 40 | 6.95 / 13.06 | 1.78 / 5.62 | 0.26 / 1.22 | 1.68 / 7.68 |
| 120 | 1.14 / 4.92 | 1.78 / 5.50 | 0.26 / 1.09 | 0.48 / 3.65 |

* Honest conclusion: for the number alone, the GP beats the hybrid PINN at every label count (about 4x at
  NL 10-40, 2x at NL 120). The PINN head is always biased high when few labels are used. The check (shape from
  the network + GBT steps) is independent of the labels: with 0 labels, 5 steps match GP-40 in mean and are better
  in max; 12 steps (0.28 / 1.4 %) beat GP-120 in max (3.65 %).
* Needed before writing: (1) the check started without the network (shared basis / nearest training design) to
  show what the network adds; (2) clean timing of the check; (3) the head fix idea (validation only).
* Deleted the NL 40 restart file. NL 0 for the other 7 experts is running.

### 9 Nov 2026, 19:30: what does the network add to the check? Start baselines (crossply / symmetry)

`scripts/check_start_baselines.py` (new): the same check (block LOBPCG, energy metric, R = 8, seed 1, same 40 test
designs) started from (a) random vectors, (b) the GBT modes of the nearest labelled training design in scaled
(b_s, beta) among the first 10 or 120, copied into the test design's coordinates, (c) the trained networks.
`results/check_baselines/symmetry_crossply_SS_v2.csv` (every step for every design).

| start | 3 steps | 5 steps | 8 steps | 12 steps | 20 steps | 30 steps |
|---|---|---|---|---|---|---|
| random (no labels, no network) | 11.10 / 27.77 | 4.25 / 11.78 | 1.39 / 6.19 | 0.39 / 1.83 | 0.09 / 0.40 | 0.04 / 0.18 |
| network, physics only (L0) | 4.09 / 11.97 | 1.92 / 6.45 | 0.73 / 2.87 | 0.28 / 1.39 | 0.09 / 0.37 | 0.05 / 0.20 |
| network, hybrid 120 labels | 3.79 / 9.40 | 1.78 / 5.50 | 0.67 / 2.90 | 0.26 / 1.09 | 0.08 / 0.32 | 0.04 / 0.19 |
| nearest of 10 labelled designs | 2.81 / 8.97 | 1.18 / 4.16 | 0.41 / 1.53 | 0.14 / 0.54 | 0.05 / 0.19 | 0.03 / 0.13 |
| nearest of 120 labelled designs | 1.28 / 2.70 | 0.51 / 1.72 | 0.21 / 1.08 | 0.08 / 0.32 | 0.03 / 0.12 | 0.02 / 0.10 |

(mean / max |error| % vs full GBT.)

* Honest reading: the network start saves a few steps against a random start (5 steps: 1.9 vs 4.3 %), but copying
  the modes of the nearest of only 10 labelled designs is better than every network at every step, and the
  nearest of 120 is about 3x better. As trained now, the network's mode shapes are not good enough to justify it
  for the check (its own Ritz value is 22-25 % above full GBT, i.e. the predicted modes are far from the
  eigenvectors). Together with the GP result for the number, the current PINN does not beat simple baselines.
* Mistake (mine): the baseline ran beside the training, the memory limit was hit and the training of periodic
  crossply NL 0 was killed at round 6 (resumed from the round-5 checkpoint, cache 3 instead of 8; the cache only
  changes speed).
* Formulation check (Aya: "are all our formulations correct?"): `manuscript/verify_math.py` passes again (paper's
  modal matrices = code to 1e-15 for all lay-ups, every mode periodic at the side nodes, classical limit -0.001 %).
  I re-derived the equilibrium operators (Eq. Kop) and geometric tensors by integration by parts: signs agree.
  Physics point: the periodic condition fixes the cell width. Standard periodic conditions allow a constant jump
  (macroscopic strain, Xia, Zhang and Ellyin 2003), and Bloch-wave buckling theory (Geymonat, Mueller and
  Triantafyllidis 1993; Triantafyllidis and Schnaidt 1993) computes the pre-buckling state under the macroscopic
  load (uniaxial compression = zero mean transverse force). The paper's condition is uniaxial load with zero mean
  transverse strain. Wording point: the Ritz value bounds the Ritz-space load, not full GBT (crossply: 70 of 152
  labels are up to 0.04 % below full GBT).

### 9 Nov 2026, 19:55: why the network start is weak; combined starts; Aya: stop the label study, fix the shapes

* Aya (19:45): stop the label study after the running run (periodic crossply NL 0), use the machine for clean
  timing and better mode shapes. The chain loop was stopped; the running training continues alone.
* `scripts/check_mode_space.py`: Ritz value of the start itself (no steps), crossply / symmetry, 40 test designs,
  mean / max: network L0 23.7 / 59.3 %, L120 24.9 / 49.3 %; nearest of 10 labelled 62.7 / 221.8 %, of 120
  19.2 / 90.1 %; 2 / 4 nearest of 120 together 11.4 / 38.8 % and 6.8 / 15.6 %; POD of the 960 training label
  modes (whitened coordinates) with 8 / 16 / 32 / 48 vectors 67.3 / 35.2 / 17.9 / 13.3 % (POD keeps only 71 % of
  the energy with 64 vectors, so a common POD basis of label modes is not compact). (A first try with 128 / 256
  POD vectors ran out of memory; only my script was killed.)
* Combined starts (`check_start_baselines.py`, same 40 designs), mean / max after 5 / 12 steps: 4 nearest of 10
  1.04 / 5.30 and 0.17 / 1.70 %; 2 nearest of 120 0.50 / 1.84 and 0.08 / 0.33 %; nearest of 10 + network L0
  1.46 / 5.83 and 0.21 / 0.85 %; nearest of 10 + network L120 1.32 / 5.18 and 0.19 / 0.80 %.
* Reading: the network's modes have a good Ritz value (24 %, better than a copied neighbour: 63 %) but the check
  converges more slowly from them; adding them to the neighbour's modes makes the check worse, because the first
  Rayleigh-Ritz step keeps the network's vectors (lower Ritz value) and drops useful neighbour vectors. The network
  learns smooth shapes with a low energy quotient, not the detail the check needs; its training loss (fit to one
  physics step / Ritz value) is not the quantity that matters for the check. And with labels, it does not even
  fit the training modes (Ritz 19 % on labelled training designs).
* Running: the same baselines for the level-2 network (38 inputs, `scripts/check_start_baselines_l2.py`), where
  nearest neighbours are far apart and a network may help more.

### 9 Nov 2026, 20:20: start baselines for the level-2 network (38 inputs, 45 industrial test designs)

`scripts/check_start_baselines_l2.py`; nearest labelled design in the network's standardised 38 inputs.
`results/check_baselines/level2_industrial_hybrid_L518.csv`. Mean / max |error| % vs full GBT:

| start | 3 steps | 5 steps | 8 steps | 12 steps | 30 steps |
|---|---|---|---|---|---|
| random | 5.84 / 20.76 | 1.90 / 7.11 | 0.57 / 3.39 | 0.24 / 2.10 | 0.07 / 0.21 |
| network (hybrid, 518 labels) | 1.50 / 3.98 | 0.70 / 2.79 | 0.31 / 2.13 | 0.14 / 0.74 | 0.06 / 0.19 |
| nearest of 50 labelled | 1.80 / 6.09 | 0.62 / 2.10 | 0.23 / 1.15 | 0.11 / 0.48 | 0.06 / 0.20 |
| nearest of 518 labelled | 1.59 / 5.21 | 0.59 / 2.23 | 0.21 / 1.15 | 0.10 / 0.48 | 0.06 / 0.20 |
| 2 nearest of 518 | 1.28 / 4.28 | 0.49 / 2.10 | 0.21 / 1.27 | 0.11 / 0.54 | 0.06 / 0.23 |

* With 38 inputs the network is about as good as copying a neighbour: best worst case at 3 steps (3.98 %), a bit
  better mean than one neighbour at 3 steps, but worse from 5 steps on. No clear win either way. The floor at
  30 steps (0.06 %) is the Ritz space vs full GBT.

### 9 Nov 2026, 20:40: neighbour + correction network (nnres); timing script; night chain

* `epinn.NeighbourResidualPINN` (`--arch nnres`): output = GBT modes of the nearest labelled training design
  (bank of the labelled designs, fixed layout, whitened) + a FunctionDecoderPINN correction whose last layer starts
  at zero; head = neighbour's log N + correction. Trained exactly as before (ngfull, hybrid: physics step on every
  design + label modes on labelled ones), so it stays a PINN. A design never copies itself (exact matches are
  skipped), so labelled training designs are predicted from their neighbour, as test designs are. Unit test:
  untrained output = neighbour copy, reload from checkpoint identical. `epinn.net_from_ckpt` now rebuilds every
  network for evaluation (bank size from the weights). `logs/run_expert_nnres.sh` = run_expert.sh + ARCH
  (run_expert.sh was left unchanged because the running job was still reading it).
* `scripts/time_check.py`: clean timing of full GBT, the check's parts (operators of a new design, energy space,
  steps), the network and the GP. A first test on ONE design while training was running (not clean, not kept):
  full GBT 13.5 s; check = operators 5.6 s + energy space 3.4 s + 0.27 s per step, i.e. 10.7 s for 5 steps.
  Warning: if the clean numbers are similar, the check is NOT much cheaper than full GBT (my earlier "about 6 s
  vs 18 s" was too optimistic); the operators and the energy space dominate, not the steps. This makes the
  planned offline / online split (point 3 of the plan) essential.
* Night chain `logs/run_night_9nov.sh`: wait for periodic crossply NL 0, then clean timing (10 designs), then
  nnres crossply / symmetry NL 10 and NL 120.

### 9 Nov 2026, 22:30: clean timing; periodic crossply NL 0; nnres L10 half-way

* **Clean timing** (`results/timing/check_timing.txt`, 10 crossply / symmetry test designs, 2 CPU cores, 7.8 GB RAM,
  nothing else running), mean per design:

  | step | time |
  |---|---|
  | full GBT (reference N_cr) | 9.76 s |
  | check: operators of the new design | 4.28 s |
  | check: energy space + metric | 1.92 s |
  | check: one step | 0.17 s (first 0.27 s) |
  | check with 5 steps, total | 7.14 s |
  | check with 12 steps, total | 8.31 s |
  | network forward pass | 2.7 ms |
  | GP prediction | 0.5 ms (fit 0.2 s) |

  Honest result: the check as it is now saves only 15–27 % against full GBT (7.1–8.3 s vs 9.8 s). 87 % of the
  check time is building the new design's operators and energy space, not the steps. The network and GP are
  ~3000x faster than full GBT but give no bound. A real speed-up of the check needs the offline / online split
  (operators interpolated or assembled from precomputed parts per expert); this is now a must before writing.
* **Periodic crossply, physics only (NL 0)** (`results/experts/periodic_crossply_SS_v2_physics_L0_cert.log`, 40 test
  designs, mean / max): head 11.52 / 24.72 %, network Ritz 23.85 / 55.65 %, check 3 steps 3.01 / 10.46 %,
  5 steps 1.14 / 5.04 %, 8 steps 0.35 / 2.52 %, 12 steps 0.14 / 1.04 %. Same picture as symmetry crossply NL 0
  (head 11.93 / 24.09, s5 1.92 / 6.45, s12 0.28 / 1.39); the periodic check converges a bit faster.
* **nnres crossply / symmetry NL 10, round 13 of 30** (killed by the 2 h limit, resumed from ckpt): validation
  Ritz of the network falls slowly from 58.4 % (round 2, ~the neighbour copy) to 52.3 % (round 12); head 20 % mean
  (worse than the plain hybrid L10 head 10.45 % and far worse than the GP, 3.32 %). The head result is expected:
  with 10 labels the nearest labelled design is far, and its log N is a poor base. Check results at the end.

### 10 Nov 2026, 00:30: nnres crossply / symmetry NL 10 finished

`results/experts/symmetry_crossply_SS_v2_nnres_L10_cert.log`, 40 test designs, |error| vs full GBT, mean / max:

| start of the check | head | 3 steps | 5 steps | 8 steps | 12 steps |
|---|---|---|---|---|---|
| plain network, hybrid L10 | 10.45 / 26.14 | – | 1.91 / 6.23 | – | 0.28 / 1.42 |
| nearest of 10 labelled (copy, no network) | – | 2.81 / 8.97 | 1.18 / 4.16 | 0.41 / 1.53 | 0.14 / 0.54 |
| nnres L10 (copy + PINN correction) | 18.06 / 53.26 | 2.90 / 10.70 | 1.26 / 5.60 | 0.44 / 2.12 | 0.15 / 0.72 |
| GP on 10 labels (number only) | 3.32 / 18.30 | | | | |

* nnres makes the check start much better than the plain network (5 steps 1.26 vs 1.91 %, 12 steps 0.15 vs 0.28 %),
  but only because it starts from the neighbour copy: the PINN correction makes it slightly WORSE than the copy
  alone (all columns). Network Ritz 47 % mean (copy start ~ 58 % on validation at round 2, falling to 47 %): the
  physics step lowers the Rayleigh quotient but not in the directions the check needs.
* Head 18 % mean: worse than the plain network and much worse than the GP. With 10 labels the neighbour's log N is
  a poor base.
* Honest conclusion for NL 10: no gain from the PINN over the copy. NL 120 running (round 13/30 at 00:27).

### 10 Nov 2026, 01:45: nnres crossply / symmetry NL 120 finished (night chain done)

Two container restarts and one 2 h kill during training; resumed from the checkpoint each time (no loss).
`results/experts/symmetry_crossply_SS_v2_nnres_L120_cert.log`, 40 test designs, |error| vs full GBT, mean / max:

| start of the check | head | 3 steps | 5 steps | 8 steps | 12 steps |
|---|---|---|---|---|---|
| plain network, hybrid L120 | 1.14 / 4.92 | 3.79 / 9.40 | 1.78 / 5.50 | 0.67 / 2.90 | 0.26 / 1.09 |
| nearest of 120 labelled (copy, no network) | – | 1.28 / 2.70 | 0.51 / 1.72 | 0.21 / 1.08 | 0.08 / 0.32 |
| nnres L120 (copy + PINN correction) | 5.62 / 12.69 | 1.30 / 2.71 | 0.51 / 1.71 | 0.21 / 1.07 | 0.08 / 0.32 |
| GP on 120 labels (number only) | 0.48 / 3.65 | | | | |

* nnres = the neighbour copy, to two decimals. The physics-trained correction hardly moves the modes (start Ritz
  19.6 % vs 19.2 % for the copy; physics fit 2.517 -> 2.503 over 30 rounds). So in the 2-input experts the
  "PINN" part adds nothing measurable to the check start; all the gain comes from the labels.
* nnres head 5.6 % mean: worse than the plain network head (1.1 %) and the GP (0.5 %).
* Conclusion (both NL 10 and NL 120): for the lay-up experts the best check start is a plain neighbour copy and the
  best number is the GP. The current PINN beats neither. Options for the paper are listed in STATUS for Aya.

### 10 Nov 2026, 08:00: Aya's decision: plan B (fast certified check = main contribution, level-2 PINN inside it)

Aya (07:57): "Go with what u suggest bec level 2 is the competitive part and also fast check is important".

**Where the check time goes** (profile, crossply / symmetry):
* operators of a new design: cross-section 0.2 s, **pre-buckling state 2.4–9.7 s** (full-GBT static solve:
  assembly of the fine spline system + PARDISO factorisation + one solve), shared section 0.2 s, K, G, P 0.1 s.
* energy space 0.1 s; **energy metric 0.14–2.35 s** (Cholesky of the per-wavelength blocks, the end-layer block up
  to 4969 wide). `time_check.py` built the metric twice (once in `prepare`, once as preconditioner).
* Full GBT = the same assembly + factorisation + pre-buckling solve + eigen solve; so the old check only replaced the
  eigen solve and could not be much faster.

**Tests (not kept as results, scratchpad scripts):**
* Pre-buckling solved in the shared basis (241 modes) instead of full GBT (570 modes): SLOWER (4.1–9.2 s, the shared
  basis is dense) although N_cr changes only 0.0002–0.007 %. Rejected.
* Pre-buckling on a coarser longitudinal mesh (same solver, hx x k), resultants interpolated (cubic) to the fine Gauss
  points; N_cr = full GBT eigen solve with this state imposed, 5 crossply designs:

  | k | time of pre-buckling | dN_cr vs full GBT |
  |---|---|---|
  | 1 | 2.4–6.3 s | 0 |
  | 2 | 1.5–3.4 s | −0.004 … −0.001 % |
  | 4 | 1.1–1.8 s | −0.007 … −0.001 % |
  | 8 | 1.1–1.4 s | −0.002 … +0.023 % |
  | 16 | 0.8–1.2 s | −0.017 … +0.032 % |

  Chosen k = 4 (< 0.01 %, ~1/3 of the time). The check then needs no factorisation of the fine system at all.
  Note: this is a (tiny) change of the reference problem; the check's values are upper bounds of the Ritz-space load
  WITH the coarse state, and errors are always reported vs full GBT with the solver's own state.

**Code** (pinnbuck only, solver untouched; defaults unchanged so all earlier data stay valid):
`families._Omega.pre(coarsen)`, `build_operator(..., pre_coarsen=1)` (stored in `op.info`), `epinn.prepare` keeps
the metric as `d.em` (used as the preconditioner, built once). New `scripts/fast_check.py` (per design: full GBT
time, fast check time and error from the nn120 copy and from the expert network after 3 / 5 / 8 / 12 steps) and
`logs/run_fast_check.sh` (all 8 experts, 40 test designs each). Smoke test, 2 designs: full GBT 9.8 s; check with
nn120 start 4.1 s (5 steps, 0.12 % error) / 5.4 s (12 steps, 0.03 %).
