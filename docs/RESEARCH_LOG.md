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
