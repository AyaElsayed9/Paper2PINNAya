# Status of Paper 2 (machine-learning part) — 3 November 2026

Read this first when continuing the work. Every try and decision: `docs/RESEARCH_LOG.md`; plan: `docs/ROADMAP.md`; code: `src/pinnbuck/`.

## Goal
A physics-informed neural network (PINN) for buckling of the omega-stringer panel (later the member
sections), with the shear-deformable GBT energy as its physics. It fixes the two weak points of Paper 1
(Elsayed et al., Compos. Struct. 399 (2026) 120959): the reduced bending stiffness (RBS) approximation
of the plate energy, and the expensive Abaqus data. **Decision (2 Nov 2026): Paper 2 uses the GBT code,
not Abaqus.** (The README still mentions Abaqus for blind validation; the ROADMAP carries the decision note.)

## Done
| step | result |
|---|---|
| Phase 0 | frozen solver V32, package `pinnbuck`, design space, tests |
| Shared mode basis (POD) | 241 GBT modes for every omega panel, error vs full GBT 0.01–0.02 % |
| Operators | K, G and pre-buckling state of a panel in the shared basis, checked to 1e-13 |
| Functions along the length | dictionary (6 wavelengths, sin/cos + 8 end splines per end): 0.05–0.27 % |
| First physics-only training | too slow (close eigenvalues, end-localised modes, condition number 1e9–1e19) |
| Data-teacher detour | 160 GBT-labelled designs: about 9–16 % error, dropped |
| Energy whitening (key step) | energy ≈ \|y\|² in the new coordinates; eigen-iteration 0.1 % in 7–13 steps (was 150) |
| Network | 38 inputs (geometry, end/side flags, lamination parameters, closed-form loads, wave content) → 8 candidate modes + N_cr head; function decoder, 0.53 M parameters; Rayleigh–Ritz layer |
| Optimiser | `ngfull`: energy natural gradient with a Rayleigh–Ritz line search |
| Physics certificate | LOBPCG started from the predicted modes (replaces the −15 %/+3 % window of Paper 1) |

## Results (errors vs full GBT on unseen designs, no labels in training)
| case | head (instant) | Ritz of predicted modes | certificate, 12 steps |
|---|---|---|---|
| CC, symmetry: 259 training → 22 unseen | 5.0 % mean / 10.9 % max | 12.3 / 25.9 % | 0.54 / 1.72 % |
| SS, periodic: 320 training → 32 unseen | 3.65 / 9.77 % | 10.7 / 34.2 % | 0.22 / 1.24 % |

Last experiment (27 Oct): Paper 1 data term (SmoothL1, oracle on the 10 % most suspicious designs).
It corrects the designs it is applied to but does not transfer to unseen designs. Plateau: about
3.6 % (head) and 10 % (Ritz), also on training designs, so the limit is the network's representation
of the modes (capacity), not data or optimiser steps.

## Latest result (3 Nov 2026)
GBT at the 22 Abaqus points of Paper 1: 0.54 % mean / 0.97 % max (Paper 1 PINN: 2.53 % / 10.80 %),
tight pitch 124.25 mm: 0.89 % instead of 8.08 % mean; all 22 half-wave numbers right. The GBT engine
removes the RBS error. Paper 1 facts and Aya's plan: `docs/PAPER1_SUMMARY.md`.

## Error analysis (3 Nov 2026, afternoon)
The best network's head is biased high (27 of 32 unseen designs above GBT, mean +3.13 %; scatter after
removing the bias 2.79 % mean, 6.64 % max). Exact GBT targets are already stored in every operator file.
Revised plan: see the last entry of `docs/RESEARCH_LOG.md` (hybrid loss with exact GBT targets first).

## Three-way test (3 Nov 2026, night) — unseen designs, one seed each
| | head mean / max, bias | Ritz bound | certificate 3 / 12 steps |
|---|---|---|---|
| physics-only | 3.65 / 9.77 %, +3.13 % | 10.7 % | 1.78 / 0.22 % |
| data-only | 2.78 / 6.25 %, −1.71 % | 18.9 % | 0.80 / 0.14 % |
| hybrid | 3.94 / 11.02 %, −1.08 % | 14.1 % | 0.79 / 0.13 % |
Labels remove the bias; physics gives the best shapes; head noisy between checkpoints (±1.5 %);
average of the 3 heads 1.61 / 5.52 %. Decision 4 Nov: hybrid is the method, quality first; plan A–E in `docs/RESEARCH_LOG.md`.
## Step A (4 Nov 2026): clean training (cosine LR, EMA, choice on 32 validation designs), 288 training designs
| arm | head mean / max, bias | Ritz bound | certificate 3 / 12 steps |
|---|---|---|---|
| physics-only | 4.48 / 11.89 %, +4.47 % | 11.95 % | 1.92 / 0.22 % |
| data-only | 1.21 / 4.79 %, +0.05 % | 22.8 % | 1.10 / 0.14 % |
| hybrid | 1.42 / 4.35 %, +0.55 % | 20.1 % | 1.81 / 0.20 % |
Labels fix the head but it overfits after outer 8–10; the physics gives the best modes. Fixes: separate weight
choice for head and modes, more designs (B1), seeds (D).
Running: step B1. v2 operator sets built (periodic SS, mixed lay-ups, full GBT; test 100, train 1000 + 64 val; no failures); operator files now stored compressed (`opio`, lossless). B0 done: shared bases kept, dictionary d2 = 11 half-waves 25–1006 mm, labels with 400 iterations. End layers decided: d2 = 11 half-waves + e24s15 (24 splines over 2.5 b_s, spacing ≥ 1.5 hx; was e16s15); 24 test designs mean/max Ritz error 0.053/0.116 % (industrial), 0.078/0.141 % (general), none below the label. Exact targets (`--iters 400`): test 100 done (Ritz error mean 0.080 %, 3 general lay-ups with small b_s above 0.2 %, worst 1.53 %, being probed); end layers revised to e24s15 (worst test design 1.53 → 0.62 %, normal designs unchanged); test targets done (mean 0.064 %, max 0.62 %, none worse than e16s15); train targets done (1064, mean 0.061 % / max 0.49 %, 7 tiny values below the label, ≥ −0.02 %); B2 pilot (hybrid, 223 + 32 val labelled designs) finished: test head mean 7.75 % / max 29.8 %, Ritz bound 34 % — much worse than step A and overfits (train head 1.2 %); certificate 12 steps mean 0.21 % / max 1.30 % (all upper bounds), same level as step A. Network now uses a fixed output layout (n differs between d2 designs) and resumable training (`--ckpt`).
Running (7 Nov): step B2 label study on the 1000 v2 train designs (physics on all, GBT labels on N = 1000 / 0 / 50 / 100 / 200 / 400, then active selection; `logs/run_b2_labels.sh`); first run NL = 1000 (all labels), outer 21 of 30: validation head stalls at about 4.7 % mean / 24 % max (train 1.3 %).
Option A (7 Nov): plain data model for N_cr (GP / MLP, same 38 inputs, 1000 labels): test 5.0–5.2 % mean / 28 % max = same as the PINN head, so the limit is the general random lay-ups, not the network. Industrial lay-ups: 0.77 % mean / 4.2 % max (GP); general: 8.7 % / 28 %. Decision (Aya, 7 Nov): **lay-up experts** — one small PINN per realistic lay-up (crossply = Paper 1, quasiiso, zerorich, pm45) and side condition, inputs b_s and β only, joined by a selector + physics mix (Ritz over all experts' modes). One-network run paused at outer 21 (baseline). Expert 1 crossply / symmetry / SS done: test (40) head 1.14 % mean / 4.92 % max, certificate 12 steps 0.26 % / 1.09 %. Expert 2 crossply / periodic done: head 1.18 % / 4.75 %, certificate 12 steps 0.11 % / 0.79 % (cross-ply: periodic N_cr = symmetry N_cr). Expert 3 quasiiso / symmetry done: head 0.55 % / 1.74 %, certificate 12 steps 0.28 % / 1.02 %. Expert 4 quasiiso / periodic done: head 0.56 % / 1.49 %, certificate 12 steps 0.13 % / 0.56 %. Expert 5 zerorich / symmetry done: head 0.75 % / 3.10 %, certificate 12 steps 0.33 % / 1.09 %. Expert 6 zerorich / periodic done: head 0.67 % / 2.14 %, certificate 12 steps 0.11 % / 0.38 %. Expert 7 pm45 / symmetry done: head 0.24 % / 0.87 %, certificate 12 steps 0.14 % / 0.44 %. Expert 8 pm45 / periodic done: head 0.29 % / 0.97 %, certificate 12 steps 0.14 % / 0.44 %. ALL 8 EXPERTS DONE (9 Nov 04:08). Level 2 (industrial lay-ups, NL 518): stopped at outer 19 (validation best at round 8, 1.60 / 4.91 %, then overfitting); Aya 9 Nov: finished with the round-8 weights: test (45) head 1.35 / 8.40 % (GP 0.77 / 4.21 %), certificate 12 steps 0.14 / 0.74 %. Running (9 Nov 09:45, `logs/run_label_study.sh`): label study, crossply / symmetry NL 0, 10, 20, 40, then NL 0 for the other 7 experts. Result NL 0 (crossply / symmetry, 40 test designs, mean / max vs full GBT): head 11.93 / 24.09 % (always above) vs 1.14 / 4.92 % with 120 labels; but the physics check after the network is almost the same: 5 steps 1.92 / 6.45 % (L120: 1.78 / 5.50), 12 steps 0.28 / 1.39 % (L120: 0.26 / 1.09). So labels fix the number the network outputs; the mode shapes (used by the check) are learnt from physics alone. NL 10: head 10.45 / 26.14 % (GP 10 labels: 3.32 / 18.30 %), check 5 steps 1.91 / 6.23 %, 12 steps 0.28 / 1.42 %: 10 labels barely move the head (the 110 unlabelled designs dominate it). NL 20: head 9.34 / 26.31 % (GP 20: 2.29 / 14.90 %), check 5 steps 1.90 / 7.29 %, 12 steps 0.27 / 1.24 %. NL 40: head 6.95 / 13.06 % (GP 40: 1.68 / 7.68 %), check 5 steps 1.78 / 5.62 %, 12 steps 0.26 / 1.22 %. Crossply part done; running: NL 0 for the other 7 experts. **Start baselines (9 Nov evening, results/check_baselines/):** the same physics check started without the network: from the GBT modes of the nearest labelled training design (nn10, nn120) or at random. Crossply / symmetry, 40 test designs, mean / max after 5 / 12 steps: random 4.25 / 11.78 and 0.39 / 1.83 %; nn10 1.18 / 4.16 and 0.14 / 0.54 %; nn120 0.51 / 1.72 and 0.08 / 0.32 %; network L0 1.92 / 6.45 and 0.28 / 1.39 %; network L120 1.78 / 5.50 and 0.26 / 1.09 %. So the network start is better than random but worse than copying the nearest of only 10 labelled designs: as trained now, the network adds little to the check. Formulation check: paper equations = code (verify_math.py, 1e-15), all modes periodic, classical limit -0.001 %; operators re-derived by parts, signs agree. Open physics point: fixed cell width (see log). GP baseline (same labels, 2 inputs): with 120 labels the GP is as good as or better than the PINN head (e.g. crossply / sym 0.48 / 3.65 vs 1.14 / 4.92 %); case for the PINN = few labels (to be shown) + modes + physics check. Literature plan: `docs/LITERATURE_PLAN.md`. Decision 8 Nov: no unseen lay-up; only the 4 trained lay-ups; keep all data. Decision 8 Nov: level 2 = one PINN for all industrial lay-ups (38 inputs, periodic SS, reuses B2 data), after the experts (`logs/run_level2.sh`). Planned after all levels: GBT-call study as Paper 1 (training labels + 1000-point sweep with tiered check: network → ensemble trigger → physics check → full GBT). Planned (Aya 9 Nov, "add it to the plan"): **cheaper physics check** — (a) time the check cleanly (matrix build vs steps); (b) profile the about 2.5 s GBT matrix build; (c) prepare what does not depend on the design once per expert (offline) so that each new design only needs a cheap online step, as in certified reduced-basis methods (literature topic 2); (d) fewer steps where 5 are enough. Goal: the check much cheaper than full GBT for every design, not only for flagged ones; risk: b_s and β change the geometry, so it may only partly work.
Next: design space as Paper 1 (b_s from 124.25 mm, β to 3.5, symmetry + periodic), ~1000 designs, GBT-call
count study. All runs indexed in `results/EXPERIMENTS.md` (`scripts/collect_results.py`).

## Where we stopped — open questions
1. Break the plateau: more candidates (Ritz mixture of experts), experts per mode family or geometry range,
   larger decoder, longer training on the laptop.
2. The projection step (network following the physics step) is slow: fit loss 0.98 → 0.30 in 20 steps.
3. Cost: the certificate needs the operators of the design (5–10 s), about the cost of one GBT solve
   (about 10 s). The instant head is the fast answer; the certificate is a check, not a speed-up.
4. Validation without Abaqus: GBT itself is validated against Abaqus in the periodic-GBT paper.
5. Later phases: surrogate and uncertainty (3), Koiter imperfection sensitivity in GBT (4), figures and manuscript (5).
