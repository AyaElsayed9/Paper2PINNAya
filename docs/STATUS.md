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
Running: step B1, building the v2 operator sets (periodic SS, mixed lay-ups, full GBT; test 100, train 1000 + 64 val), `logs/run_b1_build.sh`. B0 done: shared bases kept (max 0.14 %), dictionary d2 = 11 half-waves 25–1006 mm, labels with 400 iterations. Open: general lay-ups may need longer end layers (one smoke design 0.36 %); decided from the exact targets of the B1 sets.
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
