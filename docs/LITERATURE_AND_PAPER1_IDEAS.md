# Ideas from Paper 1 and the literature: used, not needed, still to try (3 Nov 2026)

Literature list as recorded in *Paper2_status_report_simple.pdf* (§6) and *Paper2_PINN_details.pdf* (§12);
the full scan notes of October are not in the repo.

## From Paper 1 (Eq. 18 and §2.5)
| Paper 1 element | role there | Paper 2 status |
|---|---|---|
| Rayleigh quotient λ_cr as main loss | physics | **used**, generalised to Rayleigh–Ritz over 8 candidates (Ky Fan sum of 3 modes) |
| interface penalties α (ΔW, ΔW', ΔM, ΔV), α doubling every 100 epochs | continuity of sub-plates | **not needed**: GBT modes are continuous by construction |
| Loss_sym (zero slope at the symmetry line) | boundary condition | **not needed**: side conditions built into the modes (symmetry / periodic) |
| η (RMS_W − 1)² amplitude normaliser | avoid w = 0 | **not needed**: Rayleigh–Ritz is scale-free |
| γ·SmoothL1(λ/N_FEA, 1), γ = 1e5 | active-learning data term | **tried 27 Oct** with GBT oracle: fixes the flagged designs, no transfer to unseen ones |
| Adam (lr 5e-3, 800 epochs) → L-BFGS (strong Wolfe) | two-stage optimiser | ngfull → Adam polish tried (Adam degrades the head); **L-BFGS in whitened coordinates not yet tried** (the 7.8 % L-BFGS failure of 22 Oct was in raw coordinates) |
| transfer initialisation from the neighbour panel (median 5.8 → 1.7 %) | warm start | analogues: certificate starts from the prediction (done); **curriculum over lay-ups / geometry not yet tried** |
| deep ensemble E = 10, CoV > 5 % | epistemic uncertainty | ensemble exists only in the teacher network (`net.py`); **not yet in the PINN**. With GBT the model-form error is ≈ 0.5 %, so the CoV becomes meaningful (in Paper 1 it missed 35 % vs Abaqus) |
| −15 %/+3 % neighbour window | trigger | **replaced** by the physics certificate (Ritz upper bound + residual) and the head–Ritz gap |
| separate MLP surrogate (b_s, β) | fast tool | replaced by the PINN itself (< 1 ms, lay-up as input) |

## From the literature
| idea | source | status |
|---|---|---|
| energy natural gradient (metric = energy Hessian) | Müller & Zeinhofer, ICML 2023 | **used** (energy whitening + ngfull) |
| curvature-aware / second-order optimisers for PINNs | Jnini et al. 2026 (incl. Karniadakis) | to try (with L-BFGS polish) |
| parameter encoder + coordinate encoder | Cho et al., P2INN, ICML 2024 | **used** (FunctionDecoderPINN) |
| n-width barrier of linear reduced bases | Lee & Carlberg, JCP 2020 | diagnosis that led to the decoder |
| predict a slightly larger subspace + Rayleigh–Ritz (30 % → 2 % in their tests) | Fanaskov et al. 2025; Li, Sun, Zhang 2025 | partly: 8 candidates for 3 loads; **R = 16 to try** |
| mixture of experts / domain decomposition | POU-PINNs (arXiv 2412.06842); APINN (arXiv 2211.08939) | **to try**: the Rayleigh–Ritz layer selects candidates itself |
| spectral inference networks (trace / Ky Fan loss) | Pfau et al. | used (Ky Fan loss) |
| close and crossing eigenvalues: snapshots with several modes; GP strong with little data | Boffi et al. 2024; Bertrand, Boffi, Halim | used (3 modes, subspace loss); GP was the baseline (9.3 %) |
| physics-reinforced reduced models (physics residual + projected data) | Chen, Wang, Hesthaven, Zhang | tried 25 Oct: no gain at 160 designs |
| multi-fidelity (many cheap + few accurate analyses) | Pollet et al. 2026; Yasar & Gürses 2024 | not tried (strip or closed form as low fidelity) |
| block preconditioned eigen-iteration | Knyazev, LOBPCG 2001 | **used** (physics step, certificate) |

## Diagnosis that decides the order (27 Oct)
Plateau ≈ 3.6 % (head) also on the training designs → the network cannot yet represent the modes well
enough (capacity), not a lack of data. Ideas that add capacity come first.
