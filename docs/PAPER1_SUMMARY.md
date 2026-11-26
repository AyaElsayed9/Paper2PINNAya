# Paper 1 in short (what Paper 2 builds on)

Elsayed, Kharghani, El Yaakoubi-Mesbah, Mittelstedt. *A hybrid physics-informed neural network framework for
local buckling prediction of stiffened composite panels with epistemic uncertainty quantification.*
Composite Structures 399 (2026) 120959. (PDF uploaded by Aya: `1-s2_0-S0263822326009244-main.pdf`.)

| part | Paper 1 |
|---|---|
| panel | omega-stringer panel, symmetry line (half panel), SS loaded edges |
| lay-up | **one**: skin [0/90/0/90]s, stringer [0/90/0/90/0/90/0] (cross-ply) |
| geometry | b_s 124.25–285.5 mm, β = a/b_s up to 3.5; b3 = 27.75 mm, b1 = 0.1 (b_s − 2 b3) |
| physics engine | 1D deep energy method (CLPT, Kirchhoff), two sub-plates (foot, skin bay) |
| coupling | foot is unsymmetric → **RBS**: D* = D − B A⁻¹ B (Eq. 3) |
| training | one PINN **per panel**: Adam → L-BFGS, 800 epochs, transfer initialisation from the neighbour |
| logic check (Eq. 16) | accept if 0.85·N_prev ≤ N ≤ 1.03·N_prev (previous panel of the sweep) |
| if violated | deep ensemble (E = 10), CoV; CoV > 5 % → Abaqus (≈ 41 s per panel) |
| data term | γ·SmoothL1(N_cr/N_FEA, 1) when Abaqus was called |
| labels of 144 panels | 84 % fast PINN, 10.4 % re-solved PINN, 5.6 % Abaqus |
| accuracy vs Abaqus | 16 regular/transition points: mean 1.34 %, max 3.54 %; tight pitch b_s = 124.25: **10.3–10.8 %** (RBS / 1D limit); Table 1 (beta sweep, Fig. 6): b_s 124.25 mean 10.30 %, **max 14.54 % at beta 0.58**; b_s 285.5 max 7.11 % at beta 0.58 (Paper 1's overall largest error is 14.54 %, not 10.80 %) |
| ensemble check | missed triggers 11 % vs closed form, **35 % vs Abaqus** (cannot see the RBS model error) |
| cost | 3.31 h for 144 panels (all-Abaqus 1.64 h) |
| fast tool | separate MLP surrogate (b_s, β → N_cr), 0.818 ms, R² 0.9998, all held-out points within 5 % |

## Aya's plan for Paper 2 (3 Nov 2026, her words summarised)
1. Same idea as Paper 1, but the physics engine is GBT → removes the RBS error.
2. The GBT code replaces Abaqus as the solver/oracle.
3. Replace the −15 %/+3 % neighbour logic check by a reasonable one.
4. Lay-up as an input (generalisation over lay-ups), not one lay-up.
5. The tool must be fast.

## Post-buckling in Paper 1 (added 17 Nov 2026; was missing from this summary)
Paper 1 also ran Abaqus Riks: the PINN's predicted mode 1 was mapped onto the Abaqus skin mesh as a mode-shaped
imperfection (scale 0.01, `generate_abaqus_matched_imperfection`, `build_abaqus_riks_case2.py` in the notebook) and the
load / end-shortening path was compared with the Riks path seeded by Abaqus's own eigenmode (Case 1 m = 3; Case 2
b_s = 205 mm; figures Fig_Riks_Case1, Fig_11_Riks_Case2). So Paper 1's mode shapes were used for post-buckling.
