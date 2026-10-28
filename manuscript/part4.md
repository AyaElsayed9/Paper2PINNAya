
# Verification and discussion

## Reference models

*Abaqus shell models.* The reference solutions are linear buckling analyses in Abaqus with S4R and S3 shell elements and a seed of 3.7 mm. The stringer is tied to the skin with \*TIE and the true offset of the mid-surfaces, and the loads are \*SHELL EDGE LOAD of 1 N/mm (2 N/mm on the feet, Section 2.2) at both loaded edges. SS ends prescribe $U_2 = 0$ on every end edge of skin and stringer; CC ends prescribe $U_1 = U_2 = UR_1 = UR_2 = UR_3 = 0$ at both ends and $U_3 = 0$ at $x = 0$. Symmetry side edges prescribe $U_1 = UR_2 = UR_3 = 0$. For periodic side edges, every skin node of the edge $X = -W$ is paired with the node of $X = 0$ at the same $x$, and for each pair one constraint equation (\*EQUATION) $\mathcal U_\alpha(-W,x) - \mathcal U_\alpha(0,x) = 0$ is written for each of the six degrees of freedom, as in Eq. (@pbc); at the two end corners the equation for $U_2$ is omitted because $U_2 = 0$ is already prescribed there. The rigid-body motions left by periodicity are removed by $U_1 = U_3 = 0$ at the corner $X = 0$, $x = 0$. The Abaqus database covers 60 cases with symmetry side edges (five lay-ups, SS and CC, two pitches, three aspect ratios) and the 30 SS cases with periodic side edges. Periodic CC cases have not been computed in Abaqus.

The GBT model imposes the same physical conditions (Table 4 and Section 3.6). Two differences of form do not affect the eigenvalue: Abaqus ties the nodes of its mesh, whereas the GBT modes satisfy the periodic condition at every $x$; and the drilling rotation $UR_2$, which Abaqus ties, is not a degree of freedom of the plate theory (in S4R it carries only a small artificial stiffness).

*Spline finite strip model.* As a second, independent reference, in particular for the periodic CC cases, an in-house spline finite strip model {{cheung1976,lau1986}} is used. It discretises the same geometry, laminates, loads and end conditions with flat laminated Mindlin strips of width $b_s/50$ (twice as fine as the GBT segments), linear interpolation across each strip, a mid-strip transverse shear strain against shear locking, and six degrees of freedom per nodal line in global axes. Along $x$ it uses the same cubic B-splines as GBT. The web roots are joined to the skin by rigid offset links (master–slave constraints), the direct analogue of the Abaqus \*TIE, and the side conditions are applied node by node, exactly as in Abaqus: the six degrees of freedom of the two side nodal lines are merged for every spline coefficient. Pre-buckling state and geometric stiffness are computed as in Sections 3.7 and 3.8. The strip model shares the plate theory with GBT but none of the cross-section discretisation: it has no deformation modes, a finer mesh, rigid links instead of link walls, and boundary conditions on nodal degrees of freedom instead of built into modes. Its agreement with GBT is therefore a genuine check of the mode construction. For $b_s = 150$ mm and $L = 225$ mm it has 37 170 unknowns (177 nodal lines × 6 × 35 splines), against 18 777 for GBT.

## Exactness of the periodic condition

Since the condition is built into the modes, every buckling mode satisfies it identically. For the Cross-ply cell with $b_s = 150$ mm, the largest difference between the side nodes over all conventional modes is $7.9\times10^{-15}$ for the translations, $1.1\times10^{-16}$ for the rotation $\theta$ and $2.3\times10^{-11}$ for the warping (unit-normalised modes). In the buckling mode of the Unsym panel with CC ends, $b_s = 150$ mm and $\beta = 1.5$, the largest side-edge differences along the whole length are $1.2\times10^{-18}$ ($U_X$), $6.9\times10^{-18}$ ($U_Y$), $1.4\times10^{-20}$ ($\theta$) and $1.1\times10^{-15}$ (warping), with modal amplitudes of order $10^{-5}$ to $10^{-1}$, while the side edges do move ($U_X, \theta \neq 0$). With symmetry side edges the same panel has $U_X = \theta = 0$ at the sides and different deflections of the two side edges.

## Symmetry side edges

Table 5 compares the critical loads for the 60 cases with symmetry side edges; the individual results are listed in Appendix A. GBT agrees with Abaqus within 0.58 % on average and 1.48 % at most, and the spline finite strip model within 0.57 % and 1.30 %. GBT and strip differ by 0.27 % on average. The accuracy is the same for SS and CC ends and does not deteriorate for the coupled lay-ups, which confirms the treatment of the laminate couplings in Eqs. (@Bmat)–(@H1mat) and the clamped end condition of Section 3.6.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 5.** Absolute deviation of the critical load in %, 60 cases with symmetry side edges (six cases per row: $b_s = 150, 285.5$ mm, $\beta = 0.75, 1.5, 3$).
:::

| Lay-up | Ends | GBT vs Abaqus, mean | max | Strip vs Abaqus, mean | max | GBT vs strip, mean | max |
|:--|:--|--:|--:|--:|--:|--:|--:|
| Isotropic | SS | 0.56 | 1.08 | 0.59 | 1.10 | 0.12 | 0.28 |
| Isotropic | CC | 0.58 | 1.08 | 0.55 | 1.05 | 0.13 | 0.28 |
| Cross-ply | SS | 0.49 | 0.91 | 0.41 | 0.69 | 0.22 | 0.35 |
| Cross-ply | CC | 0.50 | 0.92 | 0.42 | 0.87 | 0.14 | 0.25 |
| Sym1 | SS | 0.52 | 0.90 | 0.65 | 1.30 | 0.31 | 0.52 |
| Sym1 | CC | 0.89 | 1.48 | 0.57 | 1.07 | 0.67 | 1.14 |
| Sym2 | SS | 0.64 | 1.08 | 0.58 | 1.19 | 0.27 | 0.43 |
| Sym2 | CC | 0.67 | 1.23 | 0.63 | 1.23 | 0.41 | 0.74 |
| Unsym | SS | 0.46 | 0.66 | 0.66 | 1.15 | 0.25 | 0.49 |
| Unsym | CC | 0.48 | 0.97 | 0.65 | 1.21 | 0.17 | 0.33 |
| All 60 | | 0.58 | 1.48 | 0.57 | 1.30 | 0.27 | 1.14 |

The deviations depend mainly on the stringer pitch (Figure 3a). For $b_s = 285.5$ mm the GBT mean deviation is 0.25 %. For $b_s = 150$ mm it is 0.91 %, and all 30 GBT results and all 30 strip results lie below Abaqus (GBT −1.48 % to −0.48 %, strip −1.30 % to −0.65 %). This offset is therefore common to the two plate models and not caused by the GBT mode basis; the refinement study of Section 4.6 shows that it is not a discretisation error of GBT either, since finer segments lower the GBT result slightly. Its likely origins are differences between the shell elements and the plate model in the short skin bays of the narrow pitch, and the representation of the tie; it has not been resolved in this study. The largest GBT deviation, −1.48 %, occurs for Sym1, CC, $b_s = 150$ mm, $\beta = 3$ (GBT 12.871, strip 12.972, Abaqus 13.064 N/mm).

![**Figure 3.** Critical loads of GBT against Abaqus. (a) Deviation for the 60 cases with symmetry side edges (open symbols: SS, filled symbols: CC). (b) Deviation for the 30 SS cases with periodic side edges. (c) Change of the critical load caused by periodic instead of symmetry side edges, $100\,(N_{cr}^{per}-N_{cr}^{sym})/N_{cr}^{sym}$, GBT against Abaqus (30 SS cases; the line marks equality).](fig/fig3_validation.png){width=15.5cm}

## Periodic side edges

The 30 SS cases with periodic side edges are summarised in Table 6 and listed in full in Table A1. GBT agrees with Abaqus within 0.57 % on average and 1.25 % at most (strip: 0.59 % and 1.37 %). As for the symmetry edges, the mean deviation is larger for $b_s = 150$ mm (0.93 %, maximum 1.25 %) than for $b_s = 285.5$ mm (0.20 %, maximum 0.46 %). The periodic formulation is therefore as accurate as the symmetry formulation (Figure 3b).

::: {custom-style="MDPI_4.1_table_caption"}
**Table 6.** Periodic side edges, SS ends, 30 cases: absolute deviation of the critical load in %, and change of the critical load from symmetry to periodic side edges, $\Delta = 100\,(N_{cr}^{per}-N_{cr}^{sym})/N_{cr}^{sym}$, range over the six cases of each lay-up.
:::

| Lay-up | GBT vs Abaqus, mean | max | Strip vs Abaqus, mean | max | $\Delta$ Abaqus [%] | $\Delta$ GBT [%] | $\Delta$ strip [%] |
|:--|--:|--:|--:|--:|:-:|:-:|:-:|
| Isotropic | 0.56 | 1.08 | 0.59 | 1.10 | 0.00 | 0.00 | 0.00 |
| Cross-ply | 0.49 | 0.91 | 0.41 | 0.69 | 0.00 | 0.00 | 0.00 |
| Sym1 | 0.53 | 0.92 | 0.65 | 1.30 | −0.10 to +0.07 | −0.10 to +0.06 | −0.10 to +0.06 |
| Sym2 | 0.73 | 1.25 | 0.60 | 1.25 | −0.17 to +1.35 | −0.17 to +1.16 | −0.19 to +1.23 |
| Unsym | 0.53 | 0.92 | 0.71 | 1.37 | +0.19 to +2.06 | +0.14 to +2.01 | +0.15 to +2.10 |
| All 30 | 0.57 | 1.25 | 0.59 | 1.37 | | | |

A more demanding test than the absolute critical load is the *change* $\Delta$ that periodic instead of symmetry side edges produce, because the two models of each case share their discretisation and the common offset of Section 4.3 largely cancels. GBT reproduces the change found by Abaqus with a mean deviation of 0.05 percentage points (maximum 0.27 points) and a correlation coefficient of 0.995 (Figure 3c); the strip model gives 0.03 points (maximum 0.23) and 0.996.

## Effect of periodic side edges

For the isotropic and the cross-ply panel, periodic and symmetry side edges give the same critical load in all three models (Table 6 and Tables A1 and A2; the largest difference is 0.01 %). The reason is that the critical mode of an orthotropic panel without coupling is symmetric about the mid-lines of the skin bays, so that the symmetry conditions $U_X = \theta = 0$ are satisfied by the exact periodic mode. For coupled laminates this is no longer the case. A symmetry plane at the side edge implies that the neighbouring cell is the mirror image of the cell, and mirroring changes the sign of the ply angles: a cell with the angle-ply laminates Sym1 or Sym2 is mirrored into a cell with the opposite bending–twisting coupling. The symmetry model therefore represents a panel that is not the real one, while the periodic model represents the real, infinitely repeated panel. Bending–twisting and membrane–bending coupling skew the nodal lines of the buckles (Figure 5), and the symmetry conditions force the skewed mode back to zero rotation at the sides.

![**Figure 4.** Change of the critical load from symmetry to periodic side edges, $100\,(N_{cr}^{per}-N_{cr}^{sym})/N_{cr}^{sym}$, by lay-up: range of GBT for SS and CC ends (bars, six cases each) and individual results of the spline finite strip model (crosses).](fig/fig4_shift_ss_cc.png){width=10cm}

The size of the effect follows the strength of the coupling (Figure 4). For Sym1 (±45°) the change remains within ±0.1 %, although its bending–twisting coupling is not zero, which indicates that its critical mode is nearly symmetric about the bay mid-lines. For Sym2 (60°/−30°) and Unsym it reaches 1.35 % and 2.06 % in Abaqus with SS ends (GBT: 1.16 % and 2.01 %). With CC ends GBT predicts much larger changes, +5.50 % to +18.16 % for Sym2 and +8.53 % to +11.61 % for Unsym, and the strip model confirms them independently (+5.57 % to +17.97 % and +8.40 % to +11.50 %; GBT and strip differ by 0.30 % on average and 1.17 % at most over the 30 periodic CC cases, Table A2). A plausible explanation is that clamping removes the end rotations through which the skewed buckle of the symmetry model can adjust near the loaded edges, so that the difference between the periodic panel and its mirrored substitute grows; a detailed explanation requires further study. These CC results have not yet been confirmed by Abaqus. Since symmetry side edges always gave the lower value for the strongly coupled lay-ups (with one exception of −0.17 % for Sym2), a symmetry model is conservative for these panels, but by up to 18 % for clamped ends.

## Discretisation

The default discretisation (segments of at most $b_s/25$, spline spans $h_x = b_s/20$) is checked by refining each parameter in turn (Table 7). All changes are below 0.5 %, i.e. below the scatter against Abaqus. The slightly larger change with spline refinement for Sym1, CC, $b_s = 285.5$ mm comes from the steep pre-buckling gradients near the clamped ends of the wide, soft panel. Finer segments lower the critical load slightly for $b_s = 150$ mm and therefore do not remove the offset against Abaqus noted in Section 4.3.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 7.** Critical load $N_{cr}$ [N/mm] with the default discretisation and with refined spline spans or cross-section segments (symmetry side edges); in brackets the change from the default.
:::

| Case | Default | Spans $h_x/2$ | Segments $b_s/40$ |
|:--|--:|--:|--:|
| Isotropic, SS, $b_s = 150$, $\beta = 1.5$ | 122.674 | 122.592 (−0.07 %) | 122.407 (−0.22 %) |
| Sym2, SS, $b_s = 150$, $\beta = 3$ | 20.000 | 20.002 (+0.01 %) | 19.964 (−0.18 %) |
| Unsym, CC, $b_s = 150$, $\beta = 1.5$ | 27.002 | 27.000 (−0.01 %) | 26.881 (−0.45 %) |
| Sym1, CC, $b_s = 285.5$, $\beta = 3$ | 3.4722 | 3.4789 (+0.19 %) | 3.4703 (−0.05 %) |
| Sym1, CC, $b_s = 285.5$, $\beta = 1.5$ | 3.8165 | 3.8266 (+0.27 %) | 3.8133 (−0.08 %) |

For the last case, a further halving of the spans ($h_x/4$) gives 3.8348 N/mm (+0.48 %), so that the spline discretisation converges from below in this case; the remaining uncertainty is about 0.2 %.

## Buckling modes and modal participation

Figure 5 shows the skin deflection of the critical mode of the Unsym panel with SS ends, $b_s = 150$ mm and $\beta = 1.5$, for the two side conditions. In both cases the mode is a local skin buckle with three half-waves along the length; the buckles are skewed by the coupling. With symmetry side edges (Figure 5a) the half-bays at the sides are forced to zero slope across the edge and their deflection differs from that of the central bay; with periodic side edges (Figure 5b) the deflection at $X = -W$ and $X = 0$ is identical (largest difference $2.2\times10^{-16}$ of the maximum, against 0.72 with symmetry edges), and the half-bays at the two sides join into one full bay of the same shape as the central one. The critical load changes from 32.06 to 32.26 N/mm (+0.62 %).

![**Figure 5.** Skin deflection $w/|w|_{max}$ of the critical mode, Unsym, SS, $b_s = 150$ mm, $\beta = 1.5$: (a) symmetry side edges, $N_{cr} = 32.06$ N/mm; (b) periodic side edges, $N_{cr} = 32.26$ N/mm. The stringers lie at $X = -75$ mm and $X = -225$ mm.](fig/fig5_buckling_modes.png){width=15.5cm}

The modal participation, Eq. (@part), identifies the mechanism. For the periodic panel of Figure 5 it is 83.0 % local (L), 8.6 % global (G), 6.9 % distortional (D), 0.9 % shear (S), 0.6 % transverse extension (T) and below 0.1 % Mindlin (M); with symmetry edges it is almost the same (82.9 %, 8.2 %, 7.2 %, 1.2 %, 0.6 %). All 60 critical modes of this study are mainly local, with $p_L$ between 0.79 and 0.92, $p_G \le 0.15$ and $p_D \le 0.08$. Since the participation is measured by displacements and not by energy, the global share mainly reflects the non-zero mean deflection of the cross-section in a local buckle rather than a global mechanism. Such a decomposition is not available from the shell models.

## Computational cost

With the full basis of 482–569 modes, a GBT analysis (cross-section analysis, pre-buckling solution and eigenvalue problem) takes 8.4 s on average (3.1–21.3 s) on two cores of a standard workstation, against 7.3 s for the spline finite strip model with twice as many unknowns. The cross-section analysis is independent of the panel length, and the sparse factors of the geometric matrices (Section 3.10) are reused for all lengths of a parametric study.

## Limitations

The following limitations should be noted. (i) GBT and the strip model share the plate assumptions: laminated Reissner–Mindlin walls with a shear correction factor of 5/6, an offset tie represented by stiff link walls or rigid links, and linear buckling about a first-order pre-buckling state. An error common to both, such as the −1 % offset for $b_s = 150$ mm, cannot be detected by comparing them. (ii) The link walls are a penalty idealisation of the tie; the results are insensitive to the penalty factor (122.656, 122.663 and 122.660 N/mm for $k_L = 10^2$, $10^3$ and $10^4$, Isotropic, SS, $b_s = 150$ mm, $\beta = 1.5$), but a stiff twisting stiffness $D_{66}$ of the link walls would overestimate the critical load by 23 %, which is why the bending and twisting stiffnesses of the link walls are reduced. (iii) The periodic condition of Eq. (@pbc) keeps the width of the cell fixed. A periodic condition with free lateral expansion, in which the two side edges differ by a common jump in $U_X$ and $u$, is not yet implemented in GBT. (iv) A cell of width $2b_s$ admits the Bloch wavenumbers $\kappa b_s = 0$ and $\pi$ only; buckling modes of an infinite panel with a longer period across the stringers require a wider cell or complex Bloch modes. (v) The verified model uses the full basis of 482–569 modes; a reduced, hierarchical basis of a few tens of modes, as typically used in GBT, has not yet been validated for the coupled lay-ups with clamped ends. (vi) The rotation $\beta_x$ of the Mindlin modes does not enter the load term, which is consistent with the shell edge load of the reference model.

# Conclusions

A shear-deformable Generalized Beam Theory with exact periodic boundary conditions was developed for the buckling analysis of omega-stringer-stiffened composite panels. The main conclusions are:

1. Periodic side edges can be built into every family of GBT deformation modes. The conventional modes are obtained by projecting out constraint rows family by family, the warping of the global and distortional modes by a constrained Vlasov least-squares problem, and the shear and Mindlin modes by merging the nodal functions of the two side edges. Every mode, and hence any truncated basis, satisfies the periodic condition to machine precision. The rigid rotation of the cross-section is eliminated, and three additional distortional modes arise from the coupled lateral motion and rotation of the side edges.
2. The formulation retains the structure and the symbols of classical GBT. The modal matrices $\bar{\mathbf B}$, $\bar{\mathbf C}$, $\bar{\mathbf D}_1$, $\bar{\mathbf D}_2$, $\bar{\mathbf F}_1$ and $\bar{\mathbf H}_1$ were given explicitly for laminated walls with membrane–bending and anisotropic coupling, and the system reduces exactly to the classical GBT buckling equations and the classical geometric matrices for uncoupled orthotropic walls under uniform compression.
3. Against 90 Abaqus shell models, the critical loads agree within 0.58 % on average (maximum 1.48 %) with symmetry side edges and within 0.57 % (maximum 1.25 %) with periodic side edges. An independent spline finite strip model agrees with GBT within 0.27 % on average. The change caused by periodic instead of symmetry side edges is reproduced with a mean deviation of 0.05 percentage points and a correlation of 0.995.
4. For uncoupled laminates the symmetry and periodic models are equivalent. For laminates with bending–twisting or membrane–bending coupling, the symmetry model represents a mirrored panel and underestimates the critical load of the periodic panel by up to about 2 % for simply supported and, according to GBT and the strip model, by 5–18 % for clamped loaded edges.
5. GBT identifies all critical modes of the study as local skin buckling with a minor participation of global and distortional modes, at a computing time of a few seconds per case.

Future work will address the confirmation of the clamped periodic cases with shell finite element models, the periodic condition with free lateral expansion, reduced hierarchical mode bases for parametric studies and surrogate models, and the comparison with the closed-form periodic solutions of El Yaakoubi-Mesbah and Mittelstedt {{elyaakoubi2025}}.
