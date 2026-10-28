
# Results and discussion

This section verifies the formulation and discusses the buckling behaviour of omega-stringer-stiffened panels with periodic side edges. The critical loads predicted by the GBT–Ritz formulation are compared with Abaqus shell models with periodic constraint equations for simply supported loaded edges and with an independent spline finite strip model for clamped loaded edges. Particular attention is given to the critical loads, the buckling modes, the contribution of the individual deformation modes and the influence of the laminate coupling. The geometry of the cross-section and the material properties are summarised in Table 4, and the lay-ups in Table 5. They include an isotropic reference, a cross-ply laminate without coupling, two symmetric angle-ply laminates with bending–twisting coupling ($D_{16}, D_{26} \neq 0$; Sym2 and the seven-ply stringer of Sym1 also with $A_{16}, A_{26} \neq 0$), and an unsymmetric laminate with membrane–bending coupling ($\mathbf B \neq \mathbf 0$). Ply angles are measured from the stringer axis $x$. Two stringer pitches, $b_s = 150$ mm and $285.5$ mm, and aspect ratios $0.5 \le \beta \le 3$ are considered.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 4.** Geometry of the cross-section (Figure 1) and material properties.
:::

| Parameter | Value |
|:--|:--|
| stringer pitch $b_s$ | 150 mm, 285.5 mm |
| web root half-width $b_3$, web height $h$, web inclination | 27.75 mm, 30 mm, 25° from the vertical |
| crown half-width $h_t = b_3 - h\tan 25^\circ$ | 13.76 mm |
| foot width $b_1 = 0.1\,(b_s - 2b_3)$ | 9.45 mm, 23.00 mm |
| skin and stringer thickness $t_s$, $t_f$ | 1.0 mm, 0.875 mm (8 and 7 plies of 0.125 mm) |
| stringer offset $c = (t_s+t_f)/2$ | 0.9375 mm |
| CFRP ply | $E_1 = 150$ GPa, $E_2 = 9$ GPa, $G_{12} = G_{13} = G_{23} = 5.3$ GPa, $\nu_{12} = 0.33$ |
| isotropic material | $E = 200$ GPa, $\nu = 0.3$ |

::: {custom-style="MDPI_4.1_table_caption"}
**Table 5.** Lay-ups of skin and stringer (angles in degrees from the $x$-axis; subscript s: symmetric, subscript 2: repeated twice).
:::

| Name | Skin | Stringer | Coupling |
|:--|:--|:--|:--|
| Isotropic | – | – | none |
| Cross-ply | [0/90/0/90]$_\mathrm{s}$ | [0/90/0/90/0/90/0] | none |
| Sym1 | [45/−45/−45/45]$_\mathrm{s}$ | [45/−45/45/−45/45/−45/45] | $D_{16}, D_{26}$ (stringer also $A_{16}, A_{26}$) |
| Sym2 | [60/−30/−30/60]$_\mathrm{s}$ | [60/−30/60/−30/60/−30/60] | $A_{16}, A_{26}, D_{16}, D_{26}$ |
| Unsym | [0/30/60/90]$_2$ | [0/30/60/90/0/30/60] | $\mathbf A$, $\mathbf B$, $\mathbf D$ fully populated |

A uniform compressive line load of 1 N/mm acts at both loaded edges on the skin, the webs and the crowns, and 2 N/mm on a foot, which carries skin and flange; the foot load acts at the height of its resultant, $Y = c/2$. The eigenvalue $\lambda$ is therefore the critical line load $N_{cr}$ in N/mm, the quantity reported by the buckling analysis of the shell model.

The shell models are linear buckling analyses in Abaqus with S4R and S3 elements and a seed of 3.7 mm. The stringer is tied to the skin by \*TIE with the true offset of the mid-surfaces, and the loads are applied as shell edge loads. The simply supported loaded edges are modelled by $\mathcal U_2 = 0$ on every end edge of skin and stringer. For the periodic side edges, every skin node of the edge $X = -W$ is paired with the node of $X = 0$ at the same $x$, and one constraint equation $\mathcal U_\alpha(-W,x) - \mathcal U_\alpha(0,x) = 0$ is written for each pair and each of the six degrees of freedom, as in Eq. (@pbc); at the two end corners the equation for $\mathcal U_2$ is omitted because $\mathcal U_2 = 0$ is already prescribed there. The rigid-body motions left by periodicity are removed by $\mathcal U_1 = \mathcal U_3 = 0$ at the corner $X = 0$, $x = 0$. The GBT model imposes the same physical conditions (Table 2). Abaqus results with periodic side edges are available for the simply supported loaded edges.

For the clamped loaded edges, an in-house spline finite strip model (FSM) {{cheung1976,lau1986}} serves as the reference. It discretises the same geometry, laminates, loads and end conditions with flat laminated Mindlin strips of width $b_s/50$ (twice as fine as the GBT segments), linear interpolation across each strip, a mid-strip transverse shear strain against shear locking, and six degrees of freedom per nodal line in global axes. Along $x$ it uses the same cubic B-splines as the GBT model. The web roots are joined to the skin by rigid offset links, the direct analogue of the \*TIE of the shell model, and the periodic condition is applied node by node by merging the six degrees of freedom of the two side nodal lines for every spline coefficient, which is literally the constraint set of the shell model. The pre-buckling state and the geometric stiffness are computed as in Section 2.3. The FSM shares the plate theory with GBT but none of the cross-section discretisation: it has no deformation modes, a finer mesh, rigid links instead of link walls, and boundary conditions on nodal degrees of freedom instead of built into modes. Its agreement with GBT is therefore an independent check of the mode construction. With simply supported loaded edges it agrees with Abaqus within 0.59 % on average and 1.37 % at most (Table 6).

## Verification of the formulation

Since the periodic condition is built into the modes, every buckling mode satisfies it identically. For the cross-ply cell with $b_s = 150$ mm, the largest difference between the two side nodes over all conventional modes is $7.9\times10^{-15}$ for the in-plane displacements, $1.1\times10^{-16}$ for the rotation $\theta$ and $2.3\times10^{-11}$ for the warping (modes scaled to a unit maximum), and it is zero for the shear and Mindlin modes. In the buckling mode of the Unsym panel with clamped loaded edges, $b_s = 150$ mm and $\beta = 1.5$, the largest differences between the side edges along the whole length are $1.2\times10^{-18}$ ($U_X$), $6.9\times10^{-18}$ ($U_Y$), $1.4\times10^{-20}$ ($\theta$) and $1.1\times10^{-15}$ (warping), with modal amplitudes of order $10^{-5}$ to $10^{-1}$, while the side edges themselves move.

The explicit modal tensors of Eqs. (@Cik)–(@couplingT) were compared with the cross-section integral of Eq. (@dU0) evaluated numerically by the implementation for the Cross-ply, Sym2 and Unsym cells; all blocks agree to within $2\times10^{-15}$ relative to their largest entry, i.e. to round-off. The classical limit of Section 2.3 was checked for the Cross-ply cell with $b_s = 150$ mm by restricting the model to the conventional modes and imposing a uniform compressive membrane force and the classical end condition $\varphi_i(0) = \varphi_i(L) = 0$. The B-spline Ritz solution then reproduces the minimum of the signature-curve eigenvalue problem of Eq. (@signature), 33.747 N/mm with three half-waves for $L = 225$ mm, to within 0.001 %; the same value is obtained for $L = 450$ mm with six half-waves.

## Pre-buckling state

Figure 3 shows the computed first-order membrane forces of the cells with $b_s = 150$ mm and $\beta = 1.5$, normalised by the reference load $n_0$. At mid-length (Figure 3a), the axial membrane force in the skin lies between $1.00\,n_0$ (Cross-ply, Unsym) and $1.16\,n_0$ (Sym1, Sym2), and the feet, which consist of skin and flange, carry about twice this value; the end load is thus redistributed between skin and stringer within a transition zone near the loaded edges (Figure 3c). The transverse membrane force (Figure 3b) is nearly uniform across the width and amounts to $0.32\,n_0$ for the isotropic cell, $0.04\,n_0$ and $0.03\,n_0$ for the Cross-ply and Unsym cells, and $0.86$–$0.89\,n_0$ for the angle-ply cells. It is the force of a restrained Poisson expansion: the periodic condition, like the constraint equations of the shell model, keeps the width of the cell constant, and for the isotropic, Cross-ply and Sym1 skins the ratio $S_{ss}/S_{xx}$ at mid-length agrees with $A_{12}/A_{11}$ of the skin laminate to within 2 %. The transverse compression is therefore part of the loading of the periodic panel and enters $\mathbf G$ through Eq. (@Xgeo); for the angle-ply cells, whose $\pm45^\circ$ and $60^\circ/{-30}^\circ$ plies have a large Poisson ratio, it is almost as large as the axial force, which contributes to their low critical loads. With clamped loaded edges (dashed lines in Figure 3c,d), the restraint of the in-plane end displacements raises the axial force at the clamped end $x = 0$, up to $1.43\,n_0$ in the skin bay of the Sym1 cell. At the end $x = L$, where the warping is free, the membrane forces of the angle-ply cells oscillate slightly within the last spline spans, consistent with the slower spline convergence of the clamped angle-ply cells reported in Section 3.7.

![**Figure 3.** Pre-buckling membrane forces of the periodic cells ($b_s = 150$ mm, $\beta = 1.5$) per unit reference load $n_0$ (compression positive): (a) axial and (b) transverse membrane force in the skin and the feet at mid-length with simply supported loaded edges (shaded: feet); (c) axial and (d) transverse membrane force along the cell at the middle of the skin bay, $X = -b_s$, for simply supported (solid lines, open markers) and clamped (dashed lines, filled markers) loaded edges.](fig/fig03_prebuckling.png){width=15.5cm}

## Simply supported loaded edges

Figure 4 shows the critical load as a function of the aspect ratio for the two stringer pitches and all lay-ups, together with the Abaqus results, Table 6 lists all 30 cases that were computed with Abaqus, and Figure 5a shows the deviations case by case. GBT agrees with Abaqus within 0.57 % on average and 1.25 % at most; the FSM agrees within 0.59 % and 1.37 %, and GBT and FSM differ by 0.22 % on average. The accuracy is the same for the uncoupled and the coupled laminates (mean deviation 0.49–0.73 % for every lay-up), which confirms the treatment of the bending–twisting and membrane–bending coupling in Eqs. (@Cik)–(@couplingT).

![**Figure 4.** Critical load $N_{cr}$ (logarithmic scale) as a function of the aspect ratio $\beta$ of the periodic cell with simply supported loaded edges for (a) $b_s = 150$ mm and (b) $b_s = 285.5$ mm: GBT (lines, computed at the small dots) and Abaqus shell models with periodic constraint equations (open markers).](fig/fig04_ncr_ss.png){width=15.5cm}

For the isotropic, cross-ply and unsymmetric cells, the critical load varies only little with the aspect ratio, with local maxima where the length of the cell does not accommodate an integer number of buckle half-waves of the preferred length, as known from the garland curves of plates. For $b_s = 150$ mm such a maximum occurs at $\beta = 0.75$, and the Abaqus results follow it (cross-ply: 34.848 N/mm in Abaqus and 34.530 N/mm in GBT at $\beta = 0.75$, against 31.464 and 31.243 N/mm at $\beta = 1.5$). For the angle-ply cells Sym1 and Sym2, the critical load decreases monotonically with $\beta$ towards the value of the long cell. Their bending–twisting coupling skews the buckles, which therefore cannot adapt to the straight loaded edges by an integer number of half-waves, and the end conditions stiffen short cells.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 6.** Critical loads $N_{cr}$ [N/mm] with periodic side edges and simply supported loaded edges: Abaqus, GBT and spline finite strip model (FSM), and their deviations $e$ from Abaqus.
:::

@@TABSS@@

The deviations depend mainly on the stringer pitch. For $b_s = 285.5$ mm the mean deviation of GBT is 0.20 % (maximum 0.46 %). For $b_s = 150$ mm it is 0.93 % (maximum 1.25 %), and all 15 GBT results and all 15 FSM results lie below Abaqus (GBT −1.25 % to −0.63 %). This offset is therefore common to the two plate models and not caused by the GBT mode basis (Figure 5a), and the refinement study of Section 3.7 shows that it is not a discretisation error of GBT either, since finer segments lower the GBT result slightly. Its likely origins are differences between the S4R shell elements and the plate model in the short skin bays of the narrow pitch, and the representation of the tie; it has not been resolved in this study.

![**Figure 5.** Deviation $\Delta$ of the critical loads of the periodic cells: (a) simply supported loaded edges, GBT (markers in the colours of the lay-ups) and spline FSM (crosses) against Abaqus; (b) clamped loaded edges, GBT against the spline FSM. Filled markers: $b_s = 150$ mm, open markers: $b_s = 285.5$ mm; within each group the aspect ratio increases from left to right ($\beta = 0.75$, 1.5, 3). Shaded: $\pm1$ %.](fig/fig05_accuracy.png){width=15.5cm}

## Clamped loaded edges

The critical loads with clamped loaded edges are shown in Figure 6 and listed in Table 7. Clamping raises the critical load of short cells considerably, e.g. by 59 % for the isotropic cell with $b_s = 150$ mm and $\beta = 0.5$, and the curves decrease monotonically with $\beta$. For long cells the clamped values approach, and for some lay-ups fall slightly below, the simply supported ones (e.g. 118.13 against 121.41 N/mm for the isotropic cell with $\beta = 3$); the FSM gives the same trend (118.17 against 121.40 N/mm). In the present models the clamped edges also restrain the in-plane displacements of the loaded edges, which changes the pre-buckling state near the ends. GBT and FSM agree within 0.30 % on average and 1.17 % at most over the 30 cases (Figure 5b), with the largest deviations for the angle-ply lay-ups, in particular Sym1 with $b_s = 285.5$ mm, where the long skin bays and the bending–twisting coupling produce steep gradients of the mode near the clamped ends (Section 3.7). Abaqus results with periodic side edges and clamped loaded edges are not yet available, so that the clamped results are verified by the FSM only.

![**Figure 6.** Critical load $N_{cr}$ (logarithmic scale) as a function of the aspect ratio $\beta$ of the periodic cell with clamped loaded edges for (a) $b_s = 150$ mm and (b) $b_s = 285.5$ mm: GBT (lines, computed at the small dots) and spline finite strip model (open markers).](fig/fig06_ncr_cc.png){width=15.5cm}

::: {custom-style="MDPI_4.1_table_caption"}
**Table 7.** Critical loads $N_{cr}$ [N/mm] with periodic side edges and clamped loaded edges: spline finite strip model (FSM) and GBT, and the deviation $e$ of GBT from the FSM.
:::

@@TABCC@@

## Mode decomposition and buckling modes

Figure 7 shows the critical mode of the Unsym cell with $b_s = 150$ mm and $\beta = 1.5$ for both end conditions, as a view of the deformed cell and as the skin deflection. In both cases the mode is a local buckling of the skin bays between the stringers, which move only slightly, and the deflection at the side edges $X = -W$ and $X = 0$ is identical (largest difference $2.2\times10^{-16}$ of the maximum deflection), so that the half-bays at the two sides join into one full bay of the same shape as the central one. With simply supported loaded edges the skin bays buckle in three half-waves along the cell, skewed by the laminate coupling; with clamped loaded edges two half-waves form, which are concentrated away from the loaded edges and are hardly skewed.

![**Figure 7.** Critical buckling mode of the Unsym cell with periodic side edges, $b_s = 150$ mm and $\beta = 1.5$: (a, b) oblique view of the deformed cell (amplitude exaggerated, colour: displacement normal to the walls) and (c, d) skin deflection $w/|w|_{\max}$ with contour lines at $\pm0.25$, $\pm0.5$, $\pm0.75$ (dashed: $w = 0$); (a, c) simply supported and (b, d) clamped loaded edges. Dotted lines: stringer axes $X = -75$ mm and $X = -225$ mm.](fig/fig07_buckling_modes.png){width=15.5cm}

The energy-based mode decomposition of Eq. (@rho) quantifies these observations for all panels (Table 8 and Figure 8a–c). In every case the local family carries 96–101 % of the strain energy of the critical mode, the global family up to @@GMAX@@ % and the distortional family between @@DMIN@@ % and @@DMAX@@ %; the transverse-extension and shear families contribute less than @@TSMMAX@@ %, and the Mindlin family contributes nothing by construction (Section 2.4). Contributions slightly above 100 % or below zero are possible in Eq. (@rho), because $U_g$ contains the interaction energy of family $g$ with the other families; small negative distortional contributions of the angle-ply laminates indicate that the distortional modes interact with the local skin buckles. The global participation is largest for the cross-ply cell with $b_s = 285.5$ mm, i.e. the local buckles of the wide skin bays are accompanied by a small overall bending of the cell. The decomposition does not change significantly with the aspect ratio, i.e. no transition to a distortional or global buckling mode occurs in the range studied. For the Unsym cell with simply supported loaded edges and $b_s = 150$ mm (Figure 8b), the global participation alternates between about 2 % and zero with $\beta$, i.e. between buckling modes with and without a small rigid-body translation of the cross-section, which corresponds to an overall bending of the panel.

The hierarchical modes of Section 2.2.5 allow a finer decomposition in the manner of classical GBT. Figure 8d shows the participation of the individual hierarchical conventional modes in the critical mode of the Unsym cell with simply supported loaded edges, $b_s = 150$ mm and $\beta = 1.5$, computed with the complete basis expressed in hierarchical modes ($N_{cr}$ = @@HN@@ N/mm, relative difference from the nodal basis @@HDIFF@@). Four modes, @@HMODES@@, carry @@H4@@ % of the strain energy, and the first 20 modes @@H20@@ %. Their amplitude functions (Figure 8e) have the same three half-waves, and their cross-section shapes (Figure 8f) combine the bending of the skin bays with a small rotation and distortion of the stringers. Such a diagram identifies the modes that a reduced basis must contain; the next section shows, however, that the many modes with small participation still affect the critical load of coupled laminates.

![**Figure 8.** Mode decomposition of the critical modes. (a)–(c) Energy-based participation $\rho_g$ of the mode families as a function of $\beta$ (upper part: local family; lower part: other families; the Mindlin participation vanishes identically): (a) Cross-ply, $b_s = 285.5$ mm, SS; (b) Unsym, $b_s = 150$ mm, SS; (c) Unsym, $b_s = 150$ mm, CC. (d) Participation $\rho_k$ of the first 60 hierarchical conventional modes (ordered at $k_r = \pi/b_s$) in the critical mode of the Unsym cell, $b_s = 150$ mm, $\beta = 1.5$, SS; (e) amplitude functions of the four dominant modes, scaled by the largest in-plane displacement of each mode shape; (f) in-plane shapes of these modes (grey: undeformed).](fig/fig08_participation.png){width=15.5cm}

::: {custom-style="MDPI_4.1_table_caption"}
**Table 8.** Energy-based participation $\rho_g$ [%] of the mode families in the critical mode, range over $0.5 \le \beta \le 3$ (eight aspect ratios).
:::

@@TABRHO@@

## Hierarchical mode bases

The results above use the complete periodic mode basis, which is equivalent to imposing the periodic condition node by node in a folded-plate model with the same nodes. The characteristic feature of GBT is that a small number of hierarchically ordered modes suffices for many problems. This was examined with the hierarchical conventional modes of Section 2.2.5, ordered at $k_r = \pi/b_s$, the pre-buckling state of the complete basis and the end conditions of Section 2.3.5. Table 9 and Figure 9 show the deviation of the critical load when only the first $n_C$ hierarchical conventional modes are kept, together with all shear and Mindlin modes, and Table 9 also the deviation when the Mindlin modes are omitted from the complete basis.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 9.** Deviation [%] of the critical load from the complete basis for truncated mode bases ($\beta = 1.5$): the first $n_C$ hierarchical conventional modes with all shear and Mindlin modes, and all conventional and shear modes without the Mindlin modes.
:::

@@TABRED@@

@@REDTEXT@@

![**Figure 9.** Convergence of the critical load with the number $n_C$ of hierarchical conventional modes (all shear and Mindlin modes kept, $\beta = 1.5$): deviation from the complete basis, which contains $n_C = 273$ ($b_s = 150$ mm) and 231 ($b_s = 285.5$ mm) conventional modes; logarithmic scales. Lay-ups and end conditions as in Table 9.](fig/fig09_hierarchical.png){width=13.5cm}

@@HCONV@@

## Accuracy and efficiency

The default discretisation (segments of at most $b_s/25$, spline spans $h_x = b_s/20$) was checked by refining each parameter in turn (Table 10). @@CONVTEXT@@

::: {custom-style="MDPI_4.1_table_caption"}
**Table 10.** Critical load $N_{cr}$ [N/mm] with the default discretisation and with refined spline spans or cross-section segments (periodic side edges); in brackets the change from the default.
:::

@@TABCONV@@

The model sizes and computing times are compared in Table 11. With the full basis of 482–569 modes, a GBT analysis, including the cross-section analysis, the pre-buckling solution and the eigenvalue problem, takes @@TGSEQ@@ s on average for the ten cells with $b_s = 150$ mm and $\beta = 1.5$ on two cores of a standard workstation; over all 160 cases of Figures 4 and 6, computed two at a time on the same two cores, it ranged from @@TMIN@@ to @@TMAX@@ s. The GBT model has about half the unknowns of the FSM and less than half the degrees of freedom of the shell model, but with the full basis it is not faster than the FSM (@@TSSEQ@@ s): the global and distortional modes extend over the whole cross-section and couple with all other modes, so that its matrices are less sparse than those of the nodal FSM. The cross-section analysis itself takes about 0.1 s and is reused for all lengths of a parametric study. The efficiency of the GBT formulation therefore lies in its modal description, which provides the mode decomposition and is the basis for reduced mode sets, rather than in the solution time of the full basis.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 11.** Model size and computing time for $b_s = 150$ mm and $L = 225$ mm.
:::

| Model | Unknowns | Computing time per case |
|:--|:--|:--|
| Abaqus shell model (S4R/S3, seed 3.7 mm) | approx. $5\times10^{4}$ DOFs (estimated from the seed) | [to be added by the authors] |
| spline finite strip model | 37 170 (177 nodal lines × 6 × 35 splines) | @@TSSEQ@@ s on average |
| GBT–Ritz (present) | 18 777 (569 modes × 33 splines) | @@TGSEQ@@ s on average |

The following limitations should be noted. First, GBT and the FSM share the plate assumptions: laminated Reissner–Mindlin walls with a shear correction factor of 5/6, an offset tie represented by stiff link walls or rigid links, and linear buckling about a first-order pre-buckling state. An error common to both, such as the −1 % offset for $b_s = 150$ mm, cannot be detected by comparing them. Second, the periodic condition of Eq. (@pbc), like the constraint equations of the shell model, keeps the width of the cell constant. The Poisson expansion of the cell under the axial load is therefore restrained, and transverse membrane forces $S_{ss}$ arise in the pre-buckling state; for the isotropic cell with $b_s = 150$ mm and $\beta = 3$, omitting $S_{ss}$ from $\mathbf G$ would raise the critical load by 23 %. The condition thus represents a panel whose lateral expansion is restrained by the surrounding structure. A periodic condition with free lateral expansion, in which the two side edges differ by a common jump in $U_X$, requires one additional generalised coordinate and is left for future work. Third, a cell of width $2b_s$ admits the Bloch wavenumbers $\kappa b_s = 0$ and $\pi$ only; buckling modes with a longer period across the stringers require a wider cell. Fourth, compact hierarchical bases are adequate for the uncoupled laminates only (Section 3.6); for the coupled laminates, enriched bases remain to be developed.
