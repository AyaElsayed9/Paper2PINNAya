
# Conclusions

A shear-deformable GBT–Ritz formulation was presented for the buckling of omega-stringer-stiffened composite panels, in which the periodic boundary conditions of the panel cell are built into the cross-section analysis. The conventional modes are obtained by projecting the periodic constraint rows out of each mode family, the warping by a constrained least-squares problem in which the cell acts as a closed contour, and the shear and Mindlin modes by merging the nodal functions of the two side edges. Every mode, and hence every combination and truncation, satisfies the periodic condition exactly, without penalty terms or Lagrange multipliers. The modal tensors, equilibrium equations and boundary conditions were derived in the notation of Silvestre and Camotim, extended by the membrane terms of closed contours, transverse shear modes and a computed pre-buckling state, and reduce to classical GBT, which was confirmed to 0.001 %.

For 30 cells with simply supported loaded edges, the critical loads agree with Abaqus shell models with periodic constraint equations within @@SSMEAN@@ % on average and @@SSMAX@@ % at most@@CONC_CC@@, with the same accuracy for cross-ply, angle-ply and unsymmetric laminates. Cells up to three stringer pitches long buckle locally, with global and distortional participations below 4 %; in hierarchical form, a few conventional modes dominate the critical mode, as in classical GBT. Longer cells change their mode: cross-ply and unsymmetric cells buckle globally, whereas isotropic and angle-ply cells buckle in a long-wave mode with a period of two stringer pitches, driven by the transverse compression that the constant cell width induces, at critical loads up to @@FACTOR@@ times lower than without it. Since the usual periodic constraint equations of shell models impose the same constant width, this effect must be considered when periodic results are transferred to panels whose lateral expansion is not restrained.

Future work will address periodic conditions with free lateral expansion, enriched hierarchical bases for coupled laminates, the comparison with the closed-form periodic solutions of El Yaakoubi-Mesbah and Mittelstedt {{elyaakoubi2025}}, and the use of the formulation as the physics-based data generator for machine-learning surrogate models of stiffened panels {{elsayed2026}}.

::: {custom-style="MDPI_6.2_back_matter"}
**Author Contributions:** Conceptualization, A.M.E. and N.K.; methodology, A.M.E. and N.K.; software, A.M.E.; validation, A.M.E.; formal analysis, A.M.E.; investigation, A.M.E.; writing—original draft preparation, A.M.E.; writing—review and editing, N.K. and C.M.; supervision, N.K. and C.M.; project administration, C.M.; funding acquisition, C.M. All authors have read and agreed to the published version of the manuscript. [To be confirmed by the authors.]
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Funding:** [To be completed by the authors.]
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Data Availability Statement:** The data presented in this study are available on request from the corresponding author.
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Conflicts of Interest:** The authors declare no conflicts of interest.
:::
