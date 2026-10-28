
# Conclusions

This study presented a shear-deformable GBT–Ritz formulation for the buckling analysis of omega-stringer-stiffened composite panels in which the periodic boundary conditions of the panel cell are built into the cross-section analysis. The conventional deformation modes are obtained by projecting out the constraint rows of the periodic condition family by family, the warping of the global and distortional modes follows from a constrained Vlasov least-squares problem, and the shear and Mindlin modes are obtained by merging the nodal functions of the two side edges. Every mode, and therefore every combination and truncation of the modes, satisfies the periodic condition to machine precision, so that no penalty terms or Lagrange multipliers are needed in the member analysis. The modal tensors, the equilibrium equations and the boundary conditions were derived in the notation of the shear-deformable GBT of Silvestre and Camotim, extended by the membrane terms that the closed contours of the cell require, the transverse shear (Mindlin) modes and a computed, non-uniform pre-buckling state; for uncoupled walls under uniform compression the formulation reduces to the classical GBT buckling equations, which was confirmed numerically to 0.001 %. For 30 panels with simply supported loaded edges, the critical loads agree with Abaqus shell models with periodic constraint equations within 0.57 % on average and 1.25 % at most, and for 30 panels with clamped loaded edges they agree with an independent spline finite strip model within 0.30 % on average and 1.17 % at most. The accuracy is the same for cross-ply, angle-ply and unsymmetric laminates, which confirms the consistent treatment of the bending–twisting and membrane–bending coupling.

Beyond accuracy, the formulation provides insight into the buckling mechanisms of the cell. The energy-based mode decomposition identifies the critical modes of all panels studied as local skin buckling, with global and distortional contributions below 4 %; expressed in hierarchical modes, the critical mode of the unsymmetric cell is dominated by four conventional modes in the manner of classical GBT. The computed pre-buckling state shows that the periodic condition, which keeps the width of the cell constant, induces a transverse membrane compression of up to about 0.9 times the applied load in the angle-ply skins, which must be considered when periodic results are transferred to panels with free lateral edges. The variation of the critical load with the aspect ratio of the cell follows directly from the same model at a computing time of a few seconds per case. Since the periodic condition is a property of the modes, hierarchical GBT bases can be formed within the periodic subspace: for uncoupled laminates a few tens of conventional modes reproduce the complete result within about 1 %, whereas the coupled laminates need about 100 conventional modes together with the transverse shear modes, which contribute up to about 2 % to the critical load. Future work will address the confirmation of the clamped cases with shell finite element models, periodic conditions with free lateral expansion of the cell, reduced hierarchical mode bases for parametric studies and surrogate models, and the comparison with the closed-form periodic solutions of El Yaakoubi-Mesbah and Mittelstedt {{elyaakoubi2025}}.

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
