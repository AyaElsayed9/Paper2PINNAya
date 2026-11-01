# Conclusions

Periodic boundary conditions were built into the cross-section analysis of shear-deformable GBT: the conventional modes are obtained by projecting the periodic constraint rows out of each mode family, the warping by a constrained least-squares problem on the cell as a closed contour, and the shear and Mindlin modes by merging the side-node functions. Every mode, and hence every truncated basis, satisfies the periodic condition exactly, without penalty terms or Lagrange multipliers, and the formulation reduces to classical GBT (confirmed to 0.001 %).

For 30 omega-stringer cells with simply supported loaded edges, the GBT–Ritz critical loads agree with Abaqus shell models within @@SSMEAN@@ % on average and @@SSMAX@@ % at most@@CONC_CC@@, for isotropic, cross-ply, angle-ply and unsymmetric laminates, and the modes agree with a MAC of at least @@MACMIN@@. Cells up to three stringer pitches long buckle locally. Longer cross-ply and unsymmetric cells buckle globally, whereas isotropic and angle-ply cells buckle in a long-wave mode driven by the transverse compression that the constant cell width induces; with free lateral expansion, Abaqus gives @@FRMIN@@–@@FRMAX@@ times higher critical loads for the two cells examined. Since the usual periodic constraint equations of shell models impose the same constant width, this effect must be considered when periodic results are transferred to panels whose lateral expansion is not restrained.

Future work will address periodic conditions with free lateral expansion, enriched bases for coupled laminates, the comparison with the closed-form solutions of El Yaakoubi-Mesbah and Mittelstedt {{elyaakoubi2025}}, and the use of the formulation as a data generator for machine-learning surrogate models {{elsayed2026}}.

::: {custom-style="MDPI_6.2_back_matter"}
**Author Contributions:** Conceptualization, A.M.E. and N.K.; methodology, A.M.E. and N.K.; software, A.M.E.; validation, A.M.E.; formal analysis, A.M.E.; investigation, A.M.E.; writing—original draft preparation, A.M.E.; writing—review and editing, N.K. and C.M.; supervision, N.K. and C.M.; project administration, C.M.; funding acquisition, C.M. All authors have read and agreed to the published version of the manuscript. [To be confirmed.]
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Funding:** This research received no external funding. [To be confirmed.]
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Data Availability Statement:** The data presented in this study are available on request from the corresponding author. [To be confirmed.]
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Acknowledgments:** During the preparation of this work, the authors used Claude (Anthropic) [model version to be confirmed] for code development, data analysis, figure generation and drafting of the manuscript text. The authors have reviewed and edited all output and take full responsibility for the content of this publication.
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Conflicts of Interest:** The authors declare no conflicts of interest.
:::
