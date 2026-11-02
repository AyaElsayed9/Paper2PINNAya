# Conclusions

Periodic boundary conditions were built into the cross-section analysis of shear-deformable GBT: conventional modes by restricting each mode family to its periodic subspace, the warping by a constrained least-squares problem on the cell as a closed contour, and shear and Mindlin modes by merging the side-node functions. Every mode, and hence every truncated basis, is exactly periodic, and the formulation reduces to classical GBT (confirmed to 0.001 %).

For 30 omega-stringer cells with simply supported loaded edges, the GBT–Ritz critical loads agree with Abaqus shell models within 0.57 % on average and 1.25 % at most, and with clamped loaded edges within 0.64 % and 1.52 %, for all lay-ups; the mode shapes of the six cells compared agree to MAC > 0.999. Cells up to three stringer pitches long buckle locally. From about six pitches, isotropic and angle-ply cells buckle in a long-wave mode driven by the transverse compression that the constant cell width induces, and from twelve pitches cross-ply and unsymmetric cells buckle globally; with free lateral expansion, Abaqus gives 2.0–2.9 times higher critical loads for the two cells examined. Since periodic shell models impose the same constant width, this must be considered for panels whose lateral expansion is not restrained.

Future work will address periodic conditions with free lateral expansion, enriched bases for coupled laminates, the comparison with the closed-form solutions of El Yaakoubi-Mesbah and Mittelstedt {{elyaakoubi2025}}, and the use of the formulation as a data generator for machine-learning surrogate models {{elsayed2026}}.

::: {custom-style="MDPI_6.2_back_matter"}
**Author Contributions:** Conceptualization, A.M.E. and N.K.; methodology, A.M.E. and N.K.; software, A.M.E.; validation, A.M.E.; formal analysis, A.M.E.; investigation, A.M.E.; writing—original draft preparation, A.M.E.; writing—review and editing, N.K. and C.M.; supervision, N.K. and C.M.; project administration, C.M.; funding acquisition, C.M. All authors have read and agreed to the published version of the manuscript.
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Funding:** This research received no external funding.
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Data Availability Statement:** The data presented in this study are available on request from the corresponding author.
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Acknowledgments:** During the preparation of this study, the authors used Claude Opus 5.5 (Anthropic) for literature research and refinement of the manuscript text. The authors have reviewed and edited the output and take full responsibility for the content of this publication.
:::

::: {custom-style="MDPI_6.2_back_matter"}
**Conflicts of Interest:** The authors declare no conflicts of interest.
:::
