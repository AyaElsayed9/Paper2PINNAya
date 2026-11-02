# Conclusions

Periodic boundary conditions were built into the cross-section analysis of shear-deformable GBT: conventional modes by restricting each mode family to its periodic subspace, the warping by a constrained least-squares problem on the cell as a closed contour, and shear and Mindlin modes by merging the side-node functions. Every mode, and hence every truncated basis, is exactly periodic, and the formulation reduces to classical GBT (confirmed to 0.001 %).

For 30 omega-stringer cells with simply supported loaded edges, the GBT–Ritz critical loads agree with Abaqus shell models within @@SSMEAN@@ % on average and @@SSMAX@@ % at most@@CONC_CC@@, for all lay-ups; the mode shapes of the six cells compared agree to MAC > 0.999. Cells up to three stringer pitches long buckle locally. From about six pitches, isotropic and angle-ply cells buckle in a long-wave mode driven by the transverse compression that the constant cell width induces, and from twelve pitches cross-ply and unsymmetric cells buckle globally; with free lateral expansion, Abaqus gives @@FRMIN@@–@@FRMAX@@ times higher critical loads for the two cells examined. Since periodic shell models impose the same constant width, this must be considered for panels whose lateral expansion is not restrained.

Future work will address periodic conditions with free lateral expansion, which require additional generalised coordinates for the common displacement jumps of the side edges (Section 3.5), enriched bases for coupled laminates (Section 3.6), and the comparison with the closed-form solutions of El Yaakoubi-Mesbah and Mittelstedt {{elyaakoubi2025}}.

The formulation will also serve as the physics model of the next step of the authors' work on physics-informed neural networks (PINNs). In {{elsayed2026}}, a hybrid PINN based on the deep energy method predicts the local buckling load of omega-stringer-stiffened panels: its physics loss is a plate energy with a reduced bending stiffness, it is complemented by Abaqus results selected by active learning, and its epistemic uncertainty is estimated from an ensemble of networks. Replacing the plate energy by the GBT–Ritz energy of the present formulation is planned in three steps. First, the network maps the design (stringer pitch, aspect ratio, end conditions and the lamination parameters of skin and stringer) to the amplitude functions of a fixed set of periodic GBT modes, so that every prediction satisfies the periodic condition exactly (Proposition 1), without penalty terms. Second, the loss is the Rayleigh quotient of the predicted modes with the matrices $\mathbf K$ and $\mathbf G$ of Eq. (@eig), including the computed pre-buckling state with $S_{ss}$; it requires no reference solutions, gives upper bounds of the critical loads and covers local, distortional and global buckling as well as the long-wave mode of Section 3.5 in one model. Third, a few eigenvalue iterations started from the predicted modes in the same basis give a physics-based reliability check, which can replace the fixed tolerance window used in {{elsayed2026}} to decide when a reference solution is needed, so that Abaqus is required only for blind validation. The family participations of Eq. (@rho) then make the predicted modes interpretable as local, distortional or global.

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
