::: {custom-style="MDPI_1.1_article_type"}
Article
:::

::: {custom-style="MDPI_1.2_title"}
Periodic Boundary Conditions in Generalized Beam Theory: A Shear-Deformable GBT–Ritz Formulation for the Buckling of Omega-Stringer-Stiffened Composite Panels
:::

::: {custom-style="MDPI_1.3_authornames"}
Aya Medhat Elsayed ^1^, Navid Kharghani ^1,\*^ and Christian Mittelstedt ^1^
:::

::: {custom-style="MDPI_1.6_affiliation"}
^1^ Institute for Lightweight Engineering and Structural Mechanics (LSM), Department of Mechanical Engineering, Technical University of Darmstadt, Otto-Berndt-Straße 2, 64287 Darmstadt, Germany
:::

::: {custom-style="MDPI_1.6_affiliation"}
\* Correspondence: navid.kharghani\@lsm.tu-darmstadt.de [author list, affiliations and corresponding author to be confirmed by the authors]
:::

::: {custom-style="MDPI_1.7_abstract"}
**Abstract:** Stringer-stiffened composite panels are commonly analysed on one representative cell whose side edges satisfy periodic boundary conditions. Generalized Beam Theory (GBT) describes the buckling of thin-walled members by cross-section deformation modes, but periodic side edges have not been treated in its cross-section analysis. This paper builds the periodic condition into every family of GBT modes, so that every mode and every truncated basis satisfies it exactly. A shear-deformable GBT–Ritz member analysis of omega-stringer-stiffened panels follows, with laminated Reissner–Mindlin walls, the true stringer offset, a computed pre-buckling state and B-spline amplitude functions. The formulation is derived in the notation of Silvestre and Camotim and reduces to classical GBT for uncoupled walls. For 30 cells with simply supported loaded edges, the critical loads agree with Abaqus shell models with periodic constraint equations within @@SSMEAN@@ % on average and @@SSMAX@@ % at most@@ABS_CC@@. Short cells buckle locally. Longer cross-ply and unsymmetric cells buckle globally, whereas isotropic and angle-ply cells buckle in a long-wave mode that is driven by the transverse compression induced by the constant cell width; Abaqus confirms this mode and gives @@FRMIN@@–@@FRMAX@@ times higher critical loads for cells free to expand laterally.
:::

::: {custom-style="MDPI_1.8_keywords"}
**Keywords:** generalized beam theory; Ritz method; periodic boundary conditions; stiffened panels; omega stringer; laminated composites; buckling
:::

# Introduction

Stringer-stiffened panels are the standard construction of aircraft fuselage and wing covers, and buckling is often their governing design criterion {{zhang2024,mittelstedt2025book}}. The omega (hat) stringer is attractive because its closed cell provides high bending and torsional stiffness and a wide bonding area {{mittelstedtbeerhorst2009,schilling2020tws}}. Since a cover carries many identical stringers, its buckling away from the longitudinal edges is represented by one repeating cell whose side edges satisfy periodic boundary conditions {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. The principle goes back to the theory of waves in periodic structures {{mead1996}}; exact strip programs such as VIPASA and VICONOPT treat infinitely wide panels through a repeating portion {{wittrick1974,kennedy1994}}, as do finite strip methods {{cheung1976,lau1986}}.

At the Institute for Lightweight Engineering and Structural Mechanics (LSM), closed-form models for stiffened composite panels with periodic boundary conditions have been developed for blade stringers {{mittelstedt2007,mittelstedt2009periodic}} and omega stringers {{mittelstedtbeerhorst2009,mittelstedtschroder2010}}, including local and column buckling {{schilling2020tws,schilling2020rineng}}, shear and combined loading {{schilling2021shear,schilling2022interaction}} and postbuckling {{schilling2022post}}, and most recently closed-form and computational solutions for local buckling and postbuckling of omega-stringer-stiffened panels {{elyaakoubi2025,elyaakoubi2024wccm,elyaakoubi2026}}; see also {{mittelstedt2025stiff}}. These models are efficient but idealise the skin–stringer interaction by elastic restraints, treat local and global buckling separately and assume a uniform pre-buckling stress; membrane–bending coupling is usually included through a reduced bending stiffness {{ewing1988}}, as in the physics-informed neural network of the authors {{elsayed2026}}. Shell finite element models with periodic constraint equations avoid these idealisations, but they are expensive in parametric studies and give the buckling mode only as a nodal field.

Generalized Beam Theory (GBT), originated by Schardt {{schardt1989,schardt1994}} and developed further by Davies and Leach {{davies1994}}, represents the deformation of a prismatic thin-walled member by cross-section deformation modes with amplitude functions along the axis; its modal decomposition into global, distortional, local and shear deformations makes it a transparent tool for coupled instabilities {{camotim2010,mittelstedt2020amr,mittelstedt2024gbt,mittelstedt2025gbtch}}. Silvestre and Camotim formulated GBT for orthotropic and laminated members {{silvestre2002a,silvestre2002b}} including shear deformation {{silvestre2013shear,silvestre2013vib}}. The cross-section analysis was extended to branched sections {{dinis2006}}, reformulated {{goncalves2010,bebiano2015}} with shear modes {{goncalves2014shear}} and implemented in GBTUL 2.0 {{bebiano2018}}; closed sections were treated for tubes {{silvestre2007cyl,goncalves2013polygon}}, and GBT was applied to FRP columns {{silva2010}}, stiffened composite panels {{silva2010ijssd}} and stiffened sections {{taig2014}}. The Ritz method complements GBT in the member analysis {{mittelstedt2025ritz,vescovini2018,kuehn2014}}: Kharghani and Mittelstedt combined both into a GBT–Ritz approach for global, local and distortional buckling of shear-deformable composite beams {{kharghani2025jcs}}, composite box girders {{kharghani2025semc}} and beams with unsymmetric lay-ups {{kharghani2026cepa,kharghani2026jcs}}, in which GBT defines the kinematics and the Ritz method the amplitude functions.

GBT has been developed for members with free longitudinal edges or closed contours. The skin of a panel cell, however, ends at two side edges whose displacements must be equal. This condition has not been formulated in GBT, and it cannot simply be imposed on a truncated set of classical modes: classical modes are generated from the natural nodes of the open contour, so that a truncated set contains too few modes that satisfy the coupling of its two ends. In addition, the omega stringer forms closed cells with the skin, which excludes Vlasov's null-shear warping, the stringer is offset from the skin, and the pre-buckling state is not uniform near the loaded edges.

The present study closes this gap. Its novelty is the formulation of the periodic condition within the cross-section analysis of GBT, such that every deformation mode, and therefore every combination and truncation of the modes, satisfies it exactly, together with a shear-deformable GBT–Ritz member analysis of omega-stringer-stiffened panels with symmetric and unsymmetric laminates. The formulation is derived in the notation of Silvestre and Camotim {{silvestre2013shear}}, solved by the GBT–Ritz approach {{kharghani2026jcs}} and verified against Abaqus shell models with periodic constraint equations. The mode decomposition is then used to follow the change of the buckling mode with the length of the cell, which reveals, besides the transition from local to global buckling, a long-wave mode caused by the constant width that the periodic condition imposes. Section 2 develops the formulation, Section 3 presents the results, and Section 4 concludes.
