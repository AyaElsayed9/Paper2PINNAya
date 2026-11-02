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
\* Correspondence: navid.kharghani\@lsm.tu-darmstadt.de
:::

::: {custom-style="MDPI_1.6_affiliation"}
A.M.E.: ayaallah.elsayed\@stud.tu-darmstadt.de; C.M.: christian.mittelstedt\@lsm.tu-darmstadt.de
:::

::: {custom-style="MDPI_1.7_abstract"}
**Abstract**
:::

::: {custom-style="MDPI_1.7_abstract"}
Stringer-stiffened composite panels are commonly analysed on one representative cell with periodic boundary conditions. Generalized Beam Theory (GBT) describes the buckling of thin-walled members by cross-section modes, but periodic side edges have not been treated in its cross-section analysis. This paper builds the periodic condition into every family of GBT modes, so that every mode and every truncated basis satisfies it exactly. A shear-deformable GBT–Ritz member analysis of omega-stringer-stiffened panels follows, with laminated Reissner–Mindlin walls, the true stringer offset, a computed pre-buckling state and B-spline amplitude functions; it reduces to classical GBT for uncoupled walls. For 30 cells with simply supported loaded edges, the critical loads agree with periodic Abaqus shell models within 0.57 % on average and 1.25 % at most; with clamped loaded edges within 0.64 % and 1.52 %. Short cells buckle locally. From about six stringer pitches, isotropic and angle-ply cells buckle in a long-wave mode driven by the transverse compression that the constant cell width induces, whereas cross-ply and unsymmetric cells buckle globally from twelve pitches. Abaqus confirms this mode and, for two cells allowed to expand laterally, gives 2.0–2.9 times higher critical loads; the lateral restraint assumed in periodic analyses thus matters for long panels.
:::

::: {custom-style="MDPI_1.8_keywords"}
**Keywords:** generalized beam theory; Ritz method; periodic boundary conditions; stiffened panels; omega stringer; laminated composites; buckling
:::

# Introduction

Stringer-stiffened panels are the standard construction of aircraft fuselage and wing covers, and buckling is often their governing design criterion {{zhang2024,mittelstedt2025book}}. The omega (hat) stringer offers high bending and torsional stiffness and a wide bonding area {{mittelstedtbeerhorst2009,schilling2020tws}}. Away from its longitudinal edges, a cover with many identical stringers is represented by one repeating cell with periodic boundary conditions {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. The principle goes back to waves in periodic structures {{mead1996}} and is used in exact strip programs (VIPASA, VICONOPT) {{wittrick1974,kennedy1994}} and finite strip methods {{cheung1976,lau1986}}.

At the Institute for Lightweight Engineering and Structural Mechanics (LSM), closed-form models of stiffened composite panels have been developed for blade stringers {{mittelstedt2007,mittelstedt2009periodic}} and omega stringers {{mittelstedtbeerhorst2009,mittelstedtschroder2010}}, including local and column buckling {{schilling2020tws,schilling2020rineng}}, shear and combined loading {{schilling2021shear,schilling2022interaction}} and postbuckling {{schilling2022post}}, and most recently local buckling and postbuckling of omega-stringer panels including the stringer feet {{elyaakoubi2025,elyaakoubi2024wccm,elyaakoubi2026}}; see also {{mittelstedt2025stiff}}. These models reduce the panel to a repeating unit with periodic {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}} or symmetry conditions {{elyaakoubi2024wccm}}, the latter admitting only modes symmetric about the symmetry lines. They are efficient but idealise the skin–stringer interaction by elastic restraints, separate local and global buckling, assume a uniform pre-buckling stress and include membrane–bending coupling through a reduced bending stiffness {{ewing1988}}, as does the physics-informed neural network of the authors {{elsayed2026}}. Shell finite element models avoid these idealisations but are expensive and give the mode only as a nodal field.

Generalized Beam Theory (GBT), originated by Schardt {{schardt1989,schardt1994}} and developed further by Davies and Leach {{davies1994}}, describes a thin-walled member by cross-section modes with amplitude functions along the axis, which makes it transparent for coupled instabilities {{camotim2010,mittelstedt2020amr,mittelstedt2024gbt,mittelstedt2025gbtch}}. Silvestre and Camotim formulated GBT for orthotropic and laminated members {{silvestre2002a,silvestre2002b}} including shear deformation {{silvestre2013shear,silvestre2013vib}}. The cross-section analysis was extended to branched sections {{dinis2006}}, reformulated {{goncalves2010,bebiano2015}} with shear modes {{goncalves2014shear}} and implemented in GBTUL 2.0 {{bebiano2018}}; closed sections were treated for tubes {{silvestre2007cyl,goncalves2013polygon}}, and GBT was applied to fibre-reinforced polymer (FRP) columns {{silva2010}}, stiffened composite panels {{silva2010ijssd}} and stiffened sections {{taig2014}}. The Ritz method complements GBT in the member analysis {{mittelstedt2025ritz,vescovini2018,kuehn2014}}: Kharghani and Mittelstedt combined both into a GBT–Ritz approach for global, local and distortional buckling of shear-deformable composite beams {{kharghani2025jcs}}, composite box girders {{kharghani2025semc}} and beams with unsymmetric lay-ups {{kharghani2026cepa,kharghani2026jcs}}.

GBT has been developed for members with free longitudinal edges or closed contours, whereas the skin of a panel cell ends at two side edges whose displacements must be equal. This condition has not been formulated in GBT, and it cannot simply be imposed on a truncated set of classical modes: these are computed for an open contour with independent side edges, so that few of their combinations are periodic. Moreover, the omega stringer forms closed cells with the skin, which excludes Vlasov's warping, it is offset from the skin, and the pre-buckling state is not uniform.

The present study closes this gap. Its novelty is the formulation of the periodic condition within the cross-section analysis of GBT, such that every mode, and therefore every truncated basis, satisfies it exactly, together with a shear-deformable GBT–Ritz analysis {{kharghani2026jcs}} of omega-stringer-stiffened panels with symmetric and unsymmetric laminates, verified against Abaqus shell models. The mode decomposition reveals, besides the transition from local to global buckling of long cells, a long-wave mode caused by the constant width that the periodic condition imposes. Since every mode is periodic by construction, the formulation is also suited as the physics model of physics-informed neural networks for such panels, extending the authors' earlier network {{elsayed2026}} (Section 4). Section 2 derives the formulation, Section 3 verifies it against Abaqus and discusses the results, and Section 4 concludes.
