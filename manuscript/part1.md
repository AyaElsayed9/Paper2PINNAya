::: {custom-style="MDPI_1.1_article_type"}
Article
:::

::: {custom-style="MDPI_1.2_title"}
Generalized Beam Theory with Exact Periodic Boundary Conditions for the Buckling Analysis of Omega-Stringer-Stiffened Composite Panels
:::

::: {custom-style="MDPI_1.3_authornames"}
Aya Medhat Elsayed ^1^, Navid Kharghani ^2,\*^ and Christian Mittelstedt ^2^
:::

::: {custom-style="MDPI_1.6_affiliation"}
^1^ Department of Mechanical Engineering, Technical University of Darmstadt, 64287 Darmstadt, Germany
:::

::: {custom-style="MDPI_1.6_affiliation"}
^2^ Institute for Lightweight Engineering and Structural Mechanics (LSM), Department of Mechanical Engineering, Technical University of Darmstadt, Otto-Berndt-Straße 2, 64287 Darmstadt, Germany
:::

::: {custom-style="MDPI_1.6_affiliation"}
\* Correspondence: navid.kharghani\@lsm.tu-darmstadt.de [author list and corresponding author to be confirmed by the authors]
:::

::: {custom-style="MDPI_1.7_abstract"}
**Abstract:** Stringer-stiffened composite panels are commonly analysed as one representative cell of a wide panel with many identical stringers, and the side edges of that cell must then satisfy periodic boundary conditions. Closed-form models with periodic conditions rely on idealisations of the skin–stringer interaction, and shell finite element models with periodic constraints are accurate but expensive and give no modal insight. Generalized Beam Theory (GBT) offers a mechanically transparent and efficient description of thin-walled members, but its cross-section analysis has been developed for members with free longitudinal edges or with closed contours, and periodic side edges have not been treated. This study formulates a shear-deformable GBT for the omega-stringer-stiffened composite panel in which the periodic condition is built into every family of deformation modes: the conventional (global, distortional, local and transverse-extension) modes by constraint rows that are projected out family by family, the warping functions by a constrained Vlasov least-squares problem, and the shear and transverse-shear (Mindlin) modes by merging the nodal functions of the two side edges. Any combination of modes, and hence any truncated basis, therefore satisfies the periodic condition exactly. The walls are laminated Reissner–Mindlin plates with full membrane–bending and anisotropic coupling, the stringer is tied to the skin with its true offset, the pre-buckling state is the first-order solution under the actual end loads, and the amplitude functions are cubic B-splines in a Ritz sense, which covers simply supported and clamped loaded edges. The formulation reduces exactly to the classical GBT system for uncoupled orthotropic walls under uniform compression. Against 90 Abaqus shell models (five lay-ups, two stringer pitches, three aspect ratios) the critical loads agree within 0.58 % on average (maximum 1.48 %), and within 0.57 % (maximum 1.25 %) for the 30 cases with periodic side edges; the change caused by periodic instead of symmetry side edges is reproduced with a correlation of 0.995. For panels with bending–twisting or membrane–bending coupling, periodic side edges raise the critical load by up to about 2 % for simply supported and by 5–18 % for clamped loaded edges, the latter confirmed by an independent spline finite strip model.
:::

::: {custom-style="MDPI_1.8_keywords"}
**Keywords:** generalized beam theory; periodic boundary conditions; stiffened panels; omega stringer; laminated composites; buckling; shear deformation
:::

# Introduction

Stringer-stiffened panels are the standard construction of aircraft fuselage and wing covers, and their buckling behaviour is often the governing design criterion of weight-critical laminated composite structures {{zhang2024,mittelstedt2025book}}. Among the stiffener shapes, the omega (hat) stringer is attractive because its closed cell provides high bending and torsional stiffness, a wide bonding area and good damage tolerance {{mittelstedtbeerhorst2009,schilling2020tws}}. A typical cover carries many identical stringers. Its buckling behaviour away from the free or supported longitudinal edges is therefore represented by one repeating cell, provided that the side edges of this cell are given the conditions that the neighbouring cells impose, namely periodic boundary conditions (PBC) {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}.

## Periodic boundary conditions in the analysis of stiffened panels

The treatment of the repeating cell is a classical topic of the analysis of stiffened panels. In exact and finite strip methods for prismatic plate assemblies, the infinitely wide panel was analysed early through its repeating portion: the exact strip programs VIPASA and VICONOPT treat infinitely wide or long panels by repeating a portion and its supports {{wittrick1974,kennedy1994}}, and the finite strip method in general handles such assemblies with a small number of strips {{cheung1976,lau1986}}. The underlying principle is the theory of waves in periodic structures (Bloch–Floquet theory), by which the response of an infinite periodic structure follows from one cell with phase-shifted boundary conditions {{mead1996}}.

At the Institute for Lightweight Engineering and Structural Mechanics (LSM), a series of closed-form and semi-analytical models for stringer-stiffened composite panels with periodic boundary conditions has been developed. Mittelstedt derived explicit buckling loads of blade-stringer-stiffened plates considering periodic boundary conditions {{mittelstedt2007}} and extended the approach to the stiffener–plate interaction in local buckling {{mittelstedt2009periodic}}. For omega stringers, Mittelstedt and Beerhorst presented a closed-form buckling analysis {{mittelstedtbeerhorst2009}}, and Mittelstedt and Schröder treated the local postbuckling of hat-stringer-stiffened plates {{mittelstedtschroder2010}}. Schilling and Mittelstedt developed a closed-form approximate solution for the local buckling of omega-stringer-stiffened panels {{schilling2020tws}}, a column-buckling model with periodic boundary conditions {{schilling2020rineng}}, and models for shear buckling, compression–shear interaction and local postbuckling {{schilling2021shear,schilling2022interaction,schilling2022post}}. Most recently, El Yaakoubi-Mesbah and Mittelstedt presented a closed-form analytical solution for the local buckling of omega-stringer-stiffened composite panels under compression {{elyaakoubi2025}}, a computational model for local buckling {{elyaakoubi2024wccm}} and a closed-form postbuckling analysis including the stringer-foot effects {{elyaakoubi2026}}; the stability of stiffened laminates is summarised in {{mittelstedt2025stiff}}. These models are efficient and transparent, but they rest on idealisations: the stringer is replaced by elastic restraints, the skin is reduced to one or a few plate strips with assumed shape functions, local and global buckling are usually treated by separate models, and the pre-buckling stress is taken as uniform. Membrane–bending coupling of unsymmetric laminates is typically included through a reduced bending stiffness {{ewing1988}}.

The physics-informed neural network of the authors for the local buckling of omega-stiffened panels {{elsayed2026}} used such a reduced one-dimensional energy model and three-dimensional Abaqus analyses as the reference. Shell finite element models with periodic constraint equations between the nodes of the two side edges are accurate and general, but each design requires a new mesh, a large eigenvalue analysis and, in the case of coupled laminates, care in the definition of the constraint set; they also provide the buckling mode only as a nodal field, without a decomposition into mechanically meaningful deformation mechanisms.

## Generalized Beam Theory

Generalized Beam Theory (GBT), originated by Schardt {{schardt1989,schardt1994}} and developed further by Davies and Leach {{davies1994}}, describes the deformation of a prismatic thin-walled member as a superposition of cross-section deformation modes whose amplitudes vary along the member axis. A cross-section analysis determines the modes and the associated modal matrices, and a member analysis determines the amplitude functions. The modal decomposition separates global, distortional, local and shear-related deformations and has made GBT one of the most transparent tools for coupled instabilities of thin-walled members, as summarised in the state-of-the-art reports by Camotim et al. {{camotim2010}} and Mittelstedt {{mittelstedt2020amr,mittelstedt2024gbt}}. Silvestre and Camotim formulated first- and second-order GBT for arbitrary orthotropic materials {{silvestre2002a,silvestre2002b}}, and later shear-deformable GBT for composite members and its vibration counterpart {{silvestre2013shear,silvestre2013vib}}. Dinis et al. extended the cross-section analysis to arbitrarily branched open sections {{dinis2006}}. Gonçalves et al. introduced a new approach to the calculation of the deformation modes {{goncalves2010}}, and Bebiano et al. rationalised the cross-section analysis into natural, shear and transverse-extension mode families {{bebiano2015}}, together with the shear deformation modes {{goncalves2014shear}}, implemented in the code GBTUL 2.0 {{bebiano2018}}. Closed sections were treated, among others, for circular cylindrical shells and tubes {{silvestre2007cyl}} and for regular polygonal tubes {{goncalves2013polygon}}, where the closed contour leads to Fourier-type modes along the circumference. For laminated composites, GBT has been applied to FRP open-section columns {{silva2010}} and to laminated CFRP cylindrical stiffened panels {{silva2010ijssd}}, and Taig and Ranzi formulated GBT for stiffened sections {{taig2014}}.

Kharghani and Mittelstedt combined GBT with the Ritz method into a GBT–Ritz approach: for global, local and distortional buckling of shear-deformable composite beams with open cross-sections {{kharghani2025jcs}}, for composite box girders {{kharghani2025semc}}, and for thin-walled composite beams with unsymmetric layups, for which membrane–bending coupling is consistently included {{kharghani2026cepa,kharghani2026jcs}}. In the GBT–Ritz approach, GBT supplies the cross-section kinematics and the modal matrices, and the Ritz method supplies the approximation of the amplitude functions along the member, with admissible functions that satisfy the essential boundary conditions {{mittelstedt2025ritz,vescovini2018}}. The present work follows this route.

## Research gaps and motivation

Two gaps motivate the present study.

*The gap in GBT.* GBT has been developed for members: open sections with free longitudinal edges, branched sections, and closed contours. The skin of a stiffened panel cell, however, ends at two side edges that are neither free nor closed: the displacement field at one side edge must equal that at the other. This condition has not been formulated in GBT. It cannot simply be imposed afterwards on a set of classical modes. Classical mode generation starts from a few natural nodes, and the condition couples the two ends of the open contour; a truncated set of classical modes then contains too few modes that satisfy it, and the rigid-body and warping modes of the open cell are not admissible for the periodic cell. Moreover, the omega stringer forms closed cells with the skin, so that Vlasov's null membrane shear condition cannot be fulfilled by the warping functions, the stringer is offset from the skin by the half thicknesses, and the pre-buckling state is not uniform near the loaded edges, where the load enters the stringer through the eccentric foot.

*The gap in the treatment of periodic conditions.* The closed-form models with periodic conditions summarised above are efficient, but they rely on idealisations of the skin–stringer interaction, separate local and global buckling, and do not account for the non-uniform pre-buckling state. Shell finite element models with periodic constraints do not have these restrictions but are expensive for parametric studies and give no modal insight. Exact and finite strip methods treat the repeating portion efficiently, but they are formulated in nodal degrees of freedom and do not provide the decomposition into global, distortional and local mechanisms.

*Why combine them.* GBT and periodic boundary conditions complement each other. The periodic cell is prismatic, so the member analysis of GBT applies directly; the modal decomposition identifies the buckling mechanisms of the cell; and the reduced number of generalised coordinates makes parametric studies cheap. Conversely, the periodic condition removes the arbitrariness of the side-edge modelling that otherwise dominates panel models built from members. The key requirement is that the condition is satisfied by the modes themselves, so that it holds for any combination and any truncation of the basis.

## Contributions

The contributions of this study are:

1. A cross-section analysis in which the periodic condition is built into every family of GBT deformation modes: conventional modes (global, distortional, local, transverse extension) through constraint rows that are projected out family by family, the warping of the global and distortional modes through a Vlasov least-squares problem with a periodicity constraint, and the shear and transverse-shear modes through merging of the nodal functions of the two side edges (Section 3.4). Every mode, and therefore every truncated basis, satisfies the periodic condition exactly.
2. A shear-deformable GBT with laminated Reissner–Mindlin walls including all membrane–bending and anisotropic couplings, transverse-extension and shear modes, a stringer–skin tie with the true offset, a first-order pre-buckling state under the actual end loads with all membrane resultants, and amplitude functions approximated by cubic B-splines in a Ritz sense, covering simply supported (SS) and clamped (CC) loaded edges (Sections 3.2–3.8).
3. The exact reduction of the formulation to the classical GBT system and to the classical geometric matrices for uncoupled orthotropic walls under uniform compression (Section 3.9).
4. A validation against 90 Abaqus shell models with symmetry and periodic side edges and an independent spline finite strip model, including the first systematic quantification of the effect of periodic instead of symmetry side edges on omega-stiffened panels with coupled laminates (Section 4).

The paper is organised as follows. Section 2 defines the panel and the load. Section 3 develops the formulation: cross-section discretisation, constitutive law, kinematics, deformation modes, periodic conditions, energy and modal matrices, amplitude functions, pre-buckling state and eigenvalue problem, and the classical limits. Section 4 presents the verification, and Section 5 the conclusions.

# Problem definition

## Geometry and coordinates

The panel is the periodic cell of an infinitely wide panel with identical omega stringers at the pitch $b_s$ (Figure 1a). The cell contains two stringers and has the width $W = 2b_s$ and the length $L = \beta b_s$, with the aspect ratio $\beta$. The global axes are $X$ across the panel, $Y$ normal to the skin towards the stringers, and $x$ along the stringers, which is the load direction and the longitudinal axis of GBT. The skin occupies $-W \le X \le 0$, $Y = 0$, $0 \le x \le L$; the stringer centre lines are at $X = -b_s/2$ and $X = -3b_s/2$, so that the side edges $X = -W$ and $X = 0$ lie in the middle of a skin bay.

The stringer geometry is fixed: web root half-width $b_3 = 27.75$ mm, web height $h = 30$ mm, web inclination $\theta = 25^\circ$ from the vertical, crown half-width $h_t = b_3 - h\tan\theta = 13.76$ mm, and foot (flange) width $b_1 = 0.1\,(b_s - 2b_3)$. The skin thickness is $t_s$ and the stringer thickness $t_f$. The stringer mid-surface is offset from the skin mid-surface by $c = (t_s+t_f)/2$, so that the web roots lie at $Y = c$ and the crown at $Y = h + c$. In each wall, $s$ is the coordinate along the wall mid-line, $z$ the coordinate through the thickness, and $u$, $v$, $w$ are the displacements in the directions $x$, $s$ and $z$, as in classical GBT (Figure 1b).

![**Figure 1.** (a) Infinitely wide omega-stringer-stiffened panel and its periodic cell of width $W = 2b_s$. (b) GBT discretisation of the cell for $b_s = 150$ mm with natural and intermediate nodes; skin, feet (skin and flange as one laminate), webs, crowns and link walls of the offset tie; $s_L$ and $s_R$ are the side nodes on which the periodic condition acts.](fig/fig1_geometry.png){width=15.5cm}

## Materials, lay-ups and load

Two materials are considered: an isotropic material ($E = 200$ GPa, $\nu = 0.3$, $t_s = 1.0$ mm, $t_f = 0.875$ mm) and a carbon-fibre reinforced plastic (CFRP) ply with $E_1 = 150$ GPa, $E_2 = 9$ GPa, $G_{12} = G_{13} = G_{23} = 5.3$ GPa, $\nu_{12} = 0.33$ and the ply thickness $t_{ply} = 0.125$ mm. The skin has eight and the stringer seven plies (Table 1). The lay-ups include a cross-ply laminate without coupling, two symmetric angle-ply laminates with bending–twisting coupling ($D_{16}, D_{26} \neq 0$, and in-plane shear–extension coupling $A_{16}, A_{26} \neq 0$ where the laminate is not balanced), and an unsymmetric laminate with membrane–bending coupling ($\mathbf B \neq \mathbf 0$). Ply angles are measured from the stringer axis $x$.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 1.** Lay-ups of skin (8 plies) and stringer (7 plies), ply thickness 0.125 mm; ply angles in degrees from the $x$-axis, subscript s: symmetric, subscript 2: repeated twice.
:::

| Name | Skin | Stringer | Coupling |
|:--|:--|:--|:--|
| Cross-ply | [0/90/0/90]$_\mathrm{s}$ | [0/90/0/90/0/90/0] | none |
| Sym1 | [45/−45/−45/45]$_\mathrm{s}$ | [45/−45/45/−45/45/−45/45] | $D_{16}, D_{26}$ (stringer also $A_{16}, A_{26}$) |
| Sym2 | [60/−30/−30/60]$_\mathrm{s}$ | [60/−30/60/−30/60/−30/60] | $A_{16}, A_{26}, D_{16}, D_{26}$ |
| Unsym | [0/30/60/90]$_2$ | [0/30/60/90/0/30/60] | $\mathbf A, \mathbf B, \mathbf D$ fully populated |

A uniform compressive line load of 1 N/mm acts on every loaded edge at $x = 0$ and $x = L$: on the skin, the webs and the crowns, and 2 N/mm on a foot, which carries skin and flange; the foot load acts at the height of its resultant, $Y = c/2$. The eigenvalue $\lambda$ is therefore the critical line load $N_{cr} = \lambda \times 1$ N/mm, the quantity reported by the Abaqus buckling step. The design points are $b_s \in \{150, 285.5\}$ mm and $\beta \in \{0.75, 1.5, 3\}$, with SS or CC loaded edges and symmetry or periodic side edges.
