
# Theory and formulations

GBT describes the deformation of a prismatic thin-walled member by cross-section deformation modes whose amplitudes vary along the member axis. The analysis is divided into two stages {{schardt1994,silvestre2002a,bebiano2015}}: a cross-section analysis, which determines the deformation modes and the modal matrices, and a member analysis, which determines the amplitude functions. Based on the strain–displacement relations, the strain energy of the cross-section is established in terms of the modal amplitudes, and the governing equations follow from the stationarity of the total potential energy, which is evaluated here by the Ritz method, as in the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}. In the present work the cross-section analysis is extended by the periodic boundary conditions of the panel cell (Section 2.2).

## Cross-section analysis

The panel is the periodic cell of an infinitely wide panel with identical omega stringers at the pitch $b_s$ (Figure 1a). The cell contains two stringers and has the width $W = 2b_s$ and the length $L = \beta b_s$, with the aspect ratio $\beta$. The global cross-section axes are $X$ across the panel and $Y$ normal to the skin towards the stringers; $x$ is the longitudinal axis along the stringers. The skin occupies $-W \le X \le 0$, and its side edges $X = -W$ and $X = 0$ lie in the middle of a skin bay. The mid-line of the cross-section is divided into walls (skin, feet, webs, crowns) and the walls into segments of length at most $b_s/25$, with at least two segments per wall. As in classical GBT, the end points of the walls are the natural nodes and the remaining segment end points are intermediate nodes (Figure 1b) {{schardt1994,bebiano2015}}. For $b_s = 150$ mm the cell has 92 nodes, 19 walls and 93 segments. The coordinate $s$ runs along the mid-line of every wall, $z$ is the coordinate through the thickness, and $u$, $v$ and $w$ are the displacements in the directions $x$, $s$ and $z$. The two side nodes of the skin are denoted $s_L$ ($X = -W$) and $s_R$ ($X = 0$).

![**Figure 1.** (a) Infinitely wide omega-stringer-stiffened panel and its periodic cell of width $W = 2b_s$. (b) GBT discretisation of the cell for $b_s = 150$ mm with natural and intermediate nodes; skin, feet (skin and flange as one laminate), webs, crowns and link walls of the offset tie; $s_L$ and $s_R$ are the side nodes on which the periodic condition acts.](fig/fig1_geometry.png){width=15.5cm}

The analysis begins with the displacement field of shear-deformable GBT {{silvestre2013shear,kharghani2025jcs}},

$$
u(x,s) = \sum_{i}\bar u_i(s)\,\varphi_{i,x}(x) + \sum_{j}\bar u_j(s)\,\varphi_j(x), \eqn{u}
$$

$$
v(x,s) = \sum_{i}\bar v_i(s)\,\varphi_i(x), \eqn{v}
$$

$$
w(x,s) = \sum_{i}\bar w_i(s)\,\varphi_i(x), \eqn{w}
$$

which is extended here by the transverse shear rotations of the Reissner–Mindlin walls,

$$
\beta_x(x,s) = \sum_{r}\bar c_r(s)\,\chi_r(x),\qquad \beta_s(x,s) = \sum_{r}\bar b_r(s)\,\vartheta_r(x). \eqn{beta}
$$

The index $i$ refers to the conventional deformation modes, $j$ to the shear modes and $r$ to the transverse shear (Mindlin) modes. The functions $\bar u_i(s)$, $\bar v_i(s)$ and $\bar w_i(s)$ describe the displacements of a conventional mode in the axial, tangential and out-of-plane directions, $\bar u_j(s)$ is the warping of a shear mode, and $\bar c_r(s)$ and $\bar b_r(s)$ are the distributions of the shear rotations. The corresponding amplitude functions are $\varphi_i(x)$, $\varphi_j(x)$, $\chi_r(x)$ and $\vartheta_r(x)$. The uniform axial displacement is included among the shear modes, with $\bar u_j = 1$.

The rotations of the wall normal are $\phi_x = -w_{,x}+\beta_x$ and $\phi_s = -w_{,s}+\beta_s$, so that the strains at a distance $z$ from the mid-surface are

$$
\varepsilon_{xx} = u_{,x} + z\,\kappa_{xx},\qquad \varepsilon_{ss} = v_{,s} + z\,\kappa_{ss},\qquad \gamma_{xs} = u_{,s} + v_{,x} + z\,\kappa_{xs}, \eqn{strains}
$$

with the curvatures and the transverse shear strains

$$
\begin{gathered}
\kappa_{xx} = -w_{,xx} + \beta_{x,x},\qquad \kappa_{ss} = -w_{,ss} + \beta_{s,s},\qquad \kappa_{xs} = -2w_{,xs} + \beta_{x,s} + \beta_{s,x},\\
\gamma_{xz} = w_{,x} + \phi_x = \beta_x,\qquad \gamma_{sz} = w_{,s} + \phi_s = \beta_s .
\end{gathered} \eqn{curv}
$$

For $\beta_x = \beta_s = 0$, Eqs. (@strains) and (@curv) are the strain–displacement relations of the GBT–Ritz approach for thin walls, $\varepsilon_{xx} = u_{,x} - z\,w_{,xx}$, $\varepsilon_{ss} = v_{,s} - z\,w_{,ss}$ and $\gamma_{xs} = u_{,s} + v_{,x} - 2z\,w_{,xs}$ {{kharghani2026jcs}}.

In every ply $p$, Hooke's law relates the in-plane stresses to the strains by the transformed reduced stiffness matrix $\bar{\mathbf Q}^{(p)}$ ($i,j = 1,2,6$), and the transverse shear stresses to $\gamma_{xz}$, $\gamma_{sz}$ by $\bar{\mathbf Q}^{(p)}_s$. Integration through the thickness gives the laminate stiffnesses {{mindlin1951,whitney1970,reddy2004}}

$$
\left(A_{ij}, B_{ij}, D_{ij}\right) = \sum_{p=1}^{N_L}\int_{z_{p}}^{z_{p+1}}\bar Q^{(p)}_{ij}\left(1, z, z^2\right)dz,\qquad
\mathbf A_s = \frac56\sum_{p=1}^{N_L}\bar{\mathbf Q}^{(p)}_s\left(z_{p+1}-z_p\right), \eqn{abd}
$$

with the number of plies $N_L$, and the constitutive law of a wall reads

$$
\hat{\boldsymbol\sigma} = \mathbf C\,\hat{\boldsymbol\varepsilon},\qquad
\mathbf C = \begin{bmatrix} \mathbf A & \mathbf B & \mathbf 0\\ \mathbf B & \mathbf D & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \mathbf A_s\end{bmatrix}, \eqn{claw}
$$

with the stress resultants $\hat{\boldsymbol\sigma} = [N_{xx}, N_{ss}, N_{xs}, M_{xx}, M_{ss}, M_{xs}, Q_x, Q_s]^{\mathsf T}$ and the generalised strains $\hat{\boldsymbol\varepsilon} = [u_{,x}, v_{,s}, u_{,s}+v_{,x}, \kappa_{xx}, \kappa_{ss}, \kappa_{xs}, \gamma_{xz}, \gamma_{sz}]^{\mathsf T}$. All entries of $\mathbf A$, $\mathbf B$ and $\mathbf D$ are retained, including the membrane–bending coupling $\mathbf B$ of unsymmetric laminates and the anisotropic terms $A_{16}$, $A_{26}$, $D_{16}$, $D_{26}$. The transverse membrane resultant $N_{ss}$ is not condensed out, because the transverse-extension modes carry the strain $v_{,s}$ explicitly.

A stringer foot consists of the skin and the stringer flange, which are tied in the shell model. In GBT the two layers form one wall referred to the skin mid-plane. With $\mathbf A_0$, $\mathbf B_0$, $\mathbf D_0$ of the combined stack about its own mid-plane and $z_f$ the distance from the skin mid-plane to the mid-plane of the stack, the parallel-axis transformation gives

$$
\mathbf A_f = \mathbf A_0,\qquad \mathbf B_f = \mathbf B_0 + z_f\,\mathbf A_0,\qquad \mathbf D_f = \mathbf D_0 + 2z_f\,\mathbf B_0 + z_f^2\,\mathbf A_0 . \eqn{foot}
$$

In the shell model, the stringer mid-surface is tied to the skin with the offset $c = (t_s+t_f)/2$, where $t_s$ and $t_f$ are the thicknesses of skin and stringer. In GBT every web root is a node at $Y = c$, joined to the skin node below it by a link wall of height $c$. The link wall must transmit the in-plane displacements and the rotation like a line of rigid offset links, but it must not add axial or twisting stiffness. Its laminate is built from the skin stiffnesses $\mathbf A_0$, $\mathbf D_0$ and $\mathbf A_{s0}$ as

$$
\begin{gathered}
A^{L}_{11} = 10^{-3}A_{0,11},\quad A^{L}_{22} = k_L A_{0,22},\quad A^{L}_{66} = k_L A_{0,66},\\
D^{L}_{11} = 10^{-3}D_{0,11},\quad D^{L}_{22} = k_L D_{0,22},\quad D^{L}_{66} = 10^{-3}D_{0,66},\quad \mathbf A^{L}_s = k_L\mathbf A_{s0},
\end{gathered} \eqn{link}
$$

with $k_L = 10^3$ and all other entries zero; link walls carry no load. The results do not depend on $k_L$, whereas a link wall with a stiff twisting term $D^{L}_{66}$ would lock the rotation of the web root and over-predict the critical load (Section 3.5); omitting it is the consistent representation of a row of rigid links.

Within every segment $p$ of length $b_p$ between the nodes $q$ and $q'$, with $\eta = (s-S_q)/b_p \in [0,1]$, the functions $\bar u$ and $\bar v$ vary linearly, so that they are described by their nodal values. The transverse displacement $\bar w$ is the cubic Hermite interpolation of the nodal normal displacements $\bar w_q$, $\bar w_{q'}$ and the nodal rotations $\theta_q$, $\theta_{q'}$ ($\theta = \bar w_{,s}$):

$$
\begin{gathered}
\bar w(\eta) = H_1(\eta)\,\bar w_q + b_p H_2(\eta)\,\theta_q + H_3(\eta)\,\bar w_{q'} + b_p H_4(\eta)\,\theta_{q'},\\
H_1 = 1-3\eta^2+2\eta^3,\quad H_2 = \eta-2\eta^2+\eta^3,\quad H_3 = 3\eta^2-2\eta^3,\quad H_4 = -\eta^2+\eta^3 .
\end{gathered} \eqn{hermite}
$$

The rotation is shared by all walls that meet at a node (rigid fold lines). The in-plane displacement of node $q$ is the vector $\mathbf U_q = (U_{X,q}, U_{Y,q})$; on a segment with unit tangent $\mathbf t$ and unit normal $\mathbf n$, $\bar v = \mathbf t\cdot\mathbf U$ and $\bar w = \mathbf n\cdot\mathbf U$ at its end nodes. The space spanned by Eq. (@hermite) is that of the classical description by nodal transverse moments, which the GBT–Ritz approach uses with the segment bending stiffness $D_{22}$ {{kharghani2025jcs,kharghani2026jcs}}; the nodal form of Eq. (@hermite) is used here because it turns the periodic conditions of Section 2.2 into linear relations between nodal quantities.

Let $N$ be the number of nodes. The in-plane nodal displacements $\mathbf U = [\mathbf U_1;\dots;\mathbf U_N] \in \mathbb R^{2N}$ and the nodal rotations $\boldsymbol\Theta \in \mathbb R^{N}$ describe the in-plane deformation of the cross-section. With the matrix $\mathbf E_\ell$ that gives the elongation of every segment, $(\mathbf E_\ell\mathbf U)_p = \mathbf t_p\cdot(\mathbf U_{q'}-\mathbf U_q)$, and following the organisation of the mode families in {{bebiano2015,goncalves2014shear}}, the conventional modes are defined by the families

* **G (global):** rigid in-plane translations and rotation of the cross-section, as far as the side conditions admit them;
* **D (distortional):** in-plane nodal displacements without segment elongation, $\mathbf E_\ell\mathbf U = \mathbf 0$ (Vlasov's null transverse membrane strain), orthogonal to G and L;
* **L (local):** a unit normal displacement of one node that is not a wall junction, and a pure rotation of one node (the in-plane fields of G, D and T carry no nodal rotation, so that the nodal rotations are independent coordinates);
* **T (transverse extension):** segment elongation patterns, the columns of $\mathbf E_\ell^{\mathsf T}$;

and the further modes by

* **S (shear):** the nodal hat functions $h_q(s)$ of the warping $\bar u_j$, one of which is replaced by the uniform axial mode $\bar u_j = 1$;
* **M (Mindlin):** for $\beta_x$, a nodal vector $\boldsymbol\Omega_q$ in the plane of the cross-section, with $\bar c_r = \mathbf t\cdot\boldsymbol\Omega_q$ interpolated linearly on every segment (one vector at a node inside a straight wall, two at a fold line); for $\beta_s$, one nodal value $\bar b_r$ per node, interpolated linearly.

Since $\mathbb R^{2N} = \ker\mathbf E_\ell \oplus \operatorname{range}\mathbf E_\ell^{\mathsf T}$, the families G, D, L and T together span all in-plane nodal displacements, and the rotation modes make the nodal rotations independent. The conventional modes are therefore a complete nodal description of the in-plane deformation of the cross-section, organised in the GBT families.

The warping $\bar u_i$ of the G and D modes follows in classical GBT from Vlasov's condition of null membrane shear strain, $\bar u_{i,s} + \bar v_i = 0$ on every wall {{vlasov1961,schardt1994}}. The omega stringer and the skin form closed cells, in which this condition cannot be satisfied exactly. The warping is therefore determined as the stiffness-weighted least-squares solution

$$
\min_{\bar{\mathbf u}}\ \sum_p \hat A_{66,p}\,b_p\left(\frac{\bar u_{q'}-\bar u_{q}}{b_p}+\bar v_p\big|_{\eta = 1/2}\right)^2
\quad\text{subject to}\quad \mathbf c_{ax}^{\mathsf T}\bar{\mathbf u} = 0, \eqn{vlasov}
$$

where the sum runs over all segments, $\hat A_{66} = A_{66}-A_{26}^2/A_{22}$ is the membrane shear stiffness for $N_{ss} = 0$, and $c_{ax,q} = \sum_{p\ni q}\tfrac12\hat A_{11,p}\,b_p$ with $\hat A_{11} = A_{11}-A_{12}^2/A_{22}$, so that the stiffness-weighted mean warping vanishes (no axial content). For an open section, Eq. (@vlasov) reproduces the classical warping with zero residual; in the closed cells, the remaining membrane shear strain is carried by the shear modes S.

Inserting Eqs. (@u)–(@beta) into Eqs. (@strains) and (@curv), the generalised strains of a mode become products of a cross-section function and an amplitude function or one of its derivatives (Table 1). The conventional entries are those of classical GBT, the entry $v_{,s} = \bar v_{i,s}\varphi_i$ belongs to the transverse extension, and the entries of the shear and Mindlin modes are those of shear-deformable GBT.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 1.** Non-zero generalised strains of one mode, multiplied by the amplitude function ($\varphi$), its first derivative ($\varphi^{I}$) or its second derivative ($\varphi^{II}$).
:::

| Mode (amplitude) | $\varphi$ | $\varphi^{I}$ | $\varphi^{II}$ |
|:--|:--|:--|:--|
| conventional ($\varphi_i$) | $v_{,s} = \bar v_{i,s}$, $\kappa_{ss} = -\bar w_{i,ss}$ | $\gamma_{xs} = \bar u_{i,s}+\bar v_i$, $\kappa_{xs} = -2\bar w_{i,s}$ | $u_{,x} = \bar u_i$, $\kappa_{xx} = -\bar w_i$ |
| shear ($\varphi_j$) | $\gamma_{xs} = \bar u_{j,s}$ | $u_{,x} = \bar u_j$ | – |
| Mindlin $\beta_x$ ($\chi_r$) | $\kappa_{xs} = \bar c_{r,s}$, $\gamma_{xz} = \bar c_r$ | $\kappa_{xx} = \bar c_r$ | – |
| Mindlin $\beta_s$ ($\vartheta_r$) | $\kappa_{ss} = \bar b_{r,s}$, $\gamma_{sz} = \bar b_r$ | $\kappa_{xs} = \bar b_r$ | – |

## Periodic boundary conditions

The cell is one period of an infinitely wide panel. The displacement field at the side edge $X = -W$ must therefore equal that at $X = 0$, for all components and at every $x$:

$$
\mathcal U_\alpha(X=-W,\,x) = \mathcal U_\alpha(X=0,\,x),\qquad \alpha = 1,\dots,6,\qquad 0\le x\le L, \eqn{pbc}
$$

where $\mathcal U_\alpha$ are the six degrees of freedom of a shell node, three translations and three rotations. This is the condition that a shell finite element model imposes with one constraint equation per pair of side nodes and degree of freedom {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. On the skin side edges, whose normal is $Y$, the degrees of freedom correspond to the GBT quantities in Table 2. Since all of them are linear in the amplitude functions, Eq. (@pbc) holds for every $x$ and every set of amplitude functions if, and only if, every mode satisfies it. The periodic condition is therefore built into the modes, family by family.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 2.** Degrees of freedom of a skin side node and the corresponding conditions in GBT.
:::

| Shell DOF | Meaning at a skin side node | Periodic condition in GBT |
|:--|:--|:--|
| $\mathcal U_1$ | in-plane displacement $U_X$ | $U_{X,s_L} = U_{X,s_R}$, Eq. (@rows) |
| $\mathcal U_2$ | deflection $U_Y$ | $U_{Y,s_L} = U_{Y,s_R}$, Eq. (@rows) |
| $\mathcal U_3$ | axial displacement $u$ | warping of G and D modes, Eq. (@kkt); shear modes, Eq. (@shearmerge); axial mode uniform |
| $\mathcal U_4$ | rotation $\phi_x = -w_{,x}+\beta_x$ | $w_{,x}$ periodic because $U_Y$ is; $\boldsymbol\Omega_{s_L} = \boldsymbol\Omega_{s_R}$, Eq. (@mindlinmerge) |
| $\mathcal U_6$ | rotation $\theta+\beta_s$ | $\theta_{s_L} = \theta_{s_R}$, Eq. (@rows); $\bar b_{s_L} = \bar b_{s_R}$, Eq. (@mindlinmerge) |
| $\mathcal U_5$ | drilling rotation | not a degree of freedom of the plate theory |

For the conventional modes, the conditions $U_{X,s_L} = U_{X,s_R}$, $U_{Y,s_L} = U_{Y,s_R}$ and $\theta_{s_L} = \theta_{s_R}$ are the rows of a constraint matrix acting on the stacked nodal vector,

$$
\mathbf R\begin{bmatrix}\mathbf U\\ \boldsymbol\Theta\end{bmatrix} = \mathbf 0,\qquad
\mathbf R = \begin{bmatrix}\mathbf e^{\mathsf T}_{X,s_L}-\mathbf e^{\mathsf T}_{X,s_R}\\ \mathbf e^{\mathsf T}_{Y,s_L}-\mathbf e^{\mathsf T}_{Y,s_R}\\ \mathbf e^{\mathsf T}_{\theta,s_L}-\mathbf e^{\mathsf T}_{\theta,s_R}\end{bmatrix}\in\mathbb R^{3\times 3N}, \eqn{rows}
$$

where $\mathbf e_{X,q}$, $\mathbf e_{Y,q}$ and $\mathbf e_{\theta,q}$ extract $U_{X,q}$, $U_{Y,q}$ and $\theta_q$. Let $\mathbf V_g \in \mathbb R^{3N\times n_g}$ hold the nodal fields $[\mathbf U;\boldsymbol\Theta]$ of the $n_g$ modes of family $g \in \{G, D, L, T\}$ of the unconstrained cell. The constraint is imposed family by family, so that the mechanical meaning of the families is kept:

1. the columns of $\mathbf V_g$ that satisfy $\mathbf R\mathbf v = \mathbf 0$ are kept unchanged;
2. the remaining columns are replaced by a basis of their constrained combinations, $\mathbf V_g \leftarrow \mathbf V_g\mathbf Z_g$ with $\mathbf Z_g = \ker(\mathbf R\mathbf V_g)$;
3. the vectors of $\ker\left(\mathbf R\,[\mathbf V_G\ \mathbf V_D\ \mathbf V_L\ \mathbf V_T]\right)$ that are not yet spanned (combinations of several families) are added, to the family T if they elongate a segment ($\mathbf E_\ell\mathbf U \neq \mathbf 0$) and to the family D otherwise;
4. each column is scaled to a unit maximum nodal displacement or rotation.

The resulting modal matrix $\mathbf V = [\mathbf V_G\ \mathbf V_D\ \mathbf V_L\ \mathbf V_T]$ satisfies $\mathbf R\mathbf V = \mathbf 0$ and spans the whole periodic subspace $\ker\mathbf R$. Consequently, for arbitrary amplitude functions,

$$
\mathbf R\sum_{i}\begin{bmatrix}\mathbf U_i\\ \boldsymbol\Theta_i\end{bmatrix}\varphi_i(x) = \mathbf 0\qquad\text{for all } x, \eqn{rowsall}
$$

and the same holds for the derivatives of the amplitude functions and for any subset of the modes. The rigid rotation of the cross-section about the $x$-axis moves the two side nodes by different amounts ($U_{Y,s_L}-U_{Y,s_R} = -W\omega$ for a rotation $\omega$) and is eliminated: a torsional rigid-body motion of one cell is not compatible with an infinitely wide panel. The two rigid translations remain. On the other hand, the side nodes may move laterally and rotate, provided that they do so identically, which yields 13 distortional modes (Table 3).

For the warping of the G and D modes, the periodic condition $\bar u_{s_L} = \bar u_{s_R}$ is added to the constraint of Eq. (@vlasov). With the segment difference matrix $\mathbf H$, $(\mathbf H\bar{\mathbf u})_p = \bar u_{q'}-\bar u_{q}$, the weight matrix $\mathbf W = \operatorname{diag}(\hat A_{66,p}/b_p)$, the vector $\mathbf g$ of the values $b_p\,\bar v_p|_{\eta=1/2}$ and the constraint matrix $\mathbf C_u = \left[\mathbf c_{ax}^{\mathsf T};\ \mathbf e^{\mathsf T}_{s_L}-\mathbf e^{\mathsf T}_{s_R}\right]$, the warping of each mode follows from the Karush–Kuhn–Tucker system

$$
\begin{bmatrix}\mathbf H^{\mathsf T}\mathbf W\mathbf H & \mathbf C_u^{\mathsf T}\\ \mathbf C_u & \mathbf 0\end{bmatrix}
\begin{bmatrix}\bar{\mathbf u}\\ \boldsymbol\mu\end{bmatrix} =
\begin{bmatrix}-\mathbf H^{\mathsf T}\mathbf W\mathbf g\\ \mathbf 0\end{bmatrix}, \eqn{kkt}
$$

which is the stationarity condition of Eq. (@vlasov) with the Lagrange multipliers $\boldsymbol\mu$.

The shear modes are the nodal hat functions of the warping. The hat functions of the two side nodes are merged into one mode with a common warping value,

$$
h_{s_R} \leftarrow h_{s_R}+h_{s_L},\qquad h_{s_L}\ \text{removed}, \eqn{shearmerge}
$$

and in the same way the nodal functions of the Mindlin modes of $s_L$ are added to those of $s_R$,

$$
\boldsymbol\Omega_{s_R}\leftarrow\boldsymbol\Omega_{s_R}+\boldsymbol\Omega_{s_L},\qquad \bar b_{s_R}\leftarrow\bar b_{s_R}+\bar b_{s_L},\qquad \boldsymbol\Omega_{s_L},\ \bar b_{s_L}\ \text{removed}. \eqn{mindlinmerge}
$$

The uniform axial mode is periodic by construction. Equations (@rows)–(@mindlinmerge) make every mode satisfy Eq. (@pbc). Any linear combination with arbitrary amplitude functions, and in particular any truncated subset of the modes, therefore satisfies the periodic condition exactly, and no penalty terms or Lagrange multipliers appear in the member analysis. The resulting numbers of modes are listed in Table 3.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 3.** Number of deformation modes per family of the periodic cell (the numbers do not depend on the lay-up).
:::

| $b_s$ [mm] | G | D | L | T | S (incl. axial) | M | total |
|:--|--:|--:|--:|--:|--:|--:|--:|
| 150 | 2 | 13 | 166 | 92 | 91 | 205 | 569 |
| 285.5 | 2 | 13 | 138 | 78 | 77 | 174 | 482 |

The periodic condition cannot be added afterwards to a set of classical modes. Classical GBT generates the conventional modes from a few natural nodes of an open contour. The periodic condition couples the two ends of the contour, so that the cell behaves like a closed loop with respect to its in-plane kinematics: the rigid rotation is no longer admissible, the warping must take the same value on both sides of the cut, and the classical hierarchy has to be recomputed. Imposed afterwards on a truncated set of classical modes, the three conditions of Eq. (@rows) remove three degrees of freedom from a small space and leave too few admissible modes. Here, the modes are generated from a complete nodal basis, and the constraint is applied before any truncation.

The periodic cell is also related to Bloch–Floquet theory. For an infinite panel with the stringer pitch $b_s$, the buckling modes are Bloch waves, $\mathbf u(X+b_s, x) = e^{\mathrm i\kappa b_s}\,\mathbf u(X,x)$ {{mead1996}}. A periodic cell of width $2b_s$ admits exactly the two real cases $\kappa b_s = 0$, in which the mode repeats in every bay, and $\kappa b_s = \pi$, in which it alternates in sign from bay to bay.

For visualisation, the conventional modes can be ordered as in classical GBT by the generalised eigenvalue problem $\hat{\mathbf B}\,\boldsymbol\psi = \Lambda\left(k_r^2\hat{\mathbf D}_1 + k_r^4\hat{\mathbf C}\right)\boldsymbol\psi$ at a reference wavenumber $k_r$, with the modal matrices of Section 2.3; a small value of $\Lambda$ indicates a preference for long waves. The first eight such modes of the periodic cell are shown in Figure 2. The two side nodes move identically in every mode, including the lateral rigid translation (mode 2), while the rigid rotation is absent.

![**Figure 2.** First eight hierarchical conventional modes of the periodic cell ($b_s = 150$ mm, Cross-ply lay-up, $k_r = 2\pi/b_s$; in-plane deformation, grey: undeformed). The side nodes $s_L$ and $s_R$ (squares) move identically in every mode.](fig/fig2_modes.png){width=15.5cm}

## Member analysis

The variation of the strain energy is expressed in terms of the amplitude functions and their variations,

$$
\delta U = \int_L\int_b\int_t\left(\sigma_{xx}\delta\varepsilon_{xx} + \sigma_{ss}\delta\varepsilon_{ss} + \tau_{xs}\delta\gamma_{xs} + \tau_{xz}\delta\gamma_{xz} + \tau_{sz}\delta\gamma_{sz}\right)dz\,ds\,dx
= \int_0^L \delta\boldsymbol\Phi^{\mathsf T}\,\boldsymbol\Xi\,\boldsymbol\Phi\,dx, \eqn{dU}
$$

where the vector $\boldsymbol\Phi$ collects the amplitude functions of all modes and their first and second derivatives,

$$
\boldsymbol\Phi = \begin{bmatrix}\boldsymbol\varphi\\ \boldsymbol\varphi^{I}\\ \boldsymbol\varphi^{II}\end{bmatrix},\qquad
\varphi^{I} = \frac{d\varphi}{dx},\qquad \varphi^{II} = \frac{d^2\varphi}{dx^2}, \eqn{Phi}
$$

with $\boldsymbol\varphi$ containing $\varphi_i$, $\varphi_j$, $\chi_r$ and $\vartheta_r$. With the generalised strains of Table 1 written as $\hat{\boldsymbol\varepsilon} = \sum_{d}\mathbf E^{(d)}(s)\,\boldsymbol\varphi^{(d)}(x)$, $d \in \{0, I, II\}$, the GBT matrix $\boldsymbol\Xi$ consists of the blocks

$$
\boldsymbol\Xi = \begin{bmatrix}\boldsymbol\Xi^{(0,0)} & \boldsymbol\Xi^{(0,I)} & \boldsymbol\Xi^{(0,II)}\\ \boldsymbol\Xi^{(I,0)} & \boldsymbol\Xi^{(I,I)} & \boldsymbol\Xi^{(I,II)}\\ \boldsymbol\Xi^{(II,0)} & \boldsymbol\Xi^{(II,I)} & \boldsymbol\Xi^{(II,II)}\end{bmatrix},\qquad
\Xi^{(d,e)}_{kl} = \int_\Gamma \mathbf E^{(d)\mathsf T}_k\,\mathbf C\,\mathbf E^{(e)}_l\,ds, \eqn{Xi}
$$

where $\Gamma$ is the mid-line of the cross-section and $\boldsymbol\Xi^{(e,d)} = \boldsymbol\Xi^{(d,e)\mathsf T}$. As in the GBT–Ritz approach {{kharghani2026jcs}}, $\boldsymbol\Xi$ is decomposed into the components of the conventional modes, $\hat{\boldsymbol\Xi}$, of the shear modes, $\check{\boldsymbol\Xi}$, and of the interaction between them, $\tilde{\boldsymbol\Xi}$, supplemented here by the components $\boldsymbol\Xi^{M}$ that involve the Mindlin modes. For the conventional modes,

$$
\hat{\boldsymbol\Xi} = \begin{bmatrix}\hat{\mathbf B} & \hat{\mathbf F}_1 & \hat{\mathbf D}_2\\ \hat{\mathbf F}_1^{\mathsf T} & \hat{\mathbf D}_1 & \hat{\mathbf H}_1\\ \hat{\mathbf D}_2^{\mathsf T} & \hat{\mathbf H}_1^{\mathsf T} & \hat{\mathbf C}\end{bmatrix}, \eqn{Xihat}
$$

and substituting Table 1 and Eq. (@claw) into Eq. (@Xi) gives, for symmetric and unsymmetric laminates, with $\gamma_i = \bar u_{i,s}+\bar v_i$ and the conventional modes $i$ and $k$,

$$
\hat B_{ik} = \int_\Gamma\Big[A_{22}\,\bar v_{i,s}\bar v_{k,s} - B_{22}\left(\bar v_{i,s}\bar w_{k,ss}+\bar w_{i,ss}\bar v_{k,s}\right) + D_{22}\,\bar w_{i,ss}\bar w_{k,ss}\Big]ds, \eqn{Bmat}
$$

$$
\hat F_{1,ik} = \int_\Gamma\Big[A_{26}\,\bar v_{i,s}\gamma_k - B_{26}\left(2\bar v_{i,s}\bar w_{k,s}+\bar w_{i,ss}\gamma_k\right) + 2D_{26}\,\bar w_{i,ss}\bar w_{k,s}\Big]ds, \eqn{F1mat}
$$

$$
\hat D_{2,ik} = \int_\Gamma\Big[A_{12}\,\bar v_{i,s}\bar u_k - B_{12}\left(\bar v_{i,s}\bar w_k+\bar w_{i,ss}\bar u_k\right) + D_{12}\,\bar w_{i,ss}\bar w_k\Big]ds, \eqn{D2mat}
$$

$$
\hat D_{1,ik} = \int_\Gamma\Big[A_{66}\,\gamma_i\gamma_k - 2B_{66}\left(\gamma_i\bar w_{k,s}+\bar w_{i,s}\gamma_k\right) + 4D_{66}\,\bar w_{i,s}\bar w_{k,s}\Big]ds, \eqn{D1mat}
$$

$$
\hat H_{1,ik} = \int_\Gamma\Big[A_{16}\,\gamma_i\bar u_k - B_{16}\left(\gamma_i\bar w_k+2\bar w_{i,s}\bar u_k\right) + 2D_{16}\,\bar w_{i,s}\bar w_k\Big]ds, \eqn{H1mat}
$$

$$
\hat C_{ik} = \int_\Gamma\Big[A_{11}\,\bar u_i\bar u_k - B_{11}\left(\bar u_i\bar w_k+\bar w_i\bar u_k\right) + D_{11}\,\bar w_i\bar w_k\Big]ds. \eqn{Cmat}
$$

The transverse membrane stiffness enters through the $A_{22}$ term of $\hat{\mathbf B}$, the bending–twisting coupling through $\hat{\mathbf F}_1$ and $\hat{\mathbf H}_1$, and the membrane–bending coupling through all $B_{ij}$ terms. For the shear modes $j$ and $l$,

$$
\begin{gathered}
\check{\boldsymbol\Xi} = \begin{bmatrix}\check{\mathbf D}_1 & \check{\mathbf H}_1 & \mathbf 0\\ \check{\mathbf H}_1^{\mathsf T} & \check{\mathbf C} & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \mathbf 0\end{bmatrix},\\
\check D_{1,jl} = \int_\Gamma A_{66}\,\bar u_{j,s}\bar u_{l,s}\,ds,\quad
\check H_{1,jl} = \int_\Gamma A_{16}\,\bar u_{j,s}\bar u_{l}\,ds,\quad
\check C_{jl} = \int_\Gamma A_{11}\,\bar u_j\bar u_l\,ds,
\end{gathered} \eqn{Xicheck}
$$

and for the interaction between a conventional mode $i$ (rows) and a shear mode $j$ (columns),

$$
\tilde{\boldsymbol\Xi} = \begin{bmatrix}\tilde{\mathbf F} & \tilde{\mathbf D}_2 & \mathbf 0\\ \tilde{\mathbf D}_1 & \tilde{\mathbf H}_1 & \mathbf 0\\ \tilde{\mathbf H}_3 & \tilde{\mathbf C} & \mathbf 0\end{bmatrix}, \eqn{Xitilde}
$$

$$
\begin{gathered}
\tilde F_{ij} = \int_\Gamma\left(A_{26}\,\bar v_{i,s} - B_{26}\,\bar w_{i,ss}\right)\bar u_{j,s}\,ds,\qquad
\tilde D_{2,ij} = \int_\Gamma\left(A_{12}\,\bar v_{i,s} - B_{12}\,\bar w_{i,ss}\right)\bar u_{j}\,ds,\\
\tilde D_{1,ij} = \int_\Gamma\left(A_{66}\,\gamma_i - 2B_{66}\,\bar w_{i,s}\right)\bar u_{j,s}\,ds,\qquad
\tilde H_{1,ij} = \int_\Gamma\left(A_{16}\,\gamma_i - 2B_{16}\,\bar w_{i,s}\right)\bar u_{j}\,ds,\\
\tilde H_{3,ij} = \int_\Gamma\left(A_{16}\,\bar u_i - B_{16}\,\bar w_i\right)\bar u_{j,s}\,ds,\qquad
\tilde C_{ij} = \int_\Gamma\left(A_{11}\,\bar u_i - B_{11}\,\bar w_i\right)\bar u_{j}\,ds .
\end{gathered} \eqn{Xitildemat}
$$

The components $\boldsymbol\Xi^{M}$ of the Mindlin modes follow from Eq. (@Xi) and Table 1 in the same way; they contain the transverse shear stiffness, e.g. $\int_\Gamma A_{s,11}\,\bar c_r\bar c_{r'}\,ds$ in the block $(0,0)$ of two $\beta_x$ modes, and the bending stiffnesses $\mathbf D$ and $\mathbf B$ through the curvatures of Table 1. All blocks are evaluated by four-point Gauss quadrature on every segment.

The panel is loaded at both loaded edges by a compressive line load $q(s)$ per unit length, which acts at the distance $e_q(s)$ from the wall reference surface. Its external work is

$$
W_{ext} = \int_\Gamma q(s)\left[u_q(0,s)-u_q(L,s)\right]ds,\qquad u_q = u - e_q\,w_{,x}, \eqn{wext}
$$

where $u_q$ is the axial displacement of the fibre at which the load acts. The pre-buckling state is not assumed but computed: it is the first-order solution of the stiffness relation under this load and under the end and side conditions of the cell, $\mathbf K\mathbf a^0 = \mathbf f$, with $\mathbf K$ and $\mathbf a$ as defined below and the load vector $\mathbf f$ from the variation of Eq. (@wext). The membrane forces of this state are

$$
\left[S_{xx},\ S_{ss},\ S_{xs}\right]^{\mathsf T} = \mathbf A\,\boldsymbol\varepsilon^{0} + \mathbf B\,\boldsymbol\kappa^{0}, \eqn{pre}
$$

with the membrane strains $\boldsymbol\varepsilon^0$ and curvatures $\boldsymbol\kappa^0$ of the pre-buckling state. They vary along $s$ and $x$: the load is introduced at the ends and shared between skin and stringer by shear lag, the eccentric foot load bends the cell near the ends, and the restrained Poisson expansion produces transverse membrane forces $S_{ss}$.

The loading is characterised by the second-order energy of the membrane forces,

$$
V = \frac12\int_0^L\!\!\int_\Gamma \sum_{c\in\{u,v,w\}}\Big[S_{xx}\,c_{,x}^2 + S_{ss}\,c_{,s}^2 + 2S_{xs}\,c_{,x}\,c_{,s}\Big]ds\,dx
= \frac12\int_0^L \boldsymbol\Phi^{\mathsf T}\boldsymbol{\mathcal X}(x)\,\boldsymbol\Phi\,dx, \eqn{V}
$$

in which all three displacement components and all three membrane forces are kept. For a conventional mode, the displacement gradients are $u_{,x} = \bar u_i\varphi_i^{II}$, $u_{,s} = \bar u_{i,s}\varphi_i^{I}$, $v_{,x} = \bar v_i\varphi_i^{I}$, $v_{,s} = \bar v_{i,s}\varphi_i$, $w_{,x} = \bar w_i\varphi_i^{I}$ and $w_{,s} = \bar w_{i,s}\varphi_i$, and for a shear mode $u_{,x} = \bar u_j\varphi_j^{I}$ and $u_{,s} = \bar u_{j,s}\varphi_j$; the Mindlin modes do not move the mid-surface and do not contribute to $\boldsymbol{\mathcal X}$. Since the membrane forces vary along $x$, the matrix $\boldsymbol{\mathcal X}(x)$ is evaluated at every integration point along the cell. For the conventional modes and a uniform membrane force $S_{xx}$ with $S_{ss} = S_{xs} = 0$, Eq. (@V) reduces to the geometric matrix of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}},

$$
\boldsymbol{\mathcal X} = \begin{bmatrix}\mathbf 0 & \mathbf 0 & \mathbf 0\\ \mathbf 0 & \boldsymbol{\mathcal X}^{VW} & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \boldsymbol{\mathcal X}^{U}\end{bmatrix},\qquad
\mathcal X^{U}_{ik} = \int_\Gamma S_{xx}\,\bar u_i\bar u_k\,ds,\qquad
\mathcal X^{VW}_{ik} = \int_\Gamma S_{xx}\left(\bar v_i\bar v_k+\bar w_i\bar w_k\right)ds. \eqn{Xclassical}
$$

Following the Ritz method, $\boldsymbol\Phi$ is written as the product of the matrix $\boldsymbol\Psi$ of shape functions and the vector $\mathbf a$ of unknown coefficients,

$$
\boldsymbol\Phi = \boldsymbol\Psi\,\mathbf a,\qquad
\varphi(x) = \sum_{b=1}^{n_b} N_b(x)\,a_b . \eqn{Psi}
$$

In the GBT–Ritz approach for members, trigonometric series are used {{kharghani2026jcs}}. In the panel cell, the pre-buckling state varies strongly near the loaded edges, where the load enters the stringer, and a single trigonometric series represents this variation only with many terms. The shape functions $N_b$ are therefore cubic B-splines on an open, uniform knot vector with spans of length $h_x = \min(15\ \text{mm},\ \max(2\ \text{mm},\ W/40),\ L/10)$, i.e. $h_x = b_s/20$ for all panels studied. They are local, represent simply supported and clamped edges alike, and with the open knot vector $\varphi(0) = a_1$, $\varphi(L) = a_{n_b}$ and $\varphi^{I}(0) = 0 \Leftrightarrow a_2 = a_1$, so that the essential end conditions act directly on the end coefficients. Two end conditions of the loaded edges are considered, as in the shell models:

* **Simply supported (SS):** the out-of-plane displacement of every end edge vanishes, $U_Y = 0$ at every node at $x = 0$ and $x = L$, i.e. $\sum_i U_{Y,qi}\,\varphi_i(0) = \sum_i U_{Y,qi}\,\varphi_i(L) = 0$ for every node $q$. These linear relations between the end coefficients are imposed by a null-space reduction; the in-plane displacement $U_X$, the warping and the rotations of the end cross-sections remain free. Rigid-body motions are removed by fixing the axial mode and $U_X$ at the side node $s_R$ at $x = 0$.
* **Clamped (CC):** $\varphi_i = \varphi_i^{I} = 0$ at both ends for every conventional mode (in-plane displacements and the rotations $w_{,x}$ and $\theta$ vanish), $\chi_r = \vartheta_r = 0$ at both ends, and $\varphi_j(0) = 0$ for the shear modes, so that $u = 0$ on the whole edge $x = 0$; the edge $x = L$ remains free to warp under the uniform load.

In classical GBT, simply supported ends are represented by $\varphi_i = \sin(h\pi x/L)$, $h = 1, 2, \dots$, which suppresses all in-plane displacements of the end cross-sections while leaving their warping free. The present SS condition is weaker: like the shell model, it only suppresses the out-of-plane displacement $U_Y$ of the end edges.

With Eq. (@Psi), the stationarity of the total potential energy, $\delta(U + \lambda V) = 0$, leads to the eigenvalue problem

$$
\left(\mathbf K + \lambda\,\mathbf G\right)\mathbf a = \mathbf 0, \eqn{eig}
$$

where $\mathbf K$ is the stiffness matrix associated with $\boldsymbol\Xi$ and $\mathbf G$ the geometric stiffness matrix arising from the membrane forces $\boldsymbol{\mathcal X}$ of the reference load,

$$
\mathbf K = \int_0^L \boldsymbol\Psi^{\mathsf T}\left(\hat{\boldsymbol\Xi} + \check{\boldsymbol\Xi} + \tilde{\boldsymbol\Xi} + \boldsymbol\Xi^{M}\right)\boldsymbol\Psi\,dx,\qquad
\mathbf G = \int_0^L \boldsymbol\Psi^{\mathsf T}\boldsymbol{\mathcal X}\,\boldsymbol\Psi\,dx, \eqn{KG}
$$

in which every component of $\boldsymbol\Xi$ is placed at the rows and columns of its modes, and $\tilde{\boldsymbol\Xi}$ also at their transposed position. The integrals are evaluated with four Gauss points per spline span. The smallest positive eigenvalue $\lambda$ multiplies the reference load of 1 N/mm, so that the critical line load is $N_{cr} = \lambda\times 1$ N/mm, and the eigenvector $\mathbf a$ defines the buckling mode. Since $\mathbf G$ is indefinite (the pre-buckling state contains membrane shear and some transverse tension), Eq. (@eig) is solved for the largest values of $\mu = 1/\lambda$ of $-\mathbf G\,\mathbf a = \mu\,\mathbf K\,\mathbf a$ by the implicitly restarted Lanczos method {{lehoucq1998}} with a sparse factorisation of $\mathbf K$. For $b_s = 150$ mm and $L = 225$ mm the model has $569\times 33 = 18\,777$ coefficients before the end conditions are applied.

Two numerical aspects are noted. First, at some intermediate nodes of the inclined webs the two segment tangents are identical up to round-off, and the orthonormal basis of the tangents contains a spurious second Mindlin vector normal to the web. For it $\bar c_r = \mathbf t\cdot\boldsymbol\Omega \approx 10^{-15}$ on every segment, so that the mode carries neither strain energy nor displacement gradients (11 such modes for $b_s = 150$ mm and 8 for $b_s = 285.5$ mm). Such modes cannot change any eigenvalue, but they make $\mathbf K$ singular; they are identified by their vanishing diagonal of $\boldsymbol\Xi$ and their coefficients are held. Second, at a cross-section integration point only a few modes are non-zero (about 25 of 569 for $b_s = 150$ mm), so that the cross-section factors of $\boldsymbol{\mathcal X}$ are stored as sparse matrices that are computed once per cross-section; for every panel length only the membrane forces change.

The formulation contains classical GBT as a special case. For the conventional modes alone, walls without coupling ($\mathbf B = \mathbf 0$, $A_{16} = A_{26} = D_{16} = D_{26} = 0$), and a uniform membrane force $S_{xx}$, the matrices $\hat{\mathbf F}_1$ and $\hat{\mathbf H}_1$ vanish, and the Euler–Lagrange equations of $\delta(U + \lambda V) = 0$ with Eqs. (@Xihat) and (@Xclassical) are

$$
\hat{\mathbf C}\,\boldsymbol\varphi^{IIII} - \hat{\mathbf D}\,\boldsymbol\varphi^{II} + \hat{\mathbf B}\,\boldsymbol\varphi + \lambda\left(\boldsymbol{\mathcal X}^{U}\boldsymbol\varphi^{IIII} - \boldsymbol{\mathcal X}^{VW}\boldsymbol\varphi^{II}\right) = \mathbf 0,\qquad
\hat{\mathbf D} = \hat{\mathbf D}_1 - \hat{\mathbf D}_2 - \hat{\mathbf D}_2^{\mathsf T}, \eqn{classical}
$$

where the integration by parts of the $\hat{\mathbf D}_2$ terms gives the classical matrix $\hat{\mathbf D}$. Without the warping term $\boldsymbol{\mathcal X}^{U}$, Eq. (@classical) is the classical GBT buckling equation of Schardt {{schardt1994}} and of Silvestre and Camotim for orthotropic members {{silvestre2002b}}. For simply supported ends and $\boldsymbol\varphi = \hat{\boldsymbol\varphi}\sin k_h x$, $k_h = h\pi/L$, it becomes the algebraic eigenvalue problem

$$
\left[\hat{\mathbf B} + k_h^2\,\hat{\mathbf D} + k_h^4\,\hat{\mathbf C} + \lambda\left(k_h^2\,\boldsymbol{\mathcal X}^{VW} + k_h^4\,\boldsymbol{\mathcal X}^{U}\right)\right]\hat{\boldsymbol\varphi} = \mathbf 0, \eqn{signature}
$$

whose smallest eigenvalue over $h$ is the classical GBT signature curve. The periodic condition does not alter this structure: it only restricts the admissible modes of the cross-section. The present model differs from the classical one only by the terms that the panel requires: laminate coupling ($\hat{\mathbf F}_1$, $\hat{\mathbf H}_1$, $B_{ij}$), shear and Mindlin modes, a non-uniform pre-buckling state with $S_{ss}$ and $S_{xs}$, and end conditions that support the loaded edges without restraining their warping.

## Mode decomposition

The eigenvector $\mathbf a$ defines the buckling mode in terms of the GBT modes but does not directly show the contribution of the individual deformation mechanisms. As in the GBT–Ritz approach {{kharghani2026jcs}}, an energy-based mode decomposition is therefore performed. The modes are grouped into the families $g \in \{G, D, L, T, S, M\}$, and $\mathbf a_g$ denotes the vector $\mathbf a$ with all coefficients set to zero except those of family $g$. The strain energy of the buckling mode, $U = \tfrac12\int_0^L\boldsymbol\Phi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Phi\,dx = \tfrac12\mathbf a^{\mathsf T}\mathbf K\,\mathbf a$, is divided into the contributions of the families,

$$
U_g = \frac12\,\mathbf a_g^{\mathsf T}\mathbf K\,\mathbf a,\qquad
\rho_g = \frac{U_g}{U},\qquad \sum_g \rho_g = 1 . \eqn{rho}
$$

If the families are orthogonal with respect to $\mathbf K$, $U_g$ equals the energy of family $g$ alone, $\tfrac12\mathbf a_g^{\mathsf T}\mathbf K\,\mathbf a_g$, which is the form used for members {{kharghani2026jcs}}. In the panel cell, the distortional and local modes are not orthogonal: the D modes move the wall junctions without rotating the nodes, and the rotation modes of family L restore the continuity of the walls. The energies of the families alone then overestimate the families that are coupled to others and do not add up to $U$, whereas Eq. (@rho) divides the strain energy exactly. For instance, for the Unsym panel with simply supported loaded edges, $b_s = 150$ mm and $\beta = 1.5$, the energies of the families alone would assign 41 % to the distortional and 58 % to the local family, whereas Eq. (@rho) gives 1 % and 97 %.
