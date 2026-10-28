
# Theory and formulations

The formulation follows the two stages of GBT {{schardt1994,silvestre2002a,bebiano2015}}: a cross-section analysis, which determines the deformation modes and the modal mechanical properties (tensors), and a member analysis, which determines the modal amplitude functions. The notation is that of the shear-deformable GBT of Silvestre and Camotim {{silvestre2013shear}}, and the member analysis follows the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}. Three extensions are required for the panel cell: (i) the deformation modes must satisfy the periodic boundary conditions of the cell (Section 2.2), (ii) Vlasov's assumptions of null membrane shear strain and transverse extension cannot be imposed on the conventional modes of a cell that contains closed contours, and (iii) the walls are shear-deformable (Reissner–Mindlin) laminates. Each of them is introduced such that the classical GBT equations are recovered as a special case (Section 2.3.6).

## Cross-section analysis

### Periodic cell, discretisation and displacement field

The panel is the periodic cell of an infinitely wide panel with identical omega stringers at the pitch $b_s$ (Figure 1a). The cell contains two stringers and has the width $W = 2b_s$ and the length $L = \beta b_s$, with the aspect ratio $\beta$. The cross-section axes are $X$ across the panel and $Y$ normal to the skin, towards the stringers; $x$ is the member axis along the stringers. The skin occupies $-W \le X \le 0$, and its side edges $X = -W$ and $X = 0$ lie in the middle of a skin bay. As in Silvestre and Camotim {{silvestre2013shear}}, every wall has the local axes $(x, s, z)$, with $s$ along the wall mid-line and $z$ through the thickness, and the displacement components $u$, $v$ and $w$ (Figure 1b). The mid-line is divided into walls (skin, feet, webs, crowns) and the walls into segments of width at most $b_s/25$, with at least two segments per wall; the end points of the walls are the natural nodes and the remaining segment end points the intermediate nodes {{schardt1994,bebiano2015}}. For $b_s = 150$ mm the cell has $N = 92$ nodes, 19 walls and 93 segments. The two side nodes of the skin are denoted $s_L$ ($X = -W$) and $s_R$ ($X = 0$).

![**Figure 1.** (a) Infinitely wide omega-stringer-stiffened panel and its periodic cell of width $W = 2b_s$. (b) GBT discretisation of the cell for $b_s = 150$ mm with natural and intermediate nodes; skin, feet (skin and flange as one laminate), webs, crowns and link walls of the offset tie; $s_L$ and $s_R$ are the side nodes on which the periodic condition acts.](fig/fig1_geometry.png){width=15.5cm}

The displacement field is the modal representation of shear-deformable GBT {{silvestre2013shear}},

$$
u(x,s) = u_k(s)\,\phi_{k,x}(x) + u_j(s)\,\varphi_j(x),\qquad
v(x,s) = v_k(s)\,\phi_k(x),\qquad
w(x,s) = w_k(s)\,\phi_k(x), \eqn{disp}
$$

extended here by the transverse shear rotations of the Reissner–Mindlin walls,

$$
\beta_x(x,s) = \xi_r(s)\,\chi_r(x),\qquad \beta_s(x,s) = \zeta_r(s)\,\chi_r(x), \eqn{beta}
$$

where the summation convention applies to the mode indices; $k$ (and $i$ for virtual quantities) denotes the conventional modes, $j$ (and $h$) the shear modes, and $r$ the transverse shear (Mindlin) modes. $u_k(s)$, $v_k(s)$ and $w_k(s)$ are the displacement profiles of conventional mode $k$ and $\phi_k(x)$ its amplitude function; $u_j(s)$ is the warping profile of shear mode $j$ and $\varphi_j(x)$ its amplitude function; $\xi_r(s)$ and $\zeta_r(s)$ are the rotation profiles of Mindlin mode $r$ and $\chi_r(x)$ its amplitude function. The uniform axial extension is represented by a shear-type mode with $u_j = 1$, i.e. by $\varphi_j = \phi_{1,x}$ of the axial extension mode of classical GBT.

### Strain–displacement relations

The rotations of the wall normal are $\phi_x = -w_{,x}+\beta_x$ and $\phi_s = -w_{,s}+\beta_s$, and the linear strains at the distance $z$ from the mid-surface follow as

$$
\varepsilon_{xx} = u_{,x} + z\,\kappa_{xx},\qquad \varepsilon_{ss} = v_{,s} + z\,\kappa_{ss},\qquad \gamma_{xs} = u_{,s} + v_{,x} + z\,\kappa_{xs},\qquad \gamma_{xz} = \beta_x,\qquad \gamma_{sz} = \beta_s, \eqn{strain0}
$$

with the curvatures $\kappa_{xx} = -w_{,xx} + \beta_{x,x}$, $\kappa_{ss} = -w_{,ss} + \beta_{s,s}$ and $\kappa_{xs} = -2w_{,xs} + \beta_{x,s} + \beta_{s,x}$. Inserting Eqs. (@disp) and (@beta) gives

$$
\begin{gathered}
\varepsilon_{xx} = \left(u_k - z\,w_k\right)\phi_{k,xx} + u_j\,\varphi_{j,x} + z\,\xi_r\,\chi_{r,x},\qquad
\varepsilon_{ss} = \left(v_{k,s} - z\,w_{k,ss}\right)\phi_k + z\,\zeta_{r,s}\,\chi_r,\\
\gamma_{xs} = \left(\gamma_k - 2z\,w_{k,s}\right)\phi_{k,x} + u_{j,s}\,\varphi_j + z\left(\xi_{r,s}\,\chi_r + \zeta_r\,\chi_{r,x}\right),\qquad
\gamma_{xz} = \xi_r\,\chi_r,\qquad \gamma_{sz} = \zeta_r\,\chi_r,
\end{gathered} \eqn{strain}
$$

where $\gamma_k = u_{k,s} + v_k$ is the membrane shear strain profile of conventional mode $k$. In the GBT of Silvestre and Camotim {{silvestre2013shear}}, the conventional modes satisfy Vlasov's assumptions, $\gamma_k = 0$ and $v_{k,s} = 0$, and Eq. (@strain) reduces to their kinematic relations ($\xi_r = \zeta_r = 0$). Here these assumptions are not imposed, for two reasons: the transverse extension modes ($v_{k,s} \ne 0$) are part of the basis, and the omega stringer and the skin form closed cells and, with the periodic condition, the skin itself behaves as a closed contour (Section 2.2.3), so that $\gamma_k = 0$ cannot be satisfied by all distortional modes. Table 1 lists the generalised strains of a mode for a unit amplitude function and its derivatives.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 1.** Non-zero generalised (mid-surface) strains and curvatures of one mode, multiplied by the amplitude function or its derivatives; terms in brackets vanish under Vlasov's assumptions.
:::

| Mode | multiplied by $\phi_k$ / $\varphi_j$ / $\chi_r$ | multiplied by $\phi_{k,x}$ / $\varphi_{j,x}$ / $\chi_{r,x}$ | multiplied by $\phi_{k,xx}$ |
|:--|:--|:--|:--|
| conventional $k$ | $[\varepsilon^{M}_{ss} = v_{k,s}]$, $\kappa_{ss} = -w_{k,ss}$ | $[\gamma^{M}_{xs} = \gamma_k]$, $\kappa_{xs} = -2w_{k,s}$ | $\varepsilon^{M}_{xx} = u_k$, $\kappa_{xx} = -w_k$ |
| shear $j$ | $\gamma^{M}_{xs} = u_{j,s}$ | $\varepsilon^{M}_{xx} = u_j$ | – |
| Mindlin $r$ | $\kappa_{xs} = \xi_{r,s}$, $\kappa_{ss} = \zeta_{r,s}$, $\gamma_{xz} = \xi_r$, $\gamma_{sz} = \zeta_r$ | $\kappa_{xx} = \xi_r$, $\kappa_{xs} = \zeta_r$ | – |

The geometrically non-linear membrane strains needed for the geometric stiffness are

$$
\varepsilon^{NL}_{xx} = \tfrac12\left(u_{,x}^2+v_{,x}^2+w_{,x}^2\right),\qquad
\varepsilon^{NL}_{ss} = \tfrac12\left(u_{,s}^2+v_{,s}^2+w_{,s}^2\right),\qquad
\gamma^{NL}_{xs} = u_{,x}u_{,s}+v_{,x}v_{,s}+w_{,x}w_{,s}, \eqn{strainNL}
$$

of which classical GBT retains $\tfrac12\left(v_{,x}^2+w_{,x}^2\right) = \tfrac12\left(v_iv_k+w_iw_k\right)\phi_{i,x}\phi_{k,x}$ {{silvestre2002b,silvestre2013shear}}.

### Constitutive relations

In every layer $p$ of a wall, the in-plane stresses follow from the reduced stiffness components $\overline{Q}^{(p)}_{ij}$ ($i,j = 1,2,6$) referred to the local axes $(x,s)$, and the transverse shear stresses from $\overline{\mathbf Q}^{(p)}_{s}$. Integration through the thickness gives the membrane, coupling and bending stiffnesses and the transverse shear stiffness {{mindlin1951,whitney1970,reddy2004}},

$$
\left(A_{ij}, B_{ij}, D_{ij}\right) = \sum_{p=1}^{N_L}\int_{z_{p}}^{z_{p+1}}\overline{Q}^{(p)}_{ij}\left(1, z, z^2\right)dz,\qquad
\mathbf A_s = \frac56\sum_{p=1}^{N_L}\overline{\mathbf Q}^{(p)}_s\left(z_{p+1}-z_p\right), \eqn{abd}
$$

and the stress resultants of a wall are

$$
\begin{bmatrix}\mathbf N\\ \mathbf M\\ \mathbf Q\end{bmatrix} =
\begin{bmatrix} \mathbf A & \mathbf B & \mathbf 0\\ \mathbf B & \mathbf D & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \mathbf A_s\end{bmatrix}
\begin{bmatrix}\boldsymbol\varepsilon^{M}\\ \boldsymbol\kappa\\ \boldsymbol\gamma_z\end{bmatrix}, \eqn{claw}
$$

with $\mathbf N = [N_{xx}, N_{ss}, N_{xs}]^{\mathsf T}$, $\mathbf M = [M_{xx}, M_{ss}, M_{xs}]^{\mathsf T}$, $\mathbf Q = [Q_x, Q_s]^{\mathsf T}$, $\boldsymbol\varepsilon^{M} = [\varepsilon^M_{xx}, \varepsilon^M_{ss}, \gamma^M_{xs}]^{\mathsf T}$, $\boldsymbol\kappa = [\kappa_{xx}, \kappa_{ss}, \kappa_{xs}]^{\mathsf T}$ and $\boldsymbol\gamma_z = [\gamma_{xz}, \gamma_{sz}]^{\mathsf T}$. All entries of $\mathbf A$, $\mathbf B$ and $\mathbf D$ are retained: extension–bending, extension–shear and bending–twisting coupling, as in {{silvestre2013shear}}, and the membrane–bending coupling of unsymmetric laminates {{kharghani2026jcs}}. Since the transverse extension is a degree of freedom, $N_{ss}$ is not condensed out.

A stringer foot consists of the skin and the flange, which are tied in the shell model; in GBT they form one wall referred to the skin mid-plane. With $\mathbf A_0$, $\mathbf B_0$, $\mathbf D_0$ of the combined stack about its own mid-plane and $z_f$ the distance between the two mid-planes,

$$
\mathbf A_f = \mathbf A_0,\qquad \mathbf B_f = \mathbf B_0 + z_f\,\mathbf A_0,\qquad \mathbf D_f = \mathbf D_0 + 2z_f\,\mathbf B_0 + z_f^2\,\mathbf A_0 . \eqn{foot}
$$

In the shell model the stringer mid-surface is tied to the skin with the offset $c = (t_s+t_f)/2$, where $t_s$ and $t_f$ are the thicknesses of skin and stringer. In GBT every web root is a node at $Y = c$, joined to the skin node below it by a link wall of height $c$, which transmits the in-plane displacements and the rotation like a line of rigid offset links without adding axial or twisting stiffness:

$$
\begin{gathered}
A^{L}_{11} = 10^{-3}A_{0,11},\quad A^{L}_{22} = k_L A_{0,22},\quad A^{L}_{66} = k_L A_{0,66},\\
D^{L}_{11} = 10^{-3}D_{0,11},\quad D^{L}_{22} = k_L D_{0,22},\quad D^{L}_{66} = 10^{-3}D_{0,66},\quad \mathbf A^{L}_s = k_L\mathbf A_{s0},
\end{gathered} \eqn{link}
$$

with the skin stiffnesses $\mathbf A_0$, $\mathbf D_0$, $\mathbf A_{s0}$, $k_L = 10^3$ and all other entries zero. The results do not depend on $k_L$, whereas a stiff twisting term $D^{L}_{66}$ would lock the rotation of the web root (Section 3.6).

### Nodal description and deformation mode families

Within every segment $p$ of width $b_p$ between the nodes $q$ and $q'$, with $\eta = (s-S_q)/b_p \in [0,1]$, the profiles $u$ and $v$ vary linearly and are described by their nodal values. The transverse displacement is the cubic Hermite interpolation of the nodal normal displacements $w_q$, $w_{q'}$ and the nodal rotations $\theta_q$, $\theta_{q'}$ ($\theta = w_{,s}$),

$$
\begin{gathered}
w(\eta) = H_1(\eta)\,w_q + b_p H_2(\eta)\,\theta_q + H_3(\eta)\,w_{q'} + b_p H_4(\eta)\,\theta_{q'},\\
H_1 = 1-3\eta^2+2\eta^3,\quad H_2 = \eta-2\eta^2+\eta^3,\quad H_3 = 3\eta^2-2\eta^3,\quad H_4 = -\eta^2+\eta^3 .
\end{gathered} \eqn{hermite}
$$

Equation (@hermite) spans the same space as Schardt's description of the transverse bending by nodal moments {{schardt1994}}, which is used in the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}; the nodal form is used here because it turns the periodic conditions into linear relations between nodal quantities. The rotation is shared by all walls that meet at a node (rigid fold lines). With the in-plane nodal displacements $\mathbf U = [\mathbf U_1;\dots;\mathbf U_N]\in\mathbb R^{2N}$, $\mathbf U_q = (U_{X,q}, U_{Y,q})$, a segment with unit tangent $\mathbf t$ and unit normal $\mathbf n$ has $v = \mathbf t\cdot\mathbf U$ and $w = \mathbf n\cdot\mathbf U$ at its end nodes, and $\boldsymbol\Theta \in\mathbb R^N$ collects the nodal rotations. Let $\mathbf E_\ell$ be the matrix of segment elongations, $(\mathbf E_\ell\mathbf U)_p = \mathbf t_p\cdot(\mathbf U_{q'}-\mathbf U_q)$. Following the organisation of the cross-section analysis by Bebiano et al. {{bebiano2015}} and the shear modes of Gonçalves et al. {{goncalves2014shear}}, the conventional modes are generated in the families

* **G (global):** rigid in-plane translations and rotation of the cross-section;
* **D (distortional):** in-plane nodal displacements without segment elongation, $\mathbf E_\ell\mathbf U = \mathbf 0$, orthogonal to G and L;
* **L (local):** a unit normal displacement of one node that is not a wall junction, and a pure rotation of one node;
* **T (transverse extension):** the segment elongation patterns, i.e. the columns of $\mathbf E_\ell^{\mathsf T}$;

and the non-conventional modes in the families

* **S (shear):** the nodal hat functions $h_q(s)$ of the warping, one of which is replaced by the axial extension mode $u_j = 1$;
* **M (Mindlin):** for $\beta_x$, a nodal vector $\boldsymbol\Omega_q$ in the plane of the cross-section with $\xi_r = \mathbf t\cdot\boldsymbol\Omega_q$, interpolated linearly on every segment (one vector at a node inside a straight wall, two at a fold line); for $\beta_s$, a nodal value of $\zeta_r$, interpolated linearly.

Since $\mathbb R^{2N} = \ker\mathbf E_\ell \oplus \operatorname{range}\mathbf E_\ell^{\mathsf T}$, the families G, D, L and T span all in-plane nodal displacements, and the rotation modes make the nodal rotations independent; the conventional modes are therefore a complete description of the in-plane deformation of the discretised cross-section, organised in the GBT families. In classical GBT the warping profile $u_k$ of the G and D modes follows from Vlasov's condition $\gamma_k = u_{k,s} + v_k = 0$ {{vlasov1961,schardt1994}}. Since this condition cannot be met in closed contours, the warping is taken as the stiffness-weighted least-squares solution

$$
\min_{\mathbf u}\ \Pi_w(\mathbf u) = \sum_p \hat A_{66,p}\,b_p\left(\frac{u_{q'}-u_{q}}{b_p}+v_p\!\left(\eta = \tfrac12\right)\right)^2
\quad\text{subject to}\quad \mathbf c_{ax}^{\mathsf T}\mathbf u = 0, \eqn{vlasov}
$$

where $\mathbf u$ collects the nodal warping values, $\hat A_{66} = A_{66}-A_{26}^2/A_{22}$, and $c_{ax,q} = \sum_{p\ni q}\tfrac12\hat A_{11,p}\,b_p$ with $\hat A_{11} = A_{11}-A_{12}^2/A_{22}$, so that the warping carries no axial resultant. $\Pi_w$ is the membrane shear energy per unit length of a unit $\phi_{k,x}$; for an open contour its minimum is zero and Eq. (@vlasov) returns the classical warping, while in closed contours the remaining membrane shear strain $\gamma_k$ is retained in the strains of Eq. (@strain) and complemented by the shear modes.

## Periodic boundary conditions in the cross-section analysis

### Continuous condition and its modal form

The cell is one period of an infinitely wide panel. Every field must therefore take the same value at $X = -W$ and $X = 0$, for every $x$:

$$
\mathcal U_\alpha(X=-W,\,x) - \mathcal U_\alpha(X=0,\,x) = 0,\qquad \alpha = 1,\dots,6,\qquad 0\le x\le L, \eqn{pbc}
$$

where $\mathcal U_\alpha$ are the three translations and three rotations of a shell node. This is the constraint set that a shell model imposes with one equation per pair of side nodes and degree of freedom {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. On the skin side edges (normal $Y$), the shell degrees of freedom correspond to the GBT quantities of Table 2.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 2.** Degrees of freedom of a skin side node and the corresponding modal conditions.
:::

| Shell DOF | Meaning at a skin side node | Modal condition |
|:--|:--|:--|
| $\mathcal U_1$ | in-plane displacement $U_X$ | $U_{X,s_L} = U_{X,s_R}$ for every conventional mode, Eq. (@rows) |
| $\mathcal U_2$ | deflection $U_Y$ | $U_{Y,s_L} = U_{Y,s_R}$ for every conventional mode, Eq. (@rows) |
| $\mathcal U_3$ | axial displacement $u$ | $u_k(s_L) = u_k(s_R)$, Eq. (@kkt); $u_j(s_L) = u_j(s_R)$, Eq. (@merge) |
| $\mathcal U_4$ | rotation $\phi_x = -w_{,x}+\beta_x$ | follows from $\mathcal U_2$ for $w_{,x}$; $\xi_r(s_L) = \xi_r(s_R)$, Eq. (@merge) |
| $\mathcal U_6$ | rotation $\theta+\beta_s$ | $\theta_{s_L} = \theta_{s_R}$, Eq. (@rows); $\zeta_r(s_L) = \zeta_r(s_R)$, Eq. (@merge) |
| $\mathcal U_5$ | drilling rotation | not a degree of freedom of the plate theory |

Inserting Eqs. (@disp) and (@beta) into Eq. (@pbc) and collecting the terms of each amplitude function gives, for example for $\mathcal U_3$,

$$
\left[u_k(s_L)-u_k(s_R)\right]\phi_{k,x}(x) + \left[u_j(s_L)-u_j(s_R)\right]\varphi_j(x) = 0\qquad \forall\,x\in[0,L]. \eqn{pbcmodal}
$$

**Proposition.** Equation (@pbc) holds for arbitrary amplitude functions if, and only if, every deformation mode satisfies it individually, i.e. $\mathbf U_k(s_L) = \mathbf U_k(s_R)$, $\theta_k(s_L) = \theta_k(s_R)$ and $u_k(s_L) = u_k(s_R)$ for all conventional modes, $u_j(s_L) = u_j(s_R)$ for all shear modes and $\xi_r(s_L) = \xi_r(s_R)$, $\zeta_r(s_L) = \zeta_r(s_R)$ for all Mindlin modes.

*Proof.* Sufficiency is evident from Eq. (@pbcmodal) and its analogues. For necessity, choose all amplitude functions but one equal to zero: Eq. (@pbcmodal) then requires the bracket of the remaining mode to vanish, since $\phi_{k,x}$ or $\varphi_j$ is not identically zero. $\square$

The proposition is the basis of the present formulation: the periodic condition is a property of the cross-section analysis, not of the member analysis. A basis that satisfies it mode by mode carries the condition into every combination and every truncation of the modes, and the member analysis (Section 2.3) needs neither penalty terms nor Lagrange multipliers. Conversely, a classical GBT basis of an open contour generally violates it mode by mode, so that the condition would have to be imposed as a constraint between amplitude functions.

### Conventional modes: constrained family subspaces

For the conventional modes, the conditions on $\mathcal U_1$, $\mathcal U_2$ and the rotation part of $\mathcal U_6$ are three rows acting on the stacked nodal vector,

$$
\mathbf R\begin{bmatrix}\mathbf U\\ \boldsymbol\Theta\end{bmatrix} = \mathbf 0,\qquad
\mathbf R = \begin{bmatrix}\mathbf e^{\mathsf T}_{X,s_L}-\mathbf e^{\mathsf T}_{X,s_R}\\ \mathbf e^{\mathsf T}_{Y,s_L}-\mathbf e^{\mathsf T}_{Y,s_R}\\ \mathbf e^{\mathsf T}_{\theta,s_L}-\mathbf e^{\mathsf T}_{\theta,s_R}\end{bmatrix}\in\mathbb R^{3\times 3N},\qquad \operatorname{rank}\mathbf R = 3, \eqn{rows}
$$

where $\mathbf e_{X,q}$, $\mathbf e_{Y,q}$ and $\mathbf e_{\theta,q}$ extract $U_{X,q}$, $U_{Y,q}$ and $\theta_q$. The admissible in-plane deformations form the subspace $\ker\mathbf R$ of dimension $3N-3$. Let $\mathbf V_g\in\mathbb R^{3N\times n_g}$ hold the nodal fields $[\mathbf U;\boldsymbol\Theta]$ of family $g\in\{G, D, L, T\}$ of the unconstrained cell, with $\operatorname{span}[\mathbf V_G\ \mathbf V_D\ \mathbf V_L\ \mathbf V_T] = \mathbb R^{3N}$. The periodic modes are constructed family by family, so that every mode keeps the mechanical character of its family:

1. *Unaffected modes.* The columns of $\mathbf V_g$ with $\mathbf R\mathbf v = \mathbf 0$ (e.g. the local modes of nodes away from the side edges) are kept unchanged.
2. *Constrained combinations within a family.* The remaining columns $\mathbf V'_g$ are replaced by $\mathbf V'_g\mathbf Z_g$, with $\mathbf Z_g$ a basis of $\ker(\mathbf R\mathbf V'_g)$. By construction, $\operatorname{span}(\mathbf V_g\mathbf Z_g) = \operatorname{span}(\mathbf V_g)\cap\ker\mathbf R$ is the largest periodic subspace of family $g$.
3. *Combinations of several families.* The intersection of $\ker\mathbf R$ with the sum of the family spaces may contain vectors that do not lie in the sum of the intersections of step 2. Such vectors, $\mathbf P = (\mathbf I - \boldsymbol\Pi_{B})\,\ker\mathbf R$ with the orthogonal projector $\boldsymbol\Pi_{B}$ onto the span of the modes of steps 1 and 2, are added to family T if they elongate a segment ($\mathbf E_\ell\mathbf U\neq\mathbf 0$) and to family D otherwise.
4. *Normalisation.* Every mode is scaled to a unit maximum nodal displacement or rotation.

The resulting modal matrix $\mathbf V = [\mathbf V_G\ \mathbf V_D\ \mathbf V_L\ \mathbf V_T]$ satisfies $\mathbf R\mathbf V = \mathbf 0$ and $\operatorname{span}\mathbf V = \ker\mathbf R$, so that no admissible deformation is lost and none is spurious; the number of conventional modes is $n_C = 3N-3$. Two consequences illustrate the mechanics of the cell. The rigid rotation about the $x$-axis moves the side nodes by $U_{Y,s_L}-U_{Y,s_R} = -W\omega$ and is not admissible: a torsional rigid-body motion of one cell is incompatible with an infinitely wide panel, so that family G contains only the two translations. The lateral translation of the side nodes and their rotation, on the other hand, are admissible if they are identical at both sides; these combinations of local and global fields enter family D through step 3 and yield 13 distortional modes.

### Warping of the conventional modes: the periodic cell as a closed contour

For the warping, $\mathcal U_3$ requires $u_k(s_L) = u_k(s_R)$. In terms of the warping problem, the side edges of the skin are thereby joined: the periodic cell becomes a contour that is closed through its side edges, like the rectangular hollow sections and other closed cross-sections of GBT {{schardt1994,goncalves2013polygon}}. The periodic condition is added to Eq. (@vlasov) as a second constraint. With the segment difference matrix $\mathbf H$, $(\mathbf H\mathbf u)_p = u_{q'}-u_{q}$, the weight matrix $\mathbf W = \mathrm{diag}(\hat A_{66,p}/b_p)$, the vector $\mathbf g$ of the values $b_p\,v_p(\eta = 1/2)$, and the constraint matrix $\mathbf C_u = \left[\mathbf c_{ax}^{\mathsf T};\ \mathbf e^{\mathsf T}_{s_L}-\mathbf e^{\mathsf T}_{s_R}\right]$, the functional of Eq. (@vlasov) reads $\Pi_w = (\mathbf H\mathbf u + \mathbf g)^{\mathsf T}\mathbf W(\mathbf H\mathbf u + \mathbf g)$, and the stationarity of the Lagrangian $\mathcal L = \tfrac12\Pi_w + \boldsymbol\mu^{\mathsf T}\mathbf C_u\mathbf u$ gives the Karush–Kuhn–Tucker system

$$
\begin{bmatrix}\mathbf H^{\mathsf T}\mathbf W\mathbf H & \mathbf C_u^{\mathsf T}\\ \mathbf C_u & \mathbf 0\end{bmatrix}
\begin{bmatrix}\mathbf u\\ \boldsymbol\mu\end{bmatrix} =
\begin{bmatrix}-\mathbf H^{\mathsf T}\mathbf W\mathbf g\\ \mathbf 0\end{bmatrix}. \eqn{kkt}
$$

$\mathbf H^{\mathsf T}\mathbf W\mathbf H$ is the stiffness matrix of a chain of shear springs; its null space is the uniform warping, which is removed by $\mathbf c_{ax}$ (with $\mathbf c_{ax}^{\mathsf T}\mathbf 1 > 0$), so that Eq. (@kkt) has a unique solution. For a mode with non-zero circulation $\oint v\,ds$ around a closed contour (here including the contour closed through the side edges), the residual $\gamma_k$ does not vanish; it is the membrane shear flow of a closed section, carried by the term $A_{66}\gamma_k\gamma_i$ of the modal tensors (Section 2.3).

### Shear and Mindlin modes: merging of the side-edge functions

The shear and Mindlin modes are nodal, so that $\mathcal U_3$, $\mathcal U_4$ and $\mathcal U_6$ are satisfied by a change of basis. With the matrix $\mathbf S$ whose columns are the nodal functions of the unconstrained cell, the periodic basis is $\mathbf S\mathbf P_m$ with the merge matrix $\mathbf P_m$, which adds the column of $s_L$ to that of $s_R$ and removes it:

$$
h_{s_R}\leftarrow h_{s_R}+h_{s_L},\qquad
\boldsymbol\Omega_{s_R}\leftarrow\boldsymbol\Omega_{s_R}+\boldsymbol\Omega_{s_L},\qquad
\zeta_{s_R}\leftarrow\zeta_{s_R}+\zeta_{s_L},\qquad h_{s_L},\ \boldsymbol\Omega_{s_L},\ \zeta_{s_L}\ \text{removed}. \eqn{merge}
$$

A combination $\mathbf S\mathbf a$ of the unconstrained functions satisfies $a_{s_L} = a_{s_R}$ if, and only if, it lies in $\operatorname{span}(\mathbf S\mathbf P_m)$, so that Eq. (@merge) spans exactly the periodic subspace. Since the side nodes lie in the straight skin, their tangent spaces coincide and the vectors $\boldsymbol\Omega_{s_L}$ and $\boldsymbol\Omega_{s_R}$ can be paired. The axial extension mode is uniform and periodic by construction. The number of shear modes, including the axial mode, is $n_S = N-1$, and the number of Mindlin modes is $n_M = \sum_q m_q - m_{s_L} + N - 1$, with $m_q$ the number of vectors $\boldsymbol\Omega_q$ at node $q$.

### Properties of the periodic mode basis

The mode numbers are listed in Table 3. The periodic mode basis has the following properties. (i) *Exactness:* by the proposition of Section 2.2.1, every combination and every truncation of the modes satisfies Eq. (@pbc) identically. (ii) *Completeness:* the modes span exactly the periodic subspace of the nodal description, so that the full basis is equivalent to imposing Eq. (@pbc) node by node. (iii) *Mechanical meaning:* every mode belongs to one of the GBT families, so that the buckling mode can be decomposed into global, distortional, local, transverse-extension, shear and Mindlin contributions. (iv) *Relation to Bloch–Floquet theory:* for the infinite panel with the pitch $b_s$, the buckling modes are Bloch waves, $\mathbf u(X+b_s,x) = e^{\mathrm i\kappa b_s}\,\mathbf u(X,x)$ {{mead1996}}; a periodic cell of width $2b_s$ admits exactly the two real cases $\kappa b_s = 0$ and $\kappa b_s = \pi$, i.e. modes that repeat in every bay or alternate in sign from bay to bay.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 3.** Number of deformation modes per family of the periodic cell (the numbers do not depend on the lay-up).
:::

| $b_s$ [mm] | $N$ | G | D | L | T | S (incl. axial) | M | total |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 150 | 92 | 2 | 13 | 166 | 92 | 91 | 205 | 569 |
| 285.5 | 78 | 2 | 13 | 138 | 78 | 77 | 174 | 482 |

The periodic condition cannot be added to a truncated set of classical modes without losing these properties. Classical GBT generates the conventional modes from the natural nodes of the open contour and orders them hierarchically; the periodic condition couples the two ends of the contour, removes the rigid rotation, closes the contour for the warping and therefore changes the modes themselves. Imposed afterwards on a truncated classical set, the three rows of Eq. (@rows) remove three degrees of freedom from a small space and leave too few admissible modes. In the present formulation the constraint is applied before any truncation, and a hierarchical ordering can be introduced afterwards within the periodic subspace. Following classical GBT, the conventional modes are ordered by the generalised eigenvalue problem $\mathbf B\,\boldsymbol\psi = \Lambda\left(k_r^2\mathbf D^{I} + k_r^4\mathbf C\right)\boldsymbol\psi$ with the modal tensors of Section 2.3 at a reference wavenumber $k_r$; small $\Lambda$ indicates a preference for long waves. The first eight such modes of the periodic cell are shown in Figure 2; the two side nodes move identically in every mode, including the lateral rigid translation (mode 2), while the rigid rotation is absent.

![**Figure 2.** First eight hierarchical conventional modes of the periodic cell ($b_s = 150$ mm, Cross-ply lay-up, $k_r = 2\pi/b_s$; in-plane deformation, grey: undeformed). The side nodes $s_L$ and $s_R$ (squares) move identically in every mode.](fig/fig2_modes.png){width=15.5cm}

## Member analysis

### Strain energy variation and modal tensors

With Eqs. (@strain) and (@claw), the variation of the strain energy of the cell is, with the stresses $\boldsymbol\sigma = [\sigma_{xx}, \sigma_{ss}, \sigma_{xs}, \sigma_{xz}, \sigma_{sz}]^{\mathsf T}$ and the strains $\boldsymbol\varepsilon = [\varepsilon_{xx}, \varepsilon_{ss}, \gamma_{xs}, \gamma_{xz}, \gamma_{sz}]^{\mathsf T}$,

$$
\delta U = \int_L\int_b\int_t\boldsymbol\sigma^{\mathsf T}\delta\boldsymbol\varepsilon\,dz\,ds\,dx = \int_L\int_b\left(\mathbf N^{\mathsf T}\delta\boldsymbol\varepsilon^{M} + \mathbf M^{\mathsf T}\delta\boldsymbol\kappa + \mathbf Q^{\mathsf T}\delta\boldsymbol\gamma_z\right)ds\,dx, \eqn{dU0}
$$

Substituting the modal strains of Table 1 and integrating over the cross-section gives, in the form of Silvestre and Camotim {{silvestre2013shear}},

$$
\begin{gathered}
\delta U = \int_L\Big(C_{ik}\phi_{k,xx}\delta\phi_{i,xx} + H^{I}_{ik}\phi_{k,x}\delta\phi_{i,xx} + H^{I}_{ki}\phi_{k,xx}\delta\phi_{i,x} + D^{I}_{ik}\phi_{k,x}\delta\phi_{i,x}\\
+ D^{II}_{ik}\phi_{k}\delta\phi_{i,xx} + D^{II}_{ki}\phi_{k,xx}\delta\phi_{i} + F^{I}_{ik}\phi_{k,x}\delta\phi_{i} + F^{I}_{ki}\phi_{k}\delta\phi_{i,x} + B_{ik}\phi_{k}\delta\phi_{i}\\
+ C_{hj}\varphi_{j,x}\delta\varphi_{h,x} + H^{I}_{hj}\varphi_{j}\delta\varphi_{h,x} + H^{I}_{jh}\varphi_{j,x}\delta\varphi_{h} + D^{I}_{hj}\varphi_{j}\delta\varphi_{h}\\
+ C_{ij}\varphi_{j,x}\delta\phi_{i,xx} + H^{I}_{ij}\varphi_{j}\delta\phi_{i,xx} + H^{II}_{ji}\varphi_{j,x}\delta\phi_{i,x}\\
+ D^{I}_{ji}\varphi_{j}\delta\phi_{i,x} + D^{II}_{ji}\varphi_{j,x}\delta\phi_{i} + F^{II}_{ij}\varphi_{j}\delta\phi_{i}\\
+ C_{kh}\phi_{k,xx}\delta\varphi_{h,x} + H^{I}_{kh}\phi_{k,xx}\delta\varphi_{h} + H^{II}_{hk}\phi_{k,x}\delta\varphi_{h,x}\\
+ D^{I}_{hk}\phi_{k,x}\delta\varphi_{h} + D^{II}_{hk}\phi_{k}\delta\varphi_{h,x} + F^{II}_{kh}\phi_{k}\delta\varphi_{h}\Big)dx + \delta U_M,
\end{gathered} \eqn{dU}
$$

where the tensors with exchanged indices in the last line are those of the third line with the roles of conventional and shear modes exchanged (e.g. $C_{kh}$ is $C_{ij}$ with $i$ and $j$ replaced by $k$ and $h$), and $\delta U_M$ contains the terms of the Mindlin modes. The letters of the tensors follow {{silvestre2013shear}}: $C$ pairs the axial strains, $B$ the transverse strains and curvatures, $D^{I}$ the shear strains and twisting curvatures, $D^{II}$ the axial with the transverse quantities, $F$ the transverse with the shear quantities and $H$ the axial with the shear quantities. For the conventional modes,

$$
C_{ik} = \int_b\Big[A_{11}\,u_k u_i - B_{11}\left(w_k u_i + w_i u_k\right) + D_{11}\,w_k w_i\Big]ds, \eqn{Cik}
$$

$$
B_{ik} = \int_b\Big[A_{22}\,v_{k,s} v_{i,s} - B_{22}\left(v_{k,s} w_{i,ss} + w_{k,ss} v_{i,s}\right) + D_{22}\,w_{k,ss} w_{i,ss}\Big]ds, \eqn{Bik}
$$

$$
D^{I}_{ik} = \int_b\Big[A_{66}\,\gamma_k\gamma_i - 2B_{66}\left(\gamma_k w_{i,s} + w_{k,s}\gamma_i\right) + 4D_{66}\,w_{k,s} w_{i,s}\Big]ds, \eqn{DIik}
$$

$$
D^{II}_{ik} = \int_b\Big[A_{12}\,v_{k,s} u_i - B_{12}\left(v_{k,s} w_i + w_{k,ss} u_i\right) + D_{12}\,w_{k,ss} w_i\Big]ds, \eqn{DIIik}
$$

$$
F^{I}_{ik} = \int_b\Big[A_{26}\,\gamma_k v_{i,s} - B_{26}\left(\gamma_k w_{i,ss} + 2w_{k,s} v_{i,s}\right) + 2D_{26}\,w_{k,s} w_{i,ss}\Big]ds, \eqn{FIik}
$$

$$
H^{I}_{ik} = \int_b\Big[A_{16}\,\gamma_k u_i - B_{16}\left(\gamma_k w_i + 2w_{k,s} u_i\right) + 2D_{16}\,w_{k,s} w_i\Big]ds; \eqn{HIik}
$$

for the shear modes,

$$
C_{hj} = \int_b A_{11}\,u_j u_h\,ds,\qquad D^{I}_{hj} = \int_b A_{66}\,u_{j,s} u_{h,s}\,ds,\qquad H^{I}_{hj} = \int_b A_{16}\,u_{j,s} u_h\,ds; \eqn{shearT}
$$

and for the coupling between conventional and shear modes,

$$
\begin{gathered}
C_{ij} = \int_b\left(A_{11}u_i - B_{11}w_i\right)u_j\,ds,\qquad
H^{I}_{ij} = \int_b\left(A_{16}u_i - B_{16}w_i\right)u_{j,s}\,ds,\\
H^{II}_{ji} = \int_b\left(A_{16}\gamma_i - 2B_{16}w_{i,s}\right)u_j\,ds,\qquad
D^{I}_{ji} = \int_b\left(A_{66}\gamma_i - 2B_{66}w_{i,s}\right)u_{j,s}\,ds,\\
D^{II}_{ji} = \int_b\left(A_{12}v_{i,s} - B_{12}w_{i,ss}\right)u_j\,ds,\qquad
F^{II}_{ij} = \int_b\left(A_{26}v_{i,s} - B_{26}w_{i,ss}\right)u_{j,s}\,ds .
\end{gathered} \eqn{couplingT}
$$

For conventional modes that satisfy Vlasov's assumptions ($\gamma_k = 0$, $v_{k,s} = 0$), Eqs. (@Cik)–(@couplingT) are exactly the tensors of Silvestre and Camotim (Eqs. (10)–(12) of {{silvestre2013shear}}); the membrane terms with $\gamma_k$ and $v_{k,s}$ are the additions required by the closed contours and the transverse extension. The Mindlin terms follow in the same way from Table 1; they contain the transverse shear stiffness, e.g. $\int_b A_{s,11}\,\xi_r\xi_{r'}\,ds$ for two $\beta_x$ modes, and the stiffnesses $\mathbf D$ and $\mathbf B$ through the curvatures. In the matrix notation of the GBT–Ritz approach {{kharghani2026jcs}}, the tensors of the conventional modes form the GBT matrix

$$
\boldsymbol\Xi = \begin{bmatrix}\mathbf B & \mathbf F^{I} & \mathbf D^{II\mathsf T}\\ \mathbf F^{I\mathsf T} & \mathbf D^{I} & \mathbf H^{I\mathsf T}\\ \mathbf D^{II} & \mathbf H^{I} & \mathbf C\end{bmatrix},\qquad
\delta U = \int_L \delta\boldsymbol\Phi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Phi\,dx,\qquad
\boldsymbol\Phi = \begin{bmatrix}\boldsymbol\phi\\ \boldsymbol\phi_{,x}\\ \boldsymbol\phi_{,xx}\end{bmatrix}, \eqn{Xi}
$$

with the row index $i$ and the column index $k$ in every block, and it is extended by the shear, coupling and Mindlin tensors in the same arrangement for the complete vector of amplitude functions $\boldsymbol\Phi$.

### Second-order energy and geometric tensors

The pre-buckling state under the reference load is described by the membrane forces $S_{xx}(x,s)$, $S_{ss}(x,s)$ and $S_{xs}(x,s)$ (Section 2.3.3). Their second-order work on the non-linear strains of Eq. (@strainNL) is

$$
V = \int_L\int_b\left(S_{xx}\,\varepsilon^{NL}_{xx} + S_{ss}\,\varepsilon^{NL}_{ss} + S_{xs}\,\gamma^{NL}_{xs}\right)ds\,dx, \eqn{V}
$$

and its variation defines the geometric tensors. For the conventional modes the displacement gradients are $u_{,x} = u_k\phi_{k,xx}$, $u_{,s} = u_{k,s}\phi_{k,x}$, $v_{,x} = v_k\phi_{k,x}$, $v_{,s} = v_{k,s}\phi_k$, $w_{,x} = w_k\phi_{k,x}$ and $w_{,s} = w_{k,s}\phi_k$, and

$$
\begin{gathered}
\delta V = \int_L\Big(X^{U}_{ik}\,\phi_{k,xx}\delta\phi_{i,xx} + X^{VW}_{ik}\,\phi_{k,x}\delta\phi_{i,x} + X^{S}_{ik}\,\phi_{k}\delta\phi_{i}\\
+ X^{T}_{ik}\,\phi_{k,x}\delta\phi_{i} + X^{T}_{ki}\,\phi_{k}\delta\phi_{i,x} + X^{R}_{ik}\,\phi_{k,xx}\delta\phi_{i,x} + X^{R}_{ki}\,\phi_{k,x}\delta\phi_{i,xx}\Big)dx + \dots,\\
X^{U}_{ik} = \int_b S_{xx}\,u_i u_k\,ds,\qquad
X^{VW}_{ik} = \int_b\left[S_{xx}\left(v_i v_k + w_i w_k\right) + S_{ss}\,u_{i,s}u_{k,s}\right]ds,\\
X^{S}_{ik} = \int_b S_{ss}\left(v_{i,s}v_{k,s} + w_{i,s}w_{k,s}\right)ds,\qquad
X^{T}_{ik} = \int_b S_{xs}\left(v_{i,s}v_k + w_{i,s}w_k\right)ds,\\
X^{R}_{ik} = \int_b S_{xs}\,u_{i,s}u_k\,ds,
\end{gathered} \eqn{Xgeo}
$$

where the dots stand for the analogous terms of the shear modes ($u_{,x} = u_j\varphi_{j,x}$, $u_{,s} = u_{j,s}\varphi_j$) and of their coupling with the conventional modes; the Mindlin modes do not move the mid-surface and do not contribute. Since the membrane forces vary along $x$, the geometric tensors are functions of $x$. For a uniform axial membrane force $S_{xx}$ and $S_{ss} = S_{xs} = 0$, only $X^{U}_{ik}$ and $X^{VW}_{ik} = \int_b S_{xx}(v_iv_k + w_iw_k)\,ds$ remain; these are the geometric matrices $\boldsymbol{\mathcal X}^{U}$ and $\boldsymbol{\mathcal X}^{VW}$ of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}, and $X^{VW}_{ik}$ is the geometric tensor $X_{pik}$ of classical GBT {{silvestre2002b,silvestre2013shear}} for the axial resultant.

### Loads and pre-buckling state

The cell is loaded at both loaded edges by a compressive line load $n_0(s)$ per unit length of the mid-line, which acts at the distance $e_0(s)$ from the wall reference surface (Section 3). The variation of the potential of the external forces is

$$
\delta\Pi_{ext} = -\int_b n_0\left[\delta u_0(0,s) - \delta u_0(L,s)\right]ds,\qquad u_0 = u - e_0\,w_{,x}, \eqn{load0}
$$

$$
\begin{gathered}
\delta\Pi_{ext} = -\left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=0} + \left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=L},\\
q^{I}_k = \int_b n_0\left(u_k - e_0 w_k\right)ds,\qquad q^{I}_h = \int_b n_0\,u_h\,ds,
\end{gathered} \eqn{load}
$$

i.e. the end load enters through the generalised end forces $q^{I}_k$ and $q^{I}_h$ associated with the warping of the modes, as the forces $q^{I}_i$ in {{silvestre2013shear}}. The pre-buckling state is the first-order solution of $\delta U + \delta\Pi_{ext} = 0$ under the end and side conditions of the cell. Its membrane forces

$$
\left[S_{xx},\ S_{ss},\ S_{xs}\right]^{\mathsf T} = \mathbf A\,\boldsymbol\varepsilon^{M,0} + \mathbf B\,\boldsymbol\kappa^{0} \eqn{pre}
$$

vary along $s$ and $x$: the load enters the stringer through the eccentric foot and is shared between skin and stringer by shear lag, the eccentric foot load bends the cell near the ends, and the restrained Poisson expansion produces transverse membrane forces $S_{ss}$. In classical GBT the pre-buckling state is expressed by the modal resultants $W^{o}_p$ of an assumed stress distribution {{silvestre2013shear}}; here it is the computed first-order state, which is necessary because the load is introduced eccentrically at the ends of a stiffened cell.

### Equilibrium equations and boundary conditions

The buckling problem follows from the stationarity of the total potential energy of the adjacent state, $\delta U + \lambda\,\delta V = 0$, where $\lambda$ multiplies the reference load. Integrating Eq. (@dU) by parts and collecting the coefficients of $\delta\phi_i$ and $\delta\varphi_h$ gives the GBT equilibrium equations, written for the terms of the conventional and shear modes:

$$
\begin{gathered}
\mathcal K_{ik}\,\phi_k + \mathcal K_{ij}\,\varphi_j + \lambda\left(\mathcal G_{ik}\,\phi_k + \mathcal G_{ij}\,\varphi_j\right) = 0,\\
\mathcal K_{hk}\,\phi_k + \mathcal K_{hj}\,\varphi_j + \lambda\left(\mathcal G_{hk}\,\phi_k + \mathcal G_{hj}\,\varphi_j\right) = 0,
\end{gathered} \eqn{equil}
$$

with the stiffness operators

$$
\begin{gathered}
\mathcal K_{ik} = C_{ik}\frac{d^4}{dx^4} + H_{ik}\frac{d^3}{dx^3} - D_{ik}\frac{d^2}{dx^2} + F_{ik}\frac{d}{dx} + B_{ik},\\
\mathcal K_{ij} = C_{ij}\frac{d^3}{dx^3} + \left(H^{I}_{ij} - H^{II}_{ji}\right)\frac{d^2}{dx^2} - \left(D^{I}_{ji} - D^{II}_{ji}\right)\frac{d}{dx} + F^{II}_{ij},\\
\mathcal K_{hk} = C_{kh}\frac{d^3}{dx^3} + \left(H^{II}_{hk} - H^{I}_{kh}\right)\frac{d^2}{dx^2} - \left(D^{I}_{hk} - D^{II}_{hk}\right)\frac{d}{dx} - F^{II}_{kh},\\
\mathcal K_{hj} = C_{hj}\frac{d^2}{dx^2} + H_{hj}\frac{d}{dx} - D^{I}_{hj},
\end{gathered} \eqn{Kop}
$$

where the second equation of Eq. (@equil) has been multiplied by $-1$, as in {{silvestre2013shear}}, and

$$
H_{ik} = H^{I}_{ik} - H^{I}_{ki},\qquad D_{ik} = D^{I}_{ik} - D^{II}_{ik} - D^{II}_{ki},\qquad F_{ik} = F^{I}_{ik} - F^{I}_{ki},\qquad H_{hj} = H^{I}_{hj} - H^{I}_{jh} . \eqn{combined}
$$

For the uniform membrane force of classical GBT, the geometric operator of the conventional modes is $\mathcal G_{ik} = X^{U}_{ik}\,d^4/dx^4 - X^{VW}_{ik}\,d^2/dx^2$; in general it follows from Eq. (@Xgeo) with coefficients that depend on $x$. The boundary terms of the integration by parts give the boundary conditions at $x = 0$ and $x = L$,

$$
\left[\left(W^{\sigma}_i + q^{I}_i\right)\delta\phi_{i,x} + W^{\tau}_i\,\delta\phi_i + \left(W^{\sigma}_h + q^{I}_h\right)\delta\varphi_h\right]_{x=0}^{x=L} = 0, \eqn{bc}
$$

with the generalised normal stress resultants

$$
\begin{gathered}
W^{\sigma}_i = C_{ik}\phi_{k,xx} + H^{I}_{ik}\phi_{k,x} + D^{II}_{ik}\phi_{k} + C_{ij}\varphi_{j,x} + H^{I}_{ij}\varphi_{j},\\
W^{\sigma}_h = C_{hj}\varphi_{j,x} + H^{I}_{hj}\varphi_{j} + C_{kh}\phi_{k,xx} + H^{II}_{hk}\phi_{k,x} + D^{II}_{hk}\phi_k,
\end{gathered} \eqn{Wsig}
$$

and the generalised shear stress resultant $W^{\tau}_i = -W^{\sigma}_{i,x} + W^{xs}_i$ with $W^{xs}_i = H^{I}_{ki}\phi_{k,xx} + D^{I}_{ik}\phi_{k,x} + F^{I}_{ki}\phi_k + H^{II}_{ji}\varphi_{j,x} + D^{I}_{ji}\varphi_j$ (geometric terms included in the buckling problem). Equations (@equil)–(@Wsig) are the equations of Silvestre and Camotim {{silvestre2013shear}}, generalised by the membrane terms of the tensors, the Mindlin modes and the non-uniform pre-buckling state. Two features of the panel cell, however, prevent a closed-form or harmonic solution of Eq. (@equil): the geometric operators have coefficients that vary along $x$, and the end conditions of a stiffened cell involve combinations of modes (Section 2.3.5). The equations are therefore solved by the Ritz method.

### Ritz solution and end conditions

As in the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}, $\boldsymbol\Phi$ is written as the product of a matrix $\boldsymbol\Psi$ of shape functions and the vector $\mathbf a$ of unknown coefficients, $\boldsymbol\Phi = \boldsymbol\Psi\,\mathbf a$, with every amplitude function approximated as

$$
\phi_k(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{kb},\qquad \varphi_j(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{jb},\qquad \chi_r(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{rb} . \eqn{Psi}
$$

For members, trigonometric series are used {{kharghani2026jcs}}. In the panel cell, the pre-buckling state varies strongly near the loaded edges, where the load enters the stringer, and the clamped condition has to be combined with free warping at one end. The shape functions $N_b$ are therefore cubic B-splines on an open, uniform knot vector with spans of length $h_x = \min(15\ \text{mm},\ \max(2\ \text{mm},\ W/40),\ L/10)$, i.e. $h_x = b_s/20$ for all panels studied. With the open knot vector, $\phi(0) = a_1$, $\phi(L) = a_{n_b}$ and $\phi_{,x}(0) = 0 \Leftrightarrow a_2 = a_1$, so that the essential end conditions act directly on the end coefficients, and the natural conditions of Eq. (@bc) are satisfied in the weak sense. Two end conditions of the loaded edges are considered, as in the shell models:

* **Simply supported (SS):** the out-of-plane displacement of every end edge vanishes, $U_Y = 0$ at every node at $x = 0$ and $x = L$, i.e. $\sum_k U_{Y,qk}\,\phi_k(0) = \sum_k U_{Y,qk}\,\phi_k(L) = 0$ for every node $q$. These linear relations between the end coefficients of the conventional modes are imposed by a null-space reduction; the in-plane displacement $U_X$, the warping and the rotations of the end cross-sections remain free. Rigid-body motions are removed by fixing the axial mode and $U_X$ of the side node $s_R$ at $x = 0$.
* **Clamped (CC):** $\phi_k = \phi_{k,x} = 0$ at both ends for every conventional mode (in-plane displacements and the rotations $w_{,x}$ and $\theta$ vanish), $\chi_r = 0$ at both ends, and $\varphi_j(0) = 0$, so that $u = 0$ on the whole edge $x = 0$; the edge $x = L$ is free to warp under the uniform load.

The classical simply supported condition of GBT, $\phi_k = \sin(m\pi x/L)$, suppresses all in-plane displacements of the end cross-sections while leaving their warping free; the present SS condition is weaker and suppresses only the out-of-plane displacement $U_Y$ of the end edges, as in the shell model. Inserting Eq. (@Psi) into $\delta U + \lambda\,\delta V = 0$ gives the eigenvalue problem

$$
\left(\mathbf K + \lambda\,\mathbf G\right)\mathbf a = \mathbf 0,\qquad
\mathbf K = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Psi\,dx,\qquad
\mathbf G = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol{\mathcal X}(x)\,\boldsymbol\Psi\,dx, \eqn{eig}
$$

with $\boldsymbol\Xi$ the complete matrix of the modal tensors and $\boldsymbol{\mathcal X}(x)$ that of the geometric tensors of Eq. (@Xgeo). The integrals are evaluated with four Gauss points per span. The smallest positive eigenvalue multiplies the reference load of 1 N/mm, so that $N_{cr} = \lambda\times1$ N/mm, and $\mathbf a$ defines the buckling mode. Since $\mathbf G$ is indefinite (the pre-buckling state contains membrane shear and some transverse tension), the problem is solved for the largest $\mu = 1/\lambda$ of $-\mathbf G\mathbf a = \mu\mathbf K\mathbf a$ by the implicitly restarted Lanczos method {{lehoucq1998}} with a sparse factorisation of $\mathbf K$. For $b_s = 150$ mm and $L = 225$ mm the model has $569\times33 = 18\,777$ coefficients before the end conditions are applied. Two numerical aspects are noted. First, at some intermediate nodes of the inclined webs the two segment tangents coincide up to round-off, and the orthonormal basis of the tangents contains a spurious Mindlin vector normal to the web, with $\xi_r = \mathbf t\cdot\boldsymbol\Omega\approx10^{-15}$ (11 such modes for $b_s = 150$ mm and 8 for $b_s = 285.5$ mm); these modes carry neither energy nor displacement gradients and are held. Second, at a cross-section integration point only about 25 of the 569 modes are non-zero, so that the cross-section factors of $\boldsymbol{\mathcal X}$ are stored as sparse matrices that are computed once per cross-section.

### Classical limit

For the conventional modes alone, walls without coupling ($\mathbf B = \mathbf 0$, $A_{16} = A_{26} = D_{16} = D_{26} = 0$) and a uniform membrane force $S_{xx}$, the tensors $H_{ik}$ and $F_{ik}$ vanish and Eq. (@equil) becomes

$$
C_{ik}\,\phi_{k,xxxx} - D_{ik}\,\phi_{k,xx} + B_{ik}\,\phi_k + \lambda\left(X^{U}_{ik}\,\phi_{k,xxxx} - X^{VW}_{ik}\,\phi_{k,xx}\right) = 0, \eqn{classical}
$$

which, without the warping term $X^{U}_{ik}$, is the classical second-order GBT equation of Schardt {{schardt1994}} and of Silvestre and Camotim {{silvestre2002b}}. For simply supported ends and $\phi_k = \hat\phi_k\sin(m\pi x/L)$, with $k_m = m\pi/L$,

$$
\left[B_{ik} + k_m^2\,D_{ik} + k_m^4\,C_{ik} + \lambda\left(k_m^2\,X^{VW}_{ik} + k_m^4\,X^{U}_{ik}\right)\right]\hat\phi_k = 0, \eqn{signature}
$$

whose smallest eigenvalue over $m$ is the classical GBT signature curve. The periodic condition does not change this structure; it restricts the admissible cross-section modes. Equation (@signature) was used to verify the implementation (Section 3.1).

## Mode decomposition

The eigenvector $\mathbf a$ defines the buckling mode in terms of the GBT modes. As in the GBT–Ritz approach {{kharghani2026jcs}}, an energy-based mode decomposition quantifies the contribution of the deformation mechanisms. The modes are grouped into the families $g\in\{G, D, L, T, S, M\}$, and $\mathbf a_g$ denotes the vector $\mathbf a$ with all coefficients set to zero except those of family $g$. The strain energy of the buckling mode, $U = \tfrac12\mathbf a^{\mathsf T}\mathbf K\mathbf a$, is divided into the contributions of the families,

$$
U_g = \tfrac12\,\mathbf a_g^{\mathsf T}\mathbf K\,\mathbf a,\qquad \rho_g = \frac{U_g}{U},\qquad \sum_g\rho_g = 1 . \eqn{rho}
$$

If the families are orthogonal with respect to $\mathbf K$, $U_g$ equals the energy of family $g$ alone, $\tfrac12\mathbf a_g^{\mathsf T}\mathbf K\mathbf a_g$, the form used for members {{kharghani2026jcs}}. In the panel cell the distortional and local modes are not orthogonal, because the D modes move the wall junctions without rotating the nodes and the rotation modes of family L restore the continuity of the walls. The energies of the families alone then overestimate the coupled families and do not add up to $U$, whereas Eq. (@rho) divides the strain energy exactly. For the Unsym cell with simply supported loaded edges, $b_s = 150$ mm and $\beta = 1.5$, the energies of the families alone would assign 41 % to the distortional and 58 % to the local family, whereas Eq. (@rho) gives 1 % and 97 %.
