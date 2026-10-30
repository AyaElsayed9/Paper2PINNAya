
# Theory and formulations

The formulation follows the two stages of GBT {{schardt1994,silvestre2002a,bebiano2015}}: a cross-section analysis, which determines the deformation modes and their modal tensors, and a member analysis, which determines the amplitude functions. The notation is that of the shear-deformable GBT of Silvestre and Camotim {{silvestre2013shear}}, and the member analysis follows the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}. Three extensions are introduced such that classical GBT is recovered as a special case (Section 2.3.4): the modes satisfy the periodic boundary conditions of the cell (Section 2.2), Vlasov's assumptions are not imposed on a cell with closed contours, and the walls are shear-deformable laminates.

## Cross-section analysis

### Periodic cell and displacement field

The cell of the infinitely wide panel with identical omega stringers at the pitch $b_s$ contains two stringers; it has the width $W = 2b_s$ and the length $L = \beta b_s$ (Figure 1). $X$ runs across the panel, $Y$ normal to the skin and $x$ along the stringers; the skin occupies $-W \le X \le 0$, with the side edges $s_L$ ($X = -W$) and $s_R$ ($X = 0$) in the middle of a skin bay. Every wall has the local axes $(x, s, z)$ and the displacements $u$, $v$, $w$ {{silvestre2013shear}}. The walls (skin, feet, webs, crowns) are divided into segments of at most $b_s/25$, with natural nodes at the wall ends and intermediate nodes in between {{schardt1994,bebiano2015}}; for $b_s = 150$ mm the cell has $N = 92$ nodes, 19 walls and 93 segments.

![**Figure 1.** Periodic cell of the omega-stringer-stiffened panel. (a) Oblique view of the cell of width $W = 2b_s$ and length $L = \beta b_s$ with the compressive end load $n_0$ and the two side edges $s_L$ and $s_R$ that are linked by the periodic condition. (b) Infinitely wide panel with identical stringers at the pitch $b_s$ and the cell (shaded). (c) GBT discretisation of the cell for $b_s = 150$ mm with natural and intermediate nodes, wall types and local axes $(x, s, z)$ with the displacements $(u, v, w)$.](fig/fig01_geometry.png){width=15.5cm}

The displacement field of shear-deformable GBT {{silvestre2013shear}},

$$
u(x,s) = u_k(s)\,\phi_{k,x}(x) + u_j(s)\,\varphi_j(x),\qquad
v(x,s) = v_k(s)\,\phi_k(x),\qquad
w(x,s) = w_k(s)\,\phi_k(x), \eqn{disp}
$$

is extended by the transverse shear rotations of the Reissner–Mindlin walls,

$$
\beta_x(x,s) = \xi_r(s)\,\chi_r(x),\qquad \beta_s(x,s) = \zeta_r(s)\,\chi_r(x), \eqn{beta}
$$

with summation over the conventional modes $k$ (virtual: $i$), the shear modes $j$ (virtual: $h$) and the transverse shear (Mindlin) modes $r$. $u_k$, $v_k$, $w_k$ are the profiles of conventional mode $k$ with amplitude $\phi_k(x)$, $u_j$ the warping profile of shear mode $j$ with amplitude $\varphi_j(x)$, and $\xi_r$, $\zeta_r$ the rotation profiles of Mindlin mode $r$ with amplitude $\chi_r(x)$. The uniform axial extension is the shear-type mode $u_j = 1$.

### Strains and constitutive relations

With the curvatures $\kappa_{xx} = -w_{,xx} + \beta_{x,x}$, $\kappa_{ss} = -w_{,ss} + \beta_{s,s}$, $\kappa_{xs} = -2w_{,xs} + \beta_{x,s} + \beta_{s,x}$ and the transverse shear strains $\gamma_{xz} = \beta_x$, $\gamma_{sz} = \beta_s$, Eqs. (@disp) and (@beta) give the strains

$$
\begin{gathered}
\varepsilon_{xx} = \left(u_k - z\,w_k\right)\phi_{k,xx} + u_j\,\varphi_{j,x} + z\,\xi_r\,\chi_{r,x},\qquad
\varepsilon_{ss} = \left(v_{k,s} - z\,w_{k,ss}\right)\phi_k + z\,\zeta_{r,s}\,\chi_r,\\
\gamma_{xs} = \left(\gamma_k - 2z\,w_{k,s}\right)\phi_{k,x} + u_{j,s}\,\varphi_j + z\left(\xi_{r,s}\,\chi_r + \zeta_r\,\chi_{r,x}\right),\qquad
\gamma_{xz} = \xi_r\,\chi_r,\qquad \gamma_{sz} = \zeta_r\,\chi_r,
\end{gathered} \eqn{strain}
$$

with the membrane shear profile $\gamma_k = u_{k,s} + v_k$. In {{silvestre2013shear}} the conventional modes satisfy Vlasov's assumptions $\gamma_k = 0$ and $v_{k,s} = 0$, and Eq. (@strain) reduces to their kinematics ($\xi_r = \zeta_r = 0$). Here neither is imposed: transverse extension modes are part of the basis, and the omega stringers form closed cells with the skin, which with the periodic condition itself behaves as a closed contour (Section 2.2.3), so that $\gamma_k = 0$ cannot hold for all distortional modes. The non-linear membrane strains of the geometric stiffness are

$$
\varepsilon^{NL}_{xx} = \tfrac12\left(u_{,x}^2+v_{,x}^2+w_{,x}^2\right),\qquad
\varepsilon^{NL}_{ss} = \tfrac12\left(u_{,s}^2+v_{,s}^2+w_{,s}^2\right),\qquad
\gamma^{NL}_{xs} = u_{,x}u_{,s}+v_{,x}v_{,s}+w_{,x}w_{,s}, \eqn{strainNL}
$$

of which classical GBT retains $\tfrac12(v_{,x}^2+w_{,x}^2)$ {{silvestre2002b,silvestre2013shear}}. The stress resultants of a wall follow from classical lamination theory with a shear correction factor of 5/6 {{mindlin1951,whitney1970,reddy2004}},

$$
\begin{bmatrix}\mathbf N\\ \mathbf M\\ \mathbf Q\end{bmatrix} =
\begin{bmatrix} \mathbf A & \mathbf B & \mathbf 0\\ \mathbf B & \mathbf D & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \mathbf A_s\end{bmatrix}
\begin{bmatrix}\boldsymbol\varepsilon^{M}\\ \boldsymbol\kappa\\ \boldsymbol\gamma_z\end{bmatrix}, \eqn{claw}
$$

retaining all entries of $\mathbf A$, $\mathbf B$, $\mathbf D$, i.e. the extension–shear and bending–twisting coupling {{silvestre2013shear}} and the membrane–bending coupling of unsymmetric laminates {{kharghani2026jcs}}; $N_{ss}$ is not condensed out. A stringer foot (skin and flange) forms one wall referred to the skin mid-plane, with the stiffnesses $\mathbf A_0$, $\mathbf B_0$, $\mathbf D_0$ of the combined stack and the distance $z_f$ between the mid-planes,

$$
\mathbf A_f = \mathbf A_0,\qquad \mathbf B_f = \mathbf B_0 + z_f\,\mathbf A_0,\qquad \mathbf D_f = \mathbf D_0 + 2z_f\,\mathbf B_0 + z_f^2\,\mathbf A_0 . \eqn{foot}
$$

The stringer mid-surface is offset from the skin by $c = (t_s+t_f)/2$. Every web root is therefore a node at $Y = c$, joined to the skin node below it by a link wall that transmits displacements and rotation like rigid offset links without adding axial or twisting stiffness,

$$
\begin{gathered}
A^{L}_{11} = 10^{-3}A_{0,11},\quad A^{L}_{22} = k_L A_{0,22},\quad A^{L}_{66} = k_L A_{0,66},\\
D^{L}_{11} = 10^{-3}D_{0,11},\quad D^{L}_{22} = k_L D_{0,22},\quad D^{L}_{66} = 10^{-3}D_{0,66},\quad \mathbf A^{L}_s = k_L\mathbf A_{s0},
\end{gathered} \eqn{link}
$$

with the skin stiffnesses $\mathbf A_0$, $\mathbf D_0$, $\mathbf A_{s0}$, $k_L = 10^3$ and all other entries zero; the results do not depend on $k_L$ (Section 3.7).

### Nodal description and mode families

On every segment of width $b_p$ ($\eta\in[0,1]$), $u$ and $v$ vary linearly, and $w$ is the cubic Hermite interpolation of the nodal normal displacements and rotations $\theta = w_{,s}$,

$$
\begin{gathered}
w(\eta) = H_1(\eta)\,w_q + b_p H_2(\eta)\,\theta_q + H_3(\eta)\,w_{q'} + b_p H_4(\eta)\,\theta_{q'},\\
H_1 = 1-3\eta^2+2\eta^3,\quad H_2 = \eta-2\eta^2+\eta^3,\quad H_3 = 3\eta^2-2\eta^3,\quad H_4 = -\eta^2+\eta^3 .
\end{gathered} \eqn{hermite}
$$

which spans the same space as Schardt's nodal moments {{schardt1994}} but turns the periodic conditions into linear relations between nodal quantities. With the in-plane nodal displacements $\mathbf U\in\mathbb R^{2N}$, the nodal rotations $\boldsymbol\Theta\in\mathbb R^N$ (shared by the walls at a node) and the segment elongation matrix $\mathbf E_\ell$, $(\mathbf E_\ell\mathbf U)_p = \mathbf t_p\cdot(\mathbf U_{q'}-\mathbf U_q)$, the modes are generated in the families of Bebiano et al. {{bebiano2015}} and Gonçalves et al. {{goncalves2014shear}}: **G** (rigid translations and rotation), **D** (nodal displacements with $\mathbf E_\ell\mathbf U = \mathbf 0$, orthogonal to G and L), **L** (unit normal displacement of a non-junction node and unit rotation of a node), **T** (columns of $\mathbf E_\ell^{\mathsf T}$), **S** (nodal hat functions of the warping, one replaced by the axial mode $u_j = 1$) and **M** (nodal vectors $\boldsymbol\Omega_q$ in the cross-section plane with $\xi_r = \mathbf t\cdot\boldsymbol\Omega_q$, and nodal values of $\zeta_r$, both interpolated linearly). Since $\mathbb R^{2N} = \ker\mathbf E_\ell \oplus \operatorname{range}\mathbf E_\ell^{\mathsf T}$, G, D, L and T span all in-plane deformations of the discretised section. Because Vlasov's condition $\gamma_k = 0$ {{vlasov1961,schardt1994}} cannot be met in closed contours, the warping of the G and D modes is the stiffness-weighted least-squares solution

$$
\min_{\mathbf u}\ \Pi_w(\mathbf u) = \sum_p \hat A_{66,p}\,b_p\left(\frac{u_{q'}-u_{q}}{b_p}+v_p\!\left(\eta = \tfrac12\right)\right)^2
\quad\text{subject to}\quad \mathbf c_{ax}^{\mathsf T}\mathbf u = 0, \eqn{vlasov}
$$

with $\hat A_{66} = A_{66}-A_{26}^2/A_{22}$ and $c_{ax,q} = \sum_{p\ni q}\tfrac12\hat A_{11,p}b_p$, $\hat A_{11} = A_{11}-A_{12}^2/A_{22}$, so that the warping has no axial resultant. For an open contour the minimum is zero and the classical warping is recovered; in closed contours the residual $\gamma_k$ is retained in Eq. (@strain).

## Periodic boundary conditions in the cross-section analysis

### Continuous condition and its modal form

As one period of the infinite panel, the cell must satisfy, for every $x$,

$$
\mathcal U_\alpha(X=-W,\,x) - \mathcal U_\alpha(X=0,\,x) = 0,\qquad \alpha = 1,\dots,6,\qquad 0\le x\le L, \eqn{pbc}
$$

where $\mathcal U_\alpha$ are the translations and rotations of a shell node; this is the constraint set of a shell model with one equation per pair of side nodes and degree of freedom {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. Table 1 relates these degrees of freedom to the GBT quantities.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 1.** Degrees of freedom of a skin side node and the corresponding modal conditions.
:::

| Shell DOF | Meaning at a skin side node | Modal condition |
|:--|:--|:--|
| $\mathcal U_1$ | in-plane displacement $U_X$ | $U_{X,s_L} = U_{X,s_R}$ for every conventional mode, Eq. (@rows) |
| $\mathcal U_2$ | deflection $U_Y$ | $U_{Y,s_L} = U_{Y,s_R}$ for every conventional mode, Eq. (@rows) |
| $\mathcal U_3$ | axial displacement $u$ | $u_k(s_L) = u_k(s_R)$, Eq. (@kkt); $u_j(s_L) = u_j(s_R)$, Eq. (@merge) |
| $\mathcal U_4$ | rotation $\phi_x = -w_{,x}+\beta_x$ | follows from $\mathcal U_2$ for $w_{,x}$; $\xi_r(s_L) = \xi_r(s_R)$, Eq. (@merge) |
| $\mathcal U_6$ | rotation $\theta+\beta_s$ | $\theta_{s_L} = \theta_{s_R}$, Eq. (@rows); $\zeta_r(s_L) = \zeta_r(s_R)$, Eq. (@merge) |
| $\mathcal U_5$ | drilling rotation | not a degree of freedom of the plate theory |

Inserting Eqs. (@disp) and (@beta) into Eq. (@pbc) gives, e.g. for $\mathcal U_3$,

$$
\left[u_k(s_L)-u_k(s_R)\right]\phi_{k,x}(x) + \left[u_j(s_L)-u_j(s_R)\right]\varphi_j(x) = 0\qquad \forall\,x\in[0,L]. \eqn{pbcmodal}
$$

**Proposition.** Equation (@pbc) holds for arbitrary amplitude functions if, and only if, every mode satisfies it individually: $\mathbf U_k$, $\theta_k$ and $u_k$ of all conventional modes, $u_j$ of all shear modes and $\xi_r$, $\zeta_r$ of all Mindlin modes take equal values at $s_L$ and $s_R$.

*Proof.* Sufficiency follows from Eq. (@pbcmodal) and its analogues. For necessity, set all amplitude functions but one to zero; the bracket of the remaining mode must then vanish. $\square$

The periodic condition is thus a property of the cross-section analysis. A mode-wise periodic basis carries it into every combination and truncation of the modes: every trial function of the Ritz solution (Section 2.3.4) is kinematically admissible with respect to the side edges, and the member analysis needs neither penalty terms nor Lagrange multipliers, whereas a classical basis of the open contour violates the condition mode by mode.

### Conventional modes: constrained family subspaces

For the conventional modes, the conditions on $\mathcal U_1$, $\mathcal U_2$ and the rotation part of $\mathcal U_6$ are three rows acting on the stacked nodal vector,

$$
\mathbf R\begin{bmatrix}\mathbf U\\ \boldsymbol\Theta\end{bmatrix} = \mathbf 0,\qquad
\mathbf R = \begin{bmatrix}\mathbf e^{\mathsf T}_{X,s_L}-\mathbf e^{\mathsf T}_{X,s_R}\\ \mathbf e^{\mathsf T}_{Y,s_L}-\mathbf e^{\mathsf T}_{Y,s_R}\\ \mathbf e^{\mathsf T}_{\theta,s_L}-\mathbf e^{\mathsf T}_{\theta,s_R}\end{bmatrix}\in\mathbb R^{3\times 3N},\qquad \operatorname{rank}\mathbf R = 3, \eqn{rows}
$$

whose null space, of dimension $3N-3$, contains the admissible in-plane deformations. With $\mathbf V_g$ the nodal fields of family $g\in\{G,D,L,T\}$ of the unconstrained cell, the periodic modes are built family by family: (1) columns with $\mathbf R\mathbf v = \mathbf 0$ are kept; (2) the remaining columns $\mathbf V'_g$ are replaced by $\mathbf V'_g\mathbf Z_g$ with $\mathbf Z_g$ a basis of $\ker(\mathbf R\mathbf V'_g)$, which spans the largest periodic subspace of family $g$; (3) the part of $\ker\mathbf R$ not reached by steps 1–2, $(\mathbf I - \boldsymbol\Pi_B)\ker\mathbf R$, is added to T if it elongates a segment and to D otherwise; (4) every mode is scaled to a unit maximum nodal value. The modal matrix then satisfies $\mathbf R\mathbf V = \mathbf 0$ and $\operatorname{span}\mathbf V = \ker\mathbf R$, so that no admissible deformation is lost and none is spurious; $n_C = 3N-3$. The rigid rotation about $x$ moves the side nodes by $U_{Y,s_L}-U_{Y,s_R} = -W\omega$ and is excluded, since a torsional rigid-body motion of one cell is incompatible with an infinite panel; family G contains the two translations only. Equal lateral translations and rotations of the side nodes are admissible and enter family D in step 3 (13 distortional modes).

### Warping: the periodic cell as a closed contour

The condition $u_k(s_L) = u_k(s_R)$ joins the side edges of the skin, so that for the warping the cell becomes a contour closed through its side edges, like the closed sections of GBT {{schardt1994,goncalves2013polygon}}. With the segment difference matrix $\mathbf H$, the weights $\mathbf W = \mathrm{diag}(\hat A_{66,p}/b_p)$, the vector $\mathbf g$ of $b_p v_p(\eta = 1/2)$ and $\mathbf C_u = [\mathbf c_{ax}^{\mathsf T};\ \mathbf e^{\mathsf T}_{s_L}-\mathbf e^{\mathsf T}_{s_R}]$, Eq. (@vlasov) with the periodic constraint becomes the Karush–Kuhn–Tucker system

$$
\begin{bmatrix}\mathbf H^{\mathsf T}\mathbf W\mathbf H & \mathbf C_u^{\mathsf T}\\ \mathbf C_u & \mathbf 0\end{bmatrix}
\begin{bmatrix}\mathbf u\\ \boldsymbol\mu\end{bmatrix} =
\begin{bmatrix}-\mathbf H^{\mathsf T}\mathbf W\mathbf g\\ \mathbf 0\end{bmatrix}. \eqn{kkt}
$$

which has a unique solution because $\mathbf c_{ax}$ removes the uniform warping, the null space of $\mathbf H^{\mathsf T}\mathbf W\mathbf H$. For a mode with non-zero circulation $\oint v\,ds$ around a closed contour, including the one closed through the side edges, the residual $\gamma_k$ is the membrane shear flow of a closed section, carried by $A_{66}\gamma_k\gamma_i$ in the modal tensors.

### Shear and Mindlin modes: merging of the side-edge functions

The shear and Mindlin modes are nodal; the periodic basis follows by merging the nodal functions of $s_L$ into those of $s_R$,

$$
h_{s_R}\leftarrow h_{s_R}+h_{s_L},\qquad
\boldsymbol\Omega_{s_R}\leftarrow\boldsymbol\Omega_{s_R}+\boldsymbol\Omega_{s_L},\qquad
\zeta_{s_R}\leftarrow\zeta_{s_R}+\zeta_{s_L},\qquad h_{s_L},\ \boldsymbol\Omega_{s_L},\ \zeta_{s_L}\ \text{removed}. \eqn{merge}
$$

which spans exactly the periodic subspace, since a combination $\mathbf S\mathbf a$ satisfies $a_{s_L} = a_{s_R}$ if, and only if, it lies in the span of the merged functions. The side nodes lie in the straight skin, so that $\boldsymbol\Omega_{s_L}$ and $\boldsymbol\Omega_{s_R}$ can be paired. This gives $n_S = N-1$ shear modes, including the axial mode, and $n_M = \sum_q m_q - m_{s_L} + N - 1$ Mindlin modes, with $m_q$ vectors $\boldsymbol\Omega_q$ at node $q$ (Table 2).

::: {custom-style="MDPI_4.1_table_caption"}
**Table 2.** Number of deformation modes per family of the periodic cell (independent of the lay-up).
:::

| $b_s$ [mm] | $N$ | G | D | L | T | S (incl. axial) | M | total |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 150 | 92 | 2 | 13 | 166 | 92 | 91 | 205 | 569 |
| 285.5 | 78 | 2 | 13 | 138 | 78 | 77 | 174 | 482 |

### Properties and hierarchical ordering

The periodic basis is (i) *exact*: every combination and truncation satisfies Eq. (@pbc) identically; (ii) *complete*: it spans the periodic subspace of the nodal description, so that the full basis is equivalent to imposing Eq. (@pbc) node by node; (iii) *mechanically interpretable* through the GBT families; and (iv) *consistent with Bloch–Floquet theory* {{mead1996}}: a cell of width $2b_s$ admits the Bloch waves $\mathbf u(X+b_s,x) = e^{\mathrm i\kappa b_s}\mathbf u(X,x)$ with $\kappa b_s = 0$ and $\kappa b_s = \pi$, i.e. modes that repeat in every bay or alternate from bay to bay. Imposed afterwards on a truncated classical set, the three rows of Eq. (@rows) would leave too few admissible modes; here the constraint is applied before any truncation, and a hierarchical ordering is introduced within the periodic subspace by the generalised eigenvalue problem $\mathbf B\,\boldsymbol\psi = \Lambda(k_r^2\mathbf D^{I} + k_r^4\mathbf C)\boldsymbol\psi$ with the tensors of Section 2.3 at a reference wavenumber $k_r$, small $\Lambda$ indicating long-wave modes. In the first eight modes (Figure 2), modes 1 and 2 are the rigid translations; the following ones combine the bending of the skin bays with the translation, rotation and distortion of the stringers and, unlike the modes of an open member, cannot be attributed to a single family.

![**Figure 2.** First eight hierarchical conventional modes of the periodic cell ($b_s = 150$ mm, Cross-ply lay-up, ordered at $k_r = \pi/b_s$). Line position: in-plane deformation (grey: undeformed); colour: warping displacement $u_k(s)$ normalised to its largest absolute value. The side nodes $s_L$ and $s_R$ (squares) have identical in-plane displacements and identical warping in every mode.](fig/fig02_modes.png){width=15.5cm}

## Member analysis

### Strain energy and modal tensors

The variation of the strain energy of the cell,

$$
\delta U = \int_L\int_b\int_t\boldsymbol\sigma^{\mathsf T}\delta\boldsymbol\varepsilon\,dz\,ds\,dx = \int_L\int_b\left(\mathbf N^{\mathsf T}\delta\boldsymbol\varepsilon^{M} + \mathbf M^{\mathsf T}\delta\boldsymbol\kappa + \mathbf Q^{\mathsf T}\delta\boldsymbol\gamma_z\right)ds\,dx, \eqn{dU0}
$$

becomes, after substitution of Eq. (@strain) and integration over the cross-section, in the form of Silvestre and Camotim {{silvestre2013shear}},

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

where the tensors with exchanged indices in the last two lines are those of the preceding two lines with the roles of conventional and shear modes exchanged, and $\delta U_M$ contains the Mindlin terms. The letters follow {{silvestre2013shear}}: $C$ pairs the axial strains, $B$ the transverse strains and curvatures, $D^{I}$ the shear strains and twisting curvatures, $D^{II}$ axial with transverse, $F$ transverse with shear and $H$ axial with shear quantities. For the conventional modes,

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

For modes that satisfy Vlasov's assumptions, Eqs. (@Cik)–(@couplingT) are exactly the tensors of Eqs. (10)–(12) of {{silvestre2013shear}}; the terms with $\gamma_k$ and $v_{k,s}$ are the additions for closed contours and transverse extension. The Mindlin terms follow in the same way and contain the transverse shear stiffness, e.g. $\int_b A_{s,11}\xi_r\xi_{r'}\,ds$. In the matrix notation of the GBT–Ritz approach {{kharghani2026jcs}}, the conventional tensors form

$$
\boldsymbol\Xi = \begin{bmatrix}\mathbf B & \mathbf F^{I} & \mathbf D^{II\mathsf T}\\ \mathbf F^{I\mathsf T} & \mathbf D^{I} & \mathbf H^{I\mathsf T}\\ \mathbf D^{II} & \mathbf H^{I} & \mathbf C\end{bmatrix},\qquad
\delta U = \int_L \delta\boldsymbol\Phi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Phi\,dx,\qquad
\boldsymbol\Phi = \begin{bmatrix}\boldsymbol\phi\\ \boldsymbol\phi_{,x}\\ \boldsymbol\phi_{,xx}\end{bmatrix}, \eqn{Xi}
$$

which is extended by the shear, coupling and Mindlin tensors for the complete vector $\boldsymbol\Phi$ of amplitude functions.

### Geometric tensors, loads and pre-buckling state

The second-order work of the pre-buckling membrane forces $S_{xx}$, $S_{ss}$, $S_{xs}$ on the strains of Eq. (@strainNL),

$$
V = \int_L\int_b\left(S_{xx}\,\varepsilon^{NL}_{xx} + S_{ss}\,\varepsilon^{NL}_{ss} + S_{xs}\,\gamma^{NL}_{xs}\right)ds\,dx, \eqn{V}
$$

defines the geometric tensors; for the conventional modes

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

where the dots stand for the analogous shear-mode terms; the Mindlin modes do not move the mid-surface and do not contribute. For a uniform $S_{xx}$ and $S_{ss} = S_{xs} = 0$, only $X^{U}_{ik}$ and $X^{VW}_{ik}$ remain, the geometric matrices of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}, $X^{VW}_{ik}$ being the tensor $X_{pik}$ of classical GBT {{silvestre2002b,silvestre2013shear}}. The cell is loaded at both ends by a compressive line load $n_0(s)$ acting at the distance $e_0(s)$ from the wall reference surface; with $u_0 = u - e_0 w_{,x}$ its potential gives

$$
\begin{gathered}
\delta\Pi_{ext} = -\left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=0} + \left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=L},\\
q^{I}_k = \int_b n_0\left(u_k - e_0 w_k\right)ds,\qquad q^{I}_h = \int_b n_0\,u_h\,ds,
\end{gathered} \eqn{load}
$$

i.e. the load enters through the generalised end forces of the warping, as $q^{I}_i$ in {{silvestre2013shear}}. Instead of the modal resultants $W^{o}_p$ of an assumed stress distribution {{silvestre2013shear}}, the pre-buckling state is the first-order solution of $\delta U + \delta\Pi_{ext} = 0$, with the membrane forces

$$
\left[S_{xx},\ S_{ss},\ S_{xs}\right]^{\mathsf T} = \mathbf A\,\boldsymbol\varepsilon^{M,0} + \mathbf B\,\boldsymbol\kappa^{0} \eqn{pre}
$$

which vary along $s$ and $x$: the load enters the stringer through the eccentric foot, is shared between skin and stringer by shear lag, and the restrained Poisson expansion produces transverse forces $S_{ss}$.

### Equilibrium equations and boundary conditions

The stationarity $\delta U + \lambda\,\delta V = 0$, with $\lambda$ the load factor, gives after integration by parts the GBT equilibrium equations

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

where the second equation is multiplied by $-1$ as in {{silvestre2013shear}}, and

$$
H_{ik} = H^{I}_{ik} - H^{I}_{ki},\qquad D_{ik} = D^{I}_{ik} - D^{II}_{ik} - D^{II}_{ki},\qquad F_{ik} = F^{I}_{ik} - F^{I}_{ki},\qquad H_{hj} = H^{I}_{hj} - H^{I}_{jh} . \eqn{combined}
$$

For a uniform $S_{xx}$, $\mathcal G_{ik} = X^{U}_{ik}\,d^4/dx^4 - X^{VW}_{ik}\,d^2/dx^2$; in general its coefficients depend on $x$. The boundary terms give

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

and the generalised shear stress resultant $W^{\tau}_i = -W^{\sigma}_{i,x} + W^{xs}_i$, $W^{xs}_i = H^{I}_{ki}\phi_{k,xx} + D^{I}_{ik}\phi_{k,x} + F^{I}_{ki}\phi_k + H^{II}_{ji}\varphi_{j,x} + D^{I}_{ji}\varphi_j$. Equations (@equil)–(@Wsig) generalise those of {{silvestre2013shear}} by the membrane terms, the Mindlin modes and the non-uniform pre-buckling state. Because of the $x$-dependent geometric operators and the end conditions of the stiffened cell, they are solved by the Ritz method.

### Ritz solution, end conditions and classical limit

As in {{kharghani2025jcs,kharghani2026jcs}}, $\boldsymbol\Phi = \boldsymbol\Psi\mathbf a$, with every amplitude function approximated as

$$
\phi_k(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{kb},\qquad \varphi_j(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{jb},\qquad \chi_r(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{rb} . \eqn{Psi}
$$

where, instead of the trigonometric series used for members {{kharghani2026jcs}}, $N_b$ are cubic B-splines on an open uniform knot vector with spans $h_x = b_s/20$, which resolve the end zones and act directly on the end coefficients ($\phi(0) = a_1$; $\phi_{,x}(0) = 0 \Leftrightarrow a_2 = a_1$). As in the shell models, the loaded edges are either **simply supported (SS)**, with $U_Y = 0$ at every end node, i.e. $\sum_k U_{Y,qk}\phi_k = 0$ at $x = 0, L$, imposed by a null-space reduction, while $U_X$, the warping and the rotations remain free (rigid-body motions are removed at $s_R$, $x = 0$), or **clamped (CC)**, with $\phi_k = \phi_{k,x} = 0$ and $\chi_r = 0$ at both ends and $\varphi_j(0) = 0$, so that $x = 0$ is fully clamped and $x = L$ is free to warp under the load. Unlike the classical condition $\phi_k = \sin(m\pi x/L)$, the SS condition suppresses only $U_Y$ of the end edges. Equations (@dU) and (@Xgeo) give the eigenvalue problem

$$
\left(\mathbf K + \lambda\,\mathbf G\right)\mathbf a = \mathbf 0,\qquad
\mathbf K = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Psi\,dx,\qquad
\mathbf G = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol{\mathcal X}(x)\,\boldsymbol\Psi\,dx, \eqn{eig}
$$

integrated with four Gauss points per span. The smallest positive $\lambda$ multiplies the reference load of 1 N/mm, so that $N_{cr} = \lambda\times1$ N/mm. Since $\mathbf G$ is indefinite, the largest $\mu = 1/\lambda$ of $-\mathbf G\mathbf a = \mu\mathbf K\mathbf a$ is computed by the implicitly restarted Lanczos method {{lehoucq1998}}. For $b_s = 150$ mm and $L = 225$ mm the model has $569\times33 = 18\,777$ coefficients.

For uncoupled walls, a uniform $S_{xx}$ and the conventional modes only, $H_{ik} = F_{ik} = 0$ and Eq. (@equil) becomes

$$
C_{ik}\,\phi_{k,xxxx} - D_{ik}\,\phi_{k,xx} + B_{ik}\,\phi_k + \lambda\left(X^{U}_{ik}\,\phi_{k,xxxx} - X^{VW}_{ik}\,\phi_{k,xx}\right) = 0, \eqn{classical}
$$

which, without $X^{U}_{ik}$, is the classical second-order GBT equation {{schardt1994,silvestre2002b}}; with $\phi_k = \hat\phi_k\sin(k_m x)$, $k_m = m\pi/L$,

$$
\left[B_{ik} + k_m^2\,D_{ik} + k_m^4\,C_{ik} + \lambda\left(k_m^2\,X^{VW}_{ik} + k_m^4\,X^{U}_{ik}\right)\right]\hat\phi_k = 0, \eqn{signature}
$$

whose smallest eigenvalue over $m$ is the signature curve. The periodic condition does not change this structure; it restricts the admissible modes.

## Mode decomposition

The buckling mode is decomposed into the contributions of the families $g\in\{G,D,L,T,S,M\}$ by an energy partition {{kharghani2026jcs}}: with $\mathbf a_g$ the vector $\mathbf a$ restricted to family $g$,

$$
U_g = \tfrac12\,\mathbf a_g^{\mathsf T}\mathbf K\,\mathbf a,\qquad \rho_g = \frac{U_g}{U},\qquad \sum_g\rho_g = 1 . \eqn{rho}
$$

Unlike the self-energies $\tfrac12\mathbf a_g^{\mathsf T}\mathbf K\mathbf a_g$ used for members {{kharghani2026jcs}}, which do not add up when families are coupled, Eq. (@rho) divides the strain energy exactly. This matters in the cell, where the D modes move the wall junctions and the L rotation modes restore wall continuity: for the Unsym cell with SS ends, $b_s = 150$ mm and $\beta = 1.5$, the self-energies would assign 41 % to D and 58 % to L, Eq. (@rho) 1 % and 97 %. Since $\mathbf K\mathbf a = -\lambda\mathbf G\mathbf a$, $\rho_g$ is also the share of family $g$ in the second-order work; the Mindlin modes, which do not enter $\mathbf G$ (Section 2.3.2), have $\rho_M = 0$ identically, and their influence is assessed by omitting them (Section 3.6). The same partition applies to individual hierarchical modes.
