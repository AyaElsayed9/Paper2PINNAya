
# Theory and formulations

GBT analyses a thin-walled member in two stages {{schardt1994,silvestre2002a,bebiano2015}}: the cross-section analysis (Sections 2.1 and 2.2) yields the deformation modes and their stiffness constants (modal tensors); the member analysis (Section 2.3) yields the amplitude functions of the modes along the length and the critical load. The notation is that of Silvestre and Camotim {{silvestre2013shear}}, and the member analysis follows the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}. Classical GBT is extended, and recovered as a special case (Section 2.3.4), by periodic modes, the main contribution (Section 2.2), closed contours without Vlasov's assumptions, Reissner–Mindlin laminated walls and a computed pre-buckling state.

## Cross-section analysis

### Periodic cell and displacement field

The cell of the infinitely wide panel (Figure 1b) contains two stringers at the pitch $b_s$; it has the width $W = 2b_s$ and the length $L = \beta b_s$, with the aspect ratio $\beta$ (Figure 1a). Two stringers are needed because a cell of one pitch contains only modes that repeat in every bay, not those that alternate from bay to bay (Section 2.2.5). $X$ runs across the panel, $Y$ normal to the skin and $x$ along the stringers; the skin occupies $-W \le X \le 0$, with the side edges $s_L$ ($X = -W$) and $s_R$ ($X = 0$) in the flat middle of a skin bay. Every wall has local axes $(x, s, z)$ ($s$ along its mid-line, $z$ normal) with mid-surface displacements $u$, $v$, $w$ {{silvestre2013shear}}. The walls are divided into segments of at most $b_s/25$ {{schardt1994,bebiano2015}}; for $b_s = 150$ mm the cell has $N = 92$ nodes, 19 walls and 93 segments (Figure 1c).

![**Figure 1.** Periodic cell of the omega-stringer-stiffened panel. (a) Cell of width $W = 2b_s$ and length $L = \beta b_s$ with the end load $n_0$ and the side edges $s_L$, $s_R$ linked by the periodic condition. (b) Infinitely wide panel and the cell (shaded). (c) GBT discretisation for $b_s = 150$ mm with nodes, wall types and local axes $(x, s, z)$.](fig/fig01_geometry.png){width=15.5cm}

In GBT, every displacement is a sum of mode profiles (functions of $s$) times amplitude functions (functions of $x$) {{silvestre2013shear}},

$$
u(x,s) = u_k(s)\,\phi_{k,x}(x) + u_j(s)\,\varphi_j(x),\qquad
v(x,s) = v_k(s)\,\phi_k(x),\qquad
w(x,s) = w_k(s)\,\phi_k(x), \eqn{disp}
$$

with summation over repeated indices. The conventional modes $k$ deform the cross-section in its plane ($v_k$, $w_k$), with a warping $u_k$ tied to $\phi_{k,x}$ as in beam theory; the shear modes $j$ add an independent warping $u_j$, i.e. mid-surface shear, and include the axial extension $u_j = 1$. The rotations of the Reissner–Mindlin wall normals are described by transverse shear (Mindlin) modes $r$,

$$
\psi_x(x,s) = \xi_r(s)\,\chi_r(x),\qquad \psi_s(x,s) = \zeta_r(s)\,\chi_r(x), \eqn{beta}
$$

with the profiles $\xi_r$, $\zeta_r$ and amplitudes $\chi_r$; in variations, $i$ and $h$ replace $k$ and $j$.

### Strains and constitutive relations

Inserting Eqs. (@disp) and (@beta) into the Reissner–Mindlin strains, i.e. membrane strains plus $z$ times the curvatures $\kappa_{xx} = -w_{,xx} + \psi_{x,x}$, $\kappa_{ss} = -w_{,ss} + \psi_{s,s}$, $\kappa_{xs} = -2w_{,xs} + \psi_{x,s} + \psi_{s,x}$, and the transverse shear strains $\gamma_{xz} = \psi_x$, $\gamma_{sz} = \psi_s$, gives

$$
\begin{gathered}
\varepsilon_{xx} = \left(u_k - z\,w_k\right)\phi_{k,xx} + u_j\,\varphi_{j,x} + z\,\xi_r\,\chi_{r,x},\qquad
\varepsilon_{ss} = \left(v_{k,s} - z\,w_{k,ss}\right)\phi_k + z\,\zeta_{r,s}\,\chi_r,\\
\gamma_{xs} = \left(\gamma_k - 2z\,w_{k,s}\right)\phi_{k,x} + u_{j,s}\,\varphi_j + z\left(\xi_{r,s}\,\chi_r + \zeta_r\,\chi_{r,x}\right),\qquad
\gamma_{xz} = \xi_r\,\chi_r,\qquad \gamma_{sz} = \zeta_r\,\chi_r,
\end{gathered} \eqn{strain}
$$

with the membrane shear strain profile $\gamma_k = u_{k,s} + v_k$. In {{silvestre2013shear}} the conventional modes satisfy Vlasov's assumptions $\gamma_k = 0$ (no membrane shear) and $v_{k,s} = 0$ (no transverse extension) and the walls are Kirchhoff plates ($\xi_r = \zeta_r = 0$); here, transverse extension modes are included and $\gamma_k = 0$ cannot hold in closed contours (Section 2.2.3). The geometric stiffness requires the quadratic membrane strains

$$
\varepsilon^{NL}_{xx} = \tfrac12\left(u_{,x}^2+v_{,x}^2+w_{,x}^2\right),\qquad
\varepsilon^{NL}_{ss} = \tfrac12\left(u_{,s}^2+v_{,s}^2+w_{,s}^2\right),\qquad
\gamma^{NL}_{xs} = u_{,x}u_{,s}+v_{,x}v_{,s}+w_{,x}w_{,s}, \eqn{strainNL}
$$

of which classical GBT retains $\tfrac12(v_{,x}^2+w_{,x}^2)$ {{silvestre2002b,silvestre2013shear}}. The membrane forces $\mathbf N$, moments $\mathbf M$ and transverse shear forces $\mathbf Q$ follow from first-order shear deformation laminate theory with a shear correction factor of 5/6 {{mindlin1951,whitney1970,reddy2004}},

$$
\begin{bmatrix}\mathbf N\\ \mathbf M\\ \mathbf Q\end{bmatrix} =
\begin{bmatrix} \mathbf A & \mathbf B & \mathbf 0\\ \mathbf B & \mathbf D & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \mathbf A_z\end{bmatrix}
\begin{bmatrix}\boldsymbol\varepsilon^{M}\\ \boldsymbol\kappa\\ \boldsymbol\gamma_z\end{bmatrix}, \eqn{claw}
$$

with the membrane, coupling, bending and transverse shear stiffnesses $\mathbf A$, $\mathbf B$, $\mathbf D$, $\mathbf A_z$. All couplings are retained {{silvestre2013shear,kharghani2026jcs}}, and $N_{ss}$ is not set to zero as in beam theories, since the cell restrains the transverse expansion of the skin (Section 3.2).

Two walls need special treatment. A stringer foot (skin plus bonded flange) is one wall whose stiffness is referred to the skin mid-plane, so that it joins the neighbouring skin at the same height. With the stiffnesses $\tilde{\mathbf A}$, $\tilde{\mathbf B}$, $\tilde{\mathbf D}$ of the combined stack about its own mid-plane, at $z_f = t_f/2$ from the skin mid-plane towards the flange, the standard shift of the reference plane ($z = z' + z_f$ in the stiffness integrals) gives

$$
\mathbf A_f = \tilde{\mathbf A},\qquad \mathbf B_f = \tilde{\mathbf B} + z_f\,\tilde{\mathbf A},\qquad \mathbf D_f = \tilde{\mathbf D} + 2z_f\,\tilde{\mathbf B} + z_f^2\,\tilde{\mathbf A} . \eqn{foot}
$$

The second wall reproduces the offset $c = (t_s+t_f)/2$ between the stringer and skin mid-surfaces (thicknesses $t_s$, $t_f$), bridged by a tie in the shell model. Every web root is a node at $Y = c$, joined to the skin node below by a link wall that must transfer displacements and rotation rigidly without adding material; its stiffnesses are multiples of the skin stiffnesses $\mathbf A_0$, $\mathbf D_0$, $\mathbf A_{z0}$,

$$
\begin{gathered}
A^{L}_{11} = 10^{-3}A_{0,11},\quad A^{L}_{22} = k_L A_{0,22},\quad A^{L}_{66} = k_L A_{0,66},\\
D^{L}_{11} = 10^{-3}D_{0,11},\quad D^{L}_{22} = k_L D_{0,22},\quad D^{L}_{66} = 10^{-3}D_{0,66},\quad \mathbf A^{L}_z = k_L\mathbf A_{z0},
\end{gathered} \eqn{link}
$$

with all other entries zero. The factor $k_L = 10^3$ makes the link rigid in the cross-section plane; the factor $10^{-3}$ removes its axial, longitudinal bending and twisting stiffness, which a row of rigid links does not have. Both are numerical choices: changing $k_L$ from $10^2$ to $10^4$ or the small factor from $10^{-2}$ to 0 changes the critical load by less than 0.01 % (Section 3.7).

### Nodal description and mode families

As in modern GBT {{bebiano2015,goncalves2014shear}}, the modes are built from nodal quantities: on every segment $p$ (width $b_p$, unit tangent $\mathbf t_p$, end nodes $q$, $q'$, local coordinate $\eta\in[0,1]$), $u$ and $v$ vary linearly, and $w$ is the cubic Hermite interpolation of the nodal normal displacements $w_q$ and rotations $\theta_q = w_{,s}$,

$$
\begin{gathered}
w(\eta) = \mathcal H_1(\eta)\,w_q + b_p \mathcal H_2(\eta)\,\theta_q + \mathcal H_3(\eta)\,w_{q'} + b_p \mathcal H_4(\eta)\,\theta_{q'},\\
\mathcal H_1 = 1-3\eta^2+2\eta^3,\quad \mathcal H_2 = \eta-2\eta^2+\eta^3,\quad \mathcal H_3 = 3\eta^2-2\eta^3,\quad \mathcal H_4 = -\eta^2+\eta^3 .
\end{gathered} \eqn{hermite}
$$

This is equivalent to Schardt's nodal moments {{schardt1994}}, but makes the periodic conditions linear in the nodal quantities. The in-plane deformation is described by the nodal displacements $\mathbf U$ (a vector of the $2N$ components $U_X$, $U_Y$ of all nodes) and the nodal rotations $\boldsymbol\Theta$ ($N$ components), and the segment elongations by the matrix $\mathbf E_\ell$, $(\mathbf E_\ell\mathbf U)_p = \mathbf t_p\cdot(\mathbf U_{q'}-\mathbf U_q)$. The modes are generated in the usual GBT families:

- **G** (global): rigid-body translations and rotation;
- **D** (distortional): nodal displacements without segment elongation ($\mathbf E_\ell\mathbf U = \mathbf 0$), orthogonal to G and L, which move the wall junctions;
- **L** (local): a unit normal displacement of a non-junction node or a unit nodal rotation, i.e. wall bending between junctions;
- **T** (transverse extension): the columns of $\mathbf E_\ell^{\mathsf T}$, which stretch the segments;
- **S** (shear): the nodal hat functions of the warping, one replaced by the axial mode $u_j = 1$;
- **M** (Mindlin): nodal vectors $\boldsymbol\Omega_q$ with $\xi_r = \mathbf t\cdot\boldsymbol\Omega_q$, and nodal values of $\zeta_r$, interpolated linearly.

G, D and L span the null space (kernel) $\ker\mathbf E_\ell = \{\mathbf U:\ \mathbf E_\ell\mathbf U = \mathbf 0\}$ of fields that elongate no segment, and T spans the range of $\mathbf E_\ell^{\mathsf T}$ (the space spanned by its columns); by the standard decomposition $\mathbb R^{2N} = \ker\mathbf E_\ell \oplus \operatorname{range}\mathbf E_\ell^{\mathsf T}$ of linear algebra, the families thus cover every in-plane deformation, with the rotations in family L. The warping $u_k$ of the G and D modes is not free: classical GBT chooses it such that the membrane shear strain $\gamma_k$ vanishes {{vlasov1961,schardt1994}}. In a closed contour this is impossible, so $u_k$ minimises the membrane shear strain energy,

$$
\min_{\mathbf u}\ \Pi_w(\mathbf u) = \sum_p \hat A_{66,p}\,b_p\left(\frac{u_{q'}-u_{q}}{b_p}+v_p\!\left(\eta = \tfrac12\right)\right)^2
\quad\text{subject to}\quad \mathbf c_{ax}^{\mathsf T}\mathbf u = 0, \eqn{vlasov}
$$

where $\mathbf u$ collects the nodal warping values, the bracket is $\gamma_k$ at the segment mid-point, $\hat A_{66} = A_{66}-A_{26}^2/A_{22}$ is the shear stiffness for zero transverse force, and the constraint, with $c_{ax,q} = \sum_{p\ni q}\tfrac12\hat A_{11,p}b_p$ and $\hat A_{11} = A_{11}-A_{12}^2/A_{22}$, removes the axial resultant (axial extension is the separate mode $u_j = 1$). For an open contour the minimum is zero and Vlasov's warping is recovered; in a closed contour the remaining $\gamma_k$ is kept in Eq. (@strain).

## Periodic boundary conditions in the cross-section analysis

### Continuous condition and its modal form

As one period of the infinite panel, the cell has side edges that are the same line of the panel, so all displacements and rotations must be equal on them,

$$
\mathcal U_\alpha(X=-W,\,x) - \mathcal U_\alpha(X=0,\,x) = 0,\qquad \alpha = 1,\dots,6,\qquad 0\le x\le L, \eqn{pbc}
$$

where $\mathcal U_1$–$\mathcal U_3$ are the translations in $X$, $Y$, $x$ and $\mathcal U_4$–$\mathcal U_6$ the rotations about $X$, $Y$, $x$ of a shell node; shell models impose Eq. (@pbc) with one equation per pair of side nodes and degree of freedom {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. Table 1 gives the corresponding conditions on the GBT modes.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 1.** Degrees of freedom of a skin side node and the corresponding modal conditions.
:::

| Shell DOF | Meaning at a skin side node | Modal condition |
|:--|:--|:--|
| $\mathcal U_1$ | in-plane displacement $U_X$ | $U_{X,s_L} = U_{X,s_R}$ for every conventional mode, Eq. (@rows) |
| $\mathcal U_2$ | deflection $U_Y$ | $U_{Y,s_L} = U_{Y,s_R}$ for every conventional mode, Eq. (@rows) |
| $\mathcal U_3$ | axial displacement $u$ | $u_k(s_L) = u_k(s_R)$, Eq. (@kkt); $u_j(s_L) = u_j(s_R)$, Eq. (@merge) |
| $\mathcal U_4$ | rotation $-w_{,x}+\psi_x$ | follows from $\mathcal U_2$ for $w_{,x}$; $\xi_r(s_L) = \xi_r(s_R)$, Eq. (@merge) |
| $\mathcal U_6$ | rotation $\theta+\psi_s$ | $\theta_{s_L} = \theta_{s_R}$, Eq. (@rows); $\zeta_r(s_L) = \zeta_r(s_R)$, Eq. (@merge) |
| $\mathcal U_5$ | drilling rotation | not a degree of freedom of the plate theory |

Inserting Eqs. (@disp) and (@beta) into Eq. (@pbc) gives, e.g. for $\mathcal U_3$,

$$
\left[u_k(s_L)-u_k(s_R)\right]\phi_{k,x}(x) + \left[u_j(s_L)-u_j(s_R)\right]\varphi_j(x) = 0\qquad \forall\,x\in[0,L]. \eqn{pbcmodal}
$$

**Proposition 1.** Equation (@pbc) holds for arbitrary amplitude functions if, and only if, every mode satisfies it individually: $\mathbf U_k$, $\theta_k$ and $u_k$ of all conventional modes, $u_j$ of all shear modes and $\xi_r$, $\zeta_r$ of all Mindlin modes take equal values at $s_L$ and $s_R$.

*Proof.* If every bracket vanishes, Eq. (@pbcmodal) and its analogues hold for any amplitudes. Conversely, if all amplitude functions but one are set to zero, the bracket of the remaining mode must vanish. $\square$

The periodic condition can therefore be built into the cross-section analysis by modes that are periodic one by one (Sections 2.2.2–2.2.4); every combination of them, in particular every Ritz trial function, is then periodic, without penalty terms or Lagrange multipliers.

### Conventional modes: constrained family subspaces

For the conventional modes, the conditions on $\mathcal U_1$, $\mathcal U_2$ and the nodal rotation part of $\mathcal U_6$ are three linear equations,

$$
\mathbf R\begin{bmatrix}\mathbf U\\ \boldsymbol\Theta\end{bmatrix} = \mathbf 0,\qquad
\mathbf R = \begin{bmatrix}\mathbf e^{\mathsf T}_{X,s_L}-\mathbf e^{\mathsf T}_{X,s_R}\\ \mathbf e^{\mathsf T}_{Y,s_L}-\mathbf e^{\mathsf T}_{Y,s_R}\\ \mathbf e^{\mathsf T}_{\theta,s_L}-\mathbf e^{\mathsf T}_{\theta,s_R}\end{bmatrix}\in\mathbb R^{3\times 3N},\qquad \operatorname{rank}\mathbf R = 3, \eqn{rows}
$$

where the unit vectors $\mathbf e$ pick the $U_X$, $U_Y$ or $\theta$ entry of a side node and the three rows are independent ($\operatorname{rank}\mathbf R = 3$). The periodic deformations form the null space $\ker\mathbf R$ of dimension $3N-3$. To keep the GBT meaning of the modes, they are built family by family. With $\mathbf V_g$ the nodal fields of family $g\in\{G,D,L,T\}$ of the non-periodic cell: (1) columns with $\mathbf R\mathbf v = \mathbf 0$ are kept; (2) the other columns $\mathbf V'_g$ are replaced by the largest periodic subspace $\mathbf V'_g\mathbf Z_g$, with $\mathbf Z_g$ a basis of $\ker(\mathbf R\mathbf V'_g)$, i.e. of all $\mathbf z$ with $\mathbf R\mathbf V'_g\mathbf z = \mathbf 0$; (3) periodic deformations not reached by steps 1–2, $(\mathbf I - \mathbf P)\ker\mathbf R$ with $\mathbf P$ the orthogonal projector onto their span, are added to T if they elongate a segment and to D otherwise; (4) every mode is scaled to a unit maximum. The columns of the modal matrix $\mathbf V$ then span exactly $\ker\mathbf R$: no admissible deformation is lost and none is spurious ($n_C = 3N-3$). Family G keeps only the two translations, since a rigid rotation $\omega$ about $x$ would lift one side edge and lower the other ($U_{Y,s_L}-U_{Y,s_R} = -W\omega$). Equal lateral translations and rotations of the side nodes enter family D in step 3.

### Warping: the periodic cell as a closed contour

The condition $u_k(s_L) = u_k(s_R)$ joins the side edges, so that for the warping the cell is a contour closed through its side edges, like the closed sections of GBT {{schardt1994,goncalves2013polygon}}. With the segment difference matrix $\boldsymbol\Delta$, $(\boldsymbol\Delta\mathbf u)_p = u_{q'} - u_q$, the weights $\mathbf A_w = \mathrm{diag}(\hat A_{66,p}/b_p)$ and the vector $\mathbf f$ of $b_p v_p(\eta = 1/2)$, the objective of Eq. (@vlasov) is $(\boldsymbol\Delta\mathbf u + \mathbf f)^{\mathsf T}\mathbf A_w(\boldsymbol\Delta\mathbf u + \mathbf f)$; with the constraints $\mathbf C_u\mathbf u = \mathbf 0$, $\mathbf C_u = [\mathbf c_{ax}^{\mathsf T};\ \mathbf e^{\mathsf T}_{s_L}-\mathbf e^{\mathsf T}_{s_R}]$, its minimum is the solution of the standard Karush–Kuhn–Tucker system

$$
\begin{bmatrix}\boldsymbol\Delta^{\mathsf T}\mathbf A_w\boldsymbol\Delta & \mathbf C_u^{\mathsf T}\\ \mathbf C_u & \mathbf 0\end{bmatrix}
\begin{bmatrix}\mathbf u\\ \boldsymbol\tau\end{bmatrix} =
\begin{bmatrix}-\boldsymbol\Delta^{\mathsf T}\mathbf A_w\mathbf f\\ \mathbf 0\end{bmatrix}, \eqn{kkt}
$$

with the Lagrange multipliers $\boldsymbol\tau$. The solution is unique, since the axial constraint excludes uniform warping, the only warping without shear energy. For a mode with non-zero circulation $\oint v\,ds$, the remaining $\gamma_k$ is the shear flow of a closed section, carried by $A_{66}\gamma_k\gamma_i$ in the modal tensors.

### Shear and Mindlin modes: merging of the side-edge functions

The shear and Mindlin modes are defined node by node (hat function $\hat u_q$, vectors $\boldsymbol\Omega_q$ and function $\zeta_q$ of node $q$). A combination with the coefficients $y_q$ is periodic if, and only if, $y_{s_L} = y_{s_R}$, so the periodic basis replaces the separate functions of the two side nodes by their sums,

$$
\bar u = \hat u_{s_L}+\hat u_{s_R},\qquad
\bar{\boldsymbol\Omega} = \boldsymbol\Omega_{s_L}+\boldsymbol\Omega_{s_R},\qquad
\bar\zeta = \zeta_{s_L}+\zeta_{s_R}, \eqn{merge}
$$

while the functions of all other nodes are kept. This basis spans exactly the periodic subspace ($\boldsymbol\Omega_{s_L}$ and $\boldsymbol\Omega_{s_R}$ can be paired because both side nodes lie in the straight skin). This gives $n_S = N-1$ shear modes, including the axial mode, and $n_M = \sum_q m_q - m_{s_L} + N - 1$ Mindlin modes, with $m_q$ independent vectors $\boldsymbol\Omega_q$ at node $q$ ($m_q = 2$ at a kink, $m_q = 1$ between collinear segments); see Table 2.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 2.** Number of deformation modes per family of the periodic cell (independent of the lay-up).
:::

| $b_s$ [mm] | $N$ | G | D | L | T | S (incl. axial) | M | total |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 150 | 92 | 2 | 13 | 166 | 92 | 91 | 194 | 558 |
| 285.5 | 78 | 2 | 13 | 138 | 78 | 77 | 166 | 474 |

### Properties and hierarchical ordering

The periodic basis is exact (every truncation satisfies Eq. (@pbc)), complete (the full basis is equivalent to imposing Eq. (@pbc) node by node), interpretable through the GBT families and consistent with Bloch–Floquet theory {{mead1996}}: a cell of width $2b_s$ contains the Bloch waves $\mathbf u(X+b_s,x) = e^{\mathrm i k_B b_s}\mathbf u(X,x)$ with the Bloch wavenumbers $k_B b_s = 0$ and $\pi$, i.e. modes that repeat in every bay and modes that alternate. To order them from global to local, the conventional modes are recombined by $\hat{\mathbf B}\,\mathbf y = \Lambda(k_r^2\hat{\mathbf D}^{I} + k_r^4\hat{\mathbf C})\,\mathbf y$ (matrices of the modal tensors, Eq. (@Xi); reference wavenumber $k_r$); a small ratio $\Lambda$ of transverse to longitudinal stiffness indicates a global mode, and the recombined modes remain periodic. In Figure 2, modes 1 and 2 are rigid translations; the others combine skin-bay bending with stringer translation, rotation and distortion and, unlike modes of open members, belong to no single family.

![**Figure 2.** First eight hierarchical conventional modes of the periodic cell ($b_s = 150$ mm, Cross-ply lay-up, ordered at $k_r = \pi/b_s$). Line: in-plane deformation (grey: undeformed); colour: normalised warping $u_k(s)$. The side nodes $s_L$ and $s_R$ (squares) move identically in every mode.](fig/fig02_modes.png){width=15.5cm}

## Member analysis

### Strain energy and modal tensors

To obtain the equations for the amplitude functions, the strain energy is expressed in terms of the modes. Its variation

$$
\delta U = \int_L\int_b\int_t\boldsymbol\sigma^{\mathsf T}\delta\boldsymbol\varepsilon\,dz\,ds\,dx = \int_L\int_b\left(\mathbf N^{\mathsf T}\delta\boldsymbol\varepsilon^{M} + \mathbf M^{\mathsf T}\delta\boldsymbol\kappa + \mathbf Q^{\mathsf T}\delta\boldsymbol\gamma_z\right)ds\,dx, \eqn{dU0}
$$

(integration over the length $L$, contour $b$ and thickness $t$) becomes, with Eqs. (@strain) and (@claw) and integration over $s$, in the form of Silvestre and Camotim {{silvestre2013shear}},

$$
\begin{gathered}
\delta U = \int_L\Big(C_{ik}\phi_{k,xx}\delta\phi_{i,xx} + H^{I}_{ik}\phi_{k,x}\delta\phi_{i,xx} + H^{I}_{ki}\phi_{k,xx}\delta\phi_{i,x} + D^{I}_{ik}\phi_{k,x}\delta\phi_{i,x}\\
+ D^{II}_{ik}\phi_{k}\delta\phi_{i,xx} + D^{II}_{ki}\phi_{k,xx}\delta\phi_{i} + F^{I}_{ik}\phi_{k,x}\delta\phi_{i} + F^{I}_{ki}\phi_{k}\delta\phi_{i,x} + B_{ik}\phi_{k}\delta\phi_{i}\\
+ C_{hj}\varphi_{j,x}\delta\varphi_{h,x} + H^{I}_{hj}\varphi_{j}\delta\varphi_{h,x} + H^{I}_{jh}\varphi_{j,x}\delta\varphi_{h} + D^{I}_{hj}\varphi_{j}\delta\varphi_{h}\\
+ C_{ij}\varphi_{j,x}\delta\phi_{i,xx} + H^{I}_{ij}\varphi_{j}\delta\phi_{i,xx} + H^{II}_{ji}\varphi_{j,x}\delta\phi_{i,x}\\
+ D^{I}_{ji}\varphi_{j}\delta\phi_{i,x} + D^{II}_{ji}\varphi_{j,x}\delta\phi_{i} + F^{II}_{ij}\varphi_{j}\delta\phi_{i}\\
+ C_{hk}\phi_{k,xx}\delta\varphi_{h,x} + H^{II}_{kh}\phi_{k,xx}\delta\varphi_{h} + H^{I}_{hk}\phi_{k,x}\delta\varphi_{h,x}\\
+ D^{I}_{hk}\phi_{k,x}\delta\varphi_{h} + D^{II}_{hk}\phi_{k}\delta\varphi_{h,x} + F^{II}_{kh}\phi_{k}\delta\varphi_{h}\Big)dx + \delta U_M,
\end{gathered} \eqn{dU}
$$

where the first two lines contain the conventional modes, the third the shear modes, the others their coupling, and $\delta U_M$ the Mindlin terms. The tensors of the last two lines follow from those of the preceding two by exchanging conventional and shear modes ($C_{ij}\to C_{hk}$, $H^{I}_{ij}\to H^{II}_{kh}$, $H^{II}_{ji}\to H^{I}_{hk}$, $D^{I}_{ji}\to D^{I}_{hk}$, $D^{II}_{ji}\to D^{II}_{hk}$, $F^{II}_{ij}\to F^{II}_{kh}$, as in Eq. (12) of {{silvestre2013shear}}). $C$ pairs axial, $B$ transverse, $D^{I}$ shear, $D^{II}$ axial–transverse, $F$ transverse–shear and $H$ axial–shear strains. For the conventional modes,

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

For modes satisfying Vlasov's assumptions, Eqs. (@Cik)–(@couplingT) reduce to Eqs. (10)–(12) of {{silvestre2013shear}}; the terms with $\gamma_k$ and $v_{k,s}$ are the additions for closed contours and transverse extension. In the block-matrix notation of the GBT–Ritz approach {{kharghani2026jcs}}, the conventional tensors form the matrices $\hat{\mathbf B} = [B_{ik}]$, $\hat{\mathbf C} = [C_{ik}]$, etc. (the hat distinguishes them from the laminate stiffnesses $\mathbf A$, $\mathbf B$, $\mathbf D$ of Eq. (@claw)) and the block matrix

$$
\boldsymbol\Xi = \begin{bmatrix}\hat{\mathbf B} & \hat{\mathbf F}^{I} & (\hat{\mathbf D}^{II})^{\mathsf T}\\ (\hat{\mathbf F}^{I})^{\mathsf T} & \hat{\mathbf D}^{I} & (\hat{\mathbf H}^{I})^{\mathsf T}\\ \hat{\mathbf D}^{II} & \hat{\mathbf H}^{I} & \hat{\mathbf C}\end{bmatrix},\qquad
\delta U = \int_L \delta\boldsymbol\Phi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Phi\,dx,\qquad
\boldsymbol\Phi = \begin{bmatrix}\boldsymbol\phi\\ \boldsymbol\phi_{,x}\\ \boldsymbol\phi_{,xx}\end{bmatrix}, \eqn{Xi}
$$

where $\boldsymbol\phi$ collects the amplitudes $\phi_k$ (the blocks $\hat{\mathbf D}_2$ and $\hat{\mathbf H}_1$ of {{kharghani2026jcs}} are $(\hat{\mathbf D}^{II})^{\mathsf T}$ and $(\hat{\mathbf H}^{I})^{\mathsf T}$ here). For all modes, including $\delta U_M$, writing the strains of Eq. (@claw) as $[\boldsymbol\varepsilon^{M};\boldsymbol\kappa;\boldsymbol\gamma_z] = \sum_{d=0}^{2}\boldsymbol{\mathcal E}^{(d)}(s)\,\boldsymbol\Phi^{(d)}(x)$, where $\boldsymbol\Phi^{(d)}$ collects the $d$-th derivatives of all amplitude functions ($\boldsymbol\Phi = [\boldsymbol\Phi^{(0)};\boldsymbol\Phi^{(1)};\boldsymbol\Phi^{(2)}]$) and $\boldsymbol{\mathcal E}^{(d)}$ contains the corresponding strain profiles of Eq. (@strain), turns Eq. (@dU0) into

$$
\delta U = \int_L\sum_{d,e=0}^{2}\delta\boldsymbol\Phi^{(d)\mathsf T}\,\boldsymbol\Xi^{(de)}\,\boldsymbol\Phi^{(e)}\,dx,\qquad
\boldsymbol\Xi^{(de)} = \int_b\boldsymbol{\mathcal E}^{(d)\mathsf T}\,\mathbf C_{\mathrm{lam}}\,\boldsymbol{\mathcal E}^{(e)}\,ds, \eqn{Xide}
$$

with $\mathbf C_{\mathrm{lam}}$ the laminate stiffness matrix of Eq. (@claw). For a Mindlin mode $r$, Eq. (@strain) gives the non-zero entries $\zeta_{r,s}$ ($\kappa_{ss}$), $\xi_{r,s}$ ($\kappa_{xs}$), $\xi_r$ ($\gamma_{xz}$), $\zeta_r$ ($\gamma_{sz}$) of $\boldsymbol{\mathcal E}^{(0)}$ and $\xi_r$ ($\kappa_{xx}$), $\zeta_r$ ($\kappa_{xs}$) of $\boldsymbol{\mathcal E}^{(1)}$; for the conventional and shear modes, Eq. (@Xide) reproduces Eqs. (@Cik)–(@couplingT).

### Geometric tensors, loads and pre-buckling state

Buckling occurs when the work of the pre-buckling membrane forces $S_{xx}$, $S_{ss}$, $S_{xs}$ (per unit length, tension positive) on the non-linear strains of Eq. (@strainNL) balances the strain energy of the mode. This second-order work

$$
V = \int_L\int_b\left(S_{xx}\,\varepsilon^{NL}_{xx} + S_{ss}\,\varepsilon^{NL}_{ss} + S_{xs}\,\gamma^{NL}_{xs}\right)ds\,dx, \eqn{V}
$$

defines the geometric tensors; with Eq. (@disp), for the conventional modes,

$$
\begin{gathered}
\delta V = \int_L\Big(X^{U}_{ik}\,\phi_{k,xx}\delta\phi_{i,xx} + X^{VW}_{ik}\,\phi_{k,x}\delta\phi_{i,x} + X^{S}_{ik}\,\phi_{k}\delta\phi_{i}\\
+ X^{P}_{ik}\,\phi_{k,x}\delta\phi_{i} + X^{P}_{ki}\,\phi_{k}\delta\phi_{i,x} + X^{R}_{ik}\,\phi_{k,xx}\delta\phi_{i,x} + X^{R}_{ki}\,\phi_{k,x}\delta\phi_{i,xx}\Big)dx + \dots,\\
X^{U}_{ik} = \int_b S_{xx}\,u_i u_k\,ds,\qquad
X^{VW}_{ik} = \int_b\left[S_{xx}\left(v_i v_k + w_i w_k\right) + S_{ss}\,u_{i,s}u_{k,s}\right]ds,\\
X^{S}_{ik} = \int_b S_{ss}\left(v_{i,s}v_{k,s} + w_{i,s}w_{k,s}\right)ds,\qquad
X^{P}_{ik} = \int_b S_{xs}\left(v_{i,s}v_k + w_{i,s}w_k\right)ds,\\
X^{R}_{ik} = \int_b S_{xs}\,u_{i,s}u_k\,ds,
\end{gathered} \eqn{Xgeo}
$$

where the dots stand for the shear-mode terms. In the general form analogous to Eq. (@Xide), with the displacement gradients of all modes $\partial[u,v,w]^{\mathsf T}/\partial x = \sum_d\mathbf J^{(d)}_x\boldsymbol\Phi^{(d)}$ and $\partial[u,v,w]^{\mathsf T}/\partial s = \sum_d\mathbf J^{(d)}_s\boldsymbol\Phi^{(d)}$,

$$
\boldsymbol{\mathcal X}^{(de)} = \int_b\Big[S_{xx}\,\mathbf J_x^{(d)\mathsf T}\mathbf J_x^{(e)} + S_{ss}\,\mathbf J_s^{(d)\mathsf T}\mathbf J_s^{(e)} + S_{xs}\left(\mathbf J_x^{(d)\mathsf T}\mathbf J_s^{(e)} + \mathbf J_s^{(d)\mathsf T}\mathbf J_x^{(e)}\right)\Big]ds, \eqn{Xgen}
$$

which contains Eq. (@Xgeo) and, with $\partial u/\partial x = u_j\varphi_{j,x}$ and $\partial u/\partial s = u_{j,s}\varphi_j$, the shear-mode terms; the Mindlin modes do not move the mid-surface and do not contribute. For a uniform $S_{xx}$ and $S_{ss} = S_{xs} = 0$, only $X^{U}_{ik}$ and $X^{VW}_{ik}$ remain, the geometric matrices of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}; $X^{VW}_{ik}$ corresponds to $X_{pik}W^{o}_{p}$ of {{silvestre2002b,silvestre2013shear}} (there, $p$ is a mode index).

The compressive end load $n_0(s)$, acting at the distance $e_0(s)$ from the wall reference surface, enters through generalised end forces, as $q^{I}_i$ in {{silvestre2013shear}}; with the axial displacement $u_0 = u - e_0 w_{,x}$ at the load,

$$
\begin{gathered}
\delta\Pi_{ext} = -\left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=0} + \left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=L},\\
q^{I}_k = \int_b n_0\left(u_k - e_0 w_k\right)ds,\qquad q^{I}_h = \int_b n_0\,u_h\,ds.
\end{gathered} \eqn{load}
$$

Instead of assuming stresses through modal resultants $W^{o}_p$ {{silvestre2013shear}}, the pre-buckling state is computed, since the load enters through the eccentric foot, is shared by shear lag and produces transverse forces (Section 3.2). It is the linear static solution of the same GBT–Ritz model under the end load, $\delta U + \delta\Pi_{ext} = 0$, which gives the membrane strains $\boldsymbol\varepsilon^{M,0}$, the curvatures $\boldsymbol\kappa^{0}$ and the forces

$$
\left[S_{xx},\ S_{ss},\ S_{xs}\right]^{\mathsf T} = \mathbf A\,\boldsymbol\varepsilon^{M,0} + \mathbf B\,\boldsymbol\kappa^{0}, \eqn{pre}
$$

which vary along $s$ and $x$.

### Equilibrium equations and boundary conditions

At the critical load, $\delta U + \lambda\,\delta V = 0$ for all variations, with the load factor $\lambda$ scaling the pre-buckling forces. Integration by parts in $x$ gives one GBT equilibrium equation per conventional and per shear mode (the Mindlin equations follow likewise),

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
\mathcal K_{hk} = C_{hk}\frac{d^3}{dx^3} + \left(H^{I}_{hk} - H^{II}_{kh}\right)\frac{d^2}{dx^2} - \left(D^{I}_{hk} - D^{II}_{hk}\right)\frac{d}{dx} - F^{II}_{kh},\\
\mathcal K_{hj} = C_{hj}\frac{d^2}{dx^2} + H_{hj}\frac{d}{dx} - D^{I}_{hj},
\end{gathered} \eqn{Kop}
$$

where the shear-mode equation is multiplied by $-1$, as $\mathcal K_{hj}$ in {{silvestre2013shear}}, and

$$
H_{ik} = H^{I}_{ik} - H^{I}_{ki},\qquad D_{ik} = D^{I}_{ik} - D^{II}_{ik} - D^{II}_{ki},\qquad F_{ik} = F^{I}_{ik} - F^{I}_{ki},\qquad H_{hj} = H^{I}_{hj} - H^{I}_{jh} . \eqn{combined}
$$

The geometric operators $\mathcal G$ follow from Eq. (@Xgeo); for a uniform $S_{xx}$, $\mathcal G_{ik} = X^{U}_{ik}\,d^4/dx^4 - X^{VW}_{ik}\,d^2/dx^2$, in general with $x$-dependent coefficients. The boundary terms give the end conditions

$$
\left[\left(W^{\sigma}_i + q^{I}_i\right)\delta\phi_{i,x} + W^{\tau}_i\,\delta\phi_i + \left(W^{\sigma}_h + q^{I}_h\right)\delta\varphi_h\right]_{x=0}^{x=L} = 0, \eqn{bc}
$$

with the generalised normal stress resultants

$$
\begin{gathered}
W^{\sigma}_i = C_{ik}\phi_{k,xx} + H^{I}_{ik}\phi_{k,x} + D^{II}_{ik}\phi_{k} + C_{ij}\varphi_{j,x} + H^{I}_{ij}\varphi_{j},\\
W^{\sigma}_h = C_{hj}\varphi_{j,x} + H^{I}_{hj}\varphi_{j} + C_{hk}\phi_{k,xx} + H^{I}_{hk}\phi_{k,x} + D^{II}_{hk}\phi_k,
\end{gathered} \eqn{Wsig}
$$

and the generalised shear stress resultant $W^{\tau}_i = -W^{\sigma}_{i,x} + W^{xs}_i$, $W^{xs}_i = H^{I}_{ki}\phi_{k,xx} + D^{I}_{ik}\phi_{k,x} + F^{I}_{ki}\phi_k + H^{II}_{ji}\varphi_{j,x} + D^{I}_{ji}\varphi_j$. Since the coefficients vary with $x$ and the end conditions couple all modes, these equations are solved by the Ritz method, directly from the energy.

### Ritz solution, end conditions and classical limit

As in {{kharghani2025jcs,kharghani2026jcs}}, every amplitude function is a sum of basis functions with unknown coefficients,

$$
\phi_k(x) = \sum_{l=1}^{n_b} \Psi_l(x)\,a_{kl},\qquad \varphi_j(x) = \sum_{l=1}^{n_b} \Psi_l(x)\,a_{jl},\qquad \chi_r(x) = \sum_{l=1}^{n_b} \Psi_l(x)\,a_{rl}, \eqn{Psi}
$$

i.e. $\boldsymbol\Phi = \boldsymbol\Psi\mathbf a$, with $\boldsymbol\Psi$ the matrix of the basis functions $\Psi_l$ and their derivatives and $\mathbf a$ all coefficients. Instead of trigonometric series {{kharghani2025jcs,kharghani2026jcs}}, the $n_b$ functions $\Psi_l$ are cubic B-splines on an open uniform knot vector with spans $h_x = b_s/20$: they are local, so that they resolve the end zones, and their end values are the end coefficients ($\phi(0) = a_1$; $\phi_{,x}(0) = 0 \Leftrightarrow a_2 = a_1$), so that end conditions are imposed exactly. Cubic B-splines are established in spline finite strips {{lau1986}} and isogeometric analysis {{hughes2005iga}}; unlike the Hermite elements of GBT-based finite elements {{bebiano2018}}, they are twice continuously differentiable.

As in the shell models (Section 3), **simply supported (SS)** loaded edges have $U_Y = 0$ at every end node, i.e. $\sum_k U_{Y,qk}\phi_k = 0$ at $x = 0, L$, imposed by reducing $\mathbf a$ to the null space of these conditions; unlike the classical $\phi_k = \sin(m\pi x/L)$, this leaves $U_X$, warping and rotations free. **Clamped (CC)** edges have $\phi_k = \phi_{k,x} = 0$ and $\chi_r = 0$ at both ends and $\varphi_j(0) = 0$, so that $x = L$ can warp under the load. Inserting $\boldsymbol\Phi = \boldsymbol\Psi\mathbf a$ into Eqs. (@Xide) and (@Xgen) turns $\delta U + \lambda\,\delta V = 0$ into

$$
\left(\mathbf K + \lambda\,\mathbf G\right)\mathbf a = \mathbf 0,\qquad
\mathbf K = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Psi\,dx,\qquad
\mathbf G = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol{\mathcal X}(x)\,\boldsymbol\Psi\,dx, \eqn{eig}
$$

with $\boldsymbol\Xi = [\boldsymbol\Xi^{(de)}]$ and $\boldsymbol{\mathcal X} = [\boldsymbol{\mathcal X}^{(de)}]$ (four Gauss points per span). With pre-buckling forces for $n_0 = 1$ N/mm, the smallest positive $\lambda$ is $N_{cr}$ in N/mm. Since $\mathbf G$ is indefinite, the largest $\mu = 1/\lambda$ of $-\mathbf G\mathbf a = \mu\mathbf K\mathbf a$ is computed by the implicitly restarted Lanczos method {{lehoucq1998}}. For $b_s = 150$ mm and $L = 225$ mm (33 B-splines), the SS model has 18,230 unknowns.

Classical GBT is a special case, used for verification (Section 3.1): for uncoupled walls ($A_{16} = A_{26} = 0$, $\mathbf B = \mathbf 0$, $D_{16} = D_{26} = 0$), a uniform $S_{xx}$ and conventional modes only, $H_{ik} = F_{ik} = 0$ and Eq. (@equil) becomes

$$
C_{ik}\,\phi_{k,xxxx} - D_{ik}\,\phi_{k,xx} + B_{ik}\,\phi_k + \lambda\left(X^{U}_{ik}\,\phi_{k,xxxx} - X^{VW}_{ik}\,\phi_{k,xx}\right) = 0, \eqn{classical}
$$

which, without $X^{U}_{ik}$, is the classical second-order GBT equation {{schardt1994,silvestre2002b}}; with $\phi_k = \hat\phi_k\sin(k_m x)$, $k_m = m\pi/L$, for $m$ half-waves,

$$
\left[B_{ik} + k_m^2\,D_{ik} + k_m^4\,C_{ik} + \lambda\left(k_m^2\,X^{VW}_{ik} + k_m^4\,X^{U}_{ik}\right)\right]\hat\phi_k = 0, \eqn{signature}
$$

whose smallest eigenvalue over $m$ is the signature curve. The periodic condition does not change this structure; it restricts the admissible modes.

## Mode decomposition

As in {{kharghani2026jcs}}, the buckling mode is decomposed into the contributions of the families $g\in\{G,D,L,T,S,M\}$: with $\mathbf a_g$ the vector $\mathbf a$ with all coefficients outside family $g$ set to zero, the strain energy $U = \tfrac12\mathbf a^{\mathsf T}\mathbf K\mathbf a$ is divided into

$$
U_g = \tfrac12\,\mathbf a_g^{\mathsf T}\mathbf K\,\mathbf a,\qquad \rho_g = \frac{U_g}{U},\qquad \sum_g\rho_g = 1 . \eqn{rho}
$$

The participation $\rho_g$ includes half of the coupling energy between $g$ and every other family, so that the participations add up exactly to one; the normalised self-energies $\tfrac12\mathbf a_g^{\mathsf T}\mathbf K\mathbf a_g$ of {{kharghani2026jcs}} omit the coupling energy. This matters for the strongly coupled D and L modes: for the Unsym cell ($b_s = 150$ mm, $\beta = 1.5$, SS) the normalised self-energies assign 41 % to D and 58 % to L, Eq. (@rho) 1 % and 97 %. Since $\mathbf K\mathbf a = -\lambda\mathbf G\mathbf a$ and the Mindlin modes do not enter $\mathbf G$, $\rho_M = 0$; their influence is assessed by omitting them (Section 3.6). The partition also applies to individual hierarchical modes.
