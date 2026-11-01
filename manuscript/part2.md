
# Theory and formulations

GBT analyses a thin-walled member in two stages {{schardt1994,silvestre2002a,bebiano2015}}: the cross-section analysis (Sections 2.1 and 2.2) determines the deformation modes and their stiffness constants, the modal tensors; the member analysis (Section 2.3) determines the amplitude functions, i.e. how much of each mode occurs along the length, and the critical load. The notation follows the shear-deformable GBT of Silvestre and Camotim {{silvestre2013shear}}, and the member analysis follows the GBT–Ritz approach of Kharghani and Mittelstedt {{kharghani2025jcs,kharghani2026jcs}}. The formulation extends GBT in four respects, each of which reduces to classical GBT in a special case (Section 2.3.4): (i) all modes satisfy the periodic boundary conditions, the main contribution of this paper (Section 2.2); (ii) Vlasov's assumptions are dropped, because the cell contains closed contours; (iii) the walls are Reissner–Mindlin laminates; and (iv) the pre-buckling state is computed.

## Cross-section analysis

### Periodic cell and displacement field

The infinitely wide panel carries identical omega stringers at the pitch $b_s$ (Figure 1b). The analysed cell contains two stringers; it has the width $W = 2b_s$ and the length $L = \beta b_s$, with the aspect ratio $\beta$ (Figure 1a). Two stringers are needed because a cell of one pitch can only represent buckling modes that repeat in every bay, whereas a cell of two pitches also contains modes that alternate from bay to bay (Section 2.2.5). $X$ runs across the panel, $Y$ normal to the skin and $x$ along the stringers. The skin occupies $-W \le X \le 0$; its side edges $s_L$ ($X = -W$) and $s_R$ ($X = 0$) lie in the middle of a skin bay, i.e. in the flat skin. Every wall has the local axes $(x, s, z)$, with $s$ along the wall mid-line and $z$ normal to it, and the mid-surface displacements $u$ (along $x$), $v$ (along $s$) and $w$ (along $z$) {{silvestre2013shear}}. The walls (skin, feet, webs, crowns) are divided into segments of at most $b_s/25$ {{schardt1994,bebiano2015}}; for $b_s = 150$ mm the cell has $N = 92$ nodes, 19 walls and 93 segments (Figure 1c).

![**Figure 1.** Periodic cell of the omega-stringer-stiffened panel. (a) Cell of width $W = 2b_s$ and length $L = \beta b_s$ with the end load $n_0$ and the side edges $s_L$, $s_R$ linked by the periodic condition. (b) Infinitely wide panel and the cell (shaded). (c) GBT discretisation for $b_s = 150$ mm with nodes, wall types and local axes $(x, s, z)$.](fig/fig01_geometry.png){width=15.5cm}

In GBT, every displacement is a sum of mode profiles (functions of $s$) times amplitude functions (functions of $x$); in the shear-deformable GBT {{silvestre2013shear}},

$$
u(x,s) = u_k(s)\,\phi_{k,x}(x) + u_j(s)\,\varphi_j(x),\qquad
v(x,s) = v_k(s)\,\phi_k(x),\qquad
w(x,s) = w_k(s)\,\phi_k(x), \eqn{disp}
$$

with summation over repeated indices. The conventional modes $k$ deform the cross-section in its plane ($v_k$, $w_k$), with a warping $u_k$ tied to the slope $\phi_{k,x}$ as in beam theory; the shear modes $j$ add a warping $u_j$ with an independent amplitude $\varphi_j$, i.e. shear deformation of the wall mid-surfaces, and include the uniform axial extension $u_j = 1$. In Reissner–Mindlin walls, the rotations of the wall normals are independent of the slopes of $w$; they are described by transverse shear (Mindlin) modes $r$,

$$
\psi_x(x,s) = \xi_r(s)\,\chi_r(x),\qquad \psi_s(x,s) = \zeta_r(s)\,\chi_r(x), \eqn{beta}
$$

with the rotation profiles $\xi_r$, $\zeta_r$ and the amplitudes $\chi_r$. In variations, the indices $i$ (conventional) and $h$ (shear) replace $k$ and $j$.

### Strains and constitutive relations

The strains of a Reissner–Mindlin plate are the membrane strains plus $z$ times the curvatures $\kappa_{xx} = -w_{,xx} + \psi_{x,x}$, $\kappa_{ss} = -w_{,ss} + \psi_{s,s}$, $\kappa_{xs} = -2w_{,xs} + \psi_{x,s} + \psi_{s,x}$, together with the transverse shear strains $\gamma_{xz} = \psi_x$, $\gamma_{sz} = \psi_s$. Inserting Eqs. (@disp) and (@beta) gives the strains in terms of the modes,

$$
\begin{gathered}
\varepsilon_{xx} = \left(u_k - z\,w_k\right)\phi_{k,xx} + u_j\,\varphi_{j,x} + z\,\xi_r\,\chi_{r,x},\qquad
\varepsilon_{ss} = \left(v_{k,s} - z\,w_{k,ss}\right)\phi_k + z\,\zeta_{r,s}\,\chi_r,\\
\gamma_{xs} = \left(\gamma_k - 2z\,w_{k,s}\right)\phi_{k,x} + u_{j,s}\,\varphi_j + z\left(\xi_{r,s}\,\chi_r + \zeta_r\,\chi_{r,x}\right),\qquad
\gamma_{xz} = \xi_r\,\chi_r,\qquad \gamma_{sz} = \zeta_r\,\chi_r,
\end{gathered} \eqn{strain}
$$

with the membrane shear strain profile $\gamma_k = u_{k,s} + v_k$ of a conventional mode. In {{silvestre2013shear}} the conventional modes satisfy Vlasov's assumptions of zero membrane shear strain, $\gamma_k = 0$, and zero transverse membrane strain, $v_{k,s} = 0$, and the walls are Kirchhoff plates ($\xi_r = \zeta_r = 0$); Eq. (@strain) then reduces to their kinematics. Neither assumption is imposed here: modes with transverse extension are part of the basis, and the cell contains closed contours, in which $\gamma_k = 0$ cannot be satisfied (Section 2.2.3). For the geometric stiffness, the non-linear (quadratic) parts of the membrane strains are needed,

$$
\varepsilon^{NL}_{xx} = \tfrac12\left(u_{,x}^2+v_{,x}^2+w_{,x}^2\right),\qquad
\varepsilon^{NL}_{ss} = \tfrac12\left(u_{,s}^2+v_{,s}^2+w_{,s}^2\right),\qquad
\gamma^{NL}_{xs} = u_{,x}u_{,s}+v_{,x}v_{,s}+w_{,x}w_{,s}, \eqn{strainNL}
$$

of which classical GBT retains only $\tfrac12(v_{,x}^2+w_{,x}^2)$ {{silvestre2002b,silvestre2013shear}}. The membrane forces $\mathbf N$, moments $\mathbf M$ and transverse shear forces $\mathbf Q$ of a wall follow from the first-order shear deformation theory of laminates with a shear correction factor of 5/6 {{mindlin1951,whitney1970,reddy2004}},

$$
\begin{bmatrix}\mathbf N\\ \mathbf M\\ \mathbf Q\end{bmatrix} =
\begin{bmatrix} \mathbf A & \mathbf B & \mathbf 0\\ \mathbf B & \mathbf D & \mathbf 0\\ \mathbf 0 & \mathbf 0 & \mathbf A_s\end{bmatrix}
\begin{bmatrix}\boldsymbol\varepsilon^{M}\\ \boldsymbol\kappa\\ \boldsymbol\gamma_z\end{bmatrix}, \eqn{claw}
$$

with the membrane strains $\boldsymbol\varepsilon^{M}$, curvatures $\boldsymbol\kappa$ and transverse shear strains $\boldsymbol\gamma_z$, and the membrane, coupling, bending and transverse shear stiffnesses $\mathbf A$, $\mathbf B$, $\mathbf D$, $\mathbf A_s$ of the laminate. All couplings are retained {{silvestre2013shear,kharghani2026jcs}}, and the transverse membrane force $N_{ss}$ is not set to zero as in beam theories, because the periodic cell restrains the transverse expansion of the skin (Section 3.2).

Two walls need special treatment. A stringer foot (skin plus bonded flange) is one wall whose stiffness must be referred to the skin mid-plane, so that it joins the neighbouring skin at the same height. With the stiffnesses $\tilde{\mathbf A}$, $\tilde{\mathbf B}$, $\tilde{\mathbf D}$ of the combined stack about its own mid-plane, which lies at $z_f = t_f/2$ from the skin mid-plane towards the flange, the standard shift of the reference plane of a laminate (obtained by inserting $z = z' + z_f$ into the stiffness integrals) gives

$$
\mathbf A_f = \tilde{\mathbf A},\qquad \mathbf B_f = \tilde{\mathbf B} + z_f\,\tilde{\mathbf A},\qquad \mathbf D_f = \tilde{\mathbf D} + 2z_f\,\tilde{\mathbf B} + z_f^2\,\tilde{\mathbf A} . \eqn{foot}
$$

The second reproduces the offset of the stringer: since the flange lies on the skin, the stringer mid-surface is offset by $c = (t_s+t_f)/2$ (skin and stringer thicknesses $t_s$, $t_f$), which the shell model bridges by a tie constraint. Here every web root is a node at $Y = c$, joined to the skin node below by a link wall of height $c$ that must act like the tie, i.e. transfer displacements and rotation rigidly without adding material. Its stiffnesses are multiples of the skin stiffnesses $\mathbf A_0$, $\mathbf D_0$, $\mathbf A_{s0}$,

$$
\begin{gathered}
A^{L}_{11} = 10^{-3}A_{0,11},\quad A^{L}_{22} = k_L A_{0,22},\quad A^{L}_{66} = k_L A_{0,66},\\
D^{L}_{11} = 10^{-3}D_{0,11},\quad D^{L}_{22} = k_L D_{0,22},\quad D^{L}_{66} = 10^{-3}D_{0,66},\quad \mathbf A^{L}_s = k_L\mathbf A_{s0},
\end{gathered} \eqn{link}
$$

with all other entries zero. The entries multiplied by the large factor $k_L = 10^3$ make the link rigid in the cross-section plane (no stretching, sliding, bending or transverse shear of the link); the entries multiplied by $10^{-3}$ remove its axial, longitudinal bending and twisting stiffness, which a row of rigid links does not have. Both factors are numerical choices; the critical loads do not depend on them (Section 3.7).

### Nodal description and mode families

As in modern GBT {{bebiano2015,goncalves2014shear}}, the modes are built from nodal quantities. On every segment $p$ of width $b_p$, unit tangent $\mathbf t_p$ and end nodes $q$, $q'$ (local coordinate $\eta\in[0,1]$), $u$ and $v$ vary linearly, and $w$ is the cubic Hermite interpolation of the nodal normal displacements $w_q$ and rotations $\theta_q = w_{,s}$,

$$
\begin{gathered}
w(\eta) = H_1(\eta)\,w_q + b_p H_2(\eta)\,\theta_q + H_3(\eta)\,w_{q'} + b_p H_4(\eta)\,\theta_{q'},\\
H_1 = 1-3\eta^2+2\eta^3,\quad H_2 = \eta-2\eta^2+\eta^3,\quad H_3 = 3\eta^2-2\eta^3,\quad H_4 = -\eta^2+\eta^3 .
\end{gathered} \eqn{hermite}
$$

This is equivalent to Schardt's nodal moments {{schardt1994}}, but makes the periodic conditions linear in the nodal quantities (Section 2.2). The in-plane deformation of the cross-section is described by the nodal displacements $\mathbf U\in\mathbb R^{2N}$ (components $U_X$, $U_Y$ of every node) and the nodal rotations $\boldsymbol\Theta\in\mathbb R^N$ (shared by the walls meeting at a node). The segment elongations follow from the elongation matrix $\mathbf E_\ell$, $(\mathbf E_\ell\mathbf U)_p = \mathbf t_p\cdot(\mathbf U_{q'}-\mathbf U_q)$. The modes are generated in the usual GBT families:

- **G** (global): rigid-body translations and rotation of the cross-section;
- **D** (distortional): nodal displacements without segment elongation, $\mathbf E_\ell\mathbf U = \mathbf 0$, orthogonal to G and L, which change the shape of the cross-section by moving the wall junctions;
- **L** (local): a unit normal displacement of a node that is not a wall junction, or a unit rotation of a node, i.e. bending of the walls between junctions;
- **T** (transverse extension): the columns of $\mathbf E_\ell^{\mathsf T}$, which stretch the segments in the $s$ direction;
- **S** (shear): the nodal hat functions of the warping, one of them replaced by the axial mode $u_j = 1$;
- **M** (Mindlin): nodal vectors $\boldsymbol\Omega_q$ in the cross-section plane, which define the rotation profile $\xi_r = \mathbf t\cdot\boldsymbol\Omega_q$, and nodal values of $\zeta_r$, both interpolated linearly.

Together, G, D, L and T describe every in-plane deformation of the discretised cross-section. Any nodal displacement field $\mathbf U$ splits uniquely into a part that elongates no segment, i.e. lies in the null space (kernel) $\ker\mathbf E_\ell = \{\mathbf U:\ \mathbf E_\ell\mathbf U = \mathbf 0\}$, which is spanned by the displacement fields of G, D and L, and a part that elongates segments, which lies in the range of $\mathbf E_\ell^{\mathsf T}$ (the space spanned by its columns), i.e. in family T. In short, $\mathbb R^{2N} = \ker\mathbf E_\ell \oplus \operatorname{range}\mathbf E_\ell^{\mathsf T}$, a standard result of linear algebra: the null space of a matrix and the range of its transpose are orthogonal complements. The nodal rotations $\boldsymbol\Theta$ are described by the rotation modes of family L. The warping $u_k$ of the G and D modes is not free: in classical GBT it is chosen such that the membrane shear strain $\gamma_k$ vanishes {{vlasov1961,schardt1994}}. In a closed contour this is impossible, so $u_k$ is chosen to make the membrane shear strain energy as small as possible,

$$
\min_{\mathbf u}\ \Pi_w(\mathbf u) = \sum_p \hat A_{66,p}\,b_p\left(\frac{u_{q'}-u_{q}}{b_p}+v_p\!\left(\eta = \tfrac12\right)\right)^2
\quad\text{subject to}\quad \mathbf c_{ax}^{\mathsf T}\mathbf u = 0, \eqn{vlasov}
$$

where $\mathbf u$ collects the nodal warping values, the bracket is the shear strain $\gamma_k$ of segment $p$ at its mid-point, $\hat A_{66} = A_{66}-A_{26}^2/A_{22}$ is the shear stiffness for zero transverse membrane force, and the constraint with $c_{ax,q} = \sum_{p\ni q}\tfrac12\hat A_{11,p}b_p$, $\hat A_{11} = A_{11}-A_{12}^2/A_{22}$, removes the axial resultant of the warping (the uniform axial extension is the separate mode $u_j = 1$). For an open contour the minimum is zero and the classical Vlasov warping is recovered; in a closed contour the remaining $\gamma_k$ is kept in Eq. (@strain).

## Periodic boundary conditions in the cross-section analysis

### Continuous condition and its modal form

The cell is one period of the infinite panel. Its two side edges are therefore the same line of the panel, and every displacement and rotation must be equal on them for every $x$,

$$
\mathcal U_\alpha(X=-W,\,x) - \mathcal U_\alpha(X=0,\,x) = 0,\qquad \alpha = 1,\dots,6,\qquad 0\le x\le L, \eqn{pbc}
$$

where $\mathcal U_1$, $\mathcal U_2$, $\mathcal U_3$ are the translations in $X$, $Y$, $x$ and $\mathcal U_4$, $\mathcal U_5$, $\mathcal U_6$ the rotations about $X$, $Y$, $x$ of a shell node. Eq. (@pbc) is the constraint set that shell models impose with one equation per pair of side nodes and degree of freedom {{mittelstedt2007,mittelstedt2009periodic,schilling2020rineng}}. Table 1 lists what each of these degrees of freedom (DOFs) means in GBT and which condition it imposes on the modes.

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

Inserting the GBT displacements, Eqs. (@disp) and (@beta), into Eq. (@pbc) gives, for example for $\mathcal U_3$,

$$
\left[u_k(s_L)-u_k(s_R)\right]\phi_{k,x}(x) + \left[u_j(s_L)-u_j(s_R)\right]\varphi_j(x) = 0\qquad \forall\,x\in[0,L]. \eqn{pbcmodal}
$$

**Proposition.** Equation (@pbc) holds for arbitrary amplitude functions if, and only if, every mode satisfies it individually: $\mathbf U_k$, $\theta_k$ and $u_k$ of all conventional modes, $u_j$ of all shear modes and $\xi_r$, $\zeta_r$ of all Mindlin modes take equal values at $s_L$ and $s_R$.

*Proof.* If every bracket vanishes, Eq. (@pbcmodal) and its analogues hold for any amplitudes. Conversely, if all amplitude functions but one are set to zero, the bracket of the remaining mode must vanish. $\square$

The periodic condition can therefore be built into the cross-section analysis by constructing modes that are periodic one by one (Sections 2.2.2–2.2.4); every combination of them, in particular every Ritz trial function (Section 2.3.4), is then periodic, and the member analysis needs neither penalty terms nor Lagrange multipliers.

### Conventional modes: constrained family subspaces

For the conventional modes, the conditions on $\mathcal U_1$, $\mathcal U_2$ and on the nodal rotation part of $\mathcal U_6$ are three linear equations for the nodal vector,

$$
\mathbf R\begin{bmatrix}\mathbf U\\ \boldsymbol\Theta\end{bmatrix} = \mathbf 0,\qquad
\mathbf R = \begin{bmatrix}\mathbf e^{\mathsf T}_{X,s_L}-\mathbf e^{\mathsf T}_{X,s_R}\\ \mathbf e^{\mathsf T}_{Y,s_L}-\mathbf e^{\mathsf T}_{Y,s_R}\\ \mathbf e^{\mathsf T}_{\theta,s_L}-\mathbf e^{\mathsf T}_{\theta,s_R}\end{bmatrix}\in\mathbb R^{3\times 3N},\qquad \operatorname{rank}\mathbf R = 3, \eqn{rows}
$$

where the unit vectors $\mathbf e$ pick the $U_X$, $U_Y$ or $\theta$ entry of a side node; the three rows are independent ($\operatorname{rank}\mathbf R = 3$). The admissible (periodic) in-plane deformations are the solutions of Eq. (@rows), i.e. the null space $\ker\mathbf R$ of dimension $3N-3$. The periodic modes are built family by family, so that each mode keeps its GBT meaning. With $\mathbf V_g$ the nodal fields of family $g\in\{G,D,L,T\}$ of the cell without periodic condition: (1) columns that already satisfy $\mathbf R\mathbf v = \mathbf 0$ are kept; (2) the remaining columns $\mathbf V'_g$ are replaced by the combinations $\mathbf V'_g\mathbf Z_g$, with $\mathbf Z_g$ a basis of $\ker(\mathbf R\mathbf V'_g)$, i.e. of all coefficient vectors $\mathbf z$ with $\mathbf R\mathbf V'_g\mathbf z = \mathbf 0$; these combinations form the largest periodic subspace of family $g$; (3) periodic deformations not reached by steps 1–2, $(\mathbf I - \boldsymbol\Pi_B)\ker\mathbf R$ with $\boldsymbol\Pi_B$ the orthogonal projector onto the span of steps 1–2, are added to T if they elongate a segment and to D otherwise; (4) every mode is scaled to a unit maximum. The resulting modal matrix $\mathbf V$ satisfies $\mathbf R\mathbf V = \mathbf 0$, and its columns span exactly $\ker\mathbf R$: no admissible deformation is lost and none is spurious, and there are $n_C = 3N-3$ conventional modes. Family G keeps only the two translations: a rigid rotation $\omega$ of the cell about $x$ would lift one side edge and lower the other ($U_{Y,s_L}-U_{Y,s_R} = -W\omega$), which is impossible in an infinite panel. Equal lateral translations and rotations of the two side nodes are admissible and enter family D in step 3.

### Warping: the periodic cell as a closed contour

The condition $u_k(s_L) = u_k(s_R)$ joins the two side edges of the skin. For the warping, the cell therefore behaves like a contour closed through its side edges, similar to the closed cross-sections of GBT {{schardt1994,goncalves2013polygon}}. The warping of the G and D modes follows from the least-squares problem of Eq. (@vlasov) with this additional constraint. With the segment difference matrix $\boldsymbol\Delta$, $(\boldsymbol\Delta\mathbf u)_p = u_{q'} - u_q$, the diagonal weight matrix $\mathbf A_w = \mathrm{diag}(\hat A_{66,p}/b_p)$, the vector $\mathbf f$ with the entries $b_p v_p(\eta = 1/2)$ and the constraint matrix $\mathbf C_u = [\mathbf c_{ax}^{\mathsf T};\ \mathbf e^{\mathsf T}_{s_L}-\mathbf e^{\mathsf T}_{s_R}]$, the objective of Eq. (@vlasov) is $(\boldsymbol\Delta\mathbf u + \mathbf f)^{\mathsf T}\mathbf A_w(\boldsymbol\Delta\mathbf u + \mathbf f)$, and its minimum under the two constraints is the solution of the standard Karush–Kuhn–Tucker system

$$
\begin{bmatrix}\boldsymbol\Delta^{\mathsf T}\mathbf A_w\boldsymbol\Delta & \mathbf C_u^{\mathsf T}\\ \mathbf C_u & \mathbf 0\end{bmatrix}
\begin{bmatrix}\mathbf u\\ \boldsymbol\tau\end{bmatrix} =
\begin{bmatrix}-\boldsymbol\Delta^{\mathsf T}\mathbf A_w\mathbf f\\ \mathbf 0\end{bmatrix}, \eqn{kkt}
$$

with the Lagrange multipliers $\boldsymbol\tau$. The solution is unique, because the only warping without shear strain energy, a uniform $\mathbf u$, is excluded by the axial constraint. For a mode whose tangential displacements circulate around a closed contour, $\oint v\,ds \neq 0$, the remaining $\gamma_k$ is the shear flow of a closed section; its stiffness enters the modal tensors through the term $A_{66}\gamma_k\gamma_i$.

### Shear and Mindlin modes: merging of the side-edge functions

The shear and Mindlin modes are defined node by node: a warping hat function $\hat u_q$, rotation vectors $\boldsymbol\Omega_q$ and a function $\zeta_q$ belong to each node $q$. A combination of nodal functions with the coefficients $a_q$ is periodic if, and only if, $a_{s_L} = a_{s_R}$. The periodic basis therefore follows by adding the function of $s_L$ to that of $s_R$ and removing the former,

$$
\hat u_{s_R}\leftarrow\hat u_{s_R}+\hat u_{s_L},\qquad
\boldsymbol\Omega_{s_R}\leftarrow\boldsymbol\Omega_{s_R}+\boldsymbol\Omega_{s_L},\qquad
\zeta_{s_R}\leftarrow\zeta_{s_R}+\zeta_{s_L},\qquad \hat u_{s_L},\ \boldsymbol\Omega_{s_L},\ \zeta_{s_L}\ \text{removed}, \eqn{merge}
$$

which spans exactly the periodic subspace. The vectors $\boldsymbol\Omega_{s_L}$ and $\boldsymbol\Omega_{s_R}$ can be paired because both side nodes lie in the straight skin. This gives $n_S = N-1$ shear modes, including the axial mode, and $n_M = \sum_q m_q - m_{s_L} + N - 1$ Mindlin modes, with $m_q$ independent vectors $\boldsymbol\Omega_q$ at node $q$ ($m_q = 2$ at a kink of the contour, $m_q = 1$ between collinear segments). Table 2 summarises the number of modes.

::: {custom-style="MDPI_4.1_table_caption"}
**Table 2.** Number of deformation modes per family of the periodic cell (independent of the lay-up).
:::

| $b_s$ [mm] | $N$ | G | D | L | T | S (incl. axial) | M | total |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| 150 | 92 | 2 | 13 | 166 | 92 | 91 | 194 | 558 |
| 285.5 | 78 | 2 | 13 | 138 | 78 | 77 | 166 | 474 |

### Properties and hierarchical ordering

The periodic basis is (i) *exact*: every truncation satisfies Eq. (@pbc) identically; (ii) *complete*: the full basis is equivalent to imposing Eq. (@pbc) node by node; (iii) interpretable through the GBT families; and (iv) consistent with Bloch–Floquet theory of periodic structures {{mead1996}}: a cell of width $2b_s$ contains the Bloch waves $\mathbf u(X+b_s,x) = e^{\mathrm i\kappa_B b_s}\mathbf u(X,x)$ with the Bloch wavenumbers $\kappa_B b_s = 0$ and $\pi$, i.e. modes that repeat in every bay and modes that alternate from bay to bay. The modes of Section 2.2.2 are nodal-type fields without a physical order. To obtain modes ordered from global to local, as in GBT, the conventional modes are recombined by the generalised eigenvalue problem $\mathbf B\,\mathbf y = \Lambda(k_r^2\mathbf D^{I} + k_r^4\mathbf C)\,\mathbf y$, with the modal tensors $\mathbf B$, $\mathbf D^{I}$, $\mathbf C$ of Section 2.3 and a reference wavenumber $k_r$. $\Lambda$ compares the transverse (in-plane) stiffness of a mode with its longitudinal (axial and shear) stiffness at the wavelength $\pi/k_r$; a small $\Lambda$ indicates a global, long-wave mode. Since this recombination acts within the periodic subspace, every ordered mode remains periodic. In Figure 2, modes 1 and 2 are the rigid translations; the following modes combine skin-bay bending with stringer translation, rotation and distortion and, unlike the modes of an open member, belong to no single family.

![**Figure 2.** First eight hierarchical conventional modes of the periodic cell ($b_s = 150$ mm, Cross-ply lay-up, ordered at $k_r = \pi/b_s$). Line: in-plane deformation (grey: undeformed); colour: normalised warping $u_k(s)$. The side nodes $s_L$ and $s_R$ (squares) move identically in every mode.](fig/fig02_modes.png){width=15.5cm}

## Member analysis

### Strain energy and modal tensors

To obtain the equations for the amplitude functions, the strain energy of the cell is expressed in terms of the modes. Its variation,

$$
\delta U = \int_L\int_b\int_t\boldsymbol\sigma^{\mathsf T}\delta\boldsymbol\varepsilon\,dz\,ds\,dx = \int_L\int_b\left(\mathbf N^{\mathsf T}\delta\boldsymbol\varepsilon^{M} + \mathbf M^{\mathsf T}\delta\boldsymbol\kappa + \mathbf Q^{\mathsf T}\delta\boldsymbol\gamma_z\right)ds\,dx, \eqn{dU0}
$$

with integration over the length $L$, the cross-section contour $b$ and the wall thickness $t$, becomes, after inserting Eqs. (@strain) and (@claw) and integrating over $s$, a sum of products of amplitude derivatives and cross-section constants. In the form of Silvestre and Camotim {{silvestre2013shear}},

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

where the first two lines contain the conventional modes, the third line the shear modes, the remaining lines their coupling, and $\delta U_M$ the Mindlin terms. The tensors of the last two lines follow from those of the preceding two by exchanging the roles of conventional and shear modes ($C_{ij}\to C_{hk}$, $H^{I}_{ij}\to H^{II}_{kh}$, $H^{II}_{ji}\to H^{I}_{hk}$, $D^{I}_{ji}\to D^{I}_{hk}$, $D^{II}_{ji}\to D^{II}_{hk}$, $F^{II}_{ij}\to F^{II}_{kh}$, as in Eq. (12) of {{silvestre2013shear}}). The letters indicate which strains are paired: $C$ axial with axial, $B$ transverse with transverse, $D^{I}$ shear with shear, $D^{II}$ axial with transverse, $F$ transverse with shear and $H$ axial with shear. The tensors follow directly from the integration over $s$; for the conventional modes,

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

For modes that satisfy Vlasov's assumptions ($\gamma_k = 0$, $v_{k,s} = 0$), Eqs. (@Cik)–(@couplingT) reduce to Eqs. (10)–(12) of {{silvestre2013shear}}; the terms with $\gamma_k$ and $v_{k,s}$ are the additions for closed contours and transverse extension. In the block-matrix notation of the GBT–Ritz approach {{kharghani2026jcs}}, the conventional tensors form

$$
\boldsymbol\Xi = \begin{bmatrix}\mathbf B & \mathbf F^{I} & (\mathbf D^{II})^{\mathsf T}\\ (\mathbf F^{I})^{\mathsf T} & \mathbf D^{I} & (\mathbf H^{I})^{\mathsf T}\\ \mathbf D^{II} & \mathbf H^{I} & \mathbf C\end{bmatrix},\qquad
\delta U = \int_L \delta\boldsymbol\Phi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Phi\,dx,\qquad
\boldsymbol\Phi = \begin{bmatrix}\boldsymbol\phi\\ \boldsymbol\phi_{,x}\\ \boldsymbol\phi_{,xx}\end{bmatrix}, \eqn{Xi}
$$

where $\boldsymbol\phi$ collects the amplitude functions $\phi_k$; the blocks $\hat{\mathbf D}_2$ and $\hat{\mathbf H}_1$ of {{kharghani2026jcs}} are $(\mathbf D^{II})^{\mathsf T}$ and $(\mathbf H^{I})^{\mathsf T}$ in the index convention used here. The same structure holds for all modes, including the Mindlin terms $\delta U_M$. The generalised strains of Eq. (@claw) can be written as $[\boldsymbol\varepsilon^{M};\boldsymbol\kappa;\boldsymbol\gamma_z] = \sum_{d=0}^{2}\boldsymbol{\mathcal E}^{(d)}(s)\,\boldsymbol\Phi^{(d)}(x)$, where $\boldsymbol\Phi^{(d)}$ collects the $d$-th derivatives of all amplitude functions, so that $\boldsymbol\Phi = [\boldsymbol\Phi^{(0)};\boldsymbol\Phi^{(1)};\boldsymbol\Phi^{(2)}]$ as in Eq. (@Xi), and the columns of $\boldsymbol{\mathcal E}^{(d)}$ are the strain profiles of Eq. (@strain) that multiply these derivatives. Inserting this into Eq. (@dU0) gives

$$
\delta U = \int_L\sum_{d,e=0}^{2}\delta\boldsymbol\Phi^{(d)\mathsf T}\,\boldsymbol\Xi^{(de)}\,\boldsymbol\Phi^{(e)}\,dx,\qquad
\boldsymbol\Xi^{(de)} = \int_b\boldsymbol{\mathcal E}^{(d)\mathsf T}\,\hat{\mathbf C}\,\boldsymbol{\mathcal E}^{(e)}\,ds, \eqn{Xide}
$$

with $\hat{\mathbf C}$ the stiffness matrix of Eq. (@claw). For a Mindlin mode $r$, the non-zero entries of $\boldsymbol{\mathcal E}^{(0)}$ are $\zeta_{r,s}$ ($\kappa_{ss}$), $\xi_{r,s}$ ($\kappa_{xs}$), $\xi_r$ ($\gamma_{xz}$) and $\zeta_r$ ($\gamma_{sz}$), and those of $\boldsymbol{\mathcal E}^{(1)}$ are $\xi_r$ ($\kappa_{xx}$) and $\zeta_r$ ($\kappa_{xs}$), as read from Eq. (@strain); for the conventional and shear modes, Eq. (@Xide) reproduces Eqs. (@Cik)–(@couplingT).

### Geometric tensors, loads and pre-buckling state

The membrane forces present before buckling, the pre-buckling forces $S_{xx}$, $S_{ss}$, $S_{xs}$ (per unit length, tension positive), do work on the non-linear strains of Eq. (@strainNL) when the cell buckles; buckling occurs at the load at which this work balances the strain energy of the buckling mode. The second-order work,

$$
V = \int_L\int_b\left(S_{xx}\,\varepsilon^{NL}_{xx} + S_{ss}\,\varepsilon^{NL}_{ss} + S_{xs}\,\gamma^{NL}_{xs}\right)ds\,dx, \eqn{V}
$$

defines the geometric tensors. Inserting Eq. (@disp) into Eq. (@strainNL) gives, for the conventional modes,

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

where the dots stand for the corresponding shear-mode terms. All terms, including these, are contained in a general form analogous to Eq. (@Xide): with the displacement gradients of all modes written as $\partial[u,v,w]^{\mathsf T}/\partial\alpha = \sum_d\mathbf J^{(d)}_\alpha\boldsymbol\Phi^{(d)}$, $\alpha\in\{x,s\}$,

$$
\boldsymbol{\mathcal X}^{(de)} = \int_b\Big[S_{xx}\,\mathbf J_x^{(d)\mathsf T}\mathbf J_x^{(e)} + S_{ss}\,\mathbf J_s^{(d)\mathsf T}\mathbf J_s^{(e)} + S_{xs}\left(\mathbf J_x^{(d)\mathsf T}\mathbf J_s^{(e)} + \mathbf J_s^{(d)\mathsf T}\mathbf J_x^{(e)}\right)\Big]ds, \eqn{Xgen}
$$

which reproduces Eq. (@Xgeo) and, with $\partial u/\partial x = u_j\varphi_{j,x}$ and $\partial u/\partial s = u_{j,s}\varphi_j$, gives the shear-mode terms. The Mindlin modes rotate the wall normals but do not move the mid-surface, so they do not contribute to $V$. For a uniform $S_{xx}$ and $S_{ss} = S_{xs} = 0$, only $X^{U}_{ik}$ and $X^{VW}_{ik}$ remain; these are the geometric matrices of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}, and $X^{VW}_{ik}$ corresponds to $X_{pik}W^{o}_{p}$ of {{silvestre2002b,silvestre2013shear}}.

The cell is loaded at both ends by a compressive line load $n_0(s)$, acting at the distance $e_0(s)$ from the wall reference surface. Its potential, with the axial displacement $u_0 = u - e_0 w_{,x}$ at the point of load application, enters through generalised end forces, as $q^{I}_i$ in {{silvestre2013shear}},

$$
\begin{gathered}
\delta\Pi_{ext} = -\left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=0} + \left(q^{I}_k\,\delta\phi_{k,x} + q^{I}_h\,\delta\varphi_{h}\right)_{x=L},\\
q^{I}_k = \int_b n_0\left(u_k - e_0 w_k\right)ds,\qquad q^{I}_h = \int_b n_0\,u_h\,ds.
\end{gathered} \eqn{load}
$$

Instead of the modal resultants $W^{o}_p$ of an assumed stress distribution {{silvestre2013shear}}, the pre-buckling state is computed, because in the stiffened cell the load enters the stringer through the eccentric foot, is shared by shear lag and produces transverse forces $S_{ss}$ (Section 3.2). It is the linear static solution of the same GBT–Ritz model under the end load, $\delta U + \delta\Pi_{ext} = 0$, with the membrane strains $\boldsymbol\varepsilon^{M,0}$, curvatures $\boldsymbol\kappa^{0}$ and, by Eq. (@claw), the membrane forces

$$
\left[S_{xx},\ S_{ss},\ S_{xs}\right]^{\mathsf T} = \mathbf A\,\boldsymbol\varepsilon^{M,0} + \mathbf B\,\boldsymbol\kappa^{0}, \eqn{pre}
$$

which vary along $s$ and $x$.

### Equilibrium equations and boundary conditions

At the critical load, the cell is in neutral equilibrium: $\delta U + \lambda\,\delta V = 0$ for all variations, with the load factor $\lambda$ that scales the pre-buckling forces. Integration by parts in $x$ gives the GBT equilibrium equations, one for each conventional mode $i$ and one for each shear mode $h$ (the Mindlin equations follow in the same way),

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

where the shear-mode equation is multiplied by $-1$, which gives $\mathcal K_{hj}$ the form of {{silvestre2013shear}}, and

$$
H_{ik} = H^{I}_{ik} - H^{I}_{ki},\qquad D_{ik} = D^{I}_{ik} - D^{II}_{ik} - D^{II}_{ki},\qquad F_{ik} = F^{I}_{ik} - F^{I}_{ki},\qquad H_{hj} = H^{I}_{hj} - H^{I}_{jh} . \eqn{combined}
$$

The geometric operators $\mathcal G$ follow in the same way from Eq. (@Xgeo); for a uniform $S_{xx}$, $\mathcal G_{ik} = X^{U}_{ik}\,d^4/dx^4 - X^{VW}_{ik}\,d^2/dx^2$, whereas in general their coefficients depend on $x$. The boundary terms of the integration by parts give the end conditions

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

and the generalised shear stress resultant $W^{\tau}_i = -W^{\sigma}_{i,x} + W^{xs}_i$, $W^{xs}_i = H^{I}_{ki}\phi_{k,xx} + D^{I}_{ik}\phi_{k,x} + F^{I}_{ki}\phi_k + H^{II}_{ji}\varphi_{j,x} + D^{I}_{ji}\varphi_j$. Eqs. (@equil)–(@Wsig) show the structure of the problem; since their coefficients vary with $x$ and the end conditions couple all modes, they are solved by the Ritz method, which works directly with the energy.

### Ritz solution, end conditions and classical limit

As in {{kharghani2025jcs,kharghani2026jcs}}, every amplitude function is approximated by a sum of known basis functions with unknown coefficients,

$$
\phi_k(x) = \sum_{l=1}^{n_b} N_l(x)\,a_{kl},\qquad \varphi_j(x) = \sum_{l=1}^{n_b} N_l(x)\,a_{jl},\qquad \chi_r(x) = \sum_{l=1}^{n_b} N_l(x)\,a_{rl}, \eqn{Psi}
$$

so that $\boldsymbol\Phi = \boldsymbol\Psi\mathbf a$, with $\boldsymbol\Psi$ the matrix of the basis functions and their derivatives and $\mathbf a$ the vector of all coefficients. Instead of the trigonometric series of {{kharghani2025jcs,kharghani2026jcs}}, the $n_b$ basis functions $N_l$ are cubic B-splines on an open uniform knot vector with spans $h_x = b_s/20$. B-splines are local, so that they resolve the zones near the loaded edges, and their end values are given directly by the end coefficients ($\phi(0) = a_1$; $\phi_{,x}(0) = 0 \Leftrightarrow a_2 = a_1$), so that any end condition is imposed exactly. Cubic B-splines are established in spline finite strips {{lau1986}} and isogeometric analysis {{hughes2005iga}}; like the Hermite elements of GBT-based finite elements {{bebiano2018}} they are piecewise cubic, but twice continuously differentiable.

The loaded edges are supported as in the shell models (Section 3). **Simply supported (SS)** edges have $U_Y = 0$ at every end node, i.e. $\sum_k U_{Y,qk}\phi_k = 0$ at $x = 0$ and $x = L$ for every node $q$, which is imposed by reducing $\mathbf a$ to the null space of these conditions; $U_X$, the warping and the rotations remain free. Unlike the classical condition $\phi_k = \sin(m\pi x/L)$, which suppresses every in-plane mode at the ends, this condition suppresses only $U_Y$. **Clamped (CC)** edges have $\phi_k = \phi_{k,x} = 0$ and $\chi_r = 0$ at both ends and $\varphi_j(0) = 0$, so that $x = 0$ is fully clamped and $x = L$ can warp under the applied load. Inserting $\boldsymbol\Phi = \boldsymbol\Psi\mathbf a$ into Eqs. (@Xide) and (@Xgen) turns $\delta U + \lambda\,\delta V = 0$ into the matrix eigenvalue problem

$$
\left(\mathbf K + \lambda\,\mathbf G\right)\mathbf a = \mathbf 0,\qquad
\mathbf K = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol\Xi\,\boldsymbol\Psi\,dx,\qquad
\mathbf G = \int_0^L\boldsymbol\Psi^{\mathsf T}\boldsymbol{\mathcal X}(x)\,\boldsymbol\Psi\,dx, \eqn{eig}
$$

with the stiffness matrix $\mathbf K$, the geometric matrix $\mathbf G$, $\boldsymbol\Xi = [\boldsymbol\Xi^{(de)}]$ and $\boldsymbol{\mathcal X} = [\boldsymbol{\mathcal X}^{(de)}]$; the integrals are evaluated with four Gauss points per span. The pre-buckling forces are computed for a reference load of $n_0 = 1$ N/mm, so that the smallest positive eigenvalue $\lambda$ is the critical line load $N_{cr}$ in N/mm. Since $\mathbf G$ is indefinite, the eigenvalues are computed in the form $-\mathbf G\mathbf a = \mu\mathbf K\mathbf a$ with $\mu = 1/\lambda$, whose largest $\mu$ gives the smallest positive $\lambda$, by the implicitly restarted Lanczos method {{lehoucq1998}}. For $b_s = 150$ mm and $L = 225$ mm (33 B-splines), the SS model has 18,230 unknowns.

Classical GBT is a special case, used for verification (Section 3.1): for uncoupled walls ($A_{16} = A_{26} = 0$, $\mathbf B = \mathbf 0$, $D_{16} = D_{26} = 0$), a uniform $S_{xx}$ and the conventional modes only, $H_{ik} = F_{ik} = 0$ and Eq. (@equil) becomes

$$
C_{ik}\,\phi_{k,xxxx} - D_{ik}\,\phi_{k,xx} + B_{ik}\,\phi_k + \lambda\left(X^{U}_{ik}\,\phi_{k,xxxx} - X^{VW}_{ik}\,\phi_{k,xx}\right) = 0, \eqn{classical}
$$

which, without $X^{U}_{ik}$, is the classical second-order GBT equation {{schardt1994,silvestre2002b}}. With $\phi_k = \hat\phi_k\sin(k_m x)$, $k_m = m\pi/L$, it becomes an algebraic eigenvalue problem for each number $m$ of half-waves,

$$
\left[B_{ik} + k_m^2\,D_{ik} + k_m^4\,C_{ik} + \lambda\left(k_m^2\,X^{VW}_{ik} + k_m^4\,X^{U}_{ik}\right)\right]\hat\phi_k = 0, \eqn{signature}
$$

whose smallest eigenvalue over $m$ as a function of $L$ is the signature curve. The periodic condition does not change this structure; it only restricts the admissible modes.

## Mode decomposition

The buckling mode $\mathbf a$ is decomposed into the contributions of the mode families $g\in\{G,D,L,T,S,M\}$, as in {{kharghani2026jcs}}. With $\mathbf a_g$ the vector $\mathbf a$ in which all coefficients not belonging to family $g$ are set to zero, the strain energy $U = \tfrac12\mathbf a^{\mathsf T}\mathbf K\mathbf a$ of the mode is divided into

$$
U_g = \tfrac12\,\mathbf a_g^{\mathsf T}\mathbf K\,\mathbf a,\qquad \rho_g = \frac{U_g}{U},\qquad \sum_g\rho_g = 1 . \eqn{rho}
$$

$\rho_g$ is the participation of family $g$; it includes half of the coupling energy between $g$ and every other family, so that the participations add up exactly to one. Kharghani and Mittelstedt {{kharghani2026jcs}} use the self-energies $\tfrac12\mathbf a_g^{\mathsf T}\mathbf K\mathbf a_g$ normalised by their sum, which omit the coupling energy. This difference matters in the cell, where the D modes move the wall junctions and the L modes restore the continuity of the walls, so that both are strongly coupled: for the Unsym cell ($b_s = 150$ mm, $\beta = 1.5$, SS) the normalised self-energies assign 41 % to D and 58 % to L, Eq. (@rho) 1 % and 97 %. Since $\mathbf K\mathbf a = -\lambda\mathbf G\mathbf a$ and the Mindlin modes do not enter $\mathbf G$, $\rho_M = 0$ identically; the influence of the Mindlin modes is assessed by omitting them (Section 3.6). The same partition applies to individual hierarchical modes.
