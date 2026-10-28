
## Strain energy and GBT modal matrices

With Eq. (@strainop), the strain energy of the cell becomes

$$
\begin{gathered}
U = \frac12\int_0^L\!\!\int_\Gamma \hat{\boldsymbol\varepsilon}^{\mathsf T}\mathbf C\,\hat{\boldsymbol\varepsilon}\,ds\,dx
= \frac12\int_0^L \sum_{d=0}^{2}\sum_{e=0}^{2}\mathbf q^{(d)\mathsf T}\mathbf A^{(de)}\,\mathbf q^{(e)}\,dx,\\
\mathbf A^{(de)} = \int_\Gamma \mathbf E^{(d)\mathsf T}\,\mathbf C\,\mathbf E^{(e)}\,ds,
\end{gathered} \eqn{energy}
$$

where $\Gamma$ is the mid-line of the cross-section and $\mathbf A^{(ed)} = \mathbf A^{(de)\mathsf T}$. In matrix form, with $\boldsymbol\Phi = \left[\mathbf q;\ \mathbf q_{,x};\ \mathbf q_{,xx}\right]$,

$$
\begin{gathered}
U = \frac12\int_0^L \boldsymbol\Phi^{\mathsf T}\,\boldsymbol\Xi\,\boldsymbol\Phi\,dx,\\
\boldsymbol\Xi = \begin{bmatrix}\bar{\mathbf B} & \bar{\mathbf F}_1 & \bar{\mathbf D}_2\\ \bar{\mathbf F}_1^{\mathsf T} & \bar{\mathbf D}_1 & \bar{\mathbf H}_1\\ \bar{\mathbf D}_2^{\mathsf T} & \bar{\mathbf H}_1^{\mathsf T} & \bar{\mathbf C}\end{bmatrix}
= \begin{bmatrix}\mathbf A^{(00)} & \mathbf A^{(01)} & \mathbf A^{(02)}\\ \mathbf A^{(10)} & \mathbf A^{(11)} & \mathbf A^{(12)}\\ \mathbf A^{(20)} & \mathbf A^{(21)} & \mathbf A^{(22)}\end{bmatrix},
\end{gathered} \eqn{xi}
$$

which is the GBT matrix $\boldsymbol\Xi$ of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}} written for the full set of modes. For the conventional modes $k, l \in \mathcal C$, Table 2 and Eq. (@abd) give the classical GBT modal matrices of a laminated wall with coupling:

$$
\bar B_{kl} = \int_\Gamma\Big[A_{22}\,\bar v_{k,s}\bar v_{l,s} - B_{22}\left(\bar v_{k,s}\bar w_{l,ss}+\bar w_{k,ss}\bar v_{l,s}\right) + D_{22}\,\bar w_{k,ss}\bar w_{l,ss}\Big]ds, \eqn{Bmat}
$$

$$
\bar C_{kl} = \int_\Gamma\Big[A_{11}\,\bar u_k\bar u_l - B_{11}\left(\bar u_k\bar w_l+\bar w_k\bar u_l\right) + D_{11}\,\bar w_k\bar w_l\Big]ds, \eqn{Cmat}
$$

$$
\begin{gathered}
\bar D_{1,kl} = \int_\Gamma\Big[A_{66}\,\gamma_k\gamma_l - 2B_{66}\left(\gamma_k\bar w_{l,s}+\bar w_{k,s}\gamma_l\right) + 4D_{66}\,\bar w_{k,s}\bar w_{l,s}\Big]ds,\\
 \gamma_k = \bar u_{k,s}+\bar v_k,
\end{gathered} \eqn{D1mat}
$$

$$
\bar D_{2,kl} = \int_\Gamma\Big[A_{12}\,\bar v_{k,s}\bar u_l - B_{12}\left(\bar v_{k,s}\bar w_l+\bar w_{k,ss}\bar u_l\right) + D_{12}\,\bar w_{k,ss}\bar w_l\Big]ds, \eqn{D2mat}
$$

$$
\bar F_{1,kl} = \int_\Gamma\Big[A_{26}\,\bar v_{k,s}\gamma_l - B_{26}\left(2\bar v_{k,s}\bar w_{l,s}+\bar w_{k,ss}\gamma_l\right) + 2D_{26}\,\bar w_{k,ss}\bar w_{l,s}\Big]ds, \eqn{F1mat}
$$

$$
\bar H_{1,kl} = \int_\Gamma\Big[A_{16}\,\gamma_k\bar u_l - B_{16}\left(\gamma_k\bar w_l+2\bar w_{k,s}\bar u_l\right) + 2D_{16}\,\bar w_{k,s}\bar w_l\Big]ds. \eqn{H1mat}
$$

The transverse-extension stiffness enters through the $A_{22}$ term of $\bar{\mathbf B}$, the bending–twisting coupling through $\bar{\mathbf F}_1$ and $\bar{\mathbf H}_1$, and the membrane–bending coupling through all $B_{ij}$ terms. For the shear modes $j, j' \in \mathcal S$, Table 2 gives

$$
A^{(00)}_{jj'} = \int_\Gamma A_{66}\,\bar u_{j,s}\bar u_{j',s}\,ds,\qquad
A^{(01)}_{jj'} = \int_\Gamma A_{16}\,\bar u_{j,s}\bar u_{j'}\,ds,\qquad
A^{(11)}_{jj'} = \int_\Gamma A_{11}\,\bar u_{j}\bar u_{j'}\,ds, \eqn{shearmat}
$$

which correspond to the matrices $\check{\mathbf D}_1$, $\check{\mathbf H}_1$ and $\check{\mathbf C}$ of shear-deformable GBT {{silvestre2013shear,kharghani2025jcs}}, while the conventional–shear coupling blocks are the matrices $\tilde{\mathbf F}$, $\tilde{\mathbf D}$, $\tilde{\mathbf H}$ and $\tilde{\mathbf C}$ of the same references, evaluated with the strains of Table 2 (e.g. $A^{(12)}_{jk} = \int_\Gamma[A_{11}\bar u_j\bar u_k - B_{11}\bar u_j\bar w_k]ds$ for the extension of a shear mode against the warping and bending of a conventional mode). The Mindlin modes add the transverse shear energy $\int_\Gamma[\bar c_m A_{s,11}\bar c_{m'} + \bar b_m A_{s,22}\bar b_{m'} + \dots]ds$ and their bending and twisting contributions through $\mathbf D$ and $\mathbf B$. In the implementation, all blocks are evaluated by the same Gauss quadrature over the segments from the operators $\mathbf E^{(d)}$ of Table 2, so that no coupling term is omitted.

## Amplitude functions and end conditions

### Ritz approximation along the cell

As in the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}, the amplitude functions are approximated by admissible functions with unknown coefficients, $\boldsymbol\Phi = \boldsymbol\Psi\,\mathbf a$. Instead of global trigonometric series, cubic B-splines on an open, uniform knot vector are used,

$$
q_k(x) = \sum_{b=1}^{n_b} N_b(x)\,a_{kb}, \eqn{spline}
$$

with spans of length $h_x = \min(15\ \text{mm},\ \max(2\ \text{mm},\ W/40),\ L/10)$, i.e. $h_x = b_s/20$ for all cases studied. Local support is essential here, because the pre-buckling state of the cell varies strongly near the loaded edges (Section 3.7), which a single trigonometric series cannot represent efficiently, and because the B-splines represent simply supported and clamped edges alike. Inserting Eq. (@spline) into Eq. (@energy) gives the stiffness matrix as a sum of Kronecker products,

$$
\mathbf K = \sum_{d=0}^{2}\sum_{e=0}^{2}\mathbf A^{(de)}\otimes\mathbf X^{(de)},\qquad X^{(de)}_{ab} = \int_0^L \frac{d^dN_a}{dx^d}\,\frac{d^eN_b}{dx^e}\,dx, \eqn{kmat}
$$

with four Gauss points per span. For $b_s = 150$ mm and $L = 225$ mm the model has $569\times 33 = 18\,777$ coefficients. With open knots, $q_k(0) = a_{k1}$ and $q_k(L) = a_{kn_b}$, and $q_{k,x}(0) = 0$ is equivalent to $a_{k2} = a_{k1}$, so that the essential end conditions act directly on the end coefficients.

### End conditions

The loaded edges are supported as in the Abaqus models:

* **SS** (out-of-plane support of every end edge, $U_2 = 0$): $U_Y = 0$ at every node at $x = 0$ and $x = L$, i.e. $\sum_{k\in\mathcal C}U_{Y,ik}\,\varphi_k(0) = \sum_{k\in\mathcal C}U_{Y,ik}\,\varphi_k(L) = 0$ for every node $i$. These linear conditions on the end coefficients are imposed by a null-space reduction; the in-plane displacement $U_X$, the warping and the rotations remain free. Rigid-body motions are removed by $\varphi_a(0) = 0$ and, for periodic side edges, $U_X(s_R, 0) = 0$.
* **CC** (clamped loaded edges): $\varphi_k = \varphi_{k,x} = 0$ at both ends for every conventional mode (in-plane displacements and the rotations $w_{,x}$ and $\theta$ vanish), $\chi_m = \vartheta_m = 0$ at both ends, and $\varphi_a(0) = \varphi_j(0) = 0$, so that $u = 0$ on the whole edge $x = 0$; the edge $x = L$ is free to warp under the uniform load.

In classical GBT, simply supported ends are represented by $\varphi_k = \sin(m\pi x/L)$, which suppresses all in-plane displacements of the end cross-sections while leaving their warping free (rigid end diaphragms). The present SS condition is weaker: like the Abaqus model it only suppresses the out-of-plane displacement $U_Y$ of the end edges. The classical end condition is recovered by prescribing $\varphi_k(0) = \varphi_k(L) = 0$ for all conventional modes.

## Loads and pre-buckling state

The external work of the end loads, a line load $q(s)$ per unit length acting at a distance $e(s)$ from the wall reference surface, is

$$
W_{ext} = \int_\Gamma q(s)\left[u_q(s,0)-u_q(s,L)\right]ds,\qquad u_q = u - e(s)\,w_{,x}, \eqn{wext}
$$

where $u_q$ is the axial displacement of the fibre at which the load acts ($e = c/2$ on a foot, $e = 0$ elsewhere). Inserting Eqs. (@kin) and (@spline) gives the load vector

$$
\begin{aligned}
f_{\varphi_k,b} &= \left(\int_\Gamma q\,\bar u_k\,ds - \int_\Gamma q\,e\,\bar w_k\,ds\right)\left[N_{b,x}(0)-N_{b,x}(L)\right],\\
f_{\varphi_j,b} &= \int_\Gamma q\,\bar u_j\,ds\left[N_b(0)-N_b(L)\right],\qquad
f_{\varphi_a,b} = \int_\Gamma q\,ds\left[N_b(0)-N_b(L)\right]
\end{aligned} \eqn{load}
$$

The pre-buckling state is not assumed: it is the first-order solution $\mathbf K\mathbf a^0 = \mathbf f$ under the end and side conditions above. The membrane resultants of this state, compression positive, are

$$
\left[P_{xx},\ P_{ss},\ P_{xs}\right]^{\mathsf T} = -\left(\mathbf A\,\boldsymbol\varepsilon^{0}_m + \mathbf B\,\boldsymbol\kappa^{0}\right), \eqn{pre}
$$

with the membrane strains $\boldsymbol\varepsilon^0_m$ and curvatures $\boldsymbol\kappa^0$ of $\mathbf a^0$. They vary along $s$ and $x$: the load is introduced at the ends and shared between skin and stringer by shear lag, the eccentric foot load bends the cell near the ends, and the restrained Poisson expansion produces a transverse resultant $P_{ss}$. Periodic and symmetry side edges both prevent a net lateral expansion of the cell and therefore both produce $P_{ss}$; they differ in the shear and rotation of the side edges.

## Geometric stiffness and eigenvalue problem

### Geometric stiffness

The second-order membrane energy of the pre-buckling resultants is

$$
\begin{aligned}
U_G &= \frac12\int_0^L\!\!\int_\Gamma \sum_{c\in\{u,v,w\}}\Big[P_{xx}\,g_{c,x}^2 + P_{ss}\,g_{c,s}^2 + 2P_{xs}\,g_{c,x}\,g_{c,s}\Big]ds\,dx,\\
g_{c,\alpha}(x,s) &= \sum_k\sum_{d}G^{(d)}_{k,c\alpha}(s)\,q^{(d)}_k(x),
\end{aligned} \eqn{ug}
$$

where $g_{c,\alpha}$ is the derivative of the displacement component $c$ with respect to $\alpha\in\{x,s\}$. For a conventional mode, $g_{u,x} = \bar u_k\varphi_{k,xx}$, $g_{u,s} = \bar u_{k,s}\varphi_{k,x}$, $g_{v,x} = \bar v_k\varphi_{k,x}$, $g_{v,s} = \bar v_{k,s}\varphi_k$, $g_{w,x} = \bar w_k\varphi_{k,x}$ and $g_{w,s} = \bar w_{k,s}\varphi_k$. All three displacement components and all three resultants are kept. With Eq. (@ug),

$$
\begin{gathered}
U_G = \frac12\int_0^L\sum_{d,e}\mathbf q^{(d)\mathsf T}\,\boldsymbol{\mathcal X}^{(de)}(x)\,\mathbf q^{(e)}\,dx,\\
\mathcal X^{(de)}_{kl}(x) = \sum_{t}\int_\Gamma P_t(s,x)\sum_{c}G^{(d)}_{k,c\alpha_t}(s)\,G^{(e)}_{l,c\beta_t}(s)\,ds,
\end{gathered} \eqn{xmat}
$$

with $(P_t,\alpha_t,\beta_t)\in\{(P_{xx},x,x),(P_{ss},s,s),(P_{xs},x,s),(P_{xs},s,x)\}$. Since the resultants depend on $x$, the geometric matrices $\boldsymbol{\mathcal X}^{(de)}(x)$ are evaluated at every Gauss point along the length, and the geometric stiffness matrix is

$$
\mathbf K_G = \sum_{d,e}\int_0^L\left(\mathbf N^{(d)}\mathbf N^{(e)\mathsf T}\right)\otimes\boldsymbol{\mathcal X}^{(de)}(x)\,dx, \eqn{kg}
$$

with $\mathbf N^{(d)}$ the vector of the $d$-th derivatives of the B-splines. The Mindlin modes do not move the mid-surface and do not enter $\mathbf K_G$.

### Eigenvalue problem and modal participation

The critical load follows from the stationarity of $U - \lambda U_G$:

$$
\left(\mathbf K - \lambda\,\mathbf K_G\right)\mathbf a = \mathbf 0,\qquad N_{cr} = \lambda_1\times 1\ \text{N/mm}, \eqn{eig}
$$

on the coefficients that satisfy the end conditions. Since $\mathbf K_G$ is indefinite (the pre-buckling state contains membrane shear and some transverse tension), the problem is solved as $\mathbf K_G\,\mathbf a = \mu\,\mathbf K\,\mathbf a$ for the largest $\mu = 1/\lambda$, with the implicitly restarted Lanczos method {{lehoucq1998}} and a sparse direct factorisation of $\mathbf K$ (tolerance $10^{-8}$).

The participation of a mode family $g$ in the buckling mode is measured as in {{kharghani2026jcs}} by the share of the displacement field that the family produces,

$$
p_g = \frac{\int_0^L \lVert\mathbf F_g(x)\rVert\,dx}{\sum_h\int_0^L\lVert\mathbf F_h(x)\rVert\,dx}, \eqn{part}
$$

where $\lVert\mathbf F_g(x)\rVert$ is the root mean square over the cross-section of the displacements produced by family $g$ alone.

## Classical limits

The formulation contains classical GBT as a special case. For the conventional modes alone, uncoupled orthotropic walls ($\mathbf B = \mathbf 0$, $A_{16} = A_{26} = D_{16} = D_{26} = 0$) and a uniform pre-buckling resultant $P_{xx} = N^0$ with $P_{ss} = P_{xs} = 0$, the odd blocks vanish, $\bar{\mathbf F}_1 = \bar{\mathbf H}_1 = \mathbf 0$, and the geometric matrices of Eq. (@xmat) reduce to

$$
\begin{gathered}
\boldsymbol{\mathcal X}^{(11)} = \boldsymbol{\mathcal X}^{VW},\quad \mathcal X^{VW}_{kl} = \int_\Gamma N^0\left(\bar v_k\bar v_l+\bar w_k\bar w_l\right)ds,\\
\boldsymbol{\mathcal X}^{(22)} = \boldsymbol{\mathcal X}^{U},\quad \mathcal X^{U}_{kl} = \int_\Gamma N^0\,\bar u_k\bar u_l\,ds,
\end{gathered} \eqn{xclassical}
$$

which are the geometric matrices of the GBT–Ritz approach {{kharghani2025jcs,kharghani2026jcs}}. The Euler–Lagrange equations of $U - \lambda U_G$ with Eqs. (@xi) and (@xclassical) are

$$
\begin{gathered}
\bar{\mathbf C}\,\boldsymbol\varphi_{,xxxx} - \bar{\mathbf D}\,\boldsymbol\varphi_{,xx} + \bar{\mathbf B}\,\boldsymbol\varphi + \lambda\left(\boldsymbol{\mathcal X}^{VW}\boldsymbol\varphi_{,xx} - \boldsymbol{\mathcal X}^{U}\boldsymbol\varphi_{,xxxx}\right) = \mathbf 0,\\
\bar{\mathbf D} = \bar{\mathbf D}_1 - \bar{\mathbf D}_2 - \bar{\mathbf D}_2^{\mathsf T},
\end{gathered} \eqn{classical}
$$

where the integration by parts of the $\bar{\mathbf D}_2$ terms gives the classical torsion-related matrix $\bar{\mathbf D}$. Without the warping contribution $\boldsymbol{\mathcal X}^U$, Eq. (@classical) is the classical GBT buckling equation $\bar{\mathbf C}\boldsymbol\varphi'''' - \bar{\mathbf D}\boldsymbol\varphi'' + \bar{\mathbf B}\boldsymbol\varphi + \lambda\boldsymbol{\mathcal X}\boldsymbol\varphi'' = \mathbf 0$ of Schardt {{schardt1994}} and of Silvestre and Camotim for orthotropic members {{silvestre2002b}}. For simply supported ends and $\boldsymbol\varphi = \hat{\boldsymbol\varphi}\sin k_mx$, $k_m = m\pi/L$, it becomes the algebraic eigenproblem

$$
\left[\bar{\mathbf B} + k_m^2\bar{\mathbf D} + k_m^4\bar{\mathbf C} - \lambda\left(k_m^2\boldsymbol{\mathcal X}^{VW} + k_m^4\boldsymbol{\mathcal X}^{U}\right)\right]\hat{\boldsymbol\varphi} = \mathbf 0, \eqn{signature}
$$

whose lowest eigenvalue over $m$ is the classical GBT signature curve. The periodic condition of Section 3.4 does not alter this structure: it only restricts the admissible modes of the cross-section. Equation (@signature) also explains the hierarchical ordering used in Figure 2. The present model differs from the classical one only through the terms that the panel requires and that classical GBT omits: laminate coupling ($\bar{\mathbf F}_1$, $\bar{\mathbf H}_1$, $B_{ij}$), shear and Mindlin modes, a non-uniform pre-buckling state with $P_{ss}$ and $P_{xs}$, and end conditions that support the loaded edges without restraining their warping.

## Numerical aspects

*Empty modes.* The Mindlin vectors $\boldsymbol\Omega_i$ are an orthonormal basis of the span of the tangents of the segments meeting at node $i$. At some intermediate nodes of the inclined webs, the two segment tangents are identical up to round-off, and a rank detection without tolerance keeps a spurious second vector normal to the web. For it $\bar c_m = \mathbf t\cdot\boldsymbol\Omega\approx 10^{-15}$ on every segment, so that the mode carries neither strain energy ($A^{(de)}_{kk}\approx10^{-25}$ against $10^9$ for the stiffest mode) nor displacement gradients. There are 11 such modes for $b_s = 150$ mm and 8 for $b_s = 285.5$ mm. They cannot change any eigenvalue, but left free they make $\mathbf K$ singular; they are identified by $\sum_{d,e}|A^{(de)}_{kk}|\le10^{-12}\max_l\sum_{d,e}|A^{(de)}_{ll}|$ together with $\mathbf G_k\equiv\mathbf 0$ and their coefficients are held. With this treatment the result does not depend on the sparse solver (e.g. 31.2431 N/mm for Cross-ply, SS, periodic, $b_s = 150$ mm, $\beta = 1.5$, with both PARDISO and SuperLU).

*Sparse geometric stiffness.* At a cross-section Gauss point $s_p$ only the modes with support there are non-zero (about 25 of 569 for $b_s = 150$ mm). Writing Eq. (@xmat) with quadrature, $\mathcal X^{(de)}_{kl}(x) = \sum_t\sum_p P_t(s_p,x)\,M^{(de)}_{t,p,(kl)}$, the sparse factors $M^{(de)}_{t,p,(kl)} = w_p\sum_c G^{(d)}_{k,c\alpha_t}(s_p)G^{(e)}_{l,c\beta_t}(s_p)$ depend only on the cross-section and are computed once; for every panel length only $P_t(s_p,x_g)$ changes.

