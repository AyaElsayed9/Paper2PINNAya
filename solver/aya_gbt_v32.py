"""
================================================================================
  AYA-GBT  V32  --  buckling of thin-walled laminated members and of the
  omega-stiffened panel  (October 2026)
================================================================================
  TWO SOLVERS, both with the Abaqus end conditions and their own first-order
  pre-buckling state:
    GBT   : shear-deformable Generalised Beam Theory (Silvestre & Camotim 2013
            + transverse-extension T, Mindlin M and nodal-rotation R modes),
            amplitude functions in cubic B-splines.  The default.
              members c, l, z, lipc, t, i, box, hat, dhat : sections 23-24
              omega-stiffened panel                       : section 25
    strip : shell-type Mindlin spline finite strip -- the fast reference
            ("oracle"): members section 20, omega panel section 22.
  Accuracy against the Abaqus models: members (144 cases C, Z, L, T, I, RHS)
  GBT 1.07 % / strip 0.62 % mean; omega (60 cases, symmetry sides) GBT 0.58 % /
  strip 0.57 %; omega periodic SS (30 cases) GBT 0.57 % / strip 0.59 %.

  The classical part (closed-form Euler / torsion / plate loads and the
  classical GBT deformation modes, sections 3-12) only gives the reference
  loads and the deformation-mode figures; it never gives the reported P_cr.

  RUN
    python AYA-GBT-V32.py                    interactive questionnaire
    python AYA-GBT-V32.py --answers=FILE     same, answers read from FILE
                                             (one per line, empty = default)
    --parametric-c / -z / -l / -t / -i / -rhs / --parametric-sections
    --parametric-omega [--omega-periodic]    parametric studies (Abaqus plans)
    --solver=strip | --omega-solver=strip    use the strip instead of GBT
    --gbt=CSTMR  --omega-gbt-fine  --mesh-fast / --mesh-fine
    --omega-mesh-fast / --omega-mesh-fine    discretisation options
    --demo                                   defaults, no questions
  Speed: install pypardiso (Intel MKL PARDISO), ~10x faster factorisation.
================================================================================
"""

import sys                                   # to read --demo and detect a terminal
import time                                  # wall-clock runtime measurement
import platform                              # host hardware / OS description
import numpy as np                           # all linear algebra / integration
import matplotlib                            # plotting
matplotlib.use("Agg")                        # write PNG files (no screen needed)
import matplotlib.pyplot as plt              # the pyplot interface

# ---- publication-quality figure style (journal-paper look) -----------------
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,                      # crisp, print-ready
    "savefig.bbox": "tight",
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.linewidth": 1.0,
    "axes.grid": True,
    "grid.alpha": 0.30,
    "grid.linewidth": 0.6,
    "lines.linewidth": 1.8,
    "legend.fontsize": 9,
    "legend.frameon": True,
    "legend.framealpha": 0.9,
    "legend.edgecolor": "0.7",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "mathtext.fontset": "cm",
})
# a clean, colour-blind-friendly palette used consistently across figures
_GBT_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
               "#56B4E9", "#000000"]


# DEMO = run on the defaults without questions (--demo and the parametric studies).
DEMO = ("--demo" in sys.argv) or any(_a.startswith("--parametric") for _a in sys.argv)

# ==============================================================================
# 0.  TINY INTERACTIVE HELPERS
# ==============================================================================
def ask(prompt, default, cast=str, choices=None):
    """Ask one question.  Show the default in [..]; <Enter> accepts it.

    cast    : str/int/float  -> converts the typed text
    choices : optional list of allowed (lower-cased) answers
    """
    while True:                                          # loop until a valid answer
        if DEMO:                                         # non-interactive ->
            return default                               #   just take the default
        raw = input(f"{prompt} [{default}]: ").strip()   # read + trim the line
        if raw == "":                                    # empty line ->
            return default                               #   accept the default
        try:
            val = cast(raw)                              # try to convert the text
        except ValueError:                               # not a number etc. ->
            print("   ! could not read that, try again."); continue
        if choices is not None and str(val).lower() not in choices:
            print(f"   ! please choose one of {choices}"); continue   # re-ask
        return val                                       # good answer -> return it


def ask_list(prompt, default, cast=float):
    """Read a comma/space separated list, e.g. ply angles '45,-45,45,-45'."""
    if DEMO:                                             # non-interactive ->
        return list(default)                             #   use the default list
    raw = input(f"{prompt} [{','.join(str(x) for x in default)}]: ").strip()
    if raw == "":                                        # empty -> default
        return list(default)
    for ch in "[]()":                                    # tolerate brackets if typed
        raw = raw.replace(ch, " ")
    return [cast(p) for p in raw.replace(",", " ").split()]   # split + convert


# ==============================================================================
# 1.  CLASSICAL LAMINATION THEORY (CLT)  ->  one wall's A, B, D stiffness
# ==============================================================================
def clt(material, layup, t_total):
    """Return the laminate membrane(A), coupling(B) and bending(D) 3x3 matrices
    of ONE wall whose total thickness is t_total and whose ply angles are `layup`.

    Calling this separately for each wall (with that wall's thickness) is what
    lets the section have a thick stiffener and a thinner base.
    """
    if material["type"] == "isotropic":                  # isotropic = special case
        E, nu = material["E"], material["nu"]            #   of orthotropic with
        E1 = E2 = E                                      #   E1=E2=E and
        G12 = E / (2 * (1 + nu)); nu12 = nu              #   G=E/2(1+nu)
    else:                                                # orthotropic lamina
        E1, E2 = material["E1"], material["E2"]
        G12, nu12 = material["G12"], material["nu12"]

    nu21 = nu12 * E2 / E1                                 # reciprocal Poisson ratio
    Delta = 1.0 - nu12 * nu21                             # plane-stress determinant
    Q11, Q22 = E1 / Delta, E2 / Delta                    # reduced stiffnesses
    Q12, Q66 = nu12 * E2 / Delta, G12                    #   (material axes)

    def Qbar(thdeg):                                     # rotate Q to wall (x,s) axes
        th = np.radians(thdeg); c, s = np.cos(th), np.sin(th)   # angle -> cos,sin
        b11 = Q11*c**4 + 2*(Q12+2*Q66)*s**2*c**2 + Q22*s**4     # transformed Q11
        b22 = Q11*s**4 + 2*(Q12+2*Q66)*s**2*c**2 + Q22*c**4     # transformed Q22
        b12 = (Q11+Q22-4*Q66)*s**2*c**2 + Q12*(c**4+s**4)       # transformed Q12
        b66 = (Q11+Q22-2*Q12)*s**2*c**2 + Q66*(c**2-s**2)**2    # transformed Q66
        b16 = (Q11-Q12-2*Q66)*s*c**3 + (Q12-Q22+2*Q66)*s**3*c   # shear coupling
        b26 = (Q11-Q12-2*Q66)*s**3*c + (Q12-Q22+2*Q66)*s*c**3   # shear coupling
        return np.array([[b11, b12, b16], [b12, b22, b26], [b16, b26, b66]])

    nply = len(layup)                                    # number of plies
    zb = np.linspace(-t_total/2, t_total/2, nply + 1)    # ply interface z-heights
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))   # accumulate
    for k, th in enumerate(layup):                       # loop plies bottom->top
        Qb = Qbar(th); z0, z1 = zb[k], zb[k+1]           # this ply's Qbar and z's
        A += Qb * (z1 - z0)                              # A = sum Qbar  dz      (1st)
        B += Qb * (z1**2 - z0**2) / 2.0                  # B = sum Qbar  z dz    (2nd)
        D += Qb * (z1**3 - z0**3) / 3.0                  # D = sum Qbar  z^2 dz  (3rd)
    return A, B, D                                       # the three 3x3 matrices


# ==============================================================================
# 2.  GEOMETRY  ->  discretized node path + per-segment wall index
# ==============================================================================
def build_section(corner_nodes, n_int_per_wall):
    """Insert the requested intermediate nodes into each wall.

    corner_nodes   : ordered (Y,Z) natural nodes [mm] (the fold lines / tips)
    n_int_per_wall : list -> intermediate nodes inside each wall  (THE per-segment
                     variable requested at the top of the code)

    returns
       path     : (Npts,2)  every node along the developed midline
       s_path   : arc-length of every node
       s_corner : arc-length of the natural (corner) nodes  (the fold lines)
       wall_seg : (Npts-1,) which wall each little segment belongs to (0,1,2,..)
    """
    corner_nodes = np.asarray(corner_nodes, float)       # to array
    pts = [corner_nodes[0]]                              # start at first natural node
    s = [0.0]; s_corner = [0.0]; wall_seg = []           # arc-lengths + wall map
    for w in range(len(corner_nodes) - 1):               # loop over walls
        p0, p1 = corner_nodes[w], corner_nodes[w + 1]    # this wall's two ends
        nseg = n_int_per_wall[w] + 1                     # sub-segments = int+1
        for j in range(1, nseg + 1):                     # walk along the wall
            p = p0 + (p1 - p0) * j / nseg                # evenly spaced new node
            pts.append(p)                                # store the node
            s.append(s[-1] + np.linalg.norm(p - pts[-2]))# accumulate arc-length
            wall_seg.append(w)                           # this segment is in wall w
        s_corner.append(s[-1])                           # mark the fold-line s
    return (np.array(pts), np.array(s),
            np.array(s_corner), np.array(wall_seg))


# ==============================================================================
# 3.  CROSS-SECTION  ->  MODULUS-WEIGHTED open thin-walled properties
# ==============================================================================
def _beam_a11(A, B=None, D=None):
    """Axial (beam) membrane stiffness of a wall per unit width [N/mm]: the
    wall carries N_zz only (N_ss = N_zs = 0, M = 0), so E t = 1 / a11 with
    a = (A - B D^-1 B)^-1.  For an isotropic wall this is E t (A11 itself
    would be E t / (1 - nu^2), the stiffness of a laterally restrained plate)."""
    A = np.asarray(A, float)
    if B is not None and D is not None and np.max(np.abs(B)) > 0:
        A = A - np.asarray(B, float) @ np.linalg.solve(np.asarray(D, float), np.asarray(B, float))
    return 1.0 / np.linalg.inv(A)[0, 0]


def sec_props(path, a11_seg, d66_seg):
    """Open-section properties, weighting each segment by its OWN membrane
    stiffness A11 (so a thick stiffener counts more than a thin base).

    a11_seg : per-segment membrane stiffness A11  (length = Npts-1)
    d66_seg : per-segment D66  (for the St-Venant torsion constant)

    For uniform stiffness these reduce to  A11 * (geometric Iyy, Izz, Iww...),
    reproducing the validated reference values exactly.

    Returns force-weighted rigidities directly:
       EA, EIyy, EIzz  (N, N*mm^2), GJ, ECw, plus centroid & shear-centre.
    """
    P = np.asarray(path, float)                          # node coordinates
    dl = np.linalg.norm(P[1:] - P[:-1], axis=1)          # length of each segment
    Lmid = dl.sum()                                      # total midline length
    w = a11_seg                                          # membrane weight per seg
    # weighted line integrals of a NODAL field along the segments (trapezoidal):
    Ilw = lambda a, b: np.sum(w*dl*(a + b)/2)                       # int w f ds
    Iqw = lambda a1, a2, b1, b2: np.sum(w*dl*(2*a1*b1+a1*b2+a2*b1+2*a2*b2)/6)  # int w f g ds
    y, z = P[:, 0], P[:, 1]                              # node Y and Z
    EA = np.sum(w*dl)                                    # axial rigidity sum A11 len
    Yc = Ilw(y[:-1], y[1:]) / EA                         # modulus-weighted centroid Y
    Zc = Ilw(z[:-1], z[1:]) / EA                         # modulus-weighted centroid Z
    Y1, Y2 = y[:-1]-Yc, y[1:]-Yc                         # centroid-relative Y (ends)
    Z1, Z2 = z[:-1]-Zc, z[1:]-Zc                         # centroid-relative Z (ends)
    EIyy = Iqw(Z1, Z2, Z1, Z2)                           # int A11 (Z-Zc)^2 ds
    EIzz = Iqw(Y1, Y2, Y1, Y2)                           # int A11 (Y-Yc)^2 ds
    EIyz = Iqw(Y1, Y2, Z1, Z2)                           # product rigidity
    dw = Y1*Z2 - Y2*Z1                                   # 2x swept-area increments
    wn = np.concatenate([[0.0], np.cumsum(dw)])          # raw sectorial coord (nodes)
    w1, w2 = wn[:-1].copy(), wn[1:].copy()               # segment-end values
    mw = Ilw(w1, w2)/EA; w1 -= mw; w2 -= mw              # remove weighted mean ONCE
    Iwy = Iqw(w1, w2, Y1, Y2); Iwz = Iqw(w1, w2, Z1, Z2) # sectorial products
    den = EIyy*EIzz - EIyz**2                            # determinant for SC solve
    ys = ( EIzz*Iwz - EIyz*Iwy)/den                      # shear-centre Y (rel C)
    zs = (-EIyy*Iwy + EIyz*Iwz)/den                      # shear-centre Z (rel C)
    wS1 = w1 + zs*Y1 - ys*Z1                             # principal sectorial coord
    wS2 = w2 + zs*Y2 - ys*Z2                             #   about the shear centre
    mS = Ilw(wS1, wS2)/EA; wS1 -= mS; wS2 -= mS          # remove its weighted mean
    ECw = Iqw(wS1, wS2, wS1, wS2)                        # warping rigidity int A11 w^2 ds
    GJ = np.sum(4*d66_seg*dl)                            # St-Venant torsion sum 4 D66 len
    w_nodes = np.concatenate([wS1, [wS2[-1]]])           # warping value at each node
    return dict(Lmid=Lmid, Yc=Yc, Zc=Zc, EA=EA, EIyy=EIyy, EIzz=EIzz, EIyz=EIyz,
                ys=ys, zs=zs, ECw=ECw, GJ=GJ, w_nodes=w_nodes,
                y=y, z=z, Y=y-Yc, Z=z-Zc)


# ==============================================================================
# 4.  LOCAL-PLATE BUCKLING  ->  TWO solver options
# ==============================================================================
# 4-point Gauss rule (for the finite-strip element integrals)


# ==============================================================================
# 4a. RIGOROUS LOCAL PLATE BUCKLING  ->  MULTI-HARMONIC FINITE STRIP
#     Captures D16/D26 (bending-twisting) AND B16/B26 (membrane-bending) EXACTLY,
#     with NO orthotropic knockdown and NO shear fudge factor.  A single harmonic
#     drops the shear-coupling stiffnesses because int sin(m pi x/L) cos(m pi x/L)
#     dx = 0 over one half-wave; coupling opposite-parity harmonics through
#     int sin(m..)cos(n..) dx (non-zero for m+n odd) recovers them rigorously.
#     Validated: isotropic SS plate -> k = 4.000 exactly;
#                [45/-45]s D16/D26 reduction -> 0.750 (mesh-converged).
# ==============================================================================
def _mh_shape(r, b):
    L1, L2 = 1-r, r; dL1, dL2 = -1/b, 1/b
    H1 = 1-3*r**2+2*r**3; H2 = b*(r-2*r**2+r**3); H3 = 3*r**2-2*r**3; H4 = b*(-r**2+r**3)
    dH1 = (-6*r+6*r**2)/b; dH2 = 1-4*r+3*r**2; dH3 = (6*r-6*r**2)/b; dH4 = -2*r+3*r**2
    ddH1 = (-6+12*r)/b**2; ddH2 = (-4+6*r)/b; ddH3 = (6-12*r)/b**2; ddH4 = (-2+6*r)/b
    return (np.array([L1, L2]), np.array([dL1, dL2]), np.array([H1, H2, H3, H4]),
            np.array([dH1, dH2, dH3, dH4]), np.array([ddH1, ddH2, ddH3, ddH4]))

_MHGP = np.array([-.8611363116, -.3399810436, .3399810436, .8611363116])
_MHGW = np.array([.3478548451, .6521451549, .6521451549, .3478548451])

def _mh_Jc(m, n, L):
    if (m+n) % 2 == 0:
        return 0.0
    return (2.0*L/np.pi)*m/(m*m-n*n)

def _mh_blocks(b, ABD, ks):
    A, B, D = ABD[0], ABD[1], ABD[2]
    C6 = np.array([[A[0,0],A[0,1],A[0,2],B[0,0],B[0,1],B[0,2]],
                   [A[0,1],A[1,1],A[1,2],B[0,1],B[1,1],B[1,2]],
                   [A[0,2],A[1,2],A[2,2],B[0,2],B[1,2],B[2,2]],
                   [B[0,0],B[0,1],B[0,2],D[0,0],D[0,1],D[0,2]],
                   [B[0,1],B[1,1],B[1,2],D[0,1],D[1,1],D[1,2]],
                   [B[0,2],B[1,2],B[2,2],D[0,2],D[1,2],D[2,2]]], float)
    Si = [0, 1, 3, 4]; Ci = [2, 5]
    Css = C6[np.ix_(Si, Si)]; Ccc = C6[np.ix_(Ci, Ci)]; Csc = C6[np.ix_(Si, Ci)]
    M = len(ks); u = [0, 4]; v = [1, 5]; w = [2, 3, 6, 7]
    Kd = [np.zeros((8, 8)) for _ in range(M)]; Kg = [np.zeros((8, 8)) for _ in range(M)]
    KSC = [[np.zeros((8, 8)) for _ in range(M)] for _ in range(M)]
    for gp, gw in zip((_MHGP+1)/2, _MHGW/2):
        Ln, dLn, Hn, dHn, ddHn = _mh_shape(gp, b); jac = b
        S = []; Cc = []; gW = []; gV = []
        for k in ks:
            ex=np.zeros(8);ey=np.zeros(8);kx=np.zeros(8);ky=np.zeros(8);gxy=np.zeros(8);kxy=np.zeros(8)
            for a, val in zip(u, -k*Ln): ex[a] += val
            for a, val in zip(v, dLn): ey[a] += val
            for a, val in zip(w, k**2*Hn): kx[a] += val
            for a, val in zip(w, -ddHn): ky[a] += val
            for a, val in zip(u, dLn): gxy[a] += val
            for a, val in zip(v, k*Ln): gxy[a] += val
            for a, val in zip(w, -2*k*dHn): kxy[a] += val
            S.append(np.vstack([ex, ey, kx, ky])); Cc.append(np.vstack([gxy, kxy]))
            gw_=np.zeros(8); gv_=np.zeros(8)
            for a, val in zip(w, k*Hn): gw_[a] += val
            for a, val in zip(v, k*Ln): gv_[a] += val
            gW.append(gw_); gV.append(gv_)
        for p in range(M):
            Kd[p] += gw*jac*(S[p].T@Css@S[p] + Cc[p].T@Ccc@Cc[p])
            Kg[p] += gw*jac*(np.outer(gW[p], gW[p]) + np.outer(gV[p], gV[p]))
            for q in range(M):
                KSC[p][q] += gw*jac*(S[p].T@Csc@Cc[q])
    return Kd, KSC, Kg

def mh_plate_ncr(width, ABD, L, free_edge=False, nstrip=8, mmax=None):
    """Rigorous critical uniaxial membrane force Ncr (N/mm) of ONE wall of given
    width and FULL laminate ABD, at member length L.  free_edge=True models an
    outstand (one longitudinal edge free); otherwise both edges simply supported.
    Captures D16/D26 and B16/B26 exactly."""
    if mmax is None:
        mmax = int(np.clip(np.ceil(2.2*L/width), 6, 20))
    mlist = list(range(1, mmax+1)); ks = [m*np.pi/L for m in mlist]
    N = nstrip+1; n1 = 4*N; M = len(mlist); ntot = M*n1
    Ke = np.zeros((ntot, ntot)); Kg = np.zeros((ntot, ntot))
    Kd, KSC, Kgeo = _mh_blocks(width/nstrip, ABD, ks)
    for s in range(nstrip):
        i, j = s, s+1; dof8 = [4*i,4*i+1,4*i+2,4*i+3, 4*j,4*j+1,4*j+2,4*j+3]
        for p in range(M):
            gp = [p*n1+d for d in dof8]
            Ke[np.ix_(gp, gp)] += (L/2)*Kd[p]; Kg[np.ix_(gp, gp)] += (L/2)*Kgeo[p]
            for q in range(M):
                if p == q or (mlist[p]+mlist[q]) % 2 == 0:
                    continue
                blk = _mh_Jc(mlist[p],mlist[q],L)*KSC[p][q] + _mh_Jc(mlist[q],mlist[p],L)*KSC[q][p].T
                Ke[np.ix_(gp, [q*n1+d for d in dof8])] += blk
    fixed = [2] if free_edge else [2, 4*nstrip+2]           # SS transverse edge(s)
    keep = np.ones(ntot, bool)
    for p in range(M):
        for d in fixed: keep[p*n1+d] = False
    Ke = Ke[np.ix_(keep, keep)]; Kg = Kg[np.ix_(keep, keep)]
    try:
        mu = np.linalg.eigvals(np.linalg.solve(Ke, Kg)).real    # mu = 1/lambda
    except Exception:
        return np.inf
    mu = mu[np.isfinite(mu) & (mu > 1e-12)]
    return 1.0/mu.max() if len(mu) else np.inf

_MH_ASPECTS = (0.7, 1.0, 1.5, 2.0, 3.0, 4.0)             # local half-wavelength scan

def rigorous_local_path(corners, wallABD, wall_seg, L, n_modes):
    """Rigorous per-wall local buckling (N/mm) for an open unbranched section.
    LOCAL buckling is length-INDEPENDENT (the wall buckles into its natural
    half-wavelength), so for each wall we take the minimum plate load over a scan
    of aspect ratios L_ref/b in {0.7..4} -- the fully-developed local plateau.
    End walls are outstands (one free edge); interior walls are simply supported
    at both fold lines.  D16/D26 & B16/B26 are captured EXACTLY by the multi-
    harmonic solver -- no knockdown, no shear fudge.  (L is accepted for interface
    compatibility but the local plateau does not depend on member length.)"""
    corners = np.asarray(corners, float); nw = len(corners)-1
    vals = []
    for w in range(nw):
        b_w = np.linalg.norm(corners[w+1]-corners[w])
        A, B, D = wallABD[w][0], wallABD[w][1], wallABD[w][2]
        free = (w == 0 or w == nw-1)                        # outstand end walls
        ncr = min(mh_plate_ncr(b_w, (A, B, D), asp*b_w, free_edge=free, mmax=8)
                  for asp in _MH_ASPECTS)                   # local plateau
        vals.append(ncr)
    vals = np.array(sorted(v for v in vals if np.isfinite(v)))
    out = np.full(n_modes, vals[-1] if len(vals) else np.inf)
    out[:min(n_modes, len(vals))] = vals[:n_modes]
    return out


def rigorous_local_graph(G):
    """Rigorous per-wall local plate buckling (N/mm) for branched/closed sections.
    A wall is an outstand (one free edge) if either end node is a free tip
    (valence 1); otherwise it is simply supported at both fold lines.  Returns the
    sorted DISTINCT per-wall local plateaus (D16/D26 & B16/B26 captured exactly)."""
    from collections import Counter
    nodes = np.asarray(G["nodes"], float); walls = G["walls"]; wABD = G["wABD"]
    val = Counter()
    for w in walls:
        val[w["i"]] += 1; val[w["j"]] += 1
    vals = []
    for wi, w in enumerate(walls):
        b_w = np.linalg.norm(nodes[w["j"]]-nodes[w["i"]])
        A, B, D, tw = wABD[wi]
        free = (val[w["i"]] == 1 or val[w["j"]] == 1)       # touches a free tip?
        ncr = min(mh_plate_ncr(b_w, (A, B, D), asp*b_w, free_edge=free, mmax=8)
                  for asp in _MH_ASPECTS)
        vals.append(ncr)
    return np.array(sorted(set(np.round([v for v in vals if np.isfinite(v)], 4))))


# ==============================================================================
# 5.  GLOBAL BUCKLING (closed form) + Engesser shear reduction
# ==============================================================================
def eng(PE, S):
    """Engesser/Timoshenko shear-reduced load: 1/Ps = 1/PE + 1/S."""
    return PE * S / (PE + S)                              # harmonic-type combination


def trial_ratio(bc, L, n=4001):
    """Rayleigh-Ritz ratio  kappa = integral(phi'')^2 dx / integral(phi')^2 dx
    of the LONGITUDINAL trial function phi(x) that satisfies the end conditions.
    This kappa replaces (pi/L)^2 in every length-dependent global load, so the
    boundary condition enters through the trial function itself (energy method),
    not through a bolted-on effective-length factor.

        SS (pinned) : phi = sin(pi x / L)          -> phi = 0 at ends
        CC (clamped): phi = 1 - cos(2 pi x / L)    -> phi = 0 AND phi' = 0 at ends

    For a prismatic member these single terms are the EXACT buckling
    eigenfunctions, so kappa is exact:  SS -> (pi/L)^2 ,  CC -> (2 pi/L)^2.
    """
    x = np.linspace(0.0, L, n)                           # integration grid
    if bc == "CC":                                       # clamped-clamped
        w = 2*np.pi/L
        d1 = w*np.sin(w*x)                               # phi'
        d2 = w**2*np.cos(w*x)                            # phi''
    else:                                                # simply supported
        w = np.pi/L
        d1 = w*np.cos(w*x)                               # phi'
        d2 = -w**2*np.sin(w*x)                           # phi''
    I1 = np.trapezoid(d1**2, x)                          # integral (phi')^2
    I2 = np.trapezoid(d2**2, x)                          # integral (phi'')^2
    return I2/I1                                         # kappa  [1/mm^2]


def global_loads(props, L, bc, shear, shear_rig):
    """All global critical loads at length L (everything in N), with rigorous
    principal-axis transformation for asymmetric sections (like Z-sections where EIyz != 0)."""
    if bc not in ("SS", "CC"):                           # tolerate Kbc numbers
        bc = {0.5: "CC"}.get(bc, "SS")
    kap = trial_ratio(bc, L)                             # trial-function ratio (=(pi/Le)^2)
    EIyy, EIzz, EIyz, ECw, GJ = props["EIyy"], props["EIzz"], props["EIyz"], props["ECw"], props["GJ"]
    
    # =========================================================================
    # PHYSICS FIX: Principal Axis Transformation (Crucial for Z-sections)
    # =========================================================================
    EI_avg = 0.5 * (EIyy + EIzz)
    R = np.sqrt((0.5 * (EIyy - EIzz))**2 + EIyz**2)
    EI_maj = EI_avg + R                                  # Principal major rigidity
    EI_min = EI_avg - R                                  # Principal minor rigidity

    y0 = abs(props["ys"])                                # centroid->shear-centre dist
    r0sq = (EIyy + EIzz)/props["EA"] + y0**2             # polar radius^2 about SC (Invariant)
    
    P2 = EI_maj * kap                                    # major-axis flexural (Principal)
    P3 = EI_min * kap                                    # minor-axis flexural (Principal)
    Pw = ECw * kap                                       # warping part of torsion
    P4 = (GJ + Pw) / r0sq                                # torsion (St-Venant+warping)
    
    if shear:                                            # apply shear softening:
        # Use average shear rigidity for principal flexural axes softening
        S_avg = 0.5 * (shear_rig["S2"] + shear_rig["S3"])
        P2 = eng(P2, S_avg)                              # major flexural
        P3 = eng(P3, S_avg)                              # minor flexural
        Pw = eng(Pw, shear_rig["Sw"])                    # warping torsion only
        P4 = (GJ + Pw) / r0sq                            # (GJ stays: it is shear)
        
    beta = 1.0 - y0**2/r0sq                              # flex-torsion coupling factor
    roots = np.roots([beta, -(P3+P4), P3*P4]).real       # quadratic for PFT
    PFT = min(roots[roots > 0]) if np.any(roots > 0) else min(roots)   # lower root
    a3 = abs(P4 - PFT); a4 = abs(P3 - PFT)               # modal participations
    tot = a3 + a4 + 1e-30                                # normalize
    part = (a3/tot, a4/tot)                              # (flexural %, torsional %)
    return dict(P2=P2, P3=P3, P4=P4, PFT=PFT, r0sq=r0sq, part=part)


# ==============================================================================
# 6.  ASSEMBLE EVERYTHING FOR ONE CONFIGURATION
# ==============================================================================


def auto_nint_for_walls(wall_lengths, target=6.0, floor=6, cap=40):
    """Unified cross-section discretisation (intermediate nodes per wall).

    Used by EVERY section family (A: C/L/Z/lipped-C, B: T/I, C: box/RHS,
    D: hat) so that the transverse mesh is consistent everywhere.

    The local and distortional plate modes of a composite wall need a finer
    transverse mesh than an isotropic one, because the deformed shape is skewed
    by the D16/D26 bend-twist coupling.  A node about every `target` mm, with a
    floor of `floor` per wall, gives converged loads: mesh-refinement checks show
    the angle-ply cases fall by 5-10% going from 6 to ~14 nodes per wall and then
    settle.  `cap` keeps very long walls from exploding the problem size.
    """
    out = []
    for wl in wall_lengths:
        wl = float(wl)
        out.append(int(min(cap, max(floor, round(wl / target)))))
    return out


def _analyse_graph(cfg):
    """Families B (branched), C (closed), D (closed+branched): topology solver.
    Rigorous EA/EI/centroid and open+Bredt torsion; warping/shear-centre tagged."""
    nodes, walls, cells = cfg["nodes"], cfg["walls"], cfg["cells"]
    G = build_graph(nodes, walls, cfg["material"], cfg["n_int_per_wall"])
    p = section_props_g(G)                               # EA, EI, centroid
    ys, zs, ECw, wconf = warping_g(G, cells)             # warping + shear centre
    GJ, GJo, GJc = torsion_GJ_g(G, cells)                # open + Bredt torsion
    S2, S3 = shear_rig_g(G)                              # flexural shear rigidities
    # torsion warping-shear rigidity: neglected (classical Vlasov assumption ->
    # no warping-shear softening).  For closed/branched sections warping is small,
    # so this is standard and removes the former 0.5*max(S2,S3) guess.
    Sw = 1e30
    props = dict(Lmid=p["Lmid"], EA=p["EA"], EIyy=p["EIyy"], EIzz=p["EIzz"],
                 EIyz=p["EIyz"], Yc=p["Yc"], Zc=p["Zc"], ys=ys, zs=zs, ECw=ECw, GJ=GJ,
                 _t_wall=[G["wABD"][wi][3] for wi in range(len(walls))])
    shear_rig = dict(S2=S2, S3=S3, Sw=Sw)
    L0 = cfg["L_design"]; Kbc = cfg["Kbc"]; shear = cfg["shear"]
    g = global_loads(props, L0, Kbc, shear, shear_rig)
    nloc = cfg["n_local"]
    # RIGOROUS per-wall local plate buckling (multi-harmonic; full A,B,D):
    loc_ncr = rigorous_local_graph(dict(G, wABD=G["wABD_raw"]))   # RAW laminate: full B, no double count
    eta_loc, gov_w = 1.0, 0                              # coupling now inside the solver
    Ploc = np.array(sorted(set(np.round(loc_ncr*props["Lmid"], 3))))
    nloc = min(nloc, len(Ploc))                          # honest count of distinct locals
    modes = {1: ("Axial extension", None, "loaded / prebuckling mode"),
             2: ("Major-axis flexural", g["P2"], "Euler bending"),
             3: ("Minor-axis flexural", g["P3"], "Euler bending"),
             4: ("Torsion", g["P4"], f"St-Venant(open {GJo:.3g}) + Bredt(cell {GJc:.3g}) + warping")}
    for j in range(nloc):
        modes[5+j] = (f"Local plate-{j+1}", Ploc[j], "per-wall plate load (distinct)")
    cand = {"3+4 Flexural-torsional": g["PFT"], "5 Local (lowest)": Ploc[0]}
    gov_name = min(cand, key=cand.get); gov_val = cand[gov_name]
    # ---- GOVERNING LOAD: unified shear-deformable spline finite strip -------
    # (section 20) -- the whole member with its actual end conditions and its
    # own pre-buckling state; the classical loads above are for interpretation.
    P_gbt, gbt_info, _ur = _unified_governing(cfg, L0)
    gov_name, gov_val = gbt_info["kind"], P_gbt
    n_natural = len(nodes); n_total = sum(cfg["n_int_per_wall"]) + n_natural
    n_walls = len(walls); m_inter = sum(cfg["n_int_per_wall"])
    n_global = 4; n_distort = max(len(cells)*1 + (n_natural - n_walls - 1), 0)
    return dict(nodes=np.asarray(nodes, float), walls=walls, cells=cells, graph=G,
                props=props, modes=modes, PFT=g["PFT"], part=(0.5, 0.5), r0sq=g["r0sq"],
                Ploc=Ploc, mloc=_ur.get("m"), gov_name=gov_name, gov_val=gov_val,
                _unified=_ur,
                wallABD=[G["wABD"][wi][:3] for wi in range(len(walls))],
                shear_rig=shear_rig, L0=L0, Kbc=Kbc, shear=shear, GJo=GJo, GJc=GJc,
                warp_conf=wconf, family=cfg.get("family", "?"),
                eta_loc=eta_loc, gov_wall=gov_w,
                n_natural=n_natural, n_total=n_total, n_walls=n_walls, m_inter=m_inter,
                n_global=n_global, n_distort=n_distort, n_conv_total=n_total,
                n_local_available=m_inter, nloc=nloc, is_graph=True,
                P_gbt=P_gbt, gbt_info=gbt_info)


def analyse(cfg):
    """Dispatch: family A -> validated path solver; B/C/D -> topology solver."""
    if cfg.get("family", "A") == "A" and "corner_nodes" in cfg:
        r = _analyse_path(cfg); r["is_graph"] = False; return r
    return _analyse_graph(cfg)


def _analyse_path(cfg):
    """Family A (open, unbranched): the validated ordered-path solver."""
    # ---- 6.1  geometry + per-segment wall map ------------------------------
    corners = cfg["corner_nodes"]                        # natural nodes
    path, s_path, s_corner, wall_seg = build_section(corners, cfg["n_int_per_wall"])
    nseg = len(path) - 1                                 # number of little segments

    # ---- 6.2  per-WALL laminates (variable thickness -> variable stiffness) -
    layup = cfg["layup"]; mat = cfg["material"]          # shared lay-up + material
    t_wall = cfg["t_wall"]                               # ONE thickness per wall
    
    # 1. Calculate classical A, B, D matrices
    wallABD_raw = [clt(mat, layup, tw) for tw in t_wall] 
    
    # 2. Account for Asymmetric/Unsymmetric Laminates using Reduced Stiffness (A* and D*)
    wallABD = []
    for A, B, D in wallABD_raw:
        if np.max(np.abs(B)) > 1e-5:                     # If B matrix exists (Unsymmetric)
            invA = np.linalg.inv(A)
            invD = np.linalg.inv(D)
            A_star = A - B @ invD @ B                    # Reduced membrane stiffness
            D_star = D - B @ invA @ B                    # Reduced bending stiffness
            wallABD.append((A_star, B, D_star))
        else:                                            # Symmetric layup (B=0)
            wallABD.append((A, B, D))

    # spread each wall's stiffness onto its segments using the (possibly reduced) matrices:
    a11_seg = np.array([_beam_a11(wallABD[wall_seg[e]][0]) for e in range(nseg)])  # beam E t (wallABD A is already B-reduced)
    a66_seg = np.array([wallABD[wall_seg[e]][0][2, 2] for e in range(nseg)])  # A66
    
    dmats   = [ (wallABD[wall_seg[e]][2][0, 0],          # D11
                 wallABD[wall_seg[e]][2][0, 1],          # D12
                 wallABD[wall_seg[e]][2][1, 1],          # D22
                 wallABD[wall_seg[e]][2][2, 2])          # D66
                for e in range(nseg) ]
                
    d66_seg = np.array([dm[3] for dm in dmats])          # D66 per segment

    # ---- 6.3  modulus-weighted section properties --------------------------
    props = sec_props(path, a11_seg, d66_seg)            # EA, EI, GJ, ECw, SC...
    props["_spath"] = s_path                             # store for the sweep
    props["_scorner"] = s_corner[1:-1]                   # interior folds (exclude tips)
    props["_dmats"] = dmats                              # per-strip D for the FSM
    props["_t_wall"] = t_wall                            # keep wall thicknesses
    bw, bf = cfg["b"], cfg["a"]                          # web height, flange width

    # ---- 6.4  per-mode shear rigidities  Shat = int A66 (du/ds)^2 ds -------
    dl  = np.linalg.norm(path[1:]-path[:-1], axis=1)     # segment lengths
    dl_ = np.where(dl == 0, 1, dl)                       # guard against zero
    dYds = (path[1:, 0]-path[:-1, 0])/dl_                # dY/ds per segment
    dZds = (path[1:, 1]-path[:-1, 1])/dl_                # dZ/ds per segment
    S2 = np.sum(a66_seg*dZds**2*dl)                      # EIyy-mode shear rigidity
    S3 = np.sum(a66_seg*dYds**2*dl)                      # EIzz-mode shear rigidity
    wsl = (props["w_nodes"][1:] - props["w_nodes"][:-1]) # delta-omega across each seg
    rho = wsl/dl_                                        # warping slope per segment
    Sw  = np.sum(a66_seg*rho**2*dl)                      # torsion warping shear rigidity
    shear_rig = dict(S2=S2, S3=S3, Sw=Sw)                # collect

    # ---- 6.5  design-length results ----------------------------------------
    L0 = cfg["L_design"]; Kbc = cfg["Kbc"]; shear = cfg["shear"]
    g = global_loads(props, L0, Kbc, shear, shear_rig)   # global loads

    # local plate modes -- RIGOROUS multi-harmonic finite strip:
    nloc = cfg["n_local"]                                # how many local modes asked
    # LOCAL PLATE BUCKLING (rigorous).  Each wall is solved by the multi-harmonic
    # finite strip that keeps the FULL laminate [[A,B],[B,D]], so D16/D26 and
    # B16/B26 are captured EXACTLY -- no orthotropic knockdown, no shear fudge.
    # End walls are outstands (one free edge); interior walls are simply supported
    # at both fold lines (the classic, design-standard local-buckling assumption).
    loc = rigorous_local_path(corners, wallABD_raw, wall_seg, Kbc*L0, nloc)
    mloc = None
    Ploc = loc * props["Lmid"]                           # N/mm -> N  (x midline)
    eta_loc, gov_w = 1.0, 0                              # coupling now inside the solver

    # ---- 6.6  build the mode list (4 global + nloc local) ------------------
    modes = {1: ("Axial extension", None, "loaded / prebuckling mode"),
             2: ("Major-axis flexural", g["P2"], "Euler bending, strong axis"),
             3: ("Minor-axis flexural", g["P3"], "Euler bending, weak axis"),
             4: ("Torsion",             g["P4"], "St-Venant + warping torsion")}
    for j in range(nloc):                                # append the local modes
        modes[5+j] = (f"Local plate-{j+1}", Ploc[j],
                      f"plate-buckling eigenvalue #{j+1}")
    g_FT = g["PFT"]                                      # coupled flex-torsional

    # ---- 6.7  governing mode (lowest finite load incl. FT coupling) --------
    candidates = {"3+4 Flexural-torsional": g_FT, "5 Local (lowest)": Ploc[0]}
    gov_name = min(candidates, key=candidates.get)       # classical envelope ...
    gov_val = candidates[gov_name]
    # ---- 6.7b  FULL GBT (paper engine) -> the DEFINITIVE coupled load ------
    #  All conventional + shear modes and every A16/B16/D16/D26 coupling
    #  interact in one eigenproblem; this replaces the classical envelope as
    #  the governing value (the mode table above remains the interpretation).
    # ---- 6.7b  GOVERNING LOAD: unified shear-deformable spline finite strip
    #  (section 20): the whole member -- every wall a laminated Mindlin plate
    #  (full A, B, D, As), folds with full continuity, cubic splines along the
    #  length, the actual end conditions imposed node by node, and the
    #  pre-buckling state from the member's own static solution.  Local,
    #  distortional and global buckling and their interaction come out of one
    #  eigenproblem.  The classical envelope above is kept for interpretation.
    P_grip = None; grip_valid = False
    P_gbt, gbt_info, _ur = _unified_governing(cfg, L0)
    gov_name, gov_val = gbt_info["kind"], P_gbt
    mloc = _ur.get("m")

    n_natural = len(corners)                             # natural nodes
    n_total = len(path)                                  # total nodes
    n_walls = n_natural - 1                              # number of walls
    m_inter = n_total - n_natural                        # intermediate nodes = m
    n_global = min(4, n_natural)                         # independent global modes
    n_distort = max(n_natural - 4, 0)                    # distortional modes
    n_local_available = m_inter                          # local modes = m
    n_conv_total = n_global + n_distort + n_local_available   # = n_total (identity)

    return dict(path=path, s_path=s_path, s_corner=s_corner, wall_seg=wall_seg,
                props=props, modes=modes, PFT=g_FT, part=g["part"], r0sq=g["r0sq"],
                Ploc=Ploc, mloc=mloc, gov_name=gov_name, gov_val=gov_val,
                a11_seg=a11_seg, a66_seg=a66_seg, dmats=dmats,
                wallABD=wallABD, shear_rig=shear_rig, L0=L0, Kbc=Kbc, shear=shear,
                n_natural=n_natural, n_total=n_total, n_walls=n_walls, m_inter=m_inter,
                n_global=n_global, n_distort=n_distort, n_conv_total=n_conv_total,
                n_local_available=n_local_available, nloc=nloc, eta_loc=eta_loc,
                gov_wall=gov_w, wallABD_raw=wallABD_raw,
                P_gbt=P_gbt, gbt_info=gbt_info, P_grip=P_grip,
                grip_valid=grip_valid, _unified=_ur)


# ==============================================================================
# 7.  LENGTH SWEEP  ->  P_cr(L) for every mode
# ==============================================================================
def sweep_length(cfg, res):
    """P_cr(L) for all modes across [Lmin,Lmax] step dL."""
    Ls = np.arange(cfg["Lmin"], cfg["Lmax"] + 1e-9, cfg["dL"])   # length grid
    props = res["props"]; Kbc = cfg["Kbc"]; shear = cfg["shear"]
    sr = res["shear_rig"]
    P2=[];P3=[];P4=[];PFT=[];PLOC=[]                     # one list per curve
    graph = res.get("is_graph", False)
    # LOCAL buckling is length-independent -> compute the plateau ONCE (fast sweep)
    if graph:
        Ploc_flat = local_buckling_g(res["graph"], res["cells"],
                                     Kbc*props["Lmid"], n_modes=1)[0]*props["Lmid"]
    else:
        Ploc_flat = rigorous_local_path(cfg["corner_nodes"],
                                        res.get("wallABD_raw", res["wallABD"]),
                                        res["wall_seg"], 0.0, 1)[0]*props["Lmid"]
    for L in Ls:                                         # classical curves
        g = global_loads(props, L, Kbc, shear, sr)
        P2.append(g["P2"]); P3.append(g["P3"]); P4.append(g["P4"]); PFT.append(g["PFT"])
        PLOC.append(Ploc_flat)
    # FULL model curve (same shear-deformable strips, longitudinal harmonics,
    # interior pre-buckling state).  The design-length value with the actual
    # end zones is res["gov_val"] (marked on the plot).
    print("\n  P_cr(L) curve of the full model:")
    PGBT = _unified_sweep(cfg, Ls); PCLP = PSHR = None
    return (Ls, np.array(P2), np.array(P3), np.array(P4),
            np.array(PFT), np.array(PLOC),
            np.array(PGBT) if PGBT is not None else None,
            np.array(PCLP) if PCLP is not None else None,
            np.array(PSHR) if PSHR is not None else None)


# ==============================================================================
# 8.  PLOTS  (all the graphs from the brief, incl. deformed mode shapes)
# ==============================================================================


def make_plots(cfg, res, outdir="out", tag="C"):
    import os; os.makedirs(outdir, exist_ok=True)
    kN = 1e3
    Ls, P2, P3, P4, PFT, PLOC, PGBT, PCLP, PSHR = sweep_length(cfg, res)
    env = PGBT if PGBT is not None else np.minimum(PFT, PLOC)

    # 1. Original standard plots
    if res.get("is_graph"):
        _plots_graph(cfg, res, Ls, P2, P3, P4, PFT, PLOC, env, outdir, tag)
    elif cfg.get("corner_nodes") is not None:            # family A (c, l, z, lipc): same figures
        _lay = list(cfg.get("layup", [0.0]))
        _cn = [tuple(map(float, p)) for p in cfg["corner_nodes"]]
        _cfgA = dict(cfg, nodes=_cn,
                     walls=[dict(i=w, j=w + 1,
                                 plies=[(float(a), float(cfg["t_wall"][w]) / len(_lay))
                                        for a in _lay])
                            for w in range(len(_cn) - 1)])
        _plots_graph(_cfgA, res, Ls, P2, P3, P4, PFT, PLOC, env, outdir, tag)

    # 2. STACKED PLOT (P_cr + participation)
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 9), sharex=True)
        ax1.plot(Ls, P2/kN, label="Mode 2 (Major flexural)", color="#3498db")
        ax1.plot(Ls, P3/kN, label="Mode 3 (Minor flexural)", color="#2ecc71")
        ax1.plot(Ls, P4/kN, label="Mode 4 (Torsion)", color="#f39c12")
        ax1.plot(Ls, PFT/kN, "r--", lw=2, label="Modes 3+4 (Flexural-Torsional)")
        ax1.plot(Ls, PLOC/kN, color="#9b59b6", label="Modes 5+ (Local plate)")
        ax1.plot(Ls, env/kN, "k", lw=2.5,
                 label="FULL model (" + ("GBT" if _gbt_opts(cfg)[0] == "gbt" else "strip") + ")")
        ax1.plot([cfg["L_design"]], [res["gov_val"]/kN], "k*", ms=13,
                 label="P_cr at design length (actual ends)")
        ax1.set_yscale("log"); ax1.set_ylabel(r"$P_{cr}$ [kN]", fontsize=12)
        ax1.set_title(f"GBT Buckling Analysis ({tag}-section)", fontsize=13, fontweight='bold')
        ax1.grid(True, which="both", ls="--", alpha=0.5)
        ax1.legend(fontsize=9, loc="upper right")

        Lp = np.linspace(cfg["Lmin"], cfg["Lmax"], 240)
        P2p = np.zeros_like(Lp); P3p = np.zeros_like(Lp); P4p = np.zeros_like(Lp)
        PLp = np.zeros_like(Lp); PFp = np.zeros_like(Lp)
        for i, L in enumerate(Lp):
            gg = global_loads(res["props"], L, cfg["Kbc"], cfg["shear"], res["shear_rig"])
            P2p[i], P3p[i], P4p[i], PFp[i] = gg["P2"], gg["P3"], gg["P4"], gg["PFT"]
            PLp[i] = np.interp(L, Ls, PLOC)
        modes = {"Local plate": PLp, "Flexural-torsional": PFp, "Minor flexural": P3p, "Major flexural": P2p}
        inv = {k: 1.0/np.maximum(v, 1e-9) for k, v in modes.items()}
        tot = sum(inv.values())
        part = {k: inv[k]/tot for k in modes}
        cols = {"Local plate": "#9b59b6", "Flexural-torsional": "#e74c3c", "Minor flexural": "#2ecc71", "Major flexural": "#3498db"}
        base = np.zeros_like(Lp)
        for k in modes:
            ax2.fill_between(Lp, base, base+part[k], color=cols[k], alpha=0.6, label=k)
            ax2.plot(Lp, base+part[k], color="white", lw=0.5)
            base = base+part[k]
        ax2.set_xlim(Lp.min(), Lp.max()); ax2.set_ylim(0, 1.0)
        ax2.set_xlabel("Length L [mm]", fontsize=12)
        ax2.set_ylabel("indicative share (1/P of the classical loads)", fontsize=10)
        ax2.grid(True, ls="--", alpha=0.5); ax2.legend(fontsize=9, loc="lower right")
        plt.tight_layout(); plt.savefig(f"{outdir}/fig_stacked_{tag}.png", dpi=300); plt.close()
    except Exception as e:
        print(f"Error plotting stacked graph: {e}")

    # 3. AMPLITUDE & STRESS PLOTS (all sections)
    try:
        L0 = cfg["L_design"]; x = np.linspace(0, L0, 400)
        plt.figure(figsize=(8.0, 4.6))
        if cfg["bc"] == "CC":
            shp = lambda m: 1-np.cos(2*np.pi*m*x/L0)
            ttl = "Illustrative amplitude functions (CC ends: 1-cos)"
        else:
            shp = lambda m: np.sin(np.pi*m*x/L0)
            ttl = "Illustrative amplitude functions (SS ends: sin)"
        for m in (1, 2): plt.plot(x, shp(m), label=f"m={m}")
        if res.get("mloc"):                                  # strip solver: half-waves known
            mg = max(int(res["mloc"]), 1)
            plt.plot(x, shp(mg), color="green", label=f"m={mg} (governing mode, strip)")
        plt.xlabel("x  [mm]"); plt.ylabel(r"amplitude  $\phi(x)$"); plt.title(ttl)
        plt.grid(True, alpha=.3); plt.legend(fontsize=8); plt.tight_layout()
        plt.savefig(f"{outdir}/fig2_amplitude_{tag}.png", dpi=130); plt.close()

        t_ref = res["props"]["_t_wall"][cfg.get("stress_wall", 0)]
        sig = env/(t_ref*res["props"]["Lmid"])
        plt.figure(figsize=(7.2, 4.6))
        plt.plot(Ls, sig, "k", lw=2); plt.xlabel("L [mm]"); plt.ylabel(r"$\sigma_{cr}$ [MPa]")
        plt.title(r"Average critical stress $P_{cr}\,/\,(t\cdot$mid-line length$)$ vs length")
        plt.grid(True, alpha=.3); plt.tight_layout()
        plt.savefig(f"{outdir}/fig4_sigma_vs_L_{tag}.png", dpi=130); plt.close()
    except Exception as e:
        pass

    try:
        _plot_unified_mode(cfg, res, outdir, tag)          # strip mode
        _plot_gbt_mode(cfg, res, outdir, tag)              # GBT mode
        _plot_gbt_participation(cfg, outdir, tag)          # participation vs L
    except Exception as e:
        print(f"critical-mode plot skipped: {e}")

    # 4. INTERACTIVE MODE SHAPES (Will appear for ALL sections)
    import sys
    # Skip interactive prompt ONLY if running a batch/parametric study
    if not any(arg in sys.argv for arg in ["--demo", "--cases", "--parametric", "--parametric-t", "--parametric-i", "--parametric-z", "--parametric-c", "--parametric-l", "--parametric-sections", "--parametric-rhs", "--paper", "--badcases"]):
        print("\n" + "="*50)
        print(" MODE SHAPE PLOTTING")
        print(" GBT cross-section deformation modes, numbered as in fig6_modeshapes")
        print("="*50)
        try:
            m1_str = input("Enter the FIRST mode number you want to plot (e.g., 1): ").strip()
            m2_str = input("Enter the SECOND mode number you want to plot (e.g., 2): ").strip()
            idx1 = max(0, int(m1_str) - 1) if m1_str else 0
            idx2 = max(0, int(m2_str) - 1) if m2_str else 1
            
            if res.get("is_graph"):
                material = cfg["material"]; nodes = np.asarray(cfg["nodes"], float); walls = cfg["walls"]
                def wallC5(w):
                    pl = w["plies"]
                    if material.get("type") == "isotropic": return material_C5_iso(material["E"], material["nu"], sum(t for _a, t in pl))
                    return material_C5_layup(material["E1"], material["E2"], material["G12"], material["nu12"], [a for a, _t in pl], sum(t for _a, t in pl), rbs=False)
                ws = [(w["i"], w["j"], wallC5(w), max(3, int(round(float(np.hypot(*(nodes[w["j"]]-nodes[w["i"]]))) / 8.0)))) for w in walls]
                sec = Section([tuple(map(float, p)) for p in nodes], ws)
            else:
                mat = cfg["material"]
                def wC5(w):
                    if mat["type"] == "isotropic": return material_C5_iso(mat["E"], mat["nu"], cfg["t_wall"][w])
                    return material_C5_layup(mat["E1"], mat["E2"], mat["G12"], mat["nu12"], list(cfg["layup"]), cfg["t_wall"][w], rbs=True)
                ws = [(w, w+1, wC5(w), int(cfg["n_int_per_wall"][w])) for w in range(len(cfg["corner_nodes"])-1)]
                sec = Section([tuple(map(float, p)) for p in cfg["corner_nodes"]], ws)
            
            ms = ModeSet(sec); npt = ms.npt
            XY = np.zeros((npt, 2)); TAN = np.zeros((npt, 2)); NRM = np.zeros((npt, 2)); WALL = np.zeros(npt, int); SPOS = np.zeros(npt)
            secnodes = np.asarray(sec.nodes, float)
            for ip, (isg, xi, wgt) in enumerate(ms.pts):
                iw, a, b, ls = sec.segs[isg]; wl = sec.walls[iw]
                pa = secnodes[wl.na]; pb = secnodes[wl.nb]
                sfrac = (wl.chain.index(a) + xi) / (wl.n_int + 1)
                XY[ip] = pa + (pb - pa) * sfrac; t = (pb - pa) / (np.linalg.norm(pb - pa) + 1e-9)
                TAN[ip] = t; NRM[ip] = np.array([-t[1], t[0]])
                WALL[ip] = iw; SPOS[ip] = sfrac
                
            fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
            fig.suptitle(f"GBT cross-section deformation modes ({tag}-section)", fontsize=13, fontweight='bold')
            axes = np.atleast_1d(axes).ravel()
            amp = 0.22 * max(np.ptp(secnodes[:, 0]), np.ptp(secnodes[:, 1]))
            
            for i, k in enumerate([idx1, idx2]):
                ax = axes[i]
                if k >= ms.M:
                    ax.text(0.5, 0.5, f"Mode {k+1} not found", ha='center'); ax.axis("off"); continue
                v = ms.V[:, k]; w = ms.W[:, k]
                disp = (v[:, None] * TAN + w[:, None] * NRM) * (amp / max(np.max(np.abs(v)), np.max(np.abs(w)), 1e-9))
                Pd = XY + disp
                for iw in range(len(sec.walls)):
                    m_ = WALL == iw
                    if m_.sum() < 2: continue
                    o = np.argsort(SPOS[m_]); xw = XY[m_][o]; pw = Pd[m_][o]
                    ax.plot(xw[:, 0], xw[:, 1], color="0.75", lw=1.6, zorder=1)
                    ax.plot(pw[:, 0], pw[:, 1], color=_GBT_COLORS[k % 7], lw=2.0, zorder=2)
                lab = ms.conv_label[k] if k < len(ms.conv_label) else f"mode {k+1}"
                ax.set_title(f"Mode {k+1}: {lab}", fontsize=11); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
                for sp in ax.spines.values(): sp.set_edgecolor("0.8")
            plt.tight_layout(); plt.savefig(f"{outdir}/fig_modes_interactive_{tag}.png", dpi=300); plt.show(); plt.close()
        except Exception as e:
            print(f"Skipping mode plotting: {e}")
            
    return Ls, env


def clt_plies(material, plies):
    """Laminate A,B,D,ttot for a wall from an explicit ply stack [(ang,thk),...]."""
    if material["type"] == "isotropic":
        E, nu = material["E"], material["nu"]; E1 = E2 = E; G12 = E/(2*(1+nu)); nu12 = nu
    else:
        E1, E2 = material["E1"], material["E2"]; G12, nu12 = material["G12"], material["nu12"]
    nu21 = nu12*E2/E1; Dl = 1-nu12*nu21
    Q11, Q22, Q12, Q66 = E1/Dl, E2/Dl, nu12*E2/Dl, G12
    def Qbar(thd):
        th = np.radians(thd); c, s = np.cos(th), np.sin(th)
        return np.array([
          [Q11*c**4+2*(Q12+2*Q66)*s**2*c**2+Q22*s**4,
           (Q11+Q22-4*Q66)*s**2*c**2+Q12*(c**4+s**4),
           (Q11-Q12-2*Q66)*s*c**3+(Q12-Q22+2*Q66)*s**3*c],
          [(Q11+Q22-4*Q66)*s**2*c**2+Q12*(c**4+s**4),
           Q11*s**4+2*(Q12+2*Q66)*s**2*c**2+Q22*c**4,
           (Q11-Q12-2*Q66)*s**3*c+(Q12-Q22+2*Q66)*s*c**3],
          [(Q11-Q12-2*Q66)*s*c**3+(Q12-Q22+2*Q66)*s**3*c,
           (Q11-Q12-2*Q66)*s**3*c+(Q12-Q22+2*Q66)*s*c**3,
           (Q11+Q22-2*Q12)*s**2*c**2+Q66*(c**2-s**2)**2]])
    ttot = sum(t for _, t in plies); z = -ttot/2
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
    for ang, t in plies:
        Qb = Qbar(ang); z0, z1 = z, z+t
        A += Qb*(z1-z0); B += Qb*(z1**2-z0**2)/2; D += Qb*(z1**3-z0**3)/3; z = z1
    return A, B, D, ttot


def make_section(shape, dims, material, layup_b, layup_s, tp_b, tp_s):
    """Return (nodes, walls, cells, family). tp_b, tp_s = per-ply thickness of
    base/stiffener; a wall's total thickness = (num plies) x (ply thickness)."""
    pb = [(a, tp_b) for a in layup_b]                    # base ply stack
    ps = [(a, tp_s) for a in layup_s]                    # stiffener ply stack
    pf = pb + ps                                         # foot = base + stiffener
    a = dims.get("a"); b = dims.get("b")
    def W(i, j, pl, part): return dict(i=i, j=j, plies=pl, part=part)
    if shape == "c":
        n = [(a, -b/2), (0, -b/2), (0, b/2), (a, b/2)]
        return n, [W(0,1,pb,"flange"), W(1,2,pb,"web"), W(2,3,pb,"flange")], [], "A"
    if shape == "l":
        return [(a,0),(0,0),(0,b)], [W(0,1,pb,"leg"), W(1,2,pb,"leg")], [], "A"
    if shape == "z":
        n = [(-a,-b/2),(0,-b/2),(0,b/2),(a,b/2)]
        return n, [W(0,1,pb,"flange"), W(1,2,pb,"web"), W(2,3,pb,"flange")], [], "A"
    if shape == "t":                                     # branch at node 1
        n = [(-a/2,b),(0,b),(a/2,b),(0,0)]
        return n, [W(0,1,pb,"flange"), W(1,2,pb,"flange"), W(1,3,pb,"web")], [], "B"
    if shape == "i":                                     # branches at 1 and 4
        n = [(-a/2,b),(0,b),(a/2,b),(-a/2,0),(0,0),(a/2,0)]
        w = [W(0,1,pb,"top"), W(1,2,pb,"top"), W(1,4,pb,"web"),
             W(3,4,pb,"bot"), W(4,5,pb,"bot")]
        return n, w, [], "B"
    if shape == "box":                                   # one closed cell
        n = [(-a/2,-b/2),(a/2,-b/2),(a/2,b/2),(-a/2,b/2)]
        return n, [W(k,(k+1)%4,pb,"wall") for k in range(4)], [[0,1,2,3]], "C"
    if shape == "hat":                                   # closed + branched (nodes 2,3)
        a1, a2 = dims["a1"], dims["a2"]
        n = [(-a/2,0),(-a2/2,0),(-a1/2,b),(a1/2,b),(a2/2,0),(a/2,0)]
        w = [W(0,1,pf,"foot"), W(1,2,ps,"leg"), W(2,3,ps,"top"),
             W(3,4,ps,"leg"), W(4,5,pf,"foot"), W(1,4,pb,"inner-base")]
        return n, w, [[1,2,3,4]], "D"
    if shape == "dhat":                                  # double hat: 2 cells, 8 natural
        a1, a2 = dims["a1"], dims["a2"]
        def hn(cx): return [(cx-a/2,0),(cx-a2/2,0),(cx-a1/2,b),(cx+a1/2,b),(cx+a2/2,0),(cx+a/2,0)]
        H1 = hn(-a/2); H2 = hn(+a/2)                     # share the middle foot tip
        n = H1 + H2[1:]; o = 5                           # H2 node k>=1 -> index o+k
        w = [W(0,1,pf,"foot"), W(1,2,ps,"leg"), W(2,3,ps,"top"), W(3,4,ps,"leg"),
             W(4,5,pf,"foot"), W(1,4,pb,"inner-base"),
             W(5,o+1,pf,"foot"), W(o+1,o+2,ps,"leg"), W(o+2,o+3,ps,"top"),
             W(o+3,o+4,ps,"leg"), W(o+4,o+5,pf,"foot"), W(o+1,o+4,pb,"inner-base")]
        return n, w, [[1,2,3,4],[o+1,o+2,o+3,o+4]], "D"
    raise ValueError("unknown shape "+shape)


def _valence(walls):
    from collections import Counter
    c = Counter()
    for w in walls:
        c[w["i"]] += 1; c[w["j"]] += 1
    return c


def build_graph(nodes, walls, material, n_int):
    """Discretize every wall into segments carrying that wall's stiffness."""
    nodes = np.asarray(nodes, float)
    
    # 1. Calculate classical A, B, D matrices for each wall
    wABD_raw = [clt_plies(material, w["plies"]) for w in walls]
    
    # 2. Account for Asymmetric/Unsymmetric Laminates using Reduced Stiffness
    wABD = []
    for A, B, D, tw in wABD_raw:
        if np.max(np.abs(B)) > 1e-5:                     # Unsymmetric
            invA = np.linalg.inv(A)
            invD = np.linalg.inv(D)
            A_star = A - B @ invD @ B
            D_star = D - B @ invA @ B
            wABD.append((A_star, B, D_star, tw))
        else:                                            # Symmetric
            wABD.append((A, B, D, tw))

    wABD_raw2 = [(A, B, D, tw) for (A, B, D, tw) in wABD_raw]
    seg = []
    for wi, w in enumerate(walls):
        p0, p1 = nodes[w["i"]], nodes[w["j"]]; A, B, D, tw = wABD[wi]
        ns = n_int[wi] + 1
        for k in range(ns):
            q0 = p0+(p1-p0)*k/ns; q1 = p0+(p1-p0)*(k+1)/ns
            # Assigning the A11, A66, and D elements from the (potentially reduced) matrices
            seg.append(dict(wall=wi, p0=q0, p1=q1, A11=_beam_a11(A), A66=A[2,2],   # A11 := beam E t
                            t=tw, D=(D[0,0], D[0,1], D[1,1], D[2,2])))
    return dict(nodes=nodes, walls=walls, wABD=wABD, wABD_raw=wABD_raw2, seg=seg)


def section_props_g(G):
    """Modulus-weighted EA, centroid, EI (rigorous for every family)."""
    seg = G["seg"]; EA = SY = SZ = Lmid = 0.0
    for s in seg:
        mid = 0.5*(s["p0"]+s["p1"]); dl = np.linalg.norm(s["p1"]-s["p0"])
        EA += s["A11"]*dl; SY += s["A11"]*mid[0]*dl; SZ += s["A11"]*mid[1]*dl; Lmid += dl
    Yc, Zc = SY/EA, SZ/EA; EIyy = EIzz = EIyz = 0.0
    for s in seg:
        y0, z0 = s["p0"]-[Yc, Zc]; y1, z1 = s["p1"]-[Yc, Zc]
        dl = np.linalg.norm(s["p1"]-s["p0"]); a = s["A11"]
        EIyy += a*dl*(z0**2+z0*z1+z1**2)/3
        EIzz += a*dl*(y0**2+y0*y1+y1**2)/3
        EIyz += a*dl*(2*y0*z0+y0*z1+y1*z0+2*y1*z1)/6
    return dict(Lmid=Lmid, EA=EA, Yc=Yc, Zc=Zc, EIyy=EIyy, EIzz=EIzz, EIyz=EIyz)


def torsion_GJ_g(G, cells):
    """GJ = open St-Venant (Sum 4 D66 L) + closed-cell Bredt (4 Omega^2/contour)."""
    nodes = G["nodes"]; walls = G["walls"]; wABD = G["wABD"]
    cell_edges = set()
    for loop in cells:
        for a, c in zip(loop, loop[1:]+loop[:1]):
            cell_edges.add(frozenset((a, c)))
    GJo = 0.0
    for wi, w in enumerate(walls):
        L = np.linalg.norm(nodes[w["i"]]-nodes[w["j"]])
        if frozenset((w["i"], w["j"])) not in cell_edges:
            GJo += 4*wABD[wi][2][2, 2]*L
    GJc = 0.0
    for loop in cells:
        P = np.array([nodes[i] for i in loop]); x, y = P[:, 0], P[:, 1]
        Om = 0.5*abs(np.dot(x, np.roll(y, -1))-np.dot(y, np.roll(x, -1)))
        contour = 0.0
        for a, c in zip(loop, loop[1:]+loop[:1]):
            L = np.linalg.norm(nodes[a]-nodes[c])
            for wi, w in enumerate(walls):
                if frozenset((w["i"], w["j"])) == frozenset((a, c)):
                    contour += L/wABD[wi][0][2, 2]; break
        if contour > 0:
            GJc += 4*Om**2/contour
    return GJo+GJc, GJo, GJc


def warping_g(G, cells):
    """Sectorial warping + shear centre; exact for open, estimate (tagged) else."""
    nodes = G["nodes"]; walls = G["walls"]; wABD = G["wABD"]
    p = section_props_g(G); Yc, Zc, EA = p["Yc"], p["Zc"], p["EA"]
    EIyy, EIzz, EIyz = p["EIyy"], p["EIzz"], p["EIyz"]; N = len(nodes)
    adj = {k: [] for k in range(N)}
    for wi, w in enumerate(walls):
        adj[w["i"]].append(w["j"]); adj[w["j"]].append(w["i"])
    om = np.full(N, np.nan); om[0] = 0.0; q = [0]; seen = {0}
    while q:
        u = q.pop(0); Pu = nodes[u]-[Yc, Zc]
        for v in adj[u]:
            if v not in seen:
                Pv = nodes[v]-[Yc, Zc]
                om[v] = om[u]+(Pu[0]*Pv[1]-Pu[1]*Pv[0]); seen.add(v); q.append(v)
    A11w = [wABD[wi][0][0, 0] for wi in range(len(walls))]
    def sint(f):
        t = 0.0
        for wi, w in enumerate(walls):
            yi, zi = nodes[w["i"]]-[Yc, Zc]; yj, zj = nodes[w["j"]]-[Yc, Zc]
            L = np.linalg.norm(nodes[w["i"]]-nodes[w["j"]])
            t += f(A11w[wi], L, om[w["i"]], om[w["j"]], yi, zi, yj, zj)
        return t
    Aw = sint(lambda a, L, oi, oj, *_: a*L*(oi+oj)/2)/EA
    Iwy = sint(lambda a, L, oi, oj, yi, zi, yj, zj:
               a*L*(2*(oi-Aw)*yi+(oi-Aw)*yj+(oj-Aw)*yi+2*(oj-Aw)*yj)/6)
    Iwz = sint(lambda a, L, oi, oj, yi, zi, yj, zj:
               a*L*(2*(oi-Aw)*zi+(oi-Aw)*zj+(oj-Aw)*zi+2*(oj-Aw)*zj)/6)
    den = EIyy*EIzz-EIyz**2
    ys = (EIzz*Iwz-EIyz*Iwy)/den; zs = (-EIyy*Iwy+EIyz*Iwz)/den
    # --- symmetry enforcement: a section mirror-symmetric about an axis has its
    #     shear centre ON that axis (exact), which the open-tree estimate misses
    #     for closed cells.  Detect mirror symmetry of the node set and zero the
    #     appropriate shear-centre component.
    def mirror_sym(axis):                                # axis 0 => Y->-Y about Yc
        pts = nodes.copy().astype(float)
        c = Yc if axis == 0 else Zc
        refl = pts.copy(); refl[:, axis] = 2*c - refl[:, axis]
        for r in refl:                                   # every reflected node present?
            if np.min(np.sum((pts-r)**2, axis=1)) > 1e-6:
                return False
        return True
    if mirror_sym(0):                                    # symmetric about Z-axis (Y->-Y)
        ys = 0.0
    if mirror_sym(1):                                    # symmetric about Y-axis (Z->-Z)
        zs = 0.0
    omS = np.array([om[k]+zs*(nodes[k][0]-Yc)-ys*(nodes[k][1]-Zc) for k in range(N)])
    mS = 0.0
    for wi, w in enumerate(walls):
        L = np.linalg.norm(nodes[w["i"]]-nodes[w["j"]])
        mS += A11w[wi]*L*(omS[w["i"]]+omS[w["j"]])/2
    mS /= EA; ECw = 0.0
    for wi, w in enumerate(walls):
        L = np.linalg.norm(nodes[w["i"]]-nodes[w["j"]])
        o0, o1 = omS[w["i"]]-mS, omS[w["j"]]-mS
        ECw += A11w[wi]*L*(o0**2+o0*o1+o1**2)/3
    val = _valence(walls); branched = max(val.values()) > 2
    conf = "Definitive" if (not cells and not branched) else "Probable"
    return ys, zs, ECw, conf


def local_buckling_g(G, cells, L, n_modes=3):
    """Per-wall plate buckling; edge held (SS) at junction/cell nodes, free at tips."""
    nodes = G["nodes"]; walls = G["walls"]; wABD = G["wABD"]
    val = _valence(walls); cellnodes = set(n for lp in cells for n in lp)
    res = []
    for wi, w in enumerate(walls):
        b = np.linalg.norm(nodes[w["i"]]-nodes[w["j"]])
        D = wABD[wi][2]; D11, D12, D22, D66 = D[0, 0], D[0, 1], D[1, 1], D[2, 2]
        hi = (val[w["i"]] > 1) or (w["i"] in cellnodes)
        hj = (val[w["j"]] > 1) or (w["j"] in cellnodes)
        if hi and hj:
            Ncr = 2*np.pi**2/b**2*(np.sqrt(D11*D22)+D12+2*D66)   # SS-SS
        else:
            Ncr = 12*D66/b**2 + np.pi**2*D11/L**2                # free outstand
        res.append(Ncr)
    res = np.sort(np.array(res)); out = np.full(n_modes, res[-1])
    out[:min(n_modes, len(res))] = res[:n_modes]; return out


def shear_rig_g(G):
    """General shear rigidities S2 (EIyy mode), S3 (EIzz mode) from direction cosines."""
    S2 = S3 = 0.0
    for s in G["seg"]:
        dl = np.linalg.norm(s["p1"]-s["p0"])
        if dl == 0:
            continue
        dY = (s["p1"][0]-s["p0"][0])/dl; dZ = (s["p1"][1]-s["p0"][1])/dl
        S2 += s["A66"]*dZ**2*dl; S3 += s["A66"]*dY**2*dl
    return S2, S3


# ---- legacy open-path generators (family A still uses the validated path code)
def shape_C(a, b):
    return [(a, -b/2), (0, -b/2), (0, b/2), (a, b/2)], ["flange", "web", "flange"]
def shape_L(a, b):
    return [(a, 0.0), (0.0, 0.0), (0.0, b)], ["leg", "leg"]
def shape_Z(a, b):
    return [(-a, -b/2), (0, -b/2), (0, b/2), (a, b/2)], ["flange", "web", "flange"]
def shape_LIPC(a, b, c):
    """Lipped channel (family A): lip - flange - web - flange - lip (lips inward)."""
    return ([(a, -b/2+c), (a, -b/2), (0, -b/2), (0, b/2), (a, b/2), (a, b/2-c)],
            ["lip", "flange", "web", "flange", "lip"])
SHAPES = {"c": shape_C, "l": shape_L, "z": shape_Z, "lipc": shape_LIPC}


def get_config():
    print("="*78)
    print("  GBT BUCKLING ANALYSIS  -  general interactive base code")
    print("  (<Enter> accepts the default in [ ];  units: N, mm, MPa, kN)")
    print("  Defaults = the Abaqus models of the parametric studies, so <Enter>")
    print("  everywhere reproduces a parametric case exactly.")
    print("="*78)

    # ---- cross-section shape (first: the defaults below depend on it) -------
    print("\nCross-section shape (families):")
    print("  A open      : c=channel   l=angle   z=Z   lipc=lipped channel")
    print("  B branched  : t=T-section i=I-section")
    print("  C closed    : box=RHS/box")
    print("  D closed+br : hat=top-hat  dhat=double-hat   (stiffener on a base)")
    print("  E panel     : omega = omega-stiffened panel (Abaqus model; GBT and/or strip)")
    shp = ask("  shape", "c", str,
              ["c", "l", "z", "lipc", "t", "i", "box", "hat", "dhat", "omega"]).lower()
    OMEGA = shp == "omega"

    # ---- material elastic constants (defaults = the Abaqus model) -----------
    #   member sections : _STUDY_MATS (E1 171400, E2 9080, G12 = G13 5290, G23 2645)
    #   omega panel     : OMEGA_CFRP  (E1 150000, E2 9000, G12 = G13 = G23 5300)
    mdef = OMEGA_CFRP if OMEGA else _STUDY_MATS["cfrp"]
    idef = OMEGA_ISO if OMEGA else _STUDY_MATS["iso"]
    print("\nMaterial type:  1 = isotropic    2 = orthotropic (laminate)")
    mnum = ask("  choose 1 or 2", 2, int, ["1", "2"])
    if mnum == 1:
        E  = ask("  Young's modulus E [MPa]", idef["E"], float)
        nu = ask("  Poisson ratio nu", idef["nu"], float)
        material = {"type": "isotropic", "E": E, "nu": nu}
        iso = True
    else:
        E1   = ask("  E1  (fibre dir.) [MPa]", mdef["E1"], float)
        E2   = ask("  E2  (transverse) [MPa]", mdef["E2"], float)
        G12  = ask("  G12 (in-plane shear) [MPa]", mdef["G12"], float)
        # G13, G23 are the THROUGH-THICKNESS shear moduli (used by the
        # shear-deformable strips; defaults = the values of the Abaqus model).
        G13  = ask("  G13 (transverse shear, x-z) [MPa]", mdef["G13"], float)
        G23  = ask("  G23 (transverse shear, s-z) [MPa]", mdef["G23"], float)
        nu12 = ask("  nu12", mdef["nu12"], float)
        material = {"type": "orthotropic", "E1": E1, "E2": E2, "G12": G12,
                    "G13": G13, "G23": G23, "nu12": nu12}
        iso = False

    def ask_layup(name, default_ply=0.25, angles=(45, -45, -45, 45), t_iso=1.0):
        """Ask a lay-up (angles + ply thickness) and report the total thickness.
        Defaults: the Sym1 lay-up of the parametric studies."""
        if iso:
            t = ask(f"  {name} wall thickness [mm]", t_iso, float)
            print(f"    -> {name}: isotropic wall, t = {t:g} mm")
            return [0.0], t                              # single 'ply', angle irrelevant
        npl = ask(f"  {name}: number of plies", len(angles), int)
        dang = (list(angles) * npl)[:npl]
        ang = ask_list(f"  {name}: ply angles [deg]", dang, float)
        ang = (ang*npl)[:npl]
        tp  = ask(f"  {name}: ply thickness [mm]", default_ply, float)
        print(f"    -> {name}: {npl} plies x {tp:g} = total {npl*tp:g} mm")
        return ang, tp

    if OMEGA:                                            # own questions + solver choice
        return _get_omega_config(material, iso, ask_layup)

    HAT = shp in ("hat", "dhat")                         # needs stiffener + base
    GRAPH = shp in ("t", "i", "box", "hat", "dhat")      # families B/C/D

    # ---- dimensions ---------------------------------------------------------
    print("\nGeometry [mm]:")
    # defaults: the parametric-study sections (flange / width 20, height 50);
    # the hat stiffeners keep their own reference geometry
    a = ask("  a  (overall width / flange / base)", 50.0 if HAT else 20.0, float)
    b = ask("  b  (height)", 100.0 if HAT else 50.0, float)
    dims = {"a": a, "b": b}
    if HAT:                                              # hat needs a1 (top), a2 (inner base)
        dims["a1"] = ask("  a1 (top width, between top corners)", 20.0, float)
        dims["a2"] = ask("  a2 (inner base, between the leg feet)", 30.0, float)
    if shp == "lipc":                                    # lipped channel needs the lip
        dims["c"] = ask("  c  (lip length)", 10.0, float)

    # ---- lay-up(s) / wall thickness -----------------------------------------
    print("\nWall thickness:" if iso else "\nLaminate lay-up:")
    if HAT:                                              # two parts, may differ
        if not iso:
            print("  (stiffener = legs+top+feet-overlay ; base = the flat base plate)")
        layup_b, tp_b = ask_layup("BASE (tb)")
        layup_s, tp_s = ask_layup("STIFFENER (ts)")
        tb = len(layup_b)*tp_b if not iso else tp_b
        ts = len(layup_s)*tp_s if not iso else tp_s
        print(f"    -> feet (stiffener on base) total = ts + tb = {ts:g} + {tb:g} = {ts+tb:g} mm")
    else:
        part_name = "section" if iso else "laminate"
        layup_b, tp_b = ask_layup(part_name)
        layup_s, tp_s = layup_b, tp_b                    # single thickness/layup
    # ---- build the section --------------------------------------------------
    if GRAPH:                                            # families B/C/D
        nodes, walls, cells, family = make_section(shp, dims, material,
                                                    layup_b, layup_s, tp_b, tp_s)
        n_walls = len(walls)
        print(f"  -> {shp.upper()}: family {family}, {len(nodes)} natural nodes, "
              f"{n_walls} walls, {len(cells)} closed cell(s)")
        corner_nodes = None; t_wall = None; layup = layup_b
    else:                                                # family A: validated path
        family = "A"
        corner_nodes, _ = (shape_LIPC(a, b, dims["c"]) if shp == "lipc"
                           else SHAPES[shp](a, b))
        n_walls = len(corner_nodes)-1
        nodes = walls = cells = None
        layup = layup_b
        print(f"  -> {shp.upper()}: family A, {len(corner_nodes)} natural nodes, "
              f"{n_walls} walls")
        # per-wall thickness (family A keeps its variable-thickness feature)
        t_uniform = len(layup_b)*tp_b
        print(f"\nWall thickness per wall ({n_walls} values; default uniform {t_uniform:g} mm):")
        t_wall = ask_list("  thickness per wall [mm]", [t_uniform]*n_walls, float)
        t_wall = (t_wall+[t_wall[-1]]*n_walls)[:n_walls]

    # ---- intermediate nodes per wall ---------------------------------------
    # ---- intermediate nodes per wall (ULTRA-FINE MESH FOR HIGHEST ACCURACY) ----------
    if GRAPH:
        nd = np.asarray(nodes, float)
        _wl = [float(np.hypot(*(nd[w["j"]] - nd[w["i"]]))) for w in walls]
    else:
        cn = np.asarray(corner_nodes, float)
        _wl = [float(np.hypot(*(cn[w+1] - cn[w]))) for w in range(n_walls)]
    auto_nint = auto_nint_for_walls(_wl)          # unified auto-mesh (all shapes)

    nint = list(auto_nint)          # nodes of the classical reference part (not the solvers)

    # ---- boundary conditions, shear, method, sweep -------------------------
    bc = ask("\nBoundary conditions (SS/CC)", "SS", str, ["ss", "cc"]).upper()
    Kbc = {"SS": 1.0, "CC": 0.5}[bc]
    
    grip = False                                         # (end restraint is modelled exactly)
    sh = True                                            # Engesser shear in the reference loads

    # ---- through-thickness PLATE THEORY  (numbered choice) -----------------
    #  [1] CLPT : classical (Kirchhoff) -- the existing, validated engine.
    #  [2] FSDT : first-order shear deformation (constant gamma_xz, k=5/6).
    #  [3] TSDT : Reddy third-order (parabolic gamma_xz, no shear factor).
    #  FSDT/TSDT are solved by the shear-deformable FINITE STRIP (Section 17);
    #  CLPT is left completely untouched.
    # Plate theory of the governing solver: FSDT (Mindlin, shear factor 5/6),
    # the same kinematics as the Abaqus S4R shell.  The classical mode table
    # (Euler / torsion / plate loads) is still printed for interpretation.
    theory = "FSDT"; method = "spline-strip"
    print("\nGoverning solver (Abaqus end conditions: SS end edges held in U1, U2;"
          " CC rigid end plates):")
    print("  1 = GBT    shear-deformable GBT (Silvestre & Camotim 2013 + T, M, R modes),")
    print("             amplitude functions in cubic B-splines   (default; 144 Abaqus cases")
    print("             C, Z, L, T, I, RHS: mean |error| 1.07 %, max 6.7 %)")
    print("  2 = strip  shell-type Mindlin spline finite strip (reference; same 144 cases:")
    print("             mean |error| 0.62 %, max 7.7 %)")
    d_s = 2 if GOVERNING_SOLVER == "strip" else 1
    if shp in ("hat", "dhat") and d_s != 2:
        d_s = 2
        print("  (hat / double-hat: default 2.  GBT works for these sections but needs about")
        print("   2 GB (hat) / 5 GB (double-hat) of memory and 2-8 min per length, and it has")
        print("   not been compared with Abaqus for them)")
    snum = ask("  choose 1 or 2", d_s, int, ["1", "2"])
    solver = "strip" if snum == 2 else "gbt"
    gbt_method = "spline"
    gbt_fam = GBT_FAMILIES
    if solver == "gbt":
        print("  GBT mode families:  CS = conventional + shear (the paper);  add")
        print("    R = nodal rotations, M = transverse shear (Mindlin), T = transverse extension")
        gbt_fam = ask("  mode families", GBT_FAMILIES, str,
                      ["cs", "csr", "csm", "csmr", "cstm", "cstmr"]).upper()
    else:
        print("\nMesh of the spline strip:")
        print("  1 = standard  strips min(2.5, 0.025 h), spans min(5, L/40)"
              "  (converged; mean error vs Abaqus 0.56 %)")
        print("  2 = fast      strips 2.5 mm, spans 10 mm"
              "                  (~3-4x faster; 0.93 %)")
        print("  3 = fine      strips 0.625 mm, spans L/80"
              "                (~2x slower; 0.51 %)")
        mch = ask("  choose 1, 2 or 3", 1, int, ["1", "2", "3"])
        globals()["MESH_STRIP"], globals()["MESH_SPAN"] = {1: (None, None), 2: (2.5, 10.0),
                                                           3: (0.625, -80.0)}[mch]
    nloc = max(1, sum(nint))                             # every classical local mode
    swall = n_walls // 2                                 # wall of the stress plot (fig4)
    print("\nLength range for the P_cr(L) sweep:")
    Lmin = ask("  Lmin [mm]", 50.0, float); Lmax = ask("  Lmax [mm]", 2000.0, float)
    dL   = ask("  dL   [mm]", 50.0, float); L0 = ask("  design length L [mm]", 1000.0, float)

    return dict(material=material, layup=layup, t_wall=t_wall, a=a, b=b,
                corner_nodes=corner_nodes, nodes=nodes, walls=walls, cells=cells,
                family=family, shape=shp, n_walls=n_walls, n_int_per_wall=nint,
                Kbc=Kbc, bc=bc, grip=grip, shear=sh, theory=theory, method=method,
                n_local=nloc, stress_wall=swall,
                Lmin=Lmin, Lmax=Lmax, dL=dL, L_design=L0,
                solver=solver, gbt_method=gbt_method, gbt_families=gbt_fam)


# ==============================================================================
# 10b.  RESTORED HELPERS  (report, graph plots, mode text, small utilities)
# ==============================================================================
from collections import Counter                       # diagnostics


def cast(x, *a, **k):
    """No-op cast kept for compatibility."""
    return x


MODE_TEXT = """
------------------------------------------------------------------------------
  THE GBT DEFORMATION MODES  (numbering of fig6_modeshapes)
------------------------------------------------------------------------------
  mode 1   axial        : uniform axial extension -- the pre-buckling (loaded)
                          state; it has no buckling load of its own.
  warp@r   (modes 2..)  : unit warping at natural node r (wall end / fold line),
                          with the in-plane motion it implies.  Together they
                          give the GLOBAL modes (flexure about both axes,
                          torsion) and the DISTORTIONAL modes (fold lines move,
                          the cross-section changes shape).
  local@q  (the rest)   : unit transverse displacement of node q (intermediate
                          node or free tip): the LOCAL modes, each wall bends
                          like a plate strip while the fold lines stay put.
  The solver combines all of them (plus shear, transverse-extension and
  Mindlin modes); the RESULT block reports the share of each family.

  WHICH ONE GOVERNS:
    Long members  -> global; the load drops quickly with length and depends
                     strongly on the end condition (SS / CC differ a lot).
    Short members -> local plate buckling; the load is nearly independent of
                     length AND of the end condition, because the walls ripple
                     in short half-waves far away from the supports.  SS and CC
                     giving almost the same value there is EXPECTED, not a
                     bug -- Abaqus shows the same behaviour.
------------------------------------------------------------------------------
"""


def report(cfg, res):
    """Concise results report: governing load, reference loads, properties."""
    kN = 1e3
    p = res["props"]
    line = "=" * 68
    gname = res.get("gov_name", "")
    is_local = "local plate" in gname

    print("\n" + line)
    print(f"  RESULT   L={cfg['L_design']:.0f} mm, BC={cfg['bc']}")
    print(line)
    print(f"  P_cr = {res['gov_val']/kN:.4f} kN   ({res['gov_val']:.1f} N)")
    print(f"  solver: {gname}")
    gi = res.get("gbt_info") or {}
    if "participation" in gi:
        print("  GBT modal participation of the critical mode:")
        for g, v in sorted(gi["participation"].items(), key=lambda t: -t[1]):
            print(f"      {_GBT_GROUP_NAMES.get(g, g):<28} {100*v:6.1f} %")
    elif "rigid_share" in gi:
        print(f"  critical mode: {gi['mode']}  (rigid cross-section share "
              f"{100*gi['rigid_share']:.0f}%),  {gi['m']} half-wave(s)")
    if len(gi.get("eigs", [])) > 1:
        print(f"  next eigenvalue: {gi['eigs'][1]/kN:.4f} kN   ({gi.get('ndof', '?')} DOF)")

    print("\n" + line)
    print("  CLASSICAL REFERENCE LOADS  (rough uncoupled estimates, shown for")
    print("  mode interpretation only -- they are NOT the reported result)")
    print(line)
    md = res.get("modes", {})
    for k, lbl in ((2, "major-axis flexural"), (3, "minor-axis flexural"),
                   (4, "torsion")):
        if k in md and md[k][1]:
            print(f"    {lbl:<28} ~ {md[k][1]/kN:10.4f} kN")
    if res.get("PFT"):
        print(f"    {'flexural-torsional (coupled)':<28} ~ {res['PFT']/kN:10.4f} kN")
    if len(res.get("Ploc", [])):
        print(f"    {'local plate (per wall)':<28} ~ {res['Ploc'][0]/kN:10.4f} kN")
    print(f"  >> REPORTED P_cr                = {res['gov_val']/kN:10.4f} kN <<")

    print("\n" + line); print("  SECTION PROPERTIES"); print(line)
    print(f"  EA   = {p['EA']:.4g} N")
    print(f"  EIyy = {p['EIyy']:.4g}   EIzz = {p['EIzz']:.4g} N*mm^2")
    ecw = p['ECw'] if abs(p['ECw']) > 1e-9 * max(p['EIyy'], p['EIzz']) else 0.0   # round-off -> 0
    print(f"  GJ   = {p['GJ']:.4g}   ECw = {ecw:.4g}" + ("   (all walls meet at one point: no warping)" if ecw == 0.0 else ""))
    print(f"  centroid (Yc,Zc) = ({p['Yc']:.2f}, {p['Zc']:.2f}) mm")

    try:
        print("\n" + line); print("  LAMINATE (per wall)"); print(line)
        for w, (A, B, D) in enumerate(res["wallABD"]):
            tag = ""
            if abs(A[0, 2]) > 1e-6:
                tag += " unbalanced(A16)"
            if np.abs(B).max() > 1e-6:
                tag += " unsymmetric(B)"
            print(f"  wall {w}: A11={A[0,0]:8.0f}  A16={A[0,2]:8.0f}  "
                  f"D16={D[0,2]:7.1f}{tag}")
    except Exception:
        pass


def _plot_section_modes(cfg, res, outdir, tag, nshow=12):
    """Cross-section deformation modes, drawn as in the GBT papers."""
    material = cfg["material"]; nodes = np.asarray(cfg["nodes"], float)
    walls = cfg["walls"]; lay = cfg.get("layup", [0])

    def wallC5(w):
        pl = w["plies"]
        if material.get("type") == "isotropic":
            return material_C5_iso(material["E"], material["nu"],
                                   sum(t for _a, t in pl))
        ang = [a for a, _t in pl]; tt = sum(t for _a, t in pl)
        return material_C5_layup(material["E1"], material["E2"], material["G12"],
                                 material["nu12"], ang, tt, rbs=False)

    ws = []
    for w in walls:
        wl = float(np.hypot(*(nodes[w["j"]] - nodes[w["i"]])))
        ws.append((w["i"], w["j"], wallC5(w), max(3, int(round(wl / 6.0)))))
    sec = Section([tuple(map(float, p)) for p in nodes], ws)
    ms = ModeSet(sec); npt = ms.npt

    XY = np.zeros((npt, 2)); TAN = np.zeros((npt, 2)); NRM = np.zeros((npt, 2))
    WALL = np.zeros(npt, int); SPOS = np.zeros(npt)
    sn = np.asarray(sec.nodes, float)
    for ip, (isg, xi, wgt) in enumerate(ms.pts):
        iw, a, b, ls = sec.segs[isg]; wl = sec.walls[iw]
        pa = sn[wl.na]; pb = sn[wl.nb]
        idx = wl.chain.index(a); sf = (idx + xi) / (wl.n_int + 1)
        XY[ip] = pa + (pb - pa) * sf
        t = (pb - pa) / (np.linalg.norm(pb - pa) + 1e-9)
        TAN[ip] = t; NRM[ip] = np.array([-t[1], t[0]])
        WALL[ip] = iw; SPOS[ip] = sf

    span = max(np.ptp(sn[:, 0]), np.ptp(sn[:, 1])); amp = 0.22 * span
    labels = ms.conv_label; M = min(nshow, ms.M)
    ncol = 4; nrow = int(np.ceil(M / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.7 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for k in range(M):
        ax = axes[k]; v = ms.V[:, k]; w_ = ms.W[:, k]
        nv = max(np.max(np.abs(v)), np.max(np.abs(w_)), 1e-9)
        Pd = XY + (v[:, None] * TAN + w_[:, None] * NRM) * (amp / nv)
        for iw in range(len(sec.walls)):
            msk = WALL == iw
            if msk.sum() < 2:
                continue
            o = np.argsort(SPOS[msk])
            xw = XY[msk][o]; pw = Pd[msk][o]
            ax.plot(xw[:, 0], xw[:, 1], color="0.75", lw=1.6, zorder=1)
            ax.plot(pw[:, 0], pw[:, 1], color=_GBT_COLORS[k % len(_GBT_COLORS)],
                    lw=2.0, zorder=2)
        ax.set_title(f"mode {k+1}: {labels[k] if k < len(labels) else ''}",
                     fontsize=9)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for k in range(M, len(axes)):
        axes[k].axis("off")
    fig.suptitle(f"GBT cross-section deformation modes  ({tag})", fontsize=13)
    plt.tight_layout()
    plt.savefig(f"{outdir}/fig6_modeshapes_{tag}.png"); plt.close()


def _plots_graph(cfg, res, Ls, P2, P3, P4, PFT, PLOC, env, outdir, tag):
    """Figures for branched / closed sections: section map + P_cr(L) curve."""
    nodes = np.asarray(cfg["nodes"], float)
    plt.figure(figsize=(5.2, 5.2))
    for w in cfg["walls"]:
        a, b = nodes[w["i"]], nodes[w["j"]]
        plt.plot([a[0], b[0]], [a[1], b[1]], "-", lw=2.4, color="#0072B2")
    plt.plot(nodes[:, 0], nodes[:, 1], "o", ms=6, color="#D55E00", zorder=3)
    for i, n in enumerate(nodes):
        plt.annotate(str(i), (n[0], n[1]), textcoords="offset points",
                     xytext=(5, 5), fontsize=9)
    plt.gca().set_aspect("equal"); plt.grid(True, alpha=.3)
    plt.xlabel("y [mm]"); plt.ylabel("z [mm]")
    plt.title(f"Cross-section ({tag})")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig5_section_{tag}.png"); plt.close()

    plt.figure(figsize=(7.2, 4.6))
    if env is not None and len(env):
        plt.plot(Ls, np.asarray(env) / 1e3, "-", lw=2.6, color="k",
                 label="governing")
    if PLOC is not None and len(PLOC):
        plt.plot(Ls, np.asarray(PLOC) / 1e3, "--", lw=1.6, color="#CC79A7",
                 label="local plate")
    if PFT is not None and len(PFT):
        plt.plot(Ls, np.asarray(PFT) / 1e3, "--", lw=1.6, color="#D55E00",
                 label="flexural-torsional")
    plt.yscale("log"); plt.xlabel("member length L [mm]")
    plt.ylabel("$P_{cr}$ [kN]"); plt.title(f"Critical load vs length ({tag})")
    plt.grid(True, which="both", alpha=.3); plt.legend()
    plt.tight_layout(); plt.savefig(f"{outdir}/fig1_Pcr_vs_L_{tag}.png"); plt.close()
    try:
        _plot_section_modes(cfg, res, outdir, tag)
    except Exception:
        pass


# ==============================================================================
# 11.  MAIN
# ==============================================================================
_INPUT_WAIT = [0.0]                                       # time spent typing answers


def _timed_input_install():
    """Count the time spent waiting for the user's answers, so that the runtime
    report shows the computing time only (not the time spent typing)."""
    import builtins
    if getattr(builtins.input, "_gbt_timed", False):
        return
    _orig = builtins.input
    def _inp(prompt=""):
        t = time.perf_counter()
        try:
            return _orig(prompt)
        finally:
            _INPUT_WAIT[0] += time.perf_counter() - t
    _inp._gbt_timed = True
    builtins.input = _inp


def main():
    _timed_input_install()
    t_start = time.perf_counter()                        # start the wall clock
    cfg = get_config()                                   # ask everything
    if cfg.get("shape") == "omega":                      # omega panel: GBT / strip / FE
        run_omega_interactive(cfg)
        _system_report(time.perf_counter() - t_start)
        return
    res = analyse(cfg)                                   # run the analysis
    print(MODE_TEXT)                                     # explain the modes
    report(cfg, res)                                     # print the results table
    tag = cfg["shape"].upper()
    save_path = "out"                                    # portable output folder
    make_plots(cfg, res, outdir=save_path, tag=tag)          # write the figures
    import glob
    print(f"\nFigures written to {save_path} :")
    import os
    for f in sorted(glob.glob(os.path.join(save_path, f"*_{tag}.png"))):
        print("  ", f)
    print("\nDone.")
    _system_report(time.perf_counter() - t_start)        # runtime + hardware


# ==============================================================================
# 12.  CLASSICAL GBT CROSS-SECTION ENGINE  (Silvestre & Camotim; Kirchhoff walls)
#      Used only for the deformation-mode figures (fig6, fig_modes_interactive)
#      and the classical reference loads; the governing solvers are 23-25 / 20-22.
# ==============================================================================
import scipy.sparse as sp
import scipy.sparse.linalg as spla

GP4 = np.array([-0.8611363115940526, -0.3399810435848563,
                 0.3399810435848563,  0.8611363115940526])
GW4 = np.array([ 0.3478548451374538,  0.6521451548625461,
                 0.6521451548625461,  0.3478548451374538])


# ==============================================================================
# 1.  CLASSICAL LAMINATION THEORY  ->  wall constitutive matrix  (paper Eq. 4)
# ==============================================================================
def ply_Q(E1, E2, G12, v12):
    """Reduced stiffness of one orthotropic layer in its material axes."""
    v21 = v12 * E2 / E1
    den = 1.0 - v12 * v21
    return np.array([[E1 / den,       v12 * E2 / den, 0.0],
                     [v12 * E2 / den, E2 / den,       0.0],
                     [0.0,            0.0,            G12]])


def Qbar(Q, theta_deg):
    """Rotate the layer stiffness to the wall axes (x, s)  -- paper Eq. (4)."""
    c, s = np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg))
    Q11, Q12, Q22, Q66 = Q[0, 0], Q[0, 1], Q[1, 1], Q[2, 2]
    Qb = np.empty((3, 3))
    Qb[0, 0] = Q11*c**4 + 2*(Q12 + 2*Q66)*c*c*s*s + Q22*s**4
    Qb[1, 1] = Q11*s**4 + 2*(Q12 + 2*Q66)*c*c*s*s + Q22*c**4
    Qb[0, 1] = Qb[1, 0] = (Q11 + Q22 - 4*Q66)*c*c*s*s + Q12*(c**4 + s**4)
    Qb[2, 2] = (Q11 + Q22 - 2*Q12 - 2*Q66)*c*c*s*s + Q66*(c**4 + s**4)
    Qb[0, 2] = Qb[2, 0] = (Q11 - Q12 - 2*Q66)*c**3*s + (Q12 - Q22 + 2*Q66)*c*s**3
    Qb[1, 2] = Qb[2, 1] = (Q11 - Q12 - 2*Q66)*c*s**3 + (Q12 - Q22 + 2*Q66)*c**3*s
    return Qb


def laminate_C6(plies):
    """plies = [(Q(3x3) in material axes, angle_deg, thickness), ...]
    Returns the full 6x6 laminate matrix [[A,B],[B,D]] (strains ordered
    [e_xx, e_ss, g_xs, k_xx, k_ss, k_xs], k_xs = -2 w,xs engineering twist)."""
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
    t_tot = sum(p[2] for p in plies)
    z = -t_tot / 2.0
    for (Q, th, tp) in plies:
        Qb = Qbar(Q, th)
        z0, z1 = z, z + tp
        A += Qb * (z1 - z0)
        B += Qb * (z1**2 - z0**2) / 2.0
        D += Qb * (z1**3 - z0**3) / 3.0
        z = z1
    C6 = np.zeros((6, 6))
    C6[:3, :3] = A; C6[:3, 3:] = B; C6[3:, :3] = B; C6[3:, 3:] = D
    return C6


def condense_Nss(C6):
    """GBT constitutive treatment (Silvestre & Camotim 2002a): the kinematic
    constraint eps_ss^M = 0 is used together with N_ss = 0; statically
    condense the transverse membrane strain (index 1) out of the 6x6 matrix.
    Remaining strain order: [e_xx, g_xs, k_xx, k_ss, k_xs]  (5x5).
    For an isotropic wall: membrane -> E t (exact), bending -> plate D."""
    keep = [0, 2, 3, 4, 5]
    Ckk = C6[np.ix_(keep, keep)]
    c1 = C6[np.ix_(keep, [1])]
    return Ckk - (c1 @ c1.T) / C6[1, 1]


# ==============================================================================
# 2.  CROSS-SECTION GEOMETRY  (natural + intermediate nodes,  paper Fig. 1)
# ==============================================================================
class Wall:
    """One flat wall between two natural nodes (fold lines / free edges)."""
    def __init__(self, na, nb, C5, n_int):
        self.na, self.nb = na, nb          # natural node ids
        self.C5 = C5                        # condensed 5x5 constitutive matrix
        self.n_int = n_int                  # number of intermediate nodes
        self.chain = []                     # all node ids  na .. int .. nb
        self.t = None; self.n = None        # unit tangent / normal (2D)
        self.L = None                       # wall width  b_i
        # convenience material scalars
        self.A11c = C5[0, 0]                # condensed membrane axial stiffn.
        self.D22c = C5[3, 3]                # condensed transverse plate bend.


class Section:
    """Open thin-walled cross-section (branchless chain of walls: C, Z, ...).

    Node bookkeeping:
       nodes[i]      = (Y, Z) coordinates
       natural ids   = 0 .. n_nat-1  (in chain order)
       intermediate  = appended afterwards, wall by wall
    """
    def __init__(self, nat_coords, walls_spec):
        self.nodes = [np.array(p, float) for p in nat_coords]
        self.n_nat = len(nat_coords)
        self.walls = []
        for (na, nb, C5, n_int) in walls_spec:
            w = Wall(na, nb, C5, n_int)
            pa, pb = self.nodes[na], self.nodes[nb]
            d = pb - pa
            w.L = float(np.hypot(*d))
            w.t = d / w.L
            w.n = np.array([-w.t[1], w.t[0]])    # +90 deg rotation of t
            w.chain = [na]
            for q in range(n_int):               # equally spaced interm. nodes
                self.nodes.append(pa + d * (q + 1) / (n_int + 1))
                w.chain.append(len(self.nodes) - 1)
            w.chain.append(nb)
            self.walls.append(w)
        self.n_nodes = len(self.nodes)
        # node -> list of (wall id, 'end'/'interior') incidences
        self.inc = [[] for _ in range(self.n_nodes)]
        for iw, w in enumerate(self.walls):
            self.inc[w.na].append(iw)
            self.inc[w.nb].append(iw)
            for q in w.chain[1:-1]:
                self.inc[q].append(iw)
        # free-end natural nodes (belong to exactly one wall)
        self.tips = [i for i in range(self.n_nat) if len(self.inc[i]) == 1]
        # sub-segments (frame members & integration cells)
        self.segs = []                            # (wall id, node a, node b, l)
        for iw, w in enumerate(self.walls):
            ls = w.L / (w.n_int + 1)
            for a, b in zip(w.chain[:-1], w.chain[1:]):
                self.segs.append((iw, a, b, ls))


# ==============================================================================
# 3.  TRANSVERSE FRAME ANALYSIS  (in-plane w(s) of the conventional modes)
#     -- the "cross-section analysis" of the paper / Silvestre-Camotim 2002a
# ==============================================================================
def _beam_k(EI, l):
    """4x4 transverse bending stiffness of a frame member, DOFs (wA,thA,wB,thB),
    th = dw/ds (CCW rotation; consistent because n = +90deg rot of t)."""
    return EI / l**3 * np.array([[ 12,    6*l,  -12,    6*l ],
                                 [ 6*l, 4*l*l, -6*l,  2*l*l ],
                                 [-12,   -6*l,   12,   -6*l ],
                                 [ 6*l, 2*l*l, -6*l,  4*l*l ]])


def frame_solve(sec, node_kind):
    """Solve the statically indeterminate transverse frame of the section.

    node_kind[i] = ('vec', d)  : node displacement VECTOR d prescribed
                                 (interior fold nodes -- from Vlasov compat.)
                   ('val', w)  : transverse displacement prescribed (node on a
                                 single wall: intermediate / free end)
                   'free'      : transverse displacement unknown

    All nodal rotations are unknown (rigid fold lines -> one rotation/node).
    Returns dict  (seg index) -> (wA, thA, wB, thB)   Hermite end values."""
    nfree = {}
    ndof = sec.n_nodes                         # rotations first
    for i in range(sec.n_nodes):
        if node_kind[i] == 'free':
            nfree[i] = ndof; ndof += 1
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    def w_of(node, iw):
        """(known?, value-or-dof) of the transverse displacement of `node`
        as seen by wall iw (projection on that wall's normal)."""
        kind = node_kind[node]
        if kind == 'free':
            return False, nfree[node]
        tag, val = kind
        if tag == 'vec':
            return True, float(np.dot(val, sec.walls[iw].n))
        return True, float(val)                # 'val'

    seg_dofs = []
    for (iw, a, b, ls) in sec.segs:
        ke = _beam_k(sec.walls[iw].D22c, ls)
        ent = []
        for (node, ipos) in ((a, 0), (b, 2)):
            known, v = w_of(node, iw)
            ent.append(('k', v) if known else ('u', v))   # w   DOF
            ent.append(('u', node))                        # theta DOF
        seg_dofs.append(ent)
        # assemble (unknowns into K, knowns into F)
        for r in range(4):
            tr, vr = ent[r]
            if tr == 'k':
                continue
            for c in range(4):
                tc, vc = ent[c]
                if tc == 'u':
                    K[vr, vc] += ke[r, c]
                else:
                    F[vr] -= ke[r, c] * vc
    # rotations of nodes not touching any member never occur; K is SPD here
    x = np.linalg.lstsq(K, F, rcond=None)[0] if ndof else np.zeros(0)

    out = {}
    for isg, ent in enumerate(seg_dofs):
        vals = []
        for (t, v) in ent:
            vals.append(v if t == 'k' else x[v])
        out[isg] = tuple(vals)
    return out


def joint_disp(sec, vwall):
    """In-plane displacement VECTOR of every interior natural node from the
    Vlasov compatibility  d . t_wall = v_wall  of all adjoining walls."""
    d = {}
    for i in range(sec.n_nat):
        wl = [iw for iw in set(sec.inc[i]) if i in (sec.walls[iw].na,
                                                    sec.walls[iw].nb)]
        if len(wl) < 2:
            continue                            # free end -> handled elsewhere
        Tm = np.array([sec.walls[iw].t for iw in wl])
        rhs = np.array([vwall[iw] for iw in wl])
        d[i], *_ = np.linalg.lstsq(Tm, rhs, rcond=None)
    return d


# ==============================================================================
# 4.  DEFORMATION MODES  (paper Figs. 2 & 3, in the equivalent nodal basis)
# ==============================================================================
def _hermite(xi, l):
    """Cubic Hermite value/1st/2nd derivative at xi in [0,1] (length l)."""
    H  = np.array([1 - 3*xi**2 + 2*xi**3, l*(xi - 2*xi**2 + xi**3),
                   3*xi**2 - 2*xi**3,     l*(-xi**2 + xi**3)])
    dH = np.array([(-6*xi + 6*xi**2)/l,   1 - 4*xi + 3*xi**2,
                   (6*xi - 6*xi**2)/l,    -2*xi + 3*xi**2])
    ddH = np.array([(-6 + 12*xi)/l**2,    (-4 + 6*xi)/l,
                    (6 - 12*xi)/l**2,     (-2 + 6*xi)/l])
    return H, dH, ddH


class ModeSet:
    """Builds the conventional (k) and shear (j) deformation modes and samples
    every field  u_k, v_k, w_k, w_k,s, w_k,ss, u_j, u_j,s  at the cross-section
    Gauss points (4 per sub-segment) for the energy integrations."""

    def __init__(self, sec):
        self.sec = sec
        # ---- Gauss points on the section -------------------------------
        pts = []                                 # (seg idx, xi, weight)
        for isg, (iw, a, b, ls) in enumerate(sec.segs):
            for gp, gw in zip((GP4 + 1) / 2, GW4 / 2):
                pts.append((isg, gp, gw * ls))
        self.pts = pts
        self.npt = len(pts)
        self.wall_of_pt = np.array([sec.segs[p[0]][0] for p in pts])

        # ---- build the conventional modes ------------------------------
        conv = []                                # list of dicts of pt-arrays
        self.conv_label = []

        # (a) axial extension  (paper mode 1):  u = 1, v = w = 0
        conv.append(self._conv_fields(u_nat=np.ones(sec.n_nat)))
        self.conv_label.append('axial')

        # (b) warping modes: unit warping at natural nodes 1..n_nat-1
        #     (node 0 dropped -- together with 'axial' they span all
        #      linear-per-wall warping fields, i.e. the paper's global +
        #      distortional conventional modes)
        for r in range(1, sec.n_nat):
            e = np.zeros(sec.n_nat); e[r] = 1.0
            conv.append(self._conv_fields(u_nat=e))
            self.conv_label.append(f'warp@{r}')

        # (c) local-plate modes: unit transverse displacement at every
        #     intermediate node and every free-end node  (u = v = 0)
        loc_nodes = [i for i in range(sec.n_nodes)
                     if i >= sec.n_nat or i in sec.tips]
        for q in loc_nodes:
            conv.append(self._conv_fields(local_node=q))
            self.conv_label.append(f'local@{q}')

        self.M = len(conv)
        self.axial_idx = 0
        # pack conventional fields:  (npt, M)
        self.U   = np.stack([c['u']   for c in conv], axis=1)
        self.V   = np.stack([c['v']   for c in conv], axis=1)
        self.W   = np.stack([c['w']   for c in conv], axis=1)
        self.Ws  = np.stack([c['ws']  for c in conv], axis=1)
        self.Wss = np.stack([c['wss'] for c in conv], axis=1)

        # ---- shear modes: hat warping at every node except node 0 ------
        #      (excluding one node removes the constant warping already
        #       represented by the axial conventional mode)
        sh_nodes = list(range(1, sec.n_nodes))
        self.S = len(sh_nodes)
        Uj  = np.zeros((self.npt, self.S))
        Ujs = np.zeros((self.npt, self.S))
        for jj, q in enumerate(sh_nodes):
            for ip, (isg, xi, wgt) in enumerate(self.pts):
                iw, a, b, ls = sec.segs[isg]
                if a == q:
                    Uj[ip, jj] = 1 - xi; Ujs[ip, jj] = -1.0 / ls
                elif b == q:
                    Uj[ip, jj] = xi;     Ujs[ip, jj] =  1.0 / ls
        self.Uj, self.Ujs = Uj, Ujs

    # ------------------------------------------------------------------
    def _conv_fields(self, u_nat=None, local_node=None):
        """Sample one conventional mode at all Gauss points."""
        sec = self.sec
        u  = np.zeros(self.npt); v  = np.zeros(self.npt)
        w  = np.zeros(self.npt); ws = np.zeros(self.npt)
        wss = np.zeros(self.npt)

        if u_nat is not None:
            # -- warping mode:  u linear per wall,  v = -u,s  (Vlasov) ----
            vwall = np.array([-(u_nat[wl.nb] - u_nat[wl.na]) / wl.L
                              for wl in sec.walls])
            djoint = joint_disp(sec, vwall)
            node_kind = []
            for i in range(sec.n_nodes):
                if i in djoint:
                    node_kind.append(('vec', djoint[i]))
                else:
                    node_kind.append('free')     # tips & intermediate nodes
            hv = frame_solve(sec, node_kind)
            vmode = vwall
        else:
            # -- local mode: unit w at `local_node`, zero at all others ---
            node_kind = []
            for i in range(sec.n_nodes):
                if i == local_node:
                    node_kind.append(('val', 1.0))
                elif len(set(sec.inc[i])) >= 2:  # interior fold node
                    node_kind.append(('vec', np.zeros(2)))
                else:
                    node_kind.append(('val', 0.0))
            hv = frame_solve(sec, node_kind)
            vmode = np.zeros(len(sec.walls))

        # sample
        for ip, (isg, xi, wgt) in enumerate(self.pts):
            iw, a, b, ls = sec.segs[isg]
            wl = sec.walls[iw]
            H, dH, ddH = _hermite(xi, ls)
            q = np.array(hv[isg])
            w[ip]  = H  @ q
            ws[ip] = dH @ q
            wss[ip] = ddH @ q
            v[ip] = vmode[iw]
            if u_nat is not None:
                # position of the point along the WHOLE wall
                idx = wl.chain.index(a)
                sfrac = (idx + xi) / (wl.n_int + 1)
                u[ip] = u_nat[wl.na] * (1 - sfrac) + u_nat[wl.nb] * sfrac
        return dict(u=u, v=v, w=w, ws=ws, wss=wss)


# ==============================================================================
# 7.  BUCKLING EIGENVALUE PROBLEM      (K - P G) d = 0
# ==============================================================================
def buckle(K, G, keep, nev=6):
    Kc = K[keep][:, keep].tocsc()
    Gc = G[keep][:, keep].tocsc()
    n = Kc.shape[0]
    # solve  G d = mu K d,  P_cr = 1 / mu_max   (K is SPD after the BCs)
    try:
        mu = spla.eigsh(Gc, k=min(nev, n - 2), M=Kc, which='LA',
                        maxiter=6000, tol=0)[0]
        mu = mu[np.isfinite(mu) & (mu > 1e-14)]
        if len(mu):
            return 1.0 / mu.max()
    except Exception:
        pass
    # dense fallback
    from scipy.linalg import eigh
    mu = eigh(Gc.toarray(), Kc.toarray(), eigvals_only=True)
    mu = mu[np.isfinite(mu) & (mu > 1e-14)]
    return 1.0 / mu.max()


# ==============================================================================
# 9.  WALL MATERIAL HELPERS OF THE CLASSICAL ENGINE
# ==============================================================================
def material_C5_iso(E, v, t):
    Q = np.array([[E/(1-v*v), v*E/(1-v*v), 0],
                  [v*E/(1-v*v), E/(1-v*v), 0],
                  [0, 0, E/(2*(1+v))]])
    return condense_Nss(laminate_C6([(Q, 0.0, t)]))


def material_C5_layup(E1, E2, G12, v12, angles, t_total, rbs=False, keep_d16=True):
    Q = ply_Q(E1, E2, G12, v12)
    tp = t_total / len(angles)
    C6 = laminate_C6([(Q, a, tp) for a in angles])
    if rbs:
        A = C6[:3, :3]; B = C6[:3, 3:]; D = C6[3:, 3:]
        if np.abs(B).max() > 1e-6:                       # unsymmetric -> reduce D
            Dstar = D - B @ np.linalg.inv(A) @ B
            
            # تطبيق تعديل الـ D16 للقطاعات المفتوحة فقط
            if keep_d16:
                Dstar[0, 2] = D[0, 2]; Dstar[2, 0] = D[2, 0]      # D16 unchanged
                Dstar[1, 2] = D[1, 2]; Dstar[2, 1] = D[2, 1]      # D26 unchanged
                
            C6 = np.block([[A, np.zeros((3, 3))],
                           [np.zeros((3, 3)), Dstar]])
    return condense_Nss(C6)


# ==============================================================================
# 14b. OMEGA-STIFFENED PANEL -- GEOMETRY, LAMINATES, CONSTANTS (the Abaqus model)
#      skin width 2 bs; omega stringers at x = -bs/2, -3bs/2; b3 = 27.75, h = 30,
#      web angle 25 deg; foot = skin + flange (one stacked laminate, reference =
#      skin mid-plane); stringer mid-surface offset c = (ts + tf)/2 (Abaqus *TIE);
#      1 N/mm on every end edge; 0-deg plies run along the stringers.
# ==============================================================================
OMEGA_B3, OMEGA_H, OMEGA_THETA = 27.75, 30.0, 25.0
OMEGA_BS = [150.0, 285.5]
# End condition "CC" of the omega panel (both solvers):
#   "edge-load"   : clamped end edges -- every end node U1 = U2 = UR1 = UR2 =
#                   UR3 = 0, the z = 0 edge also U3 = 0, and the uniform
#                   1 N/mm line load on every end edge at z = L (each node free
#                   in U3).  This is the counterpart of the SS model (same edge
#                   load, ends clamped instead of U2 = 0) and it is what your
#                   Abaqus CC results correspond to: the 30 CC cases within
#                   2.8 % (mean 0.65 %); with rigid plates up to 15 % apart.
#   "rigid-plate" : rigid end plates (kinematic coupling to a reference point,
#                   U3 common to the whole end section) -- the load then goes
#                   to the walls in proportion to their axial stiffness.
#                   Run with --omega-cc-plate.
OMEGA_CC_END = "edge-load"
# Longitudinal (side) edges of the omega panel, x = 0 and x = -W (both at the
# middle of a skin bay), in both solvers:
#   "symmetry" : U1 = UR2 = UR3 = 0 on both side edges (XSYMM) -- the Abaqus
#                validation model (BC_Symmetry); default.
#   "periodic" : the two side edges are tied, u(0, z) = u(-W, z) for all six
#                DOFs at every z (the panel is one repeating cell of an
#                infinitely wide panel); one U1 anchor removes the remaining
#                rigid translation.  Run with --omega-periodic.
OMEGA_SIDE_BC = "symmetry"
# periodic only: also let the cell take a free uniform transverse strain and
# in-plane shear (constant jumps of U1 and U3 between the two side edges),
# i.e. an infinitely wide panel that is free to expand laterally.
OMEGA_PERIODIC_MACRO = False
OMEGA_BETA = [0.75, 1.5, 3.0]
OMEGA_ISO = dict(type="isotropic", E=200000.0, nu=0.3)
OMEGA_CFRP = dict(type="orthotropic", E1=150000.0, E2=9000.0, G12=5300.0,
                  nu12=0.33, G13=5300.0, G23=5300.0)
OMEGA_LAYUPS = [
    ("Isotropic", "iso",  [], []),
    ("Cross-ply", "cfrp", [0, 90, 0, 90, 90, 0, 90, 0],         [0, 90, 0, 90, 0, 90, 0]),
    ("Sym1",      "cfrp", [45, -45, -45, 45, 45, -45, -45, 45], [45, -45, 45, -45, 45, -45, 45]),
    ("Sym2",      "cfrp", [60, -30, -30, 60, 60, -30, -30, 60], [60, -30, 60, -30, 60, -30, 60]),
    ("Unsym",     "cfrp", [0, 30, 60, 90, 0, 30, 60, 90],       [0, 30, 60, 90, 0, 30, 60]),
]


def _clt_from(material, plies, z_bottom):
    """A, B, D, t of a ply stack whose first ply starts at z = z_bottom
    (plies listed bottom -> top, +z = the strip's local w direction)."""
    ttot = sum(t for _a, t in plies)
    A0, B0, D0, _ = clt_plies(material, plies)          # about own mid-plane
    e = z_bottom + ttot / 2.0                           # shift of the mid-plane
    return A0, B0 + e * A0, D0 + 2.0 * e * B0 + e * e * A0, ttot


def make_omega_panel(bs, mat_type, skin_seq, stiff_seq,
                     fiber_ref="longitudinal", strip=8.0, spec=None):
    """Cross-section of the Abaqus omega panel as walls with explicit ABD.

    Returns dict(nodes, walls, sides, info).  Each wall carries 'ABD'
    (reference = its own mid-line), 'edges' (loaded end edges it represents:
    2 for a bonded foot) and 'yload' (height of the resultant of its edge
    loads) so that the prestress can reproduce the Abaqus line load."""
    g = spec or {}              # spec: interactive input (material, plies, b3, h, theta)
    b3 = g.get("b3", OMEGA_B3); h = g.get("h", OMEGA_H)
    th = np.radians(g.get("theta", OMEGA_THETA))
    ht = b3 - h * np.tan(th)
    b1 = 0.1 * (bs - 2.0 * b3)
    W = 2.0 * bs
    if "pskin" in g:                                     # explicit plies (0 deg = along member)
        mat = g["material"]; pskin = list(g["pskin"]); pstiff = list(g["pstiff"])
    elif mat_type == "iso":
        mat = OMEGA_ISO
        pskin, pstiff = [(0.0, 1.0)], [(0.0, 0.875)]
    else:
        mat = OMEGA_CFRP
        sh = 90.0 if fiber_ref == "abaqus_global" else 0.0
        ang = lambda a: ((a + sh + 90.0) % 180.0) - 90.0
        pskin = [(ang(a), 0.125) for a in skin_seq]
        pstiff = [(ang(a), 0.125) for a in stiff_seq]
    ts = sum(t for _a, t in pskin); tf = sum(t for _a, t in pstiff)
    clear = 0.5 * (ts + tf)                              # stringer offset (tie)
    ABD_skin = _clt_from(mat, pskin, -ts / 2.0)
    ABD_stf = _clt_from(mat, pstiff, -tf / 2.0)
    ABD_foot = _clt_from(mat, pskin + pstiff, -ts / 2.0)  # stacked, ref = skin

    # ---- nodes: skin line (left -> right, so local w = +y = Abaqus normal) ---
    centers = [-bs / 2.0, -(bs / 2.0 + bs)]
    xs = [-W, 0.0]
    for c in centers:
        xs += [c - b3 - b1, c - b3, c + b3, c + b3 + b1]
    xs = sorted(set(round(x, 9) for x in xs))
    nodes = [(x, 0.0) for x in xs]
    idx = {round(x, 9): k for k, x in enumerate(xs)}
    walls = []
    def add(i, j, ABD, part, edges, yload):
        walls.append(dict(i=i, j=j, ABD=ABD, part=part, edges=edges, yload=yload))
    feet = set()
    for c in centers:
        feet.add((round(c - b3 - b1, 9), round(c - b3, 9)))
        feet.add((round(c + b3, 9), round(c + b3 + b1, 9)))
    for k in range(len(xs) - 1):
        seg = (round(xs[k], 9), round(xs[k + 1], 9))
        if seg in feet:
            add(k, k + 1, ABD_foot, "foot", 2, 0.5 * clear)
        else:
            add(k, k + 1, ABD_skin, "skin", 1, 0.0)
    # ---- stringers: webs + crown, oriented so local w points OUT of the hat -
    cells = []
    for c in centers:
        rl, rr = idx[round(c - b3, 9)], idx[round(c + b3, 9)]
        nodes.append((c - ht, h + clear)); tl = len(nodes) - 1
        nodes.append((c + ht, h + clear)); tr = len(nodes) - 1
        add(rl, tl, ABD_stf, "web", 1, 0.5 * (clear + h + clear))
        add(tl, tr, ABD_stf, "crown", 1, h + clear)
        add(tr, rr, ABD_stf, "web", 1, 0.5 * (clear + h + clear))
        cells.append([rl, tl, tr, rr])
    sides = [idx[round(-W, 9)], idx[0.0]]
    info = dict(bs=bs, W=W, b1=b1, ht=ht, ts=ts, tf=tf, clear=clear,
                material=mat, fiber_ref=fiber_ref, cells=cells, strip=strip,
                pskin=pskin, pstiff=pstiff)
    return dict(nodes=np.asarray(nodes, float), walls=walls, sides=sides, info=info)


# ==============================================================================
# 14e. MINDLIN STRIP KERNEL (_mfs_*) and the omega strip geometry, shared by the
#      strip solvers of sections 19, 20 and 22.
# ==============================================================================
_MFS_G2 = np.array([-1.0, 1.0]) / np.sqrt(3.0)


def _mfs_long_functions(bc, L, mlist, nx=None):
    """Gauss points/weights along the member and the four function families
    S0 = Y, S1 = Y', S2 = c Y', S3 = c Y''  (M x nx), c = L/(m pi)."""
    if nx is None:
        nx = max(96, 8 * (max(mlist) + 2))
    xg, wg = np.polynomial.legendre.leggauss(nx)
    x = 0.5 * L * (xg + 1.0); w = 0.5 * L * wg
    k = np.asarray(mlist, float)[:, None] * np.pi / L; p = np.pi / L; c = 1.0 / k
    if bc == "SS":
        Y = np.sin(k * x); Y1 = k * np.cos(k * x); Y2 = -k * k * Y
    elif bc == "CC":
        a, b = np.sin(k * x), np.sin(p * x); da, db = k * np.cos(k * x), p * np.cos(p * x)
        Y = a * b; Y1 = da * b + a * db; Y2 = -k * k * a * b + 2 * da * db - p * p * a * b
    else:
        raise ValueError("omega GBT: bc must be 'SS' or 'CC'")
    return x, w, [Y, Y1, c * Y1, c * Y2]


def _mfs_gram(S, w):
    return [[(S[a] * w) @ S[b].T for b in range(4)] for a in range(4)]


# strain terms: (component, field, s-operator 0=N 1=dN/ds, x-family, coefficient)
# fields 0 u 1 v 2 w 3 wx 4 ws 5 wn ; components 0 eps_x 1 eps_s 2 gam_xs 3 kap_x
# 4 kap_s 5 kap_xs 6 gam_xn 7 gam_sn 8 drill
_MFS_TERMS = [
    (0, 0, 0, 3, 1.0), (1, 1, 1, 0, 1.0), (2, 0, 1, 2, 1.0), (2, 1, 0, 1, 1.0),
    (3, 4, 0, 3, 1.0), (4, 3, 1, 0, -1.0), (5, 4, 1, 2, 1.0), (5, 3, 0, 1, -1.0),
    (6, 2, 0, 1, 1.0), (6, 4, 0, 2, 1.0), (7, 2, 1, 0, 1.0), (7, 3, 0, 0, -1.0),
    (8, 5, 0, 2, 1.0), (8, 1, 0, 1, -0.5), (8, 0, 1, 2, 0.5),
]
_MFS_GX = [(0, 0, 3), (1, 0, 1), (2, 0, 1)]          # u,x  v,x  w,x
_MFS_GS = [(0, 1, 2), (1, 1, 0), (2, 1, 0)]          # u,s  v,s  w,s


def _mfs_svec(field, sop, r, b):
    v = np.zeros(12)
    if sop == 0:
        v[field] = 1.0 - r; v[6 + field] = r
    else:
        v[field] = -1.0 / b; v[6 + field] = 1.0 / b
    return v


def _mfs_C9(ABD, As, drill_ratio=1e-3):
    A, B, D = ABD
    C = np.zeros((9, 9))
    C[:3, :3] = A; C[:3, 3:6] = B; C[3:6, :3] = B; C[3:6, 3:6] = D
    C[6:8, 6:8] = As; C[8, 8] = drill_ratio * A[2, 2]
    return C


def _mfs_strip_blocks(b, ABD, As, P):
    """Local 12x12 blocks per pair of x-families:  K = sum kron(Kab, Gram_ab),
    Kg = sum kron(Gab, Gram_ab).  P = (Px, Ps, Pxs), compression positive."""
    C = _mfs_C9(ABD, As)
    Kab = [[np.zeros((12, 12)) for _ in range(4)] for _ in range(4)]
    Gab = [[np.zeros((12, 12)) for _ in range(4)] for _ in range(4)]
    Px, Ps, Pxs = P
    for g in _MFS_G2:
        r = 0.5 * (g + 1.0); wt = 0.5 * b
        E = [np.zeros((9, 12)) for _ in range(4)]
        for (comp, f, sop, fam, cf) in _MFS_TERMS:
            E[fam][comp] += cf * _mfs_svec(f, sop, 0.5 if comp == 7 else r, b)
        for a in range(4):
            for bb in range(4):
                Kab[a][bb] += wt * (E[a].T @ C @ E[bb])
        for (f, sop, fam) in _MFS_GX:
            vx = _mfs_svec(f, sop, r, b)
            Gab[fam][fam] += wt * Px * np.outer(vx, vx)
        for (f, sop, fam) in _MFS_GS:
            vs = _mfs_svec(f, sop, r, b)
            Gab[fam][fam] += wt * Ps * np.outer(vs, vs)
        if Pxs != 0.0:
            for (f, sop, fam), (f2, sop2, fam2) in zip(_MFS_GX, _MFS_GS):
                vx = _mfs_svec(f, sop, r, b); vs = _mfs_svec(f2, sop2, r, b)
                Gab[fam][fam2] += wt * Pxs * np.outer(vx, vs)
                Gab[fam2][fam] += wt * Pxs * np.outer(vs, vx)
    return Kab, Gab


def _mfs_T6(alpha):
    """Global [u, UX, UY, Rx, RX, RY] -> local [u, v, w, wx, ws, wn]."""
    c, s = np.cos(alpha), np.sin(alpha)
    T = np.zeros((6, 6))
    T[0, 0] = 1.0; T[1, 1], T[1, 2] = c, s; T[2, 1], T[2, 2] = -s, c
    T[3, 3] = 1.0; T[4, 4], T[4, 5] = c, s; T[5, 4], T[5, 5] = -s, c
    return T


def _mfs_link_matrix(N, nodes, links, fixed, extra=0):
    """Reduction full = C q for rigid links (slave, master) and fixed (node, dof);
    'extra' global unknowns are appended after the 6N nodal dofs."""
    slave = {s: m for s, m in links}; fixset = set(fixed)
    qidx = {}; nq = 0
    for n in range(N):
        if n in slave:
            continue
        for d in range(6):
            if (n, d) not in fixset:
                qidx[(n, d)] = nq; nq += 1
    for j in range(extra):
        qidx[("g", j)] = nq; nq += 1
    rr, cc, vv = [], [], []
    def add(row, key, v):
        if key in qidx:
            rr.append(row); cc.append(qidx[key]); vv.append(v)
    for n in range(N):
        m_ = slave.get(n, n)
        for d in range(6):
            add(6 * n + d, (m_, d), 1.0)
        if n in slave:                                   # R x r, R = (RX, RY, Rx)
            rX, rY = nodes[n] - nodes[m_]
            add(6 * n + 1, (m_, 3), -rY); add(6 * n + 2, (m_, 3), rX)
            add(6 * n + 0, (m_, 4), rY); add(6 * n + 0, (m_, 5), -rX)
    for j in range(extra):
        add(6 * N + j, ("g", j), 1.0)
    return sp.csr_matrix((vv, (rr, cc)), shape=(6 * N + extra, nq))


def _mfs_assemble(nodes, strips, L, mlist, bc, links=(), fixed=()):
    """K and Kg of the Mindlin strip model for the harmonic set mlist (all
    couplings between them included).  strips = [(iA, iB, ABD, As, P)]."""
    nodes = np.asarray(nodes, float); N = len(nodes); M = len(mlist); nd = 6 * N * M
    x, w, S = _mfs_long_functions(bc, L, mlist); G = _mfs_gram(S, w)
    rows, cols, kv, gv = [], [], [], []; cache = {}
    for (iA, iB, ABD, As, P) in strips:
        d = nodes[iB] - nodes[iA]; b = float(np.hypot(*d)); al = np.arctan2(d[1], d[0])
        key = (round(b, 10), id(ABD[0]), id(As), tuple(np.round(P, 13)))
        if key not in cache:
            cache[key] = _mfs_strip_blocks(b, ABD, As, P)
        Kab, Gab = cache[key]
        T = np.zeros((12, 12)); T[:6, :6] = _mfs_T6(al); T[6:, 6:] = _mfs_T6(al)
        Ke = np.zeros((12 * M, 12 * M)); Ge = np.zeros_like(Ke)
        for a in range(4):
            for bb in range(4):
                if np.any(Kab[a][bb]):
                    Ke += np.kron(T.T @ Kab[a][bb] @ T, G[a][bb])
                if np.any(Gab[a][bb]):
                    Ge += np.kron(T.T @ Gab[a][bb] @ T, G[a][bb])
        loc = np.concatenate([(6 * n + np.arange(6))[:, None] * M + np.arange(M)[None, :]
                              for n in (iA, iB)]).ravel()
        rows.append(np.repeat(loc, 12 * M)); cols.append(np.tile(loc, 12 * M))
        kv.append(Ke.ravel()); gv.append(Ge.ravel())
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    K = sp.coo_matrix((np.concatenate(kv), (rows, cols)), shape=(nd, nd)).tocsr()
    Kg = sp.coo_matrix((np.concatenate(gv), (rows, cols)), shape=(nd, nd)).tocsr()
    C = sp.kron(_mfs_link_matrix(N, nodes, links, fixed), sp.identity(M), format="csr")
    Kr = (C.T @ K @ C).tocsc(); Gr = (C.T @ Kg @ C).tocsc()
    return 0.5 * (Kr + Kr.T), 0.5 * (Gr + Gr.T), C


def _mfs_buckle(Kr, Gr, nev):
    """Lowest positive eigenvalues of K phi = lambda Kg phi (+ vectors)."""
    lu = spla.splu(Kr.tocsc(), permc_spec="COLAMD")
    Kinv = spla.LinearOperator(Kr.shape, matvec=lu.solve, dtype=float)
    k = max(1, min(nev, Kr.shape[0] - 2))
    mu, V = spla.eigsh(Gr, k=k, M=Kr, Minv=Kinv, which="LA",
                       ncv=min(Kr.shape[0] - 1, max(3 * k, 40)), tol=1e-10, maxiter=50000)
    ok = mu > 0; lam = 1.0 / mu[ok]; V = V[:, ok]; o = np.argsort(lam)
    return lam[o], V[:, o]


def _omega_mfs_geometry(panel, ds=2.0):
    """Strip model of the omega panel on the Abaqus geometry.  Skin and the
    stacked feet on y = 0; each web root at y = clear is a slave node rigidly
    linked to the skin node below it; webs and crown up to h + clear."""
    info = panel["info"]; nd0 = panel["nodes"]; clear = info["clear"]; mat = info["material"]
    As_s = _osfe_laminate(mat, info["pskin"])[3]
    As_t = _osfe_laminate(mat, info["pstiff"])[3]
    As_f = _osfe_laminate(mat, list(info["pskin"]) + list(info["pstiff"]))[3]
    pts = [tuple(p) for p in nd0]; links = []; strips = []; root = {}
    def root_of(i):
        if i not in root:
            pts.append((nd0[i][0], clear)); root[i] = len(pts) - 1
            links.append((root[i], i))
        return root[i]
    for w in panel["walls"]:
        i0, i1 = w["i"], w["j"]
        As = {"skin": As_s, "foot": As_f}.get(w["part"], As_t)
        if w["part"] == "web":
            if abs(nd0[i0][1]) < 1e-12: i0 = root_of(i0)
            if abs(nd0[i1][1]) < 1e-12: i1 = root_of(i1)
        P0, P1 = np.array(pts[i0]), np.array(pts[i1])
        n = max(2, int(np.ceil(np.hypot(*(P1 - P0)) / ds)))
        prev = i0
        for k in range(1, n + 1):
            if k < n:
                pts.append(tuple(P0 + (P1 - P0) * k / n)); nxt = len(pts) - 1
            else:
                nxt = i1
            strips.append((prev, nxt, w, As)); prev = nxt
    return np.asarray(pts, float), strips, links


# ==============================================================================
# 14c. LAMINATE STIFFNESS (_osfe_laminate: A, B, D, As with k = 5/6) AND THE
#      SPARSE SYMMETRIC SOLVER (_osfe_factor: PARDISO, else SuperLU).
# ==============================================================================
import scipy.sparse as _osp
import scipy.sparse.linalg as _ospla

# ---- Abaqus reference table (for PRINTING the comparison only; never used in
#      the solution).  From the omega_validation_eigs results you supplied.

# ------------------------------------------------------------------ laminate --
def _osfe_laminate(mat, plies, kshear=5.0 / 6.0):
    """A, B, D (3x3, engineering gamma_xy) and As (2x2 [xz, yz]) of a ply stack
    [(theta_deg, t), ...] listed bottom -> top along the shell normal; theta
    measured from local x (material-1 reference) towards local y."""
    if mat["type"] in ("iso", "isotropic"):
        E, nu = mat["E"], mat["nu"]; G = E / (2 * (1 + nu))
        E1 = E2 = E; nu12 = nu; G12 = G13 = G23 = G
    else:
        E1, E2, nu12, G12 = mat["E1"], mat["E2"], mat["nu12"], mat["G12"]
        G13 = mat.get("G13", G12); G23 = mat.get("G23", G12)
    nu21 = nu12 * E2 / E1; d = 1 - nu12 * nu21
    Q11, Q22, Q12, Q66 = E1 / d, E2 / d, nu12 * E2 / d, G12
    ttot = sum(t for _, t in plies); z = -ttot / 2.0
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3)); As = np.zeros((2, 2))
    for th, t in plies:
        c, s = np.cos(np.radians(th)), np.sin(np.radians(th))
        Qb = np.array([
            [Q11*c**4 + 2*(Q12+2*Q66)*s*s*c*c + Q22*s**4,
             (Q11+Q22-4*Q66)*s*s*c*c + Q12*(c**4+s**4),
             (Q11-Q12-2*Q66)*s*c**3 + (Q12-Q22+2*Q66)*s**3*c],
            [(Q11+Q22-4*Q66)*s*s*c*c + Q12*(c**4+s**4),
             Q11*s**4 + 2*(Q12+2*Q66)*s*s*c*c + Q22*c**4,
             (Q11-Q12-2*Q66)*s**3*c + (Q12-Q22+2*Q66)*s*c**3],
            [(Q11-Q12-2*Q66)*s*c**3 + (Q12-Q22+2*Q66)*s**3*c,
             (Q11-Q12-2*Q66)*s**3*c + (Q12-Q22+2*Q66)*s*c**3,
             (Q11+Q22-2*Q12)*s*s*c*c + Q66*(c*c-s*s)**2]])
        Qs = np.array([[G13*c*c + G23*s*s, (G13-G23)*c*s],
                       [(G13-G23)*c*s, G23*c*c + G13*s*s]])
        z1 = z + t
        A += Qb*(z1-z); B += Qb*(z1**2-z**2)/2; D += Qb*(z1**3-z**3)/3
        As += kshear * Qs * (z1-z); z = z1
    return A, B, D, As

# ------------------------------------------------------------------- element --


# ------------------------------------------------------------ linear solver --
def _osfe_factor(Kf):
    """Sparse SPD factorisation.  Intel PARDISO (pip install pypardiso) when
    available -- ~10x faster -- otherwise SciPy SuperLU with COLAMD ordering.
    Both are exact direct solvers; the choice changes speed only."""
    try:
        import pypardiso
        sol = pypardiso.PyPardisoSolver(); sol.set_matrix_type(2)
        Ku = _osp.triu(Kf, format="csr"); sol.factorize(Ku)
        try:                                   # direct back-substitution only:
            sol.set_iparm(8, 0); sol.set_phase(33)   # no iterative refinement,
            def _solve(b):                     # no re-check of the matrix per call
                b = np.asarray(b, float)
                x = sol._call_pardiso(Ku, np.asfortranarray(b.reshape(len(b), -1)))
                return x.reshape(b.shape)
            _solve(np.ones(Ku.shape[0]))
        except Exception:
            def _solve(b):
                return sol.solve(Ku, np.asarray(b, float))
        return _solve, (lambda: sol.free_memory(everything=True))
    except Exception:
        Kc = Kf.tocsc()
        try:                                   # SPD: symmetric ordering, no pivoting
            lu = _ospla.splu(Kc, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                             options=dict(SymmetricMode=True))   # (~20x faster than
            r = np.random.default_rng(0).standard_normal(Kc.shape[0])  # COLAMD here)
            if np.linalg.norm(Kc @ lu.solve(r) - r) <= 1e-6 * np.linalg.norm(r):
                return lu.solve, (lambda: None)
        except Exception:
            pass
        lu = _ospla.splu(Kc, permc_spec="COLAMD")      # general fallback (pivoting)
        return lu.solve, (lambda: None)

# ------------------------------------------------------------------- model ----


OMEGA_BCS = ["SS", "CC"]


def _omega_study(solver, csv_out, mesh=3.7, nmodes=6):
    """Common driver of the omega parametric studies (same layout as the other
    parametric studies: material, BC, geometry | N_cr | governing mode).
    solver = "gbt" (shear-deformable GBT, section 25), "strip" (shell-type
    Mindlin spline finite strip, section 22)."""
    title = {"gbt": f"shear-deformable GBT, {OMEGA_GBT_FAMILIES} modes",
             "strip": "shell-type Mindlin spline finite strip"}[solver]
    print("=" * 92)
    print(f" PARAMETRIC STUDY  -  OMEGA-STIFFENED PANEL  ({title})")
    print(" N_cr = critical intensity of the 1 N/mm end line load [N/mm]")
    print("   SS = ends U2 = 0 (Abaqus script)   CC = " +
          ("rigid non-rotating end plates" if OMEGA_CC_END == "rigid-plate"
           else "clamped end edges (U1=U2=UR=0), 1 N/mm edge load"))
    print("   side edges = " + ("symmetry U1 = UR2 = UR3 = 0 (Abaqus script)"
                               if OMEGA_SIDE_BC == "symmetry" else
                               "PERIODIC u(x=0) = u(x=-2bs), all 6 DOFs"
                               + (" + free lateral expansion" if OMEGA_PERIODIC_MACRO else "")))
    print("=" * 92)
    print(f"{'material':>10}{'BC':>4}{'bs':>7}{'beta':>6}{'L[mm]':>8} | "
          f"{'N_cr [N/mm]':>12}{'t [s]':>7}   governing")
    print("-" * 92)
    rows = []; t_all = time.perf_counter()
    for name, mtype, sk, st in OMEGA_LAYUPS:
        _pan = {}                                        # one panel / GBT section per bs
        for bc in OMEGA_BCS:
            for bs in OMEGA_BS:
                for beta in OMEGA_BETA:
                    L = beta * bs
                    t0 = time.perf_counter()
                    try:
                        if solver == "gbt":
                            if bs not in _pan:
                                _pan.clear(); _pan[bs] = make_omega_panel(bs, mtype, sk, st)
                            modes, inf = omega_gbt_buckling(_pan[bs], L, bc, nmodes=nmodes)
                        else:
                            modes, inf = omega_spline_eigs(make_omega_panel(bs, mtype, sk, st),
                                                           L, bc, nmodes=nmodes)
                        dt = time.perf_counter() - t0
                        N, m = modes[0]
                        gov = f"m={m}  {inf['mode_type']}"
                        rows.append((name, bc, bs, beta, L, N, m, inf["mode_type"], dt))
                        print(f"{name:>10}{bc:>4}{bs:>7g}{beta:>6g}{L:>8.1f} | "
                              f"{N:>12.4f}{dt:>7.1f}   {gov[:45]}")
                    except Exception as e:
                        rows.append((name, bc, bs, beta, L, None, None, "", time.perf_counter() - t0))
                        print(f"{name:>10}{bc:>4}{bs:>7g}{beta:>6g}{L:>8.1f} | "
                              f"{'FAILED':>12}   ({e})")
                    sys.stdout.flush()
        print("-" * 92)
    print(f" {len(rows)} cases in {time.perf_counter() - t_all:.1f} s.")
    try:
        with open(csv_out, "w") as f:
            f.write("material,BC,bs_mm,beta,L_mm,Ncr_Nmm,m,governing,time_s\n")
            for (nm, bc, bs, beta, L, N, m, kind, dt) in rows:
                f.write(f"{nm},{bc},{bs:g},{beta:g},{L:g},"
                        f"{'' if N is None else round(N, 4)},{'' if m is None else m},{kind},"
                        f"{dt:.1f}\n")
        print(f" results written to {csv_out}")
    except Exception:
        pass
    return rows


def run_parametric_omega(csv_out=None, nmodes=2):
    """Omega panel parametric study with the governing omega solver
    (OMEGA_SOLVER: "gbt" = section 25, default; "strip" = section 22):
    5 layups x SS/CC x bs {150, 285.5} x beta {0.75, 1.5, 3}."""
    solver = "strip" if OMEGA_SOLVER == "strip" else "gbt"
    if csv_out is None:
        csv_out = f"omega_panel_{'GBT' if solver == 'gbt' else 'strip'}_results" + (
            ("_periodic_free" if OMEGA_PERIODIC_MACRO else "_periodic")
            if OMEGA_SIDE_BC == "periodic" else "") + ".csv"
    return _omega_study(solver, csv_out, nmodes=nmodes)


# ------------------------------------------------------------------------------
# 14d.  OMEGA PANEL  -  INTERACTIVE ENTRY  (shape "omega" in the main menu)
#       Solver choice: GBT   = shear-deformable GBT (section 25)   [default]
#                      strip = shell-type Mindlin spline finite strip (22)
#                      FE    = shell finite element (14c); or combinations.
#       End conditions: SS (Abaqus script) or CC (clamped ends, see OMEGA_CC_END)
# ------------------------------------------------------------------------------
def _get_omega_config(material, iso, ask_layup):
    """Questions for the omega-stiffened panel (called from get_config).
    Every default is the Abaqus validation model (Sym1 lay-up, bs = 150,
    beta = 1.5), so <Enter> everywhere reproduces that parametric case."""
    print("\nOmega-stiffened panel (the Abaqus validation model): skin width 2*bs,")
    print("two omega stringers at x = -bs/2 and -3bs/2 tied to the skin, sides =")
    print("symmetry (U1 = UR2 = UR3 = 0), 1 N/mm on every end edge.")
    print("Ply angles: 0 deg runs ALONG the stringers (Abaqus AXIS_2 orientation).")
    print("Parametric cases: bs = 150 (narrow) / 285.5 (wide), L = beta*bs with")
    print("beta = 0.75 (short) / 1.5 (intermediate) / 3 (long).")
    bs = ask("  bs (stringer spacing) [mm]", 150.0, float)
    L = ask("  L  (panel length) [mm]", 1.5 * bs, float)
    b3 = ask("  b3 (half base width of the omega) [mm]", OMEGA_B3, float)
    h = ask("  h  (omega height) [mm]", OMEGA_H, float)
    th = ask("  web angle from the vertical [deg]", OMEGA_THETA, float)
    print("\nWall thickness:" if iso else "\nLaminate lay-up (default = Sym1 of the Abaqus model):")
    ang_s, tp_s = ask_layup("SKIN", 0.125, OMEGA_LAYUPS[2][2], 1.0)
    ang_t, tp_t = ask_layup("STRINGER", 0.125, OMEGA_LAYUPS[2][3], 0.875)
    pskin = [(float(a), tp_s) for a in ang_s]
    pstiff = [(float(a), tp_t) for a in ang_t]
    bc = ask("\nEnd conditions on the loaded (transverse) edges "
             "(SS = ends U2 = 0 / CC = clamped ends)", "SS", str, ["ss", "cc"]).upper()
    print("Longitudinal (side) edges x = 0 and x = -2bs (middle of a skin bay):")
    print("  1 = symmetry   U1 = UR2 = UR3 = 0             (the Abaqus validation model)")
    print("  2 = periodic   u(x=0) = u(x=-2bs), all 6 DOFs (repeating cell of an")
    print("                 infinitely wide panel)")
    sch = ask("  choose 1 or 2", 1, int, ["1", "2"])
    globals()["OMEGA_SIDE_BC"] = "symmetry" if sch == 1 else "periodic"
    globals()["OMEGA_PERIODIC_MACRO"] = False
    print("\nSolver:")
    print("  1 = GBT    shear-deformable Generalised Beam Theory          (default)")
    print("  2 = strip  shell-type Mindlin spline finite strip (fast reference)")
    print("  3 = GBT + strip (side by side)")
    snum = ask("  choose 1, 2 or 3", 2 if OMEGA_SOLVER == "strip" else 1, int, ["1", "2", "3"])
    solver = {1: "gbt", 2: "strip", 3: "gbt+strip"}[snum]
    use = solver.split("+")
    mesh = None
    pan_ = None
    if "gbt" in use:
        pan_, _gm = _omega_gbt_questions(bs, material, pskin, pstiff, b3, h, th)
    if "strip" in use:
        print("\nMesh of the spline strip:")
        print(f"  1 = standard  strips bs/50 = {bs/50:.1f} mm, spans bs/20 = {bs/20:.1f} mm"
              "   (converged: within 0.55 % of fine)")
        print(f"  2 = fast      strips bs/40 = {bs/40:.1f} mm, spans bs/15 = {bs/15:.1f} mm"
              "   (within 1 % of fine, ~1.5x faster)")
        print("  3 = fine      strips 2 mm, spans min(5, L/40) mm"
              "            (reference, ~5x slower)")
        mch = ask("  choose 1, 2 or 3", 1, int, ["1", "2", "3"])
        globals()["OMEGA_MESH"] = {1: "standard", 2: "fast", 3: "fine"}[mch]
    nmodes = ask("Number of eigenvalues to report", 4, int)
    nmodes = max(2, nmodes)
    print("\nN_cr(L) curve (same solver and mesh at every length; every length is a"
          "\nfull solve, so a curve of n lengths takes about n times longer):")
    npt = ask("  number of lengths (0 = no curve)", 0, int)
    Lmin, Lmax = round(0.5 * bs, 1), round(3.0 * bs, 1)
    if npt > 0:
        Lmin = ask("  Lmin [mm]", Lmin, float)
        Lmax = ask("  Lmax [mm]", Lmax, float)
    print("\nMode shapes to plot:")
    m1 = ask("  FIRST mode number", 1, int)
    m2 = ask("  SECOND mode number", 2, int)
    return dict(shape="omega", material=material, bs=bs, L=L, bc=bc, solver=solver, panel=pan_,
                mesh=mesh, nmodes=nmodes, Lmin=Lmin, Lmax=Lmax, npt=npt,
                plot_modes=(max(1, m1), max(1, m2)),
                spec=dict(material=material, pskin=pskin, pstiff=pstiff,
                          b3=b3, h=h, theta=th))


# ---------------------------------------------------------------- results ----
def _omega_U(M, vec, xs):
    """(N, 6, nx) node displacements of a spline-strip vector (or of a GBT
    solution, section 25) at stations xs."""
    if hasattr(M, "U_at"):
        return M.U_at(vec, xs)
    from scipy.interpolate import BSpline
    t = np.r_[[0.0] * 3, M.breaks, [M.L] * 3]
    return M.node_part(vec) @ BSpline.design_matrix(np.clip(xs, 0.0, M.L), t, 3).toarray().T


def _omega_part_share(M, strips, vec):
    """Share of the buckling displacement (sum over strips of width x |u|^2 of the
    in-plane cross-section translations, all stations) carried by the skin
    (skin + bonded feet) and by the stringer (webs + crown)."""
    xs = np.linspace(0.0, M.L, 121)
    U = _omega_U(M, vec, xs)
    tr2 = (U[:, 1, :] ** 2 + U[:, 2, :] ** 2).sum(axis=1)       # per node
    sk = st = 0.0
    for k, (iA, iB, w, _As) in enumerate(strips):
        if w["part"] == "link":
            continue
        e = 0.5 * M.S[k]["b"] * (tr2[iA] + tr2[iB])
        if w["part"] in ("skin", "foot"):
            sk += e
        else:
            st += e
    tot = max(sk + st, 1e-300)
    return sk / tot, st / tot


def _omega_section_props(M, strips):
    """Loaded edge length, area, EA and the height of the axial-stiffness centroid."""
    Lload = A = EA = EAy = 0.0
    for k, (iA, iB, w, _As) in enumerate(strips):
        if w["part"] == "link":
            continue
        b = M.S[k]["b"]; y = 0.5 * (M.nodes[iA][1] + M.nodes[iB][1])
        A11 = _beam_a11(w["ABD"][0], w["ABD"][1], w["ABD"][2]); t = float(w["ABD"][3])   # beam E t
        Lload += w["edges"] * b; A += t * b; EA += A11 * b; EAy += A11 * b * y
    return dict(Lload=Lload, A=A, EA=EA, ys=EAy / EA)


def _omega_draw_section(ax, M, U, x_i, color, amp_frac=0.18, title=None):
    """Undeformed (grey) and deformed (colour) cross-section at station x_i."""
    XY = M.nodes
    d = U[:, 1:3, x_i]
    dm = max(np.abs(d).max(), 1e-30)
    amp = amp_frac * max(np.ptp(XY[:, 0]), np.ptp(XY[:, 1])) / dm
    for (i, j, _abd, _as) in M.strips:
        ax.plot(XY[[i, j], 0], XY[[i, j], 1], color="0.78", lw=1.4, zorder=1)
        ax.plot(XY[[i, j], 0] + amp * d[[i, j], 0], XY[[i, j], 1] + amp * d[[i, j], 1],
                color=color, lw=1.9, zorder=2)
    for s_, m_ in M.links:                                # tie links (stringer feet)
        ax.plot(XY[[s_, m_], 0] + amp * d[[s_, m_], 0], XY[[s_, m_], 1] + amp * d[[s_, m_], 1],
                color=color, lw=1.0, ls=":", zorder=2)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=10)


def _omega_long_line(M, U, xs):
    """Out-of-plane skin displacement along the line through its maximum (the
    same line the Abaqus script uses to count half-waves), normalised."""
    sk = np.where(np.abs(M.nodes[:, 1]) < 1e-9)[0]
    tr = np.hypot(U[:, 1, :], U[:, 2, :])
    n_all, _ = np.unravel_index(np.argmax(tr), tr.shape)
    if len(sk):
        n0, _ = np.unravel_index(np.argmax(np.abs(U[sk, 2, :])), (len(sk), len(xs)))
        line = U[sk[n0], 2, :]
    else:
        line = U[n_all, 2, :]
    return line / max(np.abs(line).max(), 1e-30)


def _omega_report(cfg, R):
    """Result window of the omega panel, laid out like the member sections."""
    M, strips, modes, kinds, shares = R["M"], R["strips"], R["modes"], R["kinds"], R["shares"]
    N0, m0 = modes[0]
    pr = R["props"]; info = R["panel"]["info"]
    line = "=" * 74
    print("\n" + line)
    print(f"  RESULT   OMEGA PANEL   bs = {cfg['bs']:g} mm,  L = {cfg['L']:g} mm,  BC = {cfg['bc']}")
    print(f"  ends: {cfg['bc']}"
          + (" (U2 = 0)" if cfg['bc'] == 'SS' else
             (" (clamped, uniform edge load)" if OMEGA_CC_END == "edge-load"
              else " (rigid end plates)"))
          + f";  sides: {OMEGA_SIDE_BC}"
          + (" (U1 = UR2 = UR3 = 0)" if OMEGA_SIDE_BC == "symmetry"
             else " (u(x=0) = u(x=-2bs))" + (" + free lateral expansion"
                                              if OMEGA_PERIODIC_MACRO else "")))
    print(line)
    print(f"  N_cr = {N0:.4f} N/mm   (critical average line load on the end edges)")
    print(f"         total end load P_cr = N_cr x {pr['Lload']:.1f} mm = {N0 * pr['Lload'] / 1e3:.3f} kN")
    print(f"  solver: {R['solver_text']}")
    print(f"  critical mode: {kinds[0]},  {m0} half-wave(s),  "
          f"displacement share skin {100 * shares[0][0]:.0f} % / stringer {100 * shares[0][1]:.0f} %")
    if len(modes) > 1:
        print(f"  next eigenvalue: {modes[1][0]:.4f} N/mm  (m = {modes[1][1]})")
    if R.get("part"):
        print("  GBT modal participation of the critical mode: " + "  ".join(
            f"{_GBT_GROUP_NAMES.get(g, g)} {100 * v:.0f} %"
            for g, v in sorted(R["part"].items(), key=lambda t: -t[1]) if v >= 0.005))
    if R.get("mode_part"):
        top = [mp for mp in R["mode_part"] if mp[1] >= 0.02][:8]
        kinds_ = {k: nm for (k, _l, _s, nm) in omega_gbt_mode_table(R["section"], max(k for k, _v in top))}
        print("  participation of the conventional GBT modes (numbers of the mode list): "
              + ",  ".join(f"mode {k} ({kinds_[k]}) {100 * v:.0f} %" for k, v in top))
    print(f"  >> REPORTED N_cr                =   {N0:.4f} N/mm <<")
    print("\n" + line); print("  EIGENVALUES  (N_cr of each buckling mode)"); print(line)
    pn = R["pname"]; ex = R.get("extra", [])
    hdr = f"  {'mode':>4} | {pn + ' N_cr [N/mm]':>18}{'m':>4}   {'mode type':<30}"
    for (en, _em, _et) in ex:
        hdr += f"| {en + ' N_cr':>12}{'m':>4}"
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for i, ((N, m), kd) in enumerate(zip(modes, kinds)):
        row = f"  {i + 1:>4} | {N:>18.4f}{m:>4}   {kd:<30}"
        for (en, em, _et) in ex:
            row += (f"| {em[i][0]:>12.4f}{em[i][1]:>4}" if i < len(em) else "|")
        print(row)
    for (en, em, et) in ex:
        print(f"  N_cr {pn} vs {en}: {100 * (N0 - em[0][0]) / em[0][0]:+.2f} %   "
              f"({pn} {R['t_main']:.1f} s, {en} {et:.1f} s)")
    print("\n" + line); print("  PANEL GEOMETRY AND PROPERTIES"); print(line)
    Lweb = float(np.hypot(cfg["spec"]["b3"] - info["ht"], cfg["spec"]["h"]))
    print(f"  skin width W = 2 bs = {info['W']:.1f} mm   foot width b1 = {info['b1']:.2f} mm")
    print(f"  omega: b3 = {cfg['spec']['b3']:g}  h = {cfg['spec']['h']:g}  "
          f"web angle = {cfg['spec']['theta']:g} deg  crown 2ht = {2 * info['ht']:.2f}  "
          f"web = {Lweb:.2f} mm")
    print(f"  t_skin = {info['ts']:g} mm   t_stringer = {info['tf']:g} mm   "
          f"stringer offset (tie) = {info['clear']:.4f} mm")
    print(f"  loaded edge length = {pr['Lload']:.1f} mm   area A = {pr['A']:.1f} mm^2")
    print(f"  EA = {pr['EA']:.4e} N   axial-stiffness centroid y_s = {pr['ys']:.2f} mm")
    print(f"  average axial stress at N_cr: sigma_cr = {N0 * pr['Lload'] / pr['A']:.2f} MPa")
    print("\n" + line); print("  LAMINATE (per part, about its own mid-plane; foot = skin + flange)"); print(line)
    seen = set()
    for (_i, _j, w, _As) in strips:
        if w["part"] in seen or w["part"] in ("crown", "link"):
            continue
        seen.add(w["part"])
        A_, B_, D_, t_ = w["ABD"]
        e = (-0.5 * info["ts"] + 0.5 * t_) if w["part"] == "foot" else 0.0
        Bm = B_ - e * A_; Dm = D_ - 2 * e * B_ + e * e * A_        # back to own mid-plane
        nm = {"skin": "skin", "foot": "foot", "web": "stringer"}[w["part"]]
        tag = []
        if abs(A_[0, 2]) > 1e-6 * A_[0, 0]: tag.append("unbalanced(A16)")
        if np.abs(Bm).max() > 1e-6 * A_[0, 0]: tag.append("unsymmetric(B)")
        print(f"  {nm:<9} t = {t_:6.3f}  A11 = {A_[0, 0]:9.0f}  A16 = {A_[0, 2]:8.0f}  "
              f"D11 = {Dm[0, 0]:8.1f}  D16 = {Dm[0, 2]:7.1f}  {' '.join(tag)}")
    if R["sweep"]:
        print("\n  N_cr(L) curve of the panel (same solver and mesh):")
        print(f"       L [mm] |  N_cr [N/mm]    m   mode                   skin %")
        print("     " + "-" * 64)
        for (L_, N_, m_, k_, sh_, *_p) in R["sweep"]:
            star = "  <- design" if abs(L_ - cfg["L"]) < 1e-9 else ""
            print(f"     {L_:>8.1f} | {N_:>12.4f}{m_:>5}   {k_:<22}{100 * sh_[0]:>5.0f}{star}")


def _omega_plots(cfg, R, outdir="out", tag="OMEGA"):
    """The same set of figures as the member sections, for the omega panel.
    Every curve here comes from the spline-strip solution itself."""
    import os
    import matplotlib.pyplot as plt
    os.makedirs(outdir, exist_ok=True)
    M, strips, modes, vecs = R["M"], R["strips"], R["modes"], R["vecs"]
    xs = np.linspace(0.0, M.L, 401)
    Us = [_omega_U(M, v, xs) for v in vecs]
    col = ["#c0392b", "#2471a3", "#1e8449", "#b9770e", "#7d3c98", "#117a65",
           "#a04000", "#2e4053"]
    pcol = {"skin": "#2471a3", "foot": "#1b2631", "web": "#c0392b", "crown": "#d68910",
            "link": "#7f8c8d"}
    pn = R["pname"]
    N0, m0 = modes[0]

    def xmax(U):
        tr = np.hypot(U[:, 1, :], U[:, 2, :])
        return np.unravel_index(np.argmax(tr), tr.shape)[1]

    # fig5: cross-section and strip mesh
    fig, ax = plt.subplots(figsize=(11, 3.6))
    for k, (iA, iB, w, _As) in enumerate(strips):
        ax.plot(M.nodes[[iA, iB], 0], M.nodes[[iA, iB], 1], "-", color=pcol[w["part"]],
                lw=2.2, solid_capstyle="butt")
    ax.plot(M.nodes[:, 0], M.nodes[:, 1], ".", color="0.35", ms=2.5)
    for s_, m_ in M.links:
        ax.plot(M.nodes[[s_, m_], 0], M.nodes[[s_, m_], 1], ":", color="k", lw=1.0)
    for nm, c in (("skin", pcol["skin"]), ("foot (skin + flange)", pcol["foot"]),
                  ("web", pcol["web"]), ("crown", pcol["crown"])):
        ax.plot([], [], color=c, lw=2.2, label=nm)
    ax.set_aspect("equal"); ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="upper right")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title(f"Cross-section and {pn} discretisation ({tag}): bs = {cfg['bs']:g} mm, "
                 f"{len(M.strips)} {'segments' if pn == 'GBT' else 'strips'} <= {R['ds']:.2f} mm")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig5_section_{tag}.png", dpi=200); plt.close()

    sw = R["sweep"]
    if sw:
        Ls = np.array([r[0] for r in sw]); Ns = np.array([r[1] for r in sw])
        sk = np.array([r[4][0] for r in sw]); st = 1.0 - sk
    # fig1: N_cr vs L
    plt.figure(figsize=(7.4, 4.6))
    if sw:
        plt.plot(Ls, Ns, "k-o", lw=2.2, ms=4, label=f"N_cr ({pn}, every point solved)")
        for r in sw:
            plt.annotate(f"m={r[2]}", (r[0], r[1]), textcoords="offset points",
                         xytext=(0, 7), ha="center", fontsize=7)
    plt.plot([cfg["L"]], [N0], "r*", ms=14, label=f"design length: {N0:.3f} N/mm")
    plt.xlabel("panel length L [mm]"); plt.ylabel(r"$N_{cr}$ [N/mm]")
    plt.title(f"Critical line load vs length ({tag}, {cfg['bc']})")
    plt.grid(True, alpha=0.3); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(f"{outdir}/fig1_Pcr_vs_L_{tag}.png", dpi=200); plt.close()

    # fig_stacked: N_cr(L) + displacement share skin / stringer
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 8.5), sharex=True)
    if sw:
        a1.plot(Ls, Ns, "k-o", lw=2.4, ms=4, label=f"FULL model ({pn})")
    a1.plot([cfg["L"]], [N0], "r*", ms=14, label="design length")
    a1.set_ylabel(r"$N_{cr}$ [N/mm]"); a1.grid(True, alpha=0.4, ls="--"); a1.legend(fontsize=9)
    a1.set_title(f"{pn} buckling analysis ({tag}-panel, {cfg['bc']})", fontweight="bold")
    if sw:
        a2.fill_between(Ls, 0, sk, color="#2471a3", alpha=0.65, label="skin (+ feet)")
        a2.fill_between(Ls, sk, 1, color="#c0392b", alpha=0.65, label="stringer (webs + crown)")
        a2.set_xlim(Ls.min(), Ls.max())
    a2.set_ylim(0, 1); a2.set_xlabel("Length L [mm]")
    a2.set_ylabel("share of the critical-mode displacement")
    a2.grid(True, alpha=0.4, ls="--")
    if sw:
        a2.legend(fontsize=9, loc="lower right")
    else:
        a2.text(0.5, 0.5, "no N_cr(L) curve requested", ha="center", va="center",
                transform=a2.transAxes, color="0.4")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig_stacked_{tag}.png", dpi=200); plt.close()

    # fig2: amplitude along the panel of the computed modes
    plt.figure(figsize=(8.0, 4.6))
    for k in range(min(3, len(Us))):
        plt.plot(xs, _omega_long_line(M, Us[k], xs), color=col[k], lw=1.8,
                 label=f"mode {k + 1}: N = {modes[k][0]:.3f}, m = {modes[k][1]}")
    plt.xlabel("z (along the stringers) [mm]"); plt.ylabel("normalised skin deflection")
    plt.title(f"Displacement amplitude along the panel ({cfg['bc']} ends, computed modes)")
    plt.grid(True, alpha=0.3); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(f"{outdir}/fig2_amplitude_{tag}.png", dpi=150); plt.close()

    # fig4: average axial stress at buckling vs L
    pr = R["props"]
    plt.figure(figsize=(7.2, 4.6))
    if sw:
        plt.plot(Ls, Ns * pr["Lload"] / pr["A"], "k-o", lw=2, ms=4)
    plt.plot([cfg["L"]], [N0 * pr["Lload"] / pr["A"]], "r*", ms=14)
    plt.xlabel("L [mm]"); plt.ylabel(r"$\sigma_{cr}$ (average axial) [MPa]")
    plt.title("Critical average axial stress vs length")
    plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(f"{outdir}/fig4_sigma_vs_L_{tag}.png", dpi=150); plt.close()

    # fig_critical_mode
    U = Us[0]; x0 = xmax(U)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 7.0), gridspec_kw=dict(height_ratios=[1.2, 1]))
    _omega_draw_section(a1, M, U, x0, col[0], title=f"cross-section at z = {xs[x0]:.0f} mm")
    a2.plot(xs, _omega_long_line(M, U, xs), "k", lw=1.8)
    a2.set_xlabel("z [mm]"); a2.set_ylabel("normalised skin deflection")
    a2.set_title(f"along the panel (line of max. deflection), m = {m0}"); a2.grid(True, alpha=0.3)
    fig.suptitle(f"Critical mode, {pn} ({cfg['bc']}):  N_cr = {N0:.4f} N/mm,  "
                 f"{R['kinds'][0]}", fontweight="bold")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig_critical_mode_{tag}.png", dpi=200); plt.close()

    # fig6: all computed buckling modes
    n = len(Us); ncol = 2; nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(12, 2.6 * nrow))
    axs = np.atleast_1d(axs).ravel()
    for k in range(n):
        _omega_draw_section(axs[k], M, Us[k], xmax(Us[k]), col[k % len(col)],
                            title=f"mode {k + 1}: N = {modes[k][0]:.3f} N/mm, m = {modes[k][1]}, "
                                  f"{R['kinds'][k]}")
    for k in range(n, len(axs)):
        axs[k].axis("off")
    fig.suptitle(f"Buckling modes of the panel ({tag}, {cfg['bc']}, L = {cfg['L']:g} mm)",
                 fontsize=13)
    plt.tight_layout(); plt.savefig(f"{outdir}/fig6_modeshapes_{tag}.png", dpi=170); plt.close()

    # fig_modes_interactive: the two modes chosen in the questionnaire
    fig, axs = plt.subplots(2, 2, figsize=(13, 6.2), gridspec_kw=dict(height_ratios=[1.2, 1]))
    for c_, k in enumerate(cfg.get("plot_modes", (1, 2))):
        k = k - 1
        if k >= n:
            for r_ in range(2):
                axs[r_, c_].text(0.5, 0.5, f"mode {k + 1} not computed\n(raise 'number of "
                                 "eigenvalues')", ha="center"); axs[r_, c_].axis("off")
            continue
        _omega_draw_section(axs[0, c_], M, Us[k], xmax(Us[k]), col[k % len(col)],
                            title=f"mode {k + 1}: N = {modes[k][0]:.3f} N/mm, m = {modes[k][1]}")
        axs[1, c_].plot(xs, _omega_long_line(M, Us[k], xs), color=col[k % len(col)], lw=1.8)
        axs[1, c_].set_xlabel("z [mm]"); axs[1, c_].grid(True, alpha=0.3)
        axs[1, c_].set_title(R["kinds"][k], fontsize=10)
    fig.suptitle(f"Deformed buckling mode shapes ({tag}-panel, {pn})", fontsize=13,
                 fontweight="bold")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig_modes_interactive_{tag}.png", dpi=200)
    plt.close()

    # fig_gbt_mode_participation: conventional GBT modes in the critical mode
    if R.get("mode_part"):
        mp = R["mode_part"][:15]
        plt.figure(figsize=(8.5, 3.8))
        plt.bar([str(k) for k, _v in mp], [100 * v for _k, v in mp], color="#1f5f8b")
        plt.xlabel("GBT mode number (hierarchical list)"); plt.ylabel("participation [%]")
        plt.title(f"Conventional GBT modes in the critical mode ({tag}-panel, {cfg['bc']}, "
                  f"L = {cfg['L']:g} mm)")
        plt.grid(True, axis="y", alpha=0.3); plt.tight_layout()
        plt.savefig(f"{outdir}/fig_gbt_mode_participation_{tag}.png", dpi=200); plt.close()
        _omega_draw_gbt_modes(R["section"], [k for k, v in R["mode_part"][:8]],
                              f"{outdir}/fig_gbt_critical_modes_{tag}.png",
                              f"GBT modes with the largest participation in the critical mode ({tag})")

    # fig_gbt_participation: GBT modal participation of the critical mode vs L
    if sw and all(len(r) > 5 and r[5] for r in sw):
        groups = [g for g in ("G", "D", "L", "S", "T", "M") if any(r[5].get(g, 0) > 0 for r in sw)]
        gc = {"G": "#3498db", "D": "#e67e22", "L": "#9b59b6", "S": "#2ecc71", "T": "#95a5a6",
              "M": "#e74c3c"}
        base = np.zeros(len(sw))
        plt.figure(figsize=(8.5, 3.8))
        for g in groups:
            v = np.array([r[5].get(g, 0.0) for r in sw])
            plt.fill_between(Ls, base, base + v, color=gc[g], alpha=0.75,
                             label=_GBT_GROUP_NAMES[g])
            base = base + v
        plt.ylim(0, 1); plt.xlim(Ls.min(), Ls.max())
        plt.xlabel("L [mm]"); plt.ylabel("GBT modal participation")
        plt.title(f"GBT modal participation of the critical mode ({tag}-panel, {cfg['bc']})")
        plt.legend(fontsize=8, loc="lower right"); plt.grid(True, alpha=0.3); plt.tight_layout()
        plt.savefig(f"{outdir}/fig_gbt_participation_{tag}.png", dpi=200); plt.close()


def run_omega_interactive(cfg):
    """Solve one omega panel (GBT and/or strip), print the result
    window and write the figures -- the same layout as the member sections.
    The window and figures are drawn from GBT when it is selected, otherwise
    from the strip; the other solvers appear as extra columns."""
    bs, L, spec, nm, bc = cfg["bs"], cfg["L"], cfg["spec"], cfg["nmodes"], cfg.get("bc", "SS")
    pan = cfg.get("panel") or make_omega_panel(bs, None, [], [], spec=spec)
    use = cfg["solver"].split("+")
    extra = []
    primary = "gbt" if "gbt" in use else "strip"
    solve = omega_gbt_buckling if primary == "gbt" else omega_spline_eigs
    t0 = time.perf_counter()
    modes, inf = solve(pan, L, bc, nmodes=nm, return_model=True)
    t_main = time.perf_counter() - t0
    if primary == "gbt" and "strip" in use:
        t0 = time.perf_counter()
        sm, _si = omega_spline_eigs(pan, L, bc, nmodes=nm)
        extra.append(("strip", sm, time.perf_counter() - t0))
    M, strips, vecs = inf["model"], inf["strips"], inf["vecs"]
    kinds = [_omega_spline_modeinfo(M, v)[2] for v in vecs]
    shares = [_omega_part_share(M, strips, v) for v in vecs]
    sweep = []
    if cfg.get("npt", 0) > 0 and cfg["Lmax"] > cfg["Lmin"] > 0:
        Lg = sorted(set([float(x) for x in np.linspace(cfg["Lmin"], cfg["Lmax"], cfg["npt"])]
                        + [float(L)]))
        print(f"\n  N_cr(L) curve ({'GBT' if primary == 'gbt' else 'strip'}): "
              f"solving {sum(abs(x - L) > 1e-9 for x in Lg)} more lengths "
              f"(the design length L = {L:g} mm is already solved) ...", flush=True)
        for Lk in Lg:
            if abs(Lk - L) < 1e-9:
                sweep.append((Lk, modes[0][0], modes[0][1], kinds[0], shares[0],
                              inf.get("participation"))); continue
            mk, ik = solve(pan, Lk, bc, nmodes=2, return_model=True)
            sweep.append((Lk, mk[0][0], mk[0][1], ik["mode_type"],
                          _omega_part_share(ik["model"], ik["strips"], ik["vecs"][0]),
                          ik.get("participation")))
            print(f"     L = {Lk:8.1f} mm   N_cr = {mk[0][0]:10.4f} N/mm   m = {mk[0][1]}",
                  flush=True)
    if primary == "gbt":
        _sec = inf["section"]
        solver_text = (f"shear-deformable GBT, modes {OMEGA_GBT_FAMILIES}: "
                       + (f"all {inf['nm']} modes" if not OMEGA_GBT_MODES else
                          f"{inf['n_conv']} of {inf['n_conv_all']} conventional modes + their shear "
                          f"fields ({inf['nm']} modes; pre-buckling with all modes)")
                       + f"; {_sec.N} cross-section nodes, spans {inf['hx']:.2f} mm, {inf['ndof']} DOF")
        pname = "GBT"
    else:
        mesh_name = "fine" if OMEGA_MESH_FINE else OMEGA_MESH
        solver_text = (f"shell-type Mindlin spline finite strip ({mesh_name} mesh: strips "
                       f"{inf['ds']:.2f} mm, spans {inf['hx']:.2f} mm, {inf['ndof']} DOF)")
        pname = "strip"
    R = dict(M=M, strips=strips, vecs=vecs, modes=modes, kinds=kinds, shares=shares,
             sweep=sweep, panel=pan, props=_omega_section_props(M, strips),
             ds=inf["ds"], hx=inf["hx"], solver_text=solver_text, pname=pname,
             extra=extra, t_main=t_main, part=inf.get("participation"),
             mode_part=inf.get("mode_participation"), section=inf.get("section"))
    _omega_report(cfg, R)
    tag = "OMEGA"; outdir = "out"
    try:
        _omega_plots(cfg, R, outdir, tag)
        import glob, os
        print(f"\nFigures written to {outdir} :")
        for f in sorted(glob.glob(os.path.join(outdir, f"*_{tag}.png"))):
            print("  ", f)
    except Exception as e:
        print(f"  (figures skipped: {e})")
    print(f"\n  {pname} solve at the design length: {t_main:.2f} s")
    print("\nDone.")
    return R


# ==============================================================================
# 18.  RUNTIME & HARDWARE REPORT
# ==============================================================================
def _system_report(elapsed_s):
    """Print the wall-clock runtime and a description of the host machine."""
    import os
    line = "=" * 78
    print("\n" + line); print("  RUNTIME & SYSTEM"); print(line)
    wait = _INPUT_WAIT[0]
    print(f"  computing time          : {elapsed_s - wait:.3f} s")
    print(f"  wall-clock runtime      : {elapsed_s:.3f} s"
          + (f"   (incl. {wait:.1f} s waiting for your answers)" if wait > 0.05 else ""))
    proc = platform.processor() or (os.environ.get("PROCESSOR_IDENTIFIER") or "n/a")
    print(f"  processor               : {proc}")
    print(f"  machine / architecture  : {platform.machine()}")
    print(f"  system / release        : {platform.system()} {platform.release()}")
    print(f"  logical CPU cores       : {os.cpu_count()}")
    # physical RAM  (portable best-effort; silently skipped where unavailable)
    try:
        gib = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024**3)
        print(f"  physical RAM            : {gib:.1f} GiB")
    except (ValueError, AttributeError, OSError):
        pass
    # nicer CPU model name on Linux, if /proc exposes it
    try:
        with open("/proc/cpuinfo") as fh:
            for ln in fh:
                if ln.lower().startswith("model name"):
                    print(f"  CPU model               : {ln.split(':', 1)[1].strip()}")
                    break
    except OSError:
        pass
    print(f"  python                  : {platform.python_version()} "
          f"({platform.python_implementation()})")
    try:
        import scipy as _scipy
        print(f"  numpy / scipy           : {np.__version__} / {_scipy.__version__}")
    except Exception:
        print(f"  numpy                   : {np.__version__}")
    print(line)


def run_parametric_RHS(csv_out="RHS_section_GBT_results.csv"):
    """Parametric study of the RHS mirroring the Abaqus script
    (same geometry, materials, lay-ups, BCs and lengths)."""
    return _section_study("box", "RHS (height 50, width 20, t=1)", csv_out)


def run_parametric_L(csv_out="L_section_GBT_results.csv"):
    """Parametric study of the L-section (angle) mirroring the Abaqus script
    (same geometry, materials, lay-ups, BCs and lengths)."""
    return _section_study("l", "L-section (legs 50 and 20, t=1)", csv_out)


def run_parametric_Z(csv_out="Z_section_GBT_results.csv"):
    """Parametric study of the Z-section mirroring the Abaqus script
    (same geometry, materials, lay-ups, BCs and lengths)."""
    return _section_study("z", "Z-section (web 50, flange 20, t=1)", csv_out)


def run_parametric_C(csv_out="C_section_GBT_results.csv"):
    """Parametric study of the C-section mirroring the Abaqus script
    (same geometry, materials, lay-ups, BCs and lengths)."""
    return _section_study("c", "C-section (web 50, flange 20, t=1)", csv_out)


def run_parametric_T(csv_out="T_section_GBT_results.csv"):
    """Parametric study of the T-section mirroring the Abaqus script
    (same geometry, materials, lay-ups, BCs and lengths)."""
    return _section_study("t", "T-section (web 50, flange 20, t=1)", csv_out)


def run_parametric_I(csv_out="I_section_GBT_results.csv"):
    """Parametric study of the I-section mirroring the Abaqus script
    (same geometry, materials, lay-ups, BCs and lengths)."""
    return _section_study("i", "I-section (web 50, flange 20, t=1)", csv_out)


# ==============================================================================
# 19.  UNIFIED WHOLE-SECTION MINDLIN FINITE STRIP  (harmonic; P_cr(L) curve of the strip solver)
# ------------------------------------------------------------------------------
# ONE formulation for every cross-section (c, l, z, lipc, t, i, box, hat, dhat)
# and both end conditions, with NO case rules:
#   * every wall is a laminated Mindlin/FSDT plate (full [A B; B D] + As, 5/6);
#     theory="CLPT" gives the Kirchhoff limit of the SAME model (As -> rigid);
#   * the walls are joined at the fold lines with full displacement AND rotation
#     continuity (6 DOF per node), so local, distortional and global buckling
#     and all their interactions come out of one eigenproblem -- no separate
#     "local plate" estimate, no plateau, no L/h switch;
#   * the end conditions enter ONLY through the longitudinal functions
#     (SS: sin(m pi x/L);  CC: sin(m pi x/L) sin(pi x/L), warping restrained)
#     -- no effective length;
#   * the pre-buckling state is computed (generalised plane strain of the whole
#     section: axial strain + in-plane deformation + shear flows of closed cells
#     + the A16/B coupling responses), not assumed;
#   * harmonics are coupled exactly (numerical x-integration of every product).
# The kernel (_mfs_*) is the one verified in 14e against the exact Mindlin
# plate solution and against the shell FE.
# ==============================================================================
def _usec_walls(cfg):
    """Natural nodes and walls [(i, j, plies)] of the configured section, with
    the SAME wall directions (hence laminate normals) as the rest of the code."""
    mat = cfg["material"]
    if cfg.get("corner_nodes") is not None:                          # family A (ordered path)
        cn = [tuple(map(float, p)) for p in cfg["corner_nodes"]]
        lay = list(cfg.get("layup", [0.0]))
        walls = []
        for w in range(len(cn) - 1):
            t = float(cfg["t_wall"][w])
            plies = ([(0.0, t)] if mat["type"] == "isotropic"
                     else [(float(a), t / len(lay)) for a in lay])
            walls.append((w, w + 1, plies))
        return np.asarray(cn, float), walls
    nodes = np.asarray(cfg["nodes"], float)
    walls = []
    for w in cfg["walls"]:
        plies = [(float(a), float(t)) for a, t in w["plies"]]
        if mat["type"] == "isotropic":
            plies = [(0.0, sum(t for _a, t in plies))]
        walls.append((int(w["i"]), int(w["j"]), plies))
    return nodes, walls


def _usec_model(cfg, ds=None, theory="FSDT"):
    """Strip model: nodes, strips [(iA, iB, ABD, As)], wall id per strip."""
    mat = cfg["material"]
    nodes0, walls = _usec_walls(cfg)
    walls = _usec_sketch_direction(cfg, nodes0, walls)   # Abaqus shell normals
    cen = _usec_centroid(nodes0, walls)
    lens = [float(np.hypot(*(nodes0[j] - nodes0[i]))) for i, j, _p in walls]
    if ds is None:                                       # <= 1/16 of the narrowest wall
        ds = min(2.5, min(lens) / 16.0)
    pts = [tuple(p) for p in nodes0]; strips = []; wid = []
    for k, (i, j, plies) in enumerate(walls):
        A, B, D, _t = clt_plies(mat, plies)
        As = _osfe_laminate(mat, plies)[3]
        if theory == "CLPT":                             # Kirchhoff limit
            As = As * 1.0e6
        P0, P1 = nodes0[i], nodes0[j]
        n = max(2, int(np.ceil(lens[k] / ds)))
        prev = i
        for q in range(1, n + 1):
            if q < n:
                pts.append(tuple(P0 + (P1 - P0) * q / n)); nxt = len(pts) - 1
            else:
                nxt = j
            # strip axis = reverse of the sketch direction (normal = Abaqus normal)
            if _usec_normal_flip(P0, P1, cen):
                strips.append((prev, nxt, (A, B, D), As))
            else:
                strips.append((nxt, prev, (A, B, D), As))
            wid.append(k); prev = nxt
    return np.asarray(pts, float), strips, wid


def _usec_prestate(nodes, strips, bc):
    """Generalised plane-strain pre-buckling state for a unit TOTAL axial
    compression (1 N) applied through the stiffness centroid: uniform axial
    strain e0 (no member bending), full in-plane response of every wall
    (Poisson, A16/B couplings, closed-cell shear flows).  Returns the
    compression-positive (Px, Ps, Pxs) of each strip per 1 N."""
    N = len(nodes); g0 = 6 * N; nd = g0 + 1
    K = np.zeros((nd, nd)); ops = []
    for (iA, iB, ABD, As) in strips:
        d = nodes[iB] - nodes[iA]; b = float(np.hypot(*d)); al = np.arctan2(d[1], d[0])
        C = _mfs_C9(ABD, As)
        T = np.zeros((12, 12)); T[:6, :6] = _mfs_T6(al); T[6:, 6:] = _mfs_T6(al)
        dofs = list(range(6 * iA, 6 * iA + 6)) + list(range(6 * iB, 6 * iB + 6)) + [g0]
        Ke = np.zeros((13, 13)); Eg = []
        for g in _MFS_G2:
            r = 0.5 * (g + 1.0)
            E = np.zeros((9, 13))
            for (comp, fld, sop, fam, cf) in _MFS_TERMS:
                if fam in (1, 3):
                    continue
                E[comp, :12] += cf * _mfs_svec(fld, sop, 0.5 if comp == 7 else r, b)
            E[:, :12] = E[:, :12] @ T; E[0, 12] = 1.0
            Ke += 0.5 * b * (E.T @ C @ E); Eg.append(E)
        K[np.ix_(dofs, dofs)] += Ke; ops.append((dofs, Eg, C, b))
    # remove the in-plane rigid-body motions of the section and the axial one
    fixed = {0, 1, 2, 3}                                 # node 0: u, UX, UY, Rx
    far = int(np.argmax(np.hypot(*(nodes - nodes[0]).T)))
    d = nodes[far] - nodes[0]                            # second node: kill the rotation
    fixed.add(6 * far + (1 if abs(d[1]) > abs(d[0]) else 2))
    free = np.array([k for k in range(g0) if k not in fixed])
    u = np.zeros(nd); u[g0] = -1.0                       # unit shortening
    rhs = -K[np.ix_(free, [g0])] @ u[[g0]]
    u[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)
    P = []; F = 0.0
    for (dofs, Eg, C, b) in ops:
        Nm = np.mean([(C @ (E @ u[dofs]))[:3] for E in Eg], axis=0)
        P.append(-Nm); F += -Nm[0] * b
    return [tuple(p / F) for p in P]


def _usec_modeclass(nodes, vec, M):
    """Global / distortional / local character of a mode: share of the cross-
    section in-plane displacement that is a rigid-body motion of the section
    (2 translations + rotation, least squares over all nodes and harmonics)."""
    full = vec.reshape(-1, 6, M)
    X, Y = nodes[:, 0] - nodes[:, 0].mean(), nodes[:, 1] - nodes[:, 1].mean()
    tot = 0.0; rig = 0.0
    for m in range(M):
        ux, uy = full[:, 1, m], full[:, 2, m]
        A = np.zeros((2 * len(X), 3))
        A[0::2, 0] = 1.0; A[1::2, 1] = 1.0; A[0::2, 2] = -Y; A[1::2, 2] = X
        b = np.empty(2 * len(X)); b[0::2] = ux; b[1::2] = uy
        c, *_ = np.linalg.lstsq(A, b, rcond=None)
        tot += b @ b; rig += (A @ c) @ (A @ c)
    r = rig / max(tot, 1e-300)
    kind = "global" if r > 0.8 else ("local" if r < 0.2 else "distortional/mixed")
    return r, kind


def unified_section_solve(cfg, L, bc="SS", theory="FSDT", window=6, ds=None,
                          mmax=None, return_mode=False):
    """Critical axial load P_cr [N] of the member (cross-section of cfg, length
    L, end condition bc) by the unified shear-deformable finite strip."""
    nodes, strips, wid = _usec_model(cfg, ds, theory)
    P = _usec_prestate(nodes, strips, bc)
    S = [(iA, iB, ABD, As, P[k]) for k, (iA, iB, ABD, As) in enumerate(strips)]
    nat, walls = _usec_walls(cfg)
    wmin = min(float(np.hypot(*(nat[j] - nat[i]))) for i, j, _p in walls)
    if mmax is None:                                     # half-waves down to wmin/2
        mmax = int(np.clip(np.ceil(2.0 * L / wmin) + 2, 8, 160))
    def run(ml, k):
        Kr, Gr, C = _mfs_assemble(nodes, S, L, ml, bc)
        lam, V = _mfs_buckle(Kr, Gr, k)
        return lam, V, C
    # lowest single-harmonic solution: coarse (geometric) scan of m, then every
    # m between the neighbours of the best one (P(m) is smooth in m)
    single = {}
    def pick(m):
        if m not in single:
            single[m] = run([m], 1)
        return single[m][0][0]
    grid = sorted(set(np.round(np.geomspace(1, mmax, 24)).astype(int).tolist()))
    for m in grid:
        pick(m)
    for mb in sorted(single, key=lambda m: single[m][0][0])[:2]:
        k = grid.index(mb)
        for m in range(grid[max(0, k - 1)], grid[min(len(grid) - 1, k + 1)] + 1):
            pick(m)
    coupled = bc == "CC"
    if not coupled:
        for (iA, iB, ABD, As, p) in S:
            A, B, D = ABD
            if (max(abs(A[0, 2]), abs(A[1, 2]), abs(B[0, 2]), abs(B[1, 2]), abs(D[0, 2]),
                    abs(D[1, 2]), abs(As[0, 1])) > 1e-8 * max(1.0, A[0, 0])
                    or abs(p[2]) > 1e-12):
                coupled = True; break
    mstar = min(single, key=lambda m: single[m][0][0])
    if not coupled:
        lam, V, C = single[mstar]; win = [mstar]
    else:
        win = list(range(max(1, mstar - window), min(mmax, mstar + window) + 1))
        lam, V, C = run(win, 1)
    Mw = len(win); vec = C @ V[:, 0]
    full = vec.reshape(-1, Mw)
    mdom = win[int(np.argmax(np.linalg.norm(full, axis=0)))]
    r, kind = _usec_modeclass(nodes, vec, Mw)
    out = dict(P=float(lam[0]), m=mdom, kind=kind, rigid_share=r, window=win,
               theory=theory, bc=bc, ndof_per_harmonic=6 * len(nodes))
    if return_mode:
        out.update(nodes=nodes, vec=vec, mlist=win)
    return out


# ==============================================================================
# 20.  SPLINE FINITE STRIP  (shear-deformable, general end conditions)
# ------------------------------------------------------------------------------
# The unified solver for every section.  Across the section: the Mindlin strip
# of section 14e (same strain terms _MFS_TERMS, same laminate [A B; B D] + As,
# 6 DOF/node, rigid links, assumed gamma_sn).  Along the member: cubic
# B-splines on an open uniform knot vector (the spline finite strip method,
# Cheung & Fan; Lau & Hancock for thin-walled members).  Because the basis is
# local and interpolates the end values, the ACTUAL support conditions are
# imposed node by node (e.g. U1 = U2 = 0 on the end edges, rigid end plates),
# and the pre-buckling state is the strip's own linear static solution under
# the actual loads -- so the load-introduction end zones (Poisson restraint of
# gripped/clamped ends, shear-flow build-up in closed cells) are part of the
# geometric stiffness.  No effective length, no uniform-stress assumption, no
# case rules.
# ==============================================================================
_SFS_GX4, _SFS_WX4 = np.polynomial.legendre.leggauss(4)


def _sfs_basis(L, nspan, breaks=None):
    """Open cubic B-splines on [0, L] (uniform spans, or the given breakpoints):
    number of basis functions, and for every span the values / x-derivatives of
    its 4 non-zero splines at the 4 Gauss points (+ weights and points)."""
    from scipy.interpolate import BSpline
    if breaks is None:
        breaks = np.linspace(0.0, L, nspan + 1)
    breaks = np.asarray(breaks, float); nspan = len(breaks) - 1
    t = np.r_[[0.0] * 3, breaks, [L] * 3]
    nb = nspan + 3
    hj = np.diff(breaks)
    xs = breaks[:-1, None] + 0.5 * hj[:, None] * (_SFS_GX4[None, :] + 1.0)
    wx = 0.5 * hj[:, None] * _SFS_WX4[None, :]
    Bv = np.zeros((nspan, 4, 4)); Bd = np.zeros((nspan, 4, 4))
    for i in range(nb):
        c = np.zeros(nb); c[i] = 1.0
        spl = BSpline(t, c, 3); dspl = spl.derivative()
        for j in range(max(0, i - 3), min(nspan, i + 1)):
            Bv[j, :, i - j] = spl(xs[j]); Bd[j, :, i - j] = dspl(xs[j])
    return nb, xs, wx, Bv, Bd


def _sfs_breaks(L, h_int, h_end, zone):
    """Span breakpoints: spans of h_end within 'zone' of each end (the
    load-introduction / support zones, where the pre-buckling state varies and
    end-localised modes live), h_int in the interior."""
    L = float(L); zone = min(float(zone), 0.5 * L)
    ne = max(1, int(np.ceil(zone / h_end)))
    ni = int(np.ceil((L - 2 * zone) / h_int)) if L - 2 * zone > 1e-9 else 0
    b = np.r_[np.linspace(0.0, zone, ne + 1)]
    if ni:
        b = np.r_[b, np.linspace(zone, L - zone, ni + 1)[1:]]
    b = np.r_[b, np.linspace(L - zone, L, ne + 1)[1:]]
    b = np.unique(np.round(b, 12))
    b[0], b[-1] = 0.0, L            # exact ends: rounding must not move them past L
    return b


def _sfs_strip_sparts(b, ABD, As, al, bubble=True):
    """s-parts of one strip for the two s-Gauss points: strain operators split
    by x-operator (0 = value, 1 = x-derivative) Es[sg][o] (9 x nloc) and
    geometric operators Gs[sg][c][o][o'] (nloc x nloc) for unit Px (c=0),
    Ps (c=1), Pxs (c=2).  Element fields: the 12 nodal DOFs (global axes,
    via T) and, with bubble=True, two strip-internal membrane modes u_b, v_b
    with the hierarchical quadratic shape phi = 4 r (1 - r) (zero at both
    fold lines, so the strips stay conforming).  They let eps_s and gamma_xs
    vary linearly across the strip, which removes the Poisson / in-plane
    bending over-stiffness of linear membrane strips (largest for +-45 plies)."""
    nl = 14 if bubble else 12
    T = np.zeros((12, 12)); T[:6, :6] = _mfs_T6(al); T[6:, 6:] = _mfs_T6(al)
    xop = {0: 0, 1: 1, 2: 0, 3: 1}                          # family -> x-operator
    Es = []; Gs = []
    for g in _MFS_G2:
        r = 0.5 * (g + 1.0)
        E0 = [np.zeros((9, 12)), np.zeros((9, 12))]
        for (comp, f, sop, fam, cf) in _MFS_TERMS:
            E0[xop[fam]][comp] += cf * _mfs_svec(f, sop, 0.5 if comp == 7 else r, b)
        E = [np.zeros((9, nl)), np.zeros((9, nl))]
        E[0][:, :12] = E0[0] @ T; E[1][:, :12] = E0[1] @ T
        ph, phs = 4.0 * r * (1.0 - r), (4.0 - 8.0 * r) / b
        if bubble:                                           # u_b = col 12, v_b = col 13
            E[1][0, 12] += ph                                # eps_x   = u,x
            E[0][1, 13] += phs                               # eps_s   = v,s
            E[0][2, 12] += phs; E[1][2, 13] += ph            # gam_xs  = u,s + v,x
            E[0][8, 12] += 0.5 * phs; E[1][8, 13] -= 0.5 * ph   # drilling penalty
        Es.append(E)
        G = [[[np.zeros((nl, nl)) for _ in range(2)] for _ in range(2)] for _ in range(3)]
        for f in (0, 1, 2):                                  # u, v, w
            vx = np.zeros(nl); vs = np.zeros(nl)
            vx[:12] = _mfs_svec(f, 0, r, b) @ T              # d/dx  -> x-op 1
            vs[:12] = _mfs_svec(f, 1, r, b) @ T              # d/ds  -> x-op 0
            if bubble and f < 2:
                vx[12 + f] = ph; vs[12 + f] = phs
            G[0][1][1] += np.outer(vx, vx)
            G[1][0][0] += np.outer(vs, vs)
            G[2][1][0] += np.outer(vx, vs); G[2][0][1] += np.outer(vs, vx)
        Gs.append(G)
    return Es, Gs


class SplineStripModel:
    """Prismatic folded-plate member discretised by Mindlin strips (across) and
    cubic B-splines (along).  nodes (N,2) cross-section; strips [(iA, iB, ABD,
    As)]; links [(slave, master)] rigid offset links; L length; nspan spans
    (or breaks = span breakpoints).  DOF numbering: field * nb + spline, with
    fields 6 n + d for the node DOFs (d = u, UX, UY, R about x, RX, RY) and
    6 N + 2 k + j for the membrane bubbles of strip k (j = 0 u_b, 1 v_b).
    Constraint keys: (node, d, spline) and ("B", k, j, spline)."""

    def __init__(self, nodes, strips, L, nspan, links=(), breaks=None, bubble=True):
        self.nodes = np.asarray(nodes, float); self.N = len(self.nodes)
        self.strips = strips; self.L = float(L); self.bubble = bool(bubble)
        if breaks is None:
            breaks = np.linspace(0.0, float(L), int(nspan) + 1)
        self.breaks = np.asarray(breaks, float); self.nspan = len(self.breaks) - 1
        self.links = list(links)
        self.nb, self.xs, self.wx, self.Bv, self.Bd = _sfs_basis(L, self.nspan, self.breaks)
        self.nl = 14 if self.bubble else 12
        self.nf = 6 * self.N + (2 * len(strips) if self.bubble else 0)
        self.ndof = self.nf * self.nb
        self.S = []
        for k, (iA, iB, ABD, As) in enumerate(strips):
            d = self.nodes[iB] - self.nodes[iA]
            b = float(np.hypot(*d)); al = np.arctan2(d[1], d[0])
            Es, Gs = _sfs_strip_sparts(b, ABD, As, al, self.bubble)
            C9 = _mfs_C9(ABD, As)
            loc = [6 * n + j for n in (iA, iB) for j in range(6)]
            if self.bubble:
                loc += [6 * self.N + 2 * k, 6 * self.N + 2 * k + 1]
            self.S.append(dict(b=b, al=al, Es=Es, Gs=Gs, C9=C9, loc=np.array(loc),
                               ABD=ABD, nodes=(iA, iB)))
        Xo = (self.Bv, self.Bd)
        self.Xg = [[np.einsum("jg,jga,jgb->jab", self.wx, Xo[o], Xo[p]) for p in range(2)]
                   for o in range(2)]
        self._sp = np.arange(self.nspan)[:, None] + np.arange(4)[None, :]

    def node_part(self, vec):
        """(N, 6, nb) view of the node DOFs of a full vector."""
        return np.asarray(vec)[:6 * self.N * self.nb].reshape(self.N, 6, self.nb)

    def _edofs(self, st):
        return (st["loc"][:, None] * self.nb)[None, :, :] + self._sp[:, None, :]   # (nspan,nl,4)

    def _assemble(self, blocks_fn):
        """Sum of the per-strip span blocks, assembled in chunks (memory)."""
        K = None; rows, cols, vals = [], [], []; cnt = 0; ne = 4 * self.nl
        def flush():
            nonlocal K, rows, cols, vals, cnt
            if not vals:
                return
            A = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(self.ndof, self.ndof)).tocsr()
            K = A if K is None else K + A
            rows, cols, vals, cnt = [], [], [], 0
        for st in self.S:
            Ke = blocks_fn(st)                                    # (nspan, ne, ne)
            ed = self._edofs(st).reshape(self.nspan, ne).astype(np.int32)
            rows.append(np.repeat(ed, ne, axis=1).ravel())
            cols.append(np.tile(ed, (1, ne)).ravel()); vals.append(Ke.ravel())
            cnt += Ke.size; del Ke
            if cnt > 8_000_000:
                flush()
        flush()
        return K

    def _xband(self, X):
        """nb x nb banded matrix from the span blocks X (nspan, 4, 4)."""
        r = np.broadcast_to(self._sp[:, :, None], X.shape)
        c = np.broadcast_to(self._sp[:, None, :], X.shape)
        return sp.coo_matrix((X.ravel(), (r.ravel(), c.ravel())),
                             shape=(self.nb, self.nb)).tocsr()

    def stiffness(self):
        """K = sum_(o,p) Kcs_op (x) X_op : the laminate is uniform along the
        member, so the stiffness is the Kronecker product of a cross-section
        matrix (per pair of x-operators) and a longitudinal spline Gram."""
        K = None; nl = self.nl
        for o in range(2):
            for p in range(2):
                r, c, v = [], [], []
                for st in self.S:
                    A = np.zeros((nl, nl))
                    for sg in range(2):
                        A += 0.5 * st["b"] * (st["Es"][sg][o].T @ st["C9"] @ st["Es"][sg][p])
                    if np.any(A):
                        r.append(np.repeat(st["loc"], nl)); c.append(np.tile(st["loc"], nl))
                        v.append(A.ravel())
                if not v:
                    continue
                Kcs = sp.coo_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))),
                                    shape=(self.nf, self.nf)).tocsr()
                T = sp.kron(Kcs, self._xband(self.Xg[o][p]), format="csr")
                K = T if K is None else K + T
        return K

    def resultants(self, u):
        """N = A eps + B kappa (tension +) at every Gauss point:
        list per strip of (nspan, 2 s-gp, 4 x-gp, 3)."""
        out = []; Xo = (self.Bv, self.Bd)
        for st in self.S:
            ue = u[self._edofs(st)]                               # (nspan, nl, 4)
            N = np.zeros((self.nspan, 2, 4, 3))
            for sg in range(2):
                eps = np.zeros((self.nspan, 4, 9))
                for o in range(2):
                    q = np.einsum("jag,jxg->jxa", ue, Xo[o])     # (nspan, 4x, nl)
                    eps += np.einsum("ka,jxa->jxk", st["Es"][sg][o], q)
                N[:, sg] = np.einsum("ck,jxk->jxc", st["C9"][:3, :6], eps[..., :6])
            out.append(N)
        return out

    def geometric(self, Nlist):
        """Kg from the pre-buckling resultants (compression = -N)."""
        Xo = (self.Bv, self.Bd); it = iter(Nlist); ne = 4 * self.nl
        def blk(st):
            N = next(it); Kg = np.zeros((self.nspan, ne, ne))
            for sg in range(2):
                wsg = 0.5 * st["b"]
                for c in range(3):
                    P = -N[:, sg, :, c]                           # (nspan, 4x) compression +
                    if not np.any(P):
                        continue
                    for o in range(2):
                        for p in range(2):
                            G = st["Gs"][sg][c][o][p]
                            if not np.any(G):
                                continue
                            Xw = np.einsum("jg,jg,jga,jgb->jab", self.wx, P, Xo[o], Xo[p])
                            Kg += wsg * np.einsum("ab,jcd->jacbd", G, Xw).reshape(self.nspan, ne, ne)
            return Kg
        return self._assemble(blk)

    def _field(self, key):
        if key[0] == "B":
            return 6 * self.N + 2 * key[1] + key[2], key[3]
        return 6 * key[0] + key[1], key[2]

    def reduction(self, fixed, merge=()):
        """full = C q.  fixed: set of constraint keys (see class doc); merge:
        groups of keys made equal (first = master); rigid links applied."""
        slave = {s: m for s, m in self.links}
        fixed = {self._field(k) for k in fixed}
        owner = {}
        for g in merge:
            m0 = self._field(g[0])
            for k in g[1:]:
                owner[self._field(k)] = m0
        qidx = {}; nq = 0
        for f in range(self.nf):
            if f < 6 * self.N and (f // 6) in slave:
                continue
            for i in range(self.nb):
                key = (f, i)
                if key in fixed or key in owner:
                    continue
                qidx[key] = nq; nq += 1
        rr, cc, vv = [], [], []
        def add(row, key, v):
            key = owner.get(key, key)
            if key in qidx:
                rr.append(row); cc.append(qidx[key]); vv.append(v)
        for n in range(self.N):
            m_ = slave.get(n, n)
            rX = rY = 0.0
            if n in slave:
                rX, rY = self.nodes[n] - self.nodes[m_]
            for i in range(self.nb):
                for d in range(6):
                    add((6 * n + d) * self.nb + i, (6 * m_ + d, i), 1.0)
                if n in slave:
                    add((6 * n + 1) * self.nb + i, (6 * m_ + 3, i), -rY)
                    add((6 * n + 2) * self.nb + i, (6 * m_ + 3, i), rX)
                    add((6 * n + 0) * self.nb + i, (6 * m_ + 4, i), rY)
                    add((6 * n + 0) * self.nb + i, (6 * m_ + 5, i), -rX)
        for f in range(6 * self.N, self.nf):
            for i in range(self.nb):
                add(f * self.nb + i, (f, i), 1.0)
        C = sp.csr_matrix((vv, (rr, cc)), shape=(self.ndof, nq))
        # optional extra unknowns: each 'jump' is one global DOF J added to the
        # listed (field, spline) rows with the given coefficient, i.e.
        # u_slave = u_master + coeff * J  (periodic cell with a free macro strain)
        for rows in (getattr(self, "jumps", None) or []):
            r = [f * self.nb + i for f, i, _c in rows]; v = [c_ for _f, _i, c_ in rows]
            C = sp.hstack([C, sp.csr_matrix((v, (r, [0] * len(r))),
                                            shape=(self.ndof, 1))]).tocsr()
        return C

    def end_line_load(self, q_of_strip, bubble=True):
        """Consistent end loads of a compressive end line load q (N/mm) on both
        end edges of every strip (q_of_strip[k]): q b / 2 on each fold node and
        q b 2/3 on the strip's u-bubble."""
        f = np.zeros(self.ndof); nb = self.nb
        for k, st in enumerate(self.S):
            q = q_of_strip[k]
            for n in st["nodes"]:
                f[(6 * n) * nb + 0] += q * st["b"] / 2.0
                f[(6 * n) * nb + nb - 1] -= q * st["b"] / 2.0
            if self.bubble and bubble:
                fb = (6 * self.N + 2 * k) * nb
                f[fb + 0] += q * st["b"] * 2.0 / 3.0
                f[fb + nb - 1] -= q * st["b"] * 2.0 / 3.0
        return f

    def end_bubble_keys(self, which=(0, 1), strips=None):
        """Constraint keys of the strip bubbles at both ends (j in which)."""
        ks = range(len(self.S)) if strips is None else strips
        if not self.bubble:
            return set()
        return {("B", k, j, i) for k in ks for j in which for i in (0, self.nb - 1)}

    def buckle(self, fixed, f, merge=(), nev=2, return_modes=False, tol=1e-10):
        C = self.reduction(fixed, merge)
        K = self.stiffness(); Kr = (C.T @ K @ C).tocsc(); Kr = 0.5 * (Kr + Kr.T)
        solve, release = _osfe_factor(Kr)
        uq = solve(C.T @ f); u = C @ uq
        N = self.resultants(u)
        Kg = self.geometric(N); Gr = (C.T @ Kg @ C).tocsc()          # compression-positive Kg
        Gr = 0.5 * (Gr + Gr.T)
        Kinv = spla.LinearOperator(Kr.shape, matvec=solve, dtype=float)
        mu, V = spla.eigsh(Gr, k=nev, M=Kr, Minv=Kinv, which="LA",
                           ncv=max(2 * nev + 16, 20), tol=tol, maxiter=20000)
        release()
        ok = mu > 0; lam = 1.0 / mu[ok]; V = V[:, ok]; o = np.argsort(lam)
        lam, V = lam[o], V[:, o]
        out = dict(u=u, N=N)
        if return_modes:
            out["modes"] = [C @ V[:, i] for i in range(V.shape[1])]
        return lam, out


# Shell normal of the section walls: decides on which side ply 1 lies and the
# sense of the ply angles -- it matters ONLY for unsymmetric lay-ups.
#   "outward" : every wall's normal points away from the section centroid
#               (walls through the centroid: n = t x Z).  This is the side of
#               the laminate your Abaqus section models reproduce: with it the
#               30 Unsym cases of C, Z, T, I and RHS are within about 2 %
#               (mean ~1 %), e.g. the I-section -- whose result is the same for
#               any global flip, so a genuine test -- goes from +2.0 % to
#               +-0.6 %.  Please confirm on one model in Abaqus/CAE
#               (Query > Shell element normals: arrows should point outward).
#   "t x Z"   : right-hand rule with the sketch line direction t
#   "Z x t"   : the opposite of "t x Z" for every wall
# It is a statement of the model's geometry, not a fitting parameter.
SECTION_NORMAL_RULE = "outward"

# Mesh of the spline strip (mesh study on the 120 Abaqus cases):
#   default (None): strips min(2.5, 0.025 h) mm, spans min(5, L/40) mm
#                   -> mean |error| vs Abaqus 0.56 %, ~18 s per case
#   --mesh-fast   : strips 2.5 mm, spans 10 mm  -> 0.93 %, ~5 s per case
#   --mesh-fine   : strips 0.625 mm, spans L/80 mm -> 0.51 %, ~40 s per case
MESH_STRIP = None
MESH_SPAN = None


def _usec_normal_flip(P0, P1, cen):
    """True if the wall P0 -> P1 (sketch direction) must use the normal
    opposite to t x Z under SECTION_NORMAL_RULE."""
    if SECTION_NORMAL_RULE == "t x Z":
        return False
    if SECTION_NORMAL_RULE == "Z x t":
        return True
    t = (np.asarray(P1, float) - np.asarray(P0, float)); t = t / np.hypot(*t)
    n_tz = np.array([t[1], -t[0]])                          # t x Z
    d = (0.5 * (np.asarray(P0) + np.asarray(P1)) - cen) @ n_tz
    return bool(d < -1e-6)


def _usec_centroid(nodes0, walls):
    """Length-weighted centroid of the wall mid-lines."""
    L = [float(np.hypot(*(nodes0[j] - nodes0[i]))) for i, j, _p in walls]
    return sum(0.5 * (nodes0[i] + nodes0[j]) * l for (i, j, _p), l in zip(walls, L)) / sum(L)


def _usec_sketch_direction(cfg, nodes, walls):
    """Orient every wall as its line is drawn in the Abaqus parametric scripts,
    because the shell normal (n = sketch direction x Z, right-hand rule) and
    therefore the ply stacking side and the sense of the ply angle follow it:
      c, z : flanges from the web junction to the free tip, web top -> bottom
      t, i : flanges left -> right, web top -> bottom
      box  : clockwise, starting at the top-left corner
    Other shapes keep the wall direction of make_section / corner_nodes."""
    shape = str(cfg.get("shape", "")).lower()
    from collections import Counter
    val = Counter()
    for i, j, _p in walls:
        val[i] += 1; val[j] += 1
    out = []
    cx, cy = nodes[:, 0].mean(), nodes[:, 1].mean()
    for i, j, plies in walls:
        P, Q = nodes[i], nodes[j]
        vertical = abs(Q[0] - P[0]) < 1e-9
        if shape in ("c", "z"):
            if vertical:                                   # web: top -> bottom
                if P[1] < Q[1]: i, j = j, i
            elif val[i] == 1:                              # flange: junction -> tip
                i, j = j, i
        elif shape in ("t", "i"):
            if vertical:
                if nodes[i][1] < nodes[j][1]: i, j = j, i
            elif nodes[i][0] > nodes[j][0]:
                i, j = j, i
        elif shape == "box":                               # clockwise
            mid = 0.5 * (nodes[i] + nodes[j]); t = nodes[j] - nodes[i]
            r = mid - np.array([cx, cy])
            if r[0] * t[1] - r[1] * t[0] > 0:              # counter-clockwise -> flip
                i, j = j, i
        out.append((i, j, plies))
    return out


def _usec_build(cfg, L, theory="FSDT", ds=None, hx=None):
    """Spline-strip model of the member of cfg (walls oriented as sketched)."""
    mat = cfg["material"]
    nodes0, walls = _usec_walls(cfg)
    walls = _usec_sketch_direction(cfg, nodes0, walls)
    cen = _usec_centroid(nodes0, walls)
    lens = [float(np.hypot(*(nodes0[j] - nodes0[i]))) for i, j, _p in walls]
    if ds is None:                                         # ~ the Abaqus seed 0.025 h
        ds = MESH_STRIP if MESH_STRIP else min(2.5, 0.025 * max(lens))
    pts = [tuple(p) for p in nodes0]; strips = []; perim = 0.0
    for (i, j, plies), Lw in zip(walls, lens):
        A, B, D, _t = clt_plies(mat, plies)
        As = _osfe_laminate(mat, plies)[3]
        if theory == "CLPT":                               # Kirchhoff limit of the same model
            As = As * 1.0e3
        n = max(8, int(np.ceil(Lw / ds)))
        P0, P1 = nodes0[i], nodes0[j]; prev = i
        flip = _usec_normal_flip(P0, P1, cen)
        for q in range(1, n + 1):
            if q < n:
                pts.append(tuple(P0 + (P1 - P0) * q / n)); nxt = len(pts) - 1
            else:
                nxt = j
            # strip axis s = reverse of the sketch direction -> strip normal
            # x^ s = t x Z; the opposite when the rule says so (see
            # SECTION_NORMAL_RULE); ply angles measured from x towards s
            if flip:
                strips.append((prev, nxt, (A, B, D), As))
            else:
                strips.append((nxt, prev, (A, B, D), As))
            prev = nxt
        perim += Lw
    # along the member: cubic spline spans of min(5, L/40) in the interior,
    # halved within one section depth of each end (support / load-introduction
    # zones: varying pre-buckling stress, end-localised modes)
    if hx is None:
        if not MESH_SPAN:
            hx = min(5.0, L / 40.0)
        else:                          # > 0: span [mm] (<= L/20);  < 0: L / |value|
            hx = L / -MESH_SPAN if MESH_SPAN < 0 else min(MESH_SPAN, L / 20.0)
    breaks = _sfs_breaks(L, hx, 0.5 * hx, max(lens))
    return SplineStripModel(np.asarray(pts, float), strips, L, 0, breaks=breaks), perim


def _usec_supports(M, bc):
    """Support conditions of the Abaqus parametric scripts.
    SS : U1 = U2 = 0 on every end-edge node, U3 held at one node, edge load on
         both ends;  CC : rigid non-rotating end plates (kinematic coupling to a
         reference point: bottom fully fixed, top free only in U3)."""
    nb, N = M.nb, M.N
    fixed = set(); merge = []
    if bc == "SS":
        for n in range(N):
            for i in (0, nb - 1):
                fixed.add((n, 1, i)); fixed.add((n, 2, i))
        fixed |= M.end_bubble_keys(which=(1,))             # whole edge: in-plane v = 0
        fixed.add((0, 0, 0))
    elif bc == "CC":
        for n in range(N):
            for d in range(6):
                fixed.add((n, d, 0))
            for d in (1, 2, 3, 4, 5):
                fixed.add((n, d, nb - 1))
        fixed |= M.end_bubble_keys(which=(0, 1))           # rigid plates: flat, no v
        merge = [[(n, 0, nb - 1) for n in range(N)]]
    else:
        raise ValueError("bc must be 'SS' or 'CC'")
    return fixed, merge


def _usec_mode_info(M, vec):
    """Half-wave count (sign changes along x of the largest translation, as in
    the Abaqus post-processing) and the global / distortional / local character
    (share of the cross-section displacement that is rigid-body, at the section
    of maximum amplitude)."""
    from scipy.interpolate import BSpline
    t = np.r_[[0.0] * 3, M.breaks, [M.L] * 3]
    xs = np.linspace(0.0, M.L, 201)
    Bx = BSpline.design_matrix(xs, t, 3).toarray()          # (201, nb)
    U = M.node_part(vec) @ Bx.T                             # (N, 6, 201)
    tr = np.linalg.norm(U[:, 1:3, :], axis=1)               # |(UX, UY)|
    n0, x0 = np.unravel_index(np.argmax(tr), tr.shape)
    comp = 1 + int(np.argmax(np.abs(U[n0, 1:3, x0])))
    line = U[n0, comp, :]; line = line[np.abs(line) > 1e-6 * np.abs(line).max()]
    s = np.sign(line); m = int(np.sum(s[:-1] * s[1:] < 0)) + 1
    X = M.nodes[:, 0] - M.nodes[:, 0].mean(); Y = M.nodes[:, 1] - M.nodes[:, 1].mean()
    A = np.zeros((2 * M.N, 3)); A[0::2, 0] = 1.0; A[1::2, 1] = 1.0
    A[0::2, 2] = -Y; A[1::2, 2] = X
    b = np.empty(2 * M.N); b[0::2] = U[:, 1, x0]; b[1::2] = U[:, 2, x0]
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    r = float((A @ c) @ (A @ c) / max(b @ b, 1e-300))
    kind = "global" if r > 0.8 else ("local" if r < 0.3 else "distortional")
    return m, kind, r


def unified_member_buckling(cfg, L, bc="SS", theory="FSDT", nev=2, ds=None, hx=None,
                            return_mode=False):
    """Critical axial load P_cr [N] of the member (section of cfg, length L)
    with the unified shear-deformable spline finite strip, the Abaqus support
    conditions and loading, and the member's own pre-buckling state."""
    M, perim = _usec_build(cfg, L, theory, ds, hx)
    fixed, merge = _usec_supports(M, bc)
    f = M.end_line_load([1.0 / perim] * len(M.S))           # total end load = 1 N
    lam, out = M.buckle(fixed, f, merge=merge, nev=nev, return_modes=True)
    m, kind, r = _usec_mode_info(M, out["modes"][0])
    res = dict(P=float(lam[0]), eigs=[float(v) for v in lam], m=m, kind=kind,
               rigid_share=r, ndof=M.ndof, bc=bc, theory=theory,
               name=f"unified shear-deformable spline strip ({bc})")
    if return_mode:
        res.update(model=M, vec=out["modes"][0])
    return res


# ==============================================================================
# 21.  UNIFIED SOLVER  <->  REST OF THE CODE  (analyse, sweep, plots, studies)
# ------------------------------------------------------------------------------
# analyse() reports as GOVERNING load the spline-strip solution of section 20
# for every cross-section (c, l, z, lipc, t, i, box, hat, dhat) and both end
# conditions.  The classical quantities (EA, EI, GJ, Euler / torsional /
# per-wall plate loads, GBT deformation modes) are still computed and printed,
# for interpretation only.  No short-member override, no laminate-class switch,
# no acceptance threshold, no effective length, no grip on/off rule.
# ==============================================================================
def _unified_governing(cfg, L0):
    """Governing P_cr [N] at the design length and a short description."""
    bc = str(cfg.get("bc", "SS")).upper()
    solver, fam, meth, seg, hx = _gbt_opts(cfg)
    if solver == "gbt":
        r = gbt_member_buckling(cfg, float(L0), bc, return_mode=True)
        part = r["participation"]
        dom = max(part, key=part.get) if part else "?"
        info = dict(kind=f"shear-deformable GBT solved by {r['method']} finite strips "
                         f"(modes {fam}, {bc})",
                    participation=part, mode=_GBT_GROUP_NAMES.get(dom, dom),
                    eigs=r["eigs"], ndof=r["ndof"], m=r.get("m"))
        r["kind"] = info["mode"]
        return r["P"], info, r
    r = unified_member_buckling(cfg, float(L0), bc, return_mode=True)
    info = dict(kind="shell-type shear-deformable spline finite strip (FSDT)",
                m=r["m"], mode=r["kind"], rigid_share=r["rigid_share"],
                eigs=r["eigs"], ndof=r["ndof"])
    return r["P"], info, r


def _unified_sweep(cfg, Ls, nmax=10):
    """P_cr(L) curve (and GBT modal participation) at up to nmax lengths of the
    sweep grid, interpolated log-log for the curve.  GBT solver: the same GBT
    finite strips as the governing value;  strip solver: harmonic strips with
    the interior pre-buckling state (section 19)."""
    Ls = np.asarray(Ls, float); bc = str(cfg.get("bc", "SS")).upper()
    solver = _gbt_opts(cfg)[0]
    if len(Ls) <= nmax:
        Lsub = Ls.copy()
    else:
        idx = np.unique(np.round(np.geomspace(1, len(Ls), nmax)).astype(int) - 1)
        Lsub = Ls[idx]
    print(f"     {'L [mm]':>8} | {'P_cr [kN]':>12}   "
          + ("GBT participation of the mode" if solver == "gbt" else "half-waves m / mode (harmonic strip)"))
    print("     " + "-" * 64)
    Psub = []; parts = []
    sec = None
    if solver == "gbt":
        _s, fam, meth, seg, hx = _gbt_opts(cfg)
        sec = GBTSection(cfg, seg, fam, GBT_NMIN)
    for L in Lsub:
        if solver == "gbt":
            r = gbt_fsm_solve(sec, float(L), bc, method=("sine" if (_gbt_opts(cfg)[2] == "sine"
                              and bc == "SS") else "spline"), hx=min(_gbt_opts(cfg)[4], L / 20.0))
            parts.append(r["participation"])
            txt = "  ".join(f"{g} {100*v:3.0f}%" for g, v in r["participation"].items() if v > 0.005)
        else:
            r = unified_section_solve(cfg, float(L), bc)
            txt = f"m={r['m']}  {r['kind']}"
        Psub.append(r["P"])
        print(f"     {L:>8.0f} | {r['P']/1e3:>12.4f}   {txt}")
    Psub = np.asarray(Psub)
    cfg["_sweep_part"] = (Lsub, parts) if parts else None
    if len(Lsub) == 1:
        return np.full_like(Ls, Psub[0])
    return np.exp(np.interp(np.log(Ls), np.log(Lsub), np.log(Psub)))


def _plot_unified_mode(cfg, res, outdir, tag):
    """Critical buckling mode of the unified solver: cross-section at the
    position of maximum amplitude, and the amplitude along the member."""
    r = res.get("_unified")
    if not r or "model" not in r:
        return
    import matplotlib.pyplot as plt
    from scipy.interpolate import BSpline
    M, vec = r["model"], r["vec"]
    t = np.r_[[0.0] * 3, M.breaks, [M.L] * 3]
    xs = np.linspace(0.0, M.L, 401)
    U = M.node_part(vec) @ BSpline.design_matrix(xs, t, 3).toarray().T
    tr = np.hypot(U[:, 1, :], U[:, 2, :])
    n0, x0 = np.unravel_index(np.argmax(tr), tr.shape)
    XY = M.nodes; d = U[:, 1:3, x0]; d = d / max(np.abs(d).max(), 1e-30)
    amp = 0.2 * max(np.ptp(XY[:, 0]), np.ptp(XY[:, 1]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6),
                                 gridspec_kw=dict(width_ratios=[1, 1.6]))
    for st in M.strips:
        i, j = st[0], st[1]
        a1.plot(XY[[i, j], 0], XY[[i, j], 1], color="0.75", lw=1.5)
        a1.plot(XY[[i, j], 0] + amp * d[[i, j], 0], XY[[i, j], 1] + amp * d[[i, j], 1],
                color="#c0392b", lw=2.0)
    a1.set_aspect("equal"); a1.set_xticks([]); a1.set_yticks([])
    a1.set_title(f"cross-section at x = {xs[x0]:.0f} mm")
    comp = 1 + int(np.argmax(np.abs(U[n0, 1:3, x0])))
    a2.plot(xs, U[n0, comp, :] / np.abs(U[n0, comp, :]).max(), "k", lw=1.8)
    a2.set_xlabel("x [mm]"); a2.set_ylabel("normalised displacement")
    a2.set_title(f"along the member (node of max. amplitude), m = {r['m']}")
    a2.grid(True, alpha=0.3)
    fig.suptitle(f"Critical mode, unified spline strip ({cfg['bc']}):  "
                 f"P_cr = {r['P']/1e3:.4f} kN,  {r['kind']}", fontweight="bold")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig_critical_mode_{tag}.png", dpi=200)
    plt.close()


def _plot_gbt_participation(cfg, outdir, tag):
    """Modal participation diagram vs L (as the paper's Fig. 16b)."""
    sp_ = cfg.get("_sweep_part")
    if not sp_:
        return
    import matplotlib.pyplot as plt
    Lsub, parts = sp_
    groups = [g for g in ("G", "D", "L", "S", "T", "M") if any(p.get(g, 0) > 0 for p in parts)]
    cols = {"G": "#3498db", "D": "#e67e22", "L": "#9b59b6", "S": "#2ecc71", "T": "#95a5a6", "M": "#e74c3c"}
    base = np.zeros(len(Lsub))
    plt.figure(figsize=(8.5, 3.6))
    for g in groups:
        v = np.array([p.get(g, 0.0) for p in parts])
        plt.fill_between(Lsub, base, base + v, color=cols[g], alpha=0.75, step=None,
                         label=_GBT_GROUP_NAMES[g])
        base = base + v
    plt.xscale("log"); plt.ylim(0, 1); plt.xlim(Lsub.min(), Lsub.max())
    plt.xlabel("L [mm]"); plt.ylabel("GBT modal participation")
    plt.title(f"GBT modal participation of the critical mode ({tag}, {cfg['bc']})")
    plt.legend(fontsize=8, loc="lower right"); plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(f"{outdir}/fig_gbt_participation_{tag}.png", dpi=200); plt.close()


def _plot_gbt_mode(cfg, res, outdir, tag):
    """Critical GBT-FSM mode: deformed cross-section where the in-plane
    displacement is largest, and the amplitude functions of the dominant modes."""
    r = res.get("_unified")
    if not r or "coef" not in r:
        return
    import matplotlib.pyplot as plt
    sec, lb, C = r["sec"], r["lb"], r["coef"]
    nip = sec.nip
    if lb.kind == "spline":
        A = np.asarray(lb.B[0] @ C[:nip].T); xs = lb.xg
    else:
        A = lb.Bf["phi"][0] @ C[:nip].T; xs = lb.xg
    Ud = A @ sec.U.T                                          # (nx, 2N)
    mag = np.hypot(Ud[:, 0::2], Ud[:, 1::2]).max(axis=1); ix = int(np.argmax(mag))
    d = Ud[ix].reshape(-1, 2); d = d / max(np.hypot(d[:, 0], d[:, 1]).max(), 1e-30)
    amp = 0.2 * max(np.ptp(sec.X[:, 0]), np.ptp(sec.X[:, 1]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1, 1.6]))
    for s in sec.S:
        i, j = s["a"], s["b"]
        a1.plot(sec.X[[i, j], 0], sec.X[[i, j], 1], color="0.75", lw=1.5)
        a1.plot(sec.X[[i, j], 0] + amp * d[[i, j], 0], sec.X[[i, j], 1] + amp * d[[i, j], 1],
                color="#c0392b", lw=2.0)
    a1.set_aspect("equal"); a1.set_xticks([]); a1.set_yticks([])
    a1.set_title(f"cross-section at x = {xs[ix]:.0f} mm")
    I = lb.wx @ np.abs(A)
    top = np.argsort(-I)[:6]
    for k in top:
        a2.plot(xs, A[:, k] / np.abs(A).max(), lw=1.6, label=f"mode {k+1} ({sec.group[k]})")
    a2.set_xlabel("x [mm]"); a2.set_ylabel(r"amplitude $\phi_k(x)$ (normalised)")
    a2.legend(fontsize=8); a2.grid(True, alpha=0.3)
    part = r.get("participation", {})
    fig.suptitle(f"Critical mode, GBT-FSM ({cfg['bc']}): P_cr = {r['P']/1e3:.4f} kN   "
                 + "  ".join(f"{g} {100*v:.0f}%" for g, v in part.items() if v > 0.005),
                 fontweight="bold")
    plt.tight_layout(); plt.savefig(f"{outdir}/fig_critical_mode_{tag}.png", dpi=200); plt.close()


# ---------------------------------------------------------------- studies ----
_STUDY_MATS = dict(
    iso=dict(type="isotropic", E=200000.0, nu=0.3),
    cfrp=dict(type="orthotropic", E1=171400.0, E2=9080.0, G12=5290.0,
              nu12=0.32, G13=5290.0, G23=2645.0))
_STUDY_PLAN = [("Isotropic", "iso", None),
               ("Sym1", "cfrp", [45, -45, -45, 45]),
               ("Sym2", "cfrp", [60, -30, -30, 60]),
               ("Unsym", "cfrp", [0, 30, 60, 90])]
_STUDY_BCS = ["SS", "CC"]
_STUDY_LENGTHS = [200.0, 500.0, 1000.0]


def _study_cfg(shape, material, layup, H=50.0, B=20.0, TPLY=0.25):
    """Section of the Abaqus parametric scripts: height H, flange / width B,
    t = 1 (isotropic: one 1 mm layer; laminates: 4 x 0.25 mm plies)."""
    lay = layup if layup else [0]
    if shape == "l":                                     # L-section: legs B (horizontal) and H
        t = len(lay) * TPLY if layup else 1.0
        return dict(material=material, layup=lay, family="A", shape="l",
                    corner_nodes=[(B, 0.0), (0.0, 0.0), (0.0, H)],
                    t_wall=[t, t], b=H, a=B)
    if shape in ("c", "z"):
        tip = B if shape == "c" else -B
        corners = [(tip, -H / 2.0), (0.0, -H / 2.0), (0.0, H / 2.0), (B, H / 2.0)]
        t = len(lay) * TPLY if layup else 1.0
        return dict(material=material, layup=lay, corner_nodes=corners, family="A",
                    shape=shape, t_wall=[t, t, t], b=H, a=B)
    dims = dict(a=B, b=H)
    if layup is None:
        nodes, walls, cells, fam = make_section(shape, dims, material, [0], [0], 1.0, 0.0)
    else:
        nodes, walls, cells, fam = make_section(shape, dims, material, layup, layup,
                                                TPLY, TPLY)
    return dict(material=material, layup=lay, nodes=nodes, walls=walls, cells=cells,
                family=fam, shape=shape, layup_b=lay, layup_s=lay)


def _section_study(shape, title, csv_out):
    """Parametric study of one section, same plan as the Abaqus scripts:
    materials Isotropic / Sym1 / Sym2 / Unsym, BC SS / CC, L = 200/500/1000,
    with the governing solver of section 24."""
    solver, fam, meth, seg, hx = _gbt_opts({})
    if solver != "gbt":
        csv_out = csv_out.replace("_GBT_", "_strip_")
    line = "=" * 100
    print(line)
    print(f" PARAMETRIC STUDY  -  {title}")
    if solver == "gbt":
        print(f" shear-deformable GBT (Silvestre & Camotim 2013) solved by {meth} finite strips;"
              f" modes {fam}; strip {seg:g} mm, span {hx:g} mm")
        hdr = "GBT modal participation [%]  G / D / L / S" + (" / T" if "T" in fam else "") \
              + (" / M" if "M" in fam else "")
    else:
        print(" shell-type shear-deformable spline finite strip (FSDT), Abaqus end conditions")
        hdr = "m   mode    rigid%"
    print(line)
    print(f"{'material':>10}{'BC':>4}{'L[mm]':>7} | {'P_cr [N]':>12}{'time [s]':>9}   {hdr}")
    print("-" * 100)
    rows = []; t_all = time.time()
    for name, mk, layup in _STUDY_PLAN:
        cfg = _study_cfg(shape, _STUDY_MATS[mk], layup)
        sec = GBTSection(cfg, seg, fam, GBT_NMIN) if solver == "gbt" else None
        for bc in _STUDY_BCS:
            for L in _STUDY_LENGTHS:
                t0 = time.time()
                try:
                    if solver == "gbt":
                        r = gbt_fsm_solve(sec, float(L), bc,
                                          method=("sine" if (meth == "sine" and bc == "SS") else "spline"),
                                          hx=min(hx, L / 20.0))
                        pr = r["participation"]
                        txt = " / ".join(f"{100*pr.get(g, 0):3.0f}" for g in "GDLS" + ("T" if "T" in fam else "")
                                         + ("M" if "M" in fam else ""))
                        extra = ";".join(f"{g}={pr.get(g, 0):.3f}" for g in pr)
                    else:
                        r = unified_member_buckling(cfg, L, bc)
                        txt = f"{r['m']:<3} {r['kind']:<14}{100*r['rigid_share']:>4.0f}"
                        extra = f"m={r['m']};{r['kind']}"
                    dt = time.time() - t0
                    rows.append((name, bc, int(L), r["P"], extra, dt))
                    print(f"{name:>10}{bc:>4}{L:>7.0f} | {r['P']:>12.3f}{dt:>9.1f}   {txt}", flush=True)
                except Exception as e:
                    rows.append((name, bc, int(L), None, str(e), time.time() - t0))
                    print(f"{name:>10}{bc:>4}{L:>7.0f} | {'FAILED':>12}   ({e})", flush=True)
        print("-" * 100)
    print(f" total time {time.time() - t_all:.0f} s")
    print(" P_cr = first (critical) eigenvalue.  SS: every end-edge node U1=U2=0, U3 held at")
    print(" one node, uniform edge load on both ends.  CC: rigid end plates, one end fixed,")
    print(" the other free only in U3 and loaded.")
    if solver == "gbt":
        print(" G global, D distortional, L local, S shear (membrane), T transverse extension,")
        print(" M transverse shear: share of int |amplitude| dx of the critical mode (paper Fig. 16b).")
    try:
        with open(csv_out, "w") as f:
            f.write("material,BC,L_mm,Pcr_N,mode_info,time_s\n")
            for (nm, bc, L, P, extra, dt) in rows:
                f.write(f"{nm},{bc},{L},{'' if P is None else round(P, 3)},{extra},{dt:.1f}\n")
        print(f" results written to {csv_out}")
    except Exception:
        pass
    return rows


# ==============================================================================
# 22.  OMEGA PANEL  -  SPLINE FINITE STRIP  (the strip reference / oracle; GBT: section 25)
# ------------------------------------------------------------------------------
# The strips, laminates and rigid tie links of section 14e (Mindlin strips on
# the Abaqus geometry: web roots at y = clear linked to the skin/foot below),
# with the cubic-spline longitudinal basis of section 20.  So the Abaqus end
# conditions are imposed node by node and the pre-buckling state -- including
# the load-introduction zones at the ends and the side restraint -- is the
# panel's own static solution:
#   SS : ends U2 = 0 on every end node; sides U1 = UR2 = UR3 = 0; one U3
#        anchor; 1 N/mm edge load on every end edge (skin and stringer walls;
#        a bonded foot carries 2 N/mm, applied at the height of its resultant)
#   CC : clamped ends: U1 = U2 = UR1 = UR2 = UR3 = 0 at both ends, U3 held at
#        z = 0; the same 1 N/mm edge loads (OMEGA_CC_END = "edge-load"), or
#        U3 common to a rigid loaded plate at z = L ("rigid-plate")
# No harmonic window, no uniform-prestress assumption.
# ==============================================================================
# Mesh of the omega spline strip, scaled with the stringer pitch bs (the size
# of the skin bays that buckle):  strips ds = bs/50, spans hx = bs/20
# (bs = 150: 3.0 / 7.5 mm;  bs = 285.5: 5.7 / 14.3 mm).  Mesh study on the
# 20 extreme Abaqus cases (narrow-short and wide-long, 5 lay-ups, SS + CC):
# within 0.55 % (mean 0.40 %) of the fine mesh ds = 2, hx = min(5, L/40), and
# about 5x faster.  Choices (interactive menu, or flags):
#   "standard" : ds = bs/50, hx = bs/20                       (default)
#   "fast"     : ds = bs/40, hx = bs/15   within 1 % of fine, ~8x faster than fine
#   "fine"     : ds = 2 mm,  hx = min(5, L/40)                --omega-mesh-fine
OMEGA_MESH = "standard"
OMEGA_MESH_FINE = False                                   # (old switch, = "fine")


def _omega_mesh(panel, L, ds, hx):
    bs = float(panel["info"]["bs"])
    kind = "fine" if OMEGA_MESH_FINE else OMEGA_MESH
    if ds is None:
        ds = (2.0 if kind == "fine" else
              min(6.0, max(1.0, bs / (40.0 if kind == "fast" else 50.0))))
    if hx is None:
        hx = (min(5.0, L / 40.0) if kind == "fine" else
              min(15.0, max(2.0, bs / (15.0 if kind == "fast" else 20.0)), L / 10.0))
    return ds, hx


def omega_spline_model(panel, L, ds=None, hx=None, bubble=False, refine=False):
    """Spline-strip model of the omega panel.  The membrane bubbles and the
    end-zone span refinement used for the member sections change the omega
    eigenvalues by < 0.05 % (skin-dominated plate modes; measured on
    Iso SS, Sym2 CC, Unsym CC) but double the size of these wide panels, so
    they are off by default here."""
    ds, hx = _omega_mesh(panel, L, ds, hx)
    nodes, strips, links = _omega_mfs_geometry(panel, ds)
    S = [(iA, iB, (w["ABD"][0], w["ABD"][1], w["ABD"][2]), As)
         for (iA, iB, w, As) in strips]
    nd0 = panel["nodes"]                                  # end zones: widest wall
    zone = max(float(np.hypot(*(nd0[w["j"]] - nd0[w["i"]]))) for w in panel["walls"])
    breaks = _sfs_breaks(L, hx, 0.5 * hx if refine else hx, zone)
    M = SplineStripModel(nodes, S, L, 0, links=links, breaks=breaks, bubble=bubble)
    return M, strips


def _omega_spline_load(M, strips):
    """Consistent end loads: 1 N/mm per loaded edge, at the height of the
    resultant (a foot = skin + stringer flange: 2 N/mm at clear/2 above the
    skin mid-plane -> axial force + moment about the transverse axis)."""
    f = np.zeros(M.ndof); nb = M.nb
    for k, (iA, iB, w, As) in enumerate(strips):
        b = M.S[k]["b"]; q = float(w["edges"])
        e = (w["yload"] - 0.5 * (M.nodes[iA][1] + M.nodes[iB][1])) if w["part"] == "foot" else 0.0
        for n in (iA, iB):
            F = q * b / 2.0
            f[(6 * n + 0) * nb + 0] += F;        f[(6 * n + 0) * nb + nb - 1] -= F
            if e:
                f[(6 * n + 4) * nb + 0] += F * e; f[(6 * n + 4) * nb + nb - 1] -= F * e
        if M.bubble:                                     # u-bubble: q b 2/3
            fb = (6 * M.N + 2 * k) * nb
            f[fb + 0] += q * b * 2.0 / 3.0; f[fb + nb - 1] -= q * b * 2.0 / 3.0
    return f


def _omega_spline_supports(M, panel, links, bc):
    nb = M.nb; slave = {s for s, _m in links}
    ind = [n for n in range(M.N) if n not in slave]
    fixed = set(); merge = []
    sL, sR = panel["sides"]                               # x = -W, x = 0
    if OMEGA_SIDE_BC == "periodic":                       # u(-W, z) = u(0, z)
        rigid_u3 = bc == "CC" and OMEGA_CC_END == "rigid-plate"
        for d in range(6):
            for i in range(nb):
                if rigid_u3 and d == 0 and i == nb - 1:
                    continue                              # already in the plate group
                merge.append([(sR, d, i), (sL, d, i)])
        if bc != "CC":
            fixed.add((sR, 1, 0))                         # U1 anchor (rigid translation)
            if OMEGA_PERIODIC_MACRO:                      # free macro strains of the cell:
                # u1(-W) = u1(0) - W eps_xx ,  u3(-W) = u3(0) - W gam_xz  (constant
                # jumps; with CC ends both are zero anyway, U1 = U3 = 0 there)
                M.jumps = [[(6 * sL + 1, i, -1.0) for i in range(nb)],
                           [(6 * sL + 0, i, -1.0) for i in range(nb)]]
                fixed.add((sR, 1, nb - 1))                # + in-plane rigid rotation
    else:
        for s in panel["sides"]:                          # symmetry: U1 = UR2 = UR3 = 0
            for i in range(nb):
                for d in (1, 5, 3):
                    fixed.add((s, d, i))
    if bc == "CC":
        for n in ind:
            for i in (0, nb - 1):
                for d in (1, 2, 3, 4, 5):
                    fixed.add((n, d, i))
            fixed.add((n, 0, 0))
        if OMEGA_CC_END == "rigid-plate":                 # U3 common at z = L
            fixed |= M.end_bubble_keys(which=(0, 1))
            merge = merge + [[(n, 0, nb - 1) for n in ind]]
        else:                                             # clamped edges, U3 free at z = L
            if M.bubble:
                fixed |= {("B", k, j, 0) for k in range(len(M.S)) for j in (0, 1)}
                fixed |= {("B", k, 1, nb - 1) for k in range(len(M.S))}
    else:
        for n in ind:
            fixed.add((n, 2, 0)); fixed.add((n, 2, nb - 1))      # U2 = 0
        # U2 = v sin(al) + w cos(al) on the whole end edge -> v-bubble of
        # every non-horizontal strip is zero there
        fixed |= M.end_bubble_keys(which=(1,), strips=[k for k, st in enumerate(M.S)
                                                       if abs(np.sin(st["al"])) > 1e-9])
        fixed.add((panel["sides"][1], 0, 0))                     # U3 anchor (x = 0)
    return fixed, merge


def _omega_spline_modeinfo(M, vec):
    """Half-waves as in the Abaqus script (skin U2 along the line through the
    point of max |U2|) and the stringer share of the mode."""
    xs = np.linspace(0.0, M.L, 301)
    U = _omega_U(M, vec, xs)
    y = M.nodes[:, 1]; sk = np.where(np.abs(y) < 1e-9)[0]; st = np.where(y > 1e-9)[0]
    Usk = U[sk, 2, :]
    n0, _x0 = np.unravel_index(np.argmax(np.abs(Usk)), Usk.shape)
    line = Usk[n0]; line = line[np.abs(line) > 1e-6 * np.abs(line).max()]
    s = np.sign(line); m = int(np.sum(s[:-1] * s[1:] < 0)) + 1
    a_sk = np.hypot(U[sk, 1, :], U[sk, 2, :]).max()
    a_st = np.hypot(U[st, 1, :], U[st, 2, :]).max() if len(st) else 0.0
    r = a_st / max(a_sk, a_st, 1e-30)
    kind = "local skin" if r < 0.25 else ("skin+stringer (distortional)" if r < 0.8
                                          else "stringer-dominated")
    return m, r, kind


def omega_spline_eigs(panel, L, bc="SS", nmodes=3, ds=None, hx=None, bubble=False,
                      refine=False, return_model=False):
    """Eigenvalues [(N_cr_i [N/mm], m_i)] of the omega panel (spline strip)."""
    ds, hx = _omega_mesh(panel, L, ds, hx)
    M, strips = omega_spline_model(panel, L, ds, hx, bubble, refine)
    fixed, merge = _omega_spline_supports(M, panel, M.links, bc)
    f = _omega_spline_load(M, strips)
    lam, out = M.buckle(fixed, f, merge=merge, nev=max(2, nmodes), return_modes=True,
                        tol=1e-7)
    modes = []
    for l, v in list(zip(lam, out["modes"]))[:max(1, nmodes)]:
        m, _r, _k = _omega_spline_modeinfo(M, v)
        modes.append((float(l), m))
    m, r, kind = _omega_spline_modeinfo(M, out["modes"][0])
    info = dict(ndof=M.ndof, stringer_share=r, mode_type=kind, ds=ds, hx=hx)
    if return_model:
        info.update(model=M, strips=strips, vecs=out["modes"][:max(1, nmodes)],
                    panel=panel)
    return modes, info


# ==============================================================================
# 23.  SHEAR-DEFORMABLE GBT SOLVED BY THE FINITE STRIP METHOD  (GBT-FSM)
# ------------------------------------------------------------------------------
# Silvestre & Camotim, "Shear deformable generalized beam theory for the
# analysis of thin-walled composite members", J. Eng. Mech. 139(8):1010-1024
# (2013) -- with the member (longitudinal) problem solved by FINITE STRIPS
# instead of the paper's beam finite elements.
#
# 1. CROSS-SECTION ANALYSIS (paper Eq. 1, Figs. 1-3).  Natural + intermediate
#    nodes on the wall mid-lines; any topology (open, branched, closed).
#    Mode families (each with its own amplitude function along x):
#      C  conventional modes (paper): in-plane nodal displacement fields that
#         keep every wall segment transversely inextensible (eps_ss^M = 0);
#         nodal rotations from the transverse frame analysis (D22, rigid
#         joints); w(s) cubic (Hermite) in each segment; warping u_k(s)
#         linear in each segment from Vlasov's condition u_k,s + v_k = 0
#         (closed cells: least-squares / Bredt warping, min int A66 g^2 ds).
#         Split into G (3 rigid in-plane motions), D (distortional) and
#         L (local: every fold / branch node at rest).  Plus the axial mode.
#      S  shear modes (paper): nodal warping profiles u_j(s), v = w = 0,
#         independent amplitude psi_j(x)  -> membrane shear strains.
#    Optional extensions of the same GBT (NOT in the 2013 paper; off unless
#    requested):
#      T  transverse-extension modes (Silvestre, Camotim & Silva 2011): the
#         in-plane nodal fields that stretch the walls (eps_ss^M != 0);
#      M  transverse-shear (Mindlin) modes: independent wall rotations,
#         gamma_xz, gamma_sz != 0 (laminate As, shear factor 5/6).
#    Displacements:  u = u_k phi_k,x + u_axial(x) + u_j psi_j ;
#                    v = v_k phi_k ;  w = w_k phi_k   (k over C and T modes)
# 2. STRAINS / CONSTITUTIVE LAW (paper Eqs. 2-5, 10-12): Love-Kirchhoff walls
#    (+ rotations of the M modes if used); full laminate [A B; B D];
#    without T modes the N_ss = 0 condensation of Silvestre & Camotim (2002a)
#    (the paper's treatment), with T modes the full 6x6 law.  The energy is
#    integrated numerically across the section, which reproduces the tensors
#    C, D, B, H, F of Eqs. (10)-(12) term by term.
# 3. GEOMETRIC STIFFNESS: paper Eq. (3), eps_NL = (v,x^2 + w,x^2)/2 acted on
#    by N_xx^0; with T modes all membrane resultants and gradients are used.
# 4. MEMBER ANALYSIS BY FINITE STRIPS
#    method="spline": amplitude functions in cubic B-splines (spline finite
#       strip).  End conditions of the Abaqus models:
#         SS : phi = 0 at both ends (every in-plane displacement held),
#              warping and rotations free, axial anchor at x = 0;
#         CC : phi = phi,x = 0, psi = 0, rotations = 0 at both ends (rigid
#              plates), axial held at x = 0 and free (loaded) at x = L.
#       Pre-buckling state = first-order GBT analysis under the end loads
#       (as in the paper's first example): N^0 varies along the member.
#    method="sine": classical semi-analytical GBT-FSM for SS ends,
#       phi = sin(m pi x/L), psi = cos(m pi x/L), coupled harmonics
#       (numerical x-integration), uniform pre-stress of paper Eq. (8):
#       N_xx^0 = P A11(s) / int A11 ds.
# ==============================================================================


class GBTSection:
    """GBT cross-section analysis of cfg (as for analyse()); seg = target
    segment width across the section [mm]; families: subset of "CSTM"."""

    def __init__(self, cfg, seg=5.0, families="CS", nmin=1):
        import scipy.linalg as sla
        from collections import Counter
        mat = cfg["material"]; self.families = families
        nodes0, walls = _usec_walls(cfg)
        walls = _usec_sketch_direction(cfg, nodes0, walls)
        cen = _usec_centroid(nodes0, walls)
        pts = [np.asarray(p, float) for p in nodes0]
        segs = []; self.C6 = []; self.As = []
        for iw, (i, j, plies) in enumerate(walls):
            A, B, D, As = _osfe_laminate(mat, plies)
            C6 = np.zeros((6, 6)); C6[:3, :3] = A; C6[:3, 3:] = B; C6[3:, :3] = B; C6[3:, 3:] = D
            self.C6.append(C6); self.As.append(As)
            P0, P1 = pts[i], pts[j]; Lw = float(np.hypot(*(P1 - P0)))
            flip = _usec_normal_flip(P0, P1, cen)
            n = max(nmin, int(np.ceil(Lw / seg - 1e-9)))       # segments in this wall
            chain = [i]
            for q in range(1, n):
                pts.append(P0 + (P1 - P0) * q / n); chain.append(len(pts) - 1)
            chain.append(j)
            for a, b in zip(chain[:-1], chain[1:]):
                # s against the sketch direction -> n = t_sketch x Z, or the
                # opposite when SECTION_NORMAL_RULE says so (default: outward)
                segs.append((a, b, iw) if flip else (b, a, iw))
        self.X = np.array(pts); self.N = N = len(pts); self.nw = len(walls)
        val = Counter()
        for i, j, _p in walls:
            val[i] += 1; val[j] += 1
        self.junctions = sorted(n for n in val if val[n] >= 2)
        self.tips = sorted(n for n in val if val[n] == 1)
        self.S = []
        for (a, b, iw) in segs:
            d = self.X[b] - self.X[a]; L = float(np.hypot(*d)); t = d / L
            self.S.append(dict(a=a, b=b, w=iw, L=L, t=t, n=np.array([-t[1], t[0]])))
        self.perimeter = sum(s["L"] for s in self.S)
        # constitutive law used by the energy (8 x 8: exx, ess, gxs, kxx, kss, kxs, gxz, gsz)
        self.C8 = []
        for iw in range(self.nw):
            C = np.zeros((8, 8))
            if "T" in families:
                C[:6, :6] = self.C6[iw]
            else:                                             # N_ss = 0 (paper)
                C5 = condense_Nss(self.C6[iw]); k = [0, 2, 3, 4, 5]
                C[np.ix_(k, k)] = C5
            C[6:, 6:] = self.As[iw]
            self.C8.append(C)
        # ---------------- in-plane modes (C, and T if requested) -------------
        E = np.zeros((len(self.S), 2 * N))
        for r, s in enumerate(self.S):
            E[r, 2 * s["b"]:2 * s["b"] + 2] += s["t"]; E[r, 2 * s["a"]:2 * s["a"] + 2] -= s["t"]
        Z = sla.null_space(E)
        c = self.X.mean(axis=0)
        G = np.zeros((2 * N, 3)); G[0::2, 0] = 1.0; G[1::2, 1] = 1.0
        G[0::2, 2] = -(self.X[:, 1] - c[1]); G[1::2, 2] = self.X[:, 0] - c[0]
        G = sla.orth(G)
        # local modes: unit transverse (normal) displacement of one node that is
        # not a fold / branch node (intermediate nodes, free tips) -- sparse
        jset = set(self.junctions); Lm = []
        for n in range(N):
            if n in jset:
                continue
            sg = next(s_ for s_ in self.S if n in (s_["a"], s_["b"]))
            e = np.zeros(2 * N); e[2 * n:2 * n + 2] = sg["n"]; Lm.append(e)
        Lm = np.array(Lm).T if Lm else np.zeros((2 * N, 0))
        # global: rigid in-plane motions not already spanned by the local modes
        # (e.g. for a T the rotation about the junction is itself a "local" field)
        Gk = []
        for k in range(3):
            B_ = np.hstack([Lm] + Gk) if (Lm.shape[1] or Gk) else np.zeros((2 * N, 0))
            r = G[:, k] - (B_ @ np.linalg.lstsq(B_, G[:, k], rcond=None)[0] if B_.shape[1] else 0)
            if np.linalg.norm(r) > 1e-8 * np.linalg.norm(G[:, k]):
                Gk.append(G[:, [k]])
        self._Gorth = G.copy()                             # all 3 rigid motions
        self._Torth = sla.orth(E.T)                         # stretching space
        G = np.hstack(Gk)
        GL = np.hstack([G, Lm])
        Dm = Z @ sla.null_space(GL.T @ Z) if GL.shape[1] < Z.shape[1] else Z[:, :0]
        Dm = sla.orth(Dm) if Dm.shape[1] else Dm
        self._Dorth = Dm
        blocks = [("G", G), ("D", Dm), ("L", Lm)]
        if "T" in families:                               # wall-segment stretching
            blocks.append(("T", E.T.copy()))
        U = np.hstack([b for _g, b in blocks])
        grp = sum([[g] * b.shape[1] for g, b in blocks], [])
        for k in range(U.shape[1]):
            u = U[:, k].reshape(N, 2); U[:, k] /= np.hypot(u[:, 0], u[:, 1]).max()
        self.R = self._frame_R()
        Th = self.R @ U                                   # frame rotations of each mode
        if "R" in families:
            # nodal-rotation modes: unit rotation of one node, all nodes at rest
            # (w Hermite inside the adjacent segments).  They complete the
            # transverse interpolation that the frame analysis ties to the
            # nodal displacements: same converged answer, far fewer nodes.
            # with every nodal rotation free, the frame-analysis rotations of the
            # other modes add nothing to the kinematic space -> drop them (the
            # modes become local in s and the stiffness matrix stays sparse)
            U = np.hstack([U, np.zeros((2 * N, N))])
            Th = np.hstack([np.zeros_like(Th), np.eye(N)])
            grp = grp + ["L"] * N
        self.U = U; self.Th = Th; self.nip = U.shape[1]
        # warping of each in-plane mode: min sum A66 L (du/L + v)^2, int A11 u = 0
        Aw = np.zeros((len(self.S), N)); wts = np.zeros(len(self.S)); c_ax = np.zeros(N)
        for r, s in enumerate(self.S):
            Aw[r, s["b"]] += 1.0 / s["L"]; Aw[r, s["a"]] -= 1.0 / s["L"]
            C5 = condense_Nss(self.C6[s["w"]])
            wts[r] = C5[1, 1] * s["L"]
            c_ax[[s["a"], s["b"]]] += 0.5 * C5[0, 0] * s["L"]
        W = np.sqrt(wts)[:, None]; AW = Aw * W
        Msys = np.block([[AW.T @ AW, c_ax[:, None]], [c_ax[None, :], np.zeros((1, 1))]])
        self.Uw = np.zeros((N, self.nip))
        for k in range(self.nip):
            if grp[k] not in ("G", "D"):
                continue                                  # local / T / R: no warping
            Ud = U[:, k].reshape(N, 2)
            v = np.array([s["t"] @ (0.5 * (Ud[s["a"]] + Ud[s["b"]])) for s in self.S])
            self.Uw[:, k] = np.linalg.lstsq(Msys, np.r_[-AW.T @ (W[:, 0] * v), 0.0],
                                            rcond=None)[0][:N]
        self.c_ax = c_ax
        # ---------------- shear modes (paper) ---------------------------------
        if "S" in families:
            Bs = np.eye(N)[:, 1:]                         # unit nodal warping (hats)
        else:
            Bs = np.zeros((N, 0))
        self.Us = Bs; self.ns = Bs.shape[1]
        # ---------------- transverse-shear (Mindlin) modes --------------------
        if "M" in families:                                   # in-plane rotation vector,
            cols = []                                         # only the components seen
            for n in range(N):                                # by a wall (t . omega)
                ts = np.array([s_["t"] for s_ in self.S if n in (s_["a"], s_["b"])])
                for vec in sla.orth(ts.T).T:
                    e = np.zeros(2 * N); e[2 * n:2 * n + 2] = vec; cols.append(e)
            self.Om = np.array(cols).T
        else:
            self.Om = np.zeros((2 * N, 0))
        self.nmx = self.Om.shape[1]
        self.nms = N if "M" in families else 0                # rotation about x
        # mode bookkeeping:  [in-plane | axial | shear | M-x | M-s]
        self.o_ax = self.nip; self.o_s = self.nip + 1; self.o_mx = self.o_s + self.ns
        self.o_ms = self.o_mx + self.nmx; self.nm = self.o_ms + self.nms
        self.group = grp + ["A"] + ["S"] * self.ns + ["M"] * (self.nmx + self.nms)
        # amplitude family: "phi" (sin-like: in-plane, M-s) or "psi" (cos-like)
        self.fam = (["phi"] * self.nip + ["psi"] * (1 + self.ns + self.nmx)
                    + ["phi"] * self.nms)
        self._sample()

    def _frame_R(self):
        """Nodal rotations of the transverse frame (D22, rigid joints) for given
        in-plane nodal displacements:  theta = R U."""
        N = self.N; Kth = np.zeros((N, N)); KtU = np.zeros((N, 2 * N))
        for s in self.S:
            k = _beam_k(self.C6[s["w"]][4, 4], s["L"])           # D22
            ends = (s["a"], s["b"])
            for r, na in zip((1, 3), ends):
                for c, nb_ in zip((1, 3), ends):
                    Kth[na, nb_] += k[r, c]
                for c, nb_ in zip((0, 2), ends):
                    KtU[na, 2 * nb_:2 * nb_ + 2] += k[r, c] * s["n"]
        return -np.linalg.solve(Kth, KtU)

    def _sample(self):
        """Strain operators E (np, 8, nm, 3 x-orders) and displacement-gradient
        operators Gd (np, 3 comps u v w, 2 dirs x s, nm, 3 x-orders)."""
        gp, gw = np.polynomial.legendre.leggauss(4); xi = 0.5 * (gp + 1.0)
        nm, nip, N = self.nm, self.nip, self.N
        th = self.Th
        npnt = 4 * len(self.S)
        E = np.zeros((npnt, 8, nm, 3)); Gd = np.zeros((npnt, 3, 2, nm, 3))
        wp = np.zeros(npnt); C8 = np.zeros((npnt, 8, 8)); p = 0
        VW = np.zeros((npnt, 2, nip))
        for s in self.S:
            a, b, L = s["a"], s["b"], s["L"]
            Ua = self.U[2 * a:2 * a + 2]; Ub = self.U[2 * b:2 * b + 2]
            va, vb = s["t"] @ Ua, s["t"] @ Ub; wa, wb = s["n"] @ Ua, s["n"] @ Ub
            nod = np.vstack([wa, th[a], wb, th[b]])
            for q in range(4):
                H, dH, ddH = _hermite(xi[q], L)
                w, ws, wss = H @ nod, dH @ nod, ddH @ nod
                v = (1 - xi[q]) * va + xi[q] * vb; vs = (vb - va) / L
                uk = (1 - xi[q]) * self.Uw[a] + xi[q] * self.Uw[b]; uks = (self.Uw[b] - self.Uw[a]) / L
                e = E[p]; g = Gd[p]; I = slice(0, nip)
                VW[p, 0] = v; VW[p, 1] = w
                # in-plane modes: u = uk phi', v = v phi, w = w phi
                e[0, I, 2] = uk; e[1, I, 0] = vs; e[2, I, 1] = uks + v
                e[3, I, 2] = -w; e[4, I, 0] = -wss; e[5, I, 1] = -2.0 * ws
                g[0, 0, I, 2] = uk; g[0, 1, I, 1] = uks
                g[1, 0, I, 1] = v;  g[1, 1, I, 0] = vs
                g[2, 0, I, 1] = w;  g[2, 1, I, 0] = ws
                # axial: u = U1(x)
                e[0, self.o_ax, 1] = 1.0; g[0, 0, self.o_ax, 1] = 1.0
                # shear modes: u = uj psi
                if self.ns:
                    J = slice(self.o_s, self.o_s + self.ns)
                    uj = (1 - xi[q]) * self.Us[a] + xi[q] * self.Us[b]; ujs = (self.Us[b] - self.Us[a]) / L
                    e[0, J, 1] = uj; e[2, J, 0] = ujs
                    g[0, 0, J, 1] = uj; g[0, 1, J, 0] = ujs
                # Mindlin modes: psi_x = t . omega(s) chi(x);  psi_s = r(s) eta(x)
                if self.nmx:
                    ha, hb = 1 - xi[q], xi[q]
                    ta = s["t"] @ self.Om[2 * a:2 * a + 2]; tb = s["t"] @ self.Om[2 * b:2 * b + 2]
                    cx = ha * ta + hb * tb; cxs = (tb - ta) / L
                    MX = slice(self.o_mx, self.o_mx + self.nmx)
                    e[3, MX, 1] = cx; e[5, MX, 0] = cxs; e[6, MX, 0] = cx
                    cs_ = np.zeros(N); cs_[a] = ha; cs_[b] = hb
                    css = np.zeros(N); css[a] = -1.0 / L; css[b] = 1.0 / L
                    MS = slice(self.o_ms, self.o_ms + self.nms)
                    e[4, MS, 0] = css; e[5, MS, 1] = cs_; e[7, MS, 0] = cs_
                wp[p] = 0.5 * L * gw[q]; C8[p] = self.C8[s["w"]]; p += 1
        self.E, self.Gd, self.wp, self.Cp, self.VW = E, Gd, wp, C8, VW
        Pvw = np.zeros((npnt, 2, 2 * N)); p = 0; Rf = self.R
        for s in self.S:
            a, b, L = s["a"], s["b"], s["L"]
            for q in range(4):
                H, _dH, _ddH = _hermite(xi[q], L)
                Pvw[p, 0, 2 * a:2 * a + 2] += (1 - xi[q]) * s["t"]; Pvw[p, 0, 2 * b:2 * b + 2] += xi[q] * s["t"]
                Pvw[p, 1, 2 * a:2 * a + 2] += H[0] * s["n"]; Pvw[p, 1, 2 * b:2 * b + 2] += H[2] * s["n"]
                Pvw[p, 1, :] += H[1] * Rf[a] + H[3] * Rf[b]
                p += 1
        self.Pvw = Pvw
        self.A11p = np.array([condense_Nss(self.C6[s["w"]])[0, 0] for s in self.S for _q in range(4)])
        self.ukp = E[:, 0, :nip, 2]; self.ujp = E[:, 0, self.o_s:self.o_s + self.ns, 1]


class _LongSpline:
    """Cubic B-spline amplitude functions (spline finite strip)."""
    kind = "spline"

    def __init__(self, L, breaks):
        from scipy.interpolate import BSpline
        self.L = L; self.breaks = np.asarray(breaks, float); self.nspan = len(breaks) - 1
        nb, xs, wx, _Bv, _Bd = _sfs_basis(L, self.nspan, self.breaks)
        self.t = np.r_[[0.0] * 3, self.breaks, [L] * 3]
        self.nb = nb; self.xg = xs.ravel(); self.wx = wx.ravel()
        spl = [BSpline(self.t, np.eye(nb)[i], 3) for i in range(nb)]
        self.B = [sp.csr_matrix(BSpline.design_matrix(self.xg, self.t, 3))]
        for d in (1, 2):
            self.B.append(sp.csr_matrix(np.column_stack([f.derivative(d)(self.xg) for f in spl])))
        self.end = {d: [np.array([(f.derivative(d) if d else f)(x) for f in spl]) for x in (0.0, L)]
                    for d in (0, 1)}


class _LongSine:
    """Harmonic amplitude functions: phi = sin(m pi x / L), psi = cos(m pi x / L)."""
    kind = "sine"

    def __init__(self, L, harmonics, ngp=8):
        self.L = L; self.m = np.asarray(harmonics, int); self.nb = len(self.m)
        nsp = 2 * max(8, int(self.m.max()))
        gp, gw = np.polynomial.legendre.leggauss(ngp); h = L / nsp
        self.xg = (np.arange(nsp)[:, None] * h + 0.5 * h * (gp[None, :] + 1)).ravel()
        self.wx = np.tile(0.5 * h * gw, nsp)
        k = self.m * np.pi / L; s, c = np.sin(k * self.xg[:, None]), np.cos(k * self.xg[:, None])
        self.Bf = {"phi": [s, k * c, -k * k * s], "psi": [c, -k * s, -k * k * c]}


def _gbt_Kcs(sec):
    out = {}
    for d in range(3):
        for e in range(3):
            Ed, Ee = sec.E[:, :, :, d], sec.E[:, :, :, e]
            if not (np.any(Ed) and np.any(Ee)):
                continue
            Y = np.matmul(sec.Cp, Ee) * sec.wp[:, None, None]             # (np, 8, nm)
            A = Ed.reshape(-1, sec.nm).T @ Y.reshape(-1, sec.nm)
            if np.any(A):
                out[(d, e)] = A
    return out


def _gbt_kg_tensors(sec, Pxx, Pss=None, Pxs=None, full=False):
    """Cross-section geometric tensors T[(d,e)] (nx, nm, nm) for the amplitude
    derivative orders (d, e).  P** (nx, np) compression-positive resultants.
    Paper Kg: Pxx (v,x v,x + w,x w,x).  full: all membrane terms of u, v, w."""
    comps = (0, 1, 2) if full else (1, 2)
    terms = []                                     # (P, gradient dir a, dir b)
    terms.append((Pxx, 0, 0))
    if full:
        terms += [(Pss, 1, 1), (Pxs, 0, 1), (Pxs, 1, 0)]
    T = {}
    nm = sec.nm
    for d in range(3):
        for e in range(3):
            acc = None
            for (P, da, db) in terms:
                for c in comps:
                    Ga = sec.Gd[:, c, da, :, d]; Gb = sec.Gd[:, c, db, :, e]
                    ia = np.nonzero(np.abs(Ga).sum(axis=0))[0]; ib = np.nonzero(np.abs(Gb).sum(axis=0))[0]
                    if not len(ia) or not len(ib):
                        continue
                    W = P * sec.wp[None, :]                                   # (nx, np)
                    t = np.matmul((W[:, :, None] * Ga[None, :, ia]).transpose(0, 2, 1), Gb[:, ib])
                    if acc is None:
                        acc = np.zeros((P.shape[0], nm, nm))
                    acc[:, ia[:, None], ib[None, :]] += t
            if acc is not None:
                T[(d, e)] = acc
    return T


def _kg_flush(Kg, rows, cols, vals, ndof):
    """Add the collected span blocks to the sparse geometric stiffness and
    empty the lists (duplicates of overlapping spans are summed)."""
    if vals:
        T = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(ndof, ndof)).tocsr()
        Kg = T if Kg is None else Kg + T
        rows.clear(); cols.clear(); vals.clear()
    return Kg


def gbt_fsm_solve(cfg, L, bc="SS", method="spline", seg=5.0, families="CS", hx=None,
                  harmonics=None, nev=2, return_mode=False, prestress="first-order"):
    """Critical load P_cr [N] of the member: shear-deformable GBT (families
    "CS" = Silvestre & Camotim 2013; add "T", "M" for the extensions) solved
    by finite strips (method "spline": SS / CC, Abaqus ends, first-order
    pre-stress;  "sine": SS, classical harmonic strips, uniform pre-stress)."""
    sec = cfg if isinstance(cfg, GBTSection) else GBTSection(cfg, seg, families)
    bc = bc.upper()
    if method == "sine":
        if bc != "SS":
            raise ValueError("sine strips: SS ends only")
        return _gbt_sine(sec, float(L), harmonics, nev, return_mode)
    L = float(L)
    if hx is None:
        hx = min(5.0, L / 40.0)
    zone = max(np.ptp(sec.X[:, 0]), np.ptp(sec.X[:, 1]))
    lb = _LongSpline(L, _sfs_breaks(L, hx, 0.5 * hx, zone))
    nb, nm, nip = lb.nb, sec.nm, sec.nip; ndof = nm * nb
    K = None
    for (d, e), A in _gbt_Kcs(sec).items():
        X = (lb.B[d].T @ sp.diags(lb.wx) @ lb.B[e]).tocsr()
        T = sp.kron(sp.csr_matrix(A), X, format="csr")
        K = T if K is None else K + T
    # ---- end conditions ------------------------------------------------------
    fixed = set()
    def ends(m, k=1):
        for i in range(k):
            fixed.update({m * nb + i, m * nb + nb - 1 - i})
    for m in range(nip):
        ends(m, 2 if bc == "CC" else 1)                     # in-plane held (+ no warping)
    fixed.add(sec.o_ax * nb)                                # axial anchor / fixed plate
    if bc == "CC":
        for m in range(sec.o_s, sec.nm):                    # shear, rotations: none at ends
            ends(m)
    free = np.array(sorted(set(range(ndof)) - fixed))
    # ---- end loads: q = 1 / perimeter on both end edges (1 N in total) ------
    q = 1.0 / sec.perimeter; f = np.zeros(ndof)
    (B0, BL), (D0, DL) = lb.end[0], lb.end[1]
    iuk = sec.wp @ sec.ukp; iuj = sec.wp @ sec.ujp
    for m in range(nip):
        f[m * nb:(m + 1) * nb] += q * iuk[m] * (D0 - DL)
    f[sec.o_ax * nb:(sec.o_ax + 1) * nb] += q * sec.perimeter * (B0 - BL)
    for j in range(sec.ns):
        m = sec.o_s + j; f[m * nb:(m + 1) * nb] += q * iuj[j] * (B0 - BL)
    Kr = K[free][:, free].tocsc(); Kr = 0.5 * (Kr + Kr.T)
    solve, release = _osfe_factor(Kr)
    full = "T" in sec.families
    if prestress == "uniform":
        Pxx = np.tile(sec.A11p / (sec.wp @ sec.A11p), (len(lb.xg), 1)); Pss = Pxs = 0 * Pxx
    else:
        a = np.zeros(ndof); a[free] = solve(f[free]); C = a.reshape(nm, nb)
        amp = np.stack([(lb.B[d] @ C.T) for d in range(3)], axis=2)      # (nx, nm, 3)
        strain = np.einsum("pamd,xmd->xpa", sec.E, amp, optimize=True)    # (nx, np, 8)
        Nres = -np.einsum("pca,xpa->xpc", sec.Cp[:, :3, :], strain, optimize=True)  # compression +
        Pxx, Pss, Pxs = Nres[..., 0], Nres[..., 1], Nres[..., 2]
    # ---- assemble Kg span by span (cross-section tensors at the span's Gauss points)
    Bd = [lb.B[d].toarray().reshape(lb.nspan, 4, nb) for d in range(3)]
    gmask = np.abs(sec.Gd[:, (0, 1, 2) if full else (1, 2)]).sum(axis=(0, 1, 2, 4))
    act = np.nonzero(gmask)[0]
    rows, cols, vals = [], [], []; Kg = None
    for j in range(lb.nspan):
        sl = slice(4 * j, 4 * j + 4); ids = np.arange(j, j + 4)
        Tg = _gbt_kg_tensors(sec, Pxx[sl], Pss[sl] if full else None, Pxs[sl] if full else None, full)
        blk = np.zeros((len(act), len(act), 4, 4))
        for (d, e), T in Tg.items():
            Bab = np.einsum("g,ga,gb->gab", lb.wx[sl], Bd[d][j][:, ids], Bd[e][j][:, ids])
            blk += np.tensordot(T[:, act][:, :, act], Bab, axes=([0], [0]))   # (k, l, a, b)
        blk = blk.transpose(0, 2, 1, 3)
        gi = (act[:, None] * nb + ids[None, :]).ravel()
        bv = blk.reshape(gi.size, gi.size); nz = np.nonzero(bv)
        rows.append(gi[nz[0]].astype(np.int32)); cols.append(gi[nz[1]].astype(np.int32))
        vals.append(bv[nz])
        if len(vals) >= 8:                                # merge duplicates every 8 spans
            Kg = _kg_flush(Kg, rows, cols, vals, ndof)    # (keeps the memory low)
    Kg = _kg_flush(Kg, rows, cols, vals, ndof)
    Gr = Kg[free][:, free].tocsc(); Gr = 0.5 * (Gr + Gr.T)
    Kinv = spla.LinearOperator(Kr.shape, matvec=solve, dtype=float)
    mu, V = spla.eigsh(Gr, k=nev, M=Kr, Minv=Kinv, which="LA", ncv=max(2 * nev + 16, 20),
                       tol=1e-10, maxiter=20000)
    release()
    ok = mu > 0; lam = 1.0 / mu[ok]; V = V[:, ok]; o = np.argsort(lam); lam, V = lam[o], V[:, o]
    vec = np.zeros(ndof); vec[free] = V[:, 0]; Cm = vec.reshape(nm, nb)
    res = dict(P=float(lam[0]), eigs=[float(v) for v in lam],
               participation=_gbt_participation(sec, lb, Cm), ndof=len(free), nmodes=nm,
               method="spline", bc=bc, families=sec.families)
    if return_mode:
        res.update(sec=sec, lb=lb, coef=Cm)
    return res


def _gbt_participation(sec, lb, C):
    """GBT modal participation of a buckling mode (as the paper's Fig. 16b),
    per mode group: p_g = int ||field_g(x)|| dx / sum over groups, where
    field_g is the displacement field produced by the group's modes alone
    (in-plane v, w for G, D, L, T; warping u for S; wall rotations for M),
    RMS over the cross-section.  Independent of how each group's basis is
    chosen or scaled.  The axial (pre-buckling) mode is excluded."""
    if lb.kind == "spline":
        A = np.asarray(lb.B[0] @ C.T)                                  # (nx, nm)
    else:
        A = np.column_stack([lb.Bf[sec.fam[m]][0] @ C[m] for m in range(sec.nm)])
    w = sec.wp / sec.wp.sum()
    grp = np.array(sec.group)
    nrm = lambda fld: float(lb.wx @ np.sqrt(np.einsum("p,xpc->x", w, fld ** 2)))
    out = {}
    # in-plane field: rigid (G), distortional (D) and transverse-extension (T)
    # parts by projection of the nodal displacements, the rest is local (L)
    nip = sec.nip
    Utot = sec.U @ A[:, :nip].T                                         # (2N, nx)
    Ftot = np.einsum("pck,xk->xpc", sec.VW, A[:, :nip])
    UT = sec._Torth @ (sec._Torth.T @ Utot) if "T" in sec.families else 0 * Utot
    UC = Utot - UT
    # GBT split: the global and distortional parts are fixed by the fold /
    # branch-node displacements (local modes leave them at rest); a section
    # with a single branch node (T) also uses its free tips for the rigid fit
    nat = sec.junctions if len(sec.junctions) >= 2 else sec.junctions + sec.tips
    J = np.array([2 * n + c for n in nat for c in (0, 1)])
    UG = sec._Gorth @ np.linalg.lstsq(sec._Gorth[J], UC[J], rcond=None)[0]
    if sec._Dorth.shape[1]:
        UD = sec._Dorth @ np.linalg.lstsq(sec._Dorth[J], (UC - UG)[J], rcond=None)[0]
    else:
        UD = 0 * UC
    FG, FD, FT = (np.einsum("pci,ix->xpc", sec.Pvw, U_) for U_ in (UG, UD, UT))
    out["G"] = nrm(FG); out["D"] = nrm(FD); out["L"] = nrm(Ftot - FG - FD - FT)
    if "T" in sec.families:
        out["T"] = nrm(FT)
    for g, F in (("S", lambda ks: sec.E[:, 0, ks, 1][:, None, :]),
                 ("M", lambda ks: sec.E[:, 6:8, ks, 0])):
        ks = np.nonzero(grp == g)[0]
        if len(ks):
            out[g] = nrm(np.einsum("pck,xk->xpc", F(ks), A[:, ks]))
    tot = sum(out.values()) or 1.0
    return {g: v / tot for g, v in out.items()}


def _gbt_sine(sec, L, harmonics, nev, return_mode):
    """Classical semi-analytical GBT-FSM (SS ends, uniform pre-stress, Eq. 8)."""
    import scipy.linalg as sla
    nm, nip = sec.nm, sec.nip
    if harmonics is None:                 # half-waves down to half the narrowest wall
        wl = np.zeros(sec.nw)
        for s in sec.S:
            wl[s["w"]] += s["L"]
        harmonics = list(range(1, int(np.clip(np.ceil(2.0 * L / wl.min()) + 2, 6, 160)) + 1))
    Kcs = _gbt_Kcs(sec)
    Pxx = sec.A11p / (sec.wp @ sec.A11p)
    Gt = _gbt_kg_tensors(sec, Pxx[None, :], full=False)          # {(1,1): (1, nm, nm)}

    def solve_set(ms):
        lb = _LongSine(L, ms); nb = lb.nb
        F = [lb.Bf[f] for f in sec.fam]
        K = np.zeros((nm * nb, nm * nb)); Kg = np.zeros_like(K)
        for (d, e), A in Kcs.items():
            for m, n in zip(*np.nonzero(A)):
                K[m * nb:(m + 1) * nb, n * nb:(n + 1) * nb] += A[m, n] * ((F[m][d] * lb.wx[:, None]).T @ F[n][e])
        for (d, e), T in Gt.items():
            T0 = T[0]
            for m, n in zip(*np.nonzero(T0)):
                Kg[m * nb:(m + 1) * nb, n * nb:(n + 1) * nb] += T0[m, n] * ((F[m][d] * lb.wx[:, None]).T @ F[n][e])
        mu, V = sla.eigh(0.5 * (Kg + Kg.T), 0.5 * (K + K.T))
        o = np.argsort(-mu); mu, V = mu[o], V[:, o]
        return 1.0 / mu[0], V[:, 0].reshape(nm, nb), lb, [1.0 / x for x in mu[:nev] if x > 0]
    single = {m: solve_set([m]) for m in harmonics}
    mstar = min(single, key=lambda m: single[m][0])
    P, C, lb, eigs = solve_set([m for m in harmonics if abs(m - mstar) <= 6])
    res = dict(P=float(P), eigs=eigs, participation=_gbt_participation(sec, lb, C), m=mstar,
               method="sine", bc="SS", nmodes=nm, families=sec.families)
    if return_mode:
        res.update(sec=sec, lb=lb, coef=C)
    return res


# ==============================================================================
# 24.  GOVERNING SOLVER SELECTION  (GBT-FSM of section 23 / shell strip of 20)
# ------------------------------------------------------------------------------
#   GOVERNING_SOLVER = "gbt"   : shear-deformable GBT (Silvestre & Camotim 2013)
#                                solved by spline finite strips (section 23)
#                      "strip" : shell-type Mindlin spline finite strip (sec. 20)
#   GBT_FAMILIES     = "CSTMR" : default -- full hierarchy (120 Abaqus cases:
#                                mean |error| 0.94 %, ~8 s per case)
#                      "CS"    : the paper (conventional + shear modes)
#                      "CSR"   : + nodal-rotation modes (completes the transverse
#                                interpolation; same converged answer as "CS",
#                                reached with far fewer nodes)
#                      "CSMR"  : + transverse-shear (Mindlin) modes
#                      "CSTMR" : + transverse-extension modes (full kinematics)
#   GBT_SEG, GBT_HX           : strip width across the section / spline span
#                               along the member [mm]
# Command line:  --solver=strip | --gbt=CS|CSM|CSTM | --gbt-sine (SS only)
# ==============================================================================
GOVERNING_SOLVER = "gbt"
GBT_FAMILIES = "CSTMR"
GBT_METHOD = "spline"
GBT_SEG = 2.5                 # target strip width across the section [mm]
GBT_NMIN = 8                  # at least this many strips per wall
GBT_HX = 10.0                 # spline span along the member [mm] (<= L/20)

_GBT_GROUP_NAMES = {"G": "global", "D": "distortional", "L": "local", "S": "shear (membrane)",
                    "T": "transverse extension", "M": "transverse shear (Mindlin)"}


def _gbt_opts(cfg):
    return (str(cfg.get("solver", GOVERNING_SOLVER)).lower(),
            str(cfg.get("gbt_families", GBT_FAMILIES)).upper(),
            str(cfg.get("gbt_method", GBT_METHOD)).lower(),
            float(cfg.get("gbt_seg", GBT_SEG)), float(cfg.get("gbt_hx", GBT_HX)))


def gbt_member_buckling(cfg, L, bc="SS", **kw):
    """P_cr of the member by the GBT-FSM with the defaults of section 24."""
    _s, fam, meth, seg, hx = _gbt_opts(cfg)
    fam = kw.pop("families", fam); meth = kw.pop("method", meth)
    seg = kw.pop("seg", seg); hx = kw.pop("hx", hx)
    if meth == "sine" and bc.upper() != "SS":
        meth = "spline"                                   # sine strips are SS only
    sec = GBTSection(cfg, seg, fam, int(cfg.get("gbt_nmin", GBT_NMIN)))
    return gbt_fsm_solve(sec, float(L), bc, method=meth, hx=min(hx, float(L) / 20.0), **kw)


# (_unified_governing, _unified_sweep, _plot_gbt_mode, _plot_gbt_participation and
#  _section_study, which use this selection, are in section 21 above.)


# ==============================================================================
# 25.  OMEGA-STIFFENED PANEL  -  SHEAR-DEFORMABLE GBT  (governing omega solver)
# ------------------------------------------------------------------------------
# The GBT model of the omega panel (equations and verification: the note
# "GBT_omega_panel_formulation_and_PINN_physics.pdf").  In short:
#   * walls of section 14e (skin, stacked feet, webs, crowns) with their own
#     A, B, D, As; the web roots sit at the tie offset c = (ts + tf)/2 and are
#     joined to the skin by short LINK walls (stiff across the section, no
#     axial or twisting stiffness = a row of rigid offset links);
#   * GBT mode families G, D, L (+ R nodal rotations), T, A, S, M built on the
#     discretised cross-section (segments <= bs/25); the SIDE-EDGE conditions
#     (symmetry, or periodic u(x=-2bs) = u(x=0)) are built into every mode;
#   * amplitude functions in cubic B-splines, spans min(15, max(2, W/40), L/10);
#   * SS ends: U2 = 0 at every end node;  CC ends: clamped edges with the
#     uniform edge load (OMEGA_CC_END = "edge-load");  1 N/mm on every end
#     edge, 2 N/mm on a foot at the height of its resultant;
#   * first-order pre-buckling state (all membrane resultants), full Kg.
# 60 Abaqus cases (symmetry sides): mean |error| 0.58 %, max 1.48 %;
# 30 Abaqus periodic SS cases: mean 0.57 %, max 1.25 %.
# The mode-reduction tools at the end (hierarchical_modes, truncated_section,
# enriched_section) are research tools for the PINN basis; the solver itself
# always uses the full mode basis.
# ==============================================================================
import scipy.linalg as sla


class _ThisModule:
    """Attribute access to this file's own functions: the omega GBT routines
    call their helpers as mm._osfe_laminate, m29._hermite, ..."""
    def __getattr__(self, name):
        return globals()[name]


_THIS = _ThisModule()
OMEGA_SOLVER = "gbt"            # governing omega solver: "gbt" (section 25) or "strip" (22)
OMEGA_GBT_SEG = 25.0            # wall segments of at most bs / OMEGA_GBT_SEG [mm]
OMEGA_GBT_FAMILIES = "CSTMR"    # full GBT mode hierarchy
OMEGA_GBT_LINK_K = 1.0e3        # stiffness factor of the link walls (result independent)


def omega_gbt_geometry(mm, panel, link_k=1.0e3):
    """Natural nodes and walls of the panel for GBT.  Returns X (nodes),
    walls [(i, j, C6, As, edges, yload, part)] with the wall direction of
    the strip model (skin left -> right, webs/crown out of the hat), and the
    side nodes (x = -W, x = 0)."""
    info = panel["info"]; nd0 = panel["nodes"]; clear = info["clear"]; mat = info["material"]
    As_s = mm._osfe_laminate(mat, info["pskin"])[3]
    As_t = mm._osfe_laminate(mat, info["pstiff"])[3]
    As_f = mm._osfe_laminate(mat, list(info["pskin"]) + list(info["pstiff"]))[3]
    X = [np.asarray(p, float) for p in nd0]
    walls = []; root = {}
    def C6_of(ABD):
        A, B, D = ABD[0], ABD[1], ABD[2]
        C = np.zeros((6, 6)); C[:3, :3] = A; C[:3, 3:] = B; C[3:, :3] = B; C[3:, 3:] = D
        return C
    # stiff link laminate (from the skin): stiff in everything but axial
    sk = next(w for w in panel["walls"] if w["part"] == "skin")
    A0, D0 = sk["ABD"][0], sk["ABD"][2]
    CL = np.zeros((6, 6))
    CL[0, 0] = 1e-3 * A0[0, 0]; CL[1, 1] = link_k * A0[1, 1]; CL[2, 2] = link_k * A0[2, 2]
    CL[3, 3] = 1e-3 * D0[0, 0]; CL[4, 4] = link_k * D0[1, 1]; CL[5, 5] = 1e-3 * D0[2, 2]   # no twisting stiffness (a line of rigid links, not a plate)
    AsL = link_k * As_s
    def root_of(i):
        if i not in root:
            X.append(np.array([nd0[i][0], clear])); root[i] = len(X) - 1
            walls.append((i, root[i], CL, AsL, 0.0, 0.0, "link",
                          dict(part="link", ABD=(CL[:3, :3], CL[:3, 3:], CL[3:, 3:], clear),
                               edges=0.0, yload=0.0)))
        return root[i]
    for w in panel["walls"]:
        i0, i1 = w["i"], w["j"]
        As = {"skin": As_s, "foot": As_f}.get(w["part"], As_t)
        if w["part"] == "web":
            if abs(nd0[i0][1]) < 1e-12: i0 = root_of(i0)
            if abs(nd0[i1][1]) < 1e-12: i1 = root_of(i1)
        walls.append((i0, i1, C6_of(w["ABD"]), As, float(w["edges"]), float(w["yload"]), w["part"], w))
    return np.array(X), walls, list(panel["sides"])


class OmegaGBTSection:
    """GBT cross-section analysis of the omega panel (families subset of
    'CSTMR'; 'CS' = Silvestre & Camotim 2013)."""

    def __init__(self, mm, m29, panel, seg=3.0, families="CSTMR", nmin=2,
                 side="symmetry", link_k=1.0e3, nint=None):
        """nint: intermediate nodes of every panel wall (order of panel['walls']);
        None = segments of at most seg [mm] (at least nmin per wall)."""
        self.mm, self.m29 = mm, m29
        self.families = families; self.side = side
        X0, walls, sides = omega_gbt_geometry(mm, panel, link_k)
        self.walls = walls
        pts = [p for p in X0]; segs = []; self.C6 = []; self.As = []
        self.wedge = []; self.wy = []; self.wpart = []; self.wsrc = []
        for iw, (i, j, C6, As, edges, yload, part, src) in enumerate(walls):
            self.C6.append(C6); self.As.append(As); self.wsrc.append(src)
            self.wedge.append(edges); self.wy.append(yload); self.wpart.append(part)
            P0, P1 = pts[i], pts[j]; Lw = float(np.hypot(*(P1 - P0)))
            if part == "link":
                n = 1
            elif nint is not None:
                n = 1 + max(0, int(nint[next(k_ for k_, w_ in enumerate(panel["walls"]) if w_ is src)]))
            else:
                n = max(nmin, int(np.ceil(Lw / seg - 1e-9)))
            chain = [i]
            for q in range(1, n):
                pts.append(P0 + (P1 - P0) * q / n); chain.append(len(pts) - 1)
            chain.append(j)
            for a, b in zip(chain[:-1], chain[1:]):
                segs.append((a, b, iw))
        self.X = np.array(pts); self.N = N = len(pts); self.nw = len(walls)
        self.sides = sides
        from collections import Counter
        val = Counter()
        for (i, j, *_r) in walls:
            val[i] += 1; val[j] += 1
        self.junctions = sorted(n for n in val if val[n] >= 2)
        self.tips = sorted(n for n in val if val[n] == 1)
        self.S = []
        for (a, b, iw) in segs:
            d = self.X[b] - self.X[a]; L = float(np.hypot(*d)); t = d / L
            self.S.append(dict(a=a, b=b, w=iw, L=L, t=t, n=np.array([-t[1], t[0]])))
        self.perimeter = sum(s["L"] * self.wedge[s["w"]] for s in self.S)   # loaded edge length
        condense = m29.condense_Nss
        self.C8 = []
        for iw in range(self.nw):
            C = np.zeros((8, 8))
            if "T" in families:
                C[:6, :6] = self.C6[iw]
            else:
                C5 = condense(self.C6[iw]); k = [0, 2, 3, 4, 5]
                C[np.ix_(k, k)] = C5
            C[6:, 6:] = self.As[iw]
            self.C8.append(C)
        self._build_modes()
        self._sample()

    # ------------------------------------------------------------------ modes
    def _constraint_rows(self, U, Th):
        """Rows of the side-edge constraints acting on stacked [U; Th]."""
        N = self.N; sL, sR = self.sides; rows = []
        def e(k):
            v = np.zeros(3 * N); v[k] = 1.0; return v
        if self.side == "symmetry":
            for s in (sL, sR):
                rows.append(e(2 * s)); rows.append(e(2 * N + s))
        else:
            for c in (0, 1):
                rows.append(e(2 * sL + c) - e(2 * sR + c))
            rows.append(e(2 * N + sL) - e(2 * N + sR))
        return np.array(rows)

    def _build_modes(self):
        families = self.families; N = self.N; m29 = self.m29
        E = np.zeros((len(self.S), 2 * N))
        for r, s in enumerate(self.S):
            E[r, 2 * s["b"]:2 * s["b"] + 2] += s["t"]; E[r, 2 * s["a"]:2 * s["a"] + 2] -= s["t"]
        Z = sla.null_space(E)
        c = self.X.mean(axis=0)
        G = np.zeros((2 * N, 3)); G[0::2, 0] = 1.0; G[1::2, 1] = 1.0
        G[0::2, 2] = -(self.X[:, 1] - c[1]); G[1::2, 2] = self.X[:, 0] - c[0]
        G = sla.orth(G)
        jset = set(self.junctions); Lm = []
        for n in range(N):
            if n in jset:
                continue
            sg = next(s_ for s_ in self.S if n in (s_["a"], s_["b"]))
            ev = np.zeros(2 * N); ev[2 * n:2 * n + 2] = sg["n"]; Lm.append(ev)
        Lm = np.array(Lm).T if Lm else np.zeros((2 * N, 0))
        Gk = []
        for k in range(3):
            B_ = np.hstack([Lm] + Gk) if (Lm.shape[1] or Gk) else np.zeros((2 * N, 0))
            r = G[:, k] - (B_ @ np.linalg.lstsq(B_, G[:, k], rcond=None)[0] if B_.shape[1] else 0)
            if np.linalg.norm(r) > 1e-8 * np.linalg.norm(G[:, k]):
                Gk.append(G[:, [k]])
        self._Gorth = G.copy(); self._Torth = sla.orth(E.T)
        G = np.hstack(Gk)
        GL = np.hstack([G, Lm])
        Dm = Z @ sla.null_space(GL.T @ Z) if GL.shape[1] < Z.shape[1] else Z[:, :0]
        Dm = sla.orth(Dm) if Dm.shape[1] else Dm
        self._Dorth = Dm
        blocks = [("G", G), ("D", Dm), ("L", Lm)]
        if "T" in families:
            blocks.append(("T", E.T.copy()))
        self.R = self._frame_R()
        # in-plane modes with their nodal rotations, then the side constraints
        Ub = []; Tb = []; grp = []
        for g, B in blocks:
            if B.shape[1]:
                Ub.append(B); Tb.append(self.R @ B if "R" not in families else 0 * (self.R @ B))
                grp += [g] * B.shape[1]
        if "R" in families:
            Ub.append(np.zeros((2 * N, N))); Tb.append(np.eye(N)); grp += ["L"] * N
        U = np.hstack(Ub); Th = np.hstack(Tb); grp = np.array(grp)
        Cr = self._constraint_rows(U, Th)
        S3 = np.vstack([U, Th])
        # project group by group, then complete with the coupled remainder
        cols = []; glab = []
        for g in ("G", "D", "L", "T"):
            idx = np.nonzero(grp == g)[0]
            if not len(idx):
                continue
            Bg = S3[:, idx]
            CB = Cr @ Bg
            free_c = np.nonzero(np.abs(CB).sum(axis=0) < 1e-12)[0]      # untouched columns
            inv_c = np.nonzero(np.abs(CB).sum(axis=0) >= 1e-12)[0]
            if len(free_c):
                cols.append(Bg[:, free_c]); glab += [g] * len(free_c)
            if len(inv_c):
                Zg = sla.null_space(CB[:, inv_c])
                if Zg.shape[1]:
                    Bz = Bg[:, inv_c] @ Zg; cols.append(Bz); glab += [g] * Bz.shape[1]
        Bc = np.hstack(cols)
        Zall = sla.null_space(Cr @ S3)
        Vall = sla.orth(S3 @ Zall)
        Pc = Vall - Bc @ np.linalg.lstsq(Bc, Vall, rcond=None)[0]
        if np.linalg.norm(Pc) > 1e-8:
            Ex = sla.orth(Pc, rcond=1e-8)
            if Ex.shape[1]:
                cols.append(Ex)
                for k in range(Ex.shape[1]):
                    u = Ex[:2 * N, k]
                    glab.append("T" if np.linalg.norm(E @ u) > 1e-6 * max(np.linalg.norm(u), 1e-30) else "D")
        S3 = np.hstack(cols); grp = glab
        U = S3[:2 * N]; Th = S3[2 * N:]
        for k in range(U.shape[1]):
            u = U[:, k].reshape(N, 2); sc = max(np.hypot(u[:, 0], u[:, 1]).max(), np.abs(Th[:, k]).max(), 1e-30)
            U[:, k] /= sc; Th[:, k] /= sc
        self.U = U; self.Th = Th; self.nip = U.shape[1]
        # warping of G and D modes (Vlasov, least squares; periodic: u_L = u_R)
        Aw = np.zeros((len(self.S), N)); wts = np.zeros(len(self.S)); c_ax = np.zeros(N)
        for r, s in enumerate(self.S):
            Aw[r, s["b"]] += 1.0 / s["L"]; Aw[r, s["a"]] -= 1.0 / s["L"]
            C5 = m29.condense_Nss(self.C6[s["w"]])
            wts[r] = C5[1, 1] * s["L"]
            c_ax[[s["a"], s["b"]]] += 0.5 * C5[0, 0] * s["L"]
        W = np.sqrt(wts)[:, None]; AW = Aw * W
        cons = [c_ax[None, :]]
        if self.side == "periodic":
            sL, sR = self.sides; ev = np.zeros((1, N)); ev[0, sL] = 1.0; ev[0, sR] = -1.0; cons.append(ev)
        Cc = np.vstack(cons); nc = Cc.shape[0]
        Msys = np.block([[AW.T @ AW, Cc.T], [Cc, np.zeros((nc, nc))]])
        self.Uw = np.zeros((N, self.nip))
        for k in range(self.nip):
            if grp[k] not in ("G", "D"):
                continue
            Ud = U[:, k].reshape(N, 2)
            v = np.array([s["t"] @ (0.5 * (Ud[s["a"]] + Ud[s["b"]])) for s in self.S])
            self.Uw[:, k] = np.linalg.lstsq(Msys, np.r_[-AW.T @ (W[:, 0] * v), np.zeros(nc)],
                                            rcond=None)[0][:N]
        self.c_ax = c_ax
        # shear (warping hat) modes
        sL, sR = self.sides
        if "S" in families:
            Bs = np.eye(N)
            if self.side == "periodic":
                Bs[:, sR] += Bs[:, sL]; Bs = np.delete(Bs, sL, axis=1)
            drop = next(n for n in range(N) if n not in (sL, sR))
            keep = [k for k in range(Bs.shape[1]) if not (Bs[drop, k] == 1.0 and Bs[:, k].sum() == 1.0)]
            Bs = Bs[:, keep]
        else:
            Bs = np.zeros((N, 0))
        self.Us = Bs; self.ns = Bs.shape[1]
        # Mindlin rotation modes
        if "M" in families:
            cols = []; owner = []
            for n in range(N):
                ts = np.array([s_["t"] for s_ in self.S if n in (s_["a"], s_["b"])])
                for vec in sla.orth(ts.T).T:
                    ev = np.zeros(2 * N); ev[2 * n:2 * n + 2] = vec; cols.append(ev); owner.append(n)
            Om = np.array(cols).T
            Ms = np.eye(N)
            if self.side == "periodic":
                iL = [k for k, o in enumerate(owner) if o == sL]; iR = [k for k, o in enumerate(owner) if o == sR]
                for a_, b_ in zip(iL, iR):
                    Om[:, b_] += Om[:, a_]
                Om = np.delete(Om, iL, axis=1)
                Ms[:, sR] += Ms[:, sL]; Ms = np.delete(Ms, sL, axis=1)
            else:                                             # symmetry: no rotation about Z
                Ms = np.delete(Ms, [sL, sR], axis=1)
            self.Om = Om; self.Ms = Ms
        else:
            self.Om = np.zeros((2 * N, 0)); self.Ms = np.zeros((N, 0))
        self.nmx = self.Om.shape[1]; self.nms = self.Ms.shape[1]
        self.o_ax = self.nip; self.o_s = self.nip + 1; self.o_mx = self.o_s + self.ns
        self.o_ms = self.o_mx + self.nmx; self.nm = self.o_ms + self.nms
        self.group = list(grp) + ["A"] + ["S"] * self.ns + ["M"] * (self.nmx + self.nms)
        self.fam = (["phi"] * self.nip + ["psi"] * (1 + self.ns + self.nmx) + ["phi"] * self.nms)

    def _frame_R(self):
        N = self.N; Kth = np.zeros((N, N)); KtU = np.zeros((N, 2 * N))
        for s in self.S:
            k = self.m29._beam_k(self.C6[s["w"]][4, 4], s["L"])
            ends = (s["a"], s["b"])
            for r, na in zip((1, 3), ends):
                for c_, nb_ in zip((1, 3), ends):
                    Kth[na, nb_] += k[r, c_]
                for c_, nb_ in zip((0, 2), ends):
                    KtU[na, 2 * nb_:2 * nb_ + 2] += k[r, c_] * s["n"]
        if self.side == "symmetry":                          # theta_Z = 0 at the sides
            for s in self.sides:
                Kth[s, :] = 0; Kth[:, s] = 0; Kth[s, s] = 1.0; KtU[s, :] = 0
        else:                                               # theta_L = theta_R
            sL, sR = self.sides
            Kth[sR, :] += Kth[sL, :]; Kth[:, sR] += Kth[:, sL]; KtU[sR, :] += KtU[sL, :]
            Kth[sL, :] = 0; Kth[:, sL] = 0; Kth[sL, sL] = 1.0; Kth[sL, sR] = -1.0; KtU[sL, :] = 0
        return -np.linalg.solve(Kth, KtU)

    def _sample(self):
        hermite = self.m29._hermite
        gp, gw = np.polynomial.legendre.leggauss(4); xi = 0.5 * (gp + 1.0)
        nm, nip, N = self.nm, self.nip, self.N
        th = self.Th
        npnt = 4 * len(self.S)
        E = np.zeros((npnt, 8, nm, 3)); Gd = np.zeros((npnt, 3, 2, nm, 3))
        wp = np.zeros(npnt); C8 = np.zeros((npnt, 8, 8)); VW = np.zeros((npnt, 2, nip))
        qp = np.zeros(npnt); ep = np.zeros(npnt); p = 0
        for s in self.S:
            a, b, L = s["a"], s["b"], s["L"]
            Ua = self.U[2 * a:2 * a + 2]; Ub = self.U[2 * b:2 * b + 2]
            va, vb = s["t"] @ Ua, s["t"] @ Ub; wa, wb = s["n"] @ Ua, s["n"] @ Ub
            nod = np.vstack([wa, th[a], wb, th[b]])
            iw = s["w"]
            yref = 0.5 * (self.X[a][1] + self.X[b][1])
            for q in range(4):
                H, dH, ddH = hermite(xi[q], L)
                w, ws, wss = H @ nod, dH @ nod, ddH @ nod
                v = (1 - xi[q]) * va + xi[q] * vb; vs = (vb - va) / L
                uk = (1 - xi[q]) * self.Uw[a] + xi[q] * self.Uw[b]; uks = (self.Uw[b] - self.Uw[a]) / L
                e = E[p]; g = Gd[p]; I = slice(0, nip)
                VW[p, 0] = v; VW[p, 1] = w
                e[0, I, 2] = uk; e[1, I, 0] = vs; e[2, I, 1] = uks + v
                e[3, I, 2] = -w; e[4, I, 0] = -wss; e[5, I, 1] = -2.0 * ws
                g[0, 0, I, 2] = uk; g[0, 1, I, 1] = uks
                g[1, 0, I, 1] = v;  g[1, 1, I, 0] = vs
                g[2, 0, I, 1] = w;  g[2, 1, I, 0] = ws
                e[0, self.o_ax, 1] = 1.0; g[0, 0, self.o_ax, 1] = 1.0
                if self.ns:
                    J = slice(self.o_s, self.o_s + self.ns)
                    uj = (1 - xi[q]) * self.Us[a] + xi[q] * self.Us[b]; ujs = (self.Us[b] - self.Us[a]) / L
                    e[0, J, 1] = uj; e[2, J, 0] = ujs
                    g[0, 0, J, 1] = uj; g[0, 1, J, 0] = ujs
                if self.nmx:
                    ha, hb = 1 - xi[q], xi[q]
                    ta = s["t"] @ self.Om[2 * a:2 * a + 2]; tb = s["t"] @ self.Om[2 * b:2 * b + 2]
                    cx = ha * ta + hb * tb; cxs = (tb - ta) / L
                    MX = slice(self.o_mx, self.o_mx + self.nmx)
                    e[3, MX, 1] = cx; e[5, MX, 0] = cxs; e[6, MX, 0] = cx
                if self.nms:
                    ha, hb = 1 - xi[q], xi[q]
                    cs_ = ha * self.Ms[a] + hb * self.Ms[b]; css = (self.Ms[b] - self.Ms[a]) / L
                    MS = slice(self.o_ms, self.o_ms + self.nms)
                    e[4, MS, 0] = css; e[5, MS, 1] = cs_; e[7, MS, 0] = cs_
                wp[p] = 0.5 * L * gw[q]; C8[p] = self.C8[iw]
                qp[p] = self.wedge[iw]
                ep[p] = (self.wy[iw] - yref) if self.wpart[iw] == "foot" else 0.0
                p += 1
        self.E, self.Gd, self.wp, self.Cp, self.VW = E, Gd, wp, C8, VW
        self.qp, self.ep = qp, ep
        Pvw = np.zeros((npnt, 2, 2 * N)); p = 0; Rf = self.R
        for s in self.S:
            a, b, L = s["a"], s["b"], s["L"]
            for q in range(4):
                H, _dH, _ddH = hermite(xi[q], L)
                Pvw[p, 0, 2 * a:2 * a + 2] += (1 - xi[q]) * s["t"]; Pvw[p, 0, 2 * b:2 * b + 2] += xi[q] * s["t"]
                Pvw[p, 1, 2 * a:2 * a + 2] += H[0] * s["n"]; Pvw[p, 1, 2 * b:2 * b + 2] += H[2] * s["n"]
                Pvw[p, 1, :] += H[1] * Rf[a] + H[3] * Rf[b]
                p += 1
        self.Pvw = Pvw
        self.A11p = np.array([self.m29.condense_Nss(self.C6[s["w"]])[0, 0] for s in self.S for _q in range(4)])
        self.ukp = E[:, 0, :nip, 2]; self.ujp = E[:, 0, self.o_s:self.o_s + self.ns, 1]


# ----------------------------------------------------------------------------
def _omega_kg_pattern(sec, full):
    """Sparse cross-section factors of the geometric stiffness (cached on the
    section; independent of L and of the pre-stress).

    Kg_(d,e)[x, a, b] = sum_t sum_p P_t[x, p] wp[p] sum_c G[p,c,da_t,a,d] G[p,c,db_t,b,e]
    (the same terms as _gbt_kg_tensors).  At a cross-section point p only the
    modes with support there are non-zero (about 25 of several hundred), so
    the sum is stored as sparse matrices M_t^(d,e) (n_points x n_pairs) over
    the mode pairs (a, b) that share a point."""
    cache = sec.__dict__.setdefault("_kg_pattern", {})
    if full in cache:
        return cache[full]
    comps = (0, 1, 2) if full else (1, 2)
    terms = [(0, 0, 0)] + ([(1, 1, 1), (2, 0, 1), (2, 1, 0)] if full else [])  # (P_t, da, db)
    G = sec.Gd[:, comps]                                   # (np, nc, 2, nm, 3)
    npnt = G.shape[0]
    act = [np.nonzero(np.abs(G[p]).sum(axis=(0, 1, 3)))[0] for p in range(npnt)]
    if np.mean([len(x) for x in act]) > 0.3 * sec.nm:     # global modes: dense is faster
        cache[full] = None
        return None
    P_i = np.concatenate([np.full(len(a) ** 2, p) for p, a in enumerate(act)])
    A_i = np.concatenate([np.repeat(a, len(a)) for a in act])
    B_i = np.concatenate([np.tile(a, len(a)) for a in act])
    ukey, inv = np.unique(A_i.astype(np.int64) * sec.nm + B_i, return_inverse=True)
    pa, pb = ukey // sec.nm, ukey % sec.nm
    w = sec.wp[P_i]
    M = {}
    for d in range(3):
        for e in range(3):
            for (t, da, db) in terms:
                v = w * np.einsum("ic,ic->i", G[P_i, :, da, A_i, d], G[P_i, :, db, B_i, e])
                if not np.any(v):
                    continue
                Mt = sp.csr_matrix((v, (P_i, inv)), shape=(npnt, len(ukey)))
                M.setdefault((d, e), []).append((t, Mt))
    cache[full] = (pa, pb, M)
    return cache[full]


def _omega_kg_dense(sec, lb, pre, full, ndof):
    """Geometric stiffness, dense per span (the V31 assembly) -- used when the
    modes are global (reduced hierarchical basis), where the sparse pair
    pattern does not pay off."""
    m29 = sec.m29
    Pxx, Pss, Pxs = pre; nb = lb.nb
    Bd = [lb.B[d].toarray().reshape(lb.nspan, 4, nb) for d in range(3)]
    gmask = np.abs(sec.Gd[:, (0, 1, 2) if full else (1, 2)]).sum(axis=(0, 1, 2, 4))
    act = np.nonzero(gmask)[0]
    rows, cols, vals = [], [], []; Kg = None
    for j in range(lb.nspan):
        sl = slice(4 * j, 4 * j + 4); ids = np.arange(j, j + 4)
        Tg = m29._gbt_kg_tensors(sec, Pxx[sl], Pss[sl] if full else None, Pxs[sl] if full else None, full)
        blk = np.zeros((len(act), len(act), 4, 4))
        for (d, e), T in Tg.items():
            Bab = np.einsum("g,ga,gb->gab", lb.wx[sl], Bd[d][j][:, ids], Bd[e][j][:, ids])
            blk += np.tensordot(T[:, act][:, :, act], Bab, axes=([0], [0]))
        blk = blk.transpose(0, 2, 1, 3)
        gi = (act[:, None] * nb + ids[None, :]).ravel()
        bv = blk.reshape(gi.size, gi.size); nz = np.nonzero(bv)
        rows.append(gi[nz[0]].astype(np.int32)); cols.append(gi[nz[1]].astype(np.int32))
        vals.append(bv[nz])
        if len(vals) >= 8:                                # merge duplicates every 8 spans
            Kg = _kg_flush(Kg, rows, cols, vals, ndof)    # (keeps the memory low)
    return _kg_flush(Kg, rows, cols, vals, ndof)


def _omega_kg_assemble(sec, lb, pre, full, ndof, chunk=8):
    """Geometric stiffness of the GBT omega model for the first-order
    pre-stress pre = (Pxx, Pss, Pxs) at the Gauss points.  Same integrand as
    _gbt_kg_tensors with the span-wise B-spline assembly, evaluated only on
    the mode pairs that share a cross-section point (identical result,
    several times faster)."""
    pat = _omega_kg_pattern(sec, full)
    if pat is None:
        return _omega_kg_dense(sec, lb, pre, full, ndof)
    pa, pb, M = pat
    nb, ns = lb.nb, lb.nspan
    Bd = [lb.B[d].toarray().reshape(ns, 4, nb) for d in range(3)]
    jj = np.arange(ns)
    Bloc = [Bd[d][jj[:, None, None], np.arange(4)[None, :, None],
                  jj[:, None, None] + np.arange(4)[None, None, :]] for d in range(3)]  # (ns, 4, 4)
    wx = np.asarray(lb.wx).reshape(ns, 4)
    loc = np.arange(4)
    Kg = None
    for s0 in range(0, ns, chunk):
        s1 = min(ns, s0 + chunk); nsl = s1 - s0
        blk = np.zeros((nsl, len(pa), 4, 4))
        for (d, e), lst in M.items():
            T = None
            for (t, Mt) in lst:
                Tt = (Mt.T @ pre[t][4 * s0:4 * s1].T).T            # (4 nsl, n_pairs)
                T = Tt if T is None else T + Tt
            Bab = np.einsum("sg,sgi,sgj->sgij", wx[s0:s1], Bloc[d][s0:s1], Bloc[e][s0:s1])
            blk += np.einsum("sgk,sgij->skij", T.reshape(nsl, 4, -1), Bab, optimize=True)
        sidx = np.arange(s0, s1)[:, None, None, None]
        r = pa[None, :, None, None] * nb + sidx + loc[None, None, :, None]
        c = pb[None, :, None, None] * nb + sidx + loc[None, None, None, :]
        r = np.broadcast_to(r, blk.shape).ravel(); c = np.broadcast_to(c, blk.shape).ravel()
        T = sp.coo_matrix((blk.ravel(), (r, c)), shape=(ndof, ndof)).tocsr()
        Kg = T if Kg is None else Kg + T
    return Kg


def omega_gbt_solve(sec, L, bc="SS", hx=None, nev=2, keep=None, return_mode=False, pss_scale=1.0,
                    ss_xfree=True, pre=None, pre_only=False):
    """Critical line load N_cr [N/mm] of the omega panel with shear-
    deformable GBT; amplitude functions in cubic B-splines.  keep: optional
    list of retained mode indices (truncated GBT; the axial mode is always
    kept)."""
    m29 = sec.m29; mm = sec.mm
    bc = bc.upper(); L = float(L)
    if hx is None:
        hx = min(15.0, max(2.0, np.ptp(sec.X[:, 0]) / 40.0), L / 10.0)
    lb = m29._LongSpline(L, m29._sfs_breaks(L, hx, hx, L))
    nb, nm, nip = lb.nb, sec.nm, sec.nip; ndof = nm * nb
    K = None
    Kcs = m29._gbt_Kcs(sec)
    for (d, e), A in Kcs.items():
        X = (lb.B[d].T @ sp.diags(lb.wx) @ lb.B[e]).tocsr()
        T = sp.kron(sp.csr_matrix(A), X, format="csr")
        K = T if K is None else K + T
    # ---- end conditions -> reduction q = C r  (mode-major: m*nb + i) -------
    U = sec.U; N = sec.N
    rows_end = {}                                            # spline index -> rows on in-plane coefs
    fixed = set()
    def fix(m, i):
        fixed.add(m * nb + i)
    ip = np.arange(nip)
    if bc == "SS":
        Ry = U[1::2, :] if ss_xfree else U                    # UY of every node (UX too if not free)
        Ry = Ry[np.abs(Ry).sum(axis=1) > 0]
        r0 = Ry
        if sec.side == "periodic":
            r0 = np.vstack([Ry, U[2 * sec.sides[1], :][None, :]])   # U1 anchor
        rows_end[0] = r0; rows_end[nb - 1] = Ry
        fix(sec.o_ax, 0)                                      # U3 anchor
    elif bc == "CC":
        for m in range(nip):
            for i in (0, 1, nb - 2, nb - 1):
                fix(m, i)                                     # in-plane: phi = phi' = 0
        fix(sec.o_ax, 0)
        for j in range(sec.ns):
            fix(sec.o_s + j, 0)                               # warping held at z = 0 only
        for m in range(sec.o_mx, sec.nm):                     # rotations clamped
            fix(m, 0); fix(m, nb - 1)
    else:
        raise ValueError(bc)
    # empty modes: zero cross-section stiffness and zero strain gradients, so
    # no energy and no Kg (spurious Mindlin Omega vectors normal to a straight
    # web, left by round-off in the rank detection).  Holding them changes no
    # eigenvalue; left free they make K singular (PARDISO hides it by pivot
    # perturbation, SuperLU does not).
    kd = sum(np.abs(np.diag(A)) for A in Kcs.values())
    gz = np.abs(sec.Gd).sum(axis=(0, 1, 2, 4))
    for m in np.nonzero((kd <= 1e-12 * kd.max()) & (gz == 0))[0]:
        for i in range(nb):
            fix(int(m), i)
    if keep is not None:
        keepset = set(int(k) for k in keep) | {sec.o_ax}
        for m in range(nm):
            if m not in keepset:
                for i in range(nb):
                    fix(m, i)
    rr, cc, vv = [], [], []; nq = 0
    done = set()
    for i, R in rows_end.items():                             # constrained in-plane end coefs
        act = [m for m in ip if (m * nb + i) not in fixed]
        Rr = R[:, act]
        Z = sla.null_space(Rr) if Rr.size else np.eye(len(act))
        for jj in range(Z.shape[1]):
            for kk, m in enumerate(act):
                if abs(Z[kk, jj]) > 0:
                    rr.append(m * nb + i); cc.append(nq + jj); vv.append(Z[kk, jj])
        nq += Z.shape[1]
        done |= {m * nb + i for m in ip}
    for dof in range(ndof):
        if dof in fixed or dof in done:
            continue
        rr.append(dof); cc.append(nq); vv.append(1.0); nq += 1
    Cm = sp.csr_matrix((vv, (rr, cc)), shape=(ndof, nq))
    # ---- end loads (1 N/mm per edge, foot 2 N/mm at its resultant height) ---
    f = np.zeros(ndof)
    (B0, BL), (D0, DL) = lb.end[0], lb.end[1]
    w_q = sec.wp * sec.qp
    iuk = w_q @ sec.ukp                                       # int q u_k ds
    iwk = (w_q * sec.ep) @ sec.VW[:, 1, :]                    # int q e w_k ds (foot offset)
    iuj = w_q @ sec.ujp
    for m in range(nip):
        f[m * nb:(m + 1) * nb] += (iuk[m] - iwk[m]) * (D0 - DL)
    f[sec.o_ax * nb:(sec.o_ax + 1) * nb] += float(w_q.sum()) * (B0 - BL)
    for j in range(sec.ns):
        m = sec.o_s + j; f[m * nb:(m + 1) * nb] += iuj[j] * (B0 - BL)
    Ptot = float(w_q.sum())                                   # total end load for N = 1 N/mm
    Kr = (Cm.T @ K @ Cm).tocsc(); Kr = 0.5 * (Kr + Kr.T)
    solve, release = mm._osfe_factor(Kr)
    full = "T" in sec.families
    a = Cm @ solve(Cm.T @ f); C = a.reshape(nm, nb)
    amp = np.stack([(lb.B[d] @ C.T) for d in range(3)], axis=2)
    strain = np.einsum("pamd,xmd->xpa", sec.E, amp, optimize=True)
    Nres = -np.einsum("pca,xpa->xpc", sec.Cp[:, :3, :], strain, optimize=True)
    Pxx, Pss, Pxs = Nres[..., 0], pss_scale * Nres[..., 1], Nres[..., 2]
    if pre_only:                                               # first-order state only
        release()
        return dict(pre=(Pxx, Pss, Pxs))
    if pre is not None:                                        # imposed pre-buckling state
        Pxx, Pss, Pxs = pre
    Kg = _omega_kg_assemble(sec, lb, (Pxx, Pss, Pxs), full, ndof)
    Gr = (Cm.T @ Kg @ Cm).tocsc(); Gr = 0.5 * (Gr + Gr.T)
    Kinv = spla.LinearOperator(Kr.shape, matvec=solve, dtype=float)
    mu, V = spla.eigsh(Gr, k=nev, M=Kr, Minv=Kinv, which="LA", ncv=max(2 * nev + 16, 20),
                       tol=1e-8, maxiter=20000)
    release()
    ok = mu > 0; lam = 1.0 / mu[ok]; V = V[:, ok]; o = np.argsort(lam); lam, V = lam[o], V[:, o]
    # lam multiplies the reference load (total end load Ptot = 1 N/mm x edges)
    Cmode = (Cm @ V[:, 0]).reshape(nm, nb)
    Call = [(Cm @ V[:, j]).reshape(nm, nb) for j in range(V.shape[1])]
    res = dict(N=float(lam[0]), eigs=[float(v) for v in lam], ndof=nq, nmodes=nm, Ptot=Ptot,
               pre=(Pxx, Pss, Pxs),
               participation=m29._gbt_participation(sec, lb, Cmode))
    if return_mode:
        res.update(lb=lb, coef=Cmode, coefs=Call)
    return res


# ----------------------------------------------------------------------------
# Hierarchical (truncatable) GBT mode basis
# ----------------------------------------------------------------------------
def _gen_eig(Aa, Bb, rel=1e-13):
    """Aa x = lam Bb x for symmetric Aa >= 0, Bb >= 0 (Bb possibly singular):
    eigenvectors of the pencil restricted to range(Aa + Bb), ascending lam.
    Directions with Bb x = 0 get lam = inf (they come last)."""
    sa = max(np.abs(np.diag(Aa)).max(), 1e-300); sb = max(np.abs(np.diag(Bb)).max(), 1e-300)
    Aa = Aa / sa; Bb = Bb / sb                     # balance the two energies
    S = Aa + Bb; S = 0.5 * (S + S.T)
    w, V = np.linalg.eigh(S)
    k = w > rel * w.max()
    Q = V[:, k] / np.sqrt(w[k])                    # Q' S Q = I
    Ar = Q.T @ Aa @ Q; Ar = 0.5 * (Ar + Ar.T)
    nu, Y = np.linalg.eigh(Ar)                     # nu = lam / (1 + lam)
    nu = np.clip(nu, 0.0, 1.0)
    lam = np.where(nu < 1.0 - 1e-12, nu / np.maximum(1.0 - nu, 1e-300), np.inf) * (sa / sb)
    return lam, Q @ Y


def hierarchical_modes(sec, kref=None):
    """GBT cross-section modal decomposition of each family block of `sec`.

    Classical GBT diagonalises the stiffness matrices of the amplitude
    equation  C V'''' - D V'' + B V = q  family by family.  Here, for every
    block the modes solve the reference-wavelength pencil

        conventional (phi) block:  B phi = lam (k^2 D + k^4 C) phi
        shear-hat / Mindlin (psi) blocks:  A00 psi = lam (k^2 A11) psi
        Mindlin drilling (phi) block:  same as conventional

    with k = kref (default pi / pitch).  lam -> 0 are the rigid-body /
    low-energy modes; ordering by lam gives the GBT hierarchy (global ->
    distortional -> local -> higher).  Returns dict block -> (lam, T) with T
    the coefficient matrix in the block's original basis."""
    A = m29_Kcs(sec)
    if kref is None:
        kref = np.pi / (0.5 * np.ptp(sec.X[:, 0]))
    k2, k4 = kref ** 2, kref ** 4
    blocks = {"C": np.arange(0, sec.nip),
              "S": np.arange(sec.o_s, sec.o_s + sec.ns),
              "MX": np.arange(sec.o_mx, sec.o_mx + sec.nmx),
              "MS": np.arange(sec.o_ms, sec.o_ms + sec.nms)}
    Z = lambda d, e, ix: A[(d, e)][np.ix_(ix, ix)] if (d, e) in A else np.zeros((len(ix), len(ix)))
    out = {}
    for b, ix in blocks.items():
        if not len(ix):
            out[b] = (np.zeros(0), np.zeros((0, 0))); continue
        if b in ("C", "MS"):
            num = Z(0, 0, ix); den = k2 * Z(1, 1, ix) + k4 * Z(2, 2, ix)
        else:
            num = Z(0, 0, ix); den = k2 * Z(1, 1, ix)
        lam, T = _gen_eig(num, den)
        out[b] = (lam, T)
    out["_A"] = A; out["_kref"] = kref
    return out


def m29_Kcs(sec):
    return sec.m29._gbt_Kcs(sec)


def truncated_section(sec, kc, ks=None, kmx=0, kms=0, hm=None, kref=None,
                      coupled_shear=False, kfac=(1.0,), cidx=None):
    """Copy of `sec` whose mode basis is the first kc conventional, ks
    shear-hat and kmx / kms Mindlin GBT modes of each family (None = all).
    The axial mode is always kept.  Exact (all modes) -> same N_cr as sec."""
    import copy
    if hm is None:
        hm = hierarchical_modes(sec, kref)
    hm = dict(hm)
    def take(b, k):
        lam, T = hm[b]
        k = T.shape[1] if k is None else int(min(k, T.shape[1]))
        return T[:, :k], lam[:k]
    Tc, lc = take("C", kc); Ts, ls = take("S", ks); Tx, lx = take("MX", kmx); Tm, lm = take("MS", kms)
    if cidx is not None:                                    # an explicit list of conventional modes
        Tc = hm["C"][1][:, list(cidx)]; lc = hm["C"][0][list(cidx)]
    if coupled_shear and sec.ns and Tc.shape[1] < sec.nip:
        # warping that each kept conventional mode drags along: the shear-hat
        # field minimising the energy of a sin(kz) conventional / cos(kz)
        # warping pair (static condensation at the reference wavenumbers)
        A = hm.get("_A") or m29_Kcs(sec)
        iS = np.arange(sec.o_s, sec.o_s + sec.ns); iC = np.arange(0, sec.nip)
        Z = lambda d, e, a, b: A[(d, e)][np.ix_(a, b)] if (d, e) in A else 0.0
        cols = [Ts]
        for kk in (hm["_kref"] * f for f in kfac):
            Kss = Z(0, 0, iS, iS) + kk ** 2 * Z(1, 1, iS, iS)
            Ksc = kk * Z(0, 1, iS, iC) + kk ** 3 * Z(1, 2, iS, iC)
            cols.append(np.linalg.solve(Kss, Ksc @ Tc))
        Ts = sla.orth(np.hstack(cols), rcond=1e-10); ls = np.full(Ts.shape[1], np.nan)
    nm_old = sec.nm
    blocks = [(np.arange(0, sec.nip), Tc), (np.array([sec.o_ax]), np.eye(1)),
              (np.arange(sec.o_s, sec.o_s + sec.ns), Ts),
              (np.arange(sec.o_mx, sec.o_mx + sec.nmx), Tx),
              (np.arange(sec.o_ms, sec.o_ms + sec.nms), Tm)]
    ncol = sum(T.shape[1] for _, T in blocks)
    T = np.zeros((nm_old, ncol)); c = 0
    for ix, Tb in blocks:
        if Tb.shape[1]:
            T[np.ix_(ix, np.arange(c, c + Tb.shape[1]))] = Tb
        c += Tb.shape[1]
    r = copy.copy(sec)
    r.T_basis = T
    r.E = np.einsum("pcmd,mk->pckd", sec.E, T, optimize=True)
    r.Gd = np.einsum("pcamd,mk->pcakd", sec.Gd, T, optimize=True)
    r.nip = Tc.shape[1]; r.o_ax = r.nip; r.o_s = r.nip + 1; r.ns = Ts.shape[1]
    r.o_mx = r.o_s + r.ns; r.nmx = Tx.shape[1]; r.o_ms = r.o_mx + r.nmx; r.nms = Tm.shape[1]; r.nm = ncol
    r.U = sec.U @ Tc; r.Th = sec.Th @ Tc; r.Uw = sec.Uw @ Tc
    r.VW = np.einsum("pcm,mk->pck", sec.VW, Tc)
    r.ukp = sec.ukp @ Tc
    r.ujp = sec.ujp @ Ts if sec.ns else sec.ujp
    r.Us = sec.Us @ Ts if sec.ns else sec.Us
    # dominant original family of every new conventional mode (labels only)
    g0 = np.array(sec.group[:sec.nip]); lab = []
    for k in range(r.nip):
        wgt = {g: float(np.sum(Tc[g0 == g, k] ** 2)) for g in set(g0)}
        lab.append(max(wgt, key=wgt.get))
    r.group = lab + ["A"] + ["S"] * r.ns + ["M"] * (r.nmx + r.nms)
    r.fam = ["phi"] * r.nip + ["psi"] * (1 + r.ns + r.nmx) + ["phi"] * r.nms
    r.lam = dict(C=lc, S=ls, MX=lx, MS=lm)
    return r


# ----------------------------------------------------------------------------


# ---------------------------------------------------------- solver interface --
class _OmegaGBTView:
    """A GBT solution seen through the interface of the spline-strip model
    (nodes, strips, links, L, S[k]['b'], ndof) so that the omega result window
    and figures of section 14d draw GBT and strip results alike.  U_at returns
    the node displacements (N, 6, nx) in the strip order [U3, U1, U2, UR3, -, -]
    (axial, X, Y, rotation about z) at the stations xs."""

    def __init__(self, sec, lb, L, ndof):
        self.sec, self.lb, self.L, self.ndof = sec, lb, float(L), int(ndof)
        self.nodes = sec.X
        self.strips = [(s_["a"], s_["b"], None, None) for s_ in sec.S]
        self.S = [dict(b=s_["L"]) for s_ in sec.S]
        self.links = []
        self.breaks = lb.breaks

    def U_at(self, C, xs):
        from scipy.interpolate import BSpline
        xs = np.clip(np.asarray(xs, float), 0.0, self.L)
        spl = BSpline(self.lb.t, np.eye(self.lb.nb), 3)
        A0 = spl(xs) @ C.T; A1 = spl.derivative()(xs) @ C.T      # amplitudes (nx, nm)
        sec = self.sec; nip = sec.nip
        U = np.zeros((sec.N, 6, len(xs)))
        U[:, 1, :] = sec.U[0::2] @ A0[:, :nip].T                  # U1 (X)
        U[:, 2, :] = sec.U[1::2] @ A0[:, :nip].T                  # U2 (Y)
        U[:, 3, :] = sec.Th @ A0[:, :nip].T                       # rotation about z
        u = sec.Uw @ A1[:, :nip].T + A0[:, sec.o_ax][None, :]     # warping + axial
        if sec.ns:
            u = u + sec.Us @ A0[:, sec.o_s:sec.o_s + sec.ns].T    # shear modes
        U[:, 0, :] = u
        return U


_OMEGA_GBT_SECTIONS = []        # small cache of cross-section analyses
OMEGA_GBT_NINT = None           # intermediate nodes per panel wall (None = from OMEGA_GBT_SEG)
OMEGA_GBT_MODES = None          # conventional GBT modes used in the buckling analysis
                                # (None = all modes of the families = the validated model)


def omega_gbt_default_nint(panel):
    """Default intermediate nodes of every panel wall (order of panel['walls']):
    segments of at most bs / OMEGA_GBT_SEG, at least two per wall."""
    bs = float(panel["info"]["bs"]); seg = bs / OMEGA_GBT_SEG; nd = panel["nodes"]
    return [max(2, int(np.ceil(float(np.hypot(*(nd[w["j"]] - nd[w["i"]]))) / seg - 1e-9))) - 1
            for w in panel["walls"]]


def omega_gbt_section(panel):
    """GBT cross-section analysis of a panel (reused for every length)."""
    if OMEGA_SIDE_BC == "periodic" and OMEGA_PERIODIC_MACRO:
        raise ValueError("omega GBT: periodic side edges with free lateral expansion "
                         "exist in the strip solver only (--omega-periodic-free)")
    nint = tuple(OMEGA_GBT_NINT) if OMEGA_GBT_NINT is not None else None
    key = (OMEGA_SIDE_BC, OMEGA_GBT_SEG, OMEGA_GBT_FAMILIES, nint)
    for (p_, k_, sec) in _OMEGA_GBT_SECTIONS:
        if p_ is panel and k_ == key:
            return sec
    bs = float(panel["info"]["bs"])
    sec = OmegaGBTSection(_THIS, _THIS, panel, seg=bs / OMEGA_GBT_SEG, families=OMEGA_GBT_FAMILIES,
                          side=OMEGA_SIDE_BC, link_k=OMEGA_GBT_LINK_K,
                          nint=list(nint) if nint is not None else None)
    _OMEGA_GBT_SECTIONS.insert(0, (panel, key, sec))
    del _OMEGA_GBT_SECTIONS[4:]
    return sec


def omega_gbt_modes(sec):
    """GBT deformation modes of the cross-section in hierarchical order
    (section 25 / note, Sec. 8: B phi = lam (k^2 D + k^4 C) phi at k = 2 pi / bs,
    family by family).  Cached on the section.  Returns the dict of
    hierarchical_modes; hm['C'] = (lam, T) of the conventional modes."""
    if getattr(sec, "_hm", None) is None:
        bs = 0.5 * float(np.ptp(sec.X[:, 0]))
        sec._hm = hierarchical_modes(sec, kref=2.0 * np.pi / bs)
    return sec._hm


def _omega_inplane_split(sec, u):
    """Shares of global (G), distortional (D), local (L) and transverse-
    extension (T) deformation in an in-plane nodal displacement vector u (2N),
    with the split of the participation (Sec. 6.8 of the note)."""
    UT = sec._Torth @ (sec._Torth.T @ u) if "T" in sec.families else 0 * u
    UC = u - UT
    nat = sec.junctions if len(sec.junctions) >= 2 else sec.junctions + sec.tips
    J = np.array([2 * n + c for n in nat for c in (0, 1)])
    UG = sec._Gorth @ np.linalg.lstsq(sec._Gorth[J], UC[J], rcond=None)[0]
    UD = (sec._Dorth @ np.linalg.lstsq(sec._Dorth[J], (UC - UG)[J], rcond=None)[0]
          if sec._Dorth.shape[1] else 0 * UC)
    UL = UC - UG - UD
    n = {g: float(np.linalg.norm(v)) for g, v in (("G", UG), ("D", UD), ("L", UL), ("T", UT))}
    t = max(sum(n.values()), 1e-300)
    return {g: v / t for g, v in n.items()}


def omega_gbt_mode_table(sec, nshow=20):
    """[(number, lam, shares, dominant name)] of the first nshow conventional GBT modes."""
    hm = omega_gbt_modes(sec); lam, T = hm["C"]; out = []
    names = {"G": "global", "D": "distortional", "L": "local", "T": "transverse ext."}
    for k in range(min(nshow, T.shape[1])):
        sh = _omega_inplane_split(sec, sec.U @ T[:, k])
        (g1, s1), (g2, s2) = sorted(sh.items(), key=lambda t: -t[1])[:2]
        nm = names[g1] if (s1 >= 0.6 or s2 < 0.3) else f"{names[g1]}/{names[g2]}"
        out.append((k + 1, float(lam[k]), sh, nm))
    return out


def _omega_gbt_mode_participation(sec, lb, Cc, idx=None):
    """Participation of the conventional GBT modes in a buckling mode (as GBT
    papers: p_k = int |V_k(z)| dz with every mode scaled to unit maximum nodal
    displacement, normalised to 1).  Cc = coefficients of the conventional
    block (n, nb) in the hierarchical basis (columns idx of the full list)."""
    hm = omega_gbt_modes(sec); T = hm["C"][1]
    idx = np.arange(T.shape[1]) if idx is None else np.asarray(idx)
    sc = np.array([np.abs(sec.U @ T[:, k]).max() for k in idx])
    A = np.asarray(lb.B[0] @ Cc.T)                                   # (nx, n)
    I = (lb.wx @ np.abs(A)) * sc
    I = I / max(I.sum(), 1e-300)
    o = np.argsort(-I)
    return [(int(idx[k]) + 1, float(I[k])) for k in o]


def omega_gbt_buckling(panel, L, bc="SS", nmodes=3, return_model=False, hx=None):
    """Eigenvalues [(N_cr_i [N/mm], m_i)] of the omega panel by GBT -- the same
    call and return values as omega_spline_eigs (section 22).  With
    OMEGA_GBT_MODES set, only those conventional GBT modes (each with the
    warping it activates) are used in the buckling analysis; the pre-buckling
    state is then computed with all modes."""
    bc = bc.upper()
    if bc == "CC" and OMEGA_CC_END != "edge-load":
        raise ValueError("omega GBT: CC ends with rigid end plates are not implemented "
                         "(use the strip solver for --omega-cc-plate)")
    sec = omega_gbt_section(panel)
    L = float(L)
    if hx is None:
        hx = min(15.0, max(2.0, np.ptp(sec.X[:, 0]) / 40.0), L / 10.0)
    if OMEGA_GBT_MODES:
        hm = omega_gbt_modes(sec)
        idx = [k - 1 for k in OMEGA_GBT_MODES if 1 <= k <= hm["C"][1].shape[1]]
        pre = omega_gbt_solve(sec, L, bc, hx=hx, pre_only=True)["pre"]
        work = truncated_section(sec, None, 0, 0, 0, hm=hm, coupled_shear=True, cidx=idx)
        r = omega_gbt_solve(work, L, bc, hx=hx, nev=max(2, nmodes), return_mode=True, pre=pre)
        mpart = _omega_gbt_mode_participation(sec, r["lb"], r["coef"][:work.nip], idx)
    else:
        idx = None; work = sec
        r = omega_gbt_solve(sec, L, bc, hx=hx, nev=max(2, nmodes), return_mode=True)
        T = omega_gbt_modes(sec)["C"][1]
        Cc = np.linalg.lstsq(T, r["coef"][:sec.nip], rcond=None)[0]   # hierarchical coordinates
        mpart = _omega_gbt_mode_participation(sec, r["lb"], Cc)
    view = _OmegaGBTView(work, r["lb"], L, r["ndof"])
    vecs = r["coefs"][:max(1, nmodes)]
    modes = []
    for lam_, C in zip(r["eigs"], vecs):
        m, _r, _k = _omega_spline_modeinfo(view, C)
        modes.append((float(lam_), m))
    m, rr, kind = _omega_spline_modeinfo(view, vecs[0])
    info = dict(ndof=r["ndof"], stringer_share=rr, mode_type=kind,
                ds=float(panel["info"]["bs"]) / OMEGA_GBT_SEG, hx=hx, nm=work.nm,
                participation=r["participation"], mode_participation=mpart,
                n_conv=(len(idx) if idx is not None else sec.nip), n_conv_all=sec.nip)
    if return_model:
        info.update(model=view, vecs=vecs, panel=panel, section=sec,
                    strips=[(s_["a"], s_["b"], sec.wsrc[s_["w"]], sec.As[s_["w"]]) for s_ in sec.S])
    return modes, info


def _omega_draw_gbt_modes(sec, numbers, fname, title):
    """In-plane shapes of the conventional GBT modes `numbers` (1-based)."""
    import matplotlib.pyplot as plt
    hm = omega_gbt_modes(sec); lam, T = hm["C"]
    numbers = [k for k in numbers if 1 <= k <= T.shape[1]][:16]
    if not numbers:
        return
    tab = {k: (l_, sh, nm) for (k, l_, sh, nm) in omega_gbt_mode_table(sec, max(numbers))}
    ncol = 2; nrow = int(np.ceil(len(numbers) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(12, 1.9 * nrow + 0.6))
    axs = np.atleast_1d(axs).ravel(); X = sec.X; H = _hermite
    span = float(np.ptp(X[:, 0]))
    for a_, k in zip(axs, numbers):
        u = sec.U @ T[:, k - 1]; th = sec.Th @ T[:, k - 1]
        sc = 0.06 * span / max(np.abs(u).max(), 1e-30)
        for s_ in sec.S:
            A_, B_, Ls = s_["a"], s_["b"], s_["L"]; t_, n_ = s_["t"], s_["n"]
            va, vb = t_ @ u[2 * A_:2 * A_ + 2], t_ @ u[2 * B_:2 * B_ + 2]
            wa, wb = n_ @ u[2 * A_:2 * A_ + 2], n_ @ u[2 * B_:2 * B_ + 2]
            xi = np.linspace(0.0, 1.0, 7); P_ = []
            for x_ in xi:
                h_, _d, _dd = H(x_, Ls)
                w_ = h_ @ np.array([wa, th[A_], wb, th[B_]]); v_ = (1 - x_) * va + x_ * vb
                P_.append(X[A_] + x_ * Ls * t_ + sc * (v_ * t_ + w_ * n_))
            P_ = np.array(P_)
            a_.plot(X[[A_, B_], 0], X[[A_, B_], 1], color="0.82", lw=0.7)
            a_.plot(P_[:, 0], P_[:, 1], color="#1f5f8b", lw=1.0)
        l_, sh, nm = tab[k]
        a_.set_title(f"GBT mode {k}  ({nm}),  lambda = {l_:.3g}", fontsize=9)
        a_.set_aspect("equal"); a_.axis("off")
    for a_ in axs[len(numbers):]:
        a_.axis("off")
    fig.suptitle(title, fontsize=12, fontweight="bold")
    plt.tight_layout(); plt.savefig(fname, dpi=170); plt.close()


def _parse_mode_list(txt, nmax):
    """'1-4,7, 12' -> [1, 2, 3, 4, 7, 12] (within 1..nmax)."""
    out = []
    for part in str(txt).replace(";", ",").replace(" ", ",").split(","):
        if not part:
            continue
        if "-" in part:
            a_, b_ = part.split("-", 1); out += list(range(int(a_), int(b_) + 1))
        else:
            out.append(int(part))
    return sorted(set(k for k in out if 1 <= k <= nmax))


def _omega_gbt_questions(bs, material, pskin, pstiff, b3, h, th):
    """Classical GBT questions of the omega panel: cross-section discretisation,
    mode families, cross-section analysis (mode list + figure), and the modes
    used in the buckling analysis.  <Enter> everywhere = the validated model."""
    import os
    spec_ = dict(material=material, pskin=pskin, pstiff=pstiff, b3=b3, h=h, theta=th)
    pan = make_omega_panel(bs, None, [], [], spec=spec_)
    dn = omega_gbt_default_nint(pan); nd_ = pan["nodes"]
    print("\nGBT CROSS-SECTION DISCRETISATION  (natural nodes at the wall ends / fold lines;")
    print("intermediate nodes inside every wall)")
    print("   wall  part      width [mm]   intermediate nodes (default: segments <= bs/%g)" % OMEGA_GBT_SEG)
    for k_, w_ in enumerate(pan["walls"]):
        print(f"   {k_:>4}  {w_['part']:<7}{float(np.hypot(*(nd_[w_['j']] - nd_[w_['i']]))):>11.2f}   {dn[k_]:>6}")
    print("   + one link wall (offset tie, height c = (ts + tf)/2) at every web root")
    nint = ask_list("  n_int per wall", dn, int)
    nint = (list(nint) + [nint[-1]] * len(dn))[:len(dn)]
    globals()["OMEGA_GBT_NINT"] = None if list(nint) == list(dn) else nint
    print("\nGBT deformation-mode families:")
    print("  CS = conventional modes (global, distortional, local) + shear modes")
    print("       (Silvestre & Camotim 2013);   add:")
    print("  R  = nodal-rotation modes,  M = transverse-shear (Mindlin) modes,")
    print("  T  = transverse-extension modes (transverse membrane strain).")
    print("  Without T the transverse membrane force is zero (N_ss = 0), so the transverse")
    print("  pre-stress from the side restraint is lost: for this panel N_cr comes out")
    print("  25-50 % too high.  CSTMR = the validated model (0.58 % vs Abaqus).")
    fam = ask("  mode families", "CSTMR", str,
              ["cs", "csr", "csm", "csmr", "cst", "cstr", "cstm", "cstmr"]).upper()
    globals()["OMEGA_GBT_FAMILIES"] = fam
    print("\n  GBT cross-section analysis ...", flush=True)
    sec = omega_gbt_section(pan)
    tab = omega_gbt_mode_table(sec, 20)
    cnt = {g: sec.group.count(g) for g in ("G", "D", "L", "T", "A", "S", "M")}
    nat = len(set([w_[0] for w_ in sec.walls] + [w_[1] for w_ in sec.walls]))
    print("\n" + "=" * 74)
    print("  GBT CROSS-SECTION ANALYSIS")
    print("=" * 74)
    print(f"  nodes: {nat} natural + {sec.N - nat} intermediate = {sec.N};  side edges: {OMEGA_SIDE_BC}")
    print(f"  deformation modes ({fam}): {sec.nm} in total")
    print(f"     conventional (in-plane, with their warping)  {sec.nip:>4}   "
          f"[global {cnt['G']}, distortional {cnt['D']}, local {cnt['L']}"
          + (f", transverse extension {cnt['T']}" if cnt['T'] else "") + "]")
    print(f"     axial 1,  shear {cnt['S']},  transverse shear (Mindlin) {cnt['M']}")
    print("  conventional modes in GBT order (cross-section eigenproblem at the")
    print("  half-wave length bs/2; small lambda = long-wave / global modes):")
    print(f"   {'mode':>4}  {'type':<16}{'lambda':>11}   {'G %':>5}{'D %':>6}{'L %':>6}{'T %':>6}")
    for (k, l_, sh, nm) in tab:
        print(f"   {k:>4}  {nm:<16}{l_:>11.3e}   {100*sh['G']:>5.0f}{100*sh['D']:>6.0f}"
              f"{100*sh['L']:>6.0f}{100*sh['T']:>6.0f}")
    os.makedirs("out", exist_ok=True)
    try:
        _omega_draw_gbt_modes(sec, list(range(1, 13)), "out/fig_gbt_deformation_modes_OMEGA.png",
                              "GBT deformation modes of the omega panel (conventional, first 12)")
        print("  (shapes: out/fig_gbt_deformation_modes_OMEGA.png)")
    except Exception as e:
        print(f"  (mode figure skipped: {e})")
    print("\nGBT modes in the buckling analysis:")
    print(f"  1 = all {sec.nm} modes of the families above   (default; the validated model)")
    print("  2 = the first n conventional modes of the list, each with the shear (warping)")
    print("      field it activates")
    print("  3 = chosen conventional modes, e.g. 1-4,7,12")
    print("  (2 and 3 = reduced GBT, not validated against Abaqus; the pre-buckling state")
    print("   is then computed with all modes)")
    mch = ask("  choose 1, 2 or 3", 1, int, ["1", "2", "3"])
    modes = None
    if mch == 2:
        n_ = ask("  number of conventional modes n", 20, int)
        modes = list(range(1, max(1, min(n_, sec.nip)) + 1))
    elif mch == 3:
        modes = _parse_mode_list(ask("  GBT modes (list and ranges)", "1-20", str), sec.nip) or None
    if modes:
        print(f"  -> {len(modes)} conventional GBT modes: {modes if len(modes) <= 30 else str(modes[:30]) + ' ...'}")
    globals()["OMEGA_GBT_MODES"] = modes
    return pan, modes


def _answers_install(path):
    """--answers=FILE : answer the interactive questions from a text file.
    One answer per line, in the order of the questions; an empty answer
    (blank line, or a line with only a comment after spaces) = <Enter> =
    the default.  Lines starting with '#' in the first column are comments
    and are skipped; text after ' #' on an answer line is a comment."""
    import builtins
    ans = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            ans.append(ln.split("#", 1)[0].strip())
    it = iter(ans)

    def _inp(prompt=""):
        try:
            a = next(it)
        except StopIteration:
            a = ""
        print(f"{prompt}{a}")
        return a
    builtins.input = _inp


if __name__ == "__main__":
    import importlib.util as _ilu                         # fast sparse solver available?
    if _ilu.find_spec("pypardiso") is None:
        print("\n  NOTE: pypardiso (Intel MKL PARDISO) is not installed -- the code falls back\n"
              "  to SciPy SuperLU: same results, but the solves are about 2-3x slower.\n"
              "  Install it once with:   python -m pip install pypardiso\n")
    for _a in sys.argv:                                   # options (sections 24, 25)
        if _a in ("--solver=strip", "--solver=gbt"):
            GOVERNING_SOLVER = OMEGA_SOLVER = _a.split("=", 1)[1]
        elif _a.startswith("--omega-solver="):
            OMEGA_SOLVER = _a.split("=", 1)[1].lower()
        elif _a.startswith("--gbt="):
            GBT_FAMILIES = _a.split("=", 1)[1].upper()
        elif _a.startswith("--gbt-seg="):
            GBT_SEG = float(_a.split("=", 1)[1])
        elif _a.startswith("--gbt-hx="):
            GBT_HX = float(_a.split("=", 1)[1])
        elif _a == "--omega-gbt-fine":
            OMEGA_GBT_SEG = 40.0
        elif _a.startswith("--answers="):
            _answers_install(_a.split("=", 1)[1])
    if "--omega-mesh-fine" in sys.argv:
        OMEGA_MESH = "fine"
    elif "--omega-mesh-fast" in sys.argv:
        OMEGA_MESH = "fast"
    if "--omega-periodic" in sys.argv:
        OMEGA_SIDE_BC = "periodic"
    if "--mesh-fast" in sys.argv:
        MESH_STRIP, MESH_SPAN = 2.5, 10.0
    elif "--mesh-fine" in sys.argv:
        MESH_STRIP, MESH_SPAN = 0.625, -80.0
    if "--parametric-sections" in sys.argv:               # C, Z, L, T, I, RHS
        for _f in (run_parametric_C, run_parametric_Z, run_parametric_L, run_parametric_T,
                   run_parametric_I, run_parametric_RHS):
            _f()
    elif "--parametric-omega" in sys.argv:                # GBT (default) or --omega-solver=strip
        run_parametric_omega()
    elif "--parametric-c" in sys.argv:
        run_parametric_C()
    elif "--parametric-z" in sys.argv:
        run_parametric_Z()
    elif "--parametric-l" in sys.argv:
        run_parametric_L()
    elif "--parametric-t" in sys.argv:
        run_parametric_T()
    elif "--parametric-rhs" in sys.argv:
        run_parametric_RHS()
    elif "--parametric-i" in sys.argv or "--parametric" in sys.argv:
        run_parametric_I()
    else:
        main()
