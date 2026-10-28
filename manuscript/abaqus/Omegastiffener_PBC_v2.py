# -*- coding: utf-8 -*-
# =============================================================================
#  Omegastiffener_PBC_v2.py  --  periodic omega-stringer cells, SS or CC ends,
#  any list of cells, buckling mode 1 exported for the GBT-Abaqus figures.
#
#  Geometry, materials, mesh (S4R/S3, seed 3.7 mm), *TIE and the periodic
#  *EQUATIONs are those of Omegastiffener_PBC.py (your validated SS model).
#
#  End conditions (environment variable PANEL_BC):
#     SS (default) : every end edge U2 = 0, 1 N/mm ShellEdgeLoad on every end
#                    edge (identical to Omegastiffener_PBC.py)
#     CC           : z = 0  every end node U1=U2=U3=UR1=UR2=UR3 = 0
#                    z = L  every end node U1=U2=UR1=UR2=UR3 = 0, U3 free
#                    the same 1 N/mm ShellEdgeLoad on every end edge
#                    (clamped edges as in GBT; NOT rigid end plates)
#
#  Side edges (environment variable PANEL_SIDES):
#     PBC (default) : periodic, the width of the cell is constant (as in GBT)
#     PBC_FREE      : periodic with free lateral expansion and in-plane shear of
#                     the cell (common jump of U1 and U3 between the two side
#                     edges, dummy node PBC_DUMMY); SS ends only
#
#  Cells (environment variable PANEL_CASES, optional):
#     "bs:beta:layup,bs:beta:layup,..."   any pitch / aspect ratio / lay-up,
#     e.g. "150:12:Cross-ply,150:8:Isotropic"
#     not set -> the 30 standard cells (bs 150/285.5, beta 0.75/1.5/3, 5 lay-ups)
#
#  Output (folder PANEL_OUT, default ~/Downloads):
#     omega_validation_eigs_<BC>_PBC.csv     eigenvalues   (SS: ..._PBC.csv)
#     omega_validation_times_<BC>_PBC.csv    times [s]
#     fe_modes/MODE1_<cell>.csv              mode 1, all nodes: x,y,z,U1,U2,U3
#
#  RUN (PowerShell, in the folder of this file):
#     $env:PANEL_BC="CC"; abaqus cae noGUI=Omegastiffener_PBC_v2.py
#     $env:PANEL_BC="SS"; $env:PANEL_CASES="150:12:Cross-ply"; abaqus cae noGUI=Omegastiffener_PBC_v2.py
#  (Remove-Item Env:PANEL_CASES  resets to the 30 standard cells)
# =============================================================================

import os
import sys
import math

# -----------------------------------------------------------------------------
#  VALIDATION PARAMETRIC SWEEP DEFINITION (30 CASES)
# -----------------------------------------------------------------------------
BS_VALS   = [150.0, 285.5]
BETA_VALS = [0.75, 1.5, 3.0]

# Layup configs: (Name, Mat_Type, Skin_Angles(8 plies), Stiffener_Angles(7 plies))
LAYUPS = [
    ("Isotropic", "iso",  [], []),
    ("Cross-ply", "cfrp", [0, 90, 0, 90, 90, 0, 90, 0],       [0, 90, 0, 90, 0, 90, 0]),
    ("Sym1",      "cfrp", [45, -45, -45, 45, 45, -45, -45, 45], [45, -45, 45, -45, 45, -45, 45]),
    ("Sym2",      "cfrp", [60, -30, -30, 60, 60, -30, -30, 60], [60, -30, 60, -30, 60, -30, 60]),
    ("Unsym",     "cfrp", [0, 30, 60, 90, 0, 30, 60, 90],     [0, 30, 60, 90, 0, 30, 60])
]

# Generate the full list of test points
TEST_POINTS = []
for bs in BS_VALS:
    for beta in BETA_VALS:
        for lname, mtype, skin_seq, stiff_seq in LAYUPS:
            TEST_POINTS.append((bs, beta, lname, mtype, skin_seq, stiff_seq))

BC_LIST = [b.strip().upper() for b in os.environ.get("PANEL_BC", "SS").split(",") if b.strip()]
SIDES_LIST = [b.strip().upper() for b in os.environ.get("PANEL_SIDES", "PBC").split(",") if b.strip()]   # PBC | PBC_FREE (SS only)
_cases = os.environ.get("PANEL_CASES", "").strip()        # e.g. "150:12:Cross-ply,150:8:Isotropic"
if _cases:
    _lay = dict((l[0], l) for l in LAYUPS)
    TEST_POINTS = []
    for c in _cases.split(","):
        if not c.strip():
            continue
        bs_, beta_, name_ = [x.strip() for x in c.split(":")]
        if name_ not in _lay:
            raise ValueError("unknown lay-up %s (use %s)" % (name_, ", ".join(_lay)))
        l_ = _lay[name_]
        TEST_POINTS.append((float(bs_), float(beta_), l_[0], l_[1], l_[2], l_[3]))

OUTPUT_DIR = os.environ.get(
    "PANEL_OUT", os.path.join(os.path.expanduser("~"), "Downloads"))

NUM_MODES = 8          # eigenvalues requested (for ordering);  N_cr = mode 1
CSV_MODES = 6          # how many (eig,m) columns to write per point
MESH_SIZE = 3.7        # FIXED accurate mesh seed [mm]

# =============================================================================
try:
    from abaqus import *
    from abaqusConstants import *
    from caeModules import *
    import mesh
    import section
    ABAQUS_OK = True
except ImportError:
    ABAQUS_OK = False

# =============================================================================
def apply_periodic(m, a, inst, W, L, bc, sides, tol=1e-3):
    """Periodic side edges: tie every node of the skin edge x = -W to the node of
    the edge x = 0 at the same z, for all 6 DOFs (*EQUATION, the x = -W node is
    the eliminated term).  Requires the same seeding on both edges (true here:
    both are straight edges of length L, seeded with MESH_SIZE)."""
    a.regenerate()
    right = sorted([n for n in inst.nodes if abs(n.coordinates[0]) < tol],
                   key=lambda n: n.coordinates[2])
    left = sorted([n for n in inst.nodes if abs(n.coordinates[0] + W) < tol],
                  key=lambda n: n.coordinates[2])
    if len(right) != len(left):
        raise RuntimeError("PBC: %d / %d nodes on the side edges" % (len(right), len(left)))
    free = (sides == "PBC_FREE") and (bc == "SS")
    if free:
        rp = a.ReferencePoint(point=(-W / 2.0, -50.0, 0.0))           # dummy node, no element
        a.Set(referencePoints=(a.referencePoints[rp.id], ), name="PBC_DUMMY")
    n_eq = 0
    for k, (nl, nr) in enumerate(zip(left, right)):
        z = nr.coordinates[2]
        if abs(nl.coordinates[2] - z) > tol:
            raise RuntimeError("PBC: side-edge nodes do not match at z = %g" % z)
        a.Set(nodes=inst.nodes.sequenceFromLabels((nl.label, )), name="PBC_L%d" % k)
        a.Set(nodes=inst.nodes.sequenceFromLabels((nr.label, )), name="PBC_R%d" % k)
        at_end = (z < tol) or (abs(z - L) < tol)
        for dof in range(1, 7):
            if bc == "CC" and z < tol:
                continue      # CC, z = 0: all six DOFs prescribed on both corner nodes
            if bc == "CC" and abs(z - L) < tol and dof != 3:
                continue      # CC, z = L: only U3 is free on the corner nodes
            if bc == "SS" and at_end and dof == 2:
                continue      # SS: U2 = 0 is already prescribed on both corner nodes
            terms = [(1.0, "PBC_L%d" % k, dof), (-1.0, "PBC_R%d" % k, dof)]
            if free and dof in (1, 3):
                terms.append((1.0, "PBC_DUMMY", dof))                 # u_L = u_R - jump
            m.Equation(name="PBC_%d_%d" % (k, dof), terms=tuple(terms))
            n_eq += 1
    if bc == "SS":
        # rigid-body anchors on the x = 0 corner node (independent node of the equations):
        # U1 (lateral translation; the symmetry BC did this before) and U3 (axial)
        m.DisplacementBC(name="PBC_Anchor0", createStepName="Initial",
                         region=a.sets["PBC_R0"], u1=SET, u3=SET)
        if free:          # the free in-plane shear jump also frees the in-plane rotation
            m.DisplacementBC(name="PBC_AnchorL", createStepName="Initial",
                             region=a.sets["PBC_R%d" % (len(right) - 1)], u1=SET)
    print("  periodic side edges: %d node pairs, %d equations%s"
          % (len(right), n_eq, " + free macro strain" if free else ""))


def run_case(bs_val, beta_val, layup_name, mat_type, seq_skin, seq_stiff, out_dir, bc="SS",
             sides="SYM"):
    import time as _tm; _t_start = _tm.time()          # [timing] model build starts
    tag      = "bs%g_beta%g_%s" % (bs_val, beta_val, layup_name)
    if bc != "SS":
        tag += "_" + bc
    if sides != "SYM":
        tag += "_" + sides
    job_name = "val_" + tag.replace(".", "p").replace("-", "m")
    model_id = "Model_" + job_name
    run_dir  = os.path.join(out_dir, "omega_modes_runs", job_name)
    prof_dir = os.path.join(out_dir, "fe_profiles")
    for d in (run_dir, prof_dir):
        if not os.path.isdir(d):
            os.makedirs(d)
    os.chdir(run_dir)

    # ----- geometry ---------------------------------------------------------
    bs_input        = bs_val
    panel_len_input = beta_val * bs_input
    b3_fixed  = 27.75
    stiff_h   = 30.0
    theta_d   = 25.0
    panel_width = 2.0 * bs_input
    dist_origin = bs_input / 2.0
    b_val       = bs_input - (2 * b3_fixed)
    b1_calc     = 0.1 * b_val
    c1 = -dist_origin
    c2 = -(dist_origin + bs_input)
    centers = [c1, c2]
    theta_rad = math.radians(theta_d)

    if model_id in mdb.models: del mdb.models[model_id]
    if job_name in mdb.jobs:   del mdb.jobs[job_name]
    Mdb()
    mdb.models.changeKey(fromName="Model-1", toName=model_id)
    m = mdb.models[model_id]

    # ----- Material & Sections Assignment -----------------------------------
    if mat_type == "iso":
        mat = m.Material(name="ISO_MAT")
        mat.Elastic(table=((200000.0, 0.3), ))
        L_skin = [section.SectionLayer(material="ISO_MAT", thickness=1.0, numIntPts=3, plyName="S1")]
        L_stiff = [section.SectionLayer(material="ISO_MAT", thickness=0.875, numIntPts=3, plyName="St1")]
        t_skin = 1.0
        t_stiff = 0.875
    else:
        mat = m.Material(name="CFRP")
        mat.Elastic(type=LAMINA, table=((150000.0, 9000.0, 0.33, 5300.0, 5300.0, 5300.0), ))
        base_ply = 0.125
        L_skin = [section.SectionLayer(material="CFRP", thickness=base_ply, orientAngle=a, numIntPts=3, plyName="S%d"%k) for k, a in enumerate(seq_skin)]
        L_stiff = [section.SectionLayer(material="CFRP", thickness=base_ply, orientAngle=a, numIntPts=3, plyName="St%d"%k) for k, a in enumerate(seq_stiff)]
        t_skin = len(seq_skin) * base_ply
        t_stiff = len(seq_stiff) * base_ply

    m.CompositeShellSection(name="Sect_Skin", preIntegrate=ON, layup=L_skin, symmetric=False, thicknessType=UNIFORM)
    m.CompositeShellSection(name="Sect_Stiff", preIntegrate=ON, layup=L_stiff, symmetric=False, thicknessType=UNIFORM)
    clearance = (t_skin / 2.0) + (t_stiff / 2.0)

    # ----- skin part --------------------------------------------------------
    s_sk = m.ConstrainedSketch(name="Sk_Skin", sheetSize=2000.0)
    s_sk.Line(point1=(0.0, 0.0), point2=(-panel_width, 0.0))
    p_skin = m.Part(name="Skin", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_skin.BaseShellExtrude(sketch=s_sk, depth=panel_len_input)
    for cx in centers:
        for x_loc in [cx + b3_fixed + b1_calc, cx - b3_fixed - b1_calc]:
            d = p_skin.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=x_loc)
            try: p_skin.PartitionFaceByDatumPlane(datumPlane=p_skin.datums[d.id], faces=p_skin.faces)
            except Exception: pass
    reg_skin = p_skin.Set(faces=p_skin.faces, name="Set_Skin")
    p_skin.SectionAssignment(region=reg_skin, sectionName="Sect_Skin", thicknessAssignment=FROM_SECTION, offsetType=MIDDLE_SURFACE)
    p_skin.MaterialOrientation(region=reg_skin, orientationType=GLOBAL, axis=AXIS_2, additionalRotationType=ROTATION_NONE, localCsys=None, stackDirection=STACK_3)

    # ----- stringer part ----------------------------------------------------
    s_st = m.ConstrainedSketch(name="Sk_Stiff", sheetSize=2000.0)
    def draw_omega(sketch, cx):
        ht = b3_fixed - (stiff_h * math.tan(theta_rad))
        pts = [(cx + b3_fixed + b1_calc, 0), (cx + b3_fixed, 0),
               (cx + ht, stiff_h), (cx - ht, stiff_h),
               (cx - b3_fixed, 0), (cx - b3_fixed - b1_calc, 0)]
        for k in range(len(pts) - 1):
            sketch.Line(point1=pts[k], point2=pts[k + 1])
    for cx in centers:
        draw_omega(s_st, cx)
    p_stiff = m.Part(name="Stringers", dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p_stiff.BaseShellExtrude(sketch=s_st, depth=panel_len_input)
    reg_stiff = p_stiff.Set(faces=p_stiff.faces, name="Set_Stiff")
    p_stiff.SectionAssignment(region=reg_stiff, sectionName="Sect_Stiff", thicknessAssignment=FROM_SECTION, offsetType=MIDDLE_SURFACE)
    p_stiff.MaterialOrientation(region=reg_stiff, orientationType=GLOBAL, axis=AXIS_2, additionalRotationType=ROTATION_NONE, localCsys=None, stackDirection=STACK_3)

    # ----- assembly + tie ---------------------------------------------------
    a = m.rootAssembly
    i_skin  = a.Instance(name="Inst_Skin",  part=p_skin,  dependent=ON)
    i_stiff = a.Instance(name="Inst_Stiff", part=p_stiff, dependent=ON)
    a.translate(instanceList=("Inst_Stiff", ), vector=(0.0, clearance, 0.0))
    s_main = a.Surface(side1Faces=i_skin.faces,  name="Surf_Main")
    s_sec  = a.Surface(side1Faces=i_stiff.faces, name="Surf_Sec")
    m.Tie(name="Tie_Bond", main=s_main, secondary=s_sec, positionToleranceMethod=SPECIFIED, positionTolerance=1.5, adjust=ON, tieRotations=ON)

    # ----- buckling step ----------------------------------------------------
    m.BuckleStep(name="Step_Buckle", previous="Initial", numEigen=NUM_MODES, eigensolver=SUBSPACE, vectors=NUM_MODES * 2, maxIterations=1000)

    pad = 0.5
    all_e = i_skin.edges + i_stiff.edges

    # symmetry side edges
    e_sides = (i_skin.edges.getByBoundingBox(xMin=-pad, xMax=pad, yMin=-pad, yMax=pad, zMin=-pad, zMax=panel_len_input + pad) +
               i_skin.edges.getByBoundingBox(xMin=-panel_width - pad, xMax=-panel_width + pad, yMin=-pad, yMax=pad, zMin=-pad, zMax=panel_len_input + pad))
    if sides == "SYM":
        m.DisplacementBC(name="BC_Symmetry", createStepName="Step_Buckle", region=a.Set(edges=e_sides, name="Set_Sides"), u1=SET, ur2=SET, ur3=SET)

    # loaded ends: ALL end edges
    ymax = clearance + stiff_h + pad
    if bc == "SS":
        e_ends = (all_e.getByBoundingBox(xMin=-panel_width - pad, xMax=pad, yMin=-pad, yMax=ymax, zMin=-pad, zMax=pad) +
                  all_e.getByBoundingBox(xMin=-panel_width - pad, xMax=pad, yMin=-pad, yMax=ymax, zMin=panel_len_input - pad, zMax=panel_len_input + pad))
        m.DisplacementBC(name="BC_Ends_U2", createStepName="Initial", region=a.Set(edges=e_ends, name="Set_Ends"), u2=SET)

        # ORIGINAL LINE LOAD: 1.0 per unit length
        m.ShellEdgeLoad(name="Load_Comp", createStepName="Step_Buckle", region=a.Surface(side1Edges=e_ends, name="Surf_Ends"), magnitude=1.0, distributionType=UNIFORM, resultant=ON)
    else:
        # CC: clamped edges + the same uniform edge load as SS (see header).
        e_ends = (all_e.getByBoundingBox(xMin=-panel_width - pad, xMax=pad, yMin=-pad, yMax=ymax, zMin=-pad, zMax=pad) +
                  all_e.getByBoundingBox(xMin=-panel_width - pad, xMax=pad, yMin=-pad, yMax=ymax, zMin=panel_len_input - pad, zMax=panel_len_input + pad))
        m.ShellEdgeLoad(name="Load_Comp", createStepName="Step_Buckle", region=a.Surface(side1Edges=e_ends, name="Surf_Ends"), magnitude=1.0, distributionType=UNIFORM, resultant=ON)
        # the clamping BCs are node sets, built after meshing (see below)

    # ----- mesh -------------------------------------------------------------
    et1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
    et2 = mesh.ElemType(elemCode=S3,  elemLibrary=STANDARD)
    p_skin.setElementType(regions=(p_skin.faces, ),  elemTypes=(et1, et2))
    p_stiff.setElementType(regions=(p_stiff.faces, ), elemTypes=(et1, et2))
    for p in (p_skin, p_stiff):
        p.seedPart(size=MESH_SIZE)
        p.generateMesh()

    if bc == "CC":                                      # clamped end nodes (after meshing)
        a.regenerate()
        tol_z = 1e-3
        y_tie = clearance + 0.5                         # stringer nodes below: tie-secondary (flanges)
        for zc, nm in ((0.0, "0"), (panel_len_input, "L")):
            lab_s = [n.label for n in i_skin.nodes if abs(n.coordinates[2] - zc) < tol_z]
            lab_t = [n.label for n in i_stiff.nodes if abs(n.coordinates[2] - zc) < tol_z and n.coordinates[1] > y_tie]
            reg = a.SetFromNodeLabels(name="Set_Clamp" + nm,
                                      nodeLabels=(("Inst_Skin", tuple(lab_s)), ("Inst_Stiff", tuple(lab_t))))
            if nm == "0":
                m.DisplacementBC(name="BC_Clamp0", createStepName="Initial", region=reg,
                                 u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
            else:
                m.DisplacementBC(name="BC_ClampL", createStepName="Initial", region=reg,
                                 u1=SET, u2=SET, ur1=SET, ur2=SET, ur3=SET)
            print("  clamped end z = %g: %d skin + %d stringer nodes" % (zc, len(lab_s), len(lab_t)))

    if sides != "SYM":                                  # periodic side edges (after meshing)
        apply_periodic(m, a, i_skin, panel_width, panel_len_input, bc, sides)

    mdb.saveAs(pathName=os.path.join(run_dir, job_name + ".cae"))

    # ----- job --------------------------------------------------------------
    myJob = mdb.Job(name=job_name, model=model_id, type=ANALYSIS, numCpus=2, numDomains=2)
    try: os.remove(job_name + ".lck")
    except Exception: pass
    job_err = None
    try:
        _t_sub = _tm.time()                              # [timing] solver starts
        myJob.submit(); myJob.waitForCompletion()
        _t_done = _tm.time()                             # [timing] solver finished
    except Exception as e:
        job_err = str(e)

    if '_t_done' not in locals(): _t_done = _tm.time()      # [timing] job raised an error
    # ----- extract eigenvalues ----------------------------------------------
    import re
    import numpy as np
    from odbAccess import openOdb

    def eigval(desc):
        low = desc.lower()
        if "eigenvalue" in low:
            mm = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", desc[low.index("eigenvalue") + 10:])
            if mm: return float(mm.group(0))
        return None

    modes = []; crit = None; odb_note = ""
    abs_odb = os.path.join(run_dir, job_name + ".odb")
    
    if os.path.exists(abs_odb):
        try:
            odb  = openOdb(abs_odb, readOnly=True)
            step = odb.steps["Step_Buckle"]
            skin_key = None
            for key in odb.rootAssembly.instances.keys():
                if "SKIN" in key.upper():
                    skin_key = key; break
            inst = odb.rootAssembly.instances[skin_key]
            coord = {n.label: n.coordinates for n in inst.nodes}

            raw = []
            for fr in step.frames:
                ev = eigval(fr.description)
                if ev is None: continue
                try: usub = fr.fieldOutputs["U"].getSubset(region=inst)
                except Exception: continue
                xs = []; zs = []; u2 = []
                for v in usub.values:
                    c = coord.get(v.nodeLabel)
                    if c is None: continue
                    xs.append(c[0]); zs.append(c[2]); u2.append(v.data[1])
                if not u2:
                    raw.append((ev, -1, None)); continue
                xs = np.array(xs); zs = np.array(zs); u2 = np.array(u2)
                i0 = int(np.argmax(np.abs(u2)))
                x0, z0 = xs[i0], zs[i0]
                
                tolx = max(panel_width * 0.02, 1.0)
                mk = np.abs(xs - x0) < tolx
                zl = zs[mk]; ul = u2[mk]
                o = np.argsort(zl); zl = zl[o]; ul = ul[o]
                sg = np.sign(ul); sg[sg == 0] = 1.0
                m_est = int(np.sum(sg[:-1] * sg[1:] < 0)) + 1
                
                tolz = max(panel_len_input * 0.02, 1.0)
                mz = np.abs(zs - z0) < tolz
                xt = xs[mz]; ut = u2[mz]
                ox = np.argsort(xt); xt = xt[ox]; ut = ut[ox]
                raw.append((ev, m_est, (xt, ut, zl, ul)))

            # mode 1 (lowest positive eigenvalue): all nodes of skin and stringers
            try:
                fr1 = None
                for fr in step.frames:
                    ev = eigval(fr.description)
                    if ev is not None and ev > 1.0e-3 and (fr1 is None or ev < fr1[0]):
                        fr1 = (ev, fr)
                if fr1 is not None:
                    mdir = os.path.join(out_dir, "fe_modes")
                    if not os.path.isdir(mdir):
                        os.makedirs(mdir)
                    with open(os.path.join(mdir, "MODE1_%s.csv" % tag), "w") as fm:
                        fm.write("instance,node,x,y,z,U1,U2,U3,eigenvalue=%.6f\n" % fr1[0])
                        for ikey in odb.rootAssembly.instances.keys():
                            ins = odb.rootAssembly.instances[ikey]
                            cxyz = {n.label: n.coordinates for n in ins.nodes}
                            if not cxyz:
                                continue
                            for v in fr1[1].fieldOutputs["U"].getSubset(region=ins).values:
                                c = cxyz.get(v.nodeLabel)
                                if c is None:
                                    continue
                                fm.write("%s,%d,%.5f,%.5f,%.5f,%.8e,%.8e,%.8e\n" % (ikey, v.nodeLabel, c[0], c[1], c[2],
                                                                                   v.data[0], v.data[1], v.data[2]))
            except Exception as e:
                print("  mode export failed: %s" % str(e))
            odb.close()
            good = [r for r in raw if r[0] > 1.0e-3]
            good.sort(key=lambda r: r[0])
            modes = [(r[0], r[1]) for r in good]
            if good and good[0][2] is not None:
                xt, ut, zl, ul = good[0][2]
                if np.max(np.abs(ut)) > 0: ut = ut / np.max(np.abs(ut))
                if np.max(np.abs(ul)) > 0: ul = ul / np.max(np.abs(ul))
                crit = (xt, ut, zl, ul, good[0][1])
            odb_note = "modes=%d first=%s" % (len(modes), modes[:3])
        except Exception as e:
            odb_note = "ODB read failed: %s" % str(e)
    else:
        odb_note = ".odb missing"

    Ncr = modes[0][0] if modes else -1.0

    # ----- write results ----------------------------------------------------
    sfx = ("" if bc == "SS" else "_" + bc) + ("" if sides == "SYM" else "_" + sides)
    csv_path = os.path.join(out_dir, "omega_validation_eigs%s.csv" % sfx)
    new = not os.path.exists(csv_path)
    with open(csv_path, "a") as f:
        if new:
            hdr = "layup,bs,beta,L,panel_width,Ncr"
            for k in range(1, CSV_MODES + 1):
                hdr += ",m%d_eig,m%d_m" % (k, k)
            f.write(hdr + "\n")
        row = "%s,%g,%g,%g,%g,%s" % (layup_name, bs_val, beta_val,
              panel_len_input, panel_width,
              ("%.5f" % Ncr) if Ncr > 0 else "NA")
        for k in range(CSV_MODES):
            if k < len(modes): row += ",%.5f,%d" % (modes[k][0], modes[k][1])
            else: row += ",NA,NA"
        f.write(row + "\n")
    # [timing] separate file, the eigenvalue CSV above is unchanged
    _t_end = _tm.time()
    t_path = os.path.join(out_dir, "omega_validation_times%s.csv" % sfx)
    t_new = not os.path.exists(t_path)
    with open(t_path, "a") as ft:
        if t_new:
            ft.write("layup,bs,beta,L,t_model_s,t_solver_s,t_post_s,t_total_s,numCpus\n")
        ft.write("%s,%g,%g,%g,%.2f,%.2f,%.2f,%.2f,%d\n" % (layup_name, bs_val, beta_val,
                 panel_len_input, _t_sub - _t_start, _t_done - _t_sub, _t_end - _t_done,
                 _t_end - _t_start, 2))

    if crit is not None:
        xt, ut, zl, ul, mcr = crit
        with open(os.path.join(prof_dir, "FE_trans_%s.csv" % tag), "w") as f:
            f.write("x_mm,U2_norm (critical mode m=%d)\n" % mcr)
            for xx, uu in zip(xt, ut): f.write("%.4f,%.6f\n" % (xx, uu))
        with open(os.path.join(prof_dir, "FE_long_%s.csv" % tag), "w") as f:
            f.write("z_mm,U2_norm (critical mode m=%d)\n" % mcr)
            for zz, uu in zip(zl, ul): f.write("%.4f,%.6f\n" % (zz, uu))

    print("=" * 66)
    print(" %-12s bs=%g beta=%g -> Ncr=%s  modes(eig,m)=%s"
          % (layup_name, bs_val, beta_val,
             ("%.4f" % Ncr) if Ncr > 0 else "FAILED", modes[:4]))
    if job_err: print("  job note:", job_err)
    print("  " + odb_note)
    print("=" * 66)
    sys.stdout.flush()
    return Ncr

if __name__ == "__main__":
    if not ABAQUS_OK:
        sys.stderr.write("\nRun with the Abaqus kernel, e.g.:\n    abaqus cae noGUI=Omegastiffener_PBC_v2.py\n")
        sys.exit(1)
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    import time
    t0 = time.time()
    for sides in SIDES_LIST:
      for bc in BC_LIST:
        for k, (bs, beta, lname, mtype, skin_seq, stiff_seq) in enumerate(TEST_POINTS, 1):
            print(">>> [%s %s %d/%d] %s bs=%g beta=%g (elapsed %.0f s)"
                  % (bc, sides, k, len(TEST_POINTS), lname, bs, beta, time.time() - t0))
            sys.stdout.flush()
            try:
                run_case(bs, beta, lname, mtype, skin_seq, stiff_seq, OUTPUT_DIR, bc=bc,
                         sides=sides)
            except Exception as exc:
                sys.stderr.write("point %s %s bs=%g beta=%g %s failed: %s\n"
                                 % (bc, sides, bs, beta, lname, str(exc)))
            
    print(">>> DONE %d points in %.0f s -> %s" % (len(TEST_POINTS), time.time() - t0, os.path.join(OUTPUT_DIR, "omega_validation_eigs.csv")))