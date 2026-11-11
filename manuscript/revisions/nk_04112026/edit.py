import copy
import sys
from lxml import etree
sys.path.insert(0, '.')
from trk import *  # noqa

SRC, DST = sys.argv[1], sys.argv[2]
tree = etree.parse(SRC); root = tree.getroot(); body = root.find('w:body', ns)
cx = Ctx(root)
B = list(body)                      # element references BEFORE any insertion
P = {i: B[i] for i in range(len(B))}
assert vis_text(P[10]).startswith('1. Introduction') and vis_text(P[16]).startswith('2. Theory')

# ---------------------------------------------------------------- 1. Introduction (comments 7, 8, 10, 21)
INTRO = {
11: [
"Stringer-stiffened panels are the standard construction of aircraft fuselage and wing covers, and buckling is often "
"the criterion that governs their design [1,2]. In composite panels, buckling can occur locally in the skin between the "
"stringers or in the stringer walls, globally as a column of skin and stringers together, or as a coupling of these "
"forms. Since many geometries and lay-ups are compared in the early design stage, buckling models are needed that are "
"both accurate and fast.",
"Among the stringer profiles, the omega (hat) stringer is widely used. Its closed profile gives high bending and "
"torsional stiffness, and its two flanges provide a wide bonding area with the skin [3,4]. A cover usually carries many "
"identical stringers at a constant spacing, the stringer pitch. Away from the longitudinal edges of the cover, the "
"deformation therefore repeats from one stringer bay to the next, and the analysis can be reduced to one repeating cell "
"[5–7]. The two side edges of this cell are linked by periodic boundary conditions: the displacements and rotations "
"at one side edge are set equal to those at the other side edge. The cell then behaves like a part of an infinitely "
"wide panel. The same idea is used in the theory of wave propagation in periodic structures, where the response of one "
"repeating unit describes the response of the whole structure [8]. Prismatic plate assemblies such as stiffened panels "
"are also analyzed efficiently by semi-analytical strip methods. The exact strip programs VIPASA and VICONOPT assemble "
"exact stiffness matrices of the individual plate strips [9,10], and finite strip methods use polynomial functions "
"across each strip combined with harmonic [11] or spline [12] functions along the length."],
12: [
"Closed-form and semi-analytical models give fast buckling estimates for stiffened composite panels and are therefore "
"well suited to preliminary design. For blade-stringer-stiffened plates, Mittelstedt derived closed-form buckling loads "
"of a repeating unit with periodic boundary conditions [5] and explicit local buckling solutions that account for the "
"interaction between stiffener and plate [6]. For omega stringers, closed-form solutions were developed for the buckling "
"of compressively loaded plates braced by omega stringers [3] and for the local postbuckling of hat-stringer-stiffened "
"plates under transverse compression [13]. Schilling and Mittelstedt extended this work to the local buckling [4] and "
"the column buckling [7] of omega-stringer-stiffened panels, to shear buckling [14], to combined compression and shear "
"[15] and to local postbuckling [16]. Most recently, El Yaakoubi-Mesbah and Mittelstedt presented closed-form and "
"computational models for the local buckling and postbuckling of compressively loaded omega-stringer-stiffened panels "
"that include the effect of the stringer feet [17–19]. An overview of these methods is given in [20].",
"All these models reduce the panel to a repeating unit, either with periodic boundary conditions [5–7] or with "
"symmetry conditions [18]. Symmetry conditions are simpler to apply, but they admit only modes that are symmetric about "
"the symmetry lines, so antisymmetric modes are excluded. To obtain closed-form expressions, the models also need "
"several idealizations. The interaction between skin and stringer is represented by elastic restraints, and local and "
"global buckling are treated separately. The pre-buckling stress is assumed to be uniform. The membrane–bending "
"coupling of unsymmetric laminates is included only approximately, through a reduced bending stiffness whose validity "
"is limited [21]. The physics-informed neural network developed earlier by the authors for the local buckling of "
"stiffened panels uses the same reduced bending stiffness [22]. Shell finite element models avoid all these "
"idealizations. However, they are computationally expensive in parametric studies with many designs. They also give "
"the buckling mode only as a nodal displacement field, which does not show directly how local, distortional and global "
"deformations contribute to it."],
13: [
"Generalized Beam Theory (GBT) lies between these two groups of models. It was originated by Schardt [23,24] and "
"developed further by Davies and Leach [25]. GBT describes the deformation of a thin-walled member as a sum of "
"cross-section deformation modes, each multiplied by an amplitude function along the member axis. The analysis has two "
"stages. The cross-section analysis computes the deformation modes and their stiffness properties, and the member "
"analysis determines the amplitude functions and the critical load. Since every mode has a clear structural meaning, "
"such as global, distortional or local deformation, GBT shows directly how the modes interact in coupled instabilities, "
"and it needs far fewer degrees of freedom than a shell model [26–29]. Silvestre and Camotim extended first-order "
"and second-order GBT to arbitrary orthotropic materials [30,31] and later formulated it for thin-walled composite "
"members with shear deformation [32] and for their vibration [33]. The cross-section analysis was generalized to "
"arbitrarily branched open sections [34] and then reformulated as a systematic procedure that can be automated "
"[35,36]. This work also clarified the role of the shear deformation modes [37] and is implemented in the code GBTul "
"2.0 [38]. Closed sections were treated for circular cylindrical shells and tubes [39] and for regular polygonal tubes "
"[40]. For composite structures, GBT was applied to fiber-reinforced polymer (FRP) open-section columns [41] and to "
"laminated CFRP cylindrical stiffened panels [42], and Taig and Ranzi formulated it for stiffened sections [43].",
"In the member analysis, the amplitude functions can be approximated by the Ritz method, which is widely used for the "
"buckling of laminated plates and beams [44–46]. Kharghani and Mittelstedt combined GBT and the Ritz method into a "
"GBT–Ritz approach and applied it to the global, local and distortional buckling of shear-deformable composite "
"beams with open cross-sections [47], to composite box girders [48] and to thin-walled composite beams with symmetric "
"and unsymmetric lay-ups [49,50]. In this approach, the cross-section analysis provides the modes, and the Ritz method "
"provides the amplitude functions with only a few unknowns per mode."],
14: [
"However, GBT has so far been formulated for members whose cross-section is either open, with free longitudinal edges, "
"or closed. The cross-section of a periodic panel cell is neither: the skin ends at two side edges, and the periodic "
"condition requires equal displacements and rotations at both of them. This condition has not yet been formulated in "
"GBT. It also cannot simply be imposed on a truncated set of classical modes. These modes are computed for an open "
"cross-section with independent side edges, so only a few of their combinations satisfy the periodic condition. In "
"addition, the omega stringer forms closed cells with the skin, so Vlasov’s assumption of zero membrane shear "
"cannot hold; the stringer is offset from the skin; and the pre-buckling stress is not uniform over the cross-section."],
15: [
"The present study closes this gap. Its main novelty is that the periodic condition is formulated within the "
"cross-section analysis of GBT. As a result, every mode, and therefore every truncated set of modes, satisfies the "
"periodic condition exactly. This formulation is combined with a shear-deformable GBT–Ritz member analysis [50] "
"and applied to omega-stringer-stiffened panels with symmetric and unsymmetric laminates. The results are verified "
"against Abaqus shell models. The mode decomposition shows the transition from local to global buckling as the cell "
"becomes longer, and it reveals a long-wave mode that is caused by the constant width imposed by the periodic "
"condition. Since every mode is periodic by construction, the formulation is also suited as the physics model of "
"physics-informed neural networks for such panels, extending the authors’ earlier network [22] (Section 4)."],
}
for i, texts in INTRO.items():
    rewrite_paragraph(cx, P[i], texts)

# ---------------------------------------------------------------- Section 2, first paragraph: ".." left by an edit
full = vis_text(P[17]); k = full.find('critical load..'); assert k > 0
replace_span(cx, P[17], k + len('critical load.'), k + len('critical load..'), '')

# ---------------------------------------------------------------- comment 27: global axes X -> Y, Y -> Z
REN = {'X': 'Y', 'Y': 'Z'}
nren = 0
# global longitudinal axis x -> X (Aya, 10 Nov 12:48): only where x is the GLOBAL coordinate; the local wall
# coordinate x of GBT (Silvestre-Camotim: x, s, z) stays.  (paragraph, index of the visible run in the original)
XG = {20: [22], 57: [6, 14, 26], 58: [11, 23], 69: [92], 81: [13, 24], 146: [2], 148: [16, 31, 38], 152: [55]}
xruns = []
for i, ks in XG.items():
    rs = vis_runs(P[i])
    for k in ks:
        r, s_ = rs[k]
        assert r.tag == q('m:r') and s_.strip('\u2009') == 'x', (i, k, s_)
        xruns.append(r)
for i in (20, 38, 44, 57, 58, 60, 68, 69, 81, 130, 148):
    for r, s in list(vis_runs(P[i])):
        if r.tag == q('m:r') and s in REN:
            w = wrap_del(cx, r); ins = cx.mark('ins'); nr = copy.deepcopy(r)
            for t in nr.findall(q('m:t')):
                t.text = REN[s]
            ins.append(nr); w.addnext(ins); nren += 1
for r in xruns:
    w = wrap_del(cx, r); ins = cx.mark('ins'); nr = copy.deepcopy(r)
    for t in nr.findall(q('m:t')):
        t.text = t.text.replace('x', 'X')
    ins.append(nr); w.addnext(ins); nren += 1
print('axis symbols renamed:', nren)
full = vis_text(P[20]); k = full.find(',s,z ('); assert k > 0
replace_span(cx, P[20], k + 5, k + 6, '($x$ parallel to $X$, ')

# ---------------------------------------------------------------- Figure 1 caption (comment 28: sub-captions only in the caption)
cap = P[22]
replace(cx, cap, 'linked by the periodic condition.', 'linked by the periodic condition, and the global axes $X$, $Y$, $Z$.')
replace(cx, cap, 'Infinitely wide panel and the cell (shaded).',
        'Infinitely wide panel with identical stringers at the pitch $b_s$ and the cell (shaded), which contains two stringers.')
replace(cx, cap, 'GBT discretisation for ', 'GBT discretization of the cell for ')
replace(cx, cap, ' mm with nodes, wall types', ' mm (92 nodes, 19 walls and 93 segments) with natural and intermediate nodes, wall types')

# ---------------------------------------------------------------- comment 37: Eq. (9), gamma_k on the whole segment
p53 = P[53]
replace(cx, p53, 'For an open contour the minimum is zero and Vlasov’s warping is recovered; in a closed contour the remaining ',
        'Only the G and D modes receive this warping. Their segments are not elongated, so $v_p$ is constant along each '
        'segment and, with the warping linear in $s$, the bracket is the shear strain $γ_k$ along the whole segment, '
        'not only at its mid-point. For an open contour the minimum is therefore zero, $γ_k$ = 0 holds along the entire '
        'contour and Vlasov’s warping is recovered exactly; in a closed contour the remaining ')

# ---------------------------------------------------------------- comment 81: sentence before Eq. (16)
rewrite_paragraph(cx, P[88], ["Substituting the strains of Eq. (3) and the constitutive law of Eq. (5) into Eq. (15) and "
                              "integrating over the contour gives, in the notation of Silvestre and Camotim [32],"])

# ---------------------------------------------------------------- comment 66: Table 1, row U5 before row U6
tbl = P[60]; rows = tbl.findall('w:tr', ns)
r6, r5 = rows[-2], rows[-1]
assert vis_text(r6).startswith('U6') and vis_text(r5).startswith('U5'), (vis_text(r6)[:10], vis_text(r5)[:10])
new5 = copy.deepcopy(r5)
CMT = {q('w:commentRangeStart'), q('w:commentRangeEnd'), q('w:commentReference')}
for e in new5.iter():
    if q('w:id') in e.attrib and e.tag not in CMT:
        e.set(q('w:id'), str(cx.nid)); cx.nid += 1
    for a in ('{http://schemas.microsoft.com/office/word/2010/wordml}paraId',
              '{http://schemas.microsoft.com/office/word/2010/wordml}textId'):
        e.attrib.pop(a, None)
# old row: no comment anchors (they move with the new row), everything deleted
for e in list(r5.iter(*CMT)):
    par = e.getparent()
    if e.tag == q('w:commentReference'):
        par.getparent().remove(par)                 # the run holding the reference
    else:
        par.remove(e)


def row_mark(row, kind):
    trPr = row.find(q('w:trPr'))
    if trPr is None:
        trPr = etree.Element(q('w:trPr')); row.insert(0, trPr) if row[0].tag != q('w:tblPrEx') else row.insert(1, trPr)
    chg = trPr.find(q('w:trPrChange'))
    (chg.addprevious(cx.mark(kind)) if chg is not None else trPr.append(cx.mark(kind)))
    for p in row.iter(q('w:p')):
        para_mark(cx, p, kind)
    for cp in row.iter(q('m:ctrlPr')):
        mark_ctrl(cx, cp, kind)


for r_, s_ in vis_runs(r5):
    wrap_del(cx, r_)
row_mark(r5, 'del')
for r_, s_ in vis_runs(new5):
    w_ = cx.mark('ins'); par = r_.getparent(); par.replace(r_, w_); w_.append(r_)
row_mark(new5, 'ins')
r6.addprevious(new5)


# ================================================================ 11 Nov 2026: reviewer-proofing of the theory
# (Navid: "the theory is complicated"; Bloch wavenumbers; Aya: add, correct, make the idea easy to follow)
# ---- Section 2, opening: overview in plain words + Scheme 1
replace(cx, P[17], 'Classical GBT is extended, and recovered as a special case (Section 2.3.4), by periodic modes, the main '
        'contribution (Section 2.2), closed contours without Vlasov’s assumptions, Reissner–Mindlin laminated walls '
        'and a computed pre-buckling state.',
        'Scheme 1 summarizes the analysis. Compared with classical GBT, four things are different: the modes are periodic '
        '(Section 2.2, the main contribution); the contours are closed, so Vlasov’s assumptions are not used; the walls '
        'are shear-deformable (Reissner–Mindlin) laminates; and the pre-buckling state is computed instead of assumed. '
        'Classical GBT is recovered as a special case (Section 2.3.4). The idea of the periodic formulation is simple. The '
        'two side edges of the cell are the same line of the infinite panel, so the cell must deform identically at both. '
        'Shell models enforce this on the solution by constraint equations. Here, every cross-section mode is made periodic '
        'before the member analysis. Any combination of periodic modes is periodic, so the member analysis itself needs no '
        'change.')
figp, capp = P[82], P[83]
assert figp.find('.//' + q('w:drawing')) is not None and vis_text(capp).startswith('Figure 2.')
sch = copy.deepcopy(figp)
for e in list(sch.iter()):
    if e.tag in (q('w:commentRangeStart'), q('w:commentRangeEnd'), q('w:bookmarkStart'), q('w:bookmarkEnd'), q('w:lastRenderedPageBreak')):
        e.getparent().remove(e)
    for a in list(e.attrib):
        if a.endswith('}paraId') or a.endswith('}textId') or a.endswith('}anchorId') or a.endswith('}editId'):
            del e.attrib[a]
WP = 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'
A = 'http://schemas.openxmlformats.org/drawingml/2006/main'
PIC = 'http://schemas.openxmlformats.org/drawingml/2006/picture'
RNS = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
cxmax = max(int(e.get('cx')) for e in root.iter('{%s}extent' % WP))
SCH_PX = (3869, 2694)
cx_ = int(1.00 * cxmax); cy_ = int(round(cx_ * SCH_PX[1] / SCH_PX[0]))        # 0.90: fits on the page below the text
for e in sch.iter('{%s}extent' % WP, '{%s}ext' % A):
    e.set('cx', str(cx_)); e.set('cy', str(cy_))
ids = [int(e.get('id')) for e in root.iter('{%s}docPr' % WP, '{%s}cNvPr' % PIC)]
dp = sch.find('.//{%s}docPr' % WP); dp.set('id', str(max(ids) + 1)); dp.set('name', 'Scheme 1')
dp.set('descr', 'Scheme 1. Steps of the periodic GBT-Ritz analysis.')
cp_ = sch.find('.//{%s}cNvPr' % PIC); cp_.set('id', str(max(ids) + 2)); cp_.set('name', 'Scheme 1'); cp_.set('descr', 'fig/scheme01_overview.png')
sch.find('.//{%s}blip' % A).set('{%s}embed' % RNS, 'rIdScheme1')
for r_ in list(sch.findall(q('w:r'))):
    w_ = cx.mark('ins'); sch.replace(r_, w_); w_.append(r_)
para_mark(cx, sch, 'ins')
P[17].addnext(sch)
capn = new_para_after(cx, sch, capp, 'Scheme 1.')
rp_ = copy.deepcopy(first_rpr(capp)) if first_rpr(capp) is not None else etree.Element(q('w:rPr'))
for e in rp_.findall(q('w:b')) + rp_.findall(q('w:bCs')):
    rp_.remove(e)
ins_ = cx.mark('ins')
for r_ in markup_runs(' Flowchart of the periodic GBT\u2013Ritz analysis. The cross-section analysis (top) is done once per '
                      'cross-section; the frame contains the periodic modes, the contribution of this study; the member analysis '
                      '(bottom) is repeated for every length and end condition. The drawings are computed results of the Unsym cell '
                      '($b_s$ = 150 mm, $\u03b2$ = 1.5, SS; Figures 6 and 7): the cell, its three dominant periodic modes and its '
                      'critical mode.', rp_):
    ins_.append(r_)
capn.append(ins_)

# ---- 2.2.2: why the modes are built family by family (plain words before the four steps)
replace(cx, P[69], 'To keep the GBT meaning of the modes, they are built family by family. With ',
        'Any basis of this null space would give periodic modes, but such modes would no longer be global, distortional, '
        'local or transverse extension modes. To keep the GBT meaning of the modes, the periodic modes are therefore built '
        'family by family, in three ideas: a mode that is already periodic is kept; modes that are not periodic are combined '
        'within their own family so that their side-edge differences cancel; and the few periodic deformations that are '
        'then still missing are added. In detail, with ')

# ---- 2.2.3: the remaining shear strain is the shear flow of a closed section (checked in the code, see RESEARCH_LOG)
append_text(cx, P[73], ' Indeed, the stationarity condition of Eq. (9) at a node is the balance of the shear flows (shear '
            'stiffness times $γ_k$) of the segments that meet there, so the shear flow is constant between two '
            'junctions, as in closed thin-walled sections; for the cross-ply cell with $b_s$ = 150 mm the computed warping '
            'satisfies this balance to 10⁻⁹ of the largest shear flow.')

# ---- 2.2.5: which Bloch wavenumbers a cell contains (Navid, 11 Nov)
replace(cx, P[81], 'modes that repeat in every bay and modes that alternate.',
        'modes that repeat in every bay and modes that alternate. In general, a cell of $n$ pitches with real periodic modes '
        'contains only the wavenumbers $k_B$$b_s$ = 2π$m$/$n$ with integer $m$. In an infinite panel, $k_B$$b_s$ can take '
        'any value between 0 and π, and the critical load is the minimum over this range; the present cell of two '
        'pitches samples its two ends. Other wavenumbers require a wider cell or complex Bloch conditions and are examined '
        'for the long cells in Section 3.5.')

replace(cx, P[81], 'To order them from global to local, the conventional modes are recombined',
        'To order the conventional modes from global to local, they are recombined')

# ---- 2.3.2: laminates with membrane-bending coupling bend before buckling
append_text(cx, P[113], ' For laminates with membrane–bending coupling (coupling stiffness of Eq. (5) not zero), the cell '
            'bends from the beginning of loading ($κ_0$ ≠ 0), so that buckling is not a bifurcation in the strict '
            'sense. As in the linear buckling analyses of the shell models (Section 3), the critical load is here the '
            'bifurcation load of the linearized problem: the membrane forces of the first-order state enter the geometric '
            'tensors, and its deflections are neglected.')

# ---- 3.5: cells of three and four pitches (manuscript/run_npitch.py, plot_npitch.py; results/npitch*.json*)
# ---- 4: conclusions
replace(cx, P[194], 'for the two cells examined.',
        'for the two cells examined. Cells of three and four pitches, which contain further Bloch wavenumbers, give the same '
        'long-wave loads in GBT and, for four pitches, in Abaqus.')


# ---- Table 2 and Section 3.7: number of Mindlin modes corrected to the solver's section (11 Nov 2026).
# Solver: b_s = 150: M = 205, total 569 (as "nmodes": 569 in results/long_curves.jsonl); b_s = 285.5: M = 174, total 482.
# The formula of Section 2.2.4 gives these values (23 and 20 nodes with two independent tangents, incl. the link walls).
cells = list(P[79].iter(q('w:tc')))
for k_, old_, new_ in ((16, '194', '205'), (17, '558', '569'), (25, '166', '174'), (26, '474', '482')):
    assert vis_text(cells[k_]) == old_, (k_, vis_text(cells[k_]))
    replace(cx, cells[k_], old_, new_)
replace(cx, P[188], '474–558 modes', '482–569 modes')

# ================================================================ 11 Nov 2026 (evening): n-pitch study as Table 6 + Figure 10
# ---- renumbering (tracked): new Table 6 and Figure 10 in Section 3.5
assert 'finer segments lower the GBT result further, Table 6' in vis_text(P[161])
replace(cx, P[161], 'finer segments lower the GBT result further, Table 6', 'finer segments lower the GBT result further, Section 3.7')
for p_, old_, new_ in ((188, 'Table 7', 'Table 8'), (189, 'Table 7', 'Table 8'), (185, 'Table 6', 'Table 7'), (186, 'Table 6', 'Table 7'),
                       (181, 'Figure 10', 'Figure 11'), (183, 'Figure 10', 'Figure 11')):
    assert vis_text(P[p_]).count(old_) == 1, (p_, old_)
    replace(cx, P[p_], old_, new_)

# ---- Section 3.5: new paragraph after the Abaqus paragraph
assert vis_text(P[175]).startswith('Abaqus confirms both mechanisms')
np35 = new_para_after(cx, P[175], P[175],
    'The two-pitch cell contains only the Bloch wavenumbers $k_B$$b_s$ = 0 and π (Section 2.2.5). Whether another wavenumber '
    'gives a lower load for the long-wave mode was examined with cells of three and four pitches, which also contain '
    '$k_B$$b_s$ = 2π/3 and π/2, in GBT and, for four pitches, in Abaqus (the shell model of this section with the width '
    '4$b_s$). The wavenumber of every GBT mode was identified from its skin deflection shifted by one pitch. Table 6 and Figure 10 '
    'show that the modes with $k_B$$b_s$ = 2π/3 and π/2 buckle at higher loads than the alternating mode, that the '
    'four-pitch cells give the same critical loads as the two-pitch cells in both models, and that their second and third '
    'eigenvalues are equal in both models, as expected for the two Bloch waves with $k_B$$b_s$ = ±π/2. For these cells '
    'the load thus decreases steadily towards $k_B$$b_s$ = π, which the two-pitch cell contains; wavenumbers between the '
    'sampled values were not examined. GBT stays 0.3–1.1 % below Abaqus, as for the short cells.')


def caption_after(ref, like, label, text):
    cp = new_para_after(cx, ref, like, label)
    rp = copy.deepcopy(first_rpr(like)) if first_rpr(like) is not None else etree.Element(q('w:rPr'))
    for e in rp.findall(q('w:b')) + rp.findall(q('w:bCs')):
        rp.remove(e)
    ins = cx.mark('ins')
    for r in markup_runs(text, rp):
        ins.append(r)
    cp.append(ins)
    return cp


def clean_copy(el):
    el = copy.deepcopy(el)
    for e in list(el.iter()):
        if e.tag in (q('w:commentRangeStart'), q('w:commentRangeEnd'), q('w:bookmarkStart'), q('w:bookmarkEnd'),
                     q('w:lastRenderedPageBreak'), q('w:rPrChange'), q('w:pPrChange'), q('w:tcPrChange'), q('w:trPrChange')):
            e.getparent().remove(e)
            continue
        for a in list(e.attrib):
            if a.endswith('}paraId') or a.endswith('}textId') or a.endswith('}anchorId') or a.endswith('}editId'):
                del e.attrib[a]
    return el


# ---- Table 6 (after the caption of Figure 9)
assert vis_text(P[179]).startswith('Figure 9.') and vis_text(P[189]).startswith('Table 8.') and P[190].tag == q('w:tbl')
tcap = caption_after(P[179], P[189], 'Table 6.',
    ' Critical loads [N/mm] for the Bloch wavenumbers $k_B$$b_s$ contained in periodic cells of two, three and four pitches '
    '($b_s$ = 150 mm, $β$ = 8, SS; GBT with spans of 15 mm; Abaqus: shell model of the four-pitch cell, whose first '
    'eigenvalue equals that of the two-pitch cell) and deviation e of GBT from Abaqus.')
tbl = clean_copy(P[190])
rows_t = tbl.findall(q('w:tr')); hdr_t, body_t = rows_t[0], rows_t[1]
for r_ in rows_t:
    tbl.remove(r_)
grid = tbl.find(q('w:tblGrid'))
for g_ in list(grid):
    grid.remove(g_)
for w_ in (1330, 1000, 1330, 1400, 1600, 1199):
    gc = etree.SubElement(grid, q('w:gridCol')); gc.set(q('w:w'), str(w_))
import json as _json
NT = _json.load(open('/home/claude/paper2pinnaya/manuscript/results/npitch_table.json'))
KBTXT = {r'\pi': 'π', r'2\pi/3': '2π/3', r'\pi/2': 'π/2', '0': '0'}
data = [['Lay-up', '$k_B$$b_s$', 'Cells (pitches)', 'GBT', 'Abaqus', 'e [%]']]
for lay, fmt in (('Isotropic', '{:.2f}'), ('Sym1', '{:.3f}')):
    for j_, kb in enumerate((r'\pi', r'2\pi/3', r'\pi/2', '0')):
        t_ = NT[lay][kb]; ab = t_['N_abq']
        data.append([lay if j_ == 0 else '', KBTXT[kb], ', '.join(str(c) for c in t_['cells']), fmt.format(t_['N_gbt']),
                     fmt.format(ab) if ab else '–', '{:+.2f}'.format(100 * (t_['N_gbt'] / ab - 1)).replace('-', '−') if ab else '–'])


def make_row(tmpl, vals):
    r_ = copy.deepcopy(tmpl)
    tcs = r_.findall(q('w:tc'))
    while len(tcs) < len(vals):
        r_.append(copy.deepcopy(tcs[-1])); tcs = r_.findall(q('w:tc'))
    for tc, v in zip(tcs, vals):
        ps = tc.findall(q('w:p'))
        for p2 in ps[1:]:
            tc.remove(p2)
        p0 = ps[0]
        for ch in list(p0):
            if ch.tag != q('w:pPr'):
                p0.remove(ch)
        if v:
            ins = cx.mark('ins')
            for rr in markup_runs(v):
                ins.append(rr)
            p0.append(ins)
        para_mark(cx, p0, 'ins')
    trPr = r_.find(q('w:trPr'))
    if trPr is None:
        trPr = etree.Element(q('w:trPr')); r_.insert(0, trPr)
    trPr.append(cx.mark('ins'))
    return r_


for k_, vals in enumerate(data):
    tbl.append(make_row(hdr_t if k_ == 0 else body_t, vals))
tcap.addnext(tbl)

# ---- Figure 10 (after Table 6)
fig10 = clean_copy(P[82])
assert fig10.find('.//' + q('w:drawing')) is not None
F10_PX = (3756, 3017)
c10x = cxmax; c10y = int(round(c10x * F10_PX[1] / F10_PX[0]))
for e in fig10.iter('{%s}extent' % WP, '{%s}ext' % A):
    e.set('cx', str(c10x)); e.set('cy', str(c10y))
ids = [int(e.get('id')) for e in root.iter('{%s}docPr' % WP, '{%s}cNvPr' % PIC)]
dp = fig10.find('.//{%s}docPr' % WP); dp.set('id', str(max(ids) + 1)); dp.set('name', 'Figure 10'); dp.set('descr', 'Figure 10. Modes of the four-pitch cells.')
cp_ = fig10.find('.//{%s}cNvPr' % PIC); cp_.set('id', str(max(ids) + 2)); cp_.set('name', 'Figure 10'); cp_.set('descr', 'fig/fig10_npitch.png')
fig10.find('.//{%s}blip' % A).set('{%s}embed' % RNS, 'rIdNpitch')
for r_ in list(fig10.findall(q('w:r'))):
    w_ = cx.mark('ins'); fig10.replace(r_, w_); w_.append(r_)
para_mark(cx, fig10, 'ins')
tbl.addnext(fig10)
caption_after(fig10, P[179], 'Figure 10.',
    ' Buckling modes of the four-pitch cells ($b_s$ = 150 mm, $\u03b2$ = 8, SS; GBT): deformed cell (colour: displacement normal '
    'to the walls) and cross-section at the largest displacement (dotted: stringer pitches) of the (a, b) first, (c, d) second '
    'and (e, f) fourth mode; (a, c, e) isotropic, (b, d, f) Sym1 cell. The Bloch wavenumbers $k_B$$b_s$ = \u03c0, \u03c0/2 and 0 '
    'were identified from the skin deflection shifted by one pitch; the third mode equals the second (double eigenvalue). '
    'Abaqus: eigenvalues of the four-pitch shell model (Table 6).')

tree.write(DST, xml_declaration=True, encoding='UTF-8', standalone=True)
print('ok, next id', cx.nid)
