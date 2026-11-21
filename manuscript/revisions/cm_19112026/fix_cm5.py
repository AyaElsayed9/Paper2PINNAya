"""fix_cm5 (21 Nov 2026, Aya: check why the professor commented, not only reply).  Runs after fix_cm2b, fix_cm3, fix_cm4.
C31: wall count, nodes and segments explained with the numbers of the code (18 natural + 74 intermediate nodes,
     93 = 92 - 1 + 2 closed cells); Fig. 1c numbers the walls and shows a web root (caption).
C37: "hat functions" replaced by a plain definition, also in Section 2.2.4 (the professor did not know the term).
C41: grammar of the remaining statement.
C58: the dots in Eq. (28) removed; the shear-mode geometric terms written out (unnumbered display equation, OMML).
C65: the meaning of the MAC in words.
C67: how the side-node differences were obtained (and why the warping value is 1e-11, not 1e-16)."""
import copy, sys, re
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-21T12:00:00Z"
ME = trk.AUTH
t = etree.parse('u/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns)
cx = Ctx(root)


def edit(p, old, new):
    """tracked replacement of the only visible occurrence of `old` in paragraph p by `new` (markup); text inside my
    own (same-author) insertions is edited directly, the rest is deleted with a w:del mark."""
    full = vis_text(p); assert full.count(old) == 1, (old, full.count(old)); pos = full.index(old); end = pos + len(old)
    mids = []; o = 0
    for r, s in vis_runs(p):
        a, c = max(pos, o), min(end, o + len(s))
        if a < c:
            mids.append((r, a - o, c - o))
        o += len(s)
    if any(r.tag == q('m:r') for r, _, _ in mids):          # number inside an equation: one math run only
        assert len(mids) == 1 and '$' not in new, (old, mids)
        r, a, c = mids[0]; pc = split(r, a, c)
        assert pc.getparent().tag != q('w:ins'), 'math run inside an insertion'
        w_ = wrap_del(cx, pc); ins = cx.mark('ins'); nr = copy.deepcopy(pc)
        nr.find(q('m:t')).text = new; ins.append(nr); w_.addnext(ins)
        assert vis_text(p) == full[:pos] + new + full[end:]
        return
    pieces = [split(r, a, c) for r, a, c in mids]
    rpr = pieces[0].find(q('w:rPr'))
    newruns = markup_runs(new, rpr) if new else []
    for k, pc in enumerate(pieces):
        last = k == len(pieces) - 1
        par = pc.getparent()
        if par.tag == q('w:ins') and par.get(q('w:author')) == ME:
            if last:
                for nr in newruns:
                    pc.addprevious(nr)
            par.remove(pc)
            continue
        if par.tag == q('w:ins'):                      # other author's insertion: isolate the piece in a clone
            after = list(pc.itersiblings())
            mid = etree.Element(q('w:ins'), attrib=dict(par.attrib)); mid.set(q('w:id'), str(cx.nid)); cx.nid += 1
            aft = etree.Element(q('w:ins'), attrib=dict(par.attrib)); aft.set(q('w:id'), str(cx.nid)); cx.nid += 1
            par.addnext(mid); mid.append(pc); mid.addnext(aft)
            for x in after:
                aft.append(x)
            if len(aft) == 0: aft.getparent().remove(aft)
            if len(par) == 0: par.getparent().remove(par)
        w_ = wrap_del(cx, pc)
        if last and newruns:
            ins = cx.mark('ins')
            for nr in newruns: ins.append(nr)
            host = w_.getparent(); (host if host.tag == q('w:ins') else w_).addnext(ins)
    import re; exp = full[:pos] + re.sub(r'\$([^$]+)\$', lambda m: m.group(1).replace('_', ''), new) + full[end:]
    got = vis_text(p); assert got == exp, (got[max(0, pos - 60):pos + 200], exp[max(0, pos - 60):pos + 200])



def own_replace(OLD, NEW):
    ps = [p for p in b.iter(q('w:p')) if OLD in vis_text(p)]; assert len(ps) == 1, OLD
    # OLD lies entirely in my own insertion of fix_cm1: change the run text in place
    p = ps[0]; v = vis_text(p); pos = v.find(OLD); o = 0; hit = []
    for r, s_ in vis_runs(p):
        if max(pos, o) < min(pos + len(OLD), o + len(s_)):
            hit.append((r, o))
        o += len(s_)
    assert all(r.getparent().tag == q('w:ins') and r.getparent().get(q('w:author')) == ME for r, _ in hit), "not own insertion"
    r0, o0 = hit[0]; tt = r0.find(q('w:t'))
    if len(hit) == 1:
        a = pos - o0; tt.text = tt.text[:a] + NEW + tt.text[a + len(OLD):]
    else:                                    # spread over several of my runs: put NEW in the first, cut the rest
        a = pos - o0; first_len = len(tt.text) - a; tt.text = tt.text[:a] + NEW; rem = len(OLD) - first_len
        for r, oo in hit[1:]:
            t2 = r.find(q('w:t')); k = min(rem, len(t2.text)); t2.text = t2.text[k:]; rem -= k
    assert NEW in vis_text(p)


def own_replace_m(OLD, NEW):
    """own_replace with $a_b$ markup in NEW: the run holding NEW is split into plain/italic/subscript runs"""
    tag = "\u2063"                                          # invisible marker around the new text
    own_replace(OLD, tag + NEW + tag)
    for tt in b.iter(q('w:t')):
        if tt.text and tt.text.count(tag) == 2:
            r = tt.getparent(); rpr = r.find(q('w:rPr'))
            pre, mid, post = tt.text.split(tag)
            runs = ([wrun(pre, rpr)] if pre else []) + markup_runs(mid, rpr) + ([wrun(post, rpr)] if post else [])
            for nr in runs:
                r.addprevious(nr)
            r.getparent().remove(r)
            return
    raise AssertionError("marker not found")


def one(text):
    ps = [p for p in b.iter(q('w:p')) if text in vis_text(p)]; assert len(ps) == 1, (text, len(ps)); return ps[0]


# C31 (own insertion of fix_cm1 replaced)
own_replace_m("The 19 walls are five skin walls, four foot walls (skin and flange together), four webs, two crowns and four short "
            "link walls that connect each foot at its mid-plane to its web at the true offset c; the segment ends are the nodes "
            "(natural nodes at the wall junctions, intermediate nodes inside the walls).",
            "The 19 walls, numbered in Figure 1c, are five skin walls (1–5: the half bays at the two side edges, the bay between "
            "the stringers and the skin under each stringer), four feet (6–9: the parts of the skin covered by a flange, skin and "
            "flange together), four webs (10–13), two crowns (14–15) and four short link walls (16–19), which connect the skin to "
            "the web roots at the offset $c$ (Section 2.1.2). The 18 natural nodes are the 16 wall junctions and the two side nodes. "
            "Every wall is divided into equal segments of at most $b_s$/25 (at least two per wall, one per link wall), which adds 74 "
            "intermediate nodes. This gives $N$ = 92 nodes and 93 segments: an open contour with $N$ nodes has $N$ − 1 segments, and "
            "each of the two closed cells formed by a stringer and the skin adds one.")
# Fig. 1 caption
edit(one("with natural and intermediate nodes, wall types and local axes"), "with natural and intermediate nodes, wall types and local axes",
     "with natural and intermediate nodes, the numbered walls (1–5 skin, 6–9 feet, 10–13 webs, 14–15 crowns, 16–19 link "
     "walls), a detail of a web root (not to scale) and the local axes")
# C37
own_replace("the nodal hat functions of the warping (each equal to 1 at one node and 0 at all other nodes, varying linearly along "
            "the segments),",
            "one warping function per node, equal to 1 at that node, 0 at all other nodes and linear along each segment (as for "
            "the shear modes in [36,37]);")
edit(one("one replaced by the axial mode"), "one replaced by the axial mode", "one of them is replaced by the axial mode")
edit(one("(hat function "), "(hat function ", "(warping function ")
# C41
edit(one("but one is set to zero"), "but one is set to zero", "but one are set to zero")
# C65
own_replace("MAC = 1 means identical shapes (up to a factor) and MAC = 0 orthogonal shapes.",
            "In words, the MAC is the squared cosine of the angle between the two vectors of skin deflections: it does not depend "
            "on how the modes are scaled, it is 1 when the two shapes are identical and 0 when they are completely different "
            "(orthogonal).")
# C67
own_replace_m("the largest absolute differences between the values at the two side edges over all conventional modes are",
            "the largest absolute differences between the values of a mode at the two side nodes $s_L$ and $s_R$, taken over all "
            "conventional modes scaled to a unit maximum, are")
edit(one("and zero for the shear and Mindlin modes. The explicit tensors"), "and zero for the shear and Mindlin modes.",
     "and zero for the shear and Mindlin modes, whose side-node functions are merged. The first two values are round-off errors "
     "of double-precision arithmetic; the warping follows from the numerical solution of Eq. (13), whose accuracy limits it to "
     "about 10⁻¹¹.")
edit(one("Every mode is periodic to round-off:"), "periodic to round-off:", "periodic to numerical precision:")
# C58: Eq. (28) without the dots, the shear-mode terms written out
p28 = [p for p in b.iter(q('w:p')) if '+…,' in vis_text(p)]; assert len(p28) == 1
v = vis_text(p28[0]); k = v.index('+…,'); replace_span(cx, p28[0], k, k + 2, '')
pd = one("where the dots stand for the shear-mode terms.")
OLD = "where the dots stand for the shear-mode terms. "
v = vis_text(pd); assert v.startswith(OLD), v[:80]
pieces = replace_span(cx, pd, 0, len(OLD), "")          # tracked deletion of the old sentence
# new paragraph P1 (inserted mark) before pd: the moved deleted sentence + the new lead-in
P1 = etree.Element(q('w:p')); pPr = pd.find(q('w:pPr'))
if pPr is not None:
    pp = copy.deepcopy(pPr)
    for e in pp.findall(q('w:rPr')): pp.remove(e)
    P1.append(pp)
para_mark(cx, P1, 'ins')
moved = []
for ch in list(pd):                                     # everything up to the first visible old text: the deleted
    if ch.tag == q('w:pPr'):                            # sentence with its comment ranges and references
        continue
    if ch.tag == q('m:oMath') or (ch.tag == q('w:r') and ''.join(x.text or '' for x in ch.iter(q('w:t'))).strip()):
        break
    moved.append(ch)
for ch in moved:
    pd.remove(ch); P1.append(ch)
ins = cx.mark('ins')
for r in markup_runs("The shear modes and their coupling with the conventional modes add, in the same way,", first_rpr(pd)):
    ins.append(r)
P1.append(ins); pd.addprevious(P1)
# display equation (OMML generated by pandoc from LaTeX), every math run tracked as inserted
M_ = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
import subprocess, tempfile, zipfile, os
_tmp = tempfile.mkdtemp(); subprocess.run(['pandoc', 'eq/c58_shear_terms.md', '-o', os.path.join(_tmp, 'eq.docx')], check=True)
eqd = etree.fromstring(zipfile.ZipFile(os.path.join(_tmp, 'eq.docx')).read('word/document.xml')).getroottree()
for c_ in list(eqd.iter('{%s}count' % M_)):           # pandoc writes m:count after m:mcJc; the schema wants it first
    par_ = c_.getparent(); par_.remove(c_); par_.insert(0, c_)
for r_ in eqd.iter('{%s}r' % M_):                        # same font as the other equations
    if r_.find(q('w:rPr')) is None:
        rp_ = etree.Element(q('w:rPr')); f_ = etree.SubElement(rp_, q('w:rFonts'))
        f_.set(q('w:ascii'), 'Cambria Math'); f_.set(q('w:hAnsi'), 'Cambria Math')
        mr_ = r_.find('{%s}rPr' % M_)
        (mr_.addnext(rp_) if mr_ is not None else r_.insert(0, rp_))
paras = eqd.findall('.//{%s}oMathPara' % M_)
for k, omp in enumerate(paras):
    PE = etree.Element(q('w:p'))
    if pPr is not None:
        pp = copy.deepcopy(pPr)
        for e in pp.findall(q('w:rPr')): pp.remove(e)
        for e in pp.findall(q('w:ind')): pp.remove(e)
        PE.append(pp)
    para_mark(cx, PE, 'ins')
    omp = copy.deepcopy(omp)
    for r in list(omp.iter('{%s}r' % M_)):
        w_ = cx.mark('ins'); r.addprevious(w_); w_.append(r)
    for cp in omp.iter('{%s}ctrlPr' % M_):
        mark_ctrl(cx, cp, 'ins')
    PE.append(omp)
    pd.addprevious(PE)
    if k == 0:
        PW = etree.Element(q('w:p'))
        if pPr is not None:
            pp = copy.deepcopy(pPr)
            for e in pp.findall(q('w:rPr')): pp.remove(e)
            PW.append(pp)
        para_mark(cx, PW, 'ins'); ins = cx.mark('ins')
        for r in markup_runs("with, for any two conventional or shear modes a and b,", first_rpr(pd)):
            ins.append(r)
        PW.append(ins); pd.addprevious(PW)
# lead-in of the remaining paragraph
ins = cx.mark('ins')
rp0 = first_rpr(pd)
for txt, kw in (("so that ", {}), ("X", dict(italic=True)), ("U", dict(italic=True, sup=True)), (" and ", {}),
                ("X", dict(italic=True)), ("R", dict(italic=True, sup=True)),
                (" of Eq. (28) are the special case of two conventional modes; the Mindlin modes do not enter, since they do "
                 "not change the mid-surface displacements.", {})):
    ins.append(wrun(txt, rp0, **kw))
first = [ch for ch in pd if ch.tag != q('w:pPr')][0]
first.addprevious(ins)
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print("ok fix_cm5")
