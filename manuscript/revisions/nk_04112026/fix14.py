"""fix14 (13 Nov 2026): Abaqus modes 2-4 of the four-pitch cells exported by Aya and compared (check_npitch_modes.py):
wavenumbers measured from the Abaqus modes, so the "assigned by order" caveats are replaced by the result.
Based on fix13: Aya, 12 Nov 2026 ("dont make up abaqus numbers, ... make sure the paper description u wrote r correct, ... check the
all paper"): corrections found by the full audit, as tracked changes on the AME-12112026 state of Aya's own file."""
import copy, sys
sys.path.insert(0, '../nk')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-13T11:30:00Z"
ME = trk.AUTH
t = etree.parse('u/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); B = list(b)
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


KB = '$k_B$$b_s$'
E = [
 (183, "For kBbs = π/2 and 0, the Abaqus values are the second (equal to the third) and the fourth eigenvalue of the "
       "four-pitch shell model, assigned by their order and multiplicity; their mode shapes were not exported and not compared.",
       "The Abaqus modes confirm these wavenumbers: measured from their skin deflection in the same way, Abaqus mode 1 has "
       + KB + " = π, the double modes 2 and 3 have π/2 and lie entirely in the plane of the two GBT π/2 modes (subspace "
       "MAC 1.0000), and mode 4 has 0; the MAC of modes 1 and 4 with GBT is 0.9997–1.0000."),
 (188, "; for kBbs = π/2 and 0: its second (equal to the third) and fourth eigenvalue, assigned by order and multiplicity)",
       "; wavenumbers of its modes measured in the same way)"),
 (191, "Abaqus: eigenvalues 1, 2 (= 3) and 4 of the four-pitch shell model, assigned to the wavenumbers by their order (Table 6); "
       "dots in (a, b): its first mode, with the modal assurance criterion MAC of the skin deflections.",
       "Abaqus: modes 1, 2 (= 3) and 4 of the four-pitch shell model, wavenumbers measured in the same way (Table 6); dots in "
       "(a, b, e, f): Abaqus mode, with the modal assurance criterion MAC of the skin deflections; (c, d): subspace MAC = share "
       "of each Abaqus mode of the double pair in the plane of the two GBT modes."),
]
for i, old, new in E:
    e = B[i]; ps = [e] if e.tag == q('w:p') else [p for p in e.iter(q('w:p')) if old in vis_text(p)]
    assert len(ps) == 1, (i, old); edit(ps[0], old, new)
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print('ok', len(E), 'edits')
