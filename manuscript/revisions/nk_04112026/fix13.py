"""Aya, 12 Nov 2026 ("dont make up abaqus numbers, ... make sure the paper description u wrote r correct, ... check the
all paper"): corrections found by the full audit, as tracked changes on the AME-12112026 state of Aya's own file."""
import copy, sys
sys.path.insert(0, '../nk')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-12T09:00:00Z"
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
 (8, "with mode shapes matching", "and the mode shapes of the six cells compared match"),
 (15, "Schilling and Mittelstedt extended", "Schilling, Mittelstedt and co-workers extended"),
 (16, "All these models reduce", "These models reduce"),
 (16, "is included only approximately", "is neglected or included only approximately"),
 (17, " and then reformulated as a systematic procedure that can be automated [31,32]. This work also clarified the role "
      "of the shear deformation modes [33] and is implemented in the code GBTul 2.0 [34].",
      ". New approaches to the calculation of the deformation modes [31], a procedure that rationalises and automates the "
      "cross-section analysis [32] and a study of the shear deformation modes [33] followed, and the code GBTul 2.0 performs "
      "GBT-based buckling and vibration analyses of thin-walled members [34]."),
 (17, ", and Taig and Ranzi formulated it for", ". Taig and Ranzi formulated GBT for"),
 (18, " with only a few unknowns per mode", ""),
 (19, "for members whose cross-section is either open, with free longitudinal edges, or closed.",
      "for open cross-sections, whose longitudinal edges are independent of each other, and for closed cross-sections."),
 (159, "1.00", "0.98"), (159, "1.16", "1.19"),
 (159, "0.86", "0.84"), (159, "0.89", "0.90"),
 (183, "as for the short cells.", "as for the short cells. For " + KB + " = π/2 and 0, the Abaqus values are the second (equal "
       "to the third) and the fourth eigenvalue of the four-pitch shell model, assigned by their order and multiplicity; their "
       "mode shapes were not exported and not compared."),
 (188, "equals that of the two-pitch cell)", "equals that of the two-pitch cell; for " + KB + " = π/2 and 0: its second (equal "
       "to the third) and fourth eigenvalue, assigned by order and multiplicity)"),
 (191, "Abaqus: eigenvalues of the four-pitch shell model (Table 6)", "Abaqus: eigenvalues 1, 2 (= 3) and 4 of the four-pitch "
       "shell model, assigned to the wavenumbers by their order (Table 6)"),
]
for i, old, new in E:
    e = B[i]; ps = [e] if e.tag == q('w:p') else [p for p in e.iter(q('w:p')) if old in vis_text(p)]
    assert len(ps) == 1, (i, old); edit(ps[0], old, new)
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print('ok', len(E), 'edits')
