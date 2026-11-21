"""fix_cm3 (20 Nov 2026): clearer wording for C47, C49, C53 and C58 (plain-language sentences added in front of
the technical text; no equation is changed).  Based on fix_cm1: Prof. Mittelstedt's comments of 19 Nov (GBT_periodic-NK-CM-19112026) - the short,
unambiguous ones, as tracked changes under Aya's name on her file.  The longer rewrites (C22 references, C32 figure,
C47, C49, C53, C58) are drafted in replies_draft.md first.  C64 numbers: scripts/c64_uniform_strain.py."""
import copy, sys
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-20T13:30:00Z"
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



E = [  # C58 (dots) moved to fix_cm5 (21 Nov: the terms are now written out)
 # C47
 ("In detail, with ", "A nodal field v is periodic when its two side nodes move identically, i.e. when Rv = 0. In detail, with "),
 ("are replaced by the largest periodic subspace", "are replaced by all their combinations whose side-edge differences "
  "cancel, i.e. by the largest periodic subspace"),
 # C49
 ("The condition ", "In an open cross-section, the warping of every conventional mode follows from its in-plane displacements, "
  "because the membrane shear strain is set to zero. Joining the two side edges closes the contour: going once around it, "
  "the warping has to return to its starting value, so in general a constant shear flow around the closed contour remains, "
  "as in closed thin-walled sections, and the warping is found by minimising the shear energy under this closure "
  "condition. The condition "),
 # C53
 ("with Bloch–Floquet theory [8]: a cell of width", "with Bloch–Floquet theory [8], the theory of waves in periodic structures. In an infinitely wide panel, a buckling "
  "mode repeats from bay to bay up to a phase shift; a cell of two pitches therefore contains the two extreme cases, modes "
  "that repeat in every bay and modes that alternate from bay to bay. Formally, a cell of width"),
 # C58
]
done = 0
for old, new in E:
    ps = [p for p in b.iter(q('w:p')) if old in vis_text(p)]
    assert len(ps) == 1, (old, len(ps))
    edit(ps[0], old, new); done += 1
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print('ok', done, 'edits')
