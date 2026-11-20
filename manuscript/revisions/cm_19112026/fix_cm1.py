"""fix_cm1 (20 Nov 2026): Prof. Mittelstedt's comments of 19 Nov (GBT_periodic-NK-CM-19112026) - the short,
unambiguous ones, as tracked changes under Aya's name on her file.  The longer rewrites (C22 references, C32 figure,
C47, C49, C53, C58) are drafted in replies_draft.md first.  C64 numbers: scripts/c64_uniform_strain.py."""
import copy, sys
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-20T12:00:00Z"
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



E = [
 # C31
 ("(Figure 1c).", "(Figure 1c). The 19 walls are five skin walls, four foot walls (skin and flange together), four webs, two "
  "crowns and four short link walls that connect each foot at its mid-plane to its web at the true offset c; the segment "
  "ends are the nodes (natural nodes at the wall junctions, intermediate nodes inside the walls)."),
 # C37
 ("the nodal hat functions of the warping,", "the nodal hat functions of the warping (each equal to 1 at one node and 0 at "
  "all other nodes, varying linearly along the segments),"),
 # C41
 ("Proposition 1. ", ""),
 # C46
 ("of dimension 3N-3. Any basis", None),
 # C50
 ("[20,36]. With the segment difference matrix", "[20,36]. For every segment p with the end nodes q and q′, the "
  "difference of the nodal warping values is in general not zero; only the two side edges have equal warping. With the "
  "segment difference matrix"),
 # C54
 ("modes 1 and 2 are rigid translations;", "modes 1 and 2 are rigid translations of the cross-section in Z and Y, i.e. the "
  "out-of-plane and in-plane bending of the panel along X;"),
 # C61
 ("unlike the Hermite elements of GBT-based finite elements [34], they are twice continuously differentiable.",
  "unlike the cubic Hermite elements of GBT-based finite elements [34], which have a continuous slope but jumps of the "
  "curvature between elements (C1), they are twice continuously differentiable across the knots (C2)."),
 # C64
 (". The shell models are those of", ". Because the analysis is linear, the magnitude of this reference "
  "load is arbitrary; its distribution (equal line loads, as in [13,18]) is a choice. A distribution following a uniform end "
  "strain, with the line load of every wall proportional to its axial membrane stiffness, changes the critical mean force "
  "per unit width by only 0.3–0.5% for the Cross-ply and Unsym lay-ups and by 1.8–4.3% for the angle-ply Sym1 and Sym2 "
  "lay-ups (GBT, 32 cells). The shell models are those of"),
 # C65
 ("MAC = 1 means identical shapes.", "MAC = 1 means identical shapes (up to a factor) and MAC = 0 orthogonal shapes."),
 # C67
 ("the largest side-node differences of the conventional modes are", "the largest absolute differences between the values "
  "at the two side edges over all conventional modes are"),
 # C71
 ("(e) amplitude functions and (f) in-plane shapes of the four dominant modes.", "(e) amplitude functions and (f) in-plane "
  "shapes of the four dominant modes, which are different cross-section deformations; their amplitude functions have the "
  "same number of half-waves, and their signs show how they combine in the critical mode."),
]
done = 0
for old, new in E:
    ps = [p for p in b.iter(q('w:p')) if old in vis_text(p)]
    assert len(ps) == 1, (old, len(ps))
    if new is None:                                  # C46: insertion after the math, before " Any basis"
        new = old.replace(". Any basis", ", with N the number of nodes. Any basis")
        old2 = ". Any basis"; new2 = ", with N the number of nodes. Any basis"
        edit(ps[0], old2, new2)
    else:
        edit(ps[0], old, new)
    done += 1
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print('ok', done, 'edits')
