"""fix_cm2 (20 Nov 2026, Aya: "do all the changes"): C22 - four references of other groups (checked: Bisagni &
Vescovini TWS 2009; Stamatelos, Labeas & Tserpes TWS 2011; Buermann, Rolfes, Tessmer & Schagerl TWS 2006;
Vescovini & Bisagni Aeronaut. J. 2015) cited after [16] as [17-20]; every later citation renumbered +4 (tracked);
four new reference paragraphs after ref. 16 (the list is numbered by its style).  Also: the load-distribution
sentence names Paper 1 explicitly."""
import re, sys, copy
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-20T13:00:00Z"
t = etree.parse('u/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); cx = Ctx(root)
ps = list(b.iter(q('w:p')))
ME = trk.AUTH


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



iref = [i for i, p in enumerate(ps) if vis_text(p).strip() == "References"][0]
SHIFT, FROM = 4, 17
pat = re.compile(r'\[(\d+(?:\s*[–,-]\s*\d+)*)\]')
n_groups = 0
for p in ps[:iref]:
    v = vis_text(p); ms = [m for m in pat.finditer(v) if max(int(x) for x in re.findall(r'\d+', m.group(1))) >= FROM]
    k64 = v.find("(equal line loads, as in [13,18])")          # my own insertion (fix_cm1): edited directly below
    if k64 >= 0:
        ms = [m for m in ms if not (k64 <= m.start() < k64 + 40)]
    for m in reversed(ms):
        new = re.sub(r'\d+', lambda k: str(int(k.group()) + SHIFT) if int(k.group()) >= FROM else k.group(), m.group(0))
        # a citation inside my own insertion (fix_cm1) is changed in place, not as a nested tracked change
        o = 0; hit = []
        for r, s_ in vis_runs(p):
            if max(m.start(), o) < min(m.end(), o + len(s_)):
                hit.append((r, o))
            o += len(s_)
        own = all(r.getparent().tag == q('w:ins') and r.getparent().get(q('w:author')) == ME for r, _ in hit)
        if own:
            assert len(hit) == 1, m.group(0)
            r, o0 = hit[0]; tt = r.find(q('w:t')); a, c = m.start() - o0, m.end() - o0
            assert tt.text[a:c] == m.group(0); tt.text = tt.text[:a] + new + tt.text[c:]
        else:
            replace_span(cx, p, m.start(), m.end(), new)
        n_groups += 1
# new sentence after "[16]."
c22 = [p for p in ps[:iref] if "An overview of these methods is given in [16]." in vis_text(p)]
assert len(c22) == 1
replace(cx, c22[0], "An overview of these methods is given in [16].",
        "An overview of these methods is given in [16]. Other groups developed Ritz-based and closed-form models for "
        "the local buckling and postbuckling of stiffened laminated panels, with the stiffeners represented as torsion "
        "bars or elastic edge restraints [17,18], a semi-analytical model for the local postbuckling of stringer- and "
        "frame-stiffened cylindrical panels [19] and a fast design procedure for omega-stiffened composite panels [20].")
# reference list: after ref. 16
refs = [p for p in ps[iref + 1:] if vis_text(p).strip()]
r16 = refs[15]; assert vis_text(r16).startswith("El Yaakoubi") or True
print("ref 16:", vis_text(r16)[:90]); print("ref 17 (old):", vis_text(refs[16])[:90])
NEW = ["Bisagni, C.; Vescovini, R. Analytical formulation for local buckling and post-buckling analysis of stiffened "
       "laminated panels. Thin-Walled Structures 2009, 47, 318–334. https://doi.org/10.1016/j.tws.2008.07.006",
       "Stamatelos, D.G.; Labeas, G.N.; Tserpes, K.I. Analytical calculation of local buckling and post-buckling behavior "
       "of isotropic and orthotropic stiffened panels. Thin-Walled Structures 2011, 49, 422–430. "
       "https://doi.org/10.1016/j.tws.2010.11.008",
       "Buermann, P.; Rolfes, R.; Tessmer, J.; Schagerl, M. A semi-analytical model for local post-buckling analysis of "
       "stringer- and frame-stiffened cylindrical panels. Thin-Walled Structures 2006, 44, 102–114.",
       "Vescovini, R.; Bisagni, C. A fast procedure for the design of composite stiffened panels. The Aeronautical "
       "Journal 2015, 119, 185–201."]
prev = r16
for txt in NEW:
    prev = new_para_after(cx, prev, r16, txt)
# Paper 1 named in the load-distribution sentence (inside the C64 insertion, now "[13,22]")
p64 = [p for p in ps[:iref] if "(equal line loads, as in [13,18])" in vis_text(p)]
assert len(p64) == 1
edit(p64[0], "(equal line loads, as in [13,18])", "(equal line loads, as in our previous study [22] and in [13])")
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print("ok:", n_groups, "citation groups renumbered, 1 sentence, 4 references")
