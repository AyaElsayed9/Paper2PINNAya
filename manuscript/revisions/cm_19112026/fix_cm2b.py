"""fix_cm2b (21 Nov 2026, replaces fix_cm2; Aya: the other groups' references must be merged with ours, not added
as a separate block).  C22: the literature paragraph is rewritten by topic, with the four checked references of
other groups (Bisagni & Vescovini TWS 2009; Stamatelos, Labeas & Tserpes TWS 2011; Buermann, Rolfes, Tessmer &
Schagerl TWS 2006; Vescovini & Bisagni Aeronaut. J. 2015) cited where they belong; all references renumbered in the
order of first citation (tracked), reference list block 9-20 replaced accordingly.  Also: the load-distribution
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
NEWNUM = {9: 15, 10: 11, 11: 12, 12: 16, 13: 13, 14: 14, 15: 17, 16: 20}
def mp(n):
    return n if n <= 8 else NEWNUM.get(n, n + 4)
def group(nums):
    """'[a,b-c]' text for sorted numbers, ranges of three or more joined by an en dash"""
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if j == i else (f"{nums[i]},{nums[j]}" if j == i + 1 else f"{nums[i]}–{nums[j]}"))
        i = j + 1
    return "[" + ",".join(out) + "]"
def expand(g):
    nums = []
    for part in re.split(r"\s*,\s*", g.strip("[]")):
        ab = re.split(r"\s*[–-]\s*", part)
        nums += list(range(int(ab[0]), int(ab[1]) + 1)) if len(ab) == 2 else [int(ab[0])]
    return nums
# C22 paragraph: rewritten by topic (new numbers directly)
c22 = [p for p in ps[:iref] if "An overview of these methods is given in [16]." in vis_text(p)]
assert len(c22) == 1; p22 = c22[0]
NEW22 = ("Closed-form and semi-analytical models give fast buckling estimates for stiffened composite panels and are "
         "therefore well suited to preliminary design. For blade-stringer-stiffened plates, closed-form buckling loads were "
         "derived for a repeating unit with periodic boundary conditions [5] and with the interaction between stiffener and "
         "plate [6], and Ritz-based and closed-form formulations, with the stiffeners represented as torsion bars or elastic "
         "edge restraints, were developed for the local buckling and postbuckling of stiffened laminated panels [9,10]. For "
         "omega stringers, closed-form and semi-analytical solutions are available for plates braced by omega stringers [3], "
         "for local buckling [4], column buckling [7], shear buckling [11] and combined compression and shear [12], and the "
         "effect of the stringer feet on local buckling was included in [13,14]. Postbuckling was studied for "
         "hat-stringer-stiffened plates under transverse compression [15], for omega-stringer-stiffened panels [16,17] and for "
         "stringer- and frame-stiffened cylindrical panels [18]. Such fast models are also the basis of design procedures for "
         "omega-stiffened composite panels [19]; an overview of these methods is given in [20].")
rewrite_paragraph(cx, p22, [NEW22])
# renumber every other citation group
pat = re.compile(r'\[(\d+(?:\s*[–,-]\s*\d+)*)\]')
n_groups = 0
for p in ps[:iref]:
    if p is p22:
        continue
    v = vis_text(p); ms = []
    for m in pat.finditer(v):
        nums = expand(m.group(0))
        if [mp(n) for n in nums] != nums:
            ms.append(m)
    k64 = v.find("(equal line loads, as in [13,18])")          # my own insertion (fix_cm1): edited directly below
    if k64 >= 0:
        ms = [m for m in ms if not (k64 <= m.start() < k64 + 40)]
    for m in reversed(ms):
        new = group([mp(n) for n in expand(m.group(0))])
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
# reference list: old 9-16 deleted, new block 9-20 inserted in the order of first citation
refs = [p for p in ps[iref + 1:] if vis_text(p).strip()]
old = {k + 1: refs[k] for k in range(len(refs))}
assert vis_text(old[16]).startswith("Mittelstedt, C. Buckling of Stiffened Laminates"), vis_text(old[16])[:60]
assert vis_text(old[9]).startswith("Mittelstedt, C.; Schröder"), vis_text(old[9])[:60]
copies = {k: copy.deepcopy(old[k]) for k in range(9, 17)}
for k in range(9, 17):
    delete_all(cx, old[k]); para_mark(cx, old[k], 'del')
NEW = {  # (text, format) segments in the style of the list: journal italic, year bold, volume italic
    "A": [("Bisagni, C.; Vescovini, R. Analytical formulation for local buckling and post-buckling analysis of stiffened "
      "laminated panels. ", None), ("Thin-Walled Structures", "i"), (" ", None), ("2009", "b"), (", ", None), ("47", "i"),
     (", 318–334. https://doi.org/10.1016/j.tws.2008.07.006", None)],
    "B": [("Stamatelos, D.G.; Labeas, G.N.; Tserpes, K.I. Analytical calculation of local buckling and post-buckling behavior "
      "of isotropic and orthotropic stiffened panels. ", None), ("Thin-Walled Structures", "i"), (" ", None), ("2011", "b"),
     (", ", None), ("49", "i"), (", 422–430. https://doi.org/10.1016/j.tws.2010.11.008", None)],
    "C": [("Buermann, P.; Rolfes, R.; Tessmer, J.; Schagerl, M. A semi-analytical model for local post-buckling analysis of "
      "stringer- and frame-stiffened cylindrical panels. ", None), ("Thin-Walled Structures", "i"), (" ", None), ("2006", "b"),
     (", ", None), ("44", "i"), (", 102–114. https://doi.org/10.1016/j.tws.2005.08.010", None)],
    "D": [("Vescovini, R.; Bisagni, C. A fast procedure for the design of composite stiffened panels. ", None),
     ("The Aeronautical Journal", "i"), (" ", None), ("2015", "b"), (", ", None), ("119", "i"),
     (", 185–201. https://doi.org/10.1017/S0001924000010332", None)]}
def fmt_run(text, f):
    r = etree.Element(q('w:r'))
    if f:
        rp = etree.SubElement(r, q('w:rPr'))
        etree.SubElement(rp, q('w:' + f)); etree.SubElement(rp, q('w:' + f + 'Cs'))
    tt = etree.SubElement(r, q('w:t')); tt.text = text
    tt.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    return r
def as_inserted(pc):
    """copy of an old reference paragraph as a tracked insertion"""
    for el in list(pc.iter(q('w:ins'), q('w:del'))):
        par = el.getparent()
        if par.tag == q('w:rPr'):
            par.remove(el)
    runs = [r for r in pc if r.tag == q('w:r')]
    ins = cx.mark('ins')
    for r in runs:
        pc.remove(r); ins.append(r)
    pc.append(ins); para_mark(cx, pc, 'ins')
    return pc
prev = old[16]
for key in ["A", "B", 10, 11, 13, 14, 9, 12, 15, "C", "D", 16]:
    if isinstance(key, str):
        prev = new_para_after(cx, prev, old[16], None)
        ins = cx.mark('ins')
        for text, f in NEW[key]:
            ins.append(fmt_run(text, f))
        prev.append(ins)
    else:
        pc = as_inserted(copies[key]); prev.addnext(pc); prev = pc
# Paper 1 named in the load-distribution sentence (inside the C64 insertion, now "[13,22]")
p64 = [p for p in ps[:iref] if "(equal line loads, as in [13,18])" in vis_text(p)]
assert len(p64) == 1
edit(p64[0], "(equal line loads, as in [13,18])", "(equal line loads, as in our previous study [22] and in [13])")
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print("ok:", n_groups, "citation groups renumbered, C22 paragraph rewritten, reference block 9-20 replaced")
