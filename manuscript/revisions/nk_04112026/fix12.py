"""Aya, 12 Nov 2026: remove the strip-method sentence (Navid did not understand it) and its references [9]-[11];
[12] (Lau & Hancock) stays because Section 2.3.3 cites it, and moves to its place of first citation; all citations
renumbered; Scheme 1 caption shortened.  Works on Aya's own file (her edits kept)."""
import copy, re, sys
sys.path.insert(0, '/tmp/claude-0/-home-claude/f2feee3d-f406-5802-a982-08c43a6228ed/scratchpad/nk')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-12T08:00:00Z"
t = etree.parse('u/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); B = list(b)
cx = Ctx(root)
ME = trk.AUTH


def unhighlight(r):
    rp = r.find(q('w:rPr'))
    if rp is not None:
        for tag in ('w:highlight', 'w:rPrChange'):
            for e in rp.findall(q(tag)):
                rp.remove(e)


# ---- 1. strip-method sentence: it was my own insertion -> remove it (no deletion mark needed)
p14 = B[14]
hl = [r for r in p14.iter(q('w:r')) if r.find('w:rPr/w:highlight', ns) is not None]
assert ''.join(rtext(r) for r in hl).startswith('Prismatic plate assemblies') and ''.join(rtext(r) for r in hl).endswith('along the length.')
for r in hl:
    a = r.getparent(); assert a.tag == q('w:ins') and a.get(q('w:author')) == ME
    a.remove(r)
# trailing space after "[8]."
runs = vis_runs(p14); r_last, s_last = runs[-1]
if s_last.endswith(' '):
    for tt in r_last.findall(q('w:t')):
        tt.text = tt.text.rstrip()
assert vis_text(p14).endswith('describes the response of the whole structure [8].'), vis_text(p14)[-80:]

# ---- 2. Scheme 1 caption: shorter (my own insertion -> replaced in place)
p24 = B[24]
rs = [r for r in p24.iter(q('w:r')) if r.find('w:rPr/w:highlight', ns) is not None]
assert rtext(rs[0]) == 'Scheme 1.'
for r in rs:
    unhighlight(r)
ins = rs[1].getparent(); rp0 = copy.deepcopy(rs[1].find(q('w:rPr')))
for r in rs[1:]:
    r.getparent().remove(r)
for r in markup_runs(' Flowchart of the periodic GBT–Ritz analysis; the frame marks the periodic modes, the contribution of '
                     'this study (drawings: Unsym cell, $b_s$ = 150 mm, $β$ = 1.5, SS).', rp0):
    ins.append(r)

# ---- 3. reference list: delete [9]-[11], move [12] after old [54]
R = [i for i, e in enumerate(B) if vis_text(e).strip() == 'References'][-1]
refs = [e for e in B[R + 1:] if e.find('.//w:pStyle', ns) is not None and e.find('.//w:pStyle', ns).get(q('w:val')) == 'MDPI81references']
assert len(refs) == 56
assert vis_text(refs[7]).startswith('Mead') and vis_text(refs[8]).startswith('Wittrick') and vis_text(refs[9]).startswith('Kennedy') and vis_text(refs[10]).startswith('Cheung')
assert vis_text(refs[11]).startswith('Lau') and vis_text(refs[54]).strip() and vis_text(refs[12]).startswith('Mittelstedt')
lau_new = copy.deepcopy(refs[11])


def del_para(p):
    for r, _ in vis_runs(p):
        wrap_del(cx, r)
    para_mark(cx, p, 'del')


for k in (8, 9, 10, 11):                 # old [9], [10], [11], and [12] at its old place
    del_para(refs[k])
for e in list(lau_new.iter()):
    for a in list(e.attrib):
        if a.endswith('}paraId') or a.endswith('}textId'):
            del e.attrib[a]
for r in list(lau_new.iter(q('w:r'))):
    w_ = cx.mark('ins'); r.getparent().replace(r, w_); w_.append(r)
para_mark(cx, lau_new, 'ins')
refs[53].addnext(lau_new)                # after old [54]

# ---- 4. citations in the text: old -> new numbers (order of first citation)
MAP = {n: n for n in range(1, 9)}
MAP.update({n: n - 4 for n in range(13, 55)}); MAP[12] = 51; MAP[55] = 52; MAP[56] = 53
cit = re.compile(r'\[(\d+(?:\s*[,–\-]\s*\d+)*)\]')
nchg = 0
nown = [0]


def renum(c, pos, end, new):
    """replace visible characters [pos, end) of container c (inside ONE run) by `new`, keeping the run in place:
    in my own insertion -> edit the text directly; in another author's insertion -> isolate the run in a clone of that
    insertion, delete it there and insert the new text right after; in original text -> tracked deletion + insertion."""
    o = 0
    for r, s_ in vis_runs(c):
        if o <= pos and end <= o + len(s_):
            a_, b_ = pos - o, end - o
            break
        o += len(s_)
    else:
        raise RuntimeError('citation split over runs')
    par = r.getparent()
    if par.tag == q('w:ins') and par.get(q('w:author')) == ME:
        tn = q('w:t') if r.tag == q('w:r') else q('m:t')
        tt = r.find(tn); tt.text = tt.text[:a_] + new + tt.text[b_:]; nown[0] += 1
        return
    piece = split(r, a_, b_)
    par = piece.getparent()
    if par.tag == q('w:ins'):
        after = list(piece.itersiblings())
        mid = etree.Element(q('w:ins'), attrib=dict(par.attrib)); mid.set(q('w:id'), str(cx.nid)); cx.nid += 1
        aft = etree.Element(q('w:ins'), attrib=dict(par.attrib)); aft.set(q('w:id'), str(cx.nid)); cx.nid += 1
        par.addnext(mid); mid.append(piece); mid.addnext(aft)
        for x in after:
            aft.append(x)
        if len(aft) == 0:
            aft.getparent().remove(aft)
        if len(par) == 0:
            par.getparent().remove(par)
    w_ = wrap_del(cx, piece)
    ins_ = cx.mark('ins'); nr = copy.deepcopy(piece)
    if nr.tag == q('w:r'):
        for tt in nr.findall(q('w:delText')):
            tt.tag = q('w:t')
    tn = q('w:t') if nr.tag == q('w:r') else q('m:t')
    nr.find(tn).text = new; ins_.append(nr)
    host = w_.getparent()
    (host if host.tag == q('w:ins') else w_).addnext(ins_)



for e in B[:R]:
    s = vis_text(e)
    for m in reversed(list(cit.finditer(s))):
        inner = m.group(1)
        new = re.sub(r'\d+', lambda mm: str(MAP[int(mm.group(0))]), inner)
        if new != inner:
            # replace only the digits that change, token by token (from the right)
            toks = list(re.finditer(r'\d+', inner))
            for tk in reversed(toks):
                o, nw = tk.group(0), str(MAP[int(tk.group(0))])
                if o != nw:
                    a0 = m.start(1) + tk.start()
                    renum(e, a0, a0 + len(o), nw); nchg += 1
print('citation numbers changed:', nchg, 'inside my own insertions:', nown[0])
# check: all citations now in order of first appearance and within 1..53
order = []
for e in B[:R]:
    for m in cit.finditer(vis_text(e)):
        for part in re.split(r'\s*,\s*', m.group(1)):
            a = re.split(r'\s*[–\-]\s*', part)
            for n in range(int(a[0]), int(a[-1]) + 1):
                if n not in order:
                    order.append(n)
print(order)
assert order == list(range(1, 54)), order
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)

# ---- 5. my replies to comments 7 and 8
cm = open('u/word/comments.xml', encoding='utf8').read()
a1 = 'The last two lines are rewritten: the periodic condition is now explained, and [9–12] are described as strip methods for prismatic plate assemblies.'
b1 = 'The last two lines are rewritten: the periodic condition is now explained, and the sentence on strip methods was removed with its references (former [9]–[11]); the references are renumbered.'
a2 = 'The omega stringer, the periodic cell and the strip methods moved to a new second paragraph.'
b2 = 'The omega stringer and the periodic cell moved to a new second paragraph.'
assert cm.count(a1) == 1 and cm.count(a2) == 1, (cm.count(a1), cm.count(a2))
open('u/word/comments.xml', 'w', encoding='utf8').write(cm.replace(a1, b1).replace(a2, b2))
print('ok')
