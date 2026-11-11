"""Tracked-change helpers for the manuscript docx (lxml)."""
import copy
import re
from lxml import etree

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
XML = 'http://www.w3.org/XML/1998/namespace'
ns = {'w': W, 'm': M}
AUTH = "Aya Medhat Elsayed"
DATE = "2026-11-10T11:00:00Z"


def q(t):
    p, n = t.split(':')
    return '{%s}%s' % ({'w': W, 'm': M}[p], n)


class Ctx:
    def __init__(self, root):
        self.root = root
        ids = [int(v) for e in root.iter() for k, v in e.attrib.items()
               if k == q('w:id') and v.lstrip('-').isdigit()]
        self.nid = max(ids) + 1000

    def mark(self, kind):
        e = etree.Element(q('w:' + kind))
        e.set(q('w:id'), str(self.nid)); self.nid += 1
        e.set(q('w:author'), AUTH); e.set(q('w:date'), DATE)
        return e


def in_del(el):
    a = el.getparent()
    while a is not None:
        if a.tag == q('w:del'):
            return True
        a = a.getparent()
    return False


def rtext(r):
    tn = q('w:t') if r.tag == q('w:r') else q('m:t')
    return ''.join(t.text or '' for t in r.findall(tn))


def vis_runs(c):
    out = []
    for r in c.iter(q('w:r'), q('m:r')):
        if in_del(r):
            continue
        s = rtext(r)
        if s:
            out.append((r, s))
    return out


def vis_text(c):
    return ''.join(s for _, s in vis_runs(c))


def _set_text(r, s):
    tn = q('w:t') if r.tag == q('w:r') else q('m:t')
    ts = r.findall(tn)
    assert len(ts) == 1, etree.tostring(r)[:300]
    ts[0].text = s
    if r.tag == q('w:r'):
        ts[0].set('{%s}space' % XML, 'preserve')


def split(r, a, b):
    """return the copy of run r holding characters [a, b) (siblings before/after keep the rest)"""
    s = rtext(r)
    if a == 0 and b == len(s):
        return r
    mid = copy.deepcopy(r); _set_text(mid, s[a:b]); r.addnext(mid)
    if b < len(s):
        aft = copy.deepcopy(r); _set_text(aft, s[b:]); mid.addnext(aft)
    if a > 0:
        _set_text(r, s[:a])
    else:
        r.getparent().remove(r)
    return mid


def wrap_del(cx, r):
    d = cx.mark('del'); par = r.getparent(); par.replace(r, d); d.append(r)
    if r.tag == q('w:r'):
        for t in r.findall(q('w:t')):
            t.tag = q('w:delText'); t.set('{%s}space' % XML, 'preserve')
    return d


def _after_point(w):
    """insertion point after wrapper w, leaving any w:ins of another author"""
    while w.getparent().tag == q('w:ins'):
        w = w.getparent()
    return w


def wrun(text, rpr=None, italic=False, sub=False, sup=False):
    r = etree.Element(q('w:r'))
    rp = copy.deepcopy(rpr) if rpr is not None else etree.Element(q('w:rPr'))
    for k in ('w:i', 'w:vertAlign'):
        for e in rp.findall(q(k)):
            rp.remove(e)
    if italic:
        rp.append(etree.Element(q('w:i')))
    if sub or sup:
        v = etree.Element(q('w:vertAlign')); v.set(q('w:val'), 'subscript' if sub else 'superscript'); rp.append(v)
    if len(rp):
        r.append(rp)
    t = etree.SubElement(r, q('w:t')); t.text = text; t.set('{%s}space' % XML, 'preserve')
    return r


def markup_runs(s, rpr=None):
    """plain text with $a_b$ (italic a, italic subscript b) and $a$ (italic) and ~x~ (upright subscript)"""
    out = []
    for k, part in enumerate(re.split(r'(\$[^$]+\$)', s)):
        if not part:
            continue
        if part.startswith('$'):
            body = part[1:-1]
            if '_' in body:
                a, b = body.split('_', 1)
                out.append(wrun(a, rpr, italic=True)); out.append(wrun(b, rpr, italic=True, sub=True))
            else:
                out.append(wrun(body, rpr, italic=True))
        else:
            out.append(wrun(part, rpr))
    return out


def replace(cx, c, old, new, nth=0):
    """tracked replacement of the nth visible occurrence of `old` in container c by `new` (markup)"""
    full = vis_text(c); pos = -1
    for _ in range(nth + 1):
        pos = full.find(old, pos + 1)
        assert pos >= 0, (old, full[:300])
    return replace_span(cx, c, pos, pos + len(old), new)


def replace_span(cx, c, pos, end, new):
    runs = vis_runs(c)
    mids = []; o = 0
    for r, s in runs:
        a, b = max(pos, o), min(end, o + len(s))
        if a < b:
            mids.append((r, a - o, b - o))
        o += len(s)
    pieces = [split(r, a, b) for r, a, b in mids]
    first = pieces[0]
    rpr = first.find(q('w:rPr'))
    wraps = [wrap_del(cx, p) for p in pieces]
    if new:
        ins = cx.mark('ins')
        if first.tag == q('m:r'):
            nr = copy.deepcopy(first); _set_text(nr, new); ins.append(nr)
            wraps[0].addnext(ins)
        else:
            for nr in markup_runs(new, rpr):
                ins.append(nr)
            _after_point(wraps[-1]).addnext(ins)
    return pieces


def delete_all(cx, c):
    for r, _ in vis_runs(c):
        wrap_del(cx, r)
    for cp in c.iter(q('m:ctrlPr')):
        mark_ctrl(cx, cp, 'del')


def mark_ctrl(cx, cp, kind):
    if cp.find(q('w:del')) is not None or cp.find(q('w:ins')) is not None:
        return
    rp = cp.find(q('w:rPr'))
    e = cx.mark(kind)
    if rp is not None:
        cp.remove(rp); e.append(rp)
    else:
        e.append(etree.Element(q('w:rPr')))
    cp.append(e)


def para_mark(cx, p, kind):
    pPr = p.find(q('w:pPr'))
    if pPr is None:
        pPr = etree.Element(q('w:pPr')); p.insert(0, pPr)
    rPr = pPr.find(q('w:rPr'))
    if rPr is None:
        rPr = etree.SubElement(pPr, q('w:rPr'))
    rPr.insert(0, cx.mark(kind))


def first_rpr(c):
    for r in c.iter(q('w:r')):
        rp = r.find(q('w:rPr'))
        if rp is not None and rp.find(q('w:rStyle')) is None:
            return rp
        if rp is None:
            return None
    return None


def rewrite_paragraph(cx, p, texts):
    """delete the visible content of p and put texts[0] into it; texts[1:] become new inserted paragraphs after p"""
    rpr = first_rpr(p)
    delete_all(cx, p)
    ins = cx.mark('ins')
    for r in markup_runs(texts[0], rpr):
        ins.append(r)
    # put the new text after the last child that is not a comment reference/range end at the very end
    p.append(ins)
    prev = p
    for t in texts[1:]:
        np_ = etree.Element(q('w:p'))
        pPr = p.find(q('w:pPr'))
        if pPr is not None:
            pp = copy.deepcopy(pPr)
            for e in pp.findall(q('w:rPr')):
                pp.remove(e)
            np_.append(pp)
        para_mark(cx, np_, 'ins')
        i2 = cx.mark('ins')
        for r in markup_runs(t, rpr):
            i2.append(r)
        np_.append(i2)
        prev.addnext(np_); prev = np_
    return prev


def append_text(cx, p, text):
    """tracked insertion of text (markup) at the end of paragraph p"""
    ins = cx.mark('ins')
    for r in markup_runs(text, first_rpr(p)):
        ins.append(r)
    p.append(ins)
    return ins


def new_para_after(cx, ref, like, text=None):
    """new tracked paragraph after `ref` with the paragraph properties of `like`"""
    np_ = etree.Element(q('w:p'))
    pPr = like.find(q('w:pPr'))
    if pPr is not None:
        pp = copy.deepcopy(pPr)
        for e in pp.findall(q('w:rPr')):
            pp.remove(e)
        np_.append(pp)
    para_mark(cx, np_, 'ins')
    if text:
        ins = cx.mark('ins')
        for r in markup_runs(text, first_rpr(like)):
            ins.append(r)
        np_.append(ins)
    ref.addnext(np_)
    return np_
