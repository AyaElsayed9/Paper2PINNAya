"""Equation layout (comments 29, 30, 32, 33): Eqs. (1), (2), (4) one expression per line; Eq. (3) as a 2-column array.
Layout only, no symbol changed (not tracked; answered in the comment threads)."""
import copy
import sys
from lxml import etree
sys.path.insert(0, '.')
from trk import q, ns, rtext

F = sys.argv[1]
tree = etree.parse(F); body = tree.getroot().find('w:body', ns)
B = list(body)


def find_eq(num):
    for e in B:
        if e.tag == q('w:tbl') and ''.join(t.text or '' for t in e.iter(q('w:t'))).strip().endswith(f'({num})'):
            return e.find('.//m:oMath', ns)
    raise KeyError(num)


QUAD = ' '


def is_sep(c):
    return c.tag == q('m:r') and rtext(c).startswith(QUAD)


def strip_sep(c):
    t = c.find(q('m:t')); t.text = t.text.lstrip(QUAD)
    return bool(t.text)


def ctrl(om):
    cp = om.find('.//m:ctrlPr', ns)
    return copy.deepcopy(cp) if cp is not None else etree.Element(q('m:ctrlPr'))


def split_lines(children):
    lines = [[]]
    for c in children:
        if is_sep(c):
            lines.append([])
            if strip_sep(c):
                lines[-1].append(c)
        else:
            lines[-1].append(c)
    return lines


def align_eq(line):
    for c in line:                                   # alignment point at the first "="
        if c.tag == q('m:r') and rtext(c) == '=':
            rp = c.find(q('m:rPr'))
            if rp is None:
                rp = etree.Element(q('m:rPr')); c.insert(0, rp)
            rp.append(etree.Element(q('m:aln')))
            return


def to_eqarr(num, nlines):
    om = find_eq(num)
    kids = list(om)
    # trailing anchors (comment end / reference run without text) stay outside the array
    last = max(i for i, c in enumerate(kids) if c.tag in (q('m:r'), q('m:d'), q('m:sSub'), q('m:sSubSup'), q('m:f'))
               and ''.join(t.text or '' for t in c.iter(q('m:t'))))
    body_, tail = kids[:last + 1], kids[last + 1:]
    lead = []
    while body_ and body_[0].tag in (q('w:commentRangeStart'),):
        lead.append(body_.pop(0))
    lines = split_lines(body_)
    assert len(lines) == nlines, (num, len(lines))
    for c in kids:
        om.remove(c)
    ea = etree.SubElement(om, q('m:eqArr'))
    pr = etree.SubElement(ea, q("m:eqArrPr")); pr.append(copy.deepcopy(cp_template))
    for k, ln in enumerate(lines):
        align_eq(ln)
        e = etree.SubElement(ea, q('m:e'))
        if k == 0:
            for c in lead:
                e.append(c)
        for c in ln:
            e.append(c)
    for c in tail:
        om.append(c)


cp_template = copy.deepcopy(body.find('.//m:ctrlPr', ns))
for num, n in ((1, 3), (2, 2), (4, 3)):
    to_eqarr(num, n)

# Eq. (3): left column eps_xx, eps_ss, gamma_xs; right column gamma_xz, gamma_sz
om = find_eq(3); mm = om.find('m:m', ns)
r1, r2 = [mr.find('m:e', ns) for mr in mm.findall('m:mr', ns)]
A, Bq = split_lines(list(r1)); C, D, E = split_lines(list(r2))
mPr = mm.find('m:mPr', ns)
mcs = mPr.find('m:mcs', ns); mc = mcs.find('m:mc', ns)
mc.find('m:mcPr', ns).find('m:count', ns).set(q('m:val'), '2')
mc.find('m:mcPr', ns).find('m:mcJc', ns).set(q('m:val'), 'left')
for mr in mm.findall('m:mr', ns):
    mm.remove(mr)


def gap():
    r = etree.Element(q('m:r')); rp = etree.SubElement(r, q('m:rPr')); s = etree.SubElement(rp, q('m:sty')); s.set(q('m:val'), 'p')
    t = etree.SubElement(r, q('m:t')); t.text = QUAD * 2
    return r


for left, right in ((A, D), (Bq, E), (C, None)):
    mr = etree.SubElement(mm, q('m:mr'))
    e1 = etree.SubElement(mr, q('m:e'))
    for c in left:
        e1.append(c)
    e2 = etree.SubElement(mr, q('m:e'))
    if right:
        e2.append(gap())
        for c in right:
            e2.append(c)
    else:
        sp = gap(); sp.find(q('m:t')).text = QUAD; e2.append(sp)
tree.write(F, xml_declaration=True, encoding='UTF-8', standalone=True)
print('equations re-laid out')
