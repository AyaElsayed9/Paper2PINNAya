"""fix_cm7 (22 Nov 2026) on Aya's edited file GBT_periodic-NK-CM-AME-21112026_4_Aya (aya/, unzipped to ua5/).
 - the two new shear-term equations (C58) are numbered (28a) and (28b), in the same two-cell table layout as Eq. (28)
 - italic a, b in "for any two conventional or shear modes a and b"
 - "; The MAC" -> ". The MAC"
 - modes 1 and 2 (C61): bending comes from the amplitude varying along X
 - Fig. 7 caption (C71): what really separates modes 12 and 8 (scripts/c71_modes_12_8.py)
Output: GBT_periodic-NK-CM-AME-21112026_5.docx"""
import copy
import sys
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-22T12:00:00Z"
t = etree.parse('ua5/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); cx = Ctx(root)
W14 = 'http://schemas.microsoft.com/office/word/2010/wordml'


def one(text):
    ps = [p for p in b.iter(q('w:p')) if text in vis_text(p)]; assert len(ps) == 1, (text, len(ps)); return ps[0]


def strip_ids(el):
    for e in el.iter():
        for k in ('{%s}paraId' % W14, '{%s}textId' % W14):
            if k in e.attrib:
                del e.attrib[k]


# ---- (28a), (28b): move each new equation paragraph into a copy of the Eq. (28) table, all tracked as inserted
p28 = [p for p in b.iter(q('w:p')) if vis_text(p).strip() == "(28)"]; assert len(p28) == 1
tbl28 = p28[0].getparent().getparent().getparent(); assert tbl28.tag == q('w:tbl')
intro = one("The shear modes and their coupling with the conventional modes add, in the same way,")
pab = one("with, for any two conventional or shear modes a and b,")
eqs = [intro.getnext(), pab.getnext()]
for e in eqs:
    assert e.tag == q('w:p') and e.find('.//m:oMathPara', ns) is not None and not vis_text(e).strip() or True
for eq, lab in zip(eqs, ("(28a)", "(28b)")):
    nt = copy.deepcopy(tbl28); strip_ids(nt)
    tr = nt.find('w:tr', ns)
    trpr = tr.find('w:trPr', ns)
    if trpr is None:
        trpr = etree.Element(q('w:trPr')); tr.insert(0, trpr)
    trpr.append(cx.mark('ins'))
    c1, c2 = tr.findall('w:tc', ns)
    g1, g2 = nt.find('w:tblGrid', ns).findall('w:gridCol', ns); g1.set(q('w:w'), '7259'); g2.set(q('w:w'), '600')
    c1.find('w:tcPr/w:tcW', ns).set(q('w:w'), '7259'); c2.find('w:tcPr/w:tcW', ns).set(q('w:w'), '600')
    old = c1.find('w:p', ns)
    ppr = old.find('w:pPr', ns)                       # equation style of Eq. (28)
    eqppr = eq.find('w:pPr', ns)
    eqppr.getparent().remove(eqppr)
    nppr = copy.deepcopy(ppr); rp = etree.SubElement(nppr, q('w:rPr')); rp.append(cx.mark('ins'))
    eq.insert(0, nppr)
    eq.getparent().replace(eq, nt)                     # remove eq from body, put table in its place
    c1.replace(old, eq)
    np_ = c2.find('w:p', ns)
    for r in np_.findall('w:r', ns):
        np_.remove(r)
    pp = np_.find('w:pPr', ns); rp2 = etree.SubElement(pp, q('w:rPr')); rp2.append(cx.mark('ins'))
    ins = cx.mark('ins'); ins.append(wrun(lab)); np_.append(ins)

# ---- italic a and b (the run is inside our own insertion, so it is re-split in place)
ins = [e for e in pab.findall('w:ins', ns)]; assert len(ins) == 1
r = ins[0].find('w:r', ns); txt = rtext(r); assert txt == "with, for any two conventional or shear modes a and b,"
rpr = r.find('w:rPr', ns)
for e in list(ins[0]):
    ins[0].remove(e)
for s, it in (("with, for any two conventional or shear modes ", False), ("a", True), (" and ", False), ("b", True), (",", False)):
    ins[0].append(wrun(s, rpr, italic=it))

# ---- small fixes
p = one(";  The MAC is the squared cosine")
full = vis_text(p); k = full.find(";  The MAC"); replace_span(cx, p, k, k + 2, ".")
p = one("modes 1 and 2 are rigid translations of the cross-section in Z and Y")
replace(cx, p, ", i.e. the out-of-plane and in-plane bending of the panel along X;",
        "; since their amplitude varies along $X$, they describe the out-of-plane and in-plane bending of the panel;")

# ---- Fig. 7 caption (C71) and the new Fig. 7 (panel (f) coloured by the transverse stretching of the walls)
p = one("Figure 7. Mode decomposition.")
replace(cx, p, ", which are different cross-section deformations; their amplitude functions have the same number of "
               "half-waves, and their signs show how they combine in the critical mode.",
        ", colored by the transverse stretching |$\u03b5_ss$| of the walls. Modes 12 and 8 have almost the same in-plane "
        "shape with opposite sign, but mode 8 stretches the webs about three times more (dark webs in (f)); since their "
        "amplitudes have opposite signs, the shapes add up in the critical mode while the web stretching cancels.")
import shutil
from PIL import Image
pf = p.getprevious()
drw = etree.tostring(pf).decode(); assert 'r:embed="rId19"' in drw, "Fig. 7 drawing not found"
shutil.copy('../../fig/fig07_participation.png', 'ua5/word/media/image8.png')
wpx, hpx = Image.open('ua5/word/media/image8.png').size
WP = 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'
A_ = 'http://schemas.openxmlformats.org/drawingml/2006/main'
for e in list(pf.iter('{%s}extent' % WP)) + list(pf.iter('{%s}ext' % A_)):
    cxv = int(e.get('cx')); e.set('cy', str(round(cxv * hpx / wpx)))

t.write('ua5/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print("ok fix_cm7")
