"""Thin-Walled Structures (Elsevier) version of the periodic-GBT paper (21 Nov 2026, Aya: target journal TWS).

Input : the revised manuscript with all tracked changes accepted (accept.py on GBT_periodic-NK-CM-AME-21112026.docx).
Output: GBT_periodic_TWS.docx (manuscript) and GBT_periodic_TWS_highlights.docx (separate highlights file).

TWS guide for authors (checked 21 Nov 2026): abstract <= 250 words without references; 3-5 highlights of <= 85
characters each, as a separate file; 1-7 keywords; references numbered [n] in the order of citation, journal names
abbreviated (LTWA); declarations of competing interest and of generative-AI use (prescribed wording); CRediT
statement and data availability.  The scientific content is not changed."""
import copy, re, sys, os, zipfile, shutil
from lxml import etree

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
R_ = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
q = lambda t: '{%s}%s' % ({'w': W, 'm': M, 'r': R_}[t.split(':')[0]], t.split(':')[1])
XS = '{http://www.w3.org/XML/1998/namespace}space'

SRC, OUT = sys.argv[1], sys.argv[2]
tmp = OUT + '.d'; shutil.rmtree(tmp, ignore_errors=True)
zipfile.ZipFile(SRC).extractall(tmp)
P = lambda *a: os.path.join(tmp, *a)
doc = etree.parse(P('word', 'document.xml')); root = doc.getroot(); body = root.find(q('w:body'))
sty = etree.parse(P('word', 'styles.xml'))


def text(p):
    return ''.join(t.text or '' for t in p.iter(q('w:t')))


def style(p):
    s = p.find(q('w:pPr') + '/' + q('w:pStyle'))
    return s.get(q('w:val')) if s is not None else None


def run(t, b=False, i=False, sup=False, size=None):
    r = etree.Element(q('w:r')); rp = etree.SubElement(r, q('w:rPr'))
    if b: etree.SubElement(rp, q('w:b'))
    if i: etree.SubElement(rp, q('w:i'))
    if size: etree.SubElement(rp, q('w:sz')).set(q('w:val'), str(size))
    if sup: etree.SubElement(rp, q('w:vertAlign')).set(q('w:val'), 'superscript')
    tt = etree.SubElement(r, q('w:t')); tt.text = t; tt.set(XS, 'preserve')
    return r


def para(sty_id, runs, jc=None, first_line=None):
    p = etree.Element(q('w:p')); pp = etree.SubElement(p, q('w:pPr'))
    etree.SubElement(pp, q('w:pStyle')).set(q('w:val'), sty_id)
    if first_line is not None:
        ind = etree.SubElement(pp, q('w:ind')); ind.set(q('w:firstLine'), str(first_line))
    if jc:
        etree.SubElement(pp, q('w:jc')).set(q('w:val'), jc)
    for r in runs:
        p.append(r)
    return p


def set_runs(p, runs):
    for ch in list(p):
        if ch.tag != q('w:pPr'):
            p.remove(ch)
    for r in runs:
        p.append(r)


# ---------------------------------------------------------------- page layout and styles
for st in sty.getroot().findall(q('w:style')):
    ind = st.find(q('w:pPr') + '/' + q('w:ind'))
    if ind is not None and ind.get(q('w:left')) in ('2608', '2806'):
        ind.set(q('w:left'), '0')
        if ind.get(q('w:hanging')):
            del ind.attrib[q('w:hanging')]
for ti in root.iter(q('w:tblInd')):
    if ti.get(q('w:w')) == '2608':
        ti.set(q('w:w'), '0')
sp = body.find(q('w:sectPr'))
for e in list(sp):
    if e.tag in (q('w:headerReference'), q('w:footerReference'), q('w:titlePg'), q('w:bidi'), q('w:docGrid')):
        sp.remove(e)
mar = sp.find(q('w:pgMar'))
for k, v in (('top', '1417'), ('bottom', '1417'), ('left', '1417'), ('right', '1417')):
    mar.set(q('w:' + k), v)
# footer with the page number
open(P('word', 'footer1.xml'), 'w', encoding='utf-8').write(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    f'<w:ftr xmlns:w="{W}"><w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r>'
    '</w:p></w:ftr>')
rels = open(P('word', '_rels', 'document.xml.rels'), encoding='utf-8').read()
rid = re.search(r'Id="(rId\d+)"[^>]*Target="footer1\.xml"', rels) or re.search(r'Target="footer1\.xml"[^>]*Id="(rId\d+)"', rels)
fr = etree.Element(q('w:footerReference')); fr.set(q('w:type'), 'default'); fr.set(q('r:id'), rid.group(1))
sp.insert(0, fr)

# ---------------------------------------------------------------- title page
kids = list(body)
for el in kids:
    if el.tag == q('w:tbl') and 'Academic Editor' in text(el):
        body.remove(el)
for p in list(body.iter(q('w:p'))):
    s = style(p)
    if s in ('MDPI11articletype', 'MDPI63notes', 'MDPI14history', 'MDPI15academiceditor', 'MDPI61citation', 'MDPI72copyright'):
        p.getparent().remove(p)
ps = list(body.iter(q('w:p')))
pa = next(p for p in ps if style(p) == 'MDPI13authornames')
set_runs(pa, [run('Aya Medhat Elsayed'), run('a', sup=True), run(', Navid Kharghani'), run('a,', sup=True), run('*', sup=True),
              run(', Christian Mittelstedt'), run('a', sup=True)])
aff = [p for p in ps if style(p) == 'MDPI16affiliation']
set_runs(aff[0], [run('a', sup=True), run(' Institute for Lightweight Engineering and Structural Mechanics (LSM), Department of '
                                           'Mechanical Engineering, Technical University of Darmstadt, Otto-Berndt-Straße 2, '
                                           '64287 Darmstadt, Germany')])
set_runs(aff[1], [run('*', sup=True), run(' Corresponding author. E-mail address: navid.kharghani@lsm.tu-darmstadt.de '
                                           '(N. Kharghani).')])
set_runs(aff[2], [run('E-mail addresses: ayaallah.elsayed@stud.tu-darmstadt.de (A.M. Elsayed), '
                      'christian.mittelstedt@lsm.tu-darmstadt.de (C. Mittelstedt).')])
for p in aff[3:]:
    if not text(p).strip():
        p.getparent().remove(p)

# ---------------------------------------------------------------- highlights and keywords
HL = ["Periodic boundary conditions are built into the cross-section analysis of GBT",
      "Every mode, and so every truncated mode basis, is exactly periodic",
      "Shear-deformable GBT–Ritz model with laminated walls and the true stringer offset",
      "Critical loads of 30 omega-stringer cells agree with Abaqus within 0.6% on average",
      "The constant width of periodic cells causes a long-wave mode in long panels"]
assert all(len(h) <= 85 for h in HL), [len(h) for h in HL]
ab = next(p for p in body.iter(q('w:p')) if style(p) == 'MDPI17abstract' and text(p).strip() == 'Abstract')
abs_p = ab.getnext()
nw = len(text(abs_p).split()); assert nw <= 250, nw
hpos = ab
hp = para('MDPI17abstract', [run('Highlights', b=True)]); ab.addprevious(hp)
for h in HL:
    x = para('MDPI17abstract', [run('• ' + h)]); ab.addprevious(x)
kw = next(p for p in body.iter(q('w:p')) if style(p) == 'MDPI18keywords')
set_runs(kw, [run('Keywords: ', b=True), run('Generalized beam theory; Periodic boundary conditions; Stiffened panels; '
                                             'Omega stringer; Laminated composites; Buckling; Ritz method')])

# ---------------------------------------------------------------- figures: "Figure 1" -> "Fig. 1" (Elsevier style)
nfig = 0
for p in body.iter(q('w:p')):
    for tt in p.iter(q('w:t')):
        if not tt.text or 'Figure' not in tt.text:
            continue
        s0 = tt.text
        s1 = re.sub(r'Figures(\s+\d)', r'Figs.\1', s0)
        s1 = re.sub(r'Figure(\s+\d)', r'Fig.\1', s1)
        if s1 != s0:
            nfig += s0.count('Figure') - s1.count('Figure'); tt.text = s1
# a sentence must not start with an abbreviation: restore "Figure"/"Figures" at the start of a body sentence
for p in body.iter(q('w:p')):
    if style(p) not in ('MDPI31text', 'MDPI32textnoindent', 'MDPI33textspaceafter', 'MDPI35textbeforelist', 'MDPI36textafterlist'):
        continue
    ts = list(p.iter(q('w:t'))); full = ''.join(t.text or '' for t in ts); o = 0
    for t in ts:
        s0 = t.text or ''
        for m_ in reversed(list(re.finditer(r'Figs?\.(?=\s+\d)', s0))):
            before = (full[:o + m_.start()]).rstrip()
            if before == '' or before.endswith(('.', ':', '?', '!')):
                rep = 'Figures' if m_.group(0) == 'Figs.' else 'Figure'
                s0 = s0[:m_.start()] + rep + s0[m_.end():]
        t.text = s0; o += len(s0)

# ---------------------------------------------------------------- back matter (Elsevier declarations)
bm = [p for p in body.iter(q('w:p')) if style(p) == 'MDPI62backmatter']
anchor = bm[0]
BACK = [
 ("CRediT authorship contribution statement",
  "Aya Medhat Elsayed: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Writing – "
  "original draft. Navid Kharghani: Conceptualization, Methodology, Validation, Writing – review & editing, Supervision. "
  "Christian Mittelstedt: Writing – review & editing, Project administration, Funding acquisition."),
 ("Declaration of competing interest",
  "The authors declare the following financial interests which may be considered as potential competing interests: "
  "financial support of the German Research Council (DFG – Deutsche Forschungsgemeinschaft), project number 510739591."),
 ("Acknowledgements",
  "This work was funded by the German Research Council (DFG – Deutsche Forschungsgemeinschaft), project number 510739591."),
 ("Data availability",
  "The raw/processed data required to reproduce these findings cannot be shared at this time as the data also forms part "
  "of an ongoing study."),
 ("Declaration of generative AI and AI-assisted technologies in the manuscript preparation process",
  "During the preparation of this work the authors used Claude Opus 5.5 (Anthropic) in order to support the literature "
  "research and to refine the language of the manuscript. After using this tool, the authors reviewed and edited the "
  "content as needed and take full responsibility for the content of the published article.")]
for h, t in BACK:
    anchor.addprevious(para('MDPI21heading1', [run(h)]))
    anchor.addprevious(para('MDPI31text', [run(t)], first_line=0))
for p in bm:
    p.getparent().remove(p)

# ---------------------------------------------------------------- references: Elsevier style, numbered [n]
ABBR = {"Thin-Walled Structures": "Thin-Walled Struct.", "Composite Structures": "Compos. Struct.",
        "Results in Engineering": "Results Eng.", "Journal of Sound and Vibration": "J. Sound Vib.",
        "International Journal of Mechanical Sciences": "Int. J. Mech. Sci.", "The Aeronautical Journal": "Aeronaut. J.",
        "Journal of Constructional Steel Research": "J. Constr. Steel Res.", "Applied Mechanics Reviews": "Appl. Mech. Rev.",
        "Journal of Engineering Mechanics": "J. Eng. Mech.", "Computational Mechanics": "Comput. Mech.",
        "International Journal of Structural Stability and Dynamics": "Int. J. Struct. Stab. Dyn.",
        "International Journal of Steel Structures": "Int. J. Steel Struct.", "Journal of Composites Science": "J. Compos. Sci.",
        "ce/papers": "ce/papers", "Journal of Applied Mechanics": "J. Appl. Mech.",
        "Computer Methods in Applied Mechanics and Engineering": "Comput. Methods Appl. Mech. Eng."}


def authors(s):
    out = []
    for a in [x.strip() for x in s.split(';') if x.strip()]:
        if ',' in a:
            last, ini = [y.strip() for y in a.split(',', 1)]; out.append(f"{ini} {last}")
        else:
            out.append(a)
    return ', '.join(out)


def convert(s):
    s = s.replace('\xa0', ' ').strip(); doi = ''
    m = re.search(r'\s*(https://doi\.org/\S+?)\.?$', s)
    if m:
        doi = m.group(1); s = s[:m.start()].rstrip()
    m = re.match(r'^(.+?\.)\s(?=[A-Z0-9‘“])', s)          # authors end at the first ". " followed by a capital
    # the author block: "Last, I.; Last, I.J." -> split at the first ". " after the last initial
    am = re.match(r'^((?:[^;]+?,\s(?:[A-Z]\.(?:-?[A-Z]\.)*)(?:;\s)?)+)\s', s)
    au = authors(am.group(1)); rest = s[am.end():]
    tail = (' ' + doi + '.') if doi else ''
    # chapter in a book
    m = re.match(r'^(.+?)\. In (.+?); (.+?): (.+?), (\d{4}); pp\. ([\d–-]+)\.?$', rest)
    if m:
        return f"{au}, {m.group(5)}. {m.group(1)}, in: {m.group(2)}. {m.group(3)}, {m.group(4)}, pp. {m.group(6)}.{tail}"
    # proceedings
    m = re.match(r'^(.+?)\. In (Proceedings of .+?), (\d{4})\.?$', rest)
    if m:
        return f"{au}, {m.group(3)}. {m.group(1)}, in: {m.group(2)}.{tail}"
    # book
    m = re.match(r'^(.+?); (.+?): (.+?), (\d{4})\.?(.*)$', rest)
    if m:
        extra = m.group(5).strip()
        return f"{au}, {m.group(4)}. {m.group(1).rstrip('.')}. {m.group(2)}, {m.group(3)}.{(' ' + extra) if extra else ''}{tail}"
    # journal article
    m = re.match(r'^(.+?)\. (.+?) (\d{4}), (.+?), ([\w–-]+)\.?$', rest)
    if m:
        j = ABBR.get(m.group(2), m.group(2))
        return f"{au}, {m.group(3)}. {m.group(1)}. {j} {m.group(4)}, {m.group(5)}.{tail}"
    raise ValueError(s)


refs = [p for p in body.iter(q('w:p')) if style(p) == 'MDPI81references']
for st in sty.getroot().findall(q('w:style')):              # references left-aligned (no stretched DOIs)
    if st.get(q('w:styleId')) == 'MDPI81references':
        jc = st.find(q('w:pPr') + '/' + q('w:jc'))
        if jc is not None:
            jc.set(q('w:val'), 'left')
assert len(refs) == 57, len(refs)
for p in refs:
    set_runs(p, [run(convert(text(p)))])
# list label "[n]" instead of "n."
num = etree.parse(P('word', 'numbering.xml')); nr = num.getroot()
st_ref = next(s for s in sty.getroot().findall(q('w:style')) if s.get(q('w:styleId')) == 'MDPI81references')
nid = st_ref.find('.//' + q('w:numId')).get(q('w:val'))
absid = next(n.find(q('w:abstractNumId')).get(q('w:val')) for n in nr.findall(q('w:num')) if n.get(q('w:numId')) == nid)
an = next(a for a in nr.findall(q('w:abstractNum')) if a.get(q('w:abstractNumId')) == absid)
lvl = an.find(q('w:lvl')); lvl.find(q('w:lvlText')).set(q('w:val'), '[%1]')
num.write(P('word', 'numbering.xml'), xml_declaration=True, encoding='UTF-8', standalone=True)

doc.write(P('word', 'document.xml'), xml_declaration=True, encoding='UTF-8', standalone=True)
sty.write(P('word', 'styles.xml'), xml_declaration=True, encoding='UTF-8', standalone=True)
if os.path.exists(OUT):
    os.remove(OUT)
with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as z:
    for base, _, files in os.walk(tmp):
        for f in files:
            full = os.path.join(base, f); z.write(full, os.path.relpath(full, tmp))
shutil.rmtree(tmp)
print(f"ok: abstract {nw} words, {len(HL)} highlights, {nfig} figure references abbreviated, {len(refs)} references")
hl_md = os.path.join(os.path.dirname(os.path.abspath(OUT)), 'highlights.md')
open(hl_md, 'w').write('# Highlights\n\n' + ''.join(f'- {h}\n' for h in HL))
import subprocess
subprocess.run(['pandoc', hl_md, '-o', OUT.replace('.docx', '_highlights.docx')], check=True)
