"""Build the manuscript: markdown parts -> numbered markdown -> (1) MDPI docx on the template, (2) LaTeX PDF for checking."""
import re, os, subprocess, shutil, zipfile, sys
from refs import REFS


ORDER = """writeProtection view zoom removePersonalInformation removeDateAndTime doNotDisplayPageBoundaries displayBackgroundShape
printPostScriptOverText printFractionalCharacterWidth printFormsData embedTrueTypeFonts embedSystemFonts saveSubsetFonts saveFormsData
mirrorMargins alignBordersAndEdges bordersDoNotSurroundHeader bordersDoNotSurroundFooter gutterAtTop hideSpellingErrors hideGrammaticalErrors
activeWritingStyle proofState formsDesign attachedTemplate linkStyles stylePaneFormatFilter stylePaneSortMethod documentType mailMerge
revisionView trackRevisions doNotTrackMoves doNotTrackFormatting documentProtection autoFormatOverride styleLockTheme styleLockQFSet
defaultTabStop autoHyphenation consecutiveHyphenLimit hyphenationZone doNotHyphenateCaps showEnvelope summaryLength clickAndTypeStyle
defaultTableStyle evenAndOddHeaders bookFoldRevPrinting bookFoldPrinting bookFoldPrintingSheets drawingGridHorizontalSpacing
drawingGridVerticalSpacing displayHorizontalDrawingGridEvery displayVerticalDrawingGridEvery doNotUseMarginsForDrawingGridOrigin
drawingGridHorizontalOrigin drawingGridVerticalOrigin doNotShadeFormData noPunctuationKerning characterSpacingControl printTwoOnOne
strictFirstAndLastChars noLineBreaksAfter noLineBreaksBefore savePreviewPicture doNotValidateAgainstSchema saveInvalidXml ignoreMixedContent
alwaysShowPlaceholderText doNotDemarcateInvalidXml saveXmlDataOnly useXSLTWhenSaving saveThroughXslt showXMLTags alwaysMergeEmptyNamespace
updateFields hdrShapeDefaults footnotePr endnotePr compat docVars rsids mathPr attachedSchema themeFontLang clrSchemeMapping
doNotIncludeSubdocsInStats doNotAutoCompressPictures forceUpgrade captions readModeInkLockDown smartTagType schemaLibrary shapeDefaults
doNotEmbedSmartTags decimalSymbol listSeparator""".split()
def sort_settings(b):
    from lxml import etree
    root = etree.fromstring(b)
    W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
    kids = [k for k in root if isinstance(k.tag, str)]
    def key(k):
        if k.tag.startswith(W) and k.tag[len(W):] in ORDER: return ORDER.index(k.tag[len(W):])
        if k.tag.startswith('{http://schemas.openxmlformats.org/officeDocument/2006/math}'): return ORDER.index('mathPr')
        return len(ORDER)
    for k in kids: root.remove(k)
    for k in sorted(kids, key=key): root.append(k)
    return etree.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)

H = os.path.dirname(os.path.abspath(__file__))
os.chdir(H)
OUT = os.path.join(H, 'out'); os.makedirs(OUT, exist_ok=True)
NAME = 'GBT_periodic_omega_manuscript'

src = '\n\n'.join(open(f, encoding='utf-8').read() for f in ('part1.md', 'part2.md', 'part3_filled.md', 'part4.md'))

# ---------- figures at their natural size (600 dpi): all figures keep 8 pt lettering ----------
from PIL import Image as _Im
def _fw(m):
    w_px, _h = _Im.open(os.path.join(H, m.group(1))).size
    return f"]({m.group(1)}){{width={min(17.0, w_px / 600 * 2.54):.2f}cm}}"
src = re.sub(r"\]\((fig/[^)]+\.png)\)\{width=[^}]*\}", _fw, src)

# ---------- headings: manual MDPI numbering ----------
lines, c1, c2, c3, heads = [], 0, 0, 0, []
in_math = False
for ln in src.split('\n'):
    if ln.strip() == '$$':
        in_math = not in_math
    m = re.match(r'^(#{1,3}) (.*)$', ln) if not in_math else None
    if m:
        lev, txt = len(m.group(1)), m.group(2).strip()
        if txt.endswith('{-}'):
            txt = txt[:-3].strip(); lines.append('#' * lev + ' ' + txt + ' {-}'); heads.append(('', txt)); continue
        if lev == 1: c1 += 1; c2 = c3 = 0; num = f'{c1}.'
        elif lev == 2: c2 += 1; c3 = 0; num = f'{c1}.{c2}.'
        else: c3 += 1; num = f'{c1}.{c2}.{c3}.'
        heads.append((num, txt))
        lines.append('#' * lev + ' ' + num + ' ' + txt + ' {-}')
    else:
        lines.append(ln)
src = '\n'.join(lines)

# ---------- equations ----------
labels = re.findall(r'\\eqn\{([^}]+)\}', src)
assert len(labels) == len(set(labels)), 'duplicate equation labels'
eqno = {l: i + 1 for i, l in enumerate(labels)}
def eqref(m):
    l = m.group(1)
    if l not in eqno: sys.exit(f'unknown equation label {l}')
    return f'({eqno[l]})'
src = re.sub(r'\(@([A-Za-z0-9_]+)\)', eqref, src)

# ---------- citations ----------
order = []
for grp in re.findall(r'\{\{([^}]+)\}\}', src):
    for k in [x.strip() for x in grp.split(',')]:
        if k not in REFS: sys.exit(f'unknown reference {k}')
        if k not in order: order.append(k)
cno = {k: i + 1 for i, k in enumerate(order)}
def compress(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        out.append(f'{nums[i]}–{nums[j]}' if j - i >= 2 else ','.join(str(n) for n in nums[i:j + 1]))
        i = j + 1
    return '[' + ','.join(out) + ']'
src = re.sub(r'\{\{([^}]+)\}\}', lambda m: compress([cno[k.strip()] for k in m.group(1).split(',')]), src)
reflist = [REFS[k] for k in order]

unused = [k for k in REFS if k not in cno]

def variant(kind):
    s = src
    if kind == 'docx':
        s = re.sub(r'\\eqn\{([^}]+)\}', lambda m: r'\text{EQNUM%d}' % eqno[m.group(1)], s)
        s += '\n\n# References {-}\n\n' + '\n\n'.join(
            '::: {custom-style="MDPI_8.1_references"}\n' + r + '\n:::' for r in reflist) + '\n'
    else:
        s = re.sub(r'\\eqn\{([^}]+)\}', lambda m: r'\tag{%d}' % eqno[m.group(1)], s)
        s += '\n\n# References {-}\n\n' + '\n'.join(f'{i+1}. {r}' for i, r in enumerate(reflist)) + '\n'
    return s

open(os.path.join(OUT, NAME + '_docx.md'), 'w', encoding='utf-8').write(variant('docx'))
open(os.path.join(OUT, NAME + '_tex.md'), 'w', encoding='utf-8').write(variant('tex'))

# ---------- docx ----------
raw = os.path.join(OUT, 'raw.docx')
subprocess.run(['pandoc', NAME + '_docx.md', '-f', 'markdown+tex_math_dollars', '-o', 'raw.docx',
                '--reference-doc=' + os.path.join(H, 'src/template.docx'), '--resource-path=' + H], cwd=OUT, check=True)

z = zipfile.ZipFile(raw)
files = {n: z.read(n) for n in z.namelist()}
doc = files['word/document.xml'].decode('utf-8')

stylemap = {'BodyText': 'MDPI31text', 'FirstParagraph': 'MDPI31text', 'Compact': 'MDPI42tablebody',
            'Heading1': 'MDPI21heading1', 'Heading2': 'MDPI22heading2', 'Heading3': 'MDPI23heading3',
            'ImageCaption': 'MDPI51figurecaption', 'CaptionedFigure': 'MDPI52figure', 'Figure': 'MDPI52figure',
            'TableCaption': 'MDPI41tablecaption'}
for a, b in stylemap.items():
    doc = doc.replace(f'<w:pStyle w:val="{a}" />', f'<w:pStyle w:val="{b}" />').replace(f'<w:pStyle w:val="{a}"/>', f'<w:pStyle w:val="{b}"/>')

# list paragraphs: MDPI text indent, pStyle first in pPr
doc = re.sub(r'<w:pPr><w:numPr>(.*?)</w:numPr><w:pStyle w:val="[^"]*" ?/></w:pPr>',
             r'<w:pPr><w:pStyle w:val="MDPI31text"/><w:numPr>\1</w:numPr><w:ind w:left="3033" w:hanging="425"/><w:jc w:val="both"/></w:pPr>', doc)
# paragraphs inside tables: table body style
def fix_table(m):
    t = m.group(0)
    if 'EQNUM' in t: return t
    t = re.sub(r'<w:pStyle w:val="MDPI31text" ?/>', '<w:pStyle w:val="MDPI42tablebody"/>', t)
    t = re.sub(r'<w:tblStyle w:val="[^"]*" ?/>', '<w:tblStyle w:val="MDPI41threelinetable"/>', t)
    ncol = len(re.findall(r'<w:gridCol ', t))
    wide = ncol > 8
    t = re.sub(r'<w:tblW [^>]*/>', '<w:tblW w:w="%d" w:type="dxa"/><w:jc w:val="%s"/>' % ((10467, 'center') if wide else (7859, 'right')), t, count=1)
    t = t.replace('<w:tblLayout w:type="fixed" />', '').replace('<w:tblLayout w:type="fixed"/>', '')
    # header row: bottom rule
    first = re.search(r'<w:tr>.*?</w:tr>|<w:tr .*?</w:tr>', t, re.S)
    if first:
        r = first.group(0)
        B = '<w:tcBorders><w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/></w:tcBorders>'
        r2 = r.replace('<w:tcPr />', '<w:tcPr>' + B + '</w:tcPr>').replace('<w:tcPr/>', '<w:tcPr>' + B + '</w:tcPr>')
        r2 = re.sub(r'(<w:tcPr>)(<w:tcW [^>]*/>)?(?!<w:tcBorders>)', lambda m: m.group(1) + (m.group(2) or '') + B, r2)
        r2 = r2.replace(B + B, B)
        t = t.replace(r, r2, 1)
    if wide:  # small font for the wide appendix tables
        t = re.sub(r'(<w:pStyle w:val="MDPI42tablebody"/>)', r'\1', t)
        t = re.sub(r'<w:r>', '<w:r><w:rPr><w:sz w:val="14"/></w:rPr>', t)
    return t
doc = re.sub(r'<w:tbl>.*?</w:tbl>', fix_table, doc, flags=re.S)

# equations: wrap display equations with a number into the template's equation table
EQTBL = ('<w:tbl><w:tblPr><w:tblW w:w="7859" w:type="dxa"/><w:tblInd w:w="2608" w:type="dxa"/><w:tblCellMar>'
         '<w:left w:w="0" w:type="dxa"/><w:right w:w="0" w:type="dxa"/></w:tblCellMar><w:tblLook w:val="04A0" '
         'w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/></w:tblPr>'
         '<w:tblGrid><w:gridCol w:w="7428"/><w:gridCol w:w="431"/></w:tblGrid><w:tr><w:tc><w:tcPr><w:tcW w:w="7428" '
         'w:type="dxa"/></w:tcPr><w:p><w:pPr><w:pStyle w:val="MDPI39equation"/><w:ind w:left="0"/><w:jc w:val="center"/></w:pPr>{math}</w:p>'
         '</w:tc><w:tc><w:tcPr><w:tcW w:w="431" w:type="dxa"/><w:vAlign w:val="center"/></w:tcPr><w:p><w:pPr>'
         '<w:pStyle w:val="MDPI3aequationnumber"/></w:pPr><w:r><w:t>({n})</w:t></w:r></w:p></w:tc></w:tr></w:tbl>')
neq = 0
def fix_eq(m):
    global neq
    p = m.group(0)
    k = re.search(r'EQNUM(\d+)', p)
    if not k: return p
    n = k.group(1)
    p = re.sub(r'<m:r>(?:(?!</m:r>).)*?EQNUM%s(?:(?!</m:r>).)*?</m:r>' % n, '', p, flags=re.S)
    math = re.search(r'<m:oMathPara>.*</m:oMathPara>', p, re.S)
    if not math:
        sys.exit('equation %s not in its own paragraph' % n)
    rest = re.sub(r'<m:oMathPara>.*</m:oMathPara>', '', p, flags=re.S)
    if re.search(r'<w:t[ >]', rest):
        sys.exit('equation %s shares its paragraph with text' % n)
    neq += 1
    # left-aligned math inside the cell
    mth = math.group(0).replace('<m:oMathParaPr><m:jc m:val="center" /></m:oMathParaPr>', '')
    mth = mth.replace('<m:oMathPara>', '<m:oMathPara><m:oMathParaPr><m:jc m:val="center"/></m:oMathParaPr>', 1)
    return EQTBL.format(math=mth, n=n)
doc = re.sub(r'<w:p>(?:(?!</w:p>).)*?EQNUM\d+.*?</w:p>', fix_eq, doc, flags=re.S)
assert 'EQNUM' not in doc, 'unprocessed equation marker'
assert neq == len(labels), (neq, len(labels))

doc = re.sub(r'(<m:sty m:val="[a-z]+" ?/>)(<m:scr m:val="[^"]+" ?/>)', r'\2\1', doc)
doc = doc.replace('<m:nor /><m:sty m:val="p" />', '<m:nor />')
doc = re.sub(r'(<m:begChr [^>]*/>)(<m:endChr [^>]*/>)(<m:sepChr [^>]*/>)', r'\1\3\2', doc)
doc = re.sub(r'<w:jc w:val="start" ?/></w:tblPr>', '</w:tblPr>', doc)
doc = re.sub(r'(<m:mcJc [^>]*/>)(<m:count [^>]*/>)', r'\2\1', doc)
doc = re.sub(r'<w:bidi ?/>', '', doc)   # template artefact: RTL section
files['word/document.xml'] = doc.encode('utf-8')
num = files['word/numbering.xml'].decode('utf-8')
num = re.sub(r'<w:nsid w:val="([0-9A-Fa-f]+)" ?/>', lambda m: '<w:nsid w:val="%s"/>' % m.group(1).upper().zfill(8)[-8:], num)
files['word/numbering.xml'] = num.encode('utf-8')
files['word/settings.xml'] = sort_settings(files['word/settings.xml'])
ct = files['[Content_Types].xml'].decode('utf-8')
for ext, typ in (('png', 'image/png'), ('jpeg', 'image/jpeg'), ('jpg', 'image/jpeg'), ('emf', 'image/x-emf')):
    if 'Extension="%s"' % ext not in ct:
        ct = ct.replace('<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="%s" ContentType="%s"/>' % (ext, typ), 1)
files['[Content_Types].xml'] = ct.encode('utf-8')
refd = set()
for n, b in files.items():
    if n.endswith('.rels'):
        for t in re.findall(r'Target="([^"]+)"', b.decode('utf-8')):
            refd.add(os.path.basename(t))
for n in [n for n in files if n.startswith('word/media/') and os.path.basename(n) not in refd]:
    del files[n]
dst = os.path.join(OUT, NAME + '.docx')
with zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED) as zo:
    for n, b in files.items(): zo.writestr(n, b)
os.remove(raw)

# ---------- LaTeX PDF (for checking the mathematics) ----------
hdr = os.path.join(OUT, 'hdr.tex')
open(hdr, 'w').write(r'''\usepackage{amsmath,amssymb}\usepackage{booktabs}\usepackage{caption}\captionsetup{labelformat=empty}


\renewcommand{\arraystretch}{1.1}\setlength{\emergencystretch}{3em}
''')
r = subprocess.run(['pandoc', NAME + '_tex.md', '-f', 'markdown+tex_math_dollars', '-o', NAME + '_preview.pdf',
                    '--pdf-engine=xelatex', '-H', 'hdr.tex', '-V', 'fontsize=10pt', '-V', 'fontfamily=amsmath', '-V', 'mainfont=TeX Gyre Pagella', '-V', 'mathfont=Latin Modern Math', '-V', 'geometry:a4paper,margin=2.2cm', '--resource-path=' + H],
                   cwd=OUT, capture_output=True, text=True)
print(r.stderr[-3000:])
print('headings:'); [print('  ', n, t) for n, t in heads]
print('equations', len(labels), 'references', len(reflist), 'unused refs', unused)
