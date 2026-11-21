"""Accept all tracked changes of a docx and drop the comments (lxml; the markup kinds used in our revision files)."""
import sys, zipfile, shutil, os, re
from lxml import etree
W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
q = lambda t: '{%s}%s' % ({'w': W, 'm': M}[t.split(':')[0]], t.split(':')[1])


def accept(root):
    body = root.find(q('w:body'))
    # deleted content
    for d in list(root.iter(q('w:del'))):
        par = d.getparent()
        if par is None:
            continue
        if par.tag in (q('w:rPr'),):                  # paragraph-mark deletion: handled below
            continue
        if par.tag == q('m:ctrlPr') or (par.getparent() is not None and par.getparent().tag == q('m:ctrlPr')):
            par.remove(d); continue
        par.remove(d)
    for e in list(root.iter(q('w:moveFrom'))):
        e.getparent().remove(e)
    # inserted content: unwrap
    for tag in ('w:ins', 'w:moveTo'):
        for e in list(root.iter(q(tag))):
            par = e.getparent()
            if par is None:
                continue
            if par.tag == q('w:rPr') or par.tag == q('m:ctrlPr'):
                if par.tag == q('m:ctrlPr'):           # <m:ctrlPr><w:ins><w:rPr/></w:ins></m:ctrlPr>
                    for ch in list(e):
                        e.addprevious(ch)
                par.remove(e); continue
            for ch in list(e):
                e.addprevious(ch)
            par.remove(e)
    # paragraphs whose mark is deleted: merge with the next paragraph
    for p in list(root.iter(q('w:p'))):
        pPr = p.find(q('w:pPr'))
        rPr = pPr.find(q('w:rPr')) if pPr is not None else None
        dm = rPr.find(q('w:del')) if rPr is not None else None
        if dm is None:
            continue
        rPr.remove(dm)
        nxt = p.getnext()
        while nxt is not None and nxt.tag != q('w:p'):
            nxt = nxt.getnext()
        if nxt is None:
            continue
        content = [c for c in p if c.tag != q('w:pPr')]
        npPr = nxt.find(q('w:pPr'))
        for c in reversed(content):
            if npPr is not None:
                npPr.addnext(c)
            else:
                nxt.insert(0, c)
        p.getparent().remove(p)
    for tag in ('w:rPrChange', 'w:pPrChange', 'w:sectPrChange', 'w:tblPrChange', 'w:trPrChange', 'w:tcPrChange'):
        for e in list(root.iter(q(tag))):
            e.getparent().remove(e)
    # comments
    for tag in ('w:commentRangeStart', 'w:commentRangeEnd'):
        for e in list(root.iter(q(tag))):
            e.getparent().remove(e)
    for e in list(root.iter(q('w:commentReference'))):
        r = e.getparent(); r.getparent().remove(r)
    # empty paragraphs left by deletions (no text, no math, no drawing) inside the body
    return root


def main(src, dst):
    tmp = dst + '.d'; shutil.rmtree(tmp, ignore_errors=True)
    with zipfile.ZipFile(src) as z:
        z.extractall(tmp)
    doc = os.path.join(tmp, 'word', 'document.xml')
    t = etree.parse(doc); accept(t.getroot())
    t.write(doc, xml_declaration=True, encoding='UTF-8', standalone=True)
    for f in ('comments.xml', 'commentsExtended.xml', 'commentsIds.xml', 'commentsExtensible.xml', 'people.xml'):
        pth = os.path.join(tmp, 'word', f)
        if os.path.exists(pth):
            os.remove(pth)
    rel = os.path.join(tmp, 'word', '_rels', 'document.xml.rels'); s = open(rel, encoding='utf-8').read()
    s = re.sub(r'<Relationship [^>]*Target="(comments|commentsExtended|commentsIds|commentsExtensible|people)\.xml"/>', '', s)
    open(rel, 'w', encoding='utf-8').write(s)
    ct = os.path.join(tmp, '[Content_Types].xml'); s = open(ct, encoding='utf-8').read()
    s = re.sub(r'<Override PartName="/word/(comments|commentsExtended|commentsIds|commentsExtensible|people)\.xml"[^>]*/>', '', s)
    open(ct, 'w', encoding='utf-8').write(s)
    if os.path.exists(dst):
        os.remove(dst)
    with zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED) as z:
        for base, _, files in os.walk(tmp):
            for f in files:
                full = os.path.join(base, f); z.write(full, os.path.relpath(full, tmp))
    shutil.rmtree(tmp)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
