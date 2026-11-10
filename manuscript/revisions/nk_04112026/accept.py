"""Clean reading copy: accept all tracked changes (own implementation; LibreOffice drops tracked math)."""
import sys
from lxml import etree
sys.path.insert(0,'.')
from trk import q
t=etree.parse(sys.argv[1]); root=t.getroot()
for tr in list(root.iter(q('w:tr'))):
    tp=tr.find(q('w:trPr'))
    if tp is not None and tp.find(q('w:del')) is not None: tr.getparent().remove(tr)
for e in list(root.iter(q('w:del'))):
    par=e.getparent()
    if par is None: continue
    if par.tag==q('m:ctrlPr') or par.tag==q('w:rPr') and par.getparent().tag==q('w:pPr'):
        # ctrl: unwrap rPr; paragraph mark deletions only in removed rows
        for c in list(e): e.addprevious(c)
        par.remove(e)
    else:
        for c in list(e.iter(q('w:commentRangeStart'),q('w:commentRangeEnd'),q('w:bookmarkStart'),q('w:bookmarkEnd'))):
            c.getparent().remove(c); e.addprevious(c)
        par.remove(e)
for e in list(root.iter(q('w:ins'))):
    par=e.getparent()
    if par.tag in (q('w:trPr'),) or (par.tag==q('w:rPr') and par.getparent().tag==q('w:pPr')):
        par.remove(e); continue
    for c in list(e): e.addprevious(c)
    par.remove(e)
for tag in ('w:rPrChange','w:pPrChange','w:tcPrChange','w:trPrChange','w:tblPrChange','w:sectPrChange','w:tblGridChange'):
    for e in list(root.iter(q(tag))): e.getparent().remove(e)
t.write(sys.argv[1],xml_declaration=True,encoding='UTF-8',standalone=True)
t=etree.parse(sys.argv[1]); root=t.getroot()
for om in list(root.iter(q('m:oMath'))):
    if not ''.join(x.text or '' for x in om.iter(q('m:t'))) and om.find('.//'+q('w:commentReference')) is None:
        om.getparent().remove(om)
t.write(sys.argv[1],xml_declaration=True,encoding='UTF-8',standalone=True)
