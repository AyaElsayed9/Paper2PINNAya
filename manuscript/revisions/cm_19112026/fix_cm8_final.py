"""fix_cm8 (22 Nov 2026) on Aya's final file (aya/..._5_1_Aya_final.docx, unzipped to ua6/), asked by Aya after the final check:
 - abstract: "satisfies it exactly" -> "satisfies them exactly" (it = the periodic conditions)            [tracked]
 - Fig. 2 caption: "grey" -> "gray" (American spelling throughout)                                       [tracked]
 - figure alt texts: colour/normalised/discretisation -> American spelling                              [alt text]
 - reply to C24 (19 walls): "Professor you are right" -> "Professor, you are right"
 - reply to C64b (MAC): "Yes Prof. I added" -> "I agree, Professor, its meaning was not explained. I added"
 - the reply to the Fig. 1 comment is threaded under it (commentsExtended paraIdParent)
Output: GBT_periodic-NK-CM-AME-21112026_5_2.docx"""
import re
import sys
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-22T14:00:00Z"
t = etree.parse('ua6/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); cx = Ctx(root)


def one(text):
    ps = [p for p in b.iter(q('w:p')) if text in vis_text(p)]; assert len(ps) == 1, (text, len(ps)); return ps[0]


p = one("then satisfies it exactly")
replace(cx, p, "satisfies it exactly", "satisfies them exactly")
p = one("(grey: undeformed)")
replace(cx, p, "grey", "gray")
WP = '{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}docPr'
for e in root.iter(WP):
    d = e.get('descr')
    if d:
        d2 = d.replace('colour', 'color').replace('normalised', 'normalized').replace('discretisation', 'discretization') \
              .replace('grey', 'gray')
        if d2 != d:
            e.set('descr', d2)
t.write('ua6/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)

# replies (text of two replies, as Aya asked)
pc = 'ua6/word/comments.xml'; x = open(pc, encoding='utf-8').read()
for old, new in (("Professor you are right that this was unclear.", "Professor, you are right that this was unclear."),
                 ("Yes Prof. I added one sentence in words:",
                  "I agree, Professor, its meaning was not explained. I added one sentence in words:")):
    assert x.count(old) == 1, old
    x = x.replace(old, new)
open(pc, 'w', encoding='utf-8').write(x)

# thread the Fig. 1 reply under the professor's comment (last paragraph of his comment)
pe = 'ua6/word/commentsExtended.xml'; y = open(pe, encoding='utf-8').read()
old = '<w15:commentEx w15:paraId="776261CC" w15:done="0"/>'
assert y.count(old) == 1
y = y.replace(old, '<w15:commentEx w15:paraId="776261CC" w15:paraIdParent="2FAB13C9" w15:done="0"/>')
open(pe, 'w', encoding='utf-8').write(y)
print("ok fix_cm8")
