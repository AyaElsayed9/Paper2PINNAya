"""fix_cm6 (21 Nov 2026) on Aya's own edited file GBT_periodic-NK-CM-AME-21112026_3 (aya/): Section 2.2.3 made
clearer (C49/C50) with plain-language additions; no equation changed.  Output: GBT_periodic-NK-CM-AME-21112026_4.docx"""
import sys
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-21T16:00:00Z"
t = etree.parse('ua/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); cx = Ctx(root)


def one(text):
    ps = [p for p in b.iter(q('w:p')) if text in vis_text(p)]; assert len(ps) == 1, (text, len(ps)); return ps[0]


p = one("; with the constraints")
replace(cx, p, "; with the constraints", " (the shear strain energy of the mode, summed over all segments); with the constraints")
p = one(" The solution is unique, since the axial constraint")
replace(cx, p, " The solution is unique,", " In words: the warping of each mode is chosen such that its shear strains are as small as "
        "possible, under two conditions: it adds no axial force, since the axial extension is the separate mode $u_j$ = 1, and "
        "both side edges have the same warping, which is the periodic condition. The solution is unique,")
p = one("For a mode with non-zero circulation")
replace(cx, p, "For a mode with non-zero circulation", "For a mode whose in-plane displacements run around the closed cell "
        "in one direction, i.e. with non-zero circulation")
t.write('ua/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print("ok fix_cm6")
