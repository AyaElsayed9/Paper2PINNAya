import subprocess, sys, glob, re, copy
from lxml import etree
sys.path.insert(0,'.')
from trk import q, ns
SK=glob.glob('/root/.claude/skills/synced/*/docx')[0]
R={
7:"Done. The first paragraph now keeps only [1,2] for the importance of stiffened panels. The omega stringer, the periodic cell and the strip methods moved to a new second paragraph. The last two lines are rewritten: the periodic condition is now explained, and [9–12] are described as strip methods for prismatic plate assemblies.",
8:"Done. These two paragraphs are rewritten as four paragraphs (about 900 words) without any new reference; each group of references is now explained. All references are kept, so the numbering is unchanged. Ref. [33] is now cited for vibration (its topic) instead of shear deformation.",
10:"Done: split into short sentences.",
21:"Done: split into short sentences.",
27:"Done: the global axes are now X (along the stringers), Y (across the panel) and Z (normal to the skin), in the text, Eqs. (10) and (12), Table 1, Sections 2.2.2, 2.2.5, 2.3.4 and 3, and Figures 1, 3, 6 and 9. The local wall coordinates of GBT stay x, s, z (Silvestre–Camotim notation), with x parallel to X; this is now stated in Section 2.1.1. Not changed: the geometric tensors X_ik and the matrices Z_g of Section 2.2.2 (we can rename them if they now look confusing).",
28:"Done: the figure shows only (a), (b), (c); all descriptions are in the caption (the numbers of nodes, walls and segments moved to the caption).",
29:"Done: three lines. Layout only, no symbol changed (made directly, not as a tracked change).",
30:"Done: two lines. Layout only (made directly, not as a tracked change).",
32:"Done: two columns, left εxx, εss, γxs and right γxz, γsz. Layout only (made directly, not as a tracked change).",
33:"Done: three lines. Layout only (made directly, not as a tracked change).",
37:"Thank you, Navid, you are right that the sentence was too short. But here γ_k is zero along the whole segment, not only at the mid-point. We compute the warping only for the G and D modes, and these modes do not stretch any segment, so v_p is the same along the segment. The warping u_k is linear between the nodes. So γ_k = (u_q' − u_q)/b_p + v_p has one constant value on each segment, and the mid-point value is the value everywhere. I also checked it in the code: v_p at the two ends of a segment differs by only about 1e-15. I added this explanation after Eq. (9).",
66:"Done: U5 now comes before U6.",
81:"Done: rewritten.",
}
ids={}
for p,txt in R.items():
    out=subprocess.run(['python3',f'{SK}/scripts/comment.py','w2',txt,'--author','Aya Medhat Elsayed','--initials','AME','--parent',str(p)],capture_output=True,text=True).stdout
    ids[p]=int(re.search(r'id=(\d+)',out).group(1))
t=etree.parse('w2/word/document.xml'); root=t.getroot()
for p,r in ids.items():
    s=root.find(f'.//w:commentRangeStart[@w:id="{p}"]',ns); e=root.find(f'.//w:commentRangeEnd[@w:id="{p}"]',ns); ref=root.find(f'.//w:commentReference[@w:id="{p}"]',ns)
    assert s is not None and e is not None and ref is not None, p
    ns_=copy.deepcopy(s); ns_.set(q('w:id'),str(r)); s.addnext(ns_)
    ne=copy.deepcopy(e); ne.set(q('w:id'),str(r)); e.addnext(ne)
    run=ref.getparent(); nr=copy.deepcopy(run); nr.find('.//w:commentReference',ns).set(q('w:id'),str(r)); run.addnext(nr)
t.write('w2/word/document.xml',xml_declaration=True,encoding='UTF-8',standalone=True)
print(ids)
