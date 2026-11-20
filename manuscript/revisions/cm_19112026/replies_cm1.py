"""replies_cm1 (20 Nov 2026): replies to Prof. Mittelstedt's comments, as threaded replies under Aya's name
(docx skill comment.py --parent), markers placed next to the parent's markers."""
import subprocess, re, sys
from lxml import etree
D = "/root/.claude/skills/synced/89684747-8288-4727-bcdd-0f93eb064259_c99e516e-a235-41ea-88e1-b79caa917733/docx/scripts"
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"; q = lambda t: "{%s}%s" % (W, t)
R = {
 22: "Agreed. We will add closed-form and semi-analytical work of other groups (e.g. Bisagni & Vescovini, TWS 2009; "
     "Stamatelos, Labeas & Tserpes, TWS 2011; Buermann et al., TWS 2006; Vescovini & Bisagni, Aeronaut. J. 2015) to this paragraph.",
 31: "Thank you - the four short link walls (foot mid-plane to web at the true offset) are hardly visible in Fig. 1c. "
     "The wall types are now listed in the text (5 skin, 4 foot, 4 web, 2 crown, 4 link), and the segments and nodes are explained.",
 32: "Agreed; Fig. 1a will show the origin of the X, Y, Z system and Fig. 1b the pitch b_s between the stringer axes.",
 37: "Defined now: each warping function is 1 at one node, 0 at all other nodes and linear along the segments.",
 41: "Agreed; the label is removed and the statement kept as plain text.",
 46: "N is the number of nodes (Section 2.1.1); this is now repeated here.",
 47: "We will rewrite this step in plain words first (keep the modes that are already periodic; combine the others within "
     "their family so that the side-edge differences cancel; add the few missing periodic deformations), followed by the equations.",
 49: "Agreed. Section 2.2.3 will start with the physical idea (joining the side edges closes the contour, so the warping "
     "can contain one extra constant shear flow), and each equation will be followed by one sentence explaining it.",
 50: "No - it is the difference of the warping between the two end nodes q and q' of segment p, which is non-zero in "
     "general; only the two side edges have equal warping. We added a sentence to make this clear.",
 53: "Agreed; the paragraph will start with the plain idea: in an infinitely wide panel a mode repeats from bay to bay up "
     "to a phase shift, and the two-pitch cell contains the two extreme cases (repeat in every bay, alternate from bay to bay).",
 54: "Both are correct: in the cross-section, modes 1 and 2 are rigid translations; with their amplitude functions along X "
     "they are the out-of-plane and in-plane bending of the panel. The sentence now says both.",
 58: "Agreed; the shear-mode terms will be written out explicitly.",
 61: "Yes, they are; the point is the continuity between elements: cubic Hermite elements have a continuous slope but jumps "
     "of the curvature (C1), whereas cubic B-splines are C2 across the knots. The sentence now says this.",
 64: "Thank you. The magnitude is only a reference (linear analysis, lambda = Ncr); the distribution follows [13,18]. We "
     "checked the uniform-end-strain distribution (line loads proportional to the axial stiffness of each wall) for 32 cells: "
     "the critical mean force per unit width changes by 0.3-0.5% (Cross-ply, Unsym) and 1.8-4.3% (Sym1, Sym2). This is now "
     "stated in the text.",
 65: "Explained now: MAC = 1 for identical shapes (up to a factor), MAC = 0 for orthogonal shapes.",
 67: "They are the largest absolute differences between the values at the two side edges over all modes of each type; "
     "the sentence now says this.",
 71: "No - Fig. 7(f) shows that they are different cross-section deformations. Their amplitude functions look alike because "
     "all components of one buckling mode have the same number of half-waves; the signs show how they combine. The caption now says this.",
}
for pid, text in R.items():
    out = subprocess.run([sys.executable, f"{D}/comment.py", "u", text, "--author", "Aya Medhat Elsayed", "--initials", "AE",
                          "--parent", str(pid)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    cm = etree.parse("u/word/comments.xml").getroot()
    nid = str(max(int(c.get(q("id"))) for c in cm.findall(q("comment"))))      # the comment just added
    t = etree.parse("u/word/document.xml"); root = t.getroot()
    s = root.find(f".//{q('commentRangeStart')}[@{q('id')}='{pid}']")
    e = root.find(f".//{q('commentRangeEnd')}[@{q('id')}='{pid}']")
    ref = root.find(f".//{q('commentReference')}[@{q('id')}='{pid}']")
    if s is not None:
        n = etree.Element(q("commentRangeStart")); n.set(q("id"), nid); s.addnext(n)
    if e is not None:
        n = etree.Element(q("commentRangeEnd")); n.set(q("id"), nid); e.addnext(n)
    run = ref.getparent()
    while run.getparent().tag not in (q("p"), q("ins"), q("hyperlink")):     # climb out of math / smartTag wrappers
        run = run.getparent()
    nr = etree.Element(q("r")); rp = etree.SubElement(nr, q("rPr")); st = etree.SubElement(rp, q("rStyle")); st.set(q("val"), "CommentReference")
    cr = etree.SubElement(nr, q("commentReference")); cr.set(q("id"), nid); run.addnext(nr)
    t.write("u/word/document.xml", xml_declaration=True, encoding="UTF-8", standalone=True)
print("ok", len(R), "replies")
