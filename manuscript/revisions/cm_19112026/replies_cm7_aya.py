"""replies_cm7 (22 Nov 2026): Aya's replies to Prof. Mittelstedt, polished in her own voice (polite student),
typos fixed, C58/C61/C71 replies updated to the fix_cm7 edits. Works on ua5/word/comments.xml."""
import re
from xml.sax.saxutils import escape

R = {
    "4": "Thank you, Professor, this is a very good point. I added four papers from other groups (Bisagni and Vescovini; "
         "Stamatelos et al.; Buermann et al.; Vescovini and Bisagni) and rewrote the paragraph by topic, so that the works "
         "of all groups are now mixed together.",
    "69": "Thank you, Professor, you are right that this was unclear. I numbered all 19 walls in Fig. 1c and added a small "
          "detail of a web root. The text now lists them (5 skin walls, 4 feet, 4 webs, 2 crowns, 4 link walls) and "
          "explains the numbers: 18 natural nodes (the 16 wall junctions and the two side nodes) plus 74 intermediate "
          "nodes give N = 92. There are 93 segments because an open contour with 92 nodes would have 91 segments, and "
          "each of the two closed cells (stringer and skin) adds one more.",
    "72": "Thank you, Professor. I corrected both: Fig. 1a now shows the origin O of the coordinate system, and in Fig. 1b "
          "the pitch b_s is now drawn between the axes of the two stringers.",
    "96": "Thank you, Professor, I agree that this term was not clear. A “hat function” is a function that is 1 at "
          "one node, 0 at all other nodes and linear in between (its graph looks like a hat). Each shear mode is such a "
          "unit warping at one node. I replaced the term by this description, here and in Section 2.2.4.",
    "107": "Thank you, Professor, I removed it.",
    "115": "N is the number of nodes of the cross-section (N = 92 for b_s = 150 mm, Section 2.1.1). Thank you for noticing, "
           "Professor; I repeated it here so that the reader does not have to look for it.",
    "119": "Thank you for reading it so carefully, Professor. I added the definition of a periodic nodal field (both side "
           "nodes move in the same way, Rv = 0) before the steps and wrote step (2) in words: all combinations of the "
           "remaining modes of a family whose side-edge differences cancel. If it is still not clear enough, please let "
           "me know and I will rewrite the whole paragraph.",
    "126": "I agree, Professor, and I am sorry that it was so hard to follow. The section now starts with the physical idea "
           "in plain words: in an open cross-section the warping follows directly from the in-plane displacements, because "
           "the membrane shear strain is zero. When the two side edges are joined, the contour becomes closed, the warping "
           "has to come back to its starting value after one loop, and therefore a constant shear flow can remain, as in "
           "closed thin-walled sections. I also explained in words what the minimization does.",
    "132": "Thank you for checking this, Professor. It is not zero in general: q and q′ are the two end nodes of segment "
           "p, so (Δu)_p is the change of the warping along this segment. Only the two side nodes s_L and s_R have "
           "equal warping. I added a sentence that says this.",
    "144": "Thank you, Professor. Bloch–Floquet theory is the theory of waves in periodic structures. I now say this in "
           "the text and start the paragraph with the idea in plain words: in an infinitely wide panel, a buckling mode "
           "repeats from bay to bay up to a phase shift, and the cell of two pitches contains the two cases, modes that "
           "repeat in every bay and modes that alternate from bay to bay.",
    "150": "You are right, Professor, the sentence was not complete. In the plane of the cross-section, modes 1 and 2 are "
           "rigid translations (in Z and in Y), but since their amplitude varies along the panel, they describe the "
           "out-of-plane and in-plane bending of the panel. I corrected the sentence.",
    "169": "I agree, Professor, the dots were not helpful. I removed them and wrote the shear-mode terms out in two new "
           "equations, (28a) and (28b), directly after Eq. (28), together with the three geometric tensors that they need. "
           "The Mindlin modes do not appear there, because they do not change the mid-surface displacements.",
    "415": "Yes, Professor, they are, and my wording was not precise. The point is the smoothness between the elements: with "
           "cubic Hermite elements the slope is continuous, but the curvature can jump from one element to the next (C1), "
           "whereas cubic B-splines also keep the curvature continuous across the knots (C2). I corrected the sentence.",
    "429": "Thank you, Professor, for this important question. Because the analysis is linear, the value of the reference "
           "load only scales the result (λ = N_cr). The equal line loads and the Abaqus shell model are the same as in "
           "our previous study [22] and in El Yaakoubi-Mesbah and Mittelstedt [13]. Following your remark, I also loaded "
           "the cell through flat, non-rotating loaded ends (a uniform end shortening), so that the end forces follow from "
           "the longitudinal strain and the wall stiffnesses. The critical force per unit width then changes by +0.2% to "
           "+4.6% for the Cross-ply and Unsym lay-ups, but by −8.7% to −15.8% for the angle-ply Sym1 and Sym2 "
           "lay-ups (GBT). Two Abaqus models with the same flat ends confirm the GBT values within 1.3% (Cross-ply 54.7 vs "
           "55.4 N/mm, Sym1 25.4 vs 25.5 N/mm). This is now stated in the text. Since GBT and Abaqus always use the same "
           "load, the comparisons in the paper are not affected.",
    "436": "Thank you, Professor. The formula is given in the sentence, but I agree that its meaning was not explained. I "
           "added one sentence in words: the MAC is the squared cosine of the angle between the two vectors of skin "
           "deflections, so it does not depend on how the modes are scaled; it is 1 for identical shapes and 0 for "
           "completely different (orthogonal) shapes.",
    "445": "Thank you, Professor. These values come directly from the code: for every mode, scaled to a unit maximum, I took "
           "the difference between its values at the two side nodes s_L and s_R, and then the largest absolute value over "
           "all conventional modes. The values of about 1e-15 are round-off errors of double precision; the warping value "
           "(2.3e-11) is larger because the warping comes from the numerical solution of Eq. (13).",
    "456": "Thank you, Professor, this is a fair point. Their in-plane shapes are indeed almost the same with opposite sign, "
           "but they are not the same mode: mode 8 stretches the webs about three times more than mode 12. In the critical "
           "mode their amplitudes have opposite signs (Fig. 7e), so the shapes add up while the web stretching cancels. "
           "I added this to the caption of Fig. 7.",
}

p = 'ua5/word/comments.xml'
x = open(p, encoding='utf-8').read()
for cid, txt in R.items():
    m = re.search(r'(<w:comment [^>]*w:id="%s"[^>]*>)(.*?)(</w:comment>)' % cid, x, re.S)
    assert m, cid
    body = m.group(2)
    assert body.count('<w:p ') + body.count('<w:p>') == 1, cid
    ann = re.search(r'<w:r><w:rPr><w:rStyle w:val="CommentReference"/></w:rPr><w:annotationRef/></w:r>', body).group(0)
    head = body[:body.index(ann) + len(ann)]
    nb = head + '<w:r><w:t xml:space="preserve">%s</w:t></w:r></w:p>' % escape(txt)
    x = x[:m.start(2)] + nb + x[m.end(2):]
open(p, 'w', encoding='utf-8').write(x)
print("ok replies", len(R))
