"""fix_cm4 (20 Nov 2026): C64 with the exact check (scripts/c64_flat_ends.py: plane loaded ends, end force on the axial
mode, distribution from the wall stiffnesses; 8 cells: -3.4 % to +1.0 %).  The first, approximate check (line loads
proportional to A11) gave -0.3 to -4.3 % and the wrong sign for Cross-ply 150 mm; it is replaced.  Replies of
C22, C32, C47, C49, C53, C58, C64 updated to what was done."""
import sys, re
sys.path.insert(0, '.')
from lxml import etree
from trk import *
import trk
trk.DATE = "2026-11-20T14:30:00Z"
t = etree.parse('u/word/document.xml'); root = t.getroot(); b = root.find('w:body', ns); cx = Ctx(root)
ME = trk.AUTH
OLD = ("A distribution following a uniform end strain, with the line load of every wall proportional to its axial membrane "
       "stiffness, changes the critical mean force per unit width by only 0.3–0.5% for the Cross-ply and Unsym lay-ups and by "
       "1.8–4.3% for the angle-ply Sym1 and Sym2 lay-ups (GBT, 32 cells).")
NEW = ("Loading the cell instead through plane loaded ends (a rigid end plate, as when the end edges are coupled to a "
       "reference point), so that the end forces distribute themselves according to the wall stiffnesses, changes the "
       "critical mean force per unit width by −3.4% to +1.0% (GBT, all four lay-ups, bs = 150 and 285.5 mm, β = 1.5).")
ps = [p for p in b.iter(q('w:p')) if OLD in vis_text(p)]; assert len(ps) == 1
# OLD lies entirely in my own insertion of fix_cm1: change the run text in place
p = ps[0]; v = vis_text(p); pos = v.find(OLD); o = 0; hit = []
for r, s in vis_runs(p):
    if max(pos, o) < min(pos + len(OLD), o + len(s)):
        hit.append((r, o))
    o += len(s)
assert all(r.getparent().tag == q('w:ins') and r.getparent().get(q('w:author')) == ME for r, _ in hit), "not own insertion"
r0, o0 = hit[0]; tt = r0.find(q('w:t'))
if len(hit) == 1:
    a = pos - o0; tt.text = tt.text[:a] + NEW + tt.text[a + len(OLD):]
else:                                    # spread over several of my runs: put NEW in the first, cut the rest
    a = pos - o0; first_len = len(tt.text) - a; tt.text = tt.text[:a] + NEW; rem = len(OLD) - first_len
    for r, oo in hit[1:]:
        t2 = r.find(q('w:t')); k = min(rem, len(t2.text)); t2.text = t2.text[k:]; rem -= k
assert NEW in vis_text(p)
t.write('u/word/document.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
# replies
c = etree.parse('u/word/comments.xml'); W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = {"Agreed. We will add": "Agreed. We added four references of other groups after [16] (Bisagni & Vescovini, TWS 2009; "
     "Stamatelos, Labeas & Tserpes, TWS 2011; Buermann et al., TWS 2006; Vescovini & Bisagni, Aeronaut. J. 2015); all later "
     "references are renumbered.",
     "Agreed; Fig. 1a will show": "Corrected: Fig. 1a now shows the origin O (loaded end, side edge s_R, skin mid-plane) with "
     "the X, Y, Z axes, and Fig. 1b the pitch b_s between the stringer axes. The link walls are drawn more clearly in Fig. 1c.",
     "We will rewrite this step": "We added a plain definition (a nodal field is periodic when both side nodes move "
     "identically, Rv = 0) and explained step (2) in words (all combinations whose side-edge differences cancel).",
     "Agreed. Section 2.2.3 will start": "Agreed. Section 2.2.3 now starts with the physical idea: joining the side edges "
     "closes the contour, so a constant shear flow around it can remain, and the warping is the minimum of the shear energy "
     "under this closure condition.",
     "Agreed; the paragraph will start": "Agreed; the paragraph now starts with the plain idea: in an infinitely wide panel a "
     "mode repeats from bay to bay up to a phase shift, and the two-pitch cell contains the two extreme cases.",
     "Agreed; the shear-mode terms will be": "The shear-mode terms have the same form as the conventional ones; the text now "
     "states exactly which terms replace which, and they are all contained in the general form that follows.",
     "Thank you. The magnitude is only a reference": "Thank you. The magnitude is only a reference (linear analysis, lambda = "
     "Ncr); the equal line loads follow our previous study [22] and [13]. We checked loading through plane loaded ends (rigid "
     "end plate, as with the end edges coupled to a reference point), where the end forces follow the wall stiffnesses: the "
     "critical mean force per unit width changes by -3.4% to +1.0% for all four lay-ups (GBT). This is now stated in the text."}
n = 0
for cm in c.getroot().findall("{%s}comment" % W):
    ts = list(cm.iter("{%s}t" % W)); full = "".join(x.text or "" for x in ts)
    for k, v in R.items():
        if full.startswith(k):
            ts[0].text = v
            for x in ts[1:]: x.text = ""
            n += 1
c.write('u/word/comments.xml', xml_declaration=True, encoding='UTF-8', standalone=True)
print("ok: C64 sentence replaced,", n, "replies updated")
