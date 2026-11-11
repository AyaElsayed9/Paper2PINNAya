set -e
F=/home/claude/paper2pinnaya/manuscript/fig; SK=$(ls -d /root/.claude/skills/synced/*/docx)
rm -rf w2; cp -r u w2
python3 edit.py u/word/document.xml w2/word/document.xml
python3 eq.py w2/word/document.xml
cp $F/fig01_geometry.png w2/word/media/image1.png; cp $F/fig03_prebuckling.png w2/word/media/image3.png; cp $F/fig06_buckling_modes.png w2/word/media/image6.png; cp $F/fig09_long_modes.png w2/word/media/image9.png
cp $F/scheme01_overview.png w2/word/media/scheme1.png
cp $F/fig10_npitch.png w2/word/media/npitch.png
python3 - <<'PY'
p='w2/word/_rels/document.xml.rels'; s=open(p,encoding='utf8').read()
assert 'rIdScheme1' not in s
s=s.replace('</Relationships>','<Relationship Id="rIdScheme1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/scheme1.png"/><Relationship Id="rIdNpitch" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/npitch.png"/></Relationships>')
open(p,'w',encoding='utf8').write(s)
PY
python3 replies.py
sed -i 's/w:date="2026-11-10T[0-9:]*Z" w:initials="AME"/w:date="2026-11-10T11:00:00Z" w:initials="AME"/g' w2/word/comments.xml
OUT="GBT_periodic-NK-04112026_AME-11112026.docx"
rm -f "$OUT"; (cd w2 && zip -qXr -D "../$OUT" .)
python3 $SK/scripts/office/validate.py "$OUT" --original nk.docx --author "Aya Medhat Elsayed" 2>&1 | tail -2
rm -rf w3; cp -r w2 w3; python3 accept.py w3/word/document.xml
CL="GBT_periodic-NK-04112026_AME-11112026_clean.docx"
rm -f "$CL"; (cd w3 && zip -qXr -D "../$CL" .)
python3 $SK/scripts/office/validate.py "$CL" 2>&1 | tail -1
