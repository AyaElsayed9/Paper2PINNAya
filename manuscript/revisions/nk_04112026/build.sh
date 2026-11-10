set -e
F=/home/claude/paper2pinnaya/manuscript/fig; SK=$(ls -d /root/.claude/skills/synced/*/docx)
rm -rf w2; cp -r u w2
python3 edit.py u/word/document.xml w2/word/document.xml
python3 eq.py w2/word/document.xml
cp $F/fig01_geometry.png w2/word/media/image1.png; cp $F/fig03_prebuckling.png w2/word/media/image3.png; cp $F/fig06_buckling_modes.png w2/word/media/image6.png
python3 replies.py
sed -i 's/w:date="2026-11-10T[0-9:]*Z" w:initials="AME"/w:date="2026-11-10T11:00:00Z" w:initials="AME"/g' w2/word/comments.xml
OUT="GBT_periodic-NK-04112026_AME-10112026.docx"
rm -f "$OUT"; (cd w2 && zip -qXr -D "../$OUT" .)
python3 $SK/scripts/office/validate.py "$OUT" --original nk.docx --author "Aya Medhat Elsayed" 2>&1 | tail -2
rm -rf w3; cp -r w2 w3; python3 accept.py w3/word/document.xml
CL="GBT_periodic-NK-04112026_AME-10112026_clean.docx"
rm -f "$CL"; (cd w3 && zip -qXr -D "../$CL" .)
python3 $SK/scripts/office/validate.py "$CL" 2>&1 | tail -1
