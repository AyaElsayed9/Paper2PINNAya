#!/bin/bash
# rebuild GBT_periodic-NK-CM-AME-20112026b.docx from AME-20112026 (fix_cm1 + replies_cm1) with fix_cm2..4 and Fig. 1
set -e
cd "$(dirname "$0")"
rm -rf u && mkdir u && (cd u && unzip -q ../GBT_periodic-NK-CM-AME-20112026.docx && find . -type l -delete)
python3 fix_cm2.py && python3 fix_cm3.py && python3 fix_cm4.py
python3 -c "
from PIL import Image
b=Image.open('../../fig/fig01_geometry.png').convert('RGB'); c=Image.new('RGB',(3455,2234),'white'); c.paste(b,(0,8)); c.save('u/word/media/image1.png', dpi=(600,600))"
(cd u && rm -f ../GBT_periodic-NK-CM-AME-20112026b.docx && zip -qXr ../GBT_periodic-NK-CM-AME-20112026b.docx . -x '*/')
D=/root/.claude/skills/synced/89684747-8288-4727-bcdd-0f93eb064259_c99e516e-a235-41ea-88e1-b79caa917733/docx/scripts
python3 $D/office/validate.py GBT_periodic-NK-CM-AME-20112026b.docx --original GBT_periodic-NK-CM-19112026.docx --author "Aya Medhat Elsayed" | tail -2
rm -rf u __pycache__
