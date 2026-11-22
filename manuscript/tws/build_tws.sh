#!/bin/bash
# Thin-Walled Structures version: accept all tracked changes of the revised file, then reformat (tws_build.py)
set -e
cd "$(dirname "$0")"
python3 accept.py ../revisions/cm_19112026/aya/GBT_periodic-NK-CM-AME-21112026_5_1_Aya_final.docx /tmp/_accepted.docx
python3 tws_build.py /tmp/_accepted.docx GBT_periodic_TWS.docx
rm -f /tmp/_accepted.docx
