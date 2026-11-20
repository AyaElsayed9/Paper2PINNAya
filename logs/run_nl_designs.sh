#!/bin/bash
# 20 Nov 2026: non-linear GBT paths of the three flagged realistic designs (self-ending)
cd /home/claude/paper2pinnaya
for n in qi_152_20_750 pm45_152_20_750 qi_124_20_500; do
  [ -f results/koiter/nl_$n.npz ] && continue
  timeout 6500 python3 -u scripts/nl_designs.py $n >> runs/experts/nl_designs.log 2>&1 || echo "FAILED_OR_TIMEOUT $n" >> runs/experts/nl_designs.log
done
echo NL_DESIGNS_DONE >> runs/experts/nl_designs.log
