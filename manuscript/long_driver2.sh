#!/bin/bash
# remaining long cells, one process each (memory); beta <= 16
cd "$(dirname "$0")"
while pgrep -f "run_long1.py Cross-ply 14" > /dev/null; do sleep 5; done
for l in Sym1 Sym2 Unsym; do python3 run_long1.py $l 14 1 2>&1 | grep -v Warn; done
for l in Isotropic Cross-ply Sym1 Sym2 Unsym; do python3 run_long1.py $l 16 1 2>&1 | grep -v Warn; done
for b in 4 6 8 10 12 16; do for l in Isotropic Sym1 Sym2; do python3 run_long1.py $l $b 0 2>&1 | grep -v Warn; done; done
echo ALLDONE
