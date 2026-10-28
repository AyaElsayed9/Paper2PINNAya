#!/bin/bash
cd "$(dirname "$0")"
python3 probe_hx.py Unsym 12 10 2>&1 | grep -v Warn
python3 probe_hx.py Sym1 12 10 2>&1 | grep -v Warn
for l in Sym2 Unsym; do python3 run_long1.py $l 14 1 9 2>&1 | grep -v Warn; done
for l in Isotropic Cross-ply Sym1 Sym2 Unsym; do python3 run_long1.py $l 16 1 10 2>&1 | grep -v Warn; done
for b in 14 16; do for l in Isotropic Sym1 Sym2; do python3 run_long1.py $l $b 0 10 2>&1 | grep -v Warn; done; done
echo ALLDONE
