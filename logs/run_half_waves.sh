#!/bin/bash
# 17 Nov 2026: number of half-waves of the critical mode, all eight experts (scripts/half_waves.py).
cd /home/claude/paper2pinnaya
for L in crossply quasiiso zerorich pm45; do for S in symmetry periodic; do
  python3 -u scripts/half_waves.py ${S}_${L}_SS_v2 >> runs/experts/half_waves.log 2>&1 || { echo "FAILED $S $L" >> runs/experts/half_waves.log; exit 1; }
  echo "done $S $L $(date -u '+%H:%M')" >> runs/experts/half_waves.log
done; done
echo "HALF_WAVES_DONE" >> runs/experts/half_waves.log
