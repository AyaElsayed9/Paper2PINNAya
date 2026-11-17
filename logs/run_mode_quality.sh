#!/bin/bash
# 17 Nov 2026: mode quality (MAC, span of the 3 lowest GBT modes, family participation) for all eight experts.
cd /home/claude/paper2pinnaya
for L in crossply quasiiso zerorich pm45; do for S in symmetry periodic; do
  python3 -u scripts/mode_quality.py ${S}_${L}_SS_v2 >> runs/experts/mode_quality.log 2>&1 || { echo "FAILED $S $L" >> runs/experts/mode_quality.log; exit 1; }
  echo "done $S $L $(date -u '+%H:%M')" >> runs/experts/mode_quality.log
done; done
echo "MODE_QUALITY_DONE" >> runs/experts/mode_quality.log
