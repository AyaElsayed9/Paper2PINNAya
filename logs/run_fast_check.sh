#!/bin/bash
# Fast physics check vs full GBT on the 40 test designs of every lay-up expert (10 Nov 2026, plan B).  Resumable.
cd /home/claude/paper2pinnaya; mkdir -p runs/fast_check
export OMP_NUM_THREADS=2
for TAG in symmetry_crossply_SS_v2 periodic_crossply_SS_v2 symmetry_pm45_SS_v2 periodic_pm45_SS_v2 \
           symmetry_quasiiso_SS_v2 periodic_quasiiso_SS_v2 symmetry_zerorich_SS_v2 periodic_zerorich_SS_v2; do
  python3 -u scripts/fast_check.py $TAG --n 40 >> runs/fast_check/$TAG.log 2>&1 || exit 1
done
echo "FAST_CHECK_DONE $(date '+%d %b %H:%M')" >> runs/fast_check/done.log
