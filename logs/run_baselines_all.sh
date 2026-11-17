#!/bin/bash
# 17 Nov 2026: random / nn10 / nn120 starts of the physics check for the 7 experts that do not have them yet
# (crossply / symmetry has them since 9 Nov), 12 steps, same test designs, seed and procedure as eval_epinn.py.
# Resumable (rows already in results/check_baselines/<tag>.csv are skipped).  Run in segments under the 2 h limit:
#     timeout 6000 bash logs/run_baselines_all.sh
cd /home/claude/paper2pinnaya
for T in periodic_crossply symmetry_quasiiso periodic_quasiiso symmetry_zerorich periodic_zerorich symmetry_pm45 periodic_pm45; do
  python3 -u scripts/check_start_baselines.py ${T}_SS_v2 --models hybrid_L120 --nn 10,120 --steps 3,5,8,12 --starts random,nn10,nn120 --threads 2 \
    >> runs/experts/baselines_all.log 2>&1 || { echo "FAILED $T $(date '+%d %b %H:%M')" >> runs/experts/baselines_all.log; exit 1; }
  echo "done $T $(date '+%d %b %H:%M')" >> runs/experts/baselines_all.log
done
echo "BASELINES_ALL_DONE $(date '+%d %b %H:%M')" >> runs/experts/baselines_all.log
