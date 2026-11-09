#!/bin/bash
# 9 Nov 2026 night (Aya: stop the label study, fix the shapes): resumable chain
#   1. wait until the running periodic crossply NL 0 run is finished (its certificate exists)
#   2. clean timing with nothing else running (scripts/time_check.py, 10 crossply / symmetry test designs)
#   3. neighbour + correction networks (ARCH=nnres), crossply / symmetry, NL 10 then NL 120
cd /home/claude/paper2pinnaya
until [ -f results/experts/periodic_crossply_SS_v2_physics_L0_cert.log ]; do sleep 30; done
while pgrep -f "train_epinn|eval_epinn" >/dev/null; do sleep 10; done
rm -f runs/experts/periodic_crossply_SS_v2_physics_L0_ckpt.pt
if [ ! -f results/timing/check_timing.txt ]; then
  python3 -u scripts/time_check.py --n 10 --threads 2 > runs/check/time_check.log 2>&1 || exit 1
fi
for NL in 10 120; do
  ARCH=nnres NL=$NL LAY=crossply SIDE=symmetry CACHE=4 bash logs/run_expert_nnres.sh || exit 1
done
echo "NIGHT_9NOV_DONE $(date '+%d %b %H:%M')" >> runs/check/night_9nov.log
