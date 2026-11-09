#!/bin/bash
# Chain (Aya, 9 Nov 2026): (1) finish level 2 with its best-on-validation weights (round 8; stopped at outer 19
# because validation got worse after round 8) and run its certificate on the 45 industrial test designs;
# (2) label study on the experts: crossply / symmetry (Paper 1's case) with NL = 0, 10, 20, 40 labels of 120,
# then NL = 0 (physics only) for the other 7 experts.  Every step is resumable; just start the chain again.
cd /home/claude/paper2pinnaya
OUTER=19 CACHE=2 NL=518 bash logs/run_level2.sh || exit 1
for NL in 0 10 20 40; do LAY=crossply SIDE=symmetry NL=$NL bash logs/run_expert.sh || exit 1; done
for X in crossply/periodic quasiiso/symmetry quasiiso/periodic zerorich/symmetry zerorich/periodic pm45/symmetry pm45/periodic; do
  LAY=${X%/*} SIDE=${X#*/} NL=0 bash logs/run_expert.sh || exit 1
done
echo "LABEL_STUDY_CHAIN_DONE $(date '+%d %b %H:%M')" >> runs/experts/label_study.log
