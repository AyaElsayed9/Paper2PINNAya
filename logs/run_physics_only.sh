#!/bin/bash
# Physics-only training (0 GBT labels in the loss) of the six experts that do not have it yet (Aya, 16 Nov 2026).
# Same settings as the two cross-ply runs of 9 Nov (run_expert.sh with NL=0: ngfull, decoder, 120 training designs,
# 30 outer iterations, seed 0), then the physics check on the 40 test designs.  CACHE=3: memory-safe, same results.
# Note: as in the cross-ply runs, the saved weights are chosen on the 32 validation designs, whose GBT loads are
# used for this choice only (both cross-ply runs chose the last iteration).
# Resumable: start it again after the 2 h limit or a restart; finished experts are skipped.
#     bash logs/run_physics_only.sh
# Do NOT run other heavy jobs beside it (9 Nov: a second job made the training run out of memory).
# While it runs, its new files are hidden from `git status` by .git/info/exclude and the six stage logs by
# `git update-index --skip-worktree`; when the chain is done:  git update-index --no-skip-worktree runs/experts/*_SS_v2.log,
# remove the lines from .git/info/exclude, then commit runs/experts/*_physics_L0.* and results/experts/*_physics_L0_cert.log.
cd /home/claude/paper2pinnaya
echo "physics-only chain started / restarted $(date '+%d %b %H:%M')" >> runs/experts/physics_only_chain.log
for X in quasiiso/symmetry quasiiso/periodic zerorich/symmetry zerorich/periodic pm45/symmetry pm45/periodic; do
  LAY=${X%/*} SIDE=${X#*/} NL=0 CACHE=3 bash logs/run_expert.sh || { echo "FAILED $X $(date '+%d %b %H:%M')" >> runs/experts/physics_only_chain.log; exit 1; }
  echo "done $X $(date '+%d %b %H:%M')" >> runs/experts/physics_only_chain.log
done
echo "PHYSICS_ONLY_CHAIN_DONE $(date '+%d %b %H:%M')" >> runs/experts/physics_only_chain.log
