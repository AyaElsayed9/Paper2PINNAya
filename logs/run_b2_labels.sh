#!/bin/bash
# Step B2 label study (7 Nov 2026): hybrid training on all 1000 v2 periodic SS train designs (physics on all),
# GBT labels on the first NL of them (sample order = random), the 64 validation designs (last 64 of
# runs/b2/files_train1000_val64.txt) always labelled; NL = 0 -> physics only.  Step-A settings, 1 thread,
# resumable (--ckpt), so it survives the 2 h limit and machine restarts.
#     NL=936 bash logs/run_b2_labels.sh      (all train designs labelled)
cd /home/claude/paper2pinnaya
export OMP_NUM_THREADS=1
NL=${NL:?set NL}
TAG=hybrid_L$NL
MODE="--mode hybrid --n-labels $NL"
[ "$NL" = 0 ] && { TAG=physics_L0; MODE="--mode physics"; }
python3 -u scripts/train_epinn.py --train data/operators/periodic_mixed_SS_v2/train --files runs/b2/files_train1000_val64.txt \
  --test data/operators/periodic_mixed_SS_v2/test --bc SS --method ngfull --n-train 1064 --val 64 --epochs 10 \
  --arch decoder --cache 2 --n-eval-train 16 --outer 30 --batch 16 --lr 1e-3 --lr-final 1e-5 --ema 0.995 \
  --R 8 --k 3 --seed 0 --threads 1 --no-test-log --exact data/exact_targets/periodic_mixed_SS_v2/train \
  $MODE --ckpt runs/b2/${TAG}_ckpt.pt --out runs/b2/$TAG.pt >> runs/b2/$TAG.log 2>&1
echo B2_TRAIN_DONE >> runs/b2/$TAG.log
