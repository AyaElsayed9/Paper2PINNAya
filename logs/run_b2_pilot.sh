#!/bin/bash
# Step B2 pilot (6 Nov 2026): hybrid training on the v2 periodic SS designs that already have exact targets
# (runs/b2/pilot_files.txt, fixed at the start; the last 32 are validation), 1 thread, while the other CPU
# makes the remaining labels.  Same settings as step A.  Resumable (--ckpt), so it survives the 2 h limit.
cd /home/claude/paper2pinnaya
export OMP_NUM_THREADS=1
N=$(wc -l < runs/b2/pilot_files.txt)
python3 -u scripts/train_epinn.py --train data/operators/periodic_mixed_SS_v2/train --files runs/b2/pilot_files.txt \
  --test data/operators/periodic_mixed_SS_v2/test --bc SS --method ngfull --n-train $N --val 32 --epochs 10 \
  --arch decoder --cache 16 --n-eval-train 16 --outer 30 --batch 16 --lr 1e-3 --lr-final 1e-5 --ema 0.995 \
  --R 8 --k 3 --seed 0 --threads 1 --no-test-log --exact data/exact_targets/periodic_mixed_SS_v2/train \
  --mode hybrid --ckpt runs/b2/pilot_hybrid_ckpt.pt --out runs/b2/pilot_hybrid.pt >> runs/b2/pilot_hybrid.log 2>&1
echo PILOT_TRAIN_DONE >> runs/b2/pilot_hybrid.log
