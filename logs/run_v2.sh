#!/bin/bash
# Step A (4 Nov 2026): same three-way test with clean training: 288 training + 32 validation designs,
# cosine learning-rate decay 1e-3 -> 1e-5, weight average (EMA 0.995), weights chosen on validation,
# test designs not looked at during training. One run at a time (memory).
cd /home/claude/paper2pinnaya
COMMON="--train data/operators/periodic_industrial_SS/train --test data/operators/periodic_industrial_SS/test --bc SS --method ngfull --n-train 400 --val 32 --epochs 10 --arch decoder --cache 16 --n-eval-train 16 --outer 30 --batch 16 --lr 1e-3 --lr-final 1e-5 --ema 0.995 --R 8 --k 3 --seed 0 --threads 2 --no-test-log --exact data/exact_targets/periodic_industrial_SS"
for m in hybrid data physics; do
  python3 -u scripts/train_epinn.py $COMMON --mode $m --out runs/v2_$m.pt > runs/v2_$m.log 2>&1
  python3 -u scripts/eval_epinn.py runs/v2_$m.pt --test data/operators/periodic_industrial_SS/test --steps 3,5,8,12 --threads 2 > runs/v2_${m}_cert.log 2>&1
  echo "$m done $(date)"
done
echo V2_DONE
