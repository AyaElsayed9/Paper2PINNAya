#!/bin/bash
# Restart of the data-only arm (3 Nov 2026): the first try was killed by the out-of-memory killer
# (data and hybrid together > 8 GB). Now: wait for hybrid, then data alone; certificates meanwhile.
cd /home/claude/paper2pinnaya
while kill -0 1758 2>/dev/null; do sleep 30; done
COMMON="--train data/operators/periodic_industrial_SS/train --test data/operators/periodic_industrial_SS/test --bc SS --method ngfull --n-train 400 --epochs 10 --arch decoder --cache 16 --n-eval-train 16 --outer 30 --batch 16 --lr 1e-3 --R 8 --k 3 --seed 0 --threads 1 --exact data/exact_targets/periodic_industrial_SS"
python3 -u scripts/train_epinn.py $COMMON --mode data --out runs/w3_data.pt > runs/w3_data.log 2>&1 &
D=$!
for m in p1_ss_dec w3_hybrid; do
  python3 -u scripts/eval_epinn.py runs/$m.pt --test data/operators/periodic_industrial_SS/test --steps 3,5,8,12 --threads 1 > runs/${m}_cert.log 2>&1
done
wait $D
python3 -u scripts/eval_epinn.py runs/w3_data.pt --test data/operators/periodic_industrial_SS/test --steps 3,5,8,12 --threads 2 > runs/w3_data_cert.log 2>&1
echo THREEWAY_DONE
