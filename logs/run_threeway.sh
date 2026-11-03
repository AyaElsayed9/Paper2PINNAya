#!/bin/bash
cd /home/claude/paper2pinnaya
until [ $(ls data/exact_targets/periodic_industrial_SS | wc -l) -ge 320 ]; do sleep 30; done
sleep 5
COMMON="--train data/operators/periodic_industrial_SS/train --test data/operators/periodic_industrial_SS/test --bc SS --method ngfull --n-train 400 --epochs 10 --arch decoder --cache 16 --n-eval-train 16 --outer 30 --batch 16 --lr 1e-3 --R 8 --k 3 --seed 0 --threads 1 --exact data/exact_targets/periodic_industrial_SS"
python3 scripts/train_epinn.py $COMMON --mode data   --out runs/w3_data.pt   > runs/w3_data.log 2>&1 &
python3 scripts/train_epinn.py $COMMON --mode hybrid --out runs/w3_hybrid.pt > runs/w3_hybrid.log 2>&1 &
wait
for m in p1_ss_dec w3_data w3_hybrid; do
  python3 scripts/eval_epinn.py runs/$m.pt --test data/operators/periodic_industrial_SS/test --steps 3,5,8,12 --threads 2 > runs/${m}_cert.log 2>&1
done
echo THREEWAY_DONE
