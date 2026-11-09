#!/bin/bash
# Lay-up experts (Aya, 7 Nov 2026): one small PINN per realistic lay-up and side condition, inputs b_s and beta
# only, same GBT physics (ngfull, hybrid with GBT labels), same physics certificate.  Every stage is resumable
# (build and labels skip finished designs, training continues from --ckpt), so the script is simply started
# again after the 2 h limit or a machine restart.
#     LAY=crossply SIDE=symmetry bash logs/run_expert.sh
# Sets per expert: train file 152 designs (seed 11, Sobol in b_s, beta; last 32 = validation), test 40 (seed 99).
cd /home/claude/paper2pinnaya
LAY=${LAY:?set LAY}; SIDE=${SIDE:?set SIDE}; BC=${BC:-SS}
TAG=${SIDE}_${LAY}_${BC}_v2
OPS=data/operators/$TAG; EX=data/exact_targets/$TAG; RUN=runs/experts; mkdir -p $RUN results/experts
LOG=$RUN/$TAG.log
export OMP_NUM_THREADS=1
stage() { echo "=== $1  $(date '+%d %b %H:%M')" >> $LOG; }

if [ ! -f $OPS/BUILD_DONE ]; then
  stage "build operators (full GBT for every design)"
  python3 -u scripts/build_operators.py --side $SIDE --bc $BC --layups $LAY --space v2 --set test --n 40 --seed 99 --full --workers 2 >> $LOG 2>&1 &&
  python3 -u scripts/build_operators.py --side $SIDE --bc $BC --layups $LAY --space v2 --set train --n 152 --seed 11 --full --workers 2 >> $LOG 2>&1 &&
  touch $OPS/BUILD_DONE || exit 1
fi
# training-file list in sample order (the last 32 are the validation designs)
python3 - <<EOF
import sys; sys.path.insert(0, "src")
from pinnbuck.designspace import sample_design
S = sample_design("omega", 152, seed=11, categorical={"side": "$SIDE", "bc": "$BC"}, layups="$LAY", space="v2")
open("$RUN/${TAG}_files.txt", "w").write("".join(f"$OPS/train/{s.sid}.pkl\n" for s in S))
EOF
if [ ! -f $EX/train/LABELS_DONE ]; then
  stage "GBT labels (LOBPCG R 8, 400 iterations)"
  for P in 0 1; do
    python3 -u scripts/exact_targets.py $OPS/train $EX/train --iters 400 --part $P --parts 2 >> $RUN/${TAG}_labels_p$P.log 2>&1 &
  done
  wait
  [ $(ls $EX/train/*.pt | wc -l) -ge 152 ] && touch $EX/train/LABELS_DONE || exit 1
fi
NL=${NL:-120}
MTAG=${TAG}_hybrid_L$NL; MODE="--mode hybrid --n-labels $NL"
[ "$NL" = 0 ] && { MTAG=${TAG}_physics_L0; MODE="--mode physics"; }     # label study (9 Nov): physics only
if [ ! -f $RUN/$MTAG.pt ]; then
  stage "train $MTAG ($NL of 120 training designs labelled)"
  python3 -u scripts/train_epinn.py --train $OPS/train --files $RUN/${TAG}_files.txt --test $OPS/test --bc $BC \
    --method ngfull --n-train 152 --val 32 --epochs 10 --arch decoder --cache ${CACHE:-8} --n-eval-train 16 --outer 30 \
    --batch 16 --lr 1e-3 --lr-final 1e-5 --ema 0.995 --R 8 --k 3 --seed 0 --threads 2 --no-test-log \
    --inputs geom --exact $EX/train $MODE --ckpt $RUN/${MTAG}_ckpt.pt --out $RUN/$MTAG.pt >> $RUN/$MTAG.log 2>&1 || exit 1
fi
if [ ! -f results/experts/${MTAG}_cert.log ]; then
  stage "physics certificate on the 40 test designs"
  python3 -u scripts/eval_epinn.py $RUN/$MTAG.pt --test $OPS/test --steps 3,5,8,12 --threads 2 > results/experts/${MTAG}_cert.tmp 2>&1 &&
  mv results/experts/${MTAG}_cert.tmp results/experts/${MTAG}_cert.log || exit 1
fi
stage "EXPERT_DONE $MTAG"
