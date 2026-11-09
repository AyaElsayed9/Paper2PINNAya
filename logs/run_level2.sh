#!/bin/bash
# Level 2 (Aya, 8 Nov 2026): ONE PINN for all industrial lay-ups (symmetric, balanced, 0/+-45/90, >= 10 % per
# direction; 24 skin x 12 stringer sequences), periodic side, SS, v2 space.  Inputs: all 38 (b_s, beta, flags,
# lamination parameters of skin and stringer, closed-form loads, wave content).  Same GBT physics (ngfull,
# hybrid), same certificate.  Reuses the B2 operators and GBT labels (no new GBT runs): the 518 industrial
# designs of the 1000 B2 training designs + the 35 industrial of the 64 B2 validation designs (validation);
# test = the 45 industrial of the 100 v2 test designs (symlinked).  Resumable (--ckpt).
#     NL=518 bash logs/run_level2.sh      (all 518 training designs labelled; NL = 0 -> physics only)
#     OUTER=19 (9 Nov): finish an earlier stopped run at its checkpoint (no more rounds; best-on-validation weights kept)
cd /home/claude/paper2pinnaya
export OMP_NUM_THREADS=2
NL=${NL:-518}
TAG=industrial_hybrid_L$NL; MODE="--mode hybrid --n-labels $NL"
[ "$NL" = 0 ] && { TAG=industrial_physics_L0; MODE="--mode physics"; }
RUN=runs/level2; TEST=data/operators/periodic_industrial_SS_v2/test
if [ ! -f $RUN/$TAG.pt ]; then
  python3 -u scripts/train_epinn.py --train data/operators/periodic_mixed_SS_v2/train --files $RUN/files_industrial.txt \
    --test $TEST --bc SS --method ngfull --n-train 553 --val 35 --epochs 10 --arch decoder --cache ${CACHE:-8} \
    --n-eval-train 16 --outer ${OUTER:-30} --batch 16 --lr 1e-3 --lr-final 1e-5 --ema 0.995 --R 8 --k 3 --seed 0 --threads 2 \
    --no-test-log --inputs all --exact data/exact_targets/periodic_mixed_SS_v2/train $MODE \
    --ckpt $RUN/${TAG}_ckpt.pt --out $RUN/$TAG.pt >> $RUN/$TAG.log 2>&1 || exit 1
fi
if [ ! -f results/level2/${TAG}_cert.log ]; then
  mkdir -p results/level2
  python3 -u scripts/eval_epinn.py $RUN/$TAG.pt --test $TEST --steps 3,5,8,12 --threads 2 > results/level2/${TAG}_cert.tmp 2>&1 &&
  mv results/level2/${TAG}_cert.tmp results/level2/${TAG}_cert.log || exit 1
fi
echo "LEVEL2_DONE $TAG $(date '+%d %b %H:%M')" >> $RUN/$TAG.log
