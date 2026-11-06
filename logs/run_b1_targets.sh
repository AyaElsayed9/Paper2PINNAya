# Step B1: exact targets (LOBPCG, R 8, 400 iterations, dictionary d2 with e24s15 end layers since 5 Nov 22:00), resumable.
# From 6 Nov 09:30 one worker only (PARTS=1) while the pilot trained; pilot done, run with PARTS=2 again.
PARTS=${PARTS:-1}
for SET in test train; do
  for P in $(seq 0 $((PARTS - 1))); do
    python3 -u scripts/exact_targets.py data/operators/periodic_mixed_SS_v2/$SET data/exact_targets/periodic_mixed_SS_v2/$SET --iters 400 --part $P --parts $PARTS >> logs/b1_targets_${SET}_p$P.log 2>&1 &
  done
  wait
  echo "B1_TARGETS_${SET}_DONE" >> logs/b1_targets_${SET}_p0.log
done
