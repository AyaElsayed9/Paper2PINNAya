# head from the chosen weights (state), mode shapes from the final EMA / final weights
for m in hybrid data; do
  for w in state_ema state_last; do
    python3 -u scripts/eval_epinn.py runs/v2_$m.pt --modes $w --test data/operators/periodic_industrial_SS/test --steps 3,5,8,12 --threads 1 > runs/v2_${m}_cert_modes_${w}.log 2>&1
  done
done
echo SPLIT_DONE >> logs/run_split_eval.log
