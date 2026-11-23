#!/bin/bash
# error-indicator study (23 Nov 2026): all 8 experts, resumable; runs for at most $1 seconds in total
B=${1:-560}; T0=$(date +%s)
for t in symmetry_crossply_SS_v2 periodic_crossply_SS_v2 symmetry_quasiiso_SS_v2 periodic_quasiiso_SS_v2 \
         symmetry_zerorich_SS_v2 periodic_zerorich_SS_v2 symmetry_pm45_SS_v2 periodic_pm45_SS_v2; do
  left=$(( B - ($(date +%s) - T0) )); [ $left -lt 40 ] && break
  python3 scripts/error_indicator.py $t --budget $(( left - 25 )) 2>&1 | grep -v Warn | grep -c "head" | sed "s/^/$t new rows: /"
done
