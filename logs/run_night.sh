#!/bin/bash
# Night chain (9 Nov 2026): expert 8 (pm45 / periodic), then level 2 (industrial lay-ups, all labelled).
# Both stages resumable, so the script is simply started again after the 2 h limit.
cd /home/claude/paper2pinnaya
LAY=pm45 SIDE=periodic bash logs/run_expert.sh || exit 1
grep -q EXPERT_DONE runs/experts/periodic_pm45_SS_v2.log || exit 1
NL=518 bash logs/run_level2.sh
