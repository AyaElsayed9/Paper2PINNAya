# Step B1: v2 operator sets, periodic SS, mixed lay-ups, full GBT for every design (resumable)
python3 -u scripts/build_operators.py --side periodic --bc SS --layups mixed --space v2 --set test --n 100 --seed 99 --full --workers 2 >> logs/b1_build_periodic_SS_test.log 2>&1
python3 -u scripts/build_operators.py --side periodic --bc SS --layups mixed --space v2 --set train --n 1064 --seed 11 --full --workers 2 >> logs/b1_build_periodic_SS_train.log 2>&1
echo B1_BUILD_DONE >> logs/b1_build_periodic_SS_train.log
