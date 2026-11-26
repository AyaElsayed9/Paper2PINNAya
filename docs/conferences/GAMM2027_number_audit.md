# Number audit of the GAMM 2027 abstract and the progress report (26 Nov 2026)

Every number was recomputed from the result files (`scripts/audit_abstract.py`) or looked up in the log entry
that produced it. "Was" = text before this audit.

## Abstract
| Claim | Was | Verified value | Source |
|---|---|---|---|
| GBT vs Abaqus, 60 periodic cases | within 1.6 % | max 1.52 % (SS max 1.25, CC max 1.52; means 0.57 / 0.64) | GBT paper, Sec. validation |
| sample panels of the POD | 64 | 64 | RESEARCH_LOG l. 30 |
| basis error over 1064 designs | at most 0.5 % | max 0.49 %, mean 0.061 % | RESEARCH_LOG (v2 space) |
| encoder / decoder | 4 x 256 / 2 x 256 | depth 4, width 256; dec 2 x 256 (default, not overridden) | runs/experts/*.pt, src/pinnbuck/epinn.py |
| candidate modes, training designs | 8, 152 | R = 8, n_train = 152 | runs/experts/*.pt |
| full GBT after | 40 iterations | MAX_STEPS = 40 | scripts/cost_study.py |
| one network for all lay-ups | max above 20 % | 1000-label net val. max 20.9-24.5 %; pilot 29.8 %; GP/MLP/GBR 27.0-29.7 % | RESEARCH_LOG 6-7 Nov |
| unseen designs | 320 | 8 experts x 40 | report_data.json |
| 120 labels, instant | within 4.9 % | max 4.918 % | results/error_indicator/*_SS_v2.csv |
| physics only overestimates | yes | 319 of 320 designs above GBT | same files |
| estimate never below true error | 3800 checks | 3800 of 3800 (760 runs x steps 0,1,3,5,12) | same files |
| 5 % tolerance, iterations | 3 to 15 | net_L0 3-15, net_L120 3-15 (random 5-18); 0 of 320 fallbacks | results/cost_study/*_SS_v2.csv |
| no labels, max error | 2.5 % | 2.51 % | same files |
| time per design | 7 to 8 s vs 11 s | net_L0 7.36 s vs full GBT 10.71 s | results/cost_study/summary_tol5.csv |
| 22 Abaqus points, checked result | max 1.19 % | L120 1.19 %, L0 0.67 % | results_p1/pinn_p1.jsonl |
| 22 Abaqus points, Paper 1 PINN | 10.80 % | 10.80 % | Paper 1 Table 2 |

## Report (all 209 distinct numbers found in the log / result files; context checked sentence by sentence)
Corrected:
| Place | Was | Now | Why |
|---|---|---|---|
| abstract, conclusion 2 | instant answer "within 0.24-1.18 %" | 0.24-1.18 % on average, at most 4.92 % | 0.24-1.18 are expert means, not a bound |
| Sec. 3 intro | checked answer "below 0.5 % in every test" | "in every test shown in the figure"; 140 % marked as training-panel error | the 26 Oct clamped test had 0.54 % after 12 steps |
| Sec. 2 | GBT vs Abaqus 0.57-0.64 % on average | + 1.52 % at most in 60 cases | maximum added |
| basis | own 80 lowest modes 1.6 % | own 80 modes, 1.6 % on average | log says "own 80 modes: 1.6 % mean" |
| data models | about 28 % at most | 27.0-29.7 % at most | exact range (mlp5 29.7, GP 27.9, gbr 27.0) |
| industrial network | 1.35 % / 0.14 % | + maxima 8.40 % / 0.74 % | maxima added |
| mode MAC | 0.39-0.40 against 0.29 | against 0.23-0.29 | nearest of 10 = 0.23, of 120 = 0.29 |
| cost | network training "pays off after about 2400 designs" | vs full GBT 1900-4400 designs; vs random start 9000-25,000 | 2400 was the break-even vs full GBT, but the sentence compared with the random start |
| Table 2 caption | error | mean error | the columns are means |
Verified unchanged (examples): 22-point table (GBT 0.54 / 0.97, tight 0.89 / 0.97, other 0.46 / 0.95; Paper 1
2.53 / 10.80, 8.08 / 10.80, 1.29 / 3.54; half-wave number 22 / 22); cost table (all 9 rows equal to the pooled
summary_tol5.csv); expert ranges (0 labels 4.68-11.93 % mean, max 16.1-24.7 %; 120 labels 0.24-1.18 %, max
0.87-4.92 %; 12 steps 0.11-0.48 %); estimate conservative 2.9-11.1x in the median, Spearman 0.88 at 5 steps
(0.87 physics only, 0.89 with labels); 640 checks at 5 steps; fixed-width effect +2.8 / +4.0 % (cross-ply),
+21-160 % (others).
