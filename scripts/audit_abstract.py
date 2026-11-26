"""Recompute every number of the GAMM abstract from the raw result files (26 Nov 2026)."""
import csv, glob, json, collections, re, subprocess
import numpy as np, torch
out = []
def say(claim, value, source): out.append((claim, value, source)); print(f"{claim:62s} | {value:38s} | {source}")

# expert networks
ck = torch.load("runs/experts/symmetry_crossply_SS_v2_hybrid_L120.pt", weights_only=False)["args"]
say("8 candidate modes; encoder layers x width; training designs", f"R={ck['R']}, depth={ck['depth']}, width={ck['width']}, n_train={ck['n_train']}", "runs/experts/*.pt args")
import inspect, sys; sys.path.insert(0, "src"); from pinnbuck import epinn
sig = inspect.signature(epinn.FunctionDecoderPINN.__init__)
say("decoder layers x width", f"dec_depth={sig.parameters['dec_depth'].default}, dec_width={sig.parameters['dec_width'].default} (not overridden)", "src/pinnbuck/epinn.py")
D = json.load(open("results/report_gamm/report_data.json"))
say("number of experts / designs per expert", f"{len(D['experts'])} experts", "report_data.json")
mx = [e["hybrid_L120"]["head_max"] for e in D["experts"]]
say("instant answer, 120 labels, max over 320 designs", f"{max(mx):.3f} %", "report_data.json (from error_indicator CSVs)")
# head bias without labels
up = n = 0; mxp = 0
for f in glob.glob("results/error_indicator/*_SS_v2.csv"):
    for r in csv.DictReader(open(f)):
        if r["model"] == "physics_L0":
            e = float(r["N_head"]) / float(r["N_full"]) - 1; n += 1; up += e > 0; mxp = max(mxp, abs(e) * 100)
say("physics-only instant answer above GBT", f"{up} of {n} designs (max {mxp:.2f} %)", "results/error_indicator/*_SS_v2.csv")
# error estimate safety
safe = tot = 0
for f in glob.glob("results/error_indicator/*_SS_v2.csv"):
    for r in csv.DictReader(open(f)):
        for st in ("0", "1", "3", "5", "12"):
            if not r["model"].startswith("gp"):
                true = float(r[f"N{st}"]) / float(r["N_full"]) - 1; tot += 1; safe += float(r[f"eta{st}"]) >= true
say("error estimate >= true error", f"{safe} of {tot} checks", "results/error_indicator/*_SS_v2.csv")
# cost study traces
st = collections.defaultdict(list); fb = collections.Counter(); err = collections.defaultdict(list); ns = collections.Counter()
for f in glob.glob("results/cost_study/*_SS_v2.csv"):
    for r in csv.DictReader(open(f)):
        e = json.loads(r["eta"]); N = json.loads(r["N"]); k = next((i for i, v in enumerate(e) if v <= 0.05), None); ns[r["start"]] += 1
        if k is None: fb[r["start"]] += 1
        else: st[r["start"]].append(k); err[r["start"]].append(100 * (N[k] / float(r["N_full"]) - 1))
for s in ("net_L0", "net_L120", "random"):
    say(f"5 % tol, start {s}: steps, fallbacks, max error", f"{min(st[s])}-{max(st[s])} steps, {fb[s]}/{ns[s]} fallback, max {max(err[s]):.2f} %", "results/cost_study/*_SS_v2.csv")
S = collections.defaultdict(list)
for r in csv.DictReader(open("results/cost_study/summary_tol5.csv")): S[r["method"]].append(float(r["online_s"]))
for m, v in S.items(): say(f"time per design: {m}", f"{np.mean(v):.2f} s", "results/cost_study/summary_tol5.csv")
src = open("scripts/cost_study.py").read(); say("fallback after N iterations", re.search(r"TOL_STOP, MAX_STEPS = (.*)", src).group(1), "scripts/cost_study.py")
# Paper 1 points
R = [json.loads(l) for l in open("results_p1/pinn_p1.jsonl")]
for m in ("physics_L0", "hybrid_L120"):
    v = [abs(r[m]["N_eta5"] / r["N_FEA"] - 1) * 100 for r in R]; say(f"22 points, {m} + check vs Abaqus max", f"{max(v):.2f} % (n={len(v)})", "results_p1/pinn_p1.jsonl")
say("22 points, Paper 1 PINN vs Abaqus max", f"{max(abs(r['N_PINN_P1']/r['N_FEA']-1)*100 for r in R):.2f} %", "Paper 1 Table 2 (copied)")
