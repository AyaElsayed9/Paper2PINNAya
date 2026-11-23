"""Summary of results/error_indicator/*.csv (scripts/error_indicator.py).

For every expert and network:
  * answers: head (instant), check after s steps; errors vs full GBT (truth only for scoring)
  * label-free estimates of the check error: eta_s^2 (relative residual squared), eta_s^2 / relgap_s and the last
    step's change; scored by the effectivity est / true error (median, 5-95 % range), the share of designs where the
    estimate is not below the true error (safe), and the Spearman rank correlation
  * adaptive stopping: steps until the estimate is <= tol (0.5 %, 1 %), the error then reached
  * GP (same labels): error and its own predictive std (coverage of +-2 std, Spearman)
Writes results/error_indicator/summary.csv and prints a table.
"""
import csv
import glob
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

STEPS = (0, 1, 3, 5, 12)
out = []
for f in sorted(glob.glob("results/error_indicator/*_SS_v2.csv")):
    tag = Path(f).stem.replace("_SS_v2", "")
    rows = list(csv.DictReader(open(f)))
    truth = {r["sid"]: float(r["N_full"]) for r in rows if r["N_full"]}
    models = sorted({r["model"] for r in rows}, key=lambda m: (m.startswith("gp"), m))
    for m in models:
        R = [r for r in rows if r["model"] == m and r["sid"] in truth]
        if not R:
            continue
        N = np.array([truth[r["sid"]] for r in R])
        head = 100 * np.abs(np.array([float(r["N_head"]) for r in R]) / N - 1)
        o = dict(expert=tag, model=m, n=len(R), head_mean=head.mean(), head_max=head.max())
        if m.startswith("gp"):
            sd = np.array([float(r["eta0"]) for r in R])               # std of log N
            err = np.log(np.array([float(r["N_head"]) for r in R]) / N)
            o.update(gp_cover2sd=np.mean(np.abs(err) <= 2 * sd), gp_spearman=spearmanr(sd, np.abs(err))[0],
                     gp_flag_recall=np.nan)
            big = np.abs(err) > 0.02
            if big.any():                                              # flag the 20 % largest std: share of >2 % errors caught
                thr = np.quantile(sd, 0.8); o["gp_flag_recall"] = np.mean(sd[big] >= thr)
            out.append(o); continue
        Ns = {s: np.array([float(r[f"N{s}"]) for r in R]) for s in STEPS}
        eta = {s: np.array([float(r[f"eta{s}"]) for r in R]) for s in STEPS}
        gap = {s: np.array([float(r[f"gap{s}"]) for r in R]) for s in STEPS}
        for s in (3, 5, 12):
            e = Ns[s] / N - 1
            o[f"s{s}_mean"] = 100 * np.abs(e).mean(); o[f"s{s}_max"] = 100 * np.abs(e).max()
            for name, est in (("eta2", eta[s] ** 2), ("eta2gap", eta[s] ** 2 / gap[s])):
                th = est / np.maximum(np.abs(e), 1e-7)
                o[f"s{s}_{name}_eff_med"] = np.median(th)
                o[f"s{s}_{name}_eff_p5"] = np.quantile(th, 0.05); o[f"s{s}_{name}_eff_p95"] = np.quantile(th, 0.95)
                o[f"s{s}_{name}_safe"] = np.mean(est >= np.abs(e))
                o[f"s{s}_{name}_spear"] = spearmanr(est, np.abs(e))[0]
        for tol in (0.005, 0.01):                                          # adaptive stop on eta^2 / relgap
            used, err = [], []
            for i in range(len(R)):
                st = next((s for s in STEPS if eta[s][i] ** 2 / gap[s][i] <= tol), STEPS[-1])
                used.append(st); err.append(100 * abs(Ns[st][i] / N[i] - 1))
            o[f"stop{tol}_steps_mean"] = np.mean(used); o[f"stop{tol}_err_mean"] = np.mean(err); o[f"stop{tol}_err_max"] = np.max(err)
        out.append(o)

keys = sorted({k for o in out for k in o}, key=lambda k: list(out[0].keys()).index(k) if k in out[0] else 999)
with open("results/error_indicator/summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
    for o in out:
        w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in o.items()})
for o in out:
    if o["model"].startswith("gp"):
        print(f"{o['expert']:18s} {o['model']:12s} head {o['head_mean']:5.2f}/{o['head_max']:5.2f} %  GP: +-2sd covers "
              f"{100 * o['gp_cover2sd']:3.0f} %, spearman(sd,|err|) {o['gp_spearman']:+.2f}, top-20 % sd catches "
              f"{100 * o['gp_flag_recall'] if o['gp_flag_recall'] == o['gp_flag_recall'] else float('nan'):3.0f} % of >2 % errors")
    else:
        print(f"{o['expert']:18s} {o['model']:12s} head {o['head_mean']:5.2f}/{o['head_max']:5.2f} % | 5 steps "
              f"{o['s5_mean']:.2f}/{o['s5_max']:.2f} %: eta2/gap eff {o['s5_eta2gap_eff_med']:.2f} "
              f"[{o['s5_eta2gap_eff_p5']:.2f}-{o['s5_eta2gap_eff_p95']:.2f}] safe {100 * o['s5_eta2gap_safe']:3.0f} % "
              f"rho {o['s5_eta2gap_spear']:+.2f} | eta2 eff {o['s5_eta2_eff_med']:.2f} safe {100 * o['s5_eta2_safe']:3.0f} % | "
              f"12 steps {o['s12_mean']:.2f}/{o['s12_max']:.2f}: eta2/gap eff {o['s12_eta2gap_eff_med']:.2f} safe "
              f"{100 * o['s12_eta2gap_safe']:3.0f} % | stop 1 %: {o['stop0.01_steps_mean']:.1f} steps -> "
              f"{o['stop0.01_err_mean']:.2f}/{o['stop0.01_err_max']:.2f} %; stop 0.5 %: {o['stop0.005_steps_mean']:.1f} steps -> "
              f"{o['stop0.005_err_mean']:.2f}/{o['stop0.005_err_max']:.2f} %")

# ---- pooled over all experts and both networks: the plain relative residual eta as a safe-side error estimate
print("\npooled (all experts, all networks):")
pool = []
for f in sorted(glob.glob("results/error_indicator/*_SS_v2.csv")):
    rows = list(csv.DictReader(open(f))); T = {r["sid"]: float(r["N_full"]) for r in rows if r["N_full"]}
    for r in rows:
        if r["model"].startswith("gp") or not r["N_full"]:
            continue
        pool.append(([float(r[f"N{s}"]) / T[r["sid"]] - 1 for s in STEPS], [float(r[f"eta{s}"]) for s in STEPS]))
E_ = np.array([p[0] for p in pool]); H_ = np.array([p[1] for p in pool])
for j, s in enumerate(STEPS):
    e, eta = E_[:, j], H_[:, j]
    print(f"  step {s:2d}: true error median {100 * np.median(e):5.2f} % max {100 * e.max():5.2f} % | eta median "
          f"{100 * np.median(eta):5.2f} % | eta >= error in {100 * np.mean(eta >= np.abs(e)):.1f} %, eta / error median "
          f"{np.median(eta / np.abs(e)):.1f}, spearman {spearmanr(eta, np.abs(e))[0]:+.2f} | eta^2 >= error in "
          f"{100 * np.mean(eta ** 2 >= np.abs(e)):.1f} %")
for tol in (0.10, 0.05, 0.03):
    k = [next((j for j in range(len(STEPS)) if H_[i, j] <= tol), len(STEPS) - 1) for i in range(len(pool))]
    err = np.array([100 * E_[i, kk] for i, kk in enumerate(k)]); st = np.array([STEPS[kk] for kk in k])
    print(f"  stop when eta <= {100 * tol:.0f} %: steps mean {st.mean():.1f}, true error mean {err.mean():.2f} / max "
          f"{err.max():.2f} %, error <= eta for all: {bool(np.all(err <= 100 * H_[range(len(pool)), k]))}")
