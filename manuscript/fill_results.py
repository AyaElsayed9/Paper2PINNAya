"""Fill the placeholders of part1/part3/part4 from the result files (no hand-typed numbers)."""
import csv
import json
import re
from pathlib import Path

import numpy as np

H = Path(__file__).parent
cur = [json.loads(l) for l in (H / "results" / "periodic_curves.jsonl").read_text().splitlines()]
conv = json.loads((H / "results" / "periodic_conv.json").read_text())
tim = json.loads((H / "results" / "strip_timing.json").read_text())      # GBT times of the sequential runs
FD = json.loads((H / "results" / "figdata.json").read_text())             # written by figures_paper.py
m = lambda v: f"{v:+.2f}".replace("-", "−")
s = (H / "part3.md").read_text(); s1 = (H / "part1.md").read_text(); s4 = (H / "part4.md").read_text()
LAYS = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
gbt = {(c["layup"], c["bs"], c["beta"], c["bc"]): c["N"] for c in cur}
ref = list(csv.DictReader(open(H / "src" / "omega_results_V32.csv")))

# ---------------- Abaqus: SS from the results table, CC from Aya's CC run (when available)
abq = {}
for r in ref:
    if r["bc"] == "SS":
        abq[(r["layup"], float(r["bs"]), float(r["beta"]), "SS")] = float(r["abq_per"])
fcc = H / "src" / "abaqus_cc_periodic.csv"                # = omega_validation_eigs_CC_PBC.csv
if fcc.exists():
    for r in csv.DictReader(open(fcc)):
        if r["Ncr"] not in ("", "NA"):
            abq[(r["layup"], float(r["bs"]), float(r["beta"]), "CC")] = float(r["Ncr"])

# ---------------- Table 5: SS and CC side by side
rows = ["| Lay-up | $b_s$ [mm] | $\\beta$ | SS: Abaqus | SS: GBT | SS: $e$ [%] | CC: Abaqus | CC: GBT | CC: $e$ [%] |",
        "|:--|--:|--:|--:|--:|--:|--:|--:|--:|"]
ecc, ess = [], []
for lay in LAYS:
    first = True
    for bs in (150.0, 285.5):
        for beta in (0.75, 1.5, 3.0):
            cells = []
            for bc in ("SS", "CC"):
                g = gbt[(lay, bs, beta, bc)]; a = abq.get((lay, bs, beta, bc))
                if a is None:
                    cells += ["[pending]", f"{g:.3f}", "[pending]"]
                else:
                    e = 100 * (g / a - 1); (ess if bc == "SS" else ecc).append(e)
                    cells += [f"{a:.3f}", f"{g:.3f}", m(e)]
            rows.append(f"| {lay if first else ''} | {bs:g} | {beta:g} | " + " | ".join(cells) + " |")
            first = False
s = s.replace("@@TABNCR@@", "\n".join(rows))
assert len(ess) == 30
ssd = {}
for lay in LAYS:
    for bs in (150.0, 285.5):
        for beta in (0.75, 1.5, 3.0):
            ssd[(lay, bs, beta)] = 100 * (gbt[(lay, bs, beta, "SS")] / abq[(lay, bs, beta, "SS")] - 1)
laym = [np.mean([abs(v) for k, v in ssd.items() if k[0] == l]) for l in LAYS]
pit = {bs: [abs(v) for k, v in ssd.items() if k[1] == bs] for bs in (150.0, 285.5)}
assert all(v < 0 for k, v in ssd.items() if k[1] == 150.0)          # "all 15 GBT results below Abaqus"
rep = {"@@SSMEAN@@": f"{np.mean(np.abs(ess)):.2f}", "@@SSMAX@@": f"{np.max(np.abs(ess)):.2f}",
       "@@LAYMIN@@": f"{min(laym):.2f}", "@@LAYMAX@@": f"{max(laym):.2f}",
       "@@M285@@": f"{np.mean(pit[285.5]):.2f}", "@@X285@@": f"{max(pit[285.5]):.2f}",
       "@@M150@@": f"{np.mean(pit[150.0]):.2f}", "@@X150@@": f"{max(pit[150.0]):.2f}"}
for k_, v_ in rep.items():
    s = s.replace(k_, v_); s1 = s1.replace(k_, v_); s4 = s4.replace(k_, v_)
if len(ecc) == 30:
    mcc, xcc = np.mean(np.abs(ecc)), np.max(np.abs(ecc))
    s1 = s1.replace("@@ABS_CC@@", f"; with clamped loaded edges within {mcc:.2f} % and {xcc:.2f} %")
    s4 = s4.replace("@@CONC_CC@@", f", and with clamped loaded edges within {mcc:.2f} % and {xcc:.2f} %")
    s = s.replace("@@CC_TEXT@@", f"With CC loaded edges the deviations are {mcc:.2f} % on average and {xcc:.2f} % at most.")
    s = s.replace("@@FIG4CC@@", "").replace("@@FIG5PANELS@@", "(a) simply supported and (b) clamped loaded edges")
else:
    s1 = s1.replace("@@ABS_CC@@", " [clamped loaded edges: Abaqus results pending]")
    s4 = s4.replace("@@CONC_CC@@", " [clamped loaded edges: Abaqus comparison pending]")
    s = s.replace("@@CC_TEXT@@", f"[CC loaded edges: {len(ecc)} of 30 Abaqus results available; to be completed.]")
    s = s.replace("@@FIG4CC@@", " (CC: Abaqus results to be added)")
    s = s.replace("@@FIG5PANELS@@", "simply supported loaded edges (clamped loaded edges to be added)")

# ---------------- Table 7: Abaqus size and time (Aya's fe_modes/job_summary.csv)
fjs = H / "src" / "abaqus_job_summary.csv"
js = list(csv.DictReader(open(fjs))) if fjs.exists() else []
sel = [r for r in js if "bs150_beta1p5" in r["job"] and r["job"].endswith("_PBC") and "_CC_" not in r["job"] and r["wallclock_s"]]
if sel:
    var = np.mean([float(r["variables"]) for r in sel]); wall = np.mean([float(r["wallclock_s"]) for r in sel])
    s = s.replace("@@ABQROW@@", f"| Abaqus shell model (S4R/S3, seed 3.7 mm) | {var:,.0f} DOFs | {wall:.1f} s on average (wall clock) |".replace(",", " "))
    s = s.replace("@@ABQRATIO@@", f"18 777 unknowns, i.e. {100 * 18777 / var:.0f} % of")
else:
    s = s.replace("@@ABQROW@@", "| Abaqus shell model (S4R/S3, seed 3.7 mm) | [from the .dat files] | [to be added] |")
    s = s.replace("@@ABQRATIO@@", "18 777 unknowns, i.e. a fraction of [ratio to be added]")
s = s.replace("@@TGSEQ@@", f"{np.mean(tim['gbt']):.1f}")

# ---------------- mode decomposition of the short cells
g = [100 * c["rho"].get("G", 0) for c in cur]; d = [100 * c["rho"].get("D", 0) for c in cur]
L_ = [100 * c["rho"].get("L", 0) for c in cur]
tsm = max(100 * sum(c["rho"].get(k, 0) for k in "TS") for c in cur)
assert 95.5 < min(L_) and max(L_) < 101.5, (min(L_), max(L_))       # "96-101 %" in the text
s = s.replace("@@GMAX@@", f"{max(g):.1f}").replace("@@DMIN@@", f"{min(d):.1f}".replace("-", "−")).replace("@@DMAX@@", f"{max(d):.1f}")
s = s.replace("@@TSMMAX@@", f"{max(tsm, 0.1):.1f}")
fd = FD["fig07"]                                                       # hierarchical participation (Figure 7d-f)
modes = [(int(k[4:]), v) for k, v in fd.items() if k.startswith("mode")]
rel = abs(fd["N_hier"] - fd["N_full"]) / fd["N_full"]; ex = int(np.floor(np.log10(rel)))
s = s.replace("@@HN@@", f"{fd['N_full']:.2f}").replace("@@HDIFF@@", f"${rel / 10 ** ex:.1f}\\times10^{{{ex}}}$")
s = s.replace("@@HMODES@@", ", ".join(str(k) for k, _v in modes[:-1]) + f" and {modes[-1][0]}")
s = s.replace("@@H4@@", f"{sum(v for _k, v in modes):.1f}").replace("@@H20@@", f"{fd['sum20']:.1f}")

# ---------------- hierarchical bases (Figure 10 data + Mindlin omission + classical small basis)
hc = json.loads((H / "results" / "hier_conv.json").read_text())
red = json.loads((H / "results" / "reduced2.json").read_text()); red1 = json.loads((H / "results" / "reduced.json").read_text())
n1 = {k: next((n for n in sorted(map(int, v["nc"])) if v["nc"][str(n)] / v["full"] - 1 < 0.01), None) for k, v in hc.items()}
r16 = [(v["nc"]["16"] / v["full"] - 1) / (v["nc"]["20"] / v["full"] - 1) for v in hc.values()]
unc = [k for k in hc if k.split("|")[0] in ("Isotropic", "Cross-ply", "Sym1", "Sym2")]
uns = [k for k in hc if k.split("|")[0] == "Unsym"]
nM = [100 * (r["noM"][0] / r["full"][0] - 1) for r in red.values()]
iso = red1["Isotropic|285.5|0.75|SS"]["a_20"]; uss = red1["Unsym|150.0|1.5|SS"]["a_20"]
txt = (f"The convergence is stepwise: some groups of modes leave the critical load unchanged, whereas the four modes from $n_C = 17$ to 20 "
       f"reduce the deviation by a factor of {min(r16):.0f} to {max(r16):.0f}. A deviation below 1 % requires {min(n1[k] for k in unc)} to "
       f"{max(n1[k] for k in unc)} modes for the isotropic, cross-ply and angle-ply cells, but {min(n1[k] for k in uns)} to {max(n1[k] for k in uns)} "
       f"of the 273 modes for the unsymmetric cell, whose skewed buckles combine many local modes. The Mindlin modes are not negligible: "
       f"omitting them raises the critical load by {min(nM):.2f} % to {max(nM):.2f} %. A small basis in the classical sense (20 conventional "
       f"modes, the warping each of them activates at $k_r$, no Mindlin modes; {iso[1]} modes) gives {m(iso[2])} % for the isotropic cell "
       f"with $b_s = 285.5$ mm and $\\beta = 0.75$, but {m(uss[2])} % for the Unsym cell with SS edges.")
s = s.replace("@@HCONV@@", txt)

# ---------------- discretisation check (Table 6)
lab = {"Isotropic|150.0|1.5|SS": "Isotropic, SS, $b_s = 150$, $\\beta = 1.5$", "Sym2|150.0|3.0|SS": "Sym2, SS, $b_s = 150$, $\\beta = 3$",
       "Unsym|150.0|1.5|CC": "Unsym, CC, $b_s = 150$, $\\beta = 1.5$", "Sym1|285.5|3.0|CC": "Sym1, CC, $b_s = 285.5$, $\\beta = 3$",
       "Sym1|285.5|1.5|CC": "Sym1, CC, $b_s = 285.5$, $\\beta = 1.5$"}
rows = ["| Case | Default | Spans $h_x/2$ | Segments $b_s/40$ |", "|:--|--:|--:|--:|"]
for k, name in lab.items():
    r = conv[k]; d0 = r["default"]; ch = lambda v: 100 * (v - d0) / d0
    rows.append(f"| {name} | {d0:.4f} | {r['spans_half']:.4f} ({m(ch(r['spans_half']))} %) | {r['seg40']:.4f} ({m(ch(r['seg40']))} %) |")
s = s.replace("@@TABCONV@@", "\n".join(rows))
q = conv["Sym1|285.5|1.5|CC"]; qc = lambda key: 100 * (q[key] - q["default"]) / q["default"]
others = [v for k, r in conv.items() if k in lab and k != "Sym1|285.5|1.5|CC"
          for v in (100 * (r["spans_half"] - r["default"]) / r["default"], 100 * (r["seg40"] - r["default"]) / r["default"])]
lk = conv["linkk"]
txt = (f"All changes are within {max(abs(v) for v in others):.2f} %, except for the Sym1 cell with CC edges, $b_s = 285.5$ mm and $\\beta = 1.5$, "
       f"for which halving the spans repeatedly raises the critical load by {m(qc('spans_half'))}, {m(qc('spans_quarter'))} and {m(qc('spans_eighth'))} %, "
       f"so that the converged value lies about 1 % above the default; finer segments lower the critical load slightly for $b_s = 150$ mm. "
       f"The critical load does not depend on the link-wall factor ({lk['100.0']:.3f}, {lk['1000.0']:.3f} and {lk['10000.0']:.3f} N/mm for "
       f"$k_L = 10^2$, $10^3$, $10^4$; isotropic cell, SS, $b_s = 150$ mm, $\\beta = 1.5$), whereas a stiff twisting term "
       f"$D^L_{{66}} = k_L D_{{0,66}}$ would lock the web roots ({100 * (conv['link_stiffD66'] - lk['1000.0']) / lk['1000.0']:+.1f} %).")
s = s.replace("@@CONVTEXT@@", txt)

# ---------------- long cells
LONG = {}
for c in cur:
    if c["bs"] == 150.0 and c["bc"] == "SS":
        LONG[(c["layup"], c["beta"], 1.0)] = c
for l in (H / "results" / "long_curves.jsonl").read_text().splitlines():
    c = json.loads(l); LONG[(c["layup"], c["beta"], float(c.get("pss", 1.0)))] = c
def _lv(mo):
    c = LONG[(mo.group(2), float(mo.group(3)), float(mo.group(4)))]
    return f"{c['N']:.2f}" if mo.group(1) == "N" else f"{c['eigs'][1]:.2f}"
for t_ in ("s", "s1", "s4"):
    pass
s = re.sub(r"@@(N2?):([A-Za-z0-9-]+):([\d.]+):([\d.]+)@@", _lv, s)
s = re.sub(r"@@R:([A-Za-z0-9-]+):([\d.]+):([\d.]+):([A-Z])@@",
           lambda mo: f"{100 * LONG[(mo.group(1), float(mo.group(2)), float(mo.group(3)))]['rho'].get(mo.group(4), 0):.1f}", s)
s = re.sub(r"@@RATIO:([A-Za-z0-9-]+):([\d.]+):([\d.]+)@@",
           lambda mo: f"{LONG[(mo.group(1), float(mo.group(2)), 1.0)]['N'] / LONG[(mo.group(1), float(mo.group(3)), 1.0)]['N']:.2f}", s)
facs = {(l, b): LONG[(l, b, 0.0)]["N"] / LONG[(l, b, 1.0)]["N"] for (l, b, p) in LONG if p == 0.0 and (l, b, 1.0) in LONG}
(fl, fb), fac = max(facs.items(), key=lambda t: t[1])
for name in ("s", "s1", "s4"):
    globals()[name] = globals()[name].replace("@@FACTOR@@", f"{fac:.1f}").replace("@@FACTORCASE@@", f"{fl}, $\\beta = {fb:g}$")
fal = H / "src" / "abaqus_long.csv"
if fal.exists():
    al = list(csv.DictReader(open(fal)))
    txt = "; the Abaqus results (Figure 8a) give " + "; ".join(
        f"{r['layup']}, $\\beta = {float(r['beta']):g}$, {('constant width' if r['sides'] == 'PBC' else 'free expansion')}: {float(r['Ncr']):.2f} N/mm"
        for r in al)
    s = s.replace("@@ABQLONG@@", txt).replace("@@FIG9ABQ@@", "; stars: Abaqus, constant width; crosses: Abaqus, free lateral expansion")
else:
    s = s.replace("@@ABQLONG@@", " [Abaqus check with constant width and free lateral expansion pending]").replace("@@FIG9ABQ@@", " (Abaqus results to be added)")

(H / "part3_filled.md").write_text(s); (H / "part1_filled.md").write_text(s1); (H / "part4_filled.md").write_text(s4)
left = [p for t in (s, s1, s4) for p in t.split("@@")[1::2]]
print("unfilled:", left)
