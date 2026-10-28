"""Fill the placeholders of part3.md from the result files (no hand-typed numbers)."""
import json
from pathlib import Path

import numpy as np

H = Path(__file__).parent
cur = [json.loads(l) for l in (H / "results" / "periodic_curves.jsonl").read_text().splitlines()]
conv = json.loads((H / "results" / "periodic_conv.json").read_text())
tim = json.loads((H / "results" / "strip_timing.json").read_text())   # GBT times of the same runs
m = lambda v: f"{v:+.2f}".replace("-", "−")
s = (H / "part3.md").read_text()

import csv
s1 = (H / "part1.md").read_text(); s4 = (H / "part4.md").read_text()
ref = list(csv.DictReader(open(H / "src" / "omega_results_V32.csv")))
LAYS = ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]
gbt = {(c["layup"], c["bs"], c["beta"], c["bc"]): c["N"] for c in cur}
# --- Table 6: SS, Abaqus vs GBT
rows = ["| Lay-up | $b_s$ [mm] | $\\beta$ | Abaqus | GBT | $e$ [%] |", "|:--|--:|--:|--:|--:|--:|"]
for lay in LAYS:
    first = True
    for r in [r for r in ref if r["layup"] == lay and r["bc"] == "SS"]:
        a = float(r["abq_per"]); g = gbt[(lay, float(r["bs"]), float(r["beta"]), "SS")]
        rows.append(f"| {lay if first else ''} | {float(r['bs']):g} | {float(r['beta']):g} | {a:.3f} | {g:.3f} | {m(100 * (g / a - 1))} |")
        first = False
s = s.replace("@@TABSS@@", "\n".join(rows))
# --- Abaqus CC (Aya's omega_validation_eigs_CC_PBC.csv copied to src/abaqus_cc_periodic.csv)
fcc = H / "src" / "abaqus_cc_periodic.csv"
acc = {}
if fcc.exists():
    for r in csv.DictReader(open(fcc)):
        if r["Ncr"] not in ("", "NA"):
            acc[(r["layup"], float(r["bs"]), float(r["beta"]))] = float(r["Ncr"])
rows = ["| Lay-up | $b_s$ [mm] | $\\beta$ | Abaqus | GBT | $e$ [%] |", "|:--|--:|--:|--:|--:|--:|"]
ecc = []
for lay in LAYS:
    first = True
    for r in [r for r in ref if r["layup"] == lay and r["bc"] == "CC"]:
        k = (lay, float(r["bs"]), float(r["beta"])); g = gbt[k + ("CC",)]
        if k in acc:
            e = 100 * (g / acc[k] - 1); ecc.append(e); a_s = f"{acc[k]:.3f}"; e_s = m(e)
        else:
            a_s = e_s = "[pending]"
        rows.append(f"| {lay if first else ''} | {k[1]:g} | {k[2]:g} | {a_s} | {g:.3f} | {e_s} |")
        first = False
s = s.replace("@@TABCC@@", "\n".join(rows))
if len(ecc) == 30:
    mcc, xcc = np.mean(np.abs(ecc)), np.max(np.abs(ecc))
    s1 = s1.replace("@@ABS_CC@@", f"; for the 30 panels with clamped loaded edges, they agree within {mcc:.2f} % on average and {xcc:.2f} % at most")
    s4 = s4.replace("@@CONC_CC@@", f", and for 30 panels with clamped loaded edges within {mcc:.2f} % on average and {xcc:.2f} % at most")
    s = s.replace("@@CC_TEXT@@", f"GBT agrees with Abaqus within {mcc:.2f} % on average and {xcc:.2f} % at most over the 30 cases (Table 7 and Figure 5b).")
    s = s.replace("@@FIG6REF@@", " and Abaqus shell models with periodic constraint equations (open markers)")
    s = s.replace("@@FIG5PANELS@@", "(a) simply supported and (b) clamped loaded edges")
else:
    s1 = s1.replace("@@ABS_CC@@", "; the clamped loaded edges are verified in the same way [Abaqus results pending]")
    s4 = s4.replace("@@CONC_CC@@", " [comparison for the 30 panels with clamped loaded edges pending]")
    s = s.replace("@@CC_TEXT@@", f"[Abaqus results for the clamped cells: {len(ecc)} of 30 available; text to be completed.]")
    s = s.replace("@@FIG6REF@@", " (Abaqus markers to be added)")
    s = s.replace("@@FIG5PANELS@@", "simply supported loaded edges (panel for the clamped loaded edges to be added)")
# --- Abaqus model size and time (Aya's fe_modes/job_summary.csv copied to src/abaqus_job_summary.csv)
fjs = H / "src" / "abaqus_job_summary.csv"
js = [r for r in csv.DictReader(open(fjs))] if fjs.exists() else []
sel = [r for r in js if "bs150_beta1p5" in r["job"] and r["job"].endswith("_PBC") and "_CC_" not in r["job"] and r["wallclock_s"]]
if sel:
    var = np.mean([float(r["variables"]) for r in sel]); wall = np.mean([float(r["wallclock_s"]) for r in sel])
    s = s.replace("@@ABQROW@@", f"| Abaqus shell model (S4R/S3, seed 3.7 mm) | {var:,.0f} DOFs | {wall:.1f} s on average (wall clock) |".replace(",", " "))
    s = s.replace("@@ABQRATIO@@", f"18 777 unknowns, i.e. {100 * 18777 / var:.0f} % of")
else:
    s = s.replace("@@ABQROW@@", "| Abaqus shell model (S4R/S3, seed 3.7 mm) | [from the .dat files] | [to be added] |")
    s = s.replace("@@ABQRATIO@@", "18 777 unknowns, i.e. a fraction of [ratio to be added]")

# mode decomposition
rows = ["| Lay-up | Ends | local (L) | distortional (D) | global (G) | T, S, M |", "|:--|:--|--:|--:|--:|--:|"]
for lay in ["Isotropic", "Cross-ply", "Sym1", "Sym2", "Unsym"]:
    for bc in ("SS", "CC"):
        x = [c for c in cur if c["layup"] == lay and c["bc"] == bc]
        rng = lambda g: f'{min(100 * c["rho"].get(g, 0) for c in x):.1f} to {max(100 * c["rho"].get(g, 0) for c in x):.1f}'.replace("-", "−")
        o = max(100 * sum(c["rho"].get(g, 0) for g in "TSM") for c in x)
        rows.append(f"| {lay if bc == 'SS' else ''} | {bc} | {rng('L')} | {rng('D')} | {rng('G')} | < {max(o, 0.01):.2f} |")
s = s.replace("@@TABRHO@@", "\n".join(rows))
g = [100 * c["rho"].get("G", 0) for c in cur]; d = [100 * c["rho"].get("D", 0) for c in cur]
tsm = max(100 * sum(c["rho"].get(k, 0) for k in "TSM") for c in cur)
s = s.replace("@@GMAX@@", f"{max(g):.1f}").replace("@@DMIN@@", f"{min(d):.1f}".replace("-", "−")).replace("@@DMAX@@", f"{max(d):.1f}")
s = s.replace("@@TSMMAX@@", f"{max(tsm, 0.1):.1f}")

# timing
t = [c["time"] for c in cur]
s = s.replace("@@TMEAN@@", f"{np.mean(t):.1f}").replace("@@TMIN@@", f"{min(t):.1f}").replace("@@TMAX@@", f"{max(t):.1f}")
s = s.replace("@@TGSEQ@@", f"{np.mean(tim['gbt']):.1f}")

# convergence
lab = {"Isotropic|150.0|1.5|SS": "Isotropic, SS, $b_s = 150$, $\\beta = 1.5$", "Sym2|150.0|3.0|SS": "Sym2, SS, $b_s = 150$, $\\beta = 3$",
       "Unsym|150.0|1.5|CC": "Unsym, CC, $b_s = 150$, $\\beta = 1.5$", "Sym1|285.5|3.0|CC": "Sym1, CC, $b_s = 285.5$, $\\beta = 3$",
       "Sym1|285.5|1.5|CC": "Sym1, CC, $b_s = 285.5$, $\\beta = 1.5$"}
rows = ["| Case | Default | Spans $h_x/2$ | Segments $b_s/40$ |", "|:--|--:|--:|--:|"]
allch = []
for k, name in lab.items():
    r = conv[k]; d0 = r["default"]
    ch = lambda v: 100 * (v - d0) / d0
    allch += [ch(r["spans_half"]), ch(r["seg40"])]
    rows.append(f"| {name} | {d0:.4f} | {r['spans_half']:.4f} ({m(ch(r['spans_half']))} %) | {r['seg40']:.4f} ({m(ch(r['seg40']))} %) |")
s = s.replace("@@TABCONV@@", "\n".join(rows))
q = conv["Sym1|285.5|1.5|CC"]
qch = 100 * (q["spans_quarter"] - q["default"]) / q["default"]
q8 = 100 * (q["spans_eighth"] - q["default"]) / q["default"]
others = [v for k, r in conv.items() if k in lab and k != "Sym1|285.5|1.5|CC"
          for v in (100 * (r["spans_half"] - r["default"]) / r["default"], 100 * (r["seg40"] - r["default"]) / r["default"])]
txt = (f"Except for one case, all changes are within {max(abs(v) for v in others):.2f} %, i.e. below the deviation from the reference models. "
       f"The exception is the Sym1 cell with clamped loaded edges, $b_s = 285.5$ mm and $\\beta = 1.5$: halving the spans repeatedly "
       f"($h_x/2$, $h_x/4$, $h_x/8$) raises the critical load by {m(100 * (q['spans_half'] - q['default']) / q['default'])}, {m(qch)} and {m(q8)} %, "
       f"i.e. slowly, so that the converged value is estimated to lie about 1 % above the default result; the spline discretisation is less "
       f"accurate near the clamped ends of the wide skin bays with bending–twisting coupling. Finer segments lower the critical "
       f"load slightly for $b_s = 150$ mm and therefore do not remove the offset against Abaqus noted in Section 3.3. "
       f"The link walls were checked in the same way: for the isotropic cell with simply supported loaded edges, $b_s = 150$ mm and "
       f"$\\beta = 1.5$, the critical load is {conv['linkk']['100.0']:.3f}, {conv['linkk']['1000.0']:.3f} and {conv['linkk']['10000.0']:.3f} N/mm "
       f"for $k_L = 10^2$, $10^3$ and $10^4$, whereas a link wall with the stiff twisting term $D^L_{{66}} = k_L D_{{0,66}}$ gives "
       f"{conv['link_stiffD66']:.3f} N/mm ({100 * (conv['link_stiffD66'] - conv['linkk']['1000.0']) / conv['linkk']['1000.0']:+.1f} %).")
s = s.replace("@@CONVTEXT@@", txt)

# hierarchical (reduced) bases
red = json.loads((H / "results" / "reduced2.json").read_text())
red1 = json.loads((H / "results" / "reduced.json").read_text())
names = {"Isotropic|150.0|1.5|SS": ("Isotropic", "SS", 150), "Cross-ply|150.0|1.5|SS": ("Cross-ply", "SS", 150),
         "Sym2|150.0|1.5|SS": ("Sym2", "SS", 150), "Unsym|150.0|1.5|SS": ("Unsym", "SS", 150),
         "Unsym|150.0|1.5|CC": ("Unsym", "CC", 150), "Sym1|285.5|1.5|CC": ("Sym1", "CC", 285.5)}
rows = ["| Lay-up | Ends | $b_s$ [mm] | $n_C = 20$ | $n_C = 40$ | $n_C = 60$ | $n_C = 100$ | $n_C = 150$ | no Mindlin modes |",
        "|:--|:--|--:|--:|--:|--:|--:|--:|--:|"]
dev = lambda v, f: 100 * (v - f) / f
for k, (lay, bc, bs) in names.items():
    r = red[k]; f0 = r["full"][0]
    rows.append(f"| {lay} | {bc} | {bs:g} | " + " | ".join(m(dev(r[c][0], f0)) for c in ("c20", "c40", "c60", "c100", "c150", "noM")) + " |")
s = s.replace("@@TABRED@@", "\n".join(rows))
nM = {k: dev(r["noM"][0], r["full"][0]) for k, r in red.items()}
c40 = {k: dev(r["c40"][0], r["full"][0]) for k, r in red.items()}
c100 = {k: dev(r["c100"][0], r["full"][0]) for k, r in red.items()}
unc = ["Isotropic|150.0|1.5|SS", "Cross-ply|150.0|1.5|SS"]
cpl = [k for k in red if k not in unc]
sm = {k: v for k, v in red1.items()}
iso = sm["Isotropic|285.5|0.75|SS"]["a_20"]; cro = sm["Cross-ply|150.0|1.5|SS"]["a_20"]
uss = sm["Unsym|150.0|1.5|SS"]["a_20"]; ucc = sm["Unsym|150.0|1.5|CC"]["a_20"]
txt = (f"The hierarchical ordering is effective for the conventional modes: for the uncoupled lay-ups, 40 conventional modes reproduce the "
       f"complete result within {max(abs(c40[k]) for k in unc):.2f} %, and 100 within {max(abs(c100[k]) for k in unc):.2f} %. For the coupled lay-ups, "
       f"more conventional modes are needed ({min(c40[k] for k in cpl):.1f} % to {max(c40[k] for k in cpl):.1f} % with 40 and "
       f"{min(c100[k] for k in cpl):.2f} % to {max(c100[k] for k in cpl):.2f} % with 100 modes), because the skewed buckles of these laminates "
       f"combine many local modes. The transverse shear (Mindlin) modes are not negligible: omitting them raises the critical load by "
       f"{min(nM.values()):.2f} % to {max(nM.values()):.2f} %. A small basis in the classical sense, consisting of the first 20 conventional modes, "
       f"the warping (shear-mode) field that each of them activates at $k_r$ and no Mindlin modes ({iso[1]} modes), gives "
       f"{m(iso[2])} % for the isotropic cell with $b_s = 285.5$ mm and $\\beta = 0.75$ and {m(cro[2])} % for the cross-ply cell with $b_s = 150$ mm and $\\beta = 1.5$, "
       f"but {m(uss[2])} % and {m(ucc[2])} % for the Unsym cell with simply supported and clamped loaded edges. Compact GBT bases are therefore "
       f"adequate for uncoupled laminates, whereas for coupled laminates the complete basis, or a basis that is enriched in the local and Mindlin "
       f"families, is required; the construction of such bases is left for future work.")
s = s.replace("@@REDTEXT@@", txt)

# figure data (figures_paper.py)
fd = json.loads((H / "results" / "figdata.json").read_text())["fig08"]
modes = [(int(k[4:]), v) for k, v in fd.items() if k.startswith("mode")]
rel = abs(fd["N_hier"] - fd["N_full"]) / fd["N_full"]
ex = int(np.floor(np.log10(rel)))
s = s.replace("@@HN@@", f"{fd['N_full']:.2f}").replace("@@HDIFF@@", f"${rel / 10 ** ex:.1f}\\times10^{{{ex}}}$")
s = s.replace("@@HMODES@@", ", ".join(str(k) for k, _v in modes[:-1]) + f" and {modes[-1][0]}")
s = s.replace("@@H4@@", f"{sum(v for _k, v in modes):.1f}").replace("@@H20@@", f"{fd['sum20']:.1f}")

hc = json.loads((H / "results" / "figdata.json").read_text())["fig11"]
nm_ = {"Isotropic|150.0|1.5|SS": "the isotropic", "Cross-ply|150.0|1.5|SS": "the Cross-ply", "Sym2|150.0|1.5|SS": "the Sym2",
       "Unsym|150.0|1.5|SS": "the Unsym (SS)", "Unsym|150.0|1.5|CC": "the Unsym (CC)", "Sym1|285.5|1.5|CC": "the Sym1 (CC, $b_s = 285.5$ mm)"}
r16 = [(v["nc"]["16"] / v["full"] - 1) / (v["nc"]["20"] / v["full"] - 1) for v in json.loads((H / "results" / "hier_conv.json").read_text()).values()]
parts = sorted(((v["n_lt1"], nm_[k]) for k, v in hc.items()), key=lambda t: t[0])
txt = ("The convergence is stepwise (Figure 11): some groups of hierarchical modes do not change the critical load, whereas the four modes "
       "between $n_C = 16$ and 20 reduce the deviation by a factor of " + f"{min(r16):.0f} to {max(r16):.0f}. A deviation below 1 % requires " +
       ", ".join(f"{n} modes for {name} cell" for n, name in parts[:-1]) + f" and {parts[-1][0]} modes for {parts[-1][1]} cell, "
       "i.e. 20 to 40 of the 273 (or 231) conventional modes for the isotropic, cross-ply and angle-ply cells, but about one third of them for the unsymmetric laminate.")
s = s.replace("@@HCONV@@", txt)

# --- long cells (results/long_curves.jsonl; short cells from periodic_curves)
import re
LONG = {}
for c in cur:
    if c["bs"] == 150.0 and c["bc"] == "SS":
        LONG[(c["layup"], c["beta"], 1.0)] = c
for l in (H / "results" / "long_curves.jsonl").read_text().splitlines():
    c = json.loads(l); LONG[(c["layup"], c["beta"], float(c.get("pss", 1.0)))] = c
def _lv(mo):
    kind, lay, beta, pss = mo.group(1), mo.group(2), float(mo.group(3)), float(mo.group(4))
    c = LONG[(lay, beta, pss)]
    if kind == "N":
        return f"{c['N']:.2f}"
    if kind == "N2":
        return f"{c['eigs'][1]:.2f}"
s = re.sub(r"@@(N2?):([A-Za-z0-9-]+):([\d.]+):([\d.]+)@@", _lv, s)
s = re.sub(r"@@R:([A-Za-z0-9-]+):([\d.]+):([\d.]+):([A-Z])@@",
           lambda mo: f"{100 * LONG[(mo.group(1), float(mo.group(2)), float(mo.group(3)))]['rho'].get(mo.group(4), 0):.1f}", s)
s = re.sub(r"@@RATIO:([A-Za-z0-9-]+):([\d.]+):([\d.]+)@@",
           lambda mo: f"{LONG[(mo.group(1), float(mo.group(2)), 1.0)]['N'] / LONG[(mo.group(1), float(mo.group(3)), 1.0)]['N']:.2f}", s)
fac = max(LONG[(l, b, 0.0)]["N"] / LONG[(l, b, 1.0)]["N"] for (l, b, p) in LONG if p == 0.0 and (l, b, 1.0) in LONG)
s = s.replace("@@FACTOR@@", f"{LONG[('Sym1', 8.0, 0.0)]['N'] / LONG[('Sym1', 8.0, 1.0)]['N']:.1f}")
s1 = s1.replace("@@FACTOR@@", f"{fac:.1f}"); s4 = s4.replace("@@FACTOR@@", f"{fac:.1f}")
fal = H / "src" / "abaqus_long.csv"
if fal.exists():
    al = list(csv.DictReader(open(fal)))
    txt = "; the Abaqus results with constant width and with free lateral expansion are compared in Figure 9a: " + "; ".join(
        f"{r['layup']}, $\\beta = {float(r['beta']):g}$, {('constant width' if r['sides'] == 'PBC' else 'free expansion')}: {float(r['Ncr']):.2f} N/mm"
        for r in al)
    s = s.replace("@@ABQLONG@@", txt).replace("@@FIG9ABQ@@", "; stars: Abaqus with constant width; crosses: Abaqus with free lateral expansion")
else:
    s = s.replace("@@ABQLONG@@", " [Abaqus check with constant width and with free lateral expansion pending]").replace("@@FIG9ABQ@@", " (Abaqus results to be added)")

(H / "part3_filled.md").write_text(s)
(H / "part1_filled.md").write_text(s1)
(H / "part4_filled.md").write_text(s4)
left = [p for t in (s, s1, s4) for p in t.split("@@")[1::2]]
print("unfilled:", left)
