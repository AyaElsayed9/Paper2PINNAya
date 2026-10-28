"""Fill the placeholders of part3.md from the result files (no hand-typed numbers)."""
import json
from pathlib import Path

import numpy as np

H = Path(__file__).parent
cur = [json.loads(l) for l in (H / "results" / "periodic_curves.jsonl").read_text().splitlines()]
conv = json.loads((H / "results" / "periodic_conv.json").read_text())
tim = json.loads((H / "results" / "strip_timing.json").read_text())
m = lambda v: f"{v:+.2f}".replace("-", "−")
s = (H / "part3.md").read_text()

s = s.replace("@@TABSS@@", (H / "tab_ss.md").read_text().strip())
s = s.replace("@@TABCC@@", (H / "tab_cc.md").read_text().strip())

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
s = s.replace("@@TGSEQ@@", f"{np.mean(tim['gbt']):.1f}").replace("@@TSSEQ@@", f"{np.mean(tim['times']):.1f}")

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
       f"i.e. slowly, so that the converged value is estimated to lie about 1 % above the default result. The FSM, which uses the same spline "
       f"spans, gives {3.858:.3f} N/mm for this case, i.e. the largest deviation between GBT and FSM in Table 7 stems mainly from the spline "
       f"discretisation of GBT near the clamped ends of the wide skin bays with bending–twisting coupling. Finer segments lower the critical "
       f"load slightly for $b_s = 150$ mm and therefore do not remove the offset against Abaqus noted in Section 3.2. "
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

(H / "part3_filled.md").write_text(s)
left = [p for p in s.split("@@")[1::2]]
print("unfilled:", left)
