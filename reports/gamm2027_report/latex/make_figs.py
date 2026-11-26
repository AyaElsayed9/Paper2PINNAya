"""Figures of the Paper 2 progress report (same style as the GBT paper figures).
All values are read from results/report_gamm/report_data.json, except the history figure,
whose points are copied from docs/RESEARCH_LOG.md (entries 22 Oct - 23 Nov 2026)."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt

H = Path(__file__).parent
D = json.load(open(H.parents[2] / "results/report_gamm/report_data.json"))
MM = 1 / 25.4; FULL = 160 * MM
plt.rcParams.update({
    "text.usetex": True, "text.latex.preamble": r"\usepackage{mathpazo}\usepackage{amsmath}",
    "font.family": "serif", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.minor.size": 1.6, "ytick.minor.size": 1.6, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.0, "lines.markersize": 4,
    "legend.frameon": False, "legend.handlelength": 1.8, "pdf.fonttype": 42})
INK = "#1b2733"; MUTED = "#8a94a0"; GRID = "#e3e7ec"
NET = "#eb6834"; CHK = "#1baf7a"; GBT = "#2a78d6"; GP = "#8a94a0"
pct = mt.FuncFormatter(lambda v, _: (f"{v:g}" if v >= 0.1 else f"{v:g}") + r"\,\%")


def grid(ax):
    ax.grid(True, color=GRID, lw=0.4, zorder=0, which="major"); ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(H / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ---- timeline of the whole work (dates from the uploads, git history and docs/RESEARCH_LOG.md) ----
import datetime as dt
d0 = dt.date(2026, 9, 20)
dd = lambda m, d: (dt.date(2026, m, d) - d0).days
rows = [("Shared space, first network, GBT as teacher", [(dd(10, 22), dd(10, 25))], NET),
        ("Conditioning fixed, Paper 1 data term", [(dd(10, 26), dd(10, 27))], NET),
        ("Physics, data or hybrid; clean training", [(dd(11, 3), dd(11, 4))], NET),
        ("Wider design space, 1000 labels, data models", [(dd(11, 4), dd(11, 7))], NET),
        ("Lay-up experts and industrial network", [(dd(11, 7), dd(11, 9))], NET),
        ("Network against simple baselines", [(dd(11, 9), dd(11, 10)), (dd(11, 16), dd(11, 17))], NET),
        ("Physics-only experts (0 labels)", [(dd(11, 16), dd(11, 17))], NET),
        ("GBT check: speed, bounds, error estimate", [(dd(11, 10), dd(11, 10)), (dd(11, 17), dd(11, 17)), (dd(11, 23), dd(11, 23))], CHK),
        (r"Cost of a 5\,\% target", [(dd(11, 23), dd(11, 23))], CHK)]
fig, (t0, ax) = plt.subplots(2, 1, figsize=(FULL, 122 * MM), gridspec_kw=dict(height_ratios=[1.0, 1.2], hspace=0.38))
for i, (lab, segs, c) in enumerate(rows):
    y = len(rows) - 1 - i
    for s, e in segs:
        t0.barh(y, e - s + 1, left=s, height=0.55, color=c, alpha=0.75, lw=0)
    t0.text(segs[0][0] - 1, y, lab, ha="right", va="center", fontsize=6.6)
t0.set_xlim(dd(10, 22) - 26, dd(11, 25)); t0.set_ylim(-0.6, len(rows) - 0.1); t0.set_yticks([])
tk = [(10, 22), (10, 29), (11, 5), (11, 12), (11, 19)]
t0.set_xticks([dd(m, d) for m, d in tk]); t0.set_xticklabels([f"{d} {'Sep' if m == 9 else 'Oct' if m == 10 else 'Nov'}" for m, d in tk])
for sp in ("left", "right", "top"): t0.spines[sp].set_visible(False)
t0.tick_params(axis="x", top=False); t0.tick_params(axis="y", left=False, right=False)
t0.text(dd(10, 22) - 26, len(rows) + 0.1, "(a)", fontsize=8)

# ---- history of the error (from the research log) -------------------------------------------
day = lambda m, d: d - 21 if m == "Oct" else d + 10
inst = [("Oct", 22, 15, 140, "physics only, Adam\n(training panels)"), ("Oct", 25, 10, 16, "GBT as teacher"),
        ("Oct", 26, 3.65, None, "conditioning fixed"), ("Nov", 4, 1.42, None, "hybrid loss"),
        ("Nov", 6, 7.75, None, "all lay-ups"), ("Nov", 7, 4.7, None, "1000 labels"),
        ("Nov", 9, 0.24, 1.18, "8 experts, 120 labels"), ("Nov", 17, 4.68, 11.93, "8 experts, 0 labels")]
chk = [("Oct", 26, 0.22, None), ("Nov", 4, 0.20, None), ("Nov", 6, 0.21, None), ("Nov", 10, 0.11, 0.33),
       ("Nov", 17, 0.14, 0.48), ("Nov", 23, 0.45, None)]
ax.set_yscale("log"); ax.set_ylim(0.08, 250); ax.set_xlim(0, 34); grid(ax)
ax.axhline(5, color=INK, lw=0.6, ls=(0, (4, 3))); ax.text(0.4, 5.6, r"5\,\% target", fontsize=6.5, color=INK)
xs = [day(m, d) for m, d, v, h, _ in inst if h is None]; ys = [v for m, d, v, h, _ in inst if h is None]
ax.plot(xs, ys, color=NET, lw=0.6, alpha=0.4, zorder=2)
off = {"physics only, Adam\n(training panels)": (6, -2, "left"), "GBT as teacher": (6, 2, "left"),
       "conditioning fixed": (6, -7, "left"), "hybrid loss": (6, -7, "left"), "all lay-ups": (-6, 3, "right"),
       "1000 labels": (6, -7, "left"), "8 experts, 120 labels": (-6, 0, "right"), "8 experts, 0 labels": (6, 0, "left")}
for m, d, v, h, lab in inst:
    x = day(m, d)
    if h is None:
        ax.plot(x, v, "o", ms=4, color=NET, mec="white", mew=0.5, zorder=4); y = v
    else:
        ax.plot([x, x], [v, h], color=NET, lw=3, solid_capstyle="round", zorder=4); y = np.sqrt(v * h)
    dx, dy, ha = off[lab]
    ax.annotate(lab, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va="center", fontsize=6.5, color=INK)
for m, d, v, h in chk:
    x = day(m, d)
    if h is None:
        ax.plot(x, v, "s", ms=3.6, color=CHK, mec="white", mew=0.5, zorder=4)
    else:
        ax.plot([x, x], [v, h], color=CHK, lw=3, solid_capstyle="round", zorder=4)
ax.annotate(r"check, 120 labels", (day("Nov", 10), 0.19), xytext=(6, 0), textcoords="offset points", fontsize=6.5, color="#0f7a53", va="center")
ax.annotate(r"check, 0 labels", (day("Nov", 17), 0.26), xytext=(6, 0), textcoords="offset points", fontsize=6.5, color="#0f7a53", va="center")
ax.annotate(r"5\,\% target met, 0 labels", (day("Nov", 23), 0.45), xytext=(2, 9), textcoords="offset points", fontsize=6.5, color="#0f7a53", ha="right")
ax.set_xticks([1, 8, 15, 22, 29]); ax.set_xticklabels(["22 Oct", "29 Oct", "5 Nov", "12 Nov", "19 Nov"])
ax.yaxis.set_major_formatter(pct); ax.set_ylabel("mean error vs full GBT")
h1 = plt.Line2D([], [], color=NET, marker="o", ls="", ms=4, label="instant network answer")
h2 = plt.Line2D([], [], color=CHK, marker="s", ls="", ms=3.6, label="after the GBT check")
ax.legend(handles=[h1, h2], loc="upper right", ncol=2)
ax.text(-0.075, 1.04, "(b)", transform=ax.transAxes, fontsize=8)
save(fig, "fig_history")

# ---- experts and labels -----------------------------------------------------------------------
names = {"crossply": "Cross-ply", "quasiiso": "Quasi-iso.", "zerorich": r"$0^\circ$-rich", "pm45": r"$\pm45^\circ$"}
E = D["experts"]
fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 64 * MM), gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.32))
for i, e in enumerate(E):
    side, lay = e["expert"].split("_"); y = len(E) - 1 - i
    i0 = e["physics_L0"]["head_mean"]; i1 = e["hybrid_L120"]["head_mean"]
    c0 = e["physics_L0"]["s12_mean"]; c1 = e["hybrid_L120"]["s12_mean"]
    a.plot([i1, i0], [y, y], color=NET, lw=0.6, alpha=0.4)
    a.plot(i0, y, "o", mfc="white", mec=NET, mew=0.9, ms=4.2)
    a.plot(i1, y, "o", color=NET, ms=4.2)
    a.plot([c0, c1], [y, y], "s", color=CHK, ms=3.4)
a.set_yticks(range(len(E))); a.set_yticklabels([f"{names[e['expert'].split('_')[1]]}, {e['expert'].split('_')[0][:3]}." for e in E][::-1])
a.set_xscale("log"); a.set_xlim(0.07, 20); a.xaxis.set_major_formatter(pct); grid(a); a.tick_params(axis="y", length=0)
a.set_xlabel("mean error vs full GBT"); a.set_ylim(-0.7, len(E) - 0.3)
hs = [plt.Line2D([], [], marker="o", ls="", mfc="white", mec=NET, ms=4.2, label="instant, 0 labels"),
      plt.Line2D([], [], marker="o", ls="", color=NET, ms=4.2, label="instant, 120 labels"),
      plt.Line2D([], [], marker="s", ls="", color=CHK, ms=3.4, label="check, 12 steps")]
a.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, columnspacing=0.8, handletextpad=0.3)
a.text(-0.02, 1.13, "(a)", transform=a.transAxes, fontsize=8, ha="right")
L = D["labels_crossply_symmetry"]; ns = [0, 10, 20, 40, 120]
net = [L["physics_L0"]["head_mean"]] + [L[f"hybrid_L{n}"]["head_mean"] for n in ns[1:]]
ck = [L["physics_L0"]["s12_mean"]] + [L[f"hybrid_L{n}"]["s12_mean"] for n in ns[1:]]
gp = [L[f"gp_L{n}"]["head_mean"] for n in ns[1:]]
xp = np.arange(len(ns))
b.plot(xp, net, "-o", color=NET, ms=3.5, label="network, instant")
b.plot(xp[1:], gp, "--o", color=GP, ms=3.5, label="Gaussian process")
b.plot(xp, ck, "-s", color=CHK, ms=3.2, label="network start + 12-step check")
b.set_yscale("log"); b.set_ylim(0.1, 30); b.yaxis.set_major_formatter(pct); grid(b)
b.set_xticks(xp); b.set_xticklabels([str(n) for n in ns]); b.set_xlabel("number of GBT labels")
b.legend(loc="center left", bbox_to_anchor=(0.0, 0.42), ncol=1, handletextpad=0.4)
b.text(-0.02, 1.03, "(b)", transform=b.transAxes, fontsize=8, ha="right")
save(fig, "fig_experts")

# ---- convergence and error estimate ------------------------------------------------------------
CV = D["convergence"]
fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 62 * MM), gridspec_kw=dict(wspace=0.3))
it = np.arange(13)
for k, c, ls, lw, lab in [("random", GP, "--", 0.9, "random"), ("nn10", INK, ":", 0.9, "nearest of 10 labelled"),
                          ("nn120", INK, "-", 0.9, "nearest of 120 labelled"), ("net_L0", NET, ":", 1.1, "network, 0 labels"),
                          ("net_L120", NET, "-", 1.1, "network, 120 labels")]:
    a.plot(it, CV[k]["median"], ls=ls, color=c, lw=lw, label=lab)
a.axhline(5, color=INK, lw=0.6, ls=(0, (4, 3)))
a.set_yscale("log"); a.set_ylim(0.08, 100); a.set_xlim(0, 12); a.yaxis.set_major_formatter(pct); grid(a)
a.set_xlabel("iteration of the check"); a.set_ylabel("median error vs full GBT")
a.legend(loc="upper right", fontsize=6.3, handlelength=2.2)
a.text(-0.02, 1.03, "(a)", transform=a.transAxes, fontsize=8, ha="right")
P = np.array([[r[0], r[1]] for r in D["estimate_step5"]])
b.plot(P[:, 1], P[:, 0], "o", ms=1.6, color=CHK, alpha=0.45, mew=0)
z = np.array([0.05, 50]); b.plot(z, z, color=INK, lw=0.7); b.fill_between(z, z, 0.01, color="#c0392b", alpha=0.06, lw=0)
b.set_xscale("log"); b.set_yscale("log"); b.set_xlim(0.05, 20); b.set_ylim(0.5, 50)
b.xaxis.set_major_formatter(pct); b.yaxis.set_major_formatter(pct); grid(b)
b.set_xlabel("true error after 5 iterations"); b.set_ylabel(r"error estimate $\eta$")
b.text(5.5, 0.7, r"$\eta$ below the true error:" + "\n" + r"0 of 640 checks", fontsize=6.5, ha="center")
b.text(-0.02, 1.03, "(b)", transform=b.transAxes, fontsize=8, ha="right")
save(fig, "fig_check")
print("figures done; estimate/true min ratio", (P[:, 0] / P[:, 1]).min(), "n", len(P))
