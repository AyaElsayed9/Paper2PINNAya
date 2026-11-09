# Literature research plan for Paper 2 (9 Nov 2026, Aya: "plan these researches")

Done before the introduction is written; every source read (not only the abstract), notes + reference in this
file, and the numbers go into one comparison table. Existing notes: `docs/LITERATURE_AND_PAPER1_IDEAS.md`.

| # | Topic | Question it answers for Paper 2 | What to collect |
|---|---|---|---|
| 1 | Machine learning / PINN for buckling of stiffened composite panels and thin-walled members | Where our accuracy, inputs and speed stand | inputs (number, type), solver for labels (FE / GBT / analytic), number of labels, error (mean / max), speed, answer check yes / no |
| 2 | Checking a fast answer for eigenvalue problems (Rayleigh–Ritz upper bound, Kato–Temple / residual bounds, certified reduced-basis methods) | What our physics check is, what is new and what is known | type of bound, cost of the check, whether it needs the full model |
| 3 | Physics-informed vs plain data surrogates with few labels (incl. Gaussian processes for buckling) | Why a PINN if a GP is as good with 120 labels (our 9 Nov baseline) | label counts where physics helps, reported comparisons |
| 4 | Solver-in-the-loop surrogates: multi-fidelity, active learning, trust / fallback rules | The GBT-call study (network → ensemble trigger → physics check → full GBT), as in Paper 1 | share of solver calls, time saved, miss rate |
| 5 | Neural eigenvalue solvers (Ritz / subspace networks, energy natural gradient, PINNs for eigenproblems) | Method section: what our training builds on | loss, optimizer, problem size |
| 6 | GBT with machine learning; GBT–Ritz works (Kharghani; Silvestre–Camotim notation) | Whether anyone has combined GBT and ML; citations for the GBT part | method, cases, validation |
