# Manuscript: GBT with exact periodic boundary conditions (omega-stiffened panels)

Structure follows the GBT–Ritz template (Kharghani & Mittelstedt): 1 Introduction, 2 Theory and formulations
(2.1 cross-section analysis, 2.2 periodic boundary conditions, 2.3 member analysis, 2.4 mode decomposition),
3 Results and discussion, 4 Conclusions. Periodic side edges only.

- `part1.md` front matter + introduction, `part2.md` theory, `part3.md` results (placeholders `@@..@@` filled
  from the result files by `fill_results.py` -> `part3_filled.md`), `part4.md` conclusions + back matter.
- `refs.py`: references (checked against Crossref / template). `build.py`: MDPI docx on `src/template.docx` + LaTeX preview.
- Data: `results/periodic_curves.jsonl` (160 GBT runs: N_cr and energy participation over beta),
  `results/periodic_conv.json` (discretisation and link-wall checks), `results/strip_timing.json` (timing GBT vs FSM),
  `src/omega_results_V32.csv` (Abaqus / FSM reference values, read-only).
- Checks: `verify_math.py` (explicit modal matrices vs code, classical limit vs signature curve, periodicity),
  `verify_interaction.py` (conventional-shear interaction matrices).
- Figures: `figs.py` (Figs. 1-2), `figs_periodic.py` (Figs. 3-5).

Build: `python fill_results.py && python build.py` (needs pandoc and xelatex).
