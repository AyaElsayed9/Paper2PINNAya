# Manuscript: GBT with exact periodic boundary conditions (omega-stiffened panels)

- `part1.md` … `part5.md`: manuscript text (Markdown + LaTeX maths; `{{key}}` citations, `\eqn{label}` equation labels).
- `refs.py`: reference list (checked against Crossref and the MDPI template).
- `build.py`: numbers sections, equations and citations, builds the MDPI Word file on `src/template.docx` and a LaTeX preview PDF into `out/`.
- `figs.py`: Figures 1–5 from the solver; `stats.py`: all error/shift statistics from `src/omega_results_V32.csv`.

Build: `python manuscript/build.py` (needs pandoc and xelatex).
