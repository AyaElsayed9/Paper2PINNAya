# Notes for Claude (read first in every session)

1. Read `docs/STATUS.md` (where we are), `docs/RESEARCH_LOG.md` (every try, result and decision) and
   `docs/PAPER1_SUMMARY.md` (Paper 1 and Aya's plan for Paper 2) and `docs/LITERATURE_AND_PAPER1_IDEAS.md`
   before proposing anything for Paper 2. Add every new try, result and decision to `docs/RESEARCH_LOG.md`
   and update `docs/STATUS.md` in the same commit.
2. Working rules from Aya:
   * Always be 100 % honest, including about failures and limits.
   * Never modify Abaqus numbers; no hardcoding or tuning to Abaqus.
   * Simple language; Arabic (Egyptian) when she asks.
   * Paper 2 uses the GBT code, not Abaqus. The network must be fast for any design and stay a PINN.
   * Solver `solver/aya_gbt_v32.py` is frozen (bug fixes only, listed in `solver/CHANGES.md`).
3. Git: author Aya, committer Claude:
   `GIT_COMMITTER_NAME=Claude GIT_COMMITTER_EMAIL=noreply@anthropic.com git -c user.name="Aya" -c user.email="ayamedhatt93@gmail.com" commit`,
   push to origin main.
4. `manuscript/` holds the finished periodic-GBT paper (Journal of Composites Science); Dr. Navid Kharghani
   is co-author (cite his GBT–Ritz works, Silvestre–Camotim notation).
