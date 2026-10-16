# MDBuddy development notes

- `src/mdbuddy/`: the pipeline. Stages in `pipeline.py` (`STAGES`), one module per stage:
  `fetch`, `prepare`, `forcefields`, `martinize`, `membrane` (COBY), `assemble` (GROMACS folder, validation), `mdp`.
- `skills/mdbuddy/`: the agent skill (SKILL.md + references). Keep it in sync with CLI flags and spec keys.
- The spec schema lives in `spec.py` (`DEFAULTS`, `_validate`). New keys need a default,
  validation, and an entry in `skills/mdbuddy/references/spec_reference.md`.
- Pinned force-field files: `forcefields.REGISTRY`. Changing a URL means changing the commit
  hash *and* the SHA-256 together.
- Core package depends only on PyYAML. vermouth/COBY are run as external programs
  (`martinize2`, `python coby_build.py`), never imported, so builds can be replayed by hand.
- Tests: `pytest`. `tests/conftest.py` provides fake martinize2/COBY/force-field files; keep
  the real tools out of the unit tests.
