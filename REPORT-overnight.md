# Overnight report (2026-10-19)

## Status by language
- **Stata:** run for the first time (Stata 19 SE, AppStream). `rd2dnet_tests.do`: ALL TESTS
  PASSED. New `install_check.do` (net install at a commit + every help example from an
  empty folder): INSTALL CHECK PASSED. Help files render; memory left clean (documented
  frames/matrices, Mata externals `rdn_C`, `rdn_DIM`). Merged in PR #4.
- **Python:** test suite green in CI (3.10-3.13, Windows, notebook).
- **R:** not built locally (no Rtools; see DECISIONS 4). New parity CI job runs R vs Python.

## Bugs fixed (ranked)
1. Stata `net install` did not install the mock-city .dta files, so `rdnetmock` failed after
   installation (pkg `f` -> `F`). PR #4.
2. Stata `rdnetdensity` returned nothing in `r()` (missing `return add`). PR #4.
3. Stale test expectation (eval 1 now fitted in s alone). PR #4.

## Open
- Parity harness (`parity/`, `.github/workflows/parity.yml`) on branch `parity-harness`:
  needs its first CI run; Stata side (`parity/run_stata.do`) not yet run.
- Background reviews (fidelity to rd2d/rddensity, invariance, Monte Carlo, use cases,
  robustness/scale) were started but their results were not processed before the
  session's usage limit; rerun them.
- All open issues listed in the brief (directed "flow" mode, precomputed metric,
  snap_scale default, boundary points, ring density MC, panels, J=1 critical value,
  eval_points ambiguity, lbate min, pandas index, naming) are untouched.

## Decisions waiting for you
See DECISIONS-overnight.md (items 4 and 5 especially).
