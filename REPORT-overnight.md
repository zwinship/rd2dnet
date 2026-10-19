# Overnight report (2026-10-19)

Work is on two PRs: #4 (Stata first run, merged) and #5 (`parity-harness`: parity harness
plus the fixes below). Decisions waiting for you are in `DECISIONS-overnight.md`.

## Status by language

| | status |
|---|---|
| **Stata** | First run in Stata 19 SE (AppStream). `rd2dnet_tests.do`: ALL TESTS PASSED at 253b969; `install_check.do` (net install + every help example from an empty folder): INSTALL CHECK PASSED; help files render; no leftover objects beyond the documented ones. **The fixes made later in PR #5 (point nodes, 1/h² side weights, leverage clip) are written in Mata but not yet run in Stata**: the AppStream session expired and needs your sign-in. Run `rd2dnet_tests.do` (new block `t_pointnodes`) and `parity/run_stata.do` first thing. |
| **Python** | Full test suite green locally and in CI (3.10-3.13, Windows, notebook). |
| **R** | Not buildable locally (no Rtools). R CMD check --as-cran green on five platforms in CI, including the new regression tests. |
| **Parity** | New `parity/` harness: 40 random designs (grids, rivers, rings, directed, transition costs, travel-time costs, fuzzy, covariates, clusters, Conley, every bwselect and kernel). Python vs R in CI (`parity` job). Before the fixes, 6 of 40 cases disagreed; see below. |

## Bugs found and fixed (most severe first)

1. **Directed networks: wrong `s`, `d`, `o` for locations on one-way edges** (Python, R, Stata). A
   location partway along an edge was reached through either end node, as if the edge
   were two-way. Found by the invariance tester; hand-computed example in the tests. PR #5.
2. **Transition costs: penalty for the location's own edge class never charged** (all three).
   Found independently by the parity harness (5 of 40 cases) and the invariance tester.
   Fix for 1 and 2: with directed networks or `TransitionCost`, every location and
   evaluation point becomes a node of the network (its edge is split there). PR #5.
3. **Splitting an already-split network reused node labels and moved existing nodes** (Python, R).
   Reusing `coords.network` (e.g. for a placebo boundary) changed `s` by up to 1.5 km.
   Found by use-case tester B after fix 1. PR #5.
4. **Covariates/pooled with h0 != h1: side weights lacked rd2d's 1/h_side² factor** (all three).
   Point estimates differed from rd2d by ~0.4%. Found by the fidelity audit. PR #5.
5. **`node_contrast()` across boundary-split edges counted the contrast against 0** (Python, R):
   new crossing nodes had no node attributes. PR #5.
6. **HC2/HC3 with leverage 1**: clip of `1 - h_ii` at 1e-12 turned rounding into SE noise
   (~1e-3, different in R and Python); now 1e-8 as in rd2d (all three). Found by parity. PR #5.
7. **Stata `net install` did not install the mock-city datasets** (`F` lines in the pkg). PR #4.
8. **Stata `rdnetdensity` returned nothing in `r()`.** PR #4.

## Found, not yet fixed (ranked)

Wrong numbers or misleading output:
- **Density test over-rejects for round or wavy zone outlines** with the default
  `normalize="network"`: rejection 0.82-0.88 at 5% for circles, 0.44 for a wavy school zone
  (0.06-0.08 with `normalize="none"`). Selected bandwidths exceed the zone's street supply.
  (use-case tester A, `C:\Users\zwins\agent_usecase_a\density_mc2.py`). Highest priority.
- **Near-singular q-fits reported as "ok"**: exploded `tau_bc`/`se_rb` (tutorial placebo: tau_bc
  92.7 with tau -0.02), or collapsed HC1 SEs with coverage 0.37 on a one-way grid. Also a
  singular q-fit silently switches the p-fit to the `s` basis and biases `tau`.
  (invariance tester #3; use-case A #4; use-case B 1a). Needs a design decision, see DECISIONS.
- **Panels: bandwidth floor counts observations, not locations** (repeated readings): half the
  points "singular", others with h of 16-150 m. (use-case A #2)
- **Lon/lat passed as metres** goes through silently (snap depth ~4,700 km). `lonlat=True`
  does not project `eval_points` coordinates; `rd2dnet(..., lonlat=True)` gives a spurious
  side-mismatch warning. `from_lines(DataFrame, crs=, to_crs=)` ignores the CRS arguments.
- **Off-network distance is free with travel-time metrics** (`snap_scale` 0): points 1.5 km
  from their street get `s = 0`. Needs a default decision.
- Ring boundaries: evaluation points can collapse onto one node without a warning (and get
  double weight in the default WBATE); `Boundary.grid()` on a closed ring returns both 0 and L.
- Treatment labels `t` that disagree with geography are used silently.

Crashes / API:
- `res.plot(true_tau=<scalar>)` crashes; `Boundary.from_region(GeoSeries)` gives a pandas error;
  `rdbwnet` docstring column names are out of date; README quick start uses undefined `take_up`.
- Feet layers with mph speeds give travel times 5.2x off without a warning.
- `from_place` falls back to the mock city after a 111 s timeout (risky default).

Fidelity to rd2d/rddensity (audit, `C:\Users\zwins\agent_fidelity`): everything else matches
rd2d to machine precision (DPI constants, all bwselect rules, CER factor, HC0-HC3, joint
covariance, WBATE/LBATE, kernels, rddensity estimator and jackknife). Intentional or
documented differences that need your call are listed in DECISIONS (density bandwidth
selector, min_obs default, cluster CR1 factor, uniform-band seed, mass points).

## Earlier open items from the brief (not touched yet)
Flow mode for rivers, precomputed-metric splitting and OD constructor, boundary points on
the line, J=1 critical value, eval_points (500, 510) ambiguity, lbate min/abs, pandas index,
naming consistency.
