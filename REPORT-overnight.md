# Overnight report (2026-10-19)

Work is in two PRs: #4 (Stata's first run, merged) and #5 (`parity-harness`: the parity
harness plus all fixes below). Decisions waiting for you are in `DECISIONS-overnight.md`.

## Status by language

| | status |
|---|---|
| **Stata** | First run in Stata 19 SE (AppStream). `rd2dnet_tests.do`: ALL TESTS PASSED at 253b969. `install_check.do` (net install, then every help example from an empty folder): INSTALL CHECK PASSED. Help files render, and nothing is left behind except the documented objects. **The Mata changes in PR #5 have not been run in Stata.** The AppStream session expired and now needs your sign-in. Please run `rd2dnet_tests.do` (new block `t_pointnodes`) and `parity/run_stata.do` first thing. |
| **Python** | Full test suite green, locally and in CI (3.10-3.13, Windows, notebook). |
| **R** | Can't be built locally (no Rtools). In CI, R CMD check --as-cran is green on five platforms, including the new regression tests. |
| **Parity** | New `parity/` harness: 40 random designs covering grids, rivers, rings, directed networks, transition and travel-time costs, fuzzy designs, covariates, clusters, Conley, every bwselect and every kernel. Python vs R runs in CI (`parity` job). Before the fixes, 6 of 40 cases disagreed. |

## Bugs found and fixed (most severe first)

1. **Directed networks: `s`, `d` and `o` were wrong for locations on one-way edges** (Python, R, Stata).
   A location partway along an edge was reached through either end node, as if the edge
   were two-way.
2. **Transition costs: the penalty for the location's own edge class was never charged**
   (all three). The parity harness (5 of 40 cases) and the invariance tester found this
   independently. The fix for 1 and 2: with directed networks or `TransitionCost`, every
   location and evaluation point becomes a node, and its edge is split there.
3. **Splitting an already-split network reused node labels and moved existing nodes**
   (Python, R). Reusing `coords.network`, for example for a placebo boundary, changed `s`
   by up to 1.5 km.
4. **Covariates or `pooled` with h0 != h1 lacked rd2d's 1/h_side² side weights** (all
   three). The tau difference from rd2d was about 0.4%.
5. **Panels: the bandwidth floor and `min_obs` counted rows, not locations** (all three).
   - With 30 readings per site the floor fell to the second-nearest site, and half the
     points came back "singular".
   - Points with only 7 sites per side were fitted at order 2 (RMSE 69 in the Monte Carlo).
   - Both now count distinct locations, as rd2d and rdrobust do.
6. **`node_contrast()` on edges split at the boundary** (Python, R) measured the contrast
   against 0: the new nodes had no attributes. They are now interpolated.
7. **Bandwidth selector: rounding decided the basis on ill-conditioned pilots** (all three). At
   Gram condition numbers around 1e10 (below the old 1e12 cut-off) plug-in variances came out
   as ±7000 depending on rounding, so Linux and Windows Python disagreed (0.3% on an `imserd`
   bandwidth) and a 1e-12 perturbation of `y` switched the basis. Cut-off now 1e8 for the
   selector's h-scaled bases; no well-posed fit on the mock city exceeds 1e7, so reference
   results are unchanged. Found by the parity harness.
8. **HC2/HC3 with leverage 1**: 1 - h_ii was clipped at 1e-12, which turned rounding into SE
   noise of about 1e-3 (different in R and Python). It is now clipped at 1e-8, as in rd2d
   (all three).
9. **Silent misuse is now caught** (Python, R):
   - Pandas inputs with a non-default index now warn. A shuffled Series gave tau_bc -0.08
     instead of 0.93.
   - 2-D `y`/`weights` are now errors.
   - `density()` rejects negative or NaN weights (it used to drop them silently).
   - Zero-length boundaries, NaN `eval_points`, `n_eval < 1` and negative or NaN
     `snap_scale` are now errors.
10. **J = 1 uniform critical value**: it was simulated (1.88-2.04 depending on the seed);
   it is now the normal quantile (all three).
11. **lon/lat input**:
    - `network_coordinates(lonlat=True)` now projects `eval_points` given as coordinates.
    - `from_lines(table, crs=4326)` now projects to UTM instead of keeping degrees.
12. **Inferring the treated side** (all three): observations beyond the ends of the line
    through `b` decided the side and triggered spurious warnings. They are now ignored.
13. **Treatment labels that contradict geography** now trigger a warning (all three).
14. **Stata: `net install` did not install the mock-city datasets.** PR #4.
15. **Stata: `rdnetdensity` returned nothing in `r()`.** PR #4.
16. **Smaller fixes**:
    - `Boundary.grid()` on a closed ring repeated the start point.
    - `plot(true_tau=<scalar>)` crashed.
    - README examples used undefined names.

## Found, not yet fixed (ranked)

**Statistical:**

- **The density test is badly mis-sized on curved boundaries.** On a round zone with
  `normalize="network"` it rejects 98-100% of the time at 5%; on a wavy school zone, 44%.
  `"none"` gives 10-20%. Street-length normalisation is only valid when observations sit
  on the streets. The README now warns about this; the method needs your decision.
- **Near-singular q-fits are reported as "ok".**
  - The tutorial's placebo point gives tau_bc 92.7 while tau is -0.02.
  - HC1 SEs collapse on a one-way grid (coverage 0.37).
  - A singular q-fit also silently switches the p-fit to the `s` basis.
- **Clustering by spatial block with 5-6 clusters per side fails**: 0.70-0.74 pointwise,
  about 0.2 for the uniform band. The README now says so. Consider a status or warning
  when few clusters fall within the bandwidth.

**Usability:**

- Forgetting `lonlat=True` goes through silently (snap depth about 4,700 km).
- Off-network distance is free with travel-time metrics (`snap_scale` 0).
- Evaluation points that snap to the same place get double weight in the default WBATE.
  The warning only fires below 75% distinct places.
- Feet layers combined with mph speeds give travel times 5.2x off, with no warning.
- `from_place` quietly falls back to the mock city after a 111 s timeout.
- Several malformed inputs raise raw numpy or GEOS errors, and collinear covariates
  give "singular" with no hint (`C:\Users\zwins\agent_robust\repro_all.py`).
- `Boundary.from_region(GeoSeries)` raises a pandas error.
- The column names in the `rdbwnet` docstring are out of date.

## Monte Carlo summary (`C:\Users\zwins\agent_mc`, 1000 replications per cell)

| design | HC1 pointwise | HC3 pointwise | uniform band (HC3) | WBATE |
|---|---|---|---|---|
| grid, straight boundary (6k) | 0.938 | 0.951 | 0.951 | 0.94 |
| irregular network (6k) | 0.941 | 0.951 | 0.930 | 0.94 |
| ring from `from_region` (6k) | 0.925 | 0.945 | 0.906 | 0.95 |
| river, s-only fits (10k) | 0.939 | 0.951 | 0.941 | 0.95 |

- With spatial shocks: Conley 0.96-0.98, HC1 0.90-0.92, cluster by block 0.70-0.74.
- Panels clustered by unit: 0.88 with about 48 sites per side, 0.94 with about 190.
- Density test size on a ring: 0.98 (`network`), 0.10 (`none`).
- The full table is in the Monte Carlo agent's summary files there.

## Fidelity to rd2d and rddensity (`C:\Users\zwins\agent_fidelity`)

Everything else matches rd2d to machine precision: the DPI constants, every bwselect
rule, the CER factor, HC0-HC3, the joint covariance, WBATE and LBATE, the kernels, and
rddensity's estimator and jackknife. `rdbwnet` scaleregul=1 and `rdnet` scaleregul=3 are
faithful to rd2d. The remaining differences are intentional; DECISIONS lists them.

## Items from the brief not yet addressed

- a flow mode for rivers;
- splitting with a precomputed metric, and an OD constructor;
- boundary points on the line;
- the `eval_points=(500, 510)` ambiguity;
- min/abs options for `lbate`;
- naming consistency.
