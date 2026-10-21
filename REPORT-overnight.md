# Overnight report (nights of 2026-10-19 and 2026-10-20)

## Night 2 summary (start here)

- **Stata verified on everything up to and including night 1.** Run in Stata 19 at the start
  of night 2: ALL TESTS PASSED, INSTALL CHECK PASSED, and PARITY PASSED (40 random designs,
  worst difference to Python 5.8e-9). The only stale references were the transition-cost
  values, which predate the own-edge penalty fix (PR #10).
- **The open questions are decided, tested and merged** (PRs #11, #12, #15, #18). See
  `DECISIONS-overnight.md` items 19-40.
  - Sorting tests: straight-line density by default, plus the Bugni-Crippa-Restrepo
    manipulation test.
  - An identifiable local basis instead of the `s`-only fallback.
  - HC3 as the default.
  - Off-network distance charged at the street's cost rate.
  - `direction="flow"` for rivers.
  - WBATE weights per place; `lbate` min/absmax; seed 0; `from_place` fails loudly.
- **CI:** all green on the merged result: pytest on 3.10-3.13 and Windows, the notebook,
  R CMD check on 5 platforms, Python-R parity, and live OSM. R's live-OSM job had kept
  failing on Overpass 502/504 errors; R now retries on mirrors.
- **Stata still to run.** The AppStream session expired during a pause at about 03:20 and
  needs your sign-in. Stata code merged after that is written, and three static reviews
  checked it against Python (all findings fixed), but it has not been run in Stata. That
  covers:
  - snap cost; place weights, lbate kinds and seed;
  - `rdnetmanip` and normalize(euclidean); direction(flow);
  - everything from the bug hunt below, including CR2.

  The commands below run all of it, including parity against regenerated cases.
- **Bug hunt (PRs #20, #21, #22):** three rounds of review and testing found 13 sources of
  wrong numbers, all fixed in all three languages (details below). The biggest change:
  **cluster-robust errors now default to CR2 with Bell-McCaffrey degrees of freedom**.
  rd2d's CR1 covered only 0.66-0.91 in Monte Carlo runs; CR2 covers 0.955-0.97. CR1
  remains as an option. See DECISIONS item 38; it departs from rd2d on purpose, so please
  check you agree.

## Status by language

| | status |
|---|---|
| **Stata** | Run in Stata 19 at the start of night 2: ALL TESTS PASSED, INSTALL CHECK PASSED, PARITY PASSED. **Everything merged after about 03:20 on 10/21 is unrun in Stata** (sign-in needed); new test blocks `t_round2` and `t_round3`. |
| **Python** | Full test suite green, locally and in CI (3.10-3.13, Windows, notebook). |
| **R** | Can't be built locally (no Rtools). In CI, R CMD check --as-cran is green on five platforms, with regression tests for every fix; CR2 matches Python's reference values to 1e-9. |
| **Parity** | New `parity/` harness: 40 random designs covering grids, rivers, rings, directed networks, transition and travel-time costs, fuzzy designs, covariates, clusters, Conley, every bwselect and every kernel. Python vs R runs in CI (`parity` job). Before the fixes, 6 of 40 cases disagreed. |

## First thing this morning: run Stata on the merged code

All PRs (#4-#22) are merged into `main` (5983390). In Stata 19 (AppStream), run these, typing
each line in the Command window:

```stata
clear all
do https://raw.githubusercontent.com/zwinship/rd2dnet/598339034dd59f96299fbef058e26cce7b7d1adb/stata/test/fetch.do 598339034dd59f96299fbef058e26cce7b7d1adb
cd test
do rd2dnet_tests.do
```

The test suite must end with `ALL TESTS PASSED`. Failures are listed at the end. Then
run the install check and the Stata side of the parity harness:

```stata
do https://raw.githubusercontent.com/zwinship/rd2dnet/598339034dd59f96299fbef058e26cce7b7d1adb/stata/test/install_check.do 598339034dd59f96299fbef058e26cce7b7d1adb
do https://raw.githubusercontent.com/zwinship/rd2dnet/0f82a4a7566675d74fb97ff95fde068d135e14c7/parity/fetch_stata.do 0f82a4a7566675d74fb97ff95fde068d135e14c7
do parity/run_stata.do
```

`run_stata.do` uses the copy that `install_check.do` just installed. It ends with
`PARITY PASSED` or a list of differences.

The `parity-cases` branch is a scratch branch holding the generated cases for this
command. It is not for merging; delete it whenever you like.

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

## Night 2 bug hunt (PRs #20, #21, #22)

Three rounds of testing ran on the merged code:
- a statistical code review;
- use-case and adversarial testers;
- a brute-force checker (its own Dijkstra on random networks);
- a Monte Carlo study of cluster-robust errors;
- static reviews of every Mata change.

Every fix is in all three languages unless noted, with regression tests in Python and R and
a Stata test block (`t_round2`, `t_round3`).

**Wrong numbers, fixed (most severe first):**

1. **Cluster-robust standard errors were far too small.** rd2d's CR1 covered 0.66-0.91 in
   Monte Carlo runs with 6 to 100 clusters per side. The default is now CR2 (Bell-McCaffrey)
   with Satterthwaite degrees of freedom and t critical values (new `df` column), which
   covers 0.955-0.97. `cluster_vce="cr1"` (Stata: `cr1`) gives rd2d's behaviour.
   DECISIONS item 38 has the mechanism and the table.
2. **Boundary splits.** Two errors, both fixed by splitting every edge that meets the
   boundary at each crossing:
   - Edges whose stored geometry runs v -> u were split at the mirror-image point. This
     covers every `from_lines` layer with lines drawn "backwards" and OSM graphs made
     undirected. In the mock city, s was off by up to 39 m for 1,870 of 3,000 points.
   - An edge crossing the boundary twice with both ends on one side (a notch) was not split
     at all.
3. **A road with both ends on the boundary** (a loop into one side) counted as running along
   it, so its locations got s = 0.
4. **Closed rings that follow streets:** evaluation points were moved to the street at the
   ring's start (e.g. requested arc 250 -> node arc 383).
5. **Transition penalties with `direction="in"`/`"flow"`** were charged in reverse order.
6. **Pooled fits with h0 != h1** scaled the common slopes differently on each side, so a
   noiseless plane was not fitted exactly.
7. **Covariates in large units** (dollars) made every point "singular".
8. **Directed and transition-cost networks:** point nodes sat on the chord between the end
   nodes, not on the street (Python, R).
9. **Callable metrics** on attributes such as minutes copied the full cost to each piece of
   a split edge. They now warn, pointing to `extensive=`. Separately, the split cache ignored
   `extensive`.
10. **`density(eval_index, window)`** did not wrap around rings.
11. **Nodes on the boundary line's extension** changed side when `b` was reversed, so
    `rd2dnet()` estimates depended on the order of `b`.
12. **Conley covariances** could be non-PSD on networks. They are now projected to PSD,
    with a warning.
13. **`rdbwnet(vce="cluster"/"conley")`** did not reproduce `rdnet`'s bandwidths
    (Python, R).

**New statuses and warnings** (thresholds from Monte Carlo runs):
- `"weak first stage"` when |fs| < 2 fs_se. These points are left out of the default
  WBATE and of LBATE.
- HC errors with repeated locations (panels).
- Sorting tests on points on or near streets: a ratio of median snap distances, which
  catches GPS noise up to about 15 m. Before, only points exactly on streets were caught.
- `from_lines` layers whose ends mostly do not join.
- `travel_time(default_kph)` covering most edges.
- `vce="conley"` with a cutoff above the bandwidth.

**Crashes and bad inputs fixed:**
- Tuple node labels with directed networks (Python).
- `TransitionCost(travel_time())` (Python).
- `wbate()` returning NaN when any point had no estimate.
- Non-integer `p`/`q`/`deriv`/`repp`/`k_max`, `repp=0`, `k_max=0`, malformed `tangvec`.
- NaN coordinates in `rd2dnet()`.
- `lat`/`lon` columns given in that order.

**Docs:** the tutorial now says HC3 and `"euclidean"` and shows the manipulation test. Also
fixed: the `network_coordinates` direction notes, `NetworkGraph(crs=)`, the multi-q
printout, and the common-bandwidth clipping.

**Also:** R's OpenStreetMap download now retries on Overpass mirrors. The live-OSM CI job
had failed on 502/504 errors from overpass-api.de.

**Checked and correct** (no change needed):
- s, d, o, pairwise distances and the HAC against brute-force Dijkstra for undirected,
  out/in/flow, travel time with snap cost, and transition costs with connectors.
- Invariance to translation, rotation, reflection, relabelling, row order, y scale and
  shift, weights, duplicated rows, and splitting an edge into collinear pieces.
- 192 vce/cluster/Conley/fuzzy/covs/pooled/deriv combinations.
- Hand-computed HC0-3 and CR1.
- The k-means DP and the LCG.
- Identical results with 1 and 8 threads.

## Known limitations (not fixed)

- **Parallel edges between the same two nodes** (dual carriageways, OSM `MultiGraph` keys):
  only the cheapest geometry is kept for snapping. A location on the other twin snaps to
  it, about 50 m off in the test case.
- **Exact ties on perfectly regular grids:** a location exactly equidistant from two
  crossings can switch crossings with 1e-12 noise. Mass points are detected by exact
  coordinate equality.
- **Parameter choices:** `manipulation_test` with `k_max=q` is O(q^3) (24 s at q = 1000).
  `h=[a, b]` means one value per point when J = 2.
- **Clustering with very few clusters:** CR2 intervals are honest but wide. With 6 blocks
  per side the degrees of freedom are about 2.
- **Untested in Stata:** everything after 03:20 on 10/21. The Mata code had three static
  reviews against Python; all findings are fixed.

## Night 1 open items, now resolved

All the night-1 "found, not yet fixed" items are settled. See DECISIONS items 19-40:
- the density test on curved boundaries;
- near-singular q-fits;
- few clusters;
- forgotten `lonlat`;
- free off-network distance;
- WBATE double weights;
- the `from_place` fallback.

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
