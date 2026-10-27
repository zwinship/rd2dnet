# Overnight decisions (2026-10-19)

Each entry: the decision, the alternatives, and why. All are reversible.

1. **Working branch `stata-testing`** for the Stata fixes, merged through a PR once CI is green.

2. **Test expectation for evaluation point 1 changed (Stata test suite).** The suite expected
   eval 1 of the mock city to have status 4 ("no bandwidth") with automatic bandwidths.
   Since 0.2.0 that point is fitted in `s` alone (basis 2) and estimated, in Python too
   (`basis = "s"`, status ok). The check now asserts status 0 and basis 2.
   Alternative: keep the old expectation and change the code; rejected because the
   0.2.0 behaviour is intended and matches Python and R.

3. **Help example for segment networks made runnable.** `help rdnetset` used
   `import delimited roads.csv`, a file nobody has. The example now builds the segment
   frame from the mock city's streets with `frlink`/`frget` and mentions the CSV route
   in the text. Alternative: ship a `roads.csv`; rejected to keep the package small.

4. **R tooling.** Rtools is not installed on this machine, so the R package (C++ code) cannot
   be built locally. I tried a silent per-user install of Rtools 4.5 from CRAN; the
   installer asked for administrator rights and is waiting behind a Windows elevation
   prompt (processes `rtools45.exe`/`rtools45.tmp`, which I cannot stop without
   elevation). Cancel that prompt in the morning, or install Rtools if you want local R
   checks. R is verified through CI (R CMD check --as-cran on five platforms) and,
   for cross-language numbers, through a GitHub Actions job that runs R on the same data.

5. **s-only fallback triggered by pilot fits (needs your decision).** On a jittered 100 m
   grid (parity case 05), points are fitted in `s` alone although the order-p (s, o) fit
   is identified: only 2-3 streets cross the boundary within the bandwidth, so the
   cubic bias-pilot fit in `o` is rank-deficient, `rdbwnet` marks (s, o) "singular" and
   falls back. Options: keep (documented behaviour), or reduce only the pilot's degree
   in `o`. Left unchanged.

6. **Parity harness on its own branch (`parity-harness`)**; generated cases are not
   committed (`parity/.gitignore`); CI regenerates them.

7. **Exact distances on directed networks and with transition costs: split edges at the
   locations.** Alternatives: (a) state-level "arc entry" distances in the engine (no graph
   growth, but a rewrite of every distance routine in three languages, including the Conley
   kernel); (b) split only for these metrics (chosen). Undirected plain-cost networks, the
   common case, are untouched. Cost: the network grows by one node per distinct location for
   directed/transition metrics. `coords.network` therefore carries `rdnet_point_*` nodes.

8. **Covariate side weights 1/h_side² (rd2d's choice)** applied in all languages, also in the
   Conley auxiliary fit and for `pooled`. Changes estimates only when h0 != h1.

9. **Leverage clip 1e-8 (rd2d) instead of 1e-12.**

10. **Parity tolerance 1e-7** instead of 1e-8: automatic bandwidths differ by up to 3.5e-8
    between Python and R from rounding in near-ill-conditioned fits.

11. **AppStream sign-in.** The Stata session expired around 06:30 and the link now asks for an
    AWS sign-in; per the rules I did not sign in. Stata changes after commit 253b969 are
    untested in Stata. PR #5 is merged only if CI is green, on the basis that main is
    unpublished and the Stata port has never been released; please run the Stata suite
    before any release.

## Waiting for your decision (from the audits)

- **Near-singular q-fits:** (a) flag points whose order-q fit has leverage 1 or a condition
  number above a threshold (new status "unstable"); (b) when only the q-fit is singular in
  (s, o), keep the p-fit in (s, o) and compute the bias correction in s alone, instead of
  refitting everything in s alone. I recommend both; not done because it changes results.
- **Density bandwidth selector** is not rddensity's (normal-reference pilots, `comb` default).
  rd2dnet uses rd2d-style constants and `common = min(hl, hr)`. Together with the
  over-rejection on round zones, consider porting rddensity's selector.
- **min_obs/bwcheck default:** rd2d's `rd2d()` uses 50 + p + 1; rd2dnet uses 20 + #terms(q).
- **Cluster SEs:** rdnet uses rd2d's `vce="hc0"` + CR1; rd2d's default (`hc1` + cluster) is
  about 1.7% larger. Also clustering under-covers in several Monte Carlo runs; Conley does well.
- **Uniform band:** rd2d uses repp=1000, a fixed seed and type-2 quantiles; rd2dnet uses
  repp=2000, `seed=None`. Consider a fixed default seed for reproducibility.
- **Off-network snapping cost for travel-time metrics** (currently free).

## Later decisions (morning of 2026-10-19)

12. **Counting locations instead of rows** for `bwcheck` and `min_obs` (panels). rd2d and
    rdrobust count unique observations for `bwcheck`; for `min_obs` I applied the same
    rule. Results change only with repeated locations.
13. **J = 1 critical value = normal quantile.** No simulation for one point.
14. **New warnings (no change to results):** pandas index not 0..n-1; treatment labels that
    contradict the boundary side; points far from the network (median snap distance
    above the network's extent); cluster-robust SEs with fewer than 10 clusters per side
    at some point (threshold chosen from the Monte Carlo: 5-6 clusters gave 0.70 coverage).
    Alternative for the clusters: a new status that drops such points. Rejected as too
    aggressive without your input.
15. **New errors for inputs that silently gave wrong results:** 2-D `y`/`weights`, negative
    or NaN `density()` weights, zero-length boundaries, NaN `eval_points`, `n_eval < 1`,
    negative or NaN `snap_scale`.
16. **Treated-side inference** uses only observations whose nearest point on the line
    through `b` is inside it (when at least 5 treated and 5 control qualify).
17. **`Precomputed.from_od` / `precomputed_from_od()`**: added (additive API) because the
    brief listed the missing OD constructor. Missing pairs are unreachable; `symmetric`
    is off by default.
18. **Not done, deliberately (need your call):** changing the density test's default
    normalisation, decoupling the p- and q-fit bases, deduplicating evaluation points in
    the default WBATE, a default off-network cost for travel-time metrics, the
    `from_place` offline fallback default, a "flow" direction for rivers, and `lbate`
    min/abs options.

Note: the Monte Carlo agent ran `taskkill /F /IM python.exe` once early on, which
stopped every Python process on this machine at that moment. Nothing of yours was
running as far as I can tell, but check if you had Python jobs open overnight.

## Night 2 (2026-10-20/21): decisions on the open questions

You asked me to decide these. Each was tested before it was adopted; details are in
CHANGELOG.md and the PRs.

19. **Sorting tests (density over-rejection on curved boundaries).**
    - Research: Bugni, Crippa and Restrepo (2026, arXiv:2609.00350) give a manipulation
      test for boundary designs and show why signed-distance density tests can fail.
      Sorting is a statement about where units locate in the plane, so the network
      should not enter the null.
    - Monte Carlo (500 reps per cell, n = 6k and 24k), comparing three options:

      | option | size with units in the plane | power |
      |---|---|---|
      | `density(normalize="euclidean")` (new) | 0.03-0.06 on straight, ring, wavy and irregular boundaries, river and panel | the most against one-sided sorting |
      | `"network"` (old default) | 0.34-1.00 on curved boundaries | not meaningful |
      | `"none"` | 0.08-0.33 | not meaningful |

    - The BCR test held its size everywhere in the plane. It is the only test with power
      against sorting that offsets along the boundary; q = 200 was best.
    - Decision: `density()` defaults to `"euclidean"`; new `manipulation_test()` /
      `rdnetmanip` with q = 200. Both warn when observations lie on the network itself,
      where no test held its size.
    - Not done: porting rddensity's bandwidth selector. Size is already correct with the
      current selector.
20. **Near-singular fits.**
    - Decision: replace the all-or-nothing `s` fallback with an "identifiable basis".
      Columns are kept in degree order; a column is dropped if it is aliased (relative
      residual < 1e-4) or would give an observation leverage 1.
    - Effects: a singular bias-correction fit no longer biases tau (2.83 → 2.5 exactly).
      On a one-way grid with an o² term pinned by one observation, the sd of tau_bc fell
      from 1.67 to 0.44.
21. **Default vce = HC3** (was HC1, rd2d's choice).
    - HC3 gave 0.945-0.951 pointwise coverage in every Monte Carlo design; HC1 gave 0.92-0.94.
    - In a high-leverage design, HC3 gave 0.95 and HC1 0.76.
    - `vce="hc1"` gives rd2d's behaviour; reference tests pin it.
22. **Off-network distance for non-length metrics.**
    - Charged at the snapped street's own cost per unit length; it was free.
    - A point at an end node uses the node's mean rate, with a relative tolerance so
      that 0 vs 1e-15 rounding does not decide.
23. **Rivers:** new `direction="flow"` (the shorter directed path either way).
24. **WBATE default:** equal weight per distinct place. Duplicates from snapping share
    one weight.
25. **LBATE:** `kind="min"` and `kind="absmax"` (Stata `estat lbate, min|absmax`).
26. **Default seed 0** for the uniform band, so results are reproducible (rd2d fixes its seed too).
27. **Points on the boundary line count as treated with `treated_side`**, as they already
    did with `treated_region`.
28. **`from_place` raises on download failure** (`offline_fallback=True` for demos).
29. **Kept as is, deliberately:**
    - `min_obs` default 20 + #terms(q): rd2d's `rd2d()` uses 52. With HC3 the coverage is
      already near nominal, and 52 would turn many sparse-network points "insufficient".
    - Standard CR1 cluster factor: rd2d's default multiplies an extra n/(n-k). The README
      warns about clustering with few clusters, and the code warns below 10 per side.
    - Python breaks exact ties between equidistant boundary nodes differently from R and
      Stata. This only matters on perfectly regular lattices (the mock city at h = 2000).
30. **Stata after the AppStream session expired (around 03:20 on 10/21).**
    - Everything up to and including the identifiable basis and HC3 was run in Stata:
      ALL TESTS PASSED and PARITY PASSED (40 random designs, worst difference 1.25e-8).
    - Later changes are untested in Stata: snap cost, place weights / lbate kinds / seed /
      on-line points, sorting tests (`rdnetmanip`), and flow. The cases for the next run
      are on the `parity-cases` branch (commit 5e4294f).
31. **Statistical review fixes (PR #20).**
    - `pooled=True` with h0 != h1: the common slopes are scaled by one bandwidth,
      max(h0, h1), on both sides. Side-specific scaling made the "common" slopes
      differ between the sides, so even a noiseless plane was not fitted exactly.
    - Singularity is judged on the column-equilibrated Gram matrix (and that matrix is
      inverted). Covariates in dollars made every point "singular" before.
    - `vce="conley"` warns when the cutoff exceeds the bandwidth: the HAC sums then
      cover most of the local sample, whose residuals sum to zero by construction, and
      the SEs collapse.
32. **A street counts as running along the boundary only if it lies on it** (its quarter
    points within `tol`), not merely because both ends do. A road that loops into one
    side and comes back gave its locations s = 0.
33. **Conley covariance: projected onto the PSD cone**, with a warning. Bartlett weights
    are positive definite on a line but not on a network distance (minimum eigenvalue
    -0.99 found at cutoff 400). Clipping only the diagonal could give se = 0. The
    projection leaves PSD matrices untouched (tolerance 1e-10 relative), so ordinary
    results do not change.
34. **Fuzzy designs: status "weak first stage" when |fs| < 2 fs_se.** The estimates stay
    in the table; the default WBATE and LBATE leave such points out, explicit WBATE
    weights may include them. Monte Carlo: 8-9% of points had |fs| < 2 fs_se, their
    se_rb reached 30, and the default WBATE SE exceeded 5 in a fifth of replications.
    A two-standard-error rule is the simplest honest flag; a weak-IV-robust
    (Anderson-Rubin) interval would be the next step.
35. **`rdbwnet(vce="cluster"|"conley")` selects with HC0 residuals**, as `rdnet` does
    internally with those standard errors, so `rdnet`'s bandwidths can be reproduced.
    The docstrings had promised HC1, which neither did consistently.
36. **Warnings, with thresholds from Monte Carlo runs:**
    - cluster-robust errors with fewer than 50 clusters per side (was 10; coverage was
      0.78 with 26 and 0.88 with 48);
    - HC errors when observations sit at fewer than 90% distinct locations (panels:
      HC3 coverage 0.49);
    - sorting tests when the locations' median distance to the nearest street is below
      0.6 times that of evenly spread points (a 45 x 45 lattice over the same extent,
      the same in all three languages). GPS noise of 3 / 8 / 15 m gives ratios of
      0.15 / 0.33 / 0.51, with rejection rates of 0.93 / 0.85 / 0.33 for the
      manipulation test. Parcels in the mock city give 0.90; uniform points about 1.
37. **Callable metrics:** the package cannot know which attributes a function reads, so
    it checks that the pieces of the first split edge add up to the original cost and
    warns otherwise (pointing to `extensive=`). Costs that read node attributes
    (`node_contrast`) are not checked. The boundary split cache is now keyed on
    `extensive` too.
38. **Cluster-robust errors: CR2 (Bell-McCaffrey) with Satterthwaite degrees of freedom by
    default; rd2d's CR1 stays available** (`cluster_vce="cr1"`, Stata `cr1`). This departs
    from rd2d on purpose. A Monte Carlo (`C:\Users\zwins\agent_cluster`, 1000 replications per
    cell, fixed bandwidth) found:
    - **Why CR1 fails:** the jump is an extrapolated intercept, so the 3 clusters nearest the
      point carry two thirds of the weight, at leverage around 0.5. Their residuals absorb
      their own cluster effect: the expected CR1 variance is 0.36-0.83 of the truth, and the
      effective number of clusters is about 5 even with 50 in the window.
    - **Coverage of 95% intervals** (panels with 16-100 sites per side; spatial blocks with
      6 and 20 per side; independent errors clustered by block):
      - CR1, normal: 0.66-0.91.
      - CR1, t(G-1): 0.71-0.91.
      - CR3, t(G-1): about 0.95 pointwise, but the uniform band reaches only 0.88-0.95.
      - Wider-fit residuals: 0.85 with 6 blocks.
      - **CR2 with Bell-McCaffrey df: 0.955-0.97 pointwise and 0.96-0.99 for the band.**
    - **The cost** is honest width: CR2 intervals are 1.2-3.2 times as long as CR1's. With
      6 blocks per side the degrees of freedom are about 2.
    - **Implementation:** an exact low-rank form, D_g = I + U M U' with U = [X_g, W_g X_g], so
      large clusters are cheap. The uniform band maps the normal sup-quantile to each point's
      t distribution. WBATE uses the smallest df among the weighted points.
    - The few-clusters warning (below 50 per side) now applies to CR1 only.
39. **Boundary splitting** now splits every edge that meets the boundary line at each
    crossing. Before, only edges whose ends were on opposite sides were split, and only at
    the first crossing. Stretches along the boundary contribute both of their ends. Stored
    geometries are oriented from u to v first. R pre-filters with
    d(u) + d(v) <= edge length, which is exact.
40. **Points on the boundary line's extension beyond its ends are treated**, whichever
    direction the line runs. They have no side, and the old rule flipped them when `b` was
    reversed.

## Nearest-crossing bins (2026-10-27)

Design choices made while implementing the bins design (an additional method;
the default is unchanged), with the owner's agreement on the main ones:

41. **Added, not replacing.** `bin_coordinates()` returns coordinates whose evaluation
    points are the crossings, so `rdnet` and every option (bias correction, the h_R
    selector, cluster/Conley/CR2, fuzzy, covariates) work unchanged; the default design's
    numbers are bit-identical to before (checked on the mock city, 7 option sets).
42. **Branch** = the edge where the path from the crossing first passes a junction
    (owner's choice). A crossing inserted by splitting an edge has one edge per side, so
    branches start at the first node with three or more neighbours (or at a crossing
    with two or more neighbours on that side).
43. **One bin per boundary node.** A first version merged crossings joined by a street
    along the boundary; on the mock city (boundary along a street) that made one bin.
    Locations on such a street belong to the nearer node.
44. **h_S:** rule of thumb `sqrt(2 chi2_k(0.95))` plus `similarity_sensitivity()`
    (owner's choice); a data-driven selector is noted in the README as future research.
45. **Edge features** are computed on the network as given: nodes the package inserted
    (crossings, locations) with two neighbours are skipped. Omega is the covariance over
    all segments (pseudo-inverse with a 1e-10 relative eigenvalue cut, the same in all
    three languages). Attributes are averaged over parallel edges and both directions.
46. **Portable tie rules.** The shortest-path tree's parent is the smallest-index
    neighbour on a shortest path (not the Dijkstra predecessor, whose ties depend on
    settling order), and a location exactly on an end node of its segment is placed on
    the node. The parity harness found the second: Python and R snapped such points to
    different incident segments.
47. **`rdnet_fe` slopes** are pooled across bins (interacted with the crossings'
    features and S, not by branch); in the Monte Carlo it has a +0.04 bias when slopes
    really differ by branch, documented in the README.
48. **GLS weights** may be negative with strongly correlated estimates; a warning says
    so rather than truncating them.
49. **Branch slopes need `min_branch` (default 10) locations on their side**; the reference
    slope is the trunk's, or else the largest branch's. Before, which aliased branch column
    a fit dropped depended on column order (two orders gave tau -0.37 and +0.38 in a parity
    case with no trunk and a one-location branch).
50. **Parity bins runs use `min_obs = 20`.** With 10 a bin fit at leverage near 1 gave
    se_rb = 305 in both languages, equal only to 5e-6: the known near-singular HC3 case
    (DECISIONS "Near-singular q-fits"), not specific to bins.
