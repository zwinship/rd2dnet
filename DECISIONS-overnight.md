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
