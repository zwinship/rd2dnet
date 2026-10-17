# Roadmap

Status as of 2026-10-16: core package, rd2d aggregate effects and scale work done (87 tests
passing; Monte Carlo and timing results in the README). Steps 3-4 remain. Scope rule: add a feature
only if rd2d has it. Everything else is documented as a recipe, not built in.

## 1. rd2d parity: aggregate effects (done 2026-10-16)

`RDNetResult.wbate(weights)`, `RDNetResult.lbate()` and `summary(wbate=, lbate=)`, matching
rd2d's `summary(WBATE=, LBATE=)` (`rd2d/results.py`, `_append_aggregate_rows`).
`RDNetResult.cov` (was `cb_cov`) is always the J x J covariance with NaN rows for points
not estimated, as rd2d's `params_cov`; `crit` (was `cb_crit`) is stored even without bands.
The cluster covariance now uses each point's own `g/(g-1)` so its diagonal equals
`se_rb**2`. Unlike rd2d (which drops aggregate SEs if any point is missing), masked points
take weight 0. Tests and Monte Carlo results are in the README.

## 2. Engineering: scale (done 2026-10-16)

Profiled on synthetic city (50k and 200k nodes), suburb and river networks
(`benchmarks/`) with up to 5M GPS pings and 10M sensor readings; results in the README's
Scale section. Fixed: the quadratic density transform and dense Conley HAC (now
sparse, bounded local Dijkstra, threaded row blocks, one pass per `rdnet` call),
per-location coordinates for repeated measurements, vectorised metric compilation
(including `TransitionCost`), edge extraction and boundary splitting, cached
bandwidth designs, and threaded geometry and evaluation points. Remaining costs are
inherent: Conley grows with the number of pairs within the cutoff, and rd2d's
median-distance bias pilot makes `rdbwnet` linear in n per evaluation point.
Not done: the `online` OSM test.

## 3. Polish

- Tutorial notebook (`examples/`) on the mock city, covering: metrics (length, travel
  time, `TransitionCost` with ramps, `Precomputed`), coordinates, `rdnet`/`rdbwnet`,
  fuzzy, covariates, `deriv`/`tangvec`, Conley, and the density test.
- Include these **recipes** (no new API needed):
  - Donut: `y_d = np.where(abs(coords.s) < delta, np.nan, y)`. Missing outcomes are
    skipped, and `delta` is in metric units.
  - Covariate balance: pass a pre-determined characteristic as `y`.
  - Placebo boundary: build a `Boundary` at a fake line and run on one side's parcels only.
  - Difference-in-discontinuities: for a panel, run on `y_post - y_pre`. For independent
    samples, difference the estimates and use `se = sqrt(se_pre**2 + se_post**2)`.
    Overlapping samples would need cross-run covariance, which rd2d doesn't offer either.
- Fill out docstrings. Add GitHub Actions CI running `pytest` on push.

## 4. Release and ports (after 1-3)

Done already: renamed to `rd2dnet` (`rdnet` was taken on PyPI by an unrelated
deep-learning model; `rd2dnet` is free on PyPI and CRAN as of 2026-10-15), GPL-3.0
LICENSE added (required, because the bandwidth and density code are ported from rd2d
and rddensity, both GPL-3.0-only), and rd2d/rddensity attribution and references added
to the README and module docstrings.

- [ ] Contact Cattaneo, Titiunik and Yu (rdpackages) before a public release, as a
      courtesy: the package builds directly on their method and code.
- [ ] Make the GitHub repo public when ready.
- [ ] **PyPI**: no review; live within minutes. `python -m build`, `twine check dist/*`
      (both passed on 2026-10-15), then `twine upload dist/*`. Try TestPyPI first.
- [ ] Optional peer review: pyOpenSci (about 1-3 months), JOSS (about 1-3 months), or
      Journal of Statistical Software (often a year or more).
- [ ] **R port -> CRAN**: first response in days to about 2 weeks; a revision round is
      common for new packages. Needs a clean `R CMD check --as-cran`, every export
      documented, fast examples, offline tests (the bundled mock city covers this) and
      the license. Likely stack: igraph or dodgr, sf, osmdata, testthat. Estimate: one
      to a few weeks; most of the effort is numerical parity with Python.
- [ ] **Stata port -> SSC** (submit to Kit Baum): light review; posted within days to
      a couple of weeks. The Stata Journal is peer-reviewed and takes months to over a
      year. Options: native Mata, which needs a hand-written Dijkstra and user-supplied
      edge lists (weeks of work), or a wrapper calling the Python package via Stata
      16+'s Python integration (fast to build, but users need a working Python).

## Getting started on another machine

```bash
gh auth login                      # once per machine
git clone https://github.com/zwinship/rd2dnet.git && cd rd2dnet
python -m venv .venv
.venv\Scripts\activate            # Windows (source .venv/bin/activate elsewhere)
pip install -e ".[dev]"
pytest
```
