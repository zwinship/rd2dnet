# Roadmap

Status as of 2026-10-15: core package done and validated on the mock city (75 tests passing,
Monte Carlo results in the README). Three agreed next steps remain. Scope rule: add a feature
only if rd2d has it. Everything else is documented as a recipe, not built in.

## 1. rd2d parity: aggregate effects

rd2d's `summary()` (`rd2d/results.py`, `_append_aggregate_rows`) offers:

- **WBATE**: weighted boundary average treatment effect, `w' tau_q` with
  `se = sqrt(w' Sigma w)`, where `Sigma` is the covariance of the bias-corrected
  estimates across evaluation points.
- **LBATE**: largest BATE. The point estimate is `max_j tau_j`; the CI is
  `[max_j (tau_q - cval * se), max_j (tau_q + cval * se)]` with `cval` the
  uniform-band critical value.

To do in rdnet:
- Always store the joint covariance across evaluation points on `RDNetResult`. It is
  currently `cb_cov`, which is only computed when `cbands=True` and only over `ok`
  points. Make it J x J with NaN rows for points that weren't estimated, as rd2d's
  `params_cov` does.
- Add `RDNetResult.wbate(weights)` and `RDNetResult.lbate()` (or `summary(WBATE=..., LBATE=...)`).
- Tests: WBATE with one-hot weights equals that point's estimate and SE; LBATE CI
  contains the max of the pointwise estimates.

## 2. Engineering: scale

Everything has only been run on the 450-node mock grid. A real city is about 50,000 nodes.
- Build a synthetic large grid (for example 200 x 250) with a few thousand to 100k points.
- Profile `network_coordinates`, `rdnet` (with `vce="hc1"` and `conley`), `rdbwnet`,
  `density`, and `TransitionCost` (its state-expanded graph is roughly |E| x classes).
- Likely hot spots: `TransitionCost.compile` (Python loops), `hac_quadratic` (block
  size, network Dijkstra per block), `_street_length_transform` (Python loop over
  edges), and `network_coordinates` edge snapping.
- Optionally run the `online`-marked OSM test once: `pytest -m online`.

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

## Getting started on another machine

```bash
gh auth login                      # once per machine
git clone https://github.com/zwinship/rdnet.git && cd rdnet
python -m venv .venv
.venv\Scripts\activate            # Windows (source .venv/bin/activate elsewhere)
pip install -e ".[dev]"
pytest
```
