"""Generate examples/tutorial.ipynb (run from the repo root)."""
import nbformat as nbf

cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip()))  # noqa: E731
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip()))  # noqa: E731

md("""
# rd2dnet tutorial

`rd2dnet` estimates boundary discontinuity designs on networks. It carries the
location-based framework of **rd2d** (Cattaneo, Titiunik and Yu) over from Euclidean
coordinates to network coordinates and user-defined network metrics: side-specific local
polynomials, robust bias-corrected inference, uniform bands, data-driven bandwidths and
aggregate effects. Please cite the rd2d papers (see the README) when you use it.

This notebook uses the bundled **mock city**, a small synthetic street grid with a
boundary along "60th Street" (treated side to the south), a highway with ramps, and a park
on the boundary. It is only an illustration; nothing here depends on that application.
Everything runs offline in about a minute.

Contents
1. Data and network
2. Distance is a metric you choose
3. Network coordinates
4. Estimation: bandwidths, BATEC curve, uniform band
5. Aggregate effects: WBATE and LBATE
6. Covariates, fuzzy designs, derivatives
7. Standard errors: cluster and Conley
8. Density (sorting) test
9. Recipes: donut, covariate balance, placebo boundary, difference-in-discontinuities
10. Large data
""")

code("""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import LineString, box

import rd2dnet
from rd2dnet.datasets import load_mock_midtown, make_mock_midtown

city = load_mock_midtown()
obs = city.obs
y = obs["outcome"].to_numpy()
obs.head()
""")

md("""
## 1. Data and network

`city.graph` is a networkx graph whose nodes carry projected `x`/`y` coordinates (metres)
and whose edges carry `length`, `highway` class, `speed_kph` and `travel_time`.
`city.boundary` is a `rd2dnet.Boundary`: a line, the treated side, and a **mask** for the
park frontage, where one side has no parcels and no effect can be estimated.
""")

code("""
fig, ax = plt.subplots(figsize=(7, 6))
for u, v in city.graph.edges():
    a, b = city.graph.nodes[u], city.graph.nodes[v]
    ax.plot([a["x"], b["x"]], [a["y"], b["y"]], color="0.75", lw=0.6)
ax.scatter(obs["x"], obs["y"], s=2, c=np.where(obs["treated"], "tab:red", "tab:blue"))
ax.plot(*city.boundary.line.xy, color="k", lw=2, label="boundary")
ax.fill(*city.park.exterior.xy, color="tab:green", alpha=0.3, label="park (masked)")
ax.set_aspect("equal"); ax.legend(loc="upper right"); ax.set_title("mock city: treated (red) / control (blue)")
plt.show()
city.boundary
""")

md("""
## 2. Distance is a metric you choose

A `NetworkGraph` is a graph **plus a metric**. Distance means dissimilarity and does not
have to be metres. Coordinates, kernel distances and bandwidths are all in the metric's
units, while the treatment side stays geographic.
""")

code("""
G = city.graph
metrics = {
    "length (m)": "length",
    "travel time (s)": rd2dnet.travel_time(),
    # 50 extra units per unit change in floor-area ratio between adjacent nodes
    "length + land-use contrast": rd2dnet.node_contrast("far", 50.0),
    # switching road class costs 500 m; ramps are transparent connectors, so
    # highway -> ramp -> highway is cheap but highway -> ramp -> street is not
    "length + class switches": rd2dnet.TransitionCost(
        "length", class_attr="highway", default_penalty=500.0, connectors={"motorway_link"}),
}
pts = obs[["x", "y"]].to_numpy()[:5]
pd.DataFrame({name: rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G, m), city.boundary, pts, n_eval=5).s
              for name, m in metrics.items()}).round(1)
""")

md("""
Any node-by-node dissimilarity matrix works too (`Precomputed`). Observations then snap to
the nearest node, because per-edge costs are unknown.
""")

code("""
base = rd2dnet.NetworkGraph(G, "length")
D = base.distances_from(np.arange(base.n_nodes))     # or your own matrix, in base.nodes order
pre = rd2dnet.NetworkGraph(G, rd2dnet.Precomputed(np.sqrt(D)))  # e.g. a concave transform
rd2dnet.network_coordinates(pre, city.boundary, pts, n_eval=5).s.round(2)
""")

md("""
## 3. Network coordinates

`network_coordinates` snaps each observation to the street network **on its own side**
and computes, as rd2d's `(X1 - b1, X2 - b2)`:

* `s`: signed network distance to the boundary (positive on the treated side);
* at each evaluation point `b`: `d(b)`, the network distance used for kernel weights, and
  `o(b)`, the signed offset along the boundary.

Evaluation points default to an even grid along the boundary; masked ones are skipped.
""")

code("""
net = city.network("length")
coords = rd2dnet.network_coordinates(net, city.boundary, obs[["x", "y"]], n_eval=21)
print(coords)
coords.frame(j=10).head()
""")

code("""
coords.eval_table().head(8)
""")

md("""
## 4. Estimation

`rdbwnet` ports rd2d's direct plug-in bandwidth selector (`mserd` by default; also
`msetwo`, `imserd`, `imsetwo` and the coverage-error-optimal `cer*` variants). `rdnet`
selects bandwidths itself when `h` is not given, fits order `p` (default 1) for the point
estimate and order `q = p + 1` for robust bias-corrected inference, and adds a uniform band.
""")

code("""
bw = rd2dnet.rdbwnet(y, coords)
bw[["arc", "h0", "h1", "n0", "n1", "status"]].head()
""")

code("""
res = rd2dnet.rdnet(y, coords, seed=0)
print(res)
""")

code("""
ax = res.plot(true_tau=city.true_tau(res.table["arc"]))
ax.set_title("BATEC along the boundary: estimate, robust CIs, uniform band, truth")
plt.show()
""")

md("""
Bandwidths can also be given directly: a scalar, `(h_control, h_treated)`, one per
evaluation point, or a `(J, 2)` array. `weights=` multiplies the kernel weights (for
example floor area), and `pooled=True` fits a common-slope model.
""")

code("""
res_fixed = rd2dnet.rdnet(y, coords, h=600.0, weights=obs["far"], seed=0)
res_fixed.table[["arc", "tau_bc", "se_rb", "status"]].head()
""")

md("""
## 5. Aggregate effects: WBATE and LBATE

As in rd2d's `summary(WBATE=..., LBATE=True)`:

* **WBATE** `w' tau`: weights are normalised to sum to one, and the standard error
  `sqrt(w' Sigma w)` uses the joint covariance `res.cov` of the estimates across
  evaluation points. Points that were not estimated (masked) must get weight 0.
* **LBATE** `max_j tau(b_j)`: its interval uses the uniform critical value, so it
  inherits the band's joint coverage (conservative for the maximum).
""")

code("""
ok = res.table["status"].eq("ok").to_numpy()
print(res.summary(wbate=ok.astype(float), lbate=True).split("aggregate effects")[1])
""")

code("""
# Weight by frontage instead: here, by the number of parcels near each point.
w = np.where(ok, res.table["n0"] + res.table["n1"], 0.0)
res.wbate(w)
""")

code("""
# Averaging correlated neighbours helps less than independence would suggest:
se_naive = np.sqrt(((ok / ok.sum()) ** 2 * res.table["se_rb"].fillna(0) ** 2).sum())
print(f"WBATE se with covariance {res.wbate(ok.astype(float))['se_rb']:.3f}  vs  ignoring it {se_naive:.3f}")
""")

md("""
## 6. Covariates, fuzzy designs, derivatives

**Covariates** enter with coefficients common to both sides (rd2d's `covs_eff`) and are
estimated jointly with the local polynomials. In the mock city the outcome depends only
weakly on log floor-area ratio (it explains about 1% of the noise), so adjusting for it
leaves estimates and standard errors almost unchanged. A covariate that predicts the
outcome well shrinks the standard errors.
""")

code("""
res_cov = rd2dnet.rdnet(y, coords, h=600.0, covs=np.log(obs["far"]), seed=0)
res_raw = rd2dnet.rdnet(y, coords, h=600.0, seed=0)
pd.DataFrame({"tau_bc": res_raw.table["tau_bc"], "tau_bc (covs)": res_cov.table["tau_bc"],
              "se_rb": res_raw.table["se_rb"], "se_rb (covs)": res_cov.table["se_rb"]}).dropna().round(3).head()
""")

md("""
**Fuzzy designs**: pass take-up `fuzzy=D`. The estimate is the ratio of the outcome jump
to the take-up jump at each point; the first stage is reported as `fs`/`fs_se`. Below,
take-up is 75% on the treated side and 25% on the other (true first stage 0.5) and the
effect of take-up is 2, in a larger simulated city (12,000 parcels). One draw scatters
around the truth, and neighbouring evaluation points move together because they share
observations; over 200 independent take-up draws the first stage averages 0.503, its
standard error matches the Monte Carlo spread, and 95% intervals cover 93% of the time.
""")

code("""
big = make_mock_midtown(n_obs=12000, seed=11)
big_coords = rd2dnet.network_coordinates(big.network(), big.boundary, big.obs[["x", "y"]], n_eval=21)
rng = np.random.default_rng(2)
T = big_coords.treated
D_takeup = np.where(T, rng.uniform(size=len(T)) < 0.75, rng.uniform(size=len(T)) < 0.25).astype(float)
y_fuzzy = 0.0008 * (big.obs["y"] - 1600) + 2.0 * D_takeup + rng.normal(0, 0.5, len(T))
res_fz = rd2dnet.rdnet(y_fuzzy, big_coords, fuzzy=D_takeup, seed=0)
res_fz.table[["arc", "fs", "fs_se", "tau", "ci_lo", "ci_hi", "status"]].dropna().round(3).head()
""")

md("""
**Derivatives**: `deriv=(1, 0)` is the jump in the slope *across* the boundary (a kink),
`deriv=(0, 1)` the jump in the slope *along* it, and `tangvec=(t_s, t_o)` any directional
derivative. The mock outcome has the same slope on both sides, so the true kink is zero
and the estimates sit near zero.
""")

code("""
kink = rd2dnet.rdnet(y, coords, h=800.0, deriv=(1, 0), seed=0)
kink.table[["arc", "tau_bc", "ci_lo", "ci_hi", "status"]].dropna().round(5).head()
""")

md("""
## 7. Standard errors

`vce` takes `hc0`-`hc3`; `cluster=` gives cluster-robust errors; `conley_cutoff=` gives a
Bartlett spatial HAC over **network** distance (or `conley_distance="euclidean"`). The
same choice drives the uniform band and WBATE/LBATE through `res.cov`.

A caution on `cluster=`: as in rd2d, cluster-robust errors use the residuals of the local
fit. When clusters are small relative to the bandwidth (here, city blocks inside a 600 m
window), the local polynomial absorbs much of the within-cluster variation and the
cluster standard errors come out too small, even with independent errors (see the
README's known limitations for a Monte Carlo). For spatially correlated shocks, the
Conley option, whose residuals come from a wider auxiliary fit, is the validated choice.
""")

code("""
y_ = y
fits = {
    "HC3 (default)": rd2dnet.rdnet(y_, coords, h=600.0, seed=0),
    "cluster (block)": rd2dnet.rdnet(y_, coords, h=600.0, cluster=obs["block"], seed=0),
    "Conley 150 m (network)": rd2dnet.rdnet(y_, coords, h=600.0, conley_cutoff=150.0, seed=0),
    "Conley 150 m (Euclidean)": rd2dnet.rdnet(y_, coords, h=600.0, conley_cutoff=150.0,
                                              conley_distance="euclidean", seed=0),
}
pd.DataFrame({k: r.table["se_rb"] for k, r in fits.items()}).dropna().round(3).head()
""")

md("""
## 8. Density (sorting) test

The local polynomial density test of Cattaneo, Jansson and Ma (`rddensity`). By default
position is the signed straight-line distance to the boundary (`normalize="euclidean"`),
which keeps its size at curved boundaries. `normalize="network"` measures position in
cumulative street length instead (more parcels per metre of street on one side); it is
valid only at straight boundaries and over-rejects badly at curved ones.
`manipulation_test()` is the test of Bugni, Crippa and Restrepo (2026): among the `q`
locations nearest the boundary, are treated and control equally likely everywhere along
it? Both tests assume units spread over the plane, not points on the streets.
""")

code("""
print(rd2dnet.density(coords))
print(rd2dnet.manipulation_test(coords, seed=0))
""")

md("""
## 9. Recipes

These need no special API; they combine the pieces above.

### Donut
Drop observations within `delta` (in metric units) of the boundary. Missing outcomes are
skipped.
""")

code("""
delta = 40.0
y_donut = np.where(np.abs(coords.s) < delta, np.nan, y)
donut = rd2dnet.rdnet(y_donut, coords, h=700.0, seed=0)
pd.DataFrame({"full": res_fixed.table["tau_bc"], "donut 40 m": donut.table["tau_bc"]}).dropna().round(3).head()
""")

md("""
### Covariate balance
Pass a pre-determined characteristic as the outcome; there should be no jump. In the mock
city, floor-area ratio is unrelated to treatment.
""")

code("""
bal = rd2dnet.rdnet(np.log(obs["far"]), coords, seed=0)
ok_b = bal.table["status"].eq("ok").to_numpy()
print(bal.summary(wbate=ok_b.astype(float), lbate=True).split("aggregate effects")[1])
""")

md("""
### Placebo boundary
Build a `Boundary` along a fake line inside one side and use only that side's
observations. Here: a street 8 blocks north, on the control side. It runs through the
park, which is masked.
""")

code("""
from rd2dnet.datasets import BOUNDARY_Y, PARK
fake_y = BOUNDARY_Y + 8 * 80.0
north = obs["y"].to_numpy() > BOUNDARY_Y
fake = rd2dnet.Boundary(LineString([(-400, fake_y), (2600, fake_y)]),
                        treated_region=box(-1e5, -1e5, 1e5, fake_y), mask=[PARK])
pc = rd2dnet.network_coordinates(city.network(), fake, obs.loc[north, ["x", "y"]], n_eval=15)
placebo = rd2dnet.rdnet(y[north], pc, seed=0)
ok_p = placebo.table["status"].eq("ok").to_numpy()
print(placebo.summary(wbate=ok_p.astype(float)).split("aggregate effects")[1])
""")

md("""
### Difference-in-discontinuities
With a **panel** (the same units before and after), run on `y_post - y_pre`. With
**independent samples**, estimate each period with the same evaluation points and
bandwidths and difference the estimates, with `se = sqrt(se_pre**2 + se_post**2)`.
Overlapping but non-identical samples would need the cross-period covariance, which rd2d
does not provide either.
""")

code("""
rng = np.random.default_rng(5)
f = y - city.true_tau(coords.crossing_arc) * coords.treated           # outcome without the effect
y_pre = f + rng.normal(0, 0.5, len(y))                               # same parcels, earlier period
y_post = y + rng.normal(0, 0.5, len(y))
did_panel = rd2dnet.rdnet(y_post - y_pre, coords, h=700.0, seed=0)

# independent samples: two cross-sections with the same evaluation points and bandwidth
c_pre, c_post = make_mock_midtown(n_obs=3000, seed=21), make_mock_midtown(n_obs=3000, seed=22)
co_pre = rd2dnet.network_coordinates(c_pre.network(), c_pre.boundary, c_pre.obs[["x", "y"]], n_eval=21)
co_post = rd2dnet.network_coordinates(c_post.network(), c_post.boundary, c_post.obs[["x", "y"]], n_eval=21)
y0 = c_pre.obs["outcome"] - c_pre.true_tau(co_pre.crossing_arc) * co_pre.treated
r_pre = rd2dnet.rdnet(y0, co_pre, h=700.0, seed=0).table
r_post = rd2dnet.rdnet(c_post.obs["outcome"], co_post, h=700.0, seed=0).table
did = pd.DataFrame({"arc": r_post["arc"],
                    "panel": did_panel.table["tau_bc"],
                    "independent": r_post["tau_bc"] - r_pre["tau_bc"],
                    "se (independent)": np.sqrt(r_pre["se_rb"] ** 2 + r_post["se_rb"] ** 2),
                    "truth": city.true_tau(r_post["arc"])})
did.dropna().round(3).head()
""")

md("""
## 10. Large data

* **Repeated measurements** (sensor panels, many readings per site) are snapped and
  measured once per distinct location: `coords.n_locations` vs `coords.n_obs`.
* Snapping, boundary geometry, the Conley kernel and bandwidth selection use several
  threads; set the environment variable `RD2DNET_NUM_THREADS` to cap them (`1` disables).
* Conley cost grows with the number of observation pairs closer than the cutoff, which
  grows with the square of local density. With millions of nearby points, a smaller cutoff
  or clustered errors are much cheaper.
* `benchmarks/profile_scale.py` times every stage on synthetic city, suburb and river
  networks with up to millions of observations; see the README's Scale section.
""")

code("""
panel_xy = np.repeat(obs[["x", "y"]].to_numpy()[:500], 20, axis=0)    # 500 sensors x 20 readings
pcoords = rd2dnet.network_coordinates(net, city.boundary, panel_xy, n_eval=11)
pcoords.n_obs, pcoords.n_locations
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "examples/tutorial.ipynb")
print("written", len(cells), "cells")
