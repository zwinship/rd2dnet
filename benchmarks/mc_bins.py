"""Monte Carlo for the nearest-crossing bins design (README: Alternative design). Usage: python benchmarks/mc_bins.py [reps]"""
import sys, warnings, numpy as np
import pathlib; ROOT = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "tests"))
warnings.simplefilter("ignore")
import rd2dnet
from test_bins import comb, points_on
G, bnd = comb(n_streets=15, spacing=1500.0)
rng0 = np.random.default_rng(0)
xy = points_on(G, 70, rng0)
coords = rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G, "length"), bnd, xy, n_eval=3)
b = rd2dnet.bin_coordinates(coords, features=["lanes"])
M = b.n_eval
tau_m = 1.0 + 0.5 * np.sin(np.arange(M))
T = coords.treated.astype(float); s = coords.s
brs = np.unique(b.branch[b.branch >= 0]); sl = {-1: 0.002}; sl.update({x: 0.004 * np.cos(k) for k, x in enumerate(brs)})
mu = 2.0 * np.cos(np.arange(M))[b.bin] + tau_m[b.bin] * T + np.array([sl[x] for x in b.branch]) * s
R = int(sys.argv[1]) if len(sys.argv) > 1 else 300
rows = []
for r in range(R):
    rng = np.random.default_rng(100 + r)
    y = mu + rng.normal(0, 1, len(mu))
    res = rd2dnet.rdnet(y, b, h=450.0, min_obs=10)
    ok = res._ok()
    cov_pt = np.mean(((res.table.ci_lo <= tau_m) & (tau_m <= res.table.ci_hi))[ok])
    wd = res.wbate("density"); wg = res.wbate("gls")
    n = (res.table.n0 + res.table.n1).to_numpy(float) * ok
    true_d = n @ tau_m / n.sum()
    fe = rd2dnet.rdnet_fe(y, b, h=450.0)
    W = np.clip(1 - np.abs(s) / 450.0, 0, None) * b.kernel_weights
    mass = np.bincount(b.bin, weights=W, minlength=M); true_fe = mass @ tau_m / mass.sum()
    rows.append([cov_pt, wd.tau_bc - true_d, wd.se_rb, wd.ci_lo <= true_d <= wd.ci_hi,
                 wg.tau_bc - tau_m.mean(), fe.tau_bc - true_fe, fe.se_rb, fe.ci_lo <= true_fe <= fe.ci_hi, ok.sum()])
a = np.array(rows, float)
print(f"reps {R}, bins ok {a[:,8].mean():.1f}/{M}")
print(f"pointwise coverage {a[:,0].mean():.3f}")
print(f"WBATE density: bias {a[:,1].mean():+.4f} sd {a[:,1].std():.4f} mean se {a[:,2].mean():.4f} coverage {a[:,3].mean():.3f}")
print(f"GLS: bias vs mean tau {a[:,4].mean():+.4f} sd {a[:,4].std():.4f}")
print(f"FE: bias {a[:,5].mean():+.4f} sd {a[:,5].std():.4f} mean se {a[:,6].mean():.4f} coverage {a[:,7].mean():.3f}")
