"""Cross-language parity cases for rd2dnet (Python, R and Stata).

Generates random designs (networks, boundaries, observations and estimation
options), runs the Python package on each, and writes

    cases/<id>/nodes.csv, edges.csv, obs.csv, boundary.csv[, region.csv]
    cases/<id>/py_*.csv        Python results (the reference)
    cases/<id>/run.R           the same analysis in R   (writes r_*.csv)
    cases/<id>/run.do          the same analysis in Stata (compares with py_*.csv)

    python parity/make_cases.py [--n 24] [--seed 1] [--out parity/cases]

Then ``Rscript parity/run_r.R parity/cases`` and ``python parity/compare.py
parity/cases`` compare R with Python, and in Stata ``do parity/run_stata.do``
(from the repository root) runs every case and prints the largest differences.

Designs avoid exact ties between equidistant boundary nodes (node positions
are jittered and lengths are not multiples of each other), because Python
breaks such ties differently from R and Stata.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import LineString, Polygon

import rd2dnet

STATUS = {"ok": 0, "masked": 1, "insufficient": 2, "singular": 3, "no bandwidth": 4,
          "weak first stage": 5, "few clusters": 6}
BASIS = {"s,o": 1, "s": 2}
RES_COLS = ["h0", "h1", "n0", "n1", "tau", "se", "tau_bc", "se_rb", "status", "basis"]


# --------------------------------------------------------------------------
# designs
# --------------------------------------------------------------------------

def grid_network(rng, n=11, spacing=100.0, jitter=25.0, drop=0.08):
    ids = np.arange(n * n)
    gx, gy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    x = gx.ravel() * spacing + rng.uniform(-jitter, jitter, n * n)
    y = gy.ravel() * spacing + rng.uniform(-jitter, jitter, n * n)
    edges = []
    for i in range(n):
        for j in range(n):
            a = i * n + j
            if i + 1 < n:
                edges.append((a, (i + 1) * n + j))
            if j + 1 < n:
                edges.append((a, i * n + j + 1))
    keep = rng.uniform(size=len(edges)) > drop
    edges = [e for e, k in zip(edges, keep) if k]
    return pd.DataFrame({"id": ids + 1, "x": x, "y": y}), edges


def river_network(rng, n_rivers=3, span=1200.0):
    """Rivers crossing a vertical border at x = 0, with tributaries."""
    xs, ys, edges = [], [], []
    for r in range(n_rivers):
        y0 = -span / 2 + (r + 0.5) * span / n_rivers + rng.uniform(-60, 60)
        k = 14
        px = np.linspace(-700, 700, k) + rng.uniform(-20, 20, k)
        py = y0 + np.cumsum(rng.uniform(-25, 25, k))
        base = len(xs)
        xs += list(px)
        ys += list(py)
        edges += [(base + t, base + t + 1) for t in range(k - 1)]
        for _ in range(2):  # tributaries
            at = base + int(rng.integers(1, k - 1))
            tx = xs[at] + rng.uniform(-150, 150)
            ty = ys[at] + rng.choice([-1, 1]) * rng.uniform(120, 220)
            xs.append(tx)
            ys.append(ty)
            edges.append((at, len(xs) - 1))
    nodes = pd.DataFrame({"id": np.arange(len(xs)) + 1, "x": xs, "y": ys})
    return nodes, edges


def edge_frame(rng, nodes, edges, directed):
    xy = nodes[["x", "y"]].to_numpy()
    rows = []
    for a, b in edges:
        ln = float(np.hypot(*(xy[a] - xy[b]))) * rng.uniform(1.0, 1.25)
        speed = rng.choice([20.0, 30.0, 50.0])
        cls = rng.choice(["a", "b", "c"], p=[0.5, 0.35, 0.15])
        if directed and rng.uniform() < 0.3:
            pair = [(a, b)] if rng.uniform() < 0.5 else [(b, a)]
        else:
            pair = [(a, b), (b, a)] if directed else [(a, b)]
        for u, v in pair:
            rows.append({"from": u + 1, "to": v + 1, "length": ln,
                         "cost2": ln / speed * 3.6, "cls": cls})
    return pd.DataFrame(rows)


def make_design(rng, kind):
    if kind == "river":
        nodes, edges = river_network(rng)
        line = np.array([[0.0, -900.0], [0.0, 900.0]])
        bnd = {"line": line, "side": "right"}
        box = (-700, 700, -800, 800)
    elif kind == "ring":
        nodes, edges = grid_network(rng, n=13, spacing=100.0)
        c = np.array([600.0, 600.0]) + rng.uniform(-30, 30, 2)
        ang = np.linspace(0, 2 * np.pi, 13)[:-1] + rng.uniform(0, 0.3)
        rad = 330 + rng.uniform(-40, 40, 12)
        ring = np.column_stack([c[0] + rad * np.cos(ang), c[1] + rad * np.sin(ang)])
        ring = np.vstack([ring, ring[:1]])
        bnd = {"line": ring, "region": ring}
        box = (0, 1200, 0, 1200)
    else:
        nodes, edges = grid_network(rng)
        y0 = 500 + rng.uniform(-30, 30)
        tilt = rng.uniform(-0.15, 0.15)
        line = np.array([[-80.0, y0 - 580 * tilt], [1080.0, y0 + 580 * tilt]])
        if rng.uniform() < 0.5:
            mid = np.array([[500.0, y0 + rng.uniform(-60, 60)]])
            line = np.vstack([line[:1], mid, line[1:]])
        if rng.uniform() < 0.5:
            bnd = {"line": line, "side": rng.choice(["left", "right"])}
        else:
            region = np.array([[-1e5, -1e5], [1e5, -1e5], [1e5, line[-1, 1]]]
                              + [list(p) for p in line[::-1]] + [[-1e5, line[0, 1]]])
            bnd = {"line": line, "region": region}
        box = (0, 1000, 0, 1000)
    return nodes, edges, bnd, box


# --------------------------------------------------------------------------
# options
# --------------------------------------------------------------------------

def random_spec(rng, kind, case):
    directed = bool(kind == "grid" and rng.uniform() < 0.25)
    metric = rng.choice(["length", "length", "cost2", "transition"]) if kind != "river" else "length"
    spec = {
        "case": case,
        "kind": kind,
        "directed": directed,
        "direction": str(rng.choice(["both", "out", "in"])) if directed else "both",
        "metric": str(metric),
        "penalty": float(np.round(rng.uniform(50, 300), 1)),
        "n_eval": int(rng.choice([5, 7, 9])),
        "n_obs": int(rng.choice([700, 1000, 1400])),
        "runs": [],
        "density": [],
        "bw": [],
    }
    kernels = ["triangular", "epanechnikov", "uniform"]
    vces = ["hc0", "hc1", "hc2", "hc3"]
    for _ in range(3):
        r = {"kernel": str(rng.choice(kernels)), "vce": str(rng.choice(vces)), "p": 1}
        u = rng.uniform()
        if u < 0.35:
            r["h"] = float(np.round(rng.uniform(250, 450), 1))
        elif u < 0.5:
            r["h"] = [float(np.round(rng.uniform(250, 450), 1)), float(np.round(rng.uniform(250, 450), 1))]
        else:
            r["bwselect"] = str(rng.choice(["mserd", "msetwo", "imserd", "imsetwo", "cerrd", "certwo"]))
        extra = rng.choice(["none", "cluster", "conley", "fuzzy", "covs", "weights", "p2", "deriv"])
        if extra == "cluster":
            r["cluster"] = True
            r.pop("vce")
        elif extra == "conley":
            r["conley"] = float(np.round(rng.uniform(80, 200), 1))
            r.pop("vce")
        elif extra == "fuzzy":
            r["fuzzy"] = True
        elif extra == "covs":
            r["covs"] = True
        elif extra == "weights":
            r["weights"] = True
        elif extra == "p2":
            r["p"] = 2
        elif extra == "deriv" and kind != "river":
            r["deriv"] = [1, 0]
        spec["runs"].append(r)
    spec["density"].append({"normalize": str(rng.choice(["network", "none"])),
                            "bwselect": str(rng.choice(["common", "each"]))})
    spec["bw"].append({"bwselect": str(rng.choice(["mserd", "msetwo", "imserd", "cerrd"]))})
    return spec


# --------------------------------------------------------------------------
# Python
# --------------------------------------------------------------------------

def py_metric(spec):
    if spec["metric"] == "cost2":
        return rd2dnet.EdgeCost("cost2")
    if spec["metric"] == "transition":
        return rd2dnet.TransitionCost("length", class_attr="cls", default_penalty=spec["penalty"],
                                      connectors=["c"])
    return "length"


def py_boundary(bnd):
    line = LineString(bnd["line"])
    if "region" in bnd:
        return rd2dnet.Boundary(line, treated_region=Polygon(bnd["region"]))
    return rd2dnet.Boundary(line, treated_side=bnd["side"])


def res_frame(table):
    t = table.copy()
    out = pd.DataFrame({c: t[c] for c in ["h0", "h1", "n0", "n1", "tau", "se", "tau_bc", "se_rb"]})
    out["status"] = t["status"].map(STATUS)
    out["basis"] = t["basis"].map(BASIS) if "basis" in t else np.nan
    return out


def run_python(d, spec, nodes, edges, obs, bnd):
    net = rd2dnet.NetworkGraph.from_frames(nodes, edges, py_metric(spec), directed=spec["directed"],
                                           direction=spec["direction"])
    cc = rd2dnet.network_coordinates(net, py_boundary(bnd), obs[["x", "y"]].to_numpy(), n_eval=spec["n_eval"])
    return cc


def run_estimates(d, spec, obs, cc):
    for k, r in enumerate(spec["runs"], 1):
        kw = dict(kernel=r["kernel"], p=r["p"], seed=0)
        if "vce" in r:
            kw["vce"] = r["vce"]
        if "h" in r:
            kw["h"] = tuple(r["h"]) if isinstance(r["h"], list) else r["h"]
        if "bwselect" in r:
            kw["bwselect"] = r["bwselect"]
        if r.get("cluster"):
            kw["cluster"] = obs["cl"].to_numpy()
        if "conley" in r:
            kw["conley_cutoff"] = r["conley"]
        if r.get("fuzzy"):
            kw["fuzzy"] = obs["D"].to_numpy()
        if r.get("covs"):
            kw["covs"] = obs[["z1"]].to_numpy()
        if r.get("weights"):
            kw["weights"] = obs["w"].to_numpy()
        if "deriv" in r:
            kw["deriv"] = tuple(r["deriv"])
        res = rd2dnet.rdnet(obs["y_out"].to_numpy(), cc, **kw)
        res_frame(res.table).to_csv(d / f"py_res_{k}.csv", index=False)
    for k, r in enumerate(spec["density"], 1):
        dr = rd2dnet.density(cc, normalize=r["normalize"], bwselect=r["bwselect"])
        pd.DataFrame([{k2: getattr(dr, k2) for k2 in ["diff", "diff_bc", "se", "h_left", "h_right"]}]).to_csv(
            d / f"py_dens_{k}.csv", index=False)
    for k, r in enumerate(spec["bw"], 1):
        bw = rd2dnet.rdbwnet(obs["y_out"].to_numpy(), cc, bwselect=r["bwselect"])
        out = pd.DataFrame({"h0": bw["h0"], "h1": bw["h1"], "status": bw["status"].map(STATUS)})
        out.to_csv(d / f"py_bw_{k}.csv", index=False)


# --------------------------------------------------------------------------
# R and Stata scripts
# --------------------------------------------------------------------------

def r_script(spec, bnd):
    m = {"length": '"length"', "cost2": 'edge_cost("cost2")',
         "transition": f'transition_cost("length", class_attr = "cls", default_penalty = {spec["penalty"]!r}, '
                       'connectors = "c")'}[spec["metric"]]
    L = ["nodes <- read.csv(file.path(d, 'nodes.csv'))",
         "edges <- read.csv(file.path(d, 'edges.csv'), stringsAsFactors = FALSE)",
         "obs <- read.csv(file.path(d, 'obs.csv'))",
         "B <- as.matrix(read.csv(file.path(d, 'boundary.csv')))",
         f"net <- network_graph(edges, nodes, metric = {m}, directed = {str(spec['directed']).upper()}, "
         f"direction = '{spec['direction']}')"]
    if "region" in bnd:
        L.append("RG <- as.matrix(read.csv(file.path(d, 'region.csv')))")
        L.append("bnd <- boundary(B, treated_region = RG)")
    else:
        L.append(f"bnd <- boundary(B, treated_side = '{bnd['side']}')")
    L.append(f"cc <- network_coordinates(net, bnd, as.matrix(obs[, c('x', 'y')]), n_eval = {spec['n_eval']})")
    L.append("write.csv(data.frame(s = cc$s, treated = as.integer(cc$treated), arc = cc$crossing_arc), "
             "file.path(d, 'r_coords.csv'), row.names = FALSE)")
    for k, r in enumerate(spec["runs"], 1):
        a = [f"kernel = '{r['kernel']}'", f"p = {r['p']}", "seed = 0"]
        if "vce" in r:
            a.append(f"vce = '{r['vce']}'")
        if "h" in r:
            a.append(f"h = c({', '.join(map(repr, r['h']))})" if isinstance(r["h"], list) else f"h = {r['h']!r}")
        if "bwselect" in r:
            a.append(f"bwselect = '{r['bwselect']}'")
        if r.get("cluster"):
            a.append("cluster = obs$cl")
        if "conley" in r:
            a.append(f"conley_cutoff = {r['conley']!r}")
        if r.get("fuzzy"):
            a.append("fuzzy = obs$D")
        if r.get("covs"):
            a.append("covs = as.matrix(obs[, 'z1', drop = FALSE])")
        if r.get("weights"):
            a.append("weights = obs$w")
        if "deriv" in r:
            a.append(f"deriv = c({r['deriv'][0]}, {r['deriv'][1]})")
        L.append(f"res <- rdnet(obs$y_out, cc, {', '.join(a)})")
        L.append(f"write_res(res$table, file.path(d, 'r_res_{k}.csv'))")
    for k, r in enumerate(spec["density"], 1):
        L.append(f"dr <- density(cc, normalize = '{r['normalize']}', bwselect = '{r['bwselect']}')")
        L.append("write.csv(data.frame(diff = dr$diff, diff_bc = dr$diff_bc, se = dr$se, h_left = dr$h_left, "
                 f"h_right = dr$h_right), file.path(d, 'r_dens_{k}.csv'), row.names = FALSE)")
    for k, r in enumerate(spec["bw"], 1):
        L.append(f"bw <- rdbwnet(obs$y_out, cc, bwselect = '{r['bwselect']}')")
        L.append("write.csv(data.frame(h0 = bw$h0, h1 = bw$h1, status = status_code(bw$status)), "
                 f"file.path(d, 'r_bw_{k}.csv'), row.names = FALSE)")
    return "\n".join(x for x in L if x) + "\n"


def stata_script(spec, bnd, cdir):
    c = f'"{cdir}'
    L = [f"* parity case {spec['case']} ({spec['kind']})",
         f'rdn_par_load {c}"',
         ]
    opts = [f"nodes(nodes) edges(edges) boundary(BL) neval({spec['n_eval']})"]
    opts.append("treatedregion(RG)" if "region" in bnd else f"treatedside({bnd['side']})")
    if spec["directed"]:
        opts.append(f"directed direction({spec['direction']})")
    if spec["metric"] == "cost2":
        opts.append("cost(cost2)")
    if spec["metric"] == "transition":
        opts.append(f"classvar(cls) penalty({spec['penalty']!r}) connectors(c)")
    L.append(f"rdnetset x y, {' '.join(opts)}")
    L.append(f'rdn_par_coords {c}" {spec["case"]}')
    for k, r in enumerate(spec["runs"], 1):
        o = [f"kernel({r['kernel']})", f"p({r['p']})", "seed(0)"]
        w = ""
        if "vce" in r:
            o.append(f"vce({r['vce']})")
        if "h" in r:
            o.append(f"h({' '.join(map(repr, r['h']))})" if isinstance(r["h"], list) else f"h({r['h']!r})")
        if "bwselect" in r:
            o.append(f"bwselect({r['bwselect']})")
        if r.get("cluster"):
            o.append("vce(cluster cl)")
        if "conley" in r:
            o.append(f"vce(conley {r['conley']!r})")
        if r.get("fuzzy"):
            o.append("fuzzy(D)")
        if r.get("covs"):
            o.append("covs(z1)")
        if r.get("weights"):
            w = " [aw=w]"
        if "deriv" in r:
            o.append(f"deriv({r['deriv'][0]} {r['deriv'][1]})")
        L.append(f"rdnet y_out{w}, {' '.join(o)}")
        L.append(f'rdn_par_res {c}/py_res_{k}.csv" {spec["case"]} rdnet{k}')
    for k, r in enumerate(spec["density"], 1):
        L.append(f"rdnetdensity, normalize({r['normalize']}) bwselect({r['bwselect']})")
        L.append(f'rdn_par_dens {c}/py_dens_{k}.csv" {spec["case"]} density{k}')
    for k, r in enumerate(spec["bw"], 1):
        L.append(f"rdbwnet y_out, bwselect({r['bwselect']})")
        L.append(f'rdn_par_bw {c}/py_bw_{k}.csv" {spec["case"]} rdbwnet{k}')
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------

def make_case(out: Path, case: int, rng):
    kind = ["grid", "grid", "ring", "river"][case % 4]
    spec = random_spec(rng, kind, case)
    nodes, edges, bnd, box = make_design(rng, kind)
    E = edge_frame(rng, nodes, edges, spec["directed"])
    n = spec["n_obs"]
    if kind == "river":
        # observations along the rivers (monitoring stations)
        xy = nodes[["x", "y"]].to_numpy()
        pick = rng.integers(0, len(edges), n)
        t = rng.uniform(size=n)
        a = np.array([e[0] for e in edges])[pick]
        b = np.array([e[1] for e in edges])[pick]
        P = xy[a] * (1 - t)[:, None] + xy[b] * t[:, None] + rng.normal(0, 8, (n, 2))
    else:
        P = np.column_stack([rng.uniform(box[0], box[1], n), rng.uniform(box[2], box[3], n)])
    obs = pd.DataFrame({"x": P[:, 0], "y": P[:, 1]})
    d = out / f"{case:02d}"
    d.mkdir(parents=True, exist_ok=True)
    cc = run_python(d, spec, nodes, E, obs, bnd)
    T = np.asarray(cc.treated, dtype=float)
    s = np.asarray(cc.s, dtype=float)
    sm = np.nan_to_num(s)
    obs["z1"] = obs["x"] / 1000 + rng.normal(0, 0.3, n)
    obs["D"] = (rng.uniform(size=n) < 0.2 + 0.6 * T).astype(int)
    obs["y_out"] = (1 + 0.002 * sm + 0.0004 * obs["x"] + 0.5 * obs["z1"]
                    + T * (1 + 0.5 * np.sin(obs["y"] / 300)) + rng.normal(0, 1, n))
    obs["cl"] = (np.floor(obs["x"] / 200) * 10 + np.floor(obs["y"] / 200)).astype(int) + 100
    obs["w"] = rng.uniform(0.5, 2, n)
    obs["py_s"] = s
    obs["py_t"] = T
    obs["py_arc"] = np.asarray(cc.crossing_arc, dtype=float)
    nodes.to_csv(d / "nodes.csv", index=False, float_format="%.17g")
    E.to_csv(d / "edges.csv", index=False, float_format="%.17g")
    obs.to_csv(d / "obs.csv", index=False, float_format="%.17g")
    pd.DataFrame(bnd["line"], columns=["x", "y"]).to_csv(d / "boundary.csv", index=False, float_format="%.17g")
    if "region" in bnd:
        pd.DataFrame(bnd["region"], columns=["x", "y"]).to_csv(d / "region.csv", index=False, float_format="%.17g")
    # reread so that Python sees exactly the numbers written to the files
    nodes = pd.read_csv(d / "nodes.csv")
    E = pd.read_csv(d / "edges.csv")
    obs = pd.read_csv(d / "obs.csv")
    bnd = {**bnd, "line": pd.read_csv(d / "boundary.csv").to_numpy()}
    if "region" in bnd:
        bnd["region"] = pd.read_csv(d / "region.csv").to_numpy()
    cc = run_python(d, spec, nodes, E, obs, bnd)
    pd.DataFrame({"s": cc.s, "treated": np.asarray(cc.treated, dtype=int), "arc": cc.crossing_arc}).to_csv(
        d / "py_coords.csv", index=False, float_format="%.17g")
    run_estimates(d, spec, obs, cc)
    (d / "spec.json").write_text(json.dumps(spec, indent=1))
    (d / "run.R").write_text(r_script(spec, bnd))
    (d / "run.do").write_text(stata_script(spec, bnd, f"parity/cases/{case:02d}"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=str(Path(__file__).parent / "cases"))
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rng = np.random.default_rng(a.seed)
    for case in range(1, a.n + 1):
        make_case(out, case, rng)
        print(f"case {case:02d} written")
    (out / "cases.txt").write_text("\n".join(f"{c:02d}" for c in range(1, a.n + 1)) + "\n")


if __name__ == "__main__":
    main()
