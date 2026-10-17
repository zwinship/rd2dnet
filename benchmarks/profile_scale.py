"""Time each rd2dnet stage on a large synthetic network.

    python benchmarks/profile_scale.py --net grid_city --n 100000 --design parcels
    python benchmarks/profile_scale.py --net river --n 1000000 --design panel --stages coords,rdnet_hc1
    python benchmarks/profile_scale.py --net grid_city --n 100000 --profile network_coordinates

Prints one line per stage (seconds) and, with ``--profile STAGE``, the top
functions by cumulative time for that stage.
"""

from __future__ import annotations

import argparse
import cProfile
import pstats
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import synthetic  # noqa: E402

import rd2dnet  # noqa: E402

ALL = ["graph", "coords", "rdbwnet", "rdnet_hc1", "rdnet_cluster", "rdnet_conley", "density",
       "transition", "travel_time"]
# Typical bandwidth / HAC cutoff / cluster cell per network family (metres).
SCALE = {"grid_city": (400.0, 150.0, 300.0), "suburb": (600.0, 200.0, 400.0), "river": (3000.0, 1000.0, 2000.0)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="grid_city", choices=["grid_city", "suburb", "river"])
    ap.add_argument("--size", type=int, default=50_000, help="approximate number of network nodes")
    ap.add_argument("--n", type=int, default=100_000)
    ap.add_argument("--design", default="parcels", choices=["parcels", "gps", "panel"])
    ap.add_argument("--sites", type=int, default=2000)
    ap.add_argument("--n-eval", type=int, default=25)
    ap.add_argument("--stages", default=",".join(ALL))
    ap.add_argument("--profile", default=None)
    a = ap.parse_args(argv)
    stages = a.stages.split(",")
    h, cutoff, cell = SCALE[a.net]

    t0 = time.perf_counter()
    if a.net == "grid_city":
        ns = int(np.sqrt(a.size * 1.25))
        net = synthetic.grid_city(n_streets=ns, n_avenues=int(a.size / ns))
    elif a.net == "suburb":
        net = synthetic.suburb(n_nodes=a.size)
    else:
        net = synthetic.river(n_nodes=a.size)
    xy, y = synthetic.observations(net, a.n, a.design, n_sites=a.sites)
    print(f"# {a.net}: {net.graph.number_of_nodes()} nodes, {net.graph.number_of_edges()} edges; "
          f"n={len(y)} ({a.design}); generated in {time.perf_counter() - t0:.1f}s", flush=True)

    out = {}

    def stage(name, fn):
        if name not in stages:
            return None
        prof = cProfile.Profile() if a.profile == name else None
        t = time.perf_counter()
        if prof:
            prof.enable()
        r = fn()
        if prof:
            prof.disable()
        out[name] = time.perf_counter() - t
        print(f"{name:14s} {out[name]:8.2f} s", flush=True)
        if prof:
            pstats.Stats(prof).sort_stats("cumulative").print_stats(25)
        return r

    G = net.graph
    ng = stage("graph", lambda: rd2dnet.NetworkGraph(G, "length")) or rd2dnet.NetworkGraph(G, "length")
    coords = stage("coords", lambda: rd2dnet.network_coordinates(ng, net.boundary, xy, n_eval=a.n_eval))
    if coords is None:
        coords = rd2dnet.network_coordinates(ng, net.boundary, xy, n_eval=a.n_eval)
    ok = np.isfinite(coords.s)
    print(f"# reachable {ok.mean():.3f}, treated {coords.treated.mean():.3f}, masked eval {coords.eval_masked.sum()}")
    bw = stage("rdbwnet", lambda: rd2dnet.rdbwnet(y, coords))
    if bw is not None:
        good = bw["status"] == "ok"
        print(f"# bandwidths ok {int(good.sum())}/{len(bw)}, median h {np.nanmedian(bw['h0']):.0f}")
        h = bw[["h0", "h1"]].to_numpy()
    r = stage("rdnet_hc1", lambda: rd2dnet.rdnet(y, coords, h=h, seed=0))
    if r is not None:
        t = r.table[r.table["status"] == "ok"]
        if len(t):
                print(f"# estimated {len(t)}/{len(r.table)}, median n0+n1 {np.median(t['n0'] + t['n1']):.0f}, "
                  f"WBATE {r.wbate((r.table['status'] == 'ok').astype(float))['tau_bc']:.3f}")
    cl = (np.floor(xy[:, 0] / cell) * 1e6 + np.floor(xy[:, 1] / cell)).astype(np.int64)
    stage("rdnet_cluster", lambda: rd2dnet.rdnet(y, coords, h=h, cluster=cl, seed=0))
    stage("rdnet_conley", lambda: rd2dnet.rdnet(y, coords, h=h, conley_cutoff=cutoff, seed=0))
    stage("density", lambda: rd2dnet.density(coords))
    if a.net == "grid_city":
        tc = rd2dnet.TransitionCost("length", class_attr="highway", default_penalty=500.0,
                                    connectors={"motorway_link"})
        stage("transition", lambda: rd2dnet.NetworkGraph(G, tc).engine)
    stage("travel_time", lambda: rd2dnet.network_coordinates(rd2dnet.NetworkGraph(G, rd2dnet.travel_time()),
                                                              net.boundary, xy, n_eval=a.n_eval))
    return out


if __name__ == "__main__":
    main()
