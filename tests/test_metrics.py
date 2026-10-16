import networkx as nx
import numpy as np
import pytest

import rd2dnet
from rd2dnet.metrics import EdgeCost, Precomputed, TransitionCost, as_metric, node_contrast, travel_time


def test_edge_length_is_sum_of_lengths(highway_graph):
    net = rd2dnet.NetworkGraph(highway_graph, "length")
    assert net.node_distance("a0", "b0") == pytest.approx(400.0)
    assert net.node_distance("a0", "l1") == pytest.approx(300.0)


def test_callable_edge_cost(highway_graph):
    net = rd2dnet.NetworkGraph(highway_graph, lambda u, v, d: 2 * d["length"])
    assert net.node_distance("a0", "b0") == pytest.approx(800.0)


def test_travel_time_metric():
    G = nx.Graph()
    G.add_node(0, x=0, y=0)
    G.add_node(1, x=1000, y=0)
    G.add_edge(0, 1, length=1000.0, speed_kph=36.0)
    net = rd2dnet.NetworkGraph(G, travel_time())
    assert net.node_distance(0, 1) == pytest.approx(100.0)  # 1 km at 10 m/s


def test_transition_cost_highway_ramp_example(highway_graph):
    """Highway -> ramp -> highway stays close; highway -> ramp -> local is far."""
    metric = TransitionCost(
        "length",
        class_attr="highway",
        default_penalty=5000.0,
        connectors={"motorway_link"},
    )
    net = rd2dnet.NetworkGraph(highway_graph, metric)
    assert net.node_distance("a0", "b0") == pytest.approx(400.0)  # no class change
    assert net.node_distance("a0", "l1") == pytest.approx(300.0 + 5000.0)  # motorway -> residential


def test_transition_penalty_dict_and_asymmetry(highway_graph):
    metric = TransitionCost(
        "length",
        penalty={("motorway", "residential"): 1000.0, ("residential", "motorway"): 10.0},
        connectors={"motorway_link"},
    )
    net = rd2dnet.NetworkGraph(highway_graph, metric)
    assert net.node_distance("a0", "l1") == pytest.approx(1300.0)
    assert net.node_distance("l1", "a0") == pytest.approx(310.0)


def test_transition_without_connectors_matches_edge_cost(highway_graph):
    net_t = rd2dnet.NetworkGraph(highway_graph, TransitionCost("length"))
    net_e = rd2dnet.NetworkGraph(highway_graph, "length")
    for a in highway_graph:
        for b in highway_graph:
            assert net_t.node_distance(a, b) == pytest.approx(net_e.node_distance(a, b))


def test_nearest_source_matches_bruteforce(highway_graph):
    metric = TransitionCost("length", default_penalty=50.0, connectors={"motorway_link"})
    net = rd2dnet.NetworkGraph(highway_graph, metric)
    sources = [net.node_index["a0"], net.node_index["l1"]]
    dist, src = net.nearest_source(sources)
    full = net.distances_from(sources)
    np.testing.assert_allclose(dist, full.min(axis=0))
    np.testing.assert_array_equal(src, np.array(sources)[full.argmin(axis=0)])


def test_node_contrast_metric():
    G = nx.Graph()
    G.add_node(0, x=0, y=0, landuse=0.0)
    G.add_node(1, x=100, y=0, landuse=3.0)
    G.add_edge(0, 1, length=100.0)
    net = rd2dnet.NetworkGraph(G, node_contrast("landuse", weight=10.0))
    assert net.node_distance(0, 1) == pytest.approx(130.0)


def test_precomputed_metric(highway_graph):
    nodes = list(highway_graph)
    D = np.abs(np.subtract.outer(np.arange(len(nodes)), np.arange(len(nodes)))).astype(float)
    net = rd2dnet.NetworkGraph(highway_graph, Precomputed(D, nodes=nodes))
    assert net.node_distance(nodes[0], nodes[3]) == pytest.approx(3.0)


def test_negative_cost_rejected(highway_graph):
    net = rd2dnet.NetworkGraph(highway_graph, lambda u, v, d: -1.0)
    with pytest.raises(ValueError, match="invalid cost"):
        net.engine


def test_infinite_cost_blocks_edge(highway_graph):
    net = rd2dnet.NetworkGraph(highway_graph, lambda u, v, d: np.inf if d["highway"] == "residential" else d["length"])
    assert np.isinf(net.node_distance("a0", "l1"))


def test_zero_cost_edges_are_kept():
    G = nx.path_graph(3)
    for n in G:
        G.nodes[n].update(x=float(n), y=0.0)
    nx.set_edge_attributes(G, 0.0, "length")
    net = rd2dnet.NetworkGraph(G, "length")
    assert net.node_distance(0, 2) == pytest.approx(0.0, abs=1e-9)


def test_directed_graph_direction():
    G = nx.DiGraph()
    G.add_node(0, x=0, y=0)
    G.add_node(1, x=1, y=0)
    G.add_edge(0, 1, length=5.0)
    assert np.isinf(rd2dnet.NetworkGraph(G, direction="out").node_distance(1, 0))
    assert rd2dnet.NetworkGraph(G, direction="in").node_distance(1, 0) == pytest.approx(5.0)
    assert rd2dnet.NetworkGraph(G, direction="both").node_distance(1, 0) == pytest.approx(5.0)


def test_as_metric_coercion():
    assert isinstance(as_metric("length"), EdgeCost)
    assert isinstance(as_metric(np.zeros((2, 2))), Precomputed)
    with pytest.raises(TypeError):
        as_metric(3)


@pytest.mark.parametrize("kind", ["edge", "transition", "precomputed"])
def test_pairs_within_matches_full_dijkstra(city, kind):
    G = city.graph
    if kind == "edge":
        net = rd2dnet.NetworkGraph(G, "length")
    elif kind == "transition":
        net = rd2dnet.NetworkGraph(G, rd2dnet.TransitionCost("length", class_attr="highway", default_penalty=120.0,
                                                            connectors={"motorway_link"}))
    else:
        D = rd2dnet.NetworkGraph(G, "length").distances_from(np.arange(G.number_of_nodes()))
        net = rd2dnet.NetworkGraph(G, rd2dnet.Precomputed(D))
    src = np.arange(0, net.n_nodes, 7)
    full = net.distances_from(src)
    i, node, d = net.engine.pairs_within(src, 400.0)
    ref = np.argwhere(full < 400.0)
    np.testing.assert_array_equal(np.lexsort((node, i)).size, len(ref))
    got = np.zeros_like(full, dtype=bool)
    got[i, node] = True
    np.testing.assert_array_equal(got, full < 400.0)
    np.testing.assert_allclose(d, full[i, node])
