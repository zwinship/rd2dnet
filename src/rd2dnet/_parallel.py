"""Thread-parallel evaluation of vectorised geometry operations on large inputs.

Shapely (GEOS) releases the GIL in its vectorised functions and STRtree
queries, so splitting a large array of points across threads scales with
the number of cores. Set ``RD2DNET_NUM_THREADS`` to cap the thread count
(``1`` disables threading).
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

MIN_CHUNK = 25_000


def n_threads(n: int) -> int:
    env = os.environ.get("RD2DNET_NUM_THREADS")
    cap = int(env) if env else min(os.cpu_count() or 1, 16)
    return max(1, min(cap, n // MIN_CHUNK))


def chunked(fn, n: int, *arrays):
    """``fn(*arrays)`` evaluated on row chunks in threads; outputs are concatenated.

    ``fn`` must act row-wise and return an array, or a tuple of arrays, with
    one leading row per input row.
    """
    k = n_threads(n)
    if k <= 1:
        return fn(*arrays)
    edges = np.linspace(0, n, 4 * k + 1).astype(int)
    parts = [slice(a, b) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    run = lambda sl: fn(*(a[sl] for a in arrays))  # noqa: E731
    # The first chunk runs on this thread so that structures GEOS builds lazily
    # on first use (STRtree, prepared geometries) exist before threads share them.
    out = [run(parts[0])]
    with ThreadPoolExecutor(k) as ex:
        out += list(ex.map(run, parts[1:]))
    if isinstance(out[0], tuple):
        return tuple(np.concatenate(z) for z in zip(*out))
    return np.concatenate(out)
