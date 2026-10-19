"""Compare R (or another language's) parity results with the Python reference.

    python parity/compare.py [parity/cases] [--prefix r] [--tol 1e-8]

Exits with status 1 when any number differs by more than ``tol`` (relative to
max(1, |reference|)) or any status code differs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def rel(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    both_nan = np.isnan(a) & np.isnan(b)
    d = np.abs(a - b) / np.maximum(1.0, np.abs(b))
    d[both_nan] = 0.0
    d[np.isnan(d)] = np.inf
    return d


def compare_case(d: Path, prefix: str):
    rows = []
    for ref in sorted(d.glob("py_*.csv")):
        other = d / ref.name.replace("py_", f"{prefix}_", 1)
        if not other.exists():
            rows.append((ref.stem, "missing", np.inf))
            continue
        P = pd.read_csv(ref)
        O = pd.read_csv(other)
        if len(P) != len(O):
            rows.append((ref.stem, "rows", np.inf))
            continue
        if "status" in P:
            bad = (P["status"].fillna(-1).to_numpy() != O["status"].fillna(-1).to_numpy())
            if bad.any():
                rows.append((ref.stem, f"status rows {list(np.where(bad)[0] + 1)}", np.inf))
            ok = (P["status"] == 0).to_numpy()
        else:
            ok = np.ones(len(P), bool)
        for c in P.columns:
            if c == "status":
                continue
            if c not in O:
                rows.append((ref.stem, c + " missing", np.inf))
                continue
            dd = rel(O[c][ok], P[c][ok])
            rows.append((ref.stem, c, float(dd.max()) if len(dd) else 0.0))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default=str(Path(__file__).parent / "cases"))
    ap.add_argument("--prefix", default="r")
    ap.add_argument("--tol", type=float, default=1e-8)
    a = ap.parse_args()
    root = Path(a.root)
    worst = 0.0
    bad = []
    for cs in (root / "cases.txt").read_text().split():
        rows = compare_case(root / cs, a.prefix)
        m = max((r[2] for r in rows), default=0.0)
        worst = max(worst, m)
        flag = "ok" if m <= a.tol else "FAIL"
        print(f"case {cs}: max rel diff {m:.2e}  {flag}")
        for name, col, v in rows:
            if v > a.tol:
                print(f"    {name} {col}: {v:.3e}")
                bad.append((cs, name, col, v))
    print(f"worst relative difference: {worst:.2e}")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
