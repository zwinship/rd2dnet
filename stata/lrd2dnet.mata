*! lrd2dnet.mata  version 0.2.0
*! Mata library for rd2dnet: boundary discontinuity designs on networks.
*!
*! A port of the Python/R package rd2dnet (https://github.com/zwinship/rd2dnet),
*! which extends rd2d (Cattaneo, Titiunik and Yu) to network coordinates. The
*! bandwidth selector is ported from rd2d and the density test from rddensity
*! (Cattaneo, Jansson and Ma); both are GPL-3.
*!
*! This file is compiled at run time by rd2dnet_mata.ado (no .mlib is needed).
*! Conventions: indices are 1-based; "infinite" distances are 1e300 inside the
*! library and are turned into missing values before results reach Stata.

version 16
mata:
mata set matastrict off

real scalar rdn_mata_version()
{
    return(1)
}

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

real scalar rdn_inf()
{
    return(1e300)
}

// Elementwise c ? a : b (a and b may be scalars).
real matrix rdn_where(real matrix c, real matrix a, real matrix b)
{
    real matrix out
    real scalar i, j, r, k, sa, sb
    r = rows(c)
    k = cols(c)
    sa = (rows(a) == 1 & cols(a) == 1)
    sb = (rows(b) == 1 & cols(b) == 1)
    out = J(r, k, .)
    for (i = 1; i <= r; i++) {
        for (j = 1; j <= k; j++) {
            if (c[i, j] != 0 & c[i, j] < .) {
                out[i, j] = (sa ? a[1, 1] : a[i, j])
            }
            else {
                out[i, j] = (sb ? b[1, 1] : b[i, j])
            }
        }
    }
    return(out)
}

// Elementwise minimum of two column vectors (either may be a scalar).
real colvector rdn_pmin(real colvector a, real colvector b)
{
    real colvector aa, bb
    real scalar n
    n = max((rows(a), rows(b)))
    aa = (rows(a) == 1 ? J(n, 1, a[1]) : a)
    bb = (rows(b) == 1 ? J(n, 1, b[1]) : b)
    if (n == 0) return(J(0, 1, .))
    return(rowmin((aa, bb)))
}

// Large values (>= 1e299) become missing (column vectors).
real colvector rdn_inf2miss(real colvector x)
{
    real colvector out
    real scalar i
    out = x
    for (i = 1; i <= rows(x); i++) {
        if (x[i] < . & abs(x[i]) >= 1e299) out[i] = .
    }
    return(out)
}

// Finite (non-missing) indicator.
real colvector rdn_fin(real colvector x)
{
    return(x :< .)
}

// First index of the minimum.
real scalar rdn_argmin(real vector v)
{
    real scalar i, k, best
    k = 1
    best = v[1]
    for (i = 2; i <= length(v); i++) {
        if (v[i] < best) {
            best = v[i]
            k = i
        }
    }
    return(k)
}

// k-th smallest value (k clamped to the data).
real scalar rdn_kth(real colvector d, real scalar k)
{
    real colvector s
    real scalar kk
    s = sort(d, 1)
    kk = min((max((k, 1)), rows(s)))
    return(s[kk])
}

// Positions of the first entry for each distinct value of x.
real colvector rdn_first_per_value(real colvector x)
{
    real colvector o, keep
    real scalar i
    if (rows(x) == 0) return(J(0, 1, .))
    o = order((x, (1::rows(x))), (1, 2))
    keep = J(rows(x), 1, 0)
    keep[1] = 1
    for (i = 2; i <= rows(x); i++) keep[i] = (x[o[i]] != x[o[i - 1]])
    return(sort(o[selectindex(keep)], 1))
}

// Quantile with linear interpolation (numpy default, R type 7).
real scalar rdn_quantile(real colvector x, real scalar p)
{
    real colvector s
    real scalar n, h, lo, hi
    s = sort(x, 1)
    n = rows(s)
    h = (n - 1) * p + 1
    lo = floor(h)
    hi = ceil(h)
    return(s[lo] + (h - lo) * (s[hi] - s[lo]))
}

real scalar rdn_median(real colvector x)
{
    return(rdn_quantile(x, 0.5))
}

// numpy.interp on a non-decreasing grid that may contain ties.
real colvector rdn_interp(real colvector x, real colvector grid, real colvector vals)
{
    real colvector out, o
    real scalar n, m, j, k, i, xi, t
    n = rows(grid)
    m = rows(x)
    out = J(m, 1, .)
    o = order(x, 1)
    j = 0
    for (k = 1; k <= m; k++) {
        i = o[k]
        xi = x[i]
        if (xi >= .) continue
        while (j < n) {
            if (grid[j + 1] <= xi) j++
            else break
        }
        if (j == 0) out[i] = vals[1]
        else if (j >= n) out[i] = vals[n]
        else {
            t = (xi - grid[j]) / (grid[j + 1] - grid[j])
            out[i] = vals[j] + t * (vals[j + 1] - vals[j])
        }
    }
    return(out)
}

// Moore-Penrose inverse with numpy's relative cutoff.
real matrix rdn_pinv(real matrix A, real scalar rcond)
{
    real matrix U, Vt
    real colvector s
    real colvector keep
    if (hasmissing(A)) return(J(cols(A), rows(A), .))
    svd(A, U, s, Vt)
    s = vec(s)
    keep = selectindex(s :> rcond * max(s))
    if (rows(keep) == 0) return(J(cols(A), rows(A), 0))
    return(Vt'[., keep] * diag(1 :/ s[keep]) * U[., keep]')
}

// Minimum-norm least squares (numpy.linalg.lstsq with rcond=None).
real colvector rdn_lstsq(real matrix A, real colvector b)
{
    return(rdn_pinv(A, epsilon(1) * max((rows(A), cols(A)))) * b)
}

// Symmetric inverse after a conditioning check; missing matrix when singular.
real matrix rdn_inv_checked(real matrix G)
{
    real matrix Gi
    if (hasmissing(G)) return(J(rows(G), cols(G), .))
    if (!(cond(G) <= 1e12)) return(J(rows(G), cols(G), .))
    Gi = luinv(G)
    return((Gi + Gi') / 2)
}

// Inverse of a local Gram matrix through its Cholesky factor (rd2d's way).
real matrix rdn_inv_gram(real matrix A)
{
    // bases are scaled by the bandwidth, so well-posed designs have condition
    // numbers far below 1e8; above it the plug-in constants are unreliable
    if (hasmissing(A)) return(J(rows(A), cols(A), .))
    if (!(cond(A) <= 1e8)) return(J(rows(A), cols(A), .))
    return(cholinv(A))
}

// Row sums of X by integer group codes; rows for groups that occur only.
real matrix rdn_groupsum(real matrix X, real colvector g)
{
    real colvector u, pos
    real matrix out
    real scalar i
    u = uniqrows(g)
    pos = J(max(u), 1, 0)
    pos[u] = (1::rows(u))
    out = J(rows(u), cols(X), 0)
    for (i = 1; i <= rows(X); i++) {
        out[pos[g[i]], .] = out[pos[g[i]], .] + X[i, .]
    }
    return(out)
}

real scalar rdn_nunique(real colvector g)
{
    if (rows(g) == 0) return(0)
    return(rows(uniqrows(g)))
}

// numpy.isclose with default tolerances (missing -> 0).
real colvector rdn_isclose(real colvector a, real colvector b)
{
    real colvector out
    out = (abs(a - b) :<= 1e-8 :+ 1e-5 :* abs(b))
    out = out :* (a :< .) :* (b :< .)
    return(out)
}

// ---------------------------------------------------------------------------
// Kernels: 1 triangular, 2 epanechnikov, 3 uniform, 4 gaussian
// ---------------------------------------------------------------------------

real scalar rdn_kernel_code(string scalar k)
{
    string scalar s
    s = strlower(strtrim(k))
    if (s == "" | s == "tri" | s == "triangular") return(1)
    if (s == "epa" | s == "epanechnikov") return(2)
    if (s == "uni" | s == "uniform") return(3)
    if (s == "gau" | s == "gaussian") return(4)
    _error(198, "unknown kernel; choose triangular, epanechnikov, uniform or gaussian")
    return(.)
}

string scalar rdn_kernel_name(real scalar k)
{
    if (k == 1) return("triangular")
    if (k == 2) return("epanechnikov")
    if (k == 3) return("uniform")
    return("gaussian")
}

real colvector rdn_kernel(real colvector u, real scalar k)
{
    real colvector a, out
    a = editmissing(abs(u), 2)
    if (k == 4) {
        a = editmissing(abs(u), 1e300)
        out = editmissing(exp(-0.5 :* (a :^ 2)) :/ sqrt(2 * pi()), 0)
        return(out :* (a :< 1e299))
    }
    if (k == 1) out = (1 :- a) :* (a :<= 1)
    else if (k == 2) out = 0.75 :* (1 :- a :^ 2) :* (a :<= 1)
    else out = 0.5 :* (a :<= 1)
    return(out)
}

// rd2d's rule-of-thumb constants (mu_2(K)^2, R(K)).
real rowvector rdn_rot_const(real scalar k)
{
    if (k == 2) return((1 / 6, 4 / (3 * pi())))
    if (k == 3) return((1 / 4, 1 / pi()))
    if (k == 4) return((1, 1 / (4 * pi())))
    return((3 / 20, 3 / (2 * pi())))
}

// ---------------------------------------------------------------------------
// Planar geometry
// ---------------------------------------------------------------------------

struct rdn_segs {
    real colvector x0, y0, x1, y1, len, cum
    real scalar total
}

// Segments of a polyline (k x 2 vertex matrix).
struct rdn_segs scalar rdn_segments(real matrix L)
{
    struct rdn_segs scalar S
    real scalar k
    k = rows(L)
    S.x0 = L[(1::(k - 1)), 1]
    S.y0 = L[(1::(k - 1)), 2]
    S.x1 = L[(2::k), 1]
    S.y1 = L[(2::k), 2]
    S.len = sqrt((S.x1 - S.x0) :^ 2 + (S.y1 - S.y0) :^ 2)
    S.cum = runningsum(S.len) - S.len
    S.total = sum(S.len)
    return(S)
}

struct rdn_near {
    real colvector seg, dist, t
}

real scalar rdn_gridcell(real scalar x, real scalar lo, real scalar cs, real scalar n)
{
    real scalar c
    c = floor((x - lo) / cs) + 1
    if (c < 1) c = 1
    if (c > n) c = n
    return(c)
}

// Nearest segment to each point (ties go to the lowest segment index), the
// distance, and the position t in [0, 1] of the foot along the segment.
struct rdn_near scalar rdn_nearest_segment(real colvector px, real colvector py,
    real colvector x0, real colvector y0, real colvector x1, real colvector y1)
{
    struct rdn_near scalar R
    real scalar m, n, xmin, xmax, ymin, ymax, W, H, cs, nx, ny, ncell
    real scalar s, a1, b1, a2, b2, cx, cy, c, i, p, q, best, bs, bt, r, lb
    real scalar ix, iy, k, dx, dy, L2, t, qx, qy, d, maxr
    real colvector cnt, ptr, fill, items, stamp
    m = rows(px)
    n = rows(x0)
    R.seg = J(m, 1, .)
    R.dist = J(m, 1, .)
    R.t = J(m, 1, .)
    if (n == 0 | m == 0) return(R)
    xmin = min((x0 \ x1))
    xmax = max((x0 \ x1))
    ymin = min((y0 \ y1))
    ymax = max((y0 \ y1))
    W = xmax - xmin
    H = ymax - ymin
    if (W * H > 0) cs = sqrt(W * H / n)
    else cs = max((W, H)) / n
    if (!(cs > 0 & cs < .)) cs = 1
    nx = min((2048, ceil(W / cs) + 1))
    ny = min((2048, ceil(H / cs) + 1))
    cs = max((cs, W / (nx - 0.5), H / (ny - 0.5)))
    ncell = nx * ny
    cnt = J(ncell, 1, 0)
    for (s = 1; s <= n; s++) {
        a1 = rdn_gridcell(min((x0[s], x1[s])), xmin, cs, nx)
        b1 = rdn_gridcell(max((x0[s], x1[s])), xmin, cs, nx)
        a2 = rdn_gridcell(min((y0[s], y1[s])), ymin, cs, ny)
        b2 = rdn_gridcell(max((y0[s], y1[s])), ymin, cs, ny)
        for (cx = a1; cx <= b1; cx++) {
            for (cy = a2; cy <= b2; cy++) {
                c = (cy - 1) * nx + cx
                cnt[c] = cnt[c] + 1
            }
        }
    }
    ptr = 1 \ (runningsum(cnt) :+ 1)
    fill = ptr[(1::ncell)]
    items = J(ptr[ncell + 1] - 1, 1, 0)
    for (s = 1; s <= n; s++) {
        a1 = rdn_gridcell(min((x0[s], x1[s])), xmin, cs, nx)
        b1 = rdn_gridcell(max((x0[s], x1[s])), xmin, cs, nx)
        a2 = rdn_gridcell(min((y0[s], y1[s])), ymin, cs, ny)
        b2 = rdn_gridcell(max((y0[s], y1[s])), ymin, cs, ny)
        for (cx = a1; cx <= b1; cx++) {
            for (cy = a2; cy <= b2; cy++) {
                c = (cy - 1) * nx + cx
                items[fill[c]] = s
                fill[c] = fill[c] + 1
            }
        }
    }
    stamp = J(n, 1, 0)
    maxr = max((nx, ny))
    for (i = 1; i <= m; i++) {
        p = px[i]
        q = py[i]
        if (p >= . | q >= .) continue
        cx = rdn_gridcell(p, xmin, cs, nx)
        cy = rdn_gridcell(q, ymin, cs, ny)
        best = rdn_inf()
        bs = 0
        bt = 0
        for (r = 0; r <= maxr; r++) {
            lb = (r == 0 ? 0 : (r - 1) * cs)
            if (lb > best) break
            for (ix = cx - r; ix <= cx + r; ix++) {
                if (ix < 1 | ix > nx) continue
                for (iy = cy - r; iy <= cy + r; iy++) {
                    if (iy < 1 | iy > ny) continue
                    if (abs(ix - cx) != r & abs(iy - cy) != r) continue
                    c = (iy - 1) * nx + ix
                    for (k = ptr[c]; k < ptr[c + 1]; k++) {
                        s = items[k]
                        if (stamp[s] == i) continue
                        stamp[s] = i
                        dx = x1[s] - x0[s]
                        dy = y1[s] - y0[s]
                        L2 = dx * dx + dy * dy
                        t = 0
                        if (L2 > 0) {
                            t = ((p - x0[s]) * dx + (q - y0[s]) * dy) / L2
                            if (t < 0) t = 0
                            if (t > 1) t = 1
                        }
                        qx = x0[s] + t * dx - p
                        qy = y0[s] + t * dy - q
                        d = sqrt(qx * qx + qy * qy)
                        if (d < best | (d == best & s < bs)) {
                            best = d
                            bs = s
                            bt = t
                        }
                    }
                }
            }
        }
        R.seg[i] = bs
        R.dist[i] = best
        R.t[i] = bt
    }
    return(R)
}

// Ring matrices: k x 2 (one ring) or k x 3 (x, y, ring id).
real matrix rdn_rings(real matrix R)
{
    if (cols(R) == 2) return((R, J(rows(R), 1, 1)))
    return(R[., (1, 2, 3)])
}

// Even-odd point-in-polygon over all rings.
real colvector rdn_in_rings(real colvector px, real colvector py, real matrix Rg)
{
    real colvector out
    real scalar m, k, i, x, y, inside, st, en, a, b, ya, yb, xc
    m = rows(px)
    k = rows(Rg)
    out = J(m, 1, 0)
    for (i = 1; i <= m; i++) {
        x = px[i]
        y = py[i]
        if (x >= . | y >= .) continue
        inside = 0
        st = 1
        while (st <= k) {
            en = st
            while (en < k) {
                if (Rg[en + 1, 3] == Rg[st, 3]) en++
                else break
            }
            b = en
            for (a = st; a <= en; a++) {
                ya = Rg[a, 2]
                yb = Rg[b, 2]
                if ((ya > y) != (yb > y)) {
                    xc = Rg[a, 1] + (y - ya) * (Rg[b, 1] - Rg[a, 1]) / (yb - ya)
                    if (x < xc) inside = !inside
                }
                b = a
            }
            st = en + 1
        }
        out[i] = inside
    }
    return(out)
}

// Closed-ring segments of all rings, as a 4-column matrix (x0, y0, x1, y1).
real matrix rdn_ring_segments(real matrix Rg)
{
    real matrix out
    real scalar k, st, en, a, b
    k = rows(Rg)
    out = J(0, 4, .)
    st = 1
    while (st <= k) {
        en = st
        while (en < k) {
            if (Rg[en + 1, 3] == Rg[st, 3]) en++
            else break
        }
        b = en
        for (a = st; a <= en; a++) {
            out = out \ (Rg[b, 1], Rg[b, 2], Rg[a, 1], Rg[a, 2])
            b = a
        }
        st = en + 1
    }
    return(out)
}

// Distance from points to a polygon (0 inside).
real colvector rdn_polygon_distance(real matrix Rg, real matrix xy)
{
    struct rdn_near scalar R
    real matrix S
    real colvector d
    S = rdn_ring_segments(Rg)
    R = rdn_nearest_segment(xy[., 1], xy[., 2], S[., 1], S[., 2], S[., 3], S[., 4])
    d = R.dist
    d = rdn_where(rdn_in_rings(xy[., 1], xy[., 2], Rg), 0, d)
    return(d)
}

// ---------------------------------------------------------------------------
// Boundary
// ---------------------------------------------------------------------------

struct rdn_bnd {
    real matrix line, region, maskarcs, maskpoly
    struct rdn_segs scalar segs
    real scalar blen, tol, hasregion
    string scalar side
}

struct rdn_bnd scalar rdn_boundary(real matrix line, real matrix region, string scalar side,
    real matrix maskarcs, real matrix maskpoly, real scalar tol)
{
    struct rdn_bnd scalar B
    if (rows(line) < 2 | cols(line) < 2) _error(198, "the boundary needs at least two vertices (a k x 2 matrix)")
    B.line = line[., (1, 2)]
    B.segs = rdn_segments(B.line)
    B.blen = B.segs.total
    B.tol = tol
    B.hasregion = (rows(region) > 0)
    if (B.hasregion) B.region = rdn_rings(region)
    else B.region = J(0, 3, .)
    B.side = side
    if (!B.hasregion & side != "left" & side != "right") {
        _error(198, "give exactly one of treatedregion() or treatedside(left|right)")
    }
    if (rows(maskarcs)) B.maskarcs = (rowmin(maskarcs[., (1, 2)]), rowmax(maskarcs[., (1, 2)]))
    else B.maskarcs = J(0, 2, .)
    if (rows(maskpoly)) B.maskpoly = rdn_rings(maskpoly)
    else B.maskpoly = J(0, 3, .)
    return(B)
}

// Whether the boundary is a closed ring (first vertex = last vertex).
real scalar rdn_bclosed(struct rdn_bnd scalar B)
{
    real scalar k
    k = rows(B.line)
    return(k > 2 & B.line[1, 1] == B.line[k, 1] & B.line[1, 2] == B.line[k, 2])
}

real colvector rdn_bdistance(struct rdn_bnd scalar B, real matrix xy)
{
    struct rdn_near scalar R
    R = rdn_nearest_segment(xy[., 1], xy[., 2], B.segs.x0, B.segs.y0, B.segs.x1, B.segs.y1)
    return(R.dist)
}

real colvector rdn_bproject(struct rdn_bnd scalar B, real matrix xy)
{
    struct rdn_near scalar R
    real colvector out
    real scalar i
    R = rdn_nearest_segment(xy[., 1], xy[., 2], B.segs.x0, B.segs.y0, B.segs.x1, B.segs.y1)
    out = J(rows(xy), 1, .)
    for (i = 1; i <= rows(xy); i++) {
        if (R.seg[i] < .) out[i] = B.segs.cum[R.seg[i]] + R.t[i] * B.segs.len[R.seg[i]]
    }
    return(out)
}

real matrix rdn_binterp(struct rdn_bnd scalar B, real colvector arc)
{
    real matrix out
    real scalar i, a, k, j, t
    out = J(rows(arc), 2, .)
    for (i = 1; i <= rows(arc); i++) {
        a = arc[i]
        if (a >= .) continue
        if (a < 0) a = 0
        if (a > B.segs.total) a = B.segs.total
        k = 1
        for (j = 1; j <= rows(B.segs.cum); j++) {
            if (B.segs.cum[j] <= a) k = j
        }
        t = 0
        if (B.segs.len[k] > 0) {
            t = (a - B.segs.cum[k]) / B.segs.len[k]
            if (t < 0) t = 0
            if (t > 1) t = 1
        }
        out[i, 1] = B.segs.x0[k] + t * (B.segs.x1[k] - B.segs.x0[k])
        out[i, 2] = B.segs.y0[k] + t * (B.segs.y1[k] - B.segs.y0[k])
    }
    return(out)
}

real colvector rdn_btreated(struct rdn_bnd scalar B, real matrix xy)
{
    real colvector arc, left, d, on
    real matrix a, z, foot
    real scalar eps
    if (B.hasregion) {
        d = rdn_polygon_distance(B.region, xy)
        return((d :<= B.tol) :* (d :< .))
    }
    arc = rdn_bproject(B, xy)
    eps = min((1e-3, B.blen * 1e-6))
    a = rdn_binterp(B, rowmax((rowmin((arc :- eps, J(rows(arc), 1, B.blen))), J(rows(arc), 1, 0))))
    z = rdn_binterp(B, rowmax((rowmin((arc :+ eps, J(rows(arc), 1, B.blen))), J(rows(arc), 1, 0))))
    foot = rdn_binterp(B, arc)
    left = ((z[., 1] - a[., 1]) :* (xy[., 2] - foot[., 2]) - (z[., 2] - a[., 2]) :* (xy[., 1] - foot[., 1])) :> 0
    // points on the line (within tol) count as treated, as with a treated region
    on = sqrt((xy[., 1] - foot[., 1]) :^ 2 + (xy[., 2] - foot[., 2]) :^ 2) :<= B.tol
    if (B.side == "left") return(left :| on)
    return((1 :- left) :| on)
}

real colvector rdn_bmasked(struct rdn_bnd scalar B, real colvector arc)
{
    real colvector out, ids, d, fin
    real matrix pts, P
    real scalar i
    out = J(rows(arc), 1, 0)
    fin = (arc :< .)
    for (i = 1; i <= rows(B.maskarcs); i++) {
        out = out :| (fin :& (arc :>= B.maskarcs[i, 1]) :& (arc :<= B.maskarcs[i, 2]))
    }
    if (rows(B.maskpoly)) {
        pts = rdn_binterp(B, arc)
        ids = uniqrows(B.maskpoly[., 3])
        for (i = 1; i <= rows(ids); i++) {
            P = select(B.maskpoly, B.maskpoly[., 3] :== ids[i])
            d = rdn_polygon_distance(P, pts)
            out = out :| ((d :< .) :& (d :<= B.tol))
        }
    }
    return(out)
}

// Where the straight segment p -> q first meets the boundary line:
// (found, x, y, t) with t the position along p -> q.
real rowvector rdn_crossing(real rowvector p, real rowvector q, struct rdn_bnd scalar B)
{
    real scalar rx, ry, best, j, qx, qy, sx, sy, den, t, u, rr, ta, tb, lo, hi
    rx = q[1] - p[1]
    ry = q[2] - p[2]
    best = .
    for (j = 1; j <= rows(B.segs.x0); j++) {
        qx = B.segs.x0[j]
        qy = B.segs.y0[j]
        sx = B.segs.x1[j] - qx
        sy = B.segs.y1[j] - qy
        den = rx * sy - ry * sx
        if (den != 0) {
            t = ((qx - p[1]) * sy - (qy - p[2]) * sx) / den
            u = ((qx - p[1]) * ry - (qy - p[2]) * rx) / den
            if (t >= 0 & t <= 1 & u >= 0 & u <= 1) {
                if (best >= . | t < best) best = t
            }
        }
    }
    if (best >= .) {
        rr = rx * rx + ry * ry
        for (j = 1; j <= rows(B.segs.x0); j++) {
            qx = B.segs.x0[j]
            qy = B.segs.y0[j]
            sx = B.segs.x1[j] - qx
            sy = B.segs.y1[j] - qy
            den = rx * sy - ry * sx
            if (den == 0 & rr > 0 & ((qx - p[1]) * ry - (qy - p[2]) * rx) == 0) {
                ta = ((qx - p[1]) * rx + (qy - p[2]) * ry) / rr
                tb = ta + (sx * rx + sy * ry) / rr
                lo = max((0, min((ta, tb))))
                hi = min((1, max((ta, tb))))
                if (lo <= hi & (best >= . | lo < best)) best = lo
            }
        }
    }
    if (best >= .) return((0, ., ., .))
    return((1, p[1] + best * rx, p[2] + best * ry, best))
}

end

mata:

// ---------------------------------------------------------------------------
// Shortest paths on a (possibly state-expanded) graph in CSR form
// ---------------------------------------------------------------------------

struct rdn_engine {
    real colvector indptr, indices, w, start_state, state_node
    real scalar nstates, nnodes
}

// Workspace reused across searches: only touched entries are reset.
struct rdn_ws {
    real colvector dist, origin, done, touched, hk, hv
    real scalar nt, hn
}

struct rdn_engine scalar rdn_csr(real colvector tail, real colvector head, real colvector w,
    real scalar ns)
{
    struct rdn_engine scalar E
    real colvector cnt, fill
    real scalar m, k, t
    m = rows(tail)
    cnt = J(ns, 1, 0)
    for (k = 1; k <= m; k++) cnt[tail[k]] = cnt[tail[k]] + 1
    E.indptr = 1 \ (runningsum(cnt) :+ 1)
    fill = E.indptr[(1::ns)]
    E.indices = J(m, 1, 0)
    E.w = J(m, 1, 0)
    for (k = 1; k <= m; k++) {
        t = tail[k]
        E.indices[fill[t]] = head[k]
        E.w[fill[t]] = w[k]
        fill[t] = fill[t] + 1
    }
    E.nstates = ns
    E.nnodes = ns
    E.start_state = (1::ns)
    E.state_node = (1::ns)
    return(E)
}

void rdn_ws_init(struct rdn_ws scalar W, real scalar ns, real scalar cap)
{
    W.dist = J(ns, 1, rdn_inf())
    W.origin = J(ns, 1, 0)
    W.done = J(ns, 1, 0)
    W.touched = J(ns, 1, 0)
    W.nt = 0
    W.hk = J(max((cap, 16)), 1, 0)
    W.hv = J(max((cap, 16)), 1, 0)
    W.hn = 0
}

void rdn_ws_reset(struct rdn_ws scalar W)
{
    real scalar k, s
    for (k = 1; k <= W.nt; k++) {
        s = W.touched[k]
        W.dist[s] = rdn_inf()
        W.origin[s] = 0
        W.done[s] = 0
    }
    W.nt = 0
    W.hn = 0
}

void rdn_heap_push(struct rdn_ws scalar W, real scalar key, real scalar val)
{
    real scalar i, p
    W.hn = W.hn + 1
    if (W.hn > rows(W.hk)) {
        W.hk = W.hk \ J(rows(W.hk), 1, 0)
        W.hv = W.hv \ J(rows(W.hv), 1, 0)
    }
    i = W.hn
    while (i > 1) {
        p = floor(i / 2)
        if (W.hk[p] <= key) break
        W.hk[i] = W.hk[p]
        W.hv[i] = W.hv[p]
        i = p
    }
    W.hk[i] = key
    W.hv[i] = val
}

// Pops the smallest key into (key, val).
void rdn_heap_pop(struct rdn_ws scalar W, real scalar key, real scalar val)
{
    real scalar lk, lv, n, i, c
    key = W.hk[1]
    val = W.hv[1]
    lk = W.hk[W.hn]
    lv = W.hv[W.hn]
    W.hn = W.hn - 1
    n = W.hn
    if (n == 0) return
    i = 1
    while (1) {
        c = 2 * i
        if (c > n) break
        if (c + 1 <= n) {
            if (W.hk[c + 1] < W.hk[c]) c = c + 1
        }
        if (W.hk[c] >= lk) break
        W.hk[i] = W.hk[c]
        W.hv[i] = W.hv[c]
        i = c
    }
    W.hk[i] = lk
    W.hv[i] = lv
}

// Dijkstra from start states (all at distance 0) up to `limit`. Each state
// remembers the label of the source it was reached from; ties go to the
// smaller label.
void rdn_dijkstra(struct rdn_engine scalar E, real colvector starts, real colvector labels,
    real scalar limit, struct rdn_ws scalar W)
{
    real scalar k, s, d, u, ed, v, nd
    for (k = 1; k <= rows(starts); k++) {
        s = starts[k]
        if (W.dist[s] >= rdn_inf()) {
            W.nt = W.nt + 1
            W.touched[W.nt] = s
        }
        if (W.dist[s] > 0 | W.origin[s] == 0) {
            W.dist[s] = 0
            W.origin[s] = labels[k]
            rdn_heap_push(W, 0, s)
        }
    }
    while (W.hn > 0) {
        d = .
        u = .
        rdn_heap_pop(W, d, u)
        if (W.done[u]) continue
        if (d > W.dist[u]) continue
        W.done[u] = 1
        for (ed = E.indptr[u]; ed < E.indptr[u + 1]; ed++) {
            v = E.indices[ed]
            nd = d + E.w[ed]
            if (nd > limit) continue
            if (W.done[v]) continue
            if (nd < W.dist[v]) {
                if (W.dist[v] >= rdn_inf()) {
                    W.nt = W.nt + 1
                    W.touched[W.nt] = v
                }
                W.dist[v] = nd
                W.origin[v] = W.origin[u]
                rdn_heap_push(W, nd, v)
            }
            else if (nd == W.dist[v] & W.origin[u] < W.origin[v]) {
                W.origin[v] = W.origin[u]
            }
        }
    }
}

// Distances from each source node to every node (k x n), 1e300 if unreachable.
real matrix rdn_distances(struct rdn_engine scalar E, real colvector sources, real scalar limit)
{
    struct rdn_ws scalar W
    real matrix out
    real scalar i, t, s, nd
    out = J(rows(sources), E.nnodes, rdn_inf())
    rdn_ws_init(W, E.nstates, rows(E.indices) + rows(sources) + 2)
    for (i = 1; i <= rows(sources); i++) {
        rdn_dijkstra(E, E.start_state[sources[i]], 1, limit, W)
        for (t = 1; t <= W.nt; t++) {
            s = W.touched[t]
            nd = E.state_node[s]
            if (W.dist[s] < out[i, nd]) out[i, nd] = W.dist[s]
        }
        rdn_ws_reset(W)
    }
    return(out)
}

// For every node: distance to the nearest source (dist) and that source's
// node index (src, 0 if unreachable).
void rdn_nearest_source(struct rdn_engine scalar E, real colvector sources,
    real colvector dist, real colvector src)
{
    struct rdn_ws scalar W
    real scalar s, nd
    rdn_ws_init(W, E.nstates, rows(E.indices) + rows(sources) + 2)
    rdn_dijkstra(E, E.start_state[sources], (1::rows(sources)), rdn_inf(), W)
    dist = J(E.nnodes, 1, rdn_inf())
    src = J(E.nnodes, 1, 0)
    for (s = 1; s <= E.nstates; s++) {
        nd = E.state_node[s]
        if (W.dist[s] < dist[nd]) {
            dist[nd] = W.dist[s]
            src[nd] = sources[W.origin[s]]
        }
    }
}

// Every (source position, node, distance) with distance strictly below `limit`.
real matrix rdn_pairs(struct rdn_engine scalar E, real colvector sources, real scalar limit)
{
    struct rdn_ws scalar W
    real matrix buf
    real colvector best, seen
    real scalar i, t, s, nd, d, ns, cnt, k
    rdn_ws_init(W, E.nstates, rows(E.indices) + 2)
    best = J(E.nnodes, 1, rdn_inf())
    seen = J(E.nnodes, 1, 0)
    buf = J(max((1024, 4 * rows(sources))), 3, .)
    cnt = 0
    for (i = 1; i <= rows(sources); i++) {
        rdn_dijkstra(E, E.start_state[sources[i]], 1, limit, W)
        ns = 0
        for (t = 1; t <= W.nt; t++) {
            s = W.touched[t]
            d = W.dist[s]
            if (!(d < limit)) continue
            nd = E.state_node[s]
            if (best[nd] >= rdn_inf()) {
                ns++
                seen[ns] = nd
            }
            if (d < best[nd]) best[nd] = d
        }
        for (k = 1; k <= ns; k++) {
            nd = seen[k]
            cnt++
            if (cnt > rows(buf)) buf = buf \ J(rows(buf), 3, .)
            buf[cnt, .] = (i, nd, best[nd])
            best[nd] = rdn_inf()
        }
        rdn_ws_reset(W)
    }
    if (cnt == 0) return(J(0, 3, .))
    return(buf[(1::cnt), .])
}

// ---------------------------------------------------------------------------
// Networks and metrics
// ---------------------------------------------------------------------------

struct rdn_net {
    string colvector ids, eclass, connectors, penrows, pencols
    real matrix xy, penmat
    real colvector efrom, eto, ecost, elen, iscross
    real colvector fu, fv, flen, fcost, fglen
    real scalar n, directed, metric, defpen, snapscale, ncross, flow
    string scalar direction
    struct rdn_engine scalar E, E2
}

// Distances, nearest sources and close pairs on the network: with
// direction(flow), the shorter of the directed paths a -> b (engine E, edges
// as drawn) and b -> a (engine E2, edges reversed).
real matrix rdn_net_distances(struct rdn_net scalar N, real colvector sources, real scalar limit)
{
    real matrix A, B
    A = rdn_distances(N.E, sources, limit)
    if (!N.flow) return(A)
    B = rdn_distances(N.E2, sources, limit)
    return(A :* (A :<= B) + B :* (B :< A))
}

void rdn_net_nearest(struct rdn_net scalar N, real colvector sources, real colvector dist, real colvector src)
{
    real colvector d2, s2, take
    rdn_nearest_source(N.E, sources, dist, src)
    if (!N.flow) return
    d2 = .
    s2 = .
    rdn_nearest_source(N.E2, sources, d2, s2)
    take = (d2 :< dist)
    dist = rdn_where(take, d2, dist)
    src = rdn_where(take, s2, src)
}

real matrix rdn_net_pairs(struct rdn_net scalar N, real colvector sources, real scalar limit)
{
    real matrix P, o
    real colvector keep
    real scalar r
    P = rdn_pairs(N.E, sources, limit)
    if (!N.flow) return(P)
    P = P \ rdn_pairs(N.E2, sources, limit)
    if (rows(P) == 0) return(P)
    P = P[order(P, (1, 2, 3)), .]
    keep = J(rows(P), 1, 1)
    for (r = 2; r <= rows(P); r++) keep[r] = (P[r, 1] != P[r - 1, 1]) | (P[r, 2] != P[r - 1, 2])
    return(select(P, keep))
}

// Edge table: one row per undirected street (of parallel edges, the cheapest,
// which is the one shortest paths use).
void rdn_edge_frame(struct rdn_net scalar N)
{
    real colvector a, b, key, o, keep, idx, dx, dy, el
    real scalar m, i
    m = rows(N.efrom)
    a = rowmin((N.efrom, N.eto))
    b = rowmax((N.efrom, N.eto))
    key = a :* (N.n + 1) :+ b
    o = order((key, N.ecost, (1::m)), (1, 2, 3))
    keep = J(m, 1, 0)
    for (i = 1; i <= m; i++) {
        if (i == 1) keep[o[i]] = 1
        else if (key[o[i]] != key[o[i - 1]]) keep[o[i]] = 1
    }
    idx = selectindex(keep)
    N.fu = N.efrom[idx]
    N.fv = N.eto[idx]
    N.fcost = N.ecost[idx]
    dx = N.xy[N.fu, 1] - N.xy[N.fv, 1]
    dy = N.xy[N.fu, 2] - N.xy[N.fv, 2]
    N.fglen = sqrt(dx :^ 2 + dy :^ 2)
    el = N.elen[idx]
    N.flen = rdn_where(el :< ., el, N.fglen)
}

real scalar rdn_penalty(struct rdn_net scalar N, string scalar a, string scalar b)
{
    real scalar i, j, v
    if (rows(N.penmat)) {
        i = 0
        j = 0
        for (v = 1; v <= rows(N.penrows); v++) if (N.penrows[v] == a) i = v
        for (v = 1; v <= rows(N.pencols); v++) if (N.pencols[v] == b) j = v
        if (i > 0 & j > 0) {
            if (N.penmat[i, j] < .) {
                if (N.penmat[i, j] < 0) _error(198, "penalties must be >= 0")
                return(N.penmat[i, j])
            }
        }
    }
    return(a == b ? 0 : N.defpen)
}

// Compile the shortest-path engine (edge costs, or the state-expanded graph
// of transition costs).
void rdn_compile(struct rdn_net scalar N)
{
    real colvector tail, head, cost, isconn, nc, cn, cid, ncstate, connbase, state_node, state_mem
    real colvector rr, cc, vv
    string colvector cls, classes
    real matrix pen
    real scalar mm, C, nmem, a, b, k, nstates, tot, s, mem, ed, p
    struct rdn_engine scalar A, E1
    transmorphic map
    N.flow = 0
    if (N.directed & N.direction == "flow") {
        // the shorter of the directed paths a -> b and b -> a: engines for the
        // edges as drawn (E) and reversed (E2)
        N.direction = "out"
        rdn_compile(N)
        E1 = N.E
        N.direction = "in"
        rdn_compile(N)
        N.E2 = N.E
        N.E = E1
        N.direction = "flow"
        N.flow = 1
        return
    }
    if (!N.directed | N.direction == "both") {
        tail = N.efrom \ N.eto
        head = N.eto \ N.efrom
        cost = N.ecost \ N.ecost
        cls = N.eclass \ N.eclass
    }
    else if (N.direction == "out") {
        tail = N.efrom
        head = N.eto
        cost = N.ecost
        cls = N.eclass
    }
    else {
        tail = N.eto
        head = N.efrom
        cost = N.ecost
        cls = N.eclass
    }
    if (N.metric == 1) {
        N.E = rdn_csr(tail, head, cost, N.n)
        return
    }
    mm = rows(tail)
    isconn = J(mm, 1, 0)
    for (k = 1; k <= mm; k++) {
        if (rows(N.connectors)) isconn[k] = anyof(N.connectors, cls[k])
    }
    nc = selectindex(!isconn)
    cn = selectindex(isconn)
    if (rows(nc)) classes = uniqrows(cls[nc])
    else classes = J(0, 1, "")
    C = rows(classes)
    nmem = C + 1
    pen = J(nmem, nmem, 0)
    map = asarray_create()
    for (a = 1; a <= C; a++) {
        asarray(map, classes[a], a)
        for (b = 1; b <= C; b++) pen[a + 1, b + 1] = rdn_penalty(N, classes[a], classes[b])
    }
    cid = J(mm, 1, 1)
    for (k = 1; k <= rows(nc); k++) cid[nc[k]] = 1 + asarray(map, cls[nc[k]])
    ncstate = J(mm, 1, 0)
    connbase = J(mm, 1, 0)
    if (rows(nc)) ncstate[nc] = N.n :+ (1::rows(nc))
    if (rows(cn)) connbase[cn] = (N.n + rows(nc)) :+ ((1::rows(cn)) :- 1) :* nmem
    nstates = N.n + rows(nc) + rows(cn) * nmem
    state_node = (1::N.n)
    state_mem = J(N.n, 1, 1)
    if (rows(nc)) {
        state_node = state_node \ head[nc]
        state_mem = state_mem \ cid[nc]
    }
    if (rows(cn)) {
        state_node = state_node \ (head[cn] # J(nmem, 1, 1))
        state_mem = state_mem \ (J(rows(cn), 1, 1) # (1::nmem))
    }
    // out-arcs of every node
    A = rdn_csr(tail, (1::mm), J(mm, 1, 0), N.n)
    tot = 0
    for (s = 1; s <= nstates; s++) {
        a = state_node[s]
        tot = tot + A.indptr[a + 1] - A.indptr[a]
    }
    rr = J(tot, 1, 0)
    cc = J(tot, 1, 0)
    vv = J(tot, 1, 0)
    p = 0
    for (s = 1; s <= nstates; s++) {
        a = state_node[s]
        mem = state_mem[s]
        for (ed = A.indptr[a]; ed < A.indptr[a + 1]; ed++) {
            k = A.indices[ed]
            p++
            rr[p] = s
            if (isconn[k]) {
                cc[p] = connbase[k] + mem
                vv[p] = cost[k]
            }
            else {
                cc[p] = ncstate[k]
                vv[p] = cost[k] + pen[mem, cid[k]]
            }
        }
    }
    N.E = rdn_csr(rr, cc, vv, nstates)
    N.E.start_state = (1::N.n)
    N.E.state_node = state_node
    N.E.nnodes = N.n
}

// Node ids as strings (numeric ids are formatted without padding).
string colvector rdn_idstrings(string scalar var)
{
    if (st_isstrvar(var)) return(strtrim(st_sdata(., var)))
    return(strtrim(strofreal(st_data(., var), "%21.0g")))
}

// Read the network from the node and edge frames.
void rdn_read_network(struct rdn_net scalar N, string scalar fnodes, string scalar fedges,
    string scalar idv, string scalar xv, string scalar yv, string scalar fromv, string scalar tov,
    string scalar costv, string scalar lenv, string scalar classv)
{
    string scalar cur
    string colvector from, to
    real colvector fi, ti, keep
    transmorphic map
    real scalar i
    cur = st_framecurrent()
    st_framecurrent(fnodes)
    if (_st_varindex(idv) >= . | _st_varindex(xv) >= . | _st_varindex(yv) >= .) {
        st_framecurrent(cur)
        _error(111, "frame " + fnodes + " needs variables " + idv + ", " + xv + " and " + yv)
    }
    N.ids = rdn_idstrings(idv)
    N.xy = st_data(., (xv, yv))
    st_framecurrent(fedges)
    if (_st_varindex(fromv) >= . | _st_varindex(tov) >= . | _st_varindex(costv) >= .) {
        st_framecurrent(cur)
        _error(111, "frame " + fedges + " needs variables " + fromv + ", " + tov + " and " + costv)
    }
    from = rdn_idstrings(fromv)
    to = rdn_idstrings(tov)
    N.ecost = st_data(., costv)
    if (lenv != "" & _st_varindex(lenv) < .) N.elen = st_data(., lenv)
    else N.elen = J(rows(from), 1, .)
    if (classv != "" & _st_varindex(classv) < .) N.eclass = rdn_idstrings(classv)
    else N.eclass = J(rows(from), 1, "")
    st_framecurrent(cur)
    if (hasmissing(N.xy)) _error(416, "every node needs x and y coordinates")
    N.n = rows(N.ids)
    map = asarray_create()
    for (i = 1; i <= N.n; i++) {
        if (asarray_contains(map, N.ids[i])) _error(459, "node ids must be unique (" + N.ids[i] + ")")
        asarray(map, N.ids[i], i)
    }
    fi = J(rows(from), 1, .)
    ti = J(rows(from), 1, .)
    for (i = 1; i <= rows(from); i++) {
        if (!asarray_contains(map, from[i]) | !asarray_contains(map, to[i])) {
            _error(459, "edge (" + from[i] + ", " + to[i] + ") refers to a node that is not in the node frame")
        }
        fi[i] = asarray(map, from[i])
        ti[i] = asarray(map, to[i])
    }
    for (i = 1; i <= rows(from); i++) {
        if (N.ecost[i] >= . | N.ecost[i] < 0) {
            _error(459, "invalid cost for edge (" + from[i] + ", " + to[i] + "); costs must be >= 0 and not missing")
        }
    }
    keep = selectindex(fi :!= ti)
    N.efrom = fi[keep]
    N.eto = ti[keep]
    N.ecost = N.ecost[keep]
    N.elen = N.elen[keep]
    N.eclass = N.eclass[keep]
    for (i = 1; i <= rows(N.eclass); i++) if (N.eclass[i] == "") N.eclass[i] = "unclassified"
    N.iscross = J(N.n, 1, 0)
    N.ncross = 0
}

// Read a network given as line segments with endpoint coordinates (x1 y1 x2 y2)
// and no node ids: segment ends within jtol of each other become one node.
void rdn_read_lines(struct rdn_net scalar N, string scalar fedges, string rowvector ecoords,
    string scalar costv, string scalar lenv, string scalar classv, real scalar jtol)
{
    string scalar cur
    real matrix seg, ends
    real colvector par, o, lab, first, glen, keep
    real scalar m, n2, i, j, a, b, ra, rb, k
    cur = st_framecurrent()
    st_framecurrent(fedges)
    for (i = 1; i <= 4; i++) {
        if (_st_varindex(ecoords[i]) >= .) {
            st_framecurrent(cur)
            _error(111, "frame " + fedges + " needs the segment end coordinates " + invtokens(ecoords))
        }
    }
    seg = st_data(., ecoords)
    m = rows(seg)
    glen = sqrt((seg[., 3] - seg[., 1]) :^ 2 + (seg[., 4] - seg[., 2]) :^ 2)
    if (costv != "" & _st_varindex(costv) < .) N.ecost = st_data(., costv)
    else if (costv == "length") N.ecost = glen
    else {
        st_framecurrent(cur)
        _error(111, "frame " + fedges + " has no cost variable " + costv)
    }
    if (lenv != "" & _st_varindex(lenv) < .) N.elen = st_data(., lenv)
    else N.elen = glen
    if (classv != "" & _st_varindex(classv) < .) N.eclass = rdn_idstrings(classv)
    else N.eclass = J(m, 1, "")
    st_framecurrent(cur)
    if (hasmissing(seg)) _error(416, "segment end coordinates must not be missing")
    ends = seg[., (1, 2)] \ seg[., (3, 4)]
    n2 = rows(ends)
    // union-find over ends closer than jtol (sweep in x order)
    par = (1::n2)
    o = order(ends, (1, 2))
    for (i = 1; i <= n2; i++) {
        a = o[i]
        for (j = i + 1; j <= n2; j++) {
            b = o[j]
            if (ends[b, 1] - ends[a, 1] > jtol) break
            if ((ends[b, 1] - ends[a, 1])^2 + (ends[b, 2] - ends[a, 2])^2 <= jtol^2) {
                ra = a
                while (par[ra] != ra) ra = par[ra]
                rb = b
                while (par[rb] != rb) rb = par[rb]
                if (ra != rb) par[max((ra, rb))] = min((ra, rb))
            }
        }
    }
    for (i = 1; i <= n2; i++) {
        ra = i
        while (par[ra] != ra) ra = par[ra]
        par[i] = ra
    }
    // node labels in order of first appearance
    lab = J(n2, 1, 0)
    first = J(n2, 1, 0)
    k = 0
    for (i = 1; i <= n2; i++) {
        if (first[par[i]] == 0) {
            k++
            first[par[i]] = k
        }
        lab[i] = first[par[i]]
    }
    N.n = k
    N.xy = J(k, 2, .)
    for (i = n2; i >= 1; i--) N.xy[lab[i], .] = ends[i, .]
    N.ids = strofreal(1::k)
    N.efrom = lab[(1::m)]
    N.eto = lab[((m + 1)::n2)]
    for (i = 1; i <= m; i++) {
        if (N.ecost[i] >= . | N.ecost[i] < 0) {
            _error(459, "invalid cost for segment " + strofreal(i) + "; costs must be >= 0 and not missing")
        }
    }
    keep = selectindex(N.efrom :!= N.eto)
    if (rows(keep) < m) printf("{txt}note: %g zero-length or closed segments dropped\n", m - rows(keep))
    N.efrom = N.efrom[keep]
    N.eto = N.eto[keep]
    N.ecost = N.ecost[keep]
    N.elen = N.elen[keep]
    N.eclass = N.eclass[keep]
    for (i = 1; i <= rows(N.eclass); i++) if (N.eclass[i] == "") N.eclass[i] = "unclassified"
    N.iscross = J(N.n, 1, 0)
    N.ncross = 0
}

// The side of a line ("left" or "right") on which the treated observations lie.
string scalar rdn_infer_side(real matrix line, real matrix xy, real colvector T)
{
    struct rdn_bnd scalar BL
    real colvector left, arc, inner, TT
    real scalar agree, share
    if (rows(T) != rows(xy)) _error(198, "treated() is needed to infer the treated side")
    if (sum(T) == 0 | sum(T) == rows(T)) _error(198, "treated() must contain treated and control observations")
    BL = rdn_boundary(line, J(0, 2, .), "left", J(0, 2, .), J(0, 3, .), 1e-6)
    left = rdn_btreated(BL, xy)
    // observations beyond the ends of the line have no clear side: use those whose
    // nearest point on the line is inside it, when there are enough of both groups
    arc = rdn_bproject(BL, xy)
    inner = selectindex((arc :> 1e-9 * BL.blen) :& (arc :< BL.blen * (1 - 1e-9)))
    TT = T
    if (rows(inner)) {
        if (sum(T[inner]) >= 5 & sum(!T[inner]) >= 5) {
            left = left[inner]
            TT = T[inner]
        }
    }
    agree = mean(left :== TT)
    share = max((agree, 1 - agree))
    if (share < 0.9) {
        printf("{txt}warning: only %4.1f%% of observations lie on the side of the boundary that matches treated();\n", 100 * share)
        printf("{txt}         check that evalpoints() are ordered along the boundary\n")
    }
    return(agree >= 0.5 ? "left" : "right")
}

// Split every edge that crosses the boundary at the crossing point, so that
// every path between the two sides passes through a node on the boundary.
// Costs and lengths are divided in proportion to the split.
void rdn_split(struct rdn_net scalar N, struct rdn_bnd scalar B)
{
    real colvector T, on, cand, drop, keep, nf, nt, nco, nle
    string colvector ncl, newids
    real matrix newxy
    real rowvector hit
    transmorphic keys
    string scalar key
    real scalar m, c, u, v, cid, t
    T = rdn_btreated(B, N.xy)
    on = (rdn_bdistance(B, N.xy) :<= B.tol)
    m = rows(N.efrom)
    cand = selectindex((!on[N.efrom]) :& (!on[N.eto]) :& (T[N.efrom] :!= T[N.eto]))
    N.ncross = 0
    if (rows(cand) == 0) return
    keys = asarray_create()
    newxy = J(0, 2, .)
    newids = J(0, 1, "")
    drop = J(m, 1, 0)
    nf = J(0, 1, .)
    nt = J(0, 1, .)
    nco = J(0, 1, .)
    nle = J(0, 1, .)
    ncl = J(0, 1, "")
    for (c = 1; c <= rows(cand); c++) {
        u = N.efrom[cand[c]]
        v = N.eto[cand[c]]
        hit = rdn_crossing(N.xy[u, .], N.xy[v, .], B)
        if (hit[1] == 0) continue
        t = hit[4]
        key = strofreal(round(hit[2], 1e-6), "%20.6f") + " " + strofreal(round(hit[3], 1e-6), "%20.6f")
        if (!asarray_contains(keys, key)) {
            newxy = newxy \ (hit[2], hit[3])
            newids = newids \ ("rdnet_crossing_" + strofreal(rows(newxy) - 1))
            asarray(keys, key, N.n + rows(newxy))
        }
        cid = asarray(keys, key)
        nf = nf \ u \ cid
        nt = nt \ cid \ v
        nco = nco \ N.ecost[cand[c]] * t \ N.ecost[cand[c]] * (1 - t)
        nle = nle \ N.elen[cand[c]] * t \ N.elen[cand[c]] * (1 - t)
        ncl = ncl \ N.eclass[cand[c]] \ N.eclass[cand[c]]
        drop[cand[c]] = 1
    }
    if (sum(drop) == 0) return
    keep = selectindex(!drop)
    N.efrom = N.efrom[keep] \ nf
    N.eto = N.eto[keep] \ nt
    N.ecost = N.ecost[keep] \ nco
    N.elen = N.elen[keep] \ nle
    N.eclass = N.eclass[keep] \ ncl
    N.xy = N.xy \ newxy
    N.ids = N.ids \ newids
    N.iscross = N.iscross \ J(rows(newxy), 1, 1)
    N.n = rows(N.xy)
    N.ncross = sum(drop)
}

// Whether locations on edges must become nodes for exact distances. A location
// partway along an edge is normally reached through the edge's end nodes plus
// part of the edge's cost; that ignores one-way edges and the penalty for
// switching onto the location's edge class, so directed networks and
// transition costs place each location on a node instead.
real scalar rdn_needs_point_nodes(struct rdn_net scalar N)
{
    return((N.directed & N.direction != "both") | N.metric == 2)
}

// Insert a node at each location frac of the way from node u to node v,
// splitting every stored edge between u and v (both directions, parallel
// edges) with cost and length divided in proportion. Returns the node of each
// location (an end node when frac is 0 or 1).
real colvector rdn_point_nodes(struct rdn_net scalar N, real colvector u, real colvector v, real colvector frac0)
{
    real colvector frac, node, atu, atv, mid, lo, hi, t, slo, shi, st, sp, ts, ids, o, cut, chain
    real colvector drop, keep, nf, nt, nco, nle
    string colvector ncl
    transmorphic spots, pairs
    string scalar key, pkey
    real scalar n, i, k, m, ns, a, b, share, fwd, r
    n = rows(u)
    frac = rowmin((rowmax((frac0, J(n, 1, 0))), J(n, 1, 1)))
    atu = (u :== v) :| (frac :<= 1e-9)
    atv = (!atu) :& (frac :>= 1 - 1e-9)
    node = rdn_where(atu, u, rdn_where(atv, v, J(n, 1, .)))
    mid = selectindex((!atu) :& (!atv))
    if (rows(mid) == 0) return(node)
    // positions measured from the smaller node index
    lo = rowmin((u[mid], v[mid]))
    hi = rowmax((u[mid], v[mid]))
    t = rdn_where(u[mid] :== lo, frac[mid], 1 :- frac[mid])
    spots = asarray_create()
    pairs = asarray_create()
    slo = J(0, 1, .)
    shi = J(0, 1, .)
    st = J(0, 1, .)
    for (i = 1; i <= rows(mid); i++) {
        key = strofreal(lo[i]) + " " + strofreal(hi[i]) + " " + strofreal(t[i], "%21.12f")
        if (!asarray_contains(spots, key)) {
            slo = slo \ lo[i]
            shi = shi \ hi[i]
            st = st \ strtoreal(strofreal(t[i], "%21.12f"))
            asarray(spots, key, rows(st))
            pkey = strofreal(lo[i]) + " " + strofreal(hi[i])
            if (asarray_contains(pairs, pkey)) asarray(pairs, pkey, asarray(pairs, pkey) \ rows(st))
            else asarray(pairs, pkey, rows(st))
        }
        node[mid[i]] = N.n + asarray(spots, key)
    }
    ns = rows(st)
    m = rows(N.efrom)
    drop = J(m, 1, 0)
    nf = J(0, 1, .)
    nt = J(0, 1, .)
    nco = J(0, 1, .)
    nle = J(0, 1, .)
    ncl = J(0, 1, "")
    for (k = 1; k <= m; k++) {
        a = min((N.efrom[k], N.eto[k]))
        b = max((N.efrom[k], N.eto[k]))
        if (a == b) continue
        pkey = strofreal(a) + " " + strofreal(b)
        if (!asarray_contains(pairs, pkey)) continue
        sp = asarray(pairs, pkey)
        o = order(st[sp], 1)
        ts = st[sp[o]]
        ids = N.n :+ sp[o]
        fwd = (N.efrom[k] == a)
        if (!fwd) {
            ts = 1 :- ts[rows(ts)::1]
            ids = ids[rows(ids)::1]
        }
        chain = N.efrom[k] \ ids \ N.eto[k]
        cut = 0 \ ts \ 1
        for (r = 1; r < rows(chain); r++) {
            share = cut[r + 1] - cut[r]
            nf = nf \ chain[r]
            nt = nt \ chain[r + 1]
            nco = nco \ N.ecost[k] * share
            nle = nle \ N.elen[k] * share
            ncl = ncl \ N.eclass[k]
        }
        drop[k] = 1
    }
    keep = selectindex(!drop)
    N.efrom = N.efrom[keep] \ nf
    N.eto = N.eto[keep] \ nt
    N.ecost = N.ecost[keep] \ nco
    N.elen = N.elen[keep] \ nle
    N.eclass = N.eclass[keep] \ ncl
    N.xy = N.xy \ (N.xy[slo, .] + (N.xy[shi, .] - N.xy[slo, .]) :* st)
    N.ids = N.ids \ ("rdnet_point_" :+ strofreal((0::(ns - 1))))
    N.iscross = N.iscross \ J(ns, 1, 0)
    N.n = rows(N.xy)
    return(node)
}

// ---------------------------------------------------------------------------
// Snapping points to the network
// ---------------------------------------------------------------------------

struct rdn_snap {
    real colvector u, v, cu, cv, depth, rate
}

// Cost per unit length of each street in the edge table (cost / length), with
// the median for streets without a usable length; and the mean over each
// node's streets.
real colvector rdn_edge_rates(struct rdn_net scalar N)
{
    real colvector r, good
    real scalar fb
    r = N.fcost :/ N.flen
    good = (r :< .) :& (N.flen :> 0) :& (N.fcost :< rdn_inf())
    fb = (sum(good) ? rdn_median(select(r, good)) : 0)
    return(rdn_where(good, r, J(rows(r), 1, fb)))
}

real colvector rdn_node_rates(struct rdn_net scalar N, real colvector er)
{
    real colvector sm, ct
    real scalar k, fb
    sm = J(N.n, 1, 0)
    ct = J(N.n, 1, 0)
    for (k = 1; k <= rows(er); k++) {
        sm[N.fu[k]] = sm[N.fu[k]] + er[k]
        sm[N.fv[k]] = sm[N.fv[k]] + er[k]
        ct[N.fu[k]] = ct[N.fu[k]] + 1
        ct[N.fv[k]] = ct[N.fv[k]] + 1
    }
    fb = (rows(er) ? rdn_median(er) : 0)
    return(rdn_where(ct :> 0, sm :/ rowmax((ct, J(N.n, 1, 1))), J(N.n, 1, fb)))
}

// Attach points to the nearest allowed street segment (or node).
struct rdn_snap scalar rdn_snap_points(struct rdn_net scalar N, real matrix xy,
    real colvector allowed, real scalar nodemode)
{
    struct rdn_snap scalar S
    struct rdn_near scalar R
    real colvector cand, ok, ei, ed, er, nr, atu, atv, tiny
    if (nodemode) {
        cand = selectindex(allowed)
        if (rows(cand) == 0) _error(3498, "no network nodes available for snapping")
        R = rdn_nearest_segment(xy[., 1], xy[., 2], N.xy[cand, 1], N.xy[cand, 2],
            N.xy[cand, 1], N.xy[cand, 2])
        S.u = cand[R.seg]
        S.v = S.u
        S.cu = J(rows(xy), 1, 0)
        S.cv = J(rows(xy), 1, 0)
        S.depth = R.dist
        S.rate = rdn_node_rates(N, rdn_edge_rates(N))[S.u]
        return(S)
    }
    ok = allowed[N.fu] :& allowed[N.fv] :& (N.fcost :< rdn_inf())
    ei = selectindex(ok)
    if (rows(ei) == 0) _error(3498, "no network edges available for snapping")
    R = rdn_nearest_segment(xy[., 1], xy[., 2], N.xy[N.fu[ei], 1], N.xy[N.fu[ei], 2],
        N.xy[N.fv[ei], 1], N.xy[N.fv[ei], 2])
    ed = ei[R.seg]
    S.u = N.fu[ed]
    S.v = N.fv[ed]
    S.cu = R.t :* N.fcost[ed]
    S.cv = (1 :- R.t) :* N.fcost[ed]
    S.depth = R.dist
    // a point at an end node uses the node's rate (which of its edges it
    // "snapped to" is arbitrary)
    er = rdn_edge_rates(N)
    nr = rdn_node_rates(N, er)
    S.rate = er[ed]
    tiny = 1e-9 :* (S.cu + S.cv)  // projections onto an edge's end, up to rounding
    atv = (S.u :!= S.v) :& (S.cv :<= tiny) :& (S.cu :> tiny)
    atu = (S.u :== S.v) :| ((S.cu :<= tiny) :& (!atv))
    S.rate = rdn_where(atu, nr[S.u], rdn_where(atv, nr[S.v], S.rate))
    return(S)
}

end

mata:

// ---------------------------------------------------------------------------
// Network coordinates around the boundary
// ---------------------------------------------------------------------------
// For observation i and evaluation point b on the boundary:
//   s_i     signed network distance to the boundary (negative on the control side)
//   o_i(b)  signed offset along the boundary from b to where i reaches it
//   d_i(b)  network distance from b to i (kernel distance)
// Location-level arrays (L*, C*) hold one entry per distinct location; obs-level
// arrays are indexed by loc.

struct rdn_coords {
    struct rdn_net scalar N
    struct rdn_bnd scalar B
    real scalar ready, nobs, neval, nloc, nunreach
    real matrix xy, exy, NR
    real colvector treated, s, crossarc, depth, loc, firstl
    real colvector Lu, Lv, Lcu, Lcv, Ldep, Ls
    real colvector Cu, Cv, Ccu, Ccv, Carc, Cnode
    real colvector nodeT, nodeon, nodedist, nodesrc
    real colvector evalarc, eu, ev, ecu, ecv, enode, enodearc, emasked
}

// (loc, first): location id of each point and one point per location (by x,
// y and side), numbered in order of first appearance.
void rdn_unique_locations(real matrix xy, real colvector T, real scalar useT,
    real colvector loc, real colvector firstl)
{
    real colvector o, isnew, g, loc0, firsts, srt, rk
    real scalar n, i, a, b, chg
    n = rows(xy)
    if (useT) o = order((xy[., (1, 2)], T, (1::n)), (1, 2, 3, 4))
    else o = order((xy[., (1, 2)], (1::n)), (1, 2, 3))
    isnew = J(n, 1, 1)
    for (i = 2; i <= n; i++) {
        a = o[i]
        b = o[i - 1]
        chg = (xy[a, 1] != xy[b, 1]) | (xy[a, 2] != xy[b, 2])
        if (useT) chg = chg | (T[a] != T[b])
        isnew[i] = chg
    }
    g = runningsum(isnew)
    loc0 = J(n, 1, .)
    loc0[o] = g
    firsts = o[selectindex(isnew)]
    srt = order(firsts, 1)
    rk = J(rows(firsts), 1, .)
    rk[srt] = (1::rows(firsts))
    loc = rk[loc0]
    firstl = firsts[srt]
}

// Place evaluation points on boundary street segments (or boundary nodes).
void rdn_snap_eval(struct rdn_coords scalar C, real colvector bnodes, real colvector bnarc)
{
    real colvector along, ei, au, av, lo, hi, inside
    real scalar k, j, a, r, frac
    k = rows(C.evalarc)
    C.eu = J(k, 1, .)
    for (j = 1; j <= k; j++) C.eu[j] = bnodes[rdn_argmin(abs(bnarc :- C.evalarc[j]))]
    C.ev = C.eu
    C.ecu = J(k, 1, 0)
    C.ecv = J(k, 1, 0)
    along = C.nodeon[C.N.fu] :& C.nodeon[C.N.fv] :& (C.N.fcost :< rdn_inf())
    ei = selectindex(along)
    if (rows(ei) == 0) return
    au = rdn_bproject(C.B, C.N.xy[C.N.fu[ei], .])
    av = rdn_bproject(C.B, C.N.xy[C.N.fv[ei], .])
    lo = rowmin((au, av))
    hi = rowmax((au, av))
    for (j = 1; j <= k; j++) {
        a = C.evalarc[j]
        if (a >= .) continue
        inside = selectindex((lo :<= a) :& (a :<= hi) :& (hi :> lo))
        if (rows(inside) == 0) continue
        r = inside[1]
        frac = (a - au[r]) / (av[r] - au[r])
        C.eu[j] = C.N.fu[ei[r]]
        C.ev[j] = C.N.fv[ei[r]]
        C.ecu[j] = frac * C.N.fcost[ei[r]]
        C.ecv[j] = (1 - frac) * C.N.fcost[ei[r]]
    }
}

// Build the coordinates of the points in xy. Tin is the user's treatment
// indicator (empty: side of the boundary). evalarc is empty for the default
// grid of neval points.
void rdn_build_coords(struct rdn_coords scalar C, real matrix xy, real colvector Tin,
    real colvector evalarc, real scalar neval, real scalar nodemode)
{
    struct rdn_snap scalar S
    real colvector bnodes, loc, firstl, TL, u, v, cu, cv, dep, m, allowed, depc
    real colvector viau, viav, dist, xnode, onedge, unreach, safe, bnarc, tot, frac, idx
    real colvector pu, pv, pc, onends, pnode, off, rate
    struct rdn_net scalar NP
    real matrix xyl, D
    real scalar n, side, j, Jn, useT
    C.nodeon = (rdn_bdistance(C.B, C.N.xy) :<= C.B.tol)
    C.nodeT = rdn_btreated(C.B, C.N.xy) :& (!C.nodeon)
    bnodes = selectindex(C.nodeon)
    if (rows(bnodes) == 0) {
        _error(3498, "no network nodes lie on the boundary; check the coordinates or increase tol()")
    }
    useT = (rows(Tin) > 0)
    loc = .
    firstl = .
    rdn_unique_locations(xy, Tin, useT, loc, firstl)
    xyl = xy[firstl, .]
    if (useT) {
        TL = Tin[firstl]
        off = (TL :!= rdn_btreated(C.B, xyl)) :& (rdn_bdistance(C.B, xyl) :> C.B.tol)
        if (sum(off)) {
            printf("{txt}note: %g observations have a treated() value that does not match their side of the boundary;\n",
                sum(off[loc]))
            printf("{txt}      they are attached to the network on the side treated() gives\n")
        }
    }
    else TL = rdn_btreated(C.B, xyl)
    n = rows(xyl)
    u = J(n, 1, .)
    v = J(n, 1, .)
    cu = J(n, 1, .)
    cv = J(n, 1, .)
    dep = J(n, 1, .)
    rate = J(n, 1, .)
    for (side = 1; side <= 2; side++) {
        if (side == 1) {
            m = selectindex(TL)
            allowed = C.nodeT :| C.nodeon
        }
        else {
            m = selectindex(!TL)
            allowed = (!C.nodeT) :| C.nodeon
        }
        if (rows(m) == 0) continue
        S = rdn_snap_points(C.N, xyl[m, .], allowed, nodemode)
        u[m] = S.u
        v[m] = S.v
        cu[m] = S.cu
        cv[m] = S.cv
        dep[m] = S.depth
        rate[m] = S.rate
    }
    // off-network cost: snapscale(), or the snapped street's cost per unit length
    depc = dep :* (C.N.snapscale < . ? J(n, 1, C.N.snapscale) : rate)

    bnarc = rdn_bproject(C.B, C.N.xy[bnodes, .])
    if (rows(evalarc)) C.evalarc = evalarc
    else if (neval == 1) C.evalarc = min(bnarc)
    else if (rdn_bclosed(C.B)) {
        // around the ring, without repeating the start
        C.evalarc = mod(min(bnarc) :+ (0::(neval - 1)) :* (C.B.blen / neval), C.B.blen)
    }
    else C.evalarc = min(bnarc) :+ (max(bnarc) - min(bnarc)) :* ((0::(neval - 1)) :/ (neval - 1))
    rdn_snap_eval(C, bnodes, bnarc)

    if (rdn_needs_point_nodes(C.N)) {
        // one-way edges and class-change penalties: put every location and
        // evaluation point on a node, so that paths respect both
        pu = u \ C.eu
        pv = v \ C.ev
        pc = cu \ C.ecu
        tot = pc + (cv \ C.ecv)
        frac = rdn_where(tot :> 0, pc :/ tot, J(rows(tot), 1, 0))
        onends = C.nodeon[pu] :& C.nodeon[pv]
        NP = C.N
        pnode = rdn_point_nodes(NP, pu, pv, frac)
        rdn_edge_frame(NP)
        rdn_compile(NP)
        C.N = NP
        C.nodeon = (rdn_bdistance(C.B, C.N.xy) :<= C.B.tol)
        if (sum(onends)) C.nodeon[pnode[selectindex(onends)]] = J(sum(onends), 1, 1)
        C.nodeT = rdn_btreated(C.B, C.N.xy) :& (!C.nodeon)
        bnodes = selectindex(C.nodeon)
        u = pnode[(1::n)]
        v = u
        cu = J(n, 1, 0)
        cv = J(n, 1, 0)
        C.eu = pnode[((n + 1)::rows(pnode))]
        C.ev = C.eu
        C.ecu = J(rows(C.eu), 1, 0)
        C.ecv = C.ecu
    }
    ndist = .
    nsrc = .
    rdn_net_nearest(C.N, bnodes, ndist, nsrc)
    C.nodedist = ndist
    C.nodesrc = nsrc
    viau = cu + C.nodedist[u]
    viav = cv + C.nodedist[v]
    dist = depc + rowmin((viau, viav))
    xnode = rdn_where(viau :<= viav, C.nodesrc[u], C.nodesrc[v])
    onedge = C.nodeon[u] :& C.nodeon[v] :& (u :!= v)
    dist = rdn_where(onedge, depc, dist)
    xnode = rdn_where(onedge, rdn_where(cu :<= cv, u, v), xnode)
    unreach = (dist :>= 1e299) :| (xnode :== 0)
    C.nunreach = sum(unreach[loc])
    C.Ls = rdn_where(unreach, ., rdn_where(TL, dist, -dist))
    xnode = rdn_where(unreach, 0, xnode)
    safe = rdn_where(xnode :== 0, 1, xnode)
    C.Lu = u
    C.Lv = v
    C.Lcu = cu
    C.Lcv = cv
    C.Ldep = depc
    C.Cu = rdn_where(onedge, u, safe)
    C.Cv = rdn_where(onedge, v, safe)
    C.Ccu = rdn_where(onedge, cu, 0)
    C.Ccv = rdn_where(onedge, cv, 0)
    C.Carc = rdn_where(onedge, rdn_bproject(C.B, xyl), rdn_bproject(C.B, C.N.xy[safe, .]))
    C.Carc = rdn_where(unreach, ., C.Carc)
    C.Cnode = xnode

    C.emasked = rdn_bmasked(C.B, C.evalarc)
    tot = C.ecu + C.ecv
    frac = rdn_where(tot :> 0, C.ecu :/ tot, 0)
    C.exy = C.N.xy[C.eu, .] :* (1 :- frac) + C.N.xy[C.ev, .] :* frac
    C.enode = rdn_where(C.ecu :<= C.ecv, C.eu, C.ev)
    C.enodearc = rdn_bproject(C.B, C.exy)

    // node rows: metric distance from each evaluation point to every node
    Jn = rows(C.evalarc)
    D = rdn_net_distances(C.N, C.eu \ C.ev, rdn_inf())
    C.NR = J(Jn, C.N.n, .)
    for (j = 1; j <= Jn; j++) {
        C.NR[j, .] = colmin((C.ecu[j] :+ D[j, .]) \ (C.ecv[j] :+ D[Jn + j, .]))
    }

    C.xy = xy
    C.loc = loc
    C.firstl = firstl
    C.treated = TL[loc]
    C.s = C.Ls[loc]
    C.crossarc = C.Carc[loc]
    C.depth = dep[loc]
    C.nobs = rows(xy)
    C.neval = Jn
    C.nloc = n
    C.ready = 1
}

// Direct distance for locations on the same segment as evaluation point j.
real colvector rdn_along(struct rdn_coords scalar C, real colvector u, real colvector v,
    real colvector cu, real scalar j)
{
    real colvector out, same, rev
    real scalar u0, v0
    u0 = C.eu[j]
    v0 = C.ev[j]
    out = J(rows(u), 1, rdn_inf())
    if (u0 == v0) return(out)
    same = selectindex((u :== u0) :& (v :== v0))
    rev = selectindex((u :== v0) :& (v :== u0))
    if (rows(same)) out[same] = abs(cu[same] :- C.ecu[j])
    if (rows(rev)) out[rev] = abs(cu[rev] :- C.ecv[j])
    return(out)
}

// Kernel distance d and boundary offset o of every observation at point j
// (missing where infinite).
void rdn_coords_at(struct rdn_coords scalar C, real scalar j, real colvector d, real colvector o)
{
    real colvector D, dl, sg, reach, ol
    D = C.NR[j, .]'
    dl = C.Ldep + rowmin((C.Lcu + D[C.Lu], C.Lcv + D[C.Lv]))
    dl = rowmin((dl, C.Ldep + rdn_along(C, C.Lu, C.Lv, C.Lcu, j)))
    sg = C.Carc :- C.enodearc[j]
    if (rdn_bclosed(C.B)) {
        // closed ring (a zone outline): the shorter way round
        sg = mod(sg :+ C.B.blen / 2, C.B.blen) :- C.B.blen / 2
    }
    sg = sign(sg)
    reach = rowmin((C.Ccu + D[C.Cu], C.Ccv + D[C.Cv]))
    reach = rowmin((reach, rdn_along(C, C.Cu, C.Cv, C.Ccu, j)))
    ol = rdn_where(C.Cnode :== 0, ., sg :* reach)
    dl = rdn_inf2miss(dl)
    ol = rdn_inf2miss(ol)
    d = dl[C.loc]
    o = ol[C.loc]
}

// ---------------------------------------------------------------------------
// Conley (Bartlett) HAC: Psi' K Psi with K_ik = max(0, 1 - d_ik / cutoff)
// ---------------------------------------------------------------------------

real matrix rdn_hac_euclid(real matrix pts, real matrix P, real scalar cutoff)
{
    real matrix out
    real colvector cxs, cys, cnt, ptr, fill, items
    real rowvector S
    real scalar m, k, xmin, ymin, xmax, ymax, cs, span, nx, ny, i, c, reach, ix, iy, t, jj, d, kij
    m = rows(pts)
    k = cols(P)
    out = J(k, k, 0)
    if (m == 0) return(out)
    xmin = min(pts[., 1])
    ymin = min(pts[., 2])
    xmax = max(pts[., 1])
    ymax = max(pts[., 2])
    cs = cutoff
    span = max((xmax - xmin, ymax - ymin))
    if (!(cs > 0) | span / cs > 2000) cs = max((span / 2000, 1e-300))
    nx = floor((xmax - xmin) / cs) + 1
    ny = floor((ymax - ymin) / cs) + 1
    cxs = J(m, 1, .)
    cys = J(m, 1, .)
    cnt = J(nx * ny, 1, 0)
    for (i = 1; i <= m; i++) {
        cxs[i] = min((nx, floor((pts[i, 1] - xmin) / cs) + 1))
        cys[i] = min((ny, floor((pts[i, 2] - ymin) / cs) + 1))
        c = (cys[i] - 1) * nx + cxs[i]
        cnt[c] = cnt[c] + 1
    }
    ptr = 1 \ (runningsum(cnt) :+ 1)
    fill = ptr[(1::(nx * ny))]
    items = J(m, 1, 0)
    for (i = 1; i <= m; i++) {
        c = (cys[i] - 1) * nx + cxs[i]
        items[fill[c]] = i
        fill[c] = fill[c] + 1
    }
    reach = ceil(cutoff / cs)
    for (i = 1; i <= m; i++) {
        S = J(1, k, 0)
        for (ix = max((1, cxs[i] - reach)); ix <= min((nx, cxs[i] + reach)); ix++) {
            for (iy = max((1, cys[i] - reach)); iy <= min((ny, cys[i] + reach)); iy++) {
                c = (iy - 1) * nx + ix
                for (t = ptr[c]; t < ptr[c + 1]; t++) {
                    jj = items[t]
                    d = sqrt((pts[i, 1] - pts[jj, 1])^2 + (pts[i, 2] - pts[jj, 2])^2)
                    kij = 1 - d / cutoff
                    if (kij > 0) S = S + kij :* P[jj, .]
                }
            }
        }
        out = out + P[i, .]' * S
    }
    return(out)
}

// Network version: location i sits on edge (lu, lv) at costs (cu, cv) from its
// ends, after an off-network depth cost dep. Node pairs closer than the cutoff
// come as triplets (na, nb, nd) over compact node ids 1..nn.
real matrix rdn_hac_net(real colvector lu, real colvector lv, real colvector cu, real colvector cv,
    real colvector dep, real colvector na, real colvector nb, real colvector nd, real scalar nn,
    real matrix P, real scalar cutoff)
{
    real matrix out
    real colvector ptr, pb, pd, fill, lptr, lj, lc, best, touched
    real rowvector S
    real scalar m, k, t, a, b, jj, i, ed, ca, dab, r, cand, nt, d, kij
    m = rows(lu)
    k = cols(P)
    out = J(k, k, 0)
    // node pairs by first node
    ptr = J(nn + 1, 1, 0)
    for (t = 1; t <= rows(na); t++) ptr[na[t] + 1] = ptr[na[t] + 1] + 1
    ptr = runningsum(ptr) :+ 1
    pb = J(rows(na), 1, 0)
    pd = J(rows(na), 1, 0)
    fill = ptr[(1::nn)]
    for (t = 1; t <= rows(na); t++) {
        a = na[t]
        pb[fill[a]] = nb[t]
        pd[fill[a]] = nd[t]
        fill[a] = fill[a] + 1
    }
    // locations by endpoint node, with the cost of reaching that node
    lptr = J(nn + 1, 1, 0)
    for (jj = 1; jj <= m; jj++) {
        lptr[lu[jj] + 1] = lptr[lu[jj] + 1] + 1
        lptr[lv[jj] + 1] = lptr[lv[jj] + 1] + 1
    }
    lptr = runningsum(lptr) :+ 1
    lj = J(2 * m, 1, 0)
    lc = J(2 * m, 1, 0)
    fill = lptr[(1::nn)]
    for (jj = 1; jj <= m; jj++) {
        a = lu[jj]
        lj[fill[a]] = jj
        lc[fill[a]] = cu[jj]
        fill[a] = fill[a] + 1
        b = lv[jj]
        lj[fill[b]] = jj
        lc[fill[b]] = cv[jj]
        fill[b] = fill[b] + 1
    }
    best = J(m, 1, rdn_inf())
    touched = J(m, 1, 0)
    for (i = 1; i <= m; i++) {
        nt = 0
        for (ed = 1; ed <= 2; ed++) {
            a = (ed == 1 ? lu[i] : lv[i])
            ca = (ed == 1 ? cu[i] : cv[i])
            for (t = ptr[a]; t < ptr[a + 1]; t++) {
                b = pb[t]
                dab = ca + pd[t]
                for (r = lptr[b]; r < lptr[b + 1]; r++) {
                    jj = lj[r]
                    cand = dab + lc[r]
                    if (best[jj] >= rdn_inf()) {
                        nt++
                        touched[nt] = jj
                    }
                    if (cand < best[jj]) best[jj] = cand
                }
            }
        }
        S = J(1, k, 0)
        for (t = 1; t <= nt; t++) {
            jj = touched[t]
            d = best[jj]
            best[jj] = rdn_inf()
            if (lu[i] == lu[jj] & lv[i] == lv[jj]) d = min((d, abs(cu[i] - cu[jj])))
            d = d + dep[i] + dep[jj]
            if (i == jj) d = 0
            kij = 1 - d / cutoff
            if (kij > 0) S = S + kij :* P[jj, .]
        }
        out = out + P[i, .]' * S
    }
    return(out)
}

// HAC quadratic form for observations idx (rows of Psi), summing rows by location.
real matrix rdn_hac(struct rdn_coords scalar C, real colvector idx, real matrix Psi,
    real scalar cutoff, real scalar euclid)
{
    real colvector lid, locs, pos, eu, ev, nodes, npos, tgt, keep
    real matrix P, out, PR
    real scalar i, L
    lid = C.loc[idx]
    locs = uniqrows(lid)
    L = rows(locs)
    pos = J(max(locs), 1, 0)
    pos[locs] = (1::L)
    P = J(L, cols(Psi), 0)
    for (i = 1; i <= rows(idx); i++) P[pos[lid[i]], .] = P[pos[lid[i]], .] + Psi[i, .]
    if (euclid) {
        out = rdn_hac_euclid(C.xy[C.firstl[locs], .], P, cutoff)
    }
    else {
        eu = C.Lu[locs]
        ev = C.Lv[locs]
        nodes = uniqrows(eu \ ev)
        npos = J(C.N.n, 1, 0)
        npos[nodes] = (1::rows(nodes))
        PR = rdn_net_pairs(C.N, nodes, cutoff)
        if (rows(PR)) {
            tgt = npos[PR[., 2]]
            keep = selectindex(tgt :> 0)
        }
        else {
            tgt = J(0, 1, .)
            keep = J(0, 1, .)
        }
        out = rdn_hac_net(npos[eu], npos[ev], C.Lcu[locs], C.Lcv[locs], C.Ldep[locs],
            (rows(keep) ? PR[keep, 1] : J(0, 1, .)), (rows(keep) ? tgt[keep] : J(0, 1, .)),
            (rows(keep) ? PR[keep, 3] : J(0, 1, .)), rows(nodes), P, cutoff)
    }
    return((out + out') / 2)
}

end

mata:

// ---------------------------------------------------------------------------
// Local polynomial basis and weighted least squares
// ---------------------------------------------------------------------------

// Dimension of the local basis: 2 for a polynomial in (s, o), 1 for a
// polynomial in s alone (used where few boundary crossings leave no variation
// in o). Kept in the external rdn_DIM; 2 unless a fit sets it.
real scalar rdn_dim()
{
    pointer() scalar pd
    pd = findexternal("rdn_DIM")
    if (pd == NULL) return(2)
    if (*pd == 1) return(1)
    return(2)
}

void rdn_set_dim(real scalar dm)
{
    pointer() scalar pd
    pd = findexternal("rdn_DIM")
    if (pd == NULL) pd = crexternal("rdn_DIM")
    *pd = dm
}

// Exponents (a, c) of s^a * o^c for all terms of total degree <= p.
real matrix rdn_poly_exp(real scalar p)
{
    real matrix E
    real scalar k, a
    E = J(0, 2, .)
    if (rdn_dim() == 1) {
        for (a = 0; a <= p; a++) E = E \ (a, 0)
        return(E)
    }
    for (k = 0; k <= p; k++) {
        for (a = k; a >= 0; a--) E = E \ (a, k - a)
    }
    return(E)
}

real scalar rdn_nterms(real scalar p)
{
    if (rdn_dim() == 1) return(p + 1)
    return((p + 1) * (p + 2) / 2)
}

real scalar rdn_exp_index(real scalar p, real scalar a, real scalar c)
{
    real matrix E
    real scalar k
    E = rdn_poly_exp(p)
    for (k = 1; k <= rows(E); k++) {
        if (E[k, 1] == a & E[k, 2] == c) return(k)
    }
    return(.)
}

real matrix rdn_poly_basis(real colvector s, real colvector o, real scalar p)
{
    real matrix E, X
    real colvector sa, oc
    real scalar k
    E = rdn_poly_exp(p)
    X = J(rows(s), rows(E), .)
    for (k = 1; k <= rows(E); k++) {
        sa = (E[k, 1] == 0 ? J(rows(s), 1, 1) : s :^ E[k, 1])
        oc = (E[k, 2] == 0 ? J(rows(o), 1, 1) : o :^ E[k, 2])
        X[., k] = sa :* oc
    }
    return(X)
}

struct rdn_fit {
    real colvector beta, resid, w
    real matrix X, Ginv
    real scalar ok
}

struct rdn_fit scalar rdn_wls(real matrix X, real colvector y, real colvector w)
{
    struct rdn_fit scalar F
    real matrix G, Gi
    F.ok = 0
    G = cross(X, w, X)
    Gi = rdn_inv_checked(G)
    if (hasmissing(Gi)) return(F)
    F.Ginv = Gi
    F.beta = Gi * cross(X, w, y)
    F.resid = y - X * F.beta
    F.X = X
    F.w = w
    F.ok = 1
    return(F)
}

// Weights a_i with c' beta = sum_i a_i y_i.
real colvector rdn_loadings(struct rdn_fit scalar F, real colvector c)
{
    return((F.X * (F.Ginv * c)) :* F.w)
}

real colvector rdn_leverage(struct rdn_fit scalar F)
{
    return(F.w :* rowsum((F.X * F.Ginv) :* F.X))
}

// Per-observation influence psi_i with c' beta - truth ~ sum_i psi_i.
real colvector rdn_influence(struct rdn_fit scalar F, real colvector c, string scalar vce)
{
    real colvector ed
    real scalar n, k
    ed = F.resid
    n = rows(ed)
    k = rows(F.beta)
    if (vce == "hc1") ed = ed :* sqrt(n / max((n - k, 1)))
    else if (vce == "hc2") ed = ed :/ sqrt(rowmax((1 :- rdn_leverage(F), J(n, 1, 1e-8))))
    else if (vce == "hc3") ed = ed :/ rowmax((1 :- rdn_leverage(F), J(n, 1, 1e-8)))
    return(rdn_loadings(F, c) :* ed)
}

// rd2d's (n - 1) / (n - k) cluster degrees-of-freedom factor.
real scalar rdn_cluster_dof(real scalar n, real scalar k, real colvector cl)
{
    if (n > k & rdn_nunique(cl) > 1) return((n - 1) / (n - k))
    return(1)
}

// Joint local design: side-specific polynomials (or pooled slopes) plus
// covariates common to both sides.
real matrix rdn_design(real colvector us, real colvector uo, real colvector Ti, real scalar ordr,
    real matrix Zi, real scalar pooled)
{
    real matrix R, X
    R = rdn_poly_basis(us, uo, ordr)
    if (pooled) {
        if (cols(R) > 1) X = (J(rows(us), 1, 1), Ti, R[., (2..cols(R))])
        else X = (J(rows(us), 1, 1), Ti)
    }
    else X = (R :* Ti, R :* (1 :- Ti))
    if (cols(Zi) > 0 & rows(Zi) > 0) X = (X, Zi)
    return(X)
}

// Columns of X kept for the weighted fit: the largest identifiable sub-basis.
// Columns are taken in order (lower degree first; s before o within a degree);
// a column is dropped when its weighted residual after projecting on the
// columns already kept is below 1e-4 times its norm (as lm() drops aliased
// terms), or when it would give some observation leverage 1 (a coefficient
// determined by that observation alone).
real colvector rdn_identifiable(real matrix X, real colvector w)
{
    real matrix A, Q
    real colvector keep, v, u, lev
    real scalar j, k, rep, n0, n1
    A = X :* sqrt(w)
    Q = J(rows(A), 0, .)
    keep = J(cols(A), 1, 0)
    lev = J(rows(A), 1, 0)
    for (j = 1; j <= cols(A); j++) {
        v = A[., j]
        n0 = sqrt(v' * v)
        if (n0 == 0) continue
        for (rep = 1; rep <= 2; rep++) {
            for (k = 1; k <= cols(Q); k++) v = v - (Q[., k]' * v) :* Q[., k]
        }
        n1 = sqrt(v' * v)
        if (n1 <= 1e-4 * n0) continue
        u = v :/ n1
        if (max(lev + u :* u) > 1 - 1e-8) continue
        keep[j] = 1
        Q = Q, u
        lev = lev + u :* u
    }
    return(keep)
}

// Target vector over the order-`order` basis and the derivative order it picks.
real colvector rdn_deriv_vector(real scalar ordr, real rowvector deriv, real rowvector tang,
    real scalar vord)
{
    real colvector vtgt
    vtgt = J(rdn_nterms(ordr), 1, 0)
    if (cols(tang) == 2) {
        vtgt[rdn_exp_index(ordr, 1, 0)] = tang[1]
        vtgt[rdn_exp_index(ordr, 0, 1)] = tang[2]
        vord = 1
        return(vtgt)
    }
    vtgt[rdn_exp_index(ordr, deriv[1], deriv[2])] = factorial(deriv[1]) * factorial(deriv[2])
    vord = deriv[1] + deriv[2]
    return(vtgt)
}

// Critical value for a uniform band: quantile of max |Z|, Z ~ N(0, corr).
real scalar rdn_uniform_crit(real matrix Sigma, real scalar level, real scalar repp, real scalar seed)
{
    real colvector sd
    real matrix corr, X, Lm, draws
    real rowvector L
    real scalar crit
    transmorphic st
    // one point: the normal quantile itself, no simulation
    if (rows(Sigma) == 1) return(invnormal(0.5 + level / 200))
    sd = sqrt(rowmax((diagonal(Sigma), J(rows(Sigma), 1, 1e-300))))
    corr = Sigma :/ (sd * sd')
    symeigensystem(corr, X, L)
    Lm = X :* sqrt(rowmax((L', J(cols(L), 1, 0))))'
    if (seed < .) {
        st = rngstate()
        rseed(seed)
    }
    draws = rnormal(repp, rows(sd), 0, 1) * Lm'
    if (seed < .) rngstate(st)
    crit = rdn_quantile(rowmax(abs(draws)), level / 100)
    return(crit)
}

// Columns of the rdnet results matrix.
string rowvector rdn_result_names()
{
    return(("eval", "arc", "node_arc", "h0", "h1", "n0", "n1", "tau", "se", "tau_bc", "se_rb",
        "ci_lo", "ci_hi", "pvalue", "status", "h", "cb_lo", "cb_hi", "itt", "fs", "fs_se", "basis"))
}

void rdn_set_rb(real matrix R, real scalar j, real scalar seq, real scalar z)
{
    real scalar tq
    tq = R[j, 10]
    R[j, 11] = seq
    R[j, 12] = tq - z * seq
    R[j, 13] = tq + z * seq
    if (seq < . & seq > 0) R[j, 14] = 2 * normal(-abs(tq / seq))
    else R[j, 14] = .
}

// Local variance of sum(psi): HC, or cluster-robust with rd2d's CR1 factor.
real scalar rdn_local_var(real colvector psi, string scalar vce, real colvector clidx, real scalar k)
{
    real colvector sums
    real scalar g
    if (vce == "conley") return(.)
    if (vce == "cluster") {
        sums = rdn_groupsum(psi, clidx)
        g = rows(sums)
        return(sum(sums :^ 2) * (g > 1 ? g / (g - 1) : 1) * rdn_cluster_dof(rows(psi), k, clidx))
    }
    return(sum(psi :^ 2))
}

real scalar rdn_sqrt0(real scalar v)
{
    if (v >= .) return(.)
    return(sqrt(max((v, 0))))
}

// ---------------------------------------------------------------------------
// rdnet: boundary average treatment effect curve
// ---------------------------------------------------------------------------
// Status codes: 0 ok, 1 masked, 2 insufficient, 3 singular, 4 no bandwidth,
// 5 weak first stage, 6 few clusters (fewer than two on a side).
// Column 22 (basis): 1 = polynomial in (s, o), 2 = in s alone (where few
// boundary crossings leave no variation in o).

real matrix rdn_rdnet(struct rdn_coords scalar C, real colvector y, real colvector D,
    real matrix Z, real colvector om, real colvector cl, real matrix Hin, real scalar p,
    real scalar q, real rowvector deriv, real matrix tang, real scalar kcode, string scalar vce,
    real scalar cutoff, real scalar euclid, real scalar rscale, real scalar rorder,
    real scalar pooled, real scalar level, real scalar cbands, real scalar repp,
    real scalar minobs, real scalar seed, real scalar scaleregul, real scalar scalebiascrct,
    string scalar bwsel, real matrix cov, real scalar crit, real matrix bw)
{
    real matrix H, R, Zi, Xa, Xi, X, Psi, V, PsiC, G
    real colvector T, Tf, s, fin, d, o, base, hside, w, use, i1, i0, idx, zc, hs, wk, m
    real colvector cvec, vtgt, psi, betaa, raux, uni, pos, se_all, ok, keys, f, qc, qj, clu
    real rowvector tj
    real scalar nJ, z, j, conley, hasD, hasZ, hascl, nord, oi, ordr, k, vord, ti, tau, itt, fs
    real scalar basevce_is, status, kk, nc, cc, gg, fail, auxorder, haveband, dimi, gside, gmin
    real colvector dimlist, bwbasis, keep, oexp, ka
    real matrix ex
    real scalar haso, ofirst, oo
    string scalar basevce
    struct rdn_fit scalar F, Fa
    transmorphic store, hstore
    real colvector orders, jumps, hj, hname
    real matrix outp
    rdn_set_dim(2)
    nJ = C.neval
    hasD = (rows(D) > 0)
    hasZ = (cols(Z) > 0 & rows(Z) > 0)
    hascl = (rows(cl) > 0)
    conley = (vce == "conley")
    basevce = (vce == "hc1" | vce == "hc2" | vce == "hc3" ? vce : "hc0")
    auxorder = (rorder < . ? rorder : q)
    bw = J(0, 0, .)
    bwbasis = J(nJ, 1, .)
    if (rows(Hin) == 0) {
        bw = rdn_rdbwnet(C, y, p, deriv, tang, Z, kcode, bwsel, basevce, cl, D, om, minobs,
            scaleregul, scalebiascrct)
        H = bw[., (21, 22)]
        bwbasis = bw[., 24]
    }
    else H = Hin
    z = invnormal(0.5 + level / 200)
    T = C.treated
    Tf = T
    s = C.s
    fin = (s :< .) :& (y :< .)
    if (hasD) fin = fin :& (D :< .)
    if (hasZ) fin = fin :& (rowmissing(Z) :== 0)
    R = J(nJ, 22, .)
    R[., 1] = (1::nJ)
    R[., 2] = C.evalarc
    R[., 3] = C.enodearc
    R[., 4] = H[., 1]
    R[., 5] = H[., 2]
    R[., 6] = J(nJ, 1, 0)
    R[., 7] = J(nJ, 1, 0)
    R[., 15] = J(nJ, 1, 0)
    store = asarray_create("real")
    hstore = asarray_create("real")
    hj = J(0, 1, .)
    hname = J(0, 1, .)
    orders = (p == q ? p : (p \ q))
    gmin = .  // clusters per side at the estimated points (fewest side)
    for (j = 1; j <= nJ; j++) {
        if (C.emasked[j]) {
            R[j, 15] = 1
            continue
        }
        if (!(H[j, 1] < . & H[j, 2] < . & H[j, 1] > 0 & H[j, 2] > 0)) {
            // the selector's reason (insufficient, singular) when it has one
            if (rows(bw)) R[j, 15] = (bw[j, 3] > 0 ? bw[j, 3] : 4)
            else R[j, 15] = 4
            continue
        }
        d = .
        o = .
        rdn_coords_at(C, j, d, o)
        base = (d :< .) :& (o :< .) :& fin :& (om :> 0)
        hside = rdn_where(T, H[j, 2], H[j, 1])
        // each side's kernel weights scaled by 1 / h_side^2, as in rd2d: this matters
        // for the coefficients shared across sides (covariates, pooled slopes)
        w = rdn_kernel(d :/ hside, kcode) :/ (hside :^ 2) :* editmissing(om, 0)
        use = (w :> 0) :& base
        i1 = selectindex(use :& T)
        i0 = selectindex(use :& (!T))
        R[j, 7] = rows(i1)
        R[j, 6] = rows(i0)
        // count distinct locations: repeated measurements at one site add no design points
        if (min((rows(rdn_first_per_value(C.loc[i1])), rows(rdn_first_per_value(C.loc[i0])))) < minobs) {
            R[j, 15] = 2
            continue
        }
        idx = i1 \ i0
        Zi = J(rows(idx), 0, .)
        zc = J(0, 1, .)
        if (hasZ) {
            zc = selectindex((colmax(Z[idx, .]) - colmin(Z[idx, .]))' :> 1e-12)
            if (rows(zc)) Zi = Z[idx, zc']
        }
        if (hascl) {
            gside = min((rdn_nunique(cl[i1]), rdn_nunique(cl[i0])))
            if (gside >= 2) gmin = min((gmin, gside))
            if (gside < 2) {
                R[j, 15] = 6
                continue
            }
        }
        // the (s, o) basis, keeping its identifiable columns (rdn_identifiable)
        rdn_set_dim(2)
        haso = 0
        if (conley) {
            hs = hside :* rscale
            wk = rdn_kernel(d :/ hs, kcode) :/ (hs :^ 2) :* editmissing(om, 0)
            m = selectindex((wk :> 0) :& base)
            Xa = rdn_design(s[m] :/ hs[m], o[m] :/ hs[m], Tf[m], auxorder,
                (rows(zc) ? Z[m, zc'] : J(rows(m), 0, .)), 0)
            Xi = rdn_design(s[idx] :/ hs[idx], o[idx] :/ hs[idx], Tf[idx], auxorder, Zi, 0)
            ka = rdn_identifiable(Xa, wk[m])
            Xa = select(Xa, ka')
            Xi = select(Xi, ka')
        }
        status = 0
        outp = J(2, 5, .)     // rows: order p, order q; cols: tau, se, k, itt, fs
        psip = J(0, 1, .)
        psiq = J(0, 1, .)
        psifs = J(0, 1, .)
        for (oi = 1; oi <= rows(orders); oi++) {
            ordr = orders[oi]
            X = rdn_design(s[idx] :/ hside[idx], o[idx] :/ hside[idx], Tf[idx], ordr, Zi, pooled)
            cvec = J(cols(X), 1, 0)
            if (pooled) cvec[2] = 1
            else {
                vord = 0
                tj = (rows(tang) ? tang[j, .] : J(1, 0, .))
                vtgt = rdn_deriv_vector(ordr, deriv, tj, vord)
                kk = rows(vtgt)
                cvec[(1::kk)] = vtgt :/ (H[j, 2] ^ vord)
                cvec[((kk + 1)::(2 * kk))] = -vtgt :/ (H[j, 1] ^ vord)
            }
            // keep the identifiable columns (few boundary crossings leave too
            // little variation in the offset o for its higher powers, or any)
            keep = rdn_identifiable(X, w[idx])
            if (any(select(cvec, !keep) :!= 0)) {
                status = 3
                break
            }
            if (ordr == p) {
                ex = rdn_poly_exp(ordr)
                oexp = (pooled ? ex[(2::rows(ex)), 2] : (ex[., 2] \ ex[., 2]))
                ofirst = (pooled ? 2 : 0)
                for (oo = 1; oo <= rows(oexp); oo++) {
                    if (oexp[oo] > 0 & keep[ofirst + oo]) haso = 1
                }
            }
            X = select(X, keep')
            cvec = select(cvec, keep)
            jumps = J(0, 1, .)
            psi1 = J(0, 1, .)
            psi2 = J(0, 1, .)
            for (ti = 1; ti <= (hasD ? 2 : 1); ti++) {
                if (ti == 1) F = rdn_wls(X, y[idx], w[idx])
                else F = rdn_wls(X, D[idx], w[idx])
                if (!F.ok) {
                    status = 3
                    break
                }
                jumps = jumps \ (cvec' * F.beta)
                if (!conley) {
                    psi = rdn_influence(F, cvec, basevce)
                }
                else {
                    if (ti == 1) Fa = rdn_wls(Xa, y[m], wk[m])
                    else Fa = rdn_wls(Xa, D[m], wk[m])
                    if (!Fa.ok) {
                        status = 3
                        break
                    }
                    if (ti == 1) raux = y[idx] - Xi * Fa.beta
                    else raux = D[idx] - Xi * Fa.beta
                    psi = rdn_loadings(F, cvec) :* raux
                }
                if (ti == 1) psi1 = psi
                else psi2 = psi
            }
            if (status) break
            if (!hasD) {
                tau = jumps[1]
                psi = psi1
            }
            else {
                itt = jumps[1]
                fs = jumps[2]
                if (abs(fs) < 1e-8) {
                    status = 5
                    break
                }
                tau = itt / fs
                if (ordr == p) {
                    outp[1, 4] = itt
                    outp[1, 5] = fs
                    psifs = psi2
                }
                psi = (psi1 - tau :* psi2) :/ fs
            }
            k = cols(X)
            if (ordr == p) {
                outp[1, 1] = tau
                outp[1, 2] = rdn_sqrt0(rdn_local_var(psi, vce, (hascl ? cl[idx] : J(0, 1, .)), k))
                outp[1, 3] = k
                psip = psi
            }
            if (ordr == q) {
                outp[2, 1] = tau
                outp[2, 2] = rdn_sqrt0(rdn_local_var(psi, vce, (hascl ? cl[idx] : J(0, 1, .)), k))
                outp[2, 3] = k
                psiq = psi
            }
        }
        if (status) {
            R[j, 15] = status
            continue
        }
        // basis 2: no term in the offset o is identified (one crossing): s alone
        R[j, 22] = (haso ? 1 : 2)
        if (hasD) {
            R[j, 19] = outp[1, 4]
            R[j, 20] = outp[1, 5]
            R[j, 21] = rdn_sqrt0(rdn_local_var(psifs, vce, (hascl ? cl[idx] : J(0, 1, .)), outp[1, 3]))
            if (conley) {
                hj = hj \ j
                hname = hname \ 3
                asarray(hstore, rows(hj), (idx, psifs))
            }
        }
        if (conley) {
            hj = hj \ j \ j
            hname = hname \ 1 \ 2
            asarray(hstore, rows(hj) - 1, (idx, psip))
            asarray(hstore, rows(hj), (idx, psiq))
        }
        R[j, 8] = outp[1, 1]
        R[j, 9] = outp[1, 2]
        R[j, 10] = outp[2, 1]
        rdn_set_rb(R, j, outp[2, 2], z)
        asarray(store, j, (idx, psiq, J(rows(idx), 1, outp[2, 3])))
    }

    crit = .
    cov = J(0, 0, .)
    if (conley & rows(hj)) {
        uni = J(0, 1, .)
        for (cc = 1; cc <= rows(hj); cc++) {
            PsiC = asarray(hstore, cc)
            uni = uni \ PsiC[., 1]
        }
        uni = uniqrows(uni)
        pos = J(C.nobs, 1, 0)
        pos[uni] = (1::rows(uni))
        Psi = J(rows(uni), rows(hj), 0)
        for (cc = 1; cc <= rows(hj); cc++) {
            PsiC = asarray(hstore, cc)
            Psi[pos[PsiC[., 1]], cc] = PsiC[., 2]
        }
        V = rdn_hac(C, uni, Psi, cutoff, euclid)
        se_all = sqrt(rowmax((diagonal(V), J(rows(V), 1, 0))))
        for (cc = 1; cc <= rows(hj); cc++) {
            if (hname[cc] == 2) rdn_set_rb(R, hj[cc], se_all[cc], z)
            else if (hname[cc] == 1) R[hj[cc], 9] = se_all[cc]
            else R[hj[cc], 21] = se_all[cc]
        }
        qc = selectindex(hname :== 2)
        qj = hj[qc]
        cov = J(nJ, nJ, .)
        cov[qj, qj] = V[qc, qc]
    }
    else {
        keys = selectindex(R[., 15] :== 0)
        if (rows(keys)) {
            uni = J(0, 1, .)
            for (cc = 1; cc <= rows(keys); cc++) {
                PsiC = asarray(store, keys[cc])
                uni = uni \ PsiC[., 1]
            }
            uni = uniqrows(uni)
            pos = J(C.nobs, 1, 0)
            pos[uni] = (1::rows(uni))
            Psi = J(rows(uni), rows(keys), 0)
            f = J(rows(keys), 1, 1)
            for (cc = 1; cc <= rows(keys); cc++) {
                PsiC = asarray(store, keys[cc])
                Psi[pos[PsiC[., 1]], cc] = PsiC[., 2]
                if (vce == "cluster") {
                    clu = cl[PsiC[., 1]]
                    gg = rdn_nunique(clu)
                    f[cc] = sqrt((gg > 1 ? gg / (gg - 1) : 1) * rdn_cluster_dof(rows(PsiC), PsiC[1, 3], clu))
                }
            }
            if (vce == "cluster") {
                G = rdn_groupsum(Psi, cl[uni])
                V = (G' * G) :* (f * f')
            }
            else V = Psi' * Psi
            cov = J(nJ, nJ, .)
            cov[keys, keys] = V
        }
    }
    R[., 16] = rdn_where(rdn_isclose(R[., 4], R[., 5]), R[., 4], .)
    if (rows(cov)) {
        ok = selectindex(diagonal(cov) :< .)
        if (rows(ok)) {
            crit = rdn_uniform_crit(cov[ok, ok], level, repp, seed)
            if (hascl & gmin < 10) {
                printf("{txt}note: cluster-robust standard errors rest on few clusters (as few as %g per side within\n", gmin)
                printf("{txt}      the bandwidth); with fewer than about 10 per side they are too small. Consider vce(conley #).\n")
            }
            if (cbands) {
                R[., 17] = R[., 10] :- crit :* R[., 11]
                R[., 18] = R[., 10] :+ crit :* R[., 11]
            }
        }
    }
    return(R)
}

end

mata:

// ---------------------------------------------------------------------------
// Bandwidth selection (port of rd2d's rdbw2d direct plug-in)
// ---------------------------------------------------------------------------
// Result columns: 1 eval, 2 arc, 3 status, 4 dn, 5 bn0, 6 bn1, 7 n0, 8 n1,
// 9 hmin0, 10 hmin1, 11 hmax0, 12 hmax1, 13 V0, 14 B0, 15 Reg10, 16 Reg20,
// 17 V1, 18 B1, 19 Reg11, 20 Reg21, 21 h0, 22 h1, 23 h, 24 basis (1 = s,o;
// 2 = s alone, where few boundary crossings leave no variation in o).

struct rdn_side {
    real colvector s, o, d, y, om, cl
    real scalar kcode, hascl
}

string rowvector rdn_bw_names()
{
    return(("eval", "arc", "status", "dn", "bn0", "bn1", "n0", "n1", "hmin0", "hmin1", "hmax0",
        "hmax1", "V0", "B0", "Reg10", "Reg20", "V1", "B1", "Reg11", "Reg21", "h0", "h1", "h", "basis"))
}

// rd2d's bivariate normal-reference rule of thumb.
real scalar rdn_rot(real colvector s, real colvector o, real scalar kcode, real scalar neff)
{
    real rowvector kc
    real matrix S, inv
    real scalar dt, sqdt, tc
    if (rdn_dim() == 1) return(1.06 * sqrt(variance(s)) * neff^(-1 / 5))
    kc = rdn_rot_const(kcode)
    S = variance((s, o))
    inv = rdn_pinv(S, 1e-20)
    dt = det(S)
    if (dt < . & dt > 0) sqdt = sqrt(dt)
    else sqdt = sqrt(max((rowmax(symeigenvalues(S)) * rowmin(symeigenvalues(S)), 1e-300)))
    tc = 1 / (2^4 * pi() * sqdt) * (2 * trace(inv * inv) + trace(inv)^2)
    return(((2 * kc[2]) / (neff * kc[1] * tc))^(1 / 6))
}

// Kernel weights K(d/h)/h^2 * om, the basis in (s/h, o/h) and the outcome.
void rdn_side_design(struct rdn_side scalar sd, real scalar h, real scalar ordr, real colvector yalt,
    real matrix R, real colvector w, real colvector yy, real colvector clk)
{
    real colvector wa, keep
    wa = rdn_kernel(sd.d :/ h, sd.kcode) :/ (h^2) :* sd.om
    keep = selectindex(wa :> 0)
    w = wa[keep]
    R = rdn_poly_basis(sd.s[keep] :/ h, sd.o[keep] :/ h, ordr)
    if (rows(yalt)) yy = yalt[keep]
    else yy = sd.y[keep]
    if (sd.hascl) clk = sd.cl[keep]
    else clk = J(0, 1, .)
}

real matrix rdn_bw_sigma(real matrix WR, real colvector resd, real scalar h, real colvector clk)
{
    real matrix scores, summed
    real scalar n, k, g, adj
    scores = WR :* resd
    if (rows(clk) == 0) return(cross(scores, scores) :* h^2)
    summed = rdn_groupsum(scores, clk)
    n = rows(clk)
    k = cols(scores)
    g = rows(summed)
    adj = (g > 1 & n > k ? ((n - 1) / (n - k)) * (g / (g - 1)) : 1)
    return(cross(summed, summed) :* h^2 :* adj)
}

real colvector rdn_bw_adjust(real colvector resd, real matrix R, real colvector w, real matrix invG,
    string scalar vce)
{
    real colvector hii
    real scalar n, k
    n = rows(R)
    k = cols(R)
    if (vce == "hc1" & n > k) return(resd :* sqrt(n / (n - k)))
    if (vce == "hc2" | vce == "hc3") {
        hii = w :* rowsum((R * invG) :* R)
        if (vce == "hc2") return(resd :/ sqrt(rowmax((1 :- hii, J(n, 1, 1e-8)))))
        return(resd :/ rowmax((1 :- hii, J(n, 1, 1e-8))))
    }
    return(resd)
}

// Raw-scale coefficients (and their scaled covariance constant).
void rdn_side_lm(struct rdn_side scalar sd, real scalar h, real scalar ordr, string scalar vce,
    real scalar varr, real colvector yalt, real colvector beta, real matrix covm)
{
    real matrix R, invG, E
    real colvector w, yy, clk, betas, resd
    R = .
    w = .
    yy = .
    clk = .
    rdn_side_design(sd, h, ordr, yalt, R, w, yy, clk)
    invG = rdn_inv_gram(cross(R, w, R))
    betas = invG * cross(R, w, yy)
    E = rdn_poly_exp(ordr)
    beta = betas :/ (h :^ (E[., 1] + E[., 2]))
    covm = J(0, 0, .)
    if (varr) {
        resd = rdn_bw_adjust(yy - R * betas, R, w, invG, vce)
        covm = invG * rdn_bw_sigma(R :* w, resd, h, clk) * invG
    }
}

// Projection of the degree p+1 terms onto the target vec (rd2d's vec_q).
real colvector rdn_side_coeff(struct rdn_side scalar sd, real scalar dn, real scalar p, real colvector vtgt)
{
    real matrix R, Rp, S, invG
    real colvector w, yy, clk
    real scalar k
    R = .
    w = .
    yy = .
    clk = .
    rdn_side_design(sd, dn, p + 1, J(0, 1, .), R, w, yy, clk)
    k = rdn_nterms(p)
    Rp = R[., (1..k)]
    S = R[., ((k + 1)..cols(R))]
    invG = rdn_inv_gram(cross(Rp, w, Rp))
    return(J(k, 1, 0) \ (vtgt' * invG * cross(Rp, w, S))')
}

// (V, B, Reg1, Reg2); bn2 missing skips the higher-order bias term.
real rowvector rdn_side_consts(struct rdn_side scalar sd, real scalar p, real colvector vtgt, real scalar dn,
    real scalar bn1, real scalar bn2, string scalar vce, real colvector yalt)
{
    real matrix R, Rp, S, Tm, WR, invG, cov1
    real colvector w, yy, clk, vecq, betas, resd, beta1, beta2, vect
    real scalar kp, kp1, V, B, reg1, reg2
    R = .
    w = .
    yy = .
    clk = .
    rdn_side_design(sd, dn, p + (bn2 >= . ? 1 : 2), yalt, R, w, yy, clk)
    kp = rdn_nterms(p)
    kp1 = rdn_nterms(p + 1)
    Rp = R[., (1..kp)]
    S = R[., ((kp + 1)..kp1)]
    WR = Rp :* w
    invG = rdn_inv_gram(cross(WR, Rp))
    vecq = J(kp, 1, 0) \ (vtgt' * invG * cross(WR, S))'
    betas = invG * cross(WR, yy)
    resd = rdn_bw_adjust(yy - Rp * betas, Rp, w, invG, vce)
    V = vtgt' * invG * rdn_bw_sigma(WR, resd, dn, clk) * invG * vtgt
    beta1 = .
    cov1 = .
    rdn_side_lm(sd, bn1, p + 1, vce, 1, yalt, beta1, cov1)
    B = vecq' * beta1
    reg1 = (vecq' * cov1 * vecq) / bn1^(2 + 2 * (p + 1))
    // a quadratic form in a covariance matrix: negative values are rounding
    // (a degenerate bias direction), and their sign must not decide the fit
    if (reg1 < 0) reg1 = 0
    reg2 = .
    if (bn2 < .) {
        Tm = R[., ((kp1 + 1)..cols(R))]
        vect = J(kp1, 1, 0) \ (vtgt' * invG * cross(WR, Tm))'
        beta2 = .
        cov1 = .
        rdn_side_lm(sd, bn2, p + 2, vce, 0, yalt, beta2, cov1)
        reg2 = dn * (vect' * beta2)
    }
    return((V, B, reg1, reg2))
}

real scalar rdn_consts_ok(real rowvector c)
{
    return(c[1] < . & c[2] < . & c[3] < . & c[1] > 0 & c[3] >= 0)
}

// rd2d's shrinkage from MSE-optimal to coverage-error-optimal bandwidth.
real scalar rdn_cer_factor(real scalar n, real scalar p)
{
    real scalar nn
    nn = max((n, 1))
    return(nn^(1 / (2 * p + 4) - 1 / (p + 4)))
}

// Covariate coefficients common to both sides, partialling out side polynomials.
real colvector rdn_covariate_gamma(struct rdn_side scalar s0, struct rdn_side scalar s1, real matrix Z0,
    real matrix Z1, real rowvector dns, real scalar p, real colvector y0, real colvector y1)
{
    real matrix A, Bm, R, Zk, invG, uz, RtW
    real colvector w, k, W, yk
    real scalar t
    struct rdn_side scalar sd
    A = 0
    Bm = 0
    for (t = 1; t <= 2; t++) {
        if (t == 1) sd = s0
        else sd = s1
        w = rdn_kernel(sd.d :/ dns[t], sd.kcode) :/ (dns[t]^2) :* sd.om
        k = selectindex(w :> 0)
        R = rdn_poly_basis(sd.s[k] :/ dns[t], sd.o[k] :/ dns[t], p)
        W = w[k]
        if (t == 1) {
            Zk = Z0[k, .]
            yk = y0[k]
        }
        else {
            Zk = Z1[k, .]
            yk = y1[k]
        }
        RtW = (R :* W)'
        invG = rdn_pinv(RtW * R, 1e-12)
        uz = RtW * Zk
        A = A :+ (cross(Zk, W, Zk) - uz' * invG * uz)
        Bm = Bm :+ (cross(Zk, W, yk) - uz' * invG * (RtW * yk))
    }
    return(rdn_lstsq(A, Bm))
}

real matrix rdn_rdbwnet(struct rdn_coords scalar C, real colvector y, real scalar p, real rowvector deriv,
    real matrix tang, real matrix Z, real scalar kcode, string scalar bwsel, string scalar vce,
    real colvector cl, real colvector D, real colvector om, real scalar bwcheck,
    real scalar scaleregul, real scalar scalebiascrct)
{
    struct rdn_side scalar s0, s1, sd
    real matrix Out, Zs0, Zs1, g
    real colvector T, s, fin, d, o, ok, m0, m1, edv, vtgt, zc, gy, gd, D0, D1, ybw0, ybw1, vq, Bc0, Bc1
    real colvector num, den, h, good, rg, r2, okc, f0, f1
    real rowvector bnd0, bnd1, dns, thr, bc, hc0, hc1
    real colvector bnv
    real scalar nJ, j, neff, neff0, neff1, vord, nt, i10, i01, expo_b, dn, itt, fs, t, sgn, bt, expo
    real scalar hasD, hasZ, numc, denc, cer, isi, dimi
    real colvector dimlist
    string scalar base
    nJ = C.neval
    hasD = (rows(D) > 0)
    hasZ = (cols(Z) > 0 & rows(Z) > 0)
    vv = vce
    if (vv != "hc0" & vv != "hc1" & vv != "hc2" & vv != "hc3") vv = "hc1"
    T = C.treated
    s = C.s
    fin = (s :< .)
    // effective sample sizes count distinct locations, as in rd2d
    neff = rdn_nunique(select(C.loc, fin))
    neff0 = rdn_nunique(select(C.loc, fin :& (!T)))
    neff1 = rdn_nunique(select(C.loc, fin :& T))
    vord = (rows(tang) ? 1 : deriv[1] + deriv[2])
    expo_b = 1 / (2 * p + 6)
    if (deriv[2] > 0 | rows(tang)) dimlist = 2
    else dimlist = (2 \ 1)
    Out = J(nJ, 24, .)
    Out[., 1] = (1::nJ)
    Out[., 2] = C.evalarc
    Out[., 3] = J(nJ, 1, 0)
    for (j = 1; j <= nJ; j++) {
        if (C.emasked[j]) {
            Out[j, 3] = 1
            continue
        }
        // try the (s, o) basis first; where it is singular, fit in s alone
        for (dimi = 1; dimi <= rows(dimlist); dimi++) {
        rdn_set_dim(dimlist[dimi])
        Out[j, 3] = 0
        nt = rdn_nterms(p)
        edv = J(nt, 1, 0)
        edv[rdn_exp_index(p, deriv[1], deriv[2])] = factorial(deriv[1]) * factorial(deriv[2])
        i10 = rdn_exp_index(p, 1, 0)
        i01 = rdn_exp_index(p, 0, 1)
        vtgt = edv
        if (rows(tang)) {
            vtgt = J(nt, 1, 0)
            vtgt[i10] = tang[j, 1]
            vtgt[i01] = tang[j, 2]
        }
        d = .
        o = .
        rdn_coords_at(C, j, d, o)
        ok = (d :< .) :& (o :< .) :& (s :< .) :& (y :< .) :& (editmissing(om, 0) :> 0)
        if (hasD) ok = ok :& (D :< .)
        if (hasZ) ok = ok :& (rowmissing(Z) :== 0)
        m0 = selectindex(ok :& (!T))
        m1 = selectindex(ok :& T)
        // bwcheck counts distinct locations (repeated measurements count once), as
        // rdrobust and rd2d count unique observations
        f0 = rdn_first_per_value(C.loc[m0])
        f1 = rdn_first_per_value(C.loc[m1])
        if (min((rows(f0), rows(f1))) < bwcheck) {
            Out[j, 3] = 2
            continue
        }
        s0.s = s[m0]; s0.o = o[m0]; s0.d = d[m0]; s0.y = y[m0]; s0.om = om[m0]; s0.kcode = kcode
        s1.s = s[m1]; s1.o = o[m1]; s1.d = d[m1]; s1.y = y[m1]; s1.om = om[m1]; s1.kcode = kcode
        s0.hascl = (rows(cl) > 0)
        s1.hascl = (rows(cl) > 0)
        if (s0.hascl) {
            s0.cl = cl[m0]
            s1.cl = cl[m1]
        }
        okc = selectindex(ok)
        dn = rdn_rot(s[okc], o[okc], kcode, neff)
        // nudged up so the bwcheck-th observation gets positive kernel weight
        bnd0 = (rdn_kth(s0.d[f0], bwcheck) * (1 + 1e-9), max(s0.d))
        bnd1 = (rdn_kth(s1.d[f1], bwcheck) * (1 + 1e-9), max(s1.d))
        dns = (min((max((dn, bnd0[1])), bnd0[2])), min((max((dn, bnd1[1])), bnd1[2])))
        thr = (rdn_median(s0.d), rdn_median(s1.d))
        if (hasD) {
            D0 = D[m0]
            D1 = D[m1]
        }
        if (hasZ) {
            rg = (colmax(Z[okc, .]) - colmin(Z[okc, .]))'
            zc = selectindex(rg :> 1e-12)
            if (rows(zc)) {
                Zs0 = Z[m0, zc']
                Zs1 = Z[m1, zc']
                gy = rdn_covariate_gamma(s0, s1, Zs0, Zs1, dns, p, s0.y, s1.y)
                s0.y = s0.y - Zs0 * gy
                s1.y = s1.y - Zs1 * gy
                if (hasD) {
                    gd = rdn_covariate_gamma(s0, s1, Zs0, Zs1, dns, p, D0, D1)
                    D0 = D0 - Zs0 * gd
                    D1 = D1 - Zs1 * gd
                }
            }
        }
        ybw0 = J(0, 1, .)
        ybw1 = J(0, 1, .)
        if (hasD) {
            itt = 0
            fs = 0
            bnv = .
            g = .
            rdn_side_lm(s0, dns[1], p, vv, 0, J(0, 1, .), bnv, g)
            itt = itt - vtgt' * bnv
            rdn_side_lm(s1, dns[2], p, vv, 0, J(0, 1, .), bnv, g)
            itt = itt + vtgt' * bnv
            rdn_side_lm(s0, dns[1], p, vv, 0, D0, bnv, g)
            fs = fs - vtgt' * bnv
            rdn_side_lm(s1, dns[2], p, vv, 0, D1, bnv, g)
            fs = fs + vtgt' * bnv
            if (fs < . & abs(fs) > 1e-8) {
                ybw0 = s0.y :/ fs - itt / fs^2 :* D0
                ybw1 = s1.y :/ fs - itt / fs^2 :* D1
            }
        }
        isi = 0
        for (t = 1; t <= 2; t++) {
            if (t == 1) sd = s0
            else sd = s1
            vq = rdn_side_coeff(sd, dns[t], p, vtgt)
            bc = rdn_side_consts(sd, p + 1, vq, dns[t], thr[t], ., vv, (t == 1 ? ybw0 : ybw1))
            if (!rdn_consts_ok(bc)) {
                isi = 1
                break
            }
            bt = ((2 + 2 * (p + 1)) * bc[1] / (2 * (bc[2]^2 + scaleregul * bc[3])))^expo_b
            if (t == 1) bt = min((max((bt, bnd0[1])), bnd0[2]))
            else bt = min((max((bt, bnd1[1])), bnd1[2]))
            Out[j, 4 + t] = bt
            if (t == 1) {
                hc0 = rdn_side_consts(sd, p, vtgt, dns[t], bt, thr[t], vv, ybw0)
                if (!rdn_consts_ok(hc0)) {
                    isi = 1
                    break
                }
            }
            else {
                hc1 = rdn_side_consts(sd, p, vtgt, dns[t], bt, thr[t], vv, ybw1)
                if (!rdn_consts_ok(hc1)) {
                    isi = 1
                    break
                }
            }
        }
        if (isi) {
            Out[j, 3] = 3
            Out[j, (5, 6)] = (., .)
            continue
        }
        Out[j, 4] = dn
        Out[j, 7] = rows(m0)
        Out[j, 8] = rows(m1)
        Out[j, (9, 10, 11, 12)] = (bnd0[1], bnd1[1], bnd0[2], bnd1[2])
        Out[j, (13, 14, 15, 16)] = hc0
        Out[j, (17, 18, 19, 20)] = hc1
        Out[j, 24] = (dimlist[dimi] == 2 ? 1 : 2)
        break
        }
        rdn_set_dim(2)
    }
    rdn_set_dim(2)
    good = selectindex(Out[., 3] :== 0)
    if (rows(good) == 0) return(Out)
    g = Out[good, .]
    Bc0 = g[., 14] + scalebiascrct :* editmissing(g[., 16], 0)
    Bc1 = g[., 18] + scalebiascrct :* editmissing(g[., 20], 0)
    expo = 1 / (2 * p + 4)
    base = bwsel
    if (bwsel == "cerrd") base = "mserd"
    else if (bwsel == "certwo") base = "msetwo"
    else if (bwsel == "icerrd") base = "imserd"
    else if (bwsel == "icertwo") base = "imsetwo"
    cer = (base != bwsel)
    numc = 2 + 2 * vord
    denc = 2 * p + 2 - 2 * vord
    if (base == "mserd" | base == "imserd") {
        num = numc :* (g[., 13] + g[., 17])
        den = denc :* ((Bc0 - Bc1) :^ 2 + scaleregul :* (g[., 15] + g[., 19]))
        if (base == "imserd") h = J(rows(g), 1, (mean(num) / mean(den))^expo)
        else h = (num :/ den) :^ expo
        if (cer) h = h :* rdn_cer_factor(neff, p)
        h = rowmin((rowmax((h, rowmax((g[., 9], g[., 10])))), rowmax((g[., 11], g[., 12]))))
        Out[good, 21] = h
        Out[good, 22] = h
    }
    else {
        for (t = 0; t <= 1; t++) {
            num = numc :* g[., 13 + 4 * t]
            den = denc :* ((t == 0 ? Bc0 : Bc1) :^ 2 + scaleregul :* g[., 15 + 4 * t])
            if (base == "imsetwo") h = J(rows(g), 1, (mean(num) / mean(den))^expo)
            else h = (num :/ den) :^ expo
            if (cer) h = h :* rdn_cer_factor((t == 0 ? neff0 : neff1), p)
            h = rowmin((rowmax((h, g[., 9 + t])), g[., 11 + t]))
            Out[good, 21 + t] = h
        }
    }
    Out[., 23] = rdn_where(rdn_isclose(Out[., 21], Out[., 22]), Out[., 21], .)
    return(Out)
}

end

mata:

// ---------------------------------------------------------------------------
// Density (sorting) test, ported from rddensity and adapted to networks
// ---------------------------------------------------------------------------

// Local polynomial CDF fit. Columns [1_L, 1_R, u_L, u_R, u_L^2, u_R^2, ...];
// the density on each side is the coefficient on u (positions 3 and 4).
// Returns 0 when one side has too few observations in the window.
real scalar rdn_cjm(real colvector x0, real colvector om0, real scalar hl, real scalar hr,
    real scalar p, real scalar kcode, real colvector b, real matrix V, real scalar nl, real scalar nr)
{
    real colvector o, x, om, F, win, xw, Fw, omw, left, u, K, keep, Hp, col
    real matrix Xp, XpW, Sinv, tl, L, Hinv
    real scalar n, denom, m, k, c
    n = rows(x0)
    o = order((x0, (1::n)), (1, 2))
    x = x0[o]
    om = om0[o]
    denom = sum(om) - mean(om)
    F = (runningsum(om) - om) :/ denom
    win = selectindex((x :>= -hl) :& (x :<= hr))
    xw = x[win]
    Fw = F[win]
    omw = om[win]
    left = (xw :< 0)
    u = rdn_where(left, xw :/ hl, xw :/ hr)
    K = rdn_where(left, rdn_kernel(abs(u), kcode) :/ hl, rdn_kernel(abs(u), kcode) :/ hr)
    keep = selectindex(K :> 0)
    xw = xw[keep]
    Fw = Fw[keep]
    omw = omw[keep]
    left = left[keep]
    u = u[keep]
    K = K[keep]
    nl = sum(left)
    nr = rows(left) - nl
    if (min((nl, nr)) <= p + 1) return(0)
    m = rows(xw)
    Xp = J(m, 2 * p + 2, 0)
    Hp = J(2 * p + 2, 1, 0)
    for (k = 0; k <= p; k++) {
        col = (k == 0 ? J(m, 1, 1) : u :^ k)
        Xp[., 2 * k + 1] = left :* col
        Xp[., 2 * k + 2] = (1 :- left) :* col
        Hp[2 * k + 1] = hl^k
        Hp[2 * k + 2] = hr^k
    }
    XpW = Xp :* K
    Sinv = luinv(cross(XpW, Xp))
    b = (Sinv * cross(XpW, Fw)) :/ Hp
    // sum over observations k with x_k >= x_i
    tl = J(m, cols(XpW), 0)
    for (c = 1; c <= cols(XpW); c++) {
        col = runningsum(XpW[(m::1), c])
        tl[., c] = col[(m::1)]
    }
    L = tl :* (omw :/ denom)
    Hinv = diag(1 :/ Hp)
    V = Hinv * Sinv * cross(L, L) * Sinv * Hinv
    return(1)
}

// MSE-optimal bandwidth for one side's density.
real scalar rdn_side_bandwidth(real colvector x, real colvector om, real scalar right, real scalar p,
    real scalar kcode, real scalar scaleregul, real scalar nmin)
{
    real colvector a, b0, b1, xs, u, K, on, us, Ks
    real matrix V0, V1, X
    real scalar sdx, lo, hi, h0, bn, colc, Cv, gi, gam, vg, prj, Cb2, h, nl, nr, k
    if (right) a = abs(select(x, x :>= 0))
    else a = abs(select(x, x :< 0))
    if (rows(a) < nmin) _error(2001, "too few observations on one side for bandwidth selection")
    sdx = sqrt(mean((x :- mean(x)) :^ 2))
    lo = rdn_kth(a, nmin)
    hi = max(a)
    h0 = min((max((2 * sdx * rows(x)^(-1 / (2 * p + 1)), lo)), hi))
    bn = min((max((rdn_median(a), lo)), hi))
    colc = (right ? 4 : 3)
    b0 = .
    V0 = .
    nl = .
    nr = .
    if (!rdn_cjm(x, om, h0, h0, p, kcode, b0, V0, nl, nr)) {
        _error(2001, "too few observations within the bandwidth on one side")
    }
    Cv = V0[colc, colc] * h0
    b1 = .
    V1 = .
    if (!rdn_cjm(x, om, bn, bn, p + 1, kcode, b1, V1, nl, nr)) {
        _error(2001, "too few observations within the bandwidth on one side")
    }
    gi = 2 * (p + 1) + (colc - 2)
    gam = b1[gi]
    vg = V1[gi, gi]
    xs = select(x, (x :>= -h0) :& (x :<= h0))
    u = xs :/ h0
    K = rdn_kernel(abs(u), kcode) :/ h0
    on = (right ? (xs :>= 0) : (xs :< 0))
    us = select(u, on)
    Ks = select(K, on)
    X = J(rows(us), p + 1, 1)
    for (k = 1; k <= p; k++) X[., k + 1] = us :^ k
    b0 = luinv(cross(X, Ks, X)) * cross(X, Ks, us :^ (p + 1))
    prj = b0[2]
    Cb2 = (prj * gam)^2 + scaleregul * prj^2 * vg
    h = (Cv / (2 * p * max((Cb2, 1e-300))))^(1 / (2 * p + 1))
    return(min((max((h, lo)), hi)))
}

// Cumulative street length within network distance a of the boundary, per
// side: returns the grids and cumulative lengths in G (pointer row vector).
void rdn_street_length(struct rdn_coords scalar C, real scalar center, real scalar window,
    real colvector gl, real colvector cl, real colvector gr, real colvector cr)
{
    real colvector u, v, on, nd, side, cost, len, keep, su, sv, arc, exact, ks, a, b, ex, cc, mm
    real colvector ed, tt, xs, ws, o, grid, cum, fk
    real scalar med, flag, i, tot, pos, k
    u = C.N.fu
    v = C.N.fv
    on = C.nodeon
    nd = rdn_where(on, 0, C.nodedist)
    side = rdn_where(on[u], C.nodeT[v], C.nodeT[u])
    cost = C.N.fcost
    len = C.N.flen
    keep = (!(on[u] :& on[v])) :& (nd[u] :< 1e299) :& (nd[v] :< 1e299)
    su = rdn_where(C.nodesrc[u] :== 0, 1, C.nodesrc[u])
    sv = rdn_where(C.nodesrc[v] :== 0, 1, C.nodesrc[v])
    arc = (rdn_bproject(C.B, C.N.xy[su, .]) + rdn_bproject(C.B, C.N.xy[sv, .])) :/ 2
    keep = keep :& (!rdn_bmasked(C.B, arc))
    if (center < .) keep = keep :& (abs(arc :- center) :<= window)
    exact = (cost :< 1e299)
    fk = selectindex(keep)
    if (rows(fk) == 0) _error(2001, "no streets on either side of the boundary in this window")
    med = rdn_median(cost[fk])
    for (flag = 0; flag <= 1; flag++) {
        ks = selectindex(keep :& (side :== flag))
        if (rows(ks) == 0) _error(2001, "no streets on one side of the boundary in this window")
        a = nd[u[ks]]
        b = nd[v[ks]]
        ex = exact[ks]
        cc = rdn_where(ex, cost[ks], abs(b - a))
        // sample each street at mm points (more on long ones), spreading its length
        mm = rowmax((J(rows(ks), 1, 4), ceil(16 :* cc :/ max((med, 1e-12)))))
        tot = sum(mm)
        xs = J(tot, 1, .)
        ws = J(tot, 1, .)
        pos = 0
        for (i = 1; i <= rows(ks); i++) {
            for (k = 1; k <= mm[i]; k++) {
                pos++
                tt = (k - 0.5) / mm[i]
                if (ex[i]) xs[pos] = min((a[i] + tt * cc[i], b[i] + (1 - tt) * cc[i]))
                else xs[pos] = a[i] + tt * (b[i] - a[i])
                ws[pos] = len[ks[i]] / mm[i]
            }
        }
        o = order((xs, (1::tot)), (1, 2))
        grid = 0 \ xs[o]
        cum = 0 \ runningsum(ws[o])
        if (flag == 0) {
            gl = grid
            cl = cum
        }
        else {
            gr = grid
            cr = cum
        }
    }
}

// The density test. Returns a row vector:
// (diff, diff_bc, se, z, pvalue, theta, f_left, f_right, h_left, h_right,
//  n_left, n_right, n_left_h, n_right_h, p, q)
real rowvector rdn_density(struct rdn_coords scalar C, real colvector sample, real colvector om0,
    real scalar evalindex, real scalar window, real rowvector h, real scalar p, real scalar q,
    string scalar normalize, real scalar kcode, string scalar bwselect, real scalar scaleregul)
{
    real colvector T, s, om, arc, sel, gl, cl, gr, cr, x, bp, bq, reach
    real matrix Vp, Vq
    real scalar center, hl, hr, nmin, nl, nr, nl2, nr2, fl, fr, diff, diffbc, se, zs, pv, th
    T = C.treated
    om = (rows(om0) ? om0 : J(C.nobs, 1, 1))
    if (normalize == "euclidean") {
        // straight-line distance to the boundary line, where it is reached in the plane
        s = rdn_bdistance(C.B, C.xy)
        arc = rdn_bproject(C.B, C.xy)
        reach = !rdn_bmasked(C.B, arc)
        if (!rdn_bclosed(C.B)) reach = reach :& (arc :> 1e-9 * C.B.blen) :& (arc :< C.B.blen * (1 - 1e-9))
    }
    else {
        s = abs(C.s)
        arc = C.crossarc
        reach = !rdn_bmasked(C.B, editmissing(arc, -1))
    }
    sel = (s :< .) :& (editmissing(om, 0) :> 0) :& reach :& sample
    center = .
    if (evalindex < .) {
        if (window >= .) _error(198, "give window() together with evalindex()")
        center = C.enodearc[evalindex]
        sel = sel :& (arc :< .) :& (abs(arc :- center) :<= window)
    }
    if (normalize == "network") {
        gl = .
        cl = .
        gr = .
        cr = .
        rdn_street_length(C, center, window, gl, cl, gr, cr)
        s = editmissing(s, 0)
        s = rdn_where(T, rdn_interp(s, gr, cr), rdn_interp(s, gl, cl))
    }
    // controls exactly on the boundary (s = 0) stay on the control side
    x = select(rdn_where(T, s, rdn_where(s :== 0, -1e-300, -s)), sel)
    om = select(om, sel)
    if (rows(x) < 20 | sum(x :>= 0) < 10 | sum(x :< 0) < 10) {
        _error(2001, "too few observations for a density test")
    }
    if (cols(h) == 0) {
        nmin = max((10, 2 * (q + 2)))
        hl = rdn_side_bandwidth(x, om, 0, p, kcode, scaleregul, nmin)
        hr = rdn_side_bandwidth(x, om, 1, p, kcode, scaleregul, nmin)
        if (bwselect == "common") {
            hl = min((hl, hr))
            hr = hl
        }
    }
    else {
        hl = h[1]
        hr = (cols(h) > 1 ? h[2] : h[1])
    }
    bp = .
    Vp = .
    bq = .
    Vq = .
    nl = .
    nr = .
    nl2 = .
    nr2 = .
    if (!rdn_cjm(x, om, hl, hr, p, kcode, bp, Vp, nl, nr)) {
        _error(2001, "too few observations within the bandwidth on one side")
    }
    if (!rdn_cjm(x, om, hl, hr, q, kcode, bq, Vq, nl2, nr2)) {
        _error(2001, "too few observations within the bandwidth on one side")
    }
    fl = bp[3]
    fr = bp[4]
    diff = fr - fl
    diffbc = bq[4] - bq[3]
    se = sqrt(max((Vq[3, 3] + Vq[4, 4] - 2 * Vq[3, 4], 0)))
    zs = (se > 0 ? diffbc / se : .)
    pv = (zs < . ? 2 * normal(-abs(zs)) : .)
    th = (fl > 0 & fr > 0 ? ln(fr / fl) : .)
    return((diff, diffbc, se, zs, pv, th, fl, fr, hl, hr, sum(x :< 0), sum(x :>= 0), nl, nr, p, q))
}

end

mata:

// ---------------------------------------------------------------------------
// Entry points called from the ado-files (options arrive as locals)
// ---------------------------------------------------------------------------

real scalar rdn_lnum(string scalar name, real scalar dflt)
{
    string scalar s
    s = st_local(name)
    if (s == "") return(dflt)
    return(strtoreal(s))
}

real matrix rdn_lmat(string scalar name)
{
    string scalar s
    s = st_local(name)
    if (s == "") return(J(0, 0, .))
    return(st_matrix(s))
}

// rdnetset: build the network, boundary and coordinates.
void rdn_setup()
{
    external struct rdn_coords scalar rdn_C
    struct rdn_net scalar N
    struct rdn_bnd scalar B
    struct rdn_coords scalar C
    string scalar touse, pmat, costv, side
    real matrix xy, region, marcs, mpoly, line, emat, ET
    real colvector Tin, evalarc
    real scalar ss
    touse = st_local("touse")
    xy = st_data(., (st_local("xvar"), st_local("yvar")), touse)
    if (st_local("tvar") != "") Tin = (st_data(., st_local("tvar"), touse) :!= 0)
    else Tin = J(0, 1, .)
    line = rdn_lmat("bmat")
    emat = rdn_lmat("emat")
    region = rdn_lmat("rmat")
    if (rows(region) == 0) region = J(0, 2, .)
    if (rows(line) == 0) {
        // as in rd2d: the boundary is the line through the evaluation points
        if (cols(emat) != 2) _error(198, "without boundary(), evalpoints() must be a J x 2 matrix of points")
        line = emat[1, .]
        for (ss = 2; ss <= rows(emat); ss++) {
            if (emat[ss, 1] != line[rows(line), 1] | emat[ss, 2] != line[rows(line), 2]) line = line \ emat[ss, .]
        }
        if (rows(line) < 2) _error(198, "evalpoints() needs at least two distinct points")
    }
    side = st_local("side")
    if (side == "" & rows(region) == 0) {
        side = rdn_infer_side(line, xy, Tin)
        printf("{txt}  treated side inferred from treated(): %s of the boundary line\n", side)
    }
    marcs = rdn_lmat("marcs")
    if (rows(marcs) == 0) marcs = J(0, 2, .)
    mpoly = rdn_lmat("mpoly")
    if (rows(mpoly) == 0) mpoly = J(0, 3, .)
    B = rdn_boundary(line, region, side, marcs, mpoly, rdn_lnum("tol", 1e-6))

    costv = st_local("costv")
    if (st_local("fnodes") != "") {
        rdn_read_network(N, st_local("fnodes"), st_local("fedges"), st_local("nodeid"), st_local("nodex"),
            st_local("nodey"), st_local("fromv"), st_local("tov"), costv, st_local("lenv"), st_local("classv"))
    }
    else {
        rdn_read_lines(N, st_local("fedges"), tokens(st_local("ecoords")), costv, st_local("lenv"),
            st_local("classv"), rdn_lnum("jointol", 1e-6))
    }
    N.directed = rdn_lnum("directed", 0)
    N.direction = st_local("direction")
    if (N.direction == "") N.direction = "both"
    N.metric = (st_local("classv") != "" ? 2 : 1)
    N.defpen = rdn_lnum("penalty", 0)
    N.connectors = tokens(st_local("connectors"))'
    pmat = st_local("pmat")
    if (pmat != "") {
        N.penmat = st_matrix(pmat)
        stripe = st_matrixrowstripe(pmat)
        N.penrows = stripe[., 2]
        stripe = st_matrixcolstripe(pmat)
        N.pencols = stripe[., 2]
    }
    else N.penmat = J(0, 0, .)
    ss = rdn_lnum("snapscale", -1)
    // default: 1 for length; otherwise missing, meaning the snapped street's
    // cost per unit length (travel time: that street's seconds per metre)
    if (ss < 0) ss = (costv == "length" ? 1 : .)
    N.snapscale = ss
    if (!rdn_lnum("nosplit", 0)) rdn_split(N, B)
    rdn_edge_frame(N)
    rdn_compile(N)
    C.N = N
    C.B = B

    if (rows(emat) == 0) evalarc = J(0, 1, .)
    else if (cols(emat) == 2) evalarc = rdn_bproject(B, emat)
    else evalarc = vec(emat)
    rdn_build_coords(C, xy, Tin, evalarc, rdn_lnum("neval", 25), rdn_lnum("snapnode", 0))
    rdn_C = C

    st_store(., "_rdn_s", touse, C.s)
    st_store(., "_rdn_treated", touse, C.treated)
    st_store(., "_rdn_arc", touse, C.crossarc)
    st_numscalar("r(N)", C.nobs)
    st_numscalar("r(n_treated)", sum(C.treated))
    st_numscalar("r(n_locations)", C.nloc)
    st_numscalar("r(n_eval)", C.neval)
    st_numscalar("r(n_masked)", sum(C.emasked))
    st_numscalar("r(n_unreachable)", C.nunreach)
    st_numscalar("r(n_crossings)", N.ncross)
    st_numscalar("r(n_nodes)", N.n)
    st_numscalar("r(n_edges)", rows(N.efrom))
    st_numscalar("r(length)", B.blen)
    ET = ((1::C.neval), C.evalarc, C.enodearc, C.exy, C.emasked)
    st_matrix("r(eval)", ET)
    st_matrixcolstripe("r(eval)", (J(6, 1, ""), ("eval" \ "arc" \ "node_arc" \ "x" \ "y" \ "masked")))
}

// Shared: map the estimation sample onto the coordinates' observations.
real colvector rdn_sample_ids(struct rdn_coords scalar C, string scalar touse)
{
    real colvector ids
    if (C.ready != 1) _error(301, "network coordinates not found; run rdnetset first (they are kept in memory only)")
    ids = st_data(., "_rdn_id", touse)
    if (rows(ids) == 0) _error(2000, "no observations")
    if (any((ids :< 1) :| (ids :> C.nobs))) {
        _error(459, "_rdn_id does not match the stored coordinates; rerun rdnetset")
    }
    return(ids)
}

real colvector rdn_spread(real colvector ids, real colvector x, real scalar n, real scalar fillv)
{
    real colvector out
    out = J(n, 1, fillv)
    out[ids] = x
    return(out)
}

// Bandwidths as a J x 2 matrix of (control, treated).
real matrix rdn_bandwidth_matrix(real scalar nJ)
{
    real colvector hv
    real matrix hm
    hm = rdn_lmat("hmat")
    if (rows(hm)) {
        if (rows(hm) != nJ | cols(hm) != 2) _error(503, "hmatrix() must be J x 2 (one row per evaluation point)")
        return(hm)
    }
    if (st_local("hlist") == "") return(J(0, 0, .))
    hv = strtoreal(tokens(st_local("hlist")))'
    if (rows(hv) == 1) return(J(nJ, 2, hv[1]))
    // two numbers are always (h_control h_treated); with two evaluation
    // points, give one bandwidth per point through hmatrix()
    if (rows(hv) == 2) return(J(nJ, 1, 1) * hv')
    if (rows(hv) == nJ) return((hv, hv))
    _error(503, "h() must be one number, (h_control h_treated), or one value per evaluation point (or use hmatrix())")
    return(J(0, 0, .))
}

void rdn_read_common(struct rdn_coords scalar C, string scalar touse, real colvector ids, real colvector y,
    real colvector D, real matrix Z, real colvector om, real colvector cl, real rowvector deriv,
    real matrix tang)
{
    real matrix Zs, Zc
    real colvector tv
    real scalar n
    n = C.nobs
    y = rdn_spread(ids, st_data(., st_local("yvar"), touse), n, .)
    if (st_local("fuzzy") != "") D = rdn_spread(ids, st_data(., st_local("fuzzy"), touse), n, .)
    else D = J(0, 1, .)
    if (st_local("covs") != "") {
        Zs = st_data(., tokens(st_local("covs")), touse)
        if (cols(Zs) > 1) {
            Zc = select(Zs, rowmissing(Zs) :== 0)
            if (rows(Zc) > cols(Zc)) {
                Zc = Zc :- mean(Zc)
                if (rank(cross(Zc, Zc)) < cols(Zc)) {
                    printf("{txt}warning: the covariates are collinear (after removing means); every local fit that\n")
                    printf("{txt}         uses them is singular; drop the redundant variables\n")
                }
            }
        }
        Z = J(n, cols(Zs), .)
        Z[ids, .] = Zs
    }
    else Z = J(0, 0, .)
    if (st_local("wvar") != "") om = rdn_spread(ids, st_data(., st_local("wvar"), touse), n, 0)
    else om = J(n, 1, 1)
    if (st_local("clvar") != "") cl = rdn_spread(ids, st_data(., st_local("clvar"), touse), n, .)
    else cl = J(0, 1, .)
    if (st_local("deriv") != "") deriv = strtoreal(tokens(st_local("deriv")))
    else deriv = (0, 0)
    if (st_local("tangvec") != "") {
        tv = strtoreal(tokens(st_local("tangvec")))'
        tang = J(C.neval, 1, 1) * tv'
    }
    else tang = J(0, 0, .)
}

void rdn_rdnet_entry()
{
    external struct rdn_coords scalar rdn_C
    string scalar touse, vce
    real colvector ids, y, D, om, cl
    real matrix Z, tang, H, R, cov, bw
    real rowvector deriv
    real scalar p, q, crit, minobs, kcode, nJ
    touse = st_local("touse")
    ids = rdn_sample_ids(rdn_C, touse)
    y = .
    D = .
    Z = .
    om = .
    cl = .
    deriv = .
    tang = .
    rdn_read_common(rdn_C, touse, ids, y, D, Z, om, cl, deriv, tang)
    nJ = rdn_C.neval
    p = rdn_lnum("p", 1)
    q = rdn_lnum("q", p + 1)
    if (p < 0 | q < p) _error(198, "need 0 <= p <= q")
    if (cols(deriv) != 2 | min(deriv) < 0 | sum(deriv) > p) _error(198, "deriv() must be two integers a c >= 0 with a + c <= p")
    if (rows(tang)) {
        if (sum(deriv) != 0) _error(198, "give either deriv() or tangvec(), not both")
        if (p < 1) _error(198, "tangvec() needs p >= 1")
    }
    if (rdn_lnum("pooled", 0) & (sum(deriv) > 0 | rows(tang) > 0)) _error(198, "pooled estimates level jumps only")
    kcode = rdn_kernel_code(st_local("kernel"))
    vce = st_local("vcetype")
    minobs = rdn_lnum("minobs", -1)
    if (minobs < 0) minobs = 20 + rdn_nterms(q)
    H = rdn_bandwidth_matrix(nJ)
    cov = .
    crit = .
    bw = .
    R = rdn_rdnet(rdn_C, y, D, Z, om, cl, H, p, q, deriv, tang, kcode, vce,
        rdn_lnum("cutoff", .), rdn_lnum("euclid", 0), rdn_lnum("rscale", 3), rdn_lnum("rorder", .),
        rdn_lnum("pooled", 0), rdn_lnum("level", 95), rdn_lnum("cbands", 1), rdn_lnum("repp", 2000),
        minobs, rdn_lnum("seed", .), rdn_lnum("scaleregul", 3), rdn_lnum("scalebiascrct", 1),
        st_local("bwselect"), cov, crit, bw)
    st_matrix(st_local("resmat"), R)
    st_matrixcolstripe(st_local("resmat"), (J(22, 1, ""), rdn_result_names()'))
    st_matrixrowstripe(st_local("resmat"), (J(nJ, 1, ""), "eval" :+ strofreal(1::nJ)))
    if (rows(cov)) {
        st_matrix(st_local("covmat"), cov)
        st_matrixcolstripe(st_local("covmat"), (J(nJ, 1, ""), "eval" :+ strofreal(1::nJ)))
        st_matrixrowstripe(st_local("covmat"), (J(nJ, 1, ""), "eval" :+ strofreal(1::nJ)))
    }
    if (rows(bw)) {
        st_matrix(st_local("bwmat"), bw)
        st_matrixcolstripe(st_local("bwmat"), (J(24, 1, ""), rdn_bw_names()'))
    }
    st_local("crit", strofreal(crit, "%21.0g"))
    st_local("nobs", strofreal(rows(ids)))
    st_local("kname", rdn_kernel_name(kcode))
    st_local("q", strofreal(q))
    st_local("minobs", strofreal(minobs))
}

void rdn_rdbwnet_entry()
{
    external struct rdn_coords scalar rdn_C
    string scalar touse, bwsel
    real colvector ids, y, D, om, cl
    real matrix Z, tang, bw
    real rowvector deriv
    real scalar p, bwcheck
    touse = st_local("touse")
    ids = rdn_sample_ids(rdn_C, touse)
    y = .
    D = .
    Z = .
    om = .
    cl = .
    deriv = .
    tang = .
    rdn_read_common(rdn_C, touse, ids, y, D, Z, om, cl, deriv, tang)
    p = rdn_lnum("p", 1)
    if (cols(deriv) != 2 | min(deriv) < 0 | sum(deriv) > p) _error(198, "deriv() must be two integers a c >= 0 with a + c <= p")
    bwcheck = rdn_lnum("bwcheck", -1)
    if (bwcheck < 0) bwcheck = 20 + rdn_nterms(p + 1)
    bwsel = strlower(st_local("bwselect"))
    if (bwsel == "") bwsel = "mserd"
    if (!anyof(("mserd", "msetwo", "imserd", "imsetwo", "cerrd", "certwo", "icerrd", "icertwo"), bwsel)) {
        _error(198, "bwselect() must be mserd, msetwo, imserd, imsetwo, cerrd, certwo, icerrd or icertwo")
    }
    bw = rdn_rdbwnet(rdn_C, y, p, deriv, tang, Z, rdn_kernel_code(st_local("kernel")), bwsel,
        st_local("vcetype"), cl, D, om, bwcheck, rdn_lnum("scaleregul", 1), rdn_lnum("scalebiascrct", 1))
    st_matrix(st_local("bwmat"), bw)
    st_matrixcolstripe(st_local("bwmat"), (J(24, 1, ""), rdn_bw_names()'))
    st_matrixrowstripe(st_local("bwmat"), (J(rows(bw), 1, ""), "eval" :+ strofreal(1::rows(bw))))
}

void rdn_density_entry()
{
    external struct rdn_coords scalar rdn_C
    string scalar touse, norm, bws
    real colvector ids, sample, om, seen, first
    real rowvector h, out
    real scalar p, q, i
    touse = st_local("touse")
    ids = rdn_sample_ids(rdn_C, touse)
    sample = rdn_spread(ids, J(rows(ids), 1, 1), rdn_C.nobs, 0)
    // repeated measurements: keep the first observation at each location
    seen = J(rdn_C.nloc, 1, 0)
    first = J(rdn_C.nobs, 1, 0)
    for (i = 1; i <= rdn_C.nobs; i++) {
        if (sample[i] & !seen[rdn_C.loc[i]]) {
            first[i] = 1
            seen[rdn_C.loc[i]] = 1
        }
    }
    if (rdn_lnum("unique", 0)) sample = sample :* first
    else if (sum(sample :& !first) > 0) {
        printf("{txt}note: several observations share a location (repeated measurements?); each is counted.\n")
        printf("{txt}      Use option unique to test the density of locations.\n")
    }
    if (st_local("wvar") != "") om = rdn_spread(ids, st_data(., st_local("wvar"), touse), rdn_C.nobs, 0)
    else om = J(0, 1, .)
    p = rdn_lnum("p", 2)
    q = rdn_lnum("q", p + 1)
    if (p < 1 | q < p) _error(198, "need 1 <= p <= q")
    norm = st_local("normalize")
    if (norm == "") norm = "euclidean"
    if (norm != "euclidean" & norm != "network" & norm != "none") {
        _error(198, "normalize() must be euclidean, network or none")
    }
    rdn_warn_on_network(rdn_C)
    bws = st_local("bwselect")
    if (bws == "") bws = "common"
    if (bws != "common" & bws != "each") _error(198, "bwselect() must be common or each")
    if (st_local("hlist") != "") h = strtoreal(tokens(st_local("hlist")))
    else h = J(1, 0, .)
    out = rdn_density(rdn_C, sample, om, rdn_lnum("evalindex", .), rdn_lnum("window", .), h, p, q,
        norm, rdn_kernel_code(st_local("kernel")), bws, rdn_lnum("scaleregul", 1))
    st_numscalar("r(diff)", out[1])
    st_numscalar("r(diff_bc)", out[2])
    st_numscalar("r(se)", out[3])
    st_numscalar("r(z)", out[4])
    st_numscalar("r(pvalue)", out[5])
    st_numscalar("r(theta)", out[6])
    st_numscalar("r(f_left)", out[7])
    st_numscalar("r(f_right)", out[8])
    st_numscalar("r(h_left)", out[9])
    st_numscalar("r(h_right)", out[10])
    st_numscalar("r(N_left)", out[11])
    st_numscalar("r(N_right)", out[12])
    st_numscalar("r(N_h_left)", out[13])
    st_numscalar("r(N_h_right)", out[14])
    st_numscalar("r(p)", out[15])
    st_numscalar("r(q)", out[16])
    st_global("r(normalize)", norm)
}

// ---------------------------------------------------------------------------
// Manipulation (sorting) test of Bugni, Crippa and Restrepo (2026), as in the
// Python and R packages: the same clusters (exact 1-D k-means of positions
// along the boundary) and the same relabelling draws (Park-Miller stream).
// ---------------------------------------------------------------------------

// reps x q Bernoulli(1/2) labels: x <- 16807 x mod (2^31 - 1), treated when
// x <= 2^30 - 1.
real matrix rdn_lcg_labels(real scalar seed, real scalar reps, real scalar q)
{
    real matrix out
    real scalar x, r, j
    x = mod(seed, 2147483646) + 1
    out = J(reps, q, 0)
    for (r = 1; r <= reps; r++) {
        for (j = 1; j <= q; j++) {
            x = mod(x * 16807, 2147483647)
            out[r, j] = (x <= 1073741823)
        }
    }
    return(out)
}

// Optimal 1-D k-means partitions of sorted positions for k = 1..kmax: row k
// holds the indices (into pos) of the first element of each of its k clusters.
real matrix rdn_kmeans_parts(real colvector pos, real scalar kmax0)
{
    real colvector vals, start, w, s1, s2, sw, ii, cc
    real matrix D, A, F
    real scalar m, kmax, k, j, i, best, bi, kk, n, a, t
    vals = J(0, 1, .)
    start = J(0, 1, .)
    w = J(0, 1, .)
    for (i = 1; i <= rows(pos); i++) {
        if (i == 1 | pos[i] != pos[i - 1]) {
            vals = vals \ pos[i]
            start = start \ i
            w = w \ 1
        }
        else w[rows(w)] = w[rows(w)] + 1
    }
    m = rows(vals)
    kmax = min((kmax0, m))
    s1 = 0 \ runningsum(w :* vals)
    s2 = 0 \ runningsum(w :* vals :* vals)
    sw = 0 \ runningsum(w)
    D = J(kmax, m, .)
    A = J(kmax, m, 1)
    for (j = 1; j <= m; j++) {
        n = sw[j + 1] - sw[1]
        a = s1[j + 1] - s1[1]
        D[1, j] = max(((s2[j + 1] - s2[1]) - a * a / n, 0))
    }
    for (k = 2; k <= kmax; k++) {
        for (j = k; j <= m; j++) {
            ii = (k::j)
            cc = J(rows(ii), 1, .)
            for (t = 1; t <= rows(ii); t++) {
                i = ii[t]
                n = sw[j + 1] - sw[i]
                a = s1[j + 1] - s1[i]
                cc[t] = D[k - 1, i - 1] + max(((s2[j + 1] - s2[i]) - a * a / n, 0))
            }
            best = min(cc)
            bi = ii[selectindex(cc :<= best + 1e-12 * max((1, abs(best))))[1]]
            D[k, j] = cc[bi - k + 1]
            A[k, j] = bi
        }
    }
    F = J(kmax, kmax, .)
    for (k = 1; k <= kmax; k++) {
        j = m
        kk = k
        while (kk >= 1) {
            i = (kk > 1 ? A[kk, j] : 1)
            F[k, kk] = start[i]
            j = i - 1
            kk = kk - 1
        }
    }
    return(F)
}

// The statistic S for each row of labels (draws x q, in position order).
real colvector rdn_mt_stat(real matrix lab, real matrix F, real scalar q)
{
    real matrix P, t, dev
    real rowvector f, ends, n
    real colvector S
    real scalar k, c
    P = J(rows(lab), q + 1, 0)
    for (c = 1; c <= q; c++) P[., c + 1] = P[., c] + lab[., c]
    S = J(rows(lab), 1, 0)
    for (k = 1; k <= rows(F); k++) {
        f = F[k, (1..k)]
        ends = (k > 1 ? (f[(2..k)], q + 1) : (q + 1))
        n = ends - f
        t = P[., ends] - P[., f]
        dev = abs(t :/ n :- 0.5) :* sqrt(n)
        S = rowmax((S, rowmax(dev)))
    }
    return(S)
}

void rdn_manip_entry()
{
    external struct rdn_coords scalar rdn_C
    real colvector ids, sample, seen, firsts, T, dist, arc, ok, idx, o, sel, pos, srt, gaps, lab, null, S
    real rowvector qs
    real matrix F, path
    real scalar reps, seed, level, kmax, nuse, qi, qq, L, g, cut, pv, crit, km, i
    ids = rdn_sample_ids(rdn_C, st_local("touse"))
    sample = rdn_spread(ids, J(rows(ids), 1, 1), rdn_C.nobs, 0)
    rdn_warn_on_network(rdn_C)
    // one observation per location within the sample
    seen = J(rdn_C.nloc, 1, 0)
    firsts = J(0, 1, .)
    for (i = 1; i <= rdn_C.nobs; i++) {
        if (sample[i] & !seen[rdn_C.loc[i]]) {
            firsts = firsts \ i
            seen[rdn_C.loc[i]] = 1
        }
    }
    reps = rdn_lnum("reps", 10000)
    seed = rdn_lnum("seed", 0)
    level = rdn_lnum("level", 95)
    kmax = rdn_lnum("kmax", .)
    if (reps < 100) _error(198, "reps() must be at least 100")
    T = rdn_C.treated[firsts]
    dist = rdn_bdistance(rdn_C.B, rdn_C.xy[firsts, .])
    arc = rdn_bproject(rdn_C.B, rdn_C.xy[firsts, .])
    L = rdn_C.B.blen
    ok = (dist :< .) :& (dist :> rdn_C.B.tol) :& (!rdn_bmasked(rdn_C.B, arc))
    if (!rdn_bclosed(rdn_C.B)) ok = ok :& (arc :> 1e-9 * L) :& (arc :< L * (1 - 1e-9))
    idx = selectindex(ok)
    nuse = rows(idx)
    if (st_local("qlist") != "") qs = strtoreal(tokens(st_local("qlist")))
    else qs = (nuse >= 2000 ? 200 : max((floor(nuse / 10), 1)))
    path = J(cols(qs), 7, .)
    for (qi = 1; qi <= cols(qs); qi++) {
        qq = qs[qi]
        if (qq < 2 | qq > nuse) {
            _error(198, sprintf("q() must be between 2 and the number of usable locations (%g)", nuse))
        }
        if (qq < 20) printf("{txt}note: q = %g locations give a test with little power\n", qq)
        if (qq > 0.4 * nuse) {
            printf("{txt}note: q = %g is a large share of the %g usable locations; the test assumes the\n", qq, nuse)
            printf("{txt}      selected locations are close to the boundary\n")
        }
        o = order((dist[idx], idx), (1, 2))  // nearest first; ties by location order
        sel = idx[o[(1::qq)]]
        pos = arc[sel]
        if (rdn_bclosed(rdn_C.B)) {  // cut the ring at the largest gap between neighbouring positions
            srt = uniqrows(pos)
            if (rows(srt) > 1) {
                gaps = (srt[(2::rows(srt))] \ srt[1] + L) - srt
                g = rdn_argmin(-gaps)
                cut = srt[mod(g, rows(srt)) + 1]
            }
            else cut = srt[1]
            pos = mod(pos :- cut, L)
        }
        o = order((pos, sel), (1, 2))  // position order; ties by location order
        pos = pos[o]
        lab = T[sel[o]]
        km = (kmax < . ? min((kmax, qq)) : qq)
        F = rdn_kmeans_parts(pos, km)
        S = rdn_mt_stat(lab', F, qq)
        null = rdn_mt_stat(rdn_lcg_labels(seed, reps, qq), F, qq)
        pv = mean(null :>= S - 1e-12 * max((1, S)))
        crit = sort(null, 1)[ceil(level / 100 * rows(null))]
        path[qi, .] = (qq, S, pv, crit, sum(lab), max(dist[sel]), rows(F))
    }
    st_numscalar("r(q)", path[1, 1])
    st_numscalar("r(statistic)", path[1, 2])
    st_numscalar("r(pvalue)", path[1, 3])
    st_numscalar("r(crit)", path[1, 4])
    st_numscalar("r(n_treated)", path[1, 5])
    st_numscalar("r(distance)", path[1, 6])
    st_numscalar("r(k_max)", path[1, 7])
    st_numscalar("r(reps)", reps)
    st_numscalar("r(N_locations)", nuse)
    st_matrix("r(path)", path)
    st_matrixcolstripe("r(path)", (J(7, 1, ""), ("q" \ "statistic" \ "pvalue" \ "crit" \ "n_treated" \ "distance" \ "k_max")))
}

// Note when most locations lie on the network itself (sorting tests assume
// units spread over the plane).
void rdn_warn_on_network(struct rdn_coords scalar C)
{
    real colvector dep
    real scalar ext
    dep = C.depth[C.firstl]
    ext = sqrt((max(C.N.xy[., 1]) - min(C.N.xy[., 1]))^2 + (max(C.N.xy[., 2]) - min(C.N.xy[., 2]))^2)
    if (rows(dep)) {
        if (mean(dep :<= 1e-6 * max((ext, 1))) > 0.5) {
            printf("{txt}warning: most locations lie on the network itself (stations on a river, points on streets);\n")
            printf("{txt}         sorting tests assume units spread over the plane, and in Monte Carlo runs none of\n")
            printf("{txt}         them held its size for points on streets near a curved boundary\n")
        }
    }
}

// Aggregate effects from e() of rdnet: WBATE (weights from a Stata matrix)
// and LBATE.
void rdn_wbate_entry()
{
    real matrix R, cov
    real colvector w, ok, okw, bad
    real scalar nJ, z, tq, se, lev
    R = st_matrix("e(results)")
    cov = st_matrix("e(cov)")
    lev = st_numscalar("e(level)")
    nJ = rows(R)
    w = vec(st_matrix(st_local("wmat")))
    if (rows(w) != nJ) _error(503, "WBATE weights must have one entry per evaluation point (" + strofreal(nJ) + ")")
    if (hasmissing(w)) _error(504, "WBATE weights must be finite")
    if (sum(w) == 0) _error(498, "WBATE weights must have a nonzero sum")
    ok = (R[., 15] :== 0)
    bad = selectindex((!ok) :& (w :!= 0))
    if (rows(bad)) _error(498, "WBATE puts weight on evaluation points that were not estimated; give them weight 0")
    okw = selectindex(ok)
    w = w :/ sum(w)
    tq = w[okw]' * R[okw, 10]
    se = sqrt(max((w[okw]' * cov[okw, okw] * w[okw], 0)))
    z = invnormal(0.5 + lev / 200)
    st_numscalar("r(tau)", w[okw]' * R[okw, 8])
    st_numscalar("r(tau_bc)", tq)
    st_numscalar("r(se_rb)", se)
    st_numscalar("r(ci_lo)", tq - z * se)
    st_numscalar("r(ci_hi)", tq + z * se)
    st_numscalar("r(pvalue)", (se > 0 ? 2 * normal(-abs(tq / se)) : .))
}

void rdn_lbate_entry()
{
    real matrix R
    real colvector okw, lo, hi, tp, tq
    real scalar crit
    string scalar kind
    R = st_matrix("e(results)")
    crit = st_numscalar("e(crit)")
    kind = st_local("kind")
    okw = selectindex(R[., 15] :== 0)
    if (rows(okw) == 0 | crit >= .) _error(498, "LBATE needs at least one estimated evaluation point")
    tp = R[okw, 8]
    tq = R[okw, 10]
    lo = tq - crit :* R[okw, 11]
    hi = tq + crit :* R[okw, 11]
    if (kind == "min") {
        st_numscalar("r(tau)", min(tp))
        st_numscalar("r(tau_bc)", min(tq))
        st_numscalar("r(ci_lo)", min(lo))
        st_numscalar("r(ci_hi)", min(hi))
    }
    else if (kind == "absmax") {
        // the interval for max_j |tau_j| implied by the uniform band
        st_numscalar("r(tau)", max(abs(tp)))
        st_numscalar("r(tau_bc)", max(abs(tq)))
        st_numscalar("r(ci_lo)", max(rowmax((J(rows(lo), 1, 0), lo, -hi))))
        st_numscalar("r(ci_hi)", max(rowmax((abs(lo), abs(hi)))))
    }
    else {
        st_numscalar("r(tau)", max(tp))
        st_numscalar("r(tau_bc)", max(tq))
        st_numscalar("r(ci_lo)", max(lo))
        st_numscalar("r(ci_hi)", max(hi))
    }
}

// Equal weight per distinct place (node_arc, column 3) among the estimated
// points (status 0, column 15); points at the same place share its weight.
real colvector rdn_place_weights(real matrix R)
{
    real colvector w, ok, place
    real scalar j, k, c
    w = J(rows(R), 1, 0)
    ok = (R[., 15] :== 0)
    place = round(R[., 3], 1e-9)
    for (j = 1; j <= rows(R); j++) {
        if (!ok[j]) continue
        c = 0
        for (k = 1; k <= rows(R); k++) c = c + (ok[k] & place[k] == place[j])
        w[j] = 1 / c
    }
    return(w)
}

end
