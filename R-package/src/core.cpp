// Compiled core of rd2dnet: shortest paths on (state-expanded) networks,
// nearest-segment search, point-in-polygon tests and the sparse Conley kernel.
//
// Conventions: graphs arrive in compressed sparse row form with 0-based
// `indptr`/`indices`; node and segment indices returned to R are 1-based.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

const double INF = std::numeric_limits<double>::infinity();

typedef std::pair<double, int> QItem;
typedef std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem> > MinHeap;

// Shortest-path workspace reused across searches: only touched entries are reset.
struct Workspace {
  std::vector<double> dist;
  std::vector<int> origin;
  std::vector<char> done;
  std::vector<int> touched;
  explicit Workspace(int n) : dist(n, INF), origin(n, -1), done(n, 0) {}
  void reset() {
    for (size_t k = 0; k < touched.size(); ++k) {
      int s = touched[k];
      dist[s] = INF;
      origin[s] = -1;
      done[s] = 0;
    }
    touched.clear();
  }
};

// Dijkstra from the given start states (all at distance 0), stopping at `limit`.
void dijkstra(const int* indptr, const int* indices, const double* w,
              const std::vector<int>& starts, const std::vector<int>& labels,
              double limit, Workspace& ws) {
  MinHeap pq;
  for (size_t k = 0; k < starts.size(); ++k) {
    int s = starts[k];
    if (ws.dist[s] == INF) ws.touched.push_back(s);
    if (0.0 < ws.dist[s] || ws.origin[s] < 0) {
      ws.dist[s] = 0.0;
      ws.origin[s] = labels[k];
      pq.push(QItem(0.0, s));
    }
  }
  while (!pq.empty()) {
    QItem top = pq.top();
    pq.pop();
    double d = top.first;
    int u = top.second;
    if (ws.done[u] || d > ws.dist[u]) continue;
    ws.done[u] = 1;
    for (int k = indptr[u]; k < indptr[u + 1]; ++k) {
      int v = indices[k];
      double nd = d + w[k];
      if (nd > limit || ws.done[v]) continue;
      if (nd < ws.dist[v]) {
        if (ws.dist[v] == INF) ws.touched.push_back(v);
        ws.dist[v] = nd;
        ws.origin[v] = ws.origin[u];
        pq.push(QItem(nd, v));
      } else if (nd == ws.dist[v] && ws.origin[u] < ws.origin[v]) {
        ws.origin[v] = ws.origin[u]; // ties go to the first source
      }
    }
  }
}

}  // namespace

// Distances from each start state to every node (minimum over a node's states).
// [[Rcpp::export(.dijkstra_distances)]]
NumericMatrix dijkstra_distances(IntegerVector indptr, IntegerVector indices, NumericVector w,
                                 int n_states, IntegerVector starts, IntegerVector state_node,
                                 int n_nodes, double limit) {
  int k = starts.size();
  NumericMatrix out(k, n_nodes);
  std::fill(out.begin(), out.end(), INF);
  Workspace ws(n_states);
  std::vector<int> st(1), lab(1, 0);
  for (int i = 0; i < k; ++i) {
    st[0] = starts[i];
    dijkstra(indptr.begin(), indices.begin(), w.begin(), st, lab, limit, ws);
    for (size_t t = 0; t < ws.touched.size(); ++t) {
      int s = ws.touched[t];
      int node = state_node[s];
      if (ws.dist[s] < out(i, node)) out(i, node) = ws.dist[s];
    }
    ws.reset();
  }
  return out;
}

// Distance from the nearest start state to every node, and which source (1-based
// position in `starts`) it was; 0 marks unreachable nodes.
// [[Rcpp::export(.dijkstra_nearest)]]
List dijkstra_nearest(IntegerVector indptr, IntegerVector indices, NumericVector w,
                      int n_states, IntegerVector starts, IntegerVector state_node, int n_nodes) {
  Workspace ws(n_states);
  std::vector<int> st(starts.begin(), starts.end());
  std::vector<int> lab(st.size());
  for (size_t i = 0; i < lab.size(); ++i) lab[i] = static_cast<int>(i) + 1;
  dijkstra(indptr.begin(), indices.begin(), w.begin(), st, lab, INF, ws);
  NumericVector dist(n_nodes, INF);
  IntegerVector src(n_nodes, 0);
  for (int s = 0; s < n_states; ++s) {
    int node = state_node[s];
    if (ws.dist[s] < dist[node]) {
      dist[node] = ws.dist[s];
      src[node] = ws.origin[s];
    }
  }
  return List::create(Named("dist") = dist, Named("src") = src);
}

// Every (source, node) pair closer than `limit`: source position, node (both
// 1-based) and the shortest distance.
// [[Rcpp::export(.dijkstra_pairs)]]
List dijkstra_pairs(IntegerVector indptr, IntegerVector indices, NumericVector w,
                    int n_states, IntegerVector starts, IntegerVector state_node,
                    int n_nodes, double limit) {
  Workspace ws(n_states);
  std::vector<int> st(1), lab(1, 0);
  std::vector<double> best(n_nodes, INF);
  std::vector<int> seen;
  std::vector<int> oi, on;
  std::vector<double> od;
  for (int i = 0; i < starts.size(); ++i) {
    st[0] = starts[i];
    dijkstra(indptr.begin(), indices.begin(), w.begin(), st, lab, limit, ws);
    for (size_t t = 0; t < ws.touched.size(); ++t) {
      int s = ws.touched[t];
      double d = ws.dist[s];
      if (!(d < limit)) continue;
      int node = state_node[s];
      if (best[node] == INF) seen.push_back(node);
      if (d < best[node]) best[node] = d;
    }
    std::sort(seen.begin(), seen.end());
    for (size_t t = 0; t < seen.size(); ++t) {
      oi.push_back(i + 1);
      on.push_back(seen[t] + 1);
      od.push_back(best[seen[t]]);
      best[seen[t]] = INF;
    }
    seen.clear();
    ws.reset();
  }
  return List::create(Named("i") = wrap(oi), Named("node") = wrap(on), Named("dist") = wrap(od));
}

namespace {

// Uniform grid over segment bounding boxes for nearest-segment queries.
struct SegmentGrid {
  const double *x0, *y0, *x1, *y1;
  int n;
  double xmin, ymin, cs;
  int nx, ny;
  std::vector<std::vector<int> > cells;

  SegmentGrid(const double* ax0, const double* ay0, const double* ax1, const double* ay1, int an)
      : x0(ax0), y0(ay0), x1(ax1), y1(ay1), n(an) {
    xmin = INF;
    ymin = INF;
    double xmax = -INF, ymax = -INF;
    for (int s = 0; s < n; ++s) {
      xmin = std::min(xmin, std::min(x0[s], x1[s]));
      xmax = std::max(xmax, std::max(x0[s], x1[s]));
      ymin = std::min(ymin, std::min(y0[s], y1[s]));
      ymax = std::max(ymax, std::max(y0[s], y1[s]));
    }
    double W = xmax - xmin, H = ymax - ymin;
    double area = W * H;
    if (area > 0) {
      cs = std::sqrt(area / std::max(n, 1));
    } else {
      cs = std::max(W, H) / std::max(n, 1);
    }
    if (!(cs > 0)) cs = 1.0;
    nx = std::min(4096, static_cast<int>(std::ceil(W / cs)) + 1);
    ny = std::min(4096, static_cast<int>(std::ceil(H / cs)) + 1);
    cs = std::max(cs, std::max(W / (nx - 0.5), H / (ny - 0.5)));
    cells.resize(static_cast<size_t>(nx) * ny);
    for (int s = 0; s < n; ++s) {
      int ca = cellx(std::min(x0[s], x1[s])), cb = cellx(std::max(x0[s], x1[s]));
      int ra = celly(std::min(y0[s], y1[s])), rb = celly(std::max(y0[s], y1[s]));
      for (int cx = ca; cx <= cb; ++cx)
        for (int cy = ra; cy <= rb; ++cy) cells[static_cast<size_t>(cy) * nx + cx].push_back(s);
    }
  }
  int cellx(double x) const {
    int c = static_cast<int>(std::floor((x - xmin) / cs));
    return std::max(0, std::min(nx - 1, c));
  }
  int celly(double y) const {
    int c = static_cast<int>(std::floor((y - ymin) / cs));
    return std::max(0, std::min(ny - 1, c));
  }
};

inline double seg_dist(double px, double py, double ax, double ay, double bx, double by, double* t) {
  double dx = bx - ax, dy = by - ay;
  double L2 = dx * dx + dy * dy;
  double tt = 0.0;
  if (L2 > 0) {
    tt = ((px - ax) * dx + (py - ay) * dy) / L2;
    tt = std::max(0.0, std::min(1.0, tt));
  }
  *t = tt;
  double qx = ax + tt * dx - px, qy = ay + tt * dy - py;
  return std::sqrt(qx * qx + qy * qy);
}

}  // namespace

// Nearest segment to each point: segment (1-based; ties go to the lowest index),
// distance, and position t in [0, 1] of the foot along the segment.
// [[Rcpp::export(.nearest_segment)]]
List nearest_segment(NumericVector px, NumericVector py, NumericVector x0, NumericVector y0,
                     NumericVector x1, NumericVector y1) {
  int m = px.size(), n = x0.size();
  IntegerVector seg(m, NA_INTEGER);
  NumericVector dist(m, NA_REAL), tpos(m, NA_REAL);
  if (n == 0) return List::create(Named("seg") = seg, Named("dist") = dist, Named("t") = tpos);
  SegmentGrid g(x0.begin(), y0.begin(), x1.begin(), y1.begin(), n);
  std::vector<int> stamp(n, -1);
  int maxr = std::max(g.nx, g.ny);
  for (int i = 0; i < m; ++i) {
    double p = px[i], q = py[i];
    if (!std::isfinite(p) || !std::isfinite(q)) continue;
    int cx = g.cellx(p), cy = g.celly(q);
    double best = INF, bt = 0.0;
    int bs = -1;
    for (int r = 0; r <= maxr; ++r) {
      double lb = r == 0 ? 0.0 : (r - 1) * g.cs;
      if (lb > best) break;
      for (int ix = cx - r; ix <= cx + r; ++ix) {
        if (ix < 0 || ix >= g.nx) continue;
        for (int iy = cy - r; iy <= cy + r; ++iy) {
          if (iy < 0 || iy >= g.ny) continue;
          if (std::abs(ix - cx) != r && std::abs(iy - cy) != r) continue;
          const std::vector<int>& cell = g.cells[static_cast<size_t>(iy) * g.nx + ix];
          for (size_t k = 0; k < cell.size(); ++k) {
            int s = cell[k];
            if (stamp[s] == i) continue;
            stamp[s] = i;
            double t;
            double d = seg_dist(p, q, x0[s], y0[s], x1[s], y1[s], &t);
            if (d < best || (d == best && s < bs)) {
              best = d;
              bs = s;
              bt = t;
            }
          }
        }
      }
    }
    seg[i] = bs + 1;
    dist[i] = best;
    tpos[i] = bt;
  }
  return List::create(Named("seg") = seg, Named("dist") = dist, Named("t") = tpos);
}

// Even-odd point-in-polygon test over all rings (holes and multipolygons
// included). Rings are given as concatenated vertices with a ring id each.
// [[Rcpp::export(.points_in_rings)]]
LogicalVector points_in_rings(NumericVector px, NumericVector py, NumericVector rx,
                              NumericVector ry, IntegerVector ring) {
  int m = px.size(), n = rx.size();
  LogicalVector out(m, false);
  for (int i = 0; i < m; ++i) {
    double x = px[i], y = py[i];
    bool inside = false;
    int start = 0;
    while (start < n) {
      int end = start;
      while (end + 1 < n && ring[end + 1] == ring[start]) ++end;
      for (int a = start, b = end; a <= end; b = a++) {
        double ya = ry[a], yb = ry[b];
        if ((ya > y) != (yb > y)) {
          double xc = rx[a] + (y - ya) * (rx[b] - rx[a]) / (yb - ya);
          if (x < xc) inside = !inside;
        }
      }
      start = end + 1;
    }
    out[i] = inside;
  }
  return out;
}

// Conley (Bartlett) quadratic form P' K P with K_ij = max(0, 1 - d_ij / cutoff)
// over Euclidean distance; pairs come from a grid with cell size `cutoff`.
// [[Rcpp::export(.hac_euclidean)]]
NumericMatrix hac_euclidean(NumericVector x, NumericVector y, NumericMatrix P, double cutoff) {
  int m = x.size(), k = P.ncol();
  NumericMatrix out(k, k);
  if (m == 0) return out;
  double xmin = *std::min_element(x.begin(), x.end());
  double ymin = *std::min_element(y.begin(), y.end());
  double xmax = *std::max_element(x.begin(), x.end());
  double ymax = *std::max_element(y.begin(), y.end());
  double cs = cutoff;
  double span = std::max(xmax - xmin, ymax - ymin);
  if (!(cs > 0) || span / cs > 4000.0) cs = std::max(span / 4000.0, 1e-300);
  int nx = static_cast<int>(std::floor((xmax - xmin) / cs)) + 1;
  int ny = static_cast<int>(std::floor((ymax - ymin) / cs)) + 1;
  std::vector<std::vector<int> > cells(static_cast<size_t>(nx) * ny);
  std::vector<int> cxs(m), cys(m);
  for (int i = 0; i < m; ++i) {
    cxs[i] = std::min(nx - 1, static_cast<int>(std::floor((x[i] - xmin) / cs)));
    cys[i] = std::min(ny - 1, static_cast<int>(std::floor((y[i] - ymin) / cs)));
    cells[static_cast<size_t>(cys[i]) * nx + cxs[i]].push_back(i);
  }
  int reach = static_cast<int>(std::ceil(cutoff / cs));
  std::vector<double> S(k);
  for (int i = 0; i < m; ++i) {
    std::fill(S.begin(), S.end(), 0.0);
    for (int ix = std::max(0, cxs[i] - reach); ix <= std::min(nx - 1, cxs[i] + reach); ++ix) {
      for (int iy = std::max(0, cys[i] - reach); iy <= std::min(ny - 1, cys[i] + reach); ++iy) {
        const std::vector<int>& cell = cells[static_cast<size_t>(iy) * nx + ix];
        for (size_t t = 0; t < cell.size(); ++t) {
          int j = cell[t];
          double d = std::sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
          double kij = 1.0 - d / cutoff;
          if (!(kij > 0)) continue;
          for (int c = 0; c < k; ++c) S[c] += kij * P(j, c);
        }
      }
    }
    for (int a = 0; a < k; ++a) {
      double pa = P(i, a);
      if (pa == 0.0) continue;
      for (int b = 0; b < k; ++b) out(a, b) += pa * S[b];
    }
  }
  return out;
}

// Conley (Bartlett) quadratic form P' K P over network distance.
//
// Location i sits on the edge (u_i, v_i) at costs (cu_i, cv_i) from its ends,
// after an off-network depth cost dep_i. Node pairs (a, b) closer than the
// cutoff come in as triplets with 1-based compact node ids. The distance from
// i to j is the cheapest of cu/cv_i + D(a, b) + cu/cv_j, or the direct distance
// along a shared edge, plus both depths; co-located observations share one
// location, at distance 0 from itself.
// [[Rcpp::export(.hac_network)]]
NumericMatrix hac_network(IntegerVector lu, IntegerVector lv, NumericVector cu, NumericVector cv,
                          NumericVector dep, IntegerVector na, IntegerVector nb, NumericVector nd,
                          int n_nodes, NumericMatrix P, double cutoff) {
  int m = lu.size(), k = P.ncol();
  NumericMatrix out(k, k);
  // node pairs by first node
  std::vector<int> ptr(n_nodes + 1, 0);
  for (int t = 0; t < na.size(); ++t) ptr[na[t]]++;
  for (int a = 0; a < n_nodes; ++a) ptr[a + 1] += ptr[a];
  std::vector<int> pb(na.size());
  std::vector<double> pd(na.size());
  {
    std::vector<int> fill(ptr.begin(), ptr.end() - 1);
    for (int t = 0; t < na.size(); ++t) {
      int a = na[t] - 1;
      pb[fill[a]] = nb[t] - 1;
      pd[fill[a]] = nd[t];
      fill[a]++;
    }
  }
  // locations by endpoint node, with the cost of reaching it
  std::vector<int> lptr(n_nodes + 1, 0);
  for (int j = 0; j < m; ++j) {
    lptr[lu[j]]++;
    lptr[lv[j]]++;
  }
  for (int a = 0; a < n_nodes; ++a) lptr[a + 1] += lptr[a];
  std::vector<int> lj(2 * m);
  std::vector<double> lc(2 * m);
  {
    std::vector<int> fill(lptr.begin(), lptr.end() - 1);
    for (int j = 0; j < m; ++j) {
      int a = lu[j] - 1, b = lv[j] - 1;
      lj[fill[a]] = j;
      lc[fill[a]] = cu[j];
      fill[a]++;
      lj[fill[b]] = j;
      lc[fill[b]] = cv[j];
      fill[b]++;
    }
  }
  std::vector<double> best(m, INF);
  std::vector<int> touched;
  std::vector<double> S(k);
  for (int i = 0; i < m; ++i) {
    for (int e = 0; e < 2; ++e) {
      int a = (e == 0 ? lu[i] : lv[i]) - 1;
      double ca = e == 0 ? cu[i] : cv[i];
      for (int t = ptr[a]; t < ptr[a + 1]; ++t) {
        int b = pb[t];
        double dab = ca + pd[t];
        for (int r = lptr[b]; r < lptr[b + 1]; ++r) {
          int j = lj[r];
          double cand = dab + lc[r];
          if (best[j] == INF) touched.push_back(j);
          if (cand < best[j]) best[j] = cand;
        }
      }
    }
    std::fill(S.begin(), S.end(), 0.0);
    for (size_t t = 0; t < touched.size(); ++t) {
      int j = touched[t];
      double d = best[j];
      best[j] = INF;
      if (lu[i] == lu[j] && lv[i] == lv[j]) d = std::min(d, std::abs(cu[i] - cu[j]));
      d += dep[i] + dep[j];
      if (i == j) d = 0.0;
      double kij = 1.0 - d / cutoff;
      if (!(kij > 0)) continue;
      for (int c = 0; c < k; ++c) S[c] += kij * P(j, c);
    }
    touched.clear();
    for (int a = 0; a < k; ++a) {
      double pa = P(i, a);
      if (pa == 0.0) continue;
      for (int b = 0; b < k; ++b) out(a, b) += pa * S[b];
    }
  }
  return out;
}
