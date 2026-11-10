#include "pcp/structure/height_simplify.h"

#include <igl/AABB.h>
#include <igl/Hit.h>

#include "pcp/geom/boundary_topology.h"

#include <algorithm>
#include <set>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>

namespace pcp {

nlohmann::json SimplifyParams::to_json() const {
  return {{"enabled", enabled}, {"max_stretch", max_stretch}};
}

void SimplifyParams::update_from_json(const nlohmann::json& j) {
  const nlohmann::json known = to_json();
  for (auto it = j.begin(); it != j.end(); ++it)
    if (!known.contains(it.key())) throw Error("unknown simplify parameter '" + it.key() + "'");
  if (j.contains("enabled")) enabled = j.at("enabled").get<bool>();
  if (j.contains("max_stretch")) max_stretch = j.at("max_stretch").get<double>();
  if (!(max_stretch >= 1.0)) throw Error("simplify.max_stretch must be >= 1");
}

HeightGraph build_height_graph(const PatchGraph& g, const std::vector<double>& heights) {
  const size_t np = g.patches.label.size();
  if (heights.size() != np) throw Error("build_height_graph: one height per patch is required");
  HeightGraph hg;
  hg.n_patches = static_cast<int>(np);
  if (np == 0) return hg;

  // Heights equal up to round-off (e.g. area-weighted means of identical coordinates) form one
  // group. Within an axis, patches are sorted by height and consecutive ones within `tol` join.
  double lo = heights[0], hi = heights[0];
  for (double h : heights) lo = std::min(lo, h), hi = std::max(hi, h);
  const double tol = 1e-12 * std::max({hi - lo, std::abs(lo), std::abs(hi), 1e-300});

  std::vector<int> group_of(np, -1);
  std::vector<std::vector<int>> clusters;
  for (int axis = 0; axis < 3; ++axis) {
    std::vector<int> ps;
    for (size_t p = 0; p < np; ++p)
      if (g.patches.label[p] / 2 == axis) ps.push_back(static_cast<int>(p));
    std::stable_sort(ps.begin(), ps.end(), [&](int a, int b) {
      return heights[static_cast<size_t>(a)] < heights[static_cast<size_t>(b)];
    });
    for (size_t i = 0; i < ps.size(); ++i) {
      if (i == 0 || heights[static_cast<size_t>(ps[i])] - heights[static_cast<size_t>(ps[i - 1])] > tol)
        clusters.emplace_back();
      clusters.back().push_back(ps[i]);
    }
  }
  // Number groups by their smallest patch id (independent of the sort above).
  for (auto& c : clusters) std::sort(c.begin(), c.end());
  std::sort(clusters.begin(), clusters.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; });
  for (const auto& c : clusters) {
    HeightGroup grp;
    grp.axis = g.patches.label[static_cast<size_t>(c[0])] / 2;
    grp.patches = c;
    grp.height = heights[static_cast<size_t>(c[0])];
    const int id = static_cast<int>(hg.groups.size());
    for (int p : c) group_of[static_cast<size_t>(p)] = id;
    hg.groups.push_back(std::move(grp));
  }

  for (const Chain& c : g.chains) {
    if (c.lower_patch < 0 || c.upper_patch < 0) {
      ++hg.chains_without_ends;
      continue;
    }
    const auto lo_p = static_cast<size_t>(c.lower_patch), up_p = static_cast<size_t>(c.upper_patch);
    const double rest = heights[up_p] - heights[lo_p];
    if (!(rest > 0.0) || group_of[lo_p] == group_of[up_p]) {
      throw Error("build_height_graph: chain between patches " + std::to_string(c.lower_patch) +
                  " and " + std::to_string(c.upper_patch) + " is degenerate (height difference " +
                  std::to_string(rest) + ")");
    }
    hg.edges.push_back({group_of[lo_p], group_of[up_p], rest});
  }
  return hg;
}

int add_separation_edges(HeightGraph& hg, const PatchGraph& g, const std::vector<double>& heights,
                         const MatX3d& X, const MatX3i& B) {
  std::vector<int> group_of(static_cast<size_t>(hg.n_patches), -1);
  for (size_t i = 0; i < hg.groups.size(); ++i)
    for (int p : hg.groups[i].patches) group_of[static_cast<size_t>(p)] = static_cast<int>(i);
  std::set<std::pair<int, int>> existing;  // (lower patch, upper patch) already added
  const Eigen::MatrixXd V = X;
  const Eigen::MatrixXi F = B;
  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V, F);
  const double diag = (X.colwise().maxCoeff() - X.colwise().minCoeff()).norm();
  const double min_t = 1e-9 * diag;
  auto label_of = [&](Index f) {
    return g.patches.label[static_cast<size_t>(g.patches.face_patch(f))];
  };
  std::vector<igl::Hit<double>> hits;
  int added = 0;
  for (Index f = 0; f < B.rows(); ++f) {
    const int p = g.patches.face_patch(f);
    const int lab = label_of(f), axis = lab / 2;
    const Eigen::RowVector3d c = (X.row(B(f, 0)) + X.row(B(f, 1)) + X.row(B(f, 2))) / 3.0;
    for (double dirsign : {-1.0, 1.0}) {
      Eigen::RowVector3d dir = Eigen::RowVector3d::Zero();
      dir(axis) = dirsign;
      hits.clear();
      tree.intersect_ray(V, F, c, dir, hits);
      std::sort(hits.begin(), hits.end(),
                [](const igl::Hit<double>& a, const igl::Hit<double>& b) { return a.t < b.t; });
      for (const auto& h : hits) {
        if (!(h.t > min_t)) continue;
        const int hl = label_of(h.id);
        if (hl / 2 != axis) continue;           // grazing a perpendicular triangle
        if (hl == lab) break;                   // same orientation: not a facing pair
        const int q = g.patches.face_patch(h.id);
        const double hp = heights[static_cast<size_t>(p)], hq = heights[static_cast<size_t>(q)];
        if (hp == hq || group_of[static_cast<size_t>(p)] == group_of[static_cast<size_t>(q)]) break;
        const int lo = hp < hq ? p : q, up = hp < hq ? q : p;
        if (existing.insert({lo, up}).second) {
          HeightEdge e;
          e.from = group_of[static_cast<size_t>(lo)];
          e.to = group_of[static_cast<size_t>(up)];
          e.rest = heights[static_cast<size_t>(up)] - heights[static_cast<size_t>(lo)];
          e.separation = true;
          hg.edges.push_back(e);
          ++added;
        }
        break;
      }
    }
  }
  hg.separation_edges += added;
  return added;
}

nlohmann::json HeightGraph::to_json() const {
  nlohmann::json gj = nlohmann::json::array(), ej = nlohmann::json::array();
  for (const HeightGroup& grp : groups)
    gj.push_back({{"axis", grp.axis}, {"height", grp.height}, {"patches", grp.patches}});
  for (const HeightEdge& e : edges)
    ej.push_back({{"from", e.from}, {"to", e.to}, {"rest", e.rest}, {"separation", e.separation}});
  return {{"n_patches", n_patches},
          {"groups", gj},
          {"edges", ej},
          {"chains_without_ends", chains_without_ends},
          {"separation_edges", separation_edges}};
}

HeightGraph HeightGraph::from_json(const nlohmann::json& j) {
  HeightGraph g;
  g.n_patches = j.at("n_patches").get<int>();
  for (const auto& x : j.at("groups")) {
    HeightGroup grp;
    grp.axis = x.at("axis").get<int>();
    grp.height = x.at("height").get<double>();
    grp.patches = x.at("patches").get<std::vector<int>>();
    g.groups.push_back(std::move(grp));
  }
  for (const auto& x : j.at("edges")) {
    HeightEdge e;
    e.from = x.at("from").get<int>();
    e.to = x.at("to").get<int>();
    e.rest = x.at("rest").get<double>();
    e.separation = x.at("separation").get<bool>();
    g.edges.push_back(e);
  }
  g.chains_without_ends = j.at("chains_without_ends").get<int>();
  g.separation_edges = j.at("separation_edges").get<int>();
  return g;
}

int add_width_edges(HeightGraph& hg, const PatchGraph& g, const std::vector<double>& heights,
                    const MatX3d& X, const MatX3i& B) {
  std::vector<int> group_of(static_cast<size_t>(hg.n_patches), -1);
  for (size_t i = 0; i < hg.groups.size(); ++i)
    for (int p : hg.groups[i].patches) group_of[static_cast<size_t>(p)] = static_cast<int>(i);
  std::set<std::pair<int, int>> existing;
  const BoundaryTopology topo = build_boundary_topology(B);
  const size_t np = g.patches.label.size();
  struct BEdge {
    Eigen::Vector3d x0, x1;
    int nbr;
  };
  std::vector<std::vector<BEdge>> bedges(np);
  for (Index f = 0; f < B.rows(); ++f) {
    const int p = g.patches.face_patch(f);
    for (int k = 0; k < 3; ++k) {
      const int nf = topo.face_nbrs[static_cast<size_t>(f)][static_cast<size_t>(k)];
      const int q = g.patches.face_patch(nf);
      if (q == p) continue;
      bedges[static_cast<size_t>(p)].push_back(
          {X.row(B(f, k)).transpose(), X.row(B(f, (k + 1) % 3)).transpose(), q});
    }
  }
  int added = 0;
  for (Index f = 0; f < B.rows(); ++f) {
    const int p = g.patches.face_patch(f);
    const int b = g.patches.label[static_cast<size_t>(p)] / 2;
    const Eigen::Vector3d c = (X.row(B(f, 0)) + X.row(B(f, 1)) + X.row(B(f, 2))).transpose() / 3.0;
    for (int a : {(b + 1) % 3, (b + 2) % 3}) {
      const int o = 3 - a - b;  // the other in-plane axis
      double best_lo = -std::numeric_limits<double>::infinity(), best_hi = std::numeric_limits<double>::infinity();
      int q_lo = -1, q_hi = -1;
      for (const BEdge& e : bedges[static_cast<size_t>(p)]) {
        const double o0 = e.x0(o), o1 = e.x1(o);
        if ((o0 <= c(o)) == (o1 <= c(o))) continue;  // does not cross the ray's line
        const double t = (c(o) - o0) / (o1 - o0);
        const double xa = e.x0(a) + t * (e.x1(a) - e.x0(a));
        if (xa > c(a) && xa < best_hi) best_hi = xa, q_hi = e.nbr;
        if (xa < c(a) && xa > best_lo) best_lo = xa, q_lo = e.nbr;
      }
      if (q_lo < 0 || q_hi < 0) continue;
      if (g.patches.label[static_cast<size_t>(q_lo)] / 2 != a || g.patches.label[static_cast<size_t>(q_hi)] / 2 != a)
        continue;
      const double hl = heights[static_cast<size_t>(q_lo)], hh = heights[static_cast<size_t>(q_hi)];
      if (!(hh > hl) || group_of[static_cast<size_t>(q_lo)] == group_of[static_cast<size_t>(q_hi)]) continue;
      if (existing.insert({q_lo, q_hi}).second) {
        HeightEdge e;
        e.from = group_of[static_cast<size_t>(q_lo)];
        e.to = group_of[static_cast<size_t>(q_hi)];
        e.rest = hh - hl;
        e.width = true;
        hg.edges.push_back(e);
        ++added;
      }
    }
  }
  hg.width_edges += added;
  return added;
}

HeightState HeightState::initial(const HeightGraph& g) {
  HeightState st;
  st.rep.resize(g.groups.size());
  std::iota(st.rep.begin(), st.rep.end(), 0);
  st.height.resize(g.groups.size());
  for (size_t i = 0; i < g.groups.size(); ++i) st.height[i] = g.groups[i].height;
  return st;
}

namespace {

double spring(double r, double rest) { return rest / r + r / rest; }

}  // namespace

PairCost merge_cost(const HeightGraph& g, const HeightState& st, int a, int b, double max_stretch) {
  PairCost pc;
  const double s = max_stretch;
  // Incident edges: sign +1 when the merged group is the upper end (r = h - other),
  // -1 when it is the lower end (r = other - h).
  struct Inc {
    double other, rest;
    int sign;
  };
  std::vector<Inc> inc;
  double lo = -std::numeric_limits<double>::infinity(), hi = std::numeric_limits<double>::infinity();
  double old_sum = 0.0;
  for (const HeightEdge& e : g.edges) {
    const int f = st.rep[static_cast<size_t>(e.from)], t = st.rep[static_cast<size_t>(e.to)];
    const bool fin = (f == a || f == b), tin = (t == a || t == b);
    if (fin && tin) return pc;  // an edge inside the pair: the groups are ordered
    if (!fin && !tin) continue;
    const double hf = st.height[static_cast<size_t>(f)], ht = st.height[static_cast<size_t>(t)];
    old_sum += spring(ht - hf, e.rest);
    if (tin) {  // r = h - hf in [rest/s, rest*s]
      inc.push_back({hf, e.rest, +1});
      lo = std::max(lo, hf + e.rest / s);
      hi = std::min(hi, hf + e.rest * s);
    } else {  // r = ht - h
      inc.push_back({ht, e.rest, -1});
      lo = std::max(lo, ht - e.rest * s);
      hi = std::min(hi, ht - e.rest / s);
    }
  }
  if (inc.empty()) {  // unconstrained pair: no distortion change, keep the mean height
    pc.feasible = true;
    pc.height = 0.5 * (st.height[static_cast<size_t>(a)] + st.height[static_cast<size_t>(b)]);
    pc.cost = 0.0;
    return pc;
  }
  if (lo > hi) return pc;

  auto length = [](const Inc& e, double h) { return e.sign > 0 ? h - e.other : e.other - h; };
  auto derivative = [&](double h) {
    double d = 0.0;
    for (const Inc& e : inc) {
      const double r = length(e, h);
      d += e.sign * (-e.rest / (r * r) + 1.0 / e.rest);
    }
    return d;
  };
  // Convex in h on [lo, hi]: bisection on the derivative with a fixed number of steps.
  double h;
  if (derivative(lo) >= 0.0) {
    h = lo;
  } else if (derivative(hi) <= 0.0) {
    h = hi;
  } else {
    double l = lo, u = hi;
    for (int it = 0; it < 100; ++it) {
      const double m = 0.5 * (l + u);
      if (derivative(m) > 0.0) u = m;
      else l = m;
    }
    h = 0.5 * (l + u);
  }
  double new_sum = 0.0;
  for (const Inc& e : inc) new_sum += spring(length(e, h), e.rest);
  pc.feasible = true;
  pc.height = h;
  pc.cost = new_sum - old_sum;
  return pc;
}

namespace {

std::array<int, 3> count_cuts(const HeightGraph& g, const HeightState& st) {
  std::array<int, 3> c{};
  for (size_t i = 0; i < g.groups.size(); ++i)
    if (st.rep[i] == static_cast<int>(i)) ++c[static_cast<size_t>(g.groups[i].axis)];
  return c;
}

// reach[i][j] != 0 iff there is a directed path from representative i to representative j.
std::vector<std::vector<char>> reachability(const HeightGraph& g, const HeightState& st) {
  const size_t n = g.groups.size();
  std::vector<std::vector<int>> adj(n);
  for (const HeightEdge& e : g.edges)
    adj[static_cast<size_t>(st.rep[static_cast<size_t>(e.from)])].push_back(st.rep[static_cast<size_t>(e.to)]);
  std::vector<std::vector<char>> reach(n, std::vector<char>(n, 0));
  std::vector<int> stack;
  for (size_t i = 0; i < n; ++i) {
    if (st.rep[i] != static_cast<int>(i)) continue;
    std::vector<char>& r = reach[i];
    stack.assign(adj[i].begin(), adj[i].end());
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      if (r[static_cast<size_t>(v)]) continue;
      r[static_cast<size_t>(v)] = 1;
      for (int w : adj[static_cast<size_t>(v)])
        if (!r[static_cast<size_t>(w)]) stack.push_back(w);
    }
  }
  return reach;
}

}  // namespace

SimplifyResult apply_merges(const HeightGraph& g, const std::vector<MergeRecord>& merges, size_t k) {
  k = std::min(k, merges.size());
  HeightState st = HeightState::initial(g);
  SimplifyResult res;
  res.cuts_before = count_cuts(g, st);
  for (const HeightEdge& e : g.edges) {
    if (e.separation) ++res.separation_edges;
    else if (e.width) ++res.width_edges;
    else ++res.chain_edges;
  }
  const int n = static_cast<int>(g.groups.size());
  for (size_t m = 0; m < k; ++m) {
    const MergeRecord& mr = merges[m];
    for (int i = 0; i < n; ++i)
      if (st.rep[static_cast<size_t>(i)] == mr.b) st.rep[static_cast<size_t>(i)] = mr.a;
    st.height[static_cast<size_t>(mr.a)] = mr.height;
  }
  res.merges.assign(merges.begin(), merges.begin() + static_cast<std::ptrdiff_t>(k));
  res.merges_proposed = static_cast<int>(merges.size());
  res.cuts_after = count_cuts(g, st);
  res.group_rep = st.rep;
  res.group_height.resize(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
    res.group_height[static_cast<size_t>(i)] = st.height[static_cast<size_t>(st.rep[static_cast<size_t>(i)])];
  res.group_of_patch.assign(static_cast<size_t>(g.n_patches), -1);
  res.patch_height.assign(static_cast<size_t>(g.n_patches), 0.0);
  for (int i = 0; i < n; ++i)
    for (int p : g.groups[static_cast<size_t>(i)].patches) {
      res.group_of_patch[static_cast<size_t>(p)] = st.rep[static_cast<size_t>(i)];
      res.patch_height[static_cast<size_t>(p)] = res.group_height[static_cast<size_t>(i)];
    }
  for (size_t e_i = 0; e_i < g.edges.size(); ++e_i) {
    const HeightEdge& e = g.edges[e_i];
    const double r = res.group_height[static_cast<size_t>(e.to)] - res.group_height[static_cast<size_t>(e.from)];
    const double q = r / e.rest;
    if (e_i == 0) res.min_ratio = res.max_ratio = q;
    res.min_ratio = std::min(res.min_ratio, q);
    res.max_ratio = std::max(res.max_ratio, q);
  }
  return res;
}

SimplifyResult simplify_heights(const HeightGraph& g, double max_stretch) {
  if (!(max_stretch >= 1.0)) throw Error("simplify_heights: max_stretch must be >= 1");
  HeightState st = HeightState::initial(g);
  std::vector<MergeRecord> merges;
  const int n = static_cast<int>(g.groups.size());
  while (true) {
    const auto reach = reachability(g, st);
    int best_a = -1, best_b = -1;
    PairCost best;
    for (int a = 0; a < n; ++a) {
      if (st.rep[static_cast<size_t>(a)] != a) continue;
      for (int b = a + 1; b < n; ++b) {
        if (st.rep[static_cast<size_t>(b)] != b) continue;
        if (g.groups[static_cast<size_t>(a)].axis != g.groups[static_cast<size_t>(b)].axis) continue;
        if (reach[static_cast<size_t>(a)][static_cast<size_t>(b)] ||
            reach[static_cast<size_t>(b)][static_cast<size_t>(a)])
          continue;
        const PairCost pc = merge_cost(g, st, a, b, max_stretch);
        if (!pc.feasible) continue;
        if (best_a < 0 || pc.cost < best.cost) {  // ties keep the lexicographically first pair
          best = pc;
          best_a = a;
          best_b = b;
        }
      }
    }
    if (best_a < 0) break;
    for (int i = 0; i < n; ++i)
      if (st.rep[static_cast<size_t>(i)] == best_b) st.rep[static_cast<size_t>(i)] = best_a;
    st.height[static_cast<size_t>(best_a)] = best.height;
    merges.push_back({best_a, best_b, best.height, best.cost});
  }
  return apply_merges(g, merges, merges.size());
}

nlohmann::json SimplifyResult::to_json() const {
  nlohmann::json merges_j = nlohmann::json::array();
  for (const MergeRecord& m : merges)
    merges_j.push_back({{"a", m.a}, {"b", m.b}, {"height", m.height}, {"cost", m.cost}});
  return {{"cuts_before", {{"x", cuts_before[0]}, {"y", cuts_before[1]}, {"z", cuts_before[2]},
                           {"total", cuts_total_before()}}},
          {"cuts_after", {{"x", cuts_after[0]}, {"y", cuts_after[1]}, {"z", cuts_after[2]},
                          {"total", cuts_total_after()}}},
          {"merges", merges_j},
          {"merges_proposed", merges_proposed},
          {"chain_edges", chain_edges},
          {"separation_edges", separation_edges},
          {"width_edges", width_edges},
          {"min_ratio", min_ratio},
          {"max_ratio", max_ratio}};
}

SimplifiedPolycube remap_with_rollback(const MatX3d& rest, const MatX3d& X_prev, const MatX4i& T,
                                       const MatX3i& B, const PatchGraph& g, const HeightGraph& hg,
                                       const SimplifyParams& sp, FinalMapParams fp) {
  SimplifiedPolycube out;
  out.graph = hg;
  const auto t0 = std::chrono::steady_clock::now();
  const SimplifyResult greedy = simplify_heights(hg, sp.max_stretch);
  out.seconds_simplify = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  fp.heights_fixed = true;
  FinalMapParams probe = fp;
  probe.untangle = false;
  auto inverted_with = [&](size_t k) {
    ++out.probe_maps;
    return compute_final_map(rest, X_prev, T, B, g, apply_merges(hg, greedy.merges, k).patch_height, probe)
        .inverted_tets;
  };
  const size_t m = greedy.merges.size();
  size_t keep = m;
  out.inverted_reference = inverted_with(0);
  if (m > 0) {
    Index n_full = inverted_with(m);
    if (n_full <= out.inverted_reference) {
      out.inverted_kept = n_full;
    } else {
      // Longest prefix whose map has no more inversions than the reference (bisection).
      size_t lo = 0, hi = m;
      Index n_lo = out.inverted_reference;
      while (hi - lo > 1) {
        const size_t mid = lo + (hi - lo) / 2;
        const Index n_mid = inverted_with(mid);
        if (n_mid <= out.inverted_reference) lo = mid, n_lo = n_mid;
        else hi = mid;
      }
      keep = lo;
      out.inverted_kept = n_lo;
    }
  } else {
    out.inverted_kept = out.inverted_reference;
  }
  out.simplify = apply_merges(hg, greedy.merges, keep);
  out.map = compute_final_map(rest, X_prev, T, B, g, out.simplify.patch_height, fp);

  // The structure of the re-mapped polycube, recomputed from its geometry.
  BoundaryTopology topo = build_boundary_topology(B);
  Eigen::VectorXi labels(B.rows());
  for (Index f = 0; f < B.rows(); ++f)
    labels(f) = g.patches.label[static_cast<size_t>(g.patches.face_patch(f))];
  const PatchGraph g2 = build_patch_graph(out.map.X, B, topo, labels);
  out.topology = check_polycube_topology(B, g2);
  if (g2.patches.label.size() != g.patches.label.size() || g2.corners.size() != g.corners.size() ||
      g2.chains.size() != g.chains.size()) {
    out.topology.ok = false;
    out.topology.violations.push_back("patch structure changed by the simplification");
  }
  return out;
}

SimplifiedPolycube simplify_and_remap(const MatX3d& rest, const MatX3d& X_prev, const MatX4i& T,
                                      const MatX3i& B, const PatchGraph& g,
                                      const std::vector<double>& heights,
                                      const SimplifyParams& sp, FinalMapParams fp) {
  const auto t0 = std::chrono::steady_clock::now();
  HeightGraph hg = build_height_graph(g, heights);
  add_separation_edges(hg, g, heights, X_prev, B);
  add_width_edges(hg, g, heights, X_prev, B);
  const double t_graph = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  SimplifiedPolycube out = remap_with_rollback(rest, X_prev, T, B, g, hg, sp, fp);
  out.seconds_simplify += t_graph;
  return out;
}

}  // namespace pcp
