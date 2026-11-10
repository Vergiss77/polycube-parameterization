#include "pcp/structure/patch_graph.h"

#include <algorithm>
#include <map>
#include <set>

#include <Eigen/Geometry>

namespace pcp {

namespace {

long long euler_of_faces(const MatX3i& B, const std::vector<int>& faces) {
  std::set<int> verts;
  std::set<std::pair<int, int>> edges;
  for (int f : faces)
    for (int k = 0; k < 3; ++k) {
      int a = B(f, k), b = B(f, (k + 1) % 3);
      verts.insert(a);
      edges.insert({std::min(a, b), std::max(a, b)});
    }
  return static_cast<long long>(verts.size()) - static_cast<long long>(edges.size()) +
         static_cast<long long>(faces.size());
}

}  // namespace

PatchGraph build_patch_graph(const MatX3d& X, const MatX3i& B, const BoundaryTopology& topo,
                             const Eigen::VectorXi& labels) {
  PatchGraph g;
  g.patches = compute_patches(topo, labels);
  const size_t np = g.patches.label.size();
  g.height.assign(np, 0.0);
  g.area.assign(np, 0.0);
  for (size_t p = 0; p < np; ++p) {
    const int axis = label_axis(g.patches.label[p]);
    double h = 0.0, a = 0.0;
    for (int f : g.patches.faces[p]) {
      Eigen::Vector3d x0 = X.row(B(f, 0)), x1 = X.row(B(f, 1)), x2 = X.row(B(f, 2));
      double af = 0.5 * (x1 - x0).cross(x2 - x0).norm();
      h += af * (x0(axis) + x1(axis) + x2(axis)) / 3.0;
      a += af;
    }
    g.area[p] = a;
    g.height[p] = a > 0.0 ? h / a : 0.0;
  }

  // Patches incident to each vertex.
  std::map<int, std::vector<int>> vpatch;
  for (Index f = 0; f < B.rows(); ++f)
    for (int k = 0; k < 3; ++k) vpatch[B(f, k)].push_back(g.patches.face_patch(f));
  std::vector<char> is_corner(static_cast<size_t>(X.rows()), 0);
  for (auto& [v, ps] : vpatch) {
    std::sort(ps.begin(), ps.end());
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
    if (ps.size() >= 3) {
      is_corner[static_cast<size_t>(v)] = 1;
      g.corners.push_back(v);
      g.corner_patches.push_back(ps);
    }
  }

  // Sharp edges and vertex adjacency along them.
  struct SharpEdge {
    int a, b, pa, pb;
  };
  std::vector<SharpEdge> sharp;
  std::map<int, std::vector<int>> vsharp;  // vertex -> sharp edge ids
  for (const FacePair& p : topo.pairs) {
    int pa = g.patches.face_patch(p.f0), pb = g.patches.face_patch(p.f1);
    if (pa == pb) continue;
    int id = static_cast<int>(sharp.size());
    sharp.push_back({p.v0, p.v1, std::min(pa, pb), std::max(pa, pb)});
    vsharp[p.v0].push_back(id);
    vsharp[p.v1].push_back(id);
  }
  std::vector<char> used(sharp.size(), 0);

  auto trace = [&](int start_vertex, int first_edge) {
    Chain c;
    c.patch_a = sharp[static_cast<size_t>(first_edge)].pa;
    c.patch_b = sharp[static_cast<size_t>(first_edge)].pb;
    c.vertices.push_back(start_vertex);
    int v = start_vertex, e = first_edge;
    while (true) {
      used[static_cast<size_t>(e)] = 1;
      const SharpEdge& se = sharp[static_cast<size_t>(e)];
      int w = se.a == v ? se.b : se.a;
      c.length += (X.row(v) - X.row(w)).norm();
      if (w == start_vertex) {
        c.closed = true;
        break;
      }
      c.vertices.push_back(w);
      if (is_corner[static_cast<size_t>(w)]) break;
      int next = -1;
      for (int e2 : vsharp[w])
        if (!used[static_cast<size_t>(e2)]) {
          next = e2;
          break;
        }
      if (next < 0) break;
      v = w;
      e = next;
    }
    return c;
  };

  for (int v : g.corners)
    for (int e : vsharp[v])
      if (!used[static_cast<size_t>(e)]) g.chains.push_back(trace(v, e));
  for (size_t e = 0; e < sharp.size(); ++e)
    if (!used[e]) g.chains.push_back(trace(sharp[e].a, static_cast<int>(e)));

  // Axis, corners and end patches of each chain.
  for (Chain& c : g.chains) {
    int axa = label_axis(g.patches.label[static_cast<size_t>(c.patch_a)]);
    int axb = label_axis(g.patches.label[static_cast<size_t>(c.patch_b)]);
    c.axis = axa != axb ? 3 - axa - axb : -1;
    if (c.closed) continue;
    c.corner0 = c.vertices.front();
    c.corner1 = c.vertices.back();
    if (c.axis < 0) continue;
    auto end_patch = [&](int corner) {
      int best = -1;
      double best_d = 0.0;
      for (int p : vpatch[corner]) {
        if (label_axis(g.patches.label[static_cast<size_t>(p)]) != c.axis) continue;
        double d = std::abs(g.height[static_cast<size_t>(p)] - X(corner, c.axis));
        if (best < 0 || d < best_d) {
          best = p;
          best_d = d;
        }
      }
      return best;
    };
    int p0 = end_patch(c.corner0), p1 = end_patch(c.corner1);
    if (X(c.corner0, c.axis) <= X(c.corner1, c.axis)) {
      c.lower_patch = p0;
      c.upper_patch = p1;
    } else {
      c.lower_patch = p1;
      c.upper_patch = p0;
    }
  }
  return g;
}

nlohmann::json PatchGraph::to_json(const MatX3d& X) const {
  static const char* names[6] = {"+X", "-X", "+Y", "-Y", "+Z", "-Z"};
  nlohmann::json jp = nlohmann::json::array();
  for (size_t p = 0; p < patches.label.size(); ++p) {
    int l = patches.label[p];
    jp.push_back({{"id", p},
                  {"label", names[l]},
                  {"axis", label_axis(l)},
                  {"sign", label_sign(l)},
                  {"height", height[p]},
                  {"area", area[p]},
                  {"neighbors", patches.neighbors[p]},
                  {"faces", patches.faces[p]}});
  }
  nlohmann::json jc = nlohmann::json::array();
  for (size_t i = 0; i < chains.size(); ++i) {
    const Chain& c = chains[i];
    jc.push_back({{"id", i},
                  {"patches", {c.patch_a, c.patch_b}},
                  {"axis", c.axis},
                  {"closed", c.closed},
                  {"corners", {c.corner0, c.corner1}},
                  {"lower_patch", c.lower_patch},
                  {"upper_patch", c.upper_patch},
                  {"length", c.length},
                  {"vertices", c.vertices}});
  }
  nlohmann::json jk = nlohmann::json::array();
  for (size_t i = 0; i < corners.size(); ++i) {
    int v = corners[i];
    jk.push_back({{"vertex", v},
                  {"position", {X(v, 0), X(v, 1), X(v, 2)}},
                  {"patches", corner_patches[i]}});
  }
  return {{"patches", jp}, {"chains", jc}, {"corners", jk}};
}

nlohmann::json TopologyCheck::to_json() const {
  return {{"ok", ok},
          {"violations", violations},
          {"euler_surface", euler_surface},
          {"euler_complex", euler_complex},
          {"non_disk_patches", non_disk_patches}};
}

TopologyCheck check_polycube_topology(const MatX3i& B, const PatchGraph& g) {
  TopologyCheck r;
  const auto& P = g.patches;
  for (size_t p = 0; p < P.label.size(); ++p) {
    if (P.neighbors[p].size() < 3)
      r.violations.push_back("patch " + std::to_string(p) + " has only " +
                             std::to_string(P.neighbors[p].size()) + " neighbors");
    for (int q : P.neighbors[p])
      if (static_cast<size_t>(q) > p && label_axis(P.label[p]) == label_axis(P.label[static_cast<size_t>(q)]))
        r.violations.push_back("adjacent patches " + std::to_string(p) + " and " +
                               std::to_string(q) + " lie on the same axis");
  }
  long long chain_chi = 0;
  for (size_t i = 0; i < g.chains.size(); ++i) {
    if (g.chains[i].axis < 0)
      r.violations.push_back("chain " + std::to_string(i) + " separates patches on the same axis");
    chain_chi += g.chains[i].closed ? 0 : 1;
  }
  long long patch_chi = 0;
  for (const auto& faces : P.faces) {
    long long chi = euler_of_faces(B, faces);
    patch_chi += chi;
    if (chi != 1) ++r.non_disk_patches;
  }
  std::vector<int> all(static_cast<size_t>(B.rows()));
  for (size_t f = 0; f < all.size(); ++f) all[f] = static_cast<int>(f);
  r.euler_surface = euler_of_faces(B, all);
  r.euler_complex = patch_chi - chain_chi + static_cast<long long>(g.corners.size());
  if (r.euler_complex != r.euler_surface)
    r.violations.push_back("Euler characteristic of the patch structure (" +
                           std::to_string(r.euler_complex) + ") differs from the surface (" +
                           std::to_string(r.euler_surface) + ")");
  r.ok = r.violations.empty();
  return r;
}

}  // namespace pcp
