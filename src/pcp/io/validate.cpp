#include "pcp/io/validate.h"

#include <algorithm>
#include <numeric>

#include <Eigen/Geometry>

namespace pcp {

namespace {

struct UnionFind {
  std::vector<int> p;
  explicit UnionFind(int n) : p(static_cast<size_t>(n)) { std::iota(p.begin(), p.end(), 0); }
  int find(int x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  }
  void unite(int a, int b) { p[find(a)] = find(b); }
};

}  // namespace

nlohmann::json SurfaceReport::to_json() const {
  return {{"ok", ok},
          {"vertices", n_vertices},
          {"faces", n_faces},
          {"edges", n_edges},
          {"boundary_edges", boundary_edges},
          {"nonmanifold_edges", nonmanifold_edges},
          {"degenerate_faces", degenerate_faces},
          {"unused_vertices", unused_vertices},
          {"components", components},
          {"euler_characteristic", euler},
          {"genus", genus},
          {"volume", volume},
          {"area", area},
          {"errors", errors}};
}

SurfaceReport validate_surface(const SurfaceMesh& mesh, double rel_area_eps) {
  SurfaceReport r;
  const Index nv = mesh.V.rows();
  const Index nf = mesh.F.rows();
  r.n_vertices = nv;
  r.n_faces = nf;
  if (nv == 0 || nf == 0) {
    r.errors.push_back("empty mesh");
    return r;
  }

  // Undirected edges sorted to count incident faces.
  std::vector<std::pair<int, int>> edges;
  edges.reserve(static_cast<size_t>(3 * nf));
  for (Index f = 0; f < nf; ++f)
    for (int k = 0; k < 3; ++k) {
      int a = mesh.F(f, k), b = mesh.F(f, (k + 1) % 3);
      edges.emplace_back(std::min(a, b), std::max(a, b));
    }
  std::sort(edges.begin(), edges.end());
  for (size_t i = 0; i < edges.size();) {
    size_t j = i;
    while (j < edges.size() && edges[j] == edges[i]) ++j;
    size_t cnt = j - i;
    ++r.n_edges;
    if (cnt == 1) ++r.boundary_edges;
    if (cnt >= 3) ++r.nonmanifold_edges;
    i = j;
  }

  // Degenerate faces and geometric totals.
  const double diag = (mesh.V.colwise().maxCoeff() - mesh.V.colwise().minCoeff()).norm();
  const double area_eps = rel_area_eps * diag * diag;
  for (Index f = 0; f < nf; ++f) {
    Eigen::Vector3d a = mesh.V.row(mesh.F(f, 0)), b = mesh.V.row(mesh.F(f, 1)),
                    c = mesh.V.row(mesh.F(f, 2));
    Eigen::Vector3d cr = (b - a).cross(c - a);
    double area = 0.5 * cr.norm();
    r.area += area;
    r.volume += a.dot(b.cross(c)) / 6.0;
    bool repeated = mesh.F(f, 0) == mesh.F(f, 1) || mesh.F(f, 1) == mesh.F(f, 2) ||
                    mesh.F(f, 0) == mesh.F(f, 2);
    if (repeated || area <= area_eps) ++r.degenerate_faces;
  }

  // Used vertices and face-connected components.
  std::vector<char> used(static_cast<size_t>(nv), 0);
  UnionFind uf(static_cast<int>(nv));
  for (Index f = 0; f < nf; ++f) {
    for (int k = 0; k < 3; ++k) used[mesh.F(f, k)] = 1;
    uf.unite(mesh.F(f, 0), mesh.F(f, 1));
    uf.unite(mesh.F(f, 1), mesh.F(f, 2));
  }
  Index n_used = 0;
  std::vector<char> root_seen(static_cast<size_t>(nv), 0);
  for (Index v = 0; v < nv; ++v) {
    if (!used[v]) {
      ++r.unused_vertices;
      continue;
    }
    ++n_used;
    int root = uf.find(static_cast<int>(v));
    if (!root_seen[root]) {
      root_seen[root] = 1;
      ++r.components;
    }
  }
  r.euler = static_cast<long long>(n_used) - static_cast<long long>(r.n_edges) +
            static_cast<long long>(nf);

  if (r.boundary_edges > 0)
    r.errors.push_back("surface is not closed: " + std::to_string(r.boundary_edges) +
                       " boundary edges");
  if (r.nonmanifold_edges > 0)
    r.errors.push_back("surface is not edge-manifold: " + std::to_string(r.nonmanifold_edges) +
                       " non-manifold edges");
  if (r.degenerate_faces > 0)
    r.errors.push_back(std::to_string(r.degenerate_faces) + " degenerate faces");
  if (r.unused_vertices > 0)
    r.errors.push_back(std::to_string(r.unused_vertices) + " unused vertices");

  r.ok = r.errors.empty();
  if (r.boundary_edges == 0 && r.nonmanifold_edges == 0) {
    long long twice_genus = 2LL * r.components - r.euler;
    r.genus = static_cast<int>(twice_genus / 2);
  }
  if (r.ok && r.volume <= 0.0) {
    r.errors.push_back("surface encloses non-positive volume (inverted orientation?)");
    r.ok = false;
  }
  return r;
}

}  // namespace pcp
