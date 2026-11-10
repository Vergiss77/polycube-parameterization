#include "pcp/geom/boundary_topology.h"

#include <algorithm>
#include <tuple>

namespace pcp {

BoundaryTopology build_boundary_topology(const MatX3i& B) {
  const Index nf = B.rows();
  struct HalfEdge {
    int a, b, face, k;
  };
  std::vector<HalfEdge> he;
  he.reserve(static_cast<size_t>(3 * nf));
  for (Index f = 0; f < nf; ++f)
    for (int k = 0; k < 3; ++k) {
      int a = B(f, k), b = B(f, (k + 1) % 3);
      he.push_back({std::min(a, b), std::max(a, b), static_cast<int>(f), k});
    }
  std::sort(he.begin(), he.end(), [](const HalfEdge& x, const HalfEdge& y) {
    return std::tie(x.a, x.b, x.face) < std::tie(y.a, y.b, y.face);
  });

  BoundaryTopology t;
  t.face_pairs.assign(static_cast<size_t>(nf), {-1, -1, -1});
  t.face_nbrs.assign(static_cast<size_t>(nf), {-1, -1, -1});
  for (size_t i = 0; i < he.size(); i += 2) {
    if (i + 1 >= he.size() || he[i].a != he[i + 1].a || he[i].b != he[i + 1].b ||
        (i + 2 < he.size() && he[i + 2].a == he[i].a && he[i + 2].b == he[i].b))
      throw Error("boundary topology: surface is not closed and edge-manifold");
    const HalfEdge &h0 = he[i], &h1 = he[i + 1];
    FacePair p;
    p.f0 = h0.face;
    p.f1 = h1.face;
    p.v0 = B(h0.face, h0.k);
    p.v1 = B(h0.face, (h0.k + 1) % 3);
    p.nodes = {p.v0, p.v1, B(h0.face, (h0.k + 2) % 3), B(h1.face, (h1.k + 2) % 3)};
    int id = static_cast<int>(t.pairs.size());
    t.pairs.push_back(p);
    t.face_pairs[static_cast<size_t>(h0.face)][static_cast<size_t>(h0.k)] = id;
    t.face_pairs[static_cast<size_t>(h1.face)][static_cast<size_t>(h1.k)] = id;
    t.face_nbrs[static_cast<size_t>(h0.face)][static_cast<size_t>(h0.k)] = h1.face;
    t.face_nbrs[static_cast<size_t>(h1.face)][static_cast<size_t>(h1.k)] = h0.face;
  }
  return t;
}

}  // namespace pcp
