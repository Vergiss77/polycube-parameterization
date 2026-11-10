#include <doctest/doctest.h>

#include <Eigen/Geometry>

#include "helpers.h"
#include "pcp/geom/normalize.h"
#include "pcp/io/off.h"
#include "pcp/structure/final_map.h"
#include "pcp/structure/labeling.h"
#include "pcp/structure/locator.h"
#include "pcp/structure/patch_graph.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {

// Face index of the box face (axis, side) containing a point, for injecting labels.
Eigen::VectorXi labels_of_box(const SurfaceMesh& s) { return initial_labels(s.V, s.F); }

}  // namespace

TEST_CASE("labels: box faces follow normals; zigzag face is relabeled") {
  SurfaceMesh s = testutil::make_box(4);
  BoundaryTopology topo = build_boundary_topology(s.F);
  Eigen::VectorXi L = labels_of_box(s);
  for (Index f = 0; f < s.F.rows(); ++f) {
    Eigen::Vector3d n = (Eigen::Vector3d(s.V.row(s.F(f, 1))) - Eigen::Vector3d(s.V.row(s.F(f, 0))))
                            .cross(Eigen::Vector3d(s.V.row(s.F(f, 2))) - Eigen::Vector3d(s.V.row(s.F(f, 0))));
    int axis;
    n.cwiseAbs().maxCoeff(&axis);
    CHECK(L(f) == make_label(axis, n(axis) > 0 ? 1 : -1));
  }
  // A face whose three neighbors share one label but which carries another label: zigzag.
  Eigen::VectorXi Z = L;
  int f = 0;
  while (true) {
    const auto& nb = topo.face_nbrs[static_cast<size_t>(f)];
    if (L(nb[0]) == L(f) && L(nb[1]) == L(f) && L(nb[2]) == L(f)) break;
    ++f;
  }
  Z(f) = (L(f) + 2) % 6;
  CHECK(remove_zigzags(topo, Z) >= 1);
  CHECK(Z == L);
}

TEST_CASE("cleanup merges an isolated patch and a wedge patch") {
  SurfaceMesh s = testutil::make_box(8);
  BoundaryTopology topo = build_boundary_topology(s.F);
  Eigen::VectorXi L = labels_of_box(s), ref = L;
  // Region of faces on the +Z face (z = 1) next to the crease x = 1 with x in [0.75, 1],
  // y in [0.375, 0.625]: relabel to +Y. It touches +Z and +X only -> degree 2.
  // Another region in the middle of the -Z face relabeled to +X -> isolated (degree 1).
  for (Index f = 0; f < s.F.rows(); ++f) {
    Eigen::Vector3d c = (s.V.row(s.F(f, 0)) + s.V.row(s.F(f, 1)) + s.V.row(s.F(f, 2))).transpose() / 3.0;
    if (std::abs(c.z() - 1.0) < 1e-9 && c.x() > 0.75 && c.y() > 0.375 && c.y() < 0.625) L(f) = make_label(1, 1);
    if (std::abs(c.z()) < 1e-9 && std::abs(c.x() - 0.5) < 0.2 && std::abs(c.y() - 0.5) < 0.2) L(f) = make_label(0, 1);
  }
  CleanupReport rep = cleanup_labels(s.V, s.F, topo, L);
  CHECK(L == ref);
  CHECK(rep.degree1_merged + rep.degree2_merged >= 2);
  CHECK(rep.degree2_merged >= 1);
  Patches P = compute_patches(topo, L);
  CHECK(P.label.size() == 6);
  for (const auto& nb : P.neighbors) CHECK(nb.size() >= 3);
}

TEST_CASE("patch graph of a box: 6 patches, 8 corners, 12 chains, topology ok") {
  SurfaceMesh s = testutil::make_box(4, 2.0, 1.0, 0.5);
  BoundaryTopology topo = build_boundary_topology(s.F);
  PatchGraph g = build_patch_graph(s.V, s.F, topo, initial_labels(s.V, s.F));
  CHECK(g.patches.label.size() == 6);
  CHECK(g.corners.size() == 8);
  CHECK(g.chains.size() == 12);
  for (const Chain& c : g.chains) {
    CHECK_FALSE(c.closed);
    CHECK(c.axis >= 0);
    CHECK(c.lower_patch >= 0);
    CHECK(c.upper_patch >= 0);
    CHECK(g.height[static_cast<size_t>(c.lower_patch)] < g.height[static_cast<size_t>(c.upper_patch)]);
  }
  TopologyCheck tc = check_polycube_topology(s.F, g);
  CHECK(tc.ok);
  CHECK(tc.euler_complex == 2);
  auto j = g.to_json(s.V);
  CHECK(j["patches"].size() == 6);
  CHECK(j["corners"].size() == 8);
}

TEST_CASE("topology check flags adjacent patches on the same axis") {
  SurfaceMesh s = testutil::make_box(4);
  BoundaryTopology topo = build_boundary_topology(s.F);
  Eigen::VectorXi L = initial_labels(s.V, s.F);
  // Relabel the whole +Y face as -X: it becomes adjacent to the +X/-X faces on the same axis.
  for (Index f = 0; f < s.F.rows(); ++f)
    if (L(f) == make_label(1, 1)) L(f) = make_label(0, -1);
  PatchGraph g = build_patch_graph(s.V, s.F, topo, L);
  TopologyCheck tc = check_polycube_topology(s.F, g);
  CHECK_FALSE(tc.ok);
}

TEST_CASE("final map: near identity on a box, exact planarity, fixed heights") {
  SurfaceMesh s = testutil::make_box(3, 1.0, 1.5, 0.7);
  TetMesh t = tetrahedralize(s, TetParams{});
  BoundaryTopology topo = build_boundary_topology(t.B);
  PatchGraph g = build_patch_graph(t.X, t.B, topo, initial_labels(t.X, t.B));
  const double diag = (t.X.colwise().maxCoeff() - t.X.colwise().minCoeff()).norm();

  FinalMapResult r = compute_final_map(t.X, t.X, t.T, t.B, g, g.height, FinalMapParams{});
  CHECK(r.arap_mean < 1e-10);
  CHECK(r.inverted_tets == 0);
  CHECK(r.max_planarity_error <= 1e-9 * diag);

  // Perturbed start still lands on exact planes.
  MatX3d X0 = t.X + 0.02 * MatX3d::Random(t.X.rows(), 3);
  FinalMapResult r2 = compute_final_map(t.X, X0, t.T, t.B, g, g.height, FinalMapParams{});
  CHECK(r2.max_planarity_error <= 1e-9 * diag);
  CHECK(r2.arap_mean < 1e-6);

  // Fixed heights: stretch the box along Z.
  std::vector<double> h = g.height;
  for (size_t p = 0; p < h.size(); ++p)
    if (g.patches.label[p] == make_label(2, 1)) h[p] = 1.0;
  FinalMapParams fp;
  fp.heights_fixed = true;
  FinalMapResult r3 = compute_final_map(t.X, t.X, t.T, t.B, g, h, fp);
  for (size_t p = 0; p < h.size(); ++p) CHECK(r3.heights[p] == doctest::Approx(h[p]).epsilon(1e-12));
  CHECK(r3.max_planarity_error <= 1e-9 * diag);
  CHECK(r3.inverted_tets == 0);
}

TEST_CASE("locator: node round trip and outside query") {
  SurfaceMesh s = testutil::make_box(3);
  TetMesh t = tetrahedralize(s, TetParams{});
  // A deformed copy acts as the "polycube" domain; the map goes back to t.X.
  Eigen::Matrix3d A;
  A << 1.1, 0.1, 0, 0, 0.9, 0.2, 0.05, 0, 1.2;
  MatX3d domain = t.X * A.transpose();
  TetLocator loc(domain, t.T);
  const double diag = (t.X.colwise().maxCoeff() - t.X.colwise().minCoeff()).norm();
  for (Index i = 0; i < t.X.rows(); i += 3) {
    auto p = loc.map(domain.row(i).transpose(), t.X);
    REQUIRE(p.has_value());
    CHECK((*p - t.X.row(i).transpose()).norm() <= 1e-9 * diag);
  }
  Eigen::Vector3d inside = domain.colwise().mean().transpose();
  CHECK(loc.locate(inside).has_value());
  CHECK_FALSE(loc.locate(Eigen::Vector3d(10, 10, 10)).has_value());
}
