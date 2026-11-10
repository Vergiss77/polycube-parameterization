#include <doctest/doctest.h>

#include <cmath>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "helpers.h"
#include "pcp/geom/normalize.h"
#include "pcp/io/off.h"
#include "pcp/polycube/energies.h"
#include "pcp/polycube/feature_edges.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {

struct Fixture {
  MatX3d rest, X;
  DeformMesh mesh;
};

Fixture make_fixture(unsigned seed, double noise = 0.05) {
  SurfaceMesh s = testutil::make_box(2, 1.0, 1.3, 0.8);
  TetMesh t = tetrahedralize(s, TetParams{});
  Normalization n = normalization_for(t.X, t.B);
  Fixture f;
  f.rest = n.apply(t.X);
  f.mesh = DeformMesh::build(f.rest, t.T, t.B);
  std::mt19937 rng(seed);
  std::normal_distribution<double> nd(0.0, noise);
  f.X = f.rest;
  for (Index i = 0; i < f.X.rows(); ++i)
    for (int k = 0; k < 3; ++k) f.X(i, k) += nd(rng) * 0.3;
  return f;
}

// Central-difference check of add_gradient against value on random coordinates.
void check_gradient(const EnergyTerm& term, const MatX3d& X, double tol = 1e-5) {
  Eigen::VectorXd g = Eigen::VectorXd::Zero(3 * X.rows());
  term.add_gradient(X, g);
  std::mt19937 rng(3);
  std::uniform_int_distribution<Index> pick(0, 3 * X.rows() - 1);
  const double h = 1e-6;
  double max_err = 0.0, scale = g.cwiseAbs().maxCoeff() + 1e-12;
  for (int trial = 0; trial < 60; ++trial) {
    Index i = pick(rng);
    MatX3d Xp = X, Xm = X;
    Xp(i / 3, i % 3) += h;
    Xm(i / 3, i % 3) -= h;
    double fd = (term.value(Xp) - term.value(Xm)) / (2 * h);
    max_err = std::max(max_err, std::abs(fd - g(i)));
  }
  CHECK(max_err <= tol * scale);
}

double min_hessian_eigenvalue(const EnergyTerm& term, const DeformMesh& mesh, const MatX3d& X) {
  BlockSparse H = mesh.make_hessian();
  term.add_hessian(X, H);
  Eigen::MatrixXd D = Eigen::MatrixXd(H.matrix());
  CHECK((D - D.transpose()).cwiseAbs().maxCoeff() <= 1e-12 * (1.0 + D.cwiseAbs().maxCoeff()));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(D);
  return es.eigenvalues().minCoeff() / (1.0 + es.eigenvalues().cwiseAbs().maxCoeff());
}

}  // namespace

TEST_CASE("polycube error on provided models" * doctest::test_suite("models")) {
  SurfaceMesh ispar = read_off(testutil::models_dir() / "I-SPAR-1B59CD84.off");
  SurfaceMesh tail = read_off(testutil::models_dir() / "TAILSITTER.off");
  CHECK(polycube_error(ispar.V, ispar.F) < 1e-6);
  double d = polycube_error(tail.V, tail.F);
  MESSAGE("TAILSITTER initial polycube error = " << d);
  CHECK(d > 0.05);
}

TEST_CASE("E1: gradient and PSD Hessian") {
  Fixture f = make_fixture(11);
  L1NormalEnergy e(f.mesh);
  e.weight = 0.7;
  for (double eps : {1.0, 0.1, 0.01}) {
    e.eps = eps;
    e.begin_iteration(f.X);
    check_gradient(e, f.X);
    CHECK(min_hessian_eigenvalue(e, f.mesh, f.X) >= -1e-10);
  }
}

TEST_CASE("E1 tends to the l1 norm as eps goes to zero") {
  Fixture f = make_fixture(5);
  L1NormalEnergy e(f.mesh);
  e.eps = 1e-12;
  e.begin_iteration(f.X);
  double area = total_area(f.X, f.mesh.B);
  double expected = (polycube_error(f.X, f.mesh.B) + 1.0) * area / f.mesh.rest_area;
  CHECK(e.value(f.X) == doctest::Approx(expected).epsilon(1e-5));
}

TEST_CASE("ARAP: zero for rigid motions, gradient check, PSD Hessian") {
  Fixture f = make_fixture(2);
  ArapEnergy e(f.mesh);
  Eigen::Matrix3d R = Eigen::AngleAxisd(0.7, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
  MatX3d Y = (f.rest * R.transpose()).rowwise() + Eigen::RowVector3d(0.3, -1, 2);
  CHECK(e.value(Y) < 1e-20);
  check_gradient(e, f.X);
  CHECK(min_hessian_eigenvalue(e, f.mesh, f.X) >= -1e-10);
}

TEST_CASE("E_eta: only creases contribute on a box, gradient check") {
  Fixture f = make_fixture(4);
  ComplexityEnergy e(f.mesh);
  e.begin_iteration(f.rest);
  // On the undeformed box, n_i - n_j is zero for coplanar pairs and has squared norm 2 across
  // the 90-degree creases.
  double expected = 0.0;
  for (size_t i = 0; i < f.mesh.topo.pairs.size(); ++i) {
    Eigen::Vector3d n0, n1;
    Eigen::Matrix<double, 3, 12> J0, J1;
    pair_normals(f.rest, f.mesh.topo.pairs[i], n0, n1, J0, J1);
    if (std::abs(n0.dot(n1)) < 1e-9) expected += 2.0 * e.edge_areas()(static_cast<Index>(i));
  }
  CHECK(e.value(f.rest) == doctest::Approx(expected / f.mesh.rest_area));
  e.begin_iteration(f.X);
  e.weight = 2.0;
  check_gradient(e, f.X);
  CHECK(min_hessian_eigenvalue(e, f.mesh, f.X) >= -1e-10);
}

TEST_CASE("E_e: gradient check, zero weight changes nothing, file loading") {
  Fixture f = make_fixture(8);
  std::vector<FeatureEdgeEnergy::Edge> edges;
  for (int i = 0; i < 6; ++i) edges.push_back({i * 3, i % 2 == 0 ? 0.0 : 1.0});
  FeatureEdgeEnergy e(f.mesh, edges);
  e.weight = 3.0;
  check_gradient(e, f.X);
  CHECK(min_hessian_eigenvalue(e, f.mesh, f.X) >= -1e-10);
  e.weight = 0.0;
  CHECK(e.value(f.X) == 0.0);
  FeatureEdgeEnergy none(f.mesh, {});
  CHECK(none.empty());
  CHECK(none.value(f.X) == 0.0);

  const FacePair& p = f.mesh.topo.pairs[0];
  auto path = testutil::temp_dir() / "edges.txt";
  testutil::write_text(path, "# test\n" + std::to_string(p.v0) + " " + std::to_string(p.v1) +
                                 " corner\n");
  auto loaded = read_feature_edges(path, f.mesh.topo);
  REQUIRE(loaded.size() == 1);
  CHECK(loaded[0].pair == 0);
  CHECK(loaded[0].target == 0.0);
  testutil::write_text(path, std::to_string(p.v0) + " " + std::to_string(p.v1) + " bent\n");
  CHECK_THROWS_AS(read_feature_edges(path, f.mesh.topo), Error);
}

TEST_CASE("Barrier: gradient check and blow-up near foldover") {
  Fixture f = make_fixture(9);
  FoldoverBarrier b(f.mesh);
  b.weight = 1e-3;
  check_gradient(b, f.X);
  CHECK(min_hessian_eigenvalue(b, f.mesh, f.X) >= -1e-10);

  // Fold one face pair onto itself: move the opposite vertex of f1 through the shared edge
  // until the two normals become opposite.
  const FacePair& p = f.mesh.topo.pairs[0];
  MatX3d Y = f.rest;
  Eigen::RowVector3d mid = 0.5 * (Y.row(p.v0) + Y.row(p.v1));
  Eigen::RowVector3d o0 = Y.row(p.nodes[2]);
  double prev = b.value(Y);
  for (double t : {0.9, 0.99, 0.999}) {
    // Place o1 on the same side as o0 (folding f1 onto f0).
    Y.row(p.nodes[3]) = mid + t * (o0 - mid) + (1 - t) * (f.rest.row(p.nodes[3]) - mid);
    double v = b.value(Y);
    CHECK(v > prev);
    prev = v;
  }
  Y.row(p.nodes[3]) = o0 + 0.0 * mid;
  CHECK(std::isinf(b.value(Y)));
}

TEST_CASE("area constraint: zero at rest, gradient check") {
  Fixture f = make_fixture(12);
  CHECK(std::abs(area_constraint(f.rest, f.mesh.B, f.mesh.rest_area)) < 1e-12);
  Eigen::VectorXd g;
  area_constraint_gradient(f.X, f.mesh.B, g);
  const double h = 1e-6;
  for (Index i = 0; i < g.size(); i += 7) {
    MatX3d Xp = f.X, Xm = f.X;
    Xp(i / 3, i % 3) += h;
    Xm(i / 3, i % 3) -= h;
    double fd = (area_constraint(Xp, f.mesh.B, 0) - area_constraint(Xm, f.mesh.B, 0)) / (2 * h);
    CHECK(std::abs(fd - g(i)) <= 1e-6 * (1 + std::abs(g(i))));
  }
}

TEST_CASE("Tet inversion barrier: gradient check, PSD Hessian, infinite when inverted") {
  Fixture f = make_fixture(13, 0.02);
  TetInversionBarrier b(f.mesh);
  b.weight = 1e-2;
  REQUIRE(b.inverted(f.X) == 0);
  check_gradient(b, f.X);
  CHECK(min_hessian_eigenvalue(b, f.mesh, f.X) >= -1e-10);
  CHECK(std::abs(b.value(f.rest)) < 1e-14);
  MatX3d Y = f.rest;
  // Push node 0 of tet 0 through its opposite face.
  Eigen::RowVector3d opp = (Y.row(f.mesh.T(0, 1)) + Y.row(f.mesh.T(0, 2)) + Y.row(f.mesh.T(0, 3))) / 3.0;
  Y.row(f.mesh.T(0, 0)) = opp + 1.5 * (opp - Y.row(f.mesh.T(0, 0)));
  CHECK(std::isinf(b.value(Y)));
}
