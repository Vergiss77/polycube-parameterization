#include <doctest/doctest.h>

#include <random>

#include <Eigen/Geometry>
#include <Eigen/LU>

#include "helpers.h"
#include "pcp/core/block_sparse.h"
#include "pcp/geom/boundary_topology.h"
#include "pcp/geom/normalize.h"
#include "pcp/geom/tet_gradient.h"
#include "pcp/geom/triangle.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {
Eigen::Matrix<double, 9, 1> random9(std::mt19937& rng) {
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  Eigen::Matrix<double, 9, 1> x;
  for (int i = 0; i < 9; ++i) x(i) = u(rng);
  return x;
}
}  // namespace

TEST_CASE("normalization: boundary area becomes one and inverse restores") {
  SurfaceMesh s = testutil::make_box(3, 2.0, 5.0, 1.5);
  s.V.rowwise() += Eigen::RowVector3d(10, -3, 7);
  Normalization n = normalization_for(s.V, s.F);
  MatX3d Y = n.apply(s.V);
  CHECK(total_area(Y, s.F) == doctest::Approx(1.0).epsilon(1e-12));
  CHECK((n.inverse(Y) - s.V).cwiseAbs().maxCoeff() <= 1e-12 * s.V.cwiseAbs().maxCoeff());
}

TEST_CASE("tet gradients recover affine deformations") {
  SurfaceMesh s = testutil::make_box(2);
  TetMesh t = tetrahedralize(s, TetParams{});
  TetGradients g = compute_tet_gradients(t.X, t.T);
  Eigen::Matrix3d A;
  A << 1.2, 0.3, -0.1, 0.0, 0.9, 0.4, 0.2, -0.5, 1.1;
  MatX3d Y = (t.X * A.transpose()).rowwise() + Eigen::RowVector3d(1, 2, 3);
  double vol = 0.0;
  for (Index i = 0; i < t.T.rows(); ++i) {
    CHECK((deformation_gradient(g, Y, t.T, i) - A).norm() < 1e-9);
    vol += g.rest_volume(i);
  }
  CHECK(vol == doctest::Approx(1.0));
}

TEST_CASE("closest_rotation is orthogonal with det +1") {
  std::mt19937 rng(1);
  for (int k = 0; k < 20; ++k) {
    Eigen::Matrix3d M = Eigen::Map<Eigen::Matrix3d>(random9(rng).data());
    Eigen::Matrix3d R = closest_rotation(M);
    CHECK((R.transpose() * R - Eigen::Matrix3d::Identity()).norm() < 1e-12);
    CHECK(R.determinant() == doctest::Approx(1.0));
  }
  Eigen::Matrix3d Rz = Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  CHECK((closest_rotation(2.0 * Rz) - Rz).norm() < 1e-12);
}

TEST_CASE("triangle derivatives match finite differences") {
  std::mt19937 rng(7);
  const double h = 1e-6;
  for (int trial = 0; trial < 10; ++trial) {
    Eigen::Matrix<double, 9, 1> x = random9(rng);
    auto a_of = [](const Eigen::Matrix<double, 9, 1>& y) {
      return tri::area_vector(y.segment<3>(0), y.segment<3>(3), y.segment<3>(6));
    };
    auto n_of = [&](const Eigen::Matrix<double, 9, 1>& y) { return a_of(y).normalized().eval(); };
    tri::Mat39 Ja = tri::area_vector_jacobian(x.segment<3>(0), x.segment<3>(3), x.segment<3>(6));
    tri::Mat39 Jn = tri::unit_normal_jacobian(x.segment<3>(0), x.segment<3>(3), x.segment<3>(6));
    tri::Mat39 Ja_fd, Jn_fd;
    for (int i = 0; i < 9; ++i) {
      Eigen::Matrix<double, 9, 1> xp = x, xm = x;
      xp(i) += h;
      xm(i) -= h;
      Ja_fd.col(i) = (a_of(xp) - a_of(xm)) / (2 * h);
      Jn_fd.col(i) = (n_of(xp) - n_of(xm)) / (2 * h);
    }
    CHECK((Ja - Ja_fd).norm() <= 1e-6 * (1.0 + Ja.norm()));
    CHECK((Jn - Jn_fd).norm() <= 1e-6 * (1.0 + Jn.norm()));
    // Hessian of each component: derivative of the Jacobian rows.
    for (int c = 0; c < 3; ++c) {
      Eigen::Matrix<double, 9, 9> H = tri::area_vector_hessian(c), H_fd;
      for (int i = 0; i < 9; ++i) {
        Eigen::Matrix<double, 9, 1> xp = x, xm = x;
        xp(i) += h;
        xm(i) -= h;
        tri::Mat39 Jp = tri::area_vector_jacobian(xp.segment<3>(0), xp.segment<3>(3), xp.segment<3>(6));
        tri::Mat39 Jm = tri::area_vector_jacobian(xm.segment<3>(0), xm.segment<3>(3), xm.segment<3>(6));
        H_fd.col(i) = ((Jp.row(c) - Jm.row(c)) / (2 * h)).transpose();
      }
      CHECK((H - H_fd).norm() <= 1e-6);
    }
  }
}

TEST_CASE("BlockSparse pattern matches brute force") {
  SurfaceMesh s = testutil::make_box(2);
  TetMesh t = tetrahedralize(s, TetParams{});
  BoundaryTopology topo = build_boundary_topology(t.B);
  std::vector<std::vector<int>> elems;
  for (Index i = 0; i < t.T.rows(); ++i) elems.push_back({t.T(i, 0), t.T(i, 1), t.T(i, 2), t.T(i, 3)});
  for (const auto& p : topo.pairs) elems.push_back({p.nodes[0], p.nodes[1], p.nodes[2], p.nodes[3]});
  BlockSparse bs(t.X.rows(), elems);

  std::vector<Eigen::Triplet<double>> trip;
  for (const auto& e : elems)
    for (int a : e)
      for (int b : e)
        for (int p = 0; p < 3; ++p)
          for (int q = 0; q < 3; ++q) trip.emplace_back(3 * a + p, 3 * b + q, 1.0);
  for (Index i = 0; i < t.X.rows(); ++i)
    for (int p = 0; p < 3; ++p)
      for (int q = 0; q < 3; ++q) trip.emplace_back(3 * i + p, 3 * i + q, 1.0);
  Eigen::SparseMatrix<double> ref(3 * t.X.rows(), 3 * t.X.rows());
  ref.setFromTriplets(trip.begin(), trip.end());
  REQUIRE(ref.nonZeros() == bs.matrix().nonZeros());
  for (Index c = 0; c <= ref.outerSize(); ++c)
    REQUIRE(ref.outerIndexPtr()[c] == bs.matrix().outerIndexPtr()[c]);
  for (Index k = 0; k < ref.nonZeros(); ++k)
    REQUIRE(ref.innerIndexPtr()[k] == bs.matrix().innerIndexPtr()[k]);

  // Accumulation lands in the right entries.
  bs.add_block(t.T(0, 0), t.T(0, 1), Eigen::Matrix3d::Constant(2.0));
  CHECK(bs.matrix().coeff(3 * t.T(0, 0) + 1, 3 * t.T(0, 1) + 2) == 2.0);
  int far = -1;
  for (int j = 1; j < static_cast<int>(t.X.rows()) && far < 0; ++j)
    if (ref.coeff(0, 3 * j) == 0.0) far = j;
  REQUIRE(far > 0);
  CHECK_THROWS_AS(bs.add_block(0, far, Eigen::Matrix3d::Zero()), Error);
}

TEST_CASE("boundary topology pairs every edge of a closed surface") {
  SurfaceMesh s = testutil::make_box(3);
  BoundaryTopology topo = build_boundary_topology(s.F);
  CHECK(topo.pairs.size() == static_cast<size_t>(s.F.rows() * 3 / 2));
  for (const auto& nb : topo.face_nbrs)
    for (int f : nb) CHECK(f >= 0);
}
