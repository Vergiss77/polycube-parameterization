#include <doctest/doctest.h>

#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "helpers.h"
#include "pcp/core/constants.h"
#include "pcp/geom/normalize.h"
#include "pcp/io/validate.h"
#include "pcp/polycube/kkt.h"
#include "pcp/polycube/l1_solver.h"
#include "pcp/polycube/rotation.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {

struct Problem {
  MatX3d rest;
  DeformMesh mesh;
};

Problem make_problem(const SurfaceMesh& s, const TetParams& tp = {}) {
  TetMesh t = tetrahedralize(s, tp);
  Normalization n = normalization_for(t.X, t.B);
  Problem p;
  p.rest = n.apply(t.X);
  p.mesh = DeformMesh::build(p.rest, t.T, t.B);
  return p;
}

}  // namespace

TEST_CASE("KKT step solves an equality-constrained quadratic in one step") {
  const int n = 12;
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd Hd = M * M.transpose() + n * Eigen::MatrixXd::Identity(n, n);
  Eigen::SparseMatrix<double> H = Hd.sparseView();
  Eigen::VectorXd g = Eigen::VectorXd::Random(n), a = Eigen::VectorXd::Random(n);
  const double c = 0.7;
  for (auto kind : {LinearSolverOptions::Kind::Direct, LinearSolverOptions::Kind::Iterative}) {
    LinearSolverOptions o;
    o.kind = kind;
    o.cg_tolerance = 1e-14;
    KktSolver kkt(o);
    Eigen::VectorXd dx;
    double lambda = 0.0;
    // At x = 0 the constraint value is a^T 0 - c.
    REQUIRE(kkt.solve(H, 0.0, g, a, -c, dx, lambda));
    CHECK(std::abs(a.dot(dx) - c) < 1e-10);
    CHECK((Hd * dx + g + lambda * a).norm() < 1e-8);
  }
}

TEST_CASE("line search: accepts a descent step, fails when every step is infeasible") {
  auto quad = [](double t) { return (t - 0.3) * (t - 0.3); };
  double t = backtracking_line_search(quad, quad(0.0), -0.6, 30);
  CHECK(t > 0.0);
  CHECK(quad(t) < quad(0.0));
  auto infeasible = [](double t) { return t > 0 ? std::numeric_limits<double>::infinity() : 0.0; };
  CHECK(backtracking_line_search(infeasible, 0.0, -1.0, 30) == 0.0);
}

TEST_CASE("slim triangle count") {
  MatX3d X(4, 3);
  X << 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0.01, 0;
  MatX3i B(2, 3);
  B << 0, 1, 2, 0, 1, 3;
  CHECK(count_slim_triangles(X, B, 5.0) == 1);
}

TEST_CASE("solver records and logs a stall when no step can be taken") {
  Problem p = make_problem(testutil::make_icosphere(2));
  L1Params params;
  params.global_rotation = false;
  params.max_line_search = -1;  // no trial step at all: every line search fails
  std::vector<std::string> lines;
  L1Result r = solve_l1_polycube(p.rest, p.mesh, params, {},
                                 [&](const std::string& s) { lines.push_back(s); });
  REQUIRE_FALSE(r.stalls.empty());
  bool logged = false;
  for (const auto& l : lines) logged |= l.find("barrier stall") != std::string::npos;
  CHECK(logged);
}

TEST_CASE("global orientation search recovers a 30 degree rotation") {
  SurfaceMesh s = testutil::make_box(4, 1.0, 1.7, 0.6);
  Eigen::Matrix3d Rz = Eigen::AngleAxisd(30.0 * kDegToRad, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  MatX3d V = s.V * Rz.transpose();
  Eigen::Matrix3d R = global_orientation_search(boundary_area_vectors(V, s.F));
  // R * Rz must be a cube symmetry (aligns the box with the axes again).
  double residual = rotation_angle(canonicalize_cube_rotation(R * Rz)) * kRadToDeg;
  CHECK(residual < 0.5);
  // An axis-aligned input keeps the identity.
  Eigen::Matrix3d I = global_orientation_search(boundary_area_vectors(s.V, s.F));
  CHECK(rotation_angle(I) < 1e-9);
}

TEST_CASE("l1 solver: sphere converges below threshold with per-stage logs") {
  // A sphere must be squashed into a cube, so the tets at the cube corners need room to deform:
  // with the default (coarse) tet size the tet inversion barrier blocks convergence at d ~ 0.09.
  SurfaceMesh s = testutil::make_icosphere(2);
  REQUIRE(validate_surface(s).ok);
  TetParams tp;
  tp.size_factor = 1.0;
  Problem p = make_problem(s, tp);
  std::vector<std::string> lines;
  L1Result r = solve_l1_polycube(p.rest, p.mesh, L1Params{}, {},
                                 [&](const std::string& l) { lines.push_back(l); });
  if (!r.reached_threshold)
    for (const auto& l : lines) MESSAGE(l);
  CHECK(r.reached_threshold);
  CHECK(r.error < 1e-3);
  CHECK(r.area_rel_diff < 1e-6);
  CHECK(r.stalls.empty());
  int stage_lines = 0;
  for (const auto& l : lines)
    if (l.rfind("stage", 0) == 0 && l.find("alpha") != std::string::npos &&
        l.find("eps") != std::string::npos && l.find("error") != std::string::npos &&
        l.find("arap") != std::string::npos)
      ++stage_lines;
  CHECK(stage_lines == r.stages);
}

TEST_CASE("l1 solver: identical results on repeated runs") {
  Problem p = make_problem(testutil::make_icosphere(2));
  L1Params params;
  params.max_stages = 4;
  L1Result a = solve_l1_polycube(p.rest, p.mesh, params);
  L1Result b = solve_l1_polycube(p.rest, p.mesh, params);
  CHECK(a.X == b.X);
}
