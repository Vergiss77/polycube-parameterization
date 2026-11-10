// B-spline basics: knot vectors, basis functions and derivatives, tensor product volumes.
#include <doctest/doctest.h>

#include <cmath>

#include "pcp/blocks/bspline.h"

using namespace pcp;

TEST_CASE("bspline: open uniform knots, partition of unity, derivatives") {
  for (int p = 1; p <= 4; ++p)
    for (int ne : {1, 3}) {
      KnotVector kv = open_uniform_knots(p, ne);
      CHECK(kv.n_basis() == ne + p);
      CHECK(kv.n_elements() == ne);
      CHECK(kv.knots.front() == 0.0);
      CHECK(kv.knots.back() == 1.0);
      for (int i = 0; i <= 50; ++i) {
        const double u = i / 50.0;
        const int s = find_span(kv, u);
        CHECK(s >= p);
        CHECK(s < kv.n_basis());
        double N[16], dN[16];
        basis_functions(kv, s, u, N, dN);
        double sum = 0.0, dsum = 0.0;
        for (int r = 0; r <= p; ++r) {
          sum += N[r];
          dsum += dN[r];
          CHECK(N[r] >= -1e-15);
        }
        CHECK(std::abs(sum - 1.0) <= 1e-14);
        CHECK(std::abs(dsum) <= 1e-12);
        // Finite differences of every basis function (same span on both sides away from knots).
        const double h = 1e-6;
        if (u - h <= 0.0 || u + h >= 1.0 || find_span(kv, u - h) != s || find_span(kv, u + h) != s) continue;
        double Np[16], Nm[16];
        basis_functions(kv, s, u + h, Np, nullptr);
        basis_functions(kv, s, u - h, Nm, nullptr);
        for (int r = 0; r <= p; ++r) CHECK(std::abs((Np[r] - Nm[r]) / (2 * h) - dN[r]) <= 1e-7);
      }
    }
  CHECK_THROWS_AS(open_uniform_knots(0, 2), Error);
  CHECK_THROWS_AS(open_uniform_knots(2, 0), Error);
}

TEST_CASE("bspline volume: linear precision with Greville points, corner interpolation") {
  TensorBSpline3 V;
  V.basis = {open_uniform_knots(3, 2), open_uniform_knots(2, 3), open_uniform_knots(1, 1)};
  const auto gu = greville(V.basis[0]), gv = greville(V.basis[1]), gw = greville(V.basis[2]);
  Eigen::Matrix3d A;
  A << 2, 0.5, 0, -0.3, 1, 0.2, 0.1, 0, 3;
  const Eigen::Vector3d t(1, -2, 0.5);
  V.coefs.resize(size_t(V.n(0) * V.n(1) * V.n(2)));
  for (int k = 0; k < V.n(2); ++k)
    for (int j = 0; j < V.n(1); ++j)
      for (int i = 0; i < V.n(0); ++i)
        V.coefs[size_t(V.index(i, j, k))] = A * Eigen::Vector3d(gu[size_t(i)], gv[size_t(j)], gw[size_t(k)]) + t;
  for (double u : {0.0, 0.17, 0.5, 0.93, 1.0})
    for (double v : {0.0, 0.4, 1.0})
      for (double w : {0.0, 0.66, 1.0}) {
        Eigen::Vector3d x;
        Eigen::Matrix3d J;
        V.eval_jacobian(u, v, w, x, J);
        CHECK((x - (A * Eigen::Vector3d(u, v, w) + t)).norm() <= 1e-13);
        CHECK((V.eval(u, v, w) - x).norm() <= 1e-14);
        CHECK((J - A).norm() <= 1e-12);
      }
  // Corners interpolate the corner control points.
  CHECK((V.eval(0, 0, 0) - V.coefs.front()).norm() <= 1e-15);
  CHECK((V.eval(1, 1, 1) - V.coefs.back()).norm() <= 1e-14);
}
