#pragma once

// B-spline basics (design D4): open uniform knot vectors on [0, 1], basis functions with first
// derivatives (Cox-de Boor, Piegl & Tiller A2.2 / A2.3), and trivariate tensor product volumes.

#include <array>
#include <vector>

#include "pcp/core/types.h"

namespace pcp {

struct KnotVector {
  int degree = 3;
  std::vector<double> knots;  // non-decreasing, first and last value repeated degree + 1 times
  int n_basis() const { return static_cast<int>(knots.size()) - degree - 1; }
  int n_elements() const { return n_basis() - degree; }  // for open uniform vectors
};

// Open uniform knot vector on [0, 1] with n_elements equal knot spans.
KnotVector open_uniform_knots(int degree, int n_elements);

// Knot span index s (degree <= s < n_basis) with knots[s] <= u < knots[s + 1]; u = 1 belongs to
// the last non-empty span.
int find_span(const KnotVector& kv, double u);

// Values N[0..p] and first derivatives dN[0..p] (dN may be null) of the p + 1 basis functions
// that are non-zero on span s, i.e. functions s - p .. s.
void basis_functions(const KnotVector& kv, int span, double u, double* N, double* dN);

// Greville abscissae (knot averages) of the basis functions.
std::vector<double> greville(const KnotVector& kv);

// Trivariate tensor product B-spline volume in R^3. Control point (i, j, k) is stored at
// i + n_u * (j + n_v * k) (u fastest, w slowest).
struct TensorBSpline3 {
  std::array<KnotVector, 3> basis;
  std::vector<Eigen::Vector3d> coefs;

  int n(int dir) const { return basis[static_cast<size_t>(dir)].n_basis(); }
  int index(int i, int j, int k) const { return i + n(0) * (j + n(1) * k); }
  Eigen::Vector3d eval(double u, double v, double w) const;
  // Position and Jacobian J(:, d) = d x / d (u, v, w)_d.
  void eval_jacobian(double u, double v, double w, Eigen::Vector3d& x, Eigen::Matrix3d& J) const;
};

}  // namespace pcp
