#include "pcp/blocks/bspline.h"

#include <algorithm>
#include <string>

namespace pcp {

KnotVector open_uniform_knots(int degree, int n_elements) {
  if (degree < 1) throw Error("B-spline degree must be at least 1, got " + std::to_string(degree));
  if (n_elements < 1) throw Error("B-spline needs at least one element");
  KnotVector kv;
  kv.degree = degree;
  kv.knots.assign(static_cast<size_t>(degree), 0.0);
  for (int e = 0; e <= n_elements; ++e)
    kv.knots.push_back(e == n_elements ? 1.0 : static_cast<double>(e) / n_elements);
  kv.knots.insert(kv.knots.end(), static_cast<size_t>(degree), 1.0);
  return kv;
}

int find_span(const KnotVector& kv, double u) {
  const int p = kv.degree, n = kv.n_basis();
  const auto& U = kv.knots;
  if (u >= U[static_cast<size_t>(n)]) return n - 1;
  if (u <= U[static_cast<size_t>(p)]) return p;
  // Last s with U[s] <= u.
  const auto it = std::upper_bound(U.begin() + p, U.begin() + n + 1, u);
  return static_cast<int>(it - U.begin()) - 1;
}

void basis_functions(const KnotVector& kv, int s, double u, double* N, double* dN) {
  const int p = kv.degree;
  const auto& U = kv.knots;
  // ndu: upper triangle basis values, lower triangle knot differences (Piegl & Tiller A2.3).
  constexpr int kMax = 16;
  if (p + 1 > kMax) throw Error("B-spline degree too large");
  double ndu[kMax][kMax], left[kMax], right[kMax];
  ndu[0][0] = 1.0;
  for (int j = 1; j <= p; ++j) {
    left[j] = u - U[static_cast<size_t>(s + 1 - j)];
    right[j] = U[static_cast<size_t>(s + j)] - u;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];
      const double temp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu[j][j] = saved;
  }
  for (int j = 0; j <= p; ++j) N[j] = ndu[j][p];
  if (!dN) return;
  // First derivative: N'_{i,p} = p (N_{i,p-1} / (U[i+p]-U[i]) - N_{i+1,p-1} / (U[i+p+1]-U[i+1])).
  for (int r = 0; r <= p; ++r) {
    double d = 0.0;
    if (r >= 1) d += ndu[r - 1][p - 1] / ndu[p][r - 1];
    if (r <= p - 1) d -= ndu[r][p - 1] / ndu[p][r];
    dN[r] = p * d;
  }
}

std::vector<double> greville(const KnotVector& kv) {
  std::vector<double> g(static_cast<size_t>(kv.n_basis()));
  for (int i = 0; i < kv.n_basis(); ++i) {
    double s = 0.0;
    for (int k = 1; k <= kv.degree; ++k) s += kv.knots[static_cast<size_t>(i + k)];
    g[static_cast<size_t>(i)] = s / kv.degree;
  }
  return g;
}

namespace {

struct Local {
  int span[3];
  double N[3][16], dN[3][16];
};

Local local_basis(const TensorBSpline3& V, double u, double v, double w, bool derivs) {
  Local L;
  const double t[3] = {u, v, w};
  for (int d = 0; d < 3; ++d) {
    const KnotVector& kv = V.basis[static_cast<size_t>(d)];
    L.span[d] = find_span(kv, t[d]);
    basis_functions(kv, L.span[d], t[d], L.N[d], derivs ? L.dN[d] : nullptr);
  }
  return L;
}

}  // namespace

Eigen::Vector3d TensorBSpline3::eval(double u, double v, double w) const {
  const Local L = local_basis(*this, u, v, w, false);
  const int p0 = basis[0].degree, p1 = basis[1].degree, p2 = basis[2].degree;
  Eigen::Vector3d x = Eigen::Vector3d::Zero();
  for (int c = 0; c <= p2; ++c)
    for (int b = 0; b <= p1; ++b)
      for (int a = 0; a <= p0; ++a) {
        const double wgt = L.N[0][a] * L.N[1][b] * L.N[2][c];
        x += wgt * coefs[static_cast<size_t>(index(L.span[0] - p0 + a, L.span[1] - p1 + b, L.span[2] - p2 + c))];
      }
  return x;
}

void TensorBSpline3::eval_jacobian(double u, double v, double w, Eigen::Vector3d& x, Eigen::Matrix3d& J) const {
  const Local L = local_basis(*this, u, v, w, true);
  const int p0 = basis[0].degree, p1 = basis[1].degree, p2 = basis[2].degree;
  x.setZero();
  J.setZero();
  for (int c = 0; c <= p2; ++c)
    for (int b = 0; b <= p1; ++b)
      for (int a = 0; a <= p0; ++a) {
        const Eigen::Vector3d& P =
            coefs[static_cast<size_t>(index(L.span[0] - p0 + a, L.span[1] - p1 + b, L.span[2] - p2 + c))];
        x += L.N[0][a] * L.N[1][b] * L.N[2][c] * P;
        J.col(0) += L.dN[0][a] * L.N[1][b] * L.N[2][c] * P;
        J.col(1) += L.N[0][a] * L.dN[1][b] * L.N[2][c] * P;
        J.col(2) += L.N[0][a] * L.N[1][b] * L.dN[2][c] * P;
      }
}

}  // namespace pcp
