#pragma once

#include <Eigen/Core>

namespace pcp::tri {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Mat39 = Eigen::Matrix<double, 3, 9>;

// Cross-product matrix: cross(v) * w == v.cross(w).
inline Mat3 cross(const Vec3& v) {
  Mat3 m;
  m << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
  return m;
}

// Area vector a = 0.5 (x1 - x0) x (x2 - x0) = area * unit normal.
inline Vec3 area_vector(const Vec3& x0, const Vec3& x1, const Vec3& x2) {
  return 0.5 * (x1 - x0).cross(x2 - x0);
}

// Jacobian of the area vector w.r.t. the 9 stacked vertex coordinates [x0; x1; x2].
// Block i equals 0.5 * cross(x_{i-1} - x_{i+1}).
inline Mat39 area_vector_jacobian(const Vec3& x0, const Vec3& x1, const Vec3& x2) {
  Mat39 J;
  J.block<3, 3>(0, 0) = 0.5 * cross(x2 - x1);
  J.block<3, 3>(0, 3) = 0.5 * cross(x0 - x2);
  J.block<3, 3>(0, 6) = 0.5 * cross(x1 - x0);
  return J;
}

// Constant Hessian (9x9) of component c of the area vector.
inline Eigen::Matrix<double, 9, 9> area_vector_hessian(int c) {
  // a_c = 0.5 * sum_i eps(c,p,q) x_i[p] x_{i+1}[q] over cyclic i.
  Eigen::Matrix<double, 9, 9> H = Eigen::Matrix<double, 9, 9>::Zero();
  auto eps = [](int i, int j, int k) {
    return static_cast<double>((i - j) * (j - k) * (k - i)) / 2.0;
  };
  for (int i = 0; i < 3; ++i) {
    int j = (i + 1) % 3;
    for (int p = 0; p < 3; ++p)
      for (int q = 0; q < 3; ++q) {
        double e = 0.5 * eps(c, p, q);
        H(3 * i + p, 3 * j + q) += e;
        H(3 * j + q, 3 * i + p) += e;
      }
  }
  return H;
}

// Jacobian of the unit normal n = a / |a| w.r.t. the 9 vertex coordinates.
inline Mat39 unit_normal_jacobian(const Vec3& x0, const Vec3& x1, const Vec3& x2) {
  Vec3 a = area_vector(x0, x1, x2);
  double len = a.norm();
  Vec3 n = a / len;
  Mat3 P = (Mat3::Identity() - n * n.transpose()) / len;
  return P * area_vector_jacobian(x0, x1, x2);
}

}  // namespace pcp::tri
