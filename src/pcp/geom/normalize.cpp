#include "pcp/geom/normalize.h"

#include <cmath>

#include <Eigen/Geometry>

namespace pcp {

MatX3d Normalization::apply(const MatX3d& X) const {
  return ((X.rowwise() - center) * scale).eval();
}

MatX3d Normalization::inverse(const MatX3d& Y) const {
  return ((Y / scale).rowwise() + center).eval();
}

double total_area(const MatX3d& X, const MatX3i& F) {
  double a = 0.0;
  for (Index f = 0; f < F.rows(); ++f) {
    Eigen::Vector3d p0 = X.row(F(f, 0)), p1 = X.row(F(f, 1)), p2 = X.row(F(f, 2));
    a += 0.5 * (p1 - p0).cross(p2 - p0).norm();
  }
  return a;
}

Normalization normalization_for(const MatX3d& X, const MatX3i& boundary) {
  Normalization n;
  double area = 0.0;
  Eigen::RowVector3d c = Eigen::RowVector3d::Zero();
  for (Index f = 0; f < boundary.rows(); ++f) {
    Eigen::Vector3d p0 = X.row(boundary(f, 0)), p1 = X.row(boundary(f, 1)),
                    p2 = X.row(boundary(f, 2));
    double a = 0.5 * (p1 - p0).cross(p2 - p0).norm();
    area += a;
    c += a * (p0 + p1 + p2).transpose() / 3.0;
  }
  if (!(area > 0.0)) throw Error("normalization: boundary has zero area");
  n.center = c / area;
  n.scale = 1.0 / std::sqrt(area);
  return n;
}

}  // namespace pcp
