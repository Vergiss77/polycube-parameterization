#include "pcp/geom/tet_gradient.h"

#include <Eigen/LU>
#include <Eigen/SVD>

namespace pcp {

TetGradients compute_tet_gradients(const MatX3d& rest, const MatX4i& T) {
  TetGradients g;
  g.G.resize(static_cast<size_t>(T.rows()));
  g.rest_volume.resize(T.rows());
  for (Index t = 0; t < T.rows(); ++t) {
    Eigen::Vector3d x0 = rest.row(T(t, 0));
    Eigen::Matrix3d D;
    for (int k = 0; k < 3; ++k) D.col(k) = rest.row(T(t, k + 1)).transpose() - x0;
    const double det = D.determinant();
    if (!(det > 0.0)) throw Error("tet " + std::to_string(t) + " has non-positive rest volume");
    Eigen::Matrix3d Dinv = D.inverse();
    Eigen::Matrix<double, 3, 4>& G = g.G[static_cast<size_t>(t)];
    // F = [x1-x0, x2-x0, x3-x0] * Dinv  =>  gradient of node k+1 is row k of Dinv.
    for (int k = 0; k < 3; ++k) G.col(k + 1) = Dinv.row(k).transpose();
    G.col(0) = -(G.col(1) + G.col(2) + G.col(3));
    g.rest_volume(t) = det / 6.0;
  }
  return g;
}

Eigen::Matrix3d deformation_gradient(const TetGradients& g, const MatX3d& X, const MatX4i& T,
                                     Index t) {
  const Eigen::Matrix<double, 3, 4>& G = g.G[static_cast<size_t>(t)];
  Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
  for (int i = 0; i < 4; ++i) F += X.row(T(t, i)).transpose() * G.col(i).transpose();
  return F;
}

Eigen::Matrix3d closest_rotation(const Eigen::Matrix3d& M) {
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d U = svd.matrixU(), V = svd.matrixV();
  Eigen::Matrix3d R = U * V.transpose();
  if (R.determinant() < 0.0) {
    U.col(2) = -U.col(2);
    R = U * V.transpose();
  }
  return R;
}

}  // namespace pcp
