#pragma once

#include <vector>

#include <Eigen/Core>

#include "pcp/core/types.h"

namespace pcp {

// Per-tet linear gradient operator computed on a rest shape. For tet t with node positions
// collected as columns of the 3x4 matrix P, the deformation gradient is F = P * G[t]^T,
// where G[t] is 3x4 (column i = gradient of the barycentric coordinate of node i).
struct TetGradients {
  std::vector<Eigen::Matrix<double, 3, 4>> G;
  Eigen::VectorXd rest_volume;
};

TetGradients compute_tet_gradients(const MatX3d& rest, const MatX4i& T);

// Deformation gradient of tet t at positions X.
Eigen::Matrix3d deformation_gradient(const TetGradients& g, const MatX3d& X, const MatX4i& T,
                                     Index t);

// Closest rotation to M (polar decomposition via SVD, reflection-corrected, det = +1).
Eigen::Matrix3d closest_rotation(const Eigen::Matrix3d& M);

}  // namespace pcp
