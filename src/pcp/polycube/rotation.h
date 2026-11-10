#pragma once

#include <Eigen/Core>

#include "pcp/core/types.h"

namespace pcp {

// Area-weighted l1 alignment cost of boundary normals under rotation R:
//   cost(R) = sum_f |A_f| * ||R n_f||_1 / sum_f |A_f| - 1   (exact, eps = 0), or the smoothed
//   version sum_f sum_c sqrt((R a_f)_c^2 + eps |a_f|^2) / sum_f |a_f| - 1 when eps > 0.
double rotation_cost(const MatX3d& area_vectors, const Eigen::Matrix3d& R, double eps);

// Area vectors of all boundary faces.
MatX3d boundary_area_vectors(const MatX3d& X, const MatX3i& B);

// Coarse search over the fundamental zone of the cube symmetry group (ZYX Euler angles in
// [0, 90) degrees on a regular grid), followed by local refinement. Ties are broken towards
// the identity. The result is canonicalized: among the 24 rotations equivalent under cube
// symmetries, the one with the smallest rotation angle is returned.
Eigen::Matrix3d global_orientation_search(const MatX3d& area_vectors, double grid_step_deg = 6.0);

// Local Newton refinement of R (a few iterations) on the smoothed cost.
Eigen::Matrix3d refine_rotation(const MatX3d& area_vectors, const Eigen::Matrix3d& R0, double eps,
                                int max_iters = 10);

// Among the 24 rotations S * R (S a cube symmetry), return the one closest to the identity.
Eigen::Matrix3d canonicalize_cube_rotation(const Eigen::Matrix3d& R);

// Rotation angle (radians) of R.
double rotation_angle(const Eigen::Matrix3d& R);

}  // namespace pcp
