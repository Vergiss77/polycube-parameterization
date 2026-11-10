#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

namespace pcp {

using Index = Eigen::Index;
using MatX3d = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using MatX3i = Eigen::Matrix<int, Eigen::Dynamic, 3>;
using MatX4i = Eigen::Matrix<int, Eigen::Dynamic, 4>;

// Error raised by any pipeline stage; the message explains what went wrong.
class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Closed triangle surface: V is n x 3 positions, F is m x 3 vertex indices.
struct SurfaceMesh {
  MatX3d V;
  MatX3i F;
};

// Tetrahedral mesh with its boundary.
// X: node positions, T: tets (positive orientation), B: boundary triangles with
// outward orientation, B_to_input: index of the input surface face of each boundary triangle.
struct TetMesh {
  MatX3d X;
  MatX4i T;
  MatX3i B;
  Eigen::VectorXi B_to_input;
};

}  // namespace pcp
