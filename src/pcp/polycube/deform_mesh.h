#pragma once

#include <vector>

#include "pcp/core/block_sparse.h"
#include "pcp/core/types.h"
#include "pcp/geom/boundary_topology.h"
#include "pcp/geom/tet_gradient.h"

namespace pcp {

// Static data shared by all energy terms: connectivity, rest-shape operators and the Hessian
// sparsity pattern (tets + pairs of adjacent boundary triangles).
struct DeformMesh {
  MatX4i T;
  MatX3i B;
  BoundaryTopology topo;
  TetGradients grads;
  double rest_volume = 0.0;  // total rest volume
  double rest_area = 0.0;    // total rest boundary area

  // Build from the (normalized) rest positions.
  static DeformMesh build(const MatX3d& rest, const MatX4i& T, const MatX3i& B);
  BlockSparse make_hessian() const;
  Index n_nodes() const { return n_nodes_; }

 private:
  Index n_nodes_ = 0;
};

// Flattening helpers: vectors of length 3N with entry 3*i + k for coordinate k of node i.
Eigen::VectorXd flatten(const MatX3d& X);
MatX3d unflatten(const Eigen::VectorXd& x);

// Positions of the three vertices of boundary face f.
inline void face_vertices(const MatX3d& X, const MatX3i& B, Index f, Eigen::Vector3d& x0,
                          Eigen::Vector3d& x1, Eigen::Vector3d& x2) {
  x0 = X.row(B(f, 0));
  x1 = X.row(B(f, 1));
  x2 = X.row(B(f, 2));
}

}  // namespace pcp
