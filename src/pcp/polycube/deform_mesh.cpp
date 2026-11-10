#include "pcp/polycube/deform_mesh.h"

#include "pcp/geom/normalize.h"

namespace pcp {

DeformMesh DeformMesh::build(const MatX3d& rest, const MatX4i& T, const MatX3i& B) {
  DeformMesh m;
  m.T = T;
  m.B = B;
  m.n_nodes_ = rest.rows();
  m.topo = build_boundary_topology(B);
  m.grads = compute_tet_gradients(rest, T);
  m.rest_volume = m.grads.rest_volume.sum();
  m.rest_area = total_area(rest, B);
  return m;
}

BlockSparse DeformMesh::make_hessian() const {
  std::vector<std::vector<int>> elems;
  elems.reserve(static_cast<size_t>(T.rows()) + topo.pairs.size());
  for (Index t = 0; t < T.rows(); ++t) elems.push_back({T(t, 0), T(t, 1), T(t, 2), T(t, 3)});
  for (const auto& p : topo.pairs)
    elems.push_back({p.nodes[0], p.nodes[1], p.nodes[2], p.nodes[3]});
  return BlockSparse(n_nodes_, elems);
}

Eigen::VectorXd flatten(const MatX3d& X) {
  Eigen::VectorXd x(3 * X.rows());
  for (Index i = 0; i < X.rows(); ++i)
    for (int k = 0; k < 3; ++k) x(3 * i + k) = X(i, k);
  return x;
}

MatX3d unflatten(const Eigen::VectorXd& x) {
  MatX3d X(x.size() / 3, 3);
  for (Index i = 0; i < X.rows(); ++i)
    for (int k = 0; k < 3; ++k) X(i, k) = x(3 * i + k);
  return X;
}

}  // namespace pcp
