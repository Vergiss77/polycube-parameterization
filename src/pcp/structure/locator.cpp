#include "pcp/structure/locator.h"

#include <Eigen/LU>
#include <igl/AABB.h>

namespace pcp {

struct TetLocator::Impl {
  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  igl::AABB<Eigen::MatrixXd, 3> tree;
};

TetLocator::TetLocator(const MatX3d& domain, const MatX4i& T) : impl_(std::make_unique<Impl>()) {
  impl_->V = domain;
  impl_->T = T;
  impl_->tree.init(impl_->V, impl_->T);
}

TetLocator::~TetLocator() = default;

std::optional<TetLocator::Location> TetLocator::locate(const Eigen::Vector3d& q) const {
  Eigen::RowVector3d qr = q.transpose();
  std::vector<int> hits = impl_->tree.find(impl_->V, impl_->T, qr, true);
  if (hits.empty()) return std::nullopt;
  Location loc;
  loc.tet = hits.front();
  Eigen::Vector3d x0 = impl_->V.row(impl_->T(loc.tet, 0)).transpose();
  Eigen::Matrix3d D;
  for (int k = 0; k < 3; ++k) D.col(k) = impl_->V.row(impl_->T(loc.tet, k + 1)).transpose() - x0;
  Eigen::Vector3d l = D.partialPivLu().solve(q - x0);
  loc.bary << 1.0 - l.sum(), l;
  return loc;
}

std::optional<Eigen::Vector3d> TetLocator::map(const Eigen::Vector3d& q,
                                               const MatX3d& target) const {
  auto loc = locate(q);
  if (!loc) return std::nullopt;
  Eigen::Vector3d p = Eigen::Vector3d::Zero();
  for (int i = 0; i < 4; ++i) p += loc->bary(i) * target.row(impl_->T(loc->tet, i)).transpose();
  return p;
}

}  // namespace pcp
