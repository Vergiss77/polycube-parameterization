#pragma once

#include <memory>
#include <optional>

#include "pcp/core/types.h"

namespace pcp {

// Point location in a tet mesh (AABB tree) and evaluation of the piecewise-linear map from that
// mesh ("domain", e.g. the polycube) to another embedding of the same tets ("target", e.g. the
// original shape).
class TetLocator {
 public:
  struct Location {
    Index tet = -1;
    Eigen::Vector4d bary;  // barycentric coordinates w.r.t. T(tet, 0..3)
  };

  TetLocator(const MatX3d& domain, const MatX4i& T);
  ~TetLocator();
  TetLocator(const TetLocator&) = delete;
  TetLocator& operator=(const TetLocator&) = delete;

  // Tet containing q (first found), or nullopt if q is outside the domain mesh.
  std::optional<Location> locate(const Eigen::Vector3d& q) const;
  // Image of q under the map domain -> target, or nullopt if q is outside the domain.
  std::optional<Eigen::Vector3d> map(const Eigen::Vector3d& q, const MatX3d& target) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pcp
