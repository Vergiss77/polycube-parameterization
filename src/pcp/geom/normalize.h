#pragma once

#include "pcp/core/types.h"

namespace pcp {

// Similarity transform x -> (x - center) * scale. Used to make the optimization
// scale-independent: after applying it the boundary area of the mesh is exactly 1.
struct Normalization {
  Eigen::RowVector3d center = Eigen::RowVector3d::Zero();
  double scale = 1.0;

  MatX3d apply(const MatX3d& X) const;
  MatX3d inverse(const MatX3d& Y) const;
};

// Total area of the triangles F over positions X.
double total_area(const MatX3d& X, const MatX3i& F);

// Centroid (area-weighted over the boundary) translation and scale 1/sqrt(boundary area).
Normalization normalization_for(const MatX3d& X, const MatX3i& boundary);

}  // namespace pcp
