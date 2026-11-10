#pragma once

#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/core/types.h"
#include "pcp/geom/boundary_topology.h"

namespace pcp {

// Direction labels: 0:+X 1:-X 2:+Y 3:-Y 4:+Z 5:-Z.
inline int label_axis(int label) { return label / 2; }
inline int label_sign(int label) { return label % 2 == 0 ? 1 : -1; }
inline int make_label(int axis, int sign) { return 2 * axis + (sign < 0 ? 1 : 0); }

// Initial labels from the largest normal component of each boundary face.
Eigen::VectorXi initial_labels(const MatX3d& X, const MatX3i& B);

// Removes zigzags: a face whose label differs from two neighbors that share a common label
// takes that label. Repeats until stable. Returns the number of relabeled faces.
Index remove_zigzags(const BoundaryTopology& topo, Eigen::VectorXi& labels);

// Patches = edge-connected components of equally labeled faces.
struct Patches {
  Eigen::VectorXi face_patch;               // patch id of each face
  std::vector<int> label;                   // label of each patch
  std::vector<std::vector<int>> faces;      // faces of each patch
  std::vector<std::vector<int>> neighbors;  // adjacent patches (sorted, unique)
};
Patches compute_patches(const BoundaryTopology& topo, const Eigen::VectorXi& labels);

struct CleanupReport {
  Index zigzag_faces = 0;
  int degree1_merged = 0;
  int degree2_merged = 0;
  int conflict_merged = 0;  // small patches adjacent to a larger patch on the same axis
  nlohmann::json to_json() const;
};

// Zigzag removal followed by iterative merging of degree-1 patches (into their neighbor),
// degree-2 patches, and patches adjacent to a larger patch on the same axis (e.g. a nearly
// folded triangle at a corner), the latter two into the neighbor sharing the longest boundary
// (measured on X). Repeats until no patch has fewer than 3 neighbors or a same-axis neighbor.
CleanupReport cleanup_labels(const MatX3d& X, const MatX3i& B, const BoundaryTopology& topo,
                             Eigen::VectorXi& labels);

}  // namespace pcp
