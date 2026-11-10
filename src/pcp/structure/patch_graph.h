#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/structure/labeling.h"

namespace pcp {

// Maximal chain of sharp edges (edges between two different patches) between two corners, or
// a closed loop without corners.
struct Chain {
  std::vector<int> vertices;  // ordered; for a closed loop the first vertex is not repeated
  bool closed = false;
  int patch_a = -1, patch_b = -1;  // patches on both sides
  int axis = -1;                   // direction the chain runs along (third axis), -1 if invalid
  int corner0 = -1, corner1 = -1;  // end vertices (corners), -1 for closed loops
  // Patches of orientation `axis` at the end with the smaller / larger coordinate along
  // `axis` (Chen 2019 uses them as directed edges lower -> higher). -1 if none.
  int lower_patch = -1, upper_patch = -1;
  double length = 0.0;
};

struct PatchGraph {
  Patches patches;
  std::vector<double> height;  // patch coordinate along its axis (area-weighted mean)
  std::vector<double> area;
  std::vector<int> corners;    // vertex ids where three or more patches meet
  std::vector<std::vector<int>> corner_patches;
  std::vector<Chain> chains;
  nlohmann::json to_json(const MatX3d& X) const;
};

PatchGraph build_patch_graph(const MatX3d& X, const MatX3i& B, const BoundaryTopology& topo,
                             const Eigen::VectorXi& labels);

struct TopologyCheck {
  bool ok = true;
  std::vector<std::string> violations;
  long long euler_surface = 0;  // V - E + F of the surface
  long long euler_complex = 0;  // sum chi(patch) - sum chi(chain) + #corners
  int non_disk_patches = 0;     // informational: patches with chi != 1
  nlohmann::json to_json() const;
};

// Checks: adjacent patches lie on different axes, every patch has at least 3 neighbors, every
// chain has a valid axis, and the Euler characteristic recomputed from the cell structure
// (patches, chains, corners) equals that of the surface.
TopologyCheck check_polycube_topology(const MatX3i& B, const PatchGraph& g);

}  // namespace pcp
