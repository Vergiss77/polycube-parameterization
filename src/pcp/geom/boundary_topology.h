#pragma once

#include <array>
#include <vector>

#include "pcp/core/types.h"

namespace pcp {

// Pair of boundary triangles sharing an edge.
struct FacePair {
  int f0 = -1, f1 = -1;      // face indices into B
  int v0 = -1, v1 = -1;      // shared edge vertices
  std::array<int, 4> nodes;  // v0, v1, opposite vertex in f0, opposite vertex in f1
};

struct BoundaryTopology {
  std::vector<FacePair> pairs;                 // one per interior edge of the closed surface
  std::vector<std::array<int, 3>> face_pairs;  // per face: indices into `pairs` (edge k = (k,k+1))
  std::vector<std::array<int, 3>> face_nbrs;   // per face: neighbor across edge k
};

// Requires a closed edge-manifold triangle surface.
BoundaryTopology build_boundary_topology(const MatX3i& B);

}  // namespace pcp
