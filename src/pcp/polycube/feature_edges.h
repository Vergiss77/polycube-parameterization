#pragma once

#include <filesystem>
#include <vector>

#include "pcp/polycube/energies.h"

namespace pcp {

// Reads user-guided feature edges. Text format, one edge per line (comments with '#'):
//   <vertex a> <vertex b> corner|flat
// Vertex indices refer to the input surface. "corner" asks for a right dihedral angle
// (polycube edge), "flat" asks for the two faces to end up coplanar.
std::vector<FeatureEdgeEnergy::Edge> read_feature_edges(const std::filesystem::path& path,
                                                        const BoundaryTopology& topo);

}  // namespace pcp
