#pragma once

// G+Smo MultiPatch XML (design D7). Format checked against the G+Smo reader (branch "stable"):
// src/gsIO/gsXml.cpp (getInterfaces / getBoundaries), src/gsIO/gsXmlUtils.hpp
// (gsXml<gsMultiPatch>::get) and src/gsCore/gsBoundary.h (box sides).
// - Geometry ids 0..N-1, <patches type="id_range">0 N-1</patches>;
// - box sides 1..6 = u=0, u=1, v=0, v=1, w=0, w=1 (west, east, south, north, front, back);
// - each interface: patch1 side1 patch2 side2, direction map (3 ints), orientation (3 flags).
//   All blocks share the global axis orientation: the map is the identity, the tangential
//   flags are 1 and the flag of the face normal direction is 0 (G+Smo recomputes it as
//   "both sides at the same parameter value", false for a high side glued to a low side).

#include <array>
#include <filesystem>
#include <utility>
#include <vector>

#include "pcp/blocks/block_layout.h"
#include "pcp/blocks/bspline.h"

namespace pcp {

struct GismoMultiPatch {
  std::vector<TensorBSpline3> volumes;
  std::vector<std::array<int, 10>> interfaces;
  std::vector<std::pair<int, int>> boundary;  // (patch, side)
};

// Interfaces and boundary sides of the block layout in G+Smo numbering.
GismoMultiPatch gismo_topology(const BlockLayout& layout, std::vector<TensorBSpline3> volumes);

void write_gismo_xml(const std::filesystem::path& path, const GismoMultiPatch& mp);

// Reads files written by write_gismo_xml (not a general XML parser). Throws pcp::Error.
GismoMultiPatch read_gismo_xml(const std::filesystem::path& path);

}  // namespace pcp
