#pragma once

#include <filesystem>

#include "pcp/core/types.h"

namespace pcp {

// Reads a triangle-only OFF file. Throws pcp::Error with a description on any problem
// (bad header, non-triangle face, out-of-range index, truncated data).
SurfaceMesh read_off(const std::filesystem::path& path);

// Writes an OFF file with 17 significant digits so that a read-back is lossless.
void write_off(const std::filesystem::path& path, const SurfaceMesh& mesh);

}  // namespace pcp
