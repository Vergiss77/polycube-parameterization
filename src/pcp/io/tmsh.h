#pragma once

#include <filesystem>

#include "pcp/core/types.h"

namespace pcp {

// Lossless text format for tetrahedral meshes (17 significant digits), used to save
// stage results and resume the pipeline:
//   TMSH 1
//   nodes N      followed by N lines "x y z"
//   tets M       followed by M lines "a b c d"
//   boundary K   followed by K lines "a b c input_face"
void write_tmsh(const std::filesystem::path& path, const TetMesh& mesh);
TetMesh read_tmsh(const std::filesystem::path& path);

}  // namespace pcp
