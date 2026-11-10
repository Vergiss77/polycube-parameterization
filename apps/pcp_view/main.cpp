// pcp_view: optional polyscope viewer for a pcp_cli output directory.
// Shows the original surface and the final polycube surface colored by patch, and the
// polycube tet mesh with its per-tet ARAP distortion.
#include <cstdio>
#include <filesystem>

#include <nlohmann/json.hpp>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>

#include "pcp/io/text.h"
#include "pcp/io/tmsh.h"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: pcp_view <pcp_cli output directory>\n");
    return 1;
  }
  try {
    const fs::path dir = argv[1];
    pcp::TetMesh original = pcp::read_tmsh(dir / "01_tet.tmsh");
    pcp::TetMesh polycube = pcp::read_tmsh(dir / "04_polycube.tmsh");
    nlohmann::json graph = nlohmann::json::parse(pcp::text::read_file(dir / "03_patch_graph.json"));

    Eigen::VectorXi patch = Eigen::VectorXi::Zero(original.B.rows());
    for (const auto& p : graph["patches"])
      for (int f : p["faces"]) patch(f) = p["id"].get<int>();
    Eigen::VectorXd patch_d = patch.cast<double>();

    polyscope::init();
    auto* a = polyscope::registerSurfaceMesh("original", original.X, original.B);
    a->addFaceScalarQuantity("patch", patch_d)->setEnabled(true);
    // Place the polycube next to the original.
    pcp::MatX3d shifted = polycube.X;
    double width = original.X.col(0).maxCoeff() - original.X.col(0).minCoeff();
    shifted.col(0).array() += 1.3 * width;
    auto* b = polyscope::registerSurfaceMesh("polycube", shifted, polycube.B);
    b->addFaceScalarQuantity("patch", patch_d)->setEnabled(true);
    polyscope::registerTetMesh("polycube tets", shifted, polycube.T)->setEnabled(false);
    polyscope::show();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
