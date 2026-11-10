#pragma once

#include <string>

#include <nlohmann/json.hpp>

#include "pcp/core/types.h"

namespace pcp {

struct TetParams {
  double radius_edge_ratio = 1.414;  // TetGen -q radius-edge bound
  double min_dihedral_deg = 0.0;     // TetGen -q minimum dihedral bound (0 = off). With the
                                     // boundary locked (-Y) it cannot remove boundary slivers.
  // Maximum tet volume. <= 0 means automatic: volume of a regular tet whose edge is
  // size_factor * (average input edge length).
  double max_volume = -1.0;
  double size_factor = 2.0;
  // Split interior edges whose two endpoints are both on the boundary (see
  // split_boundary_spanning_edges); needed so that flattening boundary patches cannot force
  // a tet to zero or negative volume.
  bool split_boundary_edges = true;

  nlohmann::json to_json() const;
  void update_from_json(const nlohmann::json& j);  // throws pcp::Error on unknown keys
};

struct TetReport {
  Index n_nodes = 0;
  Index n_tets = 0;
  Index n_boundary = 0;
  double min_dihedral_deg = 0.0;
  double max_dihedral_deg = 0.0;
  Index n_small_dihedral = 0;  // tets with a dihedral angle below 5 degrees
  double volume = 0.0;          // sum of tet volumes
  double surface_volume = 0.0;  // volume enclosed by the input surface
  double max_volume_used = 0.0;
  Index split_edges = 0;  // interior edges split by split_boundary_spanning_edges
  std::string switches;
  nlohmann::json to_json() const;
};

// Boundary-preserving tetrahedralization of a validated closed surface. The boundary of the
// result is exactly the input surface (same vertices, same triangles, B_to_input = face id),
// every tet has positive volume. Throws pcp::Error on failure.
TetMesh tetrahedralize(const SurfaceMesh& surface, const TetParams& params,
                       TetReport* report = nullptr);

// Splits, at their midpoints, all edges of the tet mesh that are not boundary edges but whose
// two endpoints both lie on the boundary. Afterwards no tet can have all four vertices on one
// flat boundary patch, and no tet spans across a crease without containing it. The boundary
// triangles are unchanged; new nodes are interior. Returns the number of split edges.
Index split_boundary_spanning_edges(TetMesh& mesh);

// Signed volume of tet t (positive for the orientation used throughout the library).
double tet_volume(const MatX3d& X, const MatX4i& T, Index t);

struct DihedralStats {
  double min_deg = 180.0;
  double max_deg = 0.0;
  Index n_below_5deg = 0;
};
// Dihedral angle statistics over all tets.
DihedralStats dihedral_stats(const MatX3d& X, const MatX4i& T);

}  // namespace pcp
