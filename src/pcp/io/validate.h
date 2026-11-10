#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/core/types.h"

namespace pcp {

// Topological/geometric summary of an input surface.
struct SurfaceReport {
  bool ok = false;
  Index n_vertices = 0;
  Index n_faces = 0;
  Index n_edges = 0;
  Index boundary_edges = 0;      // edges with exactly one incident face
  Index nonmanifold_edges = 0;   // edges with three or more incident faces
  Index degenerate_faces = 0;    // faces with (near) zero area or repeated vertices
  Index unused_vertices = 0;
  int components = 0;            // face-connected components
  long long euler = 0;           // V - E + F (over used vertices)
  int genus = -1;                // total genus, valid only for closed manifold surfaces
  double volume = 0.0;           // signed enclosed volume
  double area = 0.0;
  std::vector<std::string> errors;

  nlohmann::json to_json() const;
};

// Checks that the surface is closed, edge-manifold and free of degenerate faces.
// A face is degenerate if its area is below rel_area_eps * (bbox diagonal)^2.
SurfaceReport validate_surface(const SurfaceMesh& mesh, double rel_area_eps = 1e-14);

}  // namespace pcp
