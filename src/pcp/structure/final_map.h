#pragma once

#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/geom/tet_gradient.h"
#include "pcp/structure/patch_graph.h"

namespace pcp {

struct FinalMapParams {
  bool heights_fixed = false;  // keep the given patch heights (used after Chen simplification)
  int max_iterations = 300;
  double rel_tolerance = 1e-7;  // stop when the relative ARAP energy change is below this
  bool laplacian_init = true;   // Huang 2014 eq. (17) initialization
  bool untangle = true;         // locally repair inverted tets after the ARAP iterations
  // Upper bound on the work of the repair, in tet energy evaluations over all passes; the repair
  // stops when it is used up (deterministic, unlike a time limit).
  double untangle_budget = 1e9;
  nlohmann::json to_json() const;
  void update_from_json(const nlohmann::json& j);
};

struct FinalMapResult {
  MatX3d X;                      // exact polycube positions
  std::vector<double> heights;   // final height of each patch
  Eigen::VectorXd arap_per_tet;  // 0.5 ||F - R||^2 per tet
  double arap_mean = 0.0;        // rest-volume weighted
  double arap_max = 0.0;
  Index inverted_tets = 0;       // tets with non-positive volume in the polycube
  Index inverted_before_untangle = 0;
  int iterations = 0;
  double max_planarity_error = 0.0;  // max |x_axis - height| over patch nodes
  double seconds = 0.0;
};

// Re-solves all node positions under the linear constraints "nodes of a patch share the
// patch's coordinate along its axis", minimizing the volumetric ARAP distortion w.r.t. `rest`
// (local-global iterations). `X0` is the starting point (the l1 result); `heights` gives the
// initial (or fixed) patch heights, indexed like g.patches.
FinalMapResult compute_final_map(const MatX3d& rest, const MatX3d& X0, const MatX4i& T,
                                 const MatX3i& B, const PatchGraph& g,
                                 const std::vector<double>& heights, const FinalMapParams& params);

}  // namespace pcp
