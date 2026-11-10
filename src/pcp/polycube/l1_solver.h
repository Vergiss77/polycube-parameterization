#pragma once

#include <functional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/polycube/energies.h"
#include "pcp/polycube/kkt.h"

namespace pcp {

// Parameters of the l1 polycube deformation (defaults follow Huang et al. 2014).
struct L1Params {
  double alpha0 = 0.1;        // initial weight of E1
  double alpha_factor = 2.0;  // alpha multiplier per stage
  double eps0 = 1.0;          // initial l1 smoothing
  double eps_factor = 0.5;    // eps multiplier per stage while eps > eps_min
  double eps_min = 0.01;
  int max_stages = 20;
  int max_iters_per_stage = 50;
  int window = 5;             // stage ends when the last `window` error changes sum ...
  double window_ratio = 0.01; // ... to less than window_ratio * (error at stage start)
  double error_threshold = 1e-3;
  double beta = 0.1;          // complexity weight (1.0 rounds off the existing creases of CAD parts)
  double omega = 0.0;         // feature-edge weight (used only when edges are given)
  double barrier_mu = 1e-6;   // anti-foldover barrier weight
  double tet_barrier_mu = 1e-4;  // tet inversion barrier weight (0 disables)
  bool global_rotation = true;
  double rotation_grid_deg = 6.0;
  double damping0 = 1e-8;     // relative Levenberg-Marquardt damping
  int max_line_search = 30;
  std::string linear_solver = "auto";  // auto | direct | cg
  double cg_tolerance = 1e-4;
  int cg_max_iterations = 100;
  // Truncated Newton: an early-stopped CG step is still a descent step.
  bool log_iterations = false;  // log one line per SQP iteration (verbose diagnostics)

  LinearSolverOptions linear_options() const;
  nlohmann::json to_json() const;
  // Overwrites the fields present in j; throws pcp::Error on unknown keys.
  void update_from_json(const nlohmann::json& j);
};

struct StageLog {
  int stage = 0;
  double alpha = 0.0, eps = 0.0;
  double error = 0.0;  // polycube error at the end of the stage
  double arap = 0.0;   // mean ARAP distortion at the end of the stage
  int iterations = 0;
  double seconds = 0.0;
  long long cg_iterations = 0;  // total CG iterations (0 with the direct solver)
  int cg_max = 0;               // largest CG iteration count of one Newton system
  std::string end_reason;
};

struct StallEvent {
  int stage = 0, iteration = 0;
  Index slim_triangles = 0;  // boundary triangles with an angle below 5 degrees
  double gradient_norm = 0.0;
};

// Wall-clock seconds spent in each part of the solver.
struct L1Timing {
  double rotation = 0.0, assemble = 0.0, factor_solve = 0.0, line_search = 0.0;
  nlohmann::json to_json() const {
    return {{"rotation", rotation}, {"assemble", assemble}, {"factor_solve", factor_solve},
            {"line_search", line_search}};
  }
};

struct L1Result {
  MatX3d X;                  // deformed positions (normalized frame, rotation applied)
  Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();  // total global rotation
  double error = 0.0;        // final polycube error
  double arap_mean = 0.0;    // volume-weighted mean ARAP distortion
  double arap_max = 0.0;
  double area_rel_diff = 0.0;
  double min_normal_slack = 0.0;  // min over adjacent boundary faces of n_i . n_j + 1
  Index inverted_tets = 0;        // tets with non-positive volume in the result
  bool reached_threshold = false;
  int stages = 0, iterations = 0;
  double seconds = 0.0;
  std::vector<StageLog> stage_logs;
  std::vector<StallEvent> stalls;
  Eigen::VectorXd arap_per_tet;
  L1Timing timing;
};

using LogFn = std::function<void(const std::string&)>;

// Runs the scheduled SQP minimization of alpha*E1 + E_delta + beta*E_eta [+ omega*E_e]
// + barrier subject to constant boundary area. `rest` are normalized rest positions (boundary
// area 1) that also serve as the starting point.
L1Result solve_l1_polycube(const MatX3d& rest, const DeformMesh& mesh, const L1Params& params,
                           const std::vector<FeatureEdgeEnergy::Edge>& feature_edges = {},
                           const LogFn& log = {});

// Backtracking line search on a merit function. Returns the accepted step length, or 0 when no
// step in max_halvings halvings satisfies the Armijo condition (or all are infeasible).
double backtracking_line_search(const std::function<double(double)>& merit, double merit0,
                                double directional, int max_halvings);

// Number of boundary triangles with a corner angle below `deg` degrees.
Index count_slim_triangles(const MatX3d& X, const MatX3i& B, double deg = 5.0);

}  // namespace pcp
