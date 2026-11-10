#pragma once

// Multi-block trivariate B-spline fitting (design D3-D6): one tensor product volume per block,
// shared control points on shared faces (C0 by construction), global least squares to samples
// of the piecewise-linear map polycube -> original tet mesh.

#include <array>
#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/blocks/block_layout.h"
#include "pcp/blocks/bspline.h"

namespace pcp {

struct FitParams {
  int degree = 3;
  double element_size = 0.05;   // target element size relative to the polycube bounding box diagonal
  int samples_per_element = 4;  // fitting samples per element and direction
  double smoothing = 1e-6;      // weight of the control net smoothing term (dimensionless)
  nlohmann::json to_json() const;
  void update_from_json(const nlohmann::json& j);  // throws pcp::Error on unknown keys / bad values
};

// max(1, round(length / h)).
int elements_for_length(double length, double h);

// Spline space shared by all blocks (D4). Along every axis, interval k has elements[k] elements
// and its control points start at the global 1D offset offset[k]; offset[k + 1] = offset[k] +
// elements[k] + degree - 1, so neighboring intervals share one layer of control points.
struct SplineSpace {
  int degree = 3;
  std::array<std::vector<int>, 3> elements, offset;
  std::array<int, 3> lattice{};             // global lattice size per axis
  int n_unknowns = 0;
  std::vector<std::vector<int>> block_dofs;  // per block: local control index (u fastest) -> unknown

  std::array<KnotVector, 3> block_basis(const Block& b) const;
};

// h is the absolute element size (polycube units).
SplineSpace build_spline_space(const BlockLayout& layout, int degree, double h);

// Fitting samples (D3): parameter points of every block and their images under the map.
struct FitSamples {
  std::vector<int> block;
  std::vector<Eigen::Vector3d> param;   // (u, v, w) in [0, 1]^3
  std::vector<Eigen::Vector3d> target;  // physical position
  std::vector<char> surface;            // the sample lies on a boundary face of its block
  int unlocated = 0;                    // samples found in no tet (nearest tet used)
};

// samples_per_element s: m = n * s + 1 parameter values per direction (ends included). X_poly
// is the axis-aligned polycube, X_target the embedding the volumes approximate (same tets T).
FitSamples sample_blocks(const BlockLayout& layout, const SplineSpace& space, int samples_per_element,
                         const MatX3d& X_poly, const MatX4i& T, const MatX3d& X_target);

struct FitSolution {
  MatX3d controls;                      // one row per unknown
  std::vector<TensorBSpline3> volumes;  // one per block
};

// Global least squares (D5): (A^T A + smoothing * kappa * S) c = A^T b, kappa = tr(A^T A) / tr(S),
// S the sum of squared second differences of every block's control net. Throws pcp::Error if the
// system cannot be factorized.
FitSolution fit_volumes(const BlockLayout& layout, const SplineSpace& space, const FitSamples& samples,
                        double smoothing);

// Evaluation grid of one block: n * 4 + 1 points per direction (element boundaries included).
std::array<int, 3> eval_grid_size(const SplineSpace& space, const Block& b);

// Quality report (D6). Distances are divided by `diag` (bounding box diagonal of the target).
struct FitQuality {
  double rms_error = 0.0, max_error = 0.0;                  // all fitting samples
  double surface_rms_error = 0.0, surface_max_error = 0.0;  // samples on boundary faces
  double surface_distance_max = 0.0, surface_distance_mean = 0.0;  // boundary faces -> surface
  double min_scaled_jacobian = 0.0, mean_scaled_jacobian = 0.0;
  long long eval_points = 0, nonpositive_jacobian = 0;
  int blocks_with_nonpositive = 0;
  double c0_max_gap = 0.0;
  nlohmann::json to_json() const;
};

FitQuality evaluate_fit(const BlockLayout& layout, const SplineSpace& space, const FitSolution& sol,
                        const FitSamples& samples, const MatX3d& surface_V, const MatX3i& surface_F,
                        double diag);

// Evaluation grids of all blocks as hexahedra (design D8): points are not shared between
// blocks; cell field `block`, point field `scaled_jacobian`.
struct SampleMesh {
  MatX3d X;
  Eigen::Matrix<int, Eigen::Dynamic, 8> H;
  Eigen::VectorXi block;           // per cell
  Eigen::VectorXd scaled_jacobian;  // per point
};
SampleMesh sample_mesh(const BlockLayout& layout, const SplineSpace& space, const FitSolution& sol);

// Scaled Jacobian det J / (|J_u| |J_v| |J_w|) (0 for a degenerate column).
double scaled_jacobian(const Eigen::Matrix3d& J);

}  // namespace pcp
