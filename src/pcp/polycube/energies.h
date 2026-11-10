#pragma once

#include <array>
#include <vector>

#include "pcp/polycube/deform_mesh.h"

namespace pcp {

// Common interface of the energy terms. Gradients and Hessians are accumulated (+=) into the
// given containers, already multiplied by `weight`. Hessians are symmetric positive
// semi-definite approximations.
class EnergyTerm {
 public:
  explicit EnergyTerm(const DeformMesh& mesh) : mesh_(mesh) {}
  virtual ~EnergyTerm() = default;

  // Refresh quantities that are lagged (held constant) during one SQP iteration.
  virtual void begin_iteration(const MatX3d& X) {}
  virtual double value(const MatX3d& X) const = 0;
  virtual void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const = 0;
  virtual void add_hessian(const MatX3d& X, BlockSparse& H) const = 0;

  double weight = 1.0;

 protected:
  const DeformMesh& mesh_;
};

// Normalized polycube error d = (sum_f A_f ||n_f||_1 - sum_f A_f) / sum_f A_f.
double polycube_error(const MatX3d& X, const MatX3i& B);

// Boundary area constraint C(X) = area(X) - target and its gradient.
double area_constraint(const MatX3d& X, const MatX3i& B, double target);
void area_constraint_gradient(const MatX3d& X, const MatX3i& B, Eigen::VectorXd& g);

// E1: smoothed l1 norm of area vectors, sum_f sum_c sqrt(a_c^2 + eps * Abar_f^2) / A_rest,
// with Abar_f the face area at the start of the iteration (lagged).
class L1NormalEnergy : public EnergyTerm {
 public:
  explicit L1NormalEnergy(const DeformMesh& mesh);
  void begin_iteration(const MatX3d& X) override;
  double value(const MatX3d& X) const override;
  void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const override;
  void add_hessian(const MatX3d& X, BlockSparse& H) const override;
  double eps = 1.0;

 private:
  Eigen::VectorXd lagged_area_;
  std::array<Eigen::Matrix<double, 9, 9>, 3> hess_pos_, hess_neg_;  // PSD parts of +/-H_c
};

// E_delta: volume-weighted ARAP distortion, sum_t V_t * 0.5 ||F_t - R_t||^2 / sum_t V_t.
class ArapEnergy : public EnergyTerm {
 public:
  using EnergyTerm::EnergyTerm;
  double value(const MatX3d& X) const override;
  void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const override;
  void add_hessian(const MatX3d& X, BlockSparse& H) const override;
  // Per-tet distortion delta_t = 0.5 ||F - R||^2 (unweighted).
  Eigen::VectorXd per_tet(const MatX3d& X) const;
};

// E_eta: complexity penalty sum_e A(e) ||n_i - n_j||^2 / A_rest (Huang 2014, eqs. 8-10),
// with edge areas A(e) lagged. Weight = beta.
class ComplexityEnergy : public EnergyTerm {
 public:
  using EnergyTerm::EnergyTerm;
  void begin_iteration(const MatX3d& X) override;
  double value(const MatX3d& X) const override;
  void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const override;
  void add_hessian(const MatX3d& X, BlockSparse& H) const override;
  const Eigen::VectorXd& edge_areas() const { return edge_area_; }

 private:
  Eigen::VectorXd edge_area_;  // per face pair
};

// E_e: user-guided edges. Each entry is a face pair (index into topo.pairs) and a target
// dot product of the two normals: 0 for a right angle (polycube edge), 1 for flat.
// Energy sum (n_i . n_j - target)^2. Weight = omega.
class FeatureEdgeEnergy : public EnergyTerm {
 public:
  struct Edge {
    int pair = -1;
    double target = 0.0;
  };
  FeatureEdgeEnergy(const DeformMesh& mesh, std::vector<Edge> edges)
      : EnergyTerm(mesh), edges_(std::move(edges)) {}
  double value(const MatX3d& X) const override;
  void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const override;
  void add_hessian(const MatX3d& X, BlockSparse& H) const override;
  bool empty() const { return edges_.empty(); }

 private:
  std::vector<Edge> edges_;
};

// Anti-foldover log barrier -sum_e log(n_i . n_j + 1) over adjacent boundary faces.
// Weight = mu. value() is +infinity when any pair is infeasible.
class FoldoverBarrier : public EnergyTerm {
 public:
  using EnergyTerm::EnergyTerm;
  double value(const MatX3d& X) const override;
  void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const override;
  void add_hessian(const MatX3d& X, BlockSparse& H) const override;
  // Minimum of n_i . n_j + 1 over all adjacent pairs.
  double min_slack(const MatX3d& X) const;
};

// Tet inversion barrier -sum_t (V_t / V) log(vol_t / V_t) with rest volumes V_t (V = total).
// Keeps every tet positively oriented during the deformation (value() is +infinity otherwise).
// Weight = mu. Hessian: the positive semi-definite Gauss-Newton part.
class TetInversionBarrier : public EnergyTerm {
 public:
  using EnergyTerm::EnergyTerm;
  double value(const MatX3d& X) const override;
  void add_gradient(const MatX3d& X, Eigen::VectorXd& g) const override;
  void add_hessian(const MatX3d& X, BlockSparse& H) const override;
  Index inverted(const MatX3d& X) const;
};

// Unit normals of the two faces of a pair and the Jacobian of (n0, n1) w.r.t. the 12
// stacked coordinates of pair.nodes. Returns false if a face is degenerate.
bool pair_normals(const MatX3d& X, const FacePair& p, Eigen::Vector3d& n0, Eigen::Vector3d& n1,
                  Eigen::Matrix<double, 3, 12>& J0, Eigen::Matrix<double, 3, 12>& J1);

}  // namespace pcp
