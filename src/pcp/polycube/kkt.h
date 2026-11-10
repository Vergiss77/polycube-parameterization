#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace pcp {

// Linear solver used for the Newton systems.
struct LinearSolverOptions {
  enum class Kind { Auto, Direct, Iterative };
  Kind kind = Kind::Auto;
  // Auto uses the direct solver below this many unknowns, preconditioned CG above. A sparse
  // Cholesky of a 3D tet mesh fills in heavily, so for ~1e5 unknowns an inexact CG solve is
  // more than an order of magnitude faster.
  long direct_max_unknowns = 30000;
  double cg_tolerance = 1e-4;  // relative residual of each inexact Newton solve
  int cg_max_iterations = 100;
};

// Solves the equality-constrained Newton system of one SQP iteration
//   [ H + mu*s*I   gc ] [dx]   [ -g ]
//   [ gc^T         0  ] [lam] = [ -C ]
// via the Schur complement: two solves with the damped H (sharing one factorization or one
// preconditioner). s is the mean diagonal of H, so `damping` is relative. For the direct
// solver the sparsity pattern of H must stay fixed between calls (symbolic analysis once).
class KktSolver {
 public:
  explicit KktSolver(LinearSolverOptions opt = {}) : opt_(opt) {}
  bool solve(const Eigen::SparseMatrix<double>& H, double damping, const Eigen::VectorXd& g,
             const Eigen::VectorXd& gc, double C, Eigen::VectorXd& dx, double& lambda);
  // CG iterations used by the last solve (0 for the direct solver).
  int last_iterations() const { return last_iterations_; }
  bool uses_direct(Eigen::Index n) const;

 private:
  LinearSolverOptions opt_;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                           Eigen::IncompleteCholesky<double>>
      cg_;
  Eigen::SparseMatrix<double> work_;
  bool analyzed_ = false;
  int last_iterations_ = 0;
};

}  // namespace pcp
