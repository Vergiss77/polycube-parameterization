#include "pcp/polycube/kkt.h"

#include <cmath>

namespace pcp {

bool KktSolver::uses_direct(Eigen::Index n) const {
  switch (opt_.kind) {
    case LinearSolverOptions::Kind::Direct:
      return true;
    case LinearSolverOptions::Kind::Iterative:
      return false;
    default:
      return n <= opt_.direct_max_unknowns;
  }
}

bool KktSolver::solve(const Eigen::SparseMatrix<double>& H, double damping,
                      const Eigen::VectorXd& g, const Eigen::VectorXd& gc, double C,
                      Eigen::VectorXd& dx, double& lambda) {
  work_ = H;
  const Eigen::Index n = work_.rows();
  double diag_mean = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) diag_mean += work_.coeff(i, i);
  diag_mean /= static_cast<double>(n);
  if (!(diag_mean > 0.0)) diag_mean = 1.0;
  const double shift = damping * diag_mean;
  for (Eigen::Index i = 0; i < n; ++i) work_.coeffRef(i, i) += shift;

  // The gradient is often dominated by a component along the constraint gradient that the
  // multiplier cancels (E1 grows with area, the area constraint pushes back). Remove it first
  // with the least-squares multiplier lam0, so that the two (possibly inexact) solves do not
  // produce a small step as the difference of two large, noisy vectors. With exact solves this
  // is the same step.
  const double gc2 = gc.squaredNorm();
  const double lam0 = gc2 > 0.0 ? -g.dot(gc) / gc2 : 0.0;
  const Eigen::VectorXd gt = g + lam0 * gc;

  Eigen::VectorXd yg, yc;
  last_iterations_ = 0;
  if (uses_direct(n)) {
    if (!analyzed_) {
      ldlt_.analyzePattern(work_);
      analyzed_ = true;
    }
    ldlt_.factorize(work_);
    if (ldlt_.info() != Eigen::Success) return false;
    yg = ldlt_.solve(gt);
    yc = ldlt_.solve(gc);
    if (ldlt_.info() != Eigen::Success) return false;
  } else {
    cg_.setTolerance(opt_.cg_tolerance);
    cg_.setMaxIterations(opt_.cg_max_iterations);
    cg_.compute(work_);
    if (cg_.info() != Eigen::Success) return false;
    yg = cg_.solve(gt);
    last_iterations_ = static_cast<int>(cg_.iterations());
    yc = cg_.solve(gc);
    last_iterations_ += static_cast<int>(cg_.iterations());
    // NoConvergence after the iteration cap still gives a usable (inexact) direction.
    if (cg_.info() == Eigen::NumericalIssue) return false;
  }
  const double denom = gc.dot(yc);
  double dlam = 0.0;
  if (std::abs(denom) > 0.0 && std::isfinite(denom)) dlam = (C - gc.dot(yg)) / denom;
  dx = -yg - dlam * yc;
  lambda = lam0 + dlam;
  return dx.allFinite();
}

}  // namespace pcp
