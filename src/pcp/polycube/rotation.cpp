#include "pcp/polycube/rotation.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/LU>

#include "pcp/core/constants.h"
#include "pcp/geom/triangle.h"

namespace pcp {

MatX3d boundary_area_vectors(const MatX3d& X, const MatX3i& B) {
  MatX3d A(B.rows(), 3);
  for (Index f = 0; f < B.rows(); ++f) {
    Eigen::Vector3d x0 = X.row(B(f, 0)), x1 = X.row(B(f, 1)), x2 = X.row(B(f, 2));
    A.row(f) = tri::area_vector(x0, x1, x2).transpose();
  }
  return A;
}

double rotation_cost(const MatX3d& av, const Eigen::Matrix3d& R, double eps) {
  double num = 0.0, den = 0.0;
  for (Index f = 0; f < av.rows(); ++f) {
    Eigen::Vector3d a = av.row(f).transpose();
    double len2 = a.squaredNorm();
    Eigen::Vector3d b = R * a;
    if (eps > 0.0) {
      for (int c = 0; c < 3; ++c) num += std::sqrt(b(c) * b(c) + eps * len2);
    } else {
      num += b.lpNorm<1>();
    }
    den += std::sqrt(len2);
  }
  return num / den - 1.0;
}

namespace {

// Gradient of the smoothed cost w.r.t. an incremental rotation exp([w]x) applied after R.
Eigen::Vector3d rotation_gradient(const MatX3d& av, const Eigen::Matrix3d& R, double eps) {
  Eigen::Vector3d g = Eigen::Vector3d::Zero();
  double den = 0.0;
  for (Index f = 0; f < av.rows(); ++f) {
    Eigen::Vector3d a = av.row(f).transpose();
    double len2 = a.squaredNorm();
    Eigen::Vector3d b = R * a;
    Eigen::Vector3d coef;
    for (int c = 0; c < 3; ++c) coef(c) = b(c) / std::sqrt(b(c) * b(c) + eps * len2);
    // d b / d w = -cross(b)
    g += -tri::cross(b).transpose() * coef;
    den += std::sqrt(len2);
  }
  return g / den;
}

Eigen::Matrix3d exp_so3(const Eigen::Vector3d& w) {
  double th = w.norm();
  if (th < 1e-300) return Eigen::Matrix3d::Identity();
  return Eigen::AngleAxisd(th, w / th).toRotationMatrix();
}

std::vector<Eigen::Matrix3d> cube_symmetries() {
  std::vector<Eigen::Matrix3d> out;
  int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (auto& p : perm)
    for (int s = 0; s < 8; ++s) {
      Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
      for (int r = 0; r < 3; ++r) M(r, p[r]) = (s >> r) & 1 ? -1.0 : 1.0;
      if (M.determinant() > 0.0) out.push_back(M);
    }
  return out;
}

}  // namespace

double rotation_angle(const Eigen::Matrix3d& R) {
  return std::acos(std::clamp((R.trace() - 1.0) / 2.0, -1.0, 1.0));
}

Eigen::Matrix3d canonicalize_cube_rotation(const Eigen::Matrix3d& R) {
  static const std::vector<Eigen::Matrix3d> sym = cube_symmetries();
  Eigen::Matrix3d best = R;
  double best_angle = rotation_angle(R);
  for (const auto& S : sym) {
    Eigen::Matrix3d C = S * R;
    double a = rotation_angle(C);
    if (a < best_angle - 1e-12) {
      best_angle = a;
      best = C;
    }
  }
  return best;
}

Eigen::Matrix3d refine_rotation(const MatX3d& av, const Eigen::Matrix3d& R0, double eps,
                                int max_iters) {
  Eigen::Matrix3d R = R0;
  double f = rotation_cost(av, R, eps);
  const double h = 1e-5;
  for (int it = 0; it < max_iters; ++it) {
    Eigen::Vector3d g = rotation_gradient(av, R, eps);
    if (g.norm() < 1e-12) break;
    // Finite-difference Hessian of the analytic gradient, made positive definite.
    Eigen::Matrix3d H;
    for (int k = 0; k < 3; ++k) {
      Eigen::Vector3d e = Eigen::Vector3d::Zero();
      e(k) = h;
      H.col(k) = (rotation_gradient(av, exp_so3(e) * R, eps) -
                  rotation_gradient(av, exp_so3(-e) * R, eps)) / (2 * h);
    }
    H = 0.5 * (H + H.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(H);
    Eigen::Vector3d ev = es.eigenvalues();
    double floor = std::max(1e-8, 1e-3 * ev.cwiseAbs().maxCoeff());
    ev = ev.cwiseMax(floor);
    Eigen::Vector3d step = -(es.eigenvectors() * ev.cwiseInverse().asDiagonal() *
                             es.eigenvectors().transpose() * g);
    double max_step = 0.2;  // radians
    if (step.norm() > max_step) step *= max_step / step.norm();
    bool improved = false;
    for (int ls = 0; ls < 20; ++ls) {
      Eigen::Matrix3d Rn = exp_so3(step) * R;
      double fn = rotation_cost(av, Rn, eps);
      if (fn < f) {
        R = Rn;
        f = fn;
        improved = true;
        break;
      }
      step *= 0.5;
    }
    if (!improved || step.norm() < 1e-10) break;
  }
  return R;
}

Eigen::Matrix3d global_orientation_search(const MatX3d& av, double grid_step_deg) {
  const int n = std::max(1, static_cast<int>(std::lround(90.0 / grid_step_deg)));
  Eigen::Matrix3d best = Eigen::Matrix3d::Identity();
  double best_cost = rotation_cost(av, best, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) {
        if (i == 0 && j == 0 && k == 0) continue;
        double az = i * 90.0 / n * kDegToRad, ay = j * 90.0 / n * kDegToRad,
               ax = k * 90.0 / n * kDegToRad;
        Eigen::Matrix3d R = (Eigen::AngleAxisd(az, Eigen::Vector3d::UnitZ()) *
                             Eigen::AngleAxisd(ay, Eigen::Vector3d::UnitY()) *
                             Eigen::AngleAxisd(ax, Eigen::Vector3d::UnitX()))
                                .toRotationMatrix();
        double c = rotation_cost(av, R, 0.0);
        // Strict improvement required so that ties keep the earlier (closer to identity) one.
        if (c < best_cost - 1e-9 * (1.0 + best_cost)) {
          best_cost = c;
          best = R;
        }
      }
  // Refine with a progressively sharper smoothing.
  for (double eps : {1e-2, 1e-3, 1e-4}) {
    Eigen::Matrix3d Rr = refine_rotation(av, best, eps, 15);
    if (rotation_cost(av, Rr, 0.0) < rotation_cost(av, best, 0.0)) best = Rr;
  }
  return canonicalize_cube_rotation(best);
}

}  // namespace pcp
