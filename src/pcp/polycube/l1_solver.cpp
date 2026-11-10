#include "pcp/polycube/l1_solver.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>

#include <Eigen/Geometry>

#include "pcp/core/constants.h"
#include "pcp/polycube/kkt.h"
#include "pcp/polycube/rotation.h"

namespace pcp {

nlohmann::json L1Params::to_json() const {
  return {{"alpha0", alpha0},
          {"alpha_factor", alpha_factor},
          {"eps0", eps0},
          {"eps_factor", eps_factor},
          {"eps_min", eps_min},
          {"max_stages", max_stages},
          {"max_iters_per_stage", max_iters_per_stage},
          {"window", window},
          {"window_ratio", window_ratio},
          {"error_threshold", error_threshold},
          {"beta", beta},
          {"omega", omega},
          {"barrier_mu", barrier_mu},
          {"tet_barrier_mu", tet_barrier_mu},
          {"global_rotation", global_rotation},
          {"rotation_grid_deg", rotation_grid_deg},
          {"damping0", damping0},
          {"max_line_search", max_line_search},
          {"linear_solver", linear_solver},
          {"cg_tolerance", cg_tolerance},
          {"cg_max_iterations", cg_max_iterations},
          {"log_iterations", log_iterations}};
}

LinearSolverOptions L1Params::linear_options() const {
  LinearSolverOptions o;
  if (linear_solver == "auto")
    o.kind = LinearSolverOptions::Kind::Auto;
  else if (linear_solver == "direct")
    o.kind = LinearSolverOptions::Kind::Direct;
  else if (linear_solver == "cg")
    o.kind = LinearSolverOptions::Kind::Iterative;
  else
    throw Error("unknown linear_solver '" + linear_solver + "' (expected auto, direct or cg)");
  o.cg_tolerance = cg_tolerance;
  o.cg_max_iterations = cg_max_iterations;
  return o;
}

void L1Params::update_from_json(const nlohmann::json& j) {
  const nlohmann::json known = to_json();
  for (auto it = j.begin(); it != j.end(); ++it)
    if (!known.contains(it.key())) throw Error("unknown l1 parameter '" + it.key() + "'");
  auto get = [&](const char* k, auto& v) {
    if (j.contains(k)) v = j.at(k).get<std::decay_t<decltype(v)>>();
  };
  get("alpha0", alpha0);
  get("alpha_factor", alpha_factor);
  get("eps0", eps0);
  get("eps_factor", eps_factor);
  get("eps_min", eps_min);
  get("max_stages", max_stages);
  get("max_iters_per_stage", max_iters_per_stage);
  get("window", window);
  get("window_ratio", window_ratio);
  get("error_threshold", error_threshold);
  get("beta", beta);
  get("omega", omega);
  get("barrier_mu", barrier_mu);
  get("tet_barrier_mu", tet_barrier_mu);
  get("global_rotation", global_rotation);
  get("rotation_grid_deg", rotation_grid_deg);
  get("damping0", damping0);
  get("max_line_search", max_line_search);
  get("linear_solver", linear_solver);
  get("cg_tolerance", cg_tolerance);
  get("cg_max_iterations", cg_max_iterations);
  get("log_iterations", log_iterations);
  linear_options();  // validates linear_solver
}

double backtracking_line_search(const std::function<double(double)>& merit, double merit0,
                                double directional, int max_halvings) {
  const double c1 = 1e-4;
  double t = 1.0;
  for (int i = 0; i <= max_halvings; ++i) {
    double m = merit(t);
    if (std::isfinite(m)) {
      bool armijo = directional < 0.0 ? m <= merit0 + c1 * t * directional : m < merit0;
      if (armijo) return t;
    }
    t *= 0.5;
  }
  return 0.0;
}

Index count_slim_triangles(const MatX3d& X, const MatX3i& B, double deg) {
  const double cos_limit = std::cos(deg * kDegToRad);
  Index n = 0;
  for (Index f = 0; f < B.rows(); ++f) {
    for (int k = 0; k < 3; ++k) {
      Eigen::Vector3d p = X.row(B(f, k)), a = X.row(B(f, (k + 1) % 3)), b = X.row(B(f, (k + 2) % 3));
      Eigen::Vector3d u = a - p, v = b - p;
      double nu = u.norm(), nv = v.norm();
      if (nu == 0.0 || nv == 0.0 || u.dot(v) / (nu * nv) > cos_limit) {
        ++n;
        break;
      }
    }
  }
  return n;
}

namespace {

Eigen::RowVector3d boundary_centroid(const MatX3d& X, const MatX3i& B) {
  Eigen::RowVector3d c = Eigen::RowVector3d::Zero();
  double area = 0.0;
  for (Index f = 0; f < B.rows(); ++f) {
    Eigen::Vector3d x0 = X.row(B(f, 0)), x1 = X.row(B(f, 1)), x2 = X.row(B(f, 2));
    double a = 0.5 * (x1 - x0).cross(x2 - x0).norm();
    c += a * (x0 + x1 + x2).transpose() / 3.0;
    area += a;
  }
  return c / area;
}

void rotate_about_centroid(MatX3d& X, const MatX3i& B, const Eigen::Matrix3d& R) {
  Eigen::RowVector3d c = boundary_centroid(X, B);
  X = ((X.rowwise() - c) * R.transpose()).rowwise() + c;
}

template <typename... Args>
std::string fmt(const char* f, Args... args) {
  char buf[512];
  std::snprintf(buf, sizeof(buf), f, args...);
  return buf;
}

}  // namespace

L1Result solve_l1_polycube(const MatX3d& rest, const DeformMesh& mesh, const L1Params& P,
                           const std::vector<FeatureEdgeEnergy::Edge>& feature_edges,
                           const LogFn& log) {
  using clock = std::chrono::steady_clock;
  const auto t_start = clock::now();
  auto say = [&](const std::string& s) {
    if (log) log(s);
  };

  L1Result res;
  MatX3d X = rest;

  L1NormalEnergy e1(mesh);
  ArapEnergy arap(mesh);
  ComplexityEnergy eta(mesh);
  FeatureEdgeEnergy feat(mesh, feature_edges);
  FoldoverBarrier barrier(mesh);
  TetInversionBarrier tet_barrier(mesh);
  tet_barrier.weight = P.tet_barrier_mu;
  eta.weight = P.beta;
  feat.weight = P.omega;
  barrier.weight = P.barrier_mu;
  std::vector<EnergyTerm*> terms = {&e1, &arap, &eta, &barrier};
  if (!feat.empty() && P.omega > 0.0) terms.push_back(&feat);
  if (P.tet_barrier_mu > 0.0) terms.push_back(&tet_barrier);

  if (P.global_rotation) {
    Eigen::Matrix3d R0 = global_orientation_search(boundary_area_vectors(X, mesh.B),
                                                   P.rotation_grid_deg);
    rotate_about_centroid(X, mesh.B, R0);
    res.rotation = R0;
    say(fmt("initial orientation: rotation angle %.3f deg, polycube error %.6g",
            rotation_angle(R0) * kRadToDeg, polycube_error(X, mesh.B)));
  }

  if (barrier.min_slack(X) <= 0.0)
    throw Error("l1 solver: the input already contains folded adjacent boundary faces");
  if (P.tet_barrier_mu > 0.0 && tet_barrier.inverted(X) > 0)
    throw Error("l1 solver: the input tet mesh contains inverted tets");

  BlockSparse H = mesh.make_hessian();
  KktSolver kkt(P.linear_options());
  Eigen::VectorXd g(3 * X.rows()), gc;
  double damping = P.damping0;
  double rho = 1.0;
  double alpha = P.alpha0, eps = P.eps0;

  auto total_value = [&](const MatX3d& Y) {
    double v = 0.0;
    for (EnergyTerm* t : terms) {
      v += t->value(Y);
      if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
    }
    return v;
  };

  auto since = [&](clock::time_point t0) {
    return std::chrono::duration<double>(clock::now() - t0).count();
  };
  if (polycube_error(X, mesh.B) < P.error_threshold) {
    res.reached_threshold = true;
    say(fmt("initial polycube error %.3g is already below the threshold; no deformation needed",
            polycube_error(X, mesh.B)));
  }

  for (int stage = 0; stage < P.max_stages && !res.reached_threshold; ++stage) {
    const auto t_stage = clock::now();
    e1.weight = alpha;
    e1.eps = eps;
    const double d0 = polycube_error(X, mesh.B);
    double d_prev = d0;
    std::vector<double> changes;
    StageLog sl;
    sl.stage = stage;
    sl.alpha = alpha;
    sl.eps = eps;
    sl.end_reason = "iteration limit";

    int it = 0;
    for (; it < P.max_iters_per_stage; ++it) {
      auto t0 = clock::now();
      if (P.global_rotation) {
        Eigen::Matrix3d Rinc = refine_rotation(boundary_area_vectors(X, mesh.B),
                                               Eigen::Matrix3d::Identity(), std::max(eps, 1e-4), 3);
        if (rotation_angle(Rinc) > 0.0) {
          rotate_about_centroid(X, mesh.B, Rinc);
          res.rotation = Rinc * res.rotation;
        }
      }
      res.timing.rotation += since(t0);
      t0 = clock::now();
      for (EnergyTerm* t : terms) t->begin_iteration(X);

      g.setZero();
      for (EnergyTerm* t : terms) t->add_gradient(X, g);
      H.set_zero();
      for (EnergyTerm* t : terms) t->add_hessian(X, H);
      const double C = area_constraint(X, mesh.B, mesh.rest_area);
      area_constraint_gradient(X, mesh.B, gc);
      const double f0 = total_value(X);
      res.timing.assemble += since(t0);

      Eigen::VectorXd dx;
      double lambda = 0.0, step = 0.0;
      bool ok = false;
      for (int attempt = 0; attempt < 4 && !ok; ++attempt) {
        t0 = clock::now();
        bool solved = kkt.solve(H.matrix(), damping, g, gc, C, dx, lambda);
        res.timing.factor_solve += since(t0);
        sl.cg_iterations += kkt.last_iterations();
        sl.cg_max = std::max(sl.cg_max, kkt.last_iterations());
        if (!solved) {
          damping *= 10.0;
          continue;
        }
        rho = std::max(rho, 2.0 * std::abs(lambda));
        const double merit0 = f0 + rho * std::abs(C);
        const double directional = g.dot(dx) - rho * std::abs(C);
        auto merit = [&](double t) {
          MatX3d Y = X + t * unflatten(dx);
          double v = total_value(Y);
          if (!std::isfinite(v)) return v;
          return v + rho * std::abs(area_constraint(Y, mesh.B, mesh.rest_area));
        };
        t0 = clock::now();
        step = backtracking_line_search(merit, merit0, directional, P.max_line_search);
        res.timing.line_search += since(t0);
        if (step > 0.0) {
          ok = true;
          damping = step == 1.0 ? std::max(damping / 3.0, 1e-12) : damping;
        } else {
          damping *= 10.0;
        }
      }
      if (!ok) {
        StallEvent ev;
        ev.stage = stage;
        ev.iteration = it;
        ev.slim_triangles = count_slim_triangles(X, mesh.B);
        ev.gradient_norm = g.norm();
        res.stalls.push_back(ev);
        say(fmt("  barrier stall at stage %d iteration %d: %lld slim triangles (<5 deg), |g| = %.3g",
                stage, it, static_cast<long long>(ev.slim_triangles), ev.gradient_norm));
        sl.end_reason = "stall";
        break;
      }
      X += step * unflatten(dx);
      ++res.iterations;
      if (P.log_iterations)
        say(fmt("    iter %3d  f %.6g  d %.6g  step %.3g  |dx| %.3g  damping %.3g  cg %d  lambda %.3g",
                it, f0, polycube_error(X, mesh.B), step, dx.cwiseAbs().maxCoeff(), damping,
                kkt.last_iterations(), lambda));

      const double d = polycube_error(X, mesh.B);
      changes.push_back(std::abs(d - d_prev));
      d_prev = d;
      if (static_cast<int>(changes.size()) >= P.window) {
        double s = 0.0;
        for (size_t k = changes.size() - static_cast<size_t>(P.window); k < changes.size(); ++k)
          s += changes[k];
        if (s < P.window_ratio * d0) {
          sl.end_reason = "error plateau";
          ++it;
          break;
        }
      }
      if ((step * dx).cwiseAbs().maxCoeff() < 1e-12) {
        sl.end_reason = "step converged";
        ++it;
        break;
      }
    }

    sl.iterations = it;
    sl.error = polycube_error(X, mesh.B);
    sl.arap = arap.value(X) / arap.weight;
    sl.seconds = std::chrono::duration<double>(clock::now() - t_stage).count();
    res.stage_logs.push_back(sl);
    res.stages = stage + 1;
    say(fmt("stage %2d  alpha %-8.4g eps %-7.4g error %-10.4g arap %-10.4g iters %3d  %5.1fs  "
            "cg %lld (max %d)  (%s)",
            stage, alpha, eps, sl.error, sl.arap, sl.iterations, sl.seconds, sl.cg_iterations,
            sl.cg_max, sl.end_reason.c_str()));
    if (sl.error < P.error_threshold) {
      res.reached_threshold = true;
      break;
    }
    alpha *= P.alpha_factor;
    if (eps > P.eps_min) eps = std::max(eps * P.eps_factor, P.eps_min);
  }

  res.X = X;
  res.error = polycube_error(X, mesh.B);
  res.arap_per_tet = arap.per_tet(X);
  res.arap_mean = res.arap_per_tet.dot(mesh.grads.rest_volume) / mesh.rest_volume;
  res.arap_max = res.arap_per_tet.size() ? res.arap_per_tet.maxCoeff() : 0.0;
  res.area_rel_diff = std::abs(area_constraint(X, mesh.B, mesh.rest_area)) / mesh.rest_area;
  res.min_normal_slack = barrier.min_slack(X);
  res.inverted_tets = tet_barrier.inverted(X);
  res.seconds = std::chrono::duration<double>(clock::now() - t_start).count();
  if (!res.reached_threshold)
    say("warning: polycube error " + std::to_string(res.error) + " did not reach the threshold " +
        std::to_string(P.error_threshold));
  return res;
}

}  // namespace pcp
