#include "pcp/structure/final_map.h"

#include <chrono>
#include <cmath>
#include <array>
#include <numeric>

#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

#include "pcp/tet/tetrahedralize.h"

namespace pcp {

nlohmann::json FinalMapParams::to_json() const {
  return {{"heights_fixed", heights_fixed},
          {"max_iterations", max_iterations},
          {"rel_tolerance", rel_tolerance},
          {"laplacian_init", laplacian_init},
          {"untangle", untangle},
          {"untangle_budget", untangle_budget}};
}

void FinalMapParams::update_from_json(const nlohmann::json& j) {
  const nlohmann::json known = to_json();
  for (auto it = j.begin(); it != j.end(); ++it)
    if (!known.contains(it.key())) throw Error("unknown final map parameter '" + it.key() + "'");
  if (j.contains("heights_fixed")) heights_fixed = j.at("heights_fixed").get<bool>();
  if (j.contains("max_iterations")) max_iterations = j.at("max_iterations").get<int>();
  if (j.contains("rel_tolerance")) rel_tolerance = j.at("rel_tolerance").get<double>();
  if (j.contains("laplacian_init")) laplacian_init = j.at("laplacian_init").get<bool>();
  if (j.contains("untangle")) untangle = j.at("untangle").get<bool>();
  if (j.contains("untangle_budget")) untangle_budget = j.at("untangle_budget").get<double>();
}

namespace {

struct UnionFind {
  std::vector<int> p;
  explicit UnionFind(size_t n) : p(n) { std::iota(p.begin(), p.end(), 0); }
  int find(int x) {
    while (p[static_cast<size_t>(x)] != x) x = p[static_cast<size_t>(x)] = p[static_cast<size_t>(p[static_cast<size_t>(x)])];
    return x;
  }
  void unite(int a, int b) { p[static_cast<size_t>(find(a))] = find(b); }
};

using SpMat = Eigen::SparseMatrix<double>;

// Variable layout of one coordinate axis: each node maps either to a free variable or to a
// height group (patches of this axis glued by shared vertices).
struct AxisLayout {
  std::vector<int> node_var;     // variable index of each node
  int n_vars = 0;
  int n_free = 0;                // free node variables come first: [0, n_free)
  std::vector<int> patch_group;  // per patch of this axis: group id (var = n_free + group), else -1
  int n_groups = 0;
};

AxisLayout make_layout(int axis, Index n_nodes, const MatX3i& B, const PatchGraph& g) {
  const size_t np = g.patches.label.size();
  AxisLayout L;
  L.patch_group.assign(np, -1);
  UnionFind uf(np);
  std::vector<int> node_patch(static_cast<size_t>(n_nodes), -1);
  for (size_t p = 0; p < np; ++p) {
    if (label_axis(g.patches.label[p]) != axis) continue;
    for (int f : g.patches.faces[p])
      for (int k = 0; k < 3; ++k) {
        int v = B(f, k);
        int& np_v = node_patch[static_cast<size_t>(v)];
        if (np_v < 0)
          np_v = static_cast<int>(p);
        else if (np_v != static_cast<int>(p))
          uf.unite(np_v, static_cast<int>(p));
      }
  }
  std::vector<int> root_group(np, -1);
  for (size_t p = 0; p < np; ++p) {
    if (label_axis(g.patches.label[p]) != axis) continue;
    int r = uf.find(static_cast<int>(p));
    if (root_group[static_cast<size_t>(r)] < 0) root_group[static_cast<size_t>(r)] = L.n_groups++;
    L.patch_group[p] = root_group[static_cast<size_t>(r)];
  }
  L.node_var.assign(static_cast<size_t>(n_nodes), -1);
  for (Index i = 0; i < n_nodes; ++i)
    if (node_patch[static_cast<size_t>(i)] < 0) L.node_var[static_cast<size_t>(i)] = L.n_free++;
  for (Index i = 0; i < n_nodes; ++i) {
    int p = node_patch[static_cast<size_t>(i)];
    if (p >= 0) L.node_var[static_cast<size_t>(i)] = L.n_free + L.patch_group[static_cast<size_t>(p)];
  }
  L.n_vars = L.n_free + L.n_groups;
  return L;
}

// Generic reduced quadratic solve: minimize 0.5 x^T A x - b^T x over node values x = S y with
// some variables fixed. `fixed` marks fixed variables with values in y (in/out).
class ReducedSolver {
 public:
  void setup(const SpMat& A, const AxisLayout& L, const std::vector<char>& fixed) {
    fixed_ = fixed;
    map_.assign(static_cast<size_t>(L.n_vars), -1);
    n_unknown_ = 0;
    for (int v = 0; v < L.n_vars; ++v)
      if (!fixed_[static_cast<size_t>(v)]) map_[static_cast<size_t>(v)] = n_unknown_++;
    // S: nodes -> vars (each node has exactly one var).
    std::vector<Eigen::Triplet<double>> trip;
    for (size_t i = 0; i < L.node_var.size(); ++i) trip.emplace_back(static_cast<int>(i), L.node_var[i], 1.0);
    S_.resize(static_cast<Index>(L.node_var.size()), L.n_vars);
    S_.setFromTriplets(trip.begin(), trip.end());
    SpMat Avar = (S_.transpose() * A * S_).pruned();
    // Split into unknown/fixed blocks.
    std::vector<Eigen::Triplet<double>> tu, tf;
    for (Index c = 0; c < Avar.outerSize(); ++c)
      for (SpMat::InnerIterator it(Avar, c); it; ++it) {
        int r = static_cast<int>(it.row()), cc = static_cast<int>(it.col());
        int ur = map_[static_cast<size_t>(r)];
        if (ur < 0) continue;
        int uc = map_[static_cast<size_t>(cc)];
        if (uc >= 0)
          tu.emplace_back(ur, uc, it.value());
        else
          tf.emplace_back(ur, cc, it.value());
      }
    Auu_.resize(n_unknown_, n_unknown_);
    Auu_.setFromTriplets(tu.begin(), tu.end());
    Auf_.resize(n_unknown_, L.n_vars);
    Auf_.setFromTriplets(tf.begin(), tf.end());
    ldlt_.compute(Auu_);
    if (ldlt_.info() != Eigen::Success) throw Error("final map: reduced system is singular");
  }

  // b: node-space right-hand side. y: variable values (fixed entries are used, others written).
  void solve(const Eigen::VectorXd& b, Eigen::VectorXd& y) const {
    Eigen::VectorXd bvar = S_.transpose() * b;
    Eigen::VectorXd yfixed = y;
    for (size_t v = 0; v < fixed_.size(); ++v)
      if (!fixed_[v]) yfixed(static_cast<Index>(v)) = 0.0;
    Eigen::VectorXd rhs(n_unknown_);
    for (size_t v = 0; v < map_.size(); ++v)
      if (map_[v] >= 0) rhs(map_[v]) = bvar(static_cast<Index>(v));
    rhs -= Auf_ * yfixed;
    Eigen::VectorXd yu = ldlt_.solve(rhs);
    for (size_t v = 0; v < map_.size(); ++v)
      if (map_[v] >= 0) y(static_cast<Index>(v)) = yu(map_[v]);
  }

  Eigen::VectorXd nodes(const Eigen::VectorXd& y) const { return S_ * y; }

 private:
  std::vector<char> fixed_;
  std::vector<int> map_;
  int n_unknown_ = 0;
  SpMat S_, Auu_, Auf_;
  Eigen::SimplicialLDLT<SpMat> ldlt_;
};

// Local repair of inverted tets: around every inverted tet (`rings`-ring of nodes) minimize
//   sum_t V_t * 0.5 ||F_t - R_t||^2 + w * sum_t V_t * (max(0, tau_t - vol_t) / tau_t)^2
// over the coordinates that the patch constraints leave free (`movable`), by gradient descent
// with Armijo backtracking. tau_t is a small fraction of the rest volume. Returns the number of
// inverted tets left.
Index untangle_inverted(const MatX3d& rest, MatX3d& X, const MatX4i& T, const TetGradients& grads,
                        const std::vector<std::array<char, 3>>& movable, int rings,
                        double& budget) {
  auto vol = [&](const MatX3d& Y, Index t) { return tet_volume(Y, T, t); };
  std::vector<Index> inverted;
  for (Index t = 0; t < T.rows(); ++t)
    if (!(vol(X, t) > 0.0)) inverted.push_back(t);
  if (inverted.empty()) return 0;

  // Node -> tets adjacency.
  std::vector<std::vector<int>> node_tets(static_cast<size_t>(X.rows()));
  for (Index t = 0; t < T.rows(); ++t)
    for (int k = 0; k < 4; ++k) node_tets[static_cast<size_t>(T(t, k))].push_back(static_cast<int>(t));
  std::vector<char> in_region(static_cast<size_t>(X.rows()), 0);
  std::vector<int> region;
  for (Index t : inverted)
    for (int k = 0; k < 4; ++k)
      if (!in_region[static_cast<size_t>(T(t, k))]) {
        in_region[static_cast<size_t>(T(t, k))] = 1;
        region.push_back(T(t, k));
      }
  for (int ring = 0; ring < rings; ++ring) {
    std::vector<int> add;
    for (int v : region)
      for (int t : node_tets[static_cast<size_t>(v)])
        for (int k = 0; k < 4; ++k)
          if (!in_region[static_cast<size_t>(T(t, k))]) {
            in_region[static_cast<size_t>(T(t, k))] = 1;
            add.push_back(T(t, k));
          }
    region.insert(region.end(), add.begin(), add.end());
  }
  std::vector<char> tet_in(static_cast<size_t>(T.rows()), 0);
  std::vector<int> tets;
  for (int v : region)
    for (int t : node_tets[static_cast<size_t>(v)])
      if (!tet_in[static_cast<size_t>(t)]) {
        tet_in[static_cast<size_t>(t)] = 1;
        tets.push_back(t);
      }

  const double w = 1e3;
  const double cost = static_cast<double>(tets.size());  // tet evaluations per energy/gradient call
  auto energy = [&](const MatX3d& Y) {
    budget -= cost;
    double e = 0.0;
    for (int t : tets) {
      const double V = grads.rest_volume(t);
      Eigen::Matrix3d F = deformation_gradient(grads, Y, T, t);
      e += V * 0.5 * (F - closest_rotation(F)).squaredNorm();
      const double tau = 0.05 * V, v = vol(Y, t);
      if (v < tau) e += w * V * ((tau - v) / tau) * ((tau - v) / tau);
    }
    return e;
  };
  auto gradient = [&](const MatX3d& Y, MatX3d& G) {
    budget -= cost;
    G.setZero(Y.rows(), 3);
    for (int t : tets) {
      const double V = grads.rest_volume(t);
      const auto& Gt = grads.G[static_cast<size_t>(t)];
      Eigen::Matrix3d F = deformation_gradient(grads, Y, T, t);
      Eigen::Matrix3d D = F - closest_rotation(F);
      for (int i = 0; i < 4; ++i) G.row(T(t, i)) += V * (D * Gt.col(i)).transpose();
      const double tau = 0.05 * V, v = vol(Y, t);
      if (v < tau) {
        // d vol / d x_i = (1/6) * outward normal of the opposite face (area-weighted).
        static const int opp[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        const double c = -2.0 * w * V * (tau - v) / (tau * tau);
        for (int i = 0; i < 4; ++i) {
          Eigen::Vector3d a = Y.row(T(t, opp[i][0])), b = Y.row(T(t, opp[i][1])),
                          d = Y.row(T(t, opp[i][2]));
          Eigen::Vector3d dv = (b - a).cross(d - a) / 6.0;
          // Sign: volume increases when x_i moves away from the opposite face, i.e. along -n.
          G.row(T(t, i)) += c * (-dv).transpose();
        }
      }
    }
    for (Index i = 0; i < Y.rows(); ++i)
      for (int k = 0; k < 3; ++k)
        if (!in_region[static_cast<size_t>(i)] || !movable[static_cast<size_t>(i)][static_cast<size_t>(k)])
          G(i, k) = 0.0;
  };

  double e = energy(X);
  MatX3d G;
  double step = 1e-3;
  for (int it = 0; it < 5000 && budget > 0.0; ++it) {
    gradient(X, G);
    const double gn2 = G.squaredNorm();
    if (gn2 < 1e-30) break;
    bool ok = false;
    for (int ls = 0; ls < 40 && budget > 0.0; ++ls) {
      MatX3d Y = X - step * G;
      double en = energy(Y);
      if (en <= e - 1e-4 * step * gn2) {
        X = Y;
        e = en;
        step *= 2.0;
        ok = true;
        break;
      }
      step *= 0.5;
    }
    if (!ok) break;
    bool any = false;
    for (Index t : inverted) any = any || !(vol(X, t) > 0.0);
    if (!any && it % 50 == 49) {
      // Recheck the whole region before stopping.
      bool region_ok = true;
      for (int t : tets) region_ok = region_ok && vol(X, t) > 0.0;
      if (region_ok) break;
    }
  }
  (void)rest;
  Index left = 0;
  for (Index t = 0; t < T.rows(); ++t)
    if (!(vol(X, t) > 0.0)) ++left;
  return left;
}

}  // namespace

FinalMapResult compute_final_map(const MatX3d& rest, const MatX3d& X0, const MatX4i& T,
                                 const MatX3i& B, const PatchGraph& g,
                                 const std::vector<double>& heights, const FinalMapParams& P) {
  const auto t_start = std::chrono::steady_clock::now();
  const Index n = rest.rows();
  const size_t np = g.patches.label.size();
  if (heights.size() != np) throw Error("final map: heights must have one entry per patch");

  TetGradients grads = compute_tet_gradients(rest, T);
  const double total_vol = grads.rest_volume.sum();

  // Cotangent-like ARAP Laplacian L = sum_t V_t G_t^T G_t (same for every axis).
  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve(static_cast<size_t>(T.rows()) * 16);
  for (Index t = 0; t < T.rows(); ++t) {
    const auto& G = grads.G[static_cast<size_t>(t)];
    const double w = grads.rest_volume(t);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) trip.emplace_back(T(t, i), T(t, j), w * G.col(i).dot(G.col(j)));
  }
  SpMat L(n, n);
  L.setFromTriplets(trip.begin(), trip.end());

  // Umbrella operator for the eq. (17) initialization.
  SpMat U;
  if (P.laplacian_init) {
    std::vector<char> on_boundary(static_cast<size_t>(n), 0);
    for (Index f = 0; f < B.rows(); ++f)
      for (int k = 0; k < 3; ++k) on_boundary[static_cast<size_t>(B(f, k))] = 1;
    std::vector<std::vector<int>> nb(static_cast<size_t>(n)), nbb(static_cast<size_t>(n));
    for (Index t = 0; t < T.rows(); ++t)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          if (i != j) nb[static_cast<size_t>(T(t, i))].push_back(T(t, j));
    for (Index f = 0; f < B.rows(); ++f)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (i != j) nbb[static_cast<size_t>(B(f, i))].push_back(B(f, j));
    Index nbd = 0;
    for (char c : on_boundary) nbd += c;
    const double wb = 1.0 / std::max<Index>(nbd, 1), wi = 1.0 / std::max<Index>(n - nbd, 1);
    std::vector<Eigen::Triplet<double>> tu;
    for (Index i = 0; i < n; ++i) {
      auto& list = on_boundary[static_cast<size_t>(i)] ? nbb[static_cast<size_t>(i)] : nb[static_cast<size_t>(i)];
      std::sort(list.begin(), list.end());
      list.erase(std::unique(list.begin(), list.end()), list.end());
      const double w = std::sqrt(on_boundary[static_cast<size_t>(i)] ? wb : wi);
      tu.emplace_back(static_cast<int>(i), static_cast<int>(i), w);
      for (int j : list) tu.emplace_back(static_cast<int>(i), j, -w / static_cast<double>(list.size()));
    }
    U.resize(n, n);
    U.setFromTriplets(tu.begin(), tu.end());
  }

  FinalMapResult res;
  MatX3d X = X0;
  std::vector<AxisLayout> layout(3);
  std::vector<ReducedSolver> arap_solver(3);
  std::vector<Eigen::VectorXd> y(3);
  for (int k = 0; k < 3; ++k) {
    AxisLayout& Lk = layout[static_cast<size_t>(k)];
    Lk = make_layout(k, n, B, g);
    // Initial variable values: free nodes from X0, groups from the given heights (area-weighted
    // mean over glued patches).
    y[static_cast<size_t>(k)] = Eigen::VectorXd::Zero(Lk.n_vars);
    for (Index i = 0; i < n; ++i) {
      int v = Lk.node_var[static_cast<size_t>(i)];
      if (v < Lk.n_free) y[static_cast<size_t>(k)](v) = X0(i, k);
    }
    std::vector<double> gsum(static_cast<size_t>(Lk.n_groups), 0.0), garea(static_cast<size_t>(Lk.n_groups), 0.0);
    for (size_t p = 0; p < np; ++p) {
      int gr = Lk.patch_group[p];
      if (gr < 0) continue;
      gsum[static_cast<size_t>(gr)] += heights[p] * g.area[p];
      garea[static_cast<size_t>(gr)] += g.area[p];
    }
    for (int gr = 0; gr < Lk.n_groups; ++gr)
      y[static_cast<size_t>(k)](Lk.n_free + gr) =
          garea[static_cast<size_t>(gr)] > 0 ? gsum[static_cast<size_t>(gr)] / garea[static_cast<size_t>(gr)] : 0.0;

    // Eq. (17): heights fixed, free nodes relaxed.
    if (P.laplacian_init && Lk.n_groups > 0) {
      std::vector<char> fixed(static_cast<size_t>(Lk.n_vars), 0);
      for (int gr = 0; gr < Lk.n_groups; ++gr) fixed[static_cast<size_t>(Lk.n_free + gr)] = 1;
      ReducedSolver rs;
      rs.setup(SpMat(U.transpose() * U), Lk, fixed);
      rs.solve(Eigen::VectorXd::Zero(n), y[static_cast<size_t>(k)]);
    }
    // ARAP global step: heights free (one variable pinned against translation) or fixed.
    std::vector<char> fixed(static_cast<size_t>(Lk.n_vars), 0);
    if (P.heights_fixed) {
      for (int gr = 0; gr < Lk.n_groups; ++gr) fixed[static_cast<size_t>(Lk.n_free + gr)] = 1;
    }
    if (!P.heights_fixed || Lk.n_groups == 0) fixed[static_cast<size_t>(Lk.n_vars - 1)] = 1;
    arap_solver[static_cast<size_t>(k)].setup(L, Lk, fixed);
    X.col(k) = arap_solver[static_cast<size_t>(k)].nodes(y[static_cast<size_t>(k)]);
  }

  // Local-global iterations.
  auto energy = [&](const MatX3d& Y, std::vector<Eigen::Matrix3d>* rots) {
    double e = 0.0;
    for (Index t = 0; t < T.rows(); ++t) {
      Eigen::Matrix3d F = deformation_gradient(grads, Y, T, t);
      Eigen::Matrix3d R = closest_rotation(F);
      if (rots) (*rots)[static_cast<size_t>(t)] = R;
      e += grads.rest_volume(t) * 0.5 * (F - R).squaredNorm();
    }
    return e / total_vol;
  };
  std::vector<Eigen::Matrix3d> R(static_cast<size_t>(T.rows()));
  double e_prev = energy(X, &R);
  for (int it = 0; it < P.max_iterations; ++it) {
    // b_i^k = sum_t V_t g_i . (row k of R_t)
    MatX3d b = MatX3d::Zero(n, 3);
    for (Index t = 0; t < T.rows(); ++t) {
      const auto& G = grads.G[static_cast<size_t>(t)];
      const double w = grads.rest_volume(t);
      const Eigen::Matrix3d& Rt = R[static_cast<size_t>(t)];
      for (int i = 0; i < 4; ++i) b.row(T(t, i)) += w * (Rt * G.col(i)).transpose();
    }
    for (int k = 0; k < 3; ++k) {
      arap_solver[static_cast<size_t>(k)].solve(b.col(k), y[static_cast<size_t>(k)]);
      X.col(k) = arap_solver[static_cast<size_t>(k)].nodes(y[static_cast<size_t>(k)]);
    }
    double e = energy(X, &R);
    res.iterations = it + 1;
    if (std::abs(e_prev - e) <= P.rel_tolerance * std::max(e, 1e-12)) break;
    e_prev = e;
  }

  if (P.untangle) {
    // Coordinates left free by the patch constraints (same as the free variables of each axis).
    std::vector<std::array<char, 3>> movable(static_cast<size_t>(n));
    for (Index i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k)
        movable[static_cast<size_t>(i)][static_cast<size_t>(k)] =
            layout[static_cast<size_t>(k)].node_var[static_cast<size_t>(i)] <
                    layout[static_cast<size_t>(k)].n_free
                ? 1
                : 0;
    for (Index t = 0; t < T.rows(); ++t)
      if (!(tet_volume(X, T, t) > 0.0)) ++res.inverted_before_untangle;
    // Repeat on the remaining inverted tets with growing neighborhoods while it helps.
    Index left = res.inverted_before_untangle;
    double budget = P.untangle_budget;
    for (int pass = 0; pass < 6 && left > 0 && budget > 0.0; ++pass) {
      Index now = untangle_inverted(rest, X, T, grads, movable, 2 + pass, budget);
      if (now >= left && pass > 0) break;
      left = now;
    }
  }

  res.X = X;
  res.heights.assign(np, 0.0);
  for (size_t p = 0; p < np; ++p) {
    int k = label_axis(g.patches.label[p]);
    res.heights[p] = y[static_cast<size_t>(k)](layout[static_cast<size_t>(k)].n_free +
                                               layout[static_cast<size_t>(k)].patch_group[p]);
    for (int f : g.patches.faces[p])
      for (int c = 0; c < 3; ++c)
        res.max_planarity_error =
            std::max(res.max_planarity_error, std::abs(X(B(f, c), k) - res.heights[p]));
  }
  res.arap_per_tet.resize(T.rows());
  for (Index t = 0; t < T.rows(); ++t) {
    Eigen::Matrix3d F = deformation_gradient(grads, X, T, t);
    res.arap_per_tet(t) = 0.5 * (F - closest_rotation(F)).squaredNorm();
    if (!(tet_volume(X, T, t) > 0.0)) ++res.inverted_tets;
  }
  res.arap_mean = res.arap_per_tet.dot(grads.rest_volume) / total_vol;
  res.arap_max = res.arap_per_tet.size() ? res.arap_per_tet.maxCoeff() : 0.0;
  res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  return res;
}

}  // namespace pcp
