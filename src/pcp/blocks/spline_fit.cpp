#include "pcp/blocks/spline_fit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
#include <igl/point_mesh_squared_distance.h>

#include "pcp/structure/locator.h"

namespace pcp {

nlohmann::json FitParams::to_json() const {
  return {{"degree", degree},
          {"element_size", element_size},
          {"samples_per_element", samples_per_element},
          {"smoothing", smoothing}};
}

void FitParams::update_from_json(const nlohmann::json& j) {
  const nlohmann::json known = to_json();
  for (auto it = j.begin(); it != j.end(); ++it)
    if (!known.contains(it.key())) throw Error("unknown fit parameter '" + it.key() + "'");
  if (j.contains("degree")) degree = j.at("degree").get<int>();
  if (j.contains("element_size")) element_size = j.at("element_size").get<double>();
  if (j.contains("samples_per_element")) samples_per_element = j.at("samples_per_element").get<int>();
  if (j.contains("smoothing")) smoothing = j.at("smoothing").get<double>();
  if (degree < 1 || degree > 10) throw Error("fit.degree must be between 1 and 10, got " + std::to_string(degree));
  if (!(element_size > 0.0)) throw Error("fit.element_size must be positive");
  if (samples_per_element < 1) throw Error("fit.samples_per_element must be positive");
  if (!(smoothing >= 0.0)) throw Error("fit.smoothing must be non-negative");
}

int elements_for_length(double length, double h) {
  return std::max(1, static_cast<int>(std::lround(length / h)));
}

std::array<KnotVector, 3> SplineSpace::block_basis(const Block& b) const {
  std::array<KnotVector, 3> kv;
  for (int d = 0; d < 3; ++d)
    kv[static_cast<size_t>(d)] = open_uniform_knots(degree, elements[static_cast<size_t>(d)][static_cast<size_t>(b.cell[static_cast<size_t>(d)])]);
  return kv;
}

SplineSpace build_spline_space(const BlockLayout& L, int degree, double h) {
  SplineSpace S;
  S.degree = degree;
  for (int a = 0; a < 3; ++a) {
    auto& el = S.elements[static_cast<size_t>(a)];
    auto& off = S.offset[static_cast<size_t>(a)];
    const auto& cu = L.cuts[static_cast<size_t>(a)];
    off.push_back(0);
    for (int k = 0; k < L.intervals(a); ++k) {
      el.push_back(elements_for_length(cu[static_cast<size_t>(k) + 1] - cu[static_cast<size_t>(k)], h));
      off.push_back(off.back() + el.back() + degree - 1);
    }
    S.lattice[static_cast<size_t>(a)] = off.back() + 1;
    off.pop_back();
  }
  const long long lx = S.lattice[0], ly = S.lattice[1], lz = S.lattice[2];
  std::vector<int> unknown(static_cast<size_t>(lx * ly * lz), -1);
  auto lattice_index = [&](const Block& b, int i, int j, int k) {
    const long long gi = S.offset[0][static_cast<size_t>(b.cell[0])] + i;
    const long long gj = S.offset[1][static_cast<size_t>(b.cell[1])] + j;
    const long long gk = S.offset[2][static_cast<size_t>(b.cell[2])] + k;
    return static_cast<size_t>(gi + lx * (gj + ly * gk));
  };
  auto for_each_local = [&](const Block& b, auto&& fn) {
    const int nu = S.elements[0][static_cast<size_t>(b.cell[0])] + degree;
    const int nv = S.elements[1][static_cast<size_t>(b.cell[1])] + degree;
    const int nw = S.elements[2][static_cast<size_t>(b.cell[2])] + degree;
    for (int k = 0; k < nw; ++k)
      for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) fn(lattice_index(b, i, j, k));
  };
  for (const Block& b : L.blocks) for_each_local(b, [&](size_t g) { unknown[g] = 0; });
  for (int& u : unknown)
    if (u == 0) u = S.n_unknowns++;
  for (const Block& b : L.blocks) {
    std::vector<int> dofs;
    for_each_local(b, [&](size_t g) { dofs.push_back(unknown[g]); });
    S.block_dofs.push_back(std::move(dofs));
  }
  return S;
}

namespace {

double diagonal(const MatX3d& X) { return (X.colwise().maxCoeff() - X.colwise().minCoeff()).norm(); }

// Image of q under the map using tet t (barycentric coordinates may be outside [0, 1]).
Eigen::Vector3d map_with_tet(const MatX3d& X, const MatX4i& T, const MatX3d& Y, Index t, const Eigen::Vector3d& q) {
  const Eigen::Vector3d x0 = X.row(T(t, 0)).transpose();
  Eigen::Matrix3d D;
  for (int k = 0; k < 3; ++k) D.col(k) = X.row(T(t, k + 1)).transpose() - x0;
  const Eigen::Vector3d l = D.fullPivLu().solve(q - x0);
  Eigen::Vector3d p = (1.0 - l.sum()) * Y.row(T(t, 0)).transpose();
  for (int k = 0; k < 3; ++k) p += l(k) * Y.row(T(t, k + 1)).transpose();
  return p;
}

}  // namespace

FitSamples sample_blocks(const BlockLayout& L, const SplineSpace& S, int s, const MatX3d& X_poly,
                         const MatX4i& T, const MatX3d& X_target) {
  FitSamples out;
  TetLocator loc(X_poly, T);
  const double eps = 1e-10 * diagonal(X_poly);
  MatX3d centroids;  // built on demand for the nearest-tet fallback
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    const Block& b = L.blocks[bi];
    const Eigen::Vector3d center = 0.5 * (b.lo + b.hi);
    int m[3];
    for (int d = 0; d < 3; ++d) m[d] = S.elements[static_cast<size_t>(d)][static_cast<size_t>(b.cell[static_cast<size_t>(d)])] * s + 1;
    for (int k = 0; k < m[2]; ++k)
      for (int j = 0; j < m[1]; ++j)
        for (int i = 0; i < m[0]; ++i) {
          const int idx[3] = {i, j, k};
          Eigen::Vector3d t;
          bool surface = false;
          for (int d = 0; d < 3; ++d) {
            t(d) = static_cast<double>(idx[d]) / (m[d] - 1);
            if (idx[d] == 0 && !b.faces[static_cast<size_t>(2 * d)].interior()) surface = true;
            if (idx[d] == m[d] - 1 && !b.faces[static_cast<size_t>(2 * d + 1)].interior()) surface = true;
          }
          Eigen::Vector3d q = b.lo + t.cwiseProduct(b.hi - b.lo);
          std::optional<Eigen::Vector3d> y = loc.map(q, X_target);
          if (!y) {  // on the boundary up to round-off: move slightly towards the block center
            const Eigen::Vector3d dc = center - q;
            if (dc.norm() > eps) q += eps * dc.normalized();
            y = loc.map(q, X_target);
          }
          if (!y) {
            if (centroids.rows() == 0) {
              centroids.resize(T.rows(), 3);
              for (Index e = 0; e < T.rows(); ++e)
                centroids.row(e) = 0.25 * (X_poly.row(T(e, 0)) + X_poly.row(T(e, 1)) + X_poly.row(T(e, 2)) + X_poly.row(T(e, 3)));
            }
            Index best = 0;
            (centroids.rowwise() - q.transpose()).rowwise().squaredNorm().minCoeff(&best);
            y = map_with_tet(X_poly, T, X_target, best, q);
            ++out.unlocated;
          }
          out.block.push_back(static_cast<int>(bi));
          out.param.push_back(t);
          out.target.push_back(*y);
          out.surface.push_back(surface ? 1 : 0);
        }
  }
  return out;
}

namespace {

// Banded local accumulator of a block: entry (l, l + du + nu * (dv + nv * dw)) with |d*| <= q.
struct Band {
  int nu = 0, nv = 0, nw = 0, q = 0, width = 0;
  std::vector<double> a;
  Band(int nu_, int nv_, int nw_, int q_) : nu(nu_), nv(nv_), nw(nw_), q(q_), width(2 * q_ + 1) {
    a.assign(static_cast<size_t>(nu) * nv * nw * width * width * width, 0.0);
  }
  double& at(int i, int j, int k, int di, int dj, int dk) {
    const size_t l = static_cast<size_t>(i + nu * (j + nv * k));
    const size_t r = static_cast<size_t>((di + q) + width * ((dj + q) + width * (dk + q)));
    return a[l * static_cast<size_t>(width * width * width) + r];
  }
  void emit(const std::vector<int>& dofs, std::vector<Eigen::Triplet<double>>& out) const {
    const int W = width * width * width;
    for (int k = 0; k < nw; ++k)
      for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
          const size_t l = static_cast<size_t>(i + nu * (j + nv * k));
          for (int r = 0; r < W; ++r) {
            const double v = a[l * static_cast<size_t>(W) + static_cast<size_t>(r)];
            if (v == 0.0) continue;
            const int di = r % width - q, dj = (r / width) % width - q, dk = r / (width * width) - q;
            const int l2 = (i + di) + nu * ((j + dj) + nv * (k + dk));
            out.emplace_back(dofs[l], dofs[static_cast<size_t>(l2)], v);
          }
        }
  }
};

}  // namespace

FitSolution fit_volumes(const BlockLayout& L, const SplineSpace& S, const FitSamples& samples, double smoothing) {
  const int p = S.degree, n = S.n_unknowns;
  const int q = std::max(p, 2);
  std::vector<Eigen::Triplet<double>> ta, ts;
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(n, 3);
  size_t si = 0;
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    const Block& b = L.blocks[bi];
    const auto kv = S.block_basis(b);
    const int nu = kv[0].n_basis(), nv = kv[1].n_basis(), nw = kv[2].n_basis();
    const std::vector<int>& dofs = S.block_dofs[bi];
    Band A(nu, nv, nw, q), Sm(nu, nv, nw, q);
    double N[3][16];
    int span[3];
    for (; si < samples.block.size() && samples.block[si] == static_cast<int>(bi); ++si) {
      for (int d = 0; d < 3; ++d) {
        span[d] = find_span(kv[static_cast<size_t>(d)], samples.param[si](d));
        basis_functions(kv[static_cast<size_t>(d)], span[d], samples.param[si](d), N[d], nullptr);
      }
      const int i0 = span[0] - p, j0 = span[1] - p, k0 = span[2] - p;
      for (int c = 0; c <= p; ++c)
        for (int bb = 0; bb <= p; ++bb)
          for (int a = 0; a <= p; ++a) {
            const double wi = N[0][a] * N[1][bb] * N[2][c];
            if (wi == 0.0) continue;
            const int li = (i0 + a) + nu * ((j0 + bb) + nv * (k0 + c));
            rhs.row(dofs[static_cast<size_t>(li)]) += wi * samples.target[si].transpose();
            for (int c2 = 0; c2 <= p; ++c2)
              for (int b2 = 0; b2 <= p; ++b2)
                for (int a2 = 0; a2 <= p; ++a2) {
                  const double wj = N[0][a2] * N[1][b2] * N[2][c2];
                  if (wj != 0.0) A.at(i0 + a, j0 + bb, k0 + c, a2 - a, b2 - bb, c2 - c) += wi * wj;
                }
          }
    }
    // Second differences of the control net along each direction.
    const int dims[3] = {nu, nv, nw};
    for (int d = 0; d < 3; ++d) {
      if (dims[d] < 3) continue;
      for (int k = 0; k < nw; ++k)
        for (int j = 0; j < nv; ++j)
          for (int i = 0; i < nu; ++i) {
            int c[3] = {i, j, k};
            if (c[d] < 1 || c[d] > dims[d] - 2) continue;
            const double coef[3] = {1.0, -2.0, 1.0};
            for (int r = 0; r < 3; ++r)
              for (int s2 = 0; s2 < 3; ++s2) {
                int x[3] = {i, j, k}, off[3] = {0, 0, 0};
                x[d] += r - 1;
                off[d] = s2 - r;
                Sm.at(x[0], x[1], x[2], off[0], off[1], off[2]) += coef[r] * coef[s2];
              }
          }
    }
    A.emit(dofs, ta);
    Sm.emit(dofs, ts);
  }
  if (si != samples.block.size()) throw Error("fit_volumes: samples are not grouped by block");

  Eigen::SparseMatrix<double> AtA(n, n), Ssp(n, n);
  AtA.setFromTriplets(ta.begin(), ta.end());
  Ssp.setFromTriplets(ts.begin(), ts.end());
  double trA = 0.0, trS = 0.0;
  for (int i = 0; i < n; ++i) trA += AtA.coeff(i, i), trS += Ssp.coeff(i, i);
  const double kappa = trS > 0.0 ? trA / trS : 0.0;
  Eigen::SparseMatrix<double> M = AtA + (smoothing * kappa) * Ssp;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(M);
  if (ldlt.info() != Eigen::Success)
    throw Error("B-spline fitting: the least-squares system could not be factorized; try a larger fit.smoothing");
  Eigen::MatrixXd c = ldlt.solve(rhs);
  if (ldlt.info() != Eigen::Success || !c.allFinite() || (ldlt.vectorD().array() <= 0.0).any())
    throw Error("B-spline fitting: the least-squares system is singular; try a larger fit.smoothing");

  FitSolution sol;
  sol.controls = c;
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    TensorBSpline3 V;
    V.basis = S.block_basis(L.blocks[bi]);
    for (int dof : S.block_dofs[bi]) V.coefs.push_back(c.row(dof).transpose());
    sol.volumes.push_back(std::move(V));
  }
  return sol;
}

std::array<int, 3> eval_grid_size(const SplineSpace& S, const Block& b) {
  std::array<int, 3> m{};
  for (int d = 0; d < 3; ++d)
    m[static_cast<size_t>(d)] = S.elements[static_cast<size_t>(d)][static_cast<size_t>(b.cell[static_cast<size_t>(d)])] * 4 + 1;
  return m;
}

double scaled_jacobian(const Eigen::Matrix3d& J) {
  const double l = J.col(0).norm() * J.col(1).norm() * J.col(2).norm();
  return l > 0.0 ? J.determinant() / l : 0.0;
}

SampleMesh sample_mesh(const BlockLayout& L, const SplineSpace& S, const FitSolution& sol) {
  SampleMesh M;
  Index np = 0, nc = 0;
  for (const Block& b : L.blocks) {
    const auto m = eval_grid_size(S, b);
    np += static_cast<Index>(m[0]) * m[1] * m[2];
    nc += static_cast<Index>(m[0] - 1) * (m[1] - 1) * (m[2] - 1);
  }
  M.X.resize(np, 3);
  M.scaled_jacobian.resize(np);
  M.H.resize(nc, 8);
  M.block.resize(nc);
  Index p0 = 0, c = 0;
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    const auto m = eval_grid_size(S, L.blocks[bi]);
    auto pid = [&](int i, int j, int k) { return static_cast<int>(p0 + i + m[0] * (j + m[1] * k)); };
    for (int k = 0; k < m[2]; ++k)
      for (int j = 0; j < m[1]; ++j)
        for (int i = 0; i < m[0]; ++i) {
          Eigen::Vector3d x;
          Eigen::Matrix3d J;
          sol.volumes[bi].eval_jacobian(double(i) / (m[0] - 1), double(j) / (m[1] - 1), double(k) / (m[2] - 1), x, J);
          M.X.row(pid(i, j, k)) = x.transpose();
          M.scaled_jacobian(pid(i, j, k)) = scaled_jacobian(J);
        }
    for (int k = 0; k + 1 < m[2]; ++k)
      for (int j = 0; j + 1 < m[1]; ++j)
        for (int i = 0; i + 1 < m[0]; ++i, ++c) {
          M.H.row(c) << pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k), pid(i, j + 1, k),
              pid(i, j, k + 1), pid(i + 1, j, k + 1), pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1);
          M.block(c) = static_cast<int>(bi);
        }
    p0 += static_cast<Index>(m[0]) * m[1] * m[2];
  }
  return M;
}

nlohmann::json FitQuality::to_json() const {
  return {{"rms_error", rms_error},
          {"max_error", max_error},
          {"surface_rms_error", surface_rms_error},
          {"surface_max_error", surface_max_error},
          {"surface_distance_max", surface_distance_max},
          {"surface_distance_mean", surface_distance_mean},
          {"min_scaled_jacobian", min_scaled_jacobian},
          {"mean_scaled_jacobian", mean_scaled_jacobian},
          {"eval_points", eval_points},
          {"nonpositive_jacobian", nonpositive_jacobian},
          {"nonpositive_fraction", eval_points > 0 ? static_cast<double>(nonpositive_jacobian) / eval_points : 0.0},
          {"blocks_with_nonpositive", blocks_with_nonpositive},
          {"c0_max_gap", c0_max_gap}};
}

FitQuality evaluate_fit(const BlockLayout& L, const SplineSpace& S, const FitSolution& sol,
                        const FitSamples& samples, const MatX3d& surface_V, const MatX3i& surface_F, double diag) {
  FitQuality Q;
  // Errors at the fitting samples.
  double sum2 = 0.0, ssum2 = 0.0;
  long long ns = 0;
  for (size_t i = 0; i < samples.block.size(); ++i) {
    const Eigen::Vector3d& t = samples.param[i];
    const double e = (sol.volumes[static_cast<size_t>(samples.block[i])].eval(t(0), t(1), t(2)) - samples.target[i]).norm();
    sum2 += e * e;
    Q.max_error = std::max(Q.max_error, e);
    if (samples.surface[i]) {
      ssum2 += e * e;
      Q.surface_max_error = std::max(Q.surface_max_error, e);
      ++ns;
    }
  }
  if (!samples.block.empty()) Q.rms_error = std::sqrt(sum2 / static_cast<double>(samples.block.size())) / diag;
  if (ns > 0) Q.surface_rms_error = std::sqrt(ssum2 / static_cast<double>(ns)) / diag;
  Q.max_error /= diag;
  Q.surface_max_error /= diag;

  // Scaled Jacobian on the evaluation grids; boundary face points for the surface distance.
  std::vector<Eigen::Vector3d> bpts;
  double sj_sum = 0.0;
  Q.min_scaled_jacobian = std::numeric_limits<double>::infinity();
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    const Block& b = L.blocks[bi];
    const auto m = eval_grid_size(S, b);
    bool bad = false;
    for (int k = 0; k < m[2]; ++k)
      for (int j = 0; j < m[1]; ++j)
        for (int i = 0; i < m[0]; ++i) {
          const int idx[3] = {i, j, k};
          double t[3];
          bool on_surface = false;
          for (int d = 0; d < 3; ++d) {
            t[d] = static_cast<double>(idx[d]) / (m[static_cast<size_t>(d)] - 1);
            if (idx[d] == 0 && !b.faces[static_cast<size_t>(2 * d)].interior()) on_surface = true;
            if (idx[d] == m[static_cast<size_t>(d)] - 1 && !b.faces[static_cast<size_t>(2 * d + 1)].interior()) on_surface = true;
          }
          Eigen::Vector3d x;
          Eigen::Matrix3d J;
          sol.volumes[bi].eval_jacobian(t[0], t[1], t[2], x, J);
          const double sj = scaled_jacobian(J);
          sj_sum += sj;
          Q.min_scaled_jacobian = std::min(Q.min_scaled_jacobian, sj);
          ++Q.eval_points;
          if (sj <= 0.0) {
            ++Q.nonpositive_jacobian;
            bad = true;
          }
          if (on_surface) bpts.push_back(x);
        }
    if (bad) ++Q.blocks_with_nonpositive;
  }
  if (Q.eval_points > 0) Q.mean_scaled_jacobian = sj_sum / static_cast<double>(Q.eval_points);
  else Q.min_scaled_jacobian = 0.0;

  if (!bpts.empty()) {
    Eigen::MatrixXd P(static_cast<Index>(bpts.size()), 3);
    for (size_t i = 0; i < bpts.size(); ++i) P.row(static_cast<Index>(i)) = bpts[i].transpose();
    const Eigen::MatrixXd V = surface_V;
    const Eigen::MatrixXi F = surface_F;
    Eigen::VectorXd d2;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::point_mesh_squared_distance(P, V, F, d2, I, C);
    const Eigen::ArrayXd d = d2.array().sqrt();
    Q.surface_distance_max = d.maxCoeff() / diag;
    Q.surface_distance_mean = d.mean() / diag;
  }

  // C0 gap across interior faces (evaluated from both sides).
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    const Block& b = L.blocks[bi];
    const auto m = eval_grid_size(S, b);
    for (int f = 0; f < 6; ++f) {
      const int nb = b.faces[static_cast<size_t>(f)].neighbor;
      if (nb <= static_cast<int>(bi)) continue;
      const int a = f / 2, d1 = (a + 1) % 3, d2 = (a + 2) % 3;
      for (int r = 0; r < m[static_cast<size_t>(d1)]; ++r)
        for (int c = 0; c < m[static_cast<size_t>(d2)]; ++c) {
          Eigen::Vector3d t0, t1;
          t0(a) = (f % 2 == 1) ? 1.0 : 0.0;
          t1(a) = 1.0 - t0(a);
          t0(d1) = t1(d1) = static_cast<double>(r) / (m[static_cast<size_t>(d1)] - 1);
          t0(d2) = t1(d2) = static_cast<double>(c) / (m[static_cast<size_t>(d2)] - 1);
          const Eigen::Vector3d x0 = sol.volumes[bi].eval(t0(0), t0(1), t0(2));
          const Eigen::Vector3d x1 = sol.volumes[static_cast<size_t>(nb)].eval(t1(0), t1(1), t1(2));
          Q.c0_max_gap = std::max(Q.c0_max_gap, (x0 - x1).norm() / diag);
        }
    }
  }
  return Q;
}

}  // namespace pcp
