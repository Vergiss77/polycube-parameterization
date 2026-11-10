#include "pcp/tet/tetrahedralize.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

#include <Eigen/Geometry>
#include <igl/copyleft/tetgen/tetrahedralize.h>

#include "pcp/core/constants.h"

namespace pcp {

nlohmann::json TetParams::to_json() const {
  return {{"radius_edge_ratio", radius_edge_ratio},
          {"min_dihedral_deg", min_dihedral_deg},
          {"max_volume", max_volume},
          {"size_factor", size_factor},
          {"split_boundary_edges", split_boundary_edges}};
}

void TetParams::update_from_json(const nlohmann::json& j) {
  const nlohmann::json known = to_json();
  for (auto it = j.begin(); it != j.end(); ++it)
    if (!known.contains(it.key())) throw Error("unknown tet parameter '" + it.key() + "'");
  if (j.contains("radius_edge_ratio")) radius_edge_ratio = j.at("radius_edge_ratio").get<double>();
  if (j.contains("min_dihedral_deg")) min_dihedral_deg = j.at("min_dihedral_deg").get<double>();
  if (j.contains("max_volume")) max_volume = j.at("max_volume").get<double>();
  if (j.contains("size_factor")) size_factor = j.at("size_factor").get<double>();
  if (j.contains("split_boundary_edges"))
    split_boundary_edges = j.at("split_boundary_edges").get<bool>();
}

nlohmann::json TetReport::to_json() const {
  return {{"nodes", n_nodes},
          {"tets", n_tets},
          {"boundary_faces", n_boundary},
          {"min_dihedral_deg", min_dihedral_deg},
          {"max_dihedral_deg", max_dihedral_deg},
          {"tets_below_5deg_dihedral", n_small_dihedral},
          {"volume", volume},
          {"surface_volume", surface_volume},
          {"max_volume", max_volume_used},
          {"split_edges", split_edges},
          {"switches", switches}};
}

double tet_volume(const MatX3d& X, const MatX4i& T, Index t) {
  Eigen::Vector3d a = X.row(T(t, 0)), b = X.row(T(t, 1)), c = X.row(T(t, 2)), d = X.row(T(t, 3));
  return (b - a).dot((c - a).cross(d - a)) / 6.0;
}

DihedralStats dihedral_stats(const MatX3d& X, const MatX4i& T) {
  static const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  DihedralStats st;
  for (Index t = 0; t < T.rows(); ++t) {
    double tmin = 180.0;
    for (const auto& e : edges) {
      // The two vertices not on edge e.
      int o[2], k = 0;
      for (int i = 0; i < 4; ++i)
        if (i != e[0] && i != e[1]) o[k++] = i;
      Eigen::Vector3d p0 = X.row(T(t, e[0])), p1 = X.row(T(t, e[1]));
      Eigen::Vector3d u = (p1 - p0).normalized();
      Eigen::Vector3d a = X.row(T(t, o[0])), b = X.row(T(t, o[1]));
      Eigen::Vector3d va = (a - p0) - (a - p0).dot(u) * u;
      Eigen::Vector3d vb = (b - p0) - (b - p0).dot(u) * u;
      double c = va.normalized().dot(vb.normalized());
      double ang = std::acos(std::clamp(c, -1.0, 1.0)) * kRadToDeg;
      tmin = std::min(tmin, ang);
      st.max_deg = std::max(st.max_deg, ang);
    }
    st.min_deg = std::min(st.min_deg, tmin);
    if (tmin < 5.0) ++st.n_below_5deg;
  }
  return st;
}

Index split_boundary_spanning_edges(TetMesh& m) {
  std::vector<char> on_boundary(static_cast<size_t>(m.X.rows()), 0);
  std::set<std::pair<int, int>> boundary_edges;
  for (Index f = 0; f < m.B.rows(); ++f)
    for (int k = 0; k < 3; ++k) {
      int a = m.B(f, k), b = m.B(f, (k + 1) % 3);
      on_boundary[static_cast<size_t>(a)] = 1;
      boundary_edges.insert({std::min(a, b), std::max(a, b)});
    }
  static const int te[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  Index total = 0;
  for (int round = 0; round < 100; ++round) {
    // Candidate edges, and a greedy set of them such that no tet contains two.
    std::map<std::pair<int, int>, std::vector<int>> edge_tets;
    for (Index t = 0; t < m.T.rows(); ++t)
      for (const auto& e : te) {
        int a = m.T(t, e[0]), b = m.T(t, e[1]);
        if (!on_boundary[static_cast<size_t>(a)] || !on_boundary[static_cast<size_t>(b)]) continue;
        std::pair<int, int> key{std::min(a, b), std::max(a, b)};
        if (boundary_edges.count(key)) continue;
        edge_tets[key].push_back(static_cast<int>(t));
      }
    if (edge_tets.empty()) break;
    std::vector<char> tet_used(static_cast<size_t>(m.T.rows()), 0);
    std::vector<std::pair<std::pair<int, int>, const std::vector<int>*>> batch;
    for (const auto& [key, tets] : edge_tets) {
      bool free = true;
      for (int t : tets) free = free && !tet_used[static_cast<size_t>(t)];
      if (!free) continue;
      for (int t : tets) tet_used[static_cast<size_t>(t)] = 1;
      batch.emplace_back(key, &tets);
    }
    // Apply the splits: each tet (a, b, c, d) around edge (a, b) becomes (m, b, c, d) and
    // (a, m, c, d), which keeps the orientation.
    Index nX = m.X.rows(), nT = m.T.rows(), extra_tets = 0;
    for (const auto& item : batch) extra_tets += static_cast<Index>(item.second->size());
    m.X.conservativeResize(nX + static_cast<Index>(batch.size()), 3);
    m.T.conservativeResize(nT + extra_tets, 4);
    on_boundary.resize(static_cast<size_t>(nX) + batch.size(), 0);
    Index next_t = nT;
    for (size_t i = 0; i < batch.size(); ++i) {
      const auto [a, b] = batch[i].first;
      const int mid = static_cast<int>(nX + static_cast<Index>(i));
      m.X.row(mid) = 0.5 * (m.X.row(a) + m.X.row(b));
      for (int t : *batch[i].second) {
        Eigen::RowVector4i tet = m.T.row(t);
        Eigen::RowVector4i t1 = tet, t2 = tet;
        for (int k = 0; k < 4; ++k) {
          if (tet(k) == a) t1(k) = mid;
          if (tet(k) == b) t2(k) = mid;
        }
        m.T.row(t) = t1;
        m.T.row(next_t++) = t2;
      }
    }
    total += static_cast<Index>(batch.size());
  }
  return total;
}

TetMesh tetrahedralize(const SurfaceMesh& surface, const TetParams& params, TetReport* report) {
  const Index nv = surface.V.rows(), nf = surface.F.rows();
  if (nv == 0 || nf == 0) throw Error("tetrahedralize: empty surface");

  double max_vol = params.max_volume;
  if (max_vol <= 0.0) {
    double sum = 0.0;
    for (Index f = 0; f < nf; ++f)
      for (int k = 0; k < 3; ++k)
        sum += (surface.V.row(surface.F(f, k)) - surface.V.row(surface.F(f, (k + 1) % 3))).norm();
    double edge = params.size_factor * sum / (3.0 * static_cast<double>(nf));
    max_vol = edge * edge * edge / (6.0 * std::sqrt(2.0));
  }

  char buf[128];
  if (params.min_dihedral_deg > 0.0)
    std::snprintf(buf, sizeof(buf), "pq%.6g/%.6ga%.9gYQ", params.radius_edge_ratio,
                  params.min_dihedral_deg, max_vol);
  else
    std::snprintf(buf, sizeof(buf), "pq%.6ga%.9gYQ", params.radius_edge_ratio, max_vol);
  const std::string switches = buf;

  Eigen::MatrixXd V = surface.V;
  Eigen::MatrixXi F = surface.F;
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT, TF;
  int rc = igl::copyleft::tetgen::tetrahedralize(V, F, switches, TV, TT, TF);
  if (rc != 0 || TT.rows() == 0) {
    throw Error("tetrahedralize: TetGen failed (code " + std::to_string(rc) +
                "); the surface may be self-intersecting");
  }

  TetMesh m;
  m.X = TV;
  m.T = TT;

  // With -Y the input vertices come first and unchanged; check it.
  if (m.X.rows() < nv || (m.X.topRows(nv) - surface.V).cwiseAbs().maxCoeff() != 0.0)
    throw Error("tetrahedralize: TetGen did not preserve the input vertices");

  // Orient all tets positively.
  for (Index t = 0; t < m.T.rows(); ++t) {
    double vol = tet_volume(m.X, m.T, t);
    if (vol < 0.0) std::swap(m.T(t, 2), m.T(t, 3));
  }

  // Boundary faces of the tet mesh must be exactly the input triangles.
  std::map<std::array<int, 3>, int> input_face;
  for (Index f = 0; f < nf; ++f) {
    std::array<int, 3> key{surface.F(f, 0), surface.F(f, 1), surface.F(f, 2)};
    std::sort(key.begin(), key.end());
    input_face[key] = static_cast<int>(f);
  }
  std::map<std::array<int, 3>, int> face_count;
  static const int tf[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
  for (Index t = 0; t < m.T.rows(); ++t)
    for (const auto& f : tf) {
      std::array<int, 3> key{m.T(t, f[0]), m.T(t, f[1]), m.T(t, f[2])};
      std::sort(key.begin(), key.end());
      ++face_count[key];
    }
  Index n_boundary = 0;
  for (const auto& [key, cnt] : face_count) {
    if (cnt == 1) {
      ++n_boundary;
      if (!input_face.count(key))
        throw Error("tetrahedralize: TetGen modified the boundary (unexpected boundary face)");
    }
  }
  if (n_boundary != nf)
    throw Error("tetrahedralize: boundary has " + std::to_string(n_boundary) +
                " faces, expected " + std::to_string(nf));

  m.B = surface.F;  // input orientation is outward (validated positive volume)
  m.B_to_input = Eigen::VectorXi::LinSpaced(nf, 0, static_cast<int>(nf) - 1);
  const Index n_split = params.split_boundary_edges ? split_boundary_spanning_edges(m) : 0;

  double vol = 0.0;
  for (Index t = 0; t < m.T.rows(); ++t) {
    double v = tet_volume(m.X, m.T, t);
    if (!(v > 0.0)) throw Error("tetrahedralize: degenerate tet " + std::to_string(t));
    vol += v;
  }
  double svol = 0.0;
  for (Index f = 0; f < nf; ++f) {
    Eigen::Vector3d a = surface.V.row(surface.F(f, 0)), b = surface.V.row(surface.F(f, 1)),
                    c = surface.V.row(surface.F(f, 2));
    svol += a.dot(b.cross(c)) / 6.0;
  }
  if (std::abs(vol - svol) > 1e-9 * std::abs(svol))
    throw Error("tetrahedralize: tet volume does not match enclosed surface volume");

  if (report) {
    report->n_nodes = m.X.rows();
    report->n_tets = m.T.rows();
    report->n_boundary = n_boundary;
    DihedralStats ds = dihedral_stats(m.X, m.T);
    report->min_dihedral_deg = ds.min_deg;
    report->max_dihedral_deg = ds.max_deg;
    report->n_small_dihedral = ds.n_below_5deg;
    report->volume = vol;
    report->surface_volume = svol;
    report->max_volume_used = max_vol;
    report->split_edges = n_split;
    report->switches = switches;
  }
  return m;
}

}  // namespace pcp
