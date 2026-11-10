#include "pcp/blocks/block_layout.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <Eigen/Geometry>
#include <igl/point_mesh_squared_distance.h>
#include <igl/winding_number.h>

namespace pcp {

int BlockLayout::interior_faces() const {
  int n = 0;
  for (const Block& b : blocks)
    for (const BlockFace& f : b.faces) n += f.interior() ? 1 : 0;
  return n / 2;
}

int BlockLayout::boundary_faces() const {
  int n = 0;
  for (const Block& b : blocks)
    for (const BlockFace& f : b.faces) n += f.interior() ? 0 : 1;
  return n;
}

double BlockLayout::blocks_volume() const {
  double v = 0.0;
  for (const Block& b : blocks) v += b.volume;
  return v;
}

nlohmann::json BlockLayout::to_json(double scale, const Eigen::RowVector3d& center) const {
  nlohmann::json cuts_j = nlohmann::json::array();
  for (int a = 0; a < 3; ++a) {
    nlohmann::json c = nlohmann::json::array();
    for (double h : cuts[static_cast<size_t>(a)]) c.push_back(h / scale + center(a));
    cuts_j.push_back(c);
  }
  static const char* side_names[6] = {"x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"};
  nlohmann::json blocks_j = nlohmann::json::array();
  for (size_t i = 0; i < blocks.size(); ++i) {
    const Block& b = blocks[i];
    nlohmann::json faces = nlohmann::json::object();
    for (int f = 0; f < 6; ++f) {
      const BlockFace& bf = b.faces[static_cast<size_t>(f)];
      faces[side_names[f]] = bf.interior() ? nlohmann::json{{"neighbor", bf.neighbor}}
                                           : nlohmann::json{{"patch", bf.patch}};
    }
    blocks_j.push_back({{"id", i},
                        {"cell", {b.cell[0], b.cell[1], b.cell[2]}},
                        {"lo", {b.lo(0), b.lo(1), b.lo(2)}},
                        {"hi", {b.hi(0), b.hi(1), b.hi(2)}},
                        {"volume", b.volume},
                        {"faces", faces}});
  }
  return {{"cuts", cuts_j},
          {"blocks", blocks_j},
          {"grid_cells", grid_cells()},
          {"interior_faces", interior_faces()},
          {"boundary_faces", boundary_faces()},
          {"note", "cuts in output units (rotated polycube frame); lo/hi/volume in normalized polycube units"}};
}

std::array<std::vector<double>, 3> cut_coordinates(const MatX3d& X, const std::vector<int>& patch_label,
                                                   const std::vector<double>& patch_height) {
  const double diag = (X.colwise().maxCoeff() - X.colwise().minCoeff()).norm();
  const double tol = 1e-12 * diag;
  std::array<std::vector<double>, 3> cuts;
  for (size_t p = 0; p < patch_label.size(); ++p)
    cuts[static_cast<size_t>(patch_label[p] / 2)].push_back(patch_height[p]);
  for (auto& c : cuts) {
    std::sort(c.begin(), c.end());
    std::vector<double> u;
    for (double h : c)
      if (u.empty() || h - u.back() > tol) u.push_back(h);
    c = u;
  }
  return cuts;
}

BlockLayout build_block_layout(const MatX3d& X, const MatX3i& B, const Eigen::VectorXi& face_patch,
                               const std::vector<int>& patch_label,
                               const std::array<std::vector<double>, 3>& cuts) {
  BlockLayout L;
  L.cuts = cuts;
  for (int a = 0; a < 3; ++a)
    if (L.intervals(a) < 1) throw Error("block layout: fewer than two cut coordinates on axis " + std::to_string(a));

  for (Index f = 0; f < B.rows(); ++f) {
    const Eigen::Vector3d v0 = X.row(B(f, 0)), v1 = X.row(B(f, 1)), v2 = X.row(B(f, 2));
    L.polycube_volume += v0.dot(v1.cross(v2)) / 6.0;
  }

  // Inside test of every grid cell center (generalized winding number).
  const int nx = L.intervals(0), ny = L.intervals(1), nz = L.intervals(2);
  const int ncells = L.grid_cells();
  auto mid = [&](int a, int i) { return 0.5 * (cuts[static_cast<size_t>(a)][static_cast<size_t>(i)] + cuts[static_cast<size_t>(a)][static_cast<size_t>(i) + 1]); };
  Eigen::MatrixXd O(ncells, 3);
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) O.row(L.cell_index(i, j, k)) << mid(0, i), mid(1, j), mid(2, k);
  const Eigen::MatrixXd V = X;
  const Eigen::MatrixXi F = B;
  Eigen::VectorXd W;
  igl::winding_number(V, F, O, W);

  L.cell_block.assign(static_cast<size_t>(ncells), -1);
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const int c = L.cell_index(i, j, k);
        if (W(c) <= 0.5) continue;
        Block b;
        b.cell = {i, j, k};
        b.lo << cuts[0][static_cast<size_t>(i)], cuts[1][static_cast<size_t>(j)], cuts[2][static_cast<size_t>(k)];
        b.hi << cuts[0][static_cast<size_t>(i) + 1], cuts[1][static_cast<size_t>(j) + 1], cuts[2][static_cast<size_t>(k) + 1];
        b.volume = (b.hi - b.lo).prod();
        L.cell_block[static_cast<size_t>(c)] = static_cast<int>(L.blocks.size());
        L.blocks.push_back(b);
      }

  // Faces: interior if the neighboring cell is a block, otherwise on the boundary.
  std::vector<std::pair<int, int>> boundary;  // (block, face)
  for (size_t bi = 0; bi < L.blocks.size(); ++bi) {
    Block& b = L.blocks[bi];
    for (int f = 0; f < 6; ++f) {
      const int a = f / 2;
      std::array<int, 3> nc = b.cell;
      nc[static_cast<size_t>(a)] += (f % 2 == 0) ? -1 : 1;
      int nb = -1;
      if (nc[static_cast<size_t>(a)] >= 0 && nc[static_cast<size_t>(a)] < L.intervals(a))
        nb = L.cell_block[static_cast<size_t>(L.cell_index(nc[0], nc[1], nc[2]))];
      if (nb >= 0) b.faces[static_cast<size_t>(f)].neighbor = nb;
      else boundary.emplace_back(static_cast<int>(bi), f);
    }
  }

  // Patch of every boundary face: nearest boundary triangle to the face center.
  if (!boundary.empty()) {
    Eigen::MatrixXd P(static_cast<Index>(boundary.size()), 3);
    for (size_t q = 0; q < boundary.size(); ++q) {
      const Block& b = L.blocks[static_cast<size_t>(boundary[q].first)];
      const int f = boundary[q].second, a = f / 2;
      Eigen::Vector3d c = 0.5 * (b.lo + b.hi);
      c(a) = (f % 2 == 0) ? b.lo(a) : b.hi(a);
      P.row(static_cast<Index>(q)) = c.transpose();
    }
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);
    const double diag = (X.colwise().maxCoeff() - X.colwise().minCoeff()).norm();
    for (size_t q = 0; q < boundary.size(); ++q) {
      const int bi = boundary[q].first, f = boundary[q].second, a = f / 2;
      const int patch = face_patch(I(static_cast<Index>(q)));
      const int want = 2 * a + (f % 2 == 0 ? 1 : 0);  // outward normal -axis on the low side
      if (patch_label[static_cast<size_t>(patch)] != want || std::sqrt(sqrD(static_cast<Index>(q))) > 1e-6 * diag)
        throw Error("block layout: boundary face " + std::to_string(f) + " of block " + std::to_string(bi) +
                    " does not lie on a patch of the matching orientation (polycube is not a union of grid cells)");
      L.blocks[static_cast<size_t>(bi)].faces[static_cast<size_t>(f)].patch = patch;
    }
  }
  return L;
}

Eigen::VectorXi assign_tets_to_blocks(const BlockLayout& L, const MatX3d& X, const MatX4i& T) {
  Eigen::VectorXi block(T.rows());
  for (Index t = 0; t < T.rows(); ++t) {
    const Eigen::Vector3d c = 0.25 * (X.row(T(t, 0)) + X.row(T(t, 1)) + X.row(T(t, 2)) + X.row(T(t, 3))).transpose();
    std::array<int, 3> cell{};
    for (int a = 0; a < 3; ++a) {
      const auto& cu = L.cuts[static_cast<size_t>(a)];
      const int i = static_cast<int>(std::upper_bound(cu.begin(), cu.end(), c(a)) - cu.begin()) - 1;
      cell[static_cast<size_t>(a)] = std::clamp(i, 0, L.intervals(a) - 1);
    }
    int b = L.cell_block[static_cast<size_t>(L.cell_index(cell[0], cell[1], cell[2]))];
    if (b < 0) {  // nearest block box
      double best = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < L.blocks.size(); ++i) {
        const Eigen::Vector3d d = (L.blocks[i].lo - c).cwiseMax(c - L.blocks[i].hi).cwiseMax(0.0);
        if (d.squaredNorm() < best) {
          best = d.squaredNorm();
          b = static_cast<int>(i);
        }
      }
    }
    block(t) = b;
  }
  return block;
}

}  // namespace pcp
