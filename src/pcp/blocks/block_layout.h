#pragma once

// Block decomposition of an axis-aligned polycube (design D1, D2): the heights of the patch
// groups of every axis are the cut coordinates; the cells of the resulting global grid that lie
// inside the polycube are the blocks. Neighboring blocks share whole faces (no T-junctions).

#include <array>
#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/core/types.h"

namespace pcp {

// Face f of a block: f = 2 * axis + (0 for the low side, 1 for the high side). The G+Smo box
// side of face f is f + 1 (west, east, south, north, front, back).
struct BlockFace {
  int neighbor = -1;  // block on the other side (interior face), or -1
  int patch = -1;     // patch containing the face (boundary face), or -1
  bool interior() const { return neighbor >= 0; }
};

struct Block {
  std::array<int, 3> cell{};  // interval index along x, y, z
  Eigen::Vector3d lo = Eigen::Vector3d::Zero(), hi = Eigen::Vector3d::Zero();  // polycube box
  double volume = 0.0;
  std::array<BlockFace, 6> faces;
};

struct BlockLayout {
  std::array<std::vector<double>, 3> cuts;  // ascending cut coordinates of every axis
  std::vector<Block> blocks;                // ordered lexicographically by (z, y, x) cell
  std::vector<int> cell_block;              // block of every grid cell, -1 if outside
  double polycube_volume = 0.0;             // enclosed by the boundary triangles

  int intervals(int axis) const { return static_cast<int>(cuts[static_cast<size_t>(axis)].size()) - 1; }
  int grid_cells() const { return intervals(0) * intervals(1) * intervals(2); }
  int cell_index(int i, int j, int k) const { return (k * intervals(1) + j) * intervals(0) + i; }
  int interior_faces() const;  // each shared face counted once
  int boundary_faces() const;
  double blocks_volume() const;
  // Cut coordinates are written as c / scale + center(axis) (output units, see Normalization).
  nlohmann::json to_json(double scale, const Eigen::RowVector3d& center) const;
};

// Groups the patch heights of every axis (labels 0..5 = +X,-X,+Y,-Y,+Z,-Z) into cut
// coordinates: equal up to 1e-12 times the bounding box diagonal of X.
std::array<std::vector<double>, 3> cut_coordinates(const MatX3d& X, const std::vector<int>& patch_label,
                                                   const std::vector<double>& patch_height);

// Builds the blocks of the axis-aligned polycube (X, B) from the given cuts. face_patch is the
// patch of every boundary triangle. Throws pcp::Error if the polycube is not a union of grid
// cells (a boundary block face lies on no patch of the matching orientation).
BlockLayout build_block_layout(const MatX3d& X, const MatX3i& B, const Eigen::VectorXi& face_patch,
                               const std::vector<int>& patch_label,
                               const std::array<std::vector<double>, 3>& cuts);

// Block of every tet: the block containing its centroid in X (the nearest block if the centroid
// falls outside all blocks, e.g. near inverted tets).
Eigen::VectorXi assign_tets_to_blocks(const BlockLayout& layout, const MatX3d& X, const MatX4i& T);

}  // namespace pcp
