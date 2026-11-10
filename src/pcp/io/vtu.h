#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "pcp/core/types.h"

namespace pcp {

// Per-cell field attached to a .vtu export.
struct CellField {
  enum class Type { Float64, Int32 };
  std::string name;
  Type type = Type::Float64;
  std::vector<double> values;  // one value per cell (ints stored exactly as doubles)

  static CellField scalar(std::string name, const Eigen::VectorXd& v);
  static CellField integer(std::string name, const Eigen::VectorXi& v);
};

// Per-point field: same layout, one value per point.
using PointField = CellField;

using MatX8i = Eigen::Matrix<int, Eigen::Dynamic, 8>;

// ASCII VTK UnstructuredGrid writers (open directly in ParaView).
void write_vtu_tets(const std::filesystem::path& path, const MatX3d& X, const MatX4i& T,
                    const std::vector<CellField>& cell_fields = {});
void write_vtu_surface(const std::filesystem::path& path, const MatX3d& V, const MatX3i& F,
                       const std::vector<CellField>& cell_fields = {});
// Hexahedra in VTK order: bottom face (0,0,0), (1,0,0), (1,1,0), (0,1,0), then the top face.
void write_vtu_hexes(const std::filesystem::path& path, const MatX3d& X, const MatX8i& H,
                     const std::vector<CellField>& cell_fields = {},
                     const std::vector<PointField>& point_fields = {});

}  // namespace pcp
