#include "pcp/io/vtu.h"

#include "pcp/io/text.h"

namespace pcp {

CellField CellField::scalar(std::string name, const Eigen::VectorXd& v) {
  CellField f;
  f.name = std::move(name);
  f.type = Type::Float64;
  f.values.assign(v.data(), v.data() + v.size());
  return f;
}

CellField CellField::integer(std::string name, const Eigen::VectorXi& v) {
  CellField f;
  f.name = std::move(name);
  f.type = Type::Int32;
  f.values.resize(static_cast<size_t>(v.size()));
  for (Index i = 0; i < v.size(); ++i) f.values[static_cast<size_t>(i)] = v(i);
  return f;
}

namespace {

template <int K>
void write_vtu(const std::filesystem::path& path, const MatX3d& X,
               const Eigen::Matrix<int, Eigen::Dynamic, K>& C, int vtk_type,
               const std::vector<CellField>& fields) {
  const Index np = X.rows(), nc = C.rows();
  for (const auto& f : fields)
    if (static_cast<Index>(f.values.size()) != nc)
      throw Error("vtu field '" + f.name + "' has " + std::to_string(f.values.size()) +
                  " values but mesh has " + std::to_string(nc) + " cells");

  std::string s;
  s.reserve(static_cast<size_t>(np * 60 + nc * (K * 8 + 16) + 1024));
  s += "<?xml version=\"1.0\"?>\n";
  s += "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  s += "<UnstructuredGrid>\n<Piece NumberOfPoints=\"";
  text::append_int(s, np);
  s += "\" NumberOfCells=\"";
  text::append_int(s, nc);
  s += "\">\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (Index i = 0; i < np; ++i) {
    for (int k = 0; k < 3; ++k) {
      if (k) s += ' ';
      text::append_double(s, X(i, k));
    }
    s += '\n';
  }
  s += "</DataArray>\n</Points>\n<Cells>\n";
  s += "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (Index c = 0; c < nc; ++c) {
    for (int k = 0; k < K; ++k) {
      if (k) s += ' ';
      text::append_int(s, C(c, k));
    }
    s += '\n';
  }
  s += "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (Index c = 0; c < nc; ++c) {
    text::append_int(s, (c + 1) * K);
    s += '\n';
  }
  s += "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (Index c = 0; c < nc; ++c) {
    text::append_int(s, vtk_type);
    s += '\n';
  }
  s += "</DataArray>\n</Cells>\n";
  if (!fields.empty()) {
    s += "<CellData>\n";
    for (const auto& f : fields) {
      s += "<DataArray type=\"";
      s += f.type == CellField::Type::Int32 ? "Int32" : "Float64";
      s += "\" Name=\"" + f.name + "\" format=\"ascii\">\n";
      for (double v : f.values) {
        if (f.type == CellField::Type::Int32)
          text::append_int(s, static_cast<long long>(v));
        else
          text::append_double(s, v);
        s += '\n';
      }
      s += "</DataArray>\n";
    }
    s += "</CellData>\n";
  }
  s += "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  text::write_file(path, s);
}

}  // namespace

void write_vtu_tets(const std::filesystem::path& path, const MatX3d& X, const MatX4i& T,
                    const std::vector<CellField>& cell_fields) {
  write_vtu<4>(path, X, T, 10 /* VTK_TETRA */, cell_fields);
}

void write_vtu_surface(const std::filesystem::path& path, const MatX3d& V, const MatX3i& F,
                       const std::vector<CellField>& cell_fields) {
  write_vtu<3>(path, V, F, 5 /* VTK_TRIANGLE */, cell_fields);
}

}  // namespace pcp
