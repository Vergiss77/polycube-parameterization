#include "pcp/io/off.h"

#include "pcp/io/text.h"

namespace pcp {

SurfaceMesh read_off(const std::filesystem::path& path) {
  const std::string data = text::read_file(path);
  text::Tokenizer tk(data);
  const std::string where = path.string() + ": ";
  try {
    std::string_view head = tk.expect("OFF header");
    if (head != "OFF") throw Error("missing 'OFF' header");
    long long nv = tk.expect_int("vertex count");
    long long nf = tk.expect_int("face count");
    tk.expect_int("edge count");
    if (nv < 0 || nf < 0) throw Error("negative element count");

    SurfaceMesh m;
    m.V.resize(nv, 3);
    for (long long i = 0; i < nv; ++i)
      for (int k = 0; k < 3; ++k) m.V(i, k) = tk.expect_double("vertex coordinate");

    m.F.resize(nf, 3);
    for (long long f = 0; f < nf; ++f) {
      long long deg = tk.expect_int("face size");
      if (deg != 3) {
        throw Error("face " + std::to_string(f) + " has " + std::to_string(deg) +
                    " vertices; only triangles are supported");
      }
      for (int k = 0; k < 3; ++k) {
        long long idx = tk.expect_int("face index");
        if (idx < 0 || idx >= nv) {
          throw Error("face " + std::to_string(f) + " references vertex " + std::to_string(idx) +
                      " out of range [0, " + std::to_string(nv) + ")");
        }
        m.F(f, k) = static_cast<int>(idx);
      }
    }
    return m;
  } catch (const Error& e) {
    throw Error(where + e.what());
  }
}

void write_off(const std::filesystem::path& path, const SurfaceMesh& mesh) {
  std::string s;
  s.reserve(static_cast<size_t>(mesh.V.rows() * 64 + mesh.F.rows() * 24 + 32));
  s += "OFF\n";
  text::append_int(s, mesh.V.rows());
  s += ' ';
  text::append_int(s, mesh.F.rows());
  s += " 0\n";
  for (Index i = 0; i < mesh.V.rows(); ++i) {
    for (int k = 0; k < 3; ++k) {
      if (k) s += ' ';
      text::append_double(s, mesh.V(i, k));
    }
    s += '\n';
  }
  for (Index f = 0; f < mesh.F.rows(); ++f) {
    s += '3';
    for (int k = 0; k < 3; ++k) {
      s += ' ';
      text::append_int(s, mesh.F(f, k));
    }
    s += '\n';
  }
  text::write_file(path, s);
}

}  // namespace pcp
