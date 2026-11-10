#include "pcp/io/tmsh.h"

#include "pcp/io/text.h"

namespace pcp {

void write_tmsh(const std::filesystem::path& path, const TetMesh& m) {
  std::string s;
  s.reserve(static_cast<size_t>(m.X.rows() * 64 + m.T.rows() * 32 + m.B.rows() * 32 + 64));
  s += "TMSH 1\nnodes ";
  text::append_int(s, m.X.rows());
  s += '\n';
  for (Index i = 0; i < m.X.rows(); ++i) {
    for (int k = 0; k < 3; ++k) {
      if (k) s += ' ';
      text::append_double(s, m.X(i, k));
    }
    s += '\n';
  }
  s += "tets ";
  text::append_int(s, m.T.rows());
  s += '\n';
  for (Index t = 0; t < m.T.rows(); ++t) {
    for (int k = 0; k < 4; ++k) {
      if (k) s += ' ';
      text::append_int(s, m.T(t, k));
    }
    s += '\n';
  }
  s += "boundary ";
  text::append_int(s, m.B.rows());
  s += '\n';
  for (Index b = 0; b < m.B.rows(); ++b) {
    for (int k = 0; k < 3; ++k) {
      text::append_int(s, m.B(b, k));
      s += ' ';
    }
    text::append_int(s, m.B_to_input.size() ? m.B_to_input(b) : -1);
    s += '\n';
  }
  text::write_file(path, s);
}

TetMesh read_tmsh(const std::filesystem::path& path) {
  const std::string data = text::read_file(path);
  text::Tokenizer tk(data);
  try {
    if (tk.expect("TMSH header") != "TMSH") throw Error("missing 'TMSH' header");
    if (tk.expect_int("version") != 1) throw Error("unsupported TMSH version");
    TetMesh m;
    if (tk.expect("nodes keyword") != "nodes") throw Error("expected 'nodes'");
    long long n = tk.expect_int("node count");
    m.X.resize(n, 3);
    for (long long i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) m.X(i, k) = tk.expect_double("node coordinate");
    if (tk.expect("tets keyword") != "tets") throw Error("expected 'tets'");
    long long nt = tk.expect_int("tet count");
    m.T.resize(nt, 4);
    for (long long t = 0; t < nt; ++t)
      for (int k = 0; k < 4; ++k) {
        long long v = tk.expect_int("tet index");
        if (v < 0 || v >= n) throw Error("tet index out of range");
        m.T(t, k) = static_cast<int>(v);
      }
    if (tk.expect("boundary keyword") != "boundary") throw Error("expected 'boundary'");
    long long nb = tk.expect_int("boundary count");
    m.B.resize(nb, 3);
    m.B_to_input.resize(nb);
    for (long long b = 0; b < nb; ++b) {
      for (int k = 0; k < 3; ++k) {
        long long v = tk.expect_int("boundary index");
        if (v < 0 || v >= n) throw Error("boundary index out of range");
        m.B(b, k) = static_cast<int>(v);
      }
      m.B_to_input(b) = static_cast<int>(tk.expect_int("input face index"));
    }
    return m;
  } catch (const Error& e) {
    throw Error(path.string() + ": " + e.what());
  }
}

}  // namespace pcp
