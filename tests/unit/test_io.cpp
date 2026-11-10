#include <doctest/doctest.h>

#include <fstream>
#include <sstream>

#include "helpers.h"
#include "pcp/io/off.h"
#include "pcp/io/tmsh.h"
#include "pcp/io/validate.h"
#include "pcp/io/vtu.h"

using namespace pcp;

TEST_CASE("read_off: provided models have header counts") {
  struct Item {
    const char* file;
    Index nv, nf;
    int genus;
  } items[] = {{"CLEVIS-D64EEBFA.off", 20840, 41684, 2},
               {"I-SPAR-1B59CD84.off", 14767, 29530, 0},
               {"TAILSITTER.off", 28326, 56668, 5}};
  for (const auto& it : items) {
    CAPTURE(it.file);
    SurfaceMesh m = read_off(testutil::models_dir() / it.file);
    CHECK(m.V.rows() == it.nv);
    CHECK(m.F.rows() == it.nf);
    SurfaceReport r = validate_surface(m);
    CHECK(r.ok);
    CHECK(r.components == 1);
    CHECK(r.genus == it.genus);
  }
}

TEST_CASE("read_off: malformed files are rejected with a reason") {
  auto dir = testutil::temp_dir();
  testutil::write_text(dir / "quad.off", "OFF\n4 1 0\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n4 0 1 2 3\n");
  testutil::write_text(dir / "range.off", "OFF\n3 1 0\n0 0 0\n1 0 0\n1 1 0\n3 0 1 7\n");
  testutil::write_text(dir / "trunc.off", "OFF\n3 1 0\n0 0 0\n1 0 0\n");
  auto message = [](const std::filesystem::path& p) {
    try {
      read_off(p);
    } catch (const Error& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  CHECK(message(dir / "quad.off").find("only triangles") != std::string::npos);
  CHECK(message(dir / "range.off").find("out of range") != std::string::npos);
  CHECK(message(dir / "trunc.off").find("unexpected end of file") != std::string::npos);
  CHECK_THROWS_AS(read_off(dir / "does_not_exist.off"), Error);
}

TEST_CASE("write_off/read_off round trip is lossless") {
  SurfaceMesh m = read_off(testutil::models_dir() / "TAILSITTER.off");
  auto p = testutil::temp_dir() / "roundtrip.off";
  write_off(p, m);
  SurfaceMesh m2 = read_off(p);
  REQUIRE(m2.V.rows() == m.V.rows());
  REQUIRE(m2.F.rows() == m.F.rows());
  CHECK(m2.F == m.F);
  double scale = m.V.cwiseAbs().maxCoeff();
  CHECK((m2.V - m.V).cwiseAbs().maxCoeff() <= 1e-12 * scale);
}

TEST_CASE("validate_surface: genus of simple shapes") {
  SurfaceReport cube = validate_surface(testutil::make_box(3));
  CHECK(cube.ok);
  CHECK(cube.genus == 0);
  CHECK(cube.volume == doctest::Approx(1.0));

  SurfaceReport torus = validate_surface(testutil::make_torus(24, 12));
  CHECK(torus.ok);
  CHECK(torus.genus == 1);
}

TEST_CASE("validate_surface: open surface reports boundary edges") {
  SurfaceMesh m = testutil::make_box(2);
  // Drop the first face pair (two triangles of one quad) to open the surface.
  SurfaceMesh open;
  open.V = m.V;
  open.F = m.F.bottomRows(m.F.rows() - 2);
  SurfaceReport r = validate_surface(open);
  CHECK_FALSE(r.ok);
  CHECK(r.boundary_edges == 4);
  REQUIRE_FALSE(r.errors.empty());
  CHECK(r.errors[0].find("boundary edges") != std::string::npos);
}

TEST_CASE("vtu export writes all cells and fields") {
  SurfaceMesh m = testutil::make_box(2);
  Eigen::VectorXi labels = Eigen::VectorXi::LinSpaced(m.F.rows(), 0, static_cast<int>(m.F.rows()) - 1);
  Eigen::VectorXd area = Eigen::VectorXd::Constant(m.F.rows(), 0.5);
  auto p = testutil::temp_dir() / "surface.vtu";
  write_vtu_surface(p, m.V, m.F, {CellField::integer("label", labels), CellField::scalar("area", area)});
  std::ifstream in(p);
  std::stringstream ss;
  ss << in.rdbuf();
  std::string s = ss.str();
  CHECK(s.find("NumberOfCells=\"" + std::to_string(m.F.rows()) + "\"") != std::string::npos);
  CHECK(s.find("Name=\"label\"") != std::string::npos);
  CHECK(s.find("Name=\"area\"") != std::string::npos);
  CHECK_THROWS_AS(write_vtu_surface(p, m.V, m.F, {CellField::scalar("bad", Eigen::VectorXd(3))}), Error);
}

TEST_CASE("tmsh round trip is bit exact") {
  TetMesh t;
  t.X.resize(4, 3);
  t.X << 0.1, 0.2, 1.0 / 3.0, 1, 0, 0, 0, 1, 0, 0, 0, 3.14159265358979;
  t.T.resize(1, 4);
  t.T << 0, 1, 2, 3;
  t.B.resize(1, 3);
  t.B << 0, 2, 1;
  t.B_to_input.resize(1);
  t.B_to_input << 7;
  auto p = testutil::temp_dir() / "mesh.tmsh";
  write_tmsh(p, t);
  TetMesh r = read_tmsh(p);
  CHECK(r.X == t.X);
  CHECK(r.T == t.T);
  CHECK(r.B == t.B);
  CHECK(r.B_to_input == t.B_to_input);
}
