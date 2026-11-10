#include <doctest/doctest.h>

#include <set>

#include "helpers.h"
#include "pcp/io/off.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {
void check_tet_mesh(const SurfaceMesh& s, const TetMesh& t) {
  REQUIRE(t.B.rows() == s.F.rows());
  for (Index b = 0; b < t.B.rows(); ++b) {
    int f = t.B_to_input(b);
    for (int k = 0; k < 3; ++k)
      REQUIRE((t.X.row(t.B(b, k)) - s.V.row(s.F(f, k))).norm() == 0.0);
  }
  for (Index i = 0; i < t.T.rows(); ++i) REQUIRE(tet_volume(t.X, t.T, i) > 0.0);
}
}  // namespace

TEST_CASE("tetrahedralize: unit box") {
  SurfaceMesh s = testutil::make_box(4);
  TetReport rep;
  TetMesh t = tetrahedralize(s, TetParams{}, &rep);
  check_tet_mesh(s, t);
  CHECK(rep.volume == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(rep.min_dihedral_deg > 0.0);
}

TEST_CASE("tetrahedralize: smaller max volume gives more tets") {
  SurfaceMesh s = testutil::make_box(4);
  TetReport a, b;
  tetrahedralize(s, TetParams{}, &a);
  TetParams p;
  p.max_volume = a.max_volume_used / 8.0;
  tetrahedralize(s, p, &b);
  CHECK(b.n_tets > a.n_tets);
}

TEST_CASE("tetrahedralize: provided models" * doctest::test_suite("models")) {
  for (const char* name : {"I-SPAR-1B59CD84.off", "CLEVIS-D64EEBFA.off", "TAILSITTER.off"}) {
    CAPTURE(name);
    SurfaceMesh s = read_off(testutil::models_dir() / name);
    TetReport rep;
    TetMesh t = tetrahedralize(s, TetParams{}, &rep);
    check_tet_mesh(s, t);
    CHECK(std::abs(rep.volume - rep.surface_volume) <= 1e-9 * rep.surface_volume);
    MESSAGE(std::string(name) << ": nodes=" << rep.n_nodes << " tets=" << rep.n_tets
                 << " min_dihedral=" << rep.min_dihedral_deg
                 << " below5deg=" << rep.n_small_dihedral);
  }
}

TEST_CASE("tetrahedralize: self-intersecting input fails with an error") {
  SurfaceMesh a = testutil::make_box(2), b = testutil::make_box(2);
  b.V.rowwise() += Eigen::RowVector3d(0.5, 0.5, 0.5);
  SurfaceMesh s;
  s.V.resize(a.V.rows() + b.V.rows(), 3);
  s.V << a.V, b.V;
  s.F.resize(a.F.rows() + b.F.rows(), 3);
  s.F << a.F, (b.F.array() + static_cast<int>(a.V.rows())).matrix();
  CHECK_THROWS_AS(tetrahedralize(s, TetParams{}), Error);
}

TEST_CASE("split_boundary_spanning_edges: no interior edge joins two boundary vertices") {
  SurfaceMesh s = testutil::make_box(3);
  TetParams p;
  p.split_boundary_edges = false;
  TetMesh t = tetrahedralize(s, p);
  const double vol_before = [&] {
    double v = 0;
    for (Index i = 0; i < t.T.rows(); ++i) v += tet_volume(t.X, t.T, i);
    return v;
  }();
  Index n = split_boundary_spanning_edges(t);
  CHECK(n > 0);
  std::vector<char> bnd(static_cast<size_t>(t.X.rows()), 0);
  std::set<std::pair<int, int>> bedges;
  for (Index f = 0; f < t.B.rows(); ++f)
    for (int k = 0; k < 3; ++k) {
      int a = t.B(f, k), b = t.B(f, (k + 1) % 3);
      bnd[static_cast<size_t>(a)] = 1;
      bedges.insert({std::min(a, b), std::max(a, b)});
    }
  double vol = 0;
  for (Index i = 0; i < t.T.rows(); ++i) {
    REQUIRE(tet_volume(t.X, t.T, i) > 0.0);
    vol += tet_volume(t.X, t.T, i);
    for (int a = 0; a < 4; ++a)
      for (int b = a + 1; b < 4; ++b) {
        int u = t.T(i, a), v = t.T(i, b);
        if (bnd[static_cast<size_t>(u)] && bnd[static_cast<size_t>(v)])
          REQUIRE(bedges.count({std::min(u, v), std::max(u, v)}) == 1);
      }
  }
  CHECK(vol == doctest::Approx(vol_before).epsilon(1e-12));
}
