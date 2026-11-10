// Block decomposition of a polycube (design D1, D2) and B-spline fitting on small polycubes.
#include <doctest/doctest.h>

#include <cmath>
#include <set>

#include <Eigen/Geometry>

#include "helpers.h"
#include "pcp/blocks/block_layout.h"
#include "pcp/blocks/spline_fit.h"
#include "pcp/geom/boundary_topology.h"
#include "pcp/io/gismo_xml.h"
#include "pcp/io/text.h"
#include "pcp/io/vtu.h"
#include "pcp/structure/labeling.h"
#include "pcp/structure/patch_graph.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {

struct Poly {
  MatX3d X;
  MatX3i B;
  PatchGraph g;
};

Poly make_poly(const SurfaceMesh& s) {
  Poly p{s.V, s.F, {}};
  BoundaryTopology topo = build_boundary_topology(s.F);
  p.g = build_patch_graph(s.V, s.F, topo, initial_labels(s.V, s.F));
  return p;
}

BlockLayout layout_of(const Poly& p) {
  return build_block_layout(p.X, p.B, p.g.patches.face_patch, p.g.patches.label,
                            cut_coordinates(p.X, p.g.patches.label, p.g.height));
}

void check_faces(const BlockLayout& L) {
  for (size_t i = 0; i < L.blocks.size(); ++i)
    for (int f = 0; f < 6; ++f) {
      const BlockFace& bf = L.blocks[i].faces[size_t(f)];
      CHECK((bf.interior() != (bf.patch >= 0)));
      if (bf.interior()) {
        // Symmetric: the neighbor's opposite face points back.
        const int opp = f ^ 1;
        CHECK(L.blocks[size_t(bf.neighbor)].faces[size_t(opp)].neighbor == int(i));
      }
    }
}

const std::vector<std::array<int, 3>> kL = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
const std::vector<std::array<int, 3>> kU = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 0, 1}, {2, 0, 1}};

}  // namespace

TEST_CASE("blocks: L-shaped polycube gives 3 blocks of total volume 3") {
  Poly p = make_poly(testutil::make_voxels(kL, 2));
  BlockLayout L = layout_of(p);
  CHECK(L.cuts[0].size() == 3);
  CHECK(L.cuts[1].size() == 3);
  CHECK(L.cuts[2].size() == 2);
  CHECK(L.grid_cells() == 4);
  REQUIRE(L.blocks.size() == 3);
  CHECK(L.blocks_volume() == doctest::Approx(3.0).epsilon(1e-12));
  CHECK(std::abs(L.polycube_volume - 3.0) < 1e-12);
  CHECK(L.interior_faces() == 2);
  CHECK(L.boundary_faces() == 3 * 6 - 4);
  check_faces(L);
  // Numbering by (z, y, x): (0,0,0), (1,0,0), (0,1,0).
  CHECK(L.blocks[0].cell == std::array<int, 3>{0, 0, 0});
  CHECK(L.blocks[1].cell == std::array<int, 3>{1, 0, 0});
  CHECK(L.blocks[2].cell == std::array<int, 3>{0, 1, 0});
}

TEST_CASE("blocks: U-shaped polycube, the gap is not a block, posts sit on the base") {
  Poly p = make_poly(testutil::make_voxels(kU, 2));
  BlockLayout L = layout_of(p);
  CHECK(L.grid_cells() == 6);
  REQUIRE(L.blocks.size() == 5);
  CHECK(L.cell_block[size_t(L.cell_index(1, 0, 1))] == -1);
  const int base_left = L.cell_block[size_t(L.cell_index(0, 0, 0))];
  const int post_left = L.cell_block[size_t(L.cell_index(0, 0, 1))];
  CHECK(L.blocks[size_t(base_left)].faces[5].neighbor == post_left);
  CHECK(L.blocks[size_t(post_left)].faces[4].neighbor == base_left);
  // The middle base block's top face is on the surface (the gap floor), a +Z patch.
  const int mid = L.cell_block[size_t(L.cell_index(1, 0, 0))];
  const BlockFace& top = L.blocks[size_t(mid)].faces[5];
  REQUIRE(top.patch >= 0);
  CHECK(p.g.patches.label[size_t(top.patch)] == 4);
  CHECK(L.interior_faces() == 4);
  CHECK(L.boundary_faces() == 22);
  check_faces(L);
}

TEST_CASE("blocks: volume matches the tet mesh, tets assigned, repeat run identical") {
  SurfaceMesh s = testutil::make_voxels(kU, 2);
  TetParams tp;
  TetMesh t = tetrahedralize(s, tp);
  Poly p = make_poly(SurfaceMesh{t.X, t.B});
  BlockLayout L = layout_of(p);
  double vol = 0.0;
  for (Index i = 0; i < t.T.rows(); ++i) vol += tet_volume(t.X, t.T, i);
  CHECK(std::abs(L.blocks_volume() - vol) <= 1e-9 * vol);
  Eigen::VectorXi tb = assign_tets_to_blocks(L, t.X, t.T);
  CHECK(tb.minCoeff() >= 0);
  CHECK(tb.maxCoeff() < int(L.blocks.size()));
  // Every block gets tets (each voxel is meshed).
  std::vector<int> count(L.blocks.size(), 0);
  for (Index i = 0; i < tb.size(); ++i) ++count[size_t(tb(i))];
  for (int c : count) CHECK(c > 0);

  BlockLayout L2 = layout_of(p);
  CHECK(L2.to_json(1.0, Eigen::RowVector3d::Zero()) == L.to_json(1.0, Eigen::RowVector3d::Zero()));
}

TEST_CASE("blocks: cut coordinates merge round-off differences only") {
  MatX3d X(2, 3);
  X << 0, 0, 0, 1, 1, 1;
  // Labels: +X, -X, +X, +Y; the +X/-X heights 1 and 1 + 1e-14 are one cut, 1 + 1e-3 another.
  std::vector<int> label = {0, 1, 0, 2};
  auto cuts = cut_coordinates(X, label, {1.0, 1.0 + 1e-14, 0.0, 0.5});
  CHECK(cuts[0] == std::vector<double>{0.0, 1.0});
  CHECK(cuts[1] == std::vector<double>{0.5});
  CHECK(cuts[2].empty());
  CHECK(cut_coordinates(X, label, {1.0, 1.0 + 1e-3, 0.0, 0.5})[0].size() == 3);
}

namespace {

// Tet mesh of a voxel polycube and an affine image of it as the fitting target.
struct Affine {
  TetMesh t;
  Poly p;
  BlockLayout L;
  Eigen::Matrix3d A;
  Eigen::RowVector3d c;
  MatX3d Y;
};

Affine make_affine(const std::vector<std::array<int, 3>>& voxels) {
  Affine a;
  a.t = tetrahedralize(testutil::make_voxels(voxels, 2), TetParams{});
  a.p = make_poly(SurfaceMesh{a.t.X, a.t.B});
  a.L = layout_of(a.p);
  a.A << 2, 0.5, 0, -0.3, 1, 0.2, 0.1, 0, 3;
  a.c << 1, -2, 0.5;
  a.Y = (a.t.X * a.A.transpose()).rowwise() + a.c;
  return a;
}

}  // namespace

TEST_CASE("spline space: element counts, shared control points, unknown count") {
  CHECK(elements_for_length(0.3, 1.0) == 1);
  CHECK(elements_for_length(1.4, 1.0) == 1);
  CHECK(elements_for_length(2.6, 1.0) == 3);
  Poly p = make_poly(testutil::make_voxels(kL, 2));
  BlockLayout L = layout_of(p);
  SplineSpace S = build_spline_space(L, 3, 10.0);  // one element per interval
  CHECK(S.lattice == std::array<int, 3>{7, 7, 4});
  CHECK(S.n_unknowns == 160);
  // Both sides of every interior face reference the same control points on that face.
  for (size_t bi = 0; bi < L.blocks.size(); ++bi)
    for (int f = 0; f < 6; ++f) {
      const int nb = L.blocks[bi].faces[size_t(f)].neighbor;
      if (nb < 0) continue;
      const int a = f / 2;
      auto face_dofs = [&](size_t b, bool high) {
        std::set<int> out;
        const auto kv = S.block_basis(L.blocks[b]);
        const int n[3] = {kv[0].n_basis(), kv[1].n_basis(), kv[2].n_basis()};
        for (int k = 0; k < n[2]; ++k)
          for (int j = 0; j < n[1]; ++j)
            for (int i = 0; i < n[0]; ++i) {
              const int idx[3] = {i, j, k};
              if (idx[a] == (high ? n[a] - 1 : 0)) out.insert(S.block_dofs[b][size_t(i + n[0] * (j + n[1] * k))]);
            }
        return out;
      };
      CHECK(face_dofs(bi, f % 2 == 1) == face_dofs(size_t(nb), f % 2 == 0));
    }
  // Finer elements: interval lengths 1 with h = 0.5 give 2 elements each.
  SplineSpace S2 = build_spline_space(L, 2, 0.5);
  CHECK(S2.elements[0] == std::vector<int>{2, 2});
  CHECK(S2.offset[0] == std::vector<int>{0, 3});
  CHECK(S2.lattice[0] == 7);
}

TEST_CASE("fit samples: affine target mapped exactly, surface samples counted") {
  Affine a = make_affine(kL);
  SplineSpace S = build_spline_space(a.L, 3, 10.0);
  FitSamples fs = sample_blocks(a.L, S, 4, a.t.X, a.t.T, a.Y);
  CHECK(fs.block.size() == 3 * 125);
  CHECK(fs.unlocated == 0);
  double err = 0.0;
  for (size_t i = 0; i < fs.block.size(); ++i) {
    const Block& b = a.L.blocks[size_t(fs.block[i])];
    const Eigen::Vector3d x = b.lo + fs.param[i].cwiseProduct(b.hi - b.lo);
    err = std::max(err, (fs.target[i] - (a.A * x + a.c.transpose())).norm());
  }
  CHECK(err <= 1e-12);
  long long surf = 0;
  for (char c : fs.surface) surf += c;
  CHECK(surf == 77 + 89 + 89);
}

TEST_CASE("fit: affine target reproduced, Jacobian and C0 checks, deterministic") {
  Affine a = make_affine(kL);
  SplineSpace S = build_spline_space(a.L, 3, 10.0);
  FitSamples fs = sample_blocks(a.L, S, 4, a.t.X, a.t.T, a.Y);
  const double diag = (a.Y.colwise().maxCoeff() - a.Y.colwise().minCoeff()).norm();
  FitSolution sol = fit_volumes(a.L, S, fs, 0.0);
  FitQuality q = evaluate_fit(a.L, S, sol, fs, a.Y, a.t.B, diag);
  CHECK(q.max_error <= 1e-9);
  CHECK(q.surface_distance_max <= 1e-9);
  CHECK(q.c0_max_gap <= 1e-12);
  CHECK(q.nonpositive_jacobian == 0);
  CHECK(std::abs(q.min_scaled_jacobian - scaled_jacobian(a.A)) <= 1e-9);
  CHECK(std::abs(q.mean_scaled_jacobian - scaled_jacobian(a.A)) <= 1e-9);

  FitSolution sol2 = fit_volumes(a.L, S, fs, 0.0);
  CHECK(sol2.controls == sol.controls);

  // A folded control point is reported.
  FitSolution bad = sol;
  TensorBSpline3& V = bad.volumes[0];
  for (int k = 0; k < V.n(2); ++k)  // second layer behind the first: d x / d u < 0 at u = 0
    for (int j = 0; j < V.n(1); ++j)
      V.coefs[size_t(V.index(1, j, k))] = V.coefs[size_t(V.index(0, j, k))] - 2.0 * a.A.col(0);
  FitQuality qb = evaluate_fit(a.L, S, bad, fs, a.Y, a.t.B, diag);
  CHECK(qb.nonpositive_jacobian > 0);
  CHECK(qb.blocks_with_nonpositive == 1);
}

TEST_CASE("fit: U shape with smoothing, finer elements, C0 kept") {
  Affine a = make_affine(kU);
  SplineSpace S = build_spline_space(a.L, 3, 0.5);
  FitSamples fs = sample_blocks(a.L, S, 3, a.t.X, a.t.T, a.Y);
  const double diag = (a.Y.colwise().maxCoeff() - a.Y.colwise().minCoeff()).norm();
  FitSolution sol = fit_volumes(a.L, S, fs, 1e-6);
  FitQuality q = evaluate_fit(a.L, S, sol, fs, a.Y, a.t.B, diag);
  CHECK(q.max_error <= 1e-4);
  CHECK(q.c0_max_gap <= 1e-12);
  CHECK(q.nonpositive_jacobian == 0);
  // Without smoothing and with a single sample per element the system is still solvable here,
  // but a degree too high for the samples makes it singular: reported as an error.
  FitSamples sparse = sample_blocks(a.L, build_spline_space(a.L, 5, 0.5), 1, a.t.X, a.t.T, a.Y);
  CHECK_THROWS_AS(fit_volumes(a.L, build_spline_space(a.L, 5, 0.5), sparse, 0.0), Error);
}

TEST_CASE("fit params: json round trip and validation") {
  FitParams p;
  p.update_from_json({{"degree", 2}, {"element_size", 0.1}});
  CHECK(p.degree == 2);
  FitParams q;
  q.update_from_json(p.to_json());
  CHECK(q.to_json() == p.to_json());
  CHECK_THROWS_AS(FitParams{}.update_from_json({{"degree", 0}}), Error);
  CHECK_THROWS_AS(FitParams{}.update_from_json({{"element_size", 0.0}}), Error);
  CHECK_THROWS_AS(FitParams{}.update_from_json({{"samples_per_element", 0}}), Error);
  CHECK_THROWS_AS(FitParams{}.update_from_json({{"smoothing", -1.0}}), Error);
  CHECK_THROWS_AS(FitParams{}.update_from_json({{"order", 3}}), Error);
}

TEST_CASE("G+Smo XML: round trip, interfaces and boundary sides") {
  Affine a = make_affine(kL);
  SplineSpace S = build_spline_space(a.L, 3, 0.5);
  FitSamples fs = sample_blocks(a.L, S, 2, a.t.X, a.t.T, a.Y);
  FitSolution sol = fit_volumes(a.L, S, fs, 1e-6);
  GismoMultiPatch mp = gismo_topology(a.L, sol.volumes);
  CHECK(int(mp.interfaces.size()) == a.L.interior_faces());
  CHECK(int(mp.boundary.size()) == a.L.boundary_faces());
  // Blocks 0 and 1 are x-neighbors: east of 0 glued to west of 1, normal flag 0.
  CHECK(mp.interfaces[0] == std::array<int, 10>{0, 2, 1, 1, 0, 1, 2, 0, 1, 1});
  const auto path = testutil::temp_dir() / "lshape_bspline.xml";
  write_gismo_xml(path, mp);
  GismoMultiPatch r = read_gismo_xml(path);
  REQUIRE(r.volumes.size() == mp.volumes.size());
  for (size_t g = 0; g < r.volumes.size(); ++g) {
    for (int d = 0; d < 3; ++d) {
      CHECK(r.volumes[g].basis[size_t(d)].degree == 3);
      CHECK(r.volumes[g].basis[size_t(d)].knots == mp.volumes[g].basis[size_t(d)].knots);
    }
    CHECK(r.volumes[g].coefs == mp.volumes[g].coefs);
  }
  CHECK(r.interfaces == mp.interfaces);
  CHECK(r.boundary == mp.boundary);
}

TEST_CASE("sampled volume mesh: hexahedra per block, fields written") {
  Affine a = make_affine(kL);
  SplineSpace S = build_spline_space(a.L, 3, 0.5);  // 2 elements per interval -> 9 points per direction
  FitSolution sol = fit_volumes(a.L, S, sample_blocks(a.L, S, 2, a.t.X, a.t.T, a.Y), 1e-6);
  SampleMesh M = sample_mesh(a.L, S, sol);
  CHECK(M.H.rows() == 3 * 8 * 8 * 8);
  CHECK(M.X.rows() == 3 * 9 * 9 * 9);
  CHECK(M.block.maxCoeff() == 2);
  CHECK(M.scaled_jacobian.minCoeff() > 0.0);
  // Positive hex orientation: the first edge triple of every cell spans a positive volume.
  for (Index c = 0; c < M.H.rows(); c += 97) {
    const Eigen::Vector3d o = M.X.row(M.H(c, 0)), e1 = M.X.row(M.H(c, 1)), e2 = M.X.row(M.H(c, 3)), e3 = M.X.row(M.H(c, 4));
    CHECK((e1 - o).cross(e2 - o).dot(e3 - o) > 0.0);
  }
  const auto path = testutil::temp_dir() / "lshape_samples.vtu";
  write_vtu_hexes(path, M.X, M.H, {CellField::integer("block", M.block)},
                  {PointField::scalar("scaled_jacobian", M.scaled_jacobian)});
  const std::string xml = text::read_file(path);
  CHECK(xml.find("NumberOfCells=\"1536\"") != std::string::npos);
  CHECK(xml.find("<PointData>") != std::string::npos);
  CHECK(xml.find("Name=\"scaled_jacobian\"") != std::string::npos);
  CHECK(xml.find("Name=\"block\"") != std::string::npos);
  CHECK_THROWS_AS(write_vtu_hexes(path, M.X, M.H, {}, {PointField::scalar("bad", Eigen::VectorXd::Zero(3))}), Error);
}
