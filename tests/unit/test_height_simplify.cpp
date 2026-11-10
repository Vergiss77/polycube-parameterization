// Height simplification (Chen 2019 section 5): mesh-independent graph tests and a small polycube.
#include <doctest/doctest.h>

#include <cmath>

#include "helpers.h"
#include "pcp/geom/boundary_topology.h"
#include "pcp/structure/height_simplify.h"
#include "pcp/structure/labeling.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

namespace {

HeightGraph make_graph(const std::vector<std::pair<int, double>>& groups,
                       const std::vector<std::pair<int, int>>& edges) {
  HeightGraph g;
  for (size_t i = 0; i < groups.size(); ++i) {
    HeightGroup grp;
    grp.axis = groups[i].first;
    grp.height = groups[i].second;
    grp.patches = {static_cast<int>(i)};
    g.groups.push_back(grp);
  }
  g.n_patches = static_cast<int>(groups.size());
  for (auto [a, b] : edges) g.edges.push_back({a, b, groups[size_t(b)].second - groups[size_t(a)].second});
  return g;
}

void check_ratios(const HeightGraph& g, const SimplifyResult& r, double s) {
  for (const HeightEdge& e : g.edges) {
    const double q = (r.group_height[size_t(e.to)] - r.group_height[size_t(e.from)]) / e.rest;
    CHECK(q >= 1.0 / s - 1e-12);
    CHECK(q <= s + 1e-12);
  }
  CHECK(r.min_ratio >= 1.0 / s - 1e-12);
  CHECK(r.max_ratio <= s + 1e-12);
}

// Base of three voxels along X with a voxel on each end (a "U"); the left top is raised by
// `raise` so that the two tops have different heights. 10 distinct heights before simplification:
// X 0,1,2,3 / Y 0,1 / Z 0,1,2,2+raise.
SurfaceMesh make_u_shape(int n, double raise) {
  SurfaceMesh s = testutil::make_voxels({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 0, 1}, {2, 0, 1}}, n);
  for (Index i = 0; i < s.V.rows(); ++i)
    if (s.V(i, 2) == 2.0 && s.V(i, 0) <= 1.0) s.V(i, 2) += raise;
  return s;
}

}  // namespace

TEST_CASE("simplify params: json round trip and validation") {
  SimplifyParams p;
  p.enabled = false;
  p.max_stretch = 3.5;
  SimplifyParams q;
  q.update_from_json(p.to_json());
  CHECK(q.enabled == false);
  CHECK(q.max_stretch == 3.5);
  CHECK_THROWS_AS(q.update_from_json({{"bogus", 1}}), Error);
  CHECK_THROWS_AS(q.update_from_json({{"max_stretch", 0.5}}), Error);
}

TEST_CASE("merge cost: optimal merged height matches the closed form") {
  // L (h=0) -> A (h=1, rest 1), L -> B (h=4, rest 4). Merging A and B at h:
  // f(h) = 1/h + h + 4/h + h/4, minimum at h = sqrt(1 * 4) = 2.
  HeightGraph g = make_graph({{0, 0.0}, {0, 1.0}, {0, 4.0}}, {{0, 1}, {0, 2}});
  HeightState st = HeightState::initial(g);
  PairCost pc = merge_cost(g, st, 1, 2, 3.0);
  REQUIRE(pc.feasible);
  CHECK(std::abs(pc.height - 2.0) <= 1e-9 * 4.0);
  CHECK(pc.cost == doctest::Approx((0.5 + 2.0 + 2.0 + 0.5) - 4.0).epsilon(1e-12));
  // With s = 1.5 the intervals [4/1.5, 1.5] ... are disjoint: infeasible.
  CHECK_FALSE(merge_cost(g, st, 1, 2, 1.5).feasible);
  // An edge inside the pair means the groups are ordered.
  CHECK_FALSE(merge_cost(g, st, 0, 1, 10.0).feasible);
}

TEST_CASE("simplify: ordered groups are never merged") {
  HeightGraph g = make_graph({{0, 0.0}, {0, 1.0}, {0, 2.0}}, {{0, 1}, {1, 2}});
  SimplifyResult r = simplify_heights(g, 10.0);
  CHECK(r.merges.empty());
  CHECK(r.cuts_total_after() == 3);
  CHECK(r.group_height[0] < r.group_height[1]);
  CHECK(r.group_height[1] < r.group_height[2]);
}

TEST_CASE("simplify: two parallel chains collapse into one, axes stay separate") {
  // P0 (0) -> P1 (1) and Q0 (0.1) -> Q1 (1.1) on X; one isolated group on Y at height 0.
  HeightGraph g = make_graph({{0, 0.0}, {0, 1.0}, {0, 0.1}, {0, 1.1}, {1, 0.0}}, {{0, 1}, {2, 3}});
  SimplifyResult r = simplify_heights(g, 2.0);
  CHECK(r.cuts_before == std::array<int, 3>{4, 1, 0});
  CHECK(r.cuts_after == std::array<int, 3>{2, 1, 0});
  CHECK(r.merges.size() == 2);
  CHECK(r.group_rep[2] == r.group_rep[0]);
  CHECK(r.group_rep[3] == r.group_rep[1]);
  CHECK(r.group_rep[4] == 4);
  check_ratios(g, r, 2.0);
}

TEST_CASE("simplify: stretch bound respected; larger bound merges at least as much") {
  // Several short chains of different lengths and offsets on one axis.
  std::vector<std::pair<int, double>> groups;
  std::vector<std::pair<int, int>> edges;
  for (int k = 0; k < 6; ++k) {
    const double base = 0.13 * k, len = 1.0 + 0.37 * k;
    groups.push_back({2, base});
    groups.push_back({2, base + len});
    edges.push_back({2 * k, 2 * k + 1});
  }
  HeightGraph g = make_graph(groups, edges);
  CHECK(simplify_heights(g, 1.0).merges.empty());
  int prev = 1 << 30;
  for (double s : {1.05, 1.3, 2.0, 4.0}) {
    SimplifyResult r = simplify_heights(g, s);
    check_ratios(g, r, s);
    CHECK(r.cuts_total_after() <= prev);
    prev = r.cuts_total_after();
  }
  CHECK(prev == 2);
  // Determinism.
  SimplifyResult a = simplify_heights(g, 1.3), b = simplify_heights(g, 1.3);
  REQUIRE(a.merges.size() == b.merges.size());
  for (size_t i = 0; i < a.merges.size(); ++i) {
    CHECK(a.merges[i].a == b.merges[i].a);
    CHECK(a.merges[i].b == b.merges[i].b);
    CHECK(a.merges[i].height == b.merges[i].height);
  }
  CHECK(a.patch_height == b.patch_height);
}

TEST_CASE("height graph of a U-shaped polycube: groups, edges and one feasible merge") {
  SurfaceMesh s = make_u_shape(2, 0.2);
  BoundaryTopology topo = build_boundary_topology(s.F);
  PatchGraph pg = build_patch_graph(s.V, s.F, topo, initial_labels(s.V, s.F));
  REQUIRE(check_polycube_topology(s.F, pg).ok);
  // Patches: X 0,1,2,3 (4), Y front/back (2), Z bottom, base top, two tops (4).
  CHECK(pg.patches.label.size() == 10);
  HeightGraph hg = build_height_graph(pg, pg.height);
  CHECK(hg.groups.size() == 10);
  CHECK(hg.edges.size() == pg.chains.size());
  CHECK(hg.chains_without_ends == 0);

  // Chain edges alone leave the slab (bottom z=0 vs base top z=1) unordered: it collapses.
  auto group_at = [&](int label, double h) {
    for (size_t i = 0; i < hg.groups.size(); ++i)
      if (pg.patches.label[size_t(hg.groups[i].patches[0])] == label && hg.groups[i].height == h)
        return static_cast<int>(i);
    return -1;
  };
  const int bottom = group_at(make_label(2, -1), 0.0), base_top = group_at(make_label(2, 1), 1.0);
  const int inner_left = group_at(make_label(0, 1), 1.0), inner_right = group_at(make_label(0, -1), 2.0);
  REQUIRE(bottom >= 0);
  REQUIRE(base_top >= 0);
  REQUIRE(inner_left >= 0);
  REQUIRE(inner_right >= 0);
  {
    SimplifyResult r0 = simplify_heights(hg, 2.0);
    CHECK(r0.group_rep[size_t(bottom)] == r0.group_rep[size_t(base_top)]);
  }
  // Separation edges: across the slab (material) and across the gap between the two posts.
  const int added = add_separation_edges(hg, pg, pg.height, s.V, s.F);
  CHECK(added == hg.separation_edges);
  bool slab = false, gap = false;
  for (const HeightEdge& e : hg.edges) {
    CHECK(hg.groups[size_t(e.from)].axis == hg.groups[size_t(e.to)].axis);
    CHECK(e.rest > 0.0);
    slab |= e.separation && e.from == bottom && e.to == base_top;
    gap |= e.separation && e.from == inner_left && e.to == inner_right;
  }
  CHECK(slab);
  CHECK(gap);
  SimplifyResult r = simplify_heights(hg, 2.0);
  CHECK(r.group_rep[size_t(bottom)] != r.group_rep[size_t(base_top)]);
  CHECK(r.group_rep[size_t(inner_left)] != r.group_rep[size_t(inner_right)]);
  CHECK(r.cuts_before == std::array<int, 3>{4, 2, 4});
  CHECK(r.cuts_after == std::array<int, 3>{4, 2, 3});  // only the two tops merge
  REQUIRE(r.merges.size() == 1);
  CHECK(r.merges[0].height > 2.0);
  CHECK(r.merges[0].height < 2.2);
  check_ratios(hg, r, 2.0);

  // Unequal heights within round-off form one group.
  std::vector<double> h = pg.height;
  h[0] = std::nextafter(h[0], 1e9);
  CHECK(build_height_graph(pg, h).groups.size() == 10);
}

TEST_CASE("simplify and remap: U shape tops are snapped, structure kept, exact planarity") {
  SurfaceMesh s = make_u_shape(2, 0.2);
  TetParams tp;
  tp.size_factor = 1.0;
  TetMesh t = tetrahedralize(s, tp);
  BoundaryTopology topo = build_boundary_topology(t.B);
  Eigen::VectorXi labels = initial_labels(t.X, t.B);
  PatchGraph g = build_patch_graph(t.X, t.B, topo, labels);
  FinalMapResult fm = compute_final_map(t.X, t.X, t.T, t.B, g, g.height, FinalMapParams{});
  REQUIRE(fm.inverted_tets == 0);

  SimplifiedPolycube sp = simplify_and_remap(t.X, fm.X, t.T, t.B, g, fm.heights, SimplifyParams{},
                                             FinalMapParams{});
  const double diag = (t.X.colwise().maxCoeff() - t.X.colwise().minCoeff()).norm();
  CHECK(sp.simplify.cuts_total_before() == 10);
  CHECK(sp.simplify.cuts_total_after() == 9);
  CHECK(sp.topology.ok);
  CHECK(sp.map.inverted_tets == 0);
  CHECK(sp.map.max_planarity_error <= 1e-9 * diag);
  for (size_t p = 0; p < g.patches.label.size(); ++p)
    CHECK(sp.map.heights[p] == doctest::Approx(sp.simplify.patch_height[p]).epsilon(1e-12));
  // The two tops now share one height.
  std::vector<double> tops;
  for (size_t p = 0; p < g.patches.label.size(); ++p)
    if (g.patches.label[p] == make_label(2, 1) && g.height[p] > 1.5) tops.push_back(sp.map.heights[p]);
  REQUIRE(tops.size() == 2);
  CHECK(tops[0] == tops[1]);
}

TEST_CASE("width edges: the strip between an inset block and the base end is kept") {
  // Base 4 x 3 x 1 voxels with a 2 x 1 x 1 block on top, inset in Y. The block's +X face
  // (x = 2) and the base's +X end (x = 4) are only separated by a strip of the base top.
  std::vector<std::array<int, 3>> vox;
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 3; ++y) vox.push_back({x, y, 0});
  vox.push_back({0, 1, 1});
  vox.push_back({1, 1, 1});
  SurfaceMesh s = testutil::make_voxels(vox, 1);
  BoundaryTopology topo = build_boundary_topology(s.F);
  PatchGraph pg = build_patch_graph(s.V, s.F, topo, initial_labels(s.V, s.F));
  REQUIRE(check_polycube_topology(s.F, pg).ok);
  HeightGraph hg = build_height_graph(pg, pg.height);
  add_separation_edges(hg, pg, pg.height, s.V, s.F);
  int block_px = -1, base_px = -1;
  for (size_t i = 0; i < hg.groups.size(); ++i) {
    if (pg.patches.label[size_t(hg.groups[i].patches[0])] != make_label(0, 1)) continue;
    if (hg.groups[i].height == 2.0) block_px = int(i);
    if (hg.groups[i].height == 4.0) base_px = int(i);
  }
  REQUIRE(block_px >= 0);
  REQUIRE(base_px >= 0);
  {
    SimplifyResult r0 = simplify_heights(hg, 2.0);  // without width edges the strip collapses
    CHECK(r0.group_rep[size_t(block_px)] == r0.group_rep[size_t(base_px)]);
  }
  const int added = add_width_edges(hg, pg, pg.height, s.V, s.F);
  CHECK(added == hg.width_edges);
  bool strip = false;
  for (const HeightEdge& e : hg.edges) {
    CHECK(e.rest > 0.0);
    strip |= e.width && e.from == block_px && e.to == base_px;
  }
  CHECK(strip);
  SimplifyResult r = simplify_heights(hg, 2.0);
  CHECK(r.group_rep[size_t(block_px)] != r.group_rep[size_t(base_px)]);
  CHECK(r.width_edges == hg.width_edges);
}

TEST_CASE("rollback keeps only the merges that do not invert tets; untangle budget") {
  SurfaceMesh s = make_u_shape(2, 0.2);
  TetParams tp;
  tp.size_factor = 1.0;
  TetMesh t = tetrahedralize(s, tp);
  BoundaryTopology topo = build_boundary_topology(t.B);
  PatchGraph g = build_patch_graph(t.X, t.B, topo, initial_labels(t.X, t.B));
  FinalMapResult fm = compute_final_map(t.X, t.X, t.T, t.B, g, g.height, FinalMapParams{});
  REQUIRE(fm.inverted_tets == 0);

  // Chain edges only: the greedy search merges the two tops, then collapses the base slab.
  const HeightGraph chains_only = build_height_graph(g, fm.heights);
  const SimplifyResult greedy = simplify_heights(chains_only, 2.0);
  REQUIRE(greedy.merges.size() == 2);
  SimplifiedPolycube sp =
      remap_with_rollback(t.X, fm.X, t.T, t.B, g, chains_only, SimplifyParams{}, FinalMapParams{});
  CHECK(sp.simplify.merges_proposed == 2);
  CHECK(sp.simplify.merges.size() == 1);
  CHECK(sp.inverted_reference == 0);
  CHECK(sp.inverted_kept == 0);
  CHECK(sp.map.inverted_tets == 0);
  CHECK(sp.simplify.cuts_total_after() == 9);

  // The collapsed slab inverts tets; with a zero budget the repair stops immediately and the
  // remaining inversions are reported as they are.
  const SimplifyResult collapsed = apply_merges(chains_only, greedy.merges, 2);
  FinalMapParams fp;
  fp.heights_fixed = true;
  fp.untangle_budget = 0.0;
  FinalMapResult r0 = compute_final_map(t.X, fm.X, t.T, t.B, g, collapsed.patch_height, fp);
  CHECK(r0.inverted_before_untangle > 0);
  CHECK(r0.inverted_tets == r0.inverted_before_untangle);
  fp.untangle_budget = FinalMapParams{}.untangle_budget;
  FinalMapResult r1 = compute_final_map(t.X, fm.X, t.T, t.B, g, collapsed.patch_height, fp);
  CHECK(r1.inverted_tets <= r1.inverted_before_untangle);
}
